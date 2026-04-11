#include "DitherAnalyser.h"
#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QDebug>

#include <fitsio.h>
#ifdef HAVE_STELLARSOLVER
#include <stellarsolver.h>
#include <parameters.h>
#endif

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

#ifdef HAVE_STELLARSOLVER
// ---------------------------------------------------------------------------
// CLI plate-solve verification (no GUI)
// ---------------------------------------------------------------------------

static QStringList findIndexFiles()
{
    QStringList indexPaths;
    QStringList searchPaths = {
        "/usr/share/astrometry",
        "/usr/local/share/astrometry",
        "/usr/local/astrometry/data",
        "/opt/homebrew/share/astrometry"
    };
    for (const QString &path : searchPaths) {
        QDir indexDir(path);
        if (indexDir.exists()) {
            QFileInfoList indexFiles = indexDir.entryInfoList({"index-*.fits"}, QDir::Files);
            if (!indexFiles.isEmpty())
                indexPaths.append(path);
        }
    }
    return indexPaths;
}

struct CLIFrameResult {
    QString filename;
    double headerRA, headerDec;
    bool solved;
    double solvedRA, solvedDec, pixscale;
    double raErr, decErr, totalErr;
};

static bool loadFITSForSolver(StellarSolver *solver, const QString &fitsFile,
                               double hintRA, double hintDec,
                               std::vector<uint8_t> &bufferOut)
{
    fitsfile *fptr = nullptr;
    int status = 0;
    int naxis, bitpix;
    long naxes[2];

    if (fits_open_file(&fptr, fitsFile.toLocal8Bit().data(), READONLY, &status)) {
        std::cerr << "  Cannot open FITS file (cfitsio status " << status << ")\n";
        return false;
    }

    solver->setSearchPositionInDegrees(hintRA, hintDec);

    status = 0;
    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) || naxis != 2) {
        std::cerr << "  Not a 2D image (naxis=" << naxis << ", status=" << status << ")\n";
        fits_close_file(fptr, &status);
        return false;
    }

    int width = static_cast<int>(naxes[0]);
    int height = static_cast<int>(naxes[1]);
    long npixels = width * height;

    std::cerr << "  Image: " << width << "x" << height
              << " bitpix=" << bitpix << " pixels=" << npixels << "\n";

    std::vector<float> imageData(npixels);
    long firstPix[2] = {1, 1};
    if (fits_read_pix(fptr, TFLOAT, firstPix, npixels, nullptr,
                     imageData.data(), nullptr, &status)) {
        std::cerr << "  Cannot read pixel data (status " << status << ")\n";
        fits_close_file(fptr, &status);
        return false;
    }

    fits_close_file(fptr, &status);

    float minVal = *std::min_element(imageData.begin(), imageData.end());
    float maxVal = *std::max_element(imageData.begin(), imageData.end());
    std::cerr << "  Pixel range: " << minVal << " .. " << maxVal << "\n";

    const int downsampleFactor = 2;
    int outputWidth = width / downsampleFactor;
    int outputHeight = height / downsampleFactor;

    bufferOut.resize(outputWidth * outputHeight);
    float range = maxVal - minVal;

    if (range > 0) {
        for (int y = 0; y < outputHeight; y++) {
            for (int x = 0; x < outputWidth; x++) {
                int srcY = y * downsampleFactor;
                int srcX = x * downsampleFactor;
                if (srcY + 1 >= height || srcX + 1 >= width) continue;

                float r  = imageData[srcY * width + srcX];
                float g1 = imageData[srcY * width + srcX + 1];
                float g2 = imageData[(srcY + 1) * width + srcX];
                float b  = imageData[(srcY + 1) * width + srcX + 1];

                float avg = (r + (g1 + g2) * 0.5f + b) / 3.0f;
                float normalized = (avg - minVal) / range;
                bufferOut[y * outputWidth + x] =
                    static_cast<uint8_t>(std::clamp(normalized * 255.0f, 0.0f, 255.0f));
            }
        }
    } else {
        std::fill(bufferOut.begin(), bufferOut.end(), 128);
    }

    FITSImage::Statistic stats{};
    stats.width = outputWidth;
    stats.height = outputHeight;
    stats.channels = 1;
    stats.dataType = TBYTE;
    stats.bytesPerPixel = 1;
    for (int i = 0; i < 3; i++) {
        stats.min[i]    = (i == 0) ? minVal : 0.0;
        stats.max[i]    = (i == 0) ? maxVal : 0.0;
        stats.mean[i]   = (i == 0) ? (minVal + maxVal) / 2.0 : 0.0;
        stats.stddev[i] = 0.0;
        stats.median[i] = (i == 0) ? (minVal + maxVal) / 2.0 : 0.0;
    }
    stats.SNR = 1.0;

    std::cerr << "  Downsampled to " << outputWidth << "x" << outputHeight << "\n";

    if (!solver->loadNewImageBuffer(stats, bufferOut.data())) {
        std::cerr << "  loadNewImageBuffer failed\n";
        return false;
    }

    std::cerr << "  Image loaded into solver OK\n";
    return true;
}

static int runCLIVerify(const QString &directory, double searchRadius, bool singleFile)
{
    QCoreApplication *app = QCoreApplication::instance();

    QStringList indexPaths = findIndexFiles();
    if (indexPaths.isEmpty()) {
        std::cerr << "ERROR: No astrometry index files found\n";
        return 1;
    }
    std::cerr << "Index paths:";
    for (const auto &p : indexPaths)
        std::cerr << " " << p.toStdString();
    std::cerr << "\n";

    // Gather FITS files
    QStringList files;
    if (singleFile) {
        QFileInfo fi(directory);
        if (!fi.exists()) {
            std::cerr << "ERROR: File not found: " << directory.toStdString() << "\n";
            return 1;
        }
        files << fi.absoluteFilePath();
    } else {
        QDir dir(directory);
        QStringList filters = {"*.fits", "*.fit", "*.fts", "*.FITS", "*.FIT"};
        for (const auto &f : dir.entryInfoList(filters, QDir::Files, QDir::Name))
            files << f.absoluteFilePath();
    }

    if (files.isEmpty()) {
        std::cerr << "ERROR: No FITS files found\n";
        return 1;
    }

    // Setup solver parameters
    QList<Parameters> profiles = StellarSolver::getBuiltInProfiles();
    if (profiles.isEmpty()) {
        std::cerr << "ERROR: No StellarSolver profiles\n";
        return 1;
    }

    Parameters params = profiles.at(0);
    params.multiAlgorithm = SSolver::MULTI_AUTO;
    params.search_radius = searchRadius;
    params.minwidth = 0.1;
    params.maxwidth = 10.0;
    params.resort = true;
    params.autoDownsample = false;
    params.downsample = 1;
    params.inParallel = true;
    params.solverTimeLimit = 120;
    params.initialKeep = 2000;
    params.keepNum = 500;
    params.r_min = 1.0;
    params.removeBrightest = 0;
    params.removeDimmest = 50;
    params.saturationLimit = 65000;
    params.minarea = 5;
    params.threshold_bg_multiple = 2.0;

    std::cerr << "Search radius: " << searchRadius << " deg\n";
    std::cerr << "Files to verify: " << files.size() << "\n\n";

    // Print CSV header
    std::cout << "File,HeaderRA,HeaderDec,Temp_C,FocusPos,Solved,SolvedRA,SolvedDec,Pixscale,RAErr_arcsec,DecErr_arcsec,TotalErr_arcsec\n";

    int solved = 0, failed = 0, skipped = 0;
    std::vector<double> errors;

    for (int i = 0; i < files.size(); ++i) {
        QString filePath = files[i];
        QString baseName = QFileInfo(filePath).fileName();

        std::cerr << "[" << (i+1) << "/" << files.size() << "] " << baseName.toStdString() << "\n";

        // Read header RA/Dec
        fitsfile *fptr = nullptr;
        int status = 0;
        double headerRA = std::nan(""), headerDec = std::nan("");
        double temperature = std::nan(""), focusPos = std::nan("");

        if (fits_open_file(&fptr, filePath.toLocal8Bit().data(), READONLY, &status)) {
            std::cerr << "  Cannot open file, skipping\n";
            skipped++;
            continue;
        }

        auto tryKey = [&](const char *key, double &val) {
            int s = 0;
            fits_read_key(fptr, TDOUBLE, key, &val, nullptr, &s);
        };

        tryKey("CRVAL1", headerRA);
        tryKey("CRVAL2", headerDec);
        if (std::isnan(headerRA)) tryKey("RA", headerRA);
        if (std::isnan(headerDec)) tryKey("DEC", headerDec);

        tryKey("CCD-TEMP", temperature);
        if (std::isnan(temperature)) tryKey("SENSOR-T", temperature);
        if (std::isnan(temperature)) tryKey("TEMPERAT", temperature);

        tryKey("FOCUSPOS", focusPos);
        if (std::isnan(focusPos)) tryKey("FOCUS-PO", focusPos);
        if (std::isnan(focusPos)) tryKey("FOCPOS", focusPos);

        fits_close_file(fptr, &status);

        auto csvVal = [](double v) -> std::string {
            return std::isnan(v) ? "" : std::to_string(v);
        };

        if (std::isnan(headerRA) || std::isnan(headerDec)) {
            std::cerr << "  No RA/Dec in header, skipping\n";
            std::cout << baseName.toStdString() << ",,,,SKIP,,,,,,\n";
            skipped++;
            continue;
        }

        std::cerr << "  Header RA=" << headerRA << " Dec=" << headerDec << "\n";

        // Create solver
        StellarSolver solver;
        solver.setProperty("ProcessType", SSolver::SOLVE);
        solver.setProperty("ExtractorType", SSolver::EXTRACTOR_INTERNAL);
        solver.setProperty("SolverType", SSolver::SOLVER_STELLARSOLVER);
        solver.setIndexFolderPaths(indexPaths);
        solver.setParameters(params);

        std::vector<uint8_t> buffer;
        if (!loadFITSForSolver(&solver, filePath, headerRA, headerDec, buffer)) {
            std::cerr << "  Failed to load image for solving\n";
            std::cout << baseName.toStdString() << ","
                      << headerRA << "," << headerDec << ","
                      << csvVal(temperature) << "," << csvVal(focusPos)
                      << ",LOAD_FAIL,,,,,,\n";
            failed++;
            continue;
        }

        // Solve synchronously
        std::cerr << "  Starting solver...\n";

        QEventLoop loop;
        QObject::connect(&solver, &StellarSolver::finished, &loop, &QEventLoop::quit);
        solver.start();
        loop.exec();

        if (solver.solvingDone() && solver.hasWCSData()) {
            FITSImage::Solution solution = solver.getSolution();
            double solvedRA = solution.ra;
            double solvedDec = solution.dec;
            double pixscale = solution.pixscale / 2.0;

            double cosDec = std::cos(headerDec * M_PI / 180.0);
            double raErr  = (headerRA - solvedRA) * cosDec * 3600.0;
            double decErr = (headerDec - solvedDec) * 3600.0;
            double totalErr = std::sqrt(raErr * raErr + decErr * decErr);

            errors.push_back(totalErr);
            solved++;

            std::cerr << "  SOLVED: RA=" << solvedRA << " Dec=" << solvedDec
                      << " scale=" << pixscale << "\"/px"
                      << " error=" << totalErr << "\"\n";

            std::cout << baseName.toStdString() << ","
                      << headerRA << "," << headerDec << ","
                      << csvVal(temperature) << "," << csvVal(focusPos) << ",OK,"
                      << solvedRA << "," << solvedDec << ","
                      << pixscale << ","
                      << raErr << "," << decErr << "," << totalErr << "\n";
        } else {
            failed++;
            std::cerr << "  FAILED to solve";
            if (solver.failed()) std::cerr << " (solver error)";
            std::cerr << "\n";

            std::cout << baseName.toStdString() << ","
                      << headerRA << "," << headerDec << ","
                      << csvVal(temperature) << "," << csvVal(focusPos)
                      << ",FAIL,,,,,,\n";
        }
    }

    // Summary to stderr
    std::cerr << "\n=== Summary ===\n";
    std::cerr << "Solved: " << solved << "  Failed: " << failed << "  Skipped: " << skipped << "\n";
    if (!errors.empty()) {
        double sum = 0, mx = 0;
        for (double e : errors) { sum += e; mx = std::max(mx, e); }
        std::sort(errors.begin(), errors.end());
        std::cerr << "Error (arcsec): mean=" << (sum / errors.size())
                  << "  median=" << errors[errors.size()/2]
                  << "  max=" << mx << "\n";
    }

    return (failed > 0 && solved == 0) ? 1 : 0;
}
#endif // HAVE_STELLARSOLVER

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
#ifdef HAVE_STELLARSOLVER
    // Quick check for --verify before creating QApplication
    bool wantCLI = false;
    for (int i = 1; i < argc; ++i) {
        if (QString(argv[i]) == "--verify") {
            wantCLI = true;
            break;
        }
    }

    if (wantCLI) {
        // CLI mode: use QCoreApplication (no display needed)
        QCoreApplication app(argc, argv);
        app.setApplicationName("DitherAnalyser");

        QCommandLineParser parser;
        parser.setApplicationDescription("Dither pattern analyser with plate-solve verification");
        parser.addHelpOption();

        QCommandLineOption verifyOpt("verify", "Run plate-solve verification on FITS files (CLI mode)");
        parser.addOption(verifyOpt);

        QCommandLineOption radiusOpt({"r", "radius"},
            "Search radius in degrees (default: 2.0)", "degrees", "2.0");
        parser.addOption(radiusOpt);

        parser.addPositionalArgument("path", "FITS file or directory of FITS files");
        parser.process(app);

        QStringList positional = parser.positionalArguments();
        if (positional.isEmpty()) {
            std::cerr << "ERROR: --verify requires a FITS file or directory path\n";
            parser.showHelp(1);
        }

        QString path = positional.first();
        double radius = parser.value(radiusOpt).toDouble();
        bool singleFile = QFileInfo(path).isFile();

        return runCLIVerify(path, radius, singleFile);
    }
#endif // HAVE_STELLARSOLVER

    // GUI mode
    QApplication app(argc, argv);

    // Use Fusion if available (Qt <6.11), otherwise fall back to default
    if (auto *style = QApplication::setStyle("Fusion"); !style)
        QApplication::setStyle("modernwindows");

    // Dark palette
    QPalette darkPalette;
    darkPalette.setColor(QPalette::Window, QColor(30, 30, 46));
    darkPalette.setColor(QPalette::WindowText, QColor(205, 214, 244));
    darkPalette.setColor(QPalette::Base, QColor(24, 24, 37));
    darkPalette.setColor(QPalette::AlternateBase, QColor(30, 30, 46));
    darkPalette.setColor(QPalette::ToolTipBase, QColor(49, 50, 68));
    darkPalette.setColor(QPalette::ToolTipText, QColor(205, 214, 244));
    darkPalette.setColor(QPalette::Text, QColor(205, 214, 244));
    darkPalette.setColor(QPalette::Button, QColor(49, 50, 68));
    darkPalette.setColor(QPalette::ButtonText, QColor(205, 214, 244));
    darkPalette.setColor(QPalette::Link, QColor(137, 180, 250));
    darkPalette.setColor(QPalette::Highlight, QColor(137, 180, 250));
    darkPalette.setColor(QPalette::HighlightedText, QColor(30, 30, 46));
    app.setPalette(darkPalette);

    DitherAnalyser w;
    w.show();

    return app.exec();
}
