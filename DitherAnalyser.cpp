#include "DitherAnalyser.h"

#include <QApplication>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPushButton>
#include <QSplitter>
#include <QHeaderView>
#include <QDir>
#include <QMessageBox>
#include <QChart>
#include <QScatterSeries>
#include <QLineSeries>
#include <QValueAxis>
#include <QBarSeries>
#include <QBarSet>
#include <QBarCategoryAxis>
#include <QToolTip>
#include <QGroupBox>
#include <QDateTime>
#include <QEventLoop>
#include <QtConcurrent/QtConcurrentMap>
#include <QFutureWatcher>

#include <fitsio.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <set>

// ---------------------------------------------------------------------------
// HF sharpness metric via DCT on a central crop
// ---------------------------------------------------------------------------

// 1D DCT-II (unnormalised) applied to row of length N
static void dct1d(const float *in, float *out, int N)
{
    for (int k = 0; k < N; ++k) {
        double sum = 0;
        for (int n = 0; n < N; ++n)
            sum += in[n] * std::cos(M_PI * k * (2*n + 1) / (2.0 * N));
        out[k] = static_cast<float>(sum);
    }
}

// Compute sharpness as ratio of high-frequency power to total power
// using a 2D DCT on a central 128x128 crop of the Bayer-debayered image.
// Returns 0..1 (higher = sharper).
static double computeSharpness(const std::vector<float> &imgData,
                               int width, int height)
{
    const int CROP = 128;
    if (width < CROP * 2 || height < CROP * 2)
        return std::nan("");

    // Extract centre crop, debayer 2x2 -> mono
    int cx = width / 2 - CROP;
    int cy = height / 2 - CROP;

    std::vector<float> crop(CROP * CROP);
    for (int y = 0; y < CROP; ++y) {
        for (int x = 0; x < CROP; ++x) {
            int sy = cy + y * 2;
            int sx = cx + x * 2;
            float r  = imgData[sy * width + sx];
            float g1 = imgData[sy * width + sx + 1];
            float g2 = imgData[(sy+1) * width + sx];
            float b  = imgData[(sy+1) * width + sx + 1];
            crop[y * CROP + x] = (r + (g1+g2)*0.5f + b) / 3.0f;
        }
    }

    // Normalise to 0..1
    float minV = *std::min_element(crop.begin(), crop.end());
    float maxV = *std::max_element(crop.begin(), crop.end());
    float range = maxV - minV;
    if (range < 1.0f)
        return std::nan("");
    for (float &v : crop)
        v = (v - minV) / range;

    // 2D DCT: apply 1D DCT to rows, then to columns
    std::vector<float> tmp(CROP * CROP);
    std::vector<float> dctCoeffs(CROP * CROP);

    // DCT on rows
    for (int y = 0; y < CROP; ++y)
        dct1d(&crop[y * CROP], &tmp[y * CROP], CROP);

    // DCT on columns (transpose, DCT rows, transpose back)
    std::vector<float> col(CROP), colOut(CROP);
    for (int x = 0; x < CROP; ++x) {
        for (int y = 0; y < CROP; ++y)
            col[y] = tmp[y * CROP + x];
        dct1d(col.data(), colOut.data(), CROP);
        for (int y = 0; y < CROP; ++y)
            dctCoeffs[y * CROP + x] = colOut[y];
    }

    // Compute power in high-frequency vs total
    // HF = coefficients where (u + v) > CROP/4
    double totalPower = 0, hfPower = 0;
    int cutoff = CROP / 4;
    for (int v = 0; v < CROP; ++v) {
        for (int u = 0; u < CROP; ++u) {
            if (u == 0 && v == 0) continue; // skip DC
            double p = (double)dctCoeffs[v * CROP + u] * dctCoeffs[v * CROP + u];
            totalPower += p;
            if (u + v > cutoff)
                hfPower += p;
        }
    }

    return (totalPower > 0) ? hfPower / totalPower : 0;
}

// Open a FITS file and compute the HF sharpness metric
static double computeSharpnessFromFile(const QString &filePath)
{
    fitsfile *fptr = nullptr;
    int status = 0;
    if (fits_open_file(&fptr, filePath.toLocal8Bit().data(), READONLY, &status))
        return std::nan("");

    int naxis, bitpix;
    long naxes[2];
    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) != 0 || naxis != 2) {
        fits_close_file(fptr, &status);
        return std::nan("");
    }

    int w = static_cast<int>(naxes[0]);
    int h = static_cast<int>(naxes[1]);
    long npix = static_cast<long>(w) * h;
    std::vector<float> pixels(npix);
    long fp[2] = {1, 1};
    status = 0;
    if (fits_read_pix(fptr, TFLOAT, fp, npix, nullptr,
                      pixels.data(), nullptr, &status) != 0) {
        fits_close_file(fptr, &status);
        return std::nan("");
    }
    fits_close_file(fptr, &status);

    return computeSharpness(pixels, w, h);
}

// ---------------------------------------------------------------------------
// FitsLoader – runs in a worker thread
// ---------------------------------------------------------------------------

FitsLoader::FitsLoader(const QString &directory) : m_directory(directory) {}

void FitsLoader::process()
{
    QDir dir(m_directory);
    QStringList filters;
    filters << "*.fits" << "*.fit" << "*.fts" << "*.FITS" << "*.FIT";
    QStringList files = dir.entryList(filters, QDir::Files, QDir::Name);

    if (files.isEmpty()) {
        emit error("No FITS files found in " + m_directory);
        return;
    }

    std::vector<FrameInfo> frames;
    frames.reserve(files.size());

    for (int i = 0; i < files.size(); ++i) {
        emit progress(i + 1, files.size());

        QString path = dir.absoluteFilePath(files[i]);
        fitsfile *fptr = nullptr;
        int status = 0;

        if (fits_open_file(&fptr, path.toStdString().c_str(), READONLY, &status)) {
            continue; // skip unreadable files
        }

        FrameInfo fi{};
        fi.filename = files[i].toStdString();

        char value[FLEN_VALUE];

        // Read WCS centre
        auto readDouble = [&](const char *key, double &out) {
            int s = 0;
            if (fits_read_key(fptr, TDOUBLE, key, &out, nullptr, &s) != 0)
                out = std::nan("");
        };
        auto readInt = [&](const char *key, int &out) {
            int s = 0;
            long v = 0;
            if (fits_read_key(fptr, TLONG, key, &v, nullptr, &s) == 0)
                out = static_cast<int>(v);
        };
        auto readString = [&](const char *key, std::string &out) {
            int s = 0;
            char buf[FLEN_VALUE];
            if (fits_read_key(fptr, TSTRING, key, buf, nullptr, &s) == 0)
                out = buf;
        };

        readDouble("CRVAL1", fi.ra);
        readDouble("CRVAL2", fi.dec);
        readDouble("CRPIX1", fi.crpix1);
        readDouble("CRPIX2", fi.crpix2);
        readDouble("CDELT1", fi.cdelt1);
        readDouble("CDELT2", fi.cdelt2);
        readDouble("CROTA2", fi.crota2);
        readDouble("EXPTIME", fi.exptime);
        readInt("NAXIS1", fi.naxis1);
        readInt("NAXIS2", fi.naxis2);
        readString("DATE-OBS", fi.dateObs);

        // Also try RA/DEC keywords as fallback
        if (std::isnan(fi.ra)) readDouble("RA", fi.ra);
        if (std::isnan(fi.dec)) readDouble("DEC", fi.dec);

        // Temperature: try common FITS keyword names
        readDouble("CCD-TEMP", fi.temperature);
        if (std::isnan(fi.temperature)) readDouble("SENSOR-T", fi.temperature);
        if (std::isnan(fi.temperature)) readDouble("TEMPERAT", fi.temperature);
        if (std::isnan(fi.temperature)) readDouble("AMBIENT",  fi.temperature);

        // Focus position
        readDouble("FOCUSPOS", fi.focusPos);
        if (std::isnan(fi.focusPos)) readDouble("FOCUS-PO", fi.focusPos);
        if (std::isnan(fi.focusPos)) readDouble("FOCPOS",   fi.focusPos);
        if (std::isnan(fi.focusPos)) readDouble("FOCUSER",  fi.focusPos);

        fits_close_file(fptr, &status);

        if (!std::isnan(fi.ra) && !std::isnan(fi.dec))
            frames.push_back(fi);
    }

    emit finished(frames);
}

#ifdef HAVE_STELLARSOLVER
// ---------------------------------------------------------------------------
// SolverVerifier – plate-solves frames to verify headers
// ---------------------------------------------------------------------------

SolverVerifier::SolverVerifier(const std::vector<FrameInfo> &frames,
                               const QString &directory,
                               double searchRadius)
    : m_frames(frames), m_directory(directory), m_searchRadius(searchRadius)
{
}

void SolverVerifier::abort()
{
    m_abort.store(true);
}

QStringList SolverVerifier::findIndexFiles()
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
            QStringList filters;
            filters << "index-*.fits";
            QFileInfoList indexFiles = indexDir.entryInfoList(filters, QDir::Files);
            if (!indexFiles.isEmpty()) {
                indexPaths.append(path);
            }
        }
    }

    return indexPaths;
}

bool SolverVerifier::loadFITSForSolver(StellarSolver *solver, const QString &fitsFile,
                                        double hintRA, double hintDec, double searchRadius)
{
    fitsfile *fptr = nullptr;
    int status = 0;
    int naxis, bitpix;
    long naxes[2];

    if (fits_open_file(&fptr, fitsFile.toLocal8Bit().data(), READONLY, &status))
        return false;

    // Set coordinate hints
    solver->setSearchPositionInDegrees(hintRA, hintDec);

    // Get image dimensions
    status = 0;
    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) || naxis != 2) {
        fits_close_file(fptr, &status);
        return false;
    }

    int width = static_cast<int>(naxes[0]);
    int height = static_cast<int>(naxes[1]);
    long npixels = width * height;

    // Read image data as float
    std::vector<float> imageData(npixels);
    long firstPix[2] = {1, 1};
    if (fits_read_pix(fptr, TFLOAT, firstPix, npixels, nullptr,
                     imageData.data(), nullptr, &status)) {
        fits_close_file(fptr, &status);
        return false;
    }

    fits_close_file(fptr, &status);

    // Calculate statistics
    float minVal = *std::min_element(imageData.begin(), imageData.end());
    float maxVal = *std::max_element(imageData.begin(), imageData.end());

    // Downsample 2x with Bayer-aware processing
    const int downsampleFactor = 2;
    int outputWidth = width / downsampleFactor;
    int outputHeight = height / downsampleFactor;

    // Store the buffer in a member we keep alive until solver finishes
    auto buffer = std::make_shared<std::vector<uint8_t>>(outputWidth * outputHeight);
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
                (*buffer)[y * outputWidth + x] =
                    static_cast<uint8_t>(std::clamp(normalized * 255.0f, 0.0f, 255.0f));
            }
        }
    } else {
        std::fill(buffer->begin(), buffer->end(), 128);
    }

    // Create statistics
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

    if (!solver->loadNewImageBuffer(stats, buffer->data()))
        return false;

    // Keep buffer alive by storing in solver's property
    solver->setProperty("imageBuffer", QVariant::fromValue(buffer));
    return true;
}

void SolverVerifier::process()
{
    QStringList indexPaths = findIndexFiles();
    if (indexPaths.isEmpty()) {
        emit error("No astrometry index files found");
        return;
    }

    // Setup common solver parameters
    QList<Parameters> profiles = StellarSolver::getBuiltInProfiles();
    if (profiles.isEmpty()) {
        emit error("No StellarSolver parameter profiles available");
        return;
    }

    Parameters params = profiles.at(0);
    params.multiAlgorithm = SSolver::MULTI_AUTO;
    params.search_radius = m_searchRadius;
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

    QDir dir(m_directory);

    for (int i = 0; i < (int)m_frames.size(); ++i) {
        if (m_abort.load())
            break;

        emit progress(i + 1, (int)m_frames.size());

        const FrameInfo &fi = m_frames[i];
        QString filePath = dir.absoluteFilePath(QString::fromStdString(fi.filename));

        StellarSolver solver;
        solver.setProperty("ProcessType", SSolver::SOLVE);
        solver.setProperty("ExtractorType", SSolver::EXTRACTOR_INTERNAL);
        solver.setProperty("SolverType", SSolver::SOLVER_STELLARSOLVER);
        solver.setIndexFolderPaths(indexPaths);
        solver.setParameters(params);

        if (!loadFITSForSolver(&solver, filePath, fi.ra, fi.dec, m_searchRadius)) {
            emit frameVerified(i, false, 0, 0, 0, 0, 0, 0);
            continue;
        }

        // Run solver synchronously using event loop
        QEventLoop loop;
        bool solverDone = false;
        connect(&solver, &StellarSolver::finished, &loop, [&]() {
            solverDone = true;
            loop.quit();
        });

        solver.start();
        loop.exec();

        if (solver.solvingDone() && solver.hasWCSData()) {
            FITSImage::Solution solution = solver.getSolution();
            double solvedRA = solution.ra;
            double solvedDec = solution.dec;
            double pixscale = solution.pixscale / 2.0; // correct for 2x downsample

            double cosDec = std::cos(fi.dec * M_PI / 180.0);
            double raErr  = (fi.ra - solvedRA) * cosDec * 3600.0;
            double decErr = (fi.dec - solvedDec) * 3600.0;
            double totalErr = std::sqrt(raErr * raErr + decErr * decErr);

            emit frameVerified(i, true, solvedRA, solvedDec, pixscale,
                               raErr, decErr, totalErr);
        } else {
            emit frameVerified(i, false, 0, 0, 0, 0, 0, 0);
        }
    }

    emit finished();
}

// ---------------------------------------------------------------------------
// FWHM extraction helper (called per frame from QtConcurrent)
// ---------------------------------------------------------------------------

static double extractFwhm(const QString &filePath, double cdelt1, double cdelt2)
{
    fitsfile *fptr = nullptr;
    int status = 0;
    if (fits_open_file(&fptr, filePath.toLocal8Bit().data(), READONLY, &status))
        return std::nan("");

    int naxis, bitpix;
    long naxes[2];
    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) != 0 || naxis != 2) {
        fits_close_file(fptr, &status);
        return std::nan("");
    }

    int width = static_cast<int>(naxes[0]);
    int height = static_cast<int>(naxes[1]);
    long npixels = static_cast<long>(width) * height;

    const int ds = 2;
    int outW = width / ds;
    int outH = height / ds;

    std::vector<float> imgData(npixels);
    long firstPix[2] = {1, 1};
    status = 0;
    if (fits_read_pix(fptr, TFLOAT, firstPix, npixels, nullptr,
                      imgData.data(), nullptr, &status) != 0) {
        fits_close_file(fptr, &status);
        return std::nan("");
    }
    fits_close_file(fptr, &status);

    float minVal = *std::min_element(imgData.begin(), imgData.end());
    float maxVal = *std::max_element(imgData.begin(), imgData.end());
    float range = maxVal - minVal;

    std::vector<uint8_t> buf(outW * outH);
    if (range > 0) {
        for (int y = 0; y < outH; y++) {
            for (int x = 0; x < outW; x++) {
                int sy = y * ds, sx = x * ds;
                if (sy + 1 >= height || sx + 1 >= width) continue;
                float r  = imgData[sy * width + sx];
                float g1 = imgData[sy * width + sx + 1];
                float g2 = imgData[(sy+1) * width + sx];
                float b  = imgData[(sy+1) * width + sx + 1];
                float avg = (r + (g1+g2)*0.5f + b) / 3.0f;
                float norm = (avg - minVal) / range;
                buf[y * outW + x] =
                    static_cast<uint8_t>(std::clamp(norm * 255.0f, 0.0f, 255.0f));
            }
        }
    } else {
        std::fill(buf.begin(), buf.end(), 128);
    }

    FITSImage::Statistic stats{};
    stats.width = outW;
    stats.height = outH;
    stats.channels = 1;
    stats.dataType = TBYTE;
    stats.bytesPerPixel = 1;
    for (int c = 0; c < 3; c++) {
        stats.min[c]    = (c == 0) ? minVal : 0.0;
        stats.max[c]    = (c == 0) ? maxVal : 0.0;
        stats.mean[c]   = (c == 0) ? (minVal + maxVal) / 2.0 : 0.0;
        stats.stddev[c] = 0.0;
        stats.median[c] = (c == 0) ? (minVal + maxVal) / 2.0 : 0.0;
    }
    stats.SNR = 1.0;

    StellarSolver solver;
    if (!solver.loadNewImageBuffer(stats, buf.data()))
        return std::nan("");

    // Synchronous extraction with HFR computation (uses all available star detection)
    solver.extract(true);

    QList<FITSImage::Star> stars = solver.getStarList();
    if (stars.isEmpty())
        return std::nan("");

    std::vector<double> fwhms;
    fwhms.reserve(stars.size());
    for (const auto &star : stars) {
        // Use HFR if available, else fall back to semi-major axis
        double fwhm = (star.HFR > 0) ? star.HFR * 2.0 : star.a * 2.35;
        if (fwhm > 0)
            fwhms.push_back(fwhm);
    }

    if (fwhms.empty())
        return std::nan("");

    std::sort(fwhms.begin(), fwhms.end());
    double medianFWHM = fwhms[fwhms.size() / 2];

    double pixscale = std::abs(cdelt1) * 3600.0;
    if (pixscale < 0.001)
        pixscale = std::abs(cdelt2) * 3600.0;

    return (pixscale > 0.001) ? medianFWHM * ds * pixscale : medianFWHM * ds;
}
#endif // HAVE_STELLARSOLVER

// ---------------------------------------------------------------------------
// DitherAnalyser – main window
// ---------------------------------------------------------------------------

DitherAnalyser::DitherAnalyser(QWidget *parent) : QMainWindow(parent)
{
    setWindowTitle("Dither Analyser");
    resize(1400, 900);
    buildUI();
}

DitherAnalyser::~DitherAnalyser()
{
    if (m_loaderThread) {
        m_loaderThread->quit();
        m_loaderThread->wait();
    }
#ifdef HAVE_STELLARSOLVER
    if (m_verifierThread) {
        if (m_verifier) m_verifier->abort();
        m_verifierThread->quit();
        m_verifierThread->wait();
    }
#endif
    if (m_sharpWatcher) {
        m_sharpWatcher->cancel();
        m_sharpWatcher->waitForFinished();
    }
#ifdef HAVE_STELLARSOLVER
    if (m_fwhmWatcher) {
        m_fwhmWatcher->cancel();
        m_fwhmWatcher->waitForFinished();
    }
#endif
}

void DitherAnalyser::buildUI()
{
    auto *central = new QWidget;
    auto *mainLayout = new QVBoxLayout(central);

    // Top bar: directory selector
    auto *topBar = new QHBoxLayout;
    auto *browseBtn = new QPushButton("Open Directory...");
    m_dirLabel = new QLabel("No directory selected");
    m_dirLabel->setStyleSheet("color: gray;");
    m_progress = new QProgressBar;
    m_progress->setVisible(false);
    m_progress->setTextVisible(true);
    topBar->addWidget(browseBtn);
    topBar->addWidget(m_dirLabel, 1);
    topBar->addWidget(m_progress);
    mainLayout->addLayout(topBar);

    connect(browseBtn, &QPushButton::clicked, this, &DitherAnalyser::browseDirectory);

    // Controls bar
    auto *verifyBar = new QHBoxLayout;

#ifdef HAVE_STELLARSOLVER
    m_verifyBtn = new QPushButton("Verify RA/Dec (Plate Solve)");
    m_verifyBtn->setEnabled(false);
    m_verifyBtn->setToolTip("Plate-solve each frame using StellarSolver to verify header RA/Dec");

    auto *radiusLabel = new QLabel("Search radius:");
    m_radiusSpin = new QDoubleSpinBox;
    m_radiusSpin->setRange(0.5, 30.0);
    m_radiusSpin->setValue(2.0);
    m_radiusSpin->setSuffix(" deg");
    m_radiusSpin->setDecimals(1);
    m_radiusSpin->setToolTip("Search radius around header RA/Dec hint");

    verifyBar->addWidget(m_verifyBtn);
    verifyBar->addWidget(radiusLabel);
    verifyBar->addWidget(m_radiusSpin);
    verifyBar->addSpacing(20);

    connect(m_verifyBtn, &QPushButton::clicked, this, &DitherAnalyser::startVerification);
#endif

    m_binBtn = new QPushButton("Bin by Quality");
    m_binBtn->setEnabled(false);
    m_binBtn->setToolTip("Sort frames into quality subfolders (sharpness, FWHM, or combined)");

    verifyBar->addWidget(m_binBtn);
    verifyBar->addStretch();
    mainLayout->addLayout(verifyBar);

    connect(m_binBtn, &QPushButton::clicked, this, &DitherAnalyser::binByQuality);

    // Summary label
    m_summaryLabel = new QLabel;
    m_summaryLabel->setWordWrap(true);
    m_summaryLabel->setStyleSheet("QLabel { background: #1e1e2e; color: #cdd6f4; padding: 10px; border-radius: 6px; font-family: monospace; font-size: 11pt; }");
    mainLayout->addWidget(m_summaryLabel);

#ifdef HAVE_STELLARSOLVER
    // Verification summary label
    m_verifySummaryLabel = new QLabel;
    m_verifySummaryLabel->setWordWrap(true);
    m_verifySummaryLabel->setStyleSheet("QLabel { background: #1e1e2e; color: #a6e3a1; padding: 10px; border-radius: 6px; font-family: monospace; font-size: 11pt; }");
    m_verifySummaryLabel->setVisible(false);
    mainLayout->addWidget(m_verifySummaryLabel);
#endif

    // Tabs for charts + table
    m_tabs = new QTabWidget;

    // Scatter plot
    m_scatterView = new QChartView;
    m_scatterView->setRenderHint(QPainter::Antialiasing);
    m_tabs->addTab(m_scatterView, "Dither Pattern");

    // Timeline
    m_timelineView = new QChartView;
    m_timelineView->setRenderHint(QPainter::Antialiasing);
    m_tabs->addTab(m_timelineView, "Offset Timeline");

    // Rose / direction histogram
    m_roseView = new QChartView;
    m_roseView->setRenderHint(QPainter::Antialiasing);
    m_tabs->addTab(m_roseView, "Direction Histogram");

    // Field rotation rate
    m_rotationView = new QChartView;
    m_rotationView->setRenderHint(QPainter::Antialiasing);
    m_tabs->addTab(m_rotationView, "Field Rotation");

    // Environment (temp, focus) — available immediately
    m_envView = new QChartView;
    m_envView->setRenderHint(QPainter::Antialiasing);
    m_tabs->addTab(m_envView, "Environment");

    // HF Sharpness — available immediately (computed during load)
    m_sharpView = new QChartView;
    m_sharpView->setRenderHint(QPainter::Antialiasing);
    m_tabs->addTab(m_sharpView, "Sharpness");

#ifdef HAVE_STELLARSOLVER
    // Star FWHM — deferred until tab is selected
    m_fwhmView = new QChartView;
    m_fwhmView->setRenderHint(QPainter::Antialiasing);
    m_tabs->addTab(m_fwhmView, "Star FWHM");
    m_envTabIndex = m_tabs->count() - 1;
#endif

    // Table
    m_table = new QTableWidget;
    m_table->setAlternatingRowColors(true);
    m_tabs->addTab(m_table, "Frame Table");

    connect(m_tabs, &QTabWidget::currentChanged, this, &DitherAnalyser::onTabChanged);

    mainLayout->addWidget(m_tabs, 1);
    setCentralWidget(central);
}

void DitherAnalyser::browseDirectory()
{
    QString dir = QFileDialog::getExistingDirectory(this, "Select FITS Directory",
        QDir::homePath() + "/Astrophotography");
    if (dir.isEmpty()) return;

    m_currentDirectory = dir;
    m_dirLabel->setText(dir);
    m_progress->setVisible(true);
    m_progress->setValue(0);
#ifdef HAVE_STELLARSOLVER
    m_verifyBtn->setEnabled(false);
    m_verifySummaryLabel->setVisible(false);
#endif
    m_binBtn->setEnabled(false);

    // Clean up previous loader
    if (m_loaderThread) {
        m_loaderThread->quit();
        m_loaderThread->wait();
        delete m_loaderThread;
    }

    m_loaderThread = new QThread;
    m_loader = new FitsLoader(dir);
    m_loader->moveToThread(m_loaderThread);

    connect(m_loaderThread, &QThread::started, m_loader, &FitsLoader::process);
    connect(m_loader, &FitsLoader::finished, this, &DitherAnalyser::onFramesLoaded);
    connect(m_loader, &FitsLoader::progress, this, &DitherAnalyser::onProgress);
    connect(m_loader, &FitsLoader::error, this, [this](const QString &msg) {
        QMessageBox::warning(this, "Error", msg);
        m_progress->setVisible(false);
    });
    connect(m_loader, &FitsLoader::finished, m_loaderThread, &QThread::quit);
    connect(m_loaderThread, &QThread::finished, m_loader, &QObject::deleteLater);

    m_loaderThread->start();
}

void DitherAnalyser::onProgress(int current, int total)
{
    m_progress->setMaximum(total);
    m_progress->setValue(current);
    m_progress->setFormat(QString("Loading %1 / %2").arg(current).arg(total));
}

void DitherAnalyser::onFramesLoaded(std::vector<FrameInfo> frames)
{
    m_progress->setVisible(false);
    m_frames = std::move(frames);

    if (m_frames.empty()) {
        QMessageBox::information(this, "Info", "No frames with valid WCS found.");
        return;
    }

#ifdef HAVE_STELLARSOLVER
    m_verifyBtn->setEnabled(true);
    m_fwhmDone = false;
#endif

    qDebug() << "onFramesLoaded: analysing" << m_frames.size() << "frames";
    analyse(m_frames);
    qDebug() << "onFramesLoaded: populateTable";
    populateTable(m_frames);
    qDebug() << "onFramesLoaded: plotDitherScatter";
    plotDitherScatter(m_frames);
    qDebug() << "onFramesLoaded: plotDitherTimeline";
    plotDitherTimeline(m_frames);
    qDebug() << "onFramesLoaded: plotDitherRose";
    plotDitherRose(m_frames);
    qDebug() << "onFramesLoaded: plotFieldRotation";
    plotFieldRotation(m_frames);
    qDebug() << "onFramesLoaded: plotEnvironment";
    plotEnvironment(m_frames);
    qDebug() << "onFramesLoaded: plotSharpness";
    plotSharpness(m_frames);
    qDebug() << "onFramesLoaded: updateSummary";
    updateSummary(m_frames);
    qDebug() << "onFramesLoaded: startSharpnessExtraction";

    // Kick off parallel sharpness extraction
    startSharpnessExtraction();
    qDebug() << "onFramesLoaded: done";
}

#ifdef HAVE_STELLARSOLVER
// ---------------------------------------------------------------------------
// Verification
// ---------------------------------------------------------------------------

void DitherAnalyser::startVerification()
{
    if (m_frames.empty() || m_currentDirectory.isEmpty()) return;

    // Clean up previous verifier
    if (m_verifierThread) {
        if (m_verifier) m_verifier->abort();
        m_verifierThread->quit();
        m_verifierThread->wait();
        delete m_verifierThread;
        m_verifierThread = nullptr;
    }

    // Reset verification state
    for (auto &f : m_frames) {
        f.verified = false;
        f.solveSuccess = false;
    }

    m_verifyBtn->setEnabled(false);
    m_verifyBtn->setText("Solving...");
    m_progress->setVisible(true);
    m_progress->setValue(0);
    m_verifySummaryLabel->setVisible(false);

    double radius = m_radiusSpin->value();

    m_verifierThread = new QThread;
    m_verifier = new SolverVerifier(m_frames, m_currentDirectory, radius);
    m_verifier->moveToThread(m_verifierThread);

    connect(m_verifierThread, &QThread::started, m_verifier, &SolverVerifier::process);
    connect(m_verifier, &SolverVerifier::progress, this, &DitherAnalyser::onVerifyProgress);
    connect(m_verifier, &SolverVerifier::frameVerified, this, &DitherAnalyser::onFrameVerified);
    connect(m_verifier, &SolverVerifier::finished, this, &DitherAnalyser::onVerifyFinished);
    connect(m_verifier, &SolverVerifier::error, this, [this](const QString &msg) {
        QMessageBox::warning(this, "Verification Error", msg);
        m_progress->setVisible(false);
        m_verifyBtn->setEnabled(true);
        m_verifyBtn->setText("Verify RA/Dec (Plate Solve)");
    });
    connect(m_verifier, &SolverVerifier::finished, m_verifierThread, &QThread::quit);
    connect(m_verifierThread, &QThread::finished, m_verifier, &QObject::deleteLater);

    m_verifierThread->start();
}

void DitherAnalyser::onVerifyProgress(int current, int total)
{
    m_progress->setMaximum(total);
    m_progress->setValue(current);
    m_progress->setFormat(QString("Plate solving %1 / %2").arg(current).arg(total));
}

void DitherAnalyser::onFrameVerified(int index, bool success, double solvedRA, double solvedDec,
                                      double pixscale, double raErr, double decErr, double totalErr)
{
    if (index < 0 || index >= (int)m_frames.size()) return;

    FrameInfo &fi = m_frames[index];
    fi.verified = true;
    fi.solveSuccess = success;

    if (success) {
        fi.solvedRA = solvedRA;
        fi.solvedDec = solvedDec;
        fi.solvedPixscale = pixscale;
        fi.raErrorArcsec = raErr;
        fi.decErrorArcsec = decErr;
        fi.totalErrorArcsec = totalErr;
    }

    // Update table row if verification columns exist
    if (m_table->columnCount() >= 16) {
        if (success) {
            m_table->setItem(index, 12, new QTableWidgetItem(QString::number(solvedRA, 'f', 6)));
            m_table->setItem(index, 13, new QTableWidgetItem(QString::number(solvedDec, 'f', 6)));
            m_table->setItem(index, 14, new QTableWidgetItem(QString::number(totalErr, 'f', 2)));

            auto *statusItem = new QTableWidgetItem("OK");
            if (totalErr > 60.0) {
                statusItem->setText("BAD");
                statusItem->setForeground(QColor(243, 139, 168));
            } else if (totalErr > 10.0) {
                statusItem->setText("WARN");
                statusItem->setForeground(QColor(249, 226, 175));
            } else {
                statusItem->setForeground(QColor(166, 227, 161));
            }
            m_table->setItem(index, 15, new QTableWidgetItem(*statusItem));
        } else {
            m_table->setItem(index, 12, new QTableWidgetItem("--"));
            m_table->setItem(index, 13, new QTableWidgetItem("--"));
            m_table->setItem(index, 14, new QTableWidgetItem("--"));
            auto *statusItem = new QTableWidgetItem("FAIL");
            statusItem->setForeground(QColor(243, 139, 168));
            m_table->setItem(index, 15, new QTableWidgetItem(*statusItem));
        }
    }
}

void DitherAnalyser::onVerifyFinished()
{
    m_progress->setVisible(false);
    m_verifyBtn->setEnabled(true);
    m_verifyBtn->setText("Verify RA/Dec (Plate Solve)");

    // Rebuild table with verification columns
    populateTable(m_frames);
    updateVerificationSummary();
}
#endif // HAVE_STELLARSOLVER

void DitherAnalyser::startSharpnessExtraction()
{
    if (m_frames.empty()) return;

    int n = (int)m_frames.size();
    m_sharpResults.assign(n, std::nan(""));

    QDir dir(m_currentDirectory);
    QStringList paths;
    for (int i = 0; i < n; ++i)
        paths << dir.absoluteFilePath(QString::fromStdString(m_frames[i].filename));

    auto indices = std::make_shared<QList<int>>();
    for (int i = 0; i < n; ++i)
        indices->append(i);

    double *results = m_sharpResults.data();
    auto *watcher = new QFutureWatcher<void>(this);
    m_sharpWatcher = watcher;

    m_progress->setVisible(true);
    connect(watcher, &QFutureWatcher<void>::finished, this, &DitherAnalyser::onSharpnessFinished);
    connect(watcher, &QFutureWatcher<void>::progressRangeChanged, m_progress, &QProgressBar::setRange);
    connect(watcher, &QFutureWatcher<void>::progressValueChanged, this, [this, n](int value) {
        m_progress->setValue(value);
        m_progress->setFormat(QString("Sharpness %1 / %2").arg(value).arg(n));
    });

    QFuture<void> future = QtConcurrent::map(*indices, [paths, results](int i) {
        results[i] = computeSharpnessFromFile(paths[i]);
    });

    watcher->setProperty("_indices", QVariant::fromValue(indices));
    watcher->setFuture(future);
}

void DitherAnalyser::onSharpnessFinished()
{
    for (int i = 0; i < (int)m_frames.size() && i < (int)m_sharpResults.size(); ++i)
        m_frames[i].sharpness = m_sharpResults[i];

    m_progress->setVisible(false);

    bool hasSharpness = false;
    for (const auto &f : m_frames)
        if (!std::isnan(f.sharpness)) { hasSharpness = true; break; }
    m_binBtn->setEnabled(hasSharpness);

    plotSharpness(m_frames);
    populateTable(m_frames);
    updateSummary(m_frames);

    m_sharpWatcher->deleteLater();
    m_sharpWatcher = nullptr;
}

void DitherAnalyser::onTabChanged(int index)
{
#ifdef HAVE_STELLARSOLVER
    if (index == m_envTabIndex && !m_fwhmDone && !m_fwhmRunning && !m_frames.empty())
        startFwhmExtraction();
#else
    Q_UNUSED(index);
#endif
}

#ifdef HAVE_STELLARSOLVER
void DitherAnalyser::startFwhmExtraction()
{
    if (m_frames.empty() || m_fwhmRunning) return;

    m_fwhmRunning = true;
    m_progress->setVisible(true);
    m_progress->setTextVisible(true);

    int n = (int)m_frames.size();
    m_fwhmResults.assign(n, std::nan(""));

    // Build file paths and WCS scales for the worker threads
    QDir dir(m_currentDirectory);
    QStringList paths;
    std::vector<double> cd1, cd2;
    for (int i = 0; i < n; ++i) {
        paths << dir.absoluteFilePath(QString::fromStdString(m_frames[i].filename));
        cd1.push_back(m_frames[i].cdelt1);
        cd2.push_back(m_frames[i].cdelt2);
    }

    // Build index list for QtConcurrent::map
    auto indices = std::make_shared<QList<int>>();
    for (int i = 0; i < n; ++i)
        indices->append(i);

    double *results = m_fwhmResults.data();
    auto *watcher = new QFutureWatcher<void>(this);
    m_fwhmWatcher = watcher;

    connect(watcher, &QFutureWatcher<void>::finished, this, &DitherAnalyser::onFwhmFinished);
    connect(watcher, &QFutureWatcher<void>::progressRangeChanged, m_progress, &QProgressBar::setRange);
    connect(watcher, &QFutureWatcher<void>::progressValueChanged, this, [this, n](int value) {
        m_progress->setValue(value);
        m_progress->setFormat(QString("Extracting FWHM %1 / %2").arg(value).arg(n));
    });

    QFuture<void> future = QtConcurrent::map(*indices, [paths, cd1, cd2, results](int i) {
        results[i] = extractFwhm(paths[i], cd1[i], cd2[i]);
    });

    // Keep indices alive until the future completes
    watcher->setProperty("_indices", QVariant::fromValue(indices));
    watcher->setFuture(future);
}

void DitherAnalyser::onFwhmFinished()
{
    // Copy results into frames
    for (int i = 0; i < (int)m_frames.size() && i < (int)m_fwhmResults.size(); ++i)
        m_frames[i].fwhm = m_fwhmResults[i];

    m_progress->setVisible(false);
    m_fwhmRunning = false;
    m_fwhmDone = true;
    m_binBtn->setEnabled(true);

    plotFwhm(m_frames);
    populateTable(m_frames);
    updateSummary(m_frames);

    m_fwhmWatcher->deleteLater();
    m_fwhmWatcher = nullptr;
}
#endif // HAVE_STELLARSOLVER

void DitherAnalyser::binByQuality()
{
    if (m_frames.empty() || m_currentDirectory.isEmpty())
        return;

    // Check which metrics are available
    bool hasSharpness = false, hasFwhm = false;
    for (const auto &f : m_frames) {
        if (!std::isnan(f.sharpness)) hasSharpness = true;
        if (!std::isnan(f.fwhm))      hasFwhm = true;
    }

    if (!hasSharpness && !hasFwhm) {
        QMessageBox::information(this, "Bin by Quality",
            "No quality metrics available.\n"
            "Sharpness is computed during loading.\n"
            "FWHM requires clicking the Star FWHM tab first.");
        return;
    }

    // Let user choose metric
    enum Metric { SHARP, FWHM_METRIC, COMBINED };
    Metric metric = SHARP;

    if (hasSharpness && hasFwhm) {
        QMessageBox box(this);
        box.setWindowTitle("Bin by Quality");
        box.setText("Which metric should be used for quality binning?");
        auto *btnSharp = box.addButton("Sharpness (HF)", QMessageBox::AcceptRole);
        auto *btnFwhm  = box.addButton("FWHM", QMessageBox::AcceptRole);
        auto *btnBoth  = box.addButton("Combined", QMessageBox::AcceptRole);
        box.addButton(QMessageBox::Cancel);
        box.exec();
        if (box.clickedButton() == btnFwhm)       metric = FWHM_METRIC;
        else if (box.clickedButton() == btnBoth)   metric = COMBINED;
        else if (box.clickedButton() == btnSharp)  metric = SHARP;
        else return;
    } else if (hasFwhm) {
        metric = FWHM_METRIC;
    }

    // Build sorted list of (index, score) — higher score = better quality
    struct IndexedScore { int index; double score; };
    std::vector<IndexedScore> scored;
    int skipped = 0;

    for (int i = 0; i < (int)m_frames.size(); ++i) {
        double score = std::nan("");
        const auto &f = m_frames[i];

        if (metric == SHARP && !std::isnan(f.sharpness)) {
            score = f.sharpness;  // higher = better
        } else if (metric == FWHM_METRIC && !std::isnan(f.fwhm)) {
            score = -f.fwhm;     // lower FWHM = better, so negate
        } else if (metric == COMBINED && !std::isnan(f.sharpness) && !std::isnan(f.fwhm)) {
            // Normalise both to z-scores, then average
            score = 0; // placeholder, normalised below
        }

        if (!std::isnan(score))
            scored.push_back({i, score});
        else
            skipped++;
    }

    // For combined metric, normalise both to z-scores
    if (metric == COMBINED && !scored.empty()) {
        double sumS = 0, sumF = 0;
        for (auto &s : scored) {
            sumS += m_frames[s.index].sharpness;
            sumF += m_frames[s.index].fwhm;
        }
        double meanS = sumS / scored.size(), meanF = sumF / scored.size();
        double varS = 0, varF = 0;
        for (auto &s : scored) {
            double ds = m_frames[s.index].sharpness - meanS;
            double df = m_frames[s.index].fwhm - meanF;
            varS += ds * ds;
            varF += df * df;
        }
        double stdS = std::sqrt(varS / scored.size());
        double stdF = std::sqrt(varF / scored.size());
        if (stdS < 1e-12) stdS = 1;
        if (stdF < 1e-12) stdF = 1;

        for (auto &s : scored) {
            double zSharp = (m_frames[s.index].sharpness - meanS) / stdS;
            double zFwhm  = -(m_frames[s.index].fwhm - meanF) / stdF; // negate: low FWHM = good
            s.score = (zSharp + zFwhm) / 2.0;
        }
    }

    if (scored.empty()) {
        QMessageBox::information(this, "Bin by Quality",
            "No frames have the selected metric.");
        return;
    }

    // Sort by score (best first)
    std::sort(scored.begin(), scored.end(),
              [](const IndexedScore &a, const IndexedScore &b) { return a.score > b.score; });

    // Assign quintile bins
    struct Bin {
        const char *name;
        std::vector<int> frameIndices;
    };
    Bin bins[] = {
        {"1_Excellent", {}},
        {"2_Good",      {}},
        {"3_Fair",      {}},
        {"4_Mediocre",  {}},
        {"5_Poor",      {}},
    };
    int nBins = 5;

    for (int i = 0; i < (int)scored.size(); ++i) {
        int bin = i * nBins / (int)scored.size();
        if (bin >= nBins) bin = nBins - 1;
        bins[bin].frameIndices.push_back(scored[i].index);
    }

    // Build detail text
    const char *metricName = (metric == SHARP) ? "Sharpness"
                           : (metric == FWHM_METRIC) ? "FWHM"
                           : "Combined";
    QString details = QString("Metric: %1\n\n").arg(metricName);
    for (int b = 0; b < nBins; ++b) {
        if (bins[b].frameIndices.empty()) continue;
        // Find score range
        double lo = 1e9, hi = -1e9;
        for (int idx : bins[b].frameIndices) {
            double v = (metric == FWHM_METRIC) ? m_frames[idx].fwhm : m_frames[idx].sharpness;
            if (!std::isnan(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
        }
        if (metric == FWHM_METRIC) {
            details += QString("%1: %2 frames (FWHM %3\" - %4\")\n")
                .arg(bins[b].name).arg(bins[b].frameIndices.size())
                .arg(lo, 0, 'f', 2).arg(hi, 0, 'f', 2);
        } else {
            details += QString("%1: %2 frames (sharpness %3 - %4)\n")
                .arg(bins[b].name).arg(bins[b].frameIndices.size())
                .arg(lo, 0, 'f', 4).arg(hi, 0, 'f', 4);
        }
    }
    if (skipped > 0)
        details += QString("\n%1 frames skipped (no data)\n").arg(skipped);

    auto answer = QMessageBox::question(this, "Bin by Quality",
        QString("Create symlinks in quality subfolders?\n\n%1\n"
                "Originals will not be moved.").arg(details),
        QMessageBox::Ok | QMessageBox::Cancel);

    if (answer != QMessageBox::Ok)
        return;

    QDir baseDir(m_currentDirectory);
    int linked = 0, errors = 0;

    for (int b = 0; b < nBins; ++b) {
        if (bins[b].frameIndices.empty()) continue;

        QString subdir = bins[b].name;
        baseDir.mkpath(subdir);

        for (int idx : bins[b].frameIndices) {
            QString srcFile = QString::fromStdString(m_frames[idx].filename);
            QString srcPath = baseDir.absoluteFilePath(srcFile);
            QString dstPath = baseDir.absoluteFilePath(subdir + "/" + srcFile);

            if (QFile::exists(dstPath))
                QFile::remove(dstPath);

            if (QFile::link(srcPath, dstPath))
                linked++;
            else
                errors++;
        }
    }

    QMessageBox::information(this, "Bin by Quality",
        QString("Done: %1 symlinks created, %2 errors.")
            .arg(linked).arg(errors));
}

#ifdef HAVE_STELLARSOLVER
void DitherAnalyser::updateVerificationSummary()
{
    int solved = 0, failed = 0;
    double sumErr = 0, maxErr = 0;
    std::vector<double> errors;

    for (const auto &f : m_frames) {
        if (!f.verified) continue;
        if (f.solveSuccess) {
            solved++;
            double err = f.totalErrorArcsec;
            sumErr += err;
            maxErr = std::max(maxErr, err);
            errors.push_back(err);
        } else {
            failed++;
        }
    }

    if (solved == 0 && failed == 0) return;

    double meanErr = (solved > 0) ? sumErr / solved : 0;
    double medianErr = 0;
    if (!errors.empty()) {
        std::sort(errors.begin(), errors.end());
        medianErr = errors[errors.size() / 2];
    }

    int badCount = 0;
    for (double e : errors)
        if (e > 60.0) badCount++;

    QString quality;
    if (meanErr < 5.0)
        quality = "EXCELLENT - headers very accurate";
    else if (meanErr < 15.0)
        quality = "GOOD - headers are reliable";
    else if (meanErr < 60.0)
        quality = "FAIR - some header inaccuracy";
    else
        quality = "POOR - headers significantly disagree with plate solutions";

    QString html = QString(
        "<b>Plate-Solve Verification Results</b><br>"
        "Solved: %1 / %2 &nbsp;&nbsp; Failed: %3<br>"
        "<br>"
        "<b>Header vs Solved Position Error:</b><br>"
        "&nbsp;&nbsp;Mean: %4\" &nbsp;&nbsp; Median: %5\" &nbsp;&nbsp; Max: %6\"<br>"
        "&nbsp;&nbsp;Frames with error &gt; 60\": %7<br>"
        "<br>"
        "<b>Verification Quality: </b>%8"
    )
    .arg(solved).arg(m_frames.size()).arg(failed)
    .arg(meanErr, 0, 'f', 2)
    .arg(medianErr, 0, 'f', 2)
    .arg(maxErr, 0, 'f', 2)
    .arg(badCount)
    .arg(quality);

    m_verifySummaryLabel->setText(html);
    m_verifySummaryLabel->setVisible(true);
}
#endif // HAVE_STELLARSOLVER

// ---------------------------------------------------------------------------
// Analysis
// ---------------------------------------------------------------------------

void DitherAnalyser::analyse(std::vector<FrameInfo> &frames)
{
    // Sort by date
    std::sort(frames.begin(), frames.end(), [](const FrameInfo &a, const FrameInfo &b) {
        return a.dateObs < b.dateObs;
    });

    // Compute mean RA/Dec
    double meanRA = 0, meanDec = 0;
    for (auto &f : frames) { meanRA += f.ra; meanDec += f.dec; }
    meanRA /= frames.size();
    meanDec /= frames.size();

    double cosDec = std::cos(meanDec * M_PI / 180.0);

    // Compute offsets from mean in arcseconds
    for (auto &f : frames) {
        f.dRA  = (f.ra - meanRA) * cosDec * 3600.0;
        f.dDec = (f.dec - meanDec) * 3600.0;
    }

    // Compute frame-to-frame dither distance
    frames[0].ditherDist = 0;
    for (size_t i = 1; i < frames.size(); ++i) {
        double dra  = frames[i].dRA  - frames[i-1].dRA;
        double ddec = frames[i].dDec - frames[i-1].dDec;
        frames[i].ditherDist = std::sqrt(dra*dra + ddec*ddec);
    }
}

void DitherAnalyser::populateTable(const std::vector<FrameInfo> &frames)
{
    m_table->clear();

    QStringList headers = {"Frame", "Date-Obs", "RA (deg)", "Dec (deg)",
                           "dRA (\")", "dDec (\")", "Dither Step (\")", "Exposure",
                           "Temp (C)", "Focus", "FWHM (\")", "Sharpness"};

#ifdef HAVE_STELLARSOLVER
    // Check if any frames have verification data
    bool hasVerification = false;
    for (const auto &f : frames) {
        if (f.verified) { hasVerification = true; break; }
    }
    if (hasVerification) {
        headers << "Solved RA" << "Solved Dec" << "Error (\")" << "Status";
    }
#endif

    m_table->setColumnCount(headers.size());
    m_table->setHorizontalHeaderLabels(headers);
    m_table->setRowCount(static_cast<int>(frames.size()));

    for (int i = 0; i < (int)frames.size(); ++i) {
        const auto &f = frames[i];
        m_table->setItem(i, 0, new QTableWidgetItem(QString::fromStdString(f.filename)));
        m_table->setItem(i, 1, new QTableWidgetItem(QString::fromStdString(f.dateObs)));
        m_table->setItem(i, 2, new QTableWidgetItem(QString::number(f.ra, 'f', 6)));
        m_table->setItem(i, 3, new QTableWidgetItem(QString::number(f.dec, 'f', 6)));
        m_table->setItem(i, 4, new QTableWidgetItem(QString::number(f.dRA, 'f', 2)));
        m_table->setItem(i, 5, new QTableWidgetItem(QString::number(f.dDec, 'f', 2)));
        m_table->setItem(i, 6, new QTableWidgetItem(QString::number(f.ditherDist, 'f', 2)));
        m_table->setItem(i, 7, new QTableWidgetItem(QString::number(f.exptime, 'f', 1)));
        m_table->setItem(i, 8, new QTableWidgetItem(
            std::isnan(f.temperature) ? "--" : QString::number(f.temperature, 'f', 1)));
        m_table->setItem(i, 9, new QTableWidgetItem(
            std::isnan(f.focusPos) ? "--" : QString::number(f.focusPos, 'f', 0)));
        m_table->setItem(i, 10, new QTableWidgetItem(
            std::isnan(f.fwhm) ? "--" : QString::number(f.fwhm, 'f', 2)));
        m_table->setItem(i, 11, new QTableWidgetItem(
            std::isnan(f.sharpness) ? "--" : QString::number(f.sharpness, 'f', 4)));

#ifdef HAVE_STELLARSOLVER
        if (hasVerification && f.verified) {
            if (f.solveSuccess) {
                m_table->setItem(i, 12, new QTableWidgetItem(QString::number(f.solvedRA, 'f', 6)));
                m_table->setItem(i, 13, new QTableWidgetItem(QString::number(f.solvedDec, 'f', 6)));
                m_table->setItem(i, 14, new QTableWidgetItem(QString::number(f.totalErrorArcsec, 'f', 2)));

                auto *statusItem = new QTableWidgetItem;
                if (f.totalErrorArcsec > 60.0) {
                    statusItem->setText("BAD");
                    statusItem->setForeground(QColor(243, 139, 168));
                } else if (f.totalErrorArcsec > 10.0) {
                    statusItem->setText("WARN");
                    statusItem->setForeground(QColor(249, 226, 175));
                } else {
                    statusItem->setText("OK");
                    statusItem->setForeground(QColor(166, 227, 161));
                }
                m_table->setItem(i, 15, statusItem);
            } else {
                m_table->setItem(i, 12, new QTableWidgetItem("--"));
                m_table->setItem(i, 13, new QTableWidgetItem("--"));
                m_table->setItem(i, 14, new QTableWidgetItem("--"));
                auto *statusItem = new QTableWidgetItem("FAIL");
                statusItem->setForeground(QColor(243, 139, 168));
                m_table->setItem(i, 15, statusItem);
            }
        }
#endif
    }

    m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
}

void DitherAnalyser::plotDitherScatter(const std::vector<FrameInfo> &frames)
{
    auto *chart = new QChart;
    chart->setTitle("Dither Pattern (offset from mean centre)");
    chart->setTheme(QChart::ChartThemeDark);

    // Scatter of all positions
    auto *scatter = new QScatterSeries;
    scatter->setName("Frame positions");
    scatter->setMarkerSize(8);
    scatter->setColor(QColor(137, 180, 250));

    for (auto &f : frames)
        scatter->append(f.dRA, f.dDec);

    chart->addSeries(scatter);

    // Connect consecutive frames with lines
    auto *path = new QLineSeries;
    path->setName("Sequence");
    path->setColor(QColor(137, 180, 250, 80));
    for (auto &f : frames)
        path->append(f.dRA, f.dDec);
    chart->addSeries(path);

    // Mark the centre
    auto *centre = new QScatterSeries;
    centre->setName("Mean centre");
    centre->setMarkerSize(14);
    centre->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
    centre->setColor(QColor(250, 179, 135));
    centre->append(0, 0);
    chart->addSeries(centre);

    // Mark first and last
    auto *endpoints = new QScatterSeries;
    endpoints->setName("First / Last");
    endpoints->setMarkerSize(12);
    endpoints->setColor(QColor(166, 227, 161));
    endpoints->append(frames.front().dRA, frames.front().dDec);
    endpoints->append(frames.back().dRA, frames.back().dDec);
    chart->addSeries(endpoints);

    chart->createDefaultAxes();
    auto axes = chart->axes();
    if (axes.size() >= 2) {
        auto *xAxis = qobject_cast<QValueAxis *>(axes[0]);
        auto *yAxis = qobject_cast<QValueAxis *>(axes[1]);
        if (xAxis) xAxis->setTitleText("dRA (\")");
        if (yAxis) yAxis->setTitleText("dDec (\")");

        // Make axes equal scale, centred on 0
        if (xAxis && yAxis) {
            double maxAbs = std::max({std::abs(xAxis->min()), std::abs(xAxis->max()),
                                      std::abs(yAxis->min()), std::abs(yAxis->max())}) * 1.1;
            xAxis->setRange(-maxAbs, maxAbs);
            yAxis->setRange(-maxAbs, maxAbs);
        }
    }

    auto *oldChart = m_scatterView->chart();


    m_scatterView->setChart(chart);


    delete oldChart;
}

void DitherAnalyser::plotDitherTimeline(const std::vector<FrameInfo> &frames)
{
    auto *chart = new QChart;
    chart->setTitle("Dither Offset Over Time");
    chart->setTheme(QChart::ChartThemeDark);

    auto *raSeries = new QLineSeries;
    raSeries->setName("dRA (\")");
    raSeries->setColor(QColor(137, 180, 250));

    auto *decSeries = new QLineSeries;
    decSeries->setName("dDec (\")");
    decSeries->setColor(QColor(245, 194, 231));

    auto *distSeries = new QLineSeries;
    distSeries->setName("Step size (\")");
    distSeries->setColor(QColor(166, 227, 161));

    for (int i = 0; i < (int)frames.size(); ++i) {
        raSeries->append(i, frames[i].dRA);
        decSeries->append(i, frames[i].dDec);
        distSeries->append(i, frames[i].ditherDist);
    }

    chart->addSeries(raSeries);
    chart->addSeries(decSeries);
    chart->addSeries(distSeries);
    chart->createDefaultAxes();

    auto axes = chart->axes();
    if (axes.size() >= 2) {
        auto *xAxis = qobject_cast<QValueAxis *>(axes[0]);
        auto *yAxis = qobject_cast<QValueAxis *>(axes[1]);
        if (xAxis) xAxis->setTitleText("Frame #");
        if (yAxis) yAxis->setTitleText("Arcseconds");
    }

    auto *oldChart = m_timelineView->chart();


    m_timelineView->setChart(chart);


    delete oldChart;
}

void DitherAnalyser::plotDitherRose(const std::vector<FrameInfo> &frames)
{
    // Histogram of dither step directions (16 bins, 22.5 deg each)
    const int nBins = 16;
    std::vector<int> bins(nBins, 0);
    int validSteps = 0;

    for (size_t i = 1; i < frames.size(); ++i) {
        double dra  = frames[i].dRA - frames[i-1].dRA;
        double ddec = frames[i].dDec - frames[i-1].dDec;
        double dist = std::sqrt(dra*dra + ddec*ddec);
        if (dist < 0.1) continue; // skip negligible moves

        double angle = std::atan2(ddec, dra) * 180.0 / M_PI;
        if (angle < 0) angle += 360.0;
        int bin = static_cast<int>(angle / (360.0 / nBins)) % nBins;
        bins[bin]++;
        validSteps++;
    }

    auto *chart = new QChart;
    chart->setTitle("Dither Direction Distribution");
    chart->setTheme(QChart::ChartThemeDark);

    auto *set = new QBarSet("Steps");
    set->setColor(QColor(137, 180, 250));
    QStringList categories;
    for (int i = 0; i < nBins; ++i) {
        *set << bins[i];
        double angle = i * (360.0 / nBins);
        categories << QString("%1%2").arg(angle, 0, 'f', 0).arg(QChar(0x00B0));
    }

    auto *series = new QBarSeries;
    series->append(set);
    chart->addSeries(series);

    auto *axisX = new QBarCategoryAxis;
    axisX->append(categories);
    axisX->setTitleText("Direction");
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    auto *axisY = new QValueAxis;
    axisY->setTitleText("Count");
    axisY->setLabelFormat("%d");
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    auto *oldChart = m_roseView->chart();


    m_roseView->setChart(chart);


    delete oldChart;
}

void DitherAnalyser::plotFieldRotation(const std::vector<FrameInfo> &frames)
{
    auto *chart = new QChart;
    chart->setTitle("Field Rotation Rate vs Time");
    chart->setTheme(QChart::ChartThemeDark);

    if (frames.size() < 2 || frames.front().dateObs.empty()) {
        auto *oldChart = m_rotationView->chart();

        m_rotationView->setChart(chart);

        delete oldChart;
        return;
    }

    QDateTime t0 = QDateTime::fromString(
        QString::fromStdString(frames.front().dateObs), Qt::ISODate);
    if (!t0.isValid()) {
        auto *oldChart = m_rotationView->chart();

        m_rotationView->setChart(chart);

        delete oldChart;
        return;
    }

    // Cumulative rotation angle series
    auto *angleSeries = new QLineSeries;
    angleSeries->setName("Rotation angle (\u00B0)");
    angleSeries->setColor(QColor(137, 180, 250));

    // Instantaneous rotation rate series (deg/min between consecutive frames)
    auto *rateSeries = new QLineSeries;
    rateSeries->setName("Rotation rate (\u00B0/min)");
    rateSeries->setColor(QColor(249, 226, 175));

    // Use CROTA2 for the rotation angle; unwrap to handle wrapping
    double prevAngle = frames[0].crota2;
    double cumAngle = 0;

    double minutesFirst = 0;
    angleSeries->append(minutesFirst, cumAngle);

    for (size_t i = 1; i < frames.size(); ++i) {
        QDateTime ti = QDateTime::fromString(
            QString::fromStdString(frames[i].dateObs), Qt::ISODate);
        if (!ti.isValid()) continue;

        double minutes = t0.msecsTo(ti) / 60000.0;

        // Unwrap angle difference to [-180, 180]
        double dAngle = frames[i].crota2 - prevAngle;
        while (dAngle > 180.0)  dAngle -= 360.0;
        while (dAngle < -180.0) dAngle += 360.0;
        cumAngle += dAngle;
        prevAngle = frames[i].crota2;

        angleSeries->append(minutes, cumAngle);

        // Rate: use elapsed time between this frame and previous
        QDateTime tPrev = QDateTime::fromString(
            QString::fromStdString(frames[i-1].dateObs), Qt::ISODate);
        if (tPrev.isValid()) {
            double dt = tPrev.msecsTo(ti) / 60000.0;
            if (dt > 0.001) {
                double rate = dAngle / dt;  // deg/min
                rateSeries->append(minutes, rate);
            }
        }
    }

    chart->addSeries(angleSeries);
    chart->addSeries(rateSeries);
    chart->createDefaultAxes();

    auto axes = chart->axes();
    if (axes.size() >= 2) {
        auto *xAxis = qobject_cast<QValueAxis *>(axes[0]);
        auto *yAxis = qobject_cast<QValueAxis *>(axes[1]);
        if (xAxis) xAxis->setTitleText("Time (minutes)");
        if (yAxis) yAxis->setTitleText("Degrees");
    }

    auto *oldChart = m_rotationView->chart();


    m_rotationView->setChart(chart);


    delete oldChart;
}

void DitherAnalyser::plotEnvironment(const std::vector<FrameInfo> &frames)
{
    auto *chart = new QChart;
    chart->setTitle("Temperature & Focus");
    chart->setTheme(QChart::ChartThemeDark);

    if (frames.empty()) {
        auto *oldChart = m_envView->chart();

        m_envView->setChart(chart);

        delete oldChart;
        return;
    }

    bool hasTemp = false, hasFocus = false;
    for (const auto &f : frames) {
        if (!std::isnan(f.temperature)) hasTemp = true;
        if (!std::isnan(f.focusPos))    hasFocus = true;
    }

    if (!hasTemp && !hasFocus) {
        chart->setTitle("Temperature & Focus (no data in headers)");
        auto *oldChart = m_envView->chart();

        m_envView->setChart(chart);

        delete oldChart;
        return;
    }

    auto *axisX = new QValueAxis;
    axisX->setTitleText("Frame #");
    axisX->setRange(0, (int)frames.size() - 1);
    chart->addAxis(axisX, Qt::AlignBottom);

    auto *axisYLeft = new QValueAxis;
    chart->addAxis(axisYLeft, Qt::AlignLeft);

    if (hasTemp) {
        auto *tempSeries = new QLineSeries;
        tempSeries->setName("Temperature (\u00B0C)");
        tempSeries->setColor(QColor(249, 226, 175));
        for (int i = 0; i < (int)frames.size(); ++i)
            if (!std::isnan(frames[i].temperature))
                tempSeries->append(i, frames[i].temperature);
        chart->addSeries(tempSeries);
        tempSeries->attachAxis(axisX);
        tempSeries->attachAxis(axisYLeft);
        axisYLeft->setTitleText("Temperature (\u00B0C)");
    }

    if (hasFocus) {
        auto *axisYRight = new QValueAxis;
        axisYRight->setTitleText("Focus position");
        chart->addAxis(axisYRight, Qt::AlignRight);

        auto *focusSeries = new QLineSeries;
        focusSeries->setName("Focus pos");
        focusSeries->setColor(QColor(137, 180, 250));
        for (int i = 0; i < (int)frames.size(); ++i)
            if (!std::isnan(frames[i].focusPos))
                focusSeries->append(i, frames[i].focusPos);
        chart->addSeries(focusSeries);
        focusSeries->attachAxis(axisX);
        focusSeries->attachAxis(axisYRight);
    }

    if (!hasTemp)
        axisYLeft->setVisible(false);

    auto *oldChart = m_envView->chart();


    m_envView->setChart(chart);


    delete oldChart;
}

void DitherAnalyser::plotSharpness(const std::vector<FrameInfo> &frames)
{
    auto *chart = new QChart;
    chart->setTitle("HF Sharpness (DCT high-frequency power ratio)");
    chart->setTheme(QChart::ChartThemeDark);

    bool hasData = false;
    for (const auto &f : frames)
        if (!std::isnan(f.sharpness)) { hasData = true; break; }

    if (!hasData) {
        chart->setTitle("Sharpness (no data)");
        auto *oldChart = m_sharpView->chart();

        m_sharpView->setChart(chart);

        delete oldChart;
        return;
    }

    auto *series = new QLineSeries;
    series->setName("Sharpness");
    series->setColor(QColor(137, 180, 250));

    auto *scatter = new QScatterSeries;
    scatter->setName("Per frame");
    scatter->setMarkerSize(6);
    scatter->setColor(QColor(137, 180, 250));

    for (int i = 0; i < (int)frames.size(); ++i) {
        if (!std::isnan(frames[i].sharpness)) {
            series->append(i, frames[i].sharpness);
            scatter->append(i, frames[i].sharpness);
        }
    }

    chart->addSeries(series);
    chart->addSeries(scatter);
    chart->createDefaultAxes();

    auto axes = chart->axes();
    if (axes.size() >= 2) {
        auto *xAxis = qobject_cast<QValueAxis *>(axes[0]);
        auto *yAxis = qobject_cast<QValueAxis *>(axes[1]);
        if (xAxis) xAxis->setTitleText("Frame #");
        if (yAxis) yAxis->setTitleText("HF Power Ratio (higher = sharper)");
    }

    auto *oldChart = m_sharpView->chart();


    m_sharpView->setChart(chart);


    delete oldChart;
}

#ifdef HAVE_STELLARSOLVER
void DitherAnalyser::plotFwhm(const std::vector<FrameInfo> &frames)
{
    auto *chart = new QChart;
    chart->setTitle("Star FWHM Over Session");
    chart->setTheme(QChart::ChartThemeDark);

    bool hasFwhm = false;
    for (const auto &f : frames)
        if (!std::isnan(f.fwhm)) { hasFwhm = true; break; }

    if (!hasFwhm) {
        chart->setTitle("Star FWHM (no data yet)");
        auto *oldChart = m_fwhmView->chart();

        m_fwhmView->setChart(chart);

        delete oldChart;
        return;
    }

    auto *fwhmSeries = new QLineSeries;
    fwhmSeries->setName("FWHM (\")");
    fwhmSeries->setColor(QColor(166, 227, 161));

    auto *fwhmScatter = new QScatterSeries;
    fwhmScatter->setName("FWHM points");
    fwhmScatter->setMarkerSize(6);
    fwhmScatter->setColor(QColor(166, 227, 161));

    for (int i = 0; i < (int)frames.size(); ++i) {
        if (!std::isnan(frames[i].fwhm)) {
            fwhmSeries->append(i, frames[i].fwhm);
            fwhmScatter->append(i, frames[i].fwhm);
        }
    }

    chart->addSeries(fwhmSeries);
    chart->addSeries(fwhmScatter);
    chart->createDefaultAxes();

    auto axes = chart->axes();
    if (axes.size() >= 2) {
        auto *xAxis = qobject_cast<QValueAxis *>(axes[0]);
        auto *yAxis = qobject_cast<QValueAxis *>(axes[1]);
        if (xAxis) xAxis->setTitleText("Frame #");
        if (yAxis) yAxis->setTitleText("FWHM (arcsec)");
    }

    auto *oldChart = m_fwhmView->chart();


    m_fwhmView->setChart(chart);


    delete oldChart;
}
#else
void DitherAnalyser::plotFwhm(const std::vector<FrameInfo> &) {}
#endif

void DitherAnalyser::updateSummary(const std::vector<FrameInfo> &frames)
{
    if (frames.empty()) return;

    // Pixel scale from first frame (arcsec/pixel)
    double pixScale = std::abs(frames[0].cdelt1) * 3600.0;
    if (pixScale < 0.001) pixScale = std::abs(frames[0].cdelt2) * 3600.0;
    if (pixScale < 0.001) pixScale = 1.0;  // fallback to avoid divide-by-zero

    // Dither step statistics (skip frame 0)
    std::vector<double> steps;
    for (size_t i = 1; i < frames.size(); ++i)
        steps.push_back(frames[i].ditherDist);

    double meanStep = 0, maxStep = 0, minStep = 0, medianStep = 0, stdStep = 0;
    if (!steps.empty()) {
        meanStep = std::accumulate(steps.begin(), steps.end(), 0.0) / steps.size();
        maxStep = *std::max_element(steps.begin(), steps.end());
        minStep = *std::min_element(steps.begin(), steps.end());

        std::vector<double> sorted = steps;
        std::sort(sorted.begin(), sorted.end());
        medianStep = sorted[sorted.size() / 2];

        double variance = 0;
        for (double s : steps) variance += (s - meanStep) * (s - meanStep);
        variance /= steps.size();
        stdStep = std::sqrt(variance);
    }

    // Total spread
    double minRA = 1e9, maxRA = -1e9, minDec = 1e9, maxDec = -1e9;
    for (auto &f : frames) {
        minRA  = std::min(minRA, f.dRA);
        maxRA  = std::max(maxRA, f.dRA);
        minDec = std::min(minDec, f.dDec);
        maxDec = std::max(maxDec, f.dDec);
    }
    double spreadRA  = maxRA - minRA;
    double spreadDec = maxDec - minDec;
    double spreadTotal = std::sqrt(spreadRA*spreadRA + spreadDec*spreadDec);

    // Count "no-dither" frames (step < 1 pixel)
    int noMove = 0;
    for (double s : steps)
        if (s < pixScale) noMove++;

    // How many unique pixel positions (rounded to nearest pixel)
    std::set<std::pair<int,int>> pixPositions;
    for (auto &f : frames) {
        int px = static_cast<int>(std::round(f.dRA / pixScale));
        int py = static_cast<int>(std::round(f.dDec / pixScale));
        pixPositions.insert({px, py});
    }

    // Dither quality rating
    QString rating;
    double stepsPerPixel = meanStep / pixScale;
    if (stepsPerPixel < 1.0)
        rating = "POOR - dither steps smaller than pixel scale";
    else if (stepsPerPixel < 3.0)
        rating = "FAIR - small dithers, some hot-pixel benefit";
    else if (stepsPerPixel < 15.0)
        rating = "GOOD - effective dithering for hot-pixel rejection";
    else if (stepsPerPixel < 50.0)
        rating = "VERY GOOD - strong dithering, good sub-pixel sampling";
    else
        rating = "EXCELLENT - large dithers, maximum coverage diversity";

    // Check randomness: are directions well-distributed?
    const int nBins = 8;
    std::vector<int> dirBins(nBins, 0);
    for (size_t i = 1; i < frames.size(); ++i) {
        double dra  = frames[i].dRA - frames[i-1].dRA;
        double ddec = frames[i].dDec - frames[i-1].dDec;
        if (std::sqrt(dra*dra + ddec*ddec) < 0.1) continue;
        double angle = std::atan2(ddec, dra) * 180.0 / M_PI;
        if (angle < 0) angle += 360.0;
        int bin = static_cast<int>(angle / (360.0 / nBins)) % nBins;
        dirBins[bin]++;
    }
    int maxDir = *std::max_element(dirBins.begin(), dirBins.end());
    int minDir = *std::min_element(dirBins.begin(), dirBins.end());
    double dirUniformity = (maxDir > 0) ? (double)minDir / maxDir : 0;

    QString dirQuality;
    if (dirUniformity > 0.7)
        dirQuality = "Well-distributed (random)";
    else if (dirUniformity > 0.4)
        dirQuality = "Moderate directional bias";
    else
        dirQuality = "Strong directional bias (may indicate drift, not dither)";

    // Throughput analysis: exposure time vs elapsed time
    double totalExposure = 0;
    for (auto &f : frames)
        totalExposure += f.exptime;

    double elapsedSeconds = 0;
    double meanGap = 0, minGap = 1e9, maxGap = 0, medianGap = 0;
    std::vector<double> gaps;
    QString throughputHtml;

    if (frames.size() >= 2 &&
        !frames.front().dateObs.empty() && !frames.back().dateObs.empty())
    {
        QDateTime t0 = QDateTime::fromString(
            QString::fromStdString(frames.front().dateObs), Qt::ISODate);
        QDateTime tN = QDateTime::fromString(
            QString::fromStdString(frames.back().dateObs), Qt::ISODate);

        if (t0.isValid() && tN.isValid()) {
            // Elapsed = start of first to end of last exposure
            elapsedSeconds = t0.msecsTo(tN) / 1000.0 + frames.back().exptime;

            // Per-frame gaps (time between end of frame i and start of frame i+1)
            for (size_t i = 1; i < frames.size(); ++i) {
                QDateTime prev = QDateTime::fromString(
                    QString::fromStdString(frames[i-1].dateObs), Qt::ISODate);
                QDateTime curr = QDateTime::fromString(
                    QString::fromStdString(frames[i].dateObs), Qt::ISODate);
                if (prev.isValid() && curr.isValid()) {
                    double gap = prev.msecsTo(curr) / 1000.0 - frames[i-1].exptime;
                    if (gap < 0) gap = 0;
                    gaps.push_back(gap);
                }
            }

            if (!gaps.empty()) {
                meanGap = std::accumulate(gaps.begin(), gaps.end(), 0.0) / gaps.size();
                minGap = *std::min_element(gaps.begin(), gaps.end());
                maxGap = *std::max_element(gaps.begin(), gaps.end());
                std::vector<double> sortedGaps = gaps;
                std::sort(sortedGaps.begin(), sortedGaps.end());
                medianGap = sortedGaps[sortedGaps.size() / 2];
            }

            double throughputPct = (elapsedSeconds > 0)
                ? 100.0 * totalExposure / elapsedSeconds : 0;

            // Format durations as h:mm:ss
            auto fmtDuration = [](double secs) -> QString {
                int h = static_cast<int>(secs) / 3600;
                int m = (static_cast<int>(secs) % 3600) / 60;
                int s = static_cast<int>(secs) % 60;
                return QString("%1:%2:%3")
                    .arg(h).arg(m, 2, 10, QChar('0')).arg(s, 2, 10, QChar('0'));
            };

            throughputHtml = QString(
                "<br>"
                "<b>Throughput:</b><br>"
                "&nbsp;&nbsp;Total exposure: %1 (%2 s)&nbsp;&nbsp; "
                "Elapsed: %3 (%4 s)<br>"
                "&nbsp;&nbsp;Efficiency: <b>%5%</b> "
                "(overhead: %6 s per frame)<br>"
                "&nbsp;&nbsp;Inter-frame gap: mean %7 s&nbsp;&nbsp; "
                "median %8 s&nbsp;&nbsp; min %9 s&nbsp;&nbsp; max %10 s"
            )
            .arg(fmtDuration(totalExposure))
            .arg(totalExposure, 0, 'f', 0)
            .arg(fmtDuration(elapsedSeconds))
            .arg(elapsedSeconds, 0, 'f', 0)
            .arg(throughputPct, 0, 'f', 1)
            .arg(meanGap, 0, 'f', 1)
            .arg(meanGap, 0, 'f', 1)
            .arg(medianGap, 0, 'f', 1)
            .arg(minGap, 0, 'f', 1)
            .arg(maxGap, 0, 'f', 1);
        }
    }

    QString html = QString(
        "<b>Session: </b>%1 frames, %2 &ndash; %3<br>"
        "<b>Pixel scale: </b>%4\"/px<br>"
        "<br>"
        "<b>Dither Step Statistics:</b><br>"
        "&nbsp;&nbsp;Mean: %5\"&nbsp;&nbsp; Median: %6\"&nbsp;&nbsp; Std: %7\"<br>"
        "&nbsp;&nbsp;Min: %8\"&nbsp;&nbsp; Max: %9\"<br>"
        "&nbsp;&nbsp;Mean step: %10 pixels<br>"
        "<br>"
        "<b>Coverage Spread:</b><br>"
        "&nbsp;&nbsp;RA range: %11\"&nbsp;&nbsp; Dec range: %12\"&nbsp;&nbsp; Diagonal: %13\"<br>"
        "&nbsp;&nbsp;Unique pixel positions: %14<br>"
        "&nbsp;&nbsp;Stationary frames (step &lt; 1px): %15 / %16 (%17%)<br>"
        "<br>"
        "<b>Direction uniformity: </b>%18 (min/max ratio: %19)<br>"
        "<b>Dither Quality Rating: </b>%20"
    )
    .arg(frames.size())
    .arg(QString::fromStdString(frames.front().dateObs),
         QString::fromStdString(frames.back().dateObs))
    .arg(pixScale, 0, 'f', 2)
    .arg(meanStep, 0, 'f', 2)
    .arg(medianStep, 0, 'f', 2)
    .arg(stdStep, 0, 'f', 2)
    .arg(minStep, 0, 'f', 2)
    .arg(maxStep, 0, 'f', 2)
    .arg(stepsPerPixel, 0, 'f', 1)
    .arg(spreadRA, 0, 'f', 2)
    .arg(spreadDec, 0, 'f', 2)
    .arg(spreadTotal, 0, 'f', 2)
    .arg(pixPositions.size())
    .arg(noMove).arg(steps.size())
    .arg(steps.size() > 0 ? 100.0 * noMove / steps.size() : 0, 0, 'f', 1)
    .arg(dirQuality)
    .arg(dirUniformity, 0, 'f', 2)
    .arg(rating);

    html += throughputHtml;

    // Environment, sharpness & FWHM summary
    {
        std::vector<double> temps, focuses, fwhms, sharps;
        for (const auto &f : frames) {
            if (!std::isnan(f.temperature)) temps.push_back(f.temperature);
            if (!std::isnan(f.focusPos))    focuses.push_back(f.focusPos);
            if (!std::isnan(f.fwhm))        fwhms.push_back(f.fwhm);
            if (!std::isnan(f.sharpness))   sharps.push_back(f.sharpness);
        }

        QString envHtml;

        if (!temps.empty()) {
            double tMin = *std::min_element(temps.begin(), temps.end());
            double tMax = *std::max_element(temps.begin(), temps.end());
            double tMean = std::accumulate(temps.begin(), temps.end(), 0.0) / temps.size();
            envHtml += QString(
                "<br><b>Temperature:</b> mean %1\u00B0C &nbsp;&nbsp; "
                "range %2\u00B0C \u2013 %3\u00B0C &nbsp;&nbsp; "
                "\u0394 %4\u00B0C")
                .arg(tMean, 0, 'f', 1).arg(tMin, 0, 'f', 1)
                .arg(tMax, 0, 'f', 1).arg(tMax - tMin, 0, 'f', 1);
        }

        if (!focuses.empty()) {
            double fMin = *std::min_element(focuses.begin(), focuses.end());
            double fMax = *std::max_element(focuses.begin(), focuses.end());
            envHtml += QString(
                "<br><b>Focus position:</b> range %1 \u2013 %2 &nbsp;&nbsp; "
                "\u0394 %3 steps")
                .arg(fMin, 0, 'f', 0).arg(fMax, 0, 'f', 0).arg(fMax - fMin, 0, 'f', 0);
        }

        if (!sharps.empty()) {
            std::sort(sharps.begin(), sharps.end());
            double sMean = std::accumulate(sharps.begin(), sharps.end(), 0.0) / sharps.size();
            double sMedian = sharps[sharps.size() / 2];
            envHtml += QString(
                "<br><b>HF Sharpness:</b> median %1 &nbsp;&nbsp; mean %2 &nbsp;&nbsp; "
                "range %3 \u2013 %4")
                .arg(sMedian, 0, 'f', 4).arg(sMean, 0, 'f', 4)
                .arg(sharps.front(), 0, 'f', 4).arg(sharps.back(), 0, 'f', 4);
        }

        if (!fwhms.empty()) {
            std::sort(fwhms.begin(), fwhms.end());
            double fMean = std::accumulate(fwhms.begin(), fwhms.end(), 0.0) / fwhms.size();
            double fMedian = fwhms[fwhms.size() / 2];
            double fMin = fwhms.front();
            double fMax = fwhms.back();
            envHtml += QString(
                "<br><b>Star FWHM:</b> median %1\" &nbsp;&nbsp; mean %2\" &nbsp;&nbsp; "
                "range %3\" \u2013 %4\"")
                .arg(fMedian, 0, 'f', 2).arg(fMean, 0, 'f', 2)
                .arg(fMin, 0, 'f', 2).arg(fMax, 0, 'f', 2);

            // Seeing quality rating
            QString seeingQuality;
            if (fMedian < 2.0)
                seeingQuality = "EXCELLENT";
            else if (fMedian < 3.5)
                seeingQuality = "GOOD";
            else if (fMedian < 5.0)
                seeingQuality = "FAIR";
            else
                seeingQuality = "POOR";
            envHtml += QString(" &nbsp;&nbsp; (%1)").arg(seeingQuality);
        }

        html += envHtml;
    }

    m_summaryLabel->setText(html);
}
