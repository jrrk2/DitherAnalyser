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

#include <fitsio.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <set>

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

        fits_close_file(fptr, &status);

        if (!std::isnan(fi.ra) && !std::isnan(fi.dec))
            frames.push_back(fi);
    }

    emit finished(frames);
}

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

    // Summary label
    m_summaryLabel = new QLabel;
    m_summaryLabel->setWordWrap(true);
    m_summaryLabel->setStyleSheet("QLabel { background: #1e1e2e; color: #cdd6f4; padding: 10px; border-radius: 6px; font-family: monospace; font-size: 11pt; }");
    mainLayout->addWidget(m_summaryLabel);

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

    // Table
    m_table = new QTableWidget;
    m_table->setAlternatingRowColors(true);
    m_tabs->addTab(m_table, "Frame Table");

    mainLayout->addWidget(m_tabs, 1);
    setCentralWidget(central);
}

void DitherAnalyser::browseDirectory()
{
    QString dir = QFileDialog::getExistingDirectory(this, "Select FITS Directory",
        QDir::homePath() + "/Astrophotography");
    if (dir.isEmpty()) return;

    m_dirLabel->setText(dir);
    m_progress->setVisible(true);
    m_progress->setValue(0);

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

    analyse(m_frames);
    populateTable(m_frames);
    plotDitherScatter(m_frames);
    plotDitherTimeline(m_frames);
    plotDitherRose(m_frames);
    updateSummary(m_frames);
}

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
                           "dRA (\")", "dDec (\")", "Dither Step (\")", "Exposure"};
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

    delete m_scatterView->chart();
    m_scatterView->setChart(chart);
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

    delete m_timelineView->chart();
    m_timelineView->setChart(chart);
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

    delete m_roseView->chart();
    m_roseView->setChart(chart);
}

void DitherAnalyser::updateSummary(const std::vector<FrameInfo> &frames)
{
    if (frames.empty()) return;

    // Pixel scale from first frame (arcsec/pixel)
    double pixScale = std::abs(frames[0].cdelt1) * 3600.0;
    if (pixScale < 0.001) pixScale = std::abs(frames[0].cdelt2) * 3600.0;

    // Dither step statistics (skip frame 0)
    std::vector<double> steps;
    for (size_t i = 1; i < frames.size(); ++i)
        steps.push_back(frames[i].ditherDist);

    double meanStep = std::accumulate(steps.begin(), steps.end(), 0.0) / steps.size();
    double maxStep = *std::max_element(steps.begin(), steps.end());
    double minStep = *std::min_element(steps.begin(), steps.end());

    std::vector<double> sorted = steps;
    std::sort(sorted.begin(), sorted.end());
    double medianStep = sorted[sorted.size() / 2];

    double variance = 0;
    for (double s : steps) variance += (s - meanStep) * (s - meanStep);
    variance /= steps.size();
    double stdStep = std::sqrt(variance);

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

    m_summaryLabel->setText(html);
}
