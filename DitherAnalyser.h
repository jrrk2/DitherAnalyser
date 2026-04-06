#pragma once

#include <QMainWindow>
#include <QTableWidget>
#include <QLabel>
#include <QProgressBar>
#include <QThread>
#include <QChartView>
#include <QScatterSeries>
#include <QLineSeries>
#include <QValueAxis>
#include <QPolarChart>
#include <vector>
#include <string>

struct FrameInfo {
    std::string filename;
    std::string dateObs;
    double ra;          // degrees
    double dec;         // degrees
    double crpix1, crpix2;
    double cdelt1, cdelt2;
    double crota2;
    double exptime;
    int naxis1, naxis2;
    // Computed offsets from mean centre (in arcseconds)
    double dRA;
    double dDec;
    double ditherDist;  // distance from previous frame in arcsec
};

// Worker that loads FITS headers off the main thread
class FitsLoader : public QObject {
    Q_OBJECT
public:
    explicit FitsLoader(const QString &directory);
signals:
    void progress(int current, int total);
    void finished(std::vector<FrameInfo> frames);
    void error(const QString &msg);
public slots:
    void process();
private:
    QString m_directory;
};

class DitherAnalyser : public QMainWindow {
    Q_OBJECT
public:
    explicit DitherAnalyser(QWidget *parent = nullptr);
    ~DitherAnalyser();

private slots:
    void browseDirectory();
    void onFramesLoaded(std::vector<FrameInfo> frames);
    void onProgress(int current, int total);

private:
    void buildUI();
    void analyse(std::vector<FrameInfo> &frames);
    void populateTable(const std::vector<FrameInfo> &frames);
    void plotDitherScatter(const std::vector<FrameInfo> &frames);
    void plotDitherTimeline(const std::vector<FrameInfo> &frames);
    void plotDitherRose(const std::vector<FrameInfo> &frames);
    void updateSummary(const std::vector<FrameInfo> &frames);

    // UI
    QLabel      *m_dirLabel;
    QProgressBar *m_progress;
    QTabWidget  *m_tabs;
    QTableWidget *m_table;
    QChartView  *m_scatterView;
    QChartView  *m_timelineView;
    QChartView  *m_roseView;
    QLabel      *m_summaryLabel;

    QThread     *m_loaderThread = nullptr;
    FitsLoader  *m_loader = nullptr;
    std::vector<FrameInfo> m_frames;
};
