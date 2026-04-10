#pragma once

#include <QMainWindow>
#include <QTableWidget>
#include <QLabel>
#include <QProgressBar>
#include <QPushButton>
#include <QDoubleSpinBox>
#include <QThread>
#include <QFutureWatcher>
#include <QChartView>
#include <QScatterSeries>
#include <QLineSeries>
#include <QValueAxis>
#include <QPolarChart>
#include <vector>
#include <string>
#include <atomic>

#include <fitsio.h>
#ifdef HAVE_STELLARSOLVER
#include <stellarsolver.h>
#include <parameters.h>
#endif

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

    // Environmental / quality metrics
    double temperature = std::nan("");   // CCD or ambient temperature (C)
    double focusPos    = std::nan("");   // focuser position (steps)
    double fwhm        = std::nan("");   // median star FWHM (arcsec)
    double sharpness   = std::nan("");   // HF power ratio (0-1, higher = sharper)

#ifdef HAVE_STELLARSOLVER
    // Plate-solve verification results
    bool verified = false;
    bool solveSuccess = false;
    double solvedRA = std::nan("");
    double solvedDec = std::nan("");
    double solvedPixscale = std::nan("");
    double raErrorArcsec = std::nan("");   // header vs solved difference
    double decErrorArcsec = std::nan("");
    double totalErrorArcsec = std::nan("");
#endif
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

#ifdef HAVE_STELLARSOLVER
// Worker that plate-solves frames to verify RA/DEC headers
class SolverVerifier : public QObject {
    Q_OBJECT
public:
    explicit SolverVerifier(const std::vector<FrameInfo> &frames,
                            const QString &directory,
                            double searchRadius);
    void abort();

signals:
    void progress(int current, int total);
    void frameVerified(int index, bool success, double solvedRA, double solvedDec,
                       double pixscale, double raErr, double decErr, double totalErr);
    void finished();
    void error(const QString &msg);

public slots:
    void process();

private:
    bool loadFITSForSolver(StellarSolver *solver, const QString &fitsFile,
                           double hintRA, double hintDec, double searchRadius);
    QStringList findIndexFiles();

    std::vector<FrameInfo> m_frames;
    QString m_directory;
    double m_searchRadius;
    std::atomic<bool> m_abort{false};
};
#endif

class DitherAnalyser : public QMainWindow {
    Q_OBJECT
public:
    explicit DitherAnalyser(QWidget *parent = nullptr);
    ~DitherAnalyser();

private slots:
    void browseDirectory();
    void onFramesLoaded(std::vector<FrameInfo> frames);
    void onProgress(int current, int total);
#ifdef HAVE_STELLARSOLVER
    void startVerification();
    void onVerifyProgress(int current, int total);
    void onFrameVerified(int index, bool success, double solvedRA, double solvedDec,
                         double pixscale, double raErr, double decErr, double totalErr);
    void onVerifyFinished();
#endif
    void onTabChanged(int index);
    void startSharpnessExtraction();
    void onSharpnessFinished();
#ifdef HAVE_STELLARSOLVER
    void startFwhmExtraction();
    void onFwhmFinished();
#endif
    void binByQuality();

private:
    void buildUI();
    void analyse(std::vector<FrameInfo> &frames);
    void populateTable(const std::vector<FrameInfo> &frames);
    void plotDitherScatter(const std::vector<FrameInfo> &frames);
    void plotDitherTimeline(const std::vector<FrameInfo> &frames);
    void plotDitherRose(const std::vector<FrameInfo> &frames);
    void plotFieldRotation(const std::vector<FrameInfo> &frames);
    void plotEnvironment(const std::vector<FrameInfo> &frames);
    void plotSharpness(const std::vector<FrameInfo> &frames);
    void plotFwhm(const std::vector<FrameInfo> &frames);
    void updateSummary(const std::vector<FrameInfo> &frames);
#ifdef HAVE_STELLARSOLVER
    void updateVerificationSummary();
#endif

    // UI
    QLabel      *m_dirLabel;
    QProgressBar *m_progress;
    QTabWidget  *m_tabs;
    QTableWidget *m_table;
    QChartView  *m_scatterView;
    QChartView  *m_timelineView;
    QChartView  *m_roseView;
    QChartView  *m_rotationView;
    QChartView  *m_envView;
    QChartView  *m_sharpView;
    QChartView  *m_fwhmView;
    QLabel      *m_summaryLabel;

#ifdef HAVE_STELLARSOLVER
    // Verification UI
    QPushButton    *m_verifyBtn;
    QDoubleSpinBox *m_radiusSpin;
    QLabel         *m_verifySummaryLabel;
#endif

    // Quality binning UI
    QPushButton    *m_binBtn;

    QThread     *m_loaderThread = nullptr;
    FitsLoader  *m_loader = nullptr;

#ifdef HAVE_STELLARSOLVER
    QThread        *m_verifierThread = nullptr;
    SolverVerifier *m_verifier = nullptr;

    QFutureWatcher<void> *m_fwhmWatcher = nullptr;
    std::vector<double>  m_fwhmResults;
    bool            m_fwhmDone = false;
    bool            m_fwhmRunning = false;
#endif

    QFutureWatcher<void> *m_sharpWatcher = nullptr;
    std::vector<double>  m_sharpResults;

    int             m_envTabIndex = -1;

    std::vector<FrameInfo> m_frames;
    QString m_currentDirectory;
};
