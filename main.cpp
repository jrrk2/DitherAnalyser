#include "DitherAnalyser.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setStyle("Fusion");

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
