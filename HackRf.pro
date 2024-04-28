QT       += core gui multimedia

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    audiooutput.cpp \
    hackrfdevice.cpp \
    main.cpp \
    mainwindow.cpp \
    plotter.cpp

HEADERS += \
    audiooutput.h \
    buffer.h \
    constants.h \
    hackrfdevice.h \
    mainwindow.h \
    plotter.h

FORMS += \
    mainwindow.ui

macos {
    message("macos enabled")
    QMAKE_INFO_PLIST = ./macos/Info.plist
    QMAKE_ASSET_CATALOGS = $$PWD/macos/Assets.xcassets
    QMAKE_ASSET_CATALOGS_APP_ICON = "AppIcon"

    # INCLUDEPATH += /opt/homebrew/Cellar/hackrf/2024.02.1/include
    # LIBS += -L/opt/homebrew/Cellar/hackrf/2024.02.1/lib -lhackrf
    INCLUDEPATH += /usr/local/include/
    LIBS += -L/usr/local/lib -lhackrf -lvolk
}

unix:!macx{
    message("linux enabled")
    INCLUDEPATH += /usr/include
    INCLUDEPATH += /usr/local/include
    INCLUDEPATH += /usr/lib
    INCLUDEPATH += /usr/local/lib
    INCLUDEPATH += /usr/lib/x86_64-linux-gnu
    LIBS += -L/usr/local/lib -lhackrf
}

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
