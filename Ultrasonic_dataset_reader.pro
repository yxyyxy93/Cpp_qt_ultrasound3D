QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    qcustomplot.cpp \
    ultrasound_cscan_process.cpp \
    utils.cpp

HEADERS += \
    fftw3.h \
    npy.hpp \
    qcustomplot.h \
    ultrasound_cscan_process.h \
    utils.h

FORMS +=

LIBS += "$$PWD/libfftw3-3.dll"

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target


