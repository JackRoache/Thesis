QT += widgets
QT += core
QT += gui

CONFIG += c++11

TARGET = Thesis
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app
#LIBS += -lQt5Concurrent
SOURCES += main.cpp \
    bessel.cpp \
    clbessel.cpp \
    types.cpp \
    tissueproperties.cpp \
    phantomfile.cpp \
    bim.cpp

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0



HEADERS += \
    bessel.h \
    clbessel.h \
    helper.h \
    types.h \
    tissueproperties.h \
    phantomfile.h \
    bim.h

macx {
macx: LIBS += -L$$PWD/../../../../../usr/local/lib/ -laf.3.5.0
macx: LIBS += -framework OpenCL

INCLUDEPATH += /usr/local/include
DEPENDPATH += /usr/local/include
}



DISTFILES += \
    besselj0.cl

#LIBS += -L/home/jack/ArrayFire/lib/ -laf
#INCLUDEPATH += '/home/jack/ArrayFire/include'
#DEPENDPATH += '/home/jack/ArrayFire/include'

#LIBS += -L'/usr/lib/x86_64-linux-gnu' -lOpenCL

#INCLUDEPATH += '/opt/opencl-headers/include'
#DEPENDPATH += '/opt/opencl-headers/include'

#LIBS += -L'/home/jack/workspace/Thesis/' -lcomplex_bessel
#INCLUDEPATH += '/home/jack/workspace/Thesis/Complex_Bessel/'
#DEPENDPATH += '/home/jack/workspace/Thesis/Complex_Bessel/'

#win32: LIBS += -L$$PWD/'../../Program Files/ArrayFire/v3/lib/' -laf

#INCLUDEPATH += $$PWD/'../../Program Files/ArrayFire/v3/include'
#DEPENDPATH += $$PWD/'../../Program Files/ArrayFire/v3/include'

#win32: LIBS += -L$$PWD/'../../Program Files (x86)/AMD APP SDK/3.0/lib/x86_64/' -lOpenCL

#INCLUDEPATH += $$PWD/'../../Program Files (x86)/AMD APP SDK/3.0/include'
#DEPENDPATH += $$PWD/'../../Program Files (x86)/AMD APP SDK/3.0/include'
