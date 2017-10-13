QT += core
QT += gui

CONFIG += c++11

TARGET = Thesis
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    bessel.cpp \
    clbessel.cpp \
    mom.cpp \
    types.cpp \
    tissueproperties.cpp \
    phantomfile.cpp

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
    mom.h \
    types.h \
    tissueproperties.h \
    phantomfile.h

macx: QMAKE_MAC_SDK = macosx10.13
macx: LIBS += -L$$PWD/../../../../../usr/local/lib/ -laf.3.5.0
macx: LIBS += -framework OpenCL

INCLUDEPATH += /usr/local/include
DEPENDPATH += /usr/local/include




DISTFILES += \
    besselj0.cl



#win32: LIBS += -L$$PWD/'../../Program Files/ArrayFire/v3/lib/' -laf

#INCLUDEPATH += $$PWD/'../../Program Files/ArrayFire/v3/include'
#DEPENDPATH += $$PWD/'../../Program Files/ArrayFire/v3/include'

#win32: LIBS += -L$$PWD/'../../Program Files (x86)/AMD APP SDK/3.0/lib/x86_64/' -lOpenCL

#INCLUDEPATH += $$PWD/'../../Program Files (x86)/AMD APP SDK/3.0/include'
#DEPENDPATH += $$PWD/'../../Program Files (x86)/AMD APP SDK/3.0/include'
