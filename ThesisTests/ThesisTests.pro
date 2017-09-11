#-------------------------------------------------
#
# Project created by QtCreator 2017-09-10T14:44:32
#
#-------------------------------------------------

QT       += testlib
QT += core
QT += gui

TARGET = tst_thesisteststest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
        tst_thesisteststest.cpp \ 
    ../bessel.cpp \
    ../clbessel.cpp \
    ../mom.cpp \
    ../types.cpp

DEFINES += SRCDIR=\\\"$$PWD/\\\"

DISTFILES += \
    ../besselj0.cl

HEADERS += \
    ../bessel.h \
    ../clbessel.h \
    ../helper.h \
    ../mom.h \
    ../types.h

macx: LIBS += -L/usr/local/lib/ -laf.3.5.0
macx: LIBS += -framework OpenCL

INCLUDEPATH += /usr/local/include
DEPENDPATH += /usr/local/include
