#-------------------------------------------------
#
# Project created by QtCreator 2017-03-01T11:02:07
#
#-------------------------------------------------

QT       += testlib

QT       -= gui

TARGET = tst_krigingtest
CONFIG   += console
CONFIG   -= app_bundle
CONFIG += c++11

TEMPLATE = app


SOURCES += tst_krigingtest.cpp \
    krig.cpp \
    powvargram.cpp

HEADERS += \
    krig.h \
    powvargram.h

DEFINES += SRCDIR=\\\"$$PWD/\\\"
