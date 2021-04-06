TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt


INCLUDEPATH += xdrfile-1.1.4/include

QMAKE_CXXFLAGS += -Ofast -march=native -std=c++11

SOURCES += main.cpp \
    xdrfile-1.1.4/src/xdrfile.c \
    xdrfile-1.1.4/src/xdrfile_c_test.c \
    xdrfile-1.1.4/src/xdrfile_trr.c \
    xdrfile-1.1.4/src/xdrfile_xtc.c

HEADERS += \
    membrane.h \
    flatmembrane.h \
    vesicle.h \
    xdrfile-1.1.4/include/xdrfile.h \
    xdrfile-1.1.4/include/xdrfile_trr.h \
    xdrfile-1.1.4/include/xdrfile_xtc.h \
    input.h

