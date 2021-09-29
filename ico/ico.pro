TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += xdrfile-1.1.4/include/

QMAKE_CXXFLAGS += -O2 -march=native -std=c++11 -Wno-unused-parameter -Wno-sign-compare

SOURCES += main.cpp \
    xdrfile-1.1.4/src/xdrfile.c \
    xdrfile-1.1.4/src/xdrfile_c_test.c \
    xdrfile-1.1.4/src/xdrfile_trr.c \
    xdrfile-1.1.4/src/xdrfile_xtc.c

HEADERS += \
    particle.h \
    icosahedron.h \
    oblatespheroid.h \
    sphere.h \
    tennisball.h \
    spherepatch.h \
    pentamer.h \
    dodecahedron.h \
    surface.h \
    chain.h \
    data.h \
    slab.h \
    atom.h \
    cow.h \
    rng.h \
    xdrfile-1.1.4/include/xdrfile.h \
    xdrfile-1.1.4/include/xdrfile_trr.h \
    xdrfile-1.1.4/include/xdrfile_xtc.h \
    atom.h \
    chain.h \
    cow.h \
    data.h \
    dodecahedron.h \
    icosahedron.h \
    oblatespheroid.h \
    particle.h \
    pentamer.h \
    rng.h \
    slab.h \
    sphere.h \
    spherepatch.h \
    surface.h \
    tennisball.h \
    xtcanalysis.h \
    welford.h

