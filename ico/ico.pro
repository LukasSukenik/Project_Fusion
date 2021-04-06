TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -O2 -march=native -std=c++11 -Wno-unused-parameter -Wno-sign-compare

SOURCES += main.cpp

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
    rng.h

