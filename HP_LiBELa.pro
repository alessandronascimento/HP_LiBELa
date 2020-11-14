TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        trunk/src/Lattice.cpp \
        trunk/src/MonteCarlo.cpp \
        trunk/src/Thermo.cpp \
        trunk/src/main.cpp

HEADERS += \
    trunk/src/Lattice.h \
    trunk/src/MonteCarlo.h \
    trunk/src/Thermo.h

INCLUDEPATH += ../../LiBELa/trunk/lib/gsl/include

LIBS += -L../../LiBELa/trunk/lib/gsl/lib/ \
    -lgsl -lgslcblas -lm
