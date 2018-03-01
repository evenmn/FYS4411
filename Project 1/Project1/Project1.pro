TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    wavefunction.cpp \
    metropolis.cpp \
    gd.cpp \
    tools.cpp

HEADERS += \
    wavefunction.h \
    metropolis.h \
    gd.h \
    tools.h

# remove possible other optimization flags
#QMAKE_CXXFLAGS_RELEASE -= -O
#QMAKE_CXXFLAGS_RELEASE -= -O1
#QMAKE_CXXFLAGS_RELEASE -= -O2

# add the desired -O3 if not present
#QMAKE_CXXFLAGS_RELEASE *= -O3
