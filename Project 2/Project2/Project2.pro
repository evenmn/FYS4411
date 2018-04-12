TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    wavefunction.cpp \
    monte_carlo.cpp \
    gradient_descent.cpp

HEADERS += \
    wavefunction.h \
    monte_carlo.h \
    gradient_descent.h
