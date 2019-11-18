TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXX += -fopenmp



SOURCES += \
        main.cpp

INCLUDEPATH += C:\armadillo-9.800.2\include
DEPENDPATH += C:\armadillo-9.800.2\include


LIBS += \
    -LC:\armadillo-9.800.2\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT \
    -fopenmp
