TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
TARGET = lbm_mpi
SOURCES += main.cpp \
    math_api.cpp

HEADERS += \
    math_api.h
INCLUDEPATH += "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
LIBS += "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\msmpi.lib"
QMAKE_CXXFLAGS+= -openmp
QMAKE_LFLAGS +=  -openmp
