QT     -= gui core

TARGET = pegSerial

QMAKE_CXXFLAGS += -march=native -O2 -c -g -Wall -fopenmp
QMAKE_LFLAGS += -fopenmp

INCLUDEPATH += /usr/include/eigen3

HEADERS += src/PEG.h \
	src/TESolver.h \
	src/mainSupport.h

SOURCES += src/PEG.cpp\
	src/TESolver.cpp \
	src/mainSupport.cpp \
	src/mainSerial.cpp