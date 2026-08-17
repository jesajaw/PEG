QT     -= gui core

TARGET = pegSerial

QMAKE_CXXFLAGS += -march=native -O2 -c -g -Wall -fopenmp
QMAKE_LFLAGS += -fopenmp

# Eigen and Boost.Odeint are header-only -- no LIBS or rpath needed.
# On Ubuntu, Eigen typically installs to /usr/include/eigen3 and needs an explicit INCLUDEPATH:
#INCLUDEPATH += /usr/include/eigen3

HEADERS += src/PEG.h \
	src/TESolver.h \
	src/mainSupport.h

SOURCES += src/PEG.cpp\
	src/TESolver.cpp \
	src/mainSupport.cpp \
	src/mainSerial.cpp