QT     -= gui core

TARGET = pegSerial

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS += -O3 -fopenmp -DNDEBUG
QMAKE_LFLAGS += -fopenmp

# Homebrew include paths -- Apple Silicon default shown; Intel Macs (brew installs under /usr/local) should use the commented-out line instead.
INCLUDEPATH += /opt/homebrew/include /opt/homebrew/include/eigen3
# INCLUDEPATH += /usr/local/include /usr/local/include/eigen3

HEADERS += src/PEG.h \
	src/TESolver.h \
	src/TMSolver.h \
	src/mainSupport.h

SOURCES += src/PEG.cpp\
	src/TESolver.cpp \
	src/TMSolver.cpp \
	src/mainSupport.cpp \
	src/mainSerial.cpp
