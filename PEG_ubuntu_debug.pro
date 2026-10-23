# qmake project for pegSerial (Linux, debug build)
QT -= gui core
CONFIG += console c++17 release
CONFIG -= app_bundle debug_and_release debug_and_release_target
TARGET = pegSerial
DESTDIR = $$PWD

INCLUDEPATH += $$PWD/src
exists($$PWD/third_party/eigen/Eigen/Dense): INCLUDEPATH += $$PWD/third_party/eigen
else: INCLUDEPATH += /usr/include/eigen3
exists($$PWD/third_party/boost/boost/numeric/odeint.hpp): INCLUDEPATH += $$PWD/third_party/boost

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS   += -fopenmp
QMAKE_CXXFLAGS += -O0 -g -Wall -Wextra
contains(CONFIG, native): QMAKE_CXXFLAGS += -march=native

HEADERS += src/PEG.h \
	src/TESolver.h \
	src/TMSolver.h \
	src/mainSupport.h \
	src/sharedSolverSupport.h

SOURCES += src/PEG.cpp\
	src/TESolver.cpp \
	src/TMSolver.cpp \
	src/mainSupport.cpp \
	src/sharedSolverSupport.cpp \
	src/mainSerial.cpp

QMAKE_CLEAN += $$DESTDIR/pegSerial
tidy.commands = $(MAKE) clean && $(DEL_FILE) Makefile .qmake.stash
QMAKE_EXTRA_TARGETS += tidy
