QT-=gui core

TARGET=pegSerial

QMAKE_CXXFLAGS += -O3 -march=native -fopenmp
QMAKE_LFLAGS += -fopenmp

INCLUDEPATH += /usr/include/eigen3

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

DESTDIR = $$PWD
OBJECTS_DIR = $$PWD/build/qt

equals(OUT_PWD, $$_PRO_FILE_PWD_) {
    error("I won't run qmake in the main folder, because that would overwrite the Makefile. Please use: mkdir -p build/qt && qmake -o build/qt/Makefile PEG_ubuntu.pro (Qt Creator does this on its own with a shadow build).")
}

# removes object files & finished programs
QMAKE_CLEAN += $$PWD/pegSerial $$PWD/pegMPI

# additionally removes the whole build folder.
tidy.commands = rm -rf $$PWD/build
QMAKE_EXTRA_TARGETS += tidy