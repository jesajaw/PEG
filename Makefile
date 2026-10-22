# Library paths: Uncomment as required:
#####################
VPATH=src
INCLUDEPATH += /usr/include/eigen3
#LIBPATH=/home/user/lib

CC=mpic++
CFLAGS= -c -O3 -march=native -g -Wall -fopenmp -MMD -MP -Isrc $(addprefix -I,$(INCLUDEPATH))
LDFLAGS=-fopenmp $(if $(LIBPATH),-L$(LIBPATH))
SOURCES = mainSerial.cpp mainMPI.cpp PEG.cpp TESolver.cpp TMSolver.cpp mainSupport.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: pegSerial pegMPI

pegSerial: mainSerial.o $(COMMON)
	$(CC) $^ $(LDFLAGS) -o $@

pegMPI: mainMPI.o $(COMMON)
	$(CC) $^ $(LDFLAGS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(OBJECTS:.o=.d) pegSerial pegMPI

.PHONY: all clean

-include $(OBJECTS:.o=.d)