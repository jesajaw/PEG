SRCDIR  := src
OBJDIR  := build/obj
THIRD   := third_party

EIGEN_INC ?= $(if $(wildcard $(THIRD)/eigen/Eigen/Dense),$(THIRD)/eigen,/usr/include/eigen3)
BOOST_INC ?= $(if $(wildcard $(THIRD)/boost/boost/numeric/odeint.hpp),$(THIRD)/boost)

CXX     ?= g++
MPICXX  ?= mpic++
NATIVE  ?= 0
DEBUG   ?= 0

CXXFLAGS += -std=c++17 -Wall -fopenmp -MMD -MP -I$(SRCDIR) -isystem $(EIGEN_INC) $(addprefix -isystem ,$(BOOST_INC))
ifeq ($(DEBUG),1)
CXXFLAGS += -O0 -g
else
CXXFLAGS += -O3 -DNDEBUG
ifeq ($(NATIVE),1)
CXXFLAGS += -march=native
endif
endif
LDFLAGS  += -fopenmp

COMMON     := PEG TESolver TMSolver mainSupport sharedSolverSupport
COMMON_OBJ := $(addprefix $(OBJDIR)/,$(addsuffix .o,$(COMMON)))

.DEFAULT_GOAL := pegSerial

pegSerial: $(OBJDIR)/mainSerial.o $(COMMON_OBJ)
	$(CXX) $^ $(LDFLAGS) -o $@

pegMPI: $(OBJDIR)/mainMPI.o $(COMMON_OBJ)
	$(MPICXX) $^ $(LDFLAGS) -o $@

mpi: pegMPI
all: pegSerial pegMPI

# only this file includes mpi.h
$(OBJDIR)/mainMPI.o: $(SRCDIR)/mainMPI.cpp | $(OBJDIR)
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $@

check: pegSerial
	PEG_MATERIAL_DB=materialDatabase ./pegSerial --mode constantIncidence --min 100 --max 110 --increment 5 --eV \
	  --incidenceAngle 87 --N 8 --gratingType sinusoidal --gratingMaterial Au --gratingPeriod 1.2 \
	  --gratingGeometry 0.05 --computeTE --computeTM --outputFile build/check.txt > /dev/null
	@test "$$(grep -cE '^[0-9]' build/check.txt)" = 6 && echo "check OK (6 result rows: TE+TM x 3 energies)" \
	  || { echo "check FAILED"; exit 1; }

clean:
	rm -rf build pegSerial pegMPI

.PHONY: all mpi check clean
-include $(wildcard $(OBJDIR)/*.d)
