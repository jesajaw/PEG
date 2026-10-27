SRCDIR  := src
OBJDIR  := build/obj
THIRD   := third_party

EIGEN_INC := $(THIRD)/eigen
BOOST_INC := $(THIRD)/boost
ifneq ($(MAKECMDGOALS),clean)
ifeq ($(wildcard $(EIGEN_INC)/Eigen/Dense),)
$(error Eigen headers not found: $(EIGEN_INC)/Eigen/Dense is missing (see third_party/LICENSES))
endif
ifeq ($(wildcard $(BOOST_INC)/boost/numeric/odeint.hpp),)
$(error Boost.Odeint headers not found: $(BOOST_INC)/boost/numeric/odeint.hpp is missing (see third_party/LICENSES))
endif
endif

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

COMMON     := PEG TESolver TMSolver mainSupport SolverSupport
COMMON_OBJ := $(addprefix $(OBJDIR)/,$(addsuffix .o,$(COMMON)))

.DEFAULT_GOAL := pegSerial

pegSerial: $(OBJDIR)/mainSerial.o $(COMMON_OBJ)
	$(CXX) $^ $(LDFLAGS) -o $@

pegMPI: $(OBJDIR)/mainMPI.o $(COMMON_OBJ)
	$(MPICXX) $^ $(LDFLAGS) -o $@

all: pegSerial pegMPI

$(OBJDIR)/mainMPI.o: $(SRCDIR)/mainMPI.cpp | $(OBJDIR)
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $@

CHECK_ARGS := --mode constantIncidence --min 100 --max 110 --increment 5 --eV --incidenceAngle 87 --N 8 --gratingType sinusoidal --gratingMaterial Au --gratingPeriod 1.2  --gratingGeometry 0.05 --TETM
MPIEXEC ?= mpiexec
MPI_NP  ?= 2

check:
	@mkdir -p build; found=0; \
	if [ -x pegSerial ]; then found=1; \
	  PEG_MATERIAL_DB=materialDatabase ./pegSerial $(CHECK_ARGS) --outputFile build/check.txt > /dev/null \
	  && test "$$(grep -cE '^[0-9]' build/check.txt)" = 6 \
	  && echo "check OK: pegSerial" \
	  || { echo "check FAILED: pegSerial"; exit 1; }; \
	fi; \
	if [ -x pegMPI ]; then found=1; \
	  PEG_MATERIAL_DB=materialDatabase $(MPIEXEC) -n $(MPI_NP) ./pegMPI $(CHECK_ARGS) --outputFile build/checkMPI.txt > /dev/null \
	  && test "$$(grep -cE '^[0-9]' build/checkMPI.txt)" = 6 \
	  && echo "check OK: pegMPI with $(MPI_NP) processes" \
	  || { echo "check FAILED: pegMPI"; exit 1; }; \
	fi; \
	test $$found = 1 || { echo "check: no valid installation (pegSerial / pegMPI not found) -- run 'make pegSerial' or 'make pegMPI' first"; exit 1; }

clean:
	rm -rf build pegSerial pegMPI

-include $(wildcard $(OBJDIR)/*.d)
