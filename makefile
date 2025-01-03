# compiler name (gcc seems to use smarter tricks since the executable is faster when using gcc instead of clang on an arch linux machine)
CC = gcc

# Cplex Location
CPLEX_HEADERS_COMPILER_FLAG := -I/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include/ilcplex -I/opt/concorde

SRC_DIR := src/
HEADERS_DIR := $(SRC_DIR)headers/

LDFLAGS1 = -L/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/ -L/opt/concorde
LDFLAGS2 = -lcplex -lm -l:concorde.a

# release build is default
OBJ_DIR = obj/release/
BIN_DIR = bin/release/
CFLAGS = -O3 -mfma -msse4 -mavx2 -march=native -mtune=native -Isrc/headers
ifeq ($(MODE),debug)
OBJ_DIR = obj/debug/
BIN_DIR = bin/debug/
CFLAGS = -Wall -g -mfma -msse4 -mavx2 -march=native -mtune=native -Isrc/headers
endif

# avx default
# PREPROCDEF = -D COMPUTATION_TYPE=0
# ifeq ($(COMPUTATION),base)
# PREPROCDEF = -D COMPUTATION_TYPE=1
# endif
# ifeq ($(COMPUTATION),matrix)
# PREPROCDEF = -D COMPUTATION_TYPE=2
# endif

# separate files into the ones that use cplex and will need the extra compiler flag to find cplex headers and the ones which don't use it
SOURCE_NAMES_NO_CPLEX = TspUtilities.c ArgParser.c TspIOUtils.c CostMatrix.c NearestNeighbor.c ExtraMileage.c 2Opt.c 2OptMultithread.c 3Opt.c 3OptMultithread.c Tabu.c VariableNeighborhood.c Genetic.c SimulatedAnnealing.c
SOURCE_NAMES_CPLEX = main.c TspCplex.c Benders.c PatchingHeuristic.c BranchAndCut.c HardFixing.c LocalBranching.c

HEADER_NAMES_NO_CPLEX = TspBase.h Tsp.h EdgeCostFunctions.h
HEADER_NAMES_CPLEX = TspCplex.h

# files list
HEADER_FILES_NO_CPLEX := $(HEADER_NAMES_NO_CPLEX:%=$(HEADERS_DIR)%)
HEADER_FILES_CPLEX := $(HEADER_NAMES_CPLEX:%=$(HEADERS_DIR)%)

OBJ_FILES_NO_CPLEX := $(SOURCE_NAMES_NO_CPLEX:%.c=$(OBJ_DIR)%.o)
OBJ_FILES_CPLEX := $(SOURCE_NAMES_CPLEX:%.c=$(OBJ_DIR)%.o)
OBJ_FILES := $(OBJ_FILES_NO_CPLEX) $(OBJ_FILES_CPLEX)

# command used to check variables value and "debug" the makefile
print:
	@echo HEADER_FILES = $(HEADER_FILES_NO_CPLEX)$(HEADER_FILES_NO_CPLEX)
	@echo SOURCE_NAMES_NO_CPLEX = $(SOURCE_NAMES_NO_CPLEX:%=$(SRC_DIR)%)
	@echo SOURCE_NAMES_CPLEX = $(SOURCE_NAMES_CPLEX:%=$(SRC_DIR)%)
	@echo OBJ_FILES = $(OBJ_FILES)

# define dependencies for c files that use cplex(so TspCplex.h needs to be added as a dependecy)
$(OBJ_FILES_CPLEX): $(HEADER_FILES_CPLEX)

# build options when debugging
build: $(BIN_DIR)main

$(BIN_DIR)main: $(OBJ_FILES)
	$(CC) $(CFLAGS) $(LDFLAGS1) $(OBJ_FILES) -o $(BIN_DIR)main $(LDFLAGS2)

$(OBJ_DIR)%.o: $(SRC_DIR)%.c $(HEADER_FILES_NO_CPLEX)
	$(CC) $(PREPROCDEF) -c $(CFLAGS) $(CPLEX_HEADERS_COMPILER_FLAG) $(SRC_DIR)$(*F).c -o $@

SRC_FILES_PATH := $(SOURCE_NAMES_NO_CPLEX:%=$(SRC_DIR)%) $(SOURCE_NAMES_CPLEX:%=$(SRC_DIR)%)

final:
	$(CC) -D COMPUTATION_TYPE=0 -O3 -ftree-loop-im -mfma -msse4 -mavx2 -march=native -mtune=native -Isrc/headers $(LDFLAGS1) $(CPLEX_HEADERS_COMPILER_FLAG) $(SRC_FILES_PATH) -o bin/release/main $(LDFLAGS2)

# delete all gcc output files
clean:
	rm -f bin/debug/main bin/release/main
	rm -f obj/debug/*.o obj/release/*.o
