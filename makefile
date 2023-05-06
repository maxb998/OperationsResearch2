# compiler name
CC = clang

# Cplex Location
CPLEX_HEADERS_COMPILER_FLAG = -I/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include/ilcplex

SRC_DIR := src/
HEADERS_DIR := src/headers/

# build debug as default
OBJ_DIR = obj/debug/
BIN_DIR = bin/debug/
CFLAGS = -Wall -g -c -mavx2 -Isrc/headers
LDFLAGS1 = -L/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/
LDFLAGS2 = -lcplex -lm

# condition to check value passed
ifeq ($(MODE),exec)
OBJ_DIR = obj/x64/
BIN_DIR = bin/x64/
CFLAGS = -O3 -c -mavx2 -march=native -Isrc/headers
LDFLAGS1 = -L/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/
LDFLAGS2 = -lcplex -lm
endif

# separate files into the ones that use cplex and will need the extra compiler flag to find cplex headers and the ones which don't use it
SOURCE_NAMES := TspBase.c TspUtilities.c ArgParser.c TspIOUtils.c CostMatrix.c NearestNeighbor.c ExtraMileage.c 2Opt.c VariableNeighborhood.c TspCplex.c Blenders.c 
#SOURCE_NAMES_CPLEX := TspCplex.c Blenders.c 

HEADER_NAMES = TspBase.h TspFunctions.h EdgeCostFunctions.h TspCplex.h
#HEADER_NAMES_CPLEX = TspCplex.h

# files list
HEADER_FILES := $(HEADER_NAMES:%=$(HEADERS_DIR)%)
#HEADER_FILES_CPLEX := $(HEADER_NAMES_CPLEX:%=$(HEADERS_DIR)%)

SRC_FILES := $(SOURCE_NAMES:%=$(SRC_DIR)%)
#SRC_FILES_CPLEX := $(SOURCE_NAMES_CPLEX:%=$(SRC_DIR)%)

OBJ_FILES := $(SOURCE_NAMES:%.c=$(OBJ_DIR)%.o)
#OBJ_FILES_CPLEX := $(SOURCE_NAMES_CPLEX:%.c=$(OBJ_DIR)%.o)

# command used to check variables value and "debug" the makefile
print:
	echo $(HEADER_FILES)
	echo $(HEADER_FILES_CPLEX)
	echo $(SRC_FILES)
	echo $(SRC_FILES_CPLEX)
	echo $(OBJ_FILES)
	echo $(OBJ_FILES_CPLEX)

# build options when debugging
build: $(BIN_DIR)main

$(BIN_DIR)main: $(OBJ_DIR)main.o $(OBJ_FILES) $(OBJ_FILES_CPLEX)
	$(CC) $(LDFLAGS1) $(OBJ_DIR)main.o $(OBJ_FILES) $(OBJ_FILES_CPLEX) -o $(BIN_DIR)main $(LDFLAGS2)

$(OBJ_DIR)main.o: $(SRC_DIR)main.c $(HEADER_FILES) $(HEADER_FILES_CPLEX)
	$(CC) $(CFLAGS) $(CPLEX_HEADERS_COMPILER_FLAG) $(SRC_DIR)main.c -o $(OBJ_DIR)main.o

$(OBJ_DIR)%.o: $(SRC_DIR)%.c $(HEADER_FILES)
	$(CC) $(CFLAGS) $(CPLEX_HEADERS_COMPILER_FLAG) $(SRC_DIR)$(*F).c -o $@

#$(OBJ_FILES_CPLEX): $(SRC_DIR)$(*F).c $(HEADER_FILES) $(HEADER_FILES_CPLEX)
#	$(CC) $(CFLAGS) $(CPLEX_HEADERS_COMPILER_FLAG) $(SRC_DIR)$(*F).c -o $@

# delete all gcc output files
clean:
	rm -f bin/debug/main bin/x64/main
	rm -f obj/debug/*.o obj/x64/*.o
