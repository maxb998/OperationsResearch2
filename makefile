# compiler name
CC = gcc

# compiler/liker flags
CFLAGS_DEBUG = -Wall -g -c -mavx2 -Isrc/headers
LDFLAGS_DEBUG = -lm
CFLAGS_EXEC = -O3 -c -mavx2 -march=native -Isrc/headers
LDFLAGS_EXEC = -lm

# directories names
OBJ_DEBUG_DIR = obj/debug
OBJ_EXEC_DIR = obj/x64
BIN_DEBUG_DIR = bin/debug
BIN_EXEC_DIR = bin/x64

# files list
OBJ_DEBUG_FILES = $(OBJ_DEBUG_DIR)/main.o $(OBJ_DEBUG_DIR)/TspBase.o $(OBJ_DEBUG_DIR)/TspUtilities.o $(OBJ_DEBUG_DIR)/TspIOUtils.o $(OBJ_DEBUG_DIR)/CostMatrix.o $(OBJ_DEBUG_DIR)/NearestNeighbor.o $(OBJ_DEBUG_DIR)/ExtraMileage.o $(OBJ_DEBUG_DIR)/2Opt.o $(OBJ_DEBUG_DIR)/VariableNeighborhood.o $(OBJ_DEBUG_DIR)/ArgParser.o
OBJ_EXEC_FILES = $(OBJ_EXEC_DIR)/main.o $(OBJ_EXEC_DIR)/TspBase.o $(OBJ_EXEC_DIR)/TspUtilities.o $(OBJ_EXEC_DIR)/TspIOUtils.o $(OBJ_EXEC_DIR)/CostMatrix.o $(OBJ_EXEC_DIR)/NearestNeighbor.o $(OBJ_EXEC_DIR)/ExtraMileage.o $(OBJ_EXEC_DIR)/2Opt.o $(OBJ_EXEC_DIR)/VariableNeighborhood.o $(OBJ_EXEC_DIR)/ArgParser.o

# list of header files
GLOBAL_HEADERS = src/headers/TspBase.h src/headers/EdgeCostFunctions.h
SPECIFIC_HEADERS = src/headers/CostMatrix.h src/headers/NearestNeighbor.h src/headers/ExtraMileage.h src/headers/2Opt.h src/headers/VariableNeighborhood.h src/headers/ArgParser.h

# build options when debugging
buildDebug: $(OBJ_DEBUG_FILES)
	$(CC) $(OBJ_DEBUG_FILES) -o $(BIN_DEBUG_DIR)/main $(LDFLAGS_DEBUG)

$(OBJ_DEBUG_DIR)/main.o: src/main.c $(GLOBAL_HEADERS) $(SPECIFIC_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/main.c -o $(OBJ_DEBUG_DIR)/main.o

$(OBJ_DEBUG_DIR)/TspBase.o: src/TspBase.c src/headers/TspBase.h
	$(CC) $(CFLAGS_DEBUG) src/TspBase.c -o $(OBJ_DEBUG_DIR)/TspBase.o

$(OBJ_DEBUG_DIR)/TspUtilities.o: src/TspUtilities.c src/headers/TspUtilities.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/TspUtilities.c -o $(OBJ_DEBUG_DIR)/TspUtilities.o

$(OBJ_DEBUG_DIR)/TspIOUtils.o: src/TspIOUtils.c src/headers/TspIOUtils.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/TspIOUtils.c -o $(OBJ_DEBUG_DIR)/TspIOUtils.o

$(OBJ_DEBUG_DIR)/CostMatrix.o: src/CostMatrix.c src/headers/CostMatrix.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/CostMatrix.c -o $(OBJ_DEBUG_DIR)/CostMatrix.o

$(OBJ_DEBUG_DIR)/NearestNeighbor.o: src/NearestNeighbor.c src/headers/NearestNeighbor.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/NearestNeighbor.c -o $(OBJ_DEBUG_DIR)/NearestNeighbor.o

$(OBJ_DEBUG_DIR)/ExtraMileage.o: src/ExtraMileage.c src/headers/ExtraMileage.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/ExtraMileage.c -o $(OBJ_DEBUG_DIR)/ExtraMileage.o

$(OBJ_DEBUG_DIR)/2Opt.o: src/2Opt.c src/headers/2Opt.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/2Opt.c -o $(OBJ_DEBUG_DIR)/2Opt.o

$(OBJ_DEBUG_DIR)/VariableNeighborhood.o: src/VariableNeighborhood.c src/headers/VariableNeighborhood.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/VariableNeighborhood.c -o $(OBJ_DEBUG_DIR)/VariableNeighborhood.o

$(OBJ_DEBUG_DIR)/ArgParser.o: src/ArgParser.c src/headers/ArgParser.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/ArgParser.c -o $(OBJ_DEBUG_DIR)/ArgParser.o

# build options when testing performance
buildExec: $(OBJ_EXEC_FILES)
	$(CC) $(OBJ_EXEC_FILES) -o $(BIN_EXEC_DIR)/main $(LDFLAGS_EXEC)

$(OBJ_EXEC_DIR)/main.o: src/main.c $(GLOBAL_HEADERS) $(SPECIFIC_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/main.c -o $(OBJ_EXEC_DIR)/main.o

$(OBJ_EXEC_DIR)/TspBase.o: src/TspBase.c src/headers/TspBase.h
	$(CC) $(CFLAGS_EXEC) src/TspBase.c -o $(OBJ_EXEC_DIR)/TspBase.o

$(OBJ_EXEC_DIR)/TspUtilities.o: src/TspUtilities.c src/headers/TspUtilities.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/TspUtilities.c -o $(OBJ_EXEC_DIR)/TspUtilities.o

$(OBJ_EXEC_DIR)/TspIOUtils.o: src/TspIOUtils.c src/headers/TspIOUtils.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/TspIOUtils.c -o $(OBJ_EXEC_DIR)/TspIOUtils.o

$(OBJ_EXEC_DIR)/CostMatrix.o: src/CostMatrix.c src/headers/CostMatrix.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/CostMatrix.c -o $(OBJ_EXEC_DIR)/CostMatrix.o

$(OBJ_EXEC_DIR)/NearestNeighbor.o: src/NearestNeighbor.c  src/headers/NearestNeighbor.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/NearestNeighbor.c -o $(OBJ_EXEC_DIR)/NearestNeighbor.o

$(OBJ_EXEC_DIR)/ExtraMileage.o: src/ExtraMileage.c src/headers/ExtraMileage.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/ExtraMileage.c -o $(OBJ_EXEC_DIR)/ExtraMileage.o

$(OBJ_EXEC_DIR)/2Opt.o: src/2Opt.c src/headers/2Opt.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/2Opt.c -o $(OBJ_EXEC_DIR)/2Opt.o

$(OBJ_EXEC_DIR)/VariableNeighborhood.o: src/VariableNeighborhood.c src/headers/VariableNeighborhood.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/VariableNeighborhood.c -o $(OBJ_EXEC_DIR)/VariableNeighborhood.o

$(OBJ_EXEC_DIR)/ArgParser.o: src/ArgParser.c src/headers/ArgParser.h $(GLOBAL_HEADERS)
	$(CC) $(CFLAGS_EXEC) src/ArgParser.c -o $(OBJ_EXEC_DIR)/ArgParser.o


build: buildDebug buildExec

# delete all gcc output files
clean:
	rm -f $(BIN_DEBUG_DIR)/main $(BIN_EXEC_DIR)/main
	rm -f $(OBJ_DEBUG_DIR)/*.o $(OBJ_EXEC_DIR)/*.o