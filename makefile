# compiler name
CC = gcc

# compiler/liker flags
CFLAGS_DEBUG = -Wall -g -c -mavx2
LDFLAGS_DEBUG =
CFLAGS_EXEC = -O3 -c -mavx2
LDFLAGS_EXEC =

# directories names
OBJ_DEBUG_DIR = obj/debug
OBJ_EXEC_DIR = obj/x64
BIN_DEBUG_DIR = bin/debug
BIN_EXEC_DIR = bin/x64

# files list
OBJ_DEBUG_FILES = $(OBJ_DEBUG_DIR)/main.o $(OBJ_DEBUG_DIR)/utilities.o $(OBJ_DEBUG_DIR)/tsp.o
OBJ_EXEC_FILES = $(OBJ_EXEC_DIR)/main.o $(OBJ_EXEC_DIR)/utilities.o $(OBJ_EXEC_DIR)/tsp.o

# list of header files
HEADERS = src/tsp.h


# build options when debugging
buildDebug: $(OBJ_DEBUG_FILES)
	$(CC) $(LDFLAGS_DEBUG) $(OBJ_DEBUG_FILES) -o $(BIN_DEBUG_DIR)/main

$(OBJ_DEBUG_DIR)/main.o: src/main.c $(HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/main.c -o $(OBJ_DEBUG_DIR)/main.o

$(OBJ_DEBUG_DIR)/utilities.o: src/utilities.c $(HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/utilities.c -o $(OBJ_DEBUG_DIR)/utilities.o

$(OBJ_DEBUG_DIR)/tsp.o: src/tsp.c $(HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/tsp.c -o $(OBJ_DEBUG_DIR)/tsp.o


# build options when testing performance
buildExec: $(OBJ_EXEC_FILES)
	$(CC) $(LDFLAGS_EXEC) $(OBJ_EXEC_FILES) -o $(BIN_EXEC_DIR)/main

$(OBJ_EXEC_DIR)/main.o: src/main.c $(HEADERS)
	$(CC) $(CFLAGS_EXEC) src/main.c -o $(OBJ_EXEC_DIR)/main.o

$(OBJ_EXEC_DIR)/utilities.o: src/utilities.c $(HEADERS)
	$(CC) $(CFLAGS_EXEC) src/utilities.c -o $(OBJ_EXEC_DIR)/utilities.o

$(OBJ_EXEC_DIR)/tsp.o: src/tsp.c $(HEADERS)
	$(CC) $(CFLAGS_EXEC) src/tsp.c -o $(OBJ_EXEC_DIR)/tsp.o


build: buildDebug buildExec

# delete all gcc output files
clean:
	rm -f $(BIN_DEBUG_DIR)/main $(BIN_EXEC_DIR)/main
	rm -f $(OBJ_DEBUG_DIR)/*.o $(OBJ_EXEC_DIR)/*.o