# compiler name
CC = gcc

# compiler/liker flags
CFLAGS_DEBUG = -Wall -g -c
LDFLAGS_DEBUG =
CFLAGS_EXEC = -O3 -c
LDFLAGS_EXEC =

# directories names
OBJ_DEBUG = obj/debug
OBJ_EXEC = obj/x64
BIN_DEBUG = bin/debug
BIN_EXEC = bin/x64

# list of header files
HEADERS = src/utilities.h src/tsp.h

# build options when debugging
buildDebug: main.o_DEBUG utilities.o_DEBUG
	$(CC) $(LDFLAGS_DEBUG) $(OBJ_DEBUG)/main.o $(OBJ_DEBUG)/utilities.o -o $(BIN_DEBUG)/main

main.o_DEBUG: src/main.c $(HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/main.c -o $(OBJ_DEBUG)/main.o

utilities.o_DEBUG: src/utilities.c src/utilities.h
	$(CC) $(CFLAGS_DEBUG) src/utilities.c -o $(OBJ_DEBUG)/utilities.o


# build options when testing performance
buildExec: main.o_EXEC utilities.o_EXEC
	$(CC) $(CFLAGS_EXEC) $(OBJ_EXEC)/main.o $(OBJ_EXEC)/utilities.o -o $(BIN_EXEC)/main

main.o_EXEC: src/main.c $(HEADERS)
	$(CC) $(CFLAGS_EXEC) src/main.c -o $(OBJ_EXEC)/main.o

utilities.o_EXEC: src/utilities.c src/utilities.h
	$(CC) $(CFLAGS_EXEC) src/utilities.c -o $(OBJ_EXEC)/utilities.o


# delete all gcc output files
clean:
	rm -f $(BIN_DEBUG)/* $(BIN_EXEC)/*
	rm -f $(OBJ_DEBUG)/*.o $(OBJ_EXEC)/*.o