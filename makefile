# compiler name
CC = gcc

# compiler/liker flags
CFLAGS_DEBUG = -Wall -g -c
LDFLAGS_DEBUG =
CFLAGS_EXEC = -O3
LDFLAGS_EXEC =

# directories names
OBJ_DEBUG = obj/debug
OBJ_EXEC = obj/x64
BIN_DEBUG = bin/debug
BIN_EXEC = bin/x64

# list of header files
HEADERS = src/logLib.h src/tsp.h

# build options when debugging
buildDebug: main.o_DEBUG logLib.o_DEBUG
	$(CC) $(LDFLAGS_DEBUG) $(OBJ_DEBUG)/main.o $(OBJ_DEBUG)/logLib.o -o $(BIN_DEBUG)/main

main.o_DEBUG: src/main.c $(HEADERS)
	$(CC) $(CFLAGS_DEBUG) src/main.c -o $(OBJ_DEBUG)/main.o

logLib.o_DEBUG: src/logLib.c src/logLib.h
	$(CC) $(CFLAGS_DEBUG) src/logLib.c -o $(OBJ_DEBUG)/logLib.o


# build options when testing performance
buildExec: main.o_EXEC logLib.o_EXEC
	$(CC) $(CFLAGS_EXEC) $(OBJ_EXEC)/main.o $(OBJ_EXEC)/logLib.o -o $(BIN_EXEC)/main

main.o_EXEC: src/main.c $(HEADERS)
	$(CC) $(CFLAGS_EXEC) src/main.c -o $(OBJ_EXEC)/main.o

logLib.o_EXEC: src/logLib.c src/logLib.h
	$(CC) $(CFLAGS_EXEC) src/logLib.c -o $(OBJ_EXEC)/logLib.o


# delete all gcc output files
clean:
	rm -f $(BIN_DEBUG)/* $(BIN_EXEC)/*
	rm -f $(OBJ_DEBUG)/*.o $(OBJ_EXEC)/*.o