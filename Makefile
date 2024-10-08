# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -O2
LIBS = -lopenblas -lgsl -lgslcblas -lm -g3

# Files
OBJ = timer.o spmv.o my_dense.o my_sparse.o
EXEC = spmv

# Default rule to build the program
all: $(EXEC)

# Rule to create the executable
$(EXEC): $(OBJ)
	$(CC) $(OBJ) $(LIBS) -o $(EXEC)

# Rule to compile timer.c
timer.o: timer.c
	@$(CC) $(CFLAGS) -c timer.c

# Rule to compile spmv.c
spmv.o: spmv.c
	@$(CC) $(CFLAGS) -c spmv.c

# Rule to compile my_dense.c
my_dense.o: my_dense.c
	@$(CC) $(CFLAGS) -c my_dense.c

# Rule to compile my_sparse.c
my_sparse.o: my_sparse.c
	@$(CC) $(CFLAGS) -c my_sparse.c

# Clean rule to remove object files and executable
clean:
	rm -f $(OBJ) $(EXEC)

# PHONY to avoid issues when file named 'clean' exists
.PHONY: clean