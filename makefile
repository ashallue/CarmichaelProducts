# Makefile for compiling test

# Compiler and flags.  Note that -I. tells compiler to start in current directory for relative paths
CXX = g++
CXXFLAGS = -O3 -lgmp -lgmpxx -I.
MPI = mpic++
PATHS = -I/usr/local/include -L/usr/local/lib

# Target executables
TARGETS = test

# Source files
SRCS = is_prime_pl.cpp unit_tests.cpp CondensedInteger.cpp ModElement.cpp many_carmichaels.cpp gen_div.cpp

# Source files compiled with mpic++

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target (build all)
all: $(TARGETS)

# Generic rule for compiling .cpp to .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule for compiling test
test:  test_main.cpp is_prime_pl.o CondensedInteger.o ModElement.o unit_tests.o many_carmichaels.o gen_div.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Clean up object files and executables
clean:
	rm -f $(OBJS) $(TARGETS)

# .PHONY to indicate these are not real files
.PHONY: all clean
