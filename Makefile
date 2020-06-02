# Compiler options
CC=g++
FLAGS_CC = -Wall -O0

# -----------------------------------------------------------------------------------------
# Variables
LIBS_HEAD = $(wildcard src/solvers/*.h) $(wildcard src/sparse/*.h) $(wildcard src/utils/*h)
LIBS_SRC  = $(wildcard src/solvers/*.cpp) $(wildcard src/sparse/*.cpp) $(wildcard src/utils/*cpp)
LIBS_OBJ  = $(LIBS_SRC:.cpp=.o)
LIBS_INCL = $(wildcard src)

LIBS      = upa.a

# -----------------------------------------------------------------------------------------
# UPA
all: $(LIBS)
	$(CC) tests/test_directSolvers.cpp $(FLAGS_CC) -o tests/test_directSolvers $(LIBS) -I$(LIBS_INCL)
	$(CC) tests/test_utils.cpp $(FLAGS_CC) -o tests/test_utils $(LIBS) -I$(LIBS_INCL)
	$(CC) tests/test_sparse.cpp $(FLAGS_CC) -o tests/test_sparse $(LIBS) -I$(LIBS_INCL)
	$(CC) upa.cc $(FLAGS_CC) -o upa $(LIBS) -I$(LIBS_INCL)

# -----------------------------------------------------------------------------------------
# Numerical libraries
$(LIBS): $(LIBS_OBJ)
	rm -f $(LIBS)
	ar rv $(LIBS) $(LIBS_OBJ)

$(LIBS_OBJ) : $(LIBS_SRC) $(LIBS_HEAD)

%.o: %.cpp $(LIBS_HEAD)
	$(CC) $(FLAGS_CC) -c -o $@ $< -I$(LIBS_INCL)

clean :
	rm src/solvers/*.o src/sparse/*.o src/utils/*.o upa.a
