
# ----------------------------------------------------------
# Compiler options
COMP_CPP = g++
COMP_F90 = gfortran

FLAGS_CPP = -Wall -O0
FLAGS_F90 = -O0 -fdefault-real-8

# ----------------------------------------------------------
# Source files 
TEST_CPP = test_utils.cpp
SRC_CPP  = $(filter-out $(TEST_CPP), $(wildcard *.cpp) $(wildcard ../solvers/linearSolvers.cpp))
HEAD_CPP = $(wildcard *.h) $(wildcard ../solvers/linearSolvers.h)
OBJS_CPP = $(SRC_CPP:.cpp=.o)


# -----------------------------------------------------------
# Compilation
all : $(OBJS_CPP)
	$(COMP_CPP) $(TEST_CPP) $(FLAGS_CPP) -o test_cpp $(HEAD_CPP) $(OBJS_CPP)


$(OBJS_CPP) : $(SRC_CPP) $(HEAD_CPP)

%.o: %.cpp $(HEAD_CPP)
	$(COMP_CPP) $(FLAGS_CPP) -c -o $@ $<

clean : 
	rm *.o ../solvers/*.o