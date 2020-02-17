
# ----------------------------------------------------------
# Compiler options
COMP_CPP = g++
COMP_F90 = gfortran

FLAGS_CPP = -Wall -O0
FLAGS_F90 = -O0 -fdefault-real-8

# ----------------------------------------------------------
# Source files 
TEST_CPP = test_linearSolvers.cpp
SRC_CPP  = $(filter-out $(TEST_CPP), $(wildcard *.cpp))
HEAD_CPP = $(wildcard *.h)
OBJS_CPP = $(SRC_CPP:.cpp=.o)

TEST_F90 = test_linearSolvers.f90
SRC_F90  = $(filter-out $(TEST_F90), $(wildcard *.f90))
OBJS_F90 = $(SRC_F90:.f90=.o)
MOD_F90  = $(SRC_F90:.f90=.mod)




# -----------------------------------------------------------
# Compilation
all : $(OBJS_CPP) $(MOD_F90)
	$(COMP_CPP) $(TEST_CPP) $(FLAGS_CPP) -o test_cpp $(HEAD_CPP) $(OBJS_CPP)
	$(COMP_F90) $(TEST_F90) $(FLAGS_F90) -o test_f90 $(MOD_F90)


$(OBJS_CPP) : $(SRC_CPP) $(HEAD_CPP)
$(OBJS_F90) : $(SRC_F90)

%.o: %.cpp $(HEAD_CPP)
	$(COMP_CPP) $(FLAGS_CPP) -c -o $@ $<

%.o: %.f90
	$(COMP_F90) $(FLAGS_F90) -c -o $@ $<

%.mod: %.f90
	$(COMP_F90) $(FLAGS_F90) -c -o $@ $<

clean : 
	rm *.o *.mod

