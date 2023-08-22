#Updated version
#Name of the project
PROJ_NAME=pyHPLiBEla

#Sources
SOURCE=$(wildcard ../src/*.cpp)

#Objects
OBJ=$(SOURCE:.cpp=.o)

#Compiler
CC=g++

#Flags for compiler
CC_FLAGS=-fPIC                       \
         -O3                         \
	 -I/usr/include/python3.10 
#Linkage
LD_FLAGS=-L/usr/lib/x86_64-linux-gnu \
         -lboost_python310           \
	 -lgsl                       \
	 -lgslcblas                  \
	 -lm 


#Command to clean target
RM = rm -rf

#Compilation and linking

all: pyLattice.so pyMonteCarlo.so pyThermo.so
	@ echo ' '


pyLattice.so: Lattice.o
	@echo 'Building Lattice dynamic library.'
	$(CC) -shared Lattice.o -o pyLattice.so $(LD_FLAGS)
	@echo ' '

pyMonteCarlo.so: MonteCarlo.o Lattice.o
	@echo 'Building MC dynamic library.'
	$(CC) -shared Lattice.o MonteCarlo.o -o pyMonteCarlo.so $(LD_FLAGS)
	@echo ' '
	
pyThermo.so : Thermo.o
	@echo 'Building Thermo dynamic library.'
	$(CC) -shared Thermo.o Lattice.o -o pyThermo.so $(LD_FLAGS)
	@echo ' '
	
%.o: ../src/%.cpp
	$(CC) $(CC_FLAGS) $< -c -o $@
	@ echo ''

clean:
	@ $(RM) obj src test $(PROJ_NAME) *~
	@ $(RM) *.so *.o *.h *.cpp

.PHONY: all clean

