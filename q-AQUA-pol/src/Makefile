#Generic makefile for fortran90
FC := ifort 
# serial
#LIBS := -qmkl=sequential
#FFLAGS := -r8 -O3
# parallel  
LIBS := -mkl=sequential
FFLAGS := -r8 -O3 -qopenmp

OBJ := \
constants.o\
model_mod.o\
ttm2f_mod.o\
ttm3f_mod.o\
nasa_mod.o\
math.o\
smear.o\
nasa.o\
potential_mod.o\
bemsa2b.o\
bemsa3b.o\
bemsa4b.o\
pes_shell.o\

OBJPOT := \
getpot.o


%.o : %.f90
	$(FC) -c $(FFLAGS) $< 

getpot.x :  $(OBJ) $(OBJPOT)
	$(FC) -o $@  $^ $(FFLAGS) $(LIBS)

clean :
	rm *.o *.mod 
