F90 = gfortran
#INC=-I/usr/include
#LIBS=-L/usr/lib -lnetcdff -lnetcdf
INC =
LIBS =
FFLAGS =

# OBJS = Constants.o Quadrature_Points_and_Weights_Class.o Lagrange_Interpolation_Class.o NodalDG_Class.o NodalDG_main.o
# OBJS = Constants.o Quadrature_Points_and_Weights_Class.o Lagrange_Interpolation_Class.o NodalDG_Class.o NodalDG_2D_Class.o NodalDG_2Dmain.o
OBJS = Constants.o Quadrature_Points_and_Weights_Class.o Lagrange_Interpolation_Class.o NodalDG_Class.o Common_Class.o DGSEM1D_Class.o DGSEM1D_main.o

wall -fbounds-check : $(OBJS)
# wall : $(OBJS)
	$(F90) $(INC) -o test.exe $(OBJS) $(LIBS)

%.o : %.f90 *.f90
	$(F90) $(INC) -c $< 

clean :	
	rm *.o *.mod *.exe

