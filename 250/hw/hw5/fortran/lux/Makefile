FC = gfortran
OBJ = LinAl.o hw5.o
FFLAGS = -fopenmp -fdefault-real-8 -fdefault-double-8 -fbounds-check -g

helloworld: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o hw5

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJ) hw5


