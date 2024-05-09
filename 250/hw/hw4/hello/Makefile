FC = mpif90
OBJ = helloworld.o
FFLAGS = 

helloworld: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o hello_mpi

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJ) hello_mpi


