FC = mpif90
OBJ = pi.o
FFLAGS = 

pi: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o pi_mpi

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJ) pi_mpi


