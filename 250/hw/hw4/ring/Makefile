FC = mpif90
OBJ = ring.o
FFLAGS = 

ring: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o ring_mpi

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJ) ring_mpi


