FC = mpif90
OBJ = ping.o
FFLAGS = 

ping: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o ping_mpi

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJ) ping_mpi


