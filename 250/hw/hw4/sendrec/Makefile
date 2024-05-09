FC = mpif90
OBJ = sendrec.o
FFLAGS = 

sendrec: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o sendrec_mpi

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJ) sendrec_mpi


