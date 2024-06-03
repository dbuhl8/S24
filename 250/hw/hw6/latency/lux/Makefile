FC = mpif90
OBJ = latency.o
FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g

latency: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o latency.exe

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJ) latency.exe


