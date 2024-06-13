FC = mpif90

OBJECTS = pgameoflife.o pdriver.o 
	
MODULES = pgameoflife.mod

OUTPUTS = output.txt, *.dat

FFLAGS = -fbounds-check -g -cpp -DHB
# optional flags are -DHB to make the code work on hummingbird and my laptop
# This makes the code work using send recv rather than isendrecv which only some versions of mpi have
# Other flag is -DDECOMP_1D or -DDECOMP_2D which tells the code to do a one d or two d decomposition

.PHONY: clean

parallel: driver.exe

%.exe: $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(FFLAGS) -o driver.exe

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(MODULES) driver.exe 
