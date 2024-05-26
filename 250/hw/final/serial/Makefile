FC = gfortran

OBJECTS = gameoflife.o driver.o 
	
MODULES = gameoflife.mod

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g -llapack -lblas #-fopenmp

.PHONY: clean

serial: driver.exe
	./driver.exe <IC.dat

%.exe: $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(FFLAGS) -o driver.exe

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(MODULES) driver.exe 
