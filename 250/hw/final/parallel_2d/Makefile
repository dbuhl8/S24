FC = mpif90

OBJECTS = pgameoflife.o pdriver.o 
	
MODULES = pgameoflife.mod

OUTPUTS = output.txt, *.dat

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g -llapack -lblas #-fopenmp

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
