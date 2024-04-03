FC = gfortran

OBJECTS = LinAl.o hw2mod.o driver.o
	
MODULES = LinAl.mod hw2mod.mod

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g #-llapack -lblas 

.PHONY: clean

driver.exe: $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(FFLAGS) -o driver.exe

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(MODULES) driver.exe 
