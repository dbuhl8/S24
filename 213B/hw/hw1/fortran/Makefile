OBJECTS = LinAl.o NumDE.o driver.o
	
MODULES = LinAl.mod NumDE.mod

FFLAGS = -fdefault-real-8 -fdefault-double-8 -llapack -lblas -fbounds-check -g

.PHONY: clean

output.txt: driver.exe
	./driver.exe > output.txt

driver.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o driver.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(DOBJECTS) $(MODULES) driver.exe 
