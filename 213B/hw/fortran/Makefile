OBJECTS = LinAl.o NumDE.o driver.o
	
MODULES = LinAl.mod NumDE.mod

OUTPUTS = *.png *.dat *.mod output.txt

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g

.PHONY: clean

output.txt: driver.exe
	./driver.exe > output.txt
	gnuplot plot.p

driver.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o driver.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(OUTPUTS) $(MODULES) driver.exe 
