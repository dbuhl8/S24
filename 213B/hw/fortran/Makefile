OBJECTS = LinAl.o NumDE.o hw2_driver.o # hw1_driver.o
	
MODULES = LinAl.mod NumDE.mod

OUTPUTS = *.png *.dat *.mod output.txt

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g

.PHONY: clean

hw2: driver.exe
	./driver.exe > output.txt
	gnuplot hw2_plot.p

hw1: driver.exe
	./driver.exe > output.txt
	gnuplot hw1_plot.p

%.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o driver.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(OUTPUTS) $(MODULES) driver.exe 
