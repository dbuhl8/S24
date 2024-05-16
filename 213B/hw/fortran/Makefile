OBJECTS = LinAl.o NumDE.o hw3_driver.o # Make sure the correct driver is in the objects variable
	
MODULES = LinAl.mod NumDE.mod

OUTPUTS = *.png *.mod output.txt #you can optionally include *.dat, but this might accidentally delete data you want to keep if you aren't careful

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g

.PHONY: clean

hw3: driver.exe
	./driver.exe > output.txt
	gnuplot hw3_plot.p

midterm: driver.exe
	./driver.exe > output.txt
	gnuplot midterm_plot.p

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
