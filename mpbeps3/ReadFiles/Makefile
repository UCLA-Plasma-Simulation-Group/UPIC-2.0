#Makefile for 3D MPI/OpenMP PIC Periodic Field file readers

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran -fopenmp
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#LOPTS = -lgomp
#LEGACY =

# Makefile Intel compiler with Mac OS X

#FC90 = ifort -qopenmp
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

#FC90 = ifort -qopenmp
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99

#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

FC90 = gfortran -fopenmp
CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

LOPTS = -lgomp
#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
LEGACY =

# Makefile PGI compiler with Linux

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#LEGACY = -Mnomain

# Makefile Cray compiler with Linux

#FC90 = ftn -h omp
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#LEGACY =

#

# Linkage rules

all : preadf3 pvreadf3

preadf3 : preadf3.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) $(LGOPTS) -o preadf3 \
    preadf3.o input3mod.o $(GOBJS)

pvreadf3 : pvreadf3.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) $(LGOPTS) -o pvreadf3 \
    pvreadf3.o input3mod.o $(GOBJS)

# Compilation rules

#OPENMP
#omplib.o : omplib.f90
#	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

input3mod.o : input3mod.f90
	$(FC90) $(OPTS90) -o input3mod.o -c input3mod.f90

preadf3.o : preadf3.f90 input3mod.o
	$(FC90) $(OPTS90) -o preadf3.o -c preadf3.f90

pvreadf3.o : pvreadf3.f90 input3mod.o
	$(FC90) $(OPTS90) -o pvreadf3.o -c pvreadf3.f90

cpexecs : cpreadf3 cpvreadf3

cpreadf3 : preadf3
	cp preadf3 ../

cpvreadf3 : pvreadf3
	cp pvreadf3 ../

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f preadf3 pvreadf3
