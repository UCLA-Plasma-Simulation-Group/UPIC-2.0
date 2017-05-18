#Makefile for 2D MPI/OpenMP PIC Periodic Field file readers

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

# Fortran Graphics

#Tektronix
#GOBJS = libt1.o plot10.o libplt10.o libgks1.o libgks2.o modgraf2.o
# Postcript printer
#GOBJS = libpsp.o libgks1.o libgks2.o modgraf2.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libygl.o libgks1.o libgks2.o modgraf2.o
#No graphics
GOBJS = nullgks1.o nullgks2.o modgraf2.o

# Linkage rules

all : preadf2 pvreadf2

preadf2 : preadf2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) $(LGOPTS) -o preadf2 \
    preadf2.o input2mod.o $(GOBJS)

pvreadf2 : pvreadf2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) $(LGOPTS) -o pvreadf2 \
    pvreadf2.o input2mod.o $(GOBJS)

# Compilation rules

#OPENMP
#omplib.o : omplib.f90
#	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

input2mod.o : input2mod.f90
	$(FC90) $(OPTS90) -o input2mod.o -c input2mod.f90

preadf2.o : preadf2.f90 input2mod.o modgraf2.o
	$(FC90) $(OPTS90) -o preadf2.o -c preadf2.f90

pvreadf2.o : pvreadf2.f90 input2mod.o modgraf2.o
	$(FC90) $(OPTS90) -o pvreadf2.o -c pvreadf2.f90

# Fortran Graphics

libt1.o : sgks2/libt1.f90
	$(FC90) $(OPTS90) -o libt1.o -c sgks2/libt1.f90

plot10.o : sgks2/plot10.f
	$(FC90) $(OPTS90) -o plot10.o -c sgks2/plot10.f

libplt10.o : sgks2/libplt10.f
	$(FC90) $(OPTS90) -o libplt10.o -c sgks2/libplt10.f

libpsp.o : sgks2/libpsp.f
	$(FC90) $(OPTS90) -o libpsp.o -c sgks2/libpsp.f

libygl.o : sgks2/libygl.f
	$(FC90) $(OPTS90) -o libygl.o -c sgks2/libygl.f

libgks1.o : sgks2/libgks1.f
	$(FC90) $(OPTS90) -o libgks1.o -c sgks2/libgks1.f

nullgks1.o : sgks2/nullgks1.f
	$(FC90) $(OPTS90) -o nullgks1.o -c sgks2/nullgks1.f

libgks2.o : sgks2/libgks2.f
	$(FC90) $(OPTS90) -o libgks2.o -c sgks2/libgks2.f

nullgks2.o : sgks2/nullgks2.f
	$(FC90) $(OPTS90) -o nullgks2.o -c sgks2/nullgks2.f

modgraf2.o : sgks2/modgraf2.f90
	$(FC90) $(OPTS90) -o modgraf2.o -c sgks2/modgraf2.f90

cpexecs : cpreadf2 cpvreadf2

cpreadf2 : preadf2
	cp preadf2 ../

cpvreadf2 : pvreadf2
	cp pvreadf2 ../

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f preadf2 pvreadf2
