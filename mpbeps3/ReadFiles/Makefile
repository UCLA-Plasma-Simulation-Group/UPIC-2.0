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

#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

FC90 = gfortran -fopenmp
CC = gcc

PFC = f2py --fcompiler=gnu95
PCC = f2py --compiler=unix

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#Python
PYOPTS90 = $(OPTS90) -fPIC
PYCCOPTS = $(CCOPTS) -fPIC

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

LOPTS = -lgomp
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

all : preadf3 pvreadf3 pflreadf3 pfvreadf3 ptreadf3 psreadf3

preadf3 : preadf3.o libcmfield3.o nulllib3.o
	$(FC90) $(OPTS90) $(LOPTS) -o preadf3 \
    preadf3.o input3mod.o cmfield3.o libcmfield3_h.o libcmfield3.o nulllib3.o

pvreadf3 : pvreadf3.o libcmfield3.o nulllib3.o
	$(FC90) $(OPTS90) $(LOPTS) -o pvreadf3 \
    pvreadf3.o input3mod.o cmfield3.o libcmfield3_h.o libcmfield3.o nulllib3.o

pflreadf3 : pflreadf3.o libcmfield3.o nulllib3.o
	$(FC90) $(OPTS90) $(LOPTS) -o pflreadf3 \
    pflreadf3.o input3mod.o cmfield3.o libcmfield3_h.o libcmfield3.o nulllib3.o

pfvreadf3 : pfvreadf3.o libcmfield3.o nulllib3.o
	$(FC90) $(OPTS90) $(LOPTS) -o pfvreadf3 \
    pfvreadf3.o input3mod.o cmfield3.o libcmfield3_h.o libcmfield3.o nulllib3.o

ptreadf3 : ptreadf3.o libcmfield3.o nulllib3.o
	$(FC90) $(OPTS90) $(LOPTS) -o ptreadf3 \
    ptreadf3.o input3mod.o cmfield3.o libcmfield3_h.o libcmfield3.o nulllib3.o

psreadf3 : psreadf3.o libcmfield3.o nulllib3.o
	$(FC90) $(OPTS90) $(LOPTS) -o psreadf3 \
    psreadf3.o input3mod.o cmfield3.o libcmfield3_h.o libcmfield3.o nulllib3.o

list:
	echo preadf3 pvreadf3 pflreadf3 pfvreadf3 ptreadf3 psreadf3 python \
	     clean clobber

# Python

python : libcmfield3.so

# Compilation rules

#OPENMP
#omplib.o : omplib.f90
#	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

input3mod.o : input3mod.f90
	$(FC90) $(OPTS90) -o input3mod.o -c input3mod.f90

preadf3.o : preadf3.f90 input3mod.o cmfield3.o
	$(FC90) $(OPTS90) -o preadf3.o -c preadf3.f90

pvreadf3.o : pvreadf3.f90 input3mod.o cmfield3.o
	$(FC90) $(OPTS90) -o pvreadf3.o -c pvreadf3.f90

pflreadf3.o : pflreadf3.f90 input3mod.o cmfield3.o
	$(FC90) $(OPTS90) -o pflreadf3.o -c pflreadf3.f90

pfvreadf3.o : pfvreadf3.f90 input3mod.o cmfield3.o
	$(FC90) $(OPTS90) -o pfvreadf3.o -c pfvreadf3.f90

ptreadf3.o : ptreadf3.f90 input3mod.o cmfield3.o
	$(FC90) $(OPTS90) -o ptreadf3.o -c ptreadf3.f90

psreadf3.o : psreadf3.f90 input3mod.o cmfield3.o
	$(FC90) $(OPTS90) -o psreadf3.o -c psreadf3.f90

libcmfield3.o : libcmfield3.f
	$(FC90) $(OPTS90) -o libcmfield3.o -c libcmfield3.f

libcmfield3_h.o : libcmfield3_h.f90
	$(FC90) $(OPTS90) -o libcmfield3_h.o -c libcmfield3_h.f90

cmfield3.o : cmfield3.f90 libcmfield3_h.o input3mod.o
	$(FC90) $(OPTS90) -o cmfield3.o -c cmfield3.f90

nulllib3.o : nulllib3.f90
	$(FC90) $(OPTS90) -o nulllib3.o -c nulllib3.f90

# Python

pymcmfield3.o : libcmfield3.f
	$(FC90) $(PYOPTS90) -o pymcmfield3.o -c libcmfield3.f

libcmfield3.so : input3mod.f90 libcmfield3_h.o cmfield3.f90 pymcmfield3.o \
                 nulllib3.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libcmfield3 \
         -c input3mod.f90 libcmfield3_h.o cmfield3.f90 pymcmfield3.o \
         nulllib3.o

cpexecs : cpreadf3 cpvreadf3 cpflreadf3 cpfvreadf3 cptreadf3 cpsreadf3

cpreadf3 : preadf3
	cp preadf3 ../

cpvreadf3 : pvreadf3
	cp pvreadf3 ../
	
cpflreadf3 : pflreadf3
	cp pflreadf3 ../

cpfvreadf3 : pfvreadf3
	cp pfvreadf3 ../

cptreadf3 : ptreadf3
	cp ptreadf3 ../

cpsreadf3 : psreadf3
	cp psreadf3 ../

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f preadf3 pvreadf3 pflreadf3 pfvreadf3 ptreadf3 psreadf3 *.so

#Python
cpyexecs : cpyreadf3 cpyvreadf3 cpyflreadf3 cpyfvreadf3 cpytreadf3 \
           cpysreadf3

cpyreadf3 : preadf3.py
	cp preadf3.py ../

cpyvreadf3 : pvreadf3.py
	cp pvreadf3.py ../
	
cpyflreadf3 : pflreadf3.py
	cp pflreadf3.py ../

cpyfvreadf3 : pfvreadf3.py
	cp pfvreadf3.py ../

cpytreadf3 : ptreadf3.py
	cp ptreadf3.py ../

cpysreadf3 : psreadf3.py
	cp psreadf3.py ../
