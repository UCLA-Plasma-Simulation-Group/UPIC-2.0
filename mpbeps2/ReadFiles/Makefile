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
#LGOPTS = sgks2/Ygl-4.2/libYgl.a -L/usr/X11R6/lib -lX11
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
#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
#LGOPTS = sgks2/Ygl-4.2/libYgl.a -L/usr/X11R6/lib -lX11
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
#modgraf1.o
# Postcript printer
#GOBJS = libpsp.o libgks1.o libgks2.o modgraf2.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libygl.o libgks1.o libgks2.o libgraf2_h.o modgraf2.o modgraf1.o
#No graphics
GOBJS = nullgks1.o nullgks2.o modgraf2.o modgraf1.o

# Python Graphics

#Tektronix
#GMODSRC = sgks2/modgraf2.f90 sgks2/modgraf1.f90 sgks2/libt1.f90
#GPYOBJS = pymlibgks1.o pymlibgks2.o pymlibplt10.o pymplot10.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GMODSRC = sgks2/modgraf2.f90 sgks2/modgraf1.f90
#GPYOBJS = pymlibgks1.o pymlibgks2.o pymlibygl.o
#No graphics
GMODSRC = sgks2/modgraf2.f90 sgks2/modgraf1.f90
GPYOBJS = pymnullgks1.o pymnullgks2.o

# Linkage rules

all : preadf2 pvreadf2 pflreadf2 pfvreadf2 ptreadf2 psreadf2

preadf2 : preadf2.o libcmfield2.o nulllib2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o preadf2 \
    preadf2.o input2mod.o cmfield2.o libcmfield2_h.o libcmfield2.o \
    nulllib2.o $(GOBJS) $(LGOPTS)

pvreadf2 : pvreadf2.o libcmfield2.o nulllib2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pvreadf2 \
    pvreadf2.o input2mod.o cmfield2.o libcmfield2_h.o libcmfield2.o \
    nulllib2.o $(GOBJS) $(LGOPTS)

pflreadf2 : pflreadf2.o libcmfield2.o nulllib2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pflreadf2 \
    pflreadf2.o input2mod.o cmfield2.o libcmfield2_h.o libcmfield2.o \
    nulllib2.o $(GOBJS) $(LGOPTS)

pfvreadf2 : pfvreadf2.o libcmfield2.o nulllib2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pfvreadf2 \
    pfvreadf2.o input2mod.o cmfield2.o libcmfield2_h.o libcmfield2.o \
    nulllib2.o $(GOBJS) $(LGOPTS)

ptreadf2 : ptreadf2.o libcmfield2.o nulllib2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o ptreadf2 \
    ptreadf2.o input2mod.o cmfield2.o libcmfield2_h.o libcmfield2.o \
    nulllib2.o $(GOBJS) $(LGOPTS)

psreadf2 : psreadf2.o libcmfield2.o nulllib2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o psreadf2 \
    psreadf2.o input2mod.o cmfield2.o libcmfield2_h.o libcmfield2.o \
    nulllib2.o $(GOBJS) $(LGOPTS)

pcreadf2 : pcreadf2.o libcmfield2.o nulllib2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pcreadf2 \
    pcreadf2.o input2mod.o cmfield2.o libcmfield2_h.o libcmfield2.o \
    nulllib2.o $(GOBJS) $(LGOPTS)

pcvreadf2 : pcvreadf2.o libcmfield2.o nulllib2.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pcvreadf2 \
    pcvreadf2.o input2mod.o cmfield2.o libcmfield2_h.o libcmfield2.o \
    nulllib2.o $(GOBJS) $(LGOPTS)
    
list:
	echo preadf2 pvreadf2 pflreadf2 pfvreadf2 ptreadf2 psreadf2 python \
	     clean clobber

# Python

python : libcmfield2.so fgraf.so

# Compilation rules

#OPENMP
#omplib.o : omplib.f90
#	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

input2mod.o : input2mod.f90
	$(FC90) $(OPTS90) -o input2mod.o -c input2mod.f90

preadf2.o : preadf2.f90 input2mod.o cmfield2.o modgraf2.o
	$(FC90) $(OPTS90) -o preadf2.o -c preadf2.f90

pvreadf2.o : pvreadf2.f90 input2mod.o cmfield2.o modgraf2.o
	$(FC90) $(OPTS90) -o pvreadf2.o -c pvreadf2.f90

pflreadf2.o : pflreadf2.f90 input2mod.o cmfield2.o modgraf2.o
	$(FC90) $(OPTS90) -o pflreadf2.o -c pflreadf2.f90

pfvreadf2.o : pfvreadf2.f90 input2mod.o cmfield2.o modgraf2.o modgraf1.o
	$(FC90) $(OPTS90) -o pfvreadf2.o -c pfvreadf2.f90

ptreadf2.o : ptreadf2.f90 input2mod.o cmfield2.o modgraf2.o modgraf1.o
	$(FC90) $(OPTS90) -o ptreadf2.o -c ptreadf2.f90

psreadf2.o : psreadf2.f90 input2mod.o cmfield2.o modgraf2.o modgraf1.o
	$(FC90) $(OPTS90) -o psreadf2.o -c psreadf2.f90

pcreadf2.o : pcreadf2.f90 input2mod.o cmfield2.o modgraf2.o
	$(FC90) $(OPTS90) -o pcreadf2.o -c pcreadf2.f90

pcvreadf2.o : pcvreadf2.f90 input2mod.o cmfield2.o modgraf2.o
	$(FC90) $(OPTS90) -o pcvreadf2.o -c pcvreadf2.f90

libcmfield2.o : libcmfield2.f
	$(FC90) $(OPTS90) -o libcmfield2.o -c libcmfield2.f

libcmfield2_h.o : libcmfield2_h.f90
	$(FC90) $(OPTS90) -o libcmfield2_h.o -c libcmfield2_h.f90

cmfield2.o : cmfield2.f90 libcmfield2_h.o input2mod.o
	$(FC90) $(OPTS90) -o cmfield2.o -c cmfield2.f90

nulllib2.o : nulllib2.f90
	$(FC90) $(OPTS90) -o nulllib2.o -c nulllib2.f90

# Python

pymcmfield2.o : libcmfield2.f
	$(FC90) $(PYOPTS90) -o pymcmfield2.o -c libcmfield2.f

libcmfield2.so : input2mod.f90 libcmfield2_h.o cmfield2.f90 pymcmfield2.o \
                 nulllib2.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libcmfield2 \
         -c input2mod.f90 libcmfield2_h.o cmfield2.f90 pymcmfield2.o \
         nulllib2.o

# Fortran Graphics

libt1.o : sgks2/libt1.f90
	$(FC90) $(OPTS90) -o libt1.o -c sgks2/libt1.f90

plot10.o : sgks2/plot10.f
	$(FC90) $(OPTS90) -o plot10.o -c sgks2/plot10.f

libplt10.o : sgks2/libplt10.f
	$(FC90) $(OPTS90) -o libplt10.o -c sgks2/libplt10.f

libpsp.o : sgks2/libpsp.f
	$(FC90) $(OPTS90) -o libpsp.o -c sgks2/libpsp.f

sgks2/Ygl-4.2/libYgl.a :
	make -C sgks2/Ygl-4.2

libygl.o : sgks2/libygl.f sgks2/Ygl-4.2/libYgl.a
	$(FC90) $(OPTS90) -o libygl.o -c sgks2/libygl.f

libgks1.o : sgks2/libgks1.f
	$(FC90) $(OPTS90) -o libgks1.o -c sgks2/libgks1.f

nullgks1.o : sgks2/nullgks1.f
	$(FC90) $(OPTS90) -o nullgks1.o -c sgks2/nullgks1.f

libgks2.o : sgks2/libgks2.f
	$(FC90) $(OPTS90) -o libgks2.o -c sgks2/libgks2.f

nullgks2.o : sgks2/nullgks2.f
	$(FC90) $(OPTS90) -o nullgks2.o -c sgks2/nullgks2.f

libgraf2_h.o : sgks2/libgraf2_h.f90
	$(FC90) $(OPTS90) -o libgraf2_h.o -c sgks2/libgraf2_h.f90

modgraf2.o : sgks2/modgraf2.f90 libgraf2_h.o
	$(FC90) $(OPTS90) -o modgraf2.o -c sgks2/modgraf2.f90

modgraf1.o : sgks2/modgraf1.f90
	$(FC90) $(OPTS90) -o modgraf1.o -c sgks2/modgraf1.f90

# Python Graphics

pymplot10.o : sgks2/plot10.f
	$(FC90) $(PYOPTS90) -o pymplot10.o -c sgks2/plot10.f

pymlibplt10.o : sgks2/libplt10.f
	$(FC90) $(PYOPTS90) -o pymlibplt10.o -c sgks2/libplt10.f

pymlibygl.o : sgks2/libygl.f sgks2/Ygl-4.2/libYgl.a
	$(FC90) $(PYOPTS90) -o pymlibygl.o -c sgks2/libygl.f

pymlibpsp.o : sgks2/libpsp.f
	$(FC90) $(PYOPTS90) -o pymlibpsp.o -c sgks2/libpsp.f

pymlibgks1.o : sgks2/libgks1.f
	$(FC90) $(PYOPTS90) -o pymlibgks1.o -c sgks2/libgks1.f

pymnullgks1.o : sgks2/nullgks1.f
	$(FC90) $(PYOPTS90) -o pymnullgks1.o -c sgks2/nullgks1.f

pymlibgks2.o : sgks2/libgks2.f
	$(FC90) $(PYOPTS90) -o pymlibgks2.o -c sgks2/libgks2.f

pymnullgks2.o : sgks2/nullgks2.f
	$(FC90) $(PYOPTS90) -o pymnullgks2.o -c sgks2/nullgks2.f

fgraf.so : libgraf2_h.o $(GMODSRC) $(GPYOBJS)
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) $(LGOPTS) -m fgraf \
        -c libgraf2_h.o $(GMODSRC) $(GPYOBJS)

cpexecs : cpreadf2 cpvreadf2 cpflreadf2 cpfvreadf2 cptreadf2 cpsreadf2

cpreadf2 : preadf2
	cp preadf2 ../

cpvreadf2 : pvreadf2
	cp pvreadf2 ../

cpflreadf2 : pflreadf2
	cp pflreadf2 ../

cpfvreadf2 : pfvreadf2
	cp pfvreadf2 ../

cptreadf2 : ptreadf2
	cp ptreadf2 ../

cpsreadf2 : psreadf2
	cp psreadf2 ../

cpcreadf2 : pcreadf2
	cp pcreadf2 ../

cpcvreadf2 : pcvreadf2
	cp pcvreadf2 ../

clean :
	rm -f *.o *.mod sgks2/Ygl-4.2/*.o

clobber: clean
	rm -f preadf2 pvreadf2 pflreadf2 pfvreadf2 ptreadf2 psreadf2 \
	      pcreadf2 pcvreadf2 *.so sgks2/Ygl-4.2/*.a

#Python
cpyexecs : cpyreadf2 cpyvreadf2 cpyflreadf2 cpyfvreadf2 cpytreadf2 \
           cpysreadf2

cpyreadf2 : preadf2.py
	cp preadf2.py ../

cpyvreadf2 : pvreadf2.py
	cp pvreadf2.py ../
	
cpyflreadf2 : pflreadf2.py
	cp pflreadf2.py ../

cpyfvreadf2 : pfvreadf2.py
	cp pfvreadf2.py ../

cpytreadf2 : ptreadf2.py
	cp ptreadf2.py ../

cpysreadf2 : psreadf2.py
	cp psreadf2.py ../
