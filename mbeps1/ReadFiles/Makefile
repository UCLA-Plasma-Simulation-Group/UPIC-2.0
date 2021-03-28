#Makefile for 1D OpenMP PIC Periodic Field file readers

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran -fopenmp
#CC = gcc

#PFC = f2py --fcompiler=gnu95
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#Python
#PYOPTS90 = $(OPTS90) -fPIC
#PYCCOPTS = $(CCOPTS) -fPIC
#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

#LOPTS = -lgomp
#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
#LGOPTS = sgks1/Ygl-4.2/libYgl.a -L/usr/X11R6/lib -lX11
#LEGACY = -std=legacy

# Makefile Intel compiler with Mac OS X

#FC90 = ifort -qopenmp
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#NOMAIN = -nofor_main

# Makefile Intel compiler with Linux

#FC90 = ifort -qopenmp
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99

#NOMAIN = -nofor_main

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
#LGOPTS = sgks1/Ygl-4.2/libYgl.a -L/usr/X11R6/lib -lX11
LEGACY = -std=legacy

# Makefile PGI compiler with Linux

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#NOMAIN = -Mnomain

# Makefile Cray compiler with Linux

#FC90 = ftn -h omp
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#

# Fortran Graphics

#Tektronix
#GOBJS = libt1.o plot10.o libplt10.o libgks1.o dlibgks1.o libgks2.o \
#libgraf1_h.o modgraf1.o libgraf2_h.o modgraf2.o
# Postcript printer
#GOBJS = libpsp.o libgks1.o dlibgks1.o libgks2.o libgraf1_h.o modgraf1.o \
#libgraf2_h.o modgraf2.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libygl.o libgks1.o dlibgks1.o libgks2.o libgraf1_h.o modgraf1.o \
#libgraf2_h.o modgraf2.o
#No graphics
GOBJS = nullgks1.o nullgks2.o modgraf1.o modgraf2.o

# Python Graphics

#Tektronix
#GMODSRC = sgks1/modgraf1.f90 sgks1/modgraf2.f90 sgks1/libt1.f90
#GPYOBJS = pymlibgks1.o pymdlibgks1.o pymlibgks2.o pymlibplt10.o pymplot10.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GMODSRC = sgks1/modgraf1.f90 sgks1/modgraf2.f90
#GPYOBJS = pymlibgks1.o pymdlibgks1.o pymlibgks2.o pymlibygl.o
#No graphics
GMODSRC = sgks1/modgraf1.f90 sgks1/modgraf2.f90
GPYOBJS = pymnullgks1.o pymnullgks2.o

# Linkage rules

all : diag1 freaders

freaders : pcreadf1 pvcreadf1 pflreadf1 pfvreadf1 ptreadf1 psreadf1

diag1 : diag1.o modpreadf1.o libcmfield1.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o diag1 \
    diag1.o input1mod.o cmfield1.o libcmfield1_h.o modpreadf1.o \
    libcmfield1.o $(GOBJS) $(LGOPTS)

pcreadf1 : pcreadf1.o libcmfield1.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pcreadf1 \
    pcreadf1.o input1mod.o cmfield1.o libcmfield1_h.o modpreadf1.o \
    libcmfield1.o $(GOBJS) $(LGOPTS)

pvcreadf1 : pvcreadf1.o libcmfield1.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pvcreadf1 \
    pvcreadf1.o input1mod.o cmfield1.o libcmfield1_h.o modpreadf1.o \
    libcmfield1.o $(GOBJS) $(LGOPTS)

pflreadf1 : pflreadf1.o libcmfield1.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pflreadf1 \
    pflreadf1.o input1mod.o cmfield1.o libcmfield1_h.o modpreadf1.o \
    libcmfield1.o $(GOBJS) $(LGOPTS)

pfvreadf1 : pfvreadf1.o libcmfield1.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o pfvreadf1 \
    pfvreadf1.o input1mod.o cmfield1.o libcmfield1_h.o modpreadf1.o \
    libcmfield1.o $(GOBJS) $(LGOPTS)

ptreadf1 : ptreadf1.o libcmfield1.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o ptreadf1 \
    ptreadf1.o input1mod.o cmfield1.o libcmfield1_h.o modpreadf1.o \
    libcmfield1.o $(GOBJS) $(LGOPTS)

psreadf1 : psreadf1.o libcmfield1.o $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o psreadf1 \
    psreadf1.o input1mod.o cmfield1.o libcmfield1_h.o modpreadf1.o \
    libcmfield1.o $(GOBJS) $(LGOPTS)

# Python

python : libcmfield1.so fgraf.so

# Compilation rules

#OPENMP
#omplib.o : omplib.f90
#	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

input1mod.o : input1mod.f90
	$(FC90) $(OPTS90) -o input1mod.o -c input1mod.f90

cmfield1.o : cmfield1.f90 libcmfield1_h.o input1mod.o
	$(FC90) $(OPTS90) -o cmfield1.o -c cmfield1.f90

libcmfield1_h.o : libcmfield1_h.f90
	$(FC90) $(OPTS90) -o libcmfield1_h.o -c libcmfield1_h.f90

modpreadf1.o : modpreadf1.f90 input1mod.o cmfield1.o modgraf1.o modgraf2.o
	$(FC90) $(OPTS90) -o modpreadf1.o -c modpreadf1.f90

diag1.o : diag1.f90 modpreadf1.o
	$(FC90) $(OPTS90) -o diag1.o -c diag1.f90

pcreadf1.o : pcreadf1.f90 input1mod.o cmfield1.o modgraf1.o modgraf2.o
	$(FC90) $(OPTS90) -o pcreadf1.o -c pcreadf1.f90

pvcreadf1.o : pvcreadf1.f90 input1mod.o cmfield1.o modgraf1.o modgraf2.o
	$(FC90) $(OPTS90) -o pvcreadf1.o -c pvcreadf1.f90

pflreadf1.o : pflreadf1.f90 input1mod.o cmfield1.o modgraf1.o
	$(FC90) $(OPTS90) -o pflreadf1.o -c pflreadf1.f90

pfvreadf1.o : pfvreadf1.f90 input1mod.o cmfield1.o modgraf1.o modgraf2.o
	$(FC90) $(OPTS90) -o pfvreadf1.o -c pfvreadf1.f90

ptreadf1.o : ptreadf1.f90 input1mod.o cmfield1.o modgraf1.o
	$(FC90) $(OPTS90) -o ptreadf1.o -c ptreadf1.f90

psreadf1.o : psreadf1.f90 input1mod.o cmfield1.o modgraf1.o modgraf2.o
	$(FC90) $(OPTS90) -o psreadf1.o -c psreadf1.f90

libcmfield1.o : libcmfield1.f
	$(FC90) $(OPTS90) -o libcmfield1.o -c libcmfield1.f

list:
	echo diag1 pcreadf1 pvcreadf1 pflreadf1 pfvreadf1 ptreadf1 psreadf1 python \
	     clean clobber

# Python

pymcmfield1.o : libcmfield1.f
	$(FC90) $(PYOPTS90) -o pymcmfield1.o -c libcmfield1.f

libcmfield1.so : input1mod.f90 libcmfield1_h.o cmfield1.f90 pymcmfield1.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libcmfield1 \
         -c input1mod.f90 libcmfield1_h.o cmfield1.f90 pymcmfield1.o

# Fortran Graphics

libt1.o : sgks1/libt1.f90
	$(FC90) $(OPTS90) -o libt1.o -c sgks1/libt1.f90

plot10.o : sgks1/plot10.f
	$(FC90) $(OPTS90) -o plot10.o -c sgks1/plot10.f

libplt10.o : sgks1/libplt10.f
	$(FC90) $(OPTS90) -o libplt10.o -c sgks1/libplt10.f

sgks1/Ygl-4.2/libYgl.a :
	make -C sgks1/Ygl-4.2

libygl.o : sgks1/libygl.f sgks1/Ygl-4.2/libYgl.a
	$(FC90) $(OPTS90) -o libygl.o -c sgks1/libygl.f

libpsp.o : sgks1/libpsp.f
	$(FC90) $(OPTS90) -o libpsp.o -c sgks1/libpsp.f

dlibgks1.o : sgks1/dlibgks1.f
	$(FC90) $(OPTS90) -o dlibgks1.o -c sgks1/dlibgks1.f

libgks1.o : sgks1/libgks1.f
	$(FC90) $(OPTS90) -o libgks1.o -c sgks1/libgks1.f

nullgks1.o : sgks1/nullgks1.f
	$(FC90) $(OPTS90) -o nullgks1.o -c sgks1/nullgks1.f

libgks2.o : sgks1/libgks2.f
	$(FC90) $(OPTS90) -o libgks2.o -c sgks1/libgks2.f

nullgks2.o : sgks1/nullgks2.f
	$(FC90) $(OPTS90) -o nullgks2.o -c sgks1/nullgks2.f

libgraf1_h.o : sgks1/libgraf1_h.f90
	$(FC90) $(OPTS90) -o libgraf1_h.o -c sgks1/libgraf1_h.f90

modgraf1.o : sgks1/modgraf1.f90 libgraf1_h.o
	$(FC90) $(OPTS90) -o modgraf1.o -c sgks1/modgraf1.f90

libgraf2_h.o : sgks1/libgraf2_h.f90
	$(FC90) $(OPTS90) -o libgraf2_h.o -c sgks1/libgraf2_h.f90

modgraf2.o : sgks1/modgraf2.f90 libgraf2_h.o
	$(FC90) $(OPTS90) -o modgraf2.o -c sgks1/modgraf2.f90

# Python Graphics

pymplot10.o : sgks1/plot10.f
	$(FC90) $(PYOPTS90) -o pymplot10.o -c sgks1/plot10.f

pymlibplt10.o : sgks1/libplt10.f
	$(FC90) $(PYOPTS90) -o pymlibplt10.o -c sgks1/libplt10.f

pymlibygl.o : sgks1/libygl.f sgks1/Ygl-4.2/libYgl.a
	$(FC90) $(PYOPTS90) -o pymlibygl.o -c sgks1/libygl.f

pymlibpsp.o : sgks1/libpsp.f
	$(FC90) $(PYOPTS90) -o pymlibpsp.o -c sgks1/libpsp.f

pymdlibgks1.o : sgks1/dlibgks1.f
	$(FC90) $(PYOPTS90) -o pymdlibgks1.o -c sgks1/dlibgks1.f

pymlibgks1.o : sgks1/libgks1.f
	$(FC90) $(PYOPTS90) -o pymlibgks1.o -c sgks1/libgks1.f

pymnullgks1.o : sgks1/nullgks1.f
	$(FC90) $(PYOPTS90) -o pymnullgks1.o -c sgks1/nullgks1.f

pymlibgks2.o : sgks1/libgks2.f
	$(FC90) $(PYOPTS90) -o pymlibgks2.o -c sgks1/libgks2.f

pymnullgks2.o : sgks1/nullgks2.f
	$(FC90) $(PYOPTS90) -o pymnullgks2.o -c sgks1/nullgks2.f

fgraf.so : libgraf1_h.o libgraf2_h.o $(GMODSRC) $(GPYOBJS)
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) $(LGOPTS) -m fgraf \
        -c libgraf1_h.o libgraf2_h.o $(GMODSRC) $(GPYOBJS)

cpexecs : cpcreadf1 cpvcreadf1 cpflreadf1 cpfvreadf1 cptreadf1 cpsreadf1

cpdiag1 : diag1
	cp diag1 ../

cpcreadf1 : pcreadf1
	cp pcreadf1 ../

cpvcreadf1 : pvcreadf1
	cp pvcreadf1 ../

cpflreadf1 : pflreadf1
	cp pflreadf1 ../

cpfvreadf1 : pfvreadf1
	cp pfvreadf1 ../

cptreadf1 : ptreadf1
	cp ptreadf1 ../

cpsreadf1 : psreadf1
	cp psreadf1 ../

clean :
	rm -f *.o *.mod *.pyf sgks1/Ygl-4.2/*.o

clobber: clean
	rm -f diag1 pcreadf1 pvcreadf1 pflreadf1 pfvreadf1 ptreadf1 psreadf1 \
	      *.so sgks1/Ygl-4.2/*.a

#Python
cpyexecs : cpycreadf1 cpyvcreadf1 cpyflreadf1 cpyfvreadf1 cpytreadf1 \
           cpysreadf1

cpycreadf1 : pcreadf1.py
	cp pcreadf1.py ../

cpyvcreadf1 : pvcreadf1.py
	cp pvcreadf1.py ../
	
cpyflreadf1 : pflreadf1.py
	cp pflreadf1.py ../

cpyfvreadf1 : pfvreadf1.py
	cp pfvreadf1.py ../

cpytreadf1 : ptreadf1.py
	cp ptreadf1.py ../

cpysreadf1 : psreadf1.py
	cp psreadf1.py ../
