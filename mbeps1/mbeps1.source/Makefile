#Makefile for 1D OpenMP PIC codes: mbeps1, mbbeps1, mdbeps1

# Makefile gfortran compiler with MacOS X

#OpenMP
#MPFC = gfortran -fopenmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

#FC90 = gfortran
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
#LGOPTS = gks1/Ygl-4.2/libYgl.a -L/usr/X11R6/lib -lX11
#LEGACY = -std=legacy

# Makefile Intel compiler with Mac OS X

#OpenMP
#MPFC = ifort -openmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = ifort
#MPCC = gcc

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#FF03 = -Tf

#NOMAIN = -nofor_main

# Makefile Intel compiler with Linux

# OpenMP
#MPFC = ifort -openmp
#MPCC = icc -openmp
#NoOpenMP
#MPFC = ifort
#MPCC = icc

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99

#FF03 = -Tf

#NOMAIN = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMP
MPFC = gfortran -fopenmp
MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

CC = gcc

PFC = f2py --fcompiler=gnu95
PCC = f2py --compiler=unix

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#Python
PYOPTS90 = $(OPTS90) -fPIC
PYCCOPTS = $(CCOPTS) -fPIC
#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

LOPTS = -lgomp
#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
#LGOPTS = gks1/Ygl-4.2/libYgl.a -L/usr/X11R6/lib -lX11
LEGACY = -std=legacy

# Makefile PGI compiler with Linux

# OpenMP
#MPFC = pgf90
#MPCC = gcc -fopenmp

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#NOMAIN = -Mnomain

# Makefile Nag compiler with Linux

# OpenMP
#MPFC = nagfor -openmp
#MPCC = gcc -fopenmp

#FC90 = nagfor
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -default_kinds:r=64
#OPTS90 = -O3 -default_kinds:r=64 -C=array -nan -w=all -f95

#CCOPTS = -O3 -Wall -std=c99

# Makefile Cray compiler with Linux

# OpenMP
#MPFC = ftn -h omp
#MPCC = cc -h omp
#NoOpenMP
#MPFC = ftn -h noomp
#MPCC = cc -h noomp

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#

ESOBJS = libminit1.o libmpush1.o libmsort1.o libmgard1.o libmfft1.o \
libmfield1.o libmdiag1.o

EMOBJS = libmbpush1.o libmcurd1.o

EDOBJS = libmdpush1.o

ESHOBJS = libminit1_h.o libmpush1_h.o libmsort1_h.o libmgard1_h.o \
libmfft1_h.o libmfield1_h.o libmdiag1_h.o

EMHOBJS = libmbpush1_h.o libmcurd1_h.o

EDHOBJS = libmdpush1_h.o

ESMODS = input1mod.o modminit1.o modmpush1.o modmsort1.o modmgard1.o \
modmfft1.o modmfield1.o modmdiag1.o

EMMODS = modmbpush1.o modmcurd1.o

EDMODS = modmdpush1.o

# Fortran Graphics

#Tektronix
#GOBJS = libt1.o plot10.o libplt10.o dlibgks1.o libgks1.o libgraf1_h.o \
#modgraf1.o
# Postcript printer
#GOBJS = libpsp.o dlibgks1.o libgks1.o  libgraf1_h.o modgraf1.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libygl.o dlibgks1.o libgks1.o libgraf1_h.o modgraf1.o
#No graphics
GOBJS = nullgks1.o libgraf1_h.o modgraf1.o

# Python

ESMODSRC = input1mod.f90 modminit1.f90 modmpush1.f90 modmsort1.f90 \
modmgard1.f90 modmfft1.f90 modmfield1.f90 modmdiag1.f90

EMMODSRC = modmbpush1.f90 modmcurd1.f90

EDMODSRC = modmdpush1.f90

ESPYOBJS = pyminit1.o pympush1.o pymsort1.o pymgard1.o pymfft1.o \
pymfield1.o pymdiag1.o

EMPYOBJS = pymbpush1.o pymcurd1.o

EDPYOBJS = pymdpush1.o

# Python Graphics

#Tektronix
#GMODSRC = gks1/modgraf1.f90 gks1/libt1.f90
#GPYOBJS = pymlibgks1.o pymdlibgks1.o pymlibplt10.o pymplot10.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GMODSRC = gks1/modgraf1.f90
#GPYOBJS = pymlibgks1.o pymdlibgks1.o pymlibygl.o
#No graphics
GMODSRC = gks1/modgraf1.f90
GPYOBJS = pymnullgks1.o

# Linkage rules

all : mbeps1 mbbeps1 mdbeps1

python : mbeps1.py mbbeps1.py mdbeps1.py

mbeps1 : mbeps1.o $(ESOBJS) $(GOBJS) dtimer.o
	$(MPFC) $(OPTS90) $(LOPTS) -o mbeps1 \
    mbeps1.o $(ESOBJS) $(ESMODS) omplib.o $(ESHOBJS) $(GOBJS) msimul1.o \
    $(LGOPTS) dtimer.o

mbbeps1 : mbbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS) dtimer.o
	$(MPFC) $(OPTS90) $(LOPTS) -o mbbeps1 \
    mbbeps1.o $(ESOBJS) $(EMOBJS) $(ESMODS) $(EMMODS) omplib.o \
    $(ESHOBJS) $(EMHOBJS) $(GOBJS) msimul1.o mbsimul1.o $(LGOPTS) dtimer.o

mdbeps1 : mdbeps1.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(GOBJS) dtimer.o
	$(MPFC) $(OPTS90) $(LOPTS) -o mdbeps1 \
    mdbeps1.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(ESMODS) $(EMMODS) $(EDMODS) \
    omplib.o $(ESHOBJS) $(EMHOBJS) $(EDHOBJS) $(GOBJS) msimul1.o mbsimul1.o \
    mdsimul1.o $(LGOPTS) dtimer.o

# Python

mbeps1.py : libmpush1.so fomplib.so dtimer.so fgraf1.so

mbbeps1.py : libmbpush1.so

mdbeps1.py : libmdpush1.so

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#OPENMP
omplib.o : omplib.f90
	$(MPFC) $(OPTS90) -o omplib.o -c omplib.f90

libminit1.o : libminit1.f
	$(MPFC) $(OPTS90) -o libminit1.o -c libminit1.f

libmpush1.o : libmpush1.f
	$(MPFC) $(OPTS90) -o libmpush1.o -c libmpush1.f

libmbpush1.o : libmbpush1.f
	$(MPFC) $(OPTS90) -o libmbpush1.o -c libmbpush1.f

libmcurd1.o : libmcurd1.f
	$(MPFC) $(OPTS90) -o libmcurd1.o -c libmcurd1.f

libmdpush1.o : libmdpush1.f
	$(MPFC) $(OPTS90) -o libmdpush1.o -c libmdpush1.f

libmsort1.o : libmsort1.f
	$(MPFC) $(OPTS90) -o libmsort1.o -c libmsort1.f

libmgard1.o : libmgard1.f
	$(MPFC) $(OPTS90) -o libmgard1.o -c libmgard1.f

libmfft1.o : libmfft1.f
	$(MPFC) $(OPTS90) -o libmfft1.o -c libmfft1.f

libmfield1.o : libmfield1.f
	$(MPFC) $(OPTS90) -o libmfield1.o -c libmfield1.f

libmdiag1.o : libmdiag1.f
	$(MPFC) $(OPTS90) -o libmdiag1.o -c libmdiag1.f

libminit1_h.o : libminit1_h.f90
	$(MPFC) $(OPTS90) -o libminit1_h.o -c libminit1_h.f90

libmpush1_h.o : libmpush1_h.f90
	$(MPFC) $(OPTS90) -o libmpush1_h.o -c libmpush1_h.f90

libmbpush1_h.o : libmbpush1_h.f90
	$(MPFC) $(OPTS90) -o libmbpush1_h.o -c libmbpush1_h.f90

libmcurd1_h.o : libmcurd1_h.f90
	$(MPFC) $(OPTS90) -o libmcurd1_h.o -c libmcurd1_h.f90

libmdpush1_h.o : libmdpush1_h.f90
	$(MPFC) $(OPTS90) -o libmdpush1_h.o -c libmdpush1_h.f90

libmsort1_h.o : libmsort1_h.f90
	$(MPFC) $(OPTS90) -o libmsort1_h.o -c libmsort1_h.f90

libmgard1_h.o : libmgard1_h.f90
	$(MPFC) $(OPTS90) -o libmgard1_h.o -c libmgard1_h.f90

libmfft1_h.o : libmfft1_h.f90
	$(MPFC) $(OPTS90) -o libmfft1_h.o -c libmfft1_h.f90

libmfield1_h.o : libmfield1_h.f90
	$(MPFC) $(OPTS90) -o libmfield1_h.o -c libmfield1_h.f90

libmdiag1_h.o : libmdiag1_h.f90
	$(MPFC) $(OPTS90) -o libmdiag1_h.o -c libmdiag1_h.f90

input1mod.o : input1mod.f90
	$(MPFC) $(OPTS90) -o input1mod.o -c input1mod.f90

modminit1.o : modminit1.f90 libminit1_h.o
	$(MPFC) $(OPTS90) -o modminit1.o -c modminit1.f90

modmpush1.o : modmpush1.f90 libmpush1_h.o
	$(MPFC) $(OPTS90) -o modmpush1.o -c modmpush1.f90

modmbpush1.o : modmbpush1.f90 libmbpush1_h.o
	$(MPFC) $(OPTS90) -o modmbpush1.o -c modmbpush1.f90

modmcurd1.o : modmcurd1.f90 libmcurd1_h.o
	$(MPFC) $(OPTS90) -o modmcurd1.o -c modmcurd1.f90

modmdpush1.o : modmdpush1.f90 libmdpush1_h.o
	$(MPFC) $(OPTS90) -o modmdpush1.o -c modmdpush1.f90

modmsort1.o : modmsort1.f90 libmsort1_h.o
	$(MPFC) $(OPTS90) -o modmsort1.o -c modmsort1.f90

modmgard1.o : modmgard1.f90 libmgard1_h.o
	$(MPFC) $(OPTS90) -o modmgard1.o -c modmgard1.f90

modmfft1.o : modmfft1.f90 libmfft1_h.o
	$(MPFC) $(OPTS90) -o modmfft1.o -c modmfft1.f90

modmfield1.o : modmfield1.f90 libmfield1_h.o
	$(MPFC) $(OPTS90) -o modmfield1.o -c modmfield1.f90

modmdiag1.o : modmdiag1.f90 libmdiag1_h.o input1mod.o
	$(MPFC) $(OPTS90) -o modmdiag1.o -c modmdiag1.f90

msimul1.o : msimul1.f03 $(ESMODS)
	$(MPFC) $(OPTS90) -o msimul1.o -c $(FF03) msimul1.f03

mbsimul1.o : mbsimul1.f03 msimul1.o $(EMMODS)
	$(MPFC) $(OPTS90) -o mbsimul1.o -c $(FF03) mbsimul1.f03

mdsimul1.o : mdsimul1.f03 mbsimul1.o $(EDMODS)
	$(MPFC) $(OPTS90) -o mdsimul1.o -c $(FF03) mdsimul1.f03

mbeps1.o : mbeps1.f90 msimul1.o modgraf1.o omplib.o
	$(MPFC) $(OPTS90) -o mbeps1.o -c mbeps1.f90

mbbeps1.o : mbbeps1.f90 msimul1.o mbsimul1.o modgraf1.o omplib.o
	$(MPFC) $(OPTS90) -o mbbeps1.o -c mbbeps1.f90

mdbeps1.o : mdbeps1.f90 msimul1.o mbsimul1.o mdsimul1.o modgraf1.o omplib.o 
	$(MPFC) $(OPTS90) -o mdbeps1.o -c mdbeps1.f90

# Fortran Graphics

libt1.o : gks1/libt1.f90
	$(MPFC) $(OPTS90) -o libt1.o -c gks1/libt1.f90

plot10.o : gks1/plot10.f
	$(MPFC) $(OPTS90) -o plot10.o -c gks1/plot10.f

libplt10.o : gks1/libplt10.f
	$(MPFC) $(OPTS90) -o libplt10.o -c gks1/libplt10.f

gks1/Ygl-4.2/libYgl.a :
	make -C gks1/Ygl-4.2

libygl.o : gks1/libygl.f gks1/Ygl-4.2/libYgl.a
	$(MPFC) $(OPTS90) -o libygl.o -c gks1/libygl.f

libpsp.o : gks1/libpsp.f
	$(MPFC) $(OPTS90) -o libpsp.o -c gks1/libpsp.f

dlibgks1.o : gks1/dlibgks1.f
	$(MPFC) $(OPTS90) -o dlibgks1.o -c gks1/dlibgks1.f

libgks1.o : gks1/libgks1.f
	$(MPFC) $(OPTS90) -o libgks1.o -c gks1/libgks1.f

nullgks1.o : gks1/nullgks1.f
	$(MPFC) $(OPTS90) -o nullgks1.o -c gks1/nullgks1.f

libgraf1_h.o : gks1/libgraf1_h.f90
	$(MPFC) $(OPTS90) -o libgraf1_h.o -c gks1/libgraf1_h.f90

modgraf1.o : gks1/modgraf1.f90 libgraf1_h.o
	$(MPFC) $(OPTS90) -o modgraf1.o -c gks1/modgraf1.f90

libgks2.o : gks1/libgks2.f
	$(MPFC) $(OPTS90) -o libgks2.o -c gks1/libgks2.f

libgraf2_h.o : gks1/libgraf2_h.f90
	$(MPFC) $(OPTS90) -o libgraf2_h.o -c gks1/libgraf2_h.f90

modgraf2.o : gks1/modgraf2.f90 libgraf2_h.o
	$(MPFC) $(OPTS90) -o modgraf2.o -c gks1/modgraf2.f90

nullgks2.o : gks1/nullgks2.f
	$(MPFC) $(OPTS90) -o nullgks2.o -c gks1/nullgks2.f

# Python

pydtimer.o : dtimer.c
	$(CC) $(PYCCOPTS) -o pydtimer.o -c dtimer.c

pyminit1.o : libminit1.f
	$(MPFC) $(PYOPTS90) -o pyminit1.o -c libminit1.f

pympush1.o : libmpush1.f
	$(MPFC) $(PYOPTS90) -o pympush1.o -c libmpush1.f

pymbpush1.o : libmbpush1.f
	$(MPFC) $(PYOPTS90) -o pymbpush1.o -c libmbpush1.f

pymcurd1.o : libmcurd1.f
	$(MPFC) $(PYOPTS90) -o pymcurd1.o -c libmcurd1.f

pymdpush1.o : libmdpush1.f
	$(MPFC) $(PYOPTS90) -o pymdpush1.o -c libmdpush1.f

pymsort1.o : libmsort1.f
	$(MPFC) $(PYOPTS90) -o pymsort1.o -c libmsort1.f

pymgard1.o : libmgard1.f
	$(MPFC) $(PYOPTS90) -o pymgard1.o -c libmgard1.f

pymfft1.o : libmfft1.f
	$(MPFC) $(PYOPTS90) -o pymfft1.o -c libmfft1.f

pymfield1.o : libmfield1.f
	$(MPFC) $(PYOPTS90) -o pymfield1.o -c libmfield1.f

pymdiag1.o : libmdiag1.f
	$(MPFC) $(PYOPTS90) -o pymdiag1.o -c libmdiag1.f

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

fomplib.so : omplib.f90
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m fomplib -c omplib.f90

libmpush1.so : $(ESHOBJS) $(ESMODSRC) $(ESPYOBJS) pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmpush1 \
        -c $(ESHOBJS) $(ESMODSRC) $(ESPYOBJS) pydtimer.o

libmbpush1.so : $(EMHOBJS) $(EMMODSRC) $(EMPYOBJS) pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmbpush1 \
        -c $(EMHOBJS) $(EMMODSRC) $(EMPYOBJS) pydtimer.o

libmdpush1.so : $(EDHOBJS) $(EDMODSRC) $(EDPYOBJS) pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmdpush1 \
        -c $(EDHOBJS) $(EDMODSRC) $(EDPYOBJS) pydtimer.o

# Python Graphics

pymplot10.o : gks1/plot10.f
	$(MPFC) $(PYOPTS90) -o pymplot10.o -c gks1/plot10.f

pymlibplt10.o : gks1/libplt10.f
	$(MPFC) $(PYOPTS90) -o pymlibplt10.o -c gks1/libplt10.f

pymlibygl.o : gks1/libygl.f gks1/Ygl-4.2/libYgl.a
	$(MPFC) $(PYOPTS90) -o pymlibygl.o -c gks1/libygl.f

pymlibpsp.o : gks1/libpsp.f
	$(MPFC) $(PYOPTS90) -o pymlibpsp.o -c gks1/libpsp.f

pymdlibgks1.o : gks1/dlibgks1.f
	$(MPFC) $(PYOPTS90) -o pymdlibgks1.o -c gks1/dlibgks1.f

pymlibgks1.o : gks1/libgks1.f
	$(MPFC) $(PYOPTS90) -o pymlibgks1.o -c gks1/libgks1.f

pymnullgks1.o : gks1/nullgks1.f
	$(MPFC) $(PYOPTS90) -o pymnullgks1.o -c gks1/nullgks1.f

fgraf1.so : libgraf1_h.o $(GMODSRC) $(GPYOBJS)
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) $(LGOPTS) -m fgraf1 \
        -c libgraf1_h.o $(GMODSRC) $(GPYOBJS)

cpexecs : cpmbeps1 cpmbbeps1 cpmdbeps1

cpmbeps1 :
	cp mbeps1 ../

cpmbbeps1 :
	cp mbbeps1 ../

cpmdbeps1 :
	cp mdbeps1 ../

clean :
	rm -f *.o *.mod *.pyf gks1/Ygl-4.2/*.o

clobber: clean
	rm -f mbeps1 mbbeps1 mdbeps1 *.so gks1/Ygl-4.2/*.a
