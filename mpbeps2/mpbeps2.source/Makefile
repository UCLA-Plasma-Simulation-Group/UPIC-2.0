#Makefile for 2D MPI/OpenMP PIC codes

# Makefile gfortran compiler with MacOS X

#OpenMPI
#MPIFC = mpif90 -fopenmp
#MPICC = mpicc -fopenmp

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

# OpenMPI
#MPIFC = mpif90 -qopenmp
#MPICC = mpicc -fopenmp

#FC90 = ifort -qopenmp
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#FF03 = -Tf

#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

# MPI
#MPIFC = mpiifort -qopenmp
#MPICC = mpiicc -qopenmp

#FC90 = ifort -qopenmp
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -axAVX -qopt-report=5 -qopt-report-phase=vec
#OPTS90 = -O3 -axAVX2 -qopt-report=5 -qopt-report-phase=vec
#OPTS90 = -O3 -axMIC-AVX512 -qopt-report=5 -qopt-report-phase=vec
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#FF03 = -Tf

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99

#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
#LGOPTS = gks2/Ygl-4.2/libYgl.a -L/usr/X11R6/lib -lX11
#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMPI
# MPIFC = mpif90 -fopenmp
# MPICC = mpicc -fopenmp

# FC90 = gfortran -fopenmp
# CC = gcc

# PFC = f2py --fcompiler=gnu95
# PCC = f2py --compiler=unix

# OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

# CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#Python
# PYOPTS90 = $(OPTS90) -fPIC
# PYCCOPTS = $(CCOPTS) -fPIC

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

# LOPTS = -lgomp
#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
#LGOPTS = gks2/Ygl-4.2/libYgl.a -L/usr/X11R6/lib -lX11
LEGACY =

# Makefile PGI compiler with Linux

# OpenMPI
#MPIFC = mpipgf90
#MPICC = mpicc -fopenmp

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#LEGACY = -Mnomain

# Makefile Cray compiler with Linux

# MPI
#MPIFC = ftn -h omp
#MPICC = cc -h omp

#FC90 = ftn -h omp
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#LEGACY =

#

####################################################
### System Branching Starts
### Makefile starts
####################################################

# legal values of compilers are either GNU or INTEL
#
compiler = GNU
#
#

hostname := $(shell uname -n)
osname := $(shell uname -s)

cori = $(shell uname -n | cut -c 1-4)
hoffman = $(shell uname -n | cut -c 1-5)

# DEBUG
# $(info $(hostname))
# $(info $(osname))
# $(info $(cori))
# $(info $(hoffman))
# DEBUG

# This is where the machine dependencies are resolved

# Case 1: OS X (Darwin)
ifeq ($(osname),Darwin)
	ifeq ($(compiler),GNU)
		MPIFC = mpif90 -fopenmp
		MPICC = mpicc -fopenmp
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
	else
		MPIFC = mpiifort -qopenmp
		MPICC = mpiicc -qopenmp
		FC90 = ifort -qopenmp
		CC = icc
		OPTS90 = -O3
		# OPTS90 = -O3 -axAVX -qopt-report=5 -qopt-report-phase=vec
		# OPTS90 = -O3 -axAVX2 -qopt-report=5 -qopt-report-phase=vec
		# OPTS90 = -O3 -axMIC-AVX512 -qopt-report=5 -qopt-report-phase=vec
		# OPTS90 = -O3 -r8
		# OPTS90 = -O3 -CB -r8 -warn all -std90
		CCOPTS = -O3 -std=c99
		#CCOPTS = -O3 -no-vec -Wall -std=c99
		FF03 = -Tf
		#LEGACY = -nofor_main
	endif
else
	ifeq ($(cori),cori)
       	sys = CORI
	endif
	ifeq ($(hoffman),login)
		MPIFC = mpiifort -qopenmp
		MPICC = mpiicc -qopenmp
		FC90 = ifort -qopenmp
		CC = icc
		OPTS90 = -O3
		# OPTS90 = -O3 -axAVX -qopt-report=5 -qopt-report-phase=vec
		# OPTS90 = -O3 -axAVX2 -qopt-report=5 -qopt-report-phase=vec
		# OPTS90 = -O3 -axMIC-AVX512 -qopt-report=5 -qopt-report-phase=vec
		# OPTS90 = -O3 -r8
		# OPTS90 = -O3 -CB -r8 -warn all -std90
		CCOPTS = -O3 -std=c99
		#CCOPTS = -O3 -no-vec -Wall -std=c99
		FF03 = -Tf
		LEGACY = -nofor_main
	endif
endif
####################################################
####################################################
####################################################
####################################################

ESOBJS = libmpinit2.o libmppush2.o libvmppush2.o libvmpsort2.o libmpsort2.o \
libmpgard2.o libvmpfft2.o libmpfield2.o libmpdiag2.o

EMOBJS = libmpbpush2.o libvmpbpush2.o libmpcurd2.o libvmpcurd2.o

EDOBJS = libmpdpush2.o libvmpdpush2.o

ESHOBJS = libmpinit2_h.o libmppush2_h.o libvmppush2_h.o libmpsort2_h.o \
libvmpsort2_h.o libmpgard2_h.o libvmpfft2_h.o libmpfield2_h.o libmpdiag2_h.o

EMHOBJS = libmpbpush2_h.o libvmpbpush2_h.o libmpcurd2_h.o libvmpcurd2_h.o

EDHOBJS = libmpdpush2_h.o libvmpdpush2_h.o

ESMODS = input2mod.o modmpinit2.o modmppush2.o modmpsort2.o modmpgard2.o \
modmpfft2.o modmpfield2.o

EMMODS = modmpbpush2.o modmpcurd2.o

EDMODS = modmpdpush2.o

# Python

ESMODSRC =  ../ReadFiles/input2mod.f90 modmpinit2.f90 modmppush2.f90 modmpsort2.f90 \
modmpgard2.f90 modmpfft2.f90 modmpfield2.f90

EMMODSRC = modmpbpush2.f90 modmpcurd2.f90

EDMODSRC = modmpdpush2.f90

MPMODSRC = noMPI.f90 mpplib2.f90 mppmod2.f90 modmpdiag2.f90 omplib.f90 \
ompplib2.f90

ESPYOBJS = pympinit2.o pymppush2.o pyvmppush2.o pympsort2.o pyvmpsort2.o \
pympgard2.o pympfft2.o pympfield2.o pympdiag2.o

EMPYOBJS = pympbpush2.o pyvmpbpush2.o pympcurd2.o pyvmpcurd2.o

EDPYOBJS = pympdpush2.o pyvmpdpush2.o

# Fortran Graphics

#Tektronix
#GOBJS = libt1.o plot10.o libplt10.o libgks1.o libgks2.o plibgks2.o \
#modpgraf2.o modpgraf1.o
#NGOBJS = libt1.o plot10.o libplt10.o libgks1.o libgks2.o nplibgks2.o \
#nmodpgraf2.o nmodpgraf1.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libygl.o libgks1.o libgks2.o plibgks2.o modpgraf2.o modpgraf1.o
#NGOBJS = libygl.o libgks1.o libgks2.o nplibgks2.o nmodpgraf2.o nmodpgraf1.o
#No graphics
GOBJS = nullgks1.o nullgks2.o nullpgks2.o modpgraf2.o modpgraf1.o
NGOBJS = nullgks1.o nullgks2.o nullpgks2.o nmodpgraf2.o nmodpgraf1.o

# Python Graphics

#Tektronix
#GMODSRC = noMPI.f90 nompplib2.f90 gks2/plibgks2.f90 gks2/modpgraf2.f90 \
#gks2/modpgraf1.f90 gks2/libt1.f90
#GPYOBJS = pymlibgks1.o pymlibgks2.o pymlibplt10.o pymplot10.o pynoMPI.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GMODSRC = noMPI.f90 nompplib2.f90 gks2/plibgks2.f90 gks2/modpgraf2.f90 \
#gks2/modpgraf1.f90
#GPYOBJS = pymlibgks1.o pymlibgks2.o pymlibygl.o pynoMPI.o
#No graphics
GMODSRC = gks2/nullpgks2.f90 gks2/modpgraf2.f90 gks2/modpgraf1.f90
GPYOBJS = pymnullgks1.o pymnullgks2.o pynoMPI.o

# Linkage rules

all : mpbeps2 mpbbeps2 mpdbeps2

nompi : mbeps2 mbbeps2 mdbeps2

python : mbeps2.py mbbeps2.py mdbeps2.py fpgraf.so 

mpbeps2 : mpbeps2.o $(ESOBJS) $(GOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpbeps2 \
    mpbeps2.o $(ESOBJS) $(ESMODS) modmpdiag2.o mpplib2.o mppmod2.o omplib.o \
    ompplib2.o $(ESHOBJS) $(GOBJS) mpsimul2.o $(LGOPTS) dtimer.o

mpbbeps2 : mpbbeps2.o $(ESOBJS) $(EMOBJS) $(GOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpbbeps2 \
    mpbbeps2.o $(ESOBJS) $(EMOBJS) $(ESMODS) $(EMMODS) modmpdiag2.o mpplib2.o \
    mppmod2.o omplib.o ompplib2.o $(ESHOBJS) $(EMHOBJS) $(GOBJS) mpsimul2.o \
    mpbsimul2.o $(LGOPTS) dtimer.o

mpdbeps2 : mpdbeps2.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(GOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpdbeps2 \
    mpdbeps2.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(ESMODS) $(EMMODS) $(EDMODS) \
    modmpdiag2.o mpplib2.o mppmod2.o omplib.o ompplib2.o $(ESHOBJS) $(EMHOBJS) \
    $(EDHOBJS) $(GOBJS) mpsimul2.o mpbsimul2.o mpdsimul2.o $(LGOPTS) dtimer.o

#noMPI
mbeps2 : mbeps2.o $(ESOBJS) noMPI.o $(NGOBJS) dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mbeps2 \
    mbeps2.o $(ESOBJS) $(ESMODS) nmodmpdiag2.o nmpplib2.o nmppmod2.o omplib.o \
    nompplib2.o $(ESHOBJS) noMPI.o $(NGOBJS) msimul2.o $(LGOPTS) dtimer.o

mbbeps2 : mbbeps2.o $(ESOBJS) $(EMOBJS) noMPI.o $(NGOBJS) dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mbbeps2 \
    mbbeps2.o $(ESOBJS) $(EMOBJS) $(ESMODS) $(EMMODS) nmodmpdiag2.o nmpplib2.o \
    nmppmod2.o omplib.o nompplib2.o $(ESHOBJS) $(EMHOBJS) noMPI.o $(NGOBJS) \
    msimul2.o mbsimul2.o $(LGOPTS) dtimer.o

mdbeps2 : mdbeps2.o $(ESOBJS) $(EMOBJS) $(EDOBJS) noMPI.o $(NGOBJS) dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mdbeps2 \
    mdbeps2.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(ESMODS) $(EMMODS) $(EDMODS) \
    nmodmpdiag2.o nmpplib2.o nmppmod2.o omplib.o nompplib2.o $(ESHOBJS) \
    $(EMHOBJS) $(EDHOBJS) noMPI.o $(NGOBJS) msimul2.o mbsimul2.o mdsimul2.o \
    $(LGOPTS) dtimer.o

# Python

mbeps2.py : libmpush2.so dtimer.so

mbbeps2.py : libmbpush2.so

mdbeps2.py : libmdpush2.so

list:
	echo mpbeps2 mpbbeps2 mpdbeps2 nompi mbeps2 mbbeps2 mdbeps2 python \
	     clean clobber

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#noMPI
noMPImod.o : noMPI.f90
	$(FC90) $(OPTS90) -o noMPImod.o -c noMPI.f90

noMPI.o : noMPI.f
	$(FC90) $(OPTS90) -o noMPI.o -c noMPI.f

nmpplib2.o : mpplib2.f90 noMPImod.o
	$(FC90) $(OPTS90) -o nmpplib2.o -c mpplib2.f90

nmppmod2.o : mppmod2.f90 nmpplib2.o
	$(FC90) $(OPTS90) -o nmppmod2.o -c mppmod2.f90

nompplib2.o : ompplib2.f90 modmpsort2.o modmpfft2.o modmpgard2.o nmppmod2.o
	$(FC90) $(OPTS90) -o nompplib2.o -c ompplib2.f90

#OPENMP
omplib.o : omplib.f90
	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

mpplib2.o : mpplib2.f90
	$(MPIFC) $(OPTS90) -o mpplib2.o -c mpplib2.f90

mppmod2.o : mppmod2.f90 mpplib2.o
	$(FC90) $(OPTS90) -o mppmod2.o -c mppmod2.f90

libmpinit2.o : libmpinit2.f
	$(FC90) $(OPTS90) -o libmpinit2.o -c libmpinit2.f

libmppush2.o : libmppush2.f
	$(FC90) $(OPTS90) -o libmppush2.o -c libmppush2.f

libmpbpush2.o : libmpbpush2.f
	$(FC90) $(OPTS90) -o libmpbpush2.o -c libmpbpush2.f

libmpcurd2.o : libmpcurd2.f
	$(FC90) $(OPTS90) -o libmpcurd2.o -c libmpcurd2.f

libmpdpush2.o : libmpdpush2.f
	$(FC90) $(OPTS90) -o libmpdpush2.o -c libmpdpush2.f

libvmppush2.o : libvmppush2.f
	$(FC90) $(OPTS90) -o libvmppush2.o -c libvmppush2.f

libvmpbpush2.o : libvmpbpush2.f
	$(FC90) $(OPTS90) -o libvmpbpush2.o -c libvmpbpush2.f

libvmpcurd2.o : libvmpcurd2.f
	$(FC90) $(OPTS90) -o libvmpcurd2.o -c libvmpcurd2.f

libvmpdpush2.o : libvmpdpush2.f
	$(FC90) $(OPTS90) -o libvmpdpush2.o -c libvmpdpush2.f

libmpsort2.o : libmpsort2.f
	$(FC90) $(OPTS90) -o libmpsort2.o -c libmpsort2.f

libvmpsort2.o : libvmpsort2.f
	$(FC90) $(OPTS90) -o libvmpsort2.o -c libvmpsort2.f

libmpgard2.o : libmpgard2.f
	$(FC90) $(OPTS90) -o libmpgard2.o -c libmpgard2.f

libvmpfft2.o : libvmpfft2.f
	$(FC90) $(OPTS90) -o libvmpfft2.o -c libvmpfft2.f

libmpfield2.o : libmpfield2.f
	$(FC90) $(OPTS90) -o libmpfield2.o -c libmpfield2.f

libmpdiag2.o : libmpdiag2.f
	$(FC90) $(OPTS90) -o libmpdiag2.o -c libmpdiag2.f

libmpinit2_h.o : libmpinit2_h.f90
	$(FC90) $(OPTS90) -o libmpinit2_h.o -c libmpinit2_h.f90

libmppush2_h.o : libmppush2_h.f90
	$(FC90) $(OPTS90) -o libmppush2_h.o -c libmppush2_h.f90

libmpbpush2_h.o : libmpbpush2_h.f90
	$(FC90) $(OPTS90) -o libmpbpush2_h.o -c libmpbpush2_h.f90

libmpcurd2_h.o : libmpcurd2_h.f90
	$(FC90) $(OPTS90) -o libmpcurd2_h.o -c libmpcurd2_h.f90

libmpdpush2_h.o : libmpdpush2_h.f90
	$(FC90) $(OPTS90) -o libmpdpush2_h.o -c libmpdpush2_h.f90

libvmppush2_h.o : libvmppush2_h.f90
	$(FC90) $(OPTS90) -o libvmppush2_h.o -c libvmppush2_h.f90

libvmpbpush2_h.o : libvmpbpush2_h.f90
	$(FC90) $(OPTS90) -o libvmpbpush2_h.o -c libvmpbpush2_h.f90

libvmpcurd2_h.o : libvmpcurd2_h.f90
	$(FC90) $(OPTS90) -o libvmpcurd2_h.o -c libvmpcurd2_h.f90

libvmpdpush2_h.o : libvmpdpush2_h.f90
	$(FC90) $(OPTS90) -o libvmpdpush2_h.o -c libvmpdpush2_h.f90

libmpsort2_h.o : libmpsort2_h.f90
	$(FC90) $(OPTS90) -o libmpsort2_h.o -c libmpsort2_h.f90

libvmpsort2_h.o : libvmpsort2_h.f90
	$(FC90) $(OPTS90) -o libvmpsort2_h.o -c libvmpsort2_h.f90

libmpgard2_h.o : libmpgard2_h.f90
	$(FC90) $(OPTS90) -o libmpgard2_h.o -c libmpgard2_h.f90

libvmpfft2_h.o : libvmpfft2_h.f90
	$(FC90) $(OPTS90) -o libvmpfft2_h.o -c libvmpfft2_h.f90

libmpfield2_h.o : libmpfield2_h.f90
	$(FC90) $(OPTS90) -o libmpfield2_h.o -c libmpfield2_h.f90

libmpdiag2_h.o : libmpdiag2_h.f90
	$(FC90) $(OPTS90) -o libmpdiag2_h.o -c libmpdiag2_h.f90

input2mod.o : input2mod.f90
	$(FC90) $(OPTS90) -o input2mod.o -c input2mod.f90

modmpinit2.o : modmpinit2.f90 libmpinit2_h.o
	$(FC90) $(OPTS90) -o modmpinit2.o -c modmpinit2.f90

modmppush2.o : modmppush2.f90 libmppush2_h.o libvmppush2_h.o
	$(FC90) $(OPTS90) -o modmppush2.o -c modmppush2.f90

modmpbpush2.o : modmpbpush2.f90 libmpbpush2_h.o libvmpbpush2_h.o
	$(FC90) $(OPTS90) -o modmpbpush2.o -c modmpbpush2.f90

modmpcurd2.o : modmpcurd2.f90 libmpcurd2_h.o libvmpcurd2_h.o
	$(FC90) $(OPTS90) -o modmpcurd2.o -c modmpcurd2.f90

modmpdpush2.o : modmpdpush2.f90 libmpdpush2_h.o libvmpdpush2_h.o
	$(FC90) $(OPTS90) -o modmpdpush2.o -c modmpdpush2.f90

modmpsort2.o : modmpsort2.f90 libmpsort2_h.o libvmpsort2_h.o
	$(FC90) $(OPTS90) -o modmpsort2.o -c modmpsort2.f90

modmpgard2.o : modmpgard2.f90 libmpgard2_h.o
	$(FC90) $(OPTS90) -o modmpgard2.o -c modmpgard2.f90

modmpfft2.o : modmpfft2.f90 libvmpfft2_h.o
	$(FC90) $(OPTS90) -o modmpfft2.o -c modmpfft2.f90

modmpfield2.o : modmpfield2.f90 libmpfield2_h.o
	$(FC90) $(OPTS90) -o modmpfield2.o -c modmpfield2.f90

modmpdiag2.o : modmpdiag2.f90 mppmod2.o libmpdiag2_h.o
	$(FC90) $(OPTS90) -o modmpdiag2.o -c modmpdiag2.f90

ompplib2.o : ompplib2.f90 modmpsort2.o modmpfft2.o modmpgard2.o mppmod2.o
	$(FC90) $(OPTS90) -o ompplib2.o -c ompplib2.f90

mpsimul2.o : mpsimul2.f03 $(ESMODS) modmpdiag2.o mppmod2.o
	$(FC90) $(OPTS90) -o mpsimul2.o -c $(FF03) mpsimul2.f03

mpbsimul2.o : mpbsimul2.f03 mpsimul2.o $(EMMODS) modmpdiag2.o 
	$(FC90) $(OPTS90) -o mpbsimul2.o -c $(FF03) mpbsimul2.f03

mpdsimul2.o : mpdsimul2.f03 mpbsimul2.o $(EDMODS) modmpdiag2.o 
	$(FC90) $(OPTS90) -o mpdsimul2.o -c $(FF03) mpdsimul2.f03

mpbeps2.o : mpbeps2.f90 mpsimul2.o mppmod2.o omplib.o ompplib2.o modpgraf2.o \
            modpgraf1.o
	$(FC90) $(OPTS90) -o mpbeps2.o -c mpbeps2.f90

mpbbeps2.o : mpbbeps2.f90 mpbsimul2.o mppmod2.o omplib.o ompplib2.o \
             modpgraf2.o modpgraf1.o
	$(FC90) $(OPTS90) -o mpbbeps2.o -c mpbbeps2.f90

mpdbeps2.o : mpdbeps2.f90 mpdsimul2.o mppmod2.o omplib.o ompplib2.o \
             modpgraf2.o modpgraf1.o
	$(FC90) $(OPTS90) -o mpdbeps2.o -c mpdbeps2.f90

#noMPI
nmodmpdiag2.o : modmpdiag2.f90 nmppmod2.o libmpdiag2_h.o
	$(FC90) $(OPTS90) -o nmodmpdiag2.o -c modmpdiag2.f90

msimul2.o : mpsimul2.f03 $(ESMODS) nmodmpdiag2.o nmppmod2.o
	$(FC90) $(OPTS90) -o msimul2.o -c $(FF03) mpsimul2.f03

mbsimul2.o : mpbsimul2.f03 msimul2.o $(EMMODS) nmodmpdiag2.o 
	$(FC90) $(OPTS90) -o mbsimul2.o -c $(FF03) mpbsimul2.f03

mdsimul2.o : mpdsimul2.f03 mbsimul2.o $(EDMODS) nmodmpdiag2.o
	$(FC90) $(OPTS90) -o mdsimul2.o -c $(FF03) mpdsimul2.f03

mbeps2.o : mpbeps2.f90 msimul2.o nmppmod2.o omplib.o nompplib2.o \
           nmodpgraf2.o modpgraf1.o
	$(FC90) $(OPTS90) -o mbeps2.o -c mpbeps2.f90

mbbeps2.o : mpbbeps2.f90 mbsimul2.o nmppmod2.o omplib.o nompplib2.o \
            nmodpgraf2.o modpgraf1.o
	$(FC90) $(OPTS90) -o mbbeps2.o -c mpbbeps2.f90

mdbeps2.o : mpdbeps2.f90 mdsimul2.o nmppmod2.o omplib.o nompplib2.o \
            nmodpgraf2.o modpgraf1.o
	$(FC90) $(OPTS90) -o mdbeps2.o -c mpdbeps2.f90


# Python

pydtimer.o : dtimer.c
	$(CC) $(PYCCOPTS) -o pydtimer.o -c dtimer.c

pympinit2.o : libmpinit2.f
	$(FC90) $(PYOPTS90) -o pympinit2.o -c libmpinit2.f

pymppush2.o : libmppush2.f
	$(FC90) $(PYOPTS90) -o pymppush2.o -c libmppush2.f

pyvmppush2.o : libvmppush2.f
	$(FC90) $(PYOPTS90) -o pyvmppush2.o -c libvmppush2.f

pympbpush2.o : libmpbpush2.f
	$(FC90) $(PYOPTS90) -o pympbpush2.o -c libmpbpush2.f

pyvmpbpush2.o : libvmpbpush2.f
	$(FC90) $(PYOPTS90) -o pyvmpbpush2.o -c libvmpbpush2.f

pympcurd2.o : libmpcurd2.f
	$(FC90) $(PYOPTS90) -o pympcurd2.o -c libmpcurd2.f

pyvmpcurd2.o : libvmpcurd2.f
	$(FC90) $(PYOPTS90) -o pyvmpcurd2.o -c libvmpcurd2.f

pympdpush2.o : libmpdpush2.f
	$(FC90) $(PYOPTS90) -o pympdpush2.o -c libmpdpush2.f

pyvmpdpush2.o : libvmpdpush2.f
	$(FC90) $(PYOPTS90) -o pyvmpdpush2.o -c libvmpdpush2.f

pympsort2.o : libmpsort2.f
	$(FC90) $(PYOPTS90) -o pympsort2.o -c libmpsort2.f

pyvmpsort2.o : libvmpsort2.f
	$(FC90) $(PYOPTS90) -o pyvmpsort2.o -c libvmpsort2.f

pympgard2.o : libmpgard2.f
	$(FC90) $(PYOPTS90) -o pympgard2.o -c libmpgard2.f

pympfft2.o : libvmpfft2.f
	$(FC90) $(PYOPTS90) -o pympfft2.o -c libvmpfft2.f

pympfield2.o : libmpfield2.f
	$(FC90) $(PYOPTS90) -o pympfield2.o -c libmpfield2.f

pympdiag2.o : libmpdiag2.f
	$(FC90) $(PYOPTS90) -o pympdiag2.o -c libmpdiag2.f

pynoMPI.o : noMPI.f
	$(FC90) $(PYOPTS90) -o pynoMPI.o -c noMPI.f

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

libmpush2.so : $(ESHOBJS) $(ESMODSRC) $(ESPYOBJS) $(MPMODSRC) pynoMPI.o \
                pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmpush2 \
        -c $(ESHOBJS) $(ESMODSRC) $(ESPYOBJS) $(MPMODSRC) pynoMPI.o \
           pydtimer.o

libmbpush2.so : $(EMHOBJS) $(EMMODSRC) $(EMPYOBJS) pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmbpush2 \
        -c $(EMHOBJS) $(EMMODSRC) $(EMPYOBJS) pydtimer.o

libmdpush2.so : $(EDHOBJS) $(EDMODSRC) $(EDPYOBJS) pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmdpush2 \
        -c $(EDHOBJS) $(EDMODSRC) $(EDPYOBJS) pydtimer.o

# Fortran Graphics

libt1.o : gks2/libt1.f90
	$(FC90) $(OPTS90) -o libt1.o -c gks2/libt1.f90

plot10.o : gks2/plot10.f
	$(FC90) $(OPTS90) -o plot10.o -c gks2/plot10.f

libplt10.o : gks2/libplt10.f
	$(FC90) $(OPTS90) -o libplt10.o -c gks2/libplt10.f
	
gks2/Ygl-4.2/libYgl.a :
	make -C gks2/Ygl-4.2
	
libygl.o : gks2/libygl.f gks2/Ygl-4.2/libYgl.a
	$(FC90) $(OPTS90) -o libygl.o -c gks2/libygl.f

libgks1.o : gks2/libgks1.f
	$(FC90) $(OPTS90) -o libgks1.o -c gks2/libgks1.f

nullgks1.o : gks2/nullgks1.f
	$(FC90) $(OPTS90) -o nullgks1.o -c gks2/nullgks1.f

libgks2.o : gks2/libgks2.f
	$(FC90) $(OPTS90) -o libgks2.o -c gks2/libgks2.f

nullgks2.o : gks2/nullgks2.f
	$(FC90) $(OPTS90) -o nullgks2.o -c gks2/nullgks2.f

plibgks2.o : gks2/plibgks2.f90 mpplib2.o
	$(FC90) $(OPTS90) -o plibgks2.o -c gks2/plibgks2.f90

nullpgks2.o : gks2/nullpgks2.f90
	$(FC90) $(OPTS90) -o nullpgks2.o -c gks2/nullpgks2.f90

modpgraf2.o : gks2/modpgraf2.f90 plibgks2.o
	$(FC90) $(OPTS90) -o modpgraf2.o -c gks2/modpgraf2.f90

modpgraf1.o : gks2/modpgraf1.f90
	$(FC90) $(OPTS90) -o modpgraf1.o -c gks2/modpgraf1.f90

#noMPI
nplibgks2.o : gks2/plibgks2.f90 nmpplib2.o
	$(FC90) $(OPTS90) -o nplibgks2.o -c gks2/plibgks2.f90

nmodpgraf2.o : gks2/modpgraf2.f90 nplibgks2.o
	$(FC90) $(OPTS90) -o nmodpgraf2.o -c gks2/modpgraf2.f90

nmodpgraf1.o : gks2/modpgraf1.f90 nplibgks2.o
	$(FC90) $(OPTS90) -o nmodpgraf1.o -c gks2/modpgraf1.f90

# Python Graphics

pymplot10.o : gks2/plot10.f
	$(FC90) $(PYOPTS90) -o pymplot10.o -c gks2/plot10.f

pymlibplt10.o : gks2/libplt10.f
	$(FC90) $(PYOPTS90) -o pymlibplt10.o -c sgks2/libplt10.f

pymlibygl.o : gks2/libygl.f gks2/Ygl-4.2/libYgl.a
	$(FC90) $(PYOPTS90) -o pymlibygl.o -c gks2/libygl.f

pymlibgks1.o : gks2/libgks1.f
	$(FC90) $(PYOPTS90) -o pymlibgks1.o -c gks2/libgks1.f

pymnullgks1.o : gks2/nullgks1.f
	$(FC90) $(PYOPTS90) -o pymnullgks1.o -c gks2/nullgks1.f

pymlibgks2.o : gks2/libgks2.f
	$(FC90) $(PYOPTS90) -o pymlibgks2.o -c gks2/libgks2.f

pymnullgks2.o : gks2/nullgks2.f
	$(FC90) $(PYOPTS90) -o pymnullgks2.o -c gks2/nullgks2.f

pymplibgks2.o : gks2/plibgks2.f90 nmpplib2.o
	$(FC90) $(PYOPTS90) -o pymplibgks2.o -c gks2/plibgks2.f90

pymnullpgks2.o : gks2/nullpgks2.f90
	$(FC90) $(PYOPTS90) -o pymnullpgks2.o -c gks2/nullpgks2.f90

fpgraf.so : $(GMODSRC) $(GPYOBJS)
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) $(LGOPTS) -m fpgraf \
        -c $(GMODSRC) $(GPYOBJS)

cpexecs : cpmpbeps2 cpmpbbeps2 cpmpdbeps2

cpmpbeps2 :
	cp mpbeps2 ../

cpmpbbeps2 :
	cp mpbbeps2 ../

cpmpdbeps2 :
	cp mpdbeps2 ../

#noMPI
cpmexecs : cpmbeps2 cpmbbeps2 cpmdbeps2

cpmbeps2 :
	cp mbeps2 ../

cpmbbeps2 :
	cp mbbeps2 ../

cpmdbeps2 :
	cp mdbeps2 ../

clean :
	rm -f *.o *.mod gks2/Ygl-4.2/*.o *.pyf

clobber: clean
	rm -f mpbeps2 mpbbeps2 mpdbeps2 mbeps2 mbbeps2 mdbeps2 gks2/Ygl-4.2/*.a *.so
