#Makefile for 3D MPI/OpenMP PIC codes

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

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99

#FF03 = -Tf

#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMPI
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

ESOBJS = libmpinit3.o libmppush3.o libvmppush3.o libmpsort3.o \
libvmpsort3.o libmpgard3.o libvmpfft3.o libmpfield3.o libmpdiag3.o

EMOBJS = libmpbpush3.o libvmpbpush3.o libmpcurd3.o libvmpcurd3.o

EDOBJS = libmpdpush3.o libvmpdpush3.o

ESHOBJS = libmpinit3_h.o libmppush3_h.o libvmppush3_h.o libmpsort3_h.o \
libvmpsort3_h.o libmpgard3_h.o libvmpfft3_h.o libmpfield3_h.o libmpdiag3_h.o

EMHOBJS = libmpbpush3_h.o libvmpbpush3_h.o libmpcurd3_h.o libvmpcurd3_h.o

EDHOBJS = libmpdpush3_h.o libvmpdpush3_h.o

ESMODS = input3mod.o modmpinit3.o modmppush3.o modmpsort3.o modmpgard3.o \
modmpfft3.o modmpfield3.o

EMMODS = modmpbpush3.o modmpcurd3.o

EDMODS = modmpdpush3.o

# Python

ESMODSRC =  ../ReadFiles/input3mod.f90 modmpinit3.f90 modmppush3.f90 modmpsort3.f90 \
modmpgard3.f90 modmpfft3.f90 modmpfield3.f90

EMMODSRC = modmpbpush3.f90 modmpcurd3.f90

EDMODSRC = modmpdpush3.f90

MPMODSRC = noMPI.f90 mpplib3.f90 mppmod3.f90 modmpdiag3.f90 omplib.f90 \
ompplib3.f90

ESPYOBJS = pympinit3.o pymppush3.o pyvmppush3.o pympsort3.o pyvmpsort3.o \
pympgard3.o pympfft3.o pympfield3.o pympdiag3.o

EMPYOBJS = pympbpush3.o pyvmpbpush3.o pympcurd3.o pyvmpcurd3.o

EDPYOBJS = pympdpush3.o pyvmpdpush3.o

# Linkage rules

all : mpbeps3 mpbbeps3 mpdbeps3

nompi : mbeps3 mbbeps3 mdbeps3

python : mbeps3.py mbbeps3.py mdbeps3.py

mpbeps3 : mpbeps3.o $(ESOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpbeps3 \
    mpbeps3.o $(ESOBJS) $(ESMODS) modmpdiag3.o mpplib3.o mppmod3.o omplib.o \
    ompplib3.o $(ESHOBJS) mpsimul3.o dtimer.o

mpbbeps3 : mpbbeps3.o $(ESOBJS) $(EMOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpbbeps3 \
    mpbbeps3.o $(ESOBJS) $(EMOBJS) $(ESMODS) $(EMMODS) modmpdiag3.o mpplib3.o \
    mppmod3.o omplib.o ompplib3.o $(ESHOBJS) $(EMHOBJS) mpsimul3.o mpbsimul3.o \
    dtimer.o

mpdbeps3 : mpdbeps3.o $(ESOBJS) $(EMOBJS) $(EDOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpdbeps3 \
    mpdbeps3.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(ESMODS) $(EMMODS) $(EDMODS) \
    modmpdiag3.o mpplib3.o mppmod3.o omplib.o ompplib3.o $(ESHOBJS) $(EMHOBJS) \
    $(EDHOBJS) mpsimul3.o mpbsimul3.o mpdsimul3.o dtimer.o

#noMPI
mbeps3 : mbeps3.o $(ESOBJS) noMPI.o dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mbeps3 \
    mbeps3.o $(ESOBJS) $(ESMODS) nmodmpdiag3.o nmpplib3.o nmppmod3.o omplib.o \
    nompplib3.o $(ESHOBJS) noMPI.o msimul3.o dtimer.o

mbbeps3 : mbbeps3.o $(ESOBJS) $(EMOBJS) noMPI.o dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mbbeps3 \
    mbbeps3.o $(ESOBJS) $(EMOBJS) $(ESMODS) $(EMMODS) nmodmpdiag3.o nmpplib3.o \
    nmppmod3.o omplib.o nompplib3.o $(ESHOBJS) $(EMHOBJS) noMPI.o msimul3.o \
    mbsimul3.o dtimer.o

mdbeps3 : mdbeps3.o $(ESOBJS) $(EMOBJS) $(EDOBJS) noMPI.o dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mdbeps3 \
    mdbeps3.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(ESMODS) $(EMMODS) $(EDMODS) \
    nmodmpdiag3.o nmpplib3.o nmppmod3.o omplib.o nompplib3.o $(ESHOBJS) \
    $(EMHOBJS) $(EDHOBJS) noMPI.o msimul3.o mbsimul3.o mdsimul3.o dtimer.o

# Python

mbeps3.py : libmpush3.so dtimer.so

mbbeps3.py : libmbpush3.so

mdbeps3.py : libmdpush3.so

list:
	echo mpbeps3 mpbbeps3 mpdbeps3 nompi mbeps3 mbbeps3 mdbeps3 clean clobber

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#noMPI
noMPImod.o : noMPI.f90
	$(FC90) $(OPTS90) -o noMPImod.o -c noMPI.f90

noMPI.o : noMPI.f
	$(FC90) $(OPTS90) -o noMPI.o -c noMPI.f

nmpplib3.o : mpplib3.f90 noMPImod.o
	$(FC90) $(OPTS90) -o nmpplib3.o -c mpplib3.f90

nmppmod3.o : mppmod3.f90 nmpplib3.o
	$(FC90) $(OPTS90) -o nmppmod3.o -c mppmod3.f90

nompplib3.o : ompplib3.f90 modmpsort3.o modmpfft3.o modmpgard3.o nmppmod3.o
	$(FC90) $(OPTS90) -o nompplib3.o -c ompplib3.f90

#OPENMP
omplib.o : omplib.f90
	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

mpplib3.o : mpplib3.f90
	$(MPIFC) $(OPTS90) -o mpplib3.o -c mpplib3.f90

mppmod3.o : mppmod3.f90 mpplib3.o
	$(FC90) $(OPTS90) -o mppmod3.o -c mppmod3.f90

libmpinit3.o : libmpinit3.f
	$(FC90) $(OPTS90) -o libmpinit3.o -c libmpinit3.f

libmppush3.o : libmppush3.f
	$(FC90) $(OPTS90) -o libmppush3.o -c libmppush3.f

libmpbpush3.o : libmpbpush3.f
	$(FC90) $(OPTS90) -o libmpbpush3.o -c libmpbpush3.f

libmpcurd3.o : libmpcurd3.f
	$(FC90) $(OPTS90) -o libmpcurd3.o -c libmpcurd3.f

libmpdpush3.o : libmpdpush3.f
	$(FC90) $(OPTS90) -o libmpdpush3.o -c libmpdpush3.f

libvmppush3.o : libvmppush3.f
	$(FC90) $(OPTS90) -o libvmppush3.o -c libvmppush3.f

libvmpbpush3.o : libvmpbpush3.f
	$(FC90) $(OPTS90) -o libvmpbpush3.o -c libvmpbpush3.f

libvmpcurd3.o : libvmpcurd3.f
	$(FC90) $(OPTS90) -o libvmpcurd3.o -c libvmpcurd3.f

libvmpdpush3.o : libvmpdpush3.f
	$(FC90) $(OPTS90) -o libvmpdpush3.o -c libvmpdpush3.f

libmpsort3.o : libmpsort3.f
	$(FC90) $(OPTS90) -o libmpsort3.o -c libmpsort3.f

libvmpsort3.o : libvmpsort3.f
	$(FC90) $(OPTS90) -o libvmpsort3.o -c libvmpsort3.f

libmpgard3.o : libmpgard3.f
	$(FC90) $(OPTS90) -o libmpgard3.o -c libmpgard3.f

libvmpfft3.o : libvmpfft3.f
	$(FC90) $(OPTS90) -o libvmpfft3.o -c libvmpfft3.f

libmpfield3.o : libmpfield3.f
	$(FC90) $(OPTS90) -o libmpfield3.o -c libmpfield3.f

libmpdiag3.o : libmpdiag3.f
	$(FC90) $(OPTS90) -o libmpdiag3.o -c libmpdiag3.f

libmpinit3_h.o : libmpinit3_h.f90
	$(FC90) $(OPTS90) -o libmpinit3_h.o -c libmpinit3_h.f90

libmppush3_h.o : libmppush3_h.f90
	$(FC90) $(OPTS90) -o libmppush3_h.o -c libmppush3_h.f90

libmpbpush3_h.o : libmpbpush3_h.f90
	$(FC90) $(OPTS90) -o libmpbpush3_h.o -c libmpbpush3_h.f90

libmpcurd3_h.o : libmpcurd3_h.f90
	$(FC90) $(OPTS90) -o libmpcurd3_h.o -c libmpcurd3_h.f90

libmpdpush3_h.o : libmpdpush3_h.f90
	$(FC90) $(OPTS90) -o libmpdpush3_h.o -c libmpdpush3_h.f90

libvmppush3_h.o : libvmppush3_h.f90
	$(FC90) $(OPTS90) -o libvmppush3_h.o -c libvmppush3_h.f90

libvmpbpush3_h.o : libvmpbpush3_h.f90
	$(FC90) $(OPTS90) -o libvmpbpush3_h.o -c libvmpbpush3_h.f90

libvmpcurd3_h.o : libvmpcurd3_h.f90
	$(FC90) $(OPTS90) -o libvmpcurd3_h.o -c libvmpcurd3_h.f90

libvmpdpush3_h.o : libvmpdpush3_h.f90
	$(FC90) $(OPTS90) -o libvmpdpush3_h.o -c libvmpdpush3_h.f90

libmpsort3_h.o : libmpsort3_h.f90
	$(FC90) $(OPTS90) -o libmpsort3_h.o -c libmpsort3_h.f90

libvmpsort3_h.o : libvmpsort3_h.f90
	$(FC90) $(OPTS90) -o libvmpsort3_h.o -c libvmpsort3_h.f90

libmpgard3_h.o : libmpgard3_h.f90
	$(FC90) $(OPTS90) -o libmpgard3_h.o -c libmpgard3_h.f90

libvmpfft3_h.o : libvmpfft3_h.f90
	$(FC90) $(OPTS90) -o libvmpfft3_h.o -c libvmpfft3_h.f90

libmpfield3_h.o : libmpfield3_h.f90
	$(FC90) $(OPTS90) -o libmpfield3_h.o -c libmpfield3_h.f90

libmpdiag3_h.o : libmpdiag3_h.f90
	$(FC90) $(OPTS90) -o libmpdiag3_h.o -c libmpdiag3_h.f90

input3mod.o : input3mod.f90
	$(FC90) $(OPTS90) -o input3mod.o -c input3mod.f90

modmpinit3.o : modmpinit3.f90 libmpinit3_h.o
	$(FC90) $(OPTS90) -o modmpinit3.o -c modmpinit3.f90

modmppush3.o : modmppush3.f90 libmppush3_h.o libvmppush3_h.o
	$(FC90) $(OPTS90) -o modmppush3.o -c modmppush3.f90

modmpbpush3.o : modmpbpush3.f90 libmpbpush3_h.o libvmpbpush3_h.o
	$(FC90) $(OPTS90) -o modmpbpush3.o -c modmpbpush3.f90

modmpcurd3.o : modmpcurd3.f90 libmpcurd3_h.o libvmpcurd3_h.o
	$(FC90) $(OPTS90) -o modmpcurd3.o -c modmpcurd3.f90

modmpdpush3.o : modmpdpush3.f90 libmpdpush3_h.o libvmpdpush3_h.o
	$(FC90) $(OPTS90) -o modmpdpush3.o -c modmpdpush3.f90

modmpsort3.o : modmpsort3.f90 libmpsort3_h.o libvmpsort3_h.o
	$(FC90) $(OPTS90) -o modmpsort3.o -c modmpsort3.f90

modmpgard3.o : modmpgard3.f90 libmpgard3_h.o
	$(FC90) $(OPTS90) -o modmpgard3.o -c modmpgard3.f90

modmpfft3.o : modmpfft3.f90 libvmpfft3_h.o
	$(FC90) $(OPTS90) -o modmpfft3.o -c modmpfft3.f90

modmpfield3.o : modmpfield3.f90 libmpfield3_h.o
	$(FC90) $(OPTS90) -o modmpfield3.o -c modmpfield3.f90

modmpdiag3.o : modmpdiag3.f90 mppmod3.o libmpdiag3_h.o
	$(FC90) $(OPTS90) -o modmpdiag3.o -c modmpdiag3.f90

ompplib3.o : ompplib3.f90 modmpsort3.o modmpfft3.o modmpgard3.o mppmod3.o
	$(FC90) $(OPTS90) -o ompplib3.o -c ompplib3.f90

mpsimul3.o : mpsimul3.f03 $(ESMODS) modmpdiag3.o mppmod3.o
	$(FC90) $(OPTS90) -o mpsimul3.o -c $(FF03) mpsimul3.f03

mpbsimul3.o : mpbsimul3.f03 mpsimul3.o $(EMMODS) modmpdiag3.o 
	$(FC90) $(OPTS90) -o mpbsimul3.o -c $(FF03) mpbsimul3.f03

mpdsimul3.o : mpdsimul3.f03 mpbsimul3.o $(EDMODS) modmpdiag3.o 
	$(FC90) $(OPTS90) -o mpdsimul3.o -c $(FF03) mpdsimul3.f03

mpbeps3.o : mpbeps3.f90 mpsimul3.o mppmod3.o omplib.o ompplib3.o
	$(FC90) $(OPTS90) -o mpbeps3.o -c mpbeps3.f90

mpbbeps3.o : mpbbeps3.f90 mpbsimul3.o mppmod3.o omplib.o ompplib3.o
	$(FC90) $(OPTS90) -o mpbbeps3.o -c mpbbeps3.f90

mpdbeps3.o : mpdbeps3.f90 mpdsimul3.o mppmod3.o omplib.o ompplib3.o
	$(FC90) $(OPTS90) -o mpdbeps3.o -c mpdbeps3.f90

#noMPI
nmodmpdiag3.o : modmpdiag3.f90 nmppmod3.o libmpdiag3_h.o
	$(FC90) $(OPTS90) -o nmodmpdiag3.o -c modmpdiag3.f90

msimul3.o : mpsimul3.f03 $(ESMODS) nmodmpdiag3.o nmppmod3.o
	$(FC90) $(OPTS90) -o msimul3.o -c $(FF03) mpsimul3.f03

mbsimul3.o : mpbsimul3.f03 msimul3.o $(EMMODS) nmodmpdiag3.o 
	$(FC90) $(OPTS90) -o mbsimul3.o -c $(FF03) mpbsimul3.f03

mdsimul3.o : mpdsimul3.f03 mbsimul3.o $(EDMODS) nmodmpdiag3.o 
	$(FC90) $(OPTS90) -o mdsimul3.o -c $(FF03) mpdsimul3.f03

mbeps3.o : mpbeps3.f90 msimul3.o nmppmod3.o omplib.o nompplib3.o
	$(FC90) $(OPTS90) -o mbeps3.o -c mpbeps3.f90

mbbeps3.o : mpbbeps3.f90 mbsimul3.o nmppmod3.o omplib.o nompplib3.o
	$(FC90) $(OPTS90) -o mbbeps3.o -c mpbbeps3.f90

mdbeps3.o : mpdbeps3.f90 mdsimul3.o nmppmod3.o omplib.o nompplib3.o
	$(FC90) $(OPTS90) -o mdbeps3.o -c mpdbeps3.f90

# Python

pydtimer.o : dtimer.c
	$(CC) $(PYCCOPTS) -o pydtimer.o -c dtimer.c

pympinit3.o : libmpinit3.f
	$(FC90) $(PYOPTS90) -o pympinit3.o -c libmpinit3.f

pymppush3.o : libmppush3.f
	$(FC90) $(PYOPTS90) -o pymppush3.o -c libmppush3.f

pyvmppush3.o : libvmppush3.f
	$(FC90) $(PYOPTS90) -o pyvmppush3.o -c libvmppush3.f

pympbpush3.o : libmpbpush3.f
	$(FC90) $(PYOPTS90) -o pympbpush3.o -c libmpbpush3.f

pyvmpbpush3.o : libvmpbpush3.f
	$(FC90) $(PYOPTS90) -o pyvmpbpush3.o -c libvmpbpush3.f

pympcurd3.o : libmpcurd3.f
	$(FC90) $(PYOPTS90) -o pympcurd3.o -c libmpcurd3.f

pyvmpcurd3.o : libvmpcurd3.f
	$(FC90) $(PYOPTS90) -o pyvmpcurd3.o -c libvmpcurd3.f

pympdpush3.o : libmpdpush3.f
	$(FC90) $(PYOPTS90) -o pympdpush3.o -c libmpdpush3.f

pyvmpdpush3.o : libvmpdpush3.f
	$(FC90) $(PYOPTS90) -o pyvmpdpush3.o -c libvmpdpush3.f

pympsort3.o : libmpsort3.f
	$(FC90) $(PYOPTS90) -o pympsort3.o -c libmpsort3.f

pyvmpsort3.o : libvmpsort3.f
	$(FC90) $(PYOPTS90) -o pyvmpsort3.o -c libvmpsort3.f

pympgard3.o : libmpgard3.f
	$(FC90) $(PYOPTS90) -o pympgard3.o -c libmpgard3.f

pympfft3.o : libvmpfft3.f
	$(FC90) $(PYOPTS90) -o pympfft3.o -c libvmpfft3.f

pympfield3.o : libmpfield3.f
	$(FC90) $(PYOPTS90) -o pympfield3.o -c libmpfield3.f

pympdiag3.o : libmpdiag3.f
	$(FC90) $(PYOPTS90) -o pympdiag3.o -c libmpdiag3.f

pynoMPI.o : noMPI.f
	$(FC90) $(PYOPTS90) -o pynoMPI.o -c noMPI.f

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

libmpush3.so : $(ESHOBJS) $(ESMODSRC) $(ESPYOBJS) $(MPMODSRC) pynoMPI.o \
                pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmpush3 \
        -c $(ESHOBJS) $(ESMODSRC) $(ESPYOBJS) $(MPMODSRC) pynoMPI.o \
           pydtimer.o

libmbpush3.so : $(EMHOBJS) $(EMMODSRC) $(EMPYOBJS) pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmbpush3 \
        -c $(EMHOBJS) $(EMMODSRC) $(EMPYOBJS) pydtimer.o

libmdpush3.so : $(EDHOBJS) $(EDMODSRC) $(EDPYOBJS) pydtimer.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m libmdpush3 \
        -c $(EDHOBJS) $(EDMODSRC) $(EDPYOBJS) pydtimer.o

cpexecs : cpmpbeps3 cpmpbbeps3 cpmpdbeps3

cpmpbeps3 :
	cp mpbeps3 ../

cpmpbbeps3 :
	cp mpbbeps3 ../

cpmpdbeps3 :
	cp mpdbeps3 ../

#noMPI
cpmexecs : cpmbeps3 cpmbbeps3 cpmdbeps3

cpmbeps3 :
	cp mbeps3 ../

cpmbbeps3 :
	cp mbbeps3 ../

cpmdbeps3 :
	cp mdbeps3 ../

clean :
	rm -f *.o *.mod *.pyf

clobber: clean
	rm -f mpbeps3 mpbbeps3 mpdbeps3 mbeps3 mbbeps3 mdbeps3 *.so
