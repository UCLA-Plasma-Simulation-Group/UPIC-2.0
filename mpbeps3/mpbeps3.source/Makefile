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

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

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

ESOBJS = libmpinit3.o libmppush3.o libmpsort3.o libmpgard3.o libmpfft3.o \
libmpfield3.o libmpdiag3.o

EMOBJS = libmpbpush3.o libmpcurd3.o

EDOBJS = libmpdpush3.o

ESHOBJS = libmpinit3_h.o libmppush3_h.o libmpsort3_h.o libmpgard3_h.o \
libmpfft3_h.o libmpfield3_h.o libmpdiag3_h.o

EMHOBJS = libmpbpush3_h.o libmpcurd3_h.o

EDHOBJS = libmpdpush3_h.o

ESMODS = input3mod.o modmpinit3.o modmppush3.o modmpsort3.o modmpgard3.o \
modmpfft3.o modmpfield3.o modmpdiag3.o

EMMODS = modmpbpush3.o modmpcurd3.o

EDMODS = modmpdpush3.o

# Linkage rules

all : mpbeps3 mpbbeps3 mpdbeps3

nompi : mbeps3 mbbeps3 mdbeps3

mpbeps3 : mpbeps3.o $(ESOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpbeps3 \
    mpbeps3.o $(ESOBJS) $(ESMODS) mpplib3.o mppmod3.o omplib.o ompplib3.o \
    $(ESHOBJS) mpsimul3.o dtimer.o

mpbbeps3 : mpbbeps3.o $(ESOBJS) $(EMOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpbbeps3 \
    mpbbeps3.o $(ESOBJS) $(EMOBJS) $(ESMODS) $(EMMODS) mpplib3.o mppmod3.o \
    omplib.o ompplib3.o $(ESHOBJS) $(EMHOBJS) mpsimul3.o mpbsimul3.o dtimer.o

mpdbeps3 : mpdbeps3.o $(ESOBJS) $(EMOBJS) $(EDOBJS) dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mpdbeps3 \
    mpdbeps3.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(ESMODS) $(EMMODS) $(EDMODS) \
    mpplib3.o mppmod3.o omplib.o ompplib3.o $(ESHOBJS) $(EMHOBJS) $(EDHOBJS) \
    mpsimul3.o mpbsimul3.o mpdsimul3.o dtimer.o

#noMPI
mbeps3 : mbeps3.o $(ESOBJS) noMPI.o dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mbeps3 \
    mbeps3.o $(ESOBJS) $(ESMODS) nmpplib3.o nmppmod3.o omplib.o nompplib3.o \
    $(ESHOBJS) noMPI.o msimul3.o dtimer.o

mbbeps3 : mbbeps3.o $(ESOBJS) $(EMOBJS) noMPI.o dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mbbeps3 \
    mbbeps3.o $(ESOBJS) $(EMOBJS) $(ESMODS) $(EMMODS) nmpplib3.o nmppmod3.o \
    omplib.o nompplib3.o $(ESHOBJS) $(EMHOBJS) noMPI.o msimul3.o mbsimul3.o \
    dtimer.o

mdbeps3 : mdbeps3.o $(ESOBJS) $(EMOBJS) $(EDOBJS) noMPI.o dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o mdbeps3 \
    mdbeps3.o $(ESOBJS) $(EMOBJS) $(EDOBJS) $(ESMODS) $(EMMODS) $(EDMODS) \
    nmpplib3.o nmppmod3.o omplib.o nompplib3.o $(ESHOBJS) $(EMHOBJS) \
    $(EDHOBJS) noMPI.o msimul3.o mbsimul3.o mdsimul3.o dtimer.o

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

libmpsort3.o : libmpsort3.f
	$(FC90) $(OPTS90) -o libmpsort3.o -c libmpsort3.f

libmpgard3.o : libmpgard3.f
	$(FC90) $(OPTS90) -o libmpgard3.o -c libmpgard3.f

libmpfft3.o : libmpfft3.f
	$(FC90) $(OPTS90) -o libmpfft3.o -c libmpfft3.f

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

libmpsort3_h.o : libmpsort3_h.f90
	$(FC90) $(OPTS90) -o libmpsort3_h.o -c libmpsort3_h.f90

libmpgard3_h.o : libmpgard3_h.f90
	$(FC90) $(OPTS90) -o libmpgard3_h.o -c libmpgard3_h.f90

libmpfft3_h.o : libmpfft3_h.f90
	$(FC90) $(OPTS90) -o libmpfft3_h.o -c libmpfft3_h.f90

libmpfield3_h.o : libmpfield3_h.f90
	$(FC90) $(OPTS90) -o libmpfield3_h.o -c libmpfield3_h.f90

libmpdiag3_h.o : libmpdiag3_h.f90
	$(FC90) $(OPTS90) -o libmpdiag3_h.o -c libmpdiag3_h.f90

input3mod.o : input3mod.f90
	$(FC90) $(OPTS90) -o input3mod.o -c input3mod.f90

modmpinit3.o : modmpinit3.f90 libmpinit3_h.o
	$(FC90) $(OPTS90) -o modmpinit3.o -c modmpinit3.f90

modmppush3.o : modmppush3.f90 libmppush3_h.o
	$(FC90) $(OPTS90) -o modmppush3.o -c modmppush3.f90

modmpbpush3.o : modmpbpush3.f90 libmpbpush3_h.o
	$(FC90) $(OPTS90) -o modmpbpush3.o -c modmpbpush3.f90

modmpcurd3.o : modmpcurd3.f90 libmpcurd3_h.o
	$(FC90) $(OPTS90) -o modmpcurd3.o -c modmpcurd3.f90

modmpdpush3.o : modmpdpush3.f90 libmpdpush3_h.o
	$(FC90) $(OPTS90) -o modmpdpush3.o -c modmpdpush3.f90

modmpsort3.o : modmpsort3.f90 libmpsort3_h.o
	$(FC90) $(OPTS90) -o modmpsort3.o -c modmpsort3.f90

modmpgard3.o : modmpgard3.f90 libmpgard3_h.o
	$(FC90) $(OPTS90) -o modmpgard3.o -c modmpgard3.f90

modmpfft3.o : modmpfft3.f90 libmpfft3_h.o
	$(FC90) $(OPTS90) -o modmpfft3.o -c modmpfft3.f90

modmpfield3.o : modmpfield3.f90 libmpfield3_h.o
	$(FC90) $(OPTS90) -o modmpfield3.o -c modmpfield3.f90

modmpdiag3.o : modmpdiag3.f90 mppmod3.o libmpdiag3_h.o
	$(FC90) $(OPTS90) -o modmpdiag3.o -c modmpdiag3.f90

ompplib3.o : ompplib3.f90 modmpsort3.o modmpfft3.o modmpgard3.o mppmod3.o
	$(FC90) $(OPTS90) -o ompplib3.o -c ompplib3.f90

mpsimul3.o : mpsimul3.f03 $(ESMODS) mppmod3.o
	$(FC90) $(OPTS90) -o mpsimul3.o -c $(FF03) mpsimul3.f03

mpbsimul3.o : mpbsimul3.f03 mpsimul3.o $(EMMODS)
	$(FC90) $(OPTS90) -o mpbsimul3.o -c $(FF03) mpbsimul3.f03

mpdsimul3.o : mpdsimul3.f03 mpbsimul3.o $(EDMODS)
	$(FC90) $(OPTS90) -o mpdsimul3.o -c $(FF03) mpdsimul3.f03

mpbeps3.o : mpbeps3.f90 mpsimul3.o mppmod3.o omplib.o ompplib3.o
	$(FC90) $(OPTS90) -o mpbeps3.o -c mpbeps3.f90

mpbbeps3.o : mpbbeps3.f90 mpbsimul3.o mppmod3.o omplib.o ompplib3.o
	$(FC90) $(OPTS90) -o mpbbeps3.o -c mpbbeps3.f90

mpdbeps3.o : mpdbeps3.f90 mpdsimul3.o mppmod3.o omplib.o ompplib3.o
	$(FC90) $(OPTS90) -o mpdbeps3.o -c mpdbeps3.f90

#noMPI
msimul3.o : mpsimul3.f03 $(ESMODS) nmppmod3.o
	$(FC90) $(OPTS90) -o msimul3.o -c $(FF03) mpsimul3.f03

mbsimul3.o : mpbsimul3.f03 msimul3.o $(EMMODS)
	$(FC90) $(OPTS90) -o mbsimul3.o -c $(FF03) mpbsimul3.f03

mdsimul3.o : mpdsimul3.f03 mbsimul3.o $(EDMODS)
	$(FC90) $(OPTS90) -o mdsimul3.o -c $(FF03) mpdsimul3.f03

mbeps3.o : mpbeps3.f90 msimul3.o nmppmod3.o omplib.o nompplib3.o
	$(FC90) $(OPTS90) -o mbeps3.o -c mpbeps3.f90

mbbeps3.o : mpbbeps3.f90 mbsimul3.o nmppmod3.o omplib.o nompplib3.o
	$(FC90) $(OPTS90) -o mbbeps3.o -c mpbbeps3.f90

mdbeps3.o : mpdbeps3.f90 mdsimul3.o nmppmod3.o omplib.o nompplib3.o
	$(FC90) $(OPTS90) -o mdbeps3.o -c mpdbeps3.f90

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
	rm -f *.o *.mod

clobber: clean
	rm -f mpbeps3 mpbbeps3 mpdbeps3 mbeps3 mbbeps3 mdbeps3
