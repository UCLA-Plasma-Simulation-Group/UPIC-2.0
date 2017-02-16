#Makefile for 3D Periodic Darwin Field MPI/OpenMP PIC test code

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

#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
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
LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
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

# Linkage rules

all : mp_dfields3 

nompi : m_dfields3 

mp_dfields3 : mp_dfields3.o libmpfield3.o libmpfft3.o libmpinit3.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o mp_dfields3 \
    mp_dfields3.o libmpfield3.o libmpfft3.o libmpinit3.o modmpfield3.o \
    modmpfft3.o mpplib3.o mppmod3.o omplib.o omppflib3.o libmpfield3_h.o \
    libmpfft3_h.o libmpinit3_h.o dtimer.o

#noMPI
m_dfields3 : m_dfields3.o libmpfield3.o libmpfft3.o libmpinit3.o noMPI.o \
             dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o m_dfields3 \
    m_dfields3.o libmpfield3.o libmpfft3.o libmpinit3.o modmpfield3.o \
    modmpfft3.o nmpplib3.o nmppmod3.o omplib.o nomppflib3.o libmpfield3_h.o \
    libmpfft3_h.o libmpinit3_h.o noMPI.o dtimer.o

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

nomppflib3.o : omppflib3.f90 modmpfft3.o nmppmod3.o
	$(FC90) $(OPTS90) -o nomppflib3.o -c omppflib3.f90

#OPENMP
omplib.o : omplib.f90
	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

mpplib3.o : mpplib3.f90
	$(MPIFC) $(OPTS90) -o mpplib3.o -c mpplib3.f90

mppmod3.o : mppmod3.f90 mpplib3.o
	$(FC90) $(OPTS90) -o mppmod3.o -c mppmod3.f90

libmpfield3.o : libmpfield3.f
	$(FC90) $(OPTS90) -o libmpfield3.o -c libmpfield3.f

libmpfft3.o : libmpfft3.f
	$(FC90) $(OPTS90) -o libmpfft3.o -c libmpfft3.f

libmpinit3.o : libmpinit3.f
	$(FC90) $(OPTS90) -o libmpinit3.o -c libmpinit3.f

libmpfield3_h.o : libmpfield3_h.f90
	$(FC90) $(OPTS90) -o libmpfield3_h.o -c libmpfield3_h.f90

libmpfft3_h.o : libmpfft3_h.f90
	$(FC90) $(OPTS90) -o libmpfft3_h.o -c libmpfft3_h.f90

libmpinit3_h.o : libmpinit3_h.f90
	$(FC90) $(OPTS90) -o libmpinit3_h.o -c libmpinit3_h.f90

modmpfield3.o : modmpfield3.f90 libmpfield3_h.o
	$(FC90) $(OPTS90) -o modmpfield3.o -c modmpfield3.f90

modmpfft3.o : modmpfft3.f90 libmpfft3_h.o
	$(FC90) $(OPTS90) -o modmpfft3.o -c modmpfft3.f90

omppflib3.o : omppflib3.f90 modmpfft3.o mppmod3.o
	$(FC90) $(OPTS90) -o omppflib3.o -c omppflib3.f90

mp_dfields3.o : mp_dfields3.f90 modmpfield3.o mppmod3.o omplib.o omppflib3.o \
                libmpinit3_h.o
	$(FC90) $(OPTS90) -o mp_dfields3.o -c mp_dfields3.f90

#noMPI
m_dfields3.o : mp_dfields3.f90 modmpfield3.o nmppmod3.o omplib.o nomppflib3.o \
               libmpinit3_h.o
	$(FC90) $(OPTS90) -o m_dfields3.o -c mp_dfields3.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f mp_dfields3 m_dfields3
