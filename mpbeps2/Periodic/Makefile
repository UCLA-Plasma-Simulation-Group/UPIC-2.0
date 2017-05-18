#Makefile for 2-1/2D Periodic Darwin Field MPI/OpenMP PIC test code

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
#LGOPTS = -lYgl -L/usr/X11R6/lib -lX11
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

# Fortran Graphics

#Tektronix
#GOBJS = libt1.o plot10.o libplt10.o libgks1.o libgks2.o plibgks2.o \
#modpgraf2.o
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libygl.o libgks1.o libgks2.o plibgks2.o modpgraf2.o
#NGOBJS = libygl.o libgks1.o libgks2.o nplibgks2.o nmodpgraf2.o
#No graphics
GOBJS = nullgks1.o nullgks2.o nullpgks2.o modpgraf2.o
NGOBJS = nullgks1.o nullgks2.o nullpgks2.o nmodpgraf2.o

# Linkage rules

all : mp_pfields2 

nompi : m_pfields2 

mp_pfields2 : mp_pfields2.o libmpfield2.o libmpfft2.o libmpgard2.o $(GOBJS) \
              dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LGOPTS) -o mp_pfields2 \
    mp_pfields2.o libmpfield2.o libmpfft2.o libmpgard2.o modmpfield2.o \
    modmpfft2.o modmpgard2.o mpplib2.o mppmod2.o omplib.o omppflib2.o \
    libmpfield2_h.o libmpfft2_h.o libmpgard2_h.o $(GOBJS) dtimer.o

#noMPI
m_pfields2 : m_pfields2.o libmpfield2.o libmpfft2.o libmpgard2.o noMPI.o \
             $(NGOBJS) dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) $(LGOPTS) -o m_pfields2 \
    m_pfields2.o libmpfield2.o libmpfft2.o libmpgard2.o modmpfield2.o \
    modmpfft2.o modmpgard2.o nmpplib2.o nmppmod2.o omplib.o nomppflib2.o \
    libmpfield2_h.o libmpfft2_h.o libmpgard2_h.o noMPI.o $(NGOBJS) dtimer.o

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

nomppflib2.o : omppflib2.f90 modmpfft2.o modmpgard2.o nmppmod2.o
	$(FC90) $(OPTS90) -o nomppflib2.o -c omppflib2.f90

#OPENMP
omplib.o : omplib.f90
	$(FC90) $(OPTS90) -o omplib.o -c omplib.f90

mpplib2.o : mpplib2.f90
	$(MPIFC) $(OPTS90) -o mpplib2.o -c mpplib2.f90

mppmod2.o : mppmod2.f90 mpplib2.o
	$(FC90) $(OPTS90) -o mppmod2.o -c mppmod2.f90

libmpfield2.o : libmpfield2.f
	$(FC90) $(OPTS90) -o libmpfield2.o -c libmpfield2.f

libmpfft2.o : libmpfft2.f
	$(FC90) $(OPTS90) -o libmpfft2.o -c libmpfft2.f

libmpgard2.o : libmpgard2.f
	$(FC90) $(OPTS90) -o libmpgard2.o -c libmpgard2.f

libmpfield2_h.o : libmpfield2_h.f90
	$(FC90) $(OPTS90) -o libmpfield2_h.o -c libmpfield2_h.f90

libmpfft2_h.o : libmpfft2_h.f90
	$(FC90) $(OPTS90) -o libmpfft2_h.o -c libmpfft2_h.f90

libmpgard2_h.o : libmpgard2_h.f90
	$(FC90) $(OPTS90) -o libmpgard2_h.o -c libmpgard2_h.f90

modmpfield2.o : modmpfield2.f90 libmpfield2_h.o
	$(FC90) $(OPTS90) -o modmpfield2.o -c modmpfield2.f90

modmpfft2.o : modmpfft2.f90 libmpfft2_h.o
	$(FC90) $(OPTS90) -o modmpfft2.o -c modmpfft2.f90

modmpgard2.o : modmpgard2.f90 libmpgard2_h.o
	$(FC90) $(OPTS90) -o modmpgard2.o -c modmpgard2.f90

omppflib2.o : omppflib2.f90 modmpfft2.o modmpgard2.o mppmod2.o
	$(FC90) $(OPTS90) -o omppflib2.o -c omppflib2.f90

mp_pfields2.o : mp_pfields2.f90 modmpfield2.o mppmod2.o omplib.o omppflib2.o \
                modpgraf2.o
	$(FC90) $(OPTS90) -o mp_pfields2.o -c mp_pfields2.f90

#noMPI
m_pfields2.o : mp_pfields2.f90 modmpfield2.o nmppmod2.o omplib.o nomppflib2.o \
               nmodpgraf2.o
	$(FC90) $(OPTS90) -o m_pfields2.o -c mp_pfields2.f90

# Fortran Graphics

libt1.o : gks2/libt1.f90
	$(FC90) $(OPTS90) -o libt1.o -c gks2/libt1.f90

plot10.o : gks2/plot10.f
	$(FC90) $(OPTS90) -o plot10.o -c gks2/plot10.f

libplt10.o : gks2/libplt10.f
	$(FC90) $(OPTS90) -o libplt10.o -c gks2/libplt10.f
	
libygl.o : gks2/libygl.f
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

#noMPI
nplibgks2.o : gks2/plibgks2.f90 nmpplib2.o
	$(FC90) $(OPTS90) -o nplibgks2.o -c gks2/plibgks2.f90

nmodpgraf2.o : gks2/modpgraf2.f90 nplibgks2.o
	$(FC90) $(OPTS90) -o nmodpgraf2.o -c gks2/modpgraf2.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f mp_pfields2 m_pfields2
