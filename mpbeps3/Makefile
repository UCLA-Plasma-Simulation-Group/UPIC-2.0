#Makefile for 3D MPI/OpenMP PIC codes in mpbeps3.source

# Linkage rules

#all : python

all: fortran

fbuild : fmpbeps3 fmpbbeps3 fmpdbeps3

fnompi : fmbeps3 fmbbeps3 fmdbeps3

fortran : fbuild
	make cpexecs -C mpbeps3.source

mpbeps3 : fmpbeps3
	make cpmpbeps3 -C mpbeps3.source

mpbbeps3 : fmpbbeps3
	make cpmpbbeps3 -C mpbeps3.source

mpdbeps3 : fmpdbeps3
	make cpmpdbeps3 -C mpbeps3.source

fmpbeps3 :
	make mpbeps3 -C mpbeps3.source

fmpbbeps3 :
	make mpbbeps3 -C mpbeps3.source

fmpdbeps3 :
	make mpdbeps3 -C mpbeps3.source

#noMPI
nompi : fnompi
	make cpmexecs -C mpbeps3.source

mbeps3 : fmbeps3
	make cpmbeps3 -C mpbeps3.source

mbbeps3 : fmbbeps3
	make cpmbbeps3 -C mpbeps3.source

mdbeps3 : fmdbeps3
	make cpmdbeps3 -C mpbeps3.source

fmbeps3 :
	make mbeps3 -C mpbeps3.source

fmbbeps3 :
	make mbbeps3 -C mpbeps3.source

fmdbeps3 :
	make mdbeps3 -C mpbeps3.source

#python :
#	make python -C mpbeps3.source

clean :
	make clean -C mpbeps3.source

clobber: rclobber
	rm -f mpbeps3 mpbbeps3 mpdbeps3 mbeps3 mbbeps3 mdbeps3

rclobber:
	make clobber -C mpbeps3.source
