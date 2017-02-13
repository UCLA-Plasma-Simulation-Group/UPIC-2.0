#Makefile for 2D MPI/OpenMP PIC codes in mpbeps2.source

# Linkage rules

#all : python

all: fortran

fbuild : fmpbeps2 fmpbbeps2 fmpdbeps2

fnompi : fmbeps2 fmbbeps2 fmdbeps2

fortran : fbuild
	make cpexecs -C mpbeps2.source

mpbeps2 : fmpbeps2
	make cpmpbeps2 -C mpbeps2.source

mpbbeps2 : fmpbbeps2
	make cpmpbbeps2 -C mpbeps2.source

mpdbeps2 : fmpdbeps2
	make cpmpdbeps2 -C mpbeps2.source

fmpbeps2 :
	make mpbeps2 -C mpbeps2.source

fmpbbeps2 :
	make mpbbeps2 -C mpbeps2.source

fmpdbeps2 :
	make mpdbeps2 -C mpbeps2.source

#noMPI
nompi : fnompi
	make cpmexecs -C mpbeps2.source

mbeps2 : fmbeps2
	make cpmbeps2 -C mpbeps2.source

mbbeps2 : fmbbeps2
	make cpmbbeps2 -C mpbeps2.source

mdbeps2 : fmdbeps2
	make cpmdbeps2 -C mpbeps2.source

fmbeps2 :
	make mbeps2 -C mpbeps2.source

fmbbeps2 :
	make mbbeps2 -C mpbeps2.source

fmdbeps2 :
	make mdbeps2 -C mpbeps2.source

#python :
#	make python -C mpbeps2.source

clean :
	make clean -C mpbeps2.source

clobber: rclobber
	rm -f mpbeps2 mpbbeps2 mpdbeps2 mbeps2 mbbeps2 mdbeps2

rclobber:
	make clobber -C mpbeps2.source
