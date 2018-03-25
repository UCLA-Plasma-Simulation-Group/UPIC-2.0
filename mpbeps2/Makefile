#Makefile for 2D MPI/OpenMP PIC codes in mpbeps2.source

# Linkage rules

all: fortran

fortran : fbuild
	make cpexecs -C mpbeps2.source

fbuild :
	make -C mpbeps2.source

mpbeps2 : fmpbeps2
	make cpmpbeps2 -C mpbeps2.source

fmpbeps2 :
	make mpbeps2 -C mpbeps2.source

mpbbeps2 : fmpbbeps2
	make cpmpbbeps2 -C mpbeps2.source

fmpbbeps2 :
	make mpbbeps2 -C mpbeps2.source

mpdbeps2 : fmpdbeps2
	make cpmpdbeps2 -C mpbeps2.source

fmpdbeps2 :
	make mpdbeps2 -C mpbeps2.source

#noMPI
nompi : fnompi
	make cpmexecs -C mpbeps2.source

fnompi :
	make nompi -C mpbeps2.source

mbeps2 : fmbeps2
	make cpmbeps2 -C mpbeps2.source

fmbeps2 :
	make mbeps2 -C mpbeps2.source

mbbeps2 : fmbbeps2
	make cpmbbeps2 -C mpbeps2.source

fmbbeps2 :
	make mbbeps2 -C mpbeps2.source

mdbeps2 : fmdbeps2
	make cpmdbeps2 -C mpbeps2.source

fmdbeps2 :
	make mdbeps2 -C mpbeps2.source

clean :
	make clean -C mpbeps2.source

clobber: rclobber
	rm -f mpbeps2 mpbbeps2 mpdbeps2 mbeps2 mbbeps2 mdbeps2

rclobber:
	make clobber -C mpbeps2.source

#Fortran ReadFiles
readf : pfbuild
	make cpexecs -C ReadFiles

pfbuild :
	make -C ReadFiles

#Python ReadFiles
pyreadf : pybuild
	make cpyexecs -C ReadFiles

pybuild :
	make python -C ReadFiles

readfclean :
	make clean -C ReadFiles

readfclobber: readfrclobber
	rm -f preadf2 pvreadf2 pflreadf2 pfvreadf2 ptreadf2 psreadf2 \
	      preadf2.py pvreadf2.py pflreadf2.py pfvreadf2.py ptreadf2.py \
	      psreadf2.py

readfrclobber:
	make clobber -C ReadFiles

list:
	echo fortran mpbeps2 mpbbeps2 mpdbeps2 fnompi mbeps2 mbbeps2 mdbeps2 \
	     readf clean clobber pybuild pyreadf readfclean readfclobber