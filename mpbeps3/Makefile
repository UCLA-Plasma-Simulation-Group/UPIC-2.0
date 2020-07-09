#Makefile for 3D MPI/OpenMP PIC codes in mpbeps3.source

# Linkage rules

all: fortran

fbuild : fmpbeps3 fmpbbeps3 fmpdbeps3

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

fnompi : fmbeps3 fmbbeps3 fmdbeps3

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

python :
	make python -C mpbeps3.source

clean :
	make clean -C mpbeps3.source

clobber: rclobber
	rm -f mpbeps3 mpbbeps3 mpdbeps3 mbeps3 mbbeps3 mdbeps3

rclobber:
	make clobber -C mpbeps3.source

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
	rm -f preadf3 pvreadf3 pflreadf3 pfvreadf3 ptreadf3 psreadf3 \
	      preadf3.py pvreadf3.py pflreadf3.py pfvreadf3.py ptreadf3.py \
	      psreadf3.py

readfrclobber:
	make clobber -C ReadFiles

list:
	echo fortran mpbeps3 mpbbeps3 mpdbeps3 fnompi mbeps3 mbbeps3 mdbeps3 \
	     readf clean clobber pybuild pyreadf readfclean readfclobber

