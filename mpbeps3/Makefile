#Makefile for 3D MPI/OpenMP PIC codes in mpbeps3.source

# Linkage rules

all: fortran

# create all MPI/OpenMP executables
fbuild : fmpbeps3 fmpbbeps3 fmpdbeps3

# copy all MPI/OpenMP executables to mpbeps3 directory
fortran : fbuild
	make cpexecs -C mpbeps3.source

# create individual MPI/OpenMP executables
fmpbeps3 :
	make mpbeps3 -C mpbeps3.source

fmpbbeps3 :
	make mpbbeps3 -C mpbeps3.source

fmpdbeps3 :
	make mpdbeps3 -C mpbeps3.source

# copy individual MPI/OpenMP executables to mpbeps3 directory
mpbeps3 : fmpbeps3
	make cpmpbeps3 -C mpbeps3.source

mpbbeps3 : fmpbbeps3
	make cpmpbbeps3 -C mpbeps3.source

mpdbeps3 : fmpdbeps3
	make cpmpdbeps3 -C mpbeps3.source

#noMPI
# create all OpenMP executables
fnompi : fmbeps3 fmbbeps3 fmdbeps3

# copy all OpenMP executables to mpbeps3 directory
nompi : fnompi
	make cpmexecs -C mpbeps3.source

# create individual OpenMP executables
fmbeps3 :
	make mbeps3 -C mpbeps3.source

fmbbeps3 :
	make mbbeps3 -C mpbeps3.source

fmdbeps3 :
	make mdbeps3 -C mpbeps3.source

# copy individual OpenMP executables to mpbeps3 directory
mbeps3 : fmbeps3
	make cpmbeps3 -C mpbeps3.source

mbbeps3 : fmbbeps3
	make cpmbbeps3 -C mpbeps3.source

mdbeps3 : fmdbeps3
	make cpmdbeps3 -C mpbeps3.source

# create dynamic libraries:
# libmpush3.so, libmbpush3.so, libmdpush3.so, and dtimer.so
python :
	make python -C mpbeps3.source

# remove *.o, *.mod, and *.pyf files
clean :
	make clean -C mpbeps3.source

# remove executables in mpbeps3 directory
clobber: rclobber
	rm -f mpbeps3 mpbbeps3 mpdbeps3 mbeps3 mbbeps3 mdbeps3

# remove executables, dynamic libraries, *.o, *.mod, and *.pyf files
rclobber:
	make clobber -C mpbeps3.source

#Fortran ReadFiles
# create main Fortran data reader: diag3
pfbuild :
	make diag3 -C ReadFiles

# copy main Fortran data reader to mpbeps3 directory
readf :
	make cpdiag3 -C ReadFiles

# create all Fortran data readers:
# diag3 preadf3 pvreadf3 pflreadf3 pfvreadf3 ptreadf3 psreadf3
readfall :
	make -C ReadFiles

#Python ReadFiles
# create dynamic library libcmfield3.so for data readers scripts
rfpython :
	make python -C ReadFiles

# copy all data reader scripts to mpbeps3 directory: preadf3.py 
# pvreadf3.py pflreadf3.py pfvreadf3.py ptreadf3.py and psreadf3.py
rpyreadf : rfpython
	make cpyexecs -C ReadFiles

# remove *.o, *.mod for all data readers
readfclean :
	make clean -C ReadFiles

# remove Fortran data readers from mpbeps3 directory
readfclobber: rfpyclobber
	rm -f diag3 preadf3 pvreadf3 pflreadf3 pfvreadf3 ptreadf3 psreadf3

# remove Fortran executables, dynamic library, *.o, *.mod files
rfpyclobber:
	make clobber -C ReadFiles

list:
	echo fortran mpbeps3 mpbbeps3 mpdbeps3 fnompi mbeps3 mbbeps3 mdbeps3 clean \
	     clobber readf readfall rfpython readfclean readfclobber rfpyclobber
