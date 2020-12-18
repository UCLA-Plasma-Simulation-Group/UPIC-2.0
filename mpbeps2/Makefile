#Makefile for 2D MPI/OpenMP PIC codes in mpbeps2.source

# Linkage rules

all: fortran

# create all MPI/OpenMP executables
fbuild :
	make -C mpbeps2.source

# copy all MPI/OpenMP executables to mpbeps2 directory
fortran : fbuild
	make cpexecs -C mpbeps2.source

# create all MPI/OpenMP executables with X11 graphics
gfbuild :
	make -f gmake -C mpbeps2.source

# copy all MPI/OpenMP executables with X11 graphics to mpbeps2 directory
grfortran : gfbuild
	make cpexecs -C mpbeps2.source

# create individual MPI/OpenMP executables
fmpbeps2 :
	make mpbeps2 -C mpbeps2.source

fmpbbeps2 :
	make mpbbeps2 -C mpbeps2.source

fmpdbeps2 :
	make mpdbeps2 -C mpbeps2.source

# copy individual MPI/OpenMP executables to mpbeps2 directory
mpbeps2 : fmpbeps2
	make cpmpbeps2 -C mpbeps2.source

mpbbeps2 : fmpbbeps2
	make cpmpbbeps2 -C mpbeps2.source

mpdbeps2 : fmpdbeps2
	make cpmpdbeps2 -C mpbeps2.source

#noMPI
# create all OpenMP executables
fnompi :
	make nompi -C mpbeps2.source

# copy all OpenMP executables to mpbeps2 directory
nompi : fnompi
	make cpmexecs -C mpbeps2.source

# create all MPI/OpenMP executables with X11 graphics
gfnompi :
	make -f gmake nompi -C mpbeps2.source

# copy all MPI/OpenMP executables with X11 graphics to mpbeps2 directory
gnompi : gfnompi
	make cpmexecs -C mpbeps2.source

# create individual OpenMP executables
fmbeps2 :
	make mbeps2 -C mpbeps2.source

fmbbeps2 :
	make mbbeps2 -C mpbeps2.source

fmdbeps2 :
	make mdbeps2 -C mpbeps2.source

# copy individual OpenMP executables to mpbeps2 directory
mbeps2 : fmbeps2
	make cpmbeps2 -C mpbeps2.source

mbbeps2 : fmbbeps2
	make cpmbbeps2 -C mpbeps2.source

mdbeps2 : fmdbeps2
	make cpmdbeps2 -C mpbeps2.source

# create dynamic libraries:
# libmpush2.so, libmbpush2.so, libmdpush2.so, fpgraf.so, and dtimer.so
python :
	make python -C mpbeps2.source

# create dynamic libraries:
# libmpush2.so, libmbpush2.so, libmdpush2.so, fpgraf.so, and dtimer.so
# with X11 graphics
gpython :
	make -f gmake python -C mpbeps2.source

# remove *.o, *.mod, and *.pyf files
clean :
	make clean -C mpbeps2.source

# remove executables in mpbeps2 directory
clobber: rclobber
	rm -f mpbeps2 mpbbeps2 mpdbeps2 mbeps2 mbbeps2 mdbeps2

# remove executables, dynamic libraries, *.o, *.mod, and *.pyf files
rclobber:
	make clobber -C mpbeps2.source

#Fortran ReadFiles
# create main Fortran data reader: diag2
pfbuild :
	make diag2 -C ReadFiles

# copy main Fortran data reader to mpbeps2 directory
readf : pfbuild
	make cpdiag2 -C ReadFiles

# create all Fortran data readers:
# diag2 preadf2 pvreadf2 pflreadf2 pfvreadf2 ptreadf2 psreadf2
readfall :
	make -C ReadFiles

# create main Fortran data reader with X11 graphics: diag2
gpfbuild :
	make -f gmake diag2 -C ReadFiles

# copy main Fortran data reader with X11 graphics to mpbeps2 directory
greadf : gpfbuild
	make cpdiag2 -C ReadFiles

#Python ReadFiles
# create dynamic libraries libcmfield2.so and fgraf.so for data readers scripts
# with No graphics
rfpython :
	make python -C ReadFiles

# create dynamic libraries libcmfield2.so and fgraf.so for data readers scripts
# with X11 graphics
grfpython :
	make -f gmake python -C ReadFiles

# copy all data reader scripts to mpbeps2 directory: preadf2.py 
# pvreadf2.py pflreadf2.py pfvreadf2.py ptreadf2.py and psreadf2.py
rpyreadf : rfpython
	make cpyexecs -C ReadFiles

# remove *.o, *.mod for all data readers
readfclean :
	make clean -C ReadFiles

# remove Fortran data readers from mpbeps2 directory
readfclobber: rfpyclobber
	rm -f diag2 preadf2 pvreadf2 pflreadf2 pfvreadf2 ptreadf2 psreadf2

# remove Fortran executables, dynamic library, *.o, *.mod files
rfpyclobber:
	make clobber -C ReadFiles

list:
	echo fortran mpbeps2 mpbbeps2 mpdbeps2 fnompi mbeps2 mbbeps2 mdbeps2 clean \
	     clobber readf readfall rfpython readfclean readfclobber rfpyclobber \
	     grfortran gnompi gpython greadf grphython
