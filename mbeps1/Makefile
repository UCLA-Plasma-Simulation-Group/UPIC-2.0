#Makefile for 1D OpenMP PIC codes in mbeps1.source

# Linkage rules

all : fortran

# create all OpenMP executables
fbuild :
	make -C mbeps1.source

# copy all OpenMP executables to mbeps1 directory
fortran : fbuild
	make cpexecs -C mbeps1.source

# create all MPI/OpenMP executables with X11 graphics
gfbuild :
	make -f gmake -C mbeps1.source

# copy all OpenMP executables with X11 graphics to mbeps1 directory
grfortran : gfbuild
	make cpexecs -C mbeps1.source

# create individual OpenMP executables
fmbeps1 :
	make mbeps1 -C mbeps1.source

fmbbeps1 :
	make mbbeps1 -C mbeps1.source

fmdbeps1 :
	make mdbeps1 -C mbeps1.source

# copy individual OpenMP executables to mbeps1 directory
mbeps1 : fmbeps1
	make cpmbeps1 -C mbeps1.source

mbbeps1 : fmbbeps1
	make cpmbbeps1 -C mbeps1.source

mdbeps1 : fmdbeps1
	make cpmdbeps1 -C mbeps1.source

# create dynamic libraries:
# libmpush1.so, libmbpush1.so, libmdpush1.so, fgraf.so, and dtimer.so
python :
	make python -C mbeps1.source

# create dynamic libraries:
# libmpush1.so, libmbpush1.so, libmdpush1.so, fgraf.so, and dtimer.so
# with X11 graphics
gpython :
	make -f gmake python -C mbeps1.source

# remove *.o, *.mod, and *.pyf files
clean :
	make clean -C mbeps1.source

# remove executables in mbeps1 directory
clobber: rclobber
	rm -f mbeps1 mbbeps1 mdbeps1

# remove executables, dynamic libraries, *.o, *.mod, and *.pyf files
rclobber:
	make clobber -C mbeps1.source

#Fortran ReadFiles
# create main Fortran data reader: diag1
pfbuild :
	make diag1 -C ReadFiles

# copy main Fortran data reader to mbeps1 directory
readf : pfbuild
	make cpdiag1 -C ReadFiles

# create all Fortran data readers:
# diag1 pcreadf1 pvcreadf1 pflreadf1 pfvreadf1 ptreadf1 psreadf1
readfall :
	make -C ReadFiles

# create main Fortran data reader with X11 graphics: diag1
gpfbuild :
	make -f gmake diag1 -C ReadFiles

# copy main Fortran data reader with X11 graphics to mpbeps1 directory
greadf : gpfbuild
	make cpdiag1 -C ReadFiles

#Python ReadFiles
# create dynamic libraries libcmfield1.so and fgraf.so for data readers scripts
# with No graphics
rfpython :
	make python -C ReadFiles

# create dynamic libraries libcmfield1.so and fgraf.so for data readers scripts
# with X11 graphics
grfpython :
	make -f gmake python -C ReadFiles

# copy all data reader scripts to mbeps1 directory: preadf1.py 
# pvreadf1.py pflreadf1.py pfvreadf1.py ptreadf1.py and psreadf1.py
rpyreadf : rfpython
	make cpyexecs -C ReadFiles

# remove *.o, *.mod for all data readers
readfclean :
	make clean -C ReadFiles

# remove Fortran data readers from mbeps1 directory
readfclobber: rfpyclobber
	rm -f pcreadf1 pvcreadf1 pflreadf1 pcreadf1.py pvcreadf1.py pflreadf1.py

# remove Fortran executables, dynamic library, *.o, *.mod files
rfpyclobber:
	make clobber -C ReadFiles

list:
	echo fortran mbeps1 mbbeps1 mdbeps1 clean clobber readf readfall \
	     rfpython readfclean readfclobber rfpyclobber \
	     grfortran gpython greadf grfpython
