#-----------------------------------------------------------------------
# This program reads real periodic 2d vector data
# written for 2D MPI/OpenMP PIC codes
# written by Viktor K. Decyk, UCLA
import sys
import math
import numpy

sys.path.append('./ReadFiles')
from libcmfield2 import *
from fgraf import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

ns = 7
iudm = 19; iuv = 12
nplot = 1

dname = numpy.array(["LONGITUDINAL EFIELD ","ELEC CURRENT DENSITY",
                     "VECTOR POTENTIAL    ","TRANSVERSE EFIELD   ",
                     "MAGNETIC FIELD      ","RADIATIVE VPOTENTIAL",
                     "ION CURRENT DENSITY "],'S20')

# create string from idrun
idrun = int(input("enter idrun: "))
cdrun = str(idrun)
fname = "diag2." + cdrun
cmfield2.ffopen2(iudm,fname)

# nscalars = table of available diagnostics
nscalars = numpy.zeros((ns),int_type,'F')

# determine which vector diagnostics are available
cmfield2.readvdiags2(iudm,nscalars)

# select diagnostic
m = numpy.sum(nscalars)
if (m > 1):
   n = -1
   while True:
      if (n < 0):
         for i in xrange(0,ns):
            if (nscalars[i]==1):
               print "enter ", i+1," for ", numpy.str.rstrip(dname[i])
         n = int(input(""))
         if (n==0):
            exit(1)
         if ((n >= 1) and (n <= ns)):
            if (nscalars[n-1]==0):
               n = -1
         else:
            n = -1
         if (n > 0):
            break
         print "invalid entry, try again or enter 0 to quit"
elif (m==1):
   for i in xrange(0,ns):
      if (nscalars[i]==1):
         n = i + 1
         break
else:
   print "no vector diagnostic files found"
   exit(1)

print numpy.str.rstrip(dname[n-1]), " diagnostic selected"

nts = numpy.zeros((1),int_type,'F')
modesx = numpy.zeros((1),int_type,'F')
modesy = numpy.zeros((1),int_type,'F')
mrec = numpy.zeros((1),int_type,'F')
fname = numpy.array([""],'S32')

# return parameters for selected vector diagnostic
cmfield2.vdiagparams2(iudm,n,nts,modesx,modesy,mrec,fname)
nrec = mrec[0]
nplot = in2.ndim

# nx/ny = number of global grid points in x/y direction
nx = int(math.pow(2,in2.indx)); ny = int(math.pow(2,in2.indy))
# kyp = number of real grids in each field partition in y direction
kyp = int((ny - 1)/in2.nvp) + 1
# kyb = minimum number of processors in distributed array
kyb = int((ny - 1)/kyp) + 1
# nyv = third dimension of vector field array, >= ny
nyv = kyp*kyb

# allocate vector array
vfield = numpy.empty((in2.ndim,nx,nyv),float_type,'F')
dt = in2.dt*float(nts[0])

#open stream file for vector field
cmfield2.fsopen2(iuv,fname)

# nrec = number of complete records
nrec = int(nrec/kyb)
print "records found: nrec = ", nrec

ierr = numpy.zeros((1),int_type,'F')
# open graphics device
ierr[0] = graf2.open_graphs2(nplot)
# set palette to color wheel
graf2.set_palit(2)

# read and display data
for ii in xrange(0,nrec):
# read real vector field
   cmfield2.freadv2(iuv,vfield,in2.ndim,nx,nyv)
   it = nts[0]*ii
   time = dt*float(ii)
# show time
   print "it,time=",it,time
# read and display data
   fname = numpy.str.rstrip(dname[n-1])
   graf2.dvector2(vfield,fname,it,999,0,1,1,nx,ny,ierr)
   if (ierr[0]==1):
      exit(1)

cmfield2.closeff2(iudm)
cmfield2.closeff2(iuv)
graf2.close_graphs2()
