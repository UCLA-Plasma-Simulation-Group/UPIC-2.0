#-----------------------------------------------------------------------
# This program reads real periodic 3d velocity data
# written for 3D MPI/OpenMP PIC codes
# written by Viktor K. Decyk, UCLA
import sys
import math
import numpy

sys.path.append('./ReadFiles')
from libcmfield3 import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

ns = 2
iudm = 19; iuv = 12
dname = numpy.array(["electron velocity   ","ion velocity        "],
                     'S20')
sname = numpy.array(["ELECTRON","ION     "],'S8')

# create string from idrun
idrun = int(input("enter idrun: "))
cdrun = str(idrun)
fname = "diag3." + cdrun
cmfield3.ffopen3(iudm,fname)

# nscalars = table of available diagnostics
nscalars = numpy.zeros((ns),int_type,'F')

# determine which velocity diagnostics are available
cmfield3.readfvdiags3(iudm,nscalars)

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
   print "no velocity-space diagnostic files found"
   exit(1)

print numpy.str.rstrip(dname[n-1]), " diagnostic selected"

nts = numpy.zeros((1),int_type,'F')
nmv = numpy.zeros((1),int_type,'F')
nfvd = numpy.zeros((1),int_type,'F')
nfed = numpy.zeros((1),int_type,'F')
mrec = numpy.zeros((1),int_type,'F')
fname = numpy.array([""],'S32')

# return parameters for selected velocity diagnostic
cmfield3.fvdiagparams3(iudm,n,nts,nmv,nfvd,nfed,mrec,fname)
nrec = mrec[0]
nplot = 4

# nmv21 = number of velocity bins
nmv21 = 2*nmv[0] + 1; nmvf = nmv21 + 1

# allocate velocity distribution arrays
fvm = numpy.empty((in3.ndim,3),float_type,'F')
fv = numpy.empty((nmvf,nfvd[0]),float_type,'F')
fe = numpy.empty((nmvf,nfed[0]),float_type,'F')
wk = numpy.empty((1),float_type,'F')
dt = in3.dt*float(nts[0])

# open stream file for velocity field
cmfield3.fsopen3(iuv,fname)

# nrec = number of complete records
print "records found: nrec = ", nrec

# read and display data
for ii in xrange(0,nrec):
# read real velocity fields
   cmfield3.freadfv3(iuv,fvm,fv,fe,wk,in3.ndim,nmvf,nfvd,nfed)
   it = nts[0]*ii
   time = dt*float(ii)
   if (nfvd[0] > 0):
# display cylindrical distribution
      if (nfvd[0] < in3.ndim):
         print sname[n-1], " cylindrical:it,time=",it,time
# display cartesian distribution
      else:
         print sname[n-1], " cartesian:it,time=",it,time
# display energy distribution
   if (nfed[0] > 0):
      print sname[n-1], " energy:it,time=",it,time

cmfield3.closeff3(iudm)
cmfield3.closeff3(iuv)
