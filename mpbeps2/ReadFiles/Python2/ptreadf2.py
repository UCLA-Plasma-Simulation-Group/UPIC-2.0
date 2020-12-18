#-----------------------------------------------------------------------
# This program reads real 2d trajectory data
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

ns = 1
iudm = 19; iuv = 12
nplot = 1
sname = numpy.array(["ELECTRON","ION     "],'S8')

# create string from idrun
idrun = int(input("enter idrun: "))
cdrun = str(idrun)
fname = "diag2." + cdrun
cmfield2.ffopen2(iudm,fname)

# nscalars = table of available diagnostics
nscalars = numpy.zeros((ns),int_type,'F')

# determine which trajectory diagnostics are available
cmfield2.readtrdiags2(iudm,nscalars)

# select diagnostic
m = numpy.sum(nscalars)
if (m==1):
   for i in xrange(0,ns):
      if (nscalars[i]==1):
         n = i + 1
         break
else:
   print "no trajectory diagnostic files found"
   exit(1)

print " trajectory diagnostic selected"

nts = numpy.zeros((1),int_type,'F')
ndt = numpy.zeros((1),int_type,'F')
nst = numpy.zeros((1),int_type,'F')
nmv = numpy.zeros((1),int_type,'F')
ndimp = numpy.zeros((1),int_type,'F')
nprobt = numpy.zeros((1),int_type,'F')
mrec = numpy.zeros((1),int_type,'F')
fname = numpy.array([""],'S32')

# return parameters for selected trajectory diagnostic
cmfield2.trdiagparams2(iudm,n,nts,ndt,nst,nmv,ndimp,nprobt,mrec,fname)
nrec = mrec[0]
nplot = 1

print numpy.str.rstrip(sname[ndt[0]-1]), " trajectory available"

# nmvf = size of velocity distribution array
nmvf = 2*nmv[0] + 2

# allocate trajectory arrays
if ((nst[0]==1) or (nst[0]==2)):
   partt = numpy.empty((ndimp[0],nprobt[0]),float_type,'F')
   partd = numpy.empty((nrec,ndimp[0],nprobt[0]),float_type,'F')
# allocate trajectory distribution arrays
elif (nst[0]==3):
   fvmtp = numpy.empty((in2.ndim,3),float_type,'F')
   fvtp = numpy.empty((nmvf,in2.ndim),float_type,'F')
   fetp = numpy.empty((nmvf,0),float_type,'F')
   wk = numpy.empty((1),float_type,'F')
dt = in2.dt*float(nts[0])

# open stream file for trajectory field
cmfield2.fsopen2(iuv,fname)

# nrec = number of complete records
print "records found: nrec = ", nrec

ierr = numpy.zeros((1),int_type,'F')
# open graphics device
ierr[0] = graf2.open_graphs2(nplot)

# read and display data
for ii in xrange(0,nrec):
   it = nts[0]*ii
   time = dt*float(ii)
# read real trajectory fields
   if ((nst[0]==1) or (nst[0]==2)):
      cmfield2.freadtr2(iuv,partt,ndimp,nprobt)
      partd[ii,:,:] = partt
      print sname[ndt[0]-1], " trajectory:it,time=",it,time
# display cartesian distribution
   elif (nst[0]==3):
      cmfield2.freadfv2(iuv,fvmtp,fvtp,fetp,wk,in2.ndim,nmvf,in2.ndim,0)
      print sname[ndt[0]-1], " distribution:it,time=",it,time
      fname = numpy.str.rstrip(sname[ndt[0]-1])
      graf1.displayfv1(fvtp,fvmtp,fname,it,nmv,2,ierr)
      if (ierr[0]==1):
         exit(1)

# display time history of trajectories
if ((nst[0]==1) or (nst[0]==2)):
# display vx
   ix = 3
   graf1.displaytr1(partd,0.0,dt,nrec,ix,999,ierr)
   if (ierr[0]==1):
      exit(1)

cmfield2.closeff2(iudm)
cmfield2.closeff2(iuv)
graf2.close_graphs2()

