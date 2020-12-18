#-----------------------------------------------------------------------
# This program reads real 3d phase space data
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
nplot = 1
dname = numpy.array(["electron phase space","ion phase space     "],
                     'S20')
sname = numpy.array(["ELECTRON","ION     "],'S8')

# create string from idrun
idrun = int(input("enter idrun: "))
cdrun = str(idrun)
fname = "diag3." + cdrun
cmfield3.ffopen3(iudm,fname)

# nscalars = table of available diagnostics
nscalars = numpy.zeros((ns),int_type,'F')

# determine which phase space diagnostics are available
cmfield3.readpsdiags3(iudm,nscalars)

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
   print "no phase space diagnostic files found"
   exit(1)

print numpy.str.rstrip(dname[n-1]), " diagnostic selected"

nts = numpy.zeros((1),int_type,'F')
nmv = numpy.zeros((1),int_type,'F')
nsxb = numpy.zeros((1),int_type,'F')
nsyb = numpy.zeros((1),int_type,'F')
nszb = numpy.zeros((1),int_type,'F')
mrec = numpy.zeros((1),int_type,'F')
fname = numpy.array([""],'S32')

# return parameters for selected phase space diagnostic
cmfield3.psdiagparams3(iudm,n,nts,nmv,nsxb,nsyb,nszb,mrec,fname)
nrec = mrec[0]

# nx/ny/nz = number of global grid points in x/y/z direction
nx = int(math.pow(2,in3.indx)); ny = int(math.pow(2,in3.indy))
nz = int(math.pow(2,in3.indz))
# nmv21 = number of velocity bins
nmv21 = 2*nmv[0] + 1; nmvf = nmv21 + 1

# allocate velocity distribution arrays
noff = numpy.zeros((2,in3.nvpy,in3.nvpz),int_type,'F')
nyzp = numpy.zeros((2,in3.nvpy,in3.nvpz),int_type,'F')

fvs = numpy.empty((nmvf,in3.ndim,nsxb[0],nsyb[0],nszb[0]),float_type,
                  'F')
pvs = numpy.empty((nmvf,in3.ndim,max(nsxb[0],nsyb[0],nszb[0])),
                  float_type,'F')
ps = numpy.empty((nmvf,max(nsxb[0],nsyb[0],nszb[0])),float_type,'F')
psx = numpy.empty((nsxb[0],nmvf),float_type,'F')
psy = numpy.empty((nsyb[0],nmvf),float_type,'F')
psz = numpy.empty((nszb[0],nmvf),float_type,'F')
dt = in3.dt*float(nts[0])

# open stream file for phase space field
cmfield3.fsopen3(iuv,fname)

ierr = numpy.zeros((1),int_type,'F')
# reads distributed non-uniform partition information
cmfield3.freadncomp3(iuv,noff,nyzp,ierr)
print "ierr=",ierr[0]
if (ierr[0] != 0):
   exit(1)

# nrec = number of complete records
print "records found: nrec = ", nrec

# read and display data
for ii in xrange(0,nrec):
# read real velocity distribution fields
   cmfield3.freadps3(iuv,fvs,noff,nyzp,in3.ndim,nmvf,nsxb[0],in3.nvpy,
                     in3.nvpz)
   it = nts[0]*ii
   time = dt*float(ii)
   print sname[n-1], " it,time=",it,time
# sum over y and z
   cmfield3.summit2(fvs,pvs,4,5)
# select x-vx
   ps[:,:] = pvs[:,0,:nsxb[0]]
   psx[:,:] = numpy.transpose(ps)
# select x-vy
   ps[:,:] = pvs[:,1,:nsxb[0]]
   psx[:,:] = numpy.transpose(ps)
# select x-vz
   ps[:,:] = pvs[:,2,:nsxb[0]]
   psx[:,:] = numpy.transpose(ps)
# sum over x and z
   cmfield3.summit2(fvs,pvs,3,5)
# select y-vx
   ps[:,:] = pvs[:,0,:nsyb[0]]
   psy[:,:] = numpy.transpose(ps)
# select y-vy
   ps[:,:] = pvs[:,1,:nsyb[0]]
   psy[:,:] = numpy.transpose(ps)
# select y-vz
   ps[:,:] = pvs[:,2,:nsyb[0]]
   psy[:,:] = numpy.transpose(ps)
# sum over x and y
   cmfield3.summit2(fvs,pvs,3,4)
# select z-vx
   ps[:,:] = pvs[:,0,:nszb[0]]
   psz[:,:] = numpy.transpose(ps)
# select z-vy
   ps[:,:] = pvs[:,1,:nszb[0]]
   psz[:,:] = numpy.transpose(ps)
# select z-vz
   ps[:,:] = pvs[:,2,:nszb[0]]
   psz[:,:] = numpy.transpose(ps)

cmfield3.closeff3(iudm)
cmfield3.closeff3(iuv)
