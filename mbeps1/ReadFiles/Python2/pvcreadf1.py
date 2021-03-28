#-----------------------------------------------------------------------
# This program reads compressed complex periodic 1d vector data
# written for 1D OpenMP PIC codes
# written by Viktor K. Decyk, UCLA
import sys
import math
import numpy

sys.path.append('./ReadFiles')
from libcmfield1 import *
from fgraf import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

ns = 6
iudm = 19; iuv = 12
nplot = 1; iw = 800
wmin = 0.0; wmax = 8.0
dname = numpy.array(["LONGITUDINAL EFIELD ","ELEC CURRENT DENSITY",
                     "VECTOR POTENTIAL    ","TRANSVERSE EFIELD   ",
                     "MAGNETIC FIELD      ","RADIATIVE VPOTENTIAL",
                     "ION CURRENT DENSITY "],'S20')
cwk = numpy.array([" W > 0"," W < 0"],'S6')

# create string from idrun
idrun = int(input("enter idrun: "))
cdrun = str(idrun)
fname = "diag1." + cdrun
cmfield1.ffopen1(iudm,fname)

# nscalars = table of available diagnostics
nscalars = numpy.zeros((ns),int_type,'F')

# determine which vector diagnostics are available
cmfield1.readvdiags1(iudm,nscalars)

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
   print "no scalar diagnostic files found"
   exit(1)

print numpy.str.rstrip(dname[n-1]), " diagnostic selected"

nts = numpy.zeros((1),int_type,'F')
modesx = numpy.zeros((1),int_type,'F')
mrec = numpy.zeros((1),int_type,'F')
norm = numpy.zeros((1),int_type,'F')
fname = numpy.array([""],'S32')
# return parameters for selected vector diagnostic
cmfield1.vdiagparams1(iudm,n,nts,modesx,mrec,norm,fname)
nrec = mrec[0]
mdim = in1.ndim - 1
nplot = mdim

# allocate complex scalar array
nx = int(math.pow(2,in1.indx)); nxh = int(nx/2)

# allocate complex vector array
vfieldc = numpy.empty((mdim,modesx[0]),complex_type,'F')
# allocate vector array
vfield = numpy.empty((mdim,nx),float_type,'F')

# open stream file for vector field
cmfield1.fsopen1(iuv,fname)

# nrec = number of complete records
print "records found: nrec = ", nrec

# allocate and initialize data for frequency analysis
wm = numpy.empty((iw),float_type,'F')
vpkw = numpy.empty((mdim,modesx[0],iw,2),float_type,'F')
vpks = numpy.zeros((mdim,4,modesx[0],iw),double_type,'F')
vwk = numpy.empty((mdim,modesx[0],2),float_type,'F')
wm[:] = ((wmax-wmin)/float(iw))*numpy.linspace(0,iw-1,iw)
dt = in1.dt*float(nts[0])
dnx = 6.28318530717959/float(nx)
akmin = 0.0; akmax = dnx*float(modesx[0] - 1)

ierr = numpy.zeros((1),int_type,'F')
# open graphics device
ierr[0] = graf1.open_graphs(nplot)
# set palette to color wheel
graf2.set_palit(2)

# prepare fft tables for decompression
mixup = numpy.empty((nxh),int_type,'F')
sct = numpy.empty((nxh),complex_type,'F')
cmfield1.mfft1_init(mixup,sct,in1.indx)

# read complex vector data and display
for ii in xrange(0,nrec):
   cmfield1.freadvc1(iuv,vfieldc,mdim,modesx[0])
   time = dt*float(ii)
# perform incremental frequency analysis
   cmfield1.mivcspect1(vfieldc,wm,vpkw,vpks,time,0.0,nrec,iw,modesx[0],
                       nx,norm[0])
# decompress field data
   cmfield1.mwrvmodes1(vfield,vfieldc,nx,modesx[0])
# fourier transform to real space
   cmfield1.mfft1rn(vfield,1,mixup,sct,in1.indx)
# display real space data
   graf1.dvector1(vfield,numpy.str.rstrip(dname[n-1]),ii,999,0,1,nx,ierr)
   if (ierr[0]==1):
      exit(1)

graf1.reset_graphs()
graf1.reset_nplot(1,ierr)
graf1.reset_nplot(4,ierr)

# find the frequency with the maximum power for each mode
#vwk[0,:,0] = wm[numpy.argmax(vpkw[0,:,:,0],axis=1)]
#vwk[1,:,0] = wm[numpy.argmax(vpkw[1,:,:,0],axis=1)]
#vwk[0,:,1] = wm[numpy.argmax(vpkw[0,:,:,1],axis=1)]
#vwk[1,:,1] = wm[numpy.argmax(vpkw[1,:,:,1],axis=1)]
# display frequencies as a function of mode number
#fname = numpy.str.rstrip(dname[n-1])
#graf1.dmvector1(vwk,fname,nrec,999,2,2,modesx[0],cwk,ierr)
#if (ierr[0]==1):
#   exit(1)

pmin = numpy.amin(vpkw[vpkw>0.0])
vpkw[:,:,:,:] = numpy.where(vpkw>0.0,vpkw,pmin)
vpkw[:,:,:,:] = numpy.log(vpkw)

# display positive frequencies as a function of mode number
fname = numpy.str.rstrip(dname[n-1]+cwk[0])
# display positive frequencies as a function of mode number
graf2.dvectorl2(vpkw[:,:,:,0],fname,akmin,akmax,wmin,wmax,nrec,999,2,1,
                modesx[0],iw,ierr)
# display negative frequencies as a function of mode number
fname = numpy.str.rstrip(dname[n-1]+cwk[1])
graf2.dvectorl2(vpkw[:,:,:,1],fname,akmin,akmax,wmin,wmax,nrec,999,2,1,
                modesx[0],iw,ierr)

cmfield1.closeff(iudm)
cmfield1.closeff(iuv)
graf1.close_graphs()


