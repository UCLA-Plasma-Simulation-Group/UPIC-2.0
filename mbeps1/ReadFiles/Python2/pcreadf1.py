#-----------------------------------------------------------------------
# This program reads compressed complex periodic 1d scalar data
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

ns = 3
iudm = 19; ius = 11
nplot = 1; iw = 200
wmin = 0.0; wmax = 2.0
dname = numpy.array(["POTENTIAL       ","ELECTRON DENSITY",
                     "ION DENSITY     "],'S16')
cwk = numpy.array([" W > 0"," W < 0"],'S6')

# create string from idrun
idrun = int(input("enter idrun: "))
cdrun = str(idrun)
fname = "diag1." + cdrun
cmfield1.ffopen1(iudm,fname)

# nscalars = table of available diagnostics
nscalars = numpy.zeros((ns),int_type,'F')

# determine which scalar diagnostics are available
cmfield1.readsdiags1(iudm,nscalars)

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
# return parameters for selected scalar diagnostic
cmfield1.sdiagparams1(iudm,n,nts,modesx,mrec,norm,fname)
nrec = mrec[0]

# nx = number of grid points in x direction
nx = int(math.pow(2,in1.indx)); nxh = int(nx/2)

# allocate complex scalar array
sfieldc = numpy.empty((modesx[0]),complex_type,'F')
# allocate scalar array
sfield = numpy.empty((nx),float_type,'F')

# open stream file for scalar field
cmfield1.fsopen1(ius,fname)

# nrec = number of complete records
print "records found: nrec = ", nrec

# allocate and initialize data for frequency analysis
wm = numpy.empty((iw),float_type,'F')
pkw = numpy.empty((modesx[0],iw,2),float_type,'F')
pks = numpy.zeros((4,modesx[0],iw),double_type,'F')
wk = numpy.empty((modesx[0],2),float_type,'F')
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

# read complex scalar data and display
for ii in xrange(0,nrec):
   cmfield1.freadc1(ius,sfieldc,modesx[0])
   time = dt*float(ii)
# perform incremental frequency analysis
   cmfield1.micspect1(sfieldc,wm,pkw,pks,time,0.0,nrec,iw,modesx[0],nx,
                      norm[0])
# decompress field data
   cmfield1.mwrmodes1(sfield,sfieldc,nx,modesx[0])
# fourier transform to real space
   cmfield1.mfft1r(sfield,1,mixup,sct,in1.indx)
# display real space data
   graf1.dscaler1(sfield,numpy.str.rstrip(dname[n-1]),ii,999,0,nx,ierr)
   if (ierr[0]==1):
      exit(1)

graf1.reset_graphs()
graf1.reset_nplot(2,ierr)

# find the frequency with the maximum power for each mode
#wk[:,0] = wm[numpy.argmax(pkw[:,:,0],axis=1)]
#wk[:,1] = wm[numpy.argmax(pkw[:,:,1],axis=1)]
# display frequencies as a function of mode number
#fname = numpy.str.rstrip(dname[n-1])
#graf1.dmscaler1(wk,fname,nrec,999,1,modesx[0],cwk,ierr)
#if (ierr[0]==1):
#   exit(1)

pmin = numpy.amin(pkw[pkw>0.0])
pkw[:,:,:] = numpy.where(pkw>0.0,pkw,pmin)
pkw[:,:,:] = numpy.log(pkw)

# display positive frequencies as a function of mode number
fname = numpy.str.rstrip(dname[n-1]+cwk[0])
graf2.dscalerl2(pkw[:,:,0],fname,akmin,akmax,wmin,wmax,nrec,999,2,
                modesx[0],iw,ierr)
# display negative frequencies as a function of mode number
fname = numpy.str.rstrip(dname[n-1]+cwk[1])
graf2.dscalerl2(pkw[:,:,1],fname,akmin,akmax,wmin,wmax,nrec,999,2,
                modesx[0],iw,ierr)

cmfield1.closeff(iudm)
cmfield1.closeff(ius)
graf1.close_graphs()
