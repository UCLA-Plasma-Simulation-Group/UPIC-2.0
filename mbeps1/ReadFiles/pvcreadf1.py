#-----------------------------------------------------------------------
# This program reads compressed complex periodic 1d vector data
# written for 1D OpenMP PIC codes
# written by Viktor K. Decyk, UCLA
import sys
import math
import numpy

sys.path.append('./ReadFiles')
from fgraf import *
from libcmfield1 import *

if (sys.version_info.major==2):
   input = raw_input

#-----------------------------------------------------------------------
def main(idrun):
   """
   reads compressed complex periodic 1d vector data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 6
   iudm = 19; iuv = 12
   nplot = 1; iw = 800
   wmin = 0.0; wmax = 8.0
   dname = numpy.array(["ELEC CURRENT DENSITY","VECTOR POTENTIAL    ",
                        "TRANSVERSE EFIELD   ","MAGNETIC FIELD      ",
                        "RADIATIVE VPOTENTIAL","ION CURRENT DENSITY "],
                        dtype=str)
   cwk = numpy.array([" W > 0"," W < 0"],dtype=str)

# create string from idrun
   if (idrun < 0):
      cdrun = "Unknown"
      while (cdrun.isdigit() == False):
         cdrun = input("enter integer idrun: ")
      idrun = int(cdrun)
   cdrun = str(idrun)
   fname = "diag1." + cdrun
   cmfield1.ffopen1(iudm,fname)

# nscalars = table of available diagnostics
   nscalars = numpy.zeros((ns),int_type,'F')

# determine which vector diagnostics are available
   cmfield1.readvdiags1(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   modesx = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   norm = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')
   gname = numpy.array([""],'S32')
   ierr = numpy.zeros((1),int_type,'F')

# open graphics device
   mdim = in1.ndim - 1
   nplot = mdim
   ierr[0] = graf1.open_graphs(nplot)
# set palette to color wheel
   graf2.set_palit(2)

# select diagnostic
   m = numpy.sum(nscalars)
   while True:
      if (m > 0):
         n = -1
         while True:
            if (n < 0):
               for i in range(0,ns):
                  if (nscalars[i]==1):
                     print ("enter ", i+1," for", 
                            numpy.str.rstrip(dname[i]))
               print ("enter ", 0," for EXIT")
               c = input("")
               if (c.isdigit()):
                  n = int(c)
               if (n==0):
                  break
               if ((n >= 1) and (n <= ns)):
                  if (nscalars[n-1]==0):
                     n = -1
            else:
               n = -1
            if (n > 0):
               break
            print ("invalid entry, try again or enter 0 to quit")
      else:
         print ("no complex vector diagnostic files found")
         n = 0
# exit procedure
      if (n==0):
         if ("vfieldc" in globals()):
            vfieldc = None
         if ("vfield" in globals()):
            vfield = None
         if ("wm" in globals()):
            wm = None
         if ("vpkw" in globals()):
            vpkw = None
         if ("vpks" in globals()):
            vpks = None
         if ("vwk" in globals()):
            vwk = None
         if ("mixup" in globals()):
            mixup = None
         if ("sct" in globals()):
            sct = None
         cmfield1.closeff1(iudm)
         graf1.close_graphs()
         return

      tname = numpy.str.rstrip(dname[n-1])
      print (tname, " diagnostic selected")

# return parameters for selected vector diagnostic:
# nts, modesx, nrec, norm, fname
      cmfield1.vdiagparams1(iudm,n,nts,modesx,mrec,norm,fname)
      nrec = mrec[0]

# allocate complex scalar array
      nx = int(math.pow(2,in1.indx)); nxh = int(nx/2)

# allocate complex vector array
      if ("vfieldc" not in globals()):
         vfieldc = numpy.empty((mdim,modesx[0]),complex_type,'F')
# allocate vector array
      if ("vfield" not in globals()):
         vfield = numpy.empty((mdim,nx),float_type,'F')

# open stream file for vector field
      cmfield1.fsopen1(iuv,fname)

# nrec = number of complete records
      print ("records found: nrec = ", nrec)

# allocate and initialize data for frequency analysis
      if ("wm" not in globals()):
         wm = numpy.empty((iw),float_type,'F')
      if ("vpkw" not in globals()):
         vpkw = numpy.empty((mdim,modesx[0],iw,2),float_type,'F')
      if ("vpks" not in globals()):
         vpks = numpy.zeros((mdim,4,modesx[0],iw),double_type,'F')
      if ("vwk" not in globals()):
         vwk = numpy.empty((mdim,modesx[0],2),float_type,'F')
      if ("wm" not in globals()):
         wm[:] = ((wmax-wmin)/float(iw))*numpy.linspace(0,iw-1,iw)
      dt = in1.dt*float(nts[0])
      dnx = 6.28318530717959/float(nx)
      akmin = 0.0; akmax = dnx*float(modesx[0] - 1)

# prepare fft tables for decompression
      if ("mixup" not in globals()):
         mixup = numpy.empty((nxh),int_type,'F')
      if ("sct" not in globals()):
         sct = numpy.empty((nxh),complex_type,'F')
      cmfield1.mfft1_init(mixup,sct,in1.indx)

# read complex vector data and display
      for ii in range(0,nrec):
         cmfield1.freadvc1(iuv,vfieldc,mdim,modesx[0])
         it = nts*ii
         time = dt*float(ii)
# perform incremental frequency analysis
         cmfield1.mivcspect1(vfieldc,wm,vpkw,vpks,time,0.0,nrec,iw,
                             modesx[0],nx,norm[0])
# decompress field data
         cmfield1.mwrvmodes1(vfield,vfieldc,nx,modesx[0])
# fourier transform to real space
         cmfield1.mfft1rn(vfield,1,mixup,sct,in1.indx)
# display real space data
         graf1.dvector1(vfield,numpy.str.rstrip(dname[n-1]),it,999,0,1,
                        nx,ierr)
         if (ierr[0]==1):
            exit(1)

      graf1.reset_graphs()
      graf1.reset_nplot(1,ierr)
      graf1.reset_nplot(nplot,ierr)

# find the frequency with the maximum power for each mode
#      vwk[0,:,0] = wm[numpy.argmax(vpkw[0,:,:,0],axis=1)]
#      vwk[1,:,0] = wm[numpy.argmax(vpkw[1,:,:,0],axis=1)]
#      vwk[0,:,1] = wm[numpy.argmax(vpkw[0,:,:,1],axis=1)]
#      vwk[1,:,1] = wm[numpy.argmax(vpkw[1,:,:,1],axis=1)]
# display frequencies as a function of mode number
#      gname = numpy.str.rstrip(dname[n-1])
#      graf1.dmvector1(vwk,gname,nrec,999,2,2,modesx[0],cwk,ierr)
#      if (ierr[0]==1):
#         exit(1)

      pmin = numpy.amin(vpkw[vpkw>0.0])
      vpkw[:,:,:,:] = numpy.where(vpkw>0.0,vpkw,pmin)
      vpkw[:,:,:,:] = numpy.log(vpkw)

# display positive frequencies as a function of mode number
      it = nts*nrec
      gname = numpy.str.rstrip(dname[n-1]+cwk[0])
# display positive frequencies as a function of mode number
      graf2.dvectorl2(vpkw[:,:,:,0],gname,akmin,akmax,wmin,wmax,it,999,
                      2,1,modesx[0],iw,ierr)
# display negative frequencies as a function of mode number
      gname = numpy.str.rstrip(dname[n-1]+cwk[1])
      graf2.dvectorl2(vpkw[:,:,:,1],gname,akmin,akmax,wmin,wmax,it,999,
                      2,1,modesx[0],iw,ierr)

      graf1.reset_nplot(1,ierr)
      graf1.reset_nplot(nplot,ierr)
      cmfield1.closeff1(iuv)
      print()

if (__name__=="__main__"):
   main(-1)
