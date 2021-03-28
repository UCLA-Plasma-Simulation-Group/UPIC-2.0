#-----------------------------------------------------------------------
# This program reads compressed complex periodic 1d scalar data
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
   reads compressed complex periodic 1d scalar data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 4
   iudm = 19; ius = 11
   iw = 200
   wmin = 0.0; wmax = 2.0
   dname = numpy.array(["POTENTIAL           ","LONGITUDINAL EFIELD ",
                        "ELECTRON DENSITY    ","ION DENSITY         "],
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

# determine which scalar diagnostics are available
   cmfield1.readsdiags1(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   modesx = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   norm = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')
   gname = numpy.array([""],'S32')
   ierr = numpy.zeros((1),int_type,'F')

# open graphics device
   nplot = 1
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
         print ("no complex scalar diagnostic files found")
         n = 0
# exit procedure
      if (n==0):
         if ("sfieldc" in globals()):
            sfieldc = None
         if ("sfield" in globals()):
            sfield = None
         if ("wm" in globals()):
            wm = None
         if ("pkw" in globals()):
            pkw = None
         if ("pks" in globals()):
            pks = None
         if ("wk" in globals()):
            wk = None
         if ("mixup" in globals()):
            mixup = None
         if ("sct" in globals()):
            sct = None
         cmfield1.closeff1(iudm)
         graf1.close_graphs()
         return

      tname = numpy.str.rstrip(dname[n-1])
      print (tname, " diagnostic selected")

# return parameters for selected scalar diagnostic
# nts, modesx, nrec, norm, fname
      cmfield1.sdiagparams1(iudm,n,nts,modesx,mrec,norm,fname)
      nrec = mrec[0]

# nx = number of grid points in x direction
      nx = int(math.pow(2,in1.indx)); nxh = int(nx/2)

# allocate complex scalar arrays
      if ("sfieldc" not in globals()):
         sfieldc = numpy.empty((modesx[0]),complex_type,'F')
# allocate scalar array
      if ("sfield" not in globals()):
         sfield = numpy.empty((nx),float_type,'F')

# open stream file for scalar field
      cmfield1.fsopen1(ius,fname)

# nrec = number of complete records
      print ("records found: nrec = ", nrec)

# allocate and initialize data for frequency analysis
      if ("wm" not in globals()):
         wm = numpy.empty((iw),float_type,'F')
      if ("pkw" not in globals()):
         pkw = numpy.empty((modesx[0],iw,2),float_type,'F')
      if ("pks" not in globals()):
         pks = numpy.zeros((4,modesx[0],iw),double_type,'F')
      if ("wk" not in globals()):
         wk = numpy.empty((modesx[0],2),float_type,'F')
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

# read complex scalar data and display
      for ii in range(0,nrec):
         cmfield1.freadc1(ius,sfieldc,modesx[0])
         it = nts*ii
         time = dt*float(ii)
# perform incremental frequency analysis
         cmfield1.micspect1(sfieldc,wm,pkw,pks,time,0.0,nrec,iw,
                            modesx[0],nx,norm[0])
# decompress field data
         cmfield1.mwrmodes1(sfield,sfieldc,nx,modesx[0])
# fourier transform to real space
         cmfield1.mfft1r(sfield,1,mixup,sct,in1.indx)
# display real space data
         graf1.dscaler1(sfield,numpy.str.rstrip(dname[n-1]),it,999,0,nx,
                     ierr)
         if (ierr[0]==1):
            exit(1)

      graf1.reset_graphs()
      graf1.reset_nplot(2,ierr)

# find the frequency with the maximum power for each mode
#      wk[:,0] = wm[numpy.argmax(pkw[:,:,0],axis=1)]
#      wk[:,1] = wm[numpy.argmax(pkw[:,:,1],axis=1)]
# display frequencies as a function of mode number
#      gname = numpy.str.rstrip(dname[n-1])
#      graf1.dmscaler1(wk,gname,nrec,999,1,modesx[0],cwk,ierr)
#      if (ierr[0]==1):
#         exit(1)

      pmin = numpy.amin(pkw[pkw>0.0])
      pkw[:,:,:] = numpy.where(pkw>0.0,pkw,pmin)
      pkw[:,:,:] = numpy.log(pkw)

# display positive frequencies as a function of mode number
      it = nts*nrec
      gname = numpy.str.rstrip(dname[n-1]+cwk[0])
      graf2.dscalerl2(pkw[:,:,0],gname,akmin,akmax,wmin,wmax,it,999,2,
                      modesx[0],iw,ierr)
# display negative frequencies as a function of mode number
      gname = numpy.str.rstrip(dname[n-1]+cwk[1])
      graf2.dscalerl2(pkw[:,:,1],gname,akmin,akmax,wmin,wmax,it,999,2,
                      modesx[0],iw,ierr)

      graf1.reset_nplot(1,ierr)
      cmfield1.closeff1(ius)
      print()

if (__name__=="__main__"):
   main(-1)
