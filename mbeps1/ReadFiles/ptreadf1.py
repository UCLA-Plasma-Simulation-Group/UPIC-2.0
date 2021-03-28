#-----------------------------------------------------------------------
# This program reads real 1d trajectory data
# written for 1D OpenMP PIC codes
# written by Viktor K. Decyk, UCLA
from __future__ import print_function
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
   reads real periodic 1d trajectory data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 1
   iudm = 19; iuv = 12
   sname = numpy.array(["ELECTRON","ION     "],dtype=str)

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

# determine which trajectory diagnostics are available
   cmfield1.readtrdiags1(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   ndt = numpy.zeros((1),int_type,'F')
   nst = numpy.zeros((1),int_type,'F')
   nmv = numpy.zeros((1),int_type,'F')
   ndimp = numpy.zeros((1),int_type,'F')
   nprobt = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')
   ierr = numpy.zeros((1),int_type,'F')

# select diagnostic
   m = numpy.sum(nscalars)
   if (m==1):
      for i in range(0,ns):
         if (nscalars[i]==1):
            n = i + 1
            break
   else:
      print ("no trajectory diagnostic files found")
      n = 0
# exit procedure
      if (n==0):
         cmfield1.closeff1(iudm)
         return

   print (" trajectory diagnostic selected")

# open graphics device
   nplot = 1
   ierr[0] = graf1.open_graphs(nplot)

# return parameters for selected trajectory diagnostic:
# nts, ndt, nst, nmv, ndimp, nprobt, nrec, fname
   cmfield1.trdiagparams1(iudm,n,nts,ndt,nst,nmv,ndimp,nprobt,mrec,
                          fname)
   nrec = mrec[0]

   tname = numpy.str.rstrip(sname[ndt[0]-1])
   print (tname, " trajectory available")

# nmvf = size of velocity distribution array
   nmvf = 2*nmv[0] + 2

# allocate trajectory arrays
   if ((nst[0]==1) or (nst[0]==2)):
      if ("partt" not in globals()):
         partt = numpy.empty((ndimp[0],nprobt[0]),float_type,'F')
      if ("partd" not in globals()):
         partd = numpy.empty((nrec,ndimp[0],nprobt[0]),float_type,'F')
# allocate trajectory distribution array
   elif (nst[0]==3):
      if ("fvmtp" not in globals()):
         fvmtp = numpy.empty((in1.ndim,3),float_type,'F')
      if ("fvtp" not in globals()):
         fvtp = numpy.empty((nmvf,in1.ndim),float_type,'F')
      if ("fetp" not in globals()):
         fetp = numpy.empty((nmvf,0),float_type,'F')
      if ("wk" not in globals()):
         wk = numpy.empty((1),float_type,'F')
   dt = in1.dt*float(nts[0])

# open stream file for distribution field
   cmfield1.fsopen1(iuv,fname)

# nrec = number of complete records
   print ("records found: nrec = ", nrec)

# read and display data
   for ii in range(0,nrec):
      it = nts[0]*ii
      time = dt*float(ii)
# read real trajectory fields
      if ((nst[0]==1) or (nst[0]==2)):
         cmfield1.freadtr1(iuv,partt,ndimp,nprobt)
         partd[ii,:,:] = partt
      elif (nst[0]==3):
         cmfield1.freadfv1(iuv,fvmtp,fvtp,fetp,wk,in1.ndim,nmvf,
                           in1.ndim,0)
         gname = numpy.str.rstrip(sname[ndt[0]-1])
         graf1.displayfv1(fvtp,fvmtp,gname,it,nmv,2,ierr)
         if (ierr[0]==1):
            break

# display time history of trajectories
   if ((nst[0]==1) or (nst[0]==2)):
# display vx
      ix = 2
      graf1.displaytr1(partd,0.0,dt,nrec,ix,999,ierr)

   if ("partt" in globals()):
      partt = None
   if ("partd" in globals()):
      partd = None
   if ("fvmtp" in globals()):
       fvmtp = None
   if ("fvtp" in globals()):
      fvtp = None
   if ("fetp" in globals()):
      fetp = None
   if ("wk" in globals()):
      wk = None

   graf1.reset_graphs()
   cmfield1.closeff1(iudm)
   cmfield1.closeff1(iuv)
   graf1.close_graphs()
   print()

if (__name__=="__main__"):
   main(-1)
