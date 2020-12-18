#-----------------------------------------------------------------------
# This program reads real 3d trajectory data
# written for 3D MPI/OpenMP PIC codes
# written by Viktor K. Decyk, UCLA
from __future__ import print_function
import sys
import math
import numpy

sys.path.append('./ReadFiles')
from libcmfield3 import *

if (sys.version_info.major==2):
   input = raw_input

#-----------------------------------------------------------------------
def main(idrun):
   """
   reads real periodic 3d trajectory data
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
   fname = "diag3." + cdrun
   cmfield3.ffopen3(iudm,fname)

# nscalars = table of available diagnostics
   nscalars = numpy.zeros((ns),int_type,'F')

# determine which trajectory diagnostics are available
   cmfield3.readtrdiags3(iudm,nscalars)

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
         cmfield3.closeff3(iudm)
         return

   print (" trajectory diagnostic selected")

   nts = numpy.zeros((1),int_type,'F')
   ndt = numpy.zeros((1),int_type,'F')
   nst = numpy.zeros((1),int_type,'F')
   nmv = numpy.zeros((1),int_type,'F')
   ndimp = numpy.zeros((1),int_type,'F')
   nprobt = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')

# return parameters for selected trajectory diagnostic:
# nts, ndt, nst, nmv, ndimp, nprobt, nrec, fname
   cmfield3.trdiagparams3(iudm,n,nts,ndt,nst,nmv,ndimp,nprobt,mrec,fname)
   nrec = mrec[0]

   print (numpy.str.rstrip(sname[ndt[0]-1]), " trajectory available")

# nmvf = size of velocity distribution array
   nmvf = 2*nmv[0] + 2

# allocate trajectory arrays
   if ((nst[0]==1) or (nst[0]==2)):
      if ("partt" not in globals()):
         partt = numpy.empty((ndimp[0],nprobt[0]),float_type,'F')
      if ("partd" not in globals()):
         partd = numpy.empty((nrec,ndimp[0],nprobt[0]),float_type,'F')
# allocate trajectory distribution arrays
   elif (nst[0]==3):
      if ("fvmtp" not in globals()):
         fvmtp = numpy.empty((in3.ndim,3),float_type,'F')
      if ("fvtp" not in globals()):
         fvtp = numpy.empty((nmvf,in3.ndim),float_type,'F')
      if ("fetp" not in globals()):
         fetp = numpy.empty((nmvf,0),float_type,'F')
      if ("wk" not in globals()):
         wk = numpy.empty((1),float_type,'F')
   dt = in3.dt*float(nts[0])

# open stream file for trajectory field
   cmfield3.fsopen3(iuv,fname)

# nrec = number of complete records
   print ("records found: nrec = ", nrec)

# read and display data
   for ii in range(0,nrec):
      it = nts[0]*ii
      time = dt*float(ii)
# read real trajectory fields
      if ((nst[0]==1) or (nst[0]==2)):
         cmfield3.freadtr3(iuv,partt,ndimp,nprobt)
         partd[ii,:,:] = partt
         print (sname[ndt[0]-1], " trajectory:it,time=",it,time)
# display cartesian distribution
      elif (nst[0]==3):
         cmfield3.freadfv3(iuv,fvmtp,fvtp,fetp,wk,in3.ndim,nmvf,in3.ndim,0)
         print (sname[ndt[0]-1], " distribution:it,time=",it,time)

   cmfield3.closeff3(iuv)
   print()

if (__name__=="__main__"):
   main(-1)
