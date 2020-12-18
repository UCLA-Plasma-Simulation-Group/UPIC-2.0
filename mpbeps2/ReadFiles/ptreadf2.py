#-----------------------------------------------------------------------
# This program reads real 2d trajectory data
# written for 2D MPI/OpenMP PIC codes
# written by Viktor K. Decyk, UCLA
from __future__ import print_function
import sys
import math
import numpy

sys.path.append('./ReadFiles')
from fgraf import *
from libcmfield2 import *

if (sys.version_info.major==2):
   input = raw_input

#-----------------------------------------------------------------------
def main(idrun):
   """
   reads real periodic 2d trajectory data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 1
   iudm = 19; iuv = 12
   nplot = 1

   sname = numpy.array(["ELECTRON","ION     "],dtype=str)

# create string from idrun
   if (idrun < 0):
      cdrun = "Unknown"
      while (cdrun.isdigit() == False):
         cdrun = input("enter integer idrun: ")
      idrun = int(cdrun)   
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
   tname = numpy.array([""],'S32')
   gname = numpy.array([""],'S32')
   ierr = numpy.zeros((1),int_type,'F')

# open graphics device
   ierr[0] = graf2.open_graphs2(nplot)

# return parameters for selected trajectory diagnostic:
# nts, ndt, nst, nmv, ndimp, nprobt, nrec, fname
   cmfield2.trdiagparams2(iudm,n,nts,ndt,nst,nmv,ndimp,nprobt,mrec,fname)
   nrec = mrec[0]
   nplot = 1

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
# allocate trajectory distribution arrays
   elif (nst[0]==3):
      if ("fvmtp" not in globals()):
         fvmtp = numpy.empty((in2.ndim,3),float_type,'F')
      if ("fvtp" not in globals()):
         fvtp = numpy.empty((nmvf,in2.ndim),float_type,'F')
      if ("fetp" not in globals()):
         fetp = numpy.empty((nmvf,0),float_type,'F')
      if ("wk" not in globals()):
         wk = numpy.empty((1),float_type,'F')
   dt = in2.dt*float(nts[0])

# open stream file for trajectory field
   cmfield2.fsopen2(iuv,fname)

# nrec = number of complete records
   print ("records found: nrec = ", nrec)

# read and display data
   for ii in range(0,nrec):
      it = nts[0]*ii
      time = dt*float(ii)
# read real trajectory fields
      if ((nst[0]==1) or (nst[0]==2)):
         cmfield2.freadtr2(iuv,partt,ndimp,nprobt)
         partd[ii,:,:] = partt
         print (sname[ndt[0]-1], " trajectory:it,time=",it,time)
# display cartesian distribution
      elif (nst[0]==3):
         cmfield2.freadfv2(iuv,fvmtp,fvtp,fetp,wk,in2.ndim,nmvf,in2.ndim,0)
         print (sname[ndt[0]-1], " distribution:it,time=",it,time)
         gname = numpy.str.rstrip(sname[ndt[0]-1])
         graf1.displayfv1(fvtp,fvmtp,gname,it,nmv,2,ierr)
         if (ierr[0]==1):
            break

# display time history of trajectories
   if ((nst[0]==1) or (nst[0]==2)):
# display vx
      ix = 3
      graf1.displaytr1(partd,0.0,dt,nrec,ix,999,ierr)

   cmfield2.closeff2(iuv)
   graf2.close_graphs2()
   print()

if (__name__=="__main__"):
   main(-1)
