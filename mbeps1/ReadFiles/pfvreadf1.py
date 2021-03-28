#-----------------------------------------------------------------------
# This program reads real 1d velocity data
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
   reads real periodic 1d velocity data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 2
   iudm = 19; iuv = 12
   dname = numpy.array(["ELECTRON VELOCITY   ","ION VELOCITY        "],
                       dtype=str)
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

# determine which velocity diagnostics are available
   cmfield1.readfvdiags1(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   nmv = numpy.zeros((1),int_type,'F')
   nfvd = numpy.zeros((1),int_type,'F')
   nfed = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')
   gname = numpy.array([""],'S32')
   ierr = numpy.zeros((1),int_type,'F')

# open graphics device
   nplot = 4
   ierr[0] = graf1.open_graphs(nplot)

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
         print ("no velocity-space diagnostic files found")
         n = 0
# exit procedure
      if (n==0):
         if ("fvm" in globals()):
            fvm = None
         if ("fv" in globals()):
            fv = None
         if ("fe" in globals()):
            fe = None
         if ("wk" in globals()):
            wk = None
         cmfield1.closeff1(iudm)
         graf1.close_graphs()
         return

      tname = numpy.str.rstrip(dname[n-1])
      print (tname, " diagnostic selected")

# return parameters for selected velocity diagnostic:
# nts, nmv, nfvd, nfed, nrec, fname
      cmfield1.fvdiagparams1(iudm,n,nts,nmv,nfvd,nfed,mrec,fname)
      nrec = mrec[0]

# nmv21 = number of velocity bins
      nmv21 = 2*nmv[0] + 1; nmvf = nmv21 + 1

# allocate velocity distribution arrays
      if ("fvm" not in globals()):
         fvm = numpy.empty((in1.ndim,3),float_type,'F')
      if ("fv" not in globals()):
         fv = numpy.empty((nmvf,nfvd[0]),float_type,'F')
      if ("fe" not in globals()):
         fe = numpy.empty((nmvf,nfed[0]),float_type,'F')
      if ("wk" not in globals()):
         wk = numpy.empty((1),float_type,'F')
      dt = in1.dt*float(nts[0])

# open stream file for distribution field
      cmfield1.fsopen1(iuv,fname)

# nrec = number of complete records
      print ("records found: nrec = ", nrec)

# read and display data
      for ii in range(0,nrec):
# read real velocity fields
         cmfield1.freadfv1(iuv,fvm,fv,fe,wk,in1.ndim,nmvf,nfvd,nfed)
         it = nts[0]*ii
         time = dt*float(ii)
         if (nfvd[0] > 0):
# display cylindrical distribution
            if ((in1.ndim==3) and (nfvd[0] < in1.ndim)):
               gname = numpy.str.rstrip(sname[n-1])
               graf1.displayfvb1(fv,fvm,gname,it,nmv,2,ierr)
               if (ierr[0]==1):
                  break
# display cartesian distribution
            else:
               gname = numpy.str.rstrip(sname[n-1])
               graf1.displayfv1(fv,fvm,gname,it,nmv,2,ierr)
               if (ierr[0]==1):
                  break
# display energy distribution
            if (nfed[0] > 0):
               gname = numpy.str.rstrip(sname[n-1])
               graf1.displayfe1(fe,wk,gname,it,nmv,ierr)
               if (ierr[0]==1):
                  break

      graf1.reset_graphs()
      graf1.reset_nplot(4,ierr)
      cmfield1.closeff1(iuv)
      print()

if (__name__=="__main__"):
   main(-1)
