#-----------------------------------------------------------------------
# This program reads real periodic 3d velocity data
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
   reads real periodic 3d velocity data
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
   fname = "diag3." + cdrun
   cmfield3.ffopen3(iudm,fname)

# nscalars = table of available diagnostics
   nscalars = numpy.zeros((ns),int_type,'F')

# determine which velocity diagnostics are available
   cmfield3.readfvdiags3(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   nmv = numpy.zeros((1),int_type,'F')
   nfvd = numpy.zeros((1),int_type,'F')
   nfed = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')

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
         cmfield3.closeff3(iudm)
         return

      print (numpy.str.rstrip(dname[n-1]), " diagnostic selected")

# return parameters for selected velocity diagnostic:
# nts, nmv, nfvd, nfed, nrec, fname
      cmfield3.fvdiagparams3(iudm,n,nts,nmv,nfvd,nfed,mrec,fname)
      nrec = mrec[0]

# nmv21 = number of velocity bins
      nmv21 = 2*nmv[0] + 1; nmvf = nmv21 + 1

# allocate velocity distribution arrays
      if ("fvm" not in globals()):
         fvm = numpy.empty((in3.ndim,3),float_type,'F')
      if ("fv" not in globals()):
         fv = numpy.empty((nmvf,nfvd[0]),float_type,'F')
      if ("fe" not in globals()):
         fe = numpy.empty((nmvf,nfed[0]),float_type,'F')
      if ("wk" not in globals()):
         wk = numpy.empty((1),float_type,'F')
      dt = in3.dt*float(nts[0])

# open stream file for velocity field
      cmfield3.fsopen3(iuv,fname)

# nrec = number of complete records
      print ("records found: nrec = ", nrec)

# read and display data
      for ii in range(0,nrec):
# read real velocity fields
         cmfield3.freadfv3(iuv,fvm,fv,fe,wk,in3.ndim,nmvf,nfvd,nfed)
         it = nts[0]*ii
         time = dt*float(ii)
         if (nfvd[0] > 0):
# display cylindrical distribution
            if (nfvd[0] < in3.ndim):
               print (sname[n-1], " cylindrical:it,time=",it,time)
# display cartesian distribution
            else:
               print (sname[n-1], " cartesian:it,time=",it,time)
# display energy distribution
         if (nfed[0] > 0):
            print (sname[n-1], " energy:it,time=",it,time)
      cmfield3.closeff3(iuv)
      print()

if (__name__=="__main__"):
   main(-1)
