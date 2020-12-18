#-----------------------------------------------------------------------
# This program reads real periodic 3d scalar data
# written by 3D MPI/OpenMP PIC codes
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
   reads real periodic 3d scalar data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 3
   iudm = 19; ius = 11
   dname = numpy.array(["POTENTIAL       ","ELECTRON DENSITY",
                        "ION DENSITY     "],dtype=str)

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

# determine which scalar diagnostics are available
   cmfield3.readsdiags3(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   modesx = numpy.zeros((1),int_type,'F')
   modesy = numpy.zeros((1),int_type,'F')
   modesz = numpy.zeros((1),int_type,'F')
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
         print ("no scalar diagnostic files found")
         n = 0
# exit procedure
      if (n==0):
         if ("sfield" in globals()):
            sfield = None
         cmfield3.closeff3(iudm)
         return

      print (numpy.str.rstrip(dname[n-1]), " diagnostic selected")

# return parameters for selected scalar diagnostic:
# nts, modesx, modesy, modesz, nrec, fname
      cmfield3.sdiagparams3(iudm,n,nts,modesx,modesy,modesz,mrec,fname)
      nrec = mrec[0]

# nx/ny/nz = number of global grid points in x/y/z direction
      nx = int(math.pow(2,in3.indx)); ny = int(math.pow(2,in3.indy))
      nz = int(math.pow(2,in3.indz))
# kyp/kzp = number of real grids in each field partition in y/z
      kyp = int((ny - 1)/in3.nvpy) + 1; kzp = int((nz - 1)/in3.nvpz) + 1
# kyb/kzb = minimum number of processors in distributed array in y/z
      kyb = int((ny - 1)/kyp) + 1; kzb = int((nz - 1)/kzp) + 1
# nyv = second dimension of scalar field array, >= ny
# nzv = third dimension of scalar field array, >= nz
      nyv = kyp*kyb; nzv = kzp*kzb

# allocate scalar array
      if ("sfield" not in globals()):
         sfield = numpy.empty((nx,nyv,nzv),float_type,'F')
      dt = in3.dt*float(nts[0])

# open stream file for scalar field
      cmfield3.fsopen3(ius,fname)

# nrec = number of complete records
      nrec = int(nrec/(kyb*kzb))
      print ("records found: nrec = ", nrec)

# read and transpose scalar data
      for ii in range(0,nrec):
# read real scalar field
         cmfield3.fread3(ius,sfield,nx,kyp,kyb,kzp,kzb)
         it = nts[0]*ii
         time = dt*float(ii)
# show time
         print ("it,time=",it,time)
      cmfield3.closeff3(ius)
      print()

if (__name__=="__main__"):
   main(-1)
