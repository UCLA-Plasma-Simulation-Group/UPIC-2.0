#-----------------------------------------------------------------------
# This program reads real periodic 3d fluid data
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
   reads real periodic 3d fluid data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 2
   iudm = 19; iuv = 12
   dname = numpy.array(["ELECTRON FLUID MOMENTS","ION FLUID MOMENTS     "],
                       dtype=str)
   sname = numpy.array(["ELECTRON","ION     "],dtype=str)
   ename = numpy.array([" DENSITY        "," VELOCITY FIELD ",
                        " PRESSURE TENSOR"," ENERGY         ",
                        " HEAT FLUX      "],dtype=str)

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

# determine which fluid diagnostics are available
   cmfield3.readfldiags3(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   npro = numpy.zeros((1),int_type,'F')
   nprd = numpy.zeros((1),int_type,'F')
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
         print ("no fluid diagnostic files found")
         n = 0
# exit procedure
      if (n==0):
         if ("sfield" in globals()):
            sfield = None
         if ("fms" in globals()):
            sfield = None
         if ("vfield" in globals()):
            ufield = None
         if ("vfield" in globals()):
            ufield = None
         cmfield3.closeff3(iudm)
         return

      print (numpy.str.rstrip(dname[n-1]), " diagnostic selected")

# return parameters for selected fluid diagnostic:
# nts, npro, nprd, nrec, fname
      cmfield3.fldiagparams3(iudm,n,nts,npro,nprd,mrec,fname)
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

# allocate vector array
      if ("fms" not in globals()):
         fms = numpy.empty((nprd[0],nx,nyv,nzv),float_type,'F')
      if ("sfield" not in globals()):
         sfield = numpy.empty((nx,nyv,nzv),float_type,'F')
      if (in3.ndim==3):
         if ("vfield" not in globals()):
            vfield = numpy.empty((in3.ndim,nx,nyv,nzv),float_type,'F')
         if (npro[0] > 2):
            if ("ufield" not in globals()):
               ufield = numpy.empty((2*in3.ndim,nx,nyv,nzv),float_type,'F')
      dt = in3.dt*float(nts[0])

# open stream file for vector field
      cmfield3.fsopen3(iuv,fname)

# nrec = number of complete records
      nrec = int(nrec/(kyb*kzb))
      print ("records found: nrec = ", nrec)

# read and transpose vector data
      for ii in range(0,nrec):
# read real vector field
         cmfield3.freadv3(iuv,fms,nprd,nx,kyp,kyb,kzp,kzb)
         it = nts[0]*ii
         time = dt*float(ii)
# show time
         print (sname[n-1], " it,time=",it,time)
      cmfield3.closeff3(iuv)
      print()

if (__name__=="__main__"):
   main(-1)
