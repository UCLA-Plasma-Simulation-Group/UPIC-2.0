#-----------------------------------------------------------------------
# This program reads real periodic 2d scalar data
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
   reads real periodic 2d scalar data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 3
   iudm = 19; ius = 11
   nplot = 1

   dname = numpy.array(["POTENTIAL       ","ELECTRON DENSITY",
                        "ION DENSITY     "],dtype=str)

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

# determine which scalar diagnostics are available
   cmfield2.readsdiags2(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   modesx = numpy.zeros((1),int_type,'F')
   modesy = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')
   tname = numpy.array([""],'S32')
   ierr = numpy.zeros((1),int_type,'F')

# open graphics device
   ierr[0] = graf2.open_graphs2(nplot)
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
         print ("no scalar diagnostic files found")
         n = 0
# exit procedure
      if (n==0):
         if ("sfield" in globals()):
            sfield = None
         cmfield2.closeff2(iudm)
         graf2.close_graphs2()
         return

      tname = numpy.str.rstrip(dname[n-1])
      print (tname, " diagnostic selected")

# return parameters for selected scalar diagnostic:
# nts, modesx, modesy, nrec, fname
      cmfield2.sdiagparams2(iudm,n,nts,modesx,modesy,mrec,fname)
      nrec = mrec[0]

# nx/ny = number of global grid points in x/y direction
      nx = int(math.pow(2,in2.indx)); ny = int(math.pow(2,in2.indy))
# kyp = number of real grids in each field partition in y direction
      kyp = int((ny - 1)/in2.nvp) + 1
# kyb = minimum number of processors in distributed array
      kyb = int((ny - 1)/kyp) + 1
# nyv = second dimension of scalar field array, >= ny
      nyv = kyp*kyb

# allocate scalar array
      if ("sfield" not in globals()):
         sfield = numpy.empty((nx,nyv),float_type,'F')
      dt = in2.dt*float(nts[0])

# open stream file for scalar field
      cmfield2.fsopen2(ius,fname)

# nrec = number of complete records
      nrec = int(nrec/kyb)
      print ("records found: nrec = ", nrec)

# read and display data
      for ii in range(0,nrec):
# read real scalar field
         cmfield2.fread2(ius,sfield,nx,nyv)
         it = nts[0]*ii
         time = dt*float(ii)
# show time
         print ("it,time=",it,time)
# read and display data
         graf2.dscaler2(sfield,tname,it,999,0,1,nx,ny,ierr)
         if (ierr[0]==1):
            break

      cmfield2.closeff2(ius)
      print()

if (__name__=="__main__"):
   main(-1)
