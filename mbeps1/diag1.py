#-----------------------------------------------------------------------
# This program reads 1d data written by 1D OpenMP PIC codes
# written by Viktor K. Decyk, UCLA
from __future__ import print_function
import sys
import math
import numpy

sys.path.append('./ReadFiles')
from libcmfield1 import *
import pcreadf1
import pvcreadf1
import pflreadf1
import pfvreadf1
import ptreadf1
import psreadf1

if (sys.version_info.major==2):
   input = raw_input

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

ms = 6; ns = 6
iudm = 19
nplot = 1
dname = numpy.array(["COMPLEX SCALAR FIELDS","COMPLEX VECTOR FIELDS",
                     "REAL FLUID DATA      ","VELOCITY DATA      ",
                     "TRAJECTORY DATA      ","PHASE SPACE DATA   "],
                     dtype=str)

# create string from idrun
cdrun = "Unknown"
while (cdrun.isdigit() == False):
   cdrun = input("enter integer idrun: ")
idrun = int(cdrun)
fname = "diag1." + cdrun
cmfield1.ffopen1(iudm,fname)

# mscalars = table of available diagnostic types
mscalars = numpy.zeros((ms),int_type,'F')
# nscalars = table of available diagnostics
nscalars = numpy.zeros((ns),int_type,'F')

# determine which scalar diagnostics are available
cmfield1.readsdiags1(iudm,nscalars)
mscalars[0] = numpy.sum(nscalars)
# determine which vector diagnostics are available
nscalars[:] = 0
cmfield1.readvdiags1(iudm,nscalars)
mscalars[1] = numpy.sum(nscalars)
# determine which fluid diagnostics are available
nscalars[:] = 0
cmfield1.readfldiags1(iudm,nscalars)
mscalars[2] = numpy.sum(nscalars)
# determine which velocity diagnostics are available
nscalars[:] = 0
cmfield1.readfvdiags1(iudm,nscalars)
mscalars[3] = numpy.sum(nscalars)
# determine which trajectory diagnostics are available
nscalars[:] = 0
cmfield1.readtrdiags1(iudm,nscalars)
mscalars[4] = numpy.sum(nscalars)
# determine which phase space diagnostics are available
nscalars[:] = 0
cmfield1.readpsdiags1(iudm,nscalars)
mscalars[5] = numpy.sum(nscalars)

# select diagnostic
m = numpy.sum(mscalars)
while True:
   if (m > 0):
      n = -1
      while True:
         if (n < 0):
            for i in range(0,ms):
               if (mscalars[i] > 0):
                  print ("enter ", i+1," for", numpy.str.rstrip(dname[i]))
            print ("enter ", 0, " for EXIT")
            c = input("")
            if (c.isdigit()):
               n = int(c)
            if (n==0):
               exit(1)
            if ((n >= 1) and (n <= ms)):
               if (mscalars[n-1]==0):
                  n = -1
            else:
               n = -1
            if (n > 0):
               break
            print ("invalid entry, try again or enter 0 to quit")
   else:
      print ("no diagnostic data found")
      exit(1)

   print (numpy.str.rstrip(dname[n-1]), " diagnostic type selected")

   if (n==1):
      pcreadf1.main(idrun)
   elif (n==2):
      pvcreadf1.main(idrun)
   elif (n==3):
      pflreadf1.main(idrun)
   elif (n==4):
      pfvreadf1.main(idrun)
   elif (n==5):
      ptreadf1.main(idrun)
   elif (n==6):
      psreadf1.main(idrun)

cmfield1.closeff1(iudm)

