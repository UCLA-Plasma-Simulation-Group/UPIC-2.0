#-----------------------------------------------------------------------
# This program reads real periodic 1d fluid data
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
   reads real periodic 1d fluid data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 2
   iudm = 19; iuv = 12
   dname = numpy.array(["ELECTRON FLUID MOMENTS",
                        "ION FLUID MOMENTS     "],dtype=str)
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
   fname = "diag1." + cdrun
   cmfield1.ffopen1(iudm,fname)

# nscalars = table of available diagnostics
   nscalars = numpy.zeros((ns),int_type,'F')

# determine which fluid diagnostics are available
   cmfield1.readfldiags1(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   npro = numpy.zeros((1),int_type,'F')
   nprd = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')
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
         if ("ufield" in globals()):
            ufield = None
         cmfield1.closeff1(iudm)
         graf1.close_graphs()
         return

      tname = numpy.str.rstrip(dname[n-1])
      print (tname, " diagnostic selected")
   
# return parameters for selected fluid diagnostic:
# nts, npro, nprd, nrec, fname
      cmfield1.fldiagparams1(iudm,n,nts,npro,nprd,mrec,fname)
      nrec = mrec[0]

# nx = number of global grid points in x direction
      nx = int(math.pow(2,in1.indx))

# allocate vector arrays
      if ("sfield" not in globals()):
         sfield = numpy.empty((nx),float_type,'F')
      if ("fms" not in globals()):
         fms = numpy.empty((nprd[0],nx),float_type,'F')
      if (in1.ndim==3):
         if ("vfield" not in globals()):
            vfield = numpy.empty((in1.ndim,nx),float_type,'F')
         if (npro[0] > 2):
            if ("ufield" not in globals()):
               ufield = numpy.empty((2*in1.ndim,nx),float_type,'F')
      dt = in1.dt*float(nts[0])

# open stream file for vector field
      cmfield1.fsopen1(iuv,fname)

# nrec = number of complete records
      print ("records found: nrec = ", nrec)

# read and display data
      for ii in range(0,nrec):
# read real vector field
         cmfield1.freadv1(iuv,fms,nprd[0],nx)
         it = nts*ii
         time = dt*float(ii)
# electrostatic case
         if (in1.ndim==1):
# display density in real space
            if (npro[0] > 0):
               sfield[:] = fms[0,:]
               gname = numpy.str.rstrip(sname[n-1])+numpy.str.rstrip(ename[0])
               graf1.dscaler1(sfield,gname,it,999,1,nx,ierr)
               if (ierr[0]==1):
                  exit(1)
# display velocity field in real space
            if (npro[0] > 1):
               sfield[:] = fms[1,:]
               gname = numpy.str.rstrip(sname[n-1])+numpy.str.rstrip(ename[1])
               graf1.dscaler1(sfield,gname,it, 999,0,nx,ierr)
               if (ierr[0]==1):
                  exit(1)
# display pressure tensor in real space
            if (npro[0] > 2):
               sfield[:] = fms[2,:]
               gname = numpy.str.rstrip(sname[n-1])+numpy.str.rstrip(ename[2])
               graf1.dscaler1(sfield,gname,it,999,0,nx,ierr)
               if (ierr[0]==1):
                  exit(1)
# display electron heat flux in real space
            if (npro[0]==4):
               sfield[:] = fms[4,:]
               gname = numpy.str.rstrip(sname[n-1])+numpy.str.rstrip(ename[4])
               graf1.dscaler1(sfield,gname,it,999,0,nx,ierr)
               if (ierr[0]==1):
                  exit(1)
# electromagnetic case
         elif (in1.ndim==3):
# display density in real space
            if (npro[0] > 0):
               sfield[:] = fms[0,:]
               gname = numpy.str.rstrip(sname[n-1])+numpy.str.rstrip(ename[0])
               graf1.dscaler1(sfield,gname,it,999,1,nx,ierr)
               if (ierr[0]==1):
                  exit(1)
# display velocity field in real space
            if (npro[0] > 1):
               vfield[:,:] = fms[1:4,:]
               gname = numpy.str.rstrip(sname[n-1])+numpy.str.rstrip(ename[1])
               graf1.dvector1(vfield,gname,it,999,0,2,nx,ierr)
               if (ierr[0]==1):
                  exit(1)
# display pressure tensor in real space
            if (npro[0] > 2):
               ufield[:,:] = fms[4:10,:]
               gname = numpy.str.rstrip(sname[n-1])+numpy.str.rstrip(ename[2])
               graf1.dvector1(ufield,gname,it,999,0,2,nx,ierr)
               if (ierr[0]==1):
                  exit(1)
# display electron heat flux in real space
            if (npro[0]==4):
               vfield[:,:] = fms[11:14,:]
               gname = numpy.str.rstrip(sname[n-1])+numpy.str.rstrip(ename[4])
               graf1.dvector1(vfield,gname,it,999,0,2,nx,ierr)
               if (ierr[0]==1):
                  exit(1)

      graf1.reset_graphs()
      graf1.reset_nplot(4,ierr)
      cmfield1.closeff1(iuv)
      print()

if (__name__=="__main__"):
   main(-1)
