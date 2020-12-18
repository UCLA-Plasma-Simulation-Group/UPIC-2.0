#-----------------------------------------------------------------------
# This program reads real periodic 2d fluid data
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
   reads real periodic 2d fluid data
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
   fname = "diag2." + cdrun
   cmfield2.ffopen2(iudm,fname)

# nscalars = table of available diagnostics
   nscalars = numpy.zeros((ns),int_type,'F')

# determine which fluid diagnostics are available
   cmfield2.readfldiags2(iudm,nscalars)
   
   nts = numpy.zeros((1),int_type,'F')
   npro = numpy.zeros((1),int_type,'F')
   nprd = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')
   tname = numpy.array([""],'S32')
   gname = numpy.array([""],'S32')
   ierr = numpy.zeros((1),int_type,'F')

# open graphics device
   nplot = 4
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
         cmfield2.closeff2(iudm)
         graf2.close_graphs2()
         return

      tname = numpy.str.rstrip(dname[n-1])
      print (tname, " diagnostic selected")

# return parameters for selected fluid diagnostic:
# nts, npro, nprd, nrec, fname
      cmfield2.fldiagparams2(iudm,n,nts,npro,nprd,mrec,fname)
      nrec = mrec[0]

# nx/ny = number of global grid points in x/y direction
      nx = int(math.pow(2,in2.indx)); ny = int(math.pow(2,in2.indy))
# kyp = number of real grids in each field partition in y direction
      kyp = int((ny - 1)/in2.nvp) + 1
# kyb = minimum number of processors in distributed array
      kyb = int((ny - 1)/kyp) + 1
# nyv = third dimension of vector field array, >= ny
      nyv = kyp*kyb

# allocate vector array
      if ("fms" not in globals()):
         fms = numpy.empty((nprd[0],nx,nyv),float_type,'F')
      if ("sfield" not in globals()):
         sfield = numpy.empty((nx,nyv),float_type,'F')
      if (in2.ndim==3):
         if ("vfield" not in globals()):
            vfield = numpy.empty((in2.ndim,nx,nyv),float_type,'F')
         if (npro[0] > 2):
            if ("ufield" not in globals()):
               ufield = numpy.empty((2*in2.ndim,nx,nyv),float_type,'F')
      dt = in2.dt*float(nts[0])

# open stream file for vector field
      cmfield2.fsopen2(iuv,fname)

# nrec = number of complete records
      nrec = int(nrec/kyb)
      print( "records found: nrec = ", nrec)

# read and display data
      for ii in range(0,nrec):
# read real vector field
         cmfield2.freadv2(iuv,fms,nprd,nx,nyv)
         it = nts[0]*ii
         time = dt*float(ii)
# show time
         print (sname[n-1], " it,time=",it,time)
# electrostatic case
         if (in2.ndim==2):
# display density in real space
            if (npro[0] > 0):
               sfield[:] = fms[0,:]
               gname = tname + numpy.str.rstrip(ename[0])
               graf2.dscaler2(sfield,gname,it,999,1,1,nx,ny,ierr)
               if (ierr[0]==1):
                  break
# display velocity field in real space
            if (npro[0] > 1):
               sfield[:] = fms[1,:]
               gname = tname + numpy.str.rstrip(ename[1])
               graf2.dscaler2(sfield,gname,it,999,0,1,nx,ny,ierr)
               if (ierr[0]==1):
                  break
# display pressure tensor in real space
            if (npro[0] > 2):
               sfield[:] = fms[2,:]
               gname = tname + numpy.str.rstrip(ename[2])
               graf2.dscaler2(sfield,gname,it,999,0,1,nx,ny,ierr)
               if (ierr[0]==1):
                  break
# display electron heat flux in real space
            if (npro[0]==4):
               sfield[:] = fms[4,:]
               gname = tname + numpy.str.rstrip(ename[4])
               graf2.dscaler2(sfield,gname,it,999,0,1,nx,ny,ierr)
               if (ierr[0]==1):
                  break
# electromagnetic case
         elif (in2.ndim==3):
# display density in real space
            if (npro[0] > 0):
               sfield[:] = fms[0,:]
               gname = tname + numpy.str.rstrip(ename[0])
               graf2.dscaler2(sfield,gname,it,999,1,1,nx,ny,ierr)
               if (ierr[0]==1):
                  break
# display velocity field in real space
            if (npro[0] > 1):
               vfield[:,:] = fms[1:4,:]
               gname = tname + numpy.str.rstrip(ename[1])
               graf2.dvector2(vfield,gname,it,999,0,1,1,nx,ny,ierr)
               if (ierr[0]==1):
                  break
# display pressure tensor in real space
            if (npro[0] > 2):
               ufield[:,:] = fms[4:10,:]
               gname = tname + numpy.str.rstrip(ename[2])
               graf2.dvector2(ufield,gname,it,999,0,1,1,nx,ny,ierr)
               if (ierr[0]==1):
                  break
# display electron heat flux in real space
            if (npro[0]==4):
               vfield[:,:] = fms[11:14,:]
               gname = tname + numpy.str.rstrip(ename[4])
               graf2.dvector2(vfield,gname,it,999,0,1,1,nx,ny,ierr)
               if (ierr[0]==1):
                  break

      cmfield2.closeff2(iuv)
      print()

if (__name__=="__main__"):
   main(-1)

