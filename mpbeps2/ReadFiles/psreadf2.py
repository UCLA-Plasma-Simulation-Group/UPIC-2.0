#-----------------------------------------------------------------------
# This program reads real 2d phase space data
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
   reads real periodic 2d phase space data
   idrun = run identifier for current run, -1 if unknown
   """
   int_type = numpy.int32
   double_type = numpy.float64
   float_type = numpy.float32
   complex_type = numpy.complex64

   ns = 2
   iudm = 19; iuv = 12
   dname = numpy.array(["ELECTRON PHASE SPACE","ION PHASE SPACE     "],
                       dtype=str)
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

# determine which phase space diagnostics are available
   cmfield2.readpsdiags2(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   nmv = numpy.zeros((1),int_type,'F')
   nsxb = numpy.zeros((1),int_type,'F')
   nsyb = numpy.zeros((1),int_type,'F')
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
         print ("no phase space diagnostic files found")
         n = 0
# exit procedure
      if (n==0):
         if ("fvs" in globals()):
            fvs = None
         if ("pvs" in globals()):
            pvs = None
         if ("ps" in globals()):
            ps = None
         if ("psx" in globals()):
            psx = None
         if ("psy" in globals()):
            psy = None
         cmfield2.closeff2(iudm)
         graf2.close_graphs2()
         return

      tname = numpy.str.rstrip(dname[n-1])      
      print (tname, " diagnostic selected")

# return parameters for selected phase space diagnostic:
# nts, nmv, nsxb, nsyb, nrec, fname
      cmfield2.psdiagparams2(iudm,n,nts,nmv,nsxb,nsyb,mrec,fname)
      nrec = mrec[0]

# nx/ny = number of global grid points in x/y direction
      nx = int(math.pow(2,in2.indx)); ny = int(math.pow(2,in2.indy))
# nmv21 = number of velocity bins
      nmv21 = 2*nmv[0] + 1; nmvf = nmv21 + 1

# allocate velocity distribution arrays
      if ("fvs" not in globals()):
         fvs = numpy.empty((nmvf,in2.ndim,nsxb[0],nsyb[0]),float_type,
                           'F')
      if ("pvs" not in globals()):
         pvs = numpy.empty((nmvf,in2.ndim,max(nsxb[0],nsyb[0])),
                           float_type,'F')
      if ("ps" not in globals()):
         ps = numpy.empty((nmvf,max(nsxb[0],nsyb[0])),float_type,'F')
      if ("psx" not in globals()):
         psx = numpy.empty((nsxb[0],nmvf),float_type,'F')
      if ("psy" not in globals()):
         psy = numpy.empty((nsyb[0],nmvf),float_type,'F')
      dt = in2.dt*float(nts[0])

# open stream file for phase space field
      cmfield2.fsopen2(iuv,fname)

# nrec = number of complete records
      print ("records found: nrec = ", nrec)

# read and display data
      for ii in range(0,nrec):
# read real velocity distribution fields
         cmfield2.freadps2(iuv,fvs,in2.ndim,nmvf,nsxb[0],nsyb[0])
         it = nts[0]*ii
         time = dt*float(ii)
         print (sname[n-1], " it,time=",it,time)
# sum over y
         pvs = numpy.sum(fvs,axis=3)
# select x-vx
         xmax = float(nx)
         ymax = fvs[nmv21,0,0,0]; ymin = -ymax
         ps[:,:] = pvs[:,0,:nsxb[0]]
         psx[:,:] = numpy.transpose(ps)
# display real space data
         gname = numpy.str.rstrip(sname[n-1])+" X-VX"
         graf2.dscalerl2(psx,gname,0.0,xmax,ymin,ymax,it,999,1,nsxb[0],
                         nmv21,
                   ierr)
         if (ierr[0]==1):
            break
# select x-vy
         ymax = fvs[nmv21,1,0,0]; ymin = -ymax
         ps[:,:] = pvs[:,1,:nsxb[0]]
         psx[:,:] = numpy.transpose(ps)
# display real space data
         gname = numpy.str.rstrip(sname[n-1])+" X-VY"
         graf2.dscalerl2(psx,gname,0.0,xmax,ymin,ymax,it,999,1,nsxb[0],
                         nmv21,
                   ierr)
         if (ierr[0]==1):
            break
# sum over x
         pvs = numpy.sum(fvs,axis=2)
# select y-vx
         xmax = float(ny)
         ymax = fvs[nmv21,0,0,0]; ymin = -ymax
         ps[:,:] = pvs[:,0,:nsyb[0]]
         psy[:,:] = numpy.transpose(ps)
# display real space data
         gname = numpy.str.rstrip(sname[n-1])+" Y-VX"
         graf2.dscalerl2(psy,gname,0.0,xmax,ymin,ymax,it,999,1,nsxb[0],
                         nmv21,ierr)
         if (ierr[0]==1):
            break
# select y-vy
         ymax = fvs[nmv21,1,0,0]; ymin = -ymax
         ps[:,:] = pvs[:,1,:nsyb[0]]
         psy[:,:] = numpy.transpose(ps)
# display real space data
         gname = numpy.str.rstrip(sname[n-1])+" Y-VY"
         graf2.dscalerl2(psy,gname,0.0,xmax,ymin,ymax,it,999,1,nsxb[0],
                         nmv21,ierr)
         if (ierr[0]==1):
            break

      cmfield2.closeff2(iuv)
      print()

if (__name__=="__main__"):
   main(-1)
