#-----------------------------------------------------------------------
# This program reads real 1d phase space data
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
  reads real periodic 1d phase space data
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
   fname = "diag1." + cdrun
   cmfield1.ffopen1(iudm,fname)
 
# nscalars = table of available diagnostics
   nscalars = numpy.zeros((ns),int_type,'F')

# determine which phase space diagnostics are available
   cmfield1.readpsdiags1(iudm,nscalars)

   nts = numpy.zeros((1),int_type,'F')
   nmv = numpy.zeros((1),int_type,'F')
   nsxb = numpy.zeros((1),int_type,'F')
   mrec = numpy.zeros((1),int_type,'F')
   fname = numpy.array([""],'S32')
   ierr = numpy.zeros((1),int_type,'F')

# open graphics device
   nplot = in1.ndim
   ierr[0] = graf1.open_graphs(nplot)
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
         cmfield1.closeff1(iudm)
         graf1.close_graphs()
         return

      tname = numpy.str.rstrip(dname[n-1])
      print (tname, " diagnostic selected")
 
# return parameters for selected phase space diagnostic:
# nts, nmv, nsxb, nrec, fname
      cmfield1.psdiagparams1(iudm,n,nts,nmv,nsxb,mrec,fname)
      nrec = mrec[0]

# nx = number of global grid points in x direction
      nx = int(math.pow(2,in1.indx))
# nmv21 = number of velocity bins
      nmv21 = 2*nmv[0] + 1; nmvf = nmv21 + 1

# allocate velocity distribution arrays
      if ("fvs" not in globals()):
         fvs = numpy.empty((nmvf,in1.ndim,nsxb[0]),float_type,'F')
      if ("pvs" not in globals()):
         pvs = numpy.empty((nmvf,in1.ndim,nsxb[0]),float_type,'F')
      if ("ps" not in globals()):
         ps = numpy.empty((nmvf,nsxb[0]),float_type,'F')
      if ("psx" not in globals()):
         psx = numpy.empty((nsxb[0],nmvf),float_type,'F')
      dt = in1.dt*float(nts[0])

# open stream file for phase space field
      cmfield1.fsopen1(iuv,fname)

# nrec = number of complete records
      print ("records found: nrec = ", nrec)

# read and display data
      for ii in range(0,nrec):
# read real velocity distribution fields
         cmfield1.freadps1(iuv,fvs,in1.ndim,nmvf,nsxb[0])
         it = nts[0]*ii
         time = dt*float(ii)
         pvs[:,:,:] = fvs
# electrostatic case
         if (in1.ndim==1):
# select x-vx
            ps[:,:] = pvs[:,0,:nsxb[0]]
            psx[:,:] = numpy.transpose(ps)
            xmax = float(nx)
            ymax = fvs[nmv21,0,0]; ymin = -ymax
# display real space data
            gname = numpy.str.rstrip(sname[n-1])+" X-VX"
            graf2.dscalerl2(psx,gname,0.0,xmax,ymin,ymax,it,999,1,
                            nsxb[0],nmv21,ierr)
            if (ierr[0]==1):
               break
# electromagnetic case
         elif (in1.ndim==3):
# select x-vx
            ps[:,:] = pvs[:,0,:nsxb[0]]
            psx[:,:] = numpy.transpose(ps)
            xmax = float(nx)
            ymax = fvs[nmv21,0,0]; ymin = -ymax
# display real space data
            gname = numpy.str.rstrip(sname[n-1])+" X-VX"
            graf2.dscalerl2(psx,gname,0.0,xmax,ymin,ymax,it,999,1,
                            nsxb[0],nmv21,ierr)
            if (ierr[0]==1):
               break
# select x-vy
            ps[:,:] = pvs[:,1,:nsxb[0]]
            psx[:,:] = numpy.transpose(ps)
            ymax = fvs[nmv21,1,0]; ymin = -ymax
# display real space data
            gname = numpy.str.rstrip(sname[n-1])+" X-VY"
            graf2.dscalerl2(psx,gname,0.0,xmax,ymin,ymax,it,999,1,
                            nsxb[0],nmv21,ierr)
            if (ierr[0]==1):
               break
# select x-vz
            ps[:,:] = pvs[:,2,:nsxb[0]]
            psx[:,:] = numpy.transpose(ps)
            ymax = fvs[nmv21,1,0]; ymin = -ymax
# display real space data
            gname = numpy.str.rstrip(sname[n-1])+" X-VZ"
            graf2.dscalerl2(psx,gname,0.0,xmax,ymin,ymax,it,999,1,
                            nsxb[0],nmv21,ierr)
            if (ierr[0]==1):
               break

      graf1.reset_graphs()
      cmfield1.closeff1(iuv)
      print()

if (__name__=="__main__"):
   main(-1)

