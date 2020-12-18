#-----------------------------------------------------------------------
"""
High Level library for 2-1/2D Electromagnetic OpenMP PIC code

functions defined:
                          
bwrite_restart23: write out basic restart file for electromagnetic code
bread_restart23:  read in basic restart file for electromagnetic code

written by Viktor K. Decyk, UCLA
copyright 2016-2018, regents of the university of california
update: august 15, 2018
"""

import numpy

int_type = numpy.int32
complex_type = numpy.complex64

i3 = numpy.zeros((3),int_type)

#-----------------------------------------------------------------------
def bwrite_restart23(exyz,bxyz,iur):
   """
   write out basic restart file for electromagnetic code
   input:
   exyz/bxyz = transverse electric/magnetic field in fourier space
   iur = restart file descriptors
   """
   ndim = numpy.size(exyz,0); nyv = numpy.size(exyz,1)
   kxpd = numpy.size(exyz,2)
# write out electromagnetic fields
   i3[0] = ndim; i3[1] = nyv; i3[2] = kxpd
   i3.tofile(iur)
   if (ndim*nyv*kxpd > 0):
      exyz.tofile(iur)
      bxyz.tofile(iur)

#-----------------------------------------------------------------------
def bread_restart23(exyz,bxyz,iur):
   """
   read in basic restart file for electromagnetic code
   input:
   iur = restart file descriptors
   output:
   exyz/bxyz = transverse electric/magnetic field in fourier space
   """
   ndim = numpy.size(exyz,0); nyv = numpy.size(exyz,1)
   kxpd = numpy.size(exyz,2)
# read in electromagnetic fields
   i3[:] = numpy.fromfile(iur,int_type,3)
   it = i3[0]; iu = i3[1]; iv = i3[2]
   if (it != ndim):
      print "exyz/bxyz restart error, size(exyz,0)=", it, ndim
   elif (iu != nyv):
      print "exyz/bxyz restart error, size(exyz,1)=", iu, nyv
   elif (iv != kxpd):
      print "exyz/bxyz restart error, size(exyz,2)=", iv, kxpd
# read in field arrays
   il = it*iu*iv
   if (il > 0):
      exyz[:,:,:] = numpy.fromfile(iur,complex_type,il).reshape(it,iu,iv)
      bxyz[:,:,:] = numpy.fromfile(iur,complex_type,il).reshape(it,iu,iv)
