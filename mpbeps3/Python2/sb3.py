#-----------------------------------------------------------------------
"""
High Level library for 3D Electromagnetic OpenMP PIC code

functions defined:
                          
bwrite_restart3: write out basic restart file for electromagnetic code
bread_restart3:  read in basic restart file for electromagnetic code

written by Viktor K. Decyk, UCLA
copyright 2016-2018, regents of the university of california
update: june 13, 2018
"""

import numpy

int_type = numpy.int32
complex_type = numpy.complex64

i4 = numpy.zeros((4),int_type)

#-----------------------------------------------------------------------
def bwrite_brestart3(exyz,bxyz,iur):
   """
   write out basic restart file for electromagnetic code
   input:
   exyz/bxyz = transverse electric/magnetic field in fourier space
   iur = restart file descriptors
   """
   ndim = numpy.size(exyz,0); nzv = numpy.size(exyz,1)
   kxypd = numpy.size(exyz,2); kyzpd = numpy.size(exyz,3)
# write out electromagnetic fields
   i4[0] = ndim; i4[1] = nzv; i4[2] = kxypd; i4[3] = kyzpd
   i4.tofile(iur)
   if (ndim*nzv*kxypd*kyzpd > 0):
      exyz.tofile(iur)
      bxyz.tofile(iur)

#-----------------------------------------------------------------------
def bread_restart3(exyz,bxyz,iur):
   """
   read in basic restart file for electromagnetic code
   input:
   iur = restart file descriptors
   output:
   exyz/bxyz = transverse electric/magnetic field in fourier space
   """
   ndim = numpy.size(exyz,0); nzv = numpy.size(exyz,1)
   kxypd = numpy.size(exyz,2); kyzpd = numpy.size(exyz,3)
# read in electromagnetic fields
   i4[:] = numpy.fromfile(iur,int_type,4)
   it = i4[0]; iu = i4[1]; iv = i4[2]; iw = i4[3]
   if (it != ndim):
      print "exyz/bxyz restart error, size(exyz,0)=", it, ndim
   elif (iu != nzv):
      print "exyz/bxyz restart error, size(exyz,1)=", iu, nzv
   elif (iv != kxypd):
      print "exyz/bxyz restart error, size(exyz,2)=", iv, kxypd
   elif (iw != kyzpd):
      print "exyz/bxyz restart error, size(exyz,3)=", iw, kyzpd
# read in field arrays
   il = it*iu*iv*iw
   if (il > 0):
      exyz[:,:,:,:] = numpy.fromfile(iur,complex_type,il).reshape(it,iu,iv,iw)
      bxyz[:,:,:,:] = numpy.fromfile(iur,complex_type,il).reshape(it,iu,iv,iw)
