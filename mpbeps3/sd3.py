#-----------------------------------------------------------------------
"""
High Level library for 3D Darwin OpenMP PIC code

functions defined:

bwrite_drestart13: write out basic restart file for darwin code
bread_drestart13: read in basic restart file for darwin code

written by Viktor K. Decyk and Joshua Kelly, UCLA
copyright 2016-2018, regents of the university of california
update: june 27, 2018
"""

import numpy

int_type = numpy.int32
float_type = numpy.float32

i4 = numpy.zeros((4),int_type)
a2 = numpy.zeros((2),float_type)

#-----------------------------------------------------------------------
def bwrite_drestart3(cus,wpm,q2m0,iur):
   """
   write out basic restart file for darwin code
   input:
   cus = smoothed transverse electric field
   wpm = normalized total plasma frequency squared
   q2m0 = shift constant in darwin iteration
   iur = restart file descriptors
   """
   ndim = numpy.size(cus,0); nxv = numpy.size(cus,1)
   nypmx = numpy.size(cus,2); nzpmx = numpy.size(cus,3)
# write out shift constants for iteration
   a2[0] = wpm[0]; a2[1] = q2m0[0]
   a2.tofile(iur)
# write out darwin electric field
   i4[0] = ndim; i4[1] = nxv; i4[2] = nypmx; i4[3] = nzpmx
   i4.tofile(iur)
   if (ndim*nxv*nypmx*nzpmx > 0):
      cus.tofile(iur)

#-----------------------------------------------------------------------
def bread_drestart3(cus,wpm,q2m0,iur):
   """
   read in basic restart file for darwin code
   input:
   iur = restart file descriptors
   output:
   cus = smoothed transverse electric field
   wpm = normalized total plasma frequency squared
   q2m0 = shift constant in darwin iteration
   """
   ndim = numpy.size(cus,0); nxv = numpy.size(cus,1)
   nypmx = numpy.size(cus,2); nzpmx = numpy.size(cus,3)
# read in shift constant for iteration
   a2[:] = numpy.fromfile(iur,float_type,2)
   wpm[0] = a2[0]; q2m0[0] = a2[1]
# read in darwin electric field field
   i4[:] = numpy.fromfile(iur,int_type,4)
   it = i4[0]; iu = i4[1]; iv = i4[2]; iw = i4[3]
   if (it != ndim):
      print "cus restart error, size(cus,0)=", it, ndim
   elif (iu != nxv):
      print "cus restart error, size(cus,1)=", iu, nxv
   elif (iv != nypmx):
      print "cus restart error, size(cus,2)=", iv, nypmx
   elif (iw != nzpmx):
      print "cus restart error, size(cus,3)=", iw, nzpmx
# read in field array
   il = it*iu*iv*iw
   if (il > 0):
      cus[:,:,:,:] = numpy.fromfile(iur,float_type,il).reshape(it,iu,iv,iw)
