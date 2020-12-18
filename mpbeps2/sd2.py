#-----------------------------------------------------------------------
"""
High Level library for 2-1/2D Darwin OpenMP PIC code

functions defined:

bwrite_drestart33: write out basic restart file for darwin code
bread_drestart33: read in basic restart file for darwin code

written by Viktor K. Decyk and Joshua Kelly, UCLA
copyright 2016-2018, regents of the university of california
update: November 1, 2020
"""

import numpy

int_type = numpy.int32
float_type = numpy.float32

i3 = numpy.zeros((3),int_type)
a2 = numpy.zeros((2),float_type)

#-----------------------------------------------------------------------
def bwrite_drestart23(cus,wpm,q2m0,iur):
   """
   write out basic restart file for darwin code
   input:
   cus = smoothed transverse electric field
   wpm = normalized total plasma frequency squared
   q2m0 = shift constant in darwin iteration
   iur = restart file descriptors
   """
   ndim = numpy.size(cus,0); nxv = numpy.size(cus,1)
   nypmx = numpy.size(cus,2)
# write out shift constants for iteration
   a2[0] = wpm[0]; a2[1] = q2m0[0]
   a2.tofile(iur)
# write out darwin electric field
   i3[0] = ndim; i3[1] = nxv; i3[2] = nypmx
   i3.tofile(iur)
   if (ndim*nxv*nypmx > 0):
      cus.tofile(iur)

#-----------------------------------------------------------------------
def bread_drestart23(cus,wpm,q2m0,iur):
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
   nypmx = numpy.size(cus,2)
# read in shift constant for iteration
   a2[:] = numpy.fromfile(iur,float_type,2)
   wpm[0] = a2[0]; q2m0[0] = a2[1]
# read in darwin electric field field
   i3[:] = numpy.fromfile(iur,int_type,3)
   it = i3[0]; iu = i3[1]; iv = i3[2]
   if (it != ndim):
      print ("cus restart error, size(cus,0)=", it, ndim)
   elif (iu != nxv):
      print ("cus restart error, size(cus,1)=", iu, nxv)
   elif (iv != nypmx):
      print ("cus restart error, size(cus,2)=", iv, nypmx)
# read in field array
   il = it*iu*iv
   if (il > 0):
      cus[:,:,:] = numpy.fromfile(iur,float_type,il).reshape(it,iu,iv)
