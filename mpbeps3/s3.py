#-----------------------------------------------------------------------
from __future__ import print_function
"""
High Level library for 3D Electrostatic OpenMP PIC code

functions defined:

open_restart3: open reset and restart files
bwrite_restart3: write out basic restart file for electrostatic code
bread_restart3a: read in basic restart file for electrostatic code
bread_restart3b: read in basic restart file for electrostatic code
bread_restart3c: read in basic restart file for electrostatic code
close_restart3: close reset and restart files

written by Viktor K. Decyk, UCLA
copyright 2016-2018, regents of the university of california
update: November 3, 2020
"""

import sys
import numpy

# sys.path.append('./mbeps3.source')
from libmpush3 import *

int_type = numpy.int32
float_type = numpy.float32

irc = numpy.zeros((1),int_type)

#-----------------------------------------------------------------------
def open_restart3(cdrun):
   """
   open restart files
   input:
   cdrun = idrun string
   output:
   iur = restart file descriptor
   iur0 = old restart file descriptor
   """
   global i1, i2, i3
   iur = 0; iur0 = 0
# reset file
#  fname = "reset3"
#  iurr = open(fname,"wb+")
# restart file
   fname = "rstrt3." + cdrun
# start a new run from random numbers
   if (in3.nustrt==1):
      if (in3.ntr > 0):
         iur = open(fname,"wb+")
# continue a run which was interrupted
   elif (in3.nustrt==2):
      iur = open(fname,"rb+")
# start a new run with data from a previous run
   elif (in3.nustrt==0):
      if (in3.ntr > 0):
         iur = open(fname,"wb+")
      if (in3.idrun != in3.idrun0):
         cdrun0 = str(in3.idrun0)
         fname = "rstrt3." + cdrun0
         iur0 = open(fname,"rb+")
      else:
         print ("restart warning: old, new idruns identical")
# allocate scratch numpy arrays
   i1 = numpy.zeros((1),int_type)
   i2 = numpy.zeros((2),int_type)
   i3 = numpy.zeros((3),int_type)
   return iur, iur0

#-----------------------------------------------------------------------
def bwrite_restart3(part,ppart,pparti,qi,kpic,kipic,iur,ntime,ntime0,
                    irc):
   """
   write out basic restart file for electrostatic code
   input:
   ppart/pparti = tiled electron/ion particle arrays
   qi = ion charge density with guard cells
   kpic/kipic = number of electrons/ions in each tile
   iur = restart file descriptors
   ntime/ntime0 = current/initial time step
   output:
   part = particle array
   irc = maximum overflow, returned only if error occurs, when irc > 0
   """
   global i1, i2, i3
# write out basic restart file for electrostatic code
   idimp = numpy.size(part,0)
   iur.seek(0,0)
# write out current and initial time
   i2[0] = ntime; i2[1] = ntime0
   i2.tofile(iur)
# write out number of processors
   i2[0] = in3.nvpy; i2[1] = in3.nvpz
   i2.tofile(iur)
# copy ordered particles to linear array: updates part, npp, irc
   mpush3.mpcopyout3(part,ppart,kpic,i1,irc)
   npp = i1[0]
# write out size of electron array
   i2[0] = idimp; i2[1] = npp; i2.tofile(iur)
# write out electrons: updates part
   if (i1[0] > 0):
      part[:,0:npp].tofile(iur)
   i1[0] = in3.movion; i1.tofile(iur)
   if (in3.movion==1):
# copy ordered particles to linear array: updates part, nppi, irc
      mpush3.mpcopyout3(part,pparti,kipic,i1,irc)
      nppi = i1[0]
# write out size of ion array
      i2[0] = idimp; i2[1] = nppi; i2.tofile(iur)
# write out ions: updates part
      if (i1[0] > 0):
         part[:,0:nppi].tofile(iur)
# write out ion density, if ions are not moving
   else:
      nxv = numpy.size(qi,0)
      nypmx = numpy.size(qi,1); nzpmx = numpy.size(qi,2)
      i3[0] = nxv; i3[1] = nypmx; i3[2] = nzpmx
      i3.tofile(iur)
      if (nxv*nypmx*nzpmx > 0):
         qi.tofile(iur)
# write out electric field parameter
   i1[0] = in3.emf; i1.tofile(iur)

#-----------------------------------------------------------------------
def bread_restart3a(part,kpic,noff,iur,ntime,ntime0,npp,nppmx,mx1,myp1,
                    irc):
   """
   read in basic restart file for electrostatic code, first part
   input:
   noff[0] = lowermost global gridpoint in y in particle partition
   noff[1] = backmost global gridpoint in z in particle partition
   iur = restart file descriptors
   mx1 = (system length in x direction - 1)/mx + 1,
   myp1 = (partition length in y direction - 1)/my + 1, where
           mx/my = number of grids in sorting cell in x and y
   output:
   part = particle array
   kpic = number of electrons in each tile
   ntime/ntime0 = current/initial time step
   npp = number of particles in partition
   nppmx = maximum number of particles in tile
   irc = maximum overflow, returned only if error occurs, when irc > 0
   """
   global i1, i2
   irc[0] = 0
# read in current and initial time
   iur.seek(0,0)
   i2[:] = numpy.fromfile(iur,int_type,2)
   ntime[0] = i2[0]; ntime0[0] = i2[1]
   print ("restarting from ntime, idrun0=", ntime[0], in3.idrun0)
# read in number of processors
   i2[:] = numpy.fromfile(iur,int_type,2)
   it = i2[0]; iu = i2[1]
   if (it != in3.nvpy):
      print ("restart error, nvpy, nvpz=", it, in3.nvpy)
   elif (iu != in3.nvpz):
      print ("restart error, nvpy, nvpz=", iu, in3.nvpz)
# read in size of electron array
   i2[:] = numpy.fromfile(iur,int_type,2)
   ndimp = i2[0]; npp[0] = i2[1]
   if (ndimp != numpy.size(part,0)):
      print ("restart error, idimp=", ndimp, numpy.size(part,0))
# read in electrons: updates part, npp, irc
   if (npp[0] > 0):
      it = npp[0]
      il = ndimp*it
      part[:,0:it] = numpy.fromfile(iur,float_type,il).reshape(ndimp,it)
# find number of electrons in each of mx, my tiles: updates kpic, nppmx
   minit3.mpdblkp3(part,kpic,npp,noff,nppmx,in3.mx,in3.my,in3.mz,mx1,
                   myp1,irc)

#-----------------------------------------------------------------------
def bread_restart3b(part,ppart,kpic,kipic,noff,nyzp,iur,npp,nppi,nppmx,
                    nx,mx1,myp1,irc):
   """
   read in basic restart file for electrostatic code, second part
   input:
   noff[0] = lowermost global gridpoint in y in particle partition
   noff[1] = backmost global gridpoint in z in particle partition
   nyzp[0:1] = number of primary (complete) gridpoints in y/z
   iur = restart file descriptors
   npp = number of electrons in partition
   nx = system length in x direction
   mx1 = (system length in x direction - 1)/mx + 1,
   myp1 = (partition length in y direction - 1)/my + 1, where
          mx/my = number of grids in sorting cell in x and y
   output:
   part = particle array
   ppart = tiled electron particle arrays
   kpic/kipic = number of electrons/ions in each tile
   nppi = number of ions in partition
   nppmx = maximum number of particles in tile
   irc = maximum overflow, returned only if error occurs, when irc > 0
   """
   global i1, i2
   irc[0] = 0
# copy ordered electron data for OpenMP: updates ppart and kpic
   mpush3.mpmovin3p(part,ppart,kpic,npp,noff,in3.mx,in3.my,in3.mz,mx1,
                    myp1,irc)
# sanity check for electrons
   mpush3.mpcheck3(ppart,kpic,noff,nyzp,nx,in3.mx,in3.my,in3.mz,mx1,
                   myp1,irc)
# read in movion to determine if ions are moving
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in3.movion):
      print ("movion restart error, movion = ", it, in3.movion)
# ions are moving
   if (in3.movion==1):
# read in size of ion array
      i2[:] = numpy.fromfile(iur,int_type,2)
      ndimp = i2[0]; nppi[0] = i2[1]
      if (ndimp != numpy.size(part,0)):
         print ("ion restart error,idimp=",ndimp,numpy.size(part,0))
# read in ions: updates part, nppi, irc
      if (nppi[0] > 0):
         it = nppi[0]
         il = ndimp*it 
         part[:,0:it] = numpy.fromfile(iur,float_type,il).reshape(ndimp,it)
# find number of ions in each of mx, my tiles: updates kipic, nppmx
      minit3.mpdblkp3(part,kipic,nppi,noff,nppmx,in3.mx,in3.my,in3.mz,
                      mx1,myp1,irc)

#-----------------------------------------------------------------------
def bread_restart3c(part,pparti,kipic,noff,nyzp,qi,iur,ntime,ntime0,
                    nppi,nx,mx1,myp1,irc):
   """
   read in basic restart file for electrostatic code, third part
   input:
   part = particle array
   noff[0] = lowermost global gridpoint in y in particle partition
   noff[1] = backmost global gridpoint in z in particle partition
   nyzp[0:1] = number of primary (complete) gridpoints in y/z
   iur = restart file descriptors
   nppi = number of ions in partition
   nx = system length in x direction
   mx1 = (system length in x direction - 1)/mx + 1,
   myp1 = (partition length in y direction - 1)/my + 1, where
          mx/my = number of grids in sorting cell in x and y
   output:
   pparti = tiled ion particle arrays
   kipic = number of ions in each tile
   qi = ion charge density with guard cells
   ntime/ntime0 = current/initial time step
   irc = maximum overflow, returned only if error occurs, when irc > 0
   """
   global i1, i2, i3
   irc[0] = 0
# ions are moving
   if (in3.movion==1):
# copy ordered ion data for OpenMP: updates pparti and kipic
      mpush3.mpmovin3p(part,pparti,kipic,nppi,noff,in3.mx,in3.my,in3.mz,
                       mx1,myp1,irc)
# sanity check for ions
      mpush3.mpcheck3(pparti,kipic,noff,nyzp,nx,in3.mx,in3.my,in3.mz,
                      mx1,myp1,irc)
# ions are not moving, read in ion density qi
   else:
      i3[:] = numpy.fromfile(iur,int_type,3)
      it = i3[0]; iu = i3[1]; iv = i3[2]
      if (it != numpy.size(qi,0)):
         print ("qi restart error, size(qi,0)=",it,numpy.size(qi,0))
      elif (iu != numpy.size(qi,1)):
         print ("qi restart error, size(qi,1)=",iu,numpy.size(qi,1))
      elif (iv != numpy.size(qi,2)):
         print ("qi restart error, size(qi,2)=",iv,numpy.size(qi,2))
      il = it*iu*iv
      if (il > 0):
         qi[:,:,:] = numpy.fromfile(iur,float_type,il).reshape(it,iu,iv)
# read in electric field parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in3.emf):
      print ("warning: emf values differ, emf=",it,in3.emf)
   ntime0[0] = ntime0[0] + ntime[0]

#-----------------------------------------------------------------------
def close_restart3(iur,iur0):
   """
   read in basic restart file for electrostatic code
   input:
   iur, iur0 = restart, old restart file descriptors
   """
   global i1, i2, i3
# iur, iur0 = restart, old restart file descriptors
#     close(unit=iurr)
   if (in3.nustrt==1):
      if (in3.ntr > 0):
         iur.close()
   elif (in3.nustrt==2):
      iur.close()
   elif (in3.nustrt==0):
      if (in3.ntr > 0):
         iur.close()
      if (in3.idrun != in3.idrun0):
         iur0.close()
   del i1, i2, i3
