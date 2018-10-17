#-----------------------------------------------------------------------
"""
High Level library for 2D Electrostatic OpenMP PIC code

functions defined:

open_restart2: open reset and restart files
bwrite_restart2: write out basic restart file for electrostatic code
bread_restart2a: read in basic restart file for electrostatic code
bread_restart2b: read in basic restart file for electrostatic code
bread_restart2c: read in basic restart file for electrostatic code
close_restart2: close reset and restart files

written by Viktor K. Decyk, UCLA
copyright 2016-2018, regents of the university of california
update: august 15, 2018
"""

import numpy

# sys.path.append('./mbeps2.source')
from libmpush2 import *

int_type = numpy.int32
float_type = numpy.float32

irc = numpy.zeros((1),int_type)

#-----------------------------------------------------------------------
def open_restart2(cdrun):
   """
   open restart files
   input:
   cdrun = idrun string
   output:
   iur = restart file descriptor
   iur0 = old restart file descriptor
   """
   global i1, i2
   iur = 0; iur0 = 0
# reset file
#  fname = "reset2"
#  iurr = open(fname,"wb+")
# restart file
   fname = "rstrt2." + cdrun
# start a new run from random numbers
   if (in2.nustrt==1):
      if (in2.ntr > 0):
         iur = open(fname,"wb+")
# continue a run which was interrupted
   elif (in2.nustrt==2):
      iur = open(fname,"rb+")
# start a new run with data from a previous run
   elif (in2.nustrt==0):
      if (in2.ntr > 0):
         iur = open(fname,"wb+")
      if (in2.idrun != in2.idrun0):
         cdrun0 = str(in2.idrun0)
         fname = "rstrt2." + cdrun0
         iur0 = open(fname,"rb+")
      else:
         print "restart warning: old, new idruns identical"
# allocate scratch numpy arrays
   i1 = numpy.zeros((1),int_type)
   i2 = numpy.zeros((2),int_type)
   return iur, iur0

#-----------------------------------------------------------------------
def bwrite_restart2(part,ppart,pparti,qi,kpic,kipic,iur,ntime,ntime0,
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
   global i1, i2
# write out basic restart file for electrostatic code
   idimp = numpy.size(part,0)
   iur.seek(0,0)
# write out current and initial time
   i2[0] = ntime; i2[1] = ntime0
   i2.tofile(iur)
# write out number of processors
   i1[0] = in2.nvp
   i1.tofile(iur)
# copy ordered particles to linear array: updates part, npp, irc
   mpush2.mpcopyout2(part,ppart,kpic,i1,irc)
   npp = i1[0]
# write out size of electron array
   i2[0] = idimp; i2[1] = npp; i2.tofile(iur)
# write out electrons: updates part
   if (i1[0] > 0):
      part[:,0:npp].tofile(iur)
   i1[0] = in2.movion; i1.tofile(iur)
   if (in2.movion > 0):
# copy ordered particles to linear array: updates part, nppi, irc
      mpush2.mpcopyout2(part,pparti,kipic,i1,irc)
      nppi = i1[0]
# write out size of ion array
      i2[0] = idimp; i2[1] = nppi; i2.tofile(iur)
# write out ions: updates part
      if (i1[0] > 0):
         part[:,0:nppi].tofile(iur)
# write out ion density, if ions are not moving
   else:
      nxv = numpy.size(qi,0); nypmx = numpy.size(qi,1)
      i2[0] = nxv; i2[1] = nypmx
      i2.tofile(iur)
      if (nxv*nypmx > 0):
         qi.tofile(iur)
# write out electric field parameter
   i1[0] = in2.emf; i1.tofile(iur)

#-----------------------------------------------------------------------
def bread_restart2a(part,kpic,iur,ntime,ntime0,npp,nppmx,noff,mx1,irc):
   """
   read in basic restart file for electrostatic code, first part
   input:
   iur = restart file descriptors
   noff = lowermost global gridpoint in particle partition
   mx1 = (system length in x direction - 1)/mx + 1, where
         mx = number of grids in sorting cell in x
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
   print "restarting from ntime, idrun0=", ntime[0], in2.idrun0
# read in number of processors
   i1[:] = numpy.fromfile(iur,int_type,1)
   it = i1[0]
   if (it != in2.nvp):
      print "restart error, nvp=", it, in2.nvp
# read in size of electron array
   i2[:] = numpy.fromfile(iur,int_type,2)
   ndimp = i2[0]; npp[0] = i2[1]
   if (ndimp != numpy.size(part,0)):
      print "restart error, idimp=", ndimp, numpy.size(part,0)
# read in electrons: updates part, npp, irc
   if (npp[0] > 0):
      it = npp[0]
      il = ndimp*it
      part[:,0:it] = numpy.fromfile(iur,float_type,il).reshape(ndimp,it)
# find number of electrons in each of mx, my tiles: updates kpic, nppmx
   minit2.mpdblkp2(part,kpic,npp,noff,nppmx,in2.mx,in2.my,mx1,irc)

#-----------------------------------------------------------------------
def bread_restart2b(part,ppart,kpic,kipic,iur,npp,nppi,nppmx,noff,nyp,
                    nx,mx1,irc):
   """
   read in basic restart file for electrostatic code, second part
   input:
   iur = restart file descriptors
   npp = number of electrons in partition
   noff = lowermost global gridpoint in particle partition
   nyp = number of primary (complete) gridpoints in particle partition
   nx = system length in x direction
   mx1 = (system length in x direction - 1)/mx + 1, where
         mx = number of grids in sorting cell in x
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
   mpush2.mpmovin2p(part,ppart,kpic,npp,noff,in2.mx,in2.my,mx1,irc)
# sanity check for electrons
   mpush2.mpcheck2(ppart,kpic,noff,nyp,nx,in2.mx,in2.my,mx1,irc)
# read in movion to determine if ions are moving
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in2.movion):
      print "movion restart error, movion = ", it, in2.movion
# ions are moving
   if (in2.movion > 0):
# read in size of ion array
      i2[:] = numpy.fromfile(iur,int_type,2)
      ndimp = i2[0]; nppi[0] = i2[1]
      if (ndimp != numpy.size(part,0)):
         print "ion restart error,idimp=",ndimp,numpy.size(part,0)
# read in ions: updates part, nppi, irc
      if (nppi[0] > 0):
         it = nppi[0]
         il = ndimp*it 
         part[:,0:it] = numpy.fromfile(iur,float_type,il).reshape(ndimp,it)
# find number of ions in each of mx, my tiles: updates kipic, nppmx
      minit2.mpdblkp2(part,kipic,nppi,noff,nppmx,in2.mx,in2.my,mx1,irc)

#-----------------------------------------------------------------------
def bread_restart2c(part,pparti,kipic,qi,iur,ntime,ntime0,nppi,noff,nyp,
                    nx,mx1,irc):
   """
   read in basic restart file for electrostatic code, third part
   input:
   part = particle array
   iur = restart file descriptors
   nppi = number of ions in partition
   noff = lowermost global gridpoint in particle partition
   nyp = number of primary (complete) gridpoints in particle partition
   nx = system length in x direction
   mx1 = (system length in x direction - 1)/mx + 1, where
         mx = number of grids in sorting cell in x
   output:
   pparti = tiled ion particle arrays
   kipic = number of ions in each tile
   qi = ion charge density with guard cells
   ntime/ntime0 = current/initial time step
   irc = maximum overflow, returned only if error occurs, when irc > 0
   """
   global i1, i2
   irc[0] = 0
# ions are moving
   if (in2.movion > 0):
# copy ordered ion data for OpenMP: updates pparti and kipic
      mpush2.mpmovin2p(part,pparti,kipic,nppi,noff,in2.mx,in2.my,mx1,
                       irc)
# sanity check for ions
      mpush2.mpcheck2(pparti,kipic,noff,nyp,nx,in2.mx,in2.my,mx1,irc)
# ions are not moving, read in ion density qi
   else:
      i2[:] = numpy.fromfile(iur,int_type,2)
      it = i2[0]; iu = i2[1]
      if (it != numpy.size(qi,0)):
         print "qi restart error, size(qi,0)=",it,numpy.size(qi,0)
      elif (iu != numpy.size(qi,1)):
         print "qi restart error, size(qi,1)=",iu,numpy.size(qi,1)
      il = it*iu
      if (il > 0):
         qi[:,:] = numpy.fromfile(iur,float_type,il).reshape(it,iu)
# read in electric field parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in2.emf):
      print "warning: emf values differ, emf=",it,in2.emf
   ntime0[0] = ntime0[0] + ntime[0]

#-----------------------------------------------------------------------
def close_restart2(iur,iur0):
   """
   read in basic restart file for electrostatic code
   input:
   iur, iur0 = restart, old restart file descriptors
   """
   global i1, i2
# iur, iur0 = restart, old restart file descriptors
#     close(unit=iurr)
   if (in2.nustrt==1):
      if (in2.ntr > 0):
         iur.close()
   elif (in2.nustrt==2):
      iur.close()
   elif (in2.nustrt==0):
      if (in2.ntr > 0):
         iur.close()
      if (in2.idrun != in2.idrun0):
         iur0.close()
   del i1, i2
