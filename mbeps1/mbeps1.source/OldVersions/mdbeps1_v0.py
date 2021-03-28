#-----------------------------------------------------------------------
# 1-2/2D Darwin OpenMP PIC code
# written by Viktor K. Decyk and Joshua Kelly, UCLA
import math
import numpy
from libmpush1 import *
from libmbpush1 import *
from libmdpush1 import *
from fomplib import *
from dtimer import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# indx = exponent which determines grid points in x direction:
# nx = 2**indx.
indx =   9
# npx = number of electrons distributed in x direction.
npx =  18432
# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
# qme = charge on electron, in units of e.
tend = 10.0; dt = 0.1; qme = -1.0
# vtx/vty = thermal velocity of electrons in x/y direction
# vx0/vy0 = drift velocity of electrons in x/y direction.
vtx = 1.0; vty = 1.0; vx0 = 0.0; vy0 = 0.0
# vtx/vz0 = thermal/drift velocity of electrons in z direction
vtz = 1.0; vz0 = 0.0
# ax = smoothed particle size in x direction
# ci = reciprocal of velocity of light.
ax = .912871; ci = 0.1
# idimp = number of particle coordinates = 4
# ipbc = particle boundary condition: 1 = periodic
# relativity = (no,yes) = (0,1) = relativity is used
idimp = 4; ipbc = 1; relativity = 0
# omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
omx = 0.4; omy = 0.0; omz = 0.0
# ndc = number of corrections in darwin iteration
ndc = 1
# wke/we = particle kinetic/electrostatic field energy
# wf/wm/wt = magnetic field/transverse electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wf = numpy.zeros((1),float_type)
wm = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
zero = 0.0
# mx = number of grids in x in sorting tiles
mx = 32
# xtras = fraction of extra particles needed for particle management
xtras = 0.2
# list = (true,false) = list of particles leaving tiles found in push
list = True
# declare scalars for standard code
wpmax = numpy.empty((1),float_type)
wpmin = numpy.empty((1),float_type)

# declare scalars for OpenMP code
nppmx = numpy.empty((1),int_type)
irc = numpy.zeros((1),int_type)

# declare and initialize timing data
tinit = 0.0; tloop = 0.0
itime = numpy.empty((4),numpy.int32)
ltime = numpy.empty((4),numpy.int32)
tdpost = numpy.zeros((1),float_type)
tguard = numpy.zeros((1),float_type)
tfft = numpy.zeros((1),float_type)
tfield = numpy.zeros((1),float_type)
tdjpost = numpy.zeros((1),float_type)
tdcjpost = numpy.zeros((1),float_type)
tpush = numpy.zeros((1),float_type)
tsort = numpy.zeros((1),float_type)
dtime = numpy.empty((1),double_type)

# start timing initialization
dtimer(dtime,itime,-1)

# nvp = number of shared memory nodes (0=default)
nvp = 0
#nvp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
omplib.init_omp(nvp)

# initialize scalars for standard code
# np = total number of particles in simulation
# nx = number of grid points in x direction
np = npx; nx = int(math.pow(2,indx)); nxh = int(nx/2)
nxe = nx + 2; nxeh = int(nxe/2)
# mx1 = number of tiles in x direction
mx1 = int((nx - 1)/mx + 1)
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
qbme = qme
affp = float(nx)/float(np)

# check for unimplemented features
if (list & (ipbc != 1)):
   print "ipbc != 1 and list = True not yet supported"
   print "list reset to False"
   list = False

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe),float_type,'F')
# fxe = smoothed longitudinal electric field with guard cells
fxe = numpy.empty((nxe),float_type,'F')
# cue = electron current density with guard cells
cue = numpy.empty((2,nxe),float_type,'F')
# dcu = acceleration density with guard cells
dcu = numpy.empty((2,nxe),float_type,'F')
# cus = transverse electric field
cus = numpy.empty((2,nxe),float_type,'F')
# amu = momentum flux with guard cells
amu = numpy.empty((2,nxe),float_type,'F')
# exyze = smoothed electric field with guard cells
exyze = numpy.empty((3,nxe),float_type,'F')
# byze = smoothed magnetic field with guard cells
byze = numpy.empty((2,nxe),float_type,'F')
# eyz = transverse electric field in fourier space
eyz = numpy.empty((2,nxeh),complex_type,'F')
# byz = transverse magnetic field in fourier space
byz = numpy.empty((2,nxeh),complex_type,'F')
# ffc, ffe = form factor arrays for poisson solvers
ffc = numpy.empty((nxh),complex_type,'F')
ffe = numpy.empty((nxh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxh),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxh),complex_type,'F')
# arrays required for darwin initial fields
# kpic = number of particles in each tile
kpic = numpy.empty((mx1),int_type,'F')


# prepare fft tables
mfft1.mfft1_init(mixup,sct,indx)
# calculate form factor: ffc
mfield1.mpois1_init(ffc,ax,affp,nx)
# initialize electrons
minit1.mdistr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,nx,ipbc)

# find number of particles in each of mx, tiles: updates kpic, nppmx
minit1.mdblkp2(part,kpic,nppmx,mx,irc)

# allocate vector particle data
nppmx0 = int((1.0 + xtras)*nppmx)
ntmax = int(xtras*nppmx)
npbmx = int(xtras*nppmx)
# ppart = tiled particle array
ppart = numpy.empty((idimp,nppmx0,mx1),float_type,'F')
# ppbuff = buffer array for reordering tiled particle array
ppbuff = numpy.empty((idimp,npbmx,mx1),float_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((2,mx1),int_type,'F')
# ihole = location/destination of each particle departing tile
ihole = numpy.empty((2,ntmax+1,mx1),int_type,'F')

# copy ordered particle data for OpenMP: updates ppart and kpic
mpush1.mpmovin1(part,ppart,kpic,mx,irc)

# sanity check
mpush1.mcheck1(ppart,kpic,nx,mx,irc)

# find maximum and minimum initial electron density
qe.fill(0.0)
mpush1.mpost1(ppart,qe,kpic,qme,tdpost,mx)
mgard1.maguard1(qe,tguard,nx)
mdpush1.mfwpminx1(qe,qbme,wpmax,wpmin,nx)
wpm = 0.5*(wpmax[0] + wpmin[0])*affp
# accelerate convergence: update wpm
if (wpm <= 10.0):
   wpm = 0.75*wpm
print "wpm=",wpm
q2m0 = wpm/affp
# calculate form factor: ffe
mfield1.mepois1_init(ffe,ax,affp,wpm,ci,nx)

# initialize electric fields
cus.fill(0.0); fxe.fill(0.0)

# initialization time
dtimer(dtime,itime,1)
tinit = tinit + float(dtime)
# start timing loop
dtimer(dtime,ltime,-1)

print "program mdbeps1"

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit current with OpenMP: updates cue
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   dtimer(dtime,itime,1)
   tdjpost[0] = tdjpost[0] + float(dtime)
   mcurd1.wmdjpost1(ppart,cue,kpic,ncl,ihole,qme,zero,ci,tdjpost,nx,mx,
                    ipbc,relativity,False,irc)

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   dtimer(dtime,itime,1)
   tdpost[0] = tdpost[0] + float(dtime)
   mpush1.mpost1(ppart,qe,kpic,qme,tdpost,mx)

# add guard cells with standard procedure: updates qe, cue
   mgard1.maguard1(qe,tguard,nx)
   mgard1.macguard1(cue,tguard,nx)

# transform charge to fourier space with standard procedure:
# updates qe
   isign = -1
   mfft1.mfft1r(qe,isign,mixup,sct,tfft,indx)

# calculate longitudinal force/charge in fourier space with standard
# procedure: updates fxe, we
   mfield1.mpois1(qe,fxe,ffc,we,tfield,nx)

# transform longitudinal electric force to real space with standard
# procedure: updates fxe
   isign = 1
   mfft1.mfft1r(fxe,isign,mixup,sct,tfft,indx)

# transform current to fourier space with standard procedure:
# updates cue
   isign = -1
   mfft1.mfft1rn(cue,isign,mixup,sct,tfft,indx)

# calculate magnetic field in fourier space with standard procedure:
# updates byze, wm
   mfield1.mbbpois1(cue,byze,ffc,ci,wm,tfield,nx)

# transform magnetic force to real space with standard procedure:
# updates byze
   isign = 1
   mfft1.mfft1rn(byze,isign,mixup,sct,tfft,indx)

# add constant to magnetic field with standard procedure: updates byze
   mfield1.mbaddext1(byze,tfield,omy,omz,nx)

# copy guard cells with standard procedure: updates fxe, byze
   mgard1.mdguard1(fxe,tguard,nx)
   mgard1.mcguard1(byze,tguard,nx)

# add longitudinal and old transverse electric fields with standard
# procedure: updates exyze
   mfield1.maddvrfield1(exyze,cus,fxe,tfield)

# deposit electron acceleration density and momentum flux with OpenMP:
# updates dcu, amu
   dtimer(dtime,itime,-1)
   dcu.fill(0.0); amu.fill(0.0)
   dtimer(dtime,itime,1)
   tdcjpost[0] = tdcjpost[0] + float(dtime)
   mdpush1.wmgdjpost1(ppart,exyze,byze,dcu,amu,kpic,omx,qme,qbme,dt,ci,
                      tdcjpost,nx,mx,relativity)
# add old scaled electric field with standard procedure: updates dcu
   mdpush1.mascfguard1(dcu,cus,q2m0,tdcjpost,nx)

# add guard cells with standard procedure: updates dcu, amu
   mgard1.macguard1(dcu,tguard,nx)
   mgard1.macguard1(amu,tguard,nx)

# transform acceleration density and momentum flux to fourier space
# with standard procedure: updates dcu, amu
   isign = -1
   mfft1.mfft1rn(dcu,isign,mixup,sct,tfft,indx)
   mfft1.mfft1rn(amu,isign,mixup,sct,tfft,indx)

# take transverse part of time derivative of current with standard
# procedure: updates dcu
   mfield1.madcuperp1(dcu,amu,tfield,nx)

# calculate transverse electric field with standard procedure:
# updates cus, wf
   mfield1.mepois1(dcu,cus,ffe,affp,ci,wf,tfield,nx)

# transform transverse electric field to real space with standard
# procedure: updates cus
   isign = 1
   mfft1.mfft1rn(cus,isign,mixup,sct,tfft,indx)

# copy guard cells with standard procedure: updates cus
   mgard1.mcguard1(cus,tguard,nx)

# add longitudinal and transverse electric fields with standard
# procedure: exyze = cus + fxe, updates exyze
# cus needs to be retained for next time step
   mfield1.maddvrfield1(exyze,cus,fxe,tfield)

# inner iteration loop
   for k in xrange(0,ndc):

# deposit electron current and acceleration density and momentum flux
# with OpenMP: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      cue.fill(0.0); dcu.fill(0.0); amu.fill(0.0)
      dtimer(dtime,itime,1)
      tdcjpost[0] = tdcjpost[0] + float(dtime)
      mdpush1.wmgdcjpost1(ppart,exyze,byze,cue,dcu,amu,kpic,omx,qme,
                          qbme,dt,ci,tdcjpost,nx,mx,relativity)
# add scaled electric field with standard procedure: updates dcu
      mdpush1.mascfguard1(dcu,cus,q2m0,tdcjpost,nx)

# add guard cells for current, acceleration density, and momentum flux
# with standard procedure: updates cue, dcu, amu
      mgard1.macguard1(cue,tguard,nx)
      mgard1.macguard1(dcu,tguard,nx)
      mgard1.macguard1(amu,tguard,nx)

# transform current to fourier space with standard procedure:
# update cue
      isign = -1
      mfft1.mfft1rn(cue,isign,mixup,sct,tfft,indx)

# calculate magnetic field in fourier space with standard procedure:
# updates byze, wm
      mfield1.mbbpois1(cue,byze,ffc,ci,wm,tfield,nx)

# transform magnetic force to real space with standard procedure:
# updates byze
      isign = 1
      mfft1.mfft1rn(byze,isign,mixup,sct,tfft,indx)

# add constant to magnetic field with standard procedure: updates bzye
      mfield1.mbaddext1(byze,tfield,omy,omz,nx)

# transform acceleration density and momentum flux to fourier space
# with standard procedure: updates dcu, amu
      isign = -1
      mfft1.mfft1rn(dcu,isign,mixup,sct,tfft,indx)
      mfft1.mfft1rn(amu,isign,mixup,sct,tfft,indx)

# take transverse part of time derivative of current with standard
# procedure: updates dcu
      mfield1.madcuperp1(dcu,amu,tfield,nx)

# calculate transverse electric field with standard procedure:
# updates cus, wf
      mfield1.mepois1(dcu,cus,ffe,affp,ci,wf,tfield,nx)

# transform transverse electric field to real space with standard
# procedure: updates cus
      isign = 1
      mfft1.mfft1rn(cus,isign,mixup,sct,tfft,indx)

# copy guard cells with standard procedure: updates byze, cus
      mgard1.mcguard1(byze,tguard,nx)
      mgard1.mcguard1(cus,tguard,nx)

# add longitudinal and transverse electric fields with standard
# procedure: exyze = cus + fxe, updates exyze
# cus needs to be retained for next time step
      mfield1.maddvrfield1(exyze,cus,fxe,tfield)
   pass

# push particles with OpenMP:
# updates ppart and wke, and possibly ncl, ihole, irc
   wke[0] = 0.0
   mbpush1.wmbpush1(ppart,exyze,byze,kpic,ncl,ihole,omx,qbme,dt,dt,ci,
                    wke,tpush,nx,mx,ipbc,relativity,list,irc)

# reorder particles by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)

   if (ntime==0):
      wt = we + wf + wm
      print "Initial Total Field, Kinetic and Total Energies:"
      print "%14.7e %14.7e %14.7e" % (wt, wke, wke + wt)
      print "Initial Electrostatic, Transverse Electric and Magnetic " \
      "Field Energies:"
      print "%14.7e %14.7e %14.7e" % (we, wf, wm)

# loop time
   dtimer(dtime,ltime,1)
   tloop = tloop + float(dtime)

ntime = ntime + 1

# * * * end main iteration loop * * *

print "ntime, relativity, ndc = ", ntime, relativity, ndc
wt[0] = we + wm
print "Final Total Field, Kinetic and Total Energies:"
print "%14.7e %14.7e %14.7e" % (wt[0], wke, wke + wt[0])
print "Final Electrostatic, Transverse Electric and Magnetic Field " \
"Energies:"
print "%14.7e %14.7e %14.7e" % (we, wf, wm)

print ""
print "initialization time = ", tinit
print "deposit time = ", tdpost[0]
print "current deposit time = ", tdjpost[0]
print "current derivative deposit time = ", tdcjpost[0]
tdpost[0] = tdpost[0] + tdjpost[0] + tdcjpost[0]
print "total deposit time = ", tdpost[0]
print "guard time = ", tguard[0]
print "solver time = ", tfield[0]
print "fft time = ", tfft[0]
print "push time = ", tpush[0]
print "sort time = ", tsort[0]
tfield[0] = tfield[0] + tguard[0] + tfft[0]
print "total solver time = ", tfield[0]
time = tdpost[0] + tpush[0] + tsort[0]
print "total particle time = ", time
wt[0] = time + tfield[0]
tloop = tloop - wt[0]
print "total and additional time = ", wt[0], tloop
print ""

wt[0] = 1.0e+09/(float(nloop)*float(np))
print "Push Time (nsec) = ", tpush[0]*wt[0]
print "Deposit Time (nsec) = ", tdpost[0]*wt[0]
print "Sort Time (nsec) = ", tsort[0]*wt[0]
print "Total Particle Time (nsec) = ", time*wt[0]
print ""
