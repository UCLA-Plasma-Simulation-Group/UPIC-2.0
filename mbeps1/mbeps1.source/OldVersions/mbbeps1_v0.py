#-----------------------------------------------------------------------
# 1-2/2D Electromagnetic OpenMP PIC code
# written by Viktor K. Decyk and Joshua Kelly, UCLA
import math
import numpy
from libmpush1 import *
from libmbpush1 import *
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
tend = 10.0; dt = 0.05; qme = -1.0
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
idimp = 4; ipbc = 1; relativity = 1
# omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
omx = 0.0; omy = 0.0; omz = 0.0
# wke/we = particle kinetic/electrostatic field energy
# wf/wm/wt = magnetic field/transverse electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wf = numpy.zeros((1),float_type)
wm = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
# mx = number of grids in x in sorting tiles
mx = 32
# xtras = fraction of extra particles needed for particle management
xtras = 0.2
# list = (true,false) = list of particles leaving tiles found in push
list = True

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
dth = 0.0

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
# fxyze = smoothed electric field with guard cells
fxyze = numpy.empty((3,nxe),float_type,'F')
# byze = smoothed magnetic field with guard cells
byze = numpy.empty((2,nxe),float_type,'F')
# eyz = transverse electric field in fourier space
eyz = numpy.empty((2,nxeh),complex_type,'F')
# byz = transverse magnetic field in fourier space
byz = numpy.empty((2,nxeh),complex_type,'F')
# ffc = form factor array for poisson solver
ffc = numpy.empty((nxh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxh),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxh),complex_type,'F')
# kpic = number of particles in each tile
kpic = numpy.empty((mx1),int_type,'F')

# prepare fft tables
mfft1.mfft1_init(mixup,sct,indx)
# calculate form factors
mfield1.mpois1_init(ffc,ax,affp,nx)
# initialize electrons
minit1.mdistr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,nx,ipbc)

# initialize transverse electromagnetic fields
eyz.fill(numpy.complex(0.0,0.0))
byz.fill(numpy.complex(0.0,0.0))
# set magnitude of external transverse magnetic field
omt = numpy.sqrt(omy*omy + omz*omz)

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

if (dt > 0.64*ci):
   print "Warning: Courant condition may be exceeded!"

# initialization time
dtimer(dtime,itime,1)
tinit = tinit + float(dtime)
# start timing loop
dtimer(dtime,ltime,-1)

print "program mbbeps1"

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit current with OpenMP:
# updates ppart and cue, and possibly ncl, ihole, irc
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   dtimer(dtime,itime,1)
   tdjpost[0] = tdjpost[0] + float(dtime)
   mcurd1.wmdjpost1(ppart,cue,kpic,ncl,ihole,qme,dth,ci,tdjpost,nx,mx,
                    ipbc,relativity,list,irc)

# reorder particles by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   dtimer(dtime,itime,1)
   tdpost[0] = tdpost[0] + float(dtime)
   mpush1.mpost1(ppart,qe,kpic,qme,tdpost,mx)

# add guard cells with standard procedure: updates cue, qe
   mgard1.macguard1(cue,tguard,nx)
   mgard1.maguard1(qe,tguard,nx)

# transform charge to fourier space with standard procedure:
# updates qe
   isign = -1
   mfft1.mfft1r(qe,isign,mixup,sct,tfft,indx)

# transform current to fourier space with standard procedure:
# updates cue
   isign = -1
   mfft1.mfft1rn(cue,isign,mixup,sct,tfft,indx)

# calculate electromagnetic fields in fourier space with standard
# procedure: updates eyz, byz
   if (ntime==0):
# initialize electromagnetic fields from darwin fields
# calculate initial darwin magnetic field
      mfield1.mibpois1(cue,byz,ffc,ci,wm,tfield,nx)
      wf[0] = 0.0
# calculate initial darwin electric field
      amu = numpy.empty((2,nxe),float_type,'F')
      dcu = numpy.empty((2,nxe),float_type,'F')
      amu.fill(0.0)
      mcurd1.wmgmjpost1(ppart,amu,kpic,qme,ci,tdjpost,mx,relativity)
      mgard1.macguard1(amu,tguard,nx)
      isign = -1
      mfft1.mfft1rn(amu,isign,mixup,sct,tfft,indx)
      mfield1.mdcuperp1(dcu,amu,tfield,nx)
      numpy.delete(amu,amu)
      mfield1.metfield1(dcu,eyz,ffc,ci,wf,tfield,nx)
      numpy.delete(dcu,dcu)
      dth = 0.5*dt
   else:
      mfield1.mmaxwel1(eyz,byz,cue,ffc,ci,dt,wf,wm,tfield,nx)

# calculate longitudinal force/charge in fourier space with standard
# procedure: updates fxe, we
   mfield1.mpois1(qe,fxe,ffc,we,tfield,nx)

# add longitudinal and transverse electric fields with standard
# procedure: updates fxyze
   mfield1.memfield1(fxyze,fxe,eyz,ffc,tfield,nx)
# copy magnetic field with standard procedure: updates byze
   mfield1.mbmfield1(byze,byz,ffc,tfield,nx)

# transform electric force to real space with standard procedure:
# updates fxyze
   isign = 1
   mfft1.mfft1rn(fxyze,isign,mixup,sct,tfft,indx)

# transform magnetic force to real space with standard procedure:
# updates byze
   isign = 1
   mfft1.mfft1rn(byze,isign,mixup,sct,tfft,indx)

# add constant to magnetic field with OpenMP: updates bxyze
   if (omt > 0.0):
      mfield1.mbaddext1(byze,tfield,omy,omz,nx)

# copy guard cells with standard procedure: updates fxyze, byze
   mgard1.mcguard1(fxyze,tguard,nx)
   mgard1.mcguard1(byze,tguard,nx)

# push particles with OpenMP:
# updates ppart and wke, and possibly ncl, ihole, irc
   wke[0] = 0.0
   mbpush1.wmbpush1(ppart,fxyze,byze,kpic,ncl,ihole,omx,qbme,dt,dth,ci,
                    wke,tpush,nx,mx,ipbc,relativity,list,irc)

# reorder particles by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)

   if (ntime==0):
      wt[0] = we + wf + wm
      print "Initial Total Field, Kinetic and Total Energies:"
      print "%14.7e %14.7e %14.7e" % (wt[0], wke, wke + wt[0])
      print "Initial Electrostatic, Transverse Electric and Magnetic " \
      "Field Energies:"
      print "%14.7e %14.7e %14.7e" % (we, wf, wm)

# loop time
   dtimer(dtime,ltime,1)
   tloop = tloop + float(dtime)

ntime = ntime + 1

# * * * end main iteration loop * * *

print "ntime, relativity = ", ntime, relativity
wt[0] = we + wf + wm
print "Final Total Field, Kinetic and Total Energies:"
print "%14.7e %14.7e %14.7e" % (wt[0], wke, wke + wt[0])
print "Final Electrostatic, Transverse Electric and Magnetic Field " \
"Energies:"
print "%14.7e %14.7e %14.7e" % (we, wf, wm)

print ""
print "initialization time = ", tinit
print "deposit time = ", tdpost[0]
print "current deposit time = ", tdjpost[0]
tdpost[0] = tdpost[0] + tdjpost[0]
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
