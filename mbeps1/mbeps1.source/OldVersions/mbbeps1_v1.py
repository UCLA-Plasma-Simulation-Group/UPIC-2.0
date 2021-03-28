#-----------------------------------------------------------------------
# 1-2/2D Electromagnetic OpenMP PIC code
# written by Viktor K. Decyk and Joshua Kelly, UCLA
# copyright 2016, regents of the university of california
import math
import numpy
from libmpush1 import *
from libmbpush1 import *
from fomplib import *
from fgraf1 import *
from dtimer import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# idimp = number of particle coordinates = 4
# ipbc = particle boundary condition: 1 = periodic
idimp = 4; ipbc = 1
# wke/wki/we = electron/ion kinetic energies and electric field energy
# wf/wb = magnetic field/transverse electric field
wke = numpy.zeros((1),float_type)
wki = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wf = numpy.zeros((1),float_type)
wb = numpy.zeros((1),float_type)
# list = (true,false) = list of particles leaving tiles found in push
list = True

# declare scalars for standard code
npi = 0
ws = numpy.zeros((1),float_type)

# declare scalars for OpenMP code
nppmx = numpy.empty((1),int_type)
irc = numpy.zeros((1),int_type)

# declare scalars for diagnostics
iuin = 8; iudm = 19
iude = 10; iup = 11; iuel = 12
iua = 13; iuet = 14; iub = 15; iuar = 16
iudi = 20; iuji = 21
wef = numpy.zeros((1),float_type)

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
tdiag = numpy.zeros((1),float_type)
dtime = numpy.empty((1),double_type)

# start timing initialization
dtimer(dtime,itime,-1)

# override default input data
in1.emf = 1
in1.relativity = 1
# read namelist
in1.readnml1(iuin)
# override input data
in1.idcode = 2
in1.ndim = 3

# start timing initialization
dtimer(dtime,itime,-1)

# create string from idrun
cdrun = str(in1.idrun)
# text output file
fname = "output1." + cdrun
iuot = open(fname,"w")

# in1.nvp = number of shared memory nodes (0=default)
#in1.nvp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
omplib.init_omp(in1.nvp)

# open graphics device
irc[0] = graf1.open_graphs(in1.nplot)

# initialize scalars for standard code
# increase number of coordinates for particle tag
if ((in1.ntt > 0) or ((in1.nts > 0) and (in1.ntsc > 0))):
   idimp += 1
# np = total number of electrons in simulation
np = in1.npx + in1.npxb;
# nx = number of grid points in x direction
nx = int(math.pow(2,in1.indx)); nxh = int(nx/2)
# npi = total number of ions in simulation
if (in1.movion > 0):
   npi = in1.npxi + in1.npxbi
nxe = nx + 2; nxeh = int(nxe/2)
# mx1 = number of tiles in x direction
mx1 = int((nx - 1)/in1.mx + 1)
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(in1.tend/in1.dt + .0001); ntime = 0
qbme = in1.qme
affp = float(nx)/float(np)
if (in1.movion==1):
   qbmi = in1.qmi/in1.rmass
   vtxi = in1.vtx/numpy.sqrt(in1.rmass*in1.rtempxi)
   vtyi = in1.vty/numpy.sqrt(in1.rmass*in1.rtempyi)
   vtzi = in1.vtz/numpy.sqrt(in1.rmass*in1.rtempzi)
   vtdxi = in1.vtdx/numpy.sqrt(in1.rmass*in1.rtempdxi)
   vtdyi = in1.vtdy/numpy.sqrt(in1.rmass*in1.rtempdyi)
   vtdzi = in1.vtdz/numpy.sqrt(in1.rmass*in1.rtempdzi)
dth = 0.0

# check for unimplemented features
if (list):
   if (ipbc != 1):
      print "ipbc != 1 and list = True not yet supported"
      print "list reset to False"
      list = False

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe),float_type,'F')
# qi = ion charge density with guard cells
qi = numpy.empty((nxe),float_type,'F')
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
# kpic = number of electrons in each tile
kpic = numpy.empty((mx1),int_type,'F')

# prepare fft tables
mfft1.mfft1_init(mixup,sct,in1.indx)
# calculate form factors
mfield1.mpois1_init(ffc,in1.ax,affp,nx)
# initialize different ensemble of random numbers
if (in1.nextrand > 0):
   minit1.mnextran1(in1.nextrand,in1.ndim,np+npi)

# initialize electons
# background electrons
if (in1.npx > 0):
#  minit1.mudistr1(part,1,in1.npx,nx,ipbc)
   minit1.mfdistr1(part,in1.ampdx,in1.scaledx,in1.shiftdx,1,in1.npx,nx,
                   ipbc,in1.ndprof)
   minit1.wmvdistr1h(part,1,in1.vtx,in1.vty,in1.vtz,in1.vx0,in1.vy0,
                     in1.vz0,in1.npx,in1.nvdist)
# beam electrons
if (in1.npxb > 0):
   it = in1.npx + 1
#  minit1.mudistr1(part,it,in1.npxb,nx,ipbc)
   minit1.mfdistr1(part,in1.ampdx,in1.scaledx,in1.shiftdx,it,in1.npxb,
                   nx,ipbc,in1.ndprof)
   minit1.wmvdistr1h(part,it,in1.vtdx,in1.vtdy,in1.vtdz,in1.vdx,in1.vdy,
                     in1.vdz,in1.npxb,in1.nvdist)

# mark electron beam particles
if ((in1.nts > 0) and (in1.ntsc > 0)):
   mdiag1.setmbeam1(part,in1.npx)

# find number of electrons in each of mx, tiles: updates kpic, nppmx
minit1.mdblkp2(part,kpic,nppmx,in1.mx,irc)

# allocate vector electron data
nppmx0 = int((1.0 + in1.xtras)*nppmx)
ntmax = int(in1.xtras*nppmx)
npbmx = int(in1.xtras*nppmx)
# ppart = tiled electron array
ppart = numpy.empty((idimp,nppmx0,mx1),float_type,'F')
# ppbuff = buffer array for reordering tiled particle array
ppbuff = numpy.empty((idimp,npbmx,mx1),float_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((2,mx1),int_type,'F')
# ihole = location/destination of each particle departing tile
ihole = numpy.empty((2,ntmax+1,mx1),int_type,'F')
# copy ordered electron data for OpenMP: updates ppart and kpic
mpush1.mpmovin1(part,ppart,kpic,in1.mx,irc)

# sanity check for electrons
mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)
#numpy.delete(part,part)

# initialize background charge density: updates qi
if (in1.movion==0):
   qi.fill(0.0)
   qmi = -in1.qme
   mpush1.mpost1(ppart,qi,kpic,qmi,tdpost,in1.mx)
   mgard1.maguard1(qi,tguard,nx)

# initialize ions
if (in1.movion==1):
   part = numpy.empty((idimp,npi),float_type,'F')
# kipic = number of ions in each tile
   kipic = numpy.empty((mx1),int_type,'F')
# cui = ion current density with guard cells
   cui = numpy.zeros((2,nxe),float_type,'F')
   it = in1.npxi + 1
# background ions
   if (in1.npxi > 0):
#     minit1.mudistr1(part,1,in1.npxi,nx,ipbc)
      minit1.mfdistr1(part,in1.ampdxi,in1.scaledxi,in1.shiftdxi,1,
                      in1.npxi,nx,ipbc,in1.ndprofi)
      minit1.wmvdistr1h(part,1,vtxi,vtyi,vtzi,in1.vxi0,in1.vyi0,
                        in1.vzi0,in1.npxi,in1.nvdist)
# beam ions
   if (in1.npxbi > 0):
#     minit1.mudistr1(part,it,in1.npxbi,nx,ipbc)
      minit1.mfdistr1(part,in1.ampdxi,in1.scaledxi,in1.shiftdxi,it,
                      in1.npxbi,nx,ipbc,in1.ndprofi)
      minit1.wmvdistr1h(part,it,vtdxi,vtdyi,vtdzi,in1.vdxi,in1.vdyi,
                        in1.vdzi,in1.npxbi,in1.nvdist)
      
# mark ion beam particles
   if ((in1.nts > 0) and (in1.ntsc > 0)):
      mdiag1.setmbeam1(part,in1.npxi)

# find number of ions in each of mx, tiles: updates kipic, nppmx
   minit1.mdblkp2(part,kipic,nppmx,in1.mx,irc)

# allocate vector ion data
   nppmx1 = int((1.0 + in1.xtras)*nppmx)
   pparti = numpy.empty((idimp,nppmx1,mx1),float_type,'F')
# copy ordered ion data for OpenMP: updates pparti and kipic
   mpush1.mpmovin1(part,pparti,kipic,in1.mx,irc)

# sanity check for ions
   mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)
#  numpy.delete(part,part)

# initialize transverse electromagnetic fields
eyz.fill(numpy.complex(0.0,0.0))
byz.fill(numpy.complex(0.0,0.0))
cue.fill(0.0)
# set magnitude of external transverse magnetic field
omt = numpy.sqrt(in1.omy*in1.omy + in1.omz*in1.omz)

# allocate diagnostic arrays
# cwk = labels for power spectrum display
cwk = numpy.array(["   W > 0  ","   W < 0  "],'S10')
# reverse simulation at end back to start
if (in1.treverse==1):
   nloop = 2*nloop

# energy time history
if (in1.ntw > 0):
   mtw = int((nloop - 1)/in1.ntw + 1); itw = 0
# wt = energy time history array
   wt = numpy.zeros((mtw,7),float_type,'F')
   s = numpy.zeros((7),double_type,'F')

# allocate scratch arrays for scalar fields
if ((in1.ntde > 0) or (in1.ntp > 0) or (in1.ntel > 0) or (in1.ntdi > 0)):
   sfieldc = numpy.empty((nxh),complex_type,'F')
   sfield = numpy.empty((nxe),float_type,'F')

# allocate and initialize frequency array for spectral analysis
if ((in1.ntp > 0) or (in1.nta > 0) or (in1.ntet > 0) or (in1.ntar > 0)):
   iw = int((in1.wmax - in1.wmin)/in1.dw + 1.5)
   wm = numpy.empty((iw),float_type,'F')
   wm[:] = in1.wmin + in1.dw*numpy.linspace(0,iw,iw)

# allocate and initialize frequency array for ion spectral analysis
if (in1.movion==1):
   if ((in1.ntdi > 0) or (in1.ntji > 0)):
      iwi = int((in1.wimax - in1.wimin)/in1.dwi + 1.5)
      wmi = numpy.empty((iwi),float_type,'F')
      wmi[:] = in1.wimin + in1.dwi*numpy.linspace(0,iwi,iwi)

# allocate scratch arrays for vector fields
if ((in1.nta>0) or (in1.ntet>0) or (in1.ntb>0) or (in1.ntar>0) or 
    (in1.ntji>0)):
   vfieldc = numpy.empty((2,nxh),complex_type,'F')
   vfield = numpy.empty((2,nxe),float_type,'F')

# initialize electron density diagnostic
if (in1.ntde > 0):
   fdename = "denek1." + cdrun
   in1.modesxde = int(min(in1.modesxde,nxh+1))
# denet/denit = store selected fourier modes for electron density
   denet = numpy.empty((in1.modesxde),complex_type,'F')
# open file: updates nderec and possibly iude
   if (in1.nderec==0):
      mdiag1.dafopenc1(denet,iude,in1.nderec,fdename)

# initialize ion density diagnostic
if (in1.movion==1):
   if (in1.ntdi > 0):
      fdiname = "denik1." + cdrun
      in1.modesxdi = int(min(in1.modesxdi,nxh+1))
# denit = store selected fourier modes for ion density
      denit = numpy.empty((in1.modesxdi),complex_type,'F')
# open file: updates ndirec and possibly iudi
      if (in1.ndirec==0):
         mdiag1.dafopenc1(denit,iudi,in1.ndirec,fdiname)
# ion spectral analysis
      if ((in1.nddi==2) or (in1.nddi==3)):
         mtdi = int((nloop - 1)/in1.ntdi) + 1; itdi = 0
# pkwdi = power spectrum for potential
         pkwdi = numpy.empty((in1.modesxdi,iwi,2),float_type,'F')
# pksdi = accumulated complex spectrum for potential
         pksdi = numpy.zeros((4,in1.modesxdi,iwi),double_type,'F')
# wkdi = maximum frequency as a function of k for potential
         wkdi = numpy.empty((in1.modesxdi,2),float_type,'F')

# initialize potential diagnostic
if (in1.ntp > 0):
   fpname = "potk1." + cdrun
   in1.modesxp = int(min(in1.modesxp,nxh+1))
# pott = store selected fourier modes for potential
   pott = numpy.empty((in1.modesxp),complex_type,'F')
# open file: updates nprec and possibly iup
   if (in1.nprec==0):
      mdiag1.dafopenc1(pott,iup,in1.nprec,fpname)
# spectral analysis
   if ((in1.ndp==2) or (in1.ndp==3)):
      mtp = int((nloop - 1)/in1.ntp) + 1; itp = 0
# pkw = power spectrum for potential
      pkw = numpy.empty((in1.modesxp,iw,2),float_type,'F')
# pks = accumulated complex spectrum for potential
      pks = numpy.zeros((4,in1.modesxp,iw),double_type,'F')
# wk = maximum frequency as a function of k for potential
      wk = numpy.empty((in1.modesxp,2),float_type,'F')

# initialize longitudinal efield diagnostic
if (in1.ntel > 0):
   felname = "elk1." + cdrun
   in1.modesxel = int(min(in1.modesxel,nxh+1))
# elt = store selected fourier modes for longitudinal efield
   elt = numpy.empty((in1.modesxel),complex_type,'F')
# open file: updates nelrec and possibly iuel
   if (in1.nelrec==0):
      mdiag1.dafopenc1(elt,iuel,in1.nelrec,felname)

# initialize ion current density diagnostic
if (in1.movion==1):
   if (in1.ntji > 0):
      fjiname = "curik1." + cdrun
      in1.modesxji = int(min(in1.modesxji,nxh+1))
# curit = store selected fourier modes for ion current density
      curit = numpy.empty((2,in1.modesxji),complex_type,'F')
# open file: updates njirec and possibly iuji
      if (in1.njirec==0):
         mdiag1.dafopenvc1(curit,iuji,in1.njirec,fjiname)
# ion spectral analysis
      if ((in1.ndji==2) or (in1.ndji==3)):
         mtji = int((nloop - 1)/in1.ntji) + 1; itji = 0
# vpkwji = power spectrum for ion current density
         vpkwji = numpy.empty((2,in1.modesxji,iwi,2),float_type,'F')
# vpksji = accumulated complex spectrum for ion current density
         vpksji = numpy.zeros((2,4,in1.modesxji,iwi),double_type,'F')
# vwkji = maximum frequency as a function of k for ion current
         vwkji = numpy.empty((2,in1.modesxji,2),float_type,'F')

# initialize radiative vector potential diagnostic
if (in1.ntar > 0):
   farname = "vpotrk1." + cdrun
   in1.modesxar = int(min(in1.modesxar,nxh+1))
# vpotr = store selected fourier modes for radiative vector potential
   vpotr = numpy.empty((2,in1.modesxar),complex_type,'F')
# open file: updates narec and possibly iua
   if (in1.narec==0):
      mdiag1.dafopenvc1(vpotr,iuar,in1.narrec,farname)
# spectral analysis
   if ((in1.ndar==2) or (in1.ndar==3)):
      mtar = int((nloop - 1)/in1.ntar) + 1; itar = 0
# vpkwr = power spectrum for radiative vector potential
      vpkwr = numpy.empty((2,in1.modesxar,iw,2),float_type,'F')
# vpksr = accumulated complex spectrum for radiative vector potential
      vpksr = numpy.zeros((2,4,in1.modesxar,iw),double_type,'F')
# vwkr = maximum frequency as a function of k for radiative vector
#        potential
      vwkr = numpy.empty((2,in1.modesxar,2),float_type,'F')
# oldcu = previous current density with guard cells
      oldcu = numpy.zeros((2,nxe),float_type,'F')

# initialize vector potential diagnostic
if (in1.nta > 0):
   faname = "vpotk1." + cdrun
   in1.modesxa = int(min(in1.modesxa,nxh+1))
# vpott = store selected fourier modes for vector potential
   vpott = numpy.empty((2,in1.modesxa),complex_type,'F')
# open file: updates narec and possibly iua
   if (in1.narec==0):
      mdiag1.dafopenvc1(vpott,iua,in1.narec,faname)
# spectral analysis
   if ((in1.nda==2) or (in1.nda==3)):
      mta = int((nloop - 1)/in1.nta) + 1; ita = 0
# vpkw = power spectrum for vector potential
      vpkw = numpy.empty((2,in1.modesxa,iw,2),float_type,'F')
# vpks = accumulated complex spectrum for vector potential
      vpks = numpy.zeros((2,4,in1.modesxa,iw),double_type,'F')
# vwk = maximum frequency as a function of k for vector potential
      vwk = numpy.empty((2,in1.modesxa,2),float_type,'F')

# initialize transverse efield diagnostic
if (in1.ntet > 0):
   fetname = "etk1." + cdrun
   in1.modesxet = int(min(in1.modesxet,nxh+1))
# ett = store selected fourier modes for transverse efield
   ett = numpy.empty((2,in1.modesxet),complex_type,'F')
# open file: updates netrec and possibly iuet
   if (in1.netrec==0):
      mdiag1.dafopenvc1(ett,iuet,in1.netrec,fetname)
# spectral analysis
   if ((in1.ndet==2) or (in1.ndet==3)):
      mtet = int((nloop - 1)/in1.ntet) + 1; itet = 0
# vpkwet = power spectrum for transverse efield
      vpkwet = numpy.empty((2,in1.modesxet,iw,2),float_type,'F')
# vpkset = accumulated complex spectrum for transverse efield
      vpkset = numpy.zeros((2,4,in1.modesxet,iw),double_type,'F')
# vwket = maximum frequency as a function of k for transverse efield
      vwket = numpy.empty((2,in1.modesxet,2),float_type,'F')

# initialize magnetic field diagnostic
if (in1.ntb > 0):
   fbname = "bk1." + cdrun
   in1.modesxb = int(min(in1.modesxb,nxh+1))
# bt = store selected fourier modes for magnetic field
   bt = numpy.empty((2,in1.modesxb),complex_type,'F')
# open file: updates nbrec and possibly iub
   if (in1.netrec==0):
      mdiag1.dafopenvc1(bt,iub,in1.nbrec,fbname)

# initialize velocity diagnostic
if (in1.ntv > 0):
# sfv = electron velocity distribution functions in tile
   sfv = numpy.empty((2*in1.nmv+2,in1.ndim,mx1+1),float_type,'F')
# fvm = electron vdrift, vth, entropy for global distribution
   fvm = numpy.empty((in1.ndim,3),float_type,'F')
   mtv = int((nloop - 1)/in1.ntv) + 1; itv = 0
# fvtm = time history of electron vdrift, vth, and entropy
   fvtm = numpy.zeros((mtv,in1.ndim,3),float_type,'F')
   ws[0] = 2.0*max(4.0*in1.vtx+abs(in1.vx0),4.0*in1.vtdx+abs(in1.vdx))
   ws[0] = max(ws[0],2.0*max(4.0*in1.vty+abs(in1.vy0),
                             4.0*in1.vtdy+abs(in1.vdy)))
   ws[0] = max(ws[0],2.0*max(4.0*in1.vtz+abs(in1.vz0),
                             4.0*in1.vtdz+abs(in1.vdz)))
   sfv[0,0,:] = ws[0]
   sfv[0,1,:] = ws[0]
   sfv[0,2,:] = ws[0]
# ions
   if (in1.movion==1):
# sfvi = ion velocity distribution functions in tile
      sfvi = numpy.empty((2*in1.nmv+2,in1.ndim,mx1+1),float_type,'F')
# fvmi = ion vdrift, vth, entropy for global distribution
      fvmi = numpy.empty((in1.ndim,3),float_type,'F')
# fvtmi = time history of ion vdrift, vth, and entropy
      fvtmi = numpy.zeros((mtv,in1.ndim,3),float_type,'F')
      ws[0] = 2.0*max(4.0*vtxi+abs(in1.vxi0),4.0*vtdxi+abs(in1.vdxi))
      ws[0] = max(ws[0],2.0*max(4.0*vtyi+abs(in1.vyi0),
                                4.0*vtdyi+abs(in1.vdyi)))
      ws[0] = max(ws[0],2.0*max(4.0*vtzi+abs(in1.vzi0),
                                4.0*vtdzi+abs(in1.vdzi)))
      sfvi[0,0,:] = ws[0]
      sfvi[0,1,:] = ws[0]
      sfvi[0,2,:] = ws[0]

# initialize trajectory diagnostic
if (in1.ntt > 0):
# iprobt = scratch array 
   iprobt = numpy.empty((in1.nprobt),numpy.int32)
   mdiag1.setptraj1(ppart,kpic,iprobt,in1.nst,in1.vtx,in1.vtsx,in1.dvtx,
                    np,in1.nprobt)
   if (in1.nprobt > 16777215):
      print "nprobt overflow = ", in1.nprobt
      exit(1)
# partt = particle trajectories tracked
   partt = numpy.empty((idimp,in1.nprobt),float_type,'F')
   if ((in1.nst==1) or (in1.nst==2)):
      it = int((nloop - 1)/in1.ntt + 1); itt = 0
# partd = trajectory time history array
      partd = numpy.empty((it,idimp,in1.nprobt),float_type,'F')
   elif (in1.nst==3):
# fvtp = velocity distribution function for test particles
      fvtp = numpy.empty((2*in1.nmv+2,in1.ndim),float_type,'F')
# fvmtp = vdrift, vth, and entropy for test particles
      fvmtp = numpy.empty((in1.ndim,3),float_type,'F')
      fvtp[0,:] = 2.0*max(4.0*in1.vtx+abs(in1.vx0),
                          4.0*in1.vtdx+abs(in1.vdx))

# initialization time
dtimer(dtime,itime,1)
tinit = tinit + float(dtime)
# start timing loop
dtimer(dtime,ltime,-1)

if (in1.dt > 0.64*in1.ci):
   print "Warning: Courant condition may be exceeded!"

print >> iuot, "program mbbeps1"

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
   print >> iuot, "ntime = ", ntime

# save previous current in fourier space for radiative vector potential
   if (in1.ntar > 0):
      it = ntime/in1.ntar
      if (ntime==in1.ntar*it):
         oldcu[:] = numpy.copy(cue)

# deposit current with OpenMP:
# updates ppart and cue, and possibly ncl, ihole, irc
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   dtimer(dtime,itime,1)
   tdjpost[0] = tdjpost[0] + float(dtime)
   mcurd1.wmdjpost1(ppart,cue,kpic,ncl,ihole,in1.qme,dth,in1.ci,tdjpost,
                    nx,in1.mx,ipbc,in1.relativity,list,irc)
# add guard cells: updates cue
   mgard1.macguard1(cue,tguard,nx)

# reorder electrons by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,list,
                    irc)
# sanity check for electrons
   mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)

# deposit ion current with OpenMP:
   if (in1.movion==1):
# updates pparti and cui, and possibly ncl, ihole, irc
      dtimer(dtime,itime,-1)
      cui.fill(0.0)
      dtimer(dtime,itime,1)
      tdjpost[0] = tdjpost[0] + float(dtime)
      mcurd1.wmdjpost1(pparti,cui,kipic,ncl,ihole,in1.qmi,dth,in1.ci,
                       tdjpost,nx,in1.mx,ipbc,in1.relativity,list,irc)
# add guard cells: updates cui
      mgard1.macguard1(cui,tguard,nx)

# reorder ions by tile with OpenMP:
# updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
      msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,in1.mx,
                       list,irc)
# sanity check for ions
      mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)

# ion current density diagnostic
   if (in1.movion==1):
      if (in1.ntji > 0):
         it = ntime/in1.ntji
         if (ntime==in1.ntji*it):
            vfield[:] = numpy.copy(cui)
# transform ion current density to fourier space: updates vfield
            isign = -1
            mfft1.mfft1rn(vfield,isign,mixup,sct,tfft,in1.indx)
# calculate smoothed ion current in fourier space: updates vfieldc
            mfield1.msmooth13(vfield,vfieldc,ffc,tfield,nx)
# store selected fourier modes: updates curit
            mfield1.mrdvmodes1(vfieldc,curit,tfield,nx,in1.modesxji)
# write diagnostic output: updates ndirec
            mdiag1.dafwritevc1(curit,tdiag,iuji,in1.njirec,in1.modesxji)
# transform smoothed ion current to real space: updates vfield
            if ((in1.ndji==1) or (in1.ndji==3)):
               mfft1.mfft1crn(vfieldc,vfield,mixup,sct,tfft,in1.indx)
               mgard1.mcguard1(vfield,tguard,nx)
# display smoothed ion current
               graf1.dvector1(vfield,' ION CURRENT',ntime,999,0,2,nx,
                              irc)
               if (irc[0]==1):
                  break
               irc[0] = 0
# ion spectral analysis
            if ((in1.ndji==2) or (in1.ndji==3)):
               itji += 1
               ts = in1.dt*float(ntime)
               mdiag1.mivcspect1(curit,wmi,vpkwji,vpksji,ts,in1.t0,
                                 tdiag,mtji,iwi,in1.modesxji,nx,-1)
# performs frequency analysis of accumulated complex time series
               vwkji[0,:,0] = wmi[numpy.argmax(vpkwji[0,:,:,0],axis=1)]
               vwkji[1,:,0] = wmi[numpy.argmax(vpkwji[1,:,:,0],axis=1)]
               vwkji[0,:,1] = wmi[numpy.argmax(vpkwji[0,:,:,1],axis=1)]
               vwkji[1,:,1] = wmi[numpy.argmax(vpkwji[1,:,:,1],axis=1)]
# display frequency spectrum
               graf1.dmvector1(vwkji,'ION CURRENT OMEGA VS MODE',ntime,
                               999,2,2,in1.modesxji,cwk,irc)
               if (irc[0]==1):
                     break
               irc[0] = 0

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   dtimer(dtime,itime,1)
   tdpost[0] = tdpost[0] + float(dtime)
   mpush1.mpost1(ppart,qe,kpic,in1.qme,tdpost,in1.mx)
# add guard cells: updates qe
   mgard1.maguard1(qe,tguard,nx)

# electron density diagnostic
   if (in1.ntde > 0):
      it = int(ntime/in1.ntde)
      if (ntime==in1.ntde*it):
         sfield[:] = -numpy.copy(qe)
# transform electron density to fourier space: updates sfield
         isign = -1
         mfft1.mfft1r(sfield,isign,mixup,sct,tfft,in1.indx)
# calculate smoothed density in fourier space: updates sfieldc
         mfield1.msmooth1(sfield,sfieldc,ffc,tfield,nx)
# store selected fourier modes: updates denet
         mfield1.mrdmodes1(sfieldc,denet,tfield,nx,in1.modesxde)
# write diagnostic output: updates nderec
         mdiag1.dafwritec1(denet,tdiag,iude,in1.nderec,in1.modesxde)
# transform smoothed electron density to real space: updates sfield
         mfft1.mfft1cr(sfieldc,sfield,mixup,sct,tfft,in1.indx)
         mgard1.mdguard1(sfield,tguard,nx)
# display smoothed electron density
         graf1.dscaler1(sfield,' EDENSITY',ntime,999,0,nx,irc)
         if (irc[0]==1):
            break
         irc[0] = 0

# deposit ion charge with OpenMP: updates qi
   if (in1.movion==1):
      dtimer(dtime,itime,-1)
      qi.fill(0.0)
      dtimer(dtime,itime,1)
      tdpost[0] = tdpost[0] + float(dtime)
      mpush1.mpost1(pparti,qi,kipic,in1.qmi,tdpost,in1.mx)
# add guard cells: updates qi
      mgard1.maguard1(qi,tguard,nx)

# ion density diagnostic
   if (in1.movion==1):
      if (in1.ntdi > 0):
         it = int(ntime/in1.ntdi)
         if (ntime==in1.ntdi*it):
            sfield[:] = numpy.copy(qi)
# transform ion density to fourier space: updates sfield
            isign = -1
            mfft1.mfft1r(sfield,isign,mixup,sct,tfft,in1.indx)
# calculate smoothed density in fourier space: updates sfieldc
            mfield1.msmooth1(sfield,sfieldc,ffc,tfield,nx)
# store selected fourier modes: updates denit
            mfield1.mrdmodes1(sfieldc,denit,tfield,nx,in1.modesxdi)
# write diagnostic output: updates ndirec
            mdiag1.dafwritec1(denit,tdiag,iudi,in1.ndirec,in1.modesxdi)
# transform smoothed ion density to real space: updates sfield
            if ((in1.nddi==1) or (in1.nddi==3)):
               mfft1.mfft1cr(sfieldc,sfield,mixup,sct,tfft,in1.indx)
               mgard1.mdguard1(sfield,tguard,nx)
# display smoothed ion density
               graf1.dscaler1(sfield,' ION DENSITY',ntime,999,1,nx,irc)
               if (irc[0]==1):
                  break
               irc[0] = 0
# ion spectral analysis
            if ((in1.nddi==2) or (in1.nddi==3)):
               itdi += 1
               ts = in1.dt*float(ntime)
               mdiag1.micspect1(denit,wmi,pkwdi,pksdi,ts,in1.t0,tdiag,
                                mtdi,iwi,in1.modesxdi,nx,-1)
# performs frequency analysis of accumulated complex time series
               wkdi[:,0] = wmi[numpy.argmax(pkwdi[:,:,0],axis=1)]
               wkdi[:,1] = wmi[numpy.argmax(pkwdi[:,:,1],axis=1)]
# display frequency spectrum
               graf1.dmscaler1(wkdi,'ION DENSITY OMEGA VS MODE',ntime,
                               999,1,in1.modesxdi,cwk,irc)
               if (irc[0]==1):
                     break
               irc[0] = 0

# add electron and ion densities: updates qe
   mfield1.maddqei1(qe,qi,tfield,nx)

# add electron and ion current densities: updates cue
   if (in1.movion==1):
      mfield1.maddcuei1(cue,cui,tfield,nx)

# transform charge to fourier space: updates qe
   isign = -1
   mfft1.mfft1r(qe,isign,mixup,sct,tfft,in1.indx)

# transform current to fourier space: updates cue
   isign = -1
   mfft1.mfft1rn(cue,isign,mixup,sct,tfft,in1.indx)

# radiative vector potential diagnostic
   if (in1.ntar > 0):
      it = ntime/in1.ntar
      if (ntime==in1.ntar*it):
# average current: updates vfieldc = 0.5*(cue + oldcu)
         mfield1.mcuave1(vfieldc,cue,oldcu,tfield,nx)
# calculate radiative vector potential in fourier space: updates vfieldc
# vfieldc should contain averaged current on entry
         mfield1.mavrpot1(vfieldc,byz,ffc,in1.ci,tfield,nx)
# store selected fourier modes: updates vpotr
         mfield1.mrdvmodes1(vfieldc,vpotr,tfield,nx,in1.modesxar)
# write diagnostic output: updates narrec
         mdiag1.dafwritevc1(vpotr,tdiag,iuar,in1.narrec,in1.modesxar)
# transform radiative vector potential to real space: updates vfield
         if ((in1.ndp==1) or (in1.ndp==3)):
            mfft1.mfft1crn(vfieldc,vfield,mixup,sct,tfft,in1.indx)
            mgard1.mcguard1(vfield,tguard,nx)
# display radiative vector potential
            graf1.dvector1(vfield,' RADIATIVE VECTOR POTENTIAL',ntime,
                           999,0,2,nx,irc)
            if (irc[0]==1):
               break
            irc[0] = 0
# spectral analysis
         if ((in1.ndar==2) or (in1.ndar==3)):
            itar += 1
            ts = in1.dt*float(ntime)
            mdiag1.mivcspect1(vpotr,wm,vpkwr,vpksr,ts,in1.t0,tdiag,mtar,
                              iw,in1.modesxar,nx,1)
# performs frequency analysis of accumulated complex vector time series
            vwkr[0,:,0] = wm[numpy.argmax(vpkwr[0,:,:,0],axis=1)]
            vwkr[1,:,0] = wm[numpy.argmax(vpkwr[1,:,:,0],axis=1)]
            vwkr[0,:,1] = wm[numpy.argmax(vpkwr[0,:,:,1],axis=1)]
            vwkr[1,:,1] = wm[numpy.argmax(vpkwr[1,:,:,1],axis=1)]
# display frequency spectrum
            graf1.dmvector1(vwkr,'RADIATIVE VPOTENTIAL OMEGA VS MODE',
                            ntime,999,2,2,in1.modesxar,cwk,irc)
            if (irc[0]==1):
               break
            irc[0] = 0

# calculate electromagnetic fields in fourier space: updates eyz, byz
   if (ntime==0):
# initialize electromagnetic fields from darwin fields
# calculate initial darwin magnetic field
      mfield1.mibpois1(cue,byz,ffc,in1.ci,wb,tfield,nx)
      wf[0] = 0.0
# calculate initial darwin electric field
      amu = numpy.empty((2,nxe),float_type,'F')
      dcu = numpy.empty((2,nxe),float_type,'F')
      amu.fill(0.0)
      mcurd1.wmgmjpost1(ppart,amu,kpic,in1.qme,in1.ci,tdjpost,in1.mx,
                        in1.relativity)
      mgard1.macguard1(amu,tguard,nx)
      isign = -1
      mfft1.mfft1rn(amu,isign,mixup,sct,tfft,in1.indx)
      mfield1.mdcuperp1(dcu,amu,tfield,nx)
      numpy.delete(amu,amu)
      mfield1.metfield1(dcu,eyz,ffc,in1.ci,wf,tfield,nx)
      numpy.delete(dcu,dcu)
      dth = 0.5*in1.dt
   else:
      mfield1.mmaxwel1(eyz,byz,cue,ffc,in1.ci,in1.dt,wf,wb,tfield,nx)

# calculate longitudinal force/charge in fourier space:
# updates fxe, we
   mfield1.mpois1(qe,fxe,ffc,we,tfield,nx)

# add longitudinal and transverse electric fields: updates fxyze
   mfield1.memfield1(fxyze,fxe,eyz,ffc,tfield,nx)
# copy magnetic field: updates byze
   mfield1.mbmfield1(byze,byz,ffc,tfield,nx)

# transform electric force to real space: updates fxyze
   isign = 1
   mfft1.mfft1rn(fxyze,isign,mixup,sct,tfft,in1.indx)

# transform magnetic force to real space: updates byze
   isign = 1
   mfft1.mfft1rn(byze,isign,mixup,sct,tfft,in1.indx)

# add external traveling wave field
   ts = in1.dt*float(ntime)
   mfield1.meaddext13(fxyze,tfield,in1.amodex,in1.freq,ts,in1.trmp,
                      in1.toff,in1.el0,in1.er0,nx)

# add constant to magnetic field with OpenMP: updates bxyze
   if (omt > 0.0):
      mfield1.mbaddext1(byze,tfield,in1.omy,in1.omz,nx)

# copy guard cells: updates fxyze, byze
   mgard1.mcguard1(fxyze,tguard,nx)
   mgard1.mcguard1(byze,tguard,nx)

# potential diagnostic
   if (in1.ntp > 0):
      it = int(ntime/in1.ntp)
      if (ntime==in1.ntp*it):
# calculate potential in fourier space: updates sfieldc
         mfield1.mpot1(qe,sfieldc,ffc,ws,tfield,nx)
# store selected fourier modes: updates pott
         mfield1.mrdmodes1(sfieldc,pott,tfield,nx,in1.modesxp)
# write diagnostic output: updates nprec
         mdiag1.dafwritec1(pott,tdiag,iup,in1.nprec,in1.modesxp)
# transform potential to real space: updates sfield
         if ((in1.ndp==1) or (in1.ndp==3)):
            mfft1.mfft1cr(sfieldc,sfield,mixup,sct,tfft,in1.indx)
            mgard1.mdguard1(sfield,tguard,nx)
# display potential
            graf1.dscaler1(sfield,' POTENTIAL',ntime,999,0,nx,irc)
            if (irc[0]==1):
               break
            irc[0] = 0
# spectral analysis
         if ((in1.ndp==2) or (in1.ndp==3)):
            itp += 1
            ts = in1.dt*float(ntime)
            mdiag1.micspect1(pott,wm,pkw,pks,ts,in1.t0,tdiag,mtp,iw,
                             in1.modesxp,nx,1)
# performs frequency analysis of accumulated complex time series
            wk[:,0] = wm[numpy.argmax(pkw[:,:,0],axis=1)]
            wk[:,1] = wm[numpy.argmax(pkw[:,:,1],axis=1)]
# display frequency spectrum
            graf1.dmscaler1(wk,'POTENTIAL OMEGA VS MODE',ntime,999,2,
                            in1.modesxp,cwk,irc)
            if (irc[0]==1):
               break
            irc[0] = 0

# longitudinal efield diagnostic
   if (in1.ntel > 0):
      it = int(ntime/in1.ntel)
      if (ntime==in1.ntel*it):
# calculate longitudinal efield in fourier space: updates sfieldc
         mfield1.melfield1(qe,sfieldc,ffc,ws,tfield,nx)
# store selected fourier modes: updates elt
         mfield1.mrdmodes1(sfieldc,elt,tfield,nx,in1.modesxel)
# write diagnostic output: updates nelrec
         mdiag1.dafwritec1(elt,tdiag,iuel,in1.nelrec,in1.modesxel)
# transform longitudinal efield to real space: updates sfield
         if ((in1.ndel==1) or (in1.ndel==3)):
            mfft1.mfft1cr(sfieldc,sfield,mixup,sct,tfft,in1.indx)
            mgard1.mdguard1(sfield,tguard,nx)
# display longitudinal efield
            graf1.dscaler1(sfield,' ELFIELD',ntime,999,0,nx,irc)
            if (irc[0]==1):
               break
            irc[0] = 0

# vector potential diagnostic
   if (in1.nta > 0):
      it = ntime/in1.nta
      if (ntime==in1.nta*it):
# calculate vector potential in fourier space: updates vfieldc
         mfield1.mavpot1(byz,vfieldc,tfield,nx)
# store selected fourier modes: updates vpott
         mfield1.mrdvmodes1(vfieldc,vpott,tfield,nx,in1.modesxa)
# write diagnostic output: updates narec
         mdiag1.dafwritevc1(vpott,tdiag,iua,in1.narec,in1.modesxa)
# transform vector potential to real space: updates vfield
         if ((in1.nda==1) or (in1.nda==3)):
            mfft1.mfft1crn(vfieldc,vfield,mixup,sct,tfft,in1.indx)
            mgard1.mcguard1(vfield,tguard,nx)
# display vector potential
            graf1.dvector1(vfield,' VECTOR POTENTIAL',ntime,999,0,2,nx,
                           irc)
            if (irc[0]==1):
               break
            irc[0] = 0
# spectral analysis
         if ((in1.nda==2) or (in1.nda==3)):
            ita += 1
            ts = in1.dt*float(ntime)
            mdiag1.mivcspect1(vpott,wm,vpkw,vpks,ts,in1.t0,tdiag,mta,iw,
                              in1.modesxa,nx,1)
# performs frequency analysis of accumulated complex vector time series
            vwk[0,:,0] = wm[numpy.argmax(vpkw[0,:,:,0],axis=1)]
            vwk[1,:,0] = wm[numpy.argmax(vpkw[1,:,:,0],axis=1)]
            vwk[0,:,1] = wm[numpy.argmax(vpkw[0,:,:,1],axis=1)]
            vwk[1,:,1] = wm[numpy.argmax(vpkw[1,:,:,1],axis=1)]
# display frequency spectrum
            graf1.dmvector1(vwk,'VECTOR POTENTIAL OMEGA VS MODE',ntime,
                            999,2,2,in1.modesxa,cwk,irc)
            if (irc[0]==1):
               break
            irc[0] = 0

# transverse efield diagnostic
   if (in1.ntet > 0):
      it = ntime/in1.ntet
      if (ntime==in1.ntet*it):
# store selected fourier modes: updates ett
         mfield1.mrdvmodes1(eyz,ett,tfield,nx,in1.modesxet)
# write diagnostic output: updates netrec
         mdiag1.dafwritevc1(ett,tdiag,iuet,in1.netrec,in1.modesxet)
# transform transverse efield to real space: updates vfield
         if ((in1.ndet==1) or (in1.ndet==3)):
            mfft1.mfft1crn(eyz,vfield,mixup,sct,tfft,in1.indx)
            mgard1.mcguard1(vfield,tguard,nx)
# display transverse efield
            graf1.dvector1(vfield,' TRANSVERSE EFIELD',ntime,999,0,2,nx,
                           irc)
            if (irc[0]==1):
               break
            irc[0] = 0
# spectral analysis
         if ((in1.ndet==2) or (in1.ndet==3)):
            itet += 1
            ts = in1.dt*float(ntime)
            mdiag1.mivcspect1(ett,wm,vpkwet,vpkset,ts,in1.t0,tdiag,mtet,
                              iw,in1.modesxet,nx,0)
# performs frequency analysis of accumulated complex vector time series
            vwket[0,:,0] = wm[numpy.argmax(vpkwet[0,:,:,0],axis=1)]
            vwket[1,:,0] = wm[numpy.argmax(vpkwet[1,:,:,0],axis=1)]
            vwket[0,:,1] = wm[numpy.argmax(vpkwet[0,:,:,1],axis=1)]
            vwket[1,:,1] = wm[numpy.argmax(vpkwet[1,:,:,1],axis=1)]
# display frequency spectrum
            graf1.dmvector1(vwket,'TRANSVERSE EFIELD OMEGA VS MODE',
                            ntime,999,2,2,in1.modesxet,cwk,irc)
            if (irc[0]==1):
               break
            irc[0] = 0

# magnetic field diagnostic
   if (in1.ntb > 0):
      it = ntime/in1.ntb
      if (ntime==in1.ntb*it):
# store selected fourier modes: updates bt
         mfield1.mrdvmodes1(byz,bt,tfield,nx,in1.modesxb)
# write diagnostic output: updates nbrec
         mdiag1.dafwritevc1(bt,tdiag,iub,in1.nbrec,in1.modesxb)
# transform magnetic field to real space: updates vfield
         mfft1.mfft1crn(byz,vfield,mixup,sct,tfft,in1.indx)
         mgard1.mcguard1(vfield,tguard,nx)
# display magnetic field
         graf1.dvector1(vfield,' MAGNETIC FIELD',ntime,999,0,2,nx,irc)
         if (irc[0]==1):
            break
         irc[0] = 0

# velocity diagnostic
   if (in1.ntv > 0):
      it = int(ntime/in1.ntv)
      if (ntime==in1.ntv*it):
# calculate electron distribution function and moments
         mdiag1.mvpdist1(ppart,kpic,sfv,fvm,tdiag,np,in1.nmv)
# store time history electron vdrift, vth, and entropy
         fvtm[itv,:,:] = fvm
# display electron velocity distributions
         graf1.displayfv1(sfv[:,:,mx1],fvm,' ELECTRON',ntime,in1.nmv,2,
                          irc)
         if (irc[0]==1):
            break
         irc[0] = 0
# ion distribution function
         if (in1.movion==1):
            mdiag1.mvpdist1(pparti,kipic,sfvi,fvmi,tdiag,npi,in1.nmv)
# store time history ion vdrift, vth, and entropy
            fvtmi[itv,:,:] = fvmi
# display ion velocity distributions
            graf1.displayfv1(sfvi[:,:,mx1],fvmi,' ION',ntime,in1.nmv,2,
                             irc)
            if (irc[0]==1):
               break
            irc[0] = 0
         itv += 1

# trajectory diagnostic
   if (in1.ntt > 0):
      it = int(ntime/in1.ntt)
      if (ntime==in1.ntt*it):
# copies trajectories to array partt
         mdiag1.mptraj1(ppart,kpic,partt,tdiag)
         if ((in1.nst==1) or (in1.nst==2)):
            partd[itt,:,:] = partt
            itt += 1
         elif (in1.nst==3):
# calculate particle distribution function and moments
            mdiag1.mvdist1(partt,fvtp,fvmtp,tdiag,in1.nprobt,in1.nmv)
# display velocity distributions
            graf1.displayfv1(fvtp,fvmtp,' ELECTRON',ntime,in1.nmv,1,irc)
            if (irc[0]==1):
               break
            irc[0] = 0

# phase space diagnostic
   if (in1.nts > 0):
      it = int(ntime/in1.nts)
      if (ntime==in1.nts*it):
# plot electrons vx versus x
         graf1.dpmgrasp1(ppart,kpic,' ELECTRON',ntime,999,nx,2,1,
                         in1.ntsc,irc)
         if (irc[0]==1):
            break
         irc[0] = 0
# ion phase space
         if (in1.movion==1):
# plot electrons vx versus x
            graf1.dpmgrasp1(pparti,kipic,' ION',ntime,999,nx,2,1,
                            in1.ntsc,irc)
            if (irc[0]==1):
               break
            irc[0] = 0
 
# push electrons with OpenMP:
   wke[0] = 0.0
# updates ppart and wke, and possibly ncl, ihole, irc
   if (in1.mzf==0):
      mbpush1.wmbpush1(ppart,fxyze,byze,kpic,ncl,ihole,in1.omx,qbme,
                       in1.dt,dth,in1.ci,wke,tpush,nx,in1.mx,ipbc,
                       in1.relativity,list,irc)
# zero force: updates ppart, wke and possibly ncl, ihole, and irc
   else:
      mbpush1.wmpush1zf(ppart,kpic,ncl,ihole,in1.dth,in1.ci,wke,tpush,
                        in1.nx,in1.mx,ipbc,in1.relativity,list,irc)

# reorder electrons by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,list,
                    irc)

# sanity check for electrons
   mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)

# push ions with OpenMP:
   if (in1.movion==1):
      wki[0] = 0.0
# updates pparti and wki, and possibly ncl, ihole, irc
      if (in1.mzf==0):
         mbpush1.wmbpush1(pparti,fxyze,byze,kipic,ncl,ihole,in1.omx,
                          qbmi,in1.dt,dth,in1.ci,wki,tpush,nx,in1.mx,
                          ipbc,in1.relativity,list,irc)
      else:
         mpush1.wmpush1zf(pparti,kipic,ncl,ihole,in1.dth,in1.ci,wki,
                          tpush,nx,in1.mx,ipbc,in1.relativity,list,irc)
      wki[0] = wki[0]*in1.rmass

# reorder ions by tile with OpenMP:
# updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
      msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,in1.mx,
                       list,irc)

# sanity check for ions
      mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)

# start running simulation backwards:
# need to reverse time lag in leap-frog integration scheme
   if (in1.treverse==1):
      if (((ntime+1)==(nloop/2)) or ((ntime+1)==nloop)):
         in1.dt = -in1.dt
         ws[0] = 0.0
         mpush1.wmpush1zf(ppart,kpic,ncl,ihole,dth,in1.ci,ws,tpush,nx,
                          in1.mx,ipbc,in1.relativity,list,irc)
         msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,
                          list,irc)
         if (in1.movion==1):
            mpush1.wmpush1zf(pparti,kipic,ncl,ihole,dth,in1.ci,ws,tpush,
                             nx,in1.mx,ipbc,in1.relativity,list,irc)
            msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,
                             in1.mx,list,irc)
# still need to reverse field solver

# energy diagnostic
   if (in1.ntw > 0):
      it = int(ntime/in1.ntw)
      if (ntime==in1.ntw*it):
         wef[0] = we[0] + wf[0] + wb[0]
         ws[0] = wef[0] + wke[0] + wki[0]
         if (ntime==0):
            s[5] = ws[0]
         print >> iuot, "Total Field, Kinetic and Total Energies:"
         if (in1.movion==0):
            iuot.write("%14.7e %14.7e %14.7e\n" % (wef[0],wke[0],ws[0])) 
         else: 
            iuot.write("%14.7e %14.7e %14.7e %14.7e\n" % (wef[0],wke[0],
                       wki[0],ws[0])) 
         print >> iuot, "Electric(l,t) and Magnetic Energies = "
         iuot.write("%14.7e %14.7e %14.7e\n" % (we[0],wf[0],wb[0]))
         wt[itw,:] = [wef[0],wke[0],wki[0],ws[0],we[0],wf[0],wb[0]]
         itw += 1
         s[0] += we[0]
         s[1] += wke[0]
         s[2] += wf[0]
         s[3] += wb[0]
         s[4] += wki[0]
         s[5] = min(s[5],float(ws[0]))
         s[6] = max(s[6],float(ws[0]))

ntime = ntime + 1

# loop time
dtimer(dtime,ltime,1)
tloop = tloop + float(dtime)

# * * * end main iteration loop * * *

print >> iuot
print >> iuot, "ntime, relativity = ", ntime, ",", in1.relativity
if (in1.treverse==1):
   print >> iuot, "treverse = ", in1.treverse
   
if ((in1.ntw > 0) or (in1.ntt > 0)):
   graf1.reset_graphs()

# trajectory diagnostic
if (in1.ntt > 0):
   if ((in1.nst==1) or (in1.nst==2)):
      if (in1.nplot > 0):
         irc[0] = graf1.open_graphs(1)
      ts = in1.t0 + in1.dt*float(in1.ntt)
      graf1.displaytr1(partd,ts,in1.dt*float(in1.ntt),itt,2,3,irc)
      if (irc[0]==1):
         exit(1)
      graf1.reset_nplot(nplot,irc)

# energy diagnostic
if (in1.ntw > 0):
   ts = in1.t0 + in1.dt*float(in1.ntw)
   graf1.displayw1(wt,ts,in1.dt*float(in1.ntw),itw,irc)
   swe = s[0]; swke = s[1]; swf = s[2]; swb = s[3]
   s[5] = (s[6] - s[5])/wt[0,3]
   print >> iuot, "Energy Conservation = ", float(s[5])
   swe = swe/float(itw)
   print >> iuot, "Average Field Energy <WE> = ", float(swe)
   swke = swke/float(itw)
   print >> iuot, "Average Electron Kinetic Energy <WKE> = ",float(swke)
   print >> iuot, "Ratio <WE>/<WKE>= ", float(swe/swke)
   swf = swf/float(itw)
   print >> iuot, "Average Transverse EField Energy <WF> = ", float(swf)
   print >> iuot, "Ratio <WF>/<WKE>= ", float(swf/swke)
   swb = swb/float(itw)
   print >> iuot, "Average Magnetic Field Energy <WB> = ", float(swb)
   print >> iuot, "Ratio <WB>/<WKE>= ", float(swb/swke)

# velocity diagnostic
if (in1.ntv > 0):
   ts = in1.t0 + in1.dt*float(in1.ntv)
   graf1.displayfvt1(fvtm,' ELECTRON',ts,in1.dt*float(in1.ntv),itv,irc)
   if (irc[0]==1):
      exit(1)
# ions
   if (in1.movion==1):
      graf1.displayfvt1(fvtmi,' ION',ts,in1.dt*float(in1.ntv),itv,irc)
      if (irc[0]==1):
         exit(1)

# display final spectral analysis for ion density
if (in1.movion==1):
   if (in1.ntdi > 0):
      if ((in1.nddi==2) or (in1.nddi==3)):
# performs frequency analysis of accumulated complex time series
         wkdi[:,0] = wmi[numpy.argmax(pkwdi[:,:,0],axis=1)]
         wkdi[:,1] = wmi[numpy.argmax(pkwdi[:,:,1],axis=1)]
# display frequency spectrum
         graf1.dmscaler1(wkdi,'ION DENSITY OMEGA VS MODE',ntime,999,1,
                         in1.modesxdi,cwk,irc)

# display final spectral analysis for potential
if (in1.ntp > 0):
   if ((in1.ndp==2) or (in1.ndp==3)):
# performs frequency analysis of accumulated complex time series
      wk[:,0] = wm[numpy.argmax(pkw[:,:,0],axis=1)]
      wk[:,1] = wm[numpy.argmax(pkw[:,:,1],axis=1)]
# display frequency spectrum
      graf1.dmscaler1(wk,'POTENTIAL OMEGA VS MODE',ntime,999,2,
                      in1.modesxp,cwk,irc)

# display final spectral analysis for ion current density
if (in1.movion==1):
   if (in1.ntji > 0):
      if ((in1.ndji==2) or (in1.ndji==3)):
# performs frequency analysis of accumulated complex time series
         vwkji[0,:,0] = wmi[numpy.argmax(vpkwji[0,:,:,0],axis=1)]
         vwkji[1,:,0] = wmi[numpy.argmax(vpkwji[1,:,:,0],axis=1)]
         vwkji[0,:,1] = wmi[numpy.argmax(vpkwji[0,:,:,1],axis=1)]
         vwkji[1,:,1] = wmi[numpy.argmax(vpkwji[1,:,:,1],axis=1)]
# display frequency spectrum
         graf1.dmvector1(vwkji,'ION CURRENT OMEGA VS MODE',ntime,999,2,
                         2,in1.modesxji,cwk,irc)

# display final spectral analysis for radiative vector potential
if (in1.ntar > 0):
   if ((in1.ndar==2) or (in1.ndar==3)):
# performs frequency analysis of accumulated complex time series
      vwkr[0,:,0] = wm[numpy.argmax(vpkwr[0,:,:,0],axis=1)]
      vwkr[1,:,0] = wm[numpy.argmax(vpkwr[1,:,:,0],axis=1)]
      vwkr[0,:,1] = wm[numpy.argmax(vpkwr[0,:,:,1],axis=1)]
      vwkr[1,:,1] = wm[numpy.argmax(vpkwr[1,:,:,1],axis=1)]
# display frequency spectrum
      graf1.dmvector1(vwkr,'RADIATIVE VPOTENTIAL OMEGA VS MODE',ntime,
                      999,2,2,in1.modesxar,cwk,irc)

# display final spectral analysis for vector potential
if (in1.nta > 0):
   if ((in1.nda==2) or (in1.nda==3)):
# performs frequency analysis of accumulated complex time series
      vwk[0,:,0] = wm[numpy.argmax(vpkw[0,:,:,0],axis=1)]
      vwk[1,:,0] = wm[numpy.argmax(vpkw[1,:,:,0],axis=1)]
      vwk[0,:,1] = wm[numpy.argmax(vpkw[0,:,:,1],axis=1)]
      vwk[1,:,1] = wm[numpy.argmax(vpkw[1,:,:,1],axis=1)]
# display frequency spectrum
      graf1.dmvector1(vwk,'VECTOR POTENTIAL OMEGA VS MODE',ntime,999,2,
                      2,in1.modesxa,cwk,irc)

# display final spectral analysis for transverse efield
if (in1.ntet > 0):
   if ((in1.ndet==2) or (in1.ndet==3)):
# performs frequency analysis of accumulated complex time series
      vwket[0,:,0] = wm[numpy.argmax(vpkwet[0,:,:,0],axis=1)]
      vwket[1,:,0] = wm[numpy.argmax(vpkwet[1,:,:,0],axis=1)]
      vwket[0,:,1] = wm[numpy.argmax(vpkwet[0,:,:,1],axis=1)]
      vwket[1,:,1] = wm[numpy.argmax(vpkwet[1,:,:,1],axis=1)]
# display frequency spectrum
      graf1.dmvector1(vwket,'TRANSVERSE EFIELD OMEGA VS MODE',ntime,999,
                       2,2,in1.modesxet,cwk,irc)

print >> iuot
print >> iuot, "initialization time = ", tinit
print >> iuot, "deposit time = ", tdpost[0]
print >> iuot, "current deposit time = ", tdjpost[0]
tdpost[0] = tdpost[0] + tdjpost[0]
print >> iuot, "total deposit time = ", tdpost[0]
print >> iuot, "guard time = ", tguard[0]
print >> iuot, "solver time = ", tfield[0]
print >> iuot, "fft time = ", tfft[0]
print >> iuot, "push time = ", tpush[0]
print >> iuot, "sort time = ", tsort[0]
tfield[0] = tfield[0] + tguard[0] + tfft[0]
print >> iuot, "total solver time = ", tfield[0]
time = tdpost[0] + tpush[0] + tsort[0]
print >> iuot, "total particle time = ", time
print >> iuot, "total diagnostic time = ", tdiag[0]
ws[0] = time + tfield[0] + tdiag[0]
tloop = tloop - ws[0]
print >> iuot, "total and additional time = ", ws[0], ",", tloop
print >> iuot

ws[0] = 1.0e+09/(float(nloop)*float(np))
print >> iuot, "Push Time (nsec) = ", tpush[0]*ws[0]
print >> iuot, "Deposit Time (nsec) = ", tdpost[0]*ws[0]
print >> iuot, "Sort Time (nsec) = ", tsort[0]*ws[0]
print >> iuot, "Total Particle Time (nsec) = ", time*ws[0]
print >> iuot

# reset parameters for final diagnostic metafile
# electron density diagnostic
if (in1.ntde > 0):
   in1.nderec -= 1
# potential diagnostic
if (in1.ntp > 0):
   in1.nprec -= 1; ceng = affp
# longitudinal efield diagnostic
if (in1.ntel > 0):
   in1.nelrec -= 1; ceng = affp
# radiative vector potential diagnostic
if (in1.ntar > 0):
   in1.narrec -= 1; ceng = affp
# vector potential diagnostic
if (in1.nta > 0):
   in1.narec -= 1; ceng = affp
# transverse efield diagnostic
if (in1.ntet > 0):
   in1.netrec -= 1; ceng = affp
# magnetic field diagnostic
if (in1.ntb > 0):
   in1.nbrec -= 1; ceng = affp
# ion diagnostics
if (in1.movion==1):
# ion density diagnostic
   if (in1.ntdi > 0):
      in1.ndirec -= 1
# ion current diagnostic
if (in1.ntdi > 0):
   in1.njirec -= 1
# write final diagnostic metafile
in1.writnml1(iudm)
print >> iuot, " * * * q.e.d. * * *"
iuot.close()
# close graphics device
graf1.close_graphs()
