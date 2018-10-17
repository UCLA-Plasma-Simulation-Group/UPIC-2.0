#-----------------------------------------------------------------------
# 3D Darwin MPI/OpenMP PIC code
# written by Viktor K. Decyk, UCLA
import sys
import math
import numpy

sys.path.append('./mpbeps3.source')
from libmpush3 import *
from libmbpush3 import *
from libmdpush3 import *
from dtimer import *

int_type = numpy.int32
double_type = numpy.float64
#if (minit1.fprecision()==0):
float_type = numpy.float32
complex_type = numpy.complex64
#else:
#  float_type = numpy.float64
#  complex_type = numpy.complex128
#  print "using double precision"

# idimp = number of particle coordinates = 6
# ipbc = particle boundary condition: 1 = periodic
idimp = 6; ipbc = 1
# idps = number of partition boundaries = 4
# idds = dimensionality of domain decomposition = 2
idps = 4; idds =    2
# wke/wki/we = particle kinetic/electrostatic field energy
wke = numpy.zeros((1),float_type)
wki = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
# wf/wb = transverse electric field/magnetic field
wf = numpy.zeros((1),float_type)
wb = numpy.zeros((1),float_type)
zero = 0.0
# plist = (true,false) = list of particles leaving tiles found in push
#plist = True
plist = False

# declare scalars for standard code
ierr = numpy.zeros((1),int_type)
nt = numpy.empty((1),int_type)
ntime0 = numpy.zeros((1),int_type)
wpmax = numpy.zeros((1),float_type)
wpmin = numpy.zeros((1),float_type)
wpm = numpy.empty((1),float_type)
q2m0 = numpy.empty((1),float_type)
ws = numpy.zeros((1),float_type)
npi = 0.0

# declare scalars for MPI code
nvp = numpy.empty((1),int_type)
idproc = numpy.empty((1),int_type)
npp = numpy.empty((1),int_type)
nppi = numpy.empty((1),int_type)

# declare scalars for OpenMP code
nppmx = numpy.empty((1),int_type)
irc = numpy.zeros((1),int_type)
irc2 = numpy.zeros((2),int_type)

# declare scalars for diagnostics
numtp = numpy.empty((1),int_type)
wk = numpy.zeros((1),float_type)
wef = numpy.zeros((1),float_type)

# default Fortran unit numbers
iuin = 8; iuot = 18; iudm = 19
iude = 10; iup = 11; iuel = 12
iua = 13; iuet = 14; iub = 15
iuje = 21; iufe = 23; iuve = 25; iut = 28; iuse = 29
iudi = 20; iuji = 22; iufi = 24; iuvi = 26; iusi = 30

# declare arrays for standard code:
itot = numpy.empty((2),int_type)
wtot = numpy.empty((7),double_type)

# declare and initialize timing data
tinit = 0.0; tloop = 0.0
itime = numpy.empty((4),numpy.int32)
ltime = numpy.empty((4),numpy.int32)
tdpost = numpy.zeros((1),float_type)
tguard = numpy.zeros((1),float_type)
tfield = numpy.zeros((1),float_type)
tdjpost = numpy.zeros((1),float_type)
tdcjpost = numpy.zeros((1),float_type)
tpush = numpy.zeros((1),float_type)
tsort = numpy.zeros((1),float_type)
tmov = numpy.zeros((1),float_type)
tfmov = numpy.zeros((1),float_type)
tdiag = numpy.zeros((1),float_type)
tfft = numpy.zeros((2),float_type)
dtime = numpy.empty((1),double_type)

# start timing initialization
dtimer(dtime,itime,-1)
  
# nvp = number of MPI ranks
# initialize for distributed memory parallel processing
mpplib3.ppinit2(idproc,nvp)
kstrt = idproc[0] + 1

# in3.nvpp = number of shared memory nodes (0=default)
#if (kstrt==1):
#   in3.nvpp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
omplib.init_omp(in3.nvpp)

# read namelists
if (kstrt==1):
# override default input data
   in3.emf = 2
   in3.readnml3(iuin)
# override input data
   in3.idcode = 3
# create string from idrun
   cdrun = str(in3.idrun)
# text output file
   fname = "output3." + cdrun
   iuot = open(fname,"w")

# broadcast namelists to other nodes
in3.sendnmls3()

# import modules after namelist has been read
import s3
import sb3
import sd3

# increase number of coordinates for particle tag
if (in3.ntt > 0):
   idimp += 1

# initialize scalars for standard code
# np = total number of particles in simulation
npxyz =  float(in3.npx)*float(in3.npy)*float(in3.npz)
npxyzb =  float(in3.npxb)*float(in3.npyb)*float(in3.npzb)
np = npxyz + npxyzb
# npi = total number of ions in simulation
if (in3.movion > 0):
   npxyzi = float(in3.npxi)*float(in3.npyi)*float(in3.npzi)
   npxyzbi = float(in3.npxbi)*float(in3.npybi)*float(in3.npzbi)
   npi = npxyzi + npxyzbi
# nx/ny/nz = number of grid points in x/y/z direction
nx = int(math.pow(2,in3.indx)); ny = int(math.pow(2,in3.indy))
nz = int(math.pow(2,in3.indz))
nxh = int(nx/2); nyh = max(1,int(ny/2)); nzh = max(1,int(nz/2))
nxe = nx + 2; nye = ny + 2; nze = nz + 2
nxeh = nxe/2
nxyzh = max(nx,ny,nz)/2; nxhyz = max(nxh,ny,nz)
# mx1 = number of tiles in x direction
mx1 = int((nx - 1)/in3.mx + 1)
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(in3.tend/in3.dt + .0001); ntime = 0
# mdim = dimension of amu array
mdim = 2*in3.ndim
qbme = in3.qme
affp = float(nx)*float(ny)*float(nz)/np
if (in3.movion==1):
   qbmi = in3.qmi/in3.rmass
   vtxi = in3.vtx/numpy.sqrt(in3.rmass*in3.rtempxi)
   vtyi = in3.vty/numpy.sqrt(in3.rmass*in3.rtempyi)
   vtzi = in3.vtz/numpy.sqrt(in3.rmass*in3.rtempzi)
   vtdxi = in3.vtdx/numpy.sqrt(in3.rmass*in3.rtempdxi)
   vtdyi = in3.vtdy/numpy.sqrt(in3.rmass*in3.rtempdyi)
   vtdzi = in3.vtdz/numpy.sqrt(in3.rmass*in3.rtempdzi)
    
# obtain 2D partition (nvpy,nvpz) from nvp:
# nvpy/nvpz = number of processors in y/z
minit3.mpfcomp3(nvp,nx,ny,nz,in3.nvpy,in3.nvpz,ierr)
if (ierr[0] != 0):
   if (kstrt==1):
      print "mpfcomp3 error: nvp,nvpy,nvpz=",nvp[0],in3.nvpy,in3.nvpz
   mpplib3.ppexit(); exit(1)

# initialize data for MPI code
edges = numpy.empty((idps),float_type,'F')
nyzp = numpy.empty((idds),int_type,'F')
noff = numpy.empty((idds),int_type,'F')
nypmx = numpy.empty((idds),int_type,'F')
nzpmx = numpy.empty((idds),int_type,'F')
nypmn = numpy.empty((idds),int_type,'F')
nzpmn = numpy.empty((idds),int_type,'F')
# calculate partition variables:
# edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn
# edges[0:1] = lower:upper boundary of particle partition in y
# edges[2:3] = back:front boundary of particle partition in z
# nyzp[0:1] = number of primary (complete) gridpoints in y/z
# noff[0:1] = lowermost global gridpoint in y/z in particle partition
# nypmx = maximum size of particle partition in y, including guard cells
# nzpmx = maximum size of particle partition in z, including guard cells
# nypmn = minimum value of nyzp(1)
# nzpmn = minimum value of nyzp(2)
# find new 2d partition for uniform density distribution
#     call mpdcomp3(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,&
#    &nvpy,nvpz)
# find new 2d partition from initial analytic distribution function
minit3.mpfedges3(edges,nyzp,noff,in3.ampdy,in3.scaledy,in3.shiftdy,
                 in3.ampdz,in3.scaledz,in3.shiftdz,nypmx,nzpmx,nypmn,
                 nzpmn,ny,nz,kstrt,in3.nvpy,in3.nvpz,ipbc,in3.ndprof,
                 ierr)
if (ierr[0] != 0):
   mpplib3.ppexit(); exit(1)

# check for unimplemented features
if ((plist) and (ipbc != 1)):
   print "ipbc != 1 and list = True not yet supported"
   plist = False
   print "plist reset to False"

# initialize additional scalars for MPI code
# kyp = number of complex grids in each field partition in y direction
kyp = int((ny - 1)/in3.nvpy + 1)
# kzp = number of complex grids in each field partition in z direction
kzp = int((nz - 1)/in3.nvpz + 1)
# kxyp = number of complex grids in each field partition in x direction
# in transposed data
kxyp = int((nxh - 1)/in3.nvpy + 1)
# kyzp = number of complex grids in each field partition in y direction,
# in transposed data
kyzp = int((ny - 1)/in3.nvpz + 1)
# npmax/npimax = maximum number of electrons/ions in each partition
npmax = int((np/nvp)*1.25); npimax = int((npi/nvp)*1.25)
maxnp = max(npmax,npimax)
# myp1/mzp1 = number of tiles in y/z direction
myp1 = int((nyzp[0] - 1)/in3.my + 1); mzp1 = int((nyzp[1] - 1)/in3.mz + 1)
# mxzyp1 = mx1*max(max(mzp1),max(myp1))
mxzyp1 = mx1*max(int((nzpmx[0]-2)/in3.mz+1),int((nypmx[0]-2)/in3.my+1))
mxyzp1 = mx1*myp1*mzp1
# mterf = number of shifts required by field manager in y/z (0=search)
mterf = numpy.zeros((2),int_type)
# ntmax = size of iholep buffer for particles leaving node
ntmax = int(0.2*npmax)

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,maxnp),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nypmx[0],nzpmx[0]),float_type,'F')
# qi = ion charge density with guard cells
qi = numpy.empty((nxe,nypmx[0],nzpmx[0]),float_type,'F')
# fxyze = smoothed longitudinal electric field with guard cells
fxyze = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# cue = electron current density with guard cells
cue = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# dcu = electron acceleration density with guard cells
dcu = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# cus = smoothed transverse electric field with guard cells
cus = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# amu = electron momentum flux with guard cells
amu = numpy.empty((mdim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# exyze = smoothed total electric field with guard cells
exyze = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# bxyze = smoothed /magnetic field with guard cells
bxyze = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# qt = scalar charge density field array in fourier space
qt = numpy.empty((nze,kxyp,kyzp),complex_type,'F')
# exyzt = vector transverse electric field in fourier space
exyzt = numpy.empty((in3.ndim,nze,kxyp,kyzp),complex_type,'F')
# fxyzt = vector longitudinal electric field in fourier space
fxyzt = numpy.empty((in3.ndim,nze,kxyp,kyzp),complex_type,'F')
# cut = vector current density field arrays in fourier space
cut = numpy.empty((in3.ndim,nze,kxyp,kyzp),complex_type,'F')
# dcut = vector acceleration density in fourier space
dcut = numpy.empty((in3.ndim,nze,kxyp,kyzp),complex_type,'F')
# amut = tensor momentum flux in fourier space
amut = numpy.empty((in3.ndim,nze,kxyp,kyzp),complex_type,'F')
# bxyzt = vector magnetic field array in fourier space
bxyzt = numpy.empty((in3.ndim,nze,kxyp,kyzp),complex_type,'F')
# ffc, ffe = form factor arrays for poisson solvers
ffc = numpy.empty((nzh,kxyp,kyzp),complex_type,'F')
ffe = numpy.empty((nzh,kxyp,kyzp),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhyz),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyzh),complex_type,'F')
# kpic = number of electrons in each tile
kpic = numpy.empty((mxyzp1),int_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((26,mxyzp1),int_type,'F')
# iholep = location of hole left in linear particle arrays
iholep = numpy.empty((ntmax+1,2),int_type,'F')
if (in3.movion==1):
# kipic = number of ions in each tile
   kipic = numpy.empty((mxyzp1),int_type,'F')
# cui = ion current density with guard cells
   cui = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# dcui = electron/ion acceleration density with guard cells
   dcui = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# amui = electron/ion momentum flux with guard cells
   amui = numpy.empty((mdim,nxe,nypmx[0],nzpmx[0]),float_type,'F')
# define dummy arrays for restart
elif (in3.movion==0):
   kipic = numpy.empty((0),int_type,'F')
   pparti = numpy.empty((0,0,0),float_type,'F')

# prepare fft tables
mfft3.mpfft3_init(mixup,sct,in3.indx,in3.indy,in3.indz)
# calculate form factor: ffc
mfield3.mppois3_init(ffc,in3.ax,in3.ay,in3.az,affp,nx,ny,nz,kstrt,
                     in3.nvpy,in3.nvpz)
if (in3.nextrand > 0):
   minit3.mnextran3(nextrand,in3.ndim,npmax+npimax)

# open restart files
if (kstrt==1):
   iur, iur0 = s3.open_restart3(cdrun)

# new start
if (in3.nustrt==1):
# initialize electrons
   nps = 1
   npp[0] = 0
# background electrons
   if (npxyz > 0.0):
# calculates initial electron co-ordinates with various density profiles
      minit3.mpfdistr3(part,npp,in3.ampdx,in3.scaledx,in3.shiftdx,
                       in3.ampdy,in3.scaledy,in3.shiftdy,in3.ampdz,
                       in3.scaledz,in3.shiftdz,in3.npx,in3.npy,in3.npz,
                       nx,ny,nz,kstrt,in3.nvpy,in3.nvpz,ipbc,in3.ndprof,
                       ierr)
# initialize electron velocities or momenta
      if (ierr[0]==0):
         minit3.wmpvdistr3(part,nps,npp,in3.vtx,in3.vty,in3.vtz,in3.vx0,
                           in3.vy0,in3.vz0,in3.ci,in3.npx,in3.npy,
                           in3.npz,kstrt,in3.nvpy,in3.nvpz,
                           in3.relativity,ierr)
# check for background electron initialization error
      if (ierr[0] != 0):
         mpplib3.ppexit(); exit(1)
# beam electrons
   if (npxyzb > 0.0):
      nps = npp[0] + 1
# calculates initial electron co-ordinates with various density profiles
      minit3.mpfdistr3(part,npp,in3.ampdx,in3.scaledx,in3.shiftdx,
                       in3.ampdy,in3.scaledy,in3.shiftdy,in3.ampdz,
                       in3.scaledz,in3.shiftdz,in3.npxb,in3.npyb,
                       in3.npzb,nx,ny,nz,kstrt,in3.nvpy,in3.nvpz,ipbc,
                       in3.ndprof,ierr)
# initialize electron velocities or momenta
      if (ierr[0]==0):
         minit3.wmpvdistr3(part,nps,npp,in3.vtdx,in3.vtdy,in3.vtdz,
                           in3.vdx,in3.vdy,in3.vdz,in3.ci,in3.npxb,
                           in3.npyb,in3.npzb,kstrt,in3.nvpy,in3.nvpz,
                           in3.relativity,ierr)
# check for beam electron initialization error
      if (ierr[0] != 0):
         mpplib3.ppexit(); exit(1)

# check if any electrons are in the wrong node
   minit3.mpfholes3(part,edges,npp,iholep,1)
   if (iholep[0,0] < 0):
      ntmax = -iholep[0,0]
      ntmax = int(1.5*ntmax)
      if (kstrt==1):
         print "reallocating electron iholep: ntmax=", ntmax
      iholep = numpy.empty((ntmax+1,2),int_type,'F') 
      minit3.mpfholes3(part,edges,npp,iholep,1)
      if (iholep[0,0] < 0):
         if (kstrt==1):
            print "iholep overflow: ntmax=", ntmax
         mpplib3.ppexit(); exit(1)
# move electrons to correct node
      mppmod3.ipmove3(part,edges,npp,iholep,ny,nz,tinit,kstrt,in3.nvpy,
                     in3.nvpz,1,ierr)
      if (ierr[0] != 0):
         mpplib3.ppexit(); exit(1)

# find number of electrons in each of mx, my, mz tiles:
# updates kpic, nppmx
   minit3.mpdblkp3(part,kpic,npp,noff,nppmx,in3.mx,in3.my,in3.mz,mx1,
                   myp1,irc)

# allocate vector electron data
   nppmx0 = int((1.0 + in3.xtras)*nppmx)
   ntmaxp = int(in3.xtras*nppmx)
   npbmx = int(in3.xtras*nppmx)
   nbmaxp = int(0.25*mxzyp1*npbmx)
# ppart = tiled electron particle array
   ppart = numpy.empty((idimp,nppmx0,mxyzp1),float_type,'F')
# ihole = location/destination of each particle departing tile
   ihole = numpy.empty((2,ntmaxp+1,mxyzp1),int_type,'F')
# copy ordered electron data for OpenMP: updates ppart and kpic
   mpush3.mpmovin3p(part,ppart,kpic,npp,noff,in3.mx,in3.my,in3.mz,mx1,
                    myp1,irc)

# sanity check for electrons
   mpush3.mpcheck3(ppart,kpic,noff,nyzp,nx,in3.mx,in3.my,in3.mz,mx1,
                   myp1,irc)

# initialize background charge density: updates qi
   if (in3.movion==0):
      mpush3.mpset_pszero3(qi,tinit,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
      qmi = -in3.qme
      mpush3.mppost3(ppart,qi,kpic,noff,qmi,tinit,in3.mx,in3.my,in3.mz,
                     mx1,myp1,in3.popt)
      ompplib3.wmpaguard3(qi,nyzp,tinit,nx,kstrt,in3.nvpy,in3.nvpz)

# initialize ions
   if (in3.movion==1):
      mpush3.mpset_pvzero3(cui,tinit,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
      nps = 1
      nppi[0] = 0
# background ions
      if (npxyzi > 0.0):
# calculates initial ion co-ordinates with various density profiles
         minit3.mpfdistr3(part,nppi,in3.ampdxi,in3.scaledxi,
                          in3.shiftdxi,in3.ampdyi,in3.scaledyi,
                          in3.shiftdyi,in3.ampdzi,in3.scaledzi,
                          in3.shiftdzi,in3.npxi,in3.npyi,in3.npzi,nx,ny,
                          nz,kstrt,in3.nvpy,in3.nvpz,ipbc,in3.ndprofi,
                          ierr)
# initialize ion velocities or momenta
         if (ierr[0]==0):
            minit3.wmpvdistr3(part,nps,nppi,vtxi,vtyi,vtzi,in3.vxi0,
                              in3.vyi0,in3.vzi0,in3.ci,in3.npxi,
                              in3.npyi,in3.npzi,kstrt,in3.nvpy,in3.nvpz,
                              in3.relativity,ierr)
# check for background ion initialization error
         if (ierr[0] != 0):
            mpplib3.ppexit(); exit(1)
# beam ions
      if (npxyzbi > 0.0):
         nps = nppi[0] + 1
# calculates initial ion co-ordinates with various density profiles
         minit3.mpfdistr3(part,nppi,in3.ampdxi,in3.scaledxi,
                          in3.shiftdxi,in3.ampdyi,in3.scaledyi,
                          in3.shiftdyi,in3.ampdzi,in3.scaledzi,
                          in3.shiftdzi,in3.npxbi,in3.npybi,in3.npzbi,nx,
                          ny,nz,kstrt,in3.nvpy,in3.nvpz,ipbc,
                          in3.ndprofi,ierr)
# initialize ion velocities or momenta
         if (ierr[0]==0):
            minit3.wmpvdistr3(part,nps,nppi,vtdxi,vtdyi,vtdzi,in3.vdxi,
                              in3.vdyi,in3.vdzi,in3.ci,in3.npxbi,
                              in3.npybi,in3.npzbi,kstrt,in3.nvpy,
                              in3.nvpz,in3.relativity,ierr)
# check for beam ion initialization error
         if (ierr[0] != 0):
            mpplib3.ppexit(); exit(1)

# check if any ions are in the wrong node
      minit3.mpfholes3(part,edges,nppi,iholep,1)
# iholep overflow
      if (iholep[0,0] < 0):
         ntmax = -iholep[0,0]
         ntmax = int(1.5*ntmax)
         if (kstrt==1):
            print "reallocating ion iholep: ntmax=", ntmax
         iholep = numpy.empty((ntmax+1,2),int_type,'F') 
         minit3.mpfholes3(part,edges,nppi,iholep,1)
         if (iholep[0,0] < 0):
            if (kstrt==1):
               print "iholep overflow: ntmax=", ntmax
            mpplib3.ppexit(); exit(1)
# move ions to correct node
         mppmod3.ipmove3(part,edges,nppi,iholep,ny,nz,tinit,kstrt,
                         in3.nvpy,in3.nvpz,1,ierr)
         if (ierr[0] != 0):
            mpplib3.ppexit(); exit(1)

# find number of ions in each of mx, my, mz tiles:
# updates kipic, nppmx
      minit3.mpdblkp3(part,kipic,nppi,noff,nppmx,in3.mx,in3.my,in3.mz,
                      mx1,myp1,irc)

# allocate vector ion data
      nppmx1 = int((1.0 + in3.xtras)*nppmx)
# pparti = tiled ion particle array
      pparti = numpy.empty((idimp,nppmx1,mxyzp1),float_type,'F')
      if ("ihole" not in globals()):
         ntmaxp = int(in3.xtras*nppmx)
         npbmx = int(in3.xtras*nppmx)
         nbmaxp = int(0.25*mxzyp1*npbmx)
         ihole = numpy.empty((2,ntmaxp+1,mxyzp1),int_type,'F')
# copy ordered ion data for OpenMP: updates pparti and kipic
      mpush3.mpmovin3p(part,pparti,kipic,nppi,noff,in3.mx,in3.my,in3.mz,
                       mx1,myp1,irc)

# sanity check for ions
      mpush3.mpcheck3(pparti,kipic,noff,nyzp,nx,in3.mx,in3.my,in3.mz,
                      mx1,myp1,irc)

# find maximum and minimum initial electron density
   mpush3.mpset_pszero3(qe,tinit,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
   mpush3.mppost3(ppart,qe,kpic,noff,in3.qme,tinit,in3.mx,in3.my,in3.mz,
                  mx1,myp1,in3.popt)
   ompplib3.wmpaguard3(qe,nyzp,tinit,nx,kstrt,in3.nvpy,in3.nvpz)
   if (in3.movion==1):
      mpush3.mpset_pszero3(qi,tinit,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
      mpush3.mppost3(pparti,qi,kipic,noff,in3.qmi,tinit,in3.mx,in3.my,
                     in3.mz,mx1,myp1,in3.popt)
      ompplib3.wmpaguard3(qi,nyzp,tinit,nx,kstrt,in3.nvpy,in3.nvpz)
      mdpush3.mpfwptminx3(qe,qi,nyzp,qbme,qbmi,wpmax,wpmin,nx)
   else:
      mdpush3.mpfwpminx3(qe,nyzp,qbme,wpmax,wpmin,nx)
   wtot[0] = wpmax[0]
   wtot[1] = -wpmin[0]
   mppmod3.mpdmax(wtot[:2],tinit)
   wpmax[0] = wtot[0]
   wpmin[0] = -wtot[1]
   wpm[0] = 0.5*(wpmax[0] + wpmin[0])*affp
# accelerate convergence: update wpm
   if (wpm[0] <= 10.0):
      wpm[0] = 0.75*wpm[0]
   if (kstrt==1):
      print >> iuot, "wpm=", wpm[0]
   q2m0[0] = wpm[0]/affp

# initialize darwin electric field
   mpush3.mpset_pvzero3(cus,tinit,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)

# restart to continue a run which was interrupted
elif (in3.nustrt==2):
   print "nustrt = 2 not yet supported"
   mpplib3.ppexit(); exit(1)
#  nstart = ntime + 1
# start a new run with data from a previous run
elif (in3.nustrt==0):
# read in basic restart file for electrostatic code
   if (in3.movion==1):
      mpush3.mpset_pvzero3(cui,tinit,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
# read first part of data:
# updates ntime, ntime0, part, npp, kpic, nppmx, ierr
   s3.bread_restart3a(part,kpic,noff,iur0,nt,ntime0,npp,nppmx,mx1,myp1,
                      ierr)
   ntime = nt[0]
   if (ierr[0] != 0):
      mpplib3.ppexit(); exit(1)
# allocate vector electron data
   nppmx0 = int((1.0 + in3.xtras)*nppmx)
   ntmaxp = int(in3.xtras*nppmx)
   npbmx = int(in3.xtras*nppmx)
   nbmaxp = int(0.25*mxzyp1*npbmx)
   ppart = numpy.empty((idimp,nppmx0,mxyzp1),float_type,'F')
   ihole = numpy.empty((2,ntmaxp+1,mxyzp1),int_type,'F')
# read second part of data:
# updates ppart, kpic, part, nppi, kipic, nppmx, ierr
   s3.bread_restart3b(part,ppart,kpic,kipic,noff,nyzp,iur0,npp,nppi,
                      nppmx,nx,mx1,myp1,ierr)
   if (ierr[0] != 0):
      mpplib3.ppexit(); exit(1)
# allocate vector ion data
   if (in3.movion==1):
      nppmx1 = int((1.0 + in3.xtras)*nppmx)
      pparti = numpy.empty((idimp,nppmx1,mxyzp1),float_type,'F')
      if ("ihole" not in globals()):
         ntmaxp = int(in3.xtras*nppmx)
         npbmx = int(in3.xtras*nppmx)
         nbmaxp = int(0.25*mxzyp1*npbmx)
         ihole = numpy.empty((2,ntmaxp+1,mxyzp1),int_type,'F')
# read third part of data: updates pparti, kipic, qi, ierr
   s3.bread_restart3c(part,pparti,kipic,noff,nyzp,qi,iur0,nt,ntime0,
                      nppi,nx,mx1,myp1,ierr)
   ntime = nt[0]
   if (ierr[0] != 0):
      mpplib3.ppexit(); exit(1)
# read in basic restart file for darwin code: updates cus, wpm, q2m0
   sd3.bread_drestart3(cus,wpm,q2m0,iur0)

# calculate form factor: ffe
mfield3.mpepois3_init(ffe,in3.ax,in3.ay,in3.az,affp,wpm,in3.ci,nx,ny,nz,
                      kstrt,in3.nvpy,in3.nvpz)

# initialize longitudinal electric field
mpush3.mpset_pvzero3(fxyze,tinit,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)

# set magnitude of external magnetic field
omt = numpy.sqrt(in3.omx*in3.omx + in3.omy*in3.omy + in3.omz*in3.omz)

# js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
ks = int((kstrt - 1)/in3.nvpy)
js = int(kstrt - in3.nvpy*ks - 1)
# kyps/kzps = actual grids used in field partitions in y/z
kyps = int(min(kyp,max(0,ny-kyp*js)))
kzps = int(min(kzp,max(0,nz-kzp*ks)))
# kxyps/kyzps = actual grids in each transposed field partition in x/y
kxyps = int(min(kxyp,max(0,nxh-kxyp*js)))
kyzps = int(min(kyzp,max(0,ny-kyzp*ks)))

# allocate diagnostic arrays
# reverse simulation at end back to start
if (in3.treverse==1):
   nloop = 2*nloop

# energy time history
if (in3.ntw > 0):
   mtw = int((nloop - 1)/in3.ntw + 1); itw = 0
# wt = energy time history array
   wt = numpy.zeros((mtw,7),float_type,'F')
# s = scratch array for energies
   s = numpy.zeros((7),double_type,'F')

# allocate scratch arrays for scalar fields
if ((in3.ntde > 0) or (in3.ntp > 0) or (in3.ntdi > 0)):
   sfieldc = numpy.empty((nze,kxyp,kyzp),complex_type,'F')
   sfield = numpy.empty((nxe,nypmx[0],nzpmx[0]),float_type,'F')

# allocate scratch arrays for vector fields
if ((in3.ntel>0) or (in3.ntje>0) or (in3.nta>0) or (in3.ntet>0)
or (in3.ntb>0) or (in3.ntji>0)):
   vfieldc = numpy.empty((in3.ndim,nze,kxyp,kyzp),complex_type,'F')
   vfield = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,'F')

# initialize electron density diagnostic
if (in3.ntde > 0):
   in3.modesxde = int(min(in3.modesxde,nxh+1))
   in3.modesyde = int(min(in3.modesyde,nyh+1))
   in3.modeszde = int(min(in3.modeszde,nzh+1))
   modesz2de = int(min(2*in3.modeszde-1,nz))
   modesxpd = int(min(in3.modesxde,kxyp))
   if (in3.modesxde==(nxh+1)):
      modesxpd += 1
# denet = store selected fourier modes for electron density
   denet = numpy.empty((modesz2de,modesxpd,kyzp),complex_type,'F')
# open file for real data: updates nderec and possibly iude
   if (kstrt==1):
      fdename = "dener3." + cdrun
      in3.fdename[:] = fdename
      if (in3.nderec==0):
         mdiag3.dafopen3(sfield,nx,kyp,kzp,iude,in3.nderec,fdename)

# initialize ion density diagnostic
if (in3.movion==1):
   if (in3.ntdi > 0):
      in3.modesxdi = int(min(in3.modesxdi,nxh+1))
      in3.modesydi = int(min(in3.modesydi,nyh+1))
      in3.modeszdi = int(min(in3.modeszdi,nzh+1))
      modesz2di = int(min(2*in3.modeszdi-1,nz))
      modesxpd = int(min(in3.modesxdi,kxyp))
      if (in3.modesxdi==(nxh+1)):
         modesxpd += 1
# denit = store selected fourier modes for ion density
      denit = numpy.empty((modesz2di,modesxpd,kyzp),complex_type,'F')
# open file for real data: updates ndirec and possibly iudi
      if (kstrt==1):
         fdiname = "denir3." + cdrun
         in3.fdiname[:] = fdiname
         if (in3.ndirec==0):
            mdiag3.dafopen3(sfield,nx,kyp,kzp,iudi,in3.ndirec,fdiname)

# initialize potential diagnostic
if (in3.ntp > 0):
   in3.modesxp = int(min(in3.modesxp,nxh+1))
   in3.modesyp = int(min(in3.modesyp,nyh+1))
   in3.modeszp = int(min(in3.modeszp,nzh+1))
   modesz2p = int(min(2*in3.modeszp-1,nz))
   modesxpd = int(min(in3.modesxp,kxyp))
   if (in3.modesxp==(nxh+1)):
      modesxpd += 1
# pott = store selected fourier modes for potential
   pott = numpy.empty((modesz2p,modesxpd,kyzp),complex_type,'F')
# open file for real data: updates nprec and possibly iup
   if (kstrt==1):
      fpname = "potr3." + cdrun
      in3.fpname[:] = fpname
      if (in3.nprec==0):
         mdiag3.dafopen3(sfield,nx,kyp,kzp,iup,in3.nprec,fpname)

# initialize longitudinal efield diagnostic
if (in3.ntel > 0):
   in3.modesxel = int(min(in3.modesxel,nxh+1))
   in3.modesyel = int(min(in3.modesyel,nyh+1))
   in3.modeszel = int(min(in3.modeszel,nzh+1))
   modesz2el = int(min(2*in3.modeszel-1,nz))
   modesxpd = int(min(in3.modesxel,kxyp))
   if (in3.modesxel==(nxh+1)):
      modesxpd += 1
# elt = store selected fourier modes for longitudinal efield
   elt = numpy.empty((in3.ndim,modesz2el,modesxpd,kyzp),complex_type,
                     'F')
# open file for real data: updates nelrec and possibly iuel
   if (kstrt==1):
      felname = "elr3." + cdrun
      in3.felname[:] = felname
      if (in3.nelrec==0):
         mdiag3.dafopenv3(vfield,nx,kyp,kzp,iuel,in3.nelrec,felname)

# initialize electron current density diagnostic
if (in3.ntje > 0):
   in3.modesxje = int(min(in3.modesxje,nxh+1))
   in3.modesyje = int(min(in3.modesyje,nyh+1))
   in3.modeszje = int(min(in3.modeszje,nzh+1))
   modesz2je = int(min(2*in3.modeszje-1,nz))
   modesxpd = int(min(in3.modesxje,kxyp))
   if (in3.modesxje==(nxh+1)):
      modesxpd += 1
# oldcue = previous current density
   oldcue = numpy.empty((in3.ndim,nxe,nypmx[0],nzpmx[0]),float_type,
                           'F')
# curet = store selected fourier modes for electron current density
   curet = numpy.empty((in3.ndim,modesz2je,modesxpd,kyzp),complex_type,
                       'F')
# open file for real data: updates njerec and possibly iuje
   if (kstrt==1):
      fjename = "curer3." + cdrun
      in3.fjename[:] = fjename
      if (in3.njerec==0):
         mdiag3.dafopenv3(vfield,nx,kyp,kzp,iuje,in3.njerec,fjename)

# initialize ion current density diagnostic
if (in3.movion==1):
   if (in3.ntji > 0):
      in3.modesxji = int(min(in3.modesxji,nxh+1))
      in3.modesyji = int(min(in3.modesyji,nyh+1))
      in3.modeszji = int(min(in3.modeszji,nzh+1))
      modesz2ji = int(min(2*in3.modeszji-1,nz))
      modesxpd = int(min(in3.modesxji,kxyp))
      if (in3.modesxji==(nxh+1)):
         modesxpd += 1
# curit = store selected fourier modes for ion current density
      curit = numpy.empty((in3.ndim,modesz2ji,modesxpd,kyzp),
                          complex_type,'F')
# open file for real data: updates njirec and possibly iuji
      if (kstrt==1):
         fjiname = "curir3." + cdrun
         in3.fjiname[:] = fjiname
         if (in3.njirec==0):
            mdiag3.dafopenv3(vfield,nx,kyp,kzp,iuji,in3.njirec,fjiname)

# initialize vector potential diagnostic
if (in3.nta > 0):
   in3.modesxa = int(min(in3.modesxa,nxh+1))
   in3.modesya = int(min(in3.modesya,nyh+1))
   in3.modesza = int(min(in3.modesza,nzh+1))
   modesz2a = int(min(2*in3.modesza-1,nz))
   modesxpd = int(min(in3.modesxa,kxyp))
   if (in3.modesxa==(nxh+1)):
       modesxpd += 1
# vpott = store selected fourier modes for vector potential
   vpott = numpy.empty((in3.ndim,modesz2a,modesxpd,kyzp),complex_type,
                       'F')
# open file for real data: updates narec and possibly iua
   if (kstrt==1):
      faname = "vpotr3."+ cdrun
      in3.faname[:] = faname
      if (in3.narec==0):
         mdiag3.dafopenv3(vfield,nx,kyp,kzp,iua,in3.narec,faname)

# initialize transverse efield diagnostic
if (in3.ntet > 0):
   in3.modesxet = int(min(in3.modesxet,nxh+1))
   in3.modesyet = int(min(in3.modesyet,nyh+1))
   in3.modeszet = int(min(in3.modeszet,nzh+1))
   modesz2et = int(min(2*in3.modeszet-1,nz))
   modesxpd = int(min(in3.modesxet,kxyp))
   if (in3.modesxet==(nxh+1)):
      modesxpd += 1
# ett = store selected fourier modes for transverse efield
   ett = numpy.empty((in3.ndim,modesz2et,modesxpd,kyzp),complex_type,
                     'F')
# open file for real data: updates naetrec and possibly iuet
   if (kstrt==1):
      fetname = "etr3." + cdrun
      in3.fetname[:] = fetname
      if (in3.netrec==0):
         mdiag3.dafopenv3(vfield,nx,kyp,kzp,iuet,in3.netrec,fetname)

# initialize magnetic field diagnostic
if (in3.ntb > 0):
   in3.modesxb = int(min(in3.modesxb,nxh+1))
   in3.modesyb = int(min(in3.modesyb,nyh+1))
   in3.modeszb = int(min(in3.modeszb,nzh+1))
   modesz2b = int(min(2*in3.modeszb-1,nz))
   modesxpd = int(min(in3.modesxb,kxyp))
   if (in3.modesxb==(nxh+1)):
      modesxpd += 1
# bt = store selected fourier modes for magnetic field
   bt = numpy.empty((in3.ndim,modesz2b,modesxpd,kyzp),complex_type,'F')
# open file for real data: updates nbrec and possibly iub
   if (kstrt==1):
      fbname = "br3." + cdrun
      in3.fbname[:] = fbname
      if (in3.nbrec==0):
         mdiag3.dafopenv3(vfield,nx,kyp,kzp,iub,in3.nbrec,fbname)

# initialize fluid moments diagnostic
if (in3.ntfm > 0):
   in3.nprd = 0
   if (in3.npro==1):
      in3.nprd = 1
   elif (in3.npro==2):
      in3.nprd = 4
   elif (in3.npro==3):
      in3.nprd = 10
   elif (in3.npro==4):
      in3.nprd = 14
# electron moments
   if ((in3.ndfm==1) or (in3.ndfm==3)):
# fmse = electron fluid moments
      fmse = numpy.empty((in3.nprd,nxe,nypmx[0],nzpmx[0]),float_type,
                         'F')
# open file for real data: updates nferec and possibly iufe
      if (kstrt==1):
         ffename = "fmer3." + cdrun
         in3.ffename[:] = ffename
         if (in3.nferec==0):
            mdiag3.dafopenv3(fmse,nx,kyp,kzp,iufe,in3.nferec,ffename)
# ion moments
   if (in3.movion==1):
      if ((in3.ndfm==2) or (in3.ndfm==3)):
# fmsi = ion fluid moments
         fmsi = numpy.empty((in3.nprd,nxe,nypmx[0],nzpmx[0]),float_type,
                            'F')
# open file for real data: updates nfirec and possibly iufi
         if (kstrt==1):
            ffiname = "fmir3." + cdrun
            in3.ffiname[:] = ffiname
            if (in3.nfirec==0):
               mdiag3.dafopenv3(fmsi,nx,kyp,kzp,iufi,in3.nfirec,ffiname)

# initialize velocity-space diagnostic
if (in3.ntv > 0):
   in3.nfvd = 0; in3.nfed = 0
   if ((in3.nvft==1) or (in3.nvft==3)):
      in3.nfvd = in3.ndim
   if ((in3.nvft==4) or (in3.nvft==5)):
      in3.nfvd = 2
   if ((in3.nvft==2) or (in3.nvft==3) or (in3.nvft==5)):
      in3.nfed = 1
   nmv21 = 2*in3.nmv + 1
   mtv = int((nloop - 1)/in3.ntv) + 1; itv = 0
   eci = in3.ci
   if (in3.relativity==0):
      eci = 0.0
   wk = numpy.zeros((1),float_type)
# electron velocity diagnostic
   if ((in3.ndv==1) or (in3.ndv==3)):
# estimate maximum velocity or momentum
      ws[0] = 0.0
      if (npxyz > 0.0):
         ws[0] = 4.0*in3.vtx+abs(in3.vx0)
         ws[0] = max(ws[0],4.0*in3.vty+abs(in3.vy0))
         ws[0] = max(ws[0],4.0*in3.vtz+abs(in3.vz0))
      if (npxyzb > 0.0):
         ws[0] = max(ws[0],4.0*in3.vtdx+abs(in3.vdx))
         ws[0] = max(ws[0],4.0*in3.vtdy+abs(in3.vdy))
         ws[0] = max(ws[0],4.0*in3.vtdz+abs(in3.vdz))
# fvm = electron vdrift, vth, entropy for global distribution
      fvm = numpy.zeros((in3.ndim,3),float_type,'F')
# sfv = electron/ion velocity distribution functions in tile
      sfv = numpy.empty((nmv21+1,in3.ndim,mxyzp1),float_type,'F')
# fv = global electron velocity distribution functions
      fv = numpy.empty((nmv21+1,in3.nfvd),float_type,'F')
# fe = global electron energy distribution functions
      fe = numpy.empty((nmv21+1,in3.nfed),float_type,'F')
# open file for electron velocity data: updates nverec and possibly iuve
      if (kstrt==1):
         fvename = "fve3." + cdrun
         in3.fvename[:] = fvename
         if (in3.nverec==0):
            mdiag3.dafopenfv3(fvm,fv,fe,wk,iuve,in3.nverec,fvename)
# cartesian distribution
      if ((in3.nvft==1) or (in3.nvft==3)):
# fvtm = time history of electron vdrift, vth, and entropy for global
# distribution
         fvtm = numpy.zeros((mtv,in3.ndim,3),float_type,'F')
# set velocity or momentum scale
         fv[nmv21,:] = 2.0*ws[0]
# cylindrical distribution
      if ((in3.nvft==4) or (in3.nvft==5)):
# set velocity or momentum scale
         fv[nmv21,:] = 2.0*ws[0]
# energy distribution
      if ((in3.nvft==2) or (in3.nvft==3) or (in3.nvft==5)):
# set energy scale for electrons
         ws[0] = ws[0]*ws[0]
         fe[nmv21,:] = ws[0]/(1.0 + numpy.sqrt(1.0 + ws[0]*eci*eci))
# ion velocity diagnostic
   if (in3.movion==1):
      if ((in3.ndv==2) or (in3.ndv==3)):
# fvmi = ion vdrift, vth, entropy for global distribution
         fvmi = numpy.zeros((in3.ndim,3),float_type,'F')
# fvi = global ion velocity distribution functions
         fvi = numpy.empty((nmv21+1,in3.nfvd),float_type,'F')
# fei = global ion energy distribution functions
         fei = numpy.empty((nmv21+1,in3.nfed),float_type,'F')
# estimate maximum ion velocity or momentum
         ws[0] = 0.0
         if (npxyzi > 0.0):
            ws[0] = 4.0*vtxi+abs(in3.vxi0)
            ws[0] = max(ws[0],4.0*vtyi+abs(in3.vyi0))
            ws[0] = max(ws[0],4.0*vtzi+abs(in3.vzi0))
         if (npxyzbi > 0.0):
            ws[0] = max(ws[0],4.0*vtdxi+abs(in3.vdxi))
            ws[0] = max(ws[0],4.0*vtdyi+abs(in3.vdyi))
# open file for ion velocity data: updates nvirec and possibly iuvi
         if (kstrt==1):
            fviname = "fvi3." + cdrun
            in3.fviname[:] = fviname
            if (in3.nvirec==0):
               mdiag3.dafopenfv3(fvmi,fvi,fei,wk,iuvi,in3.nvirec,
                                 fviname)
# cartesian distribution
         if ((in3.nvft==1) or (in3.nvft==3)):
# fvtmi = time history of ion vdrift, vth, and entropy for global
# distribution
            fvtmi = numpy.zeros((mtv,in3.ndim,3),float_type,'F')
# set velocity or momentum scale
            fvi[nmv21,:] = 2.0*ws[0]
# cylindrical distribution
         if ((in3.nvft==4) or (in3.nvft==5)):
# set velocity or momentum scale
            fvi[nmv21,:] = 2.0*ws[0]
# energy distribution
         if ((in3.nvft==2) or (in3.nvft==3) or (in3.nvft==5)):
# set energy scale for ions
            ws[0] = ws[0]*ws[0]
            fei[nmv21,0] = ws[0]/(1.0 + numpy.sqrt(1.0 + ws[0]*eci*eci))

# initialize trajectory diagnostic
if (in3.ntt > 0):
   if ((in3.ndt==2) and (in3.movion==0)):
      in3.ndt = 0
   if ((in3.ndt==1) or (in3.ndt==2)):
# iprobt = scratch array
      iprobt = numpy.empty((in3.nprobt),numpy.int32)
# tedges[0:1] = back:front z boundaries of particle tags
# tedges[2:3] = lower:upper y boundaries of particle tags
      tedges = numpy.empty((idps),float_type,'F')
# electron trajectories
   if (in3.ndt==1):
# set particle tags: updates nprobt, tedges, ppart and possibly iprobt
      mdiag3.psetptraj3(ppart,tedges,kpic,iprobt,kstrt,in3.nst,in3.nvpy,
                        in3.nvpz,in3.vtx,in3.vtsx,in3.dvtx,np,
                        in3.nprobt,irc)
# estimate maximum electron velocity or momentum
      if (in3.nst==3):
         ws[0] = 0.0
         if (npxyz > 0.0):
            ws[0] = 4.0*in3.vtx+abs(in3.vx0)
            ws[0] = max(ws[0],4.0*in3.vty+abs(in3.vy0))
            ws[0] = max(ws[0],4.0*in3.vtz+abs(in3.vz0))
         if (npxyzb > 0.0):
            ws[0] = max(ws[0],4.0*in3.vtdx+abs(in3.vdx))
            ws[0] = max(ws[0],4.0*in3.vtdy+abs(in3.vdy))
            ws[0] = max(ws[0],4.0*in3.vtdz+abs(in3.vdz))
# ion trajectories
   elif (in3.ndt==2):
# set particle tags: updates nprobt, tedges, pparti and possibly iprobt
      mdiag3.psetptraj3(pparti,tedges,kipic,iprobt,kstrt,in3.nst,
                        in3.nvpy,in3.nvpz,vtxi,in3.vtsx,in3.dvtx,npi,
                        in3.nprobt,irc)
# estimate maximum ion velocity or momentum
      if (in3.nst==3):
         ws[0] = 0.0
         if (npxyzi > 0.0):
            ws[0] = 4.0*vtxi+abs(in3.vxi0)
            ws[0] = max(ws[0],4.0*vtyi+abs(in3.vyi0))
            ws[0] = max(ws[0],4.0*vtzi+abs(in3.vzi0))
         if (npxyzbi > 0.0):
            ws[0] = max(ws[0],4.0*vtdxi+abs(in3.vdxi))
            ws[0] = max(ws[0],4.0*vtdyi+abs(in3.vdyi))
            ws[0] = max(ws[0],4.0*vtzi+abs(in3.vzi0))
# electron or ion trajectories
   if ((in3.ndt==1) or (in3.ndt==2)):
      if (in3.nprobt > 16777215):
         print "nprobt overflow = ", in3.nprobt
         mpplib3.ppexit(); exit(1)
      in3.ndimp = idimp
# partt = particle trajectories tracked
      partt = numpy.empty((idimp,in3.nprobt),float_type,'F')
      ftname = "tr3." + cdrun
      in3.ftname[:] = ftname
      if ((in3.nst==1) or (in3.nst==2)):
         it = int((nloop - 1)/in3.ntt + 1); itt = 0
# partd = trajectory time history array
         partd = numpy.zeros((it,idimp,in3.nprobt),float_type,'F')
# open file for trajectory data: updates ntrec and possibly iut
         if (kstrt==1):
            if (in3.ntrec==0):
               mdiag3.dafopentr3(partt,iut,in3.ntrec,ftname)
      elif (in3.nst==3):
# fvtp = velocity distribution function for test particles
         fvtp = numpy.empty((2*in3.nmv+2,in3.ndim),float_type,'F')
# fvmtp = vdrift, vth, and entropy for test particles
         fvmtp = numpy.empty((in3.ndim,3),float_type,'F')
         fetp = numpy.empty((in3.ndim,0),float_type,'F')
         fvtp[2*in3.nmv+1,:] = 2.0*ws[0]
# open file for test particle diagnostic: updates ntrec and possibly iut
         if (kstrt==1):
            if (in3.ntrec==0):
               ws[0] = 0.0
               mdiag3.dafopenfv3(fvmtp,fvtp,fetp,ws,iut,in3.ntrec,ftname)

# initialize phase space diagnostic
if (in3.nts > 0):
   nmv21 = 2*in3.nmv + 1
   in3.mvx = min(in3.mvx,nx); in3.mvy = min(in3.mvy,nypmn[0])
   in3.mvz = min(in3.mvz,nzpmn[0])
   in3.nsxb = int((nx - 1)/in3.mvx + 1)
   in3.nsyb = int((ny - 1)/in3.mvy + 1)
   in3.nszb = int((nz - 1)/in3.mvz + 1)
   nyb = int(float(noff[0] - 1)/float(in3.mvy))
   nzb = int(float(noff[1] - 1)/float(in3.mvz))
   nyb = int((noff[0] + nyzp[0] - 1)/in3.mvy) - nyb
   nzb = int((noff[1] + nyzp[1] - 1)/in3.mvz) - nzb
   if (js==0):
      nyb += 1
   if (ks==0):
      nzb += 1
   itot[0] = nyb; itot[1] = nzb
   mppmod3.mpimax(itot[:2],tinit)
   nybmx = itot[0]; nzbmx = itot[1]
# electron phase space diagnostic
   if ((in3.nds==1) or (in3.nds==3)):
# estimate maximum electron velocity or momentum
      ws[0] = 0.0
      if (npxyz > 0.0):
         ws[0] = 4.0*in3.vtx+abs(in3.vx0)
         ws[0] = max(ws[0],4.0*in3.vty+abs(in3.vy0))
         ws[0] = max(ws[0],4.0*in3.vtz+abs(in3.vz0))
      if (npxyzb > 0.0):
         ws[0] = max(ws[0],4.0*in3.vtdx+abs(in3.vdx))
         ws[0] = max(ws[0],4.0*in3.vtdy+abs(in3.vdy))
         ws[0] = max(ws[0],4.0*in3.vtdz+abs(in3.vdz))
# fvs = global electron phase space distribution functions
      fvs = numpy.empty((nmv21+1,in3.ndim,in3.nsxb,nybmx+1,nzb+1),
                        float_type,'F')
      fvs[nmv21,:,0,0,0] = 1.25*ws[0]
# open file for electron phase space data:
# updates nserec and possibly iuse
# opens a new fortran unformatted stream file
      if (in3.nserec==0):
         if (kstrt==1):
            fsename = "pse3." + cdrun
            in3.fsename[:] = fsename
            iuse =  mdiag3.get_funit(iuse)
            mdiag3.fnopens3(iuse,fsename)
# writes distributed non-uniform partition information
            mppmod3.mpwrncomp3(nyb,nzb,in3.nvpy,in3.nvpz,iuse)
            in3.nserec = 1
# ion phase space diagnostic
   if (in3.movion==1):
      if ((in3.nds==2) or (in3.nds==3)):
# estimate maximum ion velocity or momentum
         ws[0] = 0.0
         if (npxyzi > 0.0):
            ws[0] = 4.0*vtxi+abs(in3.vxi0)
            ws[0] = max(ws[0],4.0*vtyi+abs(in3.vyi0))
            ws[0] = max(ws[0],4.0*vtzi+abs(in3.vzi0))
         if (npxyzbi > 0.0):
            ws[0] = max(ws[0],4.0*vtdxi+abs(in3.vdxi))
            ws[0] = max(ws[0],4.0*vtdyi+abs(in3.vdyi))
            ws[0] = max(ws[0],4.0*vtdzi+abs(in3.vdzi))
# fvsi = global ion phase space distribution functions
         fvsi = numpy.empty((nmv21+1,in3.ndim,in3.nsxb,nybmx+1,nzb+1),
                            float_type,'F')
         fvsi[nmv21,:,0,0,0] = 1.25*ws[0]
# open file for ion phase space data:
# updates nsirec and possibly iusi
# opens a new fortran unformatted stream file
         if (in3.nsirec==0):
            if (kstrt==1):
               fsiname = "psi3." + cdrun
               in3.fsiname[:] = fsiname
               iusi =  mdiag3.get_funit(iusi)
               mdiag3.fnopens3(iusi,fsiname)
# writes distributed non-uniform partition information
               mppmod3.mpwrncomp3(nyb,nzb,in3.nvpy,in3.nvpz,iusi)
               in3.nsirec = 1

# initialization time
dtimer(dtime,itime,1)
tinit = 0.0
tinit += float(dtime)
# start timing loop
dtimer(dtime,ltime,-1)

if (kstrt==1):
   print >> iuot, "program mpdbeps3"

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):

# print time step
   if (kstrt==1):
      if (in3.ntw > 0):
         it = int(ntime/in3.ntw)
         if (ntime==in3.ntw*it):
            print >> iuot, "ntime = ", ntime

# deposit electron current with OpenMP: updates cue
# updates ppart and cue, and possibly ncl, ihole, irc
   mpush3.mpset_pvzero3(cue,tdjpost,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
   mcurd3.wmpdjpost3(ppart,cue,kpic,ncl,ihole,noff,nyzp,in3.qme,zero,
                     in3.ci,tdjpost,nx,ny,nz,in3.mx,in3.my,in3.mz,mx1,
                     myp1,ipbc,in3.popt,in3.relativity,plist,irc)
# add guard cells with OpenMP: updates cue
   ompplib3.wmpnacguard3(cue,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# save electron current for electron current diagnostic later
   if (in3.ndc==0):
      if (in3.ntje > 0):
         it = int(ntime/in3.ntje)
         if (ntime==in3.ntje*it):
            oldcue[:,:,:,:] = cue


# deposit ion current with OpenMP:
   if (in3.movion==1):
# updates cui
      mpush3.mpset_pvzero3(cui,tdjpost,in3.mx,in3.my,in3.mz,mx1,myp1,
                           mzp1)
      mcurd3.wmpdjpost3(pparti,cui,kipic,ncl,ihole,noff,nyzp,in3.qmi,
                        zero,in3.ci,tdjpost,nx,ny,nz,in3.mx,in3.my,
                        in3.mz,mx1,myp1,ipbc,in3.popt,in3.relativity,
                        plist,irc)
# add guard cells with OpenMP: updates cui
      ompplib3.wmpnacguard3(cui,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# deposit electron charge with OpenMP: updates qe
   mpush3.mpset_pszero3(qe,tdpost,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
   mpush3.mppost3(ppart,qe,kpic,noff,in3.qme,tdpost,in3.mx,in3.my,
                  in3.mz,mx1,myp1,in3.popt)
# add guard cells with OpenMP: updates qe
   ompplib3.wmpaguard3(qe,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# electron density diagnostic
   if (in3.ntde > 0):
      it = int(ntime/in3.ntde)
      if (ntime==in3.ntde*it):
         sfield[:,:,:] = -numpy.copy(qe)
# transform electron density to fourier space: updates qt
# moves data to uniform partition
         isign = -1
         ompplib3.wmpfft3r(sfield,qt,noff,nyzp,isign,mixup,sct,tfft,
                           tfmov,in3.indx,in3.indy,in3.indz,kstrt,
                           in3.nvpy,in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,
                           mterf,ierr)
# calculate smoothed electron density in fourier space: updates sfieldc
         mfield3.mpsmooth3(qt,sfieldc,ffc,tfield,nx,ny,nz,kstrt,
                           in3.nvpy,in3.nvpz)
# store selected fourier modes: updates denet
         mfield3.mprdmodes3(sfieldc,denet,tfield,nx,ny,nz,in3.modesxde,
                         in3.modesyde,in3.modeszde,kstrt,in3.nvpy,
                         in3.nvpz)
# transform smoothed electron density to real space: updates sfield
         isign = 1
         mfft3.mpfft3r(sfield,sfieldc,isign,mixup,sct,tfft,in3.indx,
                       in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,kxyp,
                       kyp,kyzp,kzp)
# write real space diagnostic output: updates nderec
         mppmod3.mpwrite3(sfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iude,
                          in3.nderec)
#        mppmod3.mpread3(sfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iude,
#                        in3.nderec,irc)

# deposit ion charge with OpenMP: updates qi
   if (in3.movion==1):
      mpush3.mpset_pszero3(qi,tdpost,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
      mpush3.mppost3(pparti,qi,kipic,noff,in3.qmi,tdpost,in3.mx,in3.my,
                     in3.mz,mx1,myp1,in3.popt)
# add guard cells with OpenMP: updates qi
      ompplib3.wmpaguard3(qi,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# ion density diagnostic
   if (in3.movion==1):
      if (in3.ntdi > 0):
         it = int(ntime/in3.ntdi)
         if (ntime==in3.ntdi*it):
            sfield[:,:,:] = numpy.copy(qi)
# transform ion density to fourier space: updates qt
# moves data to uniform partition
            isign = -1
            ompplib3.wmpfft3r(sfield,qt,noff,nyzp,isign,mixup,sct,tfft,
                              tfmov,in3.indx,in3.indy,in3.indz,kstrt,
                              in3.nvpy,in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,
                              mterf,ierr)
# calculate smoothed ion density in fourier space: updates sfieldc
            mfield3.mpsmooth3(qt,sfieldc,ffc,tfield,nx,ny,nz,kstrt,
                              in3.nvpy,in3.nvpz)
# store selected fourier modes: updates denit
            mfield3.mprdmodes3(sfieldc,denit,tfield,nx,ny,nz,
                               in3.modesxdi,in3.modesydi,in3.modeszdi,
                               kstrt,in3.nvpy,in3.nvpz)
# transform smoothed ion density to real space: updates sfield
            isign = 1
            mfft3.mpfft3r(sfield,sfieldc,isign,mixup,sct,tfft,in3.indx,
                          in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,
                          kxyp,kyp,kyzp,kzp)
# write real space diagnostic output: updates ndirec
            mppmod3.mpwrite3(sfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,
                             iudi,in3.ndirec)
#           mppmod3.mpread3(sfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iudi,
#                            in3.ndirec,irc)

# add electron and ion densities: updates qe
   mfield3.mpaddqei3(qe,qi,nyzp,tfield,nx)

# add electron and ion current densities: updates cue
   if (in3.movion==1):
      mfield3.mpaddcuei3(cue,cui,nyzp,tfield,nx)

# transform charge to fourier space with OpenMP:
# moves data to uniform partition
# updates qt, mterf, and ierr, modifies qe
   isign = -1
   ompplib3.wmpfft3r(qe,qt,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                     in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,
                     kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# calculate longitudinal force/charge in fourier space with OpenMP:
# updates fxyzt, we
   mfield3.mppois3(qt,fxyzt,ffc,we,tfield,nx,ny,nz,kstrt,in3.nvpy,
                   in3.nvpz)

# transform longitudinal electric force to real space with OpenMP:
# moves data to non-uniform partition
# updates fxyze, mterf, and ierr, modifies fxyzt
   isign = 1
   ompplib3.wmpfft3rn(fxyze,fxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                      in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                      in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# transform current to fourier space with OpenMP:
# moves data to uniform partition
# updates cut, mterf, and ierr, modifies cue
   isign = -1
   ompplib3.wmpfft3rn(cue,cut,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                      in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                      in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# take transverse part of current with OpenMP: updates cut
   mfield3.mpcuperp3(cut,tfield,nx,ny,nz,kstrt,in3.nvpy,in3.nvpz)

# calculate magnetic field in fourier space with OpenMP:
# updates bxyzt, wb
   mfield3.mpbbpois3(cut,bxyzt,ffc,in3.ci,wb,tfield,nx,ny,nz,kstrt,
                     in3.nvpy,in3.nvpz)

# transform magnetic force to real space with OpenMP:
# moves data to non-uniform partition
# updates bxyze, modifies bxyzt
   isign = 1
   ompplib3.wmpfft3rn(bxyze,bxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                      in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                      in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# add constant to magnetic field with OpenMP: updates bxyze
   mfield3.mpbaddext3(bxyze,nyzp,tfield,in3.omx,in3.omy,in3.omz,nx)

# copy guard cells with OpenMP: updates fxyze, bxyze
   ompplib3.wmpncguard3(fxyze,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)
   ompplib3.wmpncguard3(bxyze,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# add longitudinal and old transverse electric fields with OpenMP:
# updates exyze
   mfield3.mpaddvrfield3(exyze,cus,fxyze,tfield)

# deposit electron acceleration density and momentum flux with OpenMP:
# updates dcu, amu
   mpush3.mpset_pvzero3(dcu,tdcjpost,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
   mpush3.mpset_pvzero3(amu,tdcjpost,in3.mx,in3.my,in3.mz,mx1,myp1,mzp1)
   mdpush3.wmpgdjpost3(ppart,exyze,bxyze,kpic,noff,nyzp,dcu,amu,in3.qme,
                       qbme,in3.dt,in3.ci,tdcjpost,nx,in3.mx,in3.my,
                       in3.mz,mx1,myp1,in3.popt,in3.relativity)
# add guard cells with OpenMP: updates dcu, amu
   ompplib3.wmpnacguard3(dcu,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)
   ompplib3.wmpnacguard3(amu,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# deposit ion acceleration density and momentum flux with OpenMP:
# updates dcui, amui
   if (in3.movion==1):
      mpush3.mpset_pvzero3(dcui,tdcjpost,in3.mx,in3.my,in3.mz,mx1,myp1,
                           mzp1)
      mpush3.mpset_pvzero3(amui,tdcjpost,in3.mx,in3.my,in3.mz,mx1,myp1,
                           mzp1)
      mdpush3.wmpgdjpost3(pparti,exyze,bxyze,kipic,noff,nyzp,dcui,amui,
                          in3.qmi,qbmi,in3.dt,in3.ci,tdcjpost,nx,in3.mx,
                          in3.my,in3.mz,mx1,myp1,in3.popt,in3.relativity)
# add guard cells with OpenMP: updates dcui, amui
      ompplib3.wmpnacguard3(dcui,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)
      ompplib3.wmpnacguard3(amui,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)
# add electron and ion momentum flux densities: updates dcu, amu
      mfield3.mpaddcuei3(dcu,dcui,nyzp,tfield,nx)
      mfield3.mpaddamui3(amu,amui,nyzp,tfield,nx)

# add old scaled electric field with standard procedure: updates dcu
   mdpush3.mpascfguard3(dcu,cus,nyzp,q2m0,tdcjpost,nx)

# transform acceleration density and momentum flux to fourier space
# with OpenMP and moves data to uniform partition
# updates dcut, amut, modifies dcu, amu
   isign = -1
   ompplib3.wmpfft3rn(dcu,dcut,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                      in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                      in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
   ompplib3.wmpfft3rn(amu,amut,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                      in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                      in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# take transverse part of time derivative of current with OpenMP:
# updates dcut
   mfield3.mpadcuperp3(dcut,amut,tfield,nx,ny,nz,kstrt,in3.nvpy,
                       in3.nvpz)

# calculate transverse electric field with OpenMP:
# updates exyzt, wf
   mfield3.mpepois3(dcut,exyzt,ffe,affp,in3.ci,wf,tfield,nx,ny,nz,kstrt,
                    in3.nvpy,in3.nvpz)

# transform transverse electric field to real space with OpenMP:
# moves data to non-uniform partition
# updates cus, modifies exyzt
   isign = 1
   ompplib3.wmpfft3rn(cus,exyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                      in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                      in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# copy guard cells with OpenMP: updates cus
   ompplib3.wmpncguard3(cus,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# add longitudinal and transverse electric fields with OpenMP:
#  exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
   mfield3.mpaddvrfield3(exyze,cus,fxyze,tfield)

# inner iteration loop
   for k in xrange(0,in3.ndc):

# deposit electron current and acceleration density and momentum flux
# with OpenMP: updates cue, dcu, amu
      mpush3.mpset_pvzero3(cue,tdjpost,in3.mx,in3.my,in3.mz,mx1,myp1,
                           mzp1)
      mpush3.mpset_pvzero3(dcu,tdcjpost,in3.mx,in3.my,in3.mz,mx1,myp1,
                           mzp1)
      mpush3.mpset_pvzero3(amu,tdcjpost,in3.mx,in3.my,in3.mz,mx1,myp1,
                           mzp1)
      mdpush3.wmpgdcjpost3(ppart,exyze,bxyze,kpic,noff,nyzp,cue,dcu,amu,
                           in3.qme,qbme,in3.dt,in3.ci,tdcjpost,nx,
                           in3.mx,in3.my,in3.mz,mx1,myp1,in3.popt,
                           in3.relativity)
# add guard cells for current, acceleration density, and momentum flux
# with OpenMP: updates cue, dcu, amu
      ompplib3.wmpnacguard3(cue,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)
      ompplib3.wmpnacguard3(dcu,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)
      ompplib3.wmpnacguard3(amu,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# save electron current for electron current diagnostic later
      if (k==in3.ndc):
         if (in3.ntje > 0):
            it = int(ntime/in3.ntje)
            if (ntime==in3.ntje*it):
               oldcue[:,:,:,:] = cue

# deposit ion current and acceleration density and momentum flux
# with OpenMP: updates cui, dcui, amui
      if (in3.movion==1):
         mpush3.mpset_pvzero3(cui,tdjpost,in3.mx,in3.my,in3.mz,mx1,myp1,
                              mzp1)
         mpush3.mpset_pvzero3(dcui,tdcjpost,in3.mx,in3.my,in3.mz,mx1,
                              myp1,mzp1)
         mpush3.mpset_pvzero3(amui,tdcjpost,in3.mx,in3.my,in3.mz,mx1,
                              myp1,mzp1)
         mdpush3.wmpgdcjpost3(pparti,exyze,bxyze,kipic,noff,nyzp,cui,
                              dcui,amui,in3.qmi,qbmi,in3.dt,in3.ci,
                              tdcjpost,nx,in3.mx,in3.my,in3.mz,mx1,myp1,
                              in3.popt,in3.relativity)
# add guard cells for current, acceleration density, and momentum flux
# with OpenMP: updates cui, dcui, amui
         ompplib3.wmpnacguard3(cui,nyzp,tguard,nx,kstrt,in3.nvpy,
                               in3.nvpz)
         ompplib3.wmpnacguard3(dcui,nyzp,tguard,nx,kstrt,in3.nvpy,
                               in3.nvpz)
         ompplib3.wmpnacguard3(amui,nyzp,tguard,nx,kstrt,in3.nvpy,
                               in3.nvpz)
# add electron and ion densities: updates cue, dcu, amu
         mfield3.mpaddcuei3(cue,cui,nyzp,tfield,nx)
         mfield3.mpaddcuei3(dcu,dcui,nyzp,tfield,nx)
         mfield3.mpaddamui3(amu,amui,nyzp,tfield,nx)

# add scaled electric field with standard procedure: updates dcu
      mdpush3.mpascfguard3(dcu,cus,nyzp,q2m0,tdcjpost,nx)

# transform current to fourier space with OpenMP:
# moves data to uniform partition
# updates cut, mterf, and ierr, modifies cue
      isign = -1
      ompplib3.wmpfft3rn(cue,cut,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                         in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                         in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# take transverse part of current with OpenMP: updates cut
      mfield3.mpcuperp3(cut,tfield,nx,ny,nz,kstrt,in3.nvpy,in3.nvpz)

# calculate magnetic field in fourier space with OpenMP:
# updates bxyzt, wb
      mfield3.mpbbpois3(cut,bxyzt,ffc,in3.ci,wb,tfield,nx,ny,nz,kstrt,
                        in3.nvpy,in3.nvpz)

# transform magnetic force to real space with OpenMP:
# moves data to non-uniform partition
# updates bxyze, modifies bxyzt
      isign = 1
      ompplib3.wmpfft3rn(bxyze,bxyzt,noff,nyzp,isign,mixup,sct,tfft,
                         tfmov,in3.indx,in3.indy,in3.indz,kstrt,
                         in3.nvpy,in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,
                         mterf,ierr)

# add constant to magnetic field with OpenMP: updates bxyze
      mfield3.mpbaddext3(bxyze,nyzp,tfield,in3.omx,in3.omy,in3.omz,nx)

# transform acceleration density and momentum flux to fourier space
# with OpenMP and moves data to uniform partition
# updates dcut, amut, modifies dcu, amu
      isign = -1
      ompplib3.wmpfft3rn(dcu,dcut,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                         in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                         in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
      ompplib3.wmpfft3rn(amu,amut,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                         in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                         in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# take transverse part of time derivative of current with OpenMP:
# updates dcut
      mfield3.mpadcuperp3(dcut,amut,tfield,nx,ny,nz,kstrt,in3.nvpy,
                          in3.nvpz)

# calculate transverse electric field with OpenMP:
# updates exyzt, wf
      mfield3.mpepois3(dcut,exyzt,ffe,affp,in3.ci,wf,tfield,nx,ny,nz,
                       kstrt,in3.nvpy,in3.nvpz)

# transform transverse electric field to real space with OpenMP:
# moves data to non-uniform partition
# updates cus, modifies exyzt
      isign = 1
      ompplib3.wmpfft3rn(cus,exyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,
                         in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
                         in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)

# copy guard cells with OpenMP: updates bxyze, cus
      ompplib3.wmpncguard3(bxyze,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)
      ompplib3.wmpncguard3(cus,nyzp,tguard,nx,kstrt,in3.nvpy,in3.nvpz)

# add longitudinal and transverse electric fields with OpenMP:
# exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
      mfield3.mpaddvrfield3(exyze,cus,fxyze,tfield)
   pass

# potential diagnostic
   if (in3.ntp > 0):
      it = int(ntime/in3.ntp)
      if (ntime==in3.ntp*it):
# calculate potential in fourier space: updates sfieldc
         mfield3.mppot3(qt,sfieldc,ffc,ws,tfield,nx,ny,nz,kstrt,
                        in3.nvpy,in3.nvpz)
# store selected fourier modes: updates pott
         mfield3.mprdmodes3(sfieldc,pott,tfield,nx,ny,nz,in3.modesxp,
                            in3.modesyp,in3.modeszp,kstrt,in3.nvpy,
                            in3.nvpz)
# transform potential to real space: updates sfield
         isign = 1
         mfft3.mpfft3r(sfield,sfieldc,isign,mixup,sct,tfft,in3.indx,
                       in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,kxyp,
                       kyp,kyzp,kzp)
# write real space diagnostic output: updates nprec
         mppmod3.mpwrite3(sfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iup,
                          in3.nprec)
#        mppmod3.mpread3(sfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iup,
#                        in3.nprec,irc)

# longitudinal efield diagnostic
   if (in3.ntel > 0):
      it = int(ntime/in3.ntel)
      if (ntime==in3.ntel*it):
# calculate longitudinal efield in fourier space: updates vfieldc
         mfield3.mpelfield3(qt,vfieldc,ffc,ws,tfield,nx,ny,nz,kstrt,
                            in3.nvpy,in3.nvpz)
# store selected fourier modes: updates elt
         mfield3.mprdvmodes3(vfieldc,elt,tfield,nx,ny,nz,in3.modesxel,
                             in3.modesyel,in3.modeszel,kstrt,in3.nvpy,
                             in3.nvpz)
# transform longitudinal efield to real space: updates vfield
         isign = 1
         mfft3.mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,in3.indx,
                        in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,kxyp,
                        kyp,kyzp,kzp)
# write real space diagnostic output: updates nelrec
         mppmod3.mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iuel,
                           in3.nelrec)
#        mppmod3.mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iuel,
#                         in3.nelrec,irc)

# vector potential diagnostic
   if (in3.nta > 0):
      it = int(ntime/in3.nta)
      if (ntime==in3.nta*it):
# calculate vector potential in fourier space: updates vfieldc
         mfield3.mpapot3(cut,vfieldc,ffc,in3.ci,ws,tfield,nx,ny,nz,
                         kstrt,in3.nvpy,in3.nvpz)
# store selected fourier modes: updates vpott
         mfield3.mprdvmodes3(vfieldc,vpott,tfield,nx,ny,nz,in3.modesxa,
                             in3.modesya,in3.modesza,kstrt,in3.nvpy,
                             in3.nvpz)
# transform vector potential to real space: updates vfield
         isign = 1
         mfft3.mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,in3.indx,
                        in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,kxyp,
                        kyp,kyzp,kzp)
# write real space diagnostic output: updates narec
         mppmod3.mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iua,
                           in3.narec)
#        mppmod3.mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iua,
#                         in3.narec,irc)

# transverse efield diagnostic
   if (in3.ntet > 0):
      it = int(ntime/in3.ntet)
      if (ntime==in3.ntet*it):
# calculate transverse efield in fourier space: updates vfieldc
         mfield3.mpetfield3(dcut,vfieldc,ffe,affp,in3.ci,ws,tfield,nx,
                            ny,nz,kstrt,in3.nvpy,in3.nvpz)
# store selected fourier modes: updates ett
         mfield3.mprdvmodes3(vfieldc,ett,tfield,nx,ny,nz,in3.modesxet,
                             in3.modesyet,in3.modeszet,kstrt,in3.nvpy,
                             in3.nvpz)
# transform transverse efield to real space: updates vfield
         isign = 1
         mfft3.mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,in3.indx,
                        in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,kxyp,
                        kyp,kyzp,kzp)
# write real space diagnostic output: updates netrec
         mppmod3.mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iuet,
                           in3.netrec)
#        mppmod3.mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iuet,
#                         in3.netrec,irc)if

# magnetic field diagnostic
   if (in3.ntb > 0):
      it = int(ntime/in3.ntb)
      if (ntime==in3.ntb*it):
# calculate magnetic field in fourier space: updates vfieldc
         mfield3.mpibpois3(cut,vfieldc,ffc,in3.ci,ws,tfield,nx,ny,nz,
                           kstrt,in3.nvpy,in3.nvpz)
# store selected fourier modes: updates bt
         mfield3.mprdvmodes3(vfieldc,bt,tfield,nx,ny,nz,in3.modesxb,
                             in3.modesyb,in3.modeszb,kstrt,in3.nvpy,
                             in3.nvpz)
# transform magnetic field to real space: updates vfield
         isign = 1
         mfft3.mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,in3.indx,
                        in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,kxyp,
                        kyp,kyzp,kzp)
# write real space diagnostic output: updates nbrec
         mppmod3.mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iub,
                           in3.nbrec)
#        mppmod3.mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iub,
#                         in3.nbrec,irc)

# electron current density diagnostic
   if (in3.ntje > 0):
      it = int(ntime/in3.ntje)
      if (ntime==in3.ntje*it):
         vfield[:,:,:,:] = oldcue
# transform electron current density to fourier space: updates cut
# moves data to uniform partition
         isign = -1
         ompplib3.wmpfft3rn(vfield,cut,noff,nyzp,isign,mixup,sct,tfft,
                            tfmov,in3.indx,in3.indy,in3.indz,kstrt,
                            in3.nvpy,in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,
                            mterf,ierr)
# calculate smoothed electron current in fourier space: updates vfieldc
         mfield3.mpsmooth33(cut,vfieldc,ffc,tfield,nx,ny,nz,kstrt,
                            in3.nvpy,in3.nvpz)
# store selected fourier modes: updates curet
         mfield3.mprdvmodes3(vfieldc,curet,tfield,nx,ny,nz,in3.modesxje,
                             in3.modesyje,in3.modeszje,kstrt,in3.nvpy,
                             in3.nvpz)
# transform smoothed electron current to real space: updates vfield
         isign = 1
         mfft3.mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,in3.indx,
                        in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,kxyp,
                        kyp,kyzp,kzp)
# write real space diagnostic output: updates njerec
         mppmod3.mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iuje,
                           in3.njerec)
#        mppmod3.mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iuje,
#                         in3.njerec,irc)

# ion current density diagnostic
   if (in3.movion==1):
      if (in3.ntji > 0):
         it = int(ntime/in3.ntji)
         if (ntime==in3.ntji*it):
            vfield[:,:,:,:] = cui
# transform ion current density to fourier space: updates cut
# moves data to uniform partition
            isign = -1
            ompplib3.wmpfft3rn(vfield,cut,noff,nyzp,isign,mixup,sct,
            tfft,tfmov,in3.indx,in3.indy,in3.indz,kstrt,in3.nvpy,
            in3.nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
# calculate smoothed ion current in fourier space: updates vfieldc
            mfield3.mpsmooth33(cut,vfieldc,ffc,tfield,nx,ny,nz,kstrt,
                               in3.nvpy,in3.nvpz)
# store selected fourier modes: updates curit
            mfield3.mprdvmodes3(vfieldc,curit,tfield,nx,ny,nz,
                                in3.modesxji,in3.modesyji,in3.modeszji,
                                kstrt,in3.nvpy,in3.nvpz)
# transform smoothed ion current to real space: updates vfield
            isign = 1
            mfft3.mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,in3.indx,
                           in3.indy,in3.indz,kstrt,in3.nvpy,in3.nvpz,
                           kxyp,kyp,kyzp,kzp)
# write real space diagnostic output: updates njirec
            mppmod3.mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,
                              iuji,in3.njirec)
#           mppmod3.mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,
#                            iuji,in3.njirec,irc)

# fluid moments diagnostic
   if (in3.ntfm > 0):
      it = int(ntime/in3.ntfm)
      if (ntime==in3.ntfm*it):
# calculate electron fluid moments
         if ((in3.ndfm==1) or (in3.ndfm==3)):
            mpush3.mpset_pvzero3(fmse,tdiag,in3.mx,in3.my,in3.mz,mx1,
                                 myp1,mzp1)
            mdiag3.wmgbprofx3(ppart,exyze,bxyze,fmse,kpic,noff,nyzp,
                              qbme,in3.dt,in3.ci,tdiag,in3.npro,nx,
                              in3.mx,in3.my,in3.mz,mx1,myp1,
                              in3.relativity)
# add guard cells with OpenMP: updates fmse
            ompplib3.wmpnacguard3(fmse,nyzp,tdiag,nx,kstrt,in3.nvpy,
                     in3.nvpz)
# moves vector grid fmse from non-uniform to uniform partition
            isign = -1
            mppmod3.mpfnmove3(fmse,noff,nyzp,isign,tdiag,kyp,kzp,ny,nz,
                              kstrt,in3.nvpy,in3.nvpz,mterf,irc)
            if (irc[0] != 0):
                  mpplib3.ppexit(); exit(1)
# calculates fluid quantities from fluid moments
            mdiag3.mpfluidqs3(fmse,tdiag,in3.npro,nx,ny,nz,kstrt,
                              in3.nvpy,kyp,kzp)
# write real space diagnostic output: updates nferec
            mppmod3.mpvwrite3(fmse,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iufe,
                              in3.nferec)
#           mppmod3.mpvread3(fmse,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,iufe,
#    &in3.nferec,irc)
# calculate ion fluid moments
         if (in3.movion==1):
            if ((in3.ndfm==2) or (in3.ndfm==3)):
               mpush3.mpset_pvzero3(fmsi,tdiag,in3.mx,in3.my,in3.mz,mx1,
                                    myp1,mzp1)
               mdiag3.wmgbprofx3(pparti,exyze,bxyze,fmsi,kipic,noff,
                                 nyzp,qbmi,in3.dt,in3.ci,tdiag,in3.npro,
                                 nx,in3.mx,in3.my,in3.mz,mx1,myp1,
                                 in3.relativity)
               fmsi[:,:,:] = in3.rmass*fmsi
# add guard cells with OpenMP: updates fmsi
               ompplib3.wmpnacguard3(fmsi,nyzp,tdiag,nx,kstrt,in3.nvpy,
                                     in3.nvpz)
# moves vector grid fmsi from non-uniform to uniform partition
               isign = -1
               mppmod3.mpfnmove3(fmsi,noff,nyzp,isign,tdiag,kyp,kzp,ny,
                                 nz,kstrt,in3.nvpy,in3.nvpz,mterf,irc)
               if (irc[0] != 0):
                  mpplib3.ppexit(); exit(1)
# calculates fluid quantities from fluid moments: updates fmsi
               mdiag3.mpfluidqs3(fmsi,tdiag,in3.npro,nx,ny,nz,kstrt,
                                 in3.nvpy,kyp,kzp)
# write real space diagnostic output: updates nfirec
               mppmod3.mpvwrite3(fmsi,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,
                                 iufi,in3.nfirec)
#              mppmod3.mpvread3(fmsi,tdiag,nx,ny,nz,kyp,kzp,in3.nvpy,
#                               iufi,in3.nfirec,irc)dif

# velocity-space diagnostic
   if (in3.ntv > 0):
      it = int(ntime/in3.ntv)
      if (ntime==in3.ntv*it):
         if ((in3.ndv==1) or (in3.ndv==3)):
# calculate electron cartesian distribution function and moments
            if ((in3.nvft==1) or (in3.nvft==3)):
               mdiag3.mpvpdist3(ppart,kpic,fv,sfv,fvm,tdiag,nvp,in3.nmv)
# store time history electron vdrift, vth, and entropy
               fvtm[itv,:,:] = fvm
               if (in3.movion==0):
                   itv += 1
# calculate electron cylindrical distribution function and moments
            if ((in3.nvft==4) or (in3.nvft==5)):
               mdiag3.mpvbpdist3(ppart,kpic,fv,sfv,fvm,in3.omx,in3.omy,
                                 in3.omz,tdiag,in3.nmv)
# electron energy distribution
            if ((in3.nvft==2) or (in3.nvft==3) or (in3.nvft==5)):
               mdiag3.mperpdist3(ppart,kpic,fe,sfv,eci,wk,tdiag,in3.nmv)
# write electron velocity-space diagnostic output: updates nverec
            if (kstrt==1):
             mdiag3.dafwritefv3(fvm,fv,fe,wk,tdiag,iuve,in3.nverec)
# ion distribution functions
         if (in3.movion==1):
             if ((in3.ndv==2) or (in3.ndv==3)):
# calculate ion cartesian distribution function and moments
                if ((in3.nvft==1) or (in3.nvft==3)):
                   mdiag3.mpvpdist3(pparti,kipic,fvi,sfv,fvmi,tdiag,nvp,
                                    in3.nmv)
# store time history of ion vdrift, vth, and entropy
                   fvtmi[itv,:,:] = fvmi
                   itv += 1
# calculate ion cylindrical distribution function and moments
                if ((in3.nvft==4) or (in3.nvft==5)):
                   mdiag3.mpvbpdist3(pparti,kipic,fvi,sfv,fvmi,in3.omx,
                                     in3.omy,in3.omz,tdiag,in3.nmv)
# ion energy distribution
                if ((in3.nvft==2) or (in3.nvft==3) or (in3.nvft==5)):
                   mdiag3.mperpdist3(pparti,kipic,fei,sfv,eci,wk,tdiag,
                                     in3.nmv)
                   wk[0] = in3.rmass*wk[0]
# write ion velocity-space diagnostic output: updates nvirec
                if (kstrt==1):
                   mdiag3.dafwritefv3(fvmi,fvi,fei,wk,tdiag,iuvi,
                                         in3.nvirec)

# trajectory diagnostic
   if (in3.ntt > 0):
      it = int(ntime/in3.ntt)
      if (ntime==in3.ntt*it):
         ierr[0] = 0
# copies tagged electrons in ppart to array partt: updates partt, numtp
         if (in3.ndt==1):
            mdiag3.mptraj3(ppart,kpic,partt,tdiag,numtp,ierr)
# copies tagged ions in ppart to array partt: updates partt, numtp
         elif (in3.ndt==2):
            if (in3.movion==1):
               mdiag3.mptraj3(pparti,kipic,partt,tdiag,numtp,ierr)
         if (ierr[0] != 0):
# partt overflow
            if (ierr[0] <  0):
               in3.nprobt = numtp[0]
               partt = numpy.empty((idimp,in3.nprobt),float_type,'F')
               ierr[0] = 0
# copies tagged electrons in ppart to array partt: updates partt, numtp
               if (in3.ndt==1):
                  mdiag3.mptraj3(ppart,kpic,partt,tdiag,numtp,ierr)
# copies tagged ions in ppart to array partt: updates partt, numtp
               elif (in3.ndt==2):
                  if (in3.movion==1):
                     mdiag3.mptraj3(pparti,kipic,partt,tdiag,numtp,ierr)
            if (ierr[0] != 0):
               mpplib3.ppexit();  exit(1)
# electron or ion trajectories
         if ((in3.ndt==1) or (in3.ndt==2)):
# reorder tagged particles
            if ((in3.nst==1) or (in3.nst==2)):
# determines list of tagged particles leaving this node: updates iholep
               minit3.mpfholes3(partt,tedges,numtp,iholep,2)
# iholep overflow
               if (iholep[0,0] < 0):
                  ntmax = -iholep[0,0]
                  ntmax = int(1.5*ntmax)
                  if (kstrt==1):
                     print "info:reallocating iholep:ntmax=", ntmax
                  iholep = numpy.empty((ntmax+1,2),int_type,'F')
                  minit3.mpfholes3(partt,tedges,numtp,iholep,2)
                  if (iholep[0,0] < 0):
                     if (kstrt==1):
                         print "iholep overflow: ntmax=", ntmax
                     mpplib3.ppexit(); exit(1)
# copies tagged particles: updates part
               mdiag3.mpcpytraj3(partt,part,tdiag,numtp)
# moves tagged electrons into original spatial region:
# updates part, numtp
               mppmod3.ipmove3(part,tedges,numtp,iholep,ny,nz,tmov,kstrt,
                               in3.nvpy,in3.nvpz,2,ierr)
               if (ierr[0] != 0):
                  mpplib3.ppexit(); exit(1)
# reorders tagged particles: updates partt
               mdiag3.mpordtraj3(part,partt,tedges,tdiag,numtp,ierr)
# collects distributed test particle data onto node 0
               mppmod3.mppartt3(partt,tdiag,numtp,in3.nvpy,in3.nvpz,ierr)
               if (ierr[0] != 0):
                  mpplib3.ppexit(); exit(1)
# write trajectory diagnostic output: updates ntrec
               if (kstrt==1):
                  mdiag3.dafwritetr3(partt,tdiag,iut,in3.ntrec)
                  partd[itt,:,:] = partt
                  itt += 1
            elif (in3.nst==3):
# calculate test particle distribution function and moments
               mdiag3.mpvdist3(partt,fvtp,fvmtp,tdiag,numtp,nvp,in3.nmv)
# write test particle diagnostic output: updates ntrec
               if (kstrt==1):
                  ws[0] = 0.0
                  mdiag3.dafwritefv3(fvmtp,fvtp,fetp,ws,tdiag,iut,
                                     in3.ntrec)

# phase space diagnostic
   if (in3.nts > 0):
      it = int(ntime/in3.nts)
      if (ntime==in3.nts*it):
# electron phase space diagnostic
         if ((in3.nds==1) or (in3.nds==3)):
# calculates velocity distribution in different regions of space:
# updates fvs
            mdiag3.mpvspdist3(ppart,kpic,fvs,noff,tdiag,in3.nmv,in3.mvx,
                              in3.mvy,in3.mvz,nyb)
# adjusts 3d velocity distribution in different regions of space:
# updates fvs
            mppmod3.mpadjfvs3(fvs,tdiag,noff,nyzp,in3.nmv,in3.mvy,
                              in3.mvz,nyb,in3.nvpy,nzbmx)
# write phase space diagnostic output: updates nserec
            if (in3.nserec > 0):
               mppmod3.mpwrfvsdata3(fvs,tdiag,nyb,nzb,nzbmx,iuse)
               in3.nserec += 1
# ion phase space
         if (in3.movion==1):
            if ((in3.nds==2) or (in3.nds==3)):
# calculates velocity distribution in different regions of space:
# updates fvsi
               mdiag3.mpvspdist3(pparti,kipic,fvsi,noff,tdiag,in3.nmv,
                                 in3.mvx,in3.mvy,in3.mvz,nyb)
# adjusts 3d velocity distribution in different regions of space:
# updates fvsi
               mppmod3.mpadjfvs3(fvsi,tdiag,noff,nyzp,in3.nmv,in3.mvy,
                                 in3.mvz,nyb,in3.nvpy,nzbmx)
# write phase space diagnostic output: updates nsirec
               if (in3.nsirec > 0):
                  mppmod3.mpwrfvsdata3(fvsi,tdiag,nyb,nzb,nzbmx,iusi)
                  in3.nsirec += 1

# push electrons with OpenMP:
# updates ppart and wke, and possibly ncl, ihole, irc
   wke[0] = 0.0
   mbpush3.wmpbpush3(ppart,exyze,bxyze,kpic,ncl,ihole,noff,nyzp,qbme,
                     in3.dt,in3.dt,in3.ci,wke,tpush,nx,ny,nz,in3.mx,
                     in3.my,in3.mz,mx1,myp1,ipbc,in3.popt,
                     in3.relativity,plist,irc)

# reorder electrons by tile with OpenMP and MPI
# updates: ppart, kpic, and irc and possibly ncl and ihole
   if (irc[0]==0):
      ompplib3.ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,in3.xtras,tsort,
                        tmov,kstrt,in3.nvpy,in3.nvpz,nx,ny,nz,in3.mx,
                        in3.my,in3.mz,npbmx,nbmaxp,mx1,myp1,mzp1,mxzyp1,
                        in3.popt,plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc[0]; irc[0] = 0

   while (irc2[0] != 0):
# ihole overflow
      if (irc2[0]==1):
         ntmaxp = int((1.0 + in3.xtras)*irc2[1])
         ihole = numpy.empty((2,ntmaxp+1,mxyzp1),int_type,'F')
         irc2[:] = 0
         ompplib3.ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,in3.xtras,
                           tsort,tmov,kstrt,in3.nvpy,in3.nvpz,nx,ny,nz,
                           in3.mx,in3.my,in3.mz,npbmx,nbmaxp,mx1,myp1,
                           mzp1,mxzyp1,in3.popt,False,irc2)
# ppart overflow
      elif (irc2[0]==4):
# restores electron coordinates from ppbuff: updates ppart, ncl
         msort3.mprstor3(ppart,ompplib3.ppbuff,ncl,ihole,tsort)
# copy ordered electrons to linear array: updates part
         mpush3.mpcopyout3(part,ppart,kpic,nt,irc)
         nppmx0 = int((1.0 + in3.xtras)*irc2[1])
         ppart = numpy.empty((idimp,nppmx0,mxyzp1),float_type,'F')
# copies unordered electrons to ordered array: updates ppart
         mpush3.mpcopyin3(part,ppart,kpic,irc)
         ompplib3.ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,in3.xtras,
                           tsort,tmov,kstrt,in3.nvpy,in3.nvpz,nx,ny,nz,
                           in3.mx,in3.my,in3.mz,npbmx,nbmaxp,mx1,myp1,
                           mzp1,mxzyp1,in3.popt,plist,irc2)

# sanity check for electrons
   if (in3.monitor > 0):
      mpush3.mpcheck3(ppart,kpic,noff,nyzp,nx,in3.mx,in3.my,in3.mz,mx1,
                      myp1,irc)

# push ions with OpenMP:
   if (in3.movion==1):
# updates pparti and wki, and possibly ncl, ihole, irc
      wki[0] = 0.0
      mbpush3.wmpbpush3(pparti,exyze,bxyze,kipic,ncl,ihole,noff,nyzp,
                        qbmi,in3.dt,in3.dt,in3.ci,wki,tpush,nx,ny,nz,
                        in3.mx,in3.my,in3.mz,mx1,myp1,ipbc,in3.popt,
                        in3.relativity,plist,irc)
      wki[0] = wki[0]*in3.rmass

# reorder ions by tile with OpenMP and MPI
# updates: pparti, kipic, and irc and possibly ncl and ihole
      if (irc[0]==0):
         ompplib3.ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,in3.xtras,
                           tsort,tmov,kstrt,in3.nvpy,in3.nvpz,nx,ny,nz,
                           in3.mx,in3.my,in3.mz,npbmx,nbmaxp,mx1,myp1,
                           mzp1,mxzyp1,in3.popt,plist,irc2)
      else:
         irc2[0] = 1; irc2[1] = irc[0]; irc[0] = 0

      while (irc2[0] != 0):
# ihole overflow
         if (irc2[0]==1):
            ntmaxp = int((1.0 + in3.xtras)*irc2[1])
            ihole = numpy.empty((2,ntmaxp+1,mxyzp1),int_type,'F')
            irc2[:] = 0
            ompplib3.ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,
                              in3.xtras,tsort,tmov,kstrt,in3.nvpy,
                              in3.nvpz,nx,ny,nz,in3.mx,in3.my,in3.mz,
                              npbmx,nbmaxp,mx1,myp1,mzp1,mxzyp1,
                              in3.popt,False,irc2)
# pparti overflow
         elif (irc2[0]==4):
# restores ion coordinates from ppbuff: updates pparti, ncl
            msort3.mprstor3(pparti,ompplib3.ppbuff,ncl,ihole,tsort)
# copy ordered electrons to linear array: updates part
            mpush3.mpcopyout3(part,pparti,kipic,nt,irc)
            nppmx1 = int((1.0 + in3.xtras)*irc2[1])
            pparti = numpy.empty((idimp,nppmx1,mxyzp1),float_type,'F')
# copies unordered ions to ordered array: updates pparti
            mpush3.mpcopyin3(part,pparti,kipic,irc)
            ompplib3.ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,
                              in3.xtras,tsort,tmov,kstrt,in3.nvpy,
                              in3.nvpz,nx,ny,nz,in3.mx,in3.my,in3.mz,
                              npbmx,nbmaxp,mx1,myp1,mzp1,mxzyp1,
                              in3.popt,plist,irc2)

# sanity check for ions
      if (in3.monitor > 0):
         mpush3.mpcheck3(pparti,kipic,noff,nyzp,nx,in3.mx,in3.my,in3.mz,
                         mx1,myp1,irc)

# start running simulation backwards:
# need to reverse time lag in leap-frog integration scheme
   if (in3.treverse==1):
      if (((ntime+1)==(nloop/2)) or ((ntime+1)==nloop)):
         in3.dt = -in3.dt
         ws[0] = 0.0
         mpush3.wmppush3zf(ppart,kpic,ncl,ihole,noff,nyzp,in3.dt,in3.ci,
                           ws,tpush,nx,ny,nz,in3.mx,in3.my,in3.mz,mx1,
                           myp1,ipbc,in3.popt,in3.relativity,plist,irc)
         ompplib3.ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,in3.xtras,
                           tsort,tmov,kstrt,in3.nvpy,in3.nvpz,nx,ny,nz,
                           in3.mx,in3.my,in3.mz,npbmx,nbmaxp,mx1,myp1,
                           mzp1,mxzyp1,in3.popt,plist,irc2)
         if (in3.movion==1):
            mpush3.wmppush3zf(pparti,kipic,ncl,ihole,noff,nyzp,in3.dt,
                              in3.ci,ws,tpush,nx,ny,nz,in3.mx,in3.my,
                              in3.mz,mx1,myp1,ipbc,in3.popt,
                              in3.relativity,plist,irc)
            ompplib3.ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,
                              in3.xtras,tsort,tmov,kstrt,in3.nvpy,
                              in3.nvpz,nx,ny,nz,in3.mx,in3.my,in3.mz,
                              npbmx,nbmaxp,mx1,myp1,mzp1,mxzyp1,
                              in3.popt,plist,irc2)

# energy diagnostic
   if (in3.ntw > 0):
      it = int(ntime/in3.ntw)
      if (ntime==in3.ntw*it):
         wef[0] = we[0] + wb[0]
         wtot[0] = wef[0]
         wtot[1] = wke[0]
         wtot[2] = wki[0]
         wtot[3] = wef[0] + wke[0]
         wtot[4] = we[0]
         wtot[5] = wf[0]
         wtot[6] = wb[0]
         mppmod3.mpdsum(wtot,tdiag)
         wke[0] = wtot[1]
         wki[0] = wtot[2]
         we[0] = wtot[4]
         wf[0] = wtot[5]
         wb[0] = wtot[6]
         wef[0] = we[0] + wb[0]
         ws[0] = wef[0] + wke[0] + wki[0]
         if (ntime==0):
            s[5] = ws[0]
         if (kstrt==1):
            print >> iuot, "Total Field, Kinetic and Total Energies:"
            if (in3.movion==0):
               iuot.write("%14.7e %14.7e %14.7e\n" % (wef[0],wke[0],
                          ws[0])) 
            else:
               iuot.write("%14.7e %14.7e %14.7e %14.7e\n" % (wef[0],
                          wke[0],wki[0],ws[0]))
            print >> iuot, "Electric(l,t) and Magnetic Energies:"
            iuot.write("%14.7e %14.7e %14.7e\n" % (we[0],wf[0],wb[0]))
# store energies in time history array
         wt[itw,:] = [wef[0],wke[0],wki[0],ws[0],we[0],wf[0],wb[0]]
         itw += 1
         s[0] += we[0]
         s[1] += wke[0]
         s[2] += wf[0]
         s[3] += wb[0]
         s[4] += wki[0]
         s[5] = min(s[5],float(ws[0]))
         s[6] = max(s[6],float(ws[0]))

# restart file
   if (in3.ntr > 0):
      n = ntime + 1
      it = int(n/in3.ntr)
      if (n==in3.ntr*it):
         dtimer(dtime,ltime,-1)
# write out basic restart file for electrostatic code
         s3.bwrite_restart3(part,ppart,pparti,qi,kpic,kipic,iur,n,
                            ntime0,irc)
# write out basic restart file for darwin code
         sd3.bwrite_drestart3(cus,wpm,q2m0,iur)
         dtimer(dtime,itime,1)
         tfield[0] += float(dtime)

ntime = ntime + 1

# loop time
dtimer(dtime,ltime,1)
tloop += float(dtime)

# * * * end main iteration loop * * *

if (kstrt==1):
   print >> iuot
   print >> iuot,"ntime,relativity,ndc=",ntime,",",in3.relativity,",",in3.ndc
   if (in3.treverse==1):
      print >> iuot, "treverse = ", in3.treverse
   print >> iuot, "MPI nodes nvpy, nvpz = ", in3.nvpy, in3.nvpz

# energy diagnostic
   if (in3.ntw > 0):
      s[5] = (s[6] - s[5])/wt[0,3]
      print >> iuot, "Energy Conservation = ", float(s[5])
      swe = s[0]; swke = s[1]; swf = s[2]; swb = s[3]
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

   print >> iuot
   print >> iuot, "initialization time = ", tinit
   print >> iuot, "deposit time = ", tdpost[0]
   print >> iuot, "current deposit time = ", tdjpost[0]
   print >> iuot, "current derivative deposit time = ", tdcjpost[0]
   tdpost[0] = tdpost[0] + tdjpost[0] +  tdcjpost[0]
   print >> iuot, "total deposit time = ", tdpost[0]
   print >> iuot, "guard time = ", tguard[0]
   print >> iuot, "solver time = ", tfield[0]
   print >> iuot, "field move time = ", tfmov[0]
   print >> iuot, "fft and transpose time = ", tfft[0], tfft[1]
   print >> iuot, "push time = ", tpush[0]
   print >> iuot, "particle move time = ", tmov[0]
   print >> iuot, "sort time = ", tsort[0]
   tfield[0] = tfield[0] + tguard[0] + tfft[0] + tfmov[0]
   print >> iuot, "total solver time = ", tfield[0]
   tsort[0] = tsort[0] + tmov[0]
   time = tdpost[0] + tpush[0] + tsort[0]
   print >> iuot, "total particle time = ", time
   print >> iuot, "total diagnostic time = ", tdiag[0]
   ws[0] = time + tfield[0] + tdiag[0]
   tloop = tloop - ws[0]
   print >> iuot, "total and additional time = ", ws[0], ",", tloop
   print >> iuot

   ws[0] = 1.0e+09/(float(nloop)*float(np+npi))
   print >> iuot, "Push Time (nsec) = ", tpush[0]*ws[0]
   print >> iuot, "Deposit Time (nsec) = ", tdpost[0]*ws[0]
   print >> iuot, "Sort Time (nsec) = ", tsort[0]*ws[0]
   print >> iuot, "Total Particle Time (nsec) = ", time*ws[0]
   print >> iuot

# reset parameters for final diagnostic metafile
if (in3.ntde > 0):
   in3.nderec -= 1
# potential diagnostic
if (in3.ntp > 0):
   in3.nprec -= 1; in3.ceng = affp
# longitudinal efield diagnostic
if (in3.ntel > 0):
   in3.nelrec -= 1; in3.ceng = affp
# electron current diagnostic
if (in3.ntje > 0):
   in3.njerec -= 1
# vector potential diagnostic
if (in3.nta > 0):
   in3.narec -= 1; in3.ceng = affp
# transverse efield diagnostic
if (in3.ntet > 0):
   in3.netrec -= 1; in3.ceng = affp
# magnetic field diagnostic
if (in3.ntb > 0):
   in3.nbrec -= 1; in3.ceng = affp
# fluid moments diagnostic
if (in3.ntfm > 0):
   if ((in3.ndfm==1) or (in3.ndfm==3)):
      in3.nferec -= 1
   if (in3.movion==1):
      if ((in3.ndfm==2) or (in3.ndfm==3)):
         in3.nfirec -= 1
   in3.ceng = affp
# velocity-space diagnostic
if (in3.ntv > 0):
   if ((in3.ndv==1) or (in3.ndv==3)):
      in3.nverec -= 1
   if (in3.movion==1):
      if ((in3.ndv==2) or (in3.ndv==3)):
         in3.nvirec -= 1
# trajectory diagnostic
if (in3.ntt > 0) :
   in3.ntrec -= 1
# phase space diagnostic
if (in3.nts > 0):
   if ((in3.nds==1) or (in3.nds==3)):
      in3.nserec -= 1
   if (in3.movion==1):
      if ((in3.nds==2) or (in3.nds==3)):
         in3.nsirec -= 1
if (in3.movion==1):
# ion density diagnostic
   if (in3.ntdi > 0):
      in3.ndirec -= 1
# ion current diagnostic
   if (in3.ntji > 0):
      in3.njirec -= 1
# write final diagnostic metafile
in3.writnml3(iudm)
#close(unit=iudm)
# close restart files
s3.close_restart3(iur,iur0)
# close output file
print >> iuot, " * * * q.e.d. * * *"
iuot.close()

mpplib3.ppexit()
