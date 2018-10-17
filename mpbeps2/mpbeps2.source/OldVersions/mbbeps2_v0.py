#-----------------------------------------------------------------------
# 2-1/2D Electromagnetic MPI/OpenMP PIC code
# written by Viktor K. Decyk, UCLA
import sys
import math
import numpy

#sys.path.append('./mpbeps2.source')
from libmpush2 import *
from libmbpush2 import *
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

# idimp = dimension of phase space = 5
# ipbc = particle boundary condition: 1 = periodic
idimp = 5; ipbc = 1
# idps = number of partition boundaries
idps = 2
# wke/wki/we = particle kinetic/electrostatic field energy
wke = numpy.zeros((1),float_type)
wki = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
# wf/wb = transverse electric field/magnetic field
wf = numpy.zeros((1),float_type)
wb = numpy.zeros((1),float_type)
zero = 0.0
# plist = (true,false) = list of particles leaving tiles found in push
plist = True

# declare scalars for standard code
ierr = numpy.zeros((1),int_type)
nt = numpy.empty((1),int_type)
ntime0 = numpy.zeros((1),int_type)
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
iur = 17; iur0 = 27; iscr = 99
iua = 13; iuet = 14; iub = 15; iuar = 16
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
mpplib2.ppinit2(idproc,nvp)
kstrt = idproc[0] + 1

# in2.nvpp = number of shared memory nodes (0=default)
#if (kstrt==1):
#   in2.nvpp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
omplib.init_omp(in2.nvpp)

# read namelists
if (kstrt==1):
# override default input data
   in2.emf = 1
   in2.relativity = 1
   in2.readnml2(iuin)
# override input data
   in2.idcode = 2
   in2.ndim = 3
# create string from idrun
   cdrun = str(in2.idrun)
# text output file
   fname = "output2." + cdrun
   iuot = open(fname,"w")

# broadcast namelists to other nodes
in2.sendnmls2()

# open graphics device
#     call IPLTCOMM(nplot)
#     if (kstrt==1) then
#        irc = open_pgraphs(nplot)
# set palette to color wheel
#        call STPALIT(2)
#     endif

# increase number of coordinates for particle tag
if (in2.ntt > 0):
   idimp += 1

# initialize scalars for standard code
# np = total number of particles in simulation
npxy =  float(in2.npx)*float(in2.npy)
npxyb =  float(in2.npxb)*float(in2.npyb)
np = npxy + npxyb
# npi = total number of ions in simulation
if (in2.movion > 0):
   npxyi = float(in2.npxi)*float(in2.npyi)
   npxybi = float(in2.npxbi)*float(in2.npybi)
   npi = npxyi + npxybi
# nx/ny = number of grid points in x/y direction
nx = int(math.pow(2,in2.indx)); ny = int(math.pow(2,in2.indy))
nxh = int(nx/2); nyh = max(1,int(ny/2))
nxe = nx + 2; nye = ny + 2; nxeh = nxe/2; nnxe = in2.ndim*nxe
nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
# mx1 = number of tiles in x direction
mx1 = int((nx - 1)/in2.mx + 1)
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(in2.tend/in2.dt + .0001); ntime = 0
qbme = in2.qme
affp = float(nx)*float(ny)/np
if (in2.movion==1):
   qbmi = in2.qmi/in2.rmass
   vtxi = in2.vtx/numpy.sqrt(in2.rmass*in2.rtempxi)
   vtyi = in2.vty/numpy.sqrt(in2.rmass*in2.rtempyi)
   vtzi = in2.vtz/numpy.sqrt(in2.rmass*in2.rtempzi)
   vtdxi = in2.vtdx/numpy.sqrt(in2.rmass*in2.rtempdxi)
   vtdyi = in2.vtdy/numpy.sqrt(in2.rmass*in2.rtempdyi)
   vtdzi = in2.vtdz/numpy.sqrt(in2.rmass*in2.rtempdzi)
dth = 0.0

# check if too many processors
if (nvp[0] > ny):
   if (kstrt==1):
      print "Too many nodes requested: ny, nvp=", ny, nvp[0]
   mpplib2.ppexit(); exit(1)

# initialize data for MPI code
edges = numpy.empty((idps),float_type,'F')
nyp = numpy.empty((1),int_type,'F')
noff = numpy.empty((1),int_type,'F')
nypmx = numpy.empty((1),int_type,'F')
nypmn = numpy.empty((1),int_type,'F')
# calculate partition variables: edges, nyp, noff, nypmx
# edges[0:1] = lower:upper boundary of particle partition
# nyp = number of primary (complete) gridpoints in particle partition
# noff = lowermost global gridpoint in particle partition
# nypmx = maximum size of particle partition, including guard cells
# nypmn = minimum value of nyp
# find new 1d partition for uniform density distribution
# minit2.mpdcomp2(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp)
# find new 1d partition from initial analytic distribution function
minit2.mpfedges2(edges,nyp,noff,in2.ampdy,in2.scaledy,in2.shiftdy,nypmx,
                 nypmn,ny,kstrt,nvp,ipbc,in2.ndprof,ierr)
# find new 1d partition from initial analytic distribution function
# where initial plasma density is confined to some subregion
# xmin = in2.dxmin*float(nx); xmax = in2.dxmax*float(nx)
# ymin = in2.dymin*float(ny); ymax = in2.dymax*float(ny)
#minit2.mpgfedges2(edges,nyp,noff,in2.ampdy,in2.scaledy,in2.shiftdy,
#                  ymin,ymax,nypmx,nypmn,ny,kstrt,nvp,in2.ndprof,ierr)
if (ierr[0] != 0):
   mpplib2.ppexit(); exit(1)

# check for unimplemented features
if ((plist) and (ipbc != 1)):
   print "ipbc != 1 and list = True not yet supported"
   plist = False
   print "plist reset to False"

# initialize additional scalars for MPI code
# kxp = number of complex grids in each field partition in x direction
kxp = int((nxh - 1)/nvp[0] + 1)
# kyp = number of complex grids in each field partition in y direction
kyp = int((ny - 1)/nvp[0] + 1)
# npmax/npimax = maximum number of electrons/ions in each partition
npmax = int((np/nvp[0])*1.25); npimax = int((npi/nvp[0])*1.25)
maxnp = max(npmax,npimax)
# myp1 = number of tiles in y direction
myp1 = int((nyp[0] - 1)/in2.my + 1); mxyp1 = mx1*myp1
# nterf = number of shifts required by field manager (0=search)
nterf = numpy.zeros((1),int_type)
# ntmax = size of iholep buffer for particles leaving node
ntmax = int(0.2*npmax)

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,maxnp),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nypmx[0]),float_type,'F')
# qi = ion charge density with guard cells
qi = numpy.empty((nxe,nypmx[0]),float_type,'F')
# cue = electron current density with guard cells
cue = numpy.empty((in2.ndim,nxe,nypmx[0]),float_type,'F')
# fxyze = smoothed electric field with guard cells
fxyze = numpy.empty((in2.ndim,nxe,nypmx[0]),float_type,'F')
# bxyze = smoothed /magnetic field with guard cells
bxyze = numpy.empty((in2.ndim,nxe,nypmx[0]),float_type,'F')
# qt = scalar charge density field array in fourier space
qt = numpy.empty((nye,kxp),complex_type,'F')
# cut = vector current density field arrays in fourier space
cut = numpy.empty((in2.ndim,nye,kxp),complex_type,'F')
# fxyt = vector electric field array in fourier space
fxyt = numpy.empty((in2.ndim,nye,kxp),complex_type,'F')
# bxyt = vector magnetic field array in fourier space
bxyt = numpy.empty((in2.ndim,nye,kxp),complex_type,'F')
# exyz = transverse electric field in fourier space
exyz = numpy.empty((in2.ndim,nye,kxp),complex_type,'F')
# bxyz = transverse magnetic field in fourier space
bxyz = numpy.empty((in2.ndim,nye,kxp),complex_type,'F')
# ffc = form factor array for poisson solver
ffc = numpy.empty((nyh,kxp),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhy),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyh),complex_type,'F')
# kpic = number of electrons in each tile
kpic = numpy.empty((mxyp1),int_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((8,mxyp1),int_type,'F')
# iholep = location of hole left in linear particle arrays
iholep = numpy.empty((ntmax+1),int_type,'F')
# kipic = number of ions in each tile
if (in2.movion==1):
# kipic = number of ions in each tile
   kipic = numpy.empty((mxyp1),int_type,'F')
# cui = ion current density with guard cells
   cui = numpy.empty((in2.ndim,nxe,nypmx[0]),float_type,'F')
# define dummy arrays for restart
elif (in2.movion==0):
   kipic = numpy.empty((0),int_type,'F')
   pparti = numpy.empty((0,0,0),float_type,'F')

# prepare fft tables
mfft2.mpfft2_init(mixup,sct,in2.indx,in2.indy)
# calculate form factors
mfield2.mppois2_init(ffc,in2.ax,in2.ay,affp,nx,ny,kstrt)
# initialize different ensemble of random numbers
if (in2.nextrand > 0):
   minit2.mnextran2(nextrand,in2.ndim,npmax+npimax)

# open restart files
if (kstrt==1):
   f2.open_restart2(nt,iur0,cdrun); iur = nt[0]

# new start
if (in2.nustrt==1):
# initialize electrons
   nps = 1
   npp[0] = 0
# background electrons
   if (npxy > 0.0):
# calculates initial electron co-ordinates with various density profiles
      minit2.mpfdistr2(part,npp,in2.ampdx,in2.scaledx,in2.shiftdx,
                       in2.ampdy,in2.scaledy,in2.shiftdy,in2.npx,
                       in2.npy,nx,ny,kstrt,nvp,ipbc,in2.ndprof,ierr)
# where initial plasma density is confined to some subregion
#     minit2.mpgfdistr2(part,npp,in2.ampdx,in2.scaledx,in2.shiftdx,
#                       in2.ampdy,in2.scaledy,in2.shiftdy,xmin,xmax,
#                       ymin,ymax,in2.npx,in2.npy,nx,ny,kstrt,nvp,
#                       in2.ndprof,ierr)
# initialize electron velocities or momenta
      if (ierr[0]==0):
# special cases
         if (in2.nvdist==3):
            minit2.mpvbdistr2h(part,nps,npp,in2.vtx,in2.vtz,in2.vx0,
                               in2.vz0,in2.omx,in2.omy,in2.omz,in2.npx,
                               in2.npy,kstrt,nvp,ierr)
# general cases
         else:
            minit2.wmpvdistr2h(part,nps,npp,in2.vtx,in2.vty,in2.vtz,
                               in2.vx0,in2.vy0,in2.vz0,in2.ci,in2.npx,
                               in2.npy,kstrt,nvp,in2.relativity,ierr)
# check for background electron initialization error
      if (ierr[0] != 0):
         mpplib2.ppexit(); exit(1)
# beam electrons
   if (npxyb > 0.0):
      nps = npp[0] + 1
# calculates initial electron co-ordinates with various density profiles
      minit2.mpfdistr2(part,npp,in2.ampdx,in2.scaledx,in2.shiftdx,
                       in2.ampdy,in2.scaledy,in2.shiftdy,in2.npxb,
                       in2.npyb,nx,ny,kstrt,nvp,ipbc,in2.ndprof,ierr)
# where initial plasma density is confined to some subregion
#     minit2.mpgfdistr2(part,npp,in2.ampdx,in2.scaledx,in2.shiftdx,
#                       in2.ampdy,in2.scaledy,in2.shiftdy,xmin,xmax,
#                       ymin,ymax,in2.npxb,in2.npyb,nx,ny,kstrt,nvp,
#                       in2.ndprof,ierr)
# initialize electron velocities or momenta
      if (ierr[0]==0):
# special cases
         if (in2.nvdist==3):
            minit2.mpvbdistr2h(part,nps,npp,in2.vtdx,in2.vtdz,in2.vdx,
                               in2.vdz,in2.omx,in2.omy,in2.omz,in2.npxb,
                               in2.npyb,kstrt,nvp,ierr)
# general cases
         else:
            minit2.wmpvdistr2h(part,nps,npp,in2.vtdx,in2.vtdy,in2.vtdz,
                               in2.vdx,in2.vdy,in2.vdz,in2.ci,in2.npxb,
                               in2.npyb,kstrt,nvp,in2.relativity,ierr)
# check for beam electron initialization error
      if (ierr[0] != 0):
         mpplib2.ppexit(); exit(1)

# check if any electrons are in the wrong node
   minit2.mpfholes2(part,edges,npp,iholep,in2.ndim,1)
# iholep overflow
   if (iholep[0] < 0):
      ntmax = -iholep[0]
      ntmax = int(1.5*ntmax)
      if (kstrt==1):
         print "reallocating electron iholep: ntmax=", ntmax
      iholep = numpy.empty((ntmax+1),int_type,'F') 
      minit2.mpfholes2(part,edges,npp,iholep,in2.ndim,1)
      if (iholep[0] < 0):
         if (kstrt==1):
            print "iholep overflow: ntmax=", ntmax
         mpplib2.ppexit(); exit(1)
# move electrons to correct node
      mppmod2.ipmove2(part,edges,npp,iholep,ny,tinit,kstrt,nvp,in2.ndim,
                      1,ierr)
      if (ierr[0] != 0):
         mpplib2.ppexit(); exit(1)

# find number of electrons in each of mx, my tiles: updates kpic, nppmx
   minit2.mpdblkp2(part,kpic,npp,noff,nppmx,in2.mx,in2.my,mx1,irc)

# allocate vector electron data
   nppmx0 = int((1.0 + in2.xtras)*nppmx)
   ntmaxp = int(in2.xtras*nppmx)
   npbmx = int(in2.xtras*nppmx)
   nbmaxp = int(0.25*mx1*npbmx)
# ppart = tiled electron particle array
   ppart = numpy.empty((idimp,nppmx0,mxyp1),float_type,'F')
# ihole = location/destination of each particle departing tile
   ihole = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
# copy ordered electron data for OpenMP
   mpush2.mpmovin2p(part,ppart,kpic,npp,noff,in2.mx,in2.my,mx1,irc)

# sanity check for electrons
   mpush2.mpcheck2(ppart,kpic,noff,nyp,nx,in2.mx,in2.my,mx1,irc)

# initialize background charge density: updates qi
   if (in2.movion==0):
      mpush2.mpset_pszero2(qi,tinit,in2.mx,in2.my,mx1,myp1)
# set background to negative of initial electron density
      if (in2.ionbkg==1):
         qmi = -in2.qme
         mpush2.mppost2(ppart,qi,kpic,noff,qmi,tinit,in2.mx,in2.my,mx1,
                        in2.popt)
         ompplib2.wmpaguard2(qi,nyp,tinit,nx,kstrt,nvp)

# initialize ions
   if (in2.movion==1):
      mpush2.mpset_pvzero2(cui,tinit,in2.mx,in2.my,mx1,myp1)
      nps = 1
      nppi[0] = 0
# background ions
      if (npxyi > 0.0):
# calculates initial ion co-ordinates with various density profiles
         minit2.mpfdistr2(part,nppi,in2.ampdxi,in2.scaledxi,
                          in2.shiftdxi,in2.ampdyi,in2.scaledyi,
                          in2.shiftdyi,in2.npxi,in2.npyi,nx,ny,kstrt,
                          nvp,ipbc,in2.ndprofi,ierr)
# where initial plasma density is confined to some subregion
#        minit2.mpgfdistr2(part,nppi,in2.ampdxi,in2.scaledxi,
#                          in2.shiftdxi,in2.ampdyi,in2.scaledyi,
#                          in2.shiftdyi,xmin,xmax,ymin,ymax,in2.npxi,
#                          in2.npyi,nx,ny,kstrt,nvp,in2.ndprofi,ierr)
# initialize ion velocities or momenta
         if (ierr[0]==0):
# special cases
            if (in2.nvdist==3):
               minit2.mpvbdistr2h(part,nps,nppi,vtxi,vtzi,in2.vxi0,
                                  in2.vzi0,in2.omx,in2.omy,in2.omz,
                                  in2.npxi,in2.npyi,kstrt,nvp,ierr)
# general cases
            else:
               minit2.wmpvdistr2h(part,nps,nppi,vtxi,vtyi,vtzi,in2.vxi0,
                                  in2.vyi0,in2.vzi0,in2.ci,in2.npxi,
                                  in2.npyi,kstrt,nvp,in2.relativity,
                                  ierr)
# check for background ion initialization error
         if (ierr[0] != 0):
            mpplib2.ppexit(); exit(1)
# beam ions
      if (npxybi > 0.0):
         nps = nppi[0] + 1
# calculates initial ion co-ordinates with various density profiles
         minit2.mpfdistr2(part,nppi,in2.ampdxi,in2.scaledxi,
                          in2.shiftdxi,in2.ampdyi,in2.scaledyi,
                          in2.shiftdyi,in2.npxbi,in2.npybi,nx,ny,kstrt,
                          nvp,ipbc,in2.ndprofi,ierr)
# where initial plasma density is confined to some subregion
#        minit2.mpgfdistr2(part,nppi,in2.ampdxi,in2.scaledxi,
#                          in2.shiftdxi,in2.ampdyi,in2.scaledyi,
#                          in2.shiftdyi,xmin,xmax,ymin,ymax,in2.npxbi,
#                          in2.npybi,nx,ny,kstrt,nvp,in2.ndprofi,ierr)
# initialize ion velocities or momenta
         if (ierr[0]==0):
# special cases
            if (in2.nvdist==3):
               minit2.mpvbdistr2h(part,nps,nppi,vtdxi,vtdzi,in2.vdxi,
                                  in2.vdzi,in2.omx,in2.omy,in2.omz,
                                  in2.npxbi,in2.npybi,kstrt,nvp,ierr)
# general cases
            else:
               minit2.wmpvdistr2h(part,nps,nppi,vtdxi,vtdyi,vtdzi,
                                  in2.vdxi,in2.vdyi,in2.vdzi,in2.ci,
                                  in2.npxbi,in2.npybi,kstrt,nvp,
                                  in2.relativity,ierr)
# check for beam ion initialization error
         if (ierr[0] != 0):
            mpplib2.ppexit(); exit(1)

# check if any ions are in the wrong node
      minit2.mpfholes2(part,edges,nppi,iholep,in2.ndim,1)
# iholep overflow
      if (iholep[0] < 0):
         ntmax = -iholep[0]
         ntmax = int(1.5*ntmax)
         if (kstrt==1):
            print "reallocating ion iholep: ntmax=", ntmax
         iholep = numpy.empty((ntmax+1),int_type,'F') 
         minit2.mpfholes2(part,edges,nppi,iholep,in2.ndim,1)
         if (iholep[0] < 0):
            if (kstrt==1):
               print "iholep overflow: ntmax=", ntmax
            mpplib2.ppexit(); exit(1)
# move ions to correct node
         mppmod2.ipmove2(part,edges,nppi,iholep,ny,tinit,kstrt,nvp,
                         in2.ndim,1,ierr)
         if (ierr[0] != 0):
            mpplib2.ppexit(); exit(1)

# find number of ions in each of mx, my tiles: updates kipic, nppmx
      minit2.mpdblkp2(part,kipic,nppi,noff,nppmx,in2.mx,in2.my,mx1,irc)

# allocate vector ion data
      nppmx1 = int((1.0 + in2.xtras)*nppmx)
      pparti = numpy.empty((idimp,nppmx1,mxyp1),float_type,'F')
      if ("ihole" not in globals()):
         ntmaxp = int(in2.xtras*nppmx)
         npbmx = int(in2.xtras*nppmx)
         nbmaxp = int(0.25*mxzy1*npbmx)
         ihole = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
# copy ordered ion data for OpenMP
      mpush2.mpmovin2p(part,pparti,kipic,nppi,noff,in2.mx,in2.my,mx1,
                       irc)

# sanity check for ions
      mpush2.mpcheck2(pparti,kipic,noff,nyp,nx,in2.mx,in2.my,mx1,irc)

# initialize transverse electromagnetic fields
   exyz.fill(numpy.complex(0.0,0.0))
   bxyz.fill(numpy.complex(0.0,0.0))
   cut.fill(numpy.complex(0.0,0.0))

# restart to continue a run which was interrupted
elif (in2.nustrt==2):
   print "nustrt = 2 not yet supported"
   mpplib2.ppexit(); exit(1)
#  if ((ntime+ntime0[0])> 0):
#     dth = 0.5*in2.dt
#  nstart = ntime + 1
# start a new run with data from a previous run
elif (in2.nustrt==0):
# read in basic restart file for electrostatic code
   if (in2.movion==1):
      mpush2.mpset_pvzero2(cui,tinit,in2.mx,in2.my,mx1,myp1)
# read first part of data:
# updates ntime, ntime0, part, npp, kpic, nppmx, ierr
   f2.bread_restart2a(part,kpic,tinit,kstrt,iur0,iscr,nt,ntime0,npp,
                      nppmx,noff,mx1,ierr)
   ntime = nt[0]
   if (ierr[0] != 0):
      mpplib2.ppexit(); exit(1)
# allocate vector electron data
   nppmx0 = int((1.0 + in2.xtras)*nppmx)
   ntmaxp = int(in2.xtras*nppmx)
   npbmx = int(in2.xtras*nppmx)
   nbmaxp = int(0.25*mxyp1*npbmx)
   ppart = numpy.empty((idimp,nppmx0,mxyp1),float_type,'F')
   ihole = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
# read second part of data:
# updates ppart, kpic, part, nppi, kipic, nppmx, ierr
   f2.bread_restart2b(part,ppart,kpic,kipic,tinit,kstrt,iur0,iscr,npp,
                      nppi,nppmx,noff,nyp,nx,mx1,ierr)
   if (ierr[0] != 0):
      mpplib2.ppexit(); exit(1)
# allocate vector ion data
   if (in2.movion==1):
      nppmx1 = int((1.0 + in2.xtras)*nppmx)
      pparti = numpy.empty((idimp,nppmx1,mxyp1),float_type,'F')
      if ("ihole" not in globals()):
         ntmaxp = int(in2.xtras*nppmx)
         npbmx = int(in2.xtras*nppmx)
         nbmaxp = int(0.25*mxyp1*npbmx)
         ihole = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
# read third part of data: updates pparti, kipic, qi, ierr
   f2.bread_restart2c(part,pparti,kipic,qi,tinit,kstrt,iur0,nt,ntime0,
                      nppi,noff,nyp,nx,mx1,ierr)
   ntime = nt[0]
   if (ierr[0] != 0):
      mpplib2.ppexit(); exit(1)
# read in basic restart file for electromagnetic code:
# updates exyz, bxyz
   fb2.bread_restart23(exyz,bxyz,tinit,kstrt,iur0,ierr)
   if (ierr[0] != 0):
      mpplib2.ppexit(); exit(1)
   if ((ntime+ntime0[0])> 0):
      dth = 0.5*in2.dt

# set magnitude of external magnetic field
omt = numpy.sqrt(in2.omx*in2.omx + in2.omy*in2.omy + in2.omz*in2.omz)

# kxps/kyps = actual grids used in field partitions in x/y direction
kxps = int(min(kxp,max(0,nxh-kxp*(kstrt-1))))
kyps = int(min(kyp,max(0,ny-kyp*(kstrt-1))))

# allocate diagnostic arrays
# reverse simulation at end back to start
if (in2.treverse==1):
   nloop = 2*nloop

# energy time history
if (in2.ntw > 0):
   mtw = int((nloop - 1)/in2.ntw + 1); itw = 0
# wt = energy time history array
   wt = numpy.zeros((mtw,7),float_type,'F')
# s = scratch array for energies
   s = numpy.zeros((7),double_type,'F')

# allocate scratch arrays for scalar fields
if ((in2.ntde > 0) or (in2.ntp > 0) or (in2.ntdi > 0)):
   sfieldc = numpy.empty((nye,kxp),complex_type,'F')
   sfield = numpy.empty((nxe,nypmx[0]),float_type,'F')

# allocate scratch arrays for vector fields
if ((in2.ntel>0) or (in2.ntje>0) or (in2.nta>0) or (in2.ntet>0)
or (in2.ntb>0) or (in2.ntar>0) or (in2.ntji>0)):
   vfieldc = numpy.empty((in2.ndim,nye,kxp),complex_type,'F')
   vfield = numpy.empty((in2.ndim,nxe,nypmx[0]),float_type,'F')

# initialize electron density diagnostic
if (in2.ntde > 0):
   in2.modesxde = int(min(in2.modesxde,nxh))
   in2.modesyde = int(min(in2.modesyde,nyh))
   modesy2de = 2*in2.modesyde - 1
   modesxpd = int(min(in2.modesxde,kxp))
# denet = store selected fourier modes for electron density
   denet = numpy.empty((modesy2de,modesxpd),complex_type,'F')
# open file
   if (kstrt==1):
# open file for complex data: updates nderec and possibly iude
#     fdename = "denek2." + cdrun
#     in2.fdename[:] = fdename
#     if (in2.nderec==0):
#        mdiag2.dafopenc2(denet,iude,in2.nderec,fdename)
# open file for real data: updates nderec and possibly iude
      fdename = "dener2." + cdrun
      in2.fdename[:] = fdename
      if (in2.nderec==0):
         mdiag2.dafopen2(sfield,nx,kyp,iude,in2.nderec,fdename)

# initialize ion density diagnostic
if (in2.movion==1):
   if (in2.ntdi > 0):
      in2.modesxdi = int(min(in2.modesxdi,nxh))
      in2.modesydi = int(min(in2.modesydi,nyh))
      modesy2di = 2*in2.modesydi - 1
      modesxpd = int(min(in2.modesxdi,kxp))
# denit = store selected fourier modes for ion density
      denit = numpy.empty((modesy2di,modesxpd),complex_type,'F')
# open file
      if (kstrt==1):
# open file for complex data: updates ndirec and possibly iudi
#        fdiname = "denik2." + cdrun
#        in2.fdiname[:] = fdiname
#        if (in2.ndirec==0):
#           mdiag2.dafopenc2(denit,iudi,in2.ndirec,fdiname)
# open file for real data: updates ndirec and possibly iudi
         fdiname = "denir2." + cdrun
         in2.fdiname[:] = fdiname
         if (in2.ndirec==0):
            mdiag2.dafopen2(sfield,nx,kyp,iudi,in2.ndirec,fdiname)

# initialize potential diagnostic
if (in2.ntp > 0):
   in2.modesxp = int(min(in2.modesxp,nxh))
   in2.modesyp = int(min(in2.modesyp,nyh))
   modesy2p = 2*in2.modesyp - 1
   modesxpd = int(min(in2.modesxp,kxp))
# pott = store selected fourier modes for potential
   pott = numpy.empty((modesy2p,modesxpd),complex_type,'F')
# open file
   if (kstrt==1):
# open file for complex data: updates nprec and possibly iup
#     fpname = "potk2." + cdrun
#     in2.fpname[:] = fpname
#     if (in2.nprec==0):
#        mdiag2.dafopenc2(pott,iup,in2.nprec,fpname)
# open file for real data: updates nprec and possibly iup
      fpname = "potr2." + cdrun
      in2.fpname[:] = fpname
      if (in2.nprec==0):
         mdiag2.dafopen2(sfield,nx,kyp,iup,in2.nprec,fpname)

# initialize longitudinal efield diagnostic
if (in2.ntel > 0):
   in2.modesxel = int(min(in2.modesxel,nxh))
   in2.modesyel = int(min(in2.modesyel,nyh))
   modesy2el = 2*in2.modesyel - 1
   modesxpd = int(min(in2.modesxel,kxp))
# elt = store selected fourier modes for longitudinal efield
   elt = numpy.empty((in2.ndim,modesy2el,modesxpd),complex_type,'F')
# open file
   if (kstrt==1):
# open file for complex data: updates nelrec and possibly iuel
#     felname = "elk2." + cdrun
#     in2.felname[:] = felname
#     if (in2.nelrec==0):
#        mdiag2.dafopenvc2(elt,iuel,in2.nelrec,felname)
# open file for real data: updates nelrec and possibly iuel
      felname = "elr2." + cdrun
      in2.felname[:] = felname
      if (in2.nelrec==0):
         mdiag2.dafopenv2(vfield,nx,kyp,iuel,in2.nelrec,felname)

# initialize electron current density diagnostic
if (in2.ntje > 0):
   in2.modesxje = int(min(in2.modesxje,nxh))
   in2.modesyje = int(min(in2.modesyje,nyh))
   modesy2je = 2*in2.modesyje - 1
   modesxpd = int(min(in2.modesxje,kxp))
# curet = store selected fourier modes for electron current density
   curet = numpy.empty((in2.ndim,modesy2je,modesxpd),complex_type,'F')
# open file
   if (kstrt==1):
# open file for complex data: updates njerec and possibly iuje
#     fjename = "curek2." + cdrun
#     in2.fjename[:] = fjename
#     if (in2.njerec==0):
#        mdiag2.dafopenvc2(curet,iuje,in2.njerec,fjename)
# open file for real data: updates nelrec and possibly iuel
      fjename = "curer2." + cdrun
      in2.fjename[:] = fjename
      if (in2.njerec==0):
         mdiag2.dafopenv2(vfield,nx,kyp,iuje,in2.njerec,fjename)

# initialize ion current density diagnostic
if (in2.movion==1):
   if (in2.ntji > 0):
      in2.modesxji = int(min(in2.modesxji,nxh))
      in2.modesyji = int(min(in2.modesyji,nyh))
      modesy2ji = 2*in2.modesyji - 1
      modesxpd = int(min(in2.modesxji,kxp))
# curit = store selected fourier modes for ion current density
      curit = numpy.empty((in2.ndim,modesy2ji,modesxpd),complex_type,
                          'F')
# open file
      if (kstrt==1):
# open file for complex data: updates njirec and possibly iuji
#     fjiname = "curik2." + cdrun
#     in2.fjiname[:] = fjiname
#     if (in2.njirec==0):
#        mdiag2.dafopenvc2(curit,iuji,in2.njirec,fjiname)
# open file for real data: updates njirec and possibly iuji
         fjiname = "curir2." + cdrun
         in2.fjiname[:] = fjiname
         if (in2.njirec==0):
            mdiag2.dafopenv2(vfield,nx,kyp,iuji,in2.njirec,fjiname)

# initialize radiative vector potential diagnostic
if (in2.ntar > 0):
   in2.modesxar = int(min(in2.modesxar,nxh))
   in2.modesyar = int(min(in2.modesyar,nyh))
   modesy2ar = 2*in2.modesyar - 1
   modesxpd = int(min(in2.modesxar,kxp))
# vpotr = store selected fourier modes for radiative vector potential
   vpotr = numpy.empty((in2.ndim,modesy2ar,modesxpd),complex_type,'F')
# oldcut = previous current density
   oldcut = numpy.empty((in2.ndim,nye,kxp),complex_type,'F')
# open file
   if (kstrt==1):
# open file for complex data: updates narrec and possibly iuar
#     farname = "vpotrk2." + cdrun
#     in2.farname[:] = farname
#     if (in2.narrec==0):
#        mdiag2.dafopenvc2(vpotr,iuar,in2.narrec,farname)
# open file for real data: updates narrec and possibly iuar
      farname = "vpotrr2." + cdrun
      in2.farname[:] = farname
      if (in2.narrec==0):
         mdiag2.dafopenv2(vfield,nx,kyp,iuar,in2.narrec,farname)

# initialize vector potential diagnostic
if (in2.nta > 0):
   in2.modesxa = int(min(in2.modesxa,nxh))
   in2.modesya = int(min(in2.modesya,nyh))
   modesy2a = 2*in2.modesya - 1
   modesxpd = int(min(in2.modesxa,kxp))
# vpott = store selected fourier modes for vector potential
   vpott = numpy.empty((in2.ndim,modesy2a,modesxpd),complex_type,'F')
# open file
   if (kstrt==1):
# open file for complex data: updates narec and possibly iua
#     faname = "vpotk2." + cdrun
#     in2.faname[:] = faname
#     if (in2.narec==0):
#        mdiag2.dafopenvc2(vpott,iua,in2.narec,faname)
# open file for real data: updates narec and possibly iua
      faname = "vpotr2." + cdrun
      in2.faname[:] = faname
      if (in2.narec==0):
         mdiag2.dafopenv2(vfield,nx,kyp,iua,in2.narec,faname)

# initialize transverse efield diagnostic
if (in2.ntet > 0):
   in2.modesxet = int(min(in2.modesxet,nxh))
   in2.modesyet = int(min(in2.modesyet,nyh))
   modesy2et = 2*in2.modesyet - 1
   modesxpd = int(min(in2.modesxet,kxp))
# ett = store selected fourier modes for transverse efield
   ett = numpy.empty((in2.ndim,modesy2et,modesxpd),complex_type,'F')
# open file
   if (kstrt==1):
# open file for complex data: updates netrec and possibly iuet
#     fetname = "etk2." + cdrun
#     in2.fetname[:] = fetname
#     if (in2.netrec==0):
#        mdiag2.dafopenvc2(ett,iuet,in2.netrec,fetname)
# open file for real data: updates netrec and possibly iuet
      fetname = "etr2." + cdrun
      in2.fetname[:] = fetname
      if (in2.netrec==0):
         mdiag2.dafopenv2(vfield,nx,kyp,iuet,in2.netrec,fetname)

# initialize magnetic field diagnostic
if (in2.ntb > 0):
   in2.modesxb = int(min(in2.modesxb,nxh))
   in2.modesyb = int(min(in2.modesyb,nyh))
   modesy2b = 2*in2.modesyb - 1
   modesxpd = int(min(in2.modesxb,kxp))
# bt = store selected fourier modes for magnetic field
   bt = numpy.empty((in2.ndim,modesy2b,modesxpd),complex_type,'F')
# open file
   if (kstrt==1):
# open file for complex data: updates nbrec and possibly iub
#     fbname = "bk2." + cdrun
#     in2.fbname[:] = fbname
#     if (in2.nbrec==0):
#        mdiag2.dafopenvc2(bt,iub,in2.nbrec,fbname)
# open file for real data: updates nbrec and possibly iub
      fbname = "br2." + cdrun
      in2.fbname[:] = fbname
      if (in2.nbrec==0):
         mdiag2.dafopenv2(vfield,nx,kyp,iub,in2.nbrec,fbname)

# initialize fluid moments diagnostic
if (in2.ntfm > 0):
   in2.nprd = 0
   if (in2.npro==1):
      in2.nprd = 1
   elif (in2.npro==2):
      in2.nprd = 4
   elif (in2.npro==3):
      in2.nprd = 10
   elif (in2.npro==4):
      in2.nprd = 14
# electron moments
   if ((in2.ndfm==1) or (in2.ndfm==3)):
      fmse = numpy.empty((in2.nprd,nxe,nypmx[0]),float_type,'F')
# open file for real data: updates nferec and possibly iufe
      if (kstrt==1):
         ffename = "fmer2." + cdrun
         in2.ffename[:] = ffename
         if (in2.nferec==0):
            mdiag2.dafopenv2(fmse,nx,kyp,iufe,in2.nferec,ffename)
# ion moments
   if (in2.movion==1):
      if ((in2.ndfm==2) or (in2.ndfm==3)):
# fmsi = ion fluid moments
         fmsi = numpy.empty((in2.nprd,nxe,nypmx[0]),float_type,'F')
# open file for real data: updates nfirec and possibly iufi
         if (kstrt==1):
            ffiname = "fmir2." + cdrun
            in2.ffiname[:] = ffiname
            if (in2.nfirec==0):
               mdiag2.dafopenv2(fmsi,nx,kyp,iufi,in2.nfirec,ffiname)

# initialize velocity-space diagnostic
if (in2.ntv > 0):
   in2.nfvd = 0; in2.nfed = 0
   if ((in2.nvft==1) or (in2.nvft==3)):
      in2.nfvd = in2.ndim
   elif ((in2.nvft==4) or (in2.nvft==5)):
      in2.nfvd = 2
   if ((in2.nvft==2) or (in2.nvft==3) or (in2.nvft==5)):
      in2.nfed = 1
   nmv21 = 2*in2.nmv + 1
   mtv = int((nloop - 1)/in2.ntv) + 1; itv = 0
   eci = in2.ci
   if (in2.relativity==0):
      eci = 0.0
   wk = numpy.zeros((1),float_type)
# electron velocity diagnostic
   if ((in2.ndv==1) or (in2.ndv==3)):
# estimate maximum velocity or momentum
      ws[0] = 0.0
      if (npxy > 0.0):
         ws[0] = 4.0*in2.vtx+abs(in2.vx0)
         ws[0] = max(ws[0],4.0*in2.vty+abs(in2.vy0))
         ws[0] = max(ws[0],4.0*in2.vtz+abs(in2.vz0))
      if (npxyb > 0.0):
         ws[0] = max(ws[0],4.0*in2.vtdx+abs(in2.vdx))
         ws[0] = max(ws[0],4.0*in2.vtdy+abs(in2.vdy))
         ws[0] = max(ws[0],4.0*in2.vtdz+abs(in2.vdz))
# fvm = electron vdrift, vth, entropy for global distribution
      fvm = numpy.zeros((in2.ndim,3),float_type,'F')
# sfv = electron/ion velocity distribution functions in tile
      sfv = numpy.empty((nmv21+1,in2.ndim,mxyp1),float_type,'F')
# fv = global electron velocity distribution functions
      fv = numpy.empty((nmv21+1,in2.nfvd),float_type,'F')
# fe = global electron energy distribution functions
      fe = numpy.empty((nmv21+1,in2.nfed),float_type,'F')
# open file for electron velocity data: updates nverec and possibly iuve
      if (kstrt==1):
         fvename = "fve2." + cdrun
         in2.fvename[:] = fvename
         if (in2.nverec==0):
            mdiag2.dafopenfv2(fvm,fv,fe,wk,iuve,in2.nverec,fvename)
# cartesian distribution
   if ((in2.nvft==1) or (in2.nvft==3)):
# fvtm = time history of electron vdrift, vth, and entropy for global
# distribution
      fvtm = numpy.zeros((mtv,in2.ndim,3),float_type,'F')
# set velocity or momentum scale
      fv[nmv21,:] = 2.0*ws[0]
# cylindrical distribution
      if ((in2.nvft==4) or (in2.nvft==5)):
# set velocity or momentum scale
         fv[nmv21,:] = 2.0*ws[0]
# energy distribution
      if ((in2.nvft==2) or (in2.nvft==3) or (in2.nvft==5)):
# set energy scale for electrons
         ws[0] = ws[0]*ws[0]
         fe[nmv21,:] = ws[0]/(1.0 + numpy.sqrt(1.0 + ws[0]*eci*eci))
# ion velocity diagnostic
   if (in2.movion==1):
      if ((in2.ndv==2) or (in2.ndv==3)):
# fvmi = ion vdrift, vth, entropy for global distribution
         fvmi = numpy.zeros((in2.ndim,3),float_type,'F')
# fvi = global ion velocity distribution functions
         fvi = numpy.empty((nmv21+1,in2.nfvd),float_type,'F')
# fei = global ion energy distribution functions
         fei = numpy.empty((nmv21+1,in2.nfed),float_type,'F')
# estimate maximum ion velocity or momentum
         ws[0] = 0.0
         if (npxyi > 0.0):
            ws[0] = 4.0*vtxi+abs(in2.vxi0)
            ws[0] = max(ws[0],4.0*vtyi+abs(in2.vyi0))
            ws[0] = max(ws[0],4.0*vtzi+abs(in2.vzi0))
         if (npxybi > 0.0):
            ws[0] = max(ws[0],4.0*vtdxi+abs(in2.vdxi))
            ws[0] = max(ws[0],4.0*vtdyi+abs(in2.vdyi))
            ws[0] = max(ws[0],4.0*vtdzi+abs(in2.vdzi))
# open file for ion velocity data: updates nvirec and possibly iuvi
         if (kstrt==1):
            fviname = "fvi2." + cdrun
            in2.fviname[:] = fviname
            if (in2.nvirec==0):
               mdiag2.dafopenfv2(fvmi,fvi,fei,wk,iuvi,in2.nvirec,
                                 fviname)
# cartesian distribution
         if ((in2.nvft==1) or (in2.nvft==3)):
# fvtmi = time history of ion vdrift, vth, and entropy for global
# distribution
            fvtmi = numpy.zeros((mtv,in2.ndim,3),float_type,'F')
# set velocity or momentum scale
            fvi[nmv21,:] = 2.0*ws[0]
# cylindrical distribution
         if ((in2.nvft==4) or (in2.nvft==5)):
# set velocity or momentum scale
            fvi[nmv21,:] = 2.0*ws[0]
# energy distribution
         if ((in2.nvft==2) or (in2.nvft==3) or (in2.nvft==5)):
# set energy scale for ions
            ws[0] = ws[0]*ws[0]
            fei[nmv21,0] = ws[0]/(1.0 + numpy.sqrt(1.0 + ws[0]*eci*eci))

# initialize trajectory diagnostic
if (in2.ntt > 0):
   if ((in2.ndt==2) and (in2.movion==0)):
      in2.ndt = 0
   if ((in2.ndt==1) or (in2.ndt==2)):
# iprobt = scratch array
      iprobt = numpy.empty((in2.nprobt),numpy.int32)
# tedges[0:1] = lower:upper y boundary of particle tags
      tedges = numpy.empty((idps),float_type,'F')
# electron trajectories
   if (in2.ndt==1):
# set electron tags: updates nprobt, tedges, ppart and possibly iprobt
      mdiag2.psetptraj2(ppart,tedges,kpic,iprobt,kstrt,in2.nst,in2.vtx,
                        in2.vtsx,in2.dvtx,np,in2.nprobt,irc)
# estimate maximum electron velocity or momentum
      if (in2.nst==3):
         ws[0] = 0.0
         if (npxy > 0.0):
            ws[0] = 4.0*in2.vtx+abs(in2.vx0)
            ws[0] = max(ws[0],4.0*in2.vty+abs(in2.vy0))
            ws[0] = max(ws[0],4.0*in2.vtz+abs(in2.vz0))
         if (npxyb > 0.0):
            ws[0] = max(ws[0],4.0*in2.vtdx+abs(in2.vdx))
            ws[0] = max(ws[0],4.0*in2.vtdy+abs(in2.vdy))
            ws[0] = max(ws[0],4.0*in2.vtdz+abs(in2.vdz))
# ion trajectories
   elif (in2.ndt==2):
# set ion tags: updates nprobt, tedges, pparti and possibly iprobt
      mdiag2.psetptraj2(pparti,tedges,kipic,iprobt,kstrt,in2.nst,
                        vtxi,in2.vtsx,in2.dvtx,npi,in2.nprobt,irc)
# estimate maximum ion velocity or momentum
      if (in2.nst==3):
         ws[0] = 0.0
         if (npxyi > 0.0):
            ws[0] = 4.0*vtxi+abs(in2.vxi0)
            ws[0] = max(ws[0],4.0*vtyi+abs(in2.vyi0))
            ws[0] = max(ws[0],4.0*vtzi+abs(in2.vzi0))
         if (npxybi > 0.0):
            ws[0] = max(ws[0],4.0*vtdxi+abs(in2.vdxi))
            ws[0] = max(ws[0],4.0*vtdyi+abs(in2.vdyi))
            ws[0] = max(ws[0],4.0*vtzi+abs(in2.vzi0))
# electron or ion trajectories
   if ((in2.ndt==1) or (in2.ndt==2)):
      if (in2.nprobt > 16777215):
         print "nprobt overflow = ", in2.nprobt
         mpplib2.ppexit(); exit(1)
      in2.ndimp = idimp
# partt = particle trajectories tracked
      partt = numpy.empty((idimp,in2.nprobt),float_type,'F')
      ftname = "tr2." + cdrun
      in2.ftname[:] = ftname
      if ((in2.nst==1) or (in2.nst==2)):
         it = int((nloop - 1)/in2.ntt + 1); itt = 0
# partd = trajectory time history array
         partd = numpy.zeros((it,idimp,in2.nprobt),float_type,'F')
# open file for trajectory data: updates ntrec and possibly iut
         if (kstrt==1):
            if (in2.ntrec==0):
               mdiag2.dafopentr2(partt,iut,in2.ntrec,ftname)
      elif (in2.nst==3):
# fvtp = velocity distribution function for test particles
         fvtp = numpy.empty((2*in2.nmv+2,in2.ndim),float_type,'F')
# fvmtp = vdrift, vth, and entropy for test particles
         fvmtp = numpy.empty((in2.ndim,3),float_type,'F')
         fetp = numpy.empty((in2.ndim,0),float_type,'F')
         fvtp[2*in2.nmv+1,:] = 2.0*ws[0]
# open file for test particle diagnostic: updates ntrec and possibly iut
         if (kstrt==1):
            if (in2.ntrec==0):
               ws[0] = 0.0
               mdiag2.dafopenfv2(fvmtp,fvtp,fetp,ws,iut,in2.ntrec,
                                 ftname)

# initialize phase space diagnostic
if (in2.nts > 0):
   nmv21 = 2*in2.nmv + 1
   in2.mvx = min(in2.mvx,nx); in2.mvy = min(in2.mvy,nypmn[0])
   in2.nsxb = int((nx - 1)/in2.mvx + 1)
   in2.nsyb = int((ny - 1)/in2.mvy + 1)
   nyb = int(float(noff[0] - 1)/float(in2.mvy))
   nyb = int((noff[0] + nyp[0] - 1)/in2.mvy) - nyb
   if (kstrt==1):
      nyb += 1
   itot[0] = nyb
   mppmod2.mpimax(itot[:1],tinit)
   nybmx = itot[0]
# electron phase space diagnostic
   if ((in2.nds==1) or (in2.nds==3)):
# estimate maximum electron velocity or momentum
      ws[0] = 0.0
      if (npxy > 0.0):
         ws[0] = 4.0*in2.vtx+abs(in2.vx0)
         ws[0] = max(ws[0],4.0*in2.vty+abs(in2.vy0))
         ws[0] = max(ws[0],4.0*in2.vtz+abs(in2.vz0))
      if (npxyb > 0.0):
         ws[0] = max(ws[0],4.0*in2.vtdx+abs(in2.vdx))
         ws[0] = max(ws[0],4.0*in2.vtdy+abs(in2.vdy))
         ws[0] = max(ws[0],4.0*in2.vtdz+abs(in2.vdz))
# fvs = global electron phase space distribution functions
      fvs = numpy.empty((nmv21+1,in2.ndim,in2.nsxb,nyb+1),float_type,
                        'F')
      fvs[nmv21,:,0,0] = 1.25*ws[0]
# open file for electron phase space data:
# updates nserec and possibly iuse
# opens a new fortran unformatted stream file
      if (in2.nserec==0):
         if (kstrt==1):
            fsename = "pse2." + cdrun
            in2.fsename[:] = fsename
            iuse =  mdiag2.get_funit(iuse)
            f2.fnopens2(iuse,fsename)
            in2.nserec = 1
# ion phase space diagnostic
   if (in2.movion==1):
      if ((in2.nds==2) or (in2.nds==3)):
# estimate maximum ion velocity or momentum
         ws[0] = 0.0
         if (npxyi > 0.0):
            ws[0] = 4.0*vtxi+abs(in2.vxi0)
            ws[0] = max(ws[0],4.0*vtyi+abs(in2.vyi0))
            ws[0] = max(ws[0],4.0*vtzi+abs(in2.vzi0))
         if (npxybi > 0.0):
            ws[0] = max(ws[0],4.0*vtdxi+abs(in2.vdxi))
            ws[0] = max(ws[0],4.0*vtdyi+abs(in2.vdyi))
            ws[0] = max(ws[0],4.0*vtdzi+abs(in2.vdzi))
# fvsi = global ion phase space distribution functions
         fvsi = numpy.empty((nmv21+1,in2.ndim,in2.nsxb,nyb+1),
                            float_type,'F')
         fvsi[nmv21,:,0,0] = 1.25*ws[0]
# open file for ion phase space data:
# updates nsirec and possibly iusi
# opens a new fortran unformatted stream file
         if (in2.nsirec==0):
            if (kstrt==1):
               fsiname = "psi2." + cdrun
               in2.fsiname[:] = fsiname
               iusi =  mdiag2.get_funit(iusi)
               f2.fnopens2(iusi,fsiname)
               in2.nsirec = 1

# initialization time
dtimer(dtime,itime,1)
tinit = 0.0
tinit += float(dtime)
# start timing loop
dtimer(dtime,ltime,-1)

if (in2.dt > 0.45*in2.ci):
   if (kstrt==1):
      print "Warning: Courant condition may be exceeded!"

if (kstrt==1):
   print >> iuot, "program mpbbeps2"

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):

# print time step
   if (kstrt==1):
      if (in2.ntw > 0):
         it = int(ntime/in2.ntw)
         if (ntime==in2.ntw*it):
            print >> iuot, "ntime = ", ntime

# save previous current in fourier space for radiative vector potential
   if (in2.ntar > 0):
      it = int(ntime/in2.ntar)
      if (ntime==in2.ntar*it):
         oldcut[:,:,:] = cut

# fluid moments diagnostic
   if (in2.ntfm > 0):
      it = int(ntime/in2.ntfm)
      if (ntime==in2.ntfm*it):
# calculate electron fluid moments
         if ((in2.ndfm==1) or (in2.ndfm==3)):
            mpush2.mpset_pvzero2(fmse,tdiag,in2.mx,in2.my,mx1,myp1)
            mdiag2.wmprofx23(ppart,fmse,kpic,noff,in2.ci,tdiag,in2.npro,
                             in2.mx,in2.my,mx1,in2.relativity)
# add guard cells with OpenMP: updates fmse
            ompplib2.wmpnacguard2(fmse,nyp,tdiag,nx,kstrt,nvp)
# moves vector grid fmse from non-uniform to uniform partition
            isign = -1
            mppmod2.mpfnmove2(fmse,noff,nyp,isign,tdiag,kyp,ny,kstrt,
                              nvp,nterf,irc)
            if (irc[0] != 0):
                  mpplib2.ppexit(); exit(1)
# calculates fluid quantities from fluid moments: updates fmsi
            mdiag2.mpfluidqs23(fmse,tdiag,in2.npro,nx,ny,kstrt,kyp)
# write real space diagnostic output: updates nferec
            mppmod2.mpvwrite2(fmse,tdiag,nx,ny,kyp,iufe,in2.nferec)
# calculate ion fluid moments
         if (in2.movion==1):
            if ((in2.ndfm==2) or (in2.ndfm==3)):
               mpush2.mpset_pvzero2(fmsi,tdiag,in2.mx,in2.my,mx1,myp1)
               mdiag2.wmprofx23(pparti,fmsi,kipic,noff,in2.ci,tdiag,
                                in2.npro,in2.mx,in2.my,mx1,
                                in2.relativity)
               fmsi[:,:,:] = in2.rmass*fmsi
# add guard cells with OpenMP: updates fmsi
               ompplib2.wmpnacguard2(fmsi,nyp,tdiag,nx,kstrt,nvp)
# moves vector grid fmsi from non-uniform to uniform partition
               isign = -1
               mppmod2.mpfnmove2(fmsi,noff,nyp,isign,tdiag,kyp,ny,kstrt,
                                 nvp,nterf,irc)
               if (irc[0] != 0):
                  mpplib2.ppexit(); exit(1)
# calculates fluid quantities from fluid moments: updates fmsi
               mdiag2.mpfluidqs23(fmsi,tdiag,in2.npro,nx,ny,kstrt,kyp)
# write real space diagnostic output: updates nfirec
               mppmod2.mpvwrite2(fmsi,tdiag,nx,ny,kyp,iufi,in2.nfirec)

# deposit electron current with OpenMP:
# updates ppart and cue, and possibly ncl, ihole, irc
   mpush2.mpset_pvzero2(cue,tdjpost,in2.mx,in2.my,mx1,myp1)
   mcurd2.wmpdjpost2(ppart,cue,kpic,ncl,ihole,noff,nyp,in2.qme,dth,
                     in2.ci,tdjpost,nx,ny,in2.mx,in2.my,mx1,ipbc,
                     in2.popt,in2.relativity,plist,irc)

# add guard cells with OpenMP: updates cue
   ompplib2.wmpnacguard2(cue,nyp,tguard,nx,kstrt,nvp)

# reorder electrons by tile with OpenMP and MPI
# updates: ppart, kpic, and irc and possibly ncl and ihole
   if (irc[0]==0):
      ompplib2.ompmove2(ppart,kpic,ncl,ihole,noff,nyp,in2.xtras,tsort,
                        tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,npbmx,nbmaxp,
                        mx1,in2.popt,plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc[0]; irc[0] = 0

   while (irc2[0] != 0):
# ihole overflow
      if (irc2[0]==1):
         ntmaxp = int((1.0 + in2.xtras)*irc2[1])
         ihole = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
         irc2[:] = 0
         ompplib2.ompmove2(ppart,kpic,ncl,ihole,noff,nyp,in2.xtras,
                           tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                           npbmx,nbmaxp,mx1,in2.popt,False,irc2)
# ppart overflow
      elif (irc2[0]==4):
# restores electron coordinates from ppbuff: updates ppart, ncl
         msort2.mprstor2(ppart,ompplib2.ppbuff,ncl,ihole,tsort)
# copy ordered electrons to linear array: updates part
         mpush2.mpcopyout2(part,ppart,kpic,nt,irc)
# part overflow
         if (irc[0] > 0):
            npmax = irc[0]
            maxnp = max(npmax,npimax)
            part = numpy.empty((idimp,maxnp),float_type,'F')
            irc[0] = 0
            mpush2.mpcopyout2(part,ppart,kpic,nt,irc)
         nppmx0 = int((1.0 + in2.xtras)*irc2[1])
         ppart = numpy.empty((idimp,nppmx0,mxyp1),float_type,'F')
# copies unordered electrons to ordered array: updates ppart
         mpush2.mpcopyin2(part,ppart,kpic,irc)
         ompplib2.ompmove2(ppart,kpic,ncl,ihole,noff,nyp,in2.xtras,
                           tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                           npbmx,nbmaxp,mx1,in2.popt,plist,irc2)

# sanity check for electrons
   if (in2.monitor > 0):
      mpush2.mpcheck2(ppart,kpic,noff,nyp,nx,in2.mx,in2.my,mx1,irc)

# electron current density diagnostic
   if (in2.ntje > 0):
      it = int(ntime/in2.ntje)
      if (ntime==in2.ntje*it):
         vfield[:,:,:] = cue
# transform electron current density to fourier space: updates cut
# moves data to uniform partition
         isign = -1
         ompplib2.wmpfft2rn(vfield,cut,noff,nyp,isign,mixup,sct,tfft,
                            tfmov,in2.indx,in2.indy,kstrt,nvp,kyp,ny,
                            nterf,ierr)
# calculate smoothed electron current in fourier space: updates vfieldc
         mfield2.mpsmooth23(cut,vfieldc,ffc,tfield,nx,ny,kstrt)
# store selected fourier modes: updates curet
         mfield2.mprdvmodes2(vfieldc,curet,tfield,nx,ny,in2.modesxje,
                             in2.modesyje,kstrt)
# write fourier space diagnostic output: updates njerec
#        mppmod2.mpvcwrite2(curit,tdiag,in2.modesxje,modesy2je,kxp,iuje,
#                           in2.njerec)
# transform smoothed electron current to real space: updates vfield
         isign = 1
         mfft2.mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,in2.indx,
                        in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates njerec
         mppmod2.mpvwrite2(vfield,tdiag,nx,ny,kyp,iuje,in2.njerec)
# move vector data to non-uniform partition, updates: vfield, ierr
         isign = 1
         mppmod2.mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,
                           nterf,ierr)
# display smoothed electron current
         ompplib2.wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
#        pgraf2.pdvector2(vfield,nyp,nvp,' ELECTRON CURRENT',ntime,999,
#                         0,ndstyle,1,nx,ny,ierr)
         if (ierr[0]==1):
            mpplib2.ppexit(); exit(1)
         ierr[0] = 0

# deposit ion current with OpenMP:
   if (in2.movion==1):
# updates pparti and cui, and possibly ncl, ihole, irc
      mpush2.mpset_pvzero2(cui,tdjpost,in2.mx,in2.my,mx1,myp1)
      mcurd2.wmpdjpost2(pparti,cui,kipic,ncl,ihole,noff,nyp,in2.qmi,dth,
                        in2.ci,tdjpost,nx,ny,in2.mx,in2.my,mx1,ipbc,
                        in2.popt,in2.relativity,plist,irc)
# add guard cells with OpenMP: updates cui
      ompplib2.wmpnacguard2(cui,nyp,tguard,nx,kstrt,nvp)

# reorder ions by tile with OpenMP and MPI
# updates: pparti, kipic, and irc and possibly ncl and ihole
      if (irc[0]==0):
         ompplib2.ompmove2(pparti,kipic,ncl,ihole,noff,nyp,in2.xtras,
                           tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                           npbmx,nbmaxp,mx1,in2.popt,plist,irc2)
      else:
         irc2[0] = 1; irc2[1] = irc[0]; irc[0] = 0

      while (irc2[0] != 0):
# ihole overflow
         if (irc2[0]==1):
            ntmaxp = int((1.0 + in2.xtras)*irc2[1])
            ihole = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
            irc2[:] = 0
            ompplib2.ompmove2(pparti,kipic,ncl,ihole,noff,nyp,in2.xtras,
                              tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                              npbmx,nbmaxp,mx1,in2.popt,False,irc2)
# pparti overflow
         elif (irc2[0]==4):
# restores ion coordinates from ppbuff: updates pparti, ncl
            msort2.mprstor2(pparti,ompplib2.ppbuff,ncl,ihole,tsort)
# copy ordered ions to linear array: updates part
            mpush2.mpcopyout2(part,pparti,kipic,nt,irc)
# part overflow
            if (irc[0] > 0):
               npimax = irc[0]
               maxnp = max(npmax,npimax)
               irc[0] = 0
               mpush2.mpcopyout2(part,pparti,kipic,nt,irc)
            nppmx1 = int((1.0 + in2.xtras)*irc2[1])
            pparti = numpy.empty((idimp,nppmx1,mxyp1),float_type,'F')
# copies unordered ions to ordered array: updates pparti
            mpush2.mpcopyin2(part,pparti,kipic,irc)
            ompplib2.ompmove2(pparti,kipic,ncl,ihole,noff,nyp,in2.xtras,
                              tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                              npbmx,nbmaxp,mx1,in2.popt,plist,irc2)

# sanity check for ions
      if (in2.monitor > 0):
         mpush2.mpcheck2(pparti,kipic,noff,nyp,nx,in2.mx,in2.my,mx1,irc)

# ion current density diagnostic
   if (in2.movion==1):
      if (in2.ntji > 0):
         it = int(ntime/in2.ntji)
         if (ntime==in2.ntji*it):
            vfield[:,:,:] = cui
# transform ion current density to fourier space: updates cut
# moves data to uniform partition
            isign = -1
            ompplib2.wmpfft2rn(vfield,cut,noff,nyp,isign,mixup,sct,tfft,
                               tfmov,in2.indx,in2.indy,kstrt,nvp,kyp,ny,
                               nterf,ierr)
# calculate smoothed ion current in fourier space: updates vfieldc
            mfield2.mpsmooth23(cut,vfieldc,ffc,tfield,nx,ny,kstrt)
# store selected fourier modes: updates curit
            mfield2.mprdvmodes2(vfieldc,curit,tfield,nx,ny,in2.modesxji,
                                in2.modesyji,kstrt)
# write fourier space diagnostic output: updates njirec
#              call mpvcwrite2(curit,tdiag,modesxji,modesy2ji,kxp,iuji, &
#    &njirec)
# transform smoothed ion current to real space: updates vfield
            isign = 1
            mfft2.mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,in2.indx,
                           in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates njirec
            mppmod2.mpvwrite2(vfield,tdiag,nx,ny,kyp,iuji,in2.njirec)
# move vector data to non-uniform partition, updates: vfield, ierr
            isign = 1
            mppmod2.mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,
                              nvp,nterf,ierr)
# display smoothed ion current
            ompplib2.wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
#           pgraf2.pdvector2(vfield,nyp,nvp,'ION CURRENT',ntime,999,0,
#                            in2.ndstyle,1,nx,ny,ierr)
            if (ierr[0]==1):
               mpplib2.ppexit(); exit(1)
            ierr[0] = 0

# deposit electron charge with OpenMP: updates qe
   mpush2.mpset_pszero2(qe,tdpost,in2.mx,in2.my,mx1,myp1)
   mpush2.mppost2(ppart,qe,kpic,noff,in2.qme,tdpost,in2.mx,in2.my,mx1,
                  in2.popt)
# add guard cells with OpenMP: updates qe
   ompplib2.wmpaguard2(qe,nyp,tguard,nx,kstrt,nvp)

# electron density diagnostic
   if (in2.ntde > 0):
      it = int(ntime/in2.ntde)
      if (ntime==in2.ntde*it):
         sfield[:,:] = -numpy.copy(qe)
# transform electron density to fourier space: updates qt
# moves data to uniform partition
         isign = -1
         ompplib2.wmpfft2r(sfield,qt,noff,nyp,isign,mixup,sct,tfft,
                           tfmov,in2.indx,in2.indy,kstrt,nvp,kyp,ny,
                           nterf,ierr)
# calculate smoothed electron density in fourier space: updates sfieldc
         mfield2.mpsmooth2(qt,sfieldc,ffc,tfield,nx,ny,kstrt)
# store selected fourier modes: updates denet
         mfield2.mprdmodes2(sfieldc,denet,tfield,nx,ny,in2.modesxde,
                            in2.modesyde,kstrt)
# write fourier space diagnostic output: updates nderec
#        mppmod2.mpcwrite2(denet,tdiag,in2.modesxde,modesy2de,kxp,iude,
#                          in2.nderec)
# transform smoothed electron density to real space: updates sfield
         isign = 1
         mfft2.mpfft2r(sfield,sfieldc,isign,mixup,sct,tfft,in2.indx,
                       in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates nderec
         mppmod2.mpwrite2(sfield,tdiag,nx,ny,kyp,iude,in2.nderec)
# move data to non-uniform partition, updates: sfield, ierr
         isign = 1
         mppmod2.mpfmove2(sfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,
                          nterf,ierr)
# display smoothed electron density
         ompplib2.wmpcguard2(sfield,nyp,tguard,nx,kstrt,nvp)
#        pgraf2.pdscaler2(sfield,nyp,nvp,'EDENSITY',ntime,999,1,
#                         in2.ndstyle,nx,ny,ierr)
         if (ierr[0]==1):
            mpplib2.ppexit(); exit(1)
         ierr[0] = 0

# deposit ion charge with OpenMP: updates qi
   if (in2.movion==1):
      mpush2.mpset_pszero2(qi,tdpost,in2.mx,in2.my,mx1,myp1)
      mpush2.mppost2(pparti,qi,kipic,noff,in2.qmi,tdpost,in2.mx,in2.my,
                     mx1,in2.popt)
# add guard cells with OpenMP: updates qi
      ompplib2.wmpaguard2(qi,nyp,tguard,nx,kstrt,nvp)

# ion density diagnostic
   if (in2.movion==1):
      if (in2.ntdi > 0):
         it = int(ntime/in2.ntdi)
         if (ntime==in2.ntdi*it):
            sfield[:,:] = numpy.copy(qi)
# transform ion density to fourier space: updates qt
# moves data to uniform partition
            isign = -1
            ompplib2.wmpfft2r(sfield,qt,noff,nyp,isign,mixup,sct,tfft,
                              tfmov,in2.indx,in2.indy,kstrt,nvp,kyp,ny,
                              nterf,ierr)
# calculate smoothed ion density in fourier space: updates sfieldc
            mfield2.mpsmooth2(qt,sfieldc,ffc,tfield,nx,ny,kstrt)
# store selected fourier modes: updates denit
            mfield2.mprdmodes2(sfieldc,denit,tfield,nx,ny,in2.modesxdi,
                               in2.modesydi,kstrt)
# write fourier space diagnostic output: updates ndirec
#           mppmod2.mpcwrite2(denit,tdiag,in2.modesxdi,modesy2di,kxp,
#                             iudi,in2.ndirec)
# transform smoothed ion density to real space: updates sfield
            isign = 1
            mfft2.mpfft2r(sfield,sfieldc,isign,mixup,sct,tfft,in2.indx,
                          in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates ndirec
            mppmod2.mpwrite2(sfield,tdiag,nx,ny,kyp,iudi,in2.ndirec)
# move data to non-uniform partition, updates: sfield, ierr
            isign = 1
            mppmod2.mpfmove2(sfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,
                             nvp,nterf,ierr)
# display smoothed ion density
            ompplib2.wmpcguard2(sfield,nyp,tguard,nx,kstrt,nvp)
#           pgraf2.pdscaler2(sfield,nyp,nvp,'ION DENSITY',ntime,999,1,
#                            in2.ndstyle,nx,ny,ierr)
            if (ierr[0]==1):
               mpplib2.ppexit(); exit(1)
            ierr[0] = 0

# add electron and ion densities: updates qe
   mfield2.mpaddqei2(qe,qi,nyp,tfield,nx)

# add electron and ion current densities: updates cue
   if (in2.movion==1):
      mfield2.mpaddcuei2(cue,cui,nyp,tfield,nx)

# transform charge to fourier space with OpenMP:
# moves data to uniform partition
# updates qt, nterf, and ierr, modifies qe
   isign = -1
   ompplib2.wmpfft2r(qe,qt,noff,nyp,isign,mixup,sct,tfft,tfmov,in2.indx,
                     in2.indy,kstrt,nvp,kyp,ny,nterf,ierr)

# transform current to fourier space with OpenMP:
# moves data to uniform partition
# updates cut, nterf, and ierr, modifies cue
   isign = -1
   ompplib2.wmpfft2rn(cue,cut,noff,nyp,isign,mixup,sct,tfft,tfmov,
                      in2.indx,in2.indy,kstrt,nvp,kyp,ny,nterf,ierr)

# take transverse part of current with OpenMP: updates cut
   mfield2.mpcuperp2(cut,tfield,nx,ny,kstrt)

# radiative vector potential diagnostic
   if (in2.ntar > 0):
      it = int(ntime/in2.ntar)
      if (ntime==in2.ntar*it):
# average current: updates vfieldc = 0.5*(cut + oldcut)
         mfield2.mcuave2(vfieldc,cut,oldcut,tfield,ny)
# calculate radiative vector potential in fourier space: updates vfieldc
# vfieldc should contain averaged current on entry
         mfield2.mpavrpot2(vfieldc,bxyz,ffc,affp,in2.ci,tfield,nx,ny,
                           kstrt)
# store selected fourier modes: updates vpotr
         mfield2.mprdvmodes2(vfieldc,vpotr,tfield,nx,ny,in2.modesxar,
                             in2.modesyar,kstrt)
# write fourier space diagnostic output: updates narrec
#        mppmod2.mpvcwrite2(vpotr,tdiag,in2.modesxar,modesy2ar,kxp,iuar,
#                           in2.narrec)
# transform radiative vector potential to real space: updates vfield
         isign = 1
         mfft2.mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,in2.indx,
                        in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates narrec
         mppmod2.mpvwrite2(vfield,tdiag,nx,ny,kyp,iuar,in2.narrec)

# move vector data to non-uniform partition, updates: vfield, ierr
         isign = 1
         mppmod2.mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,
                           nterf,ierr)
# display radiative vector potential
         ompplib2.wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
#        pgraf2.pdvector2(vfield,nyp,nvp,'RADIATIVE VPOTENTIAL',ntime,
#                         999,0,in2.ndstyle,1,nx,ny,ierr)
         if (ierr[0]==1):
            mpplib2.ppexit(); exit(1)
         ierr[0] = 0

# calculate electromagnetic fields in fourier space with OpenMP:
# updates exyz, bxyz, wf, wb
   if ((ntime+ntime0[0])==0):
# initialize electromagnetic fields from darwin fields
# calculate initial darwin magnetic field
      mfield2.mpibpois2(cut,bxyz,ffc,in2.ci,wb,tfield,nx,ny,kstrt)
      wf[0] = 0.0
# calculate initial darwin electric field
      amu = numpy.zeros((4,nxe,nypmx[0]),float_type,'F')
      amut = numpy.empty((4,nye,kxp),complex_type,'F')
      dcut = numpy.empty((in2.ndim,nye,kxp),complex_type,'F')
      mpush2.mpset_pvzero2(amu,tdjpost,in2.mx,in2.my,mx1,myp1)
      mcurd2.wmpgmjpost2(ppart,amu,kpic,noff,in2.qme,in2.ci,tdjpost,
                         in2.mx,in2.my,mx1,in2.popt,in2.relativity)
      ompplib2.wmpnacguard2(amu,nyp,tguard,nx,kstrt,nvp)
      isign = -1
      ompplib2.wmpfft2rn(amu,amut,noff,nyp,isign,mixup,sct,tfft,tfmov,
                         in2.indx,in2.indy,kstrt,nvp,kyp,ny,nterf,ierr)
#     numpy.delete(amu,amu)
      mfield2.mpdcuperp2(dcut,amut,tfield,nx,ny,kstrt)
#     numpy.delete(amut,amut)
      mfield2.mpetfield2(dcut,exyz,ffc,affp,in2.ci,wf,tfield,nx,ny,
                         kstrt)
#     numpy.delete(dcut,dcut)
      dth = 0.5*in2.dt
# update electromagnetic fields
   else:
      mfield2.mpmaxwel2(exyz,bxyz,cut,ffc,affp,in2.ci,in2.dt,wf,wb,
                        tfield,nx,ny,kstrt)

# calculate longitudinal force/charge in fourier space with OpenMP:
# updates fxyt, we
   mfield2.mppois2(qt,fxyt,ffc,we,tfield,nx,ny,kstrt)

# add longitudinal and transverse electric fields with OpenMP:
# updates fxyt
   isign = 1
   mfield2.mpemfield2(fxyt,exyz,ffc,isign,tfield,nx,ny,kstrt)
# copy magnetic field with OpenMP: updates bxyt
   isign = -1
   mfield2.mpemfield2(bxyt,bxyz,ffc,isign,tfield,nx,ny,kstrt)

# transform force to real space with OpenMP:
# moves data to non-uniform partition
# updates fxyze, nterf, and ierr, modifies fxyt
   isign = 1
   ompplib2.wmpfft2rn(fxyze,fxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,
                      in2.indx,in2.indy,kstrt,nvp,kyp,ny,nterf,ierr)

# transform magnetic field to real space with OpenMP:
# moves data to non-uniform partition
# updates bxyze, nterf, and ierr, modifies bxyt
   isign = 1
   ompplib2.wmpfft2rn(bxyze,bxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,
                      in2.indx,in2.indy,kstrt,nvp,kyp,ny,nterf,ierr)

# add constant to magnetic field with OpenMP: updates bxyze
   if (omt > 0.0):
      mfield2.mpbaddext2(bxyze,nyp,tfield,in2.omx,in2.omy,in2.omz,nx)

# copy guard cells with OpenMP: updates fxyze, bxyze
   ompplib2.wmpncguard2(fxyze,nyp,tguard,nx,kstrt,nvp)
   ompplib2.wmpncguard2(bxyze,nyp,tguard,nx,kstrt,nvp)

# potential diagnostic
   if (in2.ntp > 0):
      it = int(ntime/in2.ntp)
      if (ntime==in2.ntp*it):
# calculate potential in fourier space: updates sfieldc
         mfield2.mppot2(qt,sfieldc,ffc,ws,tfield,nx,ny,kstrt)
# store selected fourier modes: updates pott
         mfield2.mprdmodes2(sfieldc,pott,tfield,nx,ny,in2.modesxp,
                            in2.modesyp,kstrt)
# write fourier space diagnostic output: updates nprec
#        mppmod2.mpcwrite2(pott,tdiag,in2.modesxp,imodesy2p,kxp,iup,
#                          in2.nprec)
# transform potential to real space: updates sfield
         isign = 1
         mfft2.mpfft2r(sfield,sfieldc,isign,mixup,sct,tfft,in2.indx,
                       in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates nprec
         mppmod2.mpwrite2(sfield,tdiag,nx,ny,kyp,iup,in2.nprec)
# move data to non-uniform partition, updates: sfield, ierr
         isign = 1
         mppmod2.mpfmove2(sfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,
                          nterf,ierr)
# display potential
         ompplib2.wmpcguard2(sfield,nyp,tguard,nx,kstrt,nvp)
#        pgraf2.pdscaler2(sfield,nyp,nvp,'POTENTIAL',ntime,999,0,
#                         in2.ndstyle,nx,ny,ierr)
         if (ierr[0]==1):
            mpplib2.ppexit(); exit(1)
         ierr[0] = 0

# longitudinal efield diagnostic
   if (in2.ntel > 0):
      it = int(ntime/in2.ntel)
      if (ntime==in2.ntel*it):
# calculate longitudinal efield in fourier space: updates vfieldc
         mfield2.mpelfield2(qt,vfieldc,ffc,ws,tfield,nx,ny,kstrt)
# store selected fourier modes: updates elt
         mfield2.mprdvmodes2(vfieldc,elt,tfield,nx,ny,in2.modesxel,
                             in2.modesyel,kstrt)
# write fourier space diagnostic output: updates nelrec
#        mppmod2.mpvcwrite2(elt,tdiag,in2.modesxel,modesy2el,kxp,iuel,
#                           in2.nelrec)
# transform longitudinal efield to real space: updates vfield
         isign = 1
         mfft2.mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,in2.indx,
                        in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates nelrec
         mppmod2.mpvwrite2(vfield,tdiag,nx,ny,kyp,iuel,in2.nelrec)
# move data to non-uniform partition, updates: vfield, ierr
         isign = 1
         mppmod2.mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,
                           nterf,ierr)
# display longitudinal efield 
         ompplib2.wmpncguard2(vfield,nyp,tguard,nx,kstrt,nvp)
#        pgraf2.pdvector2(vfield,nyp,nvp,'ELFIELD',ntime,999,0,
#                         in2.ndstyle,1,nx,ny,ierr)
         if (ierr[0]==1):
            mpplib2.ppexit(); exit(1)
         ierr[0] = 0

# vector potential diagnostic
   if (in2.nta > 0):
      it = int(ntime/in2.nta)
      if (ntime==in2.nta*it):
# calculate vector potential in fourier space: updates vfieldc
         mfield2.mpavpot2(bxyz,vfieldc,tfield,nx,ny,kstrt)
# store selected fourier modes: updates vpott
         mfield2.mprdvmodes2(vfieldc,vpott,tfield,nx,ny,in2.modesxa,
                             in2.modesya,kstrt)
# write fourier space diagnostic output: updates narec
#        mppmod2.mpvcwrite2(vpott,tdiag,in2.modesxa,modesy2a,kxp,iua,
#                           in2.narec)
# transform vector potential to real space: updates vfield
         isign = 1
         mfft2.mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,in2.indx,
                        in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates narec
         mppmod2.mpvwrite2(vfield,tdiag,nx,ny,kyp,iua,in2.narec)
# move vector data to non-uniform partition, updates: vfield, ierr
         isign = 1
         mppmod2.mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,
                            nterf,ierr)
# display vector potential
         ompplib2.wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
#        pgraf2.pdvector2(vfield,nyp,nvp,'VECTOR POTENTIAL',ntime,999,0.
#                         in2.ndstyle,1,nx,ny,ierr)
         if (ierr[0]==1):
            mpplib2.ppexit(); exit(1)
         ierr[0] = 0

# transverse efield diagnostic
   if (in2.ntet > 0):
      it = int(ntime/in2.ntet)
      if (ntime==in2.ntet*it):
         vfieldc[:,:,:] = exyz
# store selected fourier modes: updates ett
         mfield2.mprdvmodes2(vfieldc,ett,tfield,nx,ny,in2.modesxet,
                             in2.modesyet,kstrt)
# write fourier space diagnostic output: updates netrec
#        mppmod2.mpvcwrite2(ett,tdiag,in2.modesxet,modesy2et,kxp,iuet,
#                           in2.netrec)
# transform transverse efield to real space: updates vfield
         isign = 1
         mfft2.mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,in2.indx,
                        in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates netrec
         mppmod2.mpvwrite2(vfield,tdiag,nx,ny,kyp,iuet,in2.netrec)
# move data to non-uniform partition, updates: vfield, ierr
         isign = 1
         mppmod2.mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,
                           nterf,ierr)
# display transverse efield 
         ompplib2.wmpncguard2(vfield,nyp,tguard,nx,kstrt,nvp)
#        pgraf2.pdvector2(vfield,nyp,nvp,'TRANSVERSE EFIELD',ntime,999,
#                         0,in2.ndstyle,1,nx,ny,ierr)
         if (ierr[0]==1):
            mpplib2.ppexit(); exit(1)
         ierr[0] = 0

# magnetic field diagnostic
   if (in2.ntb > 0):
      it = int(ntime/in2.ntb)
      if (ntime==in2.ntb*it):
         vfieldc[:,:,:] = bxyz
# store selected fourier modes: updates bt
         mfield2.mprdvmodes2(vfieldc,bt,tfield,nx,ny,in2.modesxb,
                             in2.modesyb,kstrt)
# write fourier space diagnostic output: updates nbrec
#        mppmod2.mpvcwrite2(bt,tdiag,in2.modesxb,modesy2b,kxp,iub,
#                           in2.nbrec)
# transform magnetic field to real space: updates vfield
         isign = 1
         mfft2.mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,in2.indx,
                        in2.indy,kstrt,nvp,kyp)
# write real space diagnostic output: updates nbrec
         mppmod2.mpvwrite2(vfield,tdiag,nx,ny,kyp,iub,in2.nbrec)
# move vector data to non-uniform partition, updates: vfield, ierr
         isign = 1
         mppmod2.mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,
                           nterf,ierr)
# display magnetic field
         ompplib2.wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
#        pgraf2.pdvector2(vfield,nyp,nvp,'MAGNETIC FIELD',ntime,999,0,
#                         in2.ndstyle,1,nx,ny,ierr)
         if (ierr[0]==1):
            mpplib2.ppexit(); exit(1)
         ierr[0] = 0

# velocity-space diagnostic
   if (in2.ntv > 0):
      it = int(ntime/in2.ntv)
      if (ntime==in2.ntv*it):
         if ((in2.ndv==1) or (in2.ndv==3)):
# calculate electron cartesian distribution function and moments
            if ((in2.nvft==1) or (in2.nvft==3)):
               mdiag2.mpvpdist2(ppart,kpic,fv,sfv,fvm,tdiag,nvp,in2.nmv)
# store time history electron vdrift, vth, and entropy
               fvtm[itv,:,:] = fvm
               if (in2.movion==0):
                   itv += 1
# display electron velocity distributions
#              pgraf1.pdisplayfv1(fv,fvm,' ELECTRON',kstrt,ntime,
#                                 in2.nmv,2,irc)
               if (irc[0]==1):
                  exit(1); irc[0] = 0
# calculate electron cylindrical distribution function and moments
            if ((in2.nvft==4) or (in2.nvft==5)):
               mdiag2.mpvbpdist2(ppart,kpic,fv,sfv,fvm,in2.omx,in2.omy,
                                 in2.omz,tdiag,in2.nmv)
# display electron velocity distributions in cylindrical co-ordinates
#              pgraf1.pdisplayfvb1(fv,fvm,' ELECTRON',kstrt,ntime,
#                                  in2.nmv,2,irc)
               if (irc[0]==1):
                  exit(1); irc[0] = 0
# electron energy distribution
            if ((in2.nvft==2) or (in2.nvft==3) or (in2.nvft==5)):
               mdiag2.mperpdist2(ppart,kpic,fe,sfv,eci,wk,tdiag,
                                 in2.ndim,in2.nmv)
# display electron energy distribution
#              pgraf1.pdisplayfe1(fe,wk,' ELECTRON',kstrt,ntime,in2.mv,
#                                  irc)
               if (irc[0]==1):
                  exit(1); irc[0] = 0
# write electron velocity-space diagnostic output: updates nverec
            if (kstrt==1):
               mdiag2.dafwritefv2(fvm,fv,fe,wk,tdiag,iuve,in2.nverec)
# ion distribution functions
         if (in2.movion==1):
             if ((in2.ndv==2) or (in2.ndv==3)):
# calculate ion cartesian distribution function and moments
                if ((in2.nvft==1) or (in2.nvft==3)):
                   mdiag2.mpvpdist2(pparti,kipic,fvi,sfv,fvmi,tdiag,nvp,
                                    in2.nmv)
# store time history of ion vdrift, vth, and entropy
                   fvtmi[itv,:,:] = fvmi
                   itv += 1
# display ion velocity distributions
#                  pgraf1.pdisplayfv1(fvi,fvmi,' ION',kstrt,ntime,
#                                     in2.nmv,2,irc)
                   if (irc[0]==1):
                      exit(1); irc[0] = 0
# calculate ion cylindrical distribution function and moments
                if ((in2.nvft==4) or (in2.nvft==5)):
                   mdiag2.mpvbpdist2(pparti,kipic,fvi,sfv,fvmi,in2.omx,
                                     in2.omy,in2.omz,tdiag,in2.nmv)
# display ion velocity distributions in cylindrical co-ordinates
#                  pgraf1.pdisplayfvb1(fvi,fvmi,' ION',kstrt,ntime,
#                                     in2.nmv,2,irc)
                   if (irc[0]==1):
                      exit(1); irc[0] = 0
# ion energy distribution
                if ((in2.nvft==2) or (in2.nvft==3) or (in2.nvft==5)):
                   mdiag2.mperpdist2(pparti,kipic,fei,sfv,eci,wk,tdiag,
                                     in2.ndim,in2.nmv)
                   wk[0] = in2.rmass*wk[0]
# display ion energy distribution
                   ts = fei[nmv21,0]
                   fei[nmv21,0] = in2.rmass*fei[nmv21,0]
#                  pgraf1.pdisplayfe1(fei,wk,' ION',kstrt,ntime,in2.nmv,
#                                     irc)
                   fei[nmv21,0] = ts
                   if (irc[0]==1):
                      exit(1); irc[0] = 0
# write ion velocity-space diagnostic output: updates nvirec
                   if (kstrt==1):
                       mdiag2.dafwritefv2(fvmi,fvi,fei,wk,tdiag,iuvi,
                                         in2.nvirec)

# trajectory diagnostic
   if (in2.ntt > 0):
      it = int(ntime/in2.ntt)
      if (ntime==in2.ntt*it):
         ierr[0] = 0
# copies tagged electrons in ppart to array partt: updates partt, numtp
         if (in2.ndt==1):
            mdiag2.mptraj2(ppart,kpic,partt,tdiag,numtp,ierr)
# copies tagged ions in ppart to array partt: updates partt, numtp
         elif (in2.ndt==2):
            if (in2.movion==1):
               mdiag2.mptraj2(pparti,kipic,partt,tdiag,numtp,ierr)
         if (ierr[0] != 0):
# partt overflow
            if (ierr[0] <  0):
               in2.nprobt = numtp[0]
               partt = numpy.empty((idimp,in2.nprobt),float_type,'F')
               ierr[0] = 0
# copies tagged electrons in ppart to array partt: updates partt, numtp
               if (in2.ndt==1):
                   mdiag2.mptraj2(ppart,kpic,partt,tdiag,numtp,ierr)
# copies tagged ions in ppart to array partt: updates partt, numtp
               elif (in2.ndt==2):
                  if (in2.movion==1):
                     mdiag2.mptraj2(pparti,kipic,partt,tdiag,numtp,ierr)
            if (ierr[0] != 0):
               mpplib2.ppexit();  exit(1)
# electron or ion trajectories
         if ((in2.ndt==1) or (in2.ndt==2)):
# reorder tagged particles
            if ((in2.nst==1) or (in2.nst==2)):
# determines list of tagged pareticles leaving this node: updates iholep
               minit2.mpfholes2(partt,tedges,numtp,iholep,in2.ndim,2)
# iholep overflow
               if (iholep[0] < 0):
                  ntmax = -iholep[0]
                  ntmax = int(1.5*ntmax)
                  if (kstrt==1):
                     print "info:reallocating iholep:ntmax=", ntmax
                  iholep = numpy.empty((ntmax+1),int_type,'F')
                  minit2.mpfholes2(partt,tedges,numtp,iholep,in2.ndim,2)
                  if (iholep[0] < 0):
                     if (kstrt==1):
                         print "iholep overflow: ntmax=", ntmax
                     mpplib2.ppexit(); exit(1)
# copies tagged particles: updates part
               mdiag2.mpcpytraj2(partt,part,tdiag,numtp)
# moves tagged electrons into original spatial region:
# updates part, numtp
               mppmod2.ipmove2(part,tedges,numtp,iholep,ny,tmov,kstrt,
                               nvp,in2.ndim,2,ierr)
               if (ierr[0] != 0):
                  mpplib2.ppexit(); exit(1)
# reorders tagged particles: updates partt
               mdiag2.mpordtraj2(part,partt,tedges,tdiag,numtp,ierr)
# collects distributed test particle data onto node 0
               mppmod2.mppartt2(partt,tdiag,numtp,ierr)
               if (ierr[0] != 0):
                  mpplib2.ppexit(); exit(1)
# write trajectory diagnostic output: updates ntrec
               if (kstrt==1):
                  mdiag2.dafwritetr2(partt,tdiag,iut,in2.ntrec)
                  partd[itt,:,:] = partt
                  itt += 1
            elif (in2.nst==3):
# calculate test particle distribution function and moments
               mdiag2.mpvdist2(partt,fvtp,fvmtp,tdiag,numtp,nvp,in2.nmv)
# write test particle diagnostic output: updates ntrec
               if (kstrt==1):
                  ws[0] = 0.0
                  mdiag2.dafwritefv2(fvmtp,fvtp,fetp,ws,tdiag,iut,
                                     in2.ntrec)
# display test particle velocity distributions
#              pgraf1.pdisplayfv1(fvtp,fvmtp,' ELECTRON',kstrt,ntime,
#                                 in2.nmv,1,irc)
               if (irc[0]==1):
                  exit(1); irc[0] = 0

# phase space diagnostic
   if (in2.nts > 0):
      it = int(ntime/in2.nts)
      if (ntime==in2.nts*it):
# electron phase space diagnostic
         if ((in2.nds==1) or (in2.nds==3)):
# calculates velocity distribution in different regions of space:
# updates fvs
            mdiag2.mpvspdist2(ppart,kpic,fvs,tdiag,noff,in2.nmv,in2.mvx,
                              in2.mvy)
# adjusts 3d velocity distribution in different regions of space:
#updates fvs
            mppmod2.mpadjfvs2(fvs,tdiag,noff,nyp,in2.nmv,in2.mvy)
# write phase space diagnostic output: updates nserec
            if (in2.nserec > 0):
               mppmod2.mpwrfvsdata2(fvs,tdiag,nyb,nybmx,iuse)
               in2.nserec += 1
# ion phase space
         if (in2.movion==1):
            if ((in2.nds==2) or (in2.nds==3)):
# calculates velocity distribution in different regions of space:
# updates fvsi
               mdiag2.mpvspdist2(pparti,kipic,fvsi,tdiag,noff,in2.nmv,
                                 in2.mvx,in2.mvy)
# adjusts 3d velocity distribution in different regions of space:
# updates fvsi
               mppmod2.mpadjfvs2(fvsi,tdiag,noff,nyp,in2.nmv,in2.mvy)
# write phase space diagnostic output: updates nsirec
               if (in2.nsirec > 0):
                  mppmod2.mpwrfvsdata2(fvsi,tdiag,nyb,nybmx,iusi)
                  in2.nsirec += 1

# push electrons with OpenMP:
# updates ppart and wke, and possibly ncl, ihole, irc
   wke[0] = 0.0
   mbpush2.wmpbpush2(ppart,fxyze,bxyze,kpic,ncl,ihole,noff,nyp,qbme,
                     in2.dt,dth,in2.ci,wke,tpush,nx,ny,in2.mx,in2.my,
                     mx1,ipbc,in2.popt,in2.relativity,plist,irc)

# reorder electrons by tile with OpenMP and MPI
# updates: ppart, kpic, and irc and possibly ncl and ihole
   if (irc[0]==0):
      ompplib2.ompmove2(ppart,kpic,ncl,ihole,noff,nyp,in2.xtras,tsort,
                        tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,npbmx,nbmaxp,
                        mx1,in2.popt,plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc[0]; irc[0] = 0

   while (irc2[0] != 0):
# ihole overflow
      if (irc2[0]==1):
         ntmaxp = int((1.0 + in2.xtras)*irc2[1])
         ihole = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
         irc2[:] = 0
         ompplib2.ompmove2(ppart,kpic,ncl,ihole,noff,nyp,in2.xtras,
                           tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                           npbmx,nbmaxp,mx1,in2.popt,False,irc2)
# ppart overflow
      elif (irc2[0]==4):
# restores electron coordinates from ppbuff: updates ppart, ncl
         msort2.mprstor2(ppart,ompplib2.ppbuff,ncl,ihole,tsort)
# copy ordered electrons to linear array: updates part
         mpush2.mpcopyout2(part,ppart,kpic,nt,irc)
# part overflow
         if (irc[0] > 0):
            npmax = irc[0]
            maxnp = max(npmax,npimax)
            part = numpy.empty((idimp,maxnp),float_type,'F')
            irc[0] = 0
            mpush2.mpcopyout2(part,ppart,kpic,nt,irc)
         nppmx0 = int((1.0 + in2.xtras)*irc2[1])
         ppart = numpy.empty((idimp,nppmx0,mxyp1),float_type,'F')
# copies unordered electrons to ordered array: updates ppart
         mpush2.mpcopyin2(part,ppart,kpic,irc)
         ompplib2.ompmove2(ppart,kpic,ncl,ihole,noff,nyp,in2.xtras,
                           tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                           npbmx,nbmaxp,mx1,in2.popt,plist,irc2)

# sanity check for electrons
   if (in2.monitor > 0):
      mpush2.mpcheck2(ppart,kpic,noff,nyp,nx,in2.mx,in2.my,mx1,irc)

# push ions with OpenMP:
   if (in2.movion==1):
# updates pparti and wki, and possibly ncl, ihole, irc
      wki[0] = 0.0
      mbpush2.wmpbpush2(pparti,fxyze,bxyze,kipic,ncl,ihole,noff,nyp,
                        qbmi,in2.dt,dth,in2.ci,wki,tpush,nx,ny,in2.mx,
                        in2.my,mx1,ipbc,in2.popt,in2.relativity,plist,
                        irc)
      wki[0] = wki[0]*in2.rmass

# reorder ions by tile with OpenMP and MPI
# updates: pparti, kipic, and irc and possibly ncl and ihole
      if (irc[0]==0):
         ompplib2.ompmove2(pparti,kipic,ncl,ihole,noff,nyp,in2.xtras,
                           tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                           npbmx,nbmaxp,mx1,in2.popt,plist,irc2)
      else:
         irc2[0] = 1; irc2[1] = irc[0]; irc[0] = 0

      while (irc2[0] != 0):
# ihole overflow
         if (irc2[0]==1):
            ntmaxp = int((1.0 + in2.xtras)*irc2[1])
            ihole = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
            irc2[:] = 0
            ompplib2.ompmove2(pparti,kipic,ncl,ihole,noff,nyp,in2.xtras,
                              tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                              npbmx,nbmaxp,mx1,in2.popt,False,irc2)
# pparti overflow
         elif (irc2[0]==4):
# restores ion coordinates from ppbuff: updates pparti, ncl
            msort2.mprstor2(pparti,ompplib2.ppbuff,ncl,ihole,tsort)
# copy ordered ions to linear array: updates part
            mpush2.mpcopyout2(part,pparti,kipic,nt,irc)
# part overflow
            if (irc[0] > 0):
               npimax = irc[0]
               maxnp = max(npmax,npimax)
               irc[0] = 0
               mpush2.mpcopyout2(part,pparti,kipic,nt,irc)
            nppmx1 = int((1.0 + in2.xtras)*irc2[1])
            pparti = numpy.empty((idimp,nppmx1,mxyp1),float_type,'F')
# copies unordered ions to ordered array: updates pparti
            mpush2.mpcopyin2(part,pparti,kipic,irc)
            ompplib2.ompmove2(pparti,kipic,ncl,ihole,noff,nyp,in2.xtras,
                              tsort,tmov,kstrt,nvp,nx,ny,in2.mx,in2.my,
                              npbmx,nbmaxp,mx1,in2.popt,plist,irc2)

# sanity check for ions
      if (in2.monitor > 0):
         mpush2.mpcheck2(pparti,kipic,noff,nyp,nx,in2.mx,in2.my,mx1,irc)

# start running simulation backwards:
# need to advance maxwell field solver one step ahead
   if (in2.treverse==1):
      if (((ntime+1)==(nloop/2)) or ((ntime+1)==nloop)):
# deposit electron current: updates cue
         mpush2.mpset_pvzero2(cue,tdjpost,in2.mx,in2.my,mx1,myp1)
         mcurd2.wmpdjpost2(ppart,cue,kpic,ncl,ihole,noff,nyp,in2.qme,
                           zero,in2.ci,tdjpost,nx,ny,in2.mx,in2.my,mx1,
                           ipbc,in2.popt,in2.relativity,plist,irc)
         ompplib2.wmpnacguard2(cue,nyp,tguard,nx,kstrt,nvp)
# deposit ion current: updates cui
         if (in2.movion==1):
            mpush2.mpset_pvzero2(cui,tdjpost,in2.mx,in2.my,mx1,myp1)
            mcurd2.wmpdjpost2(pparti,cui,kipic,ncl,ihole,noff,nyp,
                              in2.qmi,zero,in2.ci,tdjpost,nx,ny,in2.mx,
                              in2.my,mx1,ipbc,in2.popt,in2.relativity,
                              plist,irc)
            ompplib2.wmpnacguard2(cui,nyp,tguard,nx,kstrt,nvp)
            mfield2.mpaddcuei2(cue,cui,nyp,tfield,nx)
         isign = -1
         ompplib2.wmpfft2rn(cue,cut,noff,nyp,isign,mixup,sct,tfft,tfmov,
                            in2.indx,in2.indy,kstrt,nvp,kyp,ny,nterf,
                            ierr)
         mfield2.mpcuperp2(cut,tfield,nx,ny,kstrt)
# updates exyz, bxyz, wf, wb
         mfield2.mpmaxwel2(exyz,bxyz,cut,ffc,affp,in2.ci,in2.dt,wf,wb,
                           tfield,nx,ny,kstrt)
# reverse time
         in2.dt = -in2.dt; dth = -dth

# energy diagnostic
   if (in2.ntw > 0):
      it = int(ntime/in2.ntw)
      if (ntime==in2.ntw*it):
         wef[0] = we[0] + wf[0] + wb[0]
         wtot[0] = wef[0]
         wtot[1] = wke[0]
         wtot[2] = wki[0]
         wtot[3] = wef[0] + wke[0]
         wtot[4] = we[0]
         wtot[5] = wf[0]
         wtot[6] = wb[0]
         mppmod2.mpdsum(wtot,tdiag)
         wke[0] = wtot[1]
         wki[0] = wtot[2]
         we[0] = wtot[4]
         wf[0] = wtot[5]
         wb[0] = wtot[6]
         wef[0] = we[0] + wf[0] + wb[0]
         ws[0] = wef[0] + wke[0] + wki[0]
         if (ntime==0):
            s[5] = ws[0]
         if (kstrt==1):
            print >> iuot, "Total Field, Kinetic and Total Energies:"
            if (in2.movion==0):
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
   if (in2.ntr > 0):
      n = ntime + 1
      it = int(n/in2.ntr)
      if (n==in2.ntr*it):
# write out basic restart file for electrostatic code
         f2.bwrite_restart2(part,ppart,pparti,qi,kpic,kipic,tdiag,kstrt,
                            iur,iscr,n,ntime0,irc)
# write out basic restart file for electromagnetic code
         fb2.bwrite_restart23(exyz,bxyz,tdiag,kstrt,iur)

ntime = ntime + 1

# loop time
dtimer(dtime,ltime,1)
tloop += float(dtime)

# * * * end main iteration loop * * *

# reset graphs
#     if ((ntw > 0).or.(ntt > 0).or.(ntv > 0)) then
#        if (nplot > 0) call reset_pgraphs(kstrt,irc)
#     endif

if (kstrt==1):
   print >> iuot
   print >> iuot, "ntime, relativity = ", ntime, ",", in2.relativity
   if (in2.treverse==1):
      print >> iuot, "treverse = ", in2.treverse
   print >> iuot, "MPI nodes nvp = ", nvp[0]

# trajectory diagnostic
#     if (ntt > 0) then
#        if ((ndt==1).or.(ndt==2)) then
#           if ((nst==1).or.(nst==2)) then
#              ts = t0 + dt*real(ntt)
# displays time history of trajectories on node 0
#              if (kstrt==1) then
#                 if (nplot > 0) call reset_nplot(1,irc)
#                 call pdisplaytr1(partd,ts,dt*real(ntt),kstrt,itt,3,999&
#    &,irc)
#                 if (irc==1) stop
#                 if (nplot > 0) call reset_nplot(nplot,irc)
#              endif
#           endif
#        endif
#     endif

# energy diagnostic
if (in2.ntw > 0):
#        ts = t0 + dt*real(ntw)
#        call pdisplayw1(wt,ts,dt*real(ntw),kstrt,itw,irc)
   if (kstrt==1):
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

# velocity-space diagnostic
#     if (ntv > 0) then
#        ts = t0 + dt*real(ntv)
#        call pdisplayfvt1(fvtm,' ELECT',ts,dt*real(ntv),kstrt,itv,irc)
# ions
#        if (movion==1) then
#           call pdisplayfvt1(fvtmi,' ION',ts,dt*real(ntv),kstrt,itv,   &
#    &irc)
#        endif
#     endif

if (kstrt==1):
   print >> iuot
   print >> iuot, "initialization time = ", tinit
   print >> iuot, "deposit time = ", tdpost[0]
   print >> iuot, "current deposit time = ", tdjpost[0]
   tdpost[0] = tdpost[0] + tdjpost[0]
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

# reset parameters for final diagnostic metafile
# electron density diagnostic
if (in2.ntde > 0):
   in2.nderec -= 1
# potential diagnostic
if (in2.ntp > 0):
   in2.nprec -= 1; in2.ceng = affp
# longitudinal efield diagnostic
if (in2.ntel > 0):
   in2.nelrec -= 1; in2.ceng = affp
# electron current diagnostic
if (in2.ntje > 0):
   in2.njerec -= 1
# radiative vector potential diagnostic
if (in2.ntar > 0):
   in2.narrec -= 1; in2.ceng = affp
# vector potential diagnostic
if (in2.nta > 0):
   in2.narec -= 1; in2.ceng = affp
# transverse efield diagnostic
if (in2.ntet > 0):
   in2.netrec -= 1; in2.ceng = affp
# magnetic field diagnostic
if (in2.ntb > 0):
   in2.nbrec -= 1; in2.ceng = affp
# fluid moments diagnostic
if (in2.ntfm > 0):
   if ((in2.ndfm==1) or (in2.ndfm==3)):
      in2.nferec -= 1
   if (in2.movion==1):
      if ((in2.ndfm==2) or (in2.ndfm==3)):
         in2.nfirec -= 1
   in2.ceng = affp
# velocity-space diagnostic
if (in2.ntv > 0):
   if ((in2.ndv==1) or (in2.ndv==3)):
      in2.nverec -= 1
   if (in2.movion==1):
      if ((in2.ndv==2) or (in2.ndv==3)):
         in2.nvirec -= 1
# trajectory diagnostic
if (in2.ntt > 0) :
   in2.ntrec -= 1
# phase space diagnostic
if (in2.nts > 0):
   if ((in2.nds==1) or (in2.nds==3)):
      in2.nserec -= 1
   if (in2.movion==1):
      if ((in2.nds==2) or (in2.nds==3)):
         in2.nsirec -= 1
# ion diagnostics
if (in2.movion==1):
# ion density diagnostic
   if (in2.ntdi > 0):
      in2.ndirec -= 1
# ion current diagnostic
   if (in2.ntji > 0):
      in2.njirec -= 1
# write final diagnostic metafile
in2.writnml2(iudm)
#close(unit=iudm)
# close restart files
f2.close_restart2(iur,iur0)
# close output file
print >> iuot, " * * * q.e.d. * * *"
iuot.close()
# close graphics device
#pgraf1.close_pgraphs

mpplib2.ppexit()
