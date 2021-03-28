#-----------------------------------------------------------------------
"""
High Level library for 1D Electrostatic OpenMP PIC code

functions defined:

init_fields1: allocate field data for standard code
del_fields1: delete field data for standard code

init_electrons1: initialize electrons
reorder_electrons1: recover from electron buffer overflow errors
push_electrons1: push electrons with OpenMP:
del_electrons1: delete electrons

init_ions1: initialize ions
reorder_ions1: recover from ion buffer overflow errors
push_ions1: push ions with OpenMP:
del_ions1: delete ions

es_time_reverse1: start running simulation backwards

init_energy_diag1: initialize energy diagnostic
energy_diag1: energy diagnostic
print_energy1: print energy summaries
del_energy_diag1: delete energy diagnostic

init_spectrum1: allocate scratch arrays for scalar fields
del_spectrum1: delete scratch arrays for scalar fields

init_edensity_diag1: initialize electron density diagnostic
edensity_diag1: electron density diagnostic
del_edensity_diag1: delete electron density diagnostic data

init_idensity_diag1: initialize ion density diagnostic
idensity_diag1: ion density diagnostic
del_idensity_diag1: delete ion density diagnostic

init_potential_diag1: initialize potential diagnostic
potential_diag1: potential diagnostic
del_potential_diag1: delete potential diagnostic

init_elfield_diag1: initialize longitudinal efield diagnostic
elfield_diag1: longitudinal efield diagnostic
del_elfield_diag1: delete longitudinal efield diagnostic

init_efluidms_diag1: initialize electron fluid moments diagnostic
efluidms_diag1: electron fluid moments diagnostic
del_efluidms_diag1: delete electron fluid moments diagnostic

init_ifluidms_diag1: initialize ion fluid moments diagnostic
ifluidms_diag1: ion fluid moments diagnostic
del_ifluidms_diag1: delete ion fluid moments diagnostic

init_evelocity_diag1: initialize electron velocity diagnostic
evelocity_diag1: electron velocity diagnostic
del_evelocity_diag1: delete electron velocity diagnostic

init_ivelocity_diag1: initialize ion velocity diagnostic
ivelocity_diag1: ion velocity diagnostic
del_ivelocity_diag1: delete ion velocity diagnostic

init_traj_diag1: initialize trajectory diagnostic
traj_diag1: trajectory diagnostic
del_traj_diag1: delete trajectory diagnostic

print_timings1: print timing summaries

reset_diags1: reset electrostatic diagnostics
close_diags1: close diagnostics

initialize_diagnostics1: initialize all diagnostics from namelist
                         input parameters
                         
open_restart1: open reset and restart files
bwrite_restart1: write out basic restart file for electrostatic code
bread_restart1: read in basic restart file for electrostatic code
dwrite_restart1: write out restart diagnostic file for electrostatic
                 code
dread_restart1: read in restart diagnostic file for electrostatic code
close_restart1: close reset and restart files

written by Viktor K. Decyk and Joshua Kelly, UCLA
copyright 1999-2016, regents of the university of california
update: december 9, 2017
"""
import math
import numpy

# sys.path.append('./mbeps1.source')
from libmpush1 import *
from dtimer import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# idimp = number of particle coordinates = 2
# ipbc = particle boundary condition: 1 = periodic
idimp = 2; ipbc = 1
# wke/wki/we = electron/ion kinetic energies and electric field energy
wke = numpy.zeros((1),float_type)
wki = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
# plist = (true,false) = list of particles leaving tiles found in push
plist = True

# declare scalars for standard code
ntime0 = 0; npi = 0
ws = numpy.zeros((1),float_type)
nt = numpy.zeros((1),int_type)

# declare scalars for OpenMP code
irc = numpy.zeros((1),int_type)
irc2 = numpy.zeros((2),int_type)

# declare scalars for diagnostics
itw = 0; itp = 0; itv = 0; itt = 0
itdi = 0
# default Fortran unit numbers
iuin = 8; iudm = 19
iude = 10; iup = 11; iuel = 12
iufe = 23; iufi = 24
iudi = 20

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
dtime = numpy.empty((1),double_type)
tdpost = numpy.zeros((1),float_type)
tguard = numpy.zeros((1),float_type)
tfft = numpy.zeros((1),float_type)
tfield = numpy.zeros((1),float_type)
tpush = numpy.zeros((1),float_type)
tsort = numpy.zeros((1),float_type)
tdiag = numpy.zeros((1),float_type)

# create string from idrun
cdrun = str(in1.idrun)

# initialize scalars for standard code
# increase number of coordinates for particle tag
if ((in1.ntt > 0) or ((in1.nts > 0) and (in1.ntsc > 0))):
   idimp += 1
# np = total number of electrons in simulation
# nx = number of grid points in x direction
np = in1.npx + in1.npxb;
nx = int(math.pow(2,in1.indx)); nxh = int(nx/2)
# npi = total number of ions in simulation
if (in1.movion > 0):
   npi = in1.npxi + in1.npxbi
nxe = nx + 2; nxeh = nxe/2
# mx1 = number of tiles in x direction
mx1 = int((nx - 1)/in1.mx + 1)
# nloop = number of time steps in simulation
nloop = int(in1.tend/in1.dt + .0001)
qbme = in1.qme
affp = float(nx)/float(np)
if (in1.movion==1):
   qbmi = in1.qmi/in1.rmass
   vtxi = in1.vtx/numpy.sqrt(in1.rmass*in1.rtempxi)
   vtdxi = in1.vtdx/numpy.sqrt(in1.rmass*in1.rtempdxi)

# check for unimplemented features
if (ipbc != 1):
   plist = False

# cwk = labels for power spectrum display
cwk = numpy.array([" W > 0"," W < 0"],'S6')

#-----------------------------------------------------------------------
def init_fields1():
   """ allocate field data for standard code"""
   global qe, qi, fxe, ffc, mixup, sct
# qe = electron charge density with guard cells
   qe = numpy.empty((nxe),float_type,'F')
# qi = ion charge density with guard cells
   qi = numpy.empty((nxe),float_type,'F')
# fxe = smoothed electric field with guard cells
   fxe = numpy.empty((nxe),float_type,'F')
# ffc = form factor array for poisson solver
   ffc = numpy.empty((nxh),complex_type,'F')
# mixup = bit reverse table for FFT
   mixup = numpy.empty((nxh),int_type,'F')
# sct = sine/cosine table for FFT
   sct = numpy.empty((nxh),complex_type,'F')

#-----------------------------------------------------------------------
def del_fields1():
   """ delete field data for standard code """
   global qe, qi, fxe, ffc, mixup, sct
   del qe, qi, fxe, ffc, mixup, sct

#-----------------------------------------------------------------------
def init_electrons1():
   """ initialize electrons """
   global part, nppmx, kpic, ppart, ppbuff, ncl, ihole
# part = particle array
   part = numpy.empty((idimp,max(np,npi)),float_type,'F')
# background electrons
   if (in1.npx > 0):
# calculates initial particle co-ordinates with various density profiles
      minit1.mfdistr1(part,in1.ampdx,in1.scaledx,in1.shiftdx,1,in1.npx,
                      nx,ipbc,in1.ndprof,irc)
      if (irc[0] != 0): exit(1)
# initialize particle velocities
      minit1.wmvdistr1(part,1,in1.vtx,in1.vx0,in1.ci,in1.npx,in1.nvdist,
                       in1.relativity,irc)
      if (irc[0] != 0): exit(1)
# beam electrons
   if (in1.npxb > 0):
      it = in1.npx + 1
# calculates initial particle co-ordinates with various density profiles
      minit1.mfdistr1(part,in1.ampdx,in1.scaledx,in1.shiftdx,it,
                      in1.npxb,nx,ipbc,in1.ndprof,irc)
      if (irc[0] != 0): exit(1)
# initialize particle velocities
      minit1.wmvdistr1(part,it,in1.vtdx,in1.vdx,in1.ci,in1.npxb,
                       in1.nvdist,in1.relativity,irc)
      if (irc[0] != 0): exit(1)

# mark electron beam particles
   if ((in1.nts > 0) and (in1.ntsc > 0)):
      mdiag1.setmbeam1(part,in1.npx,irc)
      if (irc[0] != 0): exit(1)

# nppmx = maximum number of particles in tile
   nppmx = numpy.empty((1),int_type)
# kpic = number of electrons in each tile
   kpic = numpy.empty((mx1),int_type,'F')

# find number of electrons in each of mx, tiles: updates kpic, nppmx
   minit1.mdblkp2(part,kpic,nppmx,np,in1.mx,irc)
   if (irc[0] != 0): exit(1)

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
   if (irc[0] != 0): exit(1)

# sanity check for electrons
   mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)
   if (irc[0] != 0): exit(1)

#-----------------------------------------------------------------------
def reorder_electrons1(irc2):
   """
   recover from electron wmporder1 errors
   reallocates ihole, ppbuff, and ppart as necessary
   input/output:
   irc2 = error codes, returned only if error occurs, when irc2(1) != 0
   """
# local data
   global ntmax, npbmx, nppmx0, ihole, ppbuff, ppart
   nter = 10; iter = 0
   while (irc2[0] != 0):
      iter += 1
      if (iter > nter):
         print "reorder_electrons1: iteration exceeded"
         exit(1)
# ihole overflow
      if (irc2[0]==1):
         ntmax = int((1.0 + in1.xtras)*irc2[1])
         del ihole
         ihole = numpy.empty((2,ntmax+1,mx1),int_type,'F')
         irc2[:] = 0
         msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,
                          False,irc2)
# ppbuff overflow
      elif (irc2[0]==2):
         npbmx = int((1.0 + in1.xtras)*irc2[1])
         del ppbuff
         ppbuff = numpy.empty((idimp,npbmx,mx1),float_type,'F')
# restores initial values of ncl
         msort1.mprsncl1(ncl,tsort)
         irc2[:] = 0
         msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,
                          True,irc2)
# ppart overflow
      elif (irc2[0]==3):
# restores particle coordinates from ppbuff: updates ppart, ncl
         msort1.mprstor1(ppart,ppbuff,ncl,ihole,tsort)
# copy ordered particles to linear array: updates part
         mpush1.mpcopyout1(part,ppart,kpic,nt,irc)
         del ppart
         nppmx0 = int((1.0 + in1.xtras)*irc2[1])
         ppart = numpy.empty((idimp,nppmx0,mx1),float_type,'F')
# copies unordered particles to ordered array: updates ppart
         mpush1.mpcopyin1(part,ppart,kpic,irc)
         irc2[:] = 0
         msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,
                          plist,irc2)
# electrons not in correct tile, try again
      elif (irc2[0]==(-1)):
         irc2[:] = 0
         msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,
                          False,irc2)

# sanity check for electrons
   if (in1.monitor > 0):
      mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)
      if (irc[0] != 0): exit(1)

#-----------------------------------------------------------------------
def push_electrons1(ppart,kpic):
   """
   push electrons with OpenMP:
   input/output:
   ppart = tiled electron particle array
   kpic = number of electrons in each tile
   """
   global wke
   wke[0] = 0.0
# updates part, wke and possibly ncl, ihole, and irc
   if (in1.mzf==0):
      mpush1.wmpush1(ppart,fxe,kpic,ncl,ihole,qbme,in1.dt,in1.ci,wke,
                     tpush,nx,in1.mx,ipbc,in1.relativity,plist,irc)
# zero force: updates part, wke and possibly ncl, ihole, and irc
   else:
      mpush1.wmpush1zf(ppart,kpic,ncl,ihole,in1.dt,in1.ci,wke,tpush,nx,
                       in1.mx,ipbc,in1.relativity,plist,irc)

# reorder electrons by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   if (irc[0]==0):
      msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,
                       plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc; irc[0] = 0

# sanity check for electrons
   if (irc2[0]==0):
      if (in1.monitor > 0):
         mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)
         if (irc[0] != 0): exit(1)
# recover from wmporder1 errors: updates ppart
   elif (irc2[0] != 0):
      reorder_electrons1(irc2)

#-----------------------------------------------------------------------
def del_electrons1():
   """ delete electrons """
   global ppart, kpic
   del ppart, kpic

#-----------------------------------------------------------------------
def init_ions1():
   """ initialize ions """
   global kipic, pparti
# part = particle array
# kipic = number of ions in each tile
   kipic = numpy.empty((mx1),int_type,'F')
   it = in1.npxi + 1
# background ions
   if (in1.npxi > 0):
      minit1.mfdistr1(part,in1.ampdxi,in1.scaledxi,in1.shiftdxi,1,
                      in1.npxi,nx,ipbc,in1.ndprofi,irc)
      if (irc[0] != 0): exit(1)
      minit1.wmvdistr1(part,1,vtxi,in1.vxi0,in1.ci,in1.npxi,in1.nvdist,
                       in1.relativity,irc)
      if (irc[0] != 0): exit(1)
# beam ions
   if (in1.npxbi > 0):
      minit1.mfdistr1(part,in1.ampdxi,in1.scaledxi,in1.shiftdxi,it,
                      in1.npxbi,nx,ipbc,in1.ndprofi,irc)
      if (irc[0] != 0): exit(1)
      minit1.wmvdistr1(part,it,vtdxi,in1.vdxi,in1.ci,in1.npxbi,
                       in1.nvdist,in1.relativity,irc)
      if (irc[0] != 0): exit(1)
      
# marks ion beam particles
   if ((in1.nts > 0) and (in1.ntsc > 0)):
      mdiag1.setmbeam1(part,in1.npxi,irc)
      if (irc[0] != 0): exit(1)

# kipic = number of ions in each tile
   kipic = numpy.empty((mx1),int_type,'F')

# find number of ions in each of mx, tiles: updates kipic, nppmx
   minit1.mdblkp2(part,kipic,nppmx,npi,in1.mx,irc)
   if (irc[0] != 0): exit(1)

# allocate vector ion data
   nppmx1 = int((1.0 + in1.xtras)*nppmx)
   pparti = numpy.empty((idimp,nppmx1,mx1),float_type,'F')

# copy ordered ion data for OpenMP: updates pparti and kipic
   mpush1.mpmovin1(part,pparti,kipic,in1.mx,irc)
   if (irc[0] != 0): exit(1)

# sanity check for ions
   mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)
   if (irc[0] != 0): exit(1)

#-----------------------------------------------------------------------
def reorder_ions1(irc2):
   """
   recover from ion wmporder1 errors
   reallocates ihole, ppbuff, and pparti as necessary
   input/output:
   irc2 = error codes, returned only if error occurs, when irc2(1) != 0
   """
# local data
   global ntmax, npbmx, nppmx1, ihole, ppbuff, pparti
   nter = 10; iter = 0
   while (irc2[0] != 0):
      iter += 1
      if (iter > nter):
         print "reorder_ions1: iteration exceeded"
         exit(1)
# ihole overflow
      if (irc2[0]==1):
         ntmax = int((1.0 + in1.xtras)*irc2[1])
         del ihole
         ihole = numpy.empty((2,ntmax+1,mx1),int_type,'F')
         irc2[:] = 0
         msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,in1.mx,
                          False,irc2)
# ppbuff overflow
      elif (irc2[0]==2):
         npbmx = int((1.0 + in1.xtras)*irc2[1])
         del ppbuff
         ppbuff = numpy.empty((idimp,npbmx,mx1),float_type,'F')
# restores initial values of ncl
         msort1.mprsncl1(ncl,tsort)
         irc2[:] = 0
         msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,in1.mx,
                          True,irc2)
# pparti overflow
      elif (irc2[0]==3):
# restores particle coordinates from ppbuff: updates ppart, ncl
         msort1.mprstor1(pparti,ppbuff,ncl,ihole,tsort)
# copy ordered particles to linear array: updates part
         mpush1.mpcopyout1(part,pparti,kipic,nt,irc)
         del pparti
         nppmx1 = int((1.0 + in1.xtras)*irc2[1])
         pparti = numpy.empty((idimp,nppmx1,mx1),float_type,'F')
# copies unordered particles to ordered array: updates ppart
         mpush1.mpcopyin1(part,pparti,kipic,irc)
         irc2[:] = 0
         msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,in1.mx,
                          plist,irc2)
# ions not in correct tile, try again
      elif (irc2[0]==(-1)):
         irc2[:] = 0
         msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,in1.mx,
                          False,irc2)

# sanity check for ions
   if (in1.monitor > 0):
      mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)
      if (irc[0] != 0): exit(1)

#-----------------------------------------------------------------------
def push_ions1(pparti,kipic):
   """
   push ions with OpenMP
   input/output:
   pparti = tiled electron/ion particle arrays
   kipic = number of electrons/ions in each tile
   """
   global wki
   wki[0] = 0.0
# updates pparti, wki and possibly ncl, ihole, and irc
   if (in1.mzf==0):
      mpush1.wmpush1(pparti,fxe,kipic,ncl,ihole,qbmi,in1.dt,in1.ci,wki,
                     tpush,nx,in1.mx,ipbc,in1.relativity,plist,irc)
# zero force: updates pparti, wki and possibly ncl, ihole, and irc
   else:
      mpush1.wmpush1zf(pparti,kipic,ncl,ihole,in1.dt,in1.ci,wki,tpush,
                       nx,in1.mx,ipbc,in1.relativity,plist,irc)
   wki[0] = wki[0]*in1.rmass

# reorder ions by tile with OpenMP:
# updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
   if (irc[0]==0):
      msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,in1.mx,
                       plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc; irc[0] = 0

# sanity check for ions
   if (irc2[0]==0):
      if (in1.monitor > 0):
         mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)
         if (irc[0] != 0): exit(1)
# recover from wmporder1 errors: updates pparti
   elif (irc2[0] != 0):
      reorder_ions1(irc2)

#-----------------------------------------------------------------------
def del_ions1():
   """ delete ions """
   global pparti, kipic
   del pparti, kipic

#-----------------------------------------------------------------------
def es_time_reverse1():
   """
   start running simulation backwards:
   need to reverse time lag in leap-frog integration scheme
   """
# ppart/pparti = tiled electron/ion particle arrays
# kpic/kipic = number of electrons/ions in each tile
   global ws
   in1.dt = -in1.dt
   ws[0] = 0.0
   mpush1.wmpush1zf(ppart,kpic,ncl,ihole,in1.dt,in1.ci,ws,tpush,nx,
                    in1.mx,ipbc,in1.relativity,plist,irc)
   msort1.wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,in1.mx,plist,
                    irc2)
   if (irc2[0] != 0): exit(1)
   if (in1.movion==1):
      mpush1.wmpush1zf(pparti,kipic,ncl,ihole,in1.dt,in1.ci,ws,tpush,nx,
                       in1.mx,ipbc,in1.relativity,plist,irc)
      msort1.wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,in1.mx,
                       plist,irc2)
      if (irc2[0] != 0): exit(1)

#-----------------------------------------------------------------------
def init_energy_diag1():
   """ initialize energy diagnostic"""
# wt = energy time history array
   global mtw, itw, wt, s
   mtw = int((nloop - 1)/in1.ntw + 1); itw = 0
   wt = numpy.zeros((mtw,4),float_type,'F')
   s = numpy.zeros((4),double_type,'F')

#-----------------------------------------------------------------------
def energy_diag1(wt,ntime,iuot):
   """
   energy diagnostic
   input/output:
   wt = energy time history array
   input:
   ntime = current time step
   iuot = output file descriptor
   """
   global ws, s, itw
   ws[0] = we[0] + wke[0] + wki[0]
   if (ntime==0):
      s[2] = ws[0]
   if (in1.ndw > 0):
      print >> iuot, "Field, Kinetic and Total Energies:"
      if (in1.movion==0):
         iuot.write("%14.7e %14.7e %14.7e\n" % (we[0],wke[0],ws[0])) 
      else:
         iuot.write("%14.7e %14.7e %14.7e %14.7e\n" % (we[0],wke[0],
                    wki[0],ws[0])) 
   wt[itw,:] = [we[0],wke[0],wki[0],ws[0]]
   itw += 1
   s[0] += we[0]
   s[1] += wke[0]
   s[2] = min(s[2],float(ws[0]))
   s[3] = max(s[3],float(ws[0]))

#-----------------------------------------------------------------------
def print_energy1(wt,iuot):
   """
   print energy summaries
   input:
   wt = energy time history array
   iuot = output file descriptor
   """
   global s, itw
   swe = s[0]; swke = s[1]
   s[2] = (s[3] - s[2])/wt[0,3]
   print >> iuot, "Energy Conservation = ", float(s[2])
   swe = swe/float(itw)
   print >> iuot, "Average Field Energy <WE> = ", float(swe)
   swke = swke/float(itw)
   print >> iuot, "Average Electron Kinetic Energy <WKE> = ",float(swke)
   print >> iuot, "Ratio <WE>/<WKE>= ", float(swe/swke)
   print >> iuot

#-----------------------------------------------------------------------
def del_energy_diag1():
   """ delete energy diagnostic """
# wt = energy time history array
   global wt, s
   del wt, s

#-----------------------------------------------------------------------
def init_spectrum1():
   """ allocate scratch arrays for scalar fields """
   global sfield, sfieldc
   sfield = numpy.empty((nxe),float_type,'F')
   sfieldc = numpy.empty((nxh),complex_type,'F')
# allocate and initialize frequency array for spectral analysis
   if (in1.ntp > 0):
      global iw, wm
      iw = int((in1.wmax - in1.wmin)/in1.dw + 1.5)
      wm = numpy.empty((iw),float_type,'F')
      wm[:] = in1.wmin + in1.dw*numpy.linspace(0,iw-1,iw)
# allocate and initialize frequency array for ion spectral analysis
   if (in1.movion==1):
      if (in1.ntdi > 0):
         global iwi, wmi
         iwi = int((in1.wimax - in1.wimin)/in1.dwi + 1.5)
         wmi = numpy.empty((iwi),float_type,'F')
         wmi[:] = in1.wimin + in1.dwi*numpy.linspace(0,iwi-1,iwi)

#-----------------------------------------------------------------------
def del_spectrum1():
   """ delete scratch arrays for scalar fields """
   global sfield, sfieldc
   del sfield, sfieldc
   if (in1.ntp > 0):
      global wm
      del wm
   if (in1.movion==1):
      if (in1.ntdi > 0):
         global wmi
         del wmi

#-----------------------------------------------------------------------
def init_edensity_diag1():
   """ initialize electron density diagnostic """
   global denet, iude
   fdename = "denek1." + cdrun
   in1.fdename[:] = fdename
   in1.modesxde = int(min(in1.modesxde,nxh+1))
# denet = store selected fourier modes for electron density
   denet = numpy.empty((in1.modesxde),complex_type,'F')
# open file: updates nderec and possibly iude
   if (in1.nderec==0):
      mdiag1.dafopenc1(denet,iude,in1.nderec,fdename)

#-----------------------------------------------------------------------
def edensity_diag1(sfield):
   """
   electron density diagnostic
   input/output:
   sfield = scratch array for scalar field
   """
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

#-----------------------------------------------------------------------
def del_edensity_diag1():
   """ delete electron density diagnostic data """
   global denet
   if (in1.nderec > 0):
      in1.closeff(iude)
      in1.nderec -= 1
   del denet

#-----------------------------------------------------------------------
def init_idensity_diag1():
   """ initialize ion density diagnostic """
   global denit, iudi
   fdiname = "denik1." + cdrun
   in1.fdiname[:] = fdiname
   in1.modesxdi = int(min(in1.modesxdi,nxh+1))
# denit = store selected fourier modes for ion density
   denit = numpy.empty((in1.modesxdi),complex_type,'F')
# open file: updates ndirec and possibly iudi
   if (in1.ndirec==0):
      mdiag1.dafopenc1(denit,iudi,in1.ndirec,fdiname)
# ion spectral analysis
   global mtdi, itdi, pkwdi, pksdi, wkdi
   if ((in1.nddi==2) or (in1.nddi==3)):
      mtdi = int((nloop - 1)/in1.ntdi) + 1; itdi = 0
# pkwdi = power spectrum for potential
      pkwdi = numpy.empty((in1.modesxdi,iwi,2),float_type,'F')
# pksdi = accumulated complex spectrum for potential
      pksdi = numpy.zeros((4,in1.modesxdi,iwi),double_type,'F')
# wkdi = maximum frequency as a function of k for potential
      wkdi = numpy.empty((in1.modesxdi,2),float_type,'F')
# create dummy arrays to avoid undefined arguments later
   else:
      pkwdi = numpy.zeros((1,1,1),float_type,'F')
      wkdi = numpy.zeros((1,1),float_type,'F')

#-----------------------------------------------------------------------
def idensity_diag1(sfield,pkwdi,wkdi,ntime):
   """
   ion density diagnostic
   input/output:
   sfield = scratch array for scalar field
   pkwdi = power spectrum for ion density
   wkdi = maximum frequency as a function of k for ion density
   input:
   ntime = current time step
   """
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
# ion spectral analysis
   if ((in1.nddi==2) or (in1.nddi==3)):
      global itdi
      itdi += 1
      ts = in1.dt*float(ntime)
# performs frequency analysis of accumulated complex time series
# zero out mode 0
      denit[0] = numpy.complex(0.0,0.0)
      mdiag1.micspect1(denit,wmi,pkwdi,pksdi,ts,in1.t0,tdiag,mtdi,iwi,
                       in1.modesxdi,nx,-1)
# find frequency with maximum power for each mode
      wkdi[:,0] = wmi[numpy.argmax(pkwdi[:,:,0],axis=1)]
      wkdi[:,1] = wmi[numpy.argmax(pkwdi[:,:,1],axis=1)]

#-----------------------------------------------------------------------
def del_idensity_diag1():
   """ delete ion density diagnostic data """
   global denit
   if (in1.ndirec > 0):
      in1.closeff(iudi)
      in1.ndirec -= 1
   del denit
# spectral analysis
   if ((in1.nddi==2) or (in1.nddi==3)):
      global pkwdi, wkdi, pksdi
      del pkwdi, wkdi, pksdi

#-----------------------------------------------------------------------
def init_potential_diag1():
   """ initialize potential diagnostic """
   global pott, iup
   fpname = "potk1." + cdrun
   in1.fpname[:] = fpname
   in1.modesxp = int(min(in1.modesxp,nxh+1))
# pott = store selected fourier modes for potential
   pott = numpy.empty((in1.modesxp),complex_type,'F')
# open file: updates nprec and possibly iup
   if (in1.nprec==0):
      mdiag1.dafopenc1(pott,iup,in1.nprec,fpname)
# potential spectral analysis
   global mtp, itp, pkw, pks, wk
   if ((in1.ndp==2) or (in1.ndp==3)):
      mtp = int((nloop - 1)/in1.ntp) + 1; itp = 0
# pkw = power spectrum for potential
      pkw = numpy.empty((in1.modesxp,iw,2),float_type,'F')
# pks = accumulated complex spectrum for potential
      pks = numpy.zeros((4,in1.modesxp,iw),double_type,'F')
# wk = maximum frequency as a function of k for potential
      wk = numpy.empty((in1.modesxp,2),float_type,'F')
# create dummy arrays to avoid undefined arguments later
   else:
      pkw = numpy.zeros((1,1,1),float_type,'F')
      wk = numpy.zeros((1,1),float_type,'F')

#-----------------------------------------------------------------------
def potential_diag1(sfield,pkw,wk,ntime):
   """
   potential diagnostic
   input/output:
   sfield = scratch array for scalar field
   pkw = power spectrum for potential
   wk = maximum frequency as a function of k for potential
   input:
   ntime = current time step
   """
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
# potential spectral analysis
   if ((in1.ndp==2) or (in1.ndp==3)):
      global itp
      itp += 1
      ts = in1.dt*float(ntime)
# performs frequency analysis of accumulated complex time series
      mdiag1.micspect1(pott,wm,pkw,pks,ts,in1.t0,tdiag,mtp,iw,
                       in1.modesxp,nx,1)
# find frequency with maximum power for each mode
      wk[:,0] = wm[numpy.argmax(pkw[:,:,0],axis=1)]
      wk[:,1] = wm[numpy.argmax(pkw[:,:,1],axis=1)]

#-----------------------------------------------------------------------
def del_potential_diag1():
   """ delete potential diagnostic """
   global pott
   if (in1.nprec > 0):
      in1.closeff(iup)
      in1.nprec -= 1
   del pott
# spectral analysis
   if ((in1.ndp==2) or (in1.ndp==3)):
      global pkw, wk, pks
      del pkw, wk, pks
   in1.ceng = affp

#-----------------------------------------------------------------------
def init_elfield_diag1():
   """ initialize longitudinal efield diagnostic """
   global elt, iuel
   felname = "elk1." + cdrun
   in1.felname[:] = felname
   in1.modesxel = int(min(in1.modesxel,nxh+1))
# elt = store selected fourier modes for longitudinal efield
   elt = numpy.empty((in1.modesxel),complex_type,'F')
# open file: updates nelrec and possibly iuel
   if (in1.nelrec==0):
      mdiag1.dafopenc1(elt,iuel,in1.nelrec,felname)

#-----------------------------------------------------------------------
def elfield_diag1(sfield):
   """
   longitudinal efield diagnostic
   input/output:
   sfield = scratch array for scalar field
   """
# calculate longitudinal efield in fourier space: updates sfieldc
   mfield1.melfield1(qe,sfieldc,ffc,ws,tfield,nx)
# store selected fourier modes: updates elt
   mfield1.mrdmodes1(sfieldc,elt,tfield,nx,in1.modesxel)
# write diagnostic output: updates nelrec
   mdiag1.dafwritec1(elt,tdiag,iuel,in1.nelrec,in1.modesxel)
# transform longitudinal efield to real space: updates sfield
   mfft1.mfft1cr(sfieldc,sfield,mixup,sct,tfft,in1.indx)
   mgard1.mdguard1(sfield,tguard,nx)

#-----------------------------------------------------------------------
def del_elfield_diag1():
   """ delete longitudinal efield diagnostic """
   global elt
   if (in1.nelrec > 0):
      in1.closeff(iuel)
      in1.nelrec -= 1
   del elt
   in1.ceng = affp

#-----------------------------------------------------------------------
def init_efluidms_diag1():
   """ initialize electron fluid moments diagnostic """
   global fmse, iufe
# calculate first dimension of fluid arrays
   if (in1.npro==1):
      in1.nprd = 1
   elif (in1.npro==2):
      in1.nprd = 2
   elif (in1.npro==3):
      in1.nprd = 3
   elif (in1.npro==4):
      in1.nprd = 5
   if ((in1.ndfm==1) or (in1.ndfm==3)):
# fmse = electron fluid moments
      fmse = numpy.empty((in1.nprd,nxe),float_type,'F')
# open file for real data: updates nferec and possibly iufe
      ffename = "fmer1." + cdrun
      in1.ffename[:] = ffename
      if (in1.nferec==0):
         mdiag1.dafopenv1(fmse,nx,iufe,in1.nferec,ffename)

#-----------------------------------------------------------------------
def efluidms_diag1(fmse):
   """
   electron fluid moments diagnostic
   input/output:
   fmse = electron fluid moments
   """
# calculate electron fluid moments
   if ((in1.ndfm==1) or (in1.ndfm==3)):
      dtimer(dtime,itime,-1)
      fmse.fill(0.0)
      dtimer(dtime,itime,1)
      tdiag[0] += float(dtime)
      mdiag1.wmgprofx1(ppart,fxe,fmse,kpic,qbme,in1.dt,in1.ci,tdiag,
                       in1.npro,nx,in1.mx,in1.relativity)
# add guard cells with OpenMP: updates fmse
      mgard1.mamcguard1(fmse,tdiag,nx)
# calculates fluid quantities from fluid moments: updates fmse
      mdiag1.mfluidqs1(fmse,tdiag,in1.npro,nx)
# write real space diagnostic output: updates nferec
      mdiag1.dafwritev1(fmse,tdiag,iufe,in1.nferec,nx)

#-----------------------------------------------------------------------
def del_efluidms_diag1():
   """ delete electron fluid moments diagnostic """
   global fmse
   if (in1.nferec > 0):
      in1.closeff(iufe)
      in1.nferec -= 1
   del fmse

#-----------------------------------------------------------------------
def init_ifluidms_diag1():
   """ initialize ion fluid moments diagnostic """
   global fmsi, iufi
   if ((in1.ndfm==2) or (in1.ndfm==3)):
# fmsi = ion fluid moments
      fmsi = numpy.empty((in1.nprd,nxe),float_type,'F')
# open file for real data: updates nfirec and possibly iufi
      ffiname = "fmir1." + cdrun
      in1.ffiname[:] = ffiname
      if (in1.nfirec==0):
          mdiag1.dafopenv1(fmsi,nx,iufi,in1.nfirec,ffiname)

#-----------------------------------------------------------------------
def ifluidms_diag1(fmsi):
   """
   ion fluid moments diagnostic
   input/output:
   fmsi = ion fluid moments
   """
# calculate ion fluid moments
   if ((in1.ndfm==2) or (in1.ndfm==3)):
      dtimer(dtime,itime,-1)
      fmsi.fill(0.0)
      dtimer(dtime,itime,1)
      tdiag[0] += float(dtime)
      mdiag1.wmgprofx1(pparti,fxe,fmsi,kipic,qbmi,in1.dt,in1.ci,tdiag,
                       in1.npro,nx,in1.mx,in1.relativity)
# add guard cells with OpenMP: updates fmsi
      mgard1.mamcguard1(fmsi,tdiag,nx)
# calculates fluid quantities from fluid moments: updates fmsi
      mdiag1.mfluidqs1(fmsi,tdiag,in1.npro,nx)
      fmsi[:,:]  = in1.rmass*fmsi
# write real space diagnostic output: updates nfirec
      mdiag1.dafwritev1(fmsi,tdiag,iufi,in1.nfirec,nx)

#-----------------------------------------------------------------------
def del_ifluidms_diag1():
   """ ion fluid moments diagnostic """
   global fmsi
   if (in1.nfirec > 0):
      in1.closeff(iufi)
      in1.nfirec -= 1
   del fmsi

#-----------------------------------------------------------------------
def init_evelocity_diag1():
   """ initialize electron velocity diagnostic """
   global fv, sfv, fvm, mtv, itv, fvtm
   mtv = int((nloop - 1)/in1.ntv) + 1; itv = 0
# fv = global electron velocity distribution functions
   fv = numpy.empty((2*in1.nmv+2,in1.ndim),float_type,'F')
# sfv = electron velocity distribution functions in tile
   sfv = numpy.empty((2*in1.nmv+2,in1.ndim,mx1+1),float_type,'F')
# fvm = electron vdrift, vth, entropy for global distribution
   fvm = numpy.empty((in1.ndim,3),float_type,'F')
# fvtm = time history of electron vdrift, vth, and entropy
   fvtm = numpy.zeros((mtv,in1.ndim,3),float_type,'F')
   sfv[0,:,:] = 2.0*max(4.0*in1.vtx+abs(in1.vx0),
                           4.0*in1.vtdx+abs(in1.vdx))

#-----------------------------------------------------------------------
def evelocity_diag1(ppart,kpic,fv,fvm,fvtm):
   """
   electron velocity diagnostic
   input:
   ppart = tiled electron particle arrays
   kpic = number of electrons in each tile
   input/output:
   fv = global electron velocity distribution functions
   fvmi = electron vdrift, vth, entropy for global distribution
   fvtm = time history of electron vdrift, vth, and entropy
   """
   global itv
# calculate electron distribution function and moments
   mdiag1.mvpdist1(ppart,kpic,sfv,fvm,tdiag,np,in1.nmv)
   fv[:,:] = sfv[:,:,mx1]
# store time history electron vdrift, vth, and entropy
   fvtm[itv,:,:] = fvm
   itv += 1

#-----------------------------------------------------------------------
def del_evelocity_diag1():
   """ delete electron velocity diagnostic """
   global fv, sfv, fvm, fvtm
   del fv, sfv, fvm, fvtm

#-----------------------------------------------------------------------
def init_ivelocity_diag1():
   """ initialize ion velocity diagnostic """
   global fvi, sfvi, fvmi, fvtmi
# fvi = global ion velocity distribution functions
   fvi = numpy.empty((2*in1.nmv+2,in1.ndim),float_type,'F')
# sfvi = ion velocity distribution functions in tile
   sfvi = numpy.empty((2*in1.nmv+2,in1.ndim,mx1+1),float_type,'F')
# fvmi = ion vdrift, vth, entropy for global distribution
   fvmi = numpy.empty((in1.ndim,3),float_type,'F')
# fvtmi = time history of ion vdrift, vth, and entropy
   fvtmi = numpy.zeros((mtv,in1.ndim,3),float_type,'F')
   sfvi[0,:,:] = 2.0*max(4.0*vtxi+abs(in1.vxi0),4.0*vtdxi+abs(in1.vdxi))

#-----------------------------------------------------------------------
def ivelocity_diag1(pparti,kipic,fvi,fvmi,fvtmi):
   """
   ion velocity diagnostic
   input:
   pparti = tiled ion particle arrays
   kipic = number of ions in each tile
   input/output:
   fvi = global ion velocity distribution functions
   fvmi = ion vdrift, vth, entropy for global distribution
   fvtmi = time history of ion vdrift, vth, and entropy
   """
   mdiag1.mvpdist1(pparti,kipic,sfvi,fvmi,tdiag,npi,in1.nmv)
   fvi[:,:] = sfvi[:,:,mx1]
# update time step if electrons have not been calculated
   if (in1.ndv==2):
      s1.itv += 1
# store time history ion vdrift, vth, and entropy
   fvtmi[itv-1,:,:] = fvmi

#-----------------------------------------------------------------------
def del_ivelocity_diag1():
   """ delete electron velocity diagnostic """
   global fvi, sfvi, fvmi, fvtmi
   del fvi, sfvi, fvmi, fvtmi

#-----------------------------------------------------------------------
def init_traj_diag1(ntime):
   """
   initialize trajectory diagnostic
   ntime = current time step
   """
   global partt, fvtp, fvmtp, mtt, itt, partd
# set initial test trajectories
   if ((ntime+ntime0)==0):
# iprobt = scratch array 
      iprobt = numpy.empty((in1.nprobt),numpy.int32)
      mdiag1.setptraj1(ppart,kpic,iprobt,in1.nst,in1.vtx,in1.vtsx,
                       in1.dvtx,np,in1.nprobt,irc)
      if (irc[0] != 0): exit(1)
      if (in1.nprobt > 16777215):
         print "nprobt overflow = ", in1.nprobt
         exit(1)
      del iprobt
# find number of existing test tractories: updates nprobt
   else:
      mdiag1.mfnptraj1(ppart,kpic,in1.nprobt,irc)
      if (irc[0] != 0): exit(1)
# partt = particle trajectories tracked
   partt = numpy.empty((idimp,in1.nprobt),float_type,'F')
   if ((in1.nst==1) or (in1.nst==2)):
      mtt = int((nloop - 1)/in1.ntt + 1); itt = 0
# partd = trajectory time history array
      partd = numpy.empty((mtt,idimp,in1.nprobt),float_type,'F')
# create dummy array to avoid undefined arguments later
   else:
      partd = numpy.empty((1,1,1),float_type,'F')
   if (in1.nst==3):
# fvtp = velocity distribution function for test particles
      fvtp = numpy.empty((2*in1.nmv+2,in1.ndim),float_type,'F')
# fvmtp = vdrift, vth, and entropy for test particles
      fvmtp = numpy.empty((in1.ndim,3),float_type,'F')
      fvtp[0,:] = 2.0*max(4.0*in1.vtx+abs(in1.vx0),
                          4.0*in1.vtdx+abs(in1.vdx))
# create dummy arrays to avoid undefined arguments later
   else:
      fvtp = numpy.zeros((1,1),float_type,'F')
      fvmtp = numpy.zeros((1,1),float_type,'F')

#-----------------------------------------------------------------------
def traj_diag1(ppart,kpic,partd,fvtp,fvmtp):
   """
   trajectory diagnostic
   input:
   ppart = tiled electron particle array
   kpic = number of electrons in each tile
   input/output:
   partd = trajectory time history array
   fvtp = velocity distribution function for test particles
   fvmtp = vdrift, vth, and entropy for test particles
   """
# copies trajectories to array partt
   mdiag1.mptraj1(ppart,kpic,partt,tdiag,irc)
   if (irc[0] != 0): exit(1)
   if ((in1.nst==1) or (in1.nst==2)):
      global itt
      partd[itt,:,:] = partt
      itt += 1
   elif (in1.nst==3):
# calculate particle distribution function and moments
      mdiag1.mvdist1(partt,fvtp,fvmtp,tdiag,in1.nprobt,in1.nmv)

#-----------------------------------------------------------------------
def del_traj_diag1():
   """ delete trajectory diagnostic """
   global partt, partd, fvtp, fvmtp
   del partt
   if ((in1.nst==1) or (in1.nst==2)):
      del partd
   elif (in1.nst==3):
      del fvtp, fvmtp

#-----------------------------------------------------------------------
def print_timings1(tinit,tloop,iuot):
   """
   print timing summaries
   input:
   tinit = initialization elapsed time
   tloop = loop elapsed time
   iuot = output file descriptor
   """
   print >> iuot
   print >> iuot, "initialization time = ", tinit
   print >> iuot, "deposit time = ", tdpost[0]
   print >> iuot, "guard time = ", tguard[0]
   print >> iuot, "solver time = ", tfield[0]
   print >> iuot, "fft time = ", tfft[0]
   print >> iuot, "push time = ", tpush[0]
   print >> iuot, "sort time = ", tsort[0]
   tfield[0] += tguard[0] + tfft[0]
   print >> iuot, "total solver time = ", tfield[0]
   time = tdpost[0] + tpush[0] + tsort[0]
   print >> iuot, "total particle time = ", time
   print >> iuot, "total diagnostic time = ", tdiag[0]
   ws[0] = time + tfield[0] + tdiag[0]
   tloop = tloop - ws[0]
   print >> iuot, "total and additional time = ", ws[0], ",", tloop
   print >> iuot
# summarize particle timings
   ws[0] = 1.0e+09/(float(nloop)*float(np+npi))
   print >> iuot, "Push Time (nsec) = ", tpush[0]*ws[0]
   print >> iuot, "Deposit Time (nsec) = ", tdpost[0]*ws[0]
   print >> iuot, "Sort Time (nsec) = ", tsort[0]*ws[0]
   print >> iuot, "Total Particle Time (nsec) = ", time*ws[0]
   print >> iuot

#-----------------------------------------------------------------------
def reset_diags1():
   """ reset electrostatic diagnostics """
   if (in1.ntw > 0):
      itw = 0
      if ("wt" in globals()):
         wt.fill(0.0)
      if ("s" in globals()):
         s.fill(0.0)
   if (in1.ntde > 0):
      if (in1.nderec > 1):
         in1.nderec = 1
   if (in1.ntp > 0):
      if (in1.nprec > 1):
         in1.nprec = 1
      if ((in1.ndp==2) or (in1.ndp==3)):
         itp = 0; pks.fill(0.0)
   if (in1.movion==1):
      if (in1.ntdi > 0):
         if (in1.ndirec > 1):
            in1.ndirec = 1
         if ((in1.nddi==2) or (in1.nddi==3)):
            itdi = 0; pksdi.fill(0.0)
   if (in1.ntel > 0):
      if (in1.nelrec > 1):
         in1.nelrec = 1
   if (in1.ntfm > 0):
      if (in1.nferec > 1):
         in1.nferec = 1
      if (in1.movion==1):
         if (in1.nfirec > 0):
             in1.nfirec = 1
   if (in1.ntv > 0):
      itv = 0
      if ("fvtm" in globals()):
         fvtm.fill(0.0)
      if (in1.movion==1):
         if ("fvtmi" in globals()):
            fvtmi.fill(0.0)
   if (in1.ntt > 0):
      itt = 0

#-----------------------------------------------------------------------
def close_diags1(iudm):
   """
   close diagnostics
   delete data, close fortran files, and write out diagnostic metafile
   input:
   iudm = diagnostic metafile fortran file descriptor
   """
# electron density diagnostic
   if (in1.ntde > 0):
      del_edensity_diag1()
# potential diagnostic
   if (in1.ntp > 0):
      del_potential_diag1()
# longitudinal efield diagnostic
   if (in1.ntel > 0):
      del_elfield_diag1()
# fluid moments diagnostic
   if (in1.ntfm > 0) :
# electrons
      del_efluidms_diag1()
# ions
      if (in1.movion==1):
         del_ifluidms_diag1()
# ion density diagnostic
   if (in1.movion==1):
      if (in1.ntdi > 0):
         del_idensity_diag1()
# write final diagnostic metafile
   in1.writnml1(iudm)
   in1.closeff(iudm)
# deallocate arrays
   del_fields1()
   del_electrons1()
   global ppbuff, ncl, ihole, nppmx
   if (in1.movion==1):
      del_ions1()
   del ppbuff, ncl, ihole, nppmx
   if ((in1.ntde>0) or (in1.ntp>0) or (in1.ntel>0) or (in1.ntdi>0)):
      del_spectrum1()
   if (in1.ntw > 0):
      del_energy_diag1()
   if (in1.ntv > 0):
      del_evelocity_diag1()
      if (in1.movion==1):
         del_ivelocity_diag1()
   if (in1.ntt > 0):
      del_traj_diag1()

#-----------------------------------------------------------------------
def initialize_diagnostics1(ntime):
   """
   initialize all diagnostics from namelist input parameters
   input:
   ntime = current time step
   """
# initialize energy diagnostic: allocates wt
   if (in1.ntw > 0):
      init_energy_diag1()

# allocate and initialize scratch arrays for scalar fields:
# allocates sfield
   if ((in1.ntde>0) or (in1.ntp>0) or (in1.ntel>0) or (in1.ntdi>0)):
      init_spectrum1()

# initialize electron density diagnostic
   if (in1.ntde > 0):
      init_edensity_diag1()

# initialize ion density diagnostic: allocates pkwdi, wkdi
   if (in1.movion==1):
      if (in1.ntdi > 0):
         init_idensity_diag1()

# initialize potential diagnostic: allocates pkw, wk
   if (in1.ntp > 0):
      init_potential_diag1()

# initialize longitudinal efield diagnostic
   if (in1.ntel > 0):
      init_elfield_diag1()

# initialize fluid moments diagnostic
   if (in1.ntfm > 0):
# electrons: allocates fmse
      init_efluidms_diag1()
# ions: allocates fmsi
      if (in1.movion==1):
         init_ifluidms_diag1()

# initialize velocity diagnostic:
   if (in1.ntv > 0):
# electrons: allocates fv, fvm, fvtm
      init_evelocity_diag1()
# ions: allocates fvi, fvmi, fvtmi
      if (in1.movion==1):
         init_ivelocity_diag1()

# initialize trajectory diagnostic: allocates partd, fvtp, fvmtp
   if (in1.ntt > 0):
      init_traj_diag1(ntime)

#-----------------------------------------------------------------------
def open_restart1():
   """ open reset and restart files """
# iur, iurr = restart, reset, old restart file descriptors
   global iur, iurr, iur0, i1, i2, i3, i4, a2, fname
# reset file
   fname = "reset1"
   iurr = open(fname,"wb+")
# restart file
   fname = "rstrt1." + cdrun
# start a new run from random numbers
   if (in1.nustrt==1):
      if (in1.ntr > 0):
         iur = open(fname,"wb+")
# continue a run which was interrupted
   elif (in1.nustrt==2):
      iur = open(fname,"rb+")
# start a new run with data from a previous run
   elif (in1.nustrt==0):
      if (in1.ntr > 0):
         iur = open(fname,"wb+")
      if (in1.idrun != in1.idrun0):
         cdrun0 = str(in1.idrun0)
         fname = "rstrt1." + cdrun0
         iur0 = open(fname,"rb+")
      else:
         print "restart warning: old, new idruns identical"
# allocate scratch numpy arrays
   i1 = numpy.zeros((1),int_type)
   i2 = numpy.zeros((2),int_type)
   i3 = numpy.zeros((3),int_type)
   i4 = numpy.zeros((4),int_type)
   a2 = numpy.zeros((2),float_type)
   fname = numpy.array([""],'S32')

#-----------------------------------------------------------------------
def bwrite_restart1(iur,ntime):
   """
   write out basic restart file for electrostatic code
   input:
   iur = restart file descriptor
   ntime = current time step
   """
   global i1, i2
   idimp = numpy.size(part,0)
   iur.seek(0,0)
# write out current and initial time
   i2[0] = ntime; i2[1] = ntime0
   i2.tofile(iur)
# copy ordered particles to linear array: updates part, npp
   mpush1.mpcopyout1(part,ppart,kpic,i1,irc)
   npp = i1[0]
# write out size of electron array
   i2[0] = npp; i2[1] = idimp
   i2.tofile(iur)
# write out electrons, if non-zero
   if (npp > 0):
      part.tofile(iur)
# write out if ions are moving
   i1[0] = in1.movion; i1.tofile(iur)
   if (in1.movion==1):
# copy ordered particles to linear array: updates part, npp
      mpush1.mpcopyout1(part,pparti,kipic,i1,irc)
      npp = i1[0]
# write out size of ion array
      i2[0] = npp; i2[1] = idimp
      i2.tofile(iur)
# write out ions, if non-zero
      if (npp > 0):
         part.tofile(iur)
# write out ion density, if ions are not moving
   else:
      nxv = numpy.size(qi)
      i1[0] = nxv; i1.tofile(iur)
      if (nxv > 0):
         qi.tofile(iur)
# write out electric field parameter
   i1[0] = in1.emf
   i1.tofile(iur)

#-----------------------------------------------------------------------
def bread_restart1(iur):
   """
   read in basic restart file for electrostatic code
   input:
   iur = restart file descriptor
   """
   global part, i1, i2, ntime, ntime0, nppmx, nppmx1, np, npi
   global ppart, kpic, ppbuff, ncl, ihole, pparti, kipic
# part = linear particle array
   if ("part" not in globals()):
      part = numpy.empty((idimp,max(np,npi)),float_type,'F')

# read in current and initial time
   iur.seek(0,0)
   i2[:] = numpy.fromfile(iur,int_type,2)
   ntime = i2[0]; ntime0 = i2[1]
   if (in1.nustrt==0):
      print "restarting from ntime, idrun0 = ", ntime, in1.idrun0
   else:
      print "restarting from ntime, = ", ntime
# read in size of electron array
   i2[:] = numpy.fromfile(iur,int_type,2)
   npp = i2[0]; ndimp = i2[1]
   if (ndimp != numpy.size(part,0)):
      print "restart error, idimp=", ndimp, numpy.size(part,0)
      exit(1)
   if (npp != np):
      print "restart warning: new np/old np=", npp, np
      exit(1)
# read in electrons, if non-zero
   if (npp > 0):
      il = idimp*npp
      part[:,:] = numpy.fromfile(iur,float_type,il).reshape(idimp,npp)
# nppmx = maximum number of particles in tile
   if ("nppmx" not in globals()):
      nppmx = numpy.empty((1),int_type)
# kpic = number of electrons in each tile
   if ("kpic" not in globals()):
      kpic = numpy.empty((mx1),int_type,'F')
# find number of electrons in each of mx, tiles: updates kpic, nppmx
   minit1.mdblkp2(part,kpic,nppmx,np,in1.mx,irc)
# allocate vector electron data
   nppmx0 = int((1.0 + in1.xtras)*nppmx)
   ntmax = int(in1.xtras*nppmx)
   npbmx = int(in1.xtras*nppmx)
   if ("ppart" not in globals()):
      ppart = numpy.empty((idimp,nppmx0,mx1),float_type,'F')
   if ("ppbuff" not in globals()):
      ppbuff = numpy.empty((idimp,npbmx,mx1),float_type,'F')
   if ("ncl" not in globals()):
      ncl = numpy.empty((2,mx1),int_type,'F')
   if ("ihole" not in globals()):
      ihole = numpy.empty((2,ntmax+1,mx1),int_type,'F')
# copy ordered electron data for OpenMP: updates ppart and kpic
   mpush1.mpmovin1(part,ppart,kpic,in1.mx,irc)
# sanity check for electrons
   mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)
   if (irc[0] != 0): exit(1)
# read in to determine if ions are moving
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.movion):
      print "movion restart error, movion = ", it, in1.movion
      exit(1)
# ions are moving
   if (i1[0]==1):
# read in size of ion array
      i2[:] = numpy.fromfile(iur,int_type,2)
      npp = i2[0]; ndimp = i2[1]
      if (ndimp != numpy.size(part,0)):
         print "ion restart error, idimp=",ndimp,numpy.size(part,0)
         exit(1)
      if (npp != npi):
         print "restart warning: new npi/old npi=", npp, npi
         npi = npp
# read in ions, if non-zero
      if (npi > 0):
         il = idimp*npp
         part[:,:] = numpy.fromfile(iur,float_type,il).reshape(idimp,npp)
# kipic = number of ions in each tile
      if ("kipic" not in globals()):
         kipic = numpy.empty((mx1),int_type,'F')
# find number of ions in each of mx, tiles: updates kipic, nppmx
      minit1.mdblkp2(part,kipic,nppmx,npi,in1.mx,irc)
# allocate vector ion data
      nppmx1 = int((1.0 + in1.xtras)*nppmx)
      if ("pparti" not in globals()):
         pparti = numpy.empty((idimp,nppmx1,mx1),float_type,'F')
# copy ordered ion data for OpenMP: updates pparti and kipic
      mpush1.mpmovin1(part,pparti,kipic,in1.mx,irc)
# sanity check for ions
      mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)
      if (irc[0] != 0): exit(1)
# ions are not moving, read in ion density
   else:
      i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
      if (it > numpy.size(qi)):
         print "qi restart error, size(qi)=",it,numpy.size(qi)
         exit(1)
      if (it > 0):
         qi[:] = numpy.fromfile(iur,float_type,it)
# read in electric field parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.emf):
      print "warning: emf values differ, emf=",it,in1.emf
   ntime0 += ntime

#-----------------------------------------------------------------------
def dwrite_restart1(iur):
   """
   write out restart diagnostic file for electrostatic code
   input:
   iur = restart file descriptor
   """
   global i1, i2, i3, fname
# write out energy diagnostic parameter
   i1[0] = in1.ntw; i1.tofile(iur)
   if (in1.ntw > 0):
      i1[0] = itw; i1.tofile(iur)
# write out time history array sizes and data
      if (itw > 0):
         it = numpy.size(wt,1)
         i2[0] = numpy.size(wt,0); i2[1] = numpy.size(wt,1)
         i2.tofile(iur)
         wt[0:itw,:].tofile(iur)

# write out electron density diagnostic parameter
   i1[0] = in1.ntde; i1.tofile(iur)
# write out record location
   if (in1.ntde > 0):
      i1[0] = in1.nderec; i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.nderec > 0):
         it = mdiag1.fnrecl(in1.fdename)
         i1[0] = it; i1.tofile(iur)
         if (it > 0):
            fname[:] = ''.join(in1.fdename)
            fname.tofile(iur)

# write out potential diagnostic parameter
   i1[0] = in1.ntp; i1.tofile(iur)
# write out record location
   if (in1.ntp > 0):
      i1[0] = in1.nprec; i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.nprec > 0):
         it = mdiag1.fnrecl(in1.fpname)
         i1[0] = it; i1.tofile(iur)
         if (it > 0):
            fname[:] = ''.join(in1.fpname)
            fname.tofile(iur)
# write out spectrum flag
      if ((in1.ndp==2) or (in1.ndp==3)):
         i1[0] = itp; i1.tofile(iur)
# write out spectrum sizes and data
         if (itp > 0):
            i3[0] = numpy.size(pks,0); i3[1] = numpy.size(pks,1)
            i3[2] = numpy.size(pks,2)
            i3.tofile(iur)
            pks.tofile(iur)
      else:
         it = 0
         i1[0] = it; i1.tofile(iur)

# write out longitudinal efield diagnostic parameter
   i1[0] = in1.ntel; i1.tofile(iur)
# write out record location
   if (in1.ntel > 0):
      i1[0] = in1.nelrec; i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.nelrec > 0):
         it = mdiag1.fnrecl(in1.felname)
         i1[0] = it; i1.tofile(iur)
         if (it > 0):
            fname[:] = ''.join(in1.felname)
            fname.tofile(iur)

# write out ion density diagnostic parameter
   if (in1.movion==1):
      i1[0] = in1.ntdi; i1.tofile(iur)
# write out record location
      if (in1.ntdi > 0):
         i1[0] = in1.ndirec; i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
         if (in1.ndirec > 0):
            it = mdiag1.fnrecl(in1.fdiname)
            i1[0] = it; i1.tofile(iur)
            if (it > 0):
               fname[:] = ''.join(in1.fdiname)
               fname.tofile(iur)
# write out spectrum flag
         if ((in1.nddi==2) or (in1.nddi==3)):
            i1[0] = itdi; i1.tofile(iur)
# write out spectrum sizes and data
            if (itdi > 0):
               i3[0] = numpy.size(pksdi,0); i3[1] = numpy.size(pksdi,1)
               i3[2] = numpy.size(pksdi,2)
               i3.tofile(iur)
               pksdi.tofile(iur)
         else:
            it = 0
            i1[0] = it; i1.tofile(iur)

# write out fluid moments diagnostic parameter
   i1[0] = in1.ntfm; i1.tofile(iur)
   if (in1.ntfm > 0):
# write out electron record location
      i1[0] = in1.nferec; i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.nferec > 0):
         it = mdiag1.fnrecl(in1.ffename)
         i1[0] = it; i1.tofile(iur)
         if (it > 0):
            fname[:] = ''.join(in1.ffename)
            fname.tofile(iur)
      if (in1.movion==1):
# write out ion record location
         i1[0] = in1.nfirec; i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
         if (in1.nfirec > 0):
            it = mdiag1.fnrecl(in1.ffiname)
            i1[0] = it; i1.tofile(iur)
            if (it > 0):
               fname[:] = ''.join(in1.ffiname)
               fname.tofile(iur)

# write out velocity diagnostic parameter
   i1[0] = in1.ntv; i1.tofile(iur)
   if (in1.ntv > 0):
      i1[0] = itv; i1.tofile(iur)
# write out time history array sizes and data
      if (itv > 0):
         i3[0] = numpy.size(fvtm,0); i3[1] = numpy.size(fvtm,1)
         i3[2] = numpy.size(fvtm,2)
         i3.tofile(iur)
         fvtm[0:itv,:].tofile(iur)
         if (in1.movion==1):
            i3[0] = numpy.size(fvtmi,0); i3[1] = numpy.size(fvtmi,1)
            i3[2] = numpy.size(fvtmi,2)
            i3.tofile(iur)
            fvtmi[0:itv,:].tofile(iur)

# write out trajectory diagnostic parameter
   i1[0] = in1.ntt; i1.tofile(iur)
   if (in1.ntt > 0):
      if ((in1.nst==1) or (in1.nst==2)):
         i1[0] = itt; i1.tofile(iur)
# write out time history sizes and data
         if (itt > 0):
            i3[0] = numpy.size(partd,0); i3[1] = numpy.size(partd,1)
            i3[2] = numpy.size(partd,2)
            i3.tofile(iur)
            partd[0:itt,:,:].tofile(iur)
      else:
         it = 0
         i1[0] = it; i1.tofile(iur)

#-----------------------------------------------------------------------
def dread_restart1(iur):
   """
   read in restart diagnostic file for electrostatic code
   input:
   iur = restart file descriptor
   """
   global i1, i2, i3, fname, itw, itp, itv, itt, itdi
# read in energy diagnostic parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.ntw):
      print "restart error: read/expected ntw=", it, in1.ntw
      exit(1)
   if (in1.ntw > 0):
      i1[:] = numpy.fromfile(iur,int_type,1); itw = i1[0]
# read in time history array sizes and data
      if (itw > 0):
         i2[:] = numpy.fromfile(iur,int_type,2)
         iq = i2[0]; it = i2[1]
         if (iq != mtw):
            print "restart error: read/expected mtw=", iq, mtw
            exit(1)
         if (it > numpy.size(wt,1)):
            print "wt size error read/expected=", it, numpy.size(wt,1)
            exit(1)
         il = itw*it
         wt[0:itw,:] = numpy.fromfile(iur,float_type,il).reshape(itw,it)
# restore energy accumulations
         if ("s" in globals()):
            s.fill(0.0)
            s[0] = wt[0,0]
            s[1] = wt[0,1]
            s[2] = wt[0,3]
            s[3] = s[2]
            for it in xrange(1,itw):
               s[0] += wt[it,0]
               s[1] += wt[it,1]
               s[2] = min(s[2],wt[it,3])
               s[3] = max(s[3],wt[it,3])

# read in electron density diagnostic parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.ntde):
      print "restart error: read/expected ntde=", it, in1.ntde
      exit(1)
# read in record location
   if (in1.ntde > 0):
      i1[:] = numpy.fromfile(iur,int_type,1); in1.nderec = i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.nderec > 0):
        i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
        if (it==0):
           print "ntde zero length record error"
           exit(1)
        fname[:] = numpy.fromfile(iur,'S32',1)
        in1.fdename[:] = str(fname[0])

# read in potential diagnostic parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.ntp):
      print "restart error: read/expected ntp=", it, in1.ntp
      exit(1)
# read in record location
   if (in1.ntp > 0):
      i1[:] = numpy.fromfile(iur,int_type,1); in1.nprec = i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.nprec > 0):
         i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
         if (it==0):
            print "ntp zero length record error"
            exit(1)
         fname[:] = numpy.fromfile(iur,'S32',1)
         in1.fpname[:] = str(fname[0])
# read in spectrum flag
      i1[:] = numpy.fromfile(iur,int_type,1); itp = i1[0]
# read in spectrum sizes and data
      if (itp > 0):
         i3[:] = numpy.fromfile(iur,int_type,3)
         ir = i3[0]; it = i3[1]; iq = i3[2]
         if (ir != 4):
            print "pks size error: read/expected 4 =", ir
            exit(1)
         if (it != in1.modesxp):
            print ("pks size error: read/expected modesxp=", it,
                    in1.modesxp)
            exit(1)
         if (iq != iw):
            print "pks size error: read/expected iw=", iq, iw
            exit(1)
         il = ir*it*iq
         pks[:,:,:] = numpy.fromfile(iur,double_type,il).reshape(4,it,iq)

# read in longitudinal efield diagnostic parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.ntel):
      print "restart error: read/expected ntel=", it, in1.ntel
      exit(1)
# read in record location
   if (in1.ntel > 0):
      i1[:] = numpy.fromfile(iur,int_type,1); in1.nelrec = i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.nelrec > 0):
        i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
        if (it==0):
           print "ntel zero length record error"
           exit(1)
        fname[:] = numpy.fromfile(iur,'S32',1)
        in1.felname[:] = str(fname[0])

# read in ion density diagnostic parameter
   if (in1.movion==1):
      i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
      if (it != in1.ntdi):
         print "restart error: read/expected ntdi=", it, in1.ntdi
         exit(1)
# read in record location
      if (in1.ntdi > 0):
         i1[:] = numpy.fromfile(iur,int_type,1); in1.ndirec = i1[0]
# read in record length (zero if error) and file name (if no error)
         if (in1.ndirec > 0):
            i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
            if (it==0):
               print "ntdi zero length record error"
               exit(1)
            fname[:] = numpy.fromfile(iur,'S32',1)
            in1.fdiname[:] = str(fname[0])
# read in spectrum flag
         i1[:] = numpy.fromfile(iur,int_type,1); itdi = i1[0]
# read in spectrum sizes and data
         if (itdi > 0):
            i3[:] = numpy.fromfile(iur,int_type,3)
            ir = i3[0]; it = i3[1]; iq = i3[2]
            if (ir != 4):
               print "pksdi size error: read/expected 4 =", ir
               exit(1)
            if (it != in1.modesxdi):
               print ("pksdi size error: read/expected modesxdi=", it,
                       in1.modesxdi)
               exit(1)
            if (iq != iwi):
               print "pksdi size error: read/expected iwi=",iq,iwi
               exit(1)
            il = ir*it*iq
            pksdi[:,:,:] = (numpy.fromfile(iur,double_type,il).
                            reshape(4,it,iq))

# read in fluid moments diagnostic parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.ntfm):
      print "restart error: read/expected ntfm=", it, in1.ntfm
      exit(1)
# read in electron data
   if (in1.ntfm > 0):
# read in electron record location
      i1[:] = numpy.fromfile(iur,int_type,1); in1.nferec = i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.nferec > 0):
         i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
         if (it==0):
            print "ntfm zero length record error"
            exit(1)
         fname[:] = numpy.fromfile(iur,'S32',1)
         in1.ffename[:] = str(fname[0])
# read in ion data
      if (in1.movion==1):
# read in ion record location
         i1[:] = numpy.fromfile(iur,int_type,1); in1.nfirec = i1[0]
# read in ion record length (zero if error) and file name (if no error)
         if (in1.nfirec > 0):
            i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
            if (it==0):
               print "ntfm zero length ion record error"
               exit(1)
            fname[:] = numpy.fromfile(iur,'S32',1)
            in1.ffiname[:] = str(fname[0])

# read in velocity diagnostic parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.ntv):
      print "restart error: read/expected ntv=", it, in1.ntv
      exit(1)
   if (in1.ntv > 0):
      i1[:] = numpy.fromfile(iur,int_type,1); itv = i1[0]
# read in time history array sizes and data
      if (itv > 0):
         i3[:] = numpy.fromfile(iur,int_type,3)
         iq = i3[0]; it = i3[1]; ir = i3[2]
         if (iq != mtv):
            print "restart error: read/expected mtv=", iq, mtv
            exit(1)
         if (it != in1.ndim):
            print "fvtm size error read/expected ndim=", it, in1.ndim
            exit(1)
         if (ir != 3):
            print "fvtm size error read/expected 3=", ir
            exit(1)
         il = itv*it*ir
         fvtm[0:itv,:,:] = (numpy.fromfile(iur,float_type,il).
                            reshape(itv,it,3))
         if (in1.movion==1):
            i3[:] = numpy.fromfile(iur,int_type,3)
            iq = i3[0]; it = i3[1]; ir = i3[2]
            if (iq != mtv):
               print "ion restart error: read/expected mtv=", iq, mtv
               exit(1)
            if (it != in1.ndim):
               print "fvtmi size error read/expected ndim=",it, in1.ndim
               exit(1)
            if (ir != 3):
               print "fvtmi size error read/expected 3=", ir
               exit(1)
            il = itv*it*ir
            fvtmi[0:itv,:,:] = (numpy.fromfile(iur,float_type,il).
                                reshape(itv,it,3))

# read in trajectory diagnostic parameter
   i1[:] = numpy.fromfile(iur,int_type,1); it = i1[0]
   if (it != in1.ntt):
      print "restart error: read/expected ntt=", it, in1.ntt
      exit(1)
   if (in1.ntt > 0):
      i1[:] = numpy.fromfile(iur,int_type,1); itt = i1[0]
# read in time history sizes and data
      if (itt > 0):
         if ((in1.nst==1) or (in1.nst==2)):
            i3[:] = numpy.fromfile(iur,int_type,3)
            ir = i3[0]; it = i3[1]; iq = i3[2]
            if (ir != mtt):
               print "restart error: read/expected mtt=", ir, mtt
               exit(1)
            if (it != idimp):
               print ("partd size error read/expected idimp=", it,
                       numpy.size(partd,1))
               exit(1)
            if (iq != in1.nprobt):
               print ("partd size error read/expected nprobt=", iq,
                       in1.nprobt)
               exit(1)
            il = itt*it*iq
            partd[0:itt,:,:] = (numpy.fromfile(iur,float_type,il).
                                reshape(itt,it,iq))

#-----------------------------------------------------------------------
def close_restart1():
   """ close reset and restart files """
   global iur, iurr, iur0, i1, i2, i3, fname
# iur, iurr = restart, reset, old restart file descriptors
   iurr.close()
   if (in1.nustrt==1):
      if (in1.ntr > 0):
         iur.close()
   elif (in1.nustrt==2):
      iur.close()
   elif (in1.nustrt==0):
      if (in1.ntr > 0):
         iur.close()
      if (in1.idrun != in1.idrun0):
         iur0.close()
# deallocate scratch numpy arrays
   del i1, i2, i3, fname
