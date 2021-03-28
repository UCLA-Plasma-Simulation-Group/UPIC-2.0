#-----------------------------------------------------------------------
from __future__ import print_function
"""
High Level library for 1-2/2D Electromagnetic OpenMP PIC code

functions defined:

init_fields13: allocate electromagnetic field data for standard code
del_fields13: delete electromagnetic field data for standard code

init_detfield13: calculate initial darwin electric field with
                 approximate value
                 
init_electrons13: initialize electrons for 1-2/2d code
deposit_ecurrent13: deposit electron current with OpenMP
push_electrons13: push electrons with OpenMP

init_ions13: initialize ions for 1-2/2d code
deposit_icurrent13: deposit ion current with OpenMP
push_ions13: push ions with OpenMP

em_time_reverse1: start running simulation backwards

init_energy_diag13: initialize energy diagnostic
energy_diag13: energy diagnostic
print_energy13: print energy summaries
del_energy_diag13: delete energy diagnostic

init_spectrum13: allocate scratch arrays for vector fields
del_spectrum13: delete scratch arrays for vector fields

init_ecurrent_diag13: initialize electron current density diagnostic
ecurrent_diag13: electron current density diagnostic
del_ecurrent_diag13: delete electron current density diagnostic

init_icurrent_diag13: initialize ion current density diagnostic
icurrent_diag13: ion current density diagnostic
del_icurrent_diag13: delete ion current density diagnostic

init_vrpotential_diag13: initialize radiative vector potential
                         diagnostic
vrpotential_diag13: radiative vector potential diagnostic
del_vrpotential_diag13: delete radiative vector potential diagnostic

init_vpotential_diag13: initialize vector potential diagnostic
vpotential_diag13: vector potential diagnostic
del_vpotential_diag13: delete vector potential diagnostic

init_etfield_diag13: initialize transverse efield diagnostic
etfield_diag13: transverse efield diagnostic
del_etfield_diag13: delete transverse efield diagnostic

init_bfield_diag13: initialize magnetic field diagnostic
bfield_diag13: magnetic field diagnostic
del_bfield_diag13: delete magnetic field diagnostic

init_efluidms_diag13: initialize electron fluid moments diagnostic
efluidms_diag13: electron fluid moments diagnostic

init_ifluidms_diag13: initialize ion fluid moments diagnostic
ifluidms_diag13: ion fluid moments diagnostic

init_evelocity_diag13: initialize electron velocity diagnostic
evelocity_diag13: electron velocity diagnostic

init_ivelocity_diag13: initialize ion velocity diagnostic
ivelocity_diag13: ion velocity diagnostic

init_traj_diag13: initialize trajectory diagnostic
traj_diag13: trajectory diagnostic
del_traj_diag13: delete trajectory diagnostic

init_ephasesp_diag13: initialize electron phase space diagnostic
init_iphasesp_diag13: initialize ion phase space diagnostic

print_timings13: print timing summaries

reset_diags13: reset electrostatic/electromagnetic diagnostics
close_diags13: close diagnostics

initialize_diagnostics13: initialize all diagnostics from namelist
                          input parameters
                          
bwrite_restart13: write out basic restart file for electromagnetic
                  code
bread_restart13: read in basic restart file for electromagnetic code
dwrite_restart13: write out restart diagnostic file for
                  electromagnetic code
dread_restart13: read in restart diagnostic file for electromagnetic
                 code

written by Viktor K. Decyk, UCLA
copyright 1999-2016, regents of the university of california
update: february 2, 2021
"""
import sys
import math
import numpy

# sys.path.append('./mbeps1.source')
from libmpush1 import *
from libmbpush1 import *
import s1

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# idimp = number of particle coordinates = 4
# ipbc = particle boundary condition: 1 = periodic
idimp = 4; ipbc = 1
s1.idimp = idimp

# wf/wb = magnetic field/transverse electric field
wf = numpy.zeros((1),float_type)
wb = numpy.zeros((1),float_type)
zero = 0.0

# declare scalars for standard code
npi = 0
ws = numpy.zeros((1),float_type)

# declare scalars for OpenMP code
irc = numpy.zeros((1),int_type)
irc2 = numpy.zeros((2),int_type)

# declare scalars for diagnostics
itje = 0; ita = 0; itet = 0; itar = 0
itji = 0
# default Fortran unit numbers
iuje = 21; iua = 13; iuet = 14; iub = 15; iuar = 16
iuji = 22

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
dtime = numpy.empty((1),double_type)
tguard = numpy.zeros((1),float_type)
tfield = numpy.zeros((1),float_type)
tdjpost = numpy.zeros((1),float_type)
tpush = numpy.zeros((1),float_type)
tsort = numpy.zeros((1),float_type)
tdiag = numpy.zeros((1),float_type)

# initialize scalars for standard code
# increase number of coordinates for particle tag
if ((in1.ntt > 0) or ((in1.nts > 0) and (in1.ntsc > 0))):
   idimp += 1
   s1.idimp = idimp
# np = total number of electrons in simulation
np = s1.np
# nx = number of grid points in x direction
nx = s1.nx; nxh = int(nx/2)
# npi = total number of ions in simulation
if (in1.movion > 0):
   npi = s1.npi
nxe = nx + 2; nxeh = int(nxe/2)
# mx1 = number of tiles in x direction
mx1 = s1.mx1
# nloop = number of time steps in simulation
nloop = s1.nloop
qbme = s1.qbme
# affp = float(nx)/float(np)
affp = s1.affp
if (in1.movion==1):
   qbmi = s1.qbmi
   vtxi = s1.vtxi
   vtyi = in1.vty/numpy.sqrt(in1.rmass*in1.rtempyi)
   vtzi = in1.vtz/numpy.sqrt(in1.rmass*in1.rtempzi)
   vtdxi = s1.vtdxi
   vtdyi = in1.vtdy/numpy.sqrt(in1.rmass*in1.rtempdyi)
   vtdzi = in1.vtdz/numpy.sqrt(in1.rmass*in1.rtempdzi)
dth = 0.0

# check for unimplemented features
plist = s1.plist

#-----------------------------------------------------------------------
def init_fields13():
   """ allocate electromagnetic field data for standard code """
   global cue, fxyze, byze, eyz, byz
# allocate electrostatic field data: qe, qi, fxe, ffc, mixup, sct
   s1.init_fields1()
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

#-----------------------------------------------------------------------
def del_fields13():
   """ delete electromagnetic field data for standard code """
   global cue, fxyze, byze, eyz, byz
   del cue, fxyze, byze, eyz, byz
   s1.del_fields1()

#-----------------------------------------------------------------------
def init_detfield13():
   """
   calculate initial darwin electric field with approximate value
   approximation assumes initial accelerations are zero
   """
# calculate initial darwin electric field
   amu = numpy.zeros((2,nxe),float_type,'F')
   dcu = numpy.empty((2,nxe),float_type,'F')
   mcurd1.wmgmjpost1(s1.ppart,amu,s1.kpic,in1.qme,in1.ci,tdjpost,in1.mx,
                     in1.relativity)
   mgard1.macguard1(amu,tguard,nx)
   isign = -1
   mfft1.mfft1rn(amu,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)
   mfield1.mdcuperp1(dcu,amu,tfield,nx)
   del amu
   mfield1.metfield1(dcu,eyz,s1.ffc,in1.ci,wf,tfield,nx)
   del dcu

#-----------------------------------------------------------------------
def init_electrons13():
   """ initialize electrons for 1-2/2d code """
# part = particle array
   s1.part = numpy.empty((idimp,max(np,npi)),float_type,'F')
# background electrons
   if (in1.npx > 0):
# calculates initial particle co-ordinates with various density profiles
      minit1.mfdistr1(s1.part,in1.ampdx,in1.scaledx,in1.shiftdx,1,
                      in1.npx,nx,ipbc,in1.ndprof,irc)
      if (irc[0] != 0): exit(1)
# initialize particle velocities
# special cases
      if (in1.nvdist==3):
         minit1.mvbdistr1h(s1.part,1,in1.vtx,in1.vtz,in1.vx0,in1.vz0,
                           in1.omx,in1.omy,in1.omz,in1.npx,irc)
# general cases
      else:
         minit1.wmvdistr1h(s1.part,1,in1.vtx,in1.vty,in1.vtz,in1.vx0,
                           in1.vy0,in1.vz0,in1.ci,in1.npx,in1.nvdist,
                           in1.relativity,irc)
      if (irc[0] != 0): exit(1)
# beam electrons
   if (in1.npxb > 0):
      it = in1.npx + 1
# calculates initial particle co-ordinates with various density profiles
      minit1.mfdistr1(s1.part,in1.ampdx,in1.scaledx,in1.shiftdx,it,
                      in1.npxb,nx,ipbc,in1.ndprof,irc)
      if (irc[0] != 0): exit(1)
# initialize particle velocities
# special cases
      if (in1.nvdist==3):
         minit1.mvbdistr1h(s1.part,it,in1.vtdx,in1.vtdz,in1.vdx,in1.vdz,
                           in1.omx,in1.omy,in1.omz,in1.npxb,irc)
# general cases
      else:
         minit1.wmvdistr1h(s1.part,it,in1.vtdx,in1.vtdy,in1.vtdz,
                           in1.vdx,in1.vdy,in1.vdz,in1.ci,in1.npxb,
                           in1.nvdist,in1.relativity,irc)
      if (irc[0] != 0): exit(1)

# mark electron beam particles
   if ((in1.nts > 0) and (in1.ntsc > 0)):
      mdiag1.setmbeam1(s1.part,in1.npx,irc)
      if (irc[0] != 0): exit(1)

# nppmx = rmaximum number of particles in tile
   s1.nppmx = numpy.empty((1),int_type)
# kpic = number of electrons in each tile
   s1.kpic = numpy.empty((mx1),int_type,'F')

# find number of electrons in each of mx, tiles: updates kpic, nppmx
   minit1.mdblkp2(s1.part,s1.kpic,s1.nppmx,np,in1.mx,irc)
   if (irc[0] != 0): exit(1)

# allocate vector electron data
   nppmx0 = int((1.0 + in1.xtras)*s1.nppmx)
   ntmax = int(in1.xtras*s1.nppmx)
   npbmx = int(in1.xtras*s1.nppmx)
# ppart = tiled electron array
   s1.ppart = numpy.empty((idimp,nppmx0,mx1),float_type,'F')
# ppbuff = buffer array for reordering tiled particle array
   s1.ppbuff = numpy.empty((idimp,npbmx,mx1),float_type,'F')
# ncl = number of particles departing tile in each direction
   s1.ncl = numpy.empty((2,mx1),int_type,'F')
# ihole = location/destination of each particle departing tile
   s1.ihole = numpy.empty((2,ntmax+1,mx1),int_type,'F')

# copy ordered electron data for OpenMP: updates ppart and kpic
   mpush1.mpmovin1(s1.part,s1.ppart,s1.kpic,in1.mx,irc)
   if (irc[0] != 0): exit(1)

# sanity check for electrons
   mpush1.mcheck1(s1.ppart,s1.kpic,nx,in1.mx,irc)
   if (irc[0] != 0): exit(1)

#-----------------------------------------------------------------------
def deposit_ecurrent13(ppart,kpic):
   """
   deposit electron current with OpenMP
   input/output:
   ppart = tiled electron particle array
   kpic = number of electrons in each tile
   """
# updates ppart and cue, and possibly ncl, ihole, irc
   mcurd1.wmdjpost1(ppart,cue,kpic,s1.ncl,s1.ihole,in1.qme,dth,in1.ci,
                    tdjpost,nx,in1.mx,ipbc,in1.relativity,plist,irc)
# add guard cells: updates cue
   mgard1.macguard1(cue,tguard,nx)

# reorder electrons by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   if (irc[0]==0):
      msort1.wmporder1(ppart,s1.ppbuff,kpic,s1.ncl,s1.ihole,tsort,nx,
                       in1.mx,plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc; irc[0] = 0

# sanity check for electrons
   if (irc2[0]==0):
      if (in1.monitor > 0):
         mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)
         if (irc[0] != 0): exit(1)
# recover from wmporder1 errors: updates ppart
   elif (irc2[0] != 0):
      s1.reorder_electrons1(irc2)

#-----------------------------------------------------------------------
def push_electrons13(ppart,kpic):
   """
   push electrons with OpenMP
   input/output:
   ppart = tiled electron particle array
   kpic = number of electrons in each tile
   """
   s1.wke[0] = 0.0
# updates ppart and wke, and possibly ncl, ihole, irc
# Boris pusher
   if (in1.mzf==0):
      mbpush1.wmbpush1(ppart,fxyze,byze,kpic,s1.ncl,s1.ihole,in1.omx,
                       qbme,in1.dt,dth,in1.ci,s1.wke,tpush,nx,in1.mx,
                       ipbc,in1.relativity,plist,irc)
# Analytic Boris pusher
   elif (in1.mzf==2):
      mbpush1.wmabpush1(ppart,fxyze,byze,kpic,s1.ncl,s1.ihole,in1.omx,
                        qbme,in1.dt,dth,in1.ci,s1.wke,tpush,nx,in1.mx,
                        ipbc,in1.relativity,plist,irc)
# Exact Analytic pusher
   elif (mzf==3):
      mbpush1.wmeabpush1(ppart,fxyze,byze,kpic,s1.ncl,s1.ihole,in1.omx,
                         qbme,in1.dt,dth,in1.ci,s1.wke,tpush,nx,in1.mx,
                         ipbc,in1.relativity,plist,irc)
# zero force: updates ppart, wke and possibly ncl, ihole, and irc
   elif (in1.mzf==1):
      mbpush1.wmpush1zf(ppart,kpic,s1.ncl,s1.ihole,in1.dth,in1.ci,
                        s1.wke,tpush,in1.nx,in1.mx,ipbc,in1.relativity,
                        plist,irc)

# reorder electrons by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   if (irc[0]==0):
      msort1.wmporder1(ppart,s1.ppbuff,kpic,s1.ncl,s1.ihole,tsort,nx,
                       in1.mx,plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc; irc[0] = 0

# sanity check for electrons
   if (irc2[0]==0):
      if (in1.monitor > 0):
         mpush1.mcheck1(ppart,kpic,nx,in1.mx,irc)
         if (irc[0] != 0): exit(1)
# recover from wmporder1 errors: updates ppart
   elif (irc2[0] != 0):
      s1.reorder_electrons1(irc2)

#-----------------------------------------------------------------------
def init_ions13():
   """ initialize ions """
   global cui
# part = particle array
# kipic = number of ions in each tile
   s1.kipic = numpy.empty((mx1),int_type,'F')
# cui = ion current density with guard cells
   cui = numpy.zeros((2,nxe),float_type,'F')
# create dummy array to avoid undefined arguments later
   if (in1.movion==0):
      cui = numpy.zeros((1,1),float_type,'F')
   it = in1.npxi + 1
# background ions
   if (in1.npxi > 0):
# calculates initial particle co-ordinates with various density profiles
      minit1.mfdistr1(s1.part,in1.ampdxi,in1.scaledxi,in1.shiftdxi,1,
                      in1.npxi,nx,ipbc,in1.ndprofi,irc)
# initialize particle velocities
# special cases
      if (in1.nvdist==3):
         minit1.mvbdistr1h(s1.part,1,vtxi,vtzi,in1.vxi0,in1.vzi0,
                           in1.omx,in1.omy,in1.omz,in1.npxi,irc)
# general cases
      else:
         minit1.wmvdistr1h(s1.part,1,vtxi,vtyi,vtzi,in1.vxi0,in1.vyi0,
                           in1.vzi0,in1.ci,in1.npxi,in1.nvdist,
                           in1.relativity,irc)
# beam ions
   if (in1.npxbi > 0):
# calculates initial particle co-ordinates with various density profiles
      minit1.mfdistr1(s1.part,in1.ampdxi,in1.scaledxi,in1.shiftdxi,it,
                      in1.npxbi,nx,ipbc,in1.ndprofi,irc)
# initialize particle velocities
# special cases
      if (in1.nvdist==3):
         minit1.mvbdistr1h(s1.part,it,vtdxi,vtdzi,in1.vdxi,in1.vdzi,
                           in1.omx,in1.omy,in1.omz,in1.npxbi,irc)
# general cases
      else:
         minit1.wmvdistr1h(s1.part,it,vtdxi,vtdyi,vtdzi,in1.vdxi,
                           in1.vdyi,in1.vdzi,in1.ci,in1.npxbi,
                           in1.nvdist,in1.relativity,irc)

# mark ion beam particles
   if ((in1.nts > 0) and (in1.ntsc > 0)):
      mdiag1.setmbeam1(s1.part,in1.npxi,irc)
      if (irc[0] != 0): exit(1)

# find number of ions in each of mx, tiles: updates kipic, nppmx
   minit1.mdblkp2(s1.part,s1.kipic,s1.nppmx,npi,in1.mx,irc)
   if (irc[0] != 0): exit(1)

# allocate vector ion data
   nppmx1 = int((1.0 + in1.xtras)*s1.nppmx)
   s1.pparti = numpy.empty((idimp,nppmx1,mx1),float_type,'F')
# copy ordered ion data for OpenMP: updates pparti and kipic
   mpush1.mpmovin1(s1.part,s1.pparti,s1.kipic,in1.mx,irc)
   if (irc[0] != 0): exit(1)

# sanity check for ions
   mpush1.mcheck1(s1.pparti,s1.kipic,nx,in1.mx,irc)
   if (irc[0] != 0): exit(1)

#-----------------------------------------------------------------------
def deposit_icurrent13(pparti,kipic):
   """
   deposit ion current with OpenMP
   input/output:
   pparti = tiled electron/ion particle arrays
   kipic = number of electrons/ions in each tile
   """
# updates pparti and cui, and possibly ncl, ihole, irc
   mcurd1.wmdjpost1(pparti,cui,kipic,s1.ncl,s1.ihole,in1.qmi,dth,in1.ci,
                    tdjpost,nx,in1.mx,ipbc,in1.relativity,plist,irc)
# add guard cells: updates cui
   mgard1.macguard1(cui,tguard,nx)

# reorder ions by tile with OpenMP:
# updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
   if (irc[0]==0):
      msort1.wmporder1(pparti,s1.ppbuff,kipic,s1.ncl,s1.ihole,tsort,nx,
                       in1.mx,plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc; irc[0] = 0

# sanity check for ions
   if (irc2[0]==0):
      if (in1.monitor > 0):
         mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)
         if (irc[0] != 0): exit(1)
# recover from wmporder1 errors: updates pparti
   elif (irc2[0] != 0):
      s1.reorder_ions1(irc2)

#-----------------------------------------------------------------------
def push_ions13(pparti,kipic):
   """
   push ions with OpenMP
   input/output:
   pparti = tiled electron/ion particle arrays
   kipic = number of electrons/ions in each tile
   """
   s1.wki[0] = 0.0
# updates pparti and wki, and possibly ncl, ihole, irc
# Boris pusher
   if (in1.mzf==0):
      mbpush1.wmbpush1(pparti,fxyze,byze,kipic,s1.ncl,s1.ihole,in1.omx,
                       qbmi,in1.dt,dth,in1.ci,s1.wki,tpush,nx,in1.mx,
                       ipbc,in1.relativity,plist,irc)
# Analytic Boris pusher
   elif (in1.mzf==2):
      mbpush1.wmabpush1(pparti,fxyze,byze,kipic,s1.ncl,s1.ihole,in1.omx,
                        qbmi,in1.dt,dth,in1.ci,s1.wki,tpush,nx,in1.mx,
                        ipbc,in1.relativity,plist,irc)
# Exact Analytic pusher
   elif (in1.mzf==3):
      mbpush1.wmeabpush1(pparti,fxyze,byze,kipic,s1.ncl,s1.ihole,
                         in1.omx,qbmi,in1.dt,dth,in1.ci,s1.wki,tpush,nx,
                         in1.mx,ipbc,in1.relativity,plist,irc)
# zero force: updates pparti, wki and possibly ncl, ihole, and irc
   elif (in1.mzf==1):
      mpush1.wmpush1zf(pparti,kipic,s1.ncl,s1.ihole,in1.dth,in1.ci,
                       s1.wki,tpush,nx,in1.mx,ipbc,in1.relativity,plist,
                       irc)
   s1.wki[0] *= in1.rmass

# reorder ions by tile with OpenMP:
# updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
   if (irc[0]==0):
      msort1.wmporder1(pparti,s1.ppbuff,kipic,s1.ncl,s1.ihole,tsort,nx,
                       in1.mx,plist,irc2)
   else:
      irc2[0] = 1; irc2[1] = irc; irc[0] = 0

# sanity check for ions
   if (irc2[0]==0):
      if (in1.monitor > 0):
         mpush1.mcheck1(pparti,kipic,nx,in1.mx,irc)
         if (irc[0] != 0): exit(1)
# recover from wmporder1 errors: updates pparti
   elif (irc2[0] != 0):
      s1.reorder_ions1(irc2)

#-----------------------------------------------------------------------
def del_ions13():
   """ delete ions """
   global cui
   s1.del_ions1()
   del cui

#-----------------------------------------------------------------------
def em_time_reverse1():
   """
   start running simulation backwards
   need to advance maxwell field solver one step ahead
   """
   global dth
# deposit electron current: updates cue
   cue.fill(0.0)
   mcurd1.wmdjpost1(s1.ppart,cue,s1.kpic,s1.ncl,s1.ihole,in1.qme,zero,
                    in1.ci,tdjpost,nx,in1.mx,ipbc,in1.relativity,plist,
                    irc)
   mgard1.macguard1(cue,tguard,nx)
# deposit ion current: updates cui
   if (in1.movion==1):
      cui.fill(0.0)
      mcurd1.wmdjpost1(s1.pparti,cui,s1.kipic,s1.ncl,s1.ihole,in1.qmi,
                       zero,in1.ci,tdjpost,nx,in1.mx,ipbc,
                       in1.relativity,plist,irc)
      mgard1.macguard1(cui,tguard,nx)
   isign = -1
   mfft1.mfft1rn(cue,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)
# updates eyz, byz, wef, ws
   mfield1.mmaxwel1(eyz,byz,cue,s1.ffc,in1.ci,in1.dt,wef,ws,tfield,nx)
# reverse time
   in1.dt = -in1.dt; dth = -dth

#-----------------------------------------------------------------------
def init_energy_diag13():
   """ initialize energy diagnostic """
# wt = energy time history array
   global s, wef
   s1.mtw = int((nloop - 1)/in1.ntw + 1); s1.itw = 0
   s1.wt = numpy.zeros((s1.mtw,7),float_type,'F')
   s = numpy.zeros((7),double_type,'F')
   wef = numpy.zeros((1),float_type)

#-----------------------------------------------------------------------
def energy_diag13(wt,ntime,iuot):
   """
   energy diagnostic
   input/output:
   wt = energy time history array
   input:
   ntime = current time step
   iuot = output file descriptor
   """
   global ws, s, wef
   wef[0] = s1.we[0] + wf[0] + wb[0]
   ws[0] = wef[0] + s1.wke[0] + s1.wki[0]
   if (ntime==0):
      s[5] = ws[0]
   if (in1.ndw > 0):
      print ("Total Field, Kinetic and Total Energies:",file=iuot)
      if (in1.movion==0):
         iuot.write("%14.7e %14.7e %14.7e\n" % (wef[0],s1.wke[0],ws[0])) 
      else:
         iuot.write("%14.7e %14.7e %14.7e %14.7e\n" % (wef[0],s1.wke[0],
                    s1.wki[0],ws[0])) 
      print ("Electric(l,t) and Magnetic Energies = ",file=iuot)
      iuot.write("%14.7e %14.7e %14.7e\n" % (s1.we[0],wf[0],wb[0]))
   wt[s1.itw,:] = ([wef[0],s1.wke[0],s1.wki[0],ws[0],s1.we[0],wf[0],
                    wb[0]])
   s1.itw += 1
   s[0] += s1.we[0]
   s[1] += s1.wke[0]
   s[2] += wf[0]
   s[3] += wb[0]
   s[4] += s1.wki[0]
   s[5] = min(s[5],float(ws[0]))
   s[6] = max(s[6],float(ws[0]))

#-----------------------------------------------------------------------
def print_energy13(wt,iuot):
   """
   print energy summaries
   input:
   wt = energy time history array
   iuot = output file descriptor
   """
   global s
   swe = s[0]; swke = s[1]; swf = s[2]; swb = s[3]
   s[5] = (s[6] - s[5])/wt[0,3]
   print ("Energy Conservation = ",float(s[5]),file=iuot)
   swe = swe/float(s1.itw)
   print ("Average Field Energy <WE> = ",float(swe),file=iuot)
   swke = swke/float(s1.itw)
   print ("Average Electron Kinetic Energy <WKE> = ",float(swke),
          file=iuot)
   print ("Ratio <WE>/<WKE>= ",float(swe/swke),file=iuot)
   swf = swf/float(s1.itw)
   print ("Average Transverse EField Energy <WF> = ",float(swf),
          file=iuot)
   print ("Ratio <WF>/<WKE>= ",float(swf/swke),file=iuot)
   swb = swb/float(s1.itw)
   print ("Average Magnetic Field Energy <WB> = ",float(swb),file=iuot)
   print ("Ratio <WB>/<WKE>= ",float(swb/swke),file=iuot)
   print (file=iuot)

#-----------------------------------------------------------------------
def del_energy_diag13():
   """ delete energy diagnostic """
# wt = energy time history array
   global s, wef
   del s1.wt, s, wef

#-----------------------------------------------------------------------
def init_spectrum13():
   """ allocate scratch arrays for vector fields """
   global vfield, vfieldc
   vfield = numpy.empty((2,nxe),float_type,'F')
   vfieldc = numpy.empty((2,nxh),complex_type,'F')
# allocate and initialize high frequency array for spectral analysis
   if ((in1.nta > 0) or (in1.ntet>0) or (in1.ntar > 0)):
      global iwr, wmr
      iwr = int((in1.wrmax - in1.wrmin)/in1.dwr + 1.5)
      wmr = numpy.empty((iwr),float_type,'F')
      wmr[:] = in1.wrmin + in1.dwr*numpy.linspace(0,iwr-1,iwr)
# allocate and initialize frequency array for ion spectral analysis
   if (in1.movion==1):
      if (in1.ntji > 0):
         if ("wmi" not in globals()):
            global iwi, wmi
            iwi = int((in1.wimax - in1.wimin)/in1.dwi + 1.5)
            wmi = numpy.empty((iwi),float_type,'F')
            wmi[:] = in1.wimin + in1.dwi*numpy.linspace(0,iwi-1,iwi)

#-----------------------------------------------------------------------
def del_spectrum13():
   """ delete scratch arrays for vector fields """
   global vfield, vfieldc
   del vfield, vfieldc
   if ((in1.nta>0) or (in1.ntet>0) or (in1.ntar>0)):
      global wmr
      del wmr
   if (in1.ntji > 0):
      global wmi
      if ("wmi" in globals()):
         del wmi

#-----------------------------------------------------------------------
def init_ecurrent_diag13():
   """ initialize electron current density diagnostic """
   global curet, iuje
   fjename = "curek1." + s1.cdrun
# write filename to diagnostic metafile
   if (sys.version_info.major==3):
      in1.fjename = fjename.ljust(in1.fjename.dtype.itemsize)
   else:
      in1.fjename[:] = fjename
# open file
   in1.modesxje = int(min(in1.modesxje,nxh+1))
# curet = store selected fourier modes for electron current density
   curet = numpy.empty((2,in1.modesxje),complex_type,'F')
# open file: updates njerec and possibly iuje
   if (in1.njerec==0):
      mdiag1.dafopenvc1(curet,iuje,in1.njerec,fjename)

#-----------------------------------------------------------------------
def ecurrent_diag13(vfield):
   """
   electron current density diagnostic
   input/output:
   vfield = scratch array for vector field
   """
   vfield[:] = numpy.copy(cue)
# transform electron current density to fourier space: updates vfield
   isign = -1
   mfft1.mfft1rn(vfield,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)
# calculate smoothed electron current in fourier space: updates vfieldc
   mfield1.msmooth13(vfield,vfieldc,s1.ffc,tfield,nx)
# store selected fourier modes: updates curet
   mfield1.mrdvmodes1(vfieldc,curet,tfield,nx,in1.modesxje)
# write diagnostic output: updates njerec
   mdiag1.dafwritevc1(curet,tdiag,iuje,in1.njerec,in1.modesxje)
# transform smoothed electron current to real space: updates vfield
   mfft1.mfft1crn(vfieldc,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
   mgard1.mcguard1(vfield,tguard,nx)

#-----------------------------------------------------------------------
def del_ecurrent_diag13():
   """ delete electron current density diagnostic """
   global curet
   if (in1.njerec > 0):
      in1.closeff(iuje)
      in1.njerec -= 1
   del curet

#-----------------------------------------------------------------------
def init_icurrent_diag13():
   """ initialize ion current density diagnostic """
   global curit, iuji
   fjiname = "curik1." + s1.cdrun
# write filename to diagnostic metafile
   if (sys.version_info.major==3):
      in1.fjiname = fjiname.ljust(in1.fjiname.dtype.itemsize)
   else:
      in1.fjiname[:] = fjiname
# open file
   in1.modesxji = int(min(in1.modesxji,nxh+1))
# curit = store selected fourier modes for ion current density
   curit = numpy.empty((2,in1.modesxji),complex_type,'F')
# open file: updates njirec and possibly iuji
   if (in1.njirec==0):
      mdiag1.dafopenvc1(curit,iuji,in1.njirec,fjiname)
# ion current spectral analysis
   global mtji, itji, vpkwji, vpksji, vwkji
   if ((in1.ndji==2) or (in1.ndji==3)):
      mtji = int((nloop - 1)/in1.ntji) + 1; itji = 0
# vpkwji = power spectrum for ion current density
      vpkwji = numpy.empty((2,in1.modesxji,s1.iwi,2),float_type,'F')
# vpksji = accumulated complex spectrum for ion current density
      vpksji = numpy.zeros((2,4,in1.modesxji,s1.iwi),double_type,'F')
# vwkji = maximum frequency as a function of k for ion current
      vwkji = numpy.empty((2,in1.modesxji,2),float_type,'F')
# create dummy arrays to avoid undefined arguments later
   else:
      vpkwji = numpy.zeros((1,1,1,1),float_type,'F')
      vwkji = numpy.zeros((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def icurrent_diag13(vfield,vpkwdi,vwkdi,ntime):
   """
   ion current density diagnostic
   input/output:
   vfield = scratch array for vector field
   vpkwji = power spectrum for ion current density
   vwkji = maximum frequency as a function of k for ion current
   input:
   ntime = current time step
   """
   vfield[:] = numpy.copy(cui)
# transform ion current density to fourier space: updates vfield
   isign = -1
   mfft1.mfft1rn(vfield,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)
# calculate smoothed ion current in fourier space: updates vfieldc
   mfield1.msmooth13(vfield,vfieldc,s1.ffc,tfield,nx)
# store selected fourier modes: updates curit
   mfield1.mrdvmodes1(vfieldc,curit,tfield,nx,in1.modesxji)
# write diagnostic output: updates njirec
   mdiag1.dafwritevc1(curit,tdiag,iuji,in1.njirec,in1.modesxji)
# transform smoothed ion current to real space: updates vfield
   if ((in1.ndji==1) or (in1.ndji==3)):
      mfft1.mfft1crn(vfieldc,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
      mgard1.mcguard1(vfield,tguard,nx)
# ion current spectral analysis
   if ((in1.ndji==2) or (in1.ndji==3)):
      global itji
      itji += 1
      ts = in1.dt*float(ntime)
# performs frequency analysis of accumulated complex vector time series
# zero out mode 0
      curit[:,0] = numpy.complex(0.0,0.0)
      mdiag1.mivcspect1(curit,s1.wmi,vpkwji,vpksji,ts,in1.t0,tdiag,mtji,
                        s1.iwi,in1.modesxji,nx,-1)
# find frequency with maximum power for each mode
      vwkji[0,:,0] = s1.wmi[numpy.argmax(vpkwji[0,:,:,0],axis=1)]
      vwkji[1,:,0] = s1.wmi[numpy.argmax(vpkwji[1,:,:,0],axis=1)]
      vwkji[0,:,1] = s1.wmi[numpy.argmax(vpkwji[0,:,:,1],axis=1)]
      vwkji[1,:,1] = s1.wmi[numpy.argmax(vpkwji[1,:,:,1],axis=1)]

#-----------------------------------------------------------------------
def del_icurrent_diag13():
   """ delete ion current density diagnostic """
   global curit
   if (in1.njirec > 0):
      in1.closeff(iuji)
      in1.njirec -= 1
   del curit
# spectral analysis
   if ((in1.ndji==2) or (in1.ndji==3)):
      global vpkwji, vwkji, vpksji
      del vpkwji, vwkji, vpksji

#-----------------------------------------------------------------------
def init_vrpotential_diag13():
   """ initialize radiative vector potential diagnostic """
   global vpotr, iuar, oldcu
   farname = "vpotrk1." + s1.cdrun
# write filename to diagnostic metafile
   if (sys.version_info.major==3):
      in1.farname = farname.ljust(in1.farname.dtype.itemsize)
   else:
      in1.farname[:] = farname
# open file
   in1.modesxar = int(min(in1.modesxar,nxh+1))
# vpotr = store selected fourier modes for radiative vector potential
   vpotr = numpy.empty((2,in1.modesxar),complex_type,'F')
# open file: updates narrec and possibly iuar
   if (in1.narrec==0):
      mdiag1.dafopenvc1(vpotr,iuar,in1.narrec,farname)
# oldcu = previous current density with guard cells
   oldcu = numpy.zeros((2,nxe),float_type,'F')
# spectral analysis
   global mtar, itar, vpkwr, vpksr, vwkr
   if ((in1.ndar==2) or (in1.ndar==3)):
      mtar = int((nloop - 1)/in1.ntar) + 1; itar = 0
# vpkwr = power spectrum for radiative vector potential
      vpkwr = numpy.empty((2,in1.modesxar,iwr,2),float_type,'F')
# vpksr = accumulated complex spectrum for radiative vector potential
      vpksr = numpy.zeros((2,4,in1.modesxar,iwr),double_type,'F')
# vwkr = maximum frequency as a function of k for radiative vector
#        potential
      vwkr = numpy.empty((2,in1.modesxar,2),float_type,'F')
# create dummy arrays to avoid undefined arguments later
   else:
      vpkwr = numpy.zeros((1,1,1,1),float_type,'F')
      vwkr = numpy.zeros((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def vrpotential_diag13(vfield,vpkwr,vwkr,ntime):
   """
   radiative vector potential diagnostic
   average current: updates vfieldc = 0.5*(cue + oldcu)
   input/output:
   vfield = scratch array for vector field
   vpkwr = power spectrum for radiative vector potential
   vwkr = maximum frequency as a function of k for radiative vector
          potential
   input:
   ntime = current time step
   """
   mfield1.mcuave1(vfieldc,cue,oldcu,tfield,nx)
# calculate radiative vector potential in fourier space: updates vfieldc
# vfieldc should contain averaged current on entry
   mfield1.mavrpot1(vfieldc,byz,s1.ffc,in1.ci,tfield,nx)
# store selected fourier modes: updates vpotr
   mfield1.mrdvmodes1(vfieldc,vpotr,tfield,nx,in1.modesxar)
# write diagnostic output: updates narrec
   mdiag1.dafwritevc1(vpotr,tdiag,iuar,in1.narrec,in1.modesxar)
# transform radiative vector potential to real space: updates vfield
   if ((in1.ndar==1) or (in1.ndar==3)):
      mfft1.mfft1crn(vfieldc,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
      mgard1.mcguard1(vfield,tguard,nx)
# radiative vector potential spectral analysis
   if ((in1.ndar==2) or (in1.ndar==3)):
      global itar
      itar += 1
      ts = in1.dt*float(ntime)
# performs frequency analysis of accumulated complex vector time series
      mdiag1.mivcspect1(vpotr,wmr,vpkwr,vpksr,ts,in1.t0,tdiag,mtar,iwr,
                        in1.modesxar,nx,1)
# find frequency with maximum power for each mode
      vwkr[0,:,0] = wmr[numpy.argmax(vpkwr[0,:,:,0],axis=1)]
      vwkr[1,:,0] = wmr[numpy.argmax(vpkwr[1,:,:,0],axis=1)]
      vwkr[0,:,1] = wmr[numpy.argmax(vpkwr[0,:,:,1],axis=1)]
      vwkr[1,:,1] = wmr[numpy.argmax(vpkwr[1,:,:,1],axis=1)]

#-----------------------------------------------------------------------
def del_vrpotential_diag13():
   """ delete radiative vector potential diagnostic """
   global vpotr, oldcu
   if (in1.narrec > 0):
      in1.closeff(iuar)
      in1.narrec -= 1
   del vpotr, oldcu
# spectral analysis
   if ((in1.ndar==2) or (in1.ndar==3)):
      global vpkwr, vpksr, vwkr
      del vpkwr, vpksr, vwkr
   in1.ceng = affp

#-----------------------------------------------------------------------
def init_vpotential_diag13():
   """ initialize vector potential diagnostic """
   global vpott, iua
   faname = "vpotk1." + s1.cdrun
# write filename to diagnostic metafile
   if (sys.version_info.major==3):
      in1.faname = faname.ljust(in1.faname.dtype.itemsize)
   else:
      in1.faname[:] = faname
# open file
   in1.modesxa = int(min(in1.modesxa,nxh+1))
# vpott = store selected fourier modes for vector potential
   vpott = numpy.empty((2,in1.modesxa),complex_type,'F')
# open file: updates narec and possibly iua
   if (in1.narec==0):
      mdiag1.dafopenvc1(vpott,iua,in1.narec,faname)
# spectral analysis
   global mta, ita, vpkw, vpks, vwk
   if ((in1.nda==2) or (in1.nda==3)):
      mta = int((nloop - 1)/in1.nta) + 1; ita = 0
# vpkw = power spectrum for vector potential
      vpkw = numpy.empty((2,in1.modesxa,iwr,2),float_type,'F')
# vpks = accumulated complex spectrum for vector potential
      vpks = numpy.zeros((2,4,in1.modesxa,iwr),double_type,'F')
# vwk = maximum frequency as a function of k for vector potential
      vwk = numpy.empty((2,in1.modesxa,2),float_type,'F')
# create dummy arrays to avoid undefined arguments later
   else:
      vpkw = numpy.zeros((1,1,1,1),float_type,'F')
      vwk = numpy.zeros((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def vpotential_diag13(vfield,vpkwi,vwk,ntime):
   """
   vector potential diagnostic
   input/output:
   vfield = scratch array for vector field
   vpkw = power spectrum for vector potential
   vwk = maximum frequency as a function of k for vector potential
   input:
   ntime = current time step
   """
# calculate vector potential in fourier space: updates vfieldc
   mfield1.mavpot1(byz,vfieldc,tfield,nx)
# store selected fourier modes: updates vpott
   mfield1.mrdvmodes1(vfieldc,vpott,tfield,nx,in1.modesxa)
# write diagnostic output: updates narec
   mdiag1.dafwritevc1(vpott,tdiag,iua,in1.narec,in1.modesxa)
# transform vector potential to real space: updates vfield
   if ((in1.nda==1) or (in1.nda==3)):
      mfft1.mfft1crn(vfieldc,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
      mgard1.mcguard1(vfield,tguard,nx)
# vector potential spectral analysis
   if ((in1.nda==2) or (in1.nda==3)):
      global ita
      ita += 1
      ts = in1.dt*float(ntime)
# performs frequency analysis of accumulated complex vector time series
      mdiag1.mivcspect1(vpott,wmr,vpkw,vpks,ts,in1.t0,tdiag,mta,iwr,
                        in1.modesxa,nx,1)
# find frequency with maximum power for each mode
      vwk[0,:,0] = wmr[numpy.argmax(vpkw[0,:,:,0],axis=1)]
      vwk[1,:,0] = wmr[numpy.argmax(vpkw[1,:,:,0],axis=1)]
      vwk[0,:,1] = wmr[numpy.argmax(vpkw[0,:,:,1],axis=1)]
      vwk[1,:,1] = wmr[numpy.argmax(vpkw[1,:,:,1],axis=1)]

#-----------------------------------------------------------------------
def del_vpotential_diag13():
   """ delete vector potential diagnostic """
   global vpott
   if (in1.narec > 0):
      in1.closeff(iua)
      in1.narec -= 1
   del vpott
# spectral analysis
   if ((in1.nda==2) or (in1.nda==3)):
      global vpkw, vpks, vwk
      del vpkw, vpks, vwk
   in1.ceng = affp

#-----------------------------------------------------------------------
def init_etfield_diag13():
   """ initialize transverse efield diagnostic """
   global ett, iuet
   fetname = "etk1." + s1.cdrun
# write filename to diagnostic metafile
   if (sys.version_info.major==3):
      in1.fetname = fetname.ljust(in1.fetname.dtype.itemsize)
   else:
      in1.fetname[:] = fetname
# open file
   in1.modesxet = int(min(in1.modesxet,nxh+1))
# ett = store selected fourier modes for transverse efield
   ett = numpy.empty((2,in1.modesxet),complex_type,'F')
# open file: updates netrec and possibly iuet
   if (in1.netrec==0):
      mdiag1.dafopenvc1(ett,iuet,in1.netrec,fetname)
# spectral analysis
   global mtet, itet, vpkwet, vpkset, vwket
   if ((in1.ndet==2) or (in1.ndet==3)):
      mtet = int((nloop - 1)/in1.ntet) + 1; itet = 0
# vpkwet = power spectrum for transverse efield
      vpkwet = numpy.empty((2,in1.modesxet,iwr,2),float_type,'F')
# vpkset = accumulated complex spectrum for transverse efield
      vpkset = numpy.zeros((2,4,in1.modesxet,iwr),double_type,'F')
# vwket = maximum frequency as a function of k for transverse efield
      vwket = numpy.empty((2,in1.modesxet,2),float_type,'F')
# create dummy arrays to avoid undefined arguments later
   else:
      vpkwet = numpy.zeros((1,1,1,1),float_type,'F')
      vwket = numpy.zeros((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def etfield_diag13(vfield,vpkwet,vwket,ntime):
   """
   transverse efield diagnostic
   input/output:
   vfield = scratch array for vector field
   vpkwet = power spectrum for transverse efield
   vwket = maximum frequency as a function of k for transverse efield
   input:
   ntime = current time step
   """
# store selected fourier modes: updates ett
   mfield1.mrdvmodes1(eyz,ett,tfield,nx,in1.modesxet)
# write diagnostic output: updates netrec
   mdiag1.dafwritevc1(ett,tdiag,iuet,in1.netrec,in1.modesxet)
# transform transverse efield to real space: updates vfield
   if ((in1.ndet==1) or (in1.ndet==3)):
      mfft1.mfft1crn(eyz,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
      mgard1.mcguard1(vfield,tguard,nx)
# spectral analysis
   if ((in1.ndet==2) or (in1.ndet==3)):
      global itet
      itet += 1
      ts = in1.dt*float(ntime)
# performs frequency analysis of accumulated complex vector time series
      mdiag1.mivcspect1(ett,wmr,vpkwet,vpkset,ts,in1.t0,tdiag,mtet,iwr,
                        in1.modesxet,nx,0)
# find frequency with maximum power for each mode
      vwket[0,:,0] = wmr[numpy.argmax(vpkwet[0,:,:,0],axis=1)]
      vwket[1,:,0] = wmr[numpy.argmax(vpkwet[1,:,:,0],axis=1)]
      vwket[0,:,1] = wmr[numpy.argmax(vpkwet[0,:,:,1],axis=1)]
      vwket[1,:,1] = wmr[numpy.argmax(vpkwet[1,:,:,1],axis=1)]

#-----------------------------------------------------------------------
def del_etfield_diag13():
   """ delete transverse efield diagnostic """
   global ett
   if (in1.netrec > 0):
      in1.closeff(iuet)
      in1.netrec -= 1
   del ett
# spectral analysis
   if ((in1.ndet==2) or (in1.ndet==3)):
      global vpkwet, vpkset, vwket
      del vpkwet, vpkset, vwket
   in1.ceng = affp

#-----------------------------------------------------------------------
def init_bfield_diag13():
   """ initialize magnetic field diagnostic """
   global bt, iub
   fbname = "bk1." + s1.cdrun
# write filename to diagnostic metafile
   if (sys.version_info.major==3):
      in1.fbname = fbname.ljust(in1.fbname.dtype.itemsize)
   else:
      in1.fbname[:] = fbname
# open file
   in1.modesxb = int(min(in1.modesxb,nxh+1))
# bt = store selected fourier modes for magnetic field
   bt = numpy.empty((2,in1.modesxb),complex_type,'F')
# open file: updates nbrec and possibly iub
   if (in1.nbrec==0):
      mdiag1.dafopenvc1(bt,iub,in1.nbrec,fbname)

#-----------------------------------------------------------------------
def bfield_diag13(vfield):
   """
   magnetic field diagnostic
   input/output:
   vfield = scratch array for vector field
   """
# store selected fourier modes: updates bt
   mfield1.mrdvmodes1(byz,bt,tfield,nx,in1.modesxb)
# write diagnostic output: updates nbrec
   mdiag1.dafwritevc1(bt,tdiag,iub,in1.nbrec,in1.modesxb)
# transform magnetic field to real space: updates vfield
   mfft1.mfft1crn(byz,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
   mgard1.mcguard1(vfield,tguard,nx)

#-----------------------------------------------------------------------
def del_bfield_diag13():
   """ delete magnetic field diagnostic """
   global bt
   if (in1.nbrec > 0):
      in1.closeff(iub)
      in1.nbrec -= 1
   del bt
   in1.ceng = affp

#-----------------------------------------------------------------------
def init_efluidms_diag13():
   """ initialize electron fluid moments diagnostic """
   # calculate first dimension of fluid arrays
   if (in1.npro==1):
      in1.nprd = 1
   elif (in1.npro==2):
      in1.nprd = 4
   elif (in1.npro==3):
      in1.nprd = 10
   elif (in1.npro==4):
      in1.nprd = 14
   if ((in1.ndfm==1) or (in1.ndfm==3)):
# fmse = electron fluid moments
      s1.fmse = numpy.empty((in1.nprd,nxe),float_type,'F')
# open file for real data: updates nferec and possibly iufe
      ffename = "fmer1." + s1.cdrun
# write filename to diagnostic metafile
      if (sys.version_info.major==3):
         in1.ffename = ffename.ljust(in1.ffename.dtype.itemsize)
      else:
         in1.ffename[:] = ffename
# open file
      if (in1.nferec==0):
         mdiag1.dafopenv1(s1.fmse,nx,s1.iufe,in1.nferec,ffename)

#-----------------------------------------------------------------------
def efluidms_diag13(fmse):
   """
   electron fluid moments diagnostic
   input/output:
   fmse = electron fluid moments
   """
# calculate electron fluid moments
   if ((in1.ndfm==1) or (in1.ndfm==3)):
      s1.dtimer(dtime,itime,-1)
      fmse.fill(0.0)
      s1.dtimer(dtime,itime,1)
      tdiag[0] += float(dtime)
      mdiag1.wmprofx13(s1.ppart,fmse,s1.kpic,in1.ci,tdiag,in1.npro,
                       in1.mx,in1.relativity)
# add guard cells with OpenMP: updates fmse
      mgard1.mamcguard1(fmse,tdiag,nx)
# calculates fluid quantities from fluid moments: updates fmse
      mdiag1.mfluidqs13(fmse,tdiag,in1.npro,nx)
# write real space diagnostic output: updates nferec
      mdiag1.dafwritev1(fmse,tdiag,s1.iufe,in1.nferec,nx)

#-----------------------------------------------------------------------
def init_ifluidms_diag13():
   """ initialize ion fluid moments diagnostic """
   # calculate first dimension of fluid arrays
   if (in1.npro==1):
      in1.nprd = 1
   elif (in1.npro==2):
      in1.nprd = 4
   elif (in1.npro==3):
      in1.nprd = 10
   elif (in1.npro==4):
      in1.nprd = 14
   if ((in1.ndfm==2) or (in1.ndfm==3)):
# fmsi = ion fluid moments
      s1.fmsi = numpy.empty((in1.nprd,nxe),float_type,'F')
# open file for real data: updates nfirec and possibly iufi
      ffiname = "fmir1." + s1.cdrun
# write filename to diagnostic metafile
      if (sys.version_info.major==3):
         in1.ffiname = ffiname.ljust(in1.ffiname.dtype.itemsize)
      else:
         in1.ffiname[:] = ffiname
# open file
      if (in1.nfirec==0):
          mdiag1.dafopenv1(s1.fmsi,nx,s1.iufi,in1.nfirec,ffiname)

#-----------------------------------------------------------------------
def ifluidms_diag13(fmsi):
   """
   ion fluid moments diagnostic
   input/output:
   fmsi = ion fluid moments
   """
   if ((in1.ndfm==2) or (in1.ndfm==3)):
      s1.dtimer(dtime,itime,-1)
      fmsi.fill(0.0)
      s1.dtimer(dtime,itime,1)
      tdiag[0] += float(dtime)
      mdiag1.wmprofx13(s1.pparti,fmsi,s1.kipic,in1.ci,tdiag,in1.npro,
                       in1.mx,in1.relativity)
# add guard cells with OpenMP: updates fmsi
      mgard1.mamcguard1(fmsi,tdiag,nx)
# calculates fluid quantities from fluid moments: updates fmsi
      mdiag1.mfluidqs13(fmsi,tdiag,in1.npro,nx)
      fmsi[:,:]  = in1.rmass*fmsi
# write real space diagnostic output: updates nfirec
      mdiag1.dafwritev1(fmsi,tdiag,s1.iufi,in1.nfirec,nx)

#-----------------------------------------------------------------------
def init_evelocity_diag13():
   """ initialize electron velocity diagnostic """
   in1.nfvd = 0; in1.nfed = 0
   if ((in1.nvft==1) or (in1.nvft==3)):
      in1.nfvd = in1.ndim
   elif ((in1.nvft==4) or (in1.nvft==5)):
      in1.nfvd = 2
   if ((in1.nvft==2) or (in1.nvft==3) or (in1.nvft==5)):
      in1.nfed = 1
   s1.mtv = int((nloop - 1)/in1.ntv) + 1; s1.itv = 0
   s1.eci = in1.ci
   if (in1.relativity==0):
      s1.eci = 0.0
   s1.wkt[:] = 0.0
   if ((in1.ndv==1) or (in1.ndv==3)):
# estimate maximum electron velocity or momentum
      ws = 0.0
      if (in1.npx > 0):
         ws = 4.0*in1.vtx+abs(in1.vx0)
         ws = max(ws,4.0*in1.vty+abs(in1.vy0))
         ws = max(ws,4.0*in1.vtz+abs(in1.vz0))
      if (in1.npxb > 0):
         ws = max(ws,4.0*in1.vtdx+abs(in1.vdx))
         ws = max(ws,4.0*in1.vtdy+abs(in1.vdy))
         ws = max(ws,4.0*in1.vtdz+abs(in1.vdz))
# fv = global electron velocity distribution function
      s1.fv = numpy.empty((2*in1.nmv+2,in1.nfvd),float_type,'F')
# fvm = electron vdrift, vth, entropy for global distribution
      s1.fvm = numpy.zeros((in1.ndim,3),float_type,'F')
# fe = global electron energy distribution function
      s1.fe = numpy.empty((2*in1.nmv+2,in1.nfed),float_type,'F')
# sfv = electron velocity distribution function in tile
      s1.sfv = numpy.empty((2*in1.nmv+2,in1.ndim,mx1+1),float_type,'F')
# open file for electron velocity data: updates nverec and possibly iuve
      fvename = "fve1." + s1.cdrun
# write filename to diagnostic metafile
      if (sys.version_info.major==3):
         in1.fvename = fvename.ljust(in1.fvename.dtype.itemsize)
      else:
         in1.fvename[:] = fvename
      if (in1.nverec==0):
         mdiag1.dafopenfv1(s1.fvm,s1.fv,s1.fe,s1.wkt,s1.iuve,in1.nverec,
                           fvename)
# cartesian distribution
      if ((in1.nvft==1) or (in1.nvft==3)):
# fvtm = time history of electron vdrift, vth, and entropy
         s1.fvtm = numpy.zeros((s1.mtv,in1.ndim,3),float_type,'F')
# set velocity or momentum scale
         s1.fv[2*in1.nmv+1,:] = 2.0*ws
# create dummy array to avoid undefined arguments later
      else:
         s1.fvtm = numpy.empty((1,1,1),float_type,'F')
# cylindrical distribution
      if ((in1.nvft==4) or (in1.nvft==5)):
# set velocity or momentum scale
         s1.fv[2*in1.nmv+1,:] = 2.0*ws
# energy distribution
      if ((in1.nvft==2) or (in1.nvft==3) or (in1.nvft==5)):
# set energy scale for electrons
         ws = ws*ws
         s1.fe[2*in1.nmv+1,0] = ws/(1.0 + 
                                    numpy.sqrt(1.0 + ws*s1.eci*s1.eci))
# create dummy arrays to avoid undefined arguments later
   else:
      s1.fv = numpy.empty((1,1),float_type,'F')
      s1.fe = numpy.empty((1,1),float_type,'F')
      s1.fvm = numpy.zeros((1,1),float_type,'F')
      s1.fvtm = numpy.zeros((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def evelocity_diag13(fv,fe,fvm,fvtm,wkt):
   """
   electron velocity diagnostic
   input/output:
   fv = global electron velocity distribution functions
   fe = global electron energy distribution function
   fvm = electron vdrift, vth, entropy for global distribution
   fvtm = time history of electron vdrift, vth, and entropy
   wkt = total energy contained in distribution
   """
   if ((in1.ndv==1) or (in1.ndv==3)):
# calculate electron cartesian distribution function and moments
      if ((in1.nvft==1) or (in1.nvft==3)):
         s1.sfv[2*in1.nmv+1,:,mx1] = fv[2*in1.nmv+1,:] 
         mdiag1.mvpdist1(s1.ppart,s1.kpic,s1.sfv,fvm,tdiag,np,in1.nmv)
         fv[:,:] = s1.sfv[:,:,mx1] 
# store time history electron vdrift, vth, and entropy
         fvtm[s1.itv,:,:] = fvm
         s1.itv += 1
# calculate electron cylindrical distribution function and moments
      if ((in1.nvft==4) or (in1.nvft==5)):
         s1.sfv[2*in1.nmv+1,0:in1.nfvd,mx1] = fv[2*in1.nmv+1,0:in1.nfvd] 
         mdiag1.mvbpdist1(s1.ppart,s1.kpic,s1.sfv,fvm,in1.omx,in1.omy,
                          in1.omz,tdiag,np,in1.nmv)
         fv[:,:] = s1.sfv[:,0:in1.nfvd,mx1] 
# electron energy distribution
      if ((in1.nvft==2) or (in1.nvft==3) or (in1.nvft==5)):
         s1.sfv[2*in1.nmv+1,0,mx1] = fe[2*in1.nmv+1,0] 
         mdiag1.merpdist1(s1.ppart,s1.kpic,s1.sfv,s1.eci,wkt,tdiag,
                          in1.nmv) 
         fe[:,0] = s1.sfv[:,0,mx1] 
# write electron velocity-space diagnostic output: updates nverec
      mdiag1.dafwritefv1(fvm,fv,fe,wkt,tdiag,s1.iuve,in1.nverec)

#-----------------------------------------------------------------------
def init_ivelocity_diag13():
   """ initialize ion velocity diagnostic """
   in1.nfvd = 0; in1.nfed = 0
   if ((in1.nvft==1) or (in1.nvft==3)):
      in1.nfvd = in1.ndim
   elif ((in1.nvft==4) or (in1.nvft==5)):
      in1.nfvd = 2
   if ((in1.nvft==2) or (in1.nvft==3) or (in1.nvft==5)):
      in1.nfed = 1
   s1.mtv = int((nloop - 1)/in1.ntv) + 1; s1.itv = 0
   s1.eci = in1.ci
   if (in1.relativity==0):
      s1.eci = 0.0
   s1.wkt[:] = 0.0
   if ((in1.ndv==2) or (in1.ndv==3)):
# estimate maximum ion velocity or momentum
      ws = 0.0
      if (in1.npxi > 0):
         ws = 4.0*vtxi+abs(in1.vxi0)
         ws = max(ws,4.0*vtyi+abs(in1.vyi0))
         ws = max(ws,4.0*vtzi+abs(in1.vzi0))
      if (in1.npxbi > 0):
         ws = max(ws,4.0*vtdxi+abs(in1.vdxi))
         ws = max(ws,4.0*vtdyi+abs(in1.vdyi))
         ws = max(ws,4.0*vtdzi+abs(in1.vdzi))
# fvi = global ion velocity distribution function
      s1.fvi = numpy.empty((2*in1.nmv+2,in1.nfvd),float_type,'F')
# fvmi = ion vdrift, vth, entropy for global distribution
      s1.fvmi = numpy.zeros((in1.ndim,3),float_type,'F')
# fei = global ion energy distribution function
      s1.fei = numpy.empty((2*in1.nmv+2,in1.nfed),float_type,'F')
# sfvi = ion velocity distribution function in tile
      s1.sfvi = numpy.empty((2*in1.nmv+2,in1.ndim,mx1+1),float_type,'F')
# open file for ion velocity data: updates nvirec and possibly iuvi
      fviname = "fvi1." + s1.cdrun
# write filename to diagnostic metafile
      if (sys.version_info.major==3):
         in1.fviname = fviname.ljust(in1.fviname.dtype.itemsize)
      else:
         in1.fviname[:] = fviname
      if (in1.nvirec==0):
         mdiag1.dafopenfv1(s1.fvmi,s1.fvi,s1.fei,s1.wkt,s1.iuvi,
                           in1.nvirec,fviname)
# cartesian distribution
      if ((in1.nvft==1) or (in1.nvft==3)):
# fvtm = time history of ion vdrift, vth, and entropy
         s1.fvtmi = numpy.zeros((s1.mtv,in1.ndim,3),float_type,'F')
# set velocity or momentum scale
         s1.fvi[2*in1.nmv+1,:] = 2.0*ws
# create dummy array to avoid undefined arguments later
      else:
         s1.fvtmi = numpy.empty((1,1,1),float_type,'F')
# cylindrical distribution
      if ((in1.nvft==4) or (in1.nvft==5)):
         s1.fvi[2*in1.nmv+1,:] = 2.0*ws
# energy distribution
      if ((in1.nvft==2) or (in1.nvft==3) or (in1.nvft==5)):
# set energy scale for ions
         ws = ws*ws
         s1.fei[2*in1.nmv+1,0] = ws/(1.0 + 
                                     numpy.sqrt(1.0 + ws*s1.eci*s1.eci))
# create dummy arrays to avoid undefined arguments later
   else:
      s1.fvi = numpy.empty((1,1),float_type,'F')
      s1.fei = numpy.empty((1,1),float_type,'F')
      s1.fvmi = numpy.zeros((1,1),float_type,'F')
      s1.fvtmi = numpy.zeros((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def ivelocity_diag13(fvi,fei,fvmi,fvtmi,wkt):
   """
   ion velocity diagnostic
   input/output:
   fvi = global ion velocity distribution functions in tile
   fei = global ion energy distribution function
   fvmi = ion vdrift, vth, entropy for global distribution
   fvtmi = time history of ion vdrift, vth, and entropy
   wkt = total energy contained in distribution
   """
   if ((in1.ndv==2) or (in1.ndv==3)):
# calculate ion cartesian distribution function and moments
      if ((in1.nvft==1) or (in1.nvft==3)):
         s1.sfvi[2*in1.nmv+1,:,mx1] = fvi[2*in1.nmv+1,:]
         mdiag1.mvpdist1(s1.pparti,s1.kipic,s1.sfvi,fvmi,tdiag,npi,
                         in1.nmv)
         fvi[:,:] = s1.sfvi[:,:,mx1] 
# update time step if electrons have not been calculated
# store time history of ion vdrift, vth, and entropy
         if (in1.ndv==2):
            s1.itv += 1
         fvtmi[s1.itv-1,:,:] = fvmi
# calculate ion cylindrical distribution function and moments
      if ((in1.nvft==4) or (in1.nvft==5)):
         s1.sfvi[2*in1.nmv+1,0:in1.nfvd,mx1] = fvi[2*in1.nmv+1,0:in1.nfvd]
         mdiag1.mvbpdist1(s1.pparti,s1.kipic,s1.sfvi,fvmi,in1.omx,
                          in1.omy,in1.omz,tdiag,npi,in1.nmv)
         fvi[:,:] = s1.sfvi[:,0:in1.nfvd,mx1] 
# ion energy distribution
      if ((in1.nvft==2) or (in1.nvft==3) or (in1.nvft==5)):
         s1.sfvi[2*in1.nmv+1,0,mx1] = fei[2*in1.nmv+1,0]
         mdiag1.merpdist1(s1.pparti,s1.kipic,s1.sfvi,s1.eci,wkt,tdiag,
                          in1.nmv)
         fei[:,0] = s1.sfvi[:,0,mx1]
         wkt[:] = in1.rmass*wkt
# write ion velocity-space diagnostic output: updates nvirec
      mdiag1.dafwritefv1(fvmi,fvi,fei,wkt,tdiag,s1.iuvi,in1.nvirec)

#-----------------------------------------------------------------------
def init_traj_diag13(ntime):
   """
   initialize trajectory diagnostic
   input:
   ntime = current time step
   """
#set initial test trajectories
   if ((ntime+s1.ntime0)==0):
      if ((in1.ndt==2) and (in1.movion==0)):
         in1.ndt = 0
      if ((in1.ndt==1) or (in1.ndt==2)):
# iprobt = scratch array
         iprobt = numpy.empty((in1.nprobt),numpy.int32)
# electron trajectories
      if (in1.ndt==1):
# sets electron test charge distribution: updates ppart, iprobt, nprobt
         mdiag1.setptraj1(s1.ppart,s1.kpic,iprobt,in1.nst,in1.vtx,
                          in1.vtsx,in1.dvtx,np,in1.nprobt,irc)
         if (irc[0] != 0):
            print ("esetptraj1 error: irc=", irc[0])
            exit(1)
# estimate maximum electron velocity or momentum
         if (in1.nst==3):
            ws = 0.0
            if (in1.npx > 0):
               ws = 4.0*in1.vtx+abs(in1.vx0)
               ws = max(ws,4.0*in1.vty+abs(in1.vy0))
               ws = max(ws,4.0*in1.vtz+abs(in1.vz0))
            if (in1.npxb > 0):
               ws = max(ws,4.0*in1.vtdx+abs(in1.vdx))
               ws = max(ws,4.0*in1.vtdy+abs(in1.vdy))
               ws = max(ws,4.0*in1.vtdz+abs(in1.vdz))
# ion trajectories
      elif (in1.ndt==2):
# sets ion test charge distribution: updates pparti, iprobt, nprobt
         mdiag1.setptraj1(s1.pparti,s1.kipic,iprobt,in1.nst,vtxi,
                          in1.vtsx,in1.dvtx,npi,in1.nprobt,irc)
         if (irc[0] != 0):
            print ("isetptraj1 error: irc=", irc[0])
            exit(1)
# estimate maximum ion velocity or momentum
         if (in1.nst==3):
            ws = 0.0
            if (in1.npxi > 0):
               ws = 4.0*vtxi+abs(in1.vxi0)
               ws = max(ws,4.0*vtyi+abs(in1.vyi0))
               ws = max(ws,4.0*vtzi+abs(in1.vzi0))
            if (in1.npxbi > 0):
               ws = max(ws,4.0*vtdxi+abs(in1.vdxi))
               ws = max(ws,4.0*vtdyi+abs(in1.vdyi))
               ws = max(ws,4.0*vtdzi+abs(in1.vdzi))
      if ("iprobt" in globals()):
         del iprobt
# find number of existing test trajectories: updates nprobt
   else:
      if (in1.ndt==1):
         mdiag1.mfnptraj1(s1.ppart,s1.kpic,in1.nprobt,irc)
         if (in1.nst==3):
            ws = 0.0
            if (in1.npx > 0):
               ws = 4.0*in1.vtx+abs(in1.vx0)
               ws = max(ws,4.0*in1.vty+abs(in1.vy0))
               ws = max(ws,4.0*in1.vtz+abs(in1.vz0))
            if (in1.npxb > 0):
               ws = max(ws,4.0*in1.vtdx+abs(in1.vdx))
               ws = max(ws,4.0*in1.vtdy+abs(in1.vdy))
               ws = max(ws,4.0*in1.vtdz+abs(in1.vdz))
      elif (in1.ndt==2):
         mdiag1.mfnptraj1(s1.pparti,s1.kipic,in1.nprobt,irc)
         if (in1.nst==3):
            ws = 0.0
            if (in1.npxi > 0):
               ws = 4.0*vtxi+abs(in1.vxi0)
               ws = max(ws,4.0*vtyi+abs(in1.vyi0))
               ws = max(ws,4.0*vtzi+abs(in1.vzi0))
            if (in1.npxbi > 0):
               ws = max(ws,4.0*vtdxi+abs(in1.vdxi))
               ws = max(ws,4.0*vtdyi+abs(in1.vdyi))
               ws = max(ws,4.0*vtdzi+abs(in1.vdzi))
      if (irc[0] != 0):
         print ("mfnptraj1 error: irc=", irc[0])
         exit(1)
# electron or ion trajectories
   if ((in1.ndt==1) or (in1.ndt==2)):
      if (in1.nprobt > 16777215):
         print ("nprobt overflow = ", in1.nprobt)
         exit(1)
      in1.ndimp = idimp
# partt = particle trajectories tracked
      s1.partt = numpy.empty((idimp,in1.nprobt),float_type,'F')
      ftname = "tr1." + s1.cdrun
# write filename to diagnostic metafile
      if (sys.version_info.major==3):
         in1.ftname = ftname.ljust(in1.ftname.dtype.itemsize)
      else:
         in1.ftname[:] = ftname
# track particle trajectories
      if ((in1.nst==1) or (in1.nst==2)):
         s1.mtt = int((nloop - 1)/in1.ntt) + 1
         s1.itt = 0
# partd = trajectory time history array
         s1.partd = numpy.zeros((s1.mtt,idimp,in1.nprobt),
                                float_type,'F')
# create dummy arrays to avoid undefined arguments later
         s1.fvtp = numpy.zeros((1,1),float_type,'F')
         s1.fvmtp = numpy.zeros((1,1),float_type,'F')
         s1.fetp = numpy.zeros((1,1),float_type,'F')
# open file for trajectory data: updates ntrec and possibly iut
         if (in1.ntrec==0):
            mdiag1.dafopentr1(s1.partt,s1.iut,in1.ntrec,ftname)
# calculate test particle distribution function and moments
      elif (in1.nst==3):
# fvtp = velocity distribution function for test particles
         s1.fvtp = numpy.empty((2*in1.nmv+2,in1.ndim),float_type,'F')
# fvmtp = vdrift, vth, and entropy for test particles
         s1.fvmtp = numpy.empty((in1.ndim,3),float_type,'F')
         s1.fvtp[0,:] = 2.0*max(4.0*in1.vtx+abs(in1.vx0),
                             4.0*in1.vtdx+abs(in1.vdx))
# fetp = energy distribution function for test particles
         s1.fetp = numpy.empty((2*in1.nmv+2,0),float_type,'F')
         s1.fvtp[2*in1.nmv+1,:] = 2.0*ws
# open file for test particle diagnostic: updates ntrec and possibly iut
         if (in1.ntrec==0):
            ws = 0.0
            mdiag1.dafopenfv1(s1.fvmtp,s1.fvtp,s1.fetp,ws,s1.iut,
                              in1.ntrec,ftname)
# create dummy array to avoid undefined arguments later
         s1.partd = numpy.empty((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def traj_diag13(partd,fvtp,fvmtp):
   """
   trajectory diagnostic
   input/output:
   partd = trajectory time history array
   fvtp = velocity distribution function for test particles
   fvmtp = vdrift, vth, and entropy for test particles
   """
# copies trajectories to array partt
   if (in1.ndt==1):
      mdiag1.mptraj1(s1.ppart,s1.kpic,s1.partt,tdiag,irc)
# copies tagged ion coordinatess to array partti
   elif (in1.ndt==2):
      mdiag1.mptraj1(s1.pparti,s1.kipic,s1.partt,tdiag,irc)
   if (irc[0] != 0):
      exit(1)
# electron or ion trajectories
   if ((in1.ndt==1) or (in1.ndt==2)):
# store particle trajectories
      if ((in1.nst==1) or (in1.nst==2)):
# write trajectory diagnostic output: updates ntrec
         mdiag1.dafwritetr1(s1.partt,tdiag,s1.iut,in1.ntrec)
         partd[s1.itt,:,:] = s1.partt
         s1.itt += 1
# store cartesian distribution functions
      elif (in1.nst==3):
# calculate test particle distribution function and moments
         mdiag1.mvdist1(s1.partt,fvtp,fvmtp,tdiag,in1.nprobt,in1.nmv)
# write test particle diagnostic output: updates ntrec
         ws = 0.0
         mdiag1.dafwritefv1(fvmtp,fvtp,s1.fetp,ws,tdiag,s1.iut,
                            in1.ntrec)

#-----------------------------------------------------------------------
def init_ephasesp_diag13():
   """ initialize electron phase space diagnostic """
   in1.mvx = min(in1.mvx,nx)
   in1.nsxb = int((nx - 1)/in1.mvx + 1)
# electron phase space diagnostic
   if ((in1.nds==1) or (in1.nds==3)):
# estimate maximum electron velocity or momentum
      ws = 0.0
      if (in1.npx > 0):
         ws = 4.0*in1.vtx+abs(in1.vx0)
         ws = max(ws,4.0*in1.vty+abs(in1.vy0))
         ws = max(ws,4.0*in1.vtz+abs(in1.vz0))
      if (in1.npxb > 0):
         ws = max(ws,4.0*in1.vtdx+abs(in1.vdx))
         ws = max(ws,4.0*in1.vtdy+abs(in1.vdy))
         ws = max(ws,4.0*in1.vtdz+abs(in1.vdz))
# fvs = global electron phase space distribution function
      s1.fvs = numpy.zeros((2*in1.nmv+2,in1.ndim,in1.nsxb),float_type,'F')
      s1.fvs[2*in1.nmv+1,:,0] = 1.25*ws
# open file for electron phase space data:
# updates nserec and possibly iuse
# opens a new fortran unformatted stream file
      if (in1.nserec==0):
         fsename = "pse1." + s1.cdrun
# write filename to diagnostic metafile
         if (sys.version_info.major==3):
            in1.fsename = fsename.ljust(in1.fsename.dtype.itemsize)
         else:
            in1.fsename[:] = fsename
         s1.iuse = mdiag1.get_funit(s1.iuse)
         mdiag1.fnopens1(s1.iuse,fsename)
         in1.nserec = 1

#-----------------------------------------------------------------------
def  init_iphasesp_diag13():
# initialize ion phase space diagnostic
   in1.mvx = min(in1.mvx,nx)
   in1.nsxb = int((nx - 1)/in1.mvx + 1)
   if ((in1.nds==2) or (in1.nds==3)):
# estimate maximum ion velocity or momentum
      ws = 0.0
      if (in1.npxi > 0):
         ws = 4.0*vtxi+abs(in1.vxi0)
         ws = max(ws,4.0*vtyi+abs(in1.vyi0))
         ws = max(ws,4.0*vtzi+abs(in1.vzi0))
      if (in1.npxbi > 0):
         ws = max(ws,4.0*vtdxi+abs(in1.vdxi))
         ws = max(ws,4.0*vtdyi+abs(in1.vdyi))
         ws = max(ws,4.0*vtdzi+abs(in1.vdzi))
# fvsi = global electron phase space distribution function
      s1.fvsi = numpy.zeros((2*in1.nmv+2,in1.ndim,in1.nsxb),float_type,'F')
      s1.fvsi[2*in1.nmv+1,:,0] = 1.25*ws
# open file for ion phase space data:
# updates nsirec and possibly iusi
# opens a new fortran unformatted stream file
      if (in1.nsirec==0):
         fsiname = "psi1." + s1.cdrun
# write filename to diagnostic metafile
         if (sys.version_info.major==3):
            in1.fsiname = fsiname.ljust(in1.fsiname.dtype.itemsize)
         else:
            in1.fsiname[:] = fsiname
         s1.iusi = mdiag1.get_funit(s1.iusi)
         mdiag1.fnopens1(s1.iusi,fsiname)
         in1.nsirec = 1

#-----------------------------------------------------------------------
def print_timings13(tinit,tloop,iuot):
   """
   print timing summaries
   input:
   tinit = initialization elapsed time
   tloop = loop elapsed time
   iuot = output file descriptor
   """
   print (file=iuot)
   print ("initialization time = ",tinit,file=iuot)
   print ("deposit time = ",s1.tdpost[0],file=iuot)
   print ("current deposit time = ",tdjpost[0],file=iuot)
   s1.tdpost[0] += tdjpost[0]
   print ("total deposit time = ",s1.tdpost[0],file=iuot)
   tguard[0] += s1.tguard[0]
   print ("guard time = ",tguard[0],file=iuot)
   tfield[0] += s1.tfield[0]
   print ("solver time = ",tfield[0],file=iuot)
   print ("fft time = ",s1.tfft[0],file=iuot)
   print ("push time = ",tpush[0],file=iuot)
   print ("sort time = ",tsort[0],file=iuot)
   tfield[0] += tguard[0] + s1.tfft[0]
   print ( "total solver time = ",tfield[0],file=iuot)
   time = s1.tdpost[0] + tpush[0] + tsort[0]
   print ("total particle time = ",time,file=iuot)
   tdiag[0] += s1.tdiag[0]
   print ("total diagnostic time = ",tdiag[0],file=iuot)
   ws[0] = time + tfield[0] + tdiag[0]
   tloop = tloop - ws[0]
   print ("total and additional time = ",ws[0],",",tloop,file=iuot)
   print (file=iuot)
# summarize particle timings
   ws[0] = 1.0e+09/(float(nloop)*float(np+npi))
   print ("Push Time (nsec) = ",tpush[0]*ws[0],file=iuot)
   print ("Deposit Time (nsec) = ",s1.tdpost[0]*ws[0],file=iuot)
   print ("Sort Time (nsec) = ",tsort[0]*ws[0],file=iuot)
   print ("Total Particle Time (nsec) = ",time*ws[0],file=iuot)
   print (file=iuot)

#-----------------------------------------------------------------------
def reset_diags13():
   """ reset electrostatic/electromagnetic diagnostics """
# reset electrostatic diagnostics
   s1.reset_diags1()
   if (in1.ntw > 0):
      s1.wt.fill(0.0)
      s.fill(0.0)
# reset electromagnetic diagnostics
   if (in1.ntje > 0):
      if (in1.njerec > 1):
         in1.njerec = 1
   if (in1.ntar > 0):
      if (in1.narrec > 1):
         in1.narrec = 1
      if ((in1.ndar==2) or (in1.ndar==3)):
         itar = 0; vpksr.fill(0.0)
   if (in1.nta > 0):
      if (in1.narec > 1):
         in1.narec = 1
      if ((in1.nda==2) or (in1.nda==3)):
            ita = 0; vpks.fill(0.0)
   if (in1.ntet > 0):
      if (in1.netrec > 1):
         in1.netrec = 1
      if ((in1.ndet==2) or (in1.ndet==3)):
         itet = 0; vpkset.fill(0.0)
   if (in1.ntb > 0):
      if (in1.nbrec > 1):
         in1.nbrec = 1
   if (in1.movion==1):
      if (in1.ntji > 0):
         if (in1.njirec > 1):
            in1.njirec = 1
         if ((in1.ndji==2) or (in1.ndji==3)):
            itji = 0; vpksji.fill(0.0)

#-----------------------------------------------------------------------
def close_diags13(iudm):
   """
   close diagnostics
   delete data, close fortran files, and write out diagnostic metafile
   iudm = diagnostic metafile fortran file descriptor
   """
# electron density diagnostic
   if (in1.ntde > 0):
      s1.del_edensity_diag1()
# potential diagnostic
   if (in1.ntp > 0):
      s1.del_potential_diag1()
# longitudinal efield diagnostic
   if (in1.ntel > 0):
      s1.del_elfield_diag1()
# fluid moments diagnostic
   if (in1.ntfm > 0) :
# electrons
      s1.del_efluidms_diag1()
# ions
      if (in1.movion==1):
         s1.del_ifluidms_diag1()
# electron current diagnostic
   if (in1.ntje > 0):
      del_ecurrent_diag13()
# radiative vector potential diagnostic
   if (in1.ntar > 0):
      del_vrpotential_diag13()
# vector potential diagnostic
   if (in1.nta > 0):
      del_vpotential_diag13()
# transverse efield diagnostic
   if (in1.ntet > 0):
      del_etfield_diag13()
# magnetic field diagnostic
   if (in1.ntb > 0):
      del_bfield_diag13()
# velocity diagnostic
   if (in1.ntv > 0):
      s1.del_evelocity_diag1()
      if (in1.movion==1):
         s1.del_ivelocity_diag1()
# trajectory diagnostic
   if (in1.ntt > 0):
      s1.del_traj_diag1()
# phase space diagnostic
   if (in1.nts > 0):
      s1.del_ephasesp_diag1()
      if (in1.movion==1):
         s1.del_iphasesp_diag1()
# ion diagnostics
   if (in1.movion==1):
# ion density diagnostic
      if (in1.ntdi > 0):
         s1.del_idensity_diag1()
# ion current diagnostic
      if (in1.ntji > 0):
         del_icurrent_diag13()
# write final diagnostic metafile
   in1.writnml1(iudm)
# deallocate arrays
   del_fields13()
   s1.del_electrons1()
   if (in1.movion==1):
      del_ions13()
   del s1.part, s1.ppbuff, s1.ncl, s1.ihole, s1.nppmx
   if ((in1.ntde>0) or (in1.ntp>0) or (in1.ntel>0) or (in1.ntdi>0)):
      s1.del_spectrum1()
   if ((in1.ntje>0) or (in1.nta>0) or (in1.ntet>0) or (in1.ntb>0) or 
       (in1.ntar>0) or (in1.ntji>0)):
      del_spectrum13()
   if (in1.ntw > 0):
      del_energy_diag13()

#-----------------------------------------------------------------------
def initialize_diagnostics13(ntime):
   """
   initialize all diagnostics from namelist input parameters
   input:
   ntime = current time step
   """
# initialize energy diagnostic: allocates wt
   if (in1.ntw > 0):
      init_energy_diag13()

# allocate and initialize scratch arrays for scalar fields:
# allocates sfield
   if ((in1.ntde > 0) or (in1.ntp > 0) or (in1.ntel > 0) or (in1.ntdi > 0)):
      s1.init_spectrum1()

# allocate and initialize scratch arrays for vector fields:
# allocates vfield
   if ((in1.ntje>0) or (in1.nta>0) or (in1.ntet>0) or (in1.ntb>0) or 
       (in1.ntar>0) or (in1.ntji>0)):
      init_spectrum13()

# initialize electron density diagnostic
   if (in1.ntde > 0):
      s1.init_edensity_diag1()

# initialize ion density diagnostic: allocates pkwdi, wkdi
   if (in1.movion==1):
      if (in1.ntdi > 0):
         s1.init_idensity_diag1()

# initialize potential diagnostic: allocates pkw, wk
   if (in1.ntp > 0):
      s1.init_potential_diag1()

# initialize longitudinal efield diagnostic
   if (in1.ntel > 0):
      s1.init_elfield_diag1()

# initialize electron current density diagnostic
   if (in1.ntje > 0):
      init_ecurrent_diag13()

# initialize ion current density diagnostic: allocates vpkwji, vwkji
   if (in1.movion==1):
      if (in1.ntji > 0):
         init_icurrent_diag13()

# initialize radiative vector potential diagnostic:
# allocates vpkwr, vwkr, oldcu
   if (in1.ntar > 0):
         init_vrpotential_diag13()

# initialize vector potential diagnostic: allocates vpkw, vwk
   if (in1.nta > 0):
      init_vpotential_diag13()

# initialize transverse efield diagnostic: allocates vpkwet, vwket
   if (in1.ntet > 0):
      init_etfield_diag13()

# initialize magnetic field diagnostic
   if (in1.ntb > 0):
      init_bfield_diag13()

# initialize fluid moments diagnostic
   if (in1.ntfm > 0):
# electrons: allocates fmse
      init_efluidms_diag13()
# ions: allocates fmsi
      if (in1.movion==1):
         init_ifluidms_diag13()

# initialize velocity diagnostic
   if (in1.ntv > 0):
# electrons: allocates fv, fvm, fvtm
      init_evelocity_diag13()
# ions: allocates fvi, fvmi, fvtmi
      if (in1.movion==1):
         init_ivelocity_diag13()

# initialize trajectory diagnostic: allocates partd, fvtp, fvmtp
   if (in1.ntt > 0):
      init_traj_diag13(ntime)

# initialize phase space diagnostic:
   if (in1.nts > 0):
# electrons: allocates fvs
      init_ephasesp_diag13()
# ions: allocates fvsi
      if (in1.movion==1):
         init_iphasesp_diag13()

#-----------------------------------------------------------------------
def bwrite_restart13(iur,ntime):
   """
   write out basic restart file for electromagnetic code
   input:
   iur = restart file descriptor
   ntime = current time step
   """
# write out particles and electrostatic fields
   s1.bwrite_restart1(iur,ntime)
# write out electromagnetic fields
   nxvh = numpy.size(eyz,1)
   s1.i1[0] = nxvh; s1.i1.tofile(iur)
   if (nxvh > 0):
      eyz.tofile(iur)
      byz.tofile(iur)

#-----------------------------------------------------------------------
def bread_restart13(iur):
   """
   read in basic restart file for electromagnetic code
   input:
   iur = restart file descriptor
   """
# read in particles and electrostatic fields
   s1.bread_restart1(iur)
   global cui
# allocate ion current
   if ("cui" not in globals()):
      cui = numpy.zeros((2,nxe),float_type,'F')
# read in electromagnetic fields
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it > numpy.size(eyz,1)):
      print ("eyz/byz restart error, size(eyz)=",it,numpy.size(eyz,1))
      exit(1)
   if (it > 0):
      il = 2*it
      eyz[:,:] = numpy.fromfile(iur,complex_type,il).reshape(2,it)
      byz[:,:] = numpy.fromfile(iur,complex_type,il).reshape(2,it)

#-----------------------------------------------------------------------
def dwrite_restart13(iur):
   """
   write out restart diagnostic file for electromagnetic code
   input:
   iur = restart file descriptor
   """
# write out restart diagnostic file for electrostatic code
   s1.dwrite_restart1(iur)

# write out electron current density diagnostic parameter
   s1.i1[0] = in1.ntje; s1.i1.tofile(iur)
# write out record location
   if (in1.ntje > 0):
      s1.i1[0] = in1.njerec; s1.i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.njerec > 0):
         it = mdiag1.fnrecl(in1.fjename)
         s1.i1[0] = it; s1.i1.tofile(iur)
         if (it > 0):
            if (sys.version_info.major==3):
               s1.fname = in1.fjename
            else:
               s1.fname[:] = ''.join(in1.fjename)
            s1.fname.tofile(iur)

# write out radiative vector potential diagnostic parameter
   s1.i1[0] = in1.ntar; s1.i1.tofile(iur)
# write out record location
   if (in1.ntar > 0):
      s1.i1[0] = in1.narrec; s1.i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.narrec > 0):
         it = mdiag1.fnrecl(in1.farname)
         s1.i1[0] = it; s1.i1.tofile(iur)
         if (it > 0):
            if (sys.version_info.major==3):
               s1.fname = in1.farname
            else:
               s1.fname[:] = ''.join(in1.farname)
            s1.fname.tofile(iur)
# write out current density
      it = numpy.size(cue,1)
      s1.i1[0] = it; s1.i1.tofile(iur)
      cue.tofile(iur)
# write out spectrum flag
      if ((in1.ndar==2) or (in1.ndar==3)):
         s1.i1[0] = itar; s1.i1.tofile(iur)
# write out spectrum sizes and data
         if (itar > 0):
            s1.i4[0] = numpy.size(vpksr,0)
            s1.i4[1] = numpy.size(vpksr,1)
            s1.i4[2] = numpy.size(vpksr,2)
            s1.i4[3] = numpy.size(vpksr,3)
            s1.i4.tofile(iur)
            vpksr.tofile(iur)
      else:
         it = 0
         s1.i1[0] = it; s1.i1.tofile(iur)

# write out vector potential diagnostic parameter
   s1.i1[0] = in1.nta; s1.i1.tofile(iur)
# write out record location
   if (in1.nta > 0):
      s1.i1[0] = in1.narec; s1.i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.narec > 0):
         it = mdiag1.fnrecl(in1.faname)
         s1.i1[0] = it; s1.i1.tofile(iur)
         if (it > 0):
            if (sys.version_info.major==3):
               s1.fname = in1.faname
            else:
               s1.fname[:] = ''.join(in1.faname)
            s1.fname.tofile(iur)
# write out spectrum flag
      if ((in1.nda==2) or (in1.nda==3)):
         s1.i1[0] = ita; s1.i1.tofile(iur)
# write out spectrum sizes and data
         if (ita > 0):
            s1.i4[0] = numpy.size(vpks,0); s1.i4[1] = numpy.size(vpks,1)
            s1.i4[2] = numpy.size(vpks,2); s1.i4[3] = numpy.size(vpks,3)
            s1.i4.tofile(iur)
            vpks.tofile(iur)
      else:
         it = 0
         s1.i1[0] = it; s1.i1.tofile(iur)

# write out transverse efield diagnostic parameter
   s1.i1[0] = in1.ntet; s1.i1.tofile(iur)
# write out record location
   if (in1.ntet > 0):
      s1.i1[0] = in1.netrec; s1.i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.netrec > 0):
         it = mdiag1.fnrecl(in1.fetname)
         s1.i1[0] = it; s1.i1.tofile(iur)
         if (it > 0):
            if (sys.version_info.major==3):
               s1.fname = in1.fetname
            else:
               s1.fname[:] = ''.join(in1.fetname)
            s1.fname.tofile(iur)
# write out spectrum flag
      if ((in1.ndet==2) or (in1.ndet==3)):
         s1.i1[0] = itet; s1.i1.tofile(iur)
# write out spectrum sizes and data
         if (itet > 0):
            s1.i4[0] = numpy.size(vpkset,0)
            s1.i4[1] = numpy.size(vpkset,1)
            s1.i4[2] = numpy.size(vpkset,2)
            s1.i4[3] = numpy.size(vpkset,3)
            s1.i4.tofile(iur)
            vpkset.tofile(iur)
      else:
         it = 0
         s1.i1[0] = it; s1.i1.tofile(iur)

# write out magnetic field diagnostic diagnostic parameter
   s1.i1[0] = in1.ntb; s1.i1.tofile(iur)
# write out record location
   if (in1.ntb > 0):
      s1.i1[0] = in1.nbrec; s1.i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
      if (in1.nbrec > 0):
         it = mdiag1.fnrecl(in1.fbname)
         s1.i1[0] = it; s1.i1.tofile(iur)
         if (it > 0):
            if (sys.version_info.major==3):
               s1.fname = in1.fbname
            else:
               s1.fname[:] = ''.join(in1.fbname)
            s1.fname.tofile(iur)

# write out ion current density diagnostic parameter
   if (in1.movion==1):
      s1.i1[0] = in1.ntji; s1.i1.tofile(iur)
# write out record location
      if (in1.ntji > 0):
         s1.i1[0] = in1.njirec; s1.i1.tofile(iur)
# write out record length (zero if error) and file name (if no error)
         if (in1.njirec > 0):
            it = mdiag1.fnrecl(in1.fjiname)
            s1.i1[0] = it; s1.i1.tofile(iur)
            if (it > 0):
               if (sys.version_info.major==3):
                  s1.fname = in1.fjiname
               else:
                  s1.fname[:] = ''.join(in1.fjiname)
               s1.fname.tofile(iur)
# write out spectrum flag
         if ((in1.ndji==2) or (in1.ndji==3)):
            s1.i1[0] = itji; s1.i1.tofile(iur)
# write out spectrum sizes and data
            if (itji > 0):
               s1.i4[0] = numpy.size(vpksji,0)
               s1.i4[1] = numpy.size(vpksji,1)
               s1.i4[2] = numpy.size(vpksji,2)
               s1.i4[3] = numpy.size(vpksji,3)
               s1.i4.tofile(iur)
               vpksji.tofile(iur)
         else:
            it = 0
            s1.i1[0] = it; s1.i1.tofile(iur)

#-----------------------------------------------------------------------
def dread_restart13(iur):
   """
   read in restart diagnostic file for electromagnetic code
   input:
   iur = restart file descriptor electrostatic code
   """
   global itw, itv, itt, ita, itet, itar, itji
# read in restart diagnostic file for electrostatic code
   s1.dread_restart1(iur)
# restore energy accumulations
   if (in1.ntw > 0):
      if (s1.itw > 0):
         s[0] = s1.wt[0,4]
         s[1] = s1.wt[0,1]
         s[2] = s1.wt[0,5]
         s[3] = s1.wt[0,6]
         s[4] = s1.wt[0,2]
         s[5] = s1.wt[0,3]
         s[6] = s[5]
         for it in range(1,s1.itw):
            s[0] += s1.wt[it,4]
            s[1] += s1.wt[it,1]
            s[2] += s1.wt[it,5]
            s[3] += s1.wt[it,6]
            s[4] += s1.wt[it,2]
            s[5] = min(s[5],s1.wt[it,3])
            s[6] = max(s[6],s1.wt[it,3])

# read in electron current density diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.ntje):
      print ("restart error: read/expected ntje=",it,in1.ntje)
      exit(1)
# read in record location
   if (in1.ntje > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.njerec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.njerec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print ("ntje zero length record error")
            exit(1)
         if (sys.version_info.major==3):
            in1.fjename = numpy.fromfile(iur,'S1',32)
         else:
            s1.fname[:] = numpy.fromfile(iur,'S32',1)
            in1.fjename[:] = str(s1.fname[0])

# read in radiative vector potential diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.ntar):
      print ("restart error: read/expected ntar=",it,in1.ntar)
      exit(1)
# read in record location
   if (in1.ntar > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.narrec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.narrec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print ("ntar zero length record error")
            exit(1)
         if (sys.version_info.major==3):
            in1.farame = numpy.fromfile(iur,'S1',32)
         else:
            s1.fname[:] = numpy.fromfile(iur,'S32',1)
            in1.farname[:] = str(s1.fname[0])
# read in current density
      s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
      if (it > numpy.size(cue,1)):
         print ("cue size error: read/expected=",it,numpy.size(cue,1))
         exit(1)
      il = 2*it
      cue[:,:] = numpy.fromfile(iur,float_type,il).reshape(2,it)
# read in spectrum flag
      s1.i1[:] = numpy.fromfile(iur,int_type,1); itar = s1.i1[0]
# read in spectrum sizes and data
      if (itar > 0):
         s1.i4[:] = numpy.fromfile(iur,int_type,4)
         iq = s1.i4[0]; ir = s1.i4[1]; it = s1.i4[2]; ip = s1.i4[3]
         if (iq != 2):
            print ("vpksr size error: read/expected 2 =",iq)
            exit(1)
         if (ir != 4):
            print ("vpksr size error: read/expected 4 =",ir)
            exit(1)
         if (it != in1.modesxar):
            print ("vpksr size error: read/expected modesxar=",it,
                    in1.modesxar)
            exit(1)
         if (ip != iwr):
            print ("vpksr size error: read/expected iwr=",ip,iwr)
            exit(1)
         il = iq*ir*it*ip
         vpksr[:,:,:,:] = (numpy.fromfile(iur,double_type,il).
                           reshape(2,4,it,ip))

# read in vector potential diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.nta):
      print ("restart error: read/expected nta=",it,in1.nta)
      exit(1)
# read in record location
   if (in1.nta > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.narec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.narec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print ("nta zero length record error")
            exit(1)
         if (sys.version_info.major==3):
            in1.faname = numpy.fromfile(iur,'S1',32)
         else:
            s1.fname[:] = numpy.fromfile(iur,'S32',1)
            in1.faname[:] = str(s1.fname[0])
# read in spectrum flag
      s1.i1[:] = numpy.fromfile(iur,int_type,1); ita = s1.i1[0]
# read in spectrum sizes and data
      if (ita > 0):
         s1.i4[:] = numpy.fromfile(iur,int_type,4)
         iq = s1.i4[0]; ir = s1.i4[1]; it = s1.i4[2]; ip = s1.i4[3]
         if (iq != 2):
            print ("vpks size error: read/expected 2 =",iq)
            exit(1)
         if (ir != 4):
            print ("vpks size error: read/expected 4 =",ir)
            exit(1)
         if (it != in1.modesxa):
            print ("vpks size error: read/expected modesxa=",it,
                    in1.modesxa)
            exit(1)
         if (ip != iwr):
            print ("vpks size error: read/expected iwr=",ip,iwr)
            exit(1)
         il = iq*ir*it*ip
         vpks[:,:,:,:] = (numpy.fromfile(iur,double_type,il).
                          reshape(2,4,it,ip))

# read in transverse efield diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.ntet):
      print ("restart error: read/expected ntet=",it,in1.ntet)
      exit(1)
# read in record location
   if (in1.ntet > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.netrec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.netrec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print ("ntet zero length record error")
            exit(1)
         if (sys.version_info.major==3):
            in1.fetame = numpy.fromfile(iur,'S1',32)
         else:
            s1.fname[:] = numpy.fromfile(iur,'S32',1)
            in1.fetname[:] = str(s1.fname[0])
# read in spectrum flag
      s1.i1[:] = numpy.fromfile(iur,int_type,1); itet = s1.i1[0]
# read in spectrum sizes and data
      if (itet > 0):
         s1.i4[:] = numpy.fromfile(iur,int_type,4)
         iq = s1.i4[0]; ir = s1.i4[1]; it = s1.i4[2]; ip = s1.i4[3]
         if (iq != 2):
            print ("vpkset size error: read/expected 2 =",iq)
            exit(1)
         if (ir != 4):
            print ("vpkset size error: read/expected 4 =",ir)
            exit(1)
         if (it != in1.modesxet):
            print ("vpkset size error: read/expected modesxet=",it,
                    in1.modesxet)
            exit(1)
         if (ip != iwr):
            print ("vpkset size error: read/expected iwr=",ip,iwr)
            exit(1)
         il = iq*ir*it*ip
         vpkset[:,:,:,:] = (numpy.fromfile(iur,double_type,il).
                            reshape(2,4,it,ip))

# read in magnetic field diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.ntb):
      print ("restart error: read/expected nteb=",it,in1.ntb)
      exit(1)
# read in record location
   if (in1.ntb > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.nbrec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.nbrec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print ("ntb zero length record error")
            exit(1)
         if (sys.version_info.major==3):
            in1.fbname = numpy.fromfile(iur,'S1',32)
         else:
            s1.fname[:] = numpy.fromfile(iur,'S32',1)
            in1.fbname[:] = str(s1.fname[0])

# read in ion current density diagnostic parameter
   if (in1.movion==1):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
      if (it != in1.ntji):
         print ("restart error: read/expected ntji=",it,in1.ntji)
         exit(1)
# read in record location
      if (in1.ntji > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.njirec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
         if (in1.njirec > 0):
            s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
            if (it==0):
               print ("ntji zero length record error")
               exit(1)
            if (sys.version_info.major==3):
               in1.fjiname = numpy.fromfile(iur,'S1',32)
            else:
               s1.fname[:] = numpy.fromfile(iur,'S32',1)
               in1.fjiname[:] = str(s1.fname[0])
# read in spectrum flag
         s1.i1[:] = numpy.fromfile(iur,int_type,1); itji = s1.i1[0]
# read in spectrum sizes and data
         if (itji > 0):
            s1.i4[:] = numpy.fromfile(iur,int_type,4)
            iq = s1.i4[0]; ir = s1.i4[1]; it = s1.i4[2]; ip = s1.i4[3]
            if (iq != 2):
               print ("vpksji size error: read/expected 2 =",iq)
               exit(1)
            if (ir != 4):
               print ("vpksji size error: read/expected 4 =",ir)
               exit(1)
            if (it != in1.modesxji):
               print ("vpksji size error: read/expected modesxji=",it,
                       in1.modesxji)
               exit(1)
            if (ip != s1.iwi):
               print ("vpksji size error: read/expected iwi=",ip,s1.iwi)
               exit(1)
            il = iq*ir*it*ip
            vpksji[:,:,:,:] = (numpy.fromfile(iur,double_type,il).
                               reshape(2,4,it,ip))
