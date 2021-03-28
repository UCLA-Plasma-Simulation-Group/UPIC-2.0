#-----------------------------------------------------------------------
"""
High Level library for 1-2/2D Darwin OpenMP PIC code

functions defined:

init_dfields13: allocate darwin field data for standard code
del_dfields13: delete darwin field data for standard code

calc_shift13: calculate shift constant for iteration

dpush_electrons13: push darwin electrons with OpenMP

dpush_ions13: push darwin ions with OpenMP

d_time_reverse1: start running simulation backwards

denergy_diag13: darwin energy diagnostic

init_dspectrum13: allocate scratch arrays for darwin vector fields
del_dspectrum13: delete scratch arrays for darwin vector fields

init_edcurrent_diag13: initialize darwin electron current density
                       diagnostic
edcurrent_diag13: darwin electron current density diagnostic
del_edcurrent_diag13: delete darwin electron current density diagnostic

init_vdpotential_diag13: initialize darwin vector potential diagnostic
vdpotential_diag13: darwin vector potential diagnostic
del_vdpotential_diag13: delete darwin vector potential diagnostic

init_detfield_diag13: initialize darwin transverse efield diagnostic
detfield_diag13: darwin transverse efield diagnostic
del_detfield_diag13: delete darwin transverse efield diagnostic

init_dbfield_diag13: initialize darwin magnetic field diagnostic
dbfield_diag13: darwin magnetic field diagnostic
del_dbfield_diag13: delete darwin magnetic field diagnostic

edfluidms_diag13: darwin electron fluid moments diagnostic

idfluidms_diag13: darwin ion fluid moments diagnostic

print_dtimings13: print darwin timing summaries

reset_ddiags13: reset electrostatic/darwin diagnostics
close_ddiags13: close darwin diagnostics

initialize_ddiagnostics13: initialize all diagnostics from namelist
                           input parameters
                           
darwin_predictor13: predictor for darwin iteration
darwin_iteration13: darwin iteration loop

bwrite_drestart13: write out basic restart file for darwin code
bread_drestart13: read in basic restart file for darwin code
dwrite_drestart13: write out restart diagnostic file for darwin code
dread_drestart13: read in restart diagnostic file for darwin code

written by Viktor K. Decyk and Joshua Kelly, UCLA
copyright 1999-2016, regents of the university of california
update: december 11, 2017
"""
import math
import numpy

# sys.path.append('./mbeps1.source')
from libmpush1 import *
from libmbpush1 import *
from libmdpush1 import *
import s1
import sb1

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# idimp = number of particle coordinates = 4
# ipbc = particle boundary condition: 1 = periodic
idimp = 4; ipbc = 1
s1.idimp = idimp

# declare scalars for standard code
npi = 0
ws = numpy.zeros((1),float_type)
wpmax = numpy.empty((1),float_type)
wpmin = numpy.empty((1),float_type)

# declare scalars for OpenMP code
irc = numpy.zeros((1),int_type)
irc2 = numpy.zeros((2),int_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
dtime = numpy.empty((1),double_type)
tguard = numpy.zeros((1),float_type)
tfield = numpy.zeros((1),float_type)
tdcjpost = numpy.zeros((1),float_type)
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
# nloop = number of time steps in simulation
nloop = s1.nloop
qbme = s1.qbme
# affp = float(nx)/float(np)
affp = s1.affp
if (in1.movion==1):
   qbmi = s1.qbmi

# check for unimplemented features
plist = s1.plist

#-----------------------------------------------------------------------
def init_dfields13():
   """ allocate darwin field data for standard code """
   global dcu, cus, amu, ffe, dcui, amui
# allocate electromagnetic field data: fxyze, cue, byze, eyz, byz
# allocate electrostatic field data: qe, qi, fxe, ffc, mixup, sct
   sb1.init_fields13()
   global dcu, cus, amu, ffe
# dcu = electron acceleration density with guard cells
dcu = numpy.empty((2,nxe),float_type,'F')
# cus = transverse electric field
cus = numpy.empty((2,nxe),float_type,'F')
# amu = electron momentum flux with guard cells
amu = numpy.empty((2,nxe),float_type,'F')
# ffe = form factor arrays for poisson solvers
ffe = numpy.empty((nxh),complex_type,'F')
# dcui = ion acceleration density with guard cells
# amui = ion momentum flux with guard cells
if (in1.movion==1):
   dcui = numpy.empty((2,nxe),float_type,'F')
   amui = numpy.empty((2,nxe),float_type,'F')

#-----------------------------------------------------------------------
def del_dfields13():
   """ delete darwin field data for standard code """
   global dcu, cus, amu, ffe
   del dcu, cus, amu, ffe
   if (in1.movion==1):
      global dcui, amui
      del dcui, amui
   sb1.del_fields13()

#-----------------------------------------------------------------------
def calc_shift13(iuot):
   """
   calculate shift constant for iteration: updates q2m0
   input:
   iuot = output file descriptor
   """
   global wpm, q2m0, wpmax, wpmin
# find maximum and minimum initial electron density
   s1.qe.fill(0.0)
   mpush1.mpost1(s1.ppart,s1.qe,s1.kpic,in1.qme,s1.tdpost,in1.mx)
   mgard1.maguard1(s1.qe,s1.tguard,nx)
   if (in1.movion==1):
      s1.qi.fill(0.0)
      mpush1.mpost1(s1.pparti,s1.qi,s1.kipic,in1.qmi,s1.tdpost,in1.mx)
      mgard1.maguard1(s1.qi,s1.tguard,nx)
      mdpush1.mfwptminx1(s1.qe,s1.qi,s1.qbme,s1.qbmi,wpmax,wpmin,nx)
   else:
      mdpush1.mfwpminx1(s1.qe,s1.qbme,wpmax,wpmin,nx)
   wpm = 0.5*(wpmax[0] + wpmin[0])*s1.affp
# accelerate convergence: update wpm
   if (wpm <= 10.0):
      wpm = 0.75*wpm
   print >> iuot, "wpm=",wpm
   q2m0 = wpm/s1.affp

#-----------------------------------------------------------------------
def dpush_electrons13(ppart,kpic):
   """
   push darwin electrons with OpenMP
   input/output:
   ppart = tiled electron particle array
   kpic = number of electrons in each tile
   """
   s1.wke[0] = 0.0
# updates ppart and wke, and possibly ncl, ihole, irc
   if (in1.mzf==0):
      mbpush1.wmbpush1(ppart,sb1.fxyze,sb1.byze,kpic,s1.ncl,s1.ihole,
                       in1.omx,qbme,in1.dt,in1.dt,in1.ci,s1.wke,
                       sb1.tpush,nx,in1.mx,ipbc,in1.relativity,plist,
                       irc)
# zero force: updates ppart, wke and possibly ncl, ihole, and irc
   else:
      mbpush1.wmpush1zf(ppart,kpic,s1.ncl,s1.ihole,in1.dt,in1.ci,s1.wke,
                        sb1.tpush,in1.nx,in1.mx,ipbc,in1.relativity,
                        plist,irc)

# reorder electrons by tile with OpenMP:
# updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
   if (irc[0]==0):
      msort1.wmporder1(ppart,s1.ppbuff,kpic,s1.ncl,s1.ihole,sb1.tsort,
                       nx,in1.mx,plist,irc2)
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
def dpush_ions13(pparti,kipic):
   """
   push darwin ions with OpenMP
   input/output:
   pparti = tiled electron/ion particle arrays
   kipic = number of electrons/ions in each tile
   """
   s1.wki[0] = 0.0
# updates pparti and wki, and possibly ncl, ihole, irc
   if (in1.mzf==0):
      mbpush1.wmbpush1(pparti,sb1.fxyze,sb1.byze,kipic,s1.ncl,s1.ihole,
                       in1.omx,s1.qbmi,in1.dt,in1.dt,in1.ci,s1.wki,
                       sb1.tpush,nx,in1.mx,ipbc,in1.relativity,plist,
                       irc)
   else:
      mpush1.wmpush1zf(pparti,kipic,s1.ncl,s1.ihole,in1.dt,in1.ci,
                       s1.wki,sb1.tpush,nx,in1.mx,ipbc,in1.relativity,
                       plist,irc)
   s1.wki[0] = s1.wki[0]*in1.rmass

# reorder ions by tile with OpenMP:
# updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
   if (irc[0]==0):
      msort1.wmporder1(pparti,s1.ppbuff,kipic,s1.ncl,s1.ihole,sb1.tsort,
                       nx,in1.mx,plist,irc2)
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
def d_time_reverse1():
   """
   start running simulation backwards
   need to reverse time lag in leap-frog integration scheme
   """
   in1.dt = -in1.dt
   ws[0] = 0.0
   mpush1.wmpush1zf(sb1.ppart,sb1.kpic,sb1.ncl,sb1.ihole,in1.dt,in1.ci,
                    ws,sb1.tpush,nx,in1.mx,ipbc,in1.relativity,plist,
                    irc)
   msort1.wmporder1(sb1.ppart,sb1.ppbuff,sb1.kpic,sb1.ncl,sb1.ihole,
                    sb1.tsort,nx,in1.mx,plist,irc2)
   if (irc2[0] != 0): exit(1)
   if (in1.movion==1):
      mpush1.wmpush1zf(sb1.pparti,sb1.kipic,sb1.ncl,sb1.ihole,in1.dt,
                       in1.ci,ws,sb1.tpush,nx,in1.mx,ipbc,
                       in1.relativity,plist,irc)
      msort1.wmporder1(sb1.pparti,sb1.ppbuff,sb1.kipic,sb1.ncl,
                       sb1.ihole,sb1.tsort,nx,in1.mx,plist,irc2)
      if (irc2[0] != 0): exit(1)

#-----------------------------------------------------------------------
def denergy_diag13(wt,ntime,iuot):
   """
   darwin energy diagnostic
   input/output:
   wt = energy time history array
   input:
   ntime = current time step
   iuot = output file descriptor
   """
   global ws
   sb1.wef[0] = s1.we[0] + sb1.wb[0]
   ws[0] = sb1.wef[0] + s1.wke[0] + s1.wki[0]
   if (ntime==0):
      sb1.s[5] = ws[0]
   if (in1.ndw > 0):
      print >> iuot, "Total Field, Kinetic and Total Energies:"
      if (in1.movion==0):
         iuot.write("%14.7e %14.7e %14.7e\n" % (sb1.wef[0],s1.wke[0],
                    ws[0])) 
      else:
         iuot.write("%14.7e %14.7e %14.7e %14.7e\n" % (sb1.wef[0],
                    s1.wke[0],s1.wki[0],ws[0])) 
      print >> iuot, "Electric(l,t) and Magnetic Energies = "
      iuot.write("%14.7e %14.7e %14.7e\n" % (s1.we[0],sb1.wf[0],
                 sb1.wb[0]))
   wt[s1.itw,:] = ([sb1.wef[0],s1.wke[0],s1.wki[0],ws[0],s1.we[0],
                   sb1.wf[0],sb1.wb[0]])
   s1.itw += 1
   sb1.s[0] += s1.we[0]
   sb1.s[1] += s1.wke[0]
   sb1.s[2] += sb1.wf[0]
   sb1.s[3] += sb1.wb[0]
   sb1.s[4] += s1.wki[0]
   sb1.s[5] = min(sb1.s[5],float(ws[0]))
   sb1.s[6] = max(sb1.s[6],float(ws[0]))

#-----------------------------------------------------------------------
def init_dspectrum13():
   """ allocate scratch arrays for darwin vector fields """
   global vfield, vfieldc
   sb1.init_spectrum13()
   vfield = numpy.empty((2,nxe),float_type,'F')
   vfieldc = numpy.empty((2,nxh),complex_type,'F')
# allocate and initialize frequency array for spectral analysis
   if ((in1.nta > 0) or (in1.ntet>0)):
      global iw, wm
      iw = int((in1.wmax - in1.wmin)/in1.dw + 1.5)
      wm = numpy.empty((iw),float_type,'F')
      wm[:] = in1.wmin + in1.dw*numpy.linspace(0,iw-1,iw)

#-----------------------------------------------------------------------
def del_dspectrum13():
   """ delete scratch arrays for darwin vector fields """
   global vfield, vfieldc
   del vfield, vfieldc
   if ((in1.nta>0) or (in1.ntet>0)):
      global wm
      del wm

#-----------------------------------------------------------------------
def init_edcurrent_diag13():
   """ initialize darwin electron current density diagnostic """
   global curet, iuje, oldcue
   fjename = "curek1." + s1.cdrun
   in1.fjename[:] = fjename
   in1.modesxje = int(min(in1.modesxje,nxh+1))
# oldcue = previous current density
   oldcue = numpy.zeros((2,nxe),float_type,'F')
# curet = store selected fourier modes for electron current density
   curet = numpy.empty((2,in1.modesxje),complex_type,'F')
# open file: updates njerec and possibly iuje
   if (in1.njerec==0):
      mdiag1.dafopenvc1(curet,iuje,in1.njerec,fjename)

#-----------------------------------------------------------------------
def edcurrent_diag13(vfield):
   """
   darwin electron current density diagnostic
   input/output:
   vfield = scratch array for vector field
   """
   vfield[:] = numpy.copy(oldcue)
# transform electron current density to fourier space: updates vfield
   isign = -1
   mfft1.mfft1rn(vfield,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)
# calculate smoothed electron current in fourier space: updates vfieldc
   mfield1.msmooth13(vfield,vfieldc,s1.ffc,tfield,nx)
# store selected fourier modes: updates curet
   mfield1.mrdvmodes1(vfieldc,curet,tfield,nx,in1.modesxje)
# write diagnostic output: updates njerec
   mdiag1.dafwritevc1(curet,tdiag,sb1.iuje,in1.njerec,in1.modesxje)
# transform smoothed electron current to real space: updates vfield
   mfft1.mfft1crn(vfieldc,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
   mgard1.mcguard1(vfield,tguard,nx)

#-----------------------------------------------------------------------
def del_edcurrent_diag13():
   """ delete darwin electron current density diagnostic """
   global curet, oldcue
   if (in1.njerec > 0):
      in1.closeff(iuje)
      in1.njerec -= 1
   del curet, oldcue

#-----------------------------------------------------------------------
def init_vdpotential_diag13():
   """ initialize darwin vector potential diagnostic """
   global vpott
   faname = "vpotk1." + s1.cdrun
   in1.faname[:] = faname
   in1.modesxa = int(min(in1.modesxa,nxh+1))
# vpott = store selected fourier modes for vector potential
   vpott = numpy.empty((2,in1.modesxa),complex_type,'F')
# open file: updates narec and possibly iua
   if (in1.narec==0):
      mdiag1.dafopenvc1(vpott,sb1.iua,in1.narec,faname)
# spectral analysis
   global vpkw, vpks, vwk
   if ((in1.nda==2) or (in1.nda==3)):
      sb1.mta = int((nloop - 1)/in1.nta) + 1; sb1.ita = 0
# vpkw = power spectrum for vector potential
      vpkw = numpy.empty((2,in1.modesxa,iw,2),float_type,'F')
# vpks = accumulated complex spectrum for vector potential
      vpks = numpy.zeros((2,4,in1.modesxa,iw),double_type,'F')
# vwk = maximum frequency as a function of k for vector potential
      vwk = numpy.empty((2,in1.modesxa,2),float_type,'F')
# create dummy arrays to avoid undefined arguments later
   else:
      vpkw = numpy.zeros((1,1,1,1),float_type,'F')
      vwk = numpy.zeros((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def vdpotential_diag13(vfield,vpkwi,vwk,ntime):
   """
   vector darwin potential diagnostic
   input/output:
   vfield = scratch array for vector field
   vpkw = power spectrum for vector potential
   vwk = maximum frequency as a function of k for vector potential
   input:
   ntime = current time step
   """
# calculate vector potential in fourier space: updates vfieldc
   mfield1.mapot1(sb1.cue,vfieldc,s1.ffc,in1.ci,ws,tfield,nx)
# store selected fourier modes: updates vpott
   mfield1.mrdvmodes1(vfieldc,vpott,tfield,nx,in1.modesxa)
# write diagnostic output: updates narec
   mdiag1.dafwritevc1(vpott,tdiag,sb1.iua,in1.narec,in1.modesxa)
# transform vector potential to real space: updates vfield
   if ((in1.nda==1) or (in1.nda==3)):
      mfft1.mfft1crn(vfieldc,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
      mgard1.mcguard1(vfield,tguard,nx)
# vector potential spectral analysis
   if ((in1.nda==2) or (in1.nda==3)):
      global ita
      sb1.ita += 1
      ts = in1.dt*float(ntime)
# performs frequency analysis of accumulated complex vector time series
      mdiag1.mivcspect1(vpott,wm,vpkw,vpks,ts,in1.t0,tdiag,sb1.mta,iw,
                        in1.modesxa,nx,1)
# find frequency with maximum power for each mode
      vwk[0,:,0] = wm[numpy.argmax(vpkw[0,:,:,0],axis=1)]
      vwk[1,:,0] = wm[numpy.argmax(vpkw[1,:,:,0],axis=1)]
      vwk[0,:,1] = wm[numpy.argmax(vpkw[0,:,:,1],axis=1)]
      vwk[1,:,1] = wm[numpy.argmax(vpkw[1,:,:,1],axis=1)]

#-----------------------------------------------------------------------
def del_vdpotential_diag13():
   """ delete darwin vector potential diagnostic """
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
def init_detfield_diag13():
   """ initialize darwin transverse efield diagnostic """
   global ett
   fetname = "etk1." + s1.cdrun
   in1.fetname[:] = fetname
   in1.modesxet = int(min(in1.modesxet,nxh+1))
# ett = store selected fourier modes for transverse efield
   ett = numpy.empty((2,in1.modesxet),complex_type,'F')
# open file: updates netrec and possibly iuet
   if (in1.netrec==0):
      mdiag1.dafopenvc1(ett,sb1.iuet,in1.netrec,fetname)
# spectral analysis
   global vpkwet, vpkset, vwket
   if ((in1.ndet==2) or (in1.ndet==3)):
      sb1.mtet = int((nloop - 1)/in1.ntet) + 1; sb1.itet = 0
# vpkwet = power spectrum for transverse efield
      vpkwet = numpy.empty((2,in1.modesxet,iw,2),float_type,'F')
# vpkset = accumulated complex spectrum for transverse efield
      vpkset = numpy.zeros((2,4,in1.modesxet,iw),double_type,'F')
# vwket = maximum frequency as a function of k for transverse efield
      vwket = numpy.empty((2,in1.modesxet,2),float_type,'F')
# create dummy arrays to avoid undefined arguments later
   else:
      vpkwet = numpy.zeros((1,1,1,1),float_type,'F')
      vwket = numpy.zeros((1,1,1),float_type,'F')

#-----------------------------------------------------------------------
def detfield_diag13(vfield,vpkwet,vwket,ntime):
   """
   transverse darwin efield diagnostic
   input/output:
   vfield = scratch array for vector field
   vpkwet = power spectrum for transverse efield
   vwket = maximum frequency as a function of k for transverse efield
   input:
   ntime = current time step
   """
# calculate transverse efield in fourier space: updates vfieldc
   mfield1.metfield1(dcu,vfieldc,ffe,in1.ci,ws,tfield,nx)
# store selected fourier modes: updates ett
   mfield1.mrdvmodes1(vfieldc,ett,tfield,nx,in1.modesxet)
# write diagnostic output: updates netrec
   mdiag1.dafwritevc1(ett,tdiag,sb1.iuet,in1.netrec,in1.modesxet)
# transform transverse efield to real space: updates vfield
   if ((in1.ndet==1) or (in1.ndet==3)):
      mfft1.mfft1crn(vfieldc,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
      mgard1.mcguard1(vfield,tguard,nx)
# spectral analysis
   if ((in1.ndet==2) or (in1.ndet==3)):
      global itet
      sb1.itet += 1
      ts = in1.dt*float(ntime)
# performs frequency analysis of accumulated complex vector time series
      mdiag1.mivcspect1(ett,wm,vpkwet,vpkset,ts,in1.t0,tdiag,sb1.mtet,
                        iw,in1.modesxet,nx,0)
# find frequency with maximum power for each mode
      vwket[0,:,0] = wm[numpy.argmax(vpkwet[0,:,:,0],axis=1)]
      vwket[1,:,0] = wm[numpy.argmax(vpkwet[1,:,:,0],axis=1)]
      vwket[0,:,1] = wm[numpy.argmax(vpkwet[0,:,:,1],axis=1)]
      vwket[1,:,1] = wm[numpy.argmax(vpkwet[1,:,:,1],axis=1)]

#-----------------------------------------------------------------------
def del_detfield_diag13():
   """ delete darwin transverse efield diagnostic """
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
def init_dbfield_diag13():
   """ initialize darwin magnetic field diagnostic """
   global bt
   fbname = "bk1." + s1.cdrun
   in1.fbname[:] = fbname
   in1.modesxb = int(min(in1.modesxb,nxh+1))
# bt = store selected fourier modes for magnetic field
   bt = numpy.empty((2,in1.modesxb),complex_type,'F')
# open file: updates nbrec and possibly iub
   if (in1.nbrec==0):
      mdiag1.dafopenvc1(bt,sb1.iub,in1.nbrec,fbname)

#-----------------------------------------------------------------------
def dbfield_diag13(vfield):
   """
   darwin magnetic field diagnostic
   input/output:
   vfield = scratch array for vector field
   """
# calculate magnetic field in fourier space: updates vfieldc
   mfield1.mibpois1(sb1.cue,vfieldc,s1.ffc,in1.ci,ws,tfield,nx)
# store selected fourier modes: updates bt
   mfield1.mrdvmodes1(vfieldc,bt,tfield,nx,in1.modesxb)
# write diagnostic output: updates nbrec
   mdiag1.dafwritevc1(bt,tdiag,sb1.iub,in1.nbrec,in1.modesxb)
# transform magnetic field to real space: updates vfield
   mfft1.mfft1crn(vfieldc,vfield,s1.mixup,s1.sct,s1.tfft,in1.indx)
   mgard1.mcguard1(vfield,tguard,nx)

#-----------------------------------------------------------------------
def del_dbfield_diag13():
   """ delete darwin magnetic field diagnostic """
   global bt
   if (in1.nbrec > 0):
      in1.closeff(iub)
      in1.nbrec -= 1
   del bt
   in1.ceng = affp

#-----------------------------------------------------------------------
def edfluidms_diag13(fmse):
   """
   darwin electron fluid moments diagnostic
   input/output:
   fmse = electron fluid moments
   """
# calculate electron fluid moments
   if ((in1.ndfm==1) or (in1.ndfm==3)):
      s1.dtimer(dtime,itime,-1)
      fmse.fill(0.0)
      s1.dtimer(dtime,itime,1)
      tdiag[0] += float(dtime)
      mdiag1.wmgbprofx1(s1.ppart,sb1.fxyze,sb1.byze,fmse,s1.kpic,
                        in1.omx,qbme,in1.dt,in1.ci,tdiag,in1.npro,nx,
                        in1.mx,in1.relativity)
# add guard cells with OpenMP: updates fmse
      mgard1.mamcguard1(fmse,tdiag,nx)
# calculates fluid quantities from fluid moments: updates fmse
      mdiag1.mfluidqs13(fmse,tdiag,in1.npro,nx)
# write real space diagnostic output: updates nferec
      mdiag1.dafwritev1(fmse,tdiag,s1.iufe,in1.nferec,nx)

#-----------------------------------------------------------------------
def idfluidms_diag13(fmsi):
   """
   darwin ion fluid moments diagnostic
   input/output:
   fmsi = ion fluid moments
   """
# calculate ion fluid moments
   if ((in1.ndfm==1) or (in1.ndfm==3)):
      s1.dtimer(dtime,itime,-1)
      fmsi.fill(0.0)
      s1.dtimer(dtime,itime,1)
      tdiag[0] += float(dtime)
      mdiag1.wmgbprofx1(s1.pparti,sb1.fxyze,sb1.byze,fmsi,s1.kipic,
                        in1.omx,qbmi,in1.dt,in1.ci,tdiag,in1.npro,nx,
                        in1.mx,in1.relativity)
# add guard cells with OpenMP: updates fmsi
      mgard1.mamcguard1(fmsi,tdiag,nx)
# calculates fluid quantities from fluid moments: updates fmsi
      mdiag1.mfluidqs13(fmsi,tdiag,in1.npro,nx)
      fmsi[:,:]  = in1.rmass*fmsi
# write real space diagnostic output: updates nfirec
      mdiag1.dafwritev1(fmsi,tdiag,s1.iufi,in1.nfirec,nx)

#-----------------------------------------------------------------------
def print_dtimings13(tinit,tloop,iuot):
   """
   print darwin timing summaries
   input:
   tinit = initialization elapsed time
   tloop = loop elapsed time
   iuot = output file descriptor
   """
   print >> iuot
   print >> iuot, "initialization time = ", tinit
   print >> iuot, "deposit time = ", s1.tdpost[0]
   print >> iuot, "current deposit time = ", sb1.tdjpost[0]
   print >> iuot, "current derivative deposit time = ", tdcjpost[0]
   s1.tdpost[0] += sb1.tdjpost[0] + tdcjpost[0]
   print >> iuot, "total deposit time = ", s1.tdpost[0]
   tguard[0] += sb1.tguard[0] + s1.tguard[0]
   print >> iuot, "guard time = ", tguard[0]
   tfield[0] += sb1.tfield[0] + s1.tfield[0]
   print >> iuot, "solver time = ", tfield[0]
   print >> iuot, "fft time = ", s1.tfft[0]
   print >> iuot, "push time = ", sb1.tpush[0]
   print >> iuot, "sort time = ", sb1.tsort[0]
   tfield[0] += tguard[0] + s1.tfft[0]
   print >> iuot, "total solver time = ", tfield[0]
   time = s1.tdpost[0] + sb1.tpush[0] + sb1.tsort[0]
   print >> iuot, "total particle time = ", time
   tdiag[0] += sb1.tdiag[0] + s1.tdiag[0]
   print >> iuot, "total diagnostic time = ", tdiag[0]
   ws[0] = time + tfield[0] + tdiag[0]
   tloop = tloop - ws[0]
   print >> iuot, "total and additional time = ", ws[0], ",", tloop
   print >> iuot
# summarize particle timings
   ws[0] = 1.0e+09/(float(nloop)*float(np+npi))
   print >> iuot, "Push Time (nsec) = ", sb1.tpush[0]*ws[0]
   print >> iuot, "Deposit Time (nsec) = ", s1.tdpost[0]*ws[0]
   print >> iuot, "Sort Time (nsec) = ", sb1.tsort[0]*ws[0]
   print >> iuot, "Total Particle Time (nsec) = ", time*ws[0]
   print >> iuot

#-----------------------------------------------------------------------
def reset_ddiags13():
   """ reset electrostatic/darwin diagnostics """
# reset electrostatic diagnostics
   s1.reset_diags1()
   if (in1.ntw > 0):
      s1.wt.fill(0.0)
      sb1.s.fill(0.0)
# reset darwin diagnostics
   if (in1.ntje > 0):
      if (in1.njerec > 1):
         in1.njerec = 1
   if (in1.nta > 0):
      if (in1.narec > 1):
         in1.narec = 1
      if ((in1.nda==2) or (in1.nda==3)):
         sb1.ita = 0; vpks.fill(0.0)
   if (in1.ntet > 0):
      if (in1.netrec > 1):
         in1.netrec = 1
      if ((in1.ndet==2) or (in1.ndet==3)):
         sb1.itet = 0; vpkset.fill(0.0)
   if (in1.ntb > 0):
      if (in1.nbrec > 1):
         in1.nbrec = 1
   if (in1.movion==1):
      if (in1.ntji > 0):
         if (in1.njirec > 1):
            in1.njirec = 1
         if ((in1.ndji==2) or (in1.ndji==3)):
            sb1.itji = 0; sb1.vpksji.fill(0.0)

#-----------------------------------------------------------------------
def close_ddiags13(iudm):
   """
   close darwin diagnostics
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
# darwin fluid moments diagnostic
   if (in1.ntfm > 0) :
# electrons
      s1.del_efluidms_diag1()
# ions
      if (in1.movion==1):
         s1.del_ifluidms_diag1()
# darwin electron current diagnostic
   if (in1.ntje > 0):
      del_edcurrent_diag13()
# vector potential diagnostic
   if (in1.nta > 0):
      del_vdpotential_diag13()
# transverse efield diagnostic
   if (in1.ntet > 0):
      del_detfield_diag13()
# magnetic field diagnostic
   if (in1.ntb > 0):
      del_dbfield_diag13()
# ion diagnostics
   if (in1.movion==1):
# ion density diagnostic
      if (in1.ntdi > 0):
         s1.del_idensity_diag1()
# ion current diagnostic
      if (in1.ntji > 0):
         sb1.del_icurrent_diag13()
# write final diagnostic metafile
   in1.writnml1(iudm)
   in1.closeff(iudm)
# deallocate arrays
   del_dfields13()
   s1.del_electrons1()
   if (in1.movion==1):
      sb1.del_ions13()
   del s1.part, s1.ppbuff, s1.ncl, s1.ihole, s1.nppmx
   if ((in1.ntde>0) or (in1.ntp>0) or (in1.ntel>0) or (in1.ntdi>0)):
      s1.del_spectrum1()
   if ((in1.nta>0) or (in1.ntet>0) or (in1.ntb>0) or (in1.ntji>0)):
      del_dspectrum13()
   if (in1.ntw > 0):
      sb1.del_energy_diag13()
   if (in1.ntv > 0):
      s1.del_evelocity_diag1()
      if (in1.movion==1):
         s1.del_ivelocity_diag1()
   if (in1.ntt > 0):
      s1.del_traj_diag1()

#-----------------------------------------------------------------------
def initialize_ddiagnostics13(ntime):
   """
   initialize all diagnostics from namelist input parameters
   input:
   ntime = current time step
   """
# initialize energy diagnostic: allocates wt
   if (in1.ntw > 0):
      sb1.init_energy_diag13()

# allocate and initialize scratch arrays for scalar fields:
# allocates sfield
   if ((in1.ntde>0) or (in1.ntp>0) or (in1.ntel>0) or (in1.ntdi>0)):
      s1.init_spectrum1()

# allocate and initialize scratch arrays for vector fields:
# allocates vfield
   if ((in1.ntje>0) or (in1.nta>0) or (in1.ntet>0) or (in1.ntb>0) or 
       (in1.ntji>0)):
      init_dspectrum13()

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

# initialize darwin electron current density diagnostic
   if (in1.ntje > 0):
      init_edcurrent_diag13()

# initialize ion current density diagnostic: allocates vpkwji, vwkji
   if (in1.movion==1):
      if (in1.ntji > 0):
         sb1.init_icurrent_diag13()

# initialize darwin vector potential diagnostic: allocates vpkw, vwk
   if (in1.nta > 0):
      init_vdpotential_diag13()

# initialize darwin transverse efield diagnostic:
# allocates vpkwet, vwket
   if (in1.ntet > 0):
      init_detfield_diag13()

# initialize darwin magnetic field diagnostic
   if (in1.ntb > 0):
      init_dbfield_diag13()

# initialize fluid moments diagnostic
   if (in1.ntfm > 0):
# electrons: allocates fmse
      sb1.init_efluidms_diag13()
# ions: allocates fmsi
      if (in1.movion==1):
         sb1.init_ifluidms_diag13()

# initialize velocity diagnostic
   if (in1.ntv > 0):
# electrons: allocates fv, fvm, fvtm
      sb1.init_evelocity_diag13()
# ions: allocates fvi, fvmi, fvtmi
      if (in1.movion==1):
         sb1.init_ivelocity_diag13()

# initialize trajectory diagnostic: allocates partd, fvtp, fvmtp
   if (in1.ntt > 0):
      sb1.init_traj_diag13(ntime)

#-----------------------------------------------------------------------
def darwin_predictor13(q2m0):
   """
   predictor for darwin iteration
   input:
   q2m0 = shift constant in darwin iteration
   """
# transform current to fourier space: updates cue
   isign = -1
   mfft1.mfft1rn(sb1.cue,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# calculate magnetic field in fourier space: updates byze, wb
   mfield1.mbbpois1(sb1.cue,sb1.byze,s1.ffc,in1.ci,sb1.wb,s1.tfield,nx)

# transform magnetic force to real space: updates byze
   isign = 1
   mfft1.mfft1rn(sb1.byze,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# add constant to magnetic field: updates byze
   mfield1.mbaddext1(sb1.byze,sb1.tfield,in1.omy,in1.omz,nx)

# copy guard cells: updates byze
   mgard1.mcguard1(sb1.byze,sb1.tguard,nx)

# add longitudinal and old transverse electric fields: updates fxyze
   mfield1.maddvrfield1(sb1.fxyze,cus,s1.fxe,tfield)


# deposit electron acceleration density and momentum flux with OpenMP:
# updates dcu, amu
   mdpush1.wmgdjpost1(s1.ppart,sb1.fxyze,sb1.byze,dcu,amu,s1.kpic,
                      in1.omx,in1.qme,s1.qbme,in1.dt,in1.ci,tdcjpost,nx,
                      in1.mx,in1.relativity)
                      
# deposit ion acceleration density and momentum flux with OpenMP:
# updates dcui, amui
   if (in1.movion==1):
      mdpush1.wmgdjpost1(s1.pparti,sb1.fxyze,sb1.byze,dcui,amui,
                         s1.kipic,in1.omx,in1.qmi,s1.qbmi,in1.dt,in1.ci,
                         tdcjpost,nx,in1.mx,
                         in1.relativity)
# add electron and ion densities: updates dcu, amu
      mfield1.maddcuei1(dcu,dcui,tfield,nxe)
      mfield1.maddcuei1(amu,amui,tfield,nxe)

# add old scaled electric field: updates dcu
   mdpush1.mascfguard1(dcu,cus,q2m0,tdcjpost,nx)

# add guard cells: updates dcu, amu
   mgard1.macguard1(dcu,tguard,nx)
   mgard1.macguard1(amu,tguard,nx)

# transform acceleration density and momentum flux to fourier space:
# updates dcu, amu
   isign = -1
   mfft1.mfft1rn(dcu,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)
   mfft1.mfft1rn(amu,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# take transverse part of time derivative of current: updates dcu
   mfield1.madcuperp1(dcu,amu,tfield,nx)

# calculate transverse electric field: updates cus, wf
   mfield1.mepois1(dcu,cus,ffe,s1.affp,in1.ci,sb1.wf,tfield,nx)

# transform transverse electric field to real space: updates cus
   isign = 1
   mfft1.mfft1rn(cus,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# copy guard cells: updates cus
   mgard1.mcguard1(cus,tguard,nx)

# add longitudinal and transverse electric fields:
# fxyze = cus + fxe, updates fxyze
# cus needs to be retained for next time step
   mfield1.maddvrfield1(sb1.fxyze,cus,s1.fxe,tfield)

#-----------------------------------------------------------------------
def darwin_iteration(q2m0,ntime,k):
   """
   inner iteration loop
   input:
   q2m0 = shift constant in darwin iteration
   """
# deposit electron current and acceleration density and momentum flux
# with OpenMP: updates cue, dcu, amu
   mdpush1.wmgdcjpost1(s1.ppart,sb1.fxyze,sb1.byze,sb1.cue,dcu,amu,
                       s1.kpic,in1.omx,in1.qme,qbme,in1.dt,in1.ci,
                       tdcjpost,nx,in1.mx,in1.relativity)
# add guard cells for current, acceleration density, and momentum flux:
# updates cue, dcu, amu
   mgard1.macguard1(sb1.cue,sb1.tguard,nx)
   mgard1.macguard1(dcu,tguard,nx)
   mgard1.macguard1(amu,tguard,nx)

# save electron current for electron current diagnostic later
   if (k==(in1.ndc-1)):
      if (in1.ntje > 0):
         it = ntime/in1.ntje
         if (ntime==in1.ntje*it):
            oldcue[:] = numpy.copy(sb1.cue)
         
# deposit ion current and acceleration density and momentum flux
# with OpenMP: updates cui, dcui, amui
   if (in1.movion==1):
      mdpush1.wmgdcjpost1(s1.pparti,sb1.fxyze,sb1.byze,sb1.cui,dcui,
                          amui,s1.kipic,in1.omx,in1.qmi,qbmi,in1.dt,
                          in1.ci,tdcjpost,nx,in1.mx,in1.relativity)
# add guard cells for current, acceleration density, and momentum flux:
# updates cui, dcui, amui
      mgard1.macguard1(sb1.cui,sb1.tguard,nx)
      mgard1.macguard1(dcui,tguard,nx)
      mgard1.macguard1(amui,tguard,nx)
# add electron and ion densities: updates cue, dcu, amu
      mfield1.maddcuei1(sb1.cue,sb1.cui,sb1.tfield,nxe)
      mfield1.maddcuei1(dcu,dcui,tfield,nxe)
      mfield1.maddcuei1(amu,amui,tfield,nxe)

# add scaled electric field: updates dcu
   mdpush1.mascfguard1(dcu,cus,q2m0,tdcjpost,nx)

# transform current to fourier space: update cue
   isign = -1
   mfft1.mfft1rn(sb1.cue,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# calculate magnetic field in fourier space: updates byze, wb
   mfield1.mbbpois1(sb1.cue,sb1.byze,s1.ffc,in1.ci,sb1.wb,tfield,nx)

# transform magnetic force to real space: updates byze
   isign = 1
   mfft1.mfft1rn(sb1.byze,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# add constant to magnetic field: updates bzye
   mfield1.mbaddext1(sb1.byze,sb1.tfield,in1.omy,in1.omz,nx)

# transform acceleration density and momentum flux to fourier space:
# updates dcu, amu
   isign = -1
   mfft1.mfft1rn(dcu,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)
   mfft1.mfft1rn(amu,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# take transverse part of time derivative of current: updates dcu
   mfield1.madcuperp1(dcu,amu,tfield,nx)

# calculate transverse electric field: updates cus, wf
   mfield1.mepois1(dcu,cus,ffe,affp,in1.ci,sb1.wf,tfield,nx)

# transform transverse electric field to real space: updates cus
   isign = 1
   mfft1.mfft1rn(cus,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# copy guard cells: updates byze, cus
   mgard1.mcguard1(sb1.byze,sb1.tguard,nx)
   mgard1.mcguard1(cus,tguard,nx)

# add longitudinal and transverse electric fields:
# fxyze = cus + fxe, updates fxyze
# cus needs to be retained for next time step
   mfield1.maddvrfield1(sb1.fxyze,cus,s1.fxe,tfield)

#-----------------------------------------------------------------------
def bwrite_drestart13(iur,ntime):
   """
   write out basic restart file for darwin code
   input:
   iur = restart file descriptor
   ntime = current time step
   """
   global wpm, q2m0
# write out particles and electrostatic fields
   s1.bwrite_restart1(iur,ntime)
# write out shift constant for iteration
   s1.a2[0] = wpm; s1.a2[1] = q2m0
   s1.a2.tofile(iur)
# write out darwin electric field field
   nxv = numpy.size(cus,1)
   s1.i1[0] = nxv; s1.i1.tofile(iur)
   if (nxv > 0):
      cus.tofile(iur)

#-----------------------------------------------------------------------
def bread_drestart13(iur):
   """
   read in basic restart file for darwin code
   input:
   iur = restart file descriptor
   """
# read in particles and electrostatic fields
   s1.bread_restart1(iur)
   global cui, wpm, q2m0
# allocate ion current
   if ("sb1.cui" not in globals()):
      sb1.cui = numpy.zeros((2,nxe),float_type,'F')
# read in shift constant for iteration
   s1.a2[:] = numpy.fromfile(iur,float_type,2)
   wpm = s1.a2[0]; q2m0 = s1.a2[1]
# read in darwin electric field field
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it > numpy.size(cus,1)):
      print "cus restart error, size(cus)=",it,numpy.size(cus,1)
      exit(1)
   if (it > 0):
      il = 2*it
      cus[:,:] = numpy.fromfile(iur,float_type,il).reshape(2,it)

#-----------------------------------------------------------------------
def dwrite_drestart13(iur):
   """
   write out restart diagnostic file for darwin code
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
            s1.fname[:] = ''.join(in1.fjename)
            s1.fname.tofile(iur)

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
            s1.fname[:] = ''.join(in1.faname)
            s1.fname.tofile(iur)
# write out spectrum flag
      if ((in1.nda==2) or (in1.nda==3)):
         s1.i1[0] = sb1.ita; s1.i1.tofile(iur)
# write out spectrum sizes and data
         if (sb1.ita > 0):
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
            s1.fname[:] = ''.join(in1.fetname)
            s1.fname.tofile(iur)
# write out spectrum flag
      if ((in1.ndet==2) or (in1.ndet==3)):
         s1.i1[0] = sb1.itet; s1.i1.tofile(iur)
# write out spectrum sizes and data
         if (sb1.itet > 0):
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
               s1.fname[:] = ''.join(in1.fjiname)
               s1.fname.tofile(iur)
# write out spectrum flag
         if ((in1.ndji==2) or (in1.ndji==3)):
            s1.i1[0] = sb1.itji; s1.i1.tofile(iur)
# write out spectrum sizes and data
            if (sb1.itji > 0):
               s1.i4[0] = numpy.size(sb1.vpksji,0)
               s1.i4[1] = numpy.size(sb1.vpksji,1)
               s1.i4[2] = numpy.size(sb1.vpksji,2)
               s1.i4[3] = numpy.size(sb1.vpksji,3)
               s1.i4.tofile(iur)
               sb1.vpksji.tofile(iur)
         else:
            it = 0
            s1.i1[0] = it; s1.i1.tofile(iur)

#-----------------------------------------------------------------------
def dread_drestart13(iur):
   """
   read in restart diagnostic file for darwin code
   input:
   iur = restart file descriptor electrostatic code
   """
# read in restart diagnostic file for electrostatic code
   s1.dread_restart1(iur)
# restore energy accumulations
   if (in1.ntw > 0):
      if (s1.itw > 0):
         sb1.s[0] = s1.wt[0,4]
         sb1.s[1] = s1.wt[0,1]
         sb1.s[2] = s1.wt[0,5]
         sb1.s[3] = s1.wt[0,6]
         sb1.s[4] = s1.wt[0,2]
         sb1.s[5] = s1.wt[0,3]
         sb1.s[6] = sb1.s[5]
         for it in xrange(1,s1.itw):
            sb1.s[0] += s1.wt[it,4]
            sb1.s[1] += s1.wt[it,1]
            sb1.s[2] += s1.wt[it,5]
            sb1.s[3] += s1.wt[it,6]
            sb1.s[4] += s1.wt[it,2]
            sb1.s[5] = min(sb1.s[5],s1.wt[it,3])
            sb1.s[6] = max(sb1.s[6],s1.wt[it,3])

# read in electron current density diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.ntje):
      print "restart error: read/expected ntje=", it, in1.ntje
      exit(1)
# read in record location
   if (in1.ntje > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.njerec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.njerec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print "ntje zero length record error"
            exit(1)
         s1.fname[:] = numpy.fromfile(iur,'S32',1)
         in1.fjename[:] = str(s1.fname[0])

# read in vector potential diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.nta):
      print "restart error: read/expected nta=", it, in1.nta
      exit(1)
# read in record location
   if (in1.nta > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.narec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.narec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print "nta zero length record error"
            exit(1)
         s1.fname[:] = numpy.fromfile(iur,'S32',1)
         in1.faname[:] = str(s1.fname[0])
# read in spectrum flag
      s1.i1[:] = numpy.fromfile(iur,int_type,1); sb1.ita = s1.i1[0]
# read in spectrum sizes and data
      if (sb1.ita > 0):
         s1.i4[:] = numpy.fromfile(iur,int_type,4)
         iq = s1.i4[0]; ir = s1.i4[1]; it = s1.i4[2]; ip = s1.i4[3]
         if (iq != 2):
            print "vpks size error: read/expected 2 =", iq
            exit(1)
         if (ir != 4):
            print "vpks size error: read/expected 4 =", ir
            exit(1)
         if (it != in1.modesxa):
            print ("vpks size error: read/expected modesxa=", it,
                    in1.modesxa)
            exit(1)
         if (ip != iw):
            print "vpks size error: read/expected iw=", ip, iw
            exit(1)
         il = iq*ir*it*ip
         vpks[:,:,:,:] = (numpy.fromfile(iur,double_type,il).
                          reshape(2,4,it,ip))

# read in transverse efield diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.ntet):
      print "restart error: read/expected ntet=", it, in1.ntet
      exit(1)
# read in record location
   if (in1.ntet > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.netrec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.netrec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print "ntet zero length record error"
            exit(1)
         s1.fname[:] = numpy.fromfile(iur,'S32',1)
         in1.fetname[:] = str(s1.fname[0])
# read in spectrum flag
      s1.i1[:] = numpy.fromfile(iur,int_type,1); sb1.itet = s1.i1[0]
# read in spectrum sizes and data
      if (sb1.itet > 0):
         s1.i4[:] = numpy.fromfile(iur,int_type,4)
         iq = s1.i4[0]; ir = s1.i4[1]; it = s1.i4[2]; ip = s1.i4[3]
         if (iq != 2):
            print "vpkset size error: read/expected 2 =", iq
            exit(1)
         if (ir != 4):
            print "vpkset size error: read/expected 4 =", ir
            exit(1)
         if (it != in1.modesxet):
            print ("vpkset size error: read/expected modesxet=", it,
                    in1.modesxet)
            exit(1)
         if (ip != iw):
            print "vpkset size error: read/expected iw=", ip, iw
            exit(1)
         il = iq*ir*it*ip
         vpkset[:,:,:,:] = (numpy.fromfile(iur,double_type,il).
                            reshape(2,4,it,ip))

# read in magnetic field diagnostic parameter
   s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
   if (it != in1.ntb):
      print "restart error: read/expected nteb=", it, in1.ntb
      exit(1)
# read in record location
   if (in1.ntb > 0):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.nbrec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
      if (in1.nbrec > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
         if (it==0):
            print "ntb zero length record error"
            exit(1)
         s1.fname[:] = numpy.fromfile(iur,'S32',1)
         in1.fbname[:] = str(s1.fname[0])

# read in ion current density diagnostic parameter
   if (in1.movion==1):
      s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
      if (it != in1.ntji):
         print "restart error: read/expected ntji=", it, in1.ntji
         exit(1)
# read in record location
      if (in1.ntji > 0):
         s1.i1[:] = numpy.fromfile(iur,int_type,1); in1.njirec = s1.i1[0]
# read in record length (zero if error) and file name (if no error)
         if (in1.njirec > 0):
            s1.i1[:] = numpy.fromfile(iur,int_type,1); it = s1.i1[0]
            if (it==0):
               print "ntji zero length record error"
               exit(1)
            s1.fname[:] = numpy.fromfile(iur,'S32',1)
            in1.fjiname[:] = str(s1.fname[0])
# read in spectrum flag
         s1.i1[:] = numpy.fromfile(iur,int_type,1); sb1.itji = s1.i1[0]
# read in spectrum sizes and data
         if (sb1.itji > 0):
            s1.i4[:] = numpy.fromfile(iur,int_type,4)
            iq = s1.i4[0]; ir = s1.i4[1]; it = s1.i4[2]; ip = s1.i4[3]
            if (iq != 2):
               print "vpksji size error: read/expected 2 =", iq
               exit(1)
            if (ir != 4):
               print "vpksji size error: read/expected 4 =", ir
               exit(1)
            if (it != in1.modesxji):
               print ("vpksji size error: read/expected modesxji=",it,
                       in1.modesxji)
               exit(1)
            if (ip != s1.iwi):
               print "vpksji size error: read/expected iwi=", ip, s1.iwi
               exit(1)
            il = iq*ir*it*ip
            sb1.vpksji[:,:,:,:] = (numpy.fromfile(iur,double_type,il).
                                   reshape(2,4,it,ip))
