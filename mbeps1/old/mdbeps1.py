#-----------------------------------------------------------------------
# 1-2/2D Darwin OpenMP PIC code
# written by Viktor K. Decyk and Joshua Kelly, UCLA
# copyright 2016, regents of the university of california
from __future__ import print_function
import sys
import math
import numpy

sys.path.append('./mbeps1.source')
from libmpush1 import *
from libmbpush1 import *
from libmdpush1 import *
from fomplib import *
from fgraf1 import *
from dtimer import *

# override default input data
in1.emf = 2
# read namelist
iuin = 8
in1.readnml1(iuin)
# override input data
in1.idcode = 3
in1.ndim = 3
in1.ntar = 0

# import modules after namelist has been read
import s1
import sb1
import sd1

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# ipbc = particle boundary condition: 1 = periodic
ipbc = sd1.ipbc
zero = 0.0

# declare scalars for standard code
npi = 0
ws = numpy.zeros((1),float_type)
wpmax = numpy.empty((1),float_type)
wpmin = numpy.empty((1),float_type)

# declare scalars for OpenMP code
irc = numpy.zeros((1),int_type)
ierr = numpy.zeros((1),int_type)

# declare and initialize timing data
tinit = 0.0; tloop = 0.0
itime = numpy.empty((4),numpy.int32)
ltime = numpy.empty((4),numpy.int32)
dtime = numpy.empty((1),double_type)

# start timing initialization
dtimer(dtime,itime,-1)

# text output file
fname = "output1." + s1.cdrun
iuot = open(fname,"w")

# in1.nvp = number of shared memory nodes (0=default)
#nvp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
omplib.init_omp(in1.nvp)

# open graphics device
irc[0] = graf1.open_graphs(in1.nplot)

# initialize scalars for standard code
# np = total number of electrons in simulation
np = s1.np
# nx = number of grid points in x direction
nx = s1.nx; nxh = int(nx/2)
# npi = total number of ions in simulation
if (in1.movion > 0):
   npi = s1.npi
nxe = nx + 2; nxeh = int(nxe/2)
# nloop = number of time steps in simulation
# nstart = initial time loop index
# ntime = current time step
nloop = s1.nloop; nstart = 0; ntime = 0

# allocate field data for standard code:
# qe/qi = electron/ion charge density with guard cells
# fxe = smoothed electric field with guard cells
# ffc = form factor array for poisson solver
# mixup = bit reverse table for FFT
# sct = sine/cosine table for FFT
# cue = electron current density with guard cells
# fxyze/byze = smoothed electric/magnetic field with guard cells
# dcu/dcui = electron/ion acceleration density with guard cells
# amu/amui = electron/ion momentum flux with guard cells
# cus = transverse electric field
# ffe = form factor array for iterative poisson solver
sd1.init_dfields13()

# prepare fft tables: updates mixup, sct
mfft1.mfft1_init(s1.mixup,s1.sct,in1.indx)
# calculate form factor: updates ffc
mfield1.mpois1_init(s1.ffc,in1.ax,s1.affp,nx)
# initialize different ensemble of random numbers
if (in1.nextrand > 0):
   minit1.mnextran1(in1.nextrand,in1.ndim,np+npi)

# open reset and restart files; iur, iurr, iur0
s1.open_restart1()

# new start
if (in1.nustrt==1):
# initialize electrons: updates ppart, kpic
# ppart = tiled electron particle arrays
# kpic = number of electrons in each tile
   sb1.init_electrons13()

# initialize background charge density: updates qi
   if (in1.movion==0):
      s1.qi.fill(0.0)
      qmi = -in1.qme
      mpush1.mpost1(s1.ppart,s1.qi,s1.kpic,qmi,s1.tdpost,in1.mx)
      mgard1.maguard1(s1.qi,s1.tguard,nx)

# initialize ions: updates pparti, kipic, cui
# pparti = tiled on particle arrays
# kipic = number of ions in each tile
# cui = ion current density with guard cells
   if (in1.movion==1):
      sb1.init_ions13()

# calculate shift constant for iteration: update wpm, q2m0
   sd1.calc_shift13(iuot)

# initialize darwin electric field
   sd1.cus.fill(0.0)

# restart to continue a run which was interrupted
elif (in1.nustrt==2):
   sd1.bread_drestart13(s1.iur)
   ntime = s1.ntime
   nstart = ntime
# start a new run with data from a previous run
elif (in1.nustrt==0):
   sd1.bread_drestart13(s1.iur0)

# calculate form factor: updates ffe
mfield1.mepois1_init(sd1.ffe,in1.ax,s1.affp,sd1.wpm,in1.ci,nx)

# initialize longitudinal electric field
s1.fxe.fill(0.0)

# set magnitude of external transverse magnetic field
omt = numpy.sqrt(in1.omy*in1.omy + in1.omz*in1.omz)

# reverse simulation at end back to start
if (in1.treverse==1):
   nloop = 2*nloop
   sb1.nloop = nloop
   sd1.nloop = nloop

# initialize all diagnostics from namelist input parameters
# wt = energy time history array=
# pkw = power spectrum for potential
# pkwdi = power spectrum for ion density
# wk = maximum frequency as a function of k for potential
# wkdi = maximum frequency as a function of k for ion density
# fmse/fmsi = electron/ion fluid moments
# fv/fvi = global electron/ion velocity distribution functions
# fvm/fvmi = electron/ion vdrift, vth, entropy for global distribution
# fvtm/fvtmi = time history of electron/ion vdrift, vth, and entropy
# fvtp = velocity distribution function for test particles
# fvmtp = vdrift, vth, and entropy for test particles
# partd = trajectory time history array
# vpkwji = power spectrum for ion current density
# vwkji = maximum frequency as a function of k for ion current
# vpkw = power spectrum for vector potential
# vwk = maximum frequency as a function of k for vector potential
# vpkwet = power spectrum for transverse efield
# vwket = maximum frequency as a function of k for transverse efield
sd1.initialize_ddiagnostics13(ntime)

# read in restart diagnostic file to continue interrupted run
if (in1.nustrt==2):
   sd1.dread_drestart13(s1.iur)

# write reset file
#sd1.bwrite_drestart13(s1.iurr,ntime)

# initialization time
dtimer(dtime,itime,1)
tinit = tinit + float(dtime)
# start timing loop
dtimer(dtime,ltime,-1)

print ("program mdbeps1",file=iuot)

# * * * start main iteration loop * * *

for ntime in range(nstart,nloop):
   print ("ntime = ",ntime,file=iuot)

# debug reset
#  if (ntime==nloop/2):
#     sd1.bread_drestart13(s1.iurr)
#     sd1.reset_ddiags13()

# deposit current with OpenMP: updates cue
   dtimer(dtime,itime,-1)
   sb1.cue.fill(0.0)
   dtimer(dtime,itime,1)
   sb1.tdjpost[0] += float(dtime)
   mcurd1.wmdjpost1(s1.ppart,sb1.cue,s1.kpic,s1.ncl,s1.ihole,in1.qme,
                    zero,in1.ci,sb1.tdjpost,nx,in1.mx,ipbc,
                    in1.relativity,False,irc)
# add guard cells: updates cue
   mgard1.macguard1(sb1.cue,sb1.tguard,nx)

# save electron current for electron current diagnostic later
   if (in1.ndc==0):
      if (in1.ntje > 0):
         it = int(ntime/in1.ntje)
         if (ntime==in1.ntje*it):
            sb1.oldcue[:] = numpy.copy(sb1.cue)

# deposit ion current with OpenMP: updates cui
   if (in1.movion==1):
      dtimer(dtime,itime,-1)
      sb1.cui.fill(0.0)
      dtimer(dtime,itime,1)
      sb1.tdjpost[0] += float(dtime)
      mcurd1.wmdjpost1(s1.pparti,sb1.cui,s1.kipic,s1.ncl,s1.ihole,
                       in1.qmi,zero,in1.ci,sb1.tdjpost,nx,in1.mx,ipbc,
                       in1.relativity,list,irc)
# add guard cells: updates cui
      mgard1.macguard1(sb1.cui,sb1.tguard,nx)

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   s1.qe.fill(0.0)
   dtimer(dtime,itime,1)
   s1.tdpost[0] += float(dtime)
   mpush1.mpost1(s1.ppart,s1.qe,s1.kpic,in1.qme,s1.tdpost,in1.mx)
# add guard cells: updates qe
   mgard1.maguard1(s1.qe,s1.tguard,nx)

# electron density diagnostic: updates sfield=electron density
   if (in1.ntde > 0):
      it = int(ntime/in1.ntde)
      if (ntime==in1.ntde*it):
         s1.edensity_diag1(s1.sfield)
# display smoothed electron density
         graf1.dscaler1(s1.sfield,' EDENSITY',ntime,999,0,nx,irc)
         if (irc[0]==1): break
         irc[0] = 0

# deposit ion charge with OpenMP: updates qi
   if (in1.movion==1):
      dtimer(dtime,itime,-1)
      s1.qi.fill(0.0)
      dtimer(dtime,itime,1)
      s1.tdpost[0] += float(dtime)
      mpush1.mpost1(s1.pparti,s1.qi,s1.kipic,in1.qmi,s1.tdpost,in1.mx)
# add guard cells: updates qi
      mgard1.maguard1(s1.qi,s1.tguard,nx)

# ion density diagnostic: updates sfield=ion density, pkwdi, wkdi
   if (in1.movion==1):
      if (in1.ntdi > 0):
         it = int(ntime/in1.ntdi)
         if (ntime==in1.ntdi*it):
            s1.idensity_diag1(s1.sfield,s1.pkwdi,s1.wkdi,ntime)
            if ((in1.nddi==1) or (in1.nddi==3)):
# display smoothed ion density
               graf1.dscaler1(s1.sfield,' ION DENSITY',ntime,999,1,nx,
                              irc)
               if (irc[0]==1): break
               irc[0] = 0
# ion spectral analysis
            if ((in1.nddi==2) or (in1.nddi==3)):
# display frequency spectrum
               graf1.dmscaler1(s1.wkdi,'ION DENSITY OMEGA VS MODE',
                               ntime,999,1,in1.modesxdi,s1.cwk,irc)
               if (irc[0]==1): break
               irc[0] = 0

# add electron and ion densities: updates qe
   mfield1.maddqei1(s1.qe,s1.qi,s1.tfield,nx)

# add electron and ion current densities: updates cue
   if (in1.movion==1):
      mfield1.maddcuei1(sb1.cue,sb1.cui,sb1.tfield,nx)

# transform charge to fourier space: updates qe
   isign = -1
   mfft1.mfft1r(s1.qe,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# calculate longitudinal force/charge in fourier space:
# updates fxe, we
   mfield1.mpois1(s1.qe,s1.fxe,s1.ffc,s1.we,s1.tfield,nx)

# transform longitudinal electric force to real space: updates fxe
   isign = 1
   mfft1.mfft1r(s1.fxe,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# copy guard cells: updates fxe
   mgard1.mdguard1(s1.fxe,s1.tguard,nx)

# initialize electron deposit data
   dtimer(dtime,itime,-1)
   sd1.dcu.fill(0.0); sd1.amu.fill(0.0)
   dtimer(dtime,itime,1)
   sd1.tdcjpost[0] += float(dtime)
# initialize ion deposit data
   if (in1.movion==1):
      dtimer(dtime,itime,-1)
      sd1.dcui.fill(0.0); sd1.amui.fill(0.0)
      dtimer(dtime,itime,1)
      sd1.tdcjpost[0] += float(dtime)
# predictor for darwin iteration: updates: cue, cus, byze, fxyze
   sd1.darwin_predictor13(sd1.q2m0)

# inner iteration loop
   for k in range(0,in1.ndc):

# initialize electron deposit data
      dtimer(dtime,itime,-1)
      sb1.cue.fill(0.0); sd1.dcu.fill(0.0); sd1.amu.fill(0.0)
      dtimer(dtime,itime,1)
      sd1.tdcjpost[0] += float(dtime)
# initialize ion deposit data
      if (in1.movion==1):
         dtimer(dtime,itime,-1)
         sb1.cui.fill(0.0); sd1.dcui.fill(0.0); sd1.amui.fill(0.0)
         dtimer(dtime,itime,1)
         sd1.tdcjpost[0] += float(dtime)
# updates: dcu, cus, byze, fxyze
      sd1.darwin_iteration(sd1.q2m0,ntime,k)

   pass

# add external traveling wave field
   ts = in1.dt*float(ntime)
   mfield1.meaddext13(sb1.fxyze,sb1.tfield,in1.amodex,in1.freq,ts,
                      in1.trmp,in1.toff,in1.el0,in1.er0,in1.ey0,in1.ez0,
                      nx)

# copy guard cells: updates fxyze
   mgard1.mcguard1(sb1.fxyze,sb1.tguard,nx)

# darwin electron current density diagnostic:
# updates vfield=electron current
   if (in1.ntje > 0):
      it = int(ntime/in1.ntje)
      if (ntime==in1.ntje*it):
         sd1.edcurrent_diag13(sb1.vfield)
# display smoothed electron current
         graf1.dvector1(sb1.vfield,' ELECTRON CURRENT',ntime,999,0,2,nx,
                        irc)
         if (irc[0]==1): break
         irc[0] = 0

# ion current density diagnostic:
# updates vfield=ion current, vpkwji, vwkji
   if (in1.movion==1):
      if (in1.ntji > 0):
         it = int(ntime/in1.ntji)
         if (ntime==in1.ntji*it):
            sb1.icurrent_diag13(sb1.vfield,sb1.vpkwji,sb1.vwkji,ntime)
            if ((in1.ndji==1) or (in1.ndji==3)):
# display smoothed ion current
               graf1.dvector1(sb1.vfield,' ION CURRENT',ntime,999,0,2,
                              nx,irc)
               if (irc[0]==1): break
               irc[0] = 0
# ion spectral analysis
            if ((in1.ndji==2) or (in1.ndji==3)):
# display frequency spectrum
               graf1.dmvector1(sb1.vwkji,'ION CURRENT OMEGA VS MODE',
                               ntime,999,2,2,in1.modesxji,s1.cwk,irc)
               if (irc[0]==1): break
               irc[0] = 0

# potential diagnostic: updates sfield=potential, pkw, wk
   if (in1.ntp > 0):
      it = int(ntime/in1.ntp)
      if (ntime==in1.ntp*it):
         s1.potential_diag1(s1.sfield,s1.pkw,s1.wk,ntime)
         if ((in1.ndp==1) or (in1.ndp==3)):
# display potential
            graf1.dscaler1(s1.sfield,' POTENTIAL',ntime,999,0,nx,irc)
            if (irc[0]==1): break
            irc[0] = 0
# spectral analysis
         if ((in1.ndp==2) or (in1.ndp==3)):
# display frequency spectrum
            graf1.dmscaler1(s1.wk,'POTENTIAL OMEGA VS MODE',ntime,999,2,
                            in1.modesxp,s1.cwk,irc)
            if (irc[0]==1): break
            irc[0] = 0

# longitudinal efield diagnostic: updates sfield=longitudinal efield
   if (in1.ntel > 0):
      it = int(ntime/in1.ntel)
      if (ntime==in1.ntel*it):
         s1.elfield_diag1(s1.sfield)
# display longitudinal efield
         graf1.dscaler1(s1.sfield,' ELFIELD',ntime,999,0,nx,irc)
         if (irc[0]==1): break
         irc[0] = 0

# vector potential diagnostic:updates vfield=vector potential, vpkw, vwk
   if (in1.nta > 0):
      it = int(ntime/in1.nta)
      if (ntime==in1.nta*it):
         sd1.vdpotential_diag13(sd1.vfield,sd1.vpkw,sd1.vwk,ntime)
         if ((in1.nda==1) or (in1.nda==3)):
# display vector potential
            graf1.dvector1(sd1.vfield,' VECTOR POTENTIAL',ntime,999,0,2,
                           nx,irc)
            if (irc[0]==1): break
            irc[0] = 0
# spectral analysis
         if ((in1.nda==2) or (in1.nda==3)):
# display frequency spectrum
            graf1.dmvector1(sd1.vwk,'VECTOR POTENTIAL OMEGA VS MODE',
                            ntime,999,2,2,in1.modesxa,s1.cwk,irc)
            if (irc[0]==1): break
            irc[0] = 0

# transverse efield diagnostic:
# updates vfield=transverse efield, vpkwet, vwket
   if (in1.ntet > 0):
      it = int(ntime/in1.ntet)
      if (ntime==in1.ntet*it):
         sd1.detfield_diag13(sd1.vfield,sd1.vpkwet,sd1.vwket,ntime)
         if ((in1.ndet==1) or (in1.ndet==3)):
# display transverse efield
            graf1.dvector1(sd1.vfield,' TRANSVERSE EFIELD',ntime,999,0,
                           2,nx,irc)
            if (irc[0]==1): break
            irc[0] = 0
# spectral analysis
         if ((in1.ndet==2) or (in1.ndet==3)):
# display frequency spectrum
            graf1.dmvector1(sd1.vwket,'TRANSVERSE EFIELD OMEGA VS MODE',
                            ntime,999,2,2,in1.modesxet,s1.cwk,irc)
            if (irc[0]==1): break
            irc[0] = 0

# magnetic field diagnostic: updates vfield=bfield
   if (in1.ntb > 0):
      it = int(ntime/in1.ntb)
      if (ntime==in1.ntb*it):
         sd1.dbfield_diag13(sd1.vfield)
# display magnetic field
         graf1.dvector1(sd1.vfield,' MAGNETIC FIELD',ntime,999,0,2,nx,
                        irc)
         if (irc[0]==1): break
         irc[0] = 0

# fluid moments diagnostic
   if (in1.ntfm > 0):
      it = int(ntime/in1.ntfm)
      if (ntime==in1.ntfm*it):
# updates fmse
         sd1.edfluidms_diag13(s1.fmse)
         if (in1.movion==1):
# updates fmsi
            sd1.idfluidms_diag13(s1.fmsi)

# velocity diagnostic
   if (in1.ntv > 0):
      it = int(ntime/in1.ntv)
      if (ntime==in1.ntv*it):
# updates fv, fe, fvm, fvtm, wkt
         sb1.evelocity_diag13(s1.fv,s1.fe,s1.fvm,s1.fvtm,s1.wkt)
# display electron velocity distributions
         if ((in1.ndv==1) or (in1.ndv==3)):
            if ((in1.nvft==1) or (in1.nvft==3)):
               graf1.displayfv1(s1.fv,s1.fvm,' ELECTRON',ntime,in1.nmv,
                                2,irc)
               if (irc[0]==1): break
               irc[0] = 0
# display electron velocity distributions in cylindrical co-ordinates
            if ((in1.nvft==4) or (in1.nvft==5)):
               graf1.displayfvb1(s1.fv,s1.fvm,' ELECTRON',ntime,in1.nmv,
                                2,irc)  
               if (irc[0]==1): break
               irc[0] = 0   
# display electron energy distribution
            if ((in1.nvft==2) or (in1.nvft==3) or (in1.nvft==5)):
               graf1.displayfe1(s1.fe,s1.wkt,' ELECTRON',ntime,in1.nmv,
                                irc)
               if (irc[0]==1): break
               irc[0] = 0  
# ion distribution function
         if (in1.movion==1):
# updates fvi, fei, fvmi, fvtmi, wkt
            sb1.ivelocity_diag13(s1.fvi,s1.fei,s1.fvmi,s1.fvtmi,s1.wkt)
# display ion velocity distributions
            if ((in1.ndv==2) or (in1.ndv==3)):
               if ((in1.nvft==1) or (in1.nvft==3)):
                  graf1.displayfv1(s1.fvi,s1.fvmi,' ION',ntime,in1.nmv,
                                   2,irc)
                  if (irc[0]==1): break
                  irc[0] = 0
# display ion velocity distributions in cylindrical co-ordinates
               if ((in1.nvft==4) or (in1.nvft==5)):   
                  graf1.displayfvb1(s1.fvi,s1.fvmi,' ION',ntime,in1.nmv,
                                    2,irc) 
                  if (irc[0]==1): break
                  irc[0] = 0
# display ion energy distribution
               if ((in1.nvft==2) or (in1.nvft==3) or (in1.nvft==5)):
                  ts = s1.fei[2*in1.nmv+1,0]
                  s1.fei[2*in1.nmv+1,0] = in1.rmass*ts
                  graf1.displayfe1(s1.fei,s1.wkt,' ION',ntime,in1.nmv,
                                   irc)
                  if (irc[0]==1): break
                  irc[0] = 0
                  s1.fei[2*in1.nmv+1,0] = ts

# trajectory diagnostic: updates partd, fvtp, fvmtp
   if (in1.ntt > 0):
      it = int(ntime/in1.ntt)
      if (ntime==in1.ntt*it):
         sb1.traj_diag13(s1.partd,s1.fvtp,s1.fvmtp)
         if (in1.nst==3):
# display test particle velocity distributions
            if (in1.ndt==1):
               graf1.displayfv1(s1.fvtp,s1.fvmtp,' ELECTRON',ntime,
                                in1.nmv,2,irc)
            elif (in1.ndt==2):
               graf1.displayfv1(s1.fvtp,s1.fvmtp,' ION',ntime,in1.nmv,
                                2,irc)
               if (irc[0]==1): break
               irc[0] = 0

# phase space diagnostic
   if (in1.nts > 0):
      it = int(ntime/in1.nts)
      if (ntime==in1.nts*it):
# calculate electron phase space distribution: updates fvs
         s1.ephasesp_diag1(s1.fvs)
# plot electrons
         if ((in1.nds==1) or (in1.nds==3)):
# vx, vy, or vz versus x
            nn = in1.nsxv; ierr[0] = 0
            for i in range(0,3):
               if ((nn % 2)==1):
                  graf1.dpmbgrasp1(s1.ppart,s1.kpic,' ELECTRON',ntime,
                                   999,in1.omx,in1.omy,in1.omz,nx,i+2,1,
                                   in1.ntsc,irc)
                  if (irc[0]==1):
                     ierr[0] = 1
                     break
                  irc[0] = 0
               nn = int(nn/2)
            if (ierr[0]==1): break
# vx-vy, vx-vz or vy-vz
            nn = in1.nsvv; ierr[0] = 0
            for i in range(0,3):
               if ((nn % 2)==1):
                  graf1.dpmbgrasp1(s1.ppart,s1.kpic,' ELECTRON',ntime,
                                   999,in1.omx,in1.omy,in1.omz,nx,
                                   min(i+3,4),max(i+1,2),in1.ntsc,irc)
                  if (irc[0]==1):
                     ierr[0] = 1
                     break
                  irc[0] = 0
               nn = int(nn/2)
            if (ierr[0]==1): break
# ion phase space
         if (in1.movion==1):
# calculate ion phase space distribution: updates fvsi
            s1.iphasesp_diag1(s1.fvsi)
# plot ions
            if ((in1.nds==2) or (in1.nds==3)):
# vx, vy, or vz versus x
               nn = in1.nsxv; ierr[0] = 0
               for i in range(0,3):
                  if ((nn % 2)==1):
                     graf1.dpmbgrasp1(s1.pparti,s1.kipic,' ION',ntime,
                                      999,in1.omx,in1.omy,in1.omz,nx,
                                      i+2,1,in1.ntsc,irc)
                     if (irc[0]==1):
                        ierr[0] = 1
                        break
                     irc[0] = 0
                  nn = int(nn/2)
               if (ierr[0]==1): break
# vx-vy, vx-vz or vy-vz
               nn = in1.nsvv; ierr[0] = 0
               for i in range(0,3):
                  if ((nn % 2)==1):
                     graf1.dpmbgrasp1(s1.pparti,s1.kipic,' ION',ntime,
                                      999,in1.omx,in1.omy,in1.omz,nx,
                                      min(i+3,4),max(i+1,2),in1.ntsc,
                                      irc)
                     if (irc[0]==1):
                        ierr[0] = 1
                        break
                     irc[0] = 0
                  nn = int(nn/2)
               if (ierr[0]==1): break

# push electrons with OpenMP: updates ppart, wke, kpic
   sd1.dpush_electrons13(s1.ppart,s1.kpic)

# push ions with OpenMP: updates pparti, wki, kipic
   if (in1.movion==1):
      sd1.dpush_ions13(s1.pparti,s1.kipic)

# start running simulation backwards:
# need to reverse time lag in leap-frog integration scheme
   if (in1.treverse==1):
      if (((ntime+1)==(nloop/2)) or ((ntime+1)==nloop)):
         sd1.d_time_reverse1()

# energy diagnostic
   if (in1.ntw > 0):
      it = int(ntime/in1.ntw)
      if (ntime==in1.ntw*it):
         sd1.denergy_diag13(s1.wt,ntime,iuot)

# restart file
   if (in1.ntr > 0):
      n = ntime + 1
      it = int(n/in1.ntr)
      if (n==in1.ntr*it):
         dtimer(dtime,itime,-1)
         sd1.bwrite_drestart13(s1.iur,n)
         sd1.dwrite_drestart13(s1.iur)
         in1.writnml1(s1.iudm)
         dtimer(dtime,itime,1)
         s1.tfield[0] += float(dtime)

ntime = ntime + 1

# loop time
dtimer(dtime,ltime,1)
tloop = tloop + float(dtime)

# * * * end main iteration loop * * *

print (file=iuot)
print ("ntime,relativity,ndc=",ntime,",",in1.relativity,",",in1.ndc,
       file=iuot)
if (in1.treverse==1):
   print ("treverse = ",in1.treverse,file=iuot)

# print timing summaries
sd1.print_dtimings13(tinit,tloop,iuot)

if ((in1.ntw > 0) or (in1.ntt > 0)):
   graf1.reset_graphs()

# trajectory diagnostic
if (in1.ntt > 0):
   if ((in1.nst==1) or (in1.nst==2)):
      if (in1.nplot > 0):
         irc[0] = graf1.open_graphs(1)
      ts = in1.t0
      graf1.displaytr1(s1.partd,ts,in1.dt*float(in1.ntt),s1.itt,2,999,
                       irc)
      if (irc[0]==1):
         exit(0)
      graf1.reset_nplot(in1.nplot,irc)

# energy diagnostic
if (in1.ntw > 0):
   ts = in1.t0
# display energy histories
   graf1.displayw1(s1.wt,ts,in1.dt*float(in1.ntw),s1.itw,irc)
   if (irc[0]==1):
      exit(0)
# print energy summaries
   sb1.print_energy13(s1.wt,iuot)

# velocity diagnostic
if (in1.ntv > 0):
   ts = in1.t0
# display electron distribution time histories and entropy
   if ((in1.ndv==1) or (in1.ndv==3)):
      graf1.displayfvt1(s1.fvtm,' ELECTRON',ts,in1.dt*float(in1.ntv),
                        s1.itv,irc)
      if (irc[0]==1):
         exit(0)
# display ion distribution time histories and entropy
   if (in1.movion==1):
      if ((in1.ndv==2) or (in1.ndv==3)):
         graf1.displayfvt1(s1.fvtmi,' ION',ts,in1.dt*float(in1.ntv),
                           s1.itv,irc)
         if (irc[0]==1):
            exit(0)

# display final spectral analysis for ion density
if (in1.movion==1):
   if (in1.ntdi > 0):
      if ((in1.nddi==2) or (in1.nddi==3)):
# display frequency spectrum
         graf1.dmscaler1(s1.wkdi,'ION DENSITY OMEGA VS MODE',ntime,999,
                         1,in1.modesxdi,s1.cwk,irc)
         if (irc[0]==1):
            exit(0)

# display final spectral analysis for potential
if (in1.ntp > 0):
   if ((in1.ndp==2) or (in1.ndp==3)):
# display frequency spectrum
      graf1.dmscaler1(s1.wk,'POTENTIAL OMEGA VS MODE',ntime,999,2,
                      in1.modesxp,s1.cwk,irc)
      if (irc[0]==1):
         exit(0)

# display final spectral analysis for ion current density
if (in1.movion==1):
   if (in1.ntji > 0):
      if ((in1.ndji==2) or (in1.ndji==3)):
# display frequency spectrum
         graf1.dmvector1(sb1.vwkji,'ION CURRENT OMEGA VS MODE',ntime,
                         999,2,2,in1.modesxji,s1.cwk,irc)
         if (irc[0]==1):
            exit(0)

# display final spectral analysis for vector potential
if (in1.nta > 0):
   if ((in1.nda==2) or (in1.nda==3)):
# display frequency spectrum
      graf1.dmvector1(sd1.vwk,'VECTOR POTENTIAL OMEGA VS MODE',ntime,
                      999,2,2,in1.modesxa,s1.cwk,irc)
      if (irc[0]==1):
         exit(0)

# display final spectral analysis for transverse efield
if (in1.ntet > 0):
   if ((in1.ndet==2) or (in1.ndet==3)):
# display frequency spectrum
      graf1.dmvector1(sd1.vwket,'TRANSVERSE EFIELD OMEGA VS MODE',ntime,
                      999,2,2,in1.modesxet,s1.cwk,irc)
      if (irc[0]==1):
         exit(0)

# close diagnostics
sd1.close_ddiags13(s1.iudm)
# close reset and restart files: iur, iurr, iur0
s1.close_restart1()
# close output file
print (" * * * q.e.d. * * *",file=iuot)
iuot.close()
# close graphics device
graf1.close_graphs()
