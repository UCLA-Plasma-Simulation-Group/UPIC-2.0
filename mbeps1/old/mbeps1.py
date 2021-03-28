#-----------------------------------------------------------------------
# 1D Electrostatic OpenMP PIC code
# written by Viktor K. Decyk and Joshua Kelly, UCLA
# copyright 2016, regents of the university of california
from __future__ import print_function
import sys
import math
import numpy

sys.path.append('./mbeps1.source')
from libmpush1 import *
from fomplib import *
from fgraf1 import *
from dtimer import *

# read namelist
iuin = 8
in1.readnml1(iuin)
# override input data
in1.idcode = 1
in1.ndim = 1
if (in1.nts > 0):
   in1.nsxv = min(in1.nsxv,1); in1.nsvv = 0

# import electrostatic module after namelist has been read
import s1

int_type = numpy.int32
double_type = numpy.float64
#if (minit1.fprecision()==0):
float_type = numpy.float32
complex_type = numpy.complex64
#else:
#  float_type = numpy.float64
#  complex_type = numpy.complex128
#  print ("using double precision")

# declare scalars for standard code
npi = 0
ws = numpy.zeros((1),float_type)

# declare scalars for OpenMP code
irc = numpy.zeros((1),int_type)

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
#in1.nvp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
omplib.init_omp(in1.nvp)

# open graphics device
irc[0] = graf1.open_graphs(in1.nplot)

# initialize scalars for standard code
# np = total number of electrons in simulation
np = s1.np
# nx = number of grid points in x direction
nx = s1.nx
# npi = total number of ions in simulation
if (in1.movion > 0):
   npi = s1.npi
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
s1.init_fields1()

# prepare fft tables: updates mixup, sct
mfft1.mfft1_init(s1.mixup,s1.sct,in1.indx)
# calculate form factors: updates ffc
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
   s1.init_electrons1()

# initialize background charge density: updates qi
   if (in1.movion==0):
      s1.qi.fill(0.0)
      if (in1.ndprof > 0):
         qmi = -in1.qme
         mpush1.mpost1(s1.ppart,s1.qi,s1.kpic,qmi,s1.tdpost,in1.mx)
         mgard1.maguard1(s1.qi,s1.tguard,nx)

# initialize ions: updates pparti, kipic
# pparti = tiled on particle arrays
# kipic = number of ions in each tile
# cui = ion current density with guard cells
   if (in1.movion==1):
      s1.init_ions1()

# restart to continue a run which was interrupted
elif (in1.nustrt==2):
      s1.bread_restart1(s1.iur)
      ntime = s1.ntime
      nstart = ntime
# start a new run with data from a previous run
elif (in1.nustrt==0):
      s1.bread_restart1(s1.iur0)
      ntime = s1.ntime

# reverse simulation at end back to start
if (in1.treverse==1):
   nloop = 2*nloop
   s1.nloop = nloop

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
s1.initialize_diagnostics1(ntime)

# read in restart diagnostic file to continue interrupted run
if (in1.nustrt==2):
   s1.dread_restart1(s1.iur)

# write reset file
#s1.bwrite_restart1(s1.iurr,ntime)

# initialization time
dtimer(dtime,itime,1)
tinit = tinit + float(dtime)
# start timing loop
dtimer(dtime,ltime,-1)

print ("program mbeps1",file=iuot)

# * * * start main iteration loop * * *

for ntime in range(nstart,nloop):
   print ("ntime = ",ntime,file=iuot)

# debug reset
#  if (ntime==int(nloop/2)):
#     s1.bread_restart1(s1.iurr)
#     s1.reset_diags1()

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

# transform charge to fourier space: updates qe
   isign = -1
   mfft1.mfft1r(s1.qe,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# calculate force/charge in fourier space: updates fxe, we
   mfield1.mpois1(s1.qe,s1.fxe,s1.ffc,s1.we,s1.tfield,nx)

# transform force to real space: updates fxe
   isign = 1
   mfft1.mfft1r(s1.fxe,isign,s1.mixup,s1.sct,s1.tfft,in1.indx)

# add external traveling wave field
   ts = in1.dt*float(ntime)
   mfield1.meaddext1(s1.fxe,s1.tfield,in1.amodex,in1.freq,ts,in1.trmp,
                     in1.toff,in1.el0,in1.er0,nx)

# copy guard cells: updates fxe
   mgard1.mdguard1(s1.fxe,s1.tguard,nx)

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

# fluid moments diagnostic
   if (in1.ntfm > 0):
      it = int(ntime/in1.ntfm)
      if (ntime==in1.ntfm*it):
# updates fmse
         s1.efluidms_diag1(s1.fmse)
         if (in1.movion==1):
# updates fmsi
            s1.ifluidms_diag1(s1.fmsi)

# velocity diagnostic
   if (in1.ntv > 0):
      it = int(ntime/in1.ntv)
      if (ntime==in1.ntv*it):
# updates fv, fvm, fvtm
         s1.evelocity_diag1(s1.fv,s1.fe,s1.fvm,s1.fvtm,s1.wkt)
# display electron velocity distributions
         if ((in1.ndv==1) or (in1.ndv==3)):
            if ((in1.nvft==1) or (in1.nvft==3)):
               graf1.displayfv1(s1.fv,s1.fvm,' ELECTRON',ntime,in1.nmv,
                                1,irc)
               if (irc[0]==1): break
               irc[0] = 0
# display electron energy distribution
            if ((in1.nvft==2) or (in1.nvft==3)):
               graf1.displayfe1(s1.fe,s1.wkt,' ELECTRON',ntime,in1.nmv,
                                irc)
               if (irc[0]==1): break
               irc[0] = 0
         if (in1.movion==1):
# updates fvi, fvmi, fvtmi
            s1.ivelocity_diag1(s1.fvi,s1.fei,s1.fvmi,s1.fvtmi,s1.wkt)
# display ion velocity distributions
            if ((in1.ndv==2) or (in1.ndv==3)):
               if ((in1.nvft==1) or (in1.nvft==3)):
                  graf1.displayfv1(s1.fvi,s1.fvmi,' ION',ntime,in1.nmv,
                                   1,irc)
                  if (irc[0]==1): break
                  irc[0] = 0
# display ion energy distribution
               if ((in1.nvft==2) or (in1.nvft==3)):
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
         s1.traj_diag1(s1.partd,s1.fvtp,s1.fvmtp)
         if (in1.nst==3):
# display test particle velocity distributions
            if (in1.ndt==1):
               graf1.displayfv1(s1.fvtp,s1.fvmtp,' ELECTRON',ntime,
                                in1.nmv,1,irc)
            elif (in1.ndt==2):
               graf1.displayfv1(s1.fvtp,s1.fvmtp,' ION',ntime,in1.nmv,1,
                                irc)
            if (irc[0]==1): break
            irc[0] = 0

# phase space diagnostic
   if (in1.nts > 0):
      it = int(ntime/in1.nts)
      if (ntime==in1.nts*it):
# calculate electron phase space distribution: updates fvs
         s1.ephasesp_diag1(s1.fvs)
# plot electrons vx versus x
         if ((in1.nds==1) or (in1.nds==3)):
            graf1.dpmgrasp1(s1.ppart,s1.kpic,' ELECTRON',ntime,999,nx,2,
                            1,in1.ntsc,irc)
            if (irc[0]==1): break
            irc[0] = 0
# ion phase space
         if (in1.movion==1):
# calculate ion phase space distribution: updates fvsi
            s1.iphasesp_diag1(s1.fvsi)
# plot ions vx versus x
            if ((in1.nds==2) or (in1.nds==3)):
               graf1.dpmgrasp1(s1.pparti,s1.kipic,' ION',ntime,999,nx,2,
                               1,in1.ntsc,irc)
               if (irc[0]==1): break
               irc[0] = 0
         
# push electrons with OpenMP: updates ppart, wke, kpic
   s1.push_electrons1(s1.ppart,s1.kpic)

# push ions with OpenMP: updates pparti, wki, kipic
   if (in1.movion==1):
      s1.push_ions1(s1.pparti,s1.kipic)

# start running simulation backwards:
# need to reverse time lag in leap-frog integration scheme
   if (in1.treverse==1):
      if (((ntime+1)==int(nloop/2)) or ((ntime+1)==nloop)):
         s1.es_time_reverse1()

# energy diagnostic: updates wt
   if (in1.ntw > 0):
      it = int(ntime/in1.ntw)
      if (ntime==in1.ntw*it):
         s1.energy_diag1(s1.wt,ntime,iuot)

# restart file
   if (in1.ntr > 0):
      n = ntime + 1
      it = int(n/in1.ntr)
      if (n==in1.ntr*it):
         dtimer(dtime,itime,-1)
         s1.bwrite_restart1(s1.iur,n)
         s1.dwrite_restart1(s1.iur)
         in1.writnml1(s1.iudm)
         dtimer(dtime,itime,1)
         s1.tfield[0] += float(dtime)

ntime = ntime + 1

# loop time
dtimer(dtime,ltime,1)
tloop = tloop + float(dtime)

# * * * end main iteration loop * * *

print (file=iuot)
print ("ntime, relativity = ",ntime,",",in1.relativity,file=iuot)
if (in1.treverse==1):
   print ("treverse = ",in1.treverse,file=iuot)

# print timing summaries
s1.print_timings1(tinit,tloop,iuot)

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
   s1.print_energy1(s1.wt,iuot)

# velocity diagnostic
if (in1.ntv > 0):
   ts = in1.t0
# display electron distribution time histories and entropy
   if ((in1.ndv==1) or (in1.ndv==3)):
      graf1.displayfvt1(s1.fvtm,' ELECT',ts,in1.dt*float(in1.ntv),
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

# close diagnostics
s1.close_diags1(s1.iudm)
# close reset and restart files: iur, iurr, iur0
s1.close_restart1()
# close output file
print (" * * * q.e.d. * * *",file=iuot)
iuot.close()
# close graphics device
graf1.close_graphs()
