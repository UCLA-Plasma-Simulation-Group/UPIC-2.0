!-----------------------------------------------------------------------
! Input for 1-2/2D MPI/OpenMP PIC codes
!
      module in1
!
! input1mod.f defines namelists containing input and output variables:
! readnml1 reads namelist from unit iuin
! writnml1 writes final diagnostic metafile to unit iudm
! written by viktor k. decyk, ucla
! copyright 2011, regents of the university of california
! update: January 31, 2021
!
      implicit none
!
! Basic Input Namelist
      save
! Identification Parameters:
! idrun = run identifier for current run
! idcode = code identifier
      integer :: idrun = 1, idcode = 0
!
! Global Simulation Parameters:
! indx = exponent which determines length in x direction, nx=2**indx
      integer :: indx =  11
! psolve = type of poisson solver = (1,2,3)
!     integer :: psolve = PERIODIC_2D
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: relativity = 0
! ci = recipr0cal of velocity of light
      real :: ci = 0.1
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! ndim = number of velocity dimensions = 1 or 3
      integer :: ndim = 3
! nvdist = velocity distribution type
! nvdist = (1,2,3) = (maxwellian/juttner,waterbag,ring) distribution
! for nvdist=2, maximum velocity in x/y/z is (vtx/vty/vtz)*sqrt(3)
! for nvdist=3, x component of thermal velocity (vtx) and drift (vdx)
! are used to set radial thermal velocity (vtr)  and ring radius (vdr)
      integer :: nvdist = 1
! treverse = (0,1) = (no,yes) reverse simulation at end back to start
      integer :: treverse = 0
!
! Background Electron Parameters:
! npx = number of background electrons distributed in x direction
      integer :: npx = 409600
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
!
! Beam Electron Parameters:
! npxb = number of beam electrons in x direction
      integer :: npxb = 0
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
!
! Time Parameters:
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
      real :: tend = 45.000, dt = 0.1
!
! Numerical Parameters:
! inorder = interpolation order
! popt = particle optimization scheme
! dopt = charge deposit optimization scheme
! djopt = current deposit optimization scheme
!     integer :: inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD
!     integer :: djopt = STANDARD
! ax = half-width of particle in x direction
!     real :: ax = .816497
!     real :: ax = .866025
      real :: ax = .912871
! mx = number of grids in x in sorting tiles
      integer :: mx = 32
! nextrand = (0,N) = generate (default,Nth block) of random numbers
      integer :: nextrand = 0
!
! Initial Electron Density Parameters:
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
! for ndprof = 0, n(x) = n0
! for ndprof = 1, n(x) = n0*(1 + ampdx*(x/nx - 0.5))
! for ndprof = 2, n(x) = n0*(1 + ampdx*sin(x/scaledx - shiftdx))
! for ndprof = (3,4,5), n(x) = n0*(1 + ampdx*f((x - shiftdx)/scaledx))
! where f = (exp(-x**2/2),sech(x)**2,exp(x))
! n0 is determined by the requirement that the integral over density
! equals the total number of particles distributed by the function
      integer :: ndprof = 0
! ampdx = amplitude of density compared to uniform in x
! scaledx = scale length for spatial coordinate in x
! shiftdx = shift of spatial coordinate in x
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
!
! Zero Force Parameter:
! mzf = (0,1) = (no,yes) set forces to zero
      integer :: mzf = 0
!
! External Electrostatic Traveling Wave Driver:
! Extx(x) = e1*sin(k0*x + freq*time) + e2*cos(k0*x - freq*time)
! amodex = wave number which determines k0 = 2*pi*amodex/NX
! freq = frequency of external wave driver
! trmp = ramp-up time for external wave driver
! toff = shut-off time for external wave driver
! if toff < 0 => toff = + infinity
      real :: amodex = 0.0, freq = 0.0, trmp = 0.0, toff = 0.0
! el0/er0 = external pump amplitude for left-going/right-going wave
!     e1 = el0*(time/trmp), e2 = er0*(time/trmp), if time < trmp
!     e1 = el0,             e2 = er0,             if trmp < time < toff
!     e1 = 0,               e2 = 0,               if time > toff
      real :: el0 = 0.0, er0 = 0.0
!
! External Electromagnetic Circularly Polarized Wave Driver:
! Exty(x) = e3*cos(k0*x - freq*time)
! Extz(x) = e4*sin(k0*x - freq*time)
! ey0/ez0 = external pump amplitude for y/z circularly polarized wave
!     e3 = ey0*(time/trmp), e4 = ez0*(time/trmp), if time < trmp
!     e3 = ey0,             e4 = ez0,             if trmp < time < toff
!     e3 = 0,               e4 = 0,               if time > toff
      real :: ey0 = 0.0, ez0 = 0.0
!
! Restart Parameters:
! nustrt = type of initialization
! 0 = start a new run with data from a previous run
! 1 = start a new run from random numbers
! 2 = restart, continue a run which was interrupted
! ntr = number of time steps between restart routine
! idrun0 = run identifier for old run, used by nustrt = 0
!          if idrun0 = 0, then idrun0 is set to idrun
      integer :: nustrt = 1, ntr = 0, idrun0 = 0
!
! Energy Diagnostic Parameters:
! ntw = number of time steps between energy diagnostic
! ndw = (0,1) = (no,yes) = print energy values in output file
      integer :: ntw = 1, ndw = 1
!
! Electron Density Diagnostic Parameters:
! ntde = number of time steps between electron density diagnostic
! modesxde = number of modes in x to keep for electron density
!            diagnostic
! nderec = current record number for electron density writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntde = 0, modesxde = 41, nderec = -1
!
! Potential Diagnostic Parameters:
! ntp = number of time steps between potential diagnostic
! ndp = (0,1,2,3) = display (nothing,potential,spectrum,both)
! modesxp = number of modes in x to keep for potential diagnostic
! nprec = current record number for potential writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: ntp = 0, ndp = 1, modesxp = 41, nprec = -1
!
! Longitudinal Efield Diagnostic Parameters:
! ntel = number of time steps between longitudinal efield diagnostic
! modesxel = number of modes in x to keep for longitudinal efield
!            diagnostic
! nelrec = current record number for longitudinal efield writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntel = 0, modesxel = 41, nelrec = -1
!
! Power Spectrum Diagnostic Parameters:
! wmin/wmax = minimum/maximum frequency used in power spectrum
! dw = frequency increment used in power spectrum
      real :: wmin = 0.0, wmax = 2.0, dw = 0.01
!
! Fluid Moments Diagnostic Parameter:
! ntfm = number of time steps between fluid moments diagnostic
! ndfm = (0,1,2,3) = display (nothing,electrons,ions,both)
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy flux)
! if npro = n is selected, all profiles less than n are also calculated
      integer :: ntfm = 0, ndfm = 1, npro = 2
! nferec = current record number for electron fluid moments writes
! nfirec = current record number for ion fluid moments writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: nferec = -1, nfirec = -1
!
! Velocity-Space Diagnostic Parameter:
! ntv = number of time steps between velocity-space diagnostic
! ndv = (0,1,2,3) = display (nothing,electrons,ions,both)
! nmv = number of segments in v for velocity distribution
      integer :: ntv = 0, ndv = 3, nmv = 40
! nvft = (1,2,3,4,5) = (cartesian,energy,cartesian+energy,cylindrical,
!        (cylindrical+energy) 1d distribution functions
!        for cylindrical, z axis is along the external magnetic field
      integer :: nvft = 1
! nverec = current record number for electron velocity distribution writes
! nvirec = current record number for ion velocity distribution writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: nverec = -1, nvirec = -1
!
! Trajectory Diagnostic Parameters:
! ntt = number of time steps between trajectory diagnostic.
! ndt = (0,1,2) = process (nothing,electrons,ions)
! nst = type of test particle distribution, if ntt > 0
! 1 = uniformly distribution in real space
! 2 = uniform distribution in velocity space
! 3 = velocity slice at vtsx +- dvtx/2
! nprobt = number of test charges whose trajectories will be stored.
      integer :: ntt = 0, ndt = 1, nst = 0, nprobt = 0
! vtsx = center of velocity slice if nst = 3
! dvtx = width of velocity slice if nst = 3
      real :: vtsx = 0.0, dvtx = 0.1
! ntrec = current record number for trajectory writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntrec = -1
!
! Phase-Space Diagnostic
! nts = number of time steps between phase space diagnostic
! nds = (0,1,2,3) = display (nothing,electrons,ions,both)
! nsxv = component(s) for phase-space display(s), if nts > 0
! 1 = x-vx, 2 = x-vy, 3 = x-vx and x-vy, 4 = x-vz, 5 = x-vx and x-vz,
! 6 = x-vy and x-vz, 7 = x-vx and x-vy and x-vz
! nsvv = component(s) for phase-space display(s), if nts > 0
! 1 = vx-vy, 2 = vx-vz, 3 = vx-vy and vx-vz, 4 = vy-vz,
! 5 = vx-vy and vy-vz, 6 = vx-vy, vx-vz, and vy-vz
! if the magnetic field is at an angle to the cartesian axes, then
! vx means vperp1, vy means vperp2, and vz means vparallel
! ntsc = (0,1) = (no,yes) color beam particles
      integer :: nts = 0, nds = 3, nsxv = 1, nsvv = 0, ntsc = 0
! mvx/mvy = number of grids in x/y for phase space aggregation
      integer :: mvx = 3, mvy = 3
! nserec = current record number for electron phase space writes
! nsirec = current record number for ion phase space writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: nserec = -1, nsirec = -1
!
! ntm = number of time steps between momentum diagnostic
!     integer :: ntm = 0
!
! Ion Parameter:
! movion = (0,1) = (no,1) number of moving ions
      integer :: movion = 0
!
! Electromagnetic Field Parameter:
! emf = (0,1,2) = use (electrostatic,electromagnetic,darwin) fields
      integer :: emf = 0
!
! Display Parameter:
! nplot = maximum number of plots per page
      integer :: nplot = 4
!
! Multi-tasking Parameter:
! nvp = number of shared memory nodes (0=default)
      integer :: nvp = 0
!
! Error Processing Parameter
! monitor = (0,1,2) = (disable,normal,extended) error processing
      integer :: monitor = 0
!
! define namelist
      namelist /input1/ idrun, idcode, indx, mx, npx, npxb, qme, vtx,   &
     &vty, vtz, vx0, vy0, vz0, vdx, vdy, vdz, vtdx, vtdy, vtdz,         &
     &relativity, ci, xtras, ndim, nvdist, treverse, tend, dt, ax,      &
     &nextrand, mzf, ndprof, ampdx, scaledx, shiftdx, amodex, freq,     &
     &trmp, toff, el0, er0, ey0, ez0, ntw, ndw, ntde, modesxde, nderec, &
     &ntp, ndp, modesxp, nprec, ntel, modesxel, nelrec, wmin, wmax, dw, &
     &ntfm, ndfm, npro, nferec, nfirec, ntv, ndv, nmv, nvft, nverec,    &
     &nvirec, ntt, ndt, nst, nprobt, vtsx, dvtx, ntrec, nts, nds, nsxv, &
     &nsvv, ntsc, mvx, mvy, nserec, nsirec, movion, emf, nustrt, ntr,   &
     &idrun0, nplot, nvp, monitor
!
! Electromagnetic Namelist
! External Magnetic Field Parameters:
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
!
! Electron Current Diagnostic Parameters:
! ntje = number of time steps between ion current diagnostic
! modesxje = number of modes in x to keep for ion current diagnostic
! njerec = current record number for ion current writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntje = 0, modesxje = 41, njerec = -1
!
! Vector Potential Diagnostic Parameters:
! nta = number of time steps between vector potential diagnostic
! nda = (0,1,2,3) = display (nothing,vector potential,spectrum,both)
! modesxa = number of modes in x to keep for vector potential diagnostic
! narec = current record number for vector potential writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: nta = 0, nda = 1, modesxa = 41, narec = -1
!
! Transverse Efield Diagnostic Parameters:
! ntet = number of time steps between transverse efield diagnostic
! ndet = (0,1,2,3) = display (nothing,transverse efield,spectrum,both)
! modesxet = number of modes in x to keep for transverse efield
!            diagnostic
! netrec = current record number for ltransverse efield writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntet = 0, ndet = 1, modesxet = 41, netrec = -1
!
! Magnetic Field Diagnostic Parameters:
! ntb = number of time steps between magnetic field diagnostic
! modesxb = number of modes in x to keep for magnetic field
!           diagnostic
! nbrec = current record number for magnetic field writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: ntb = 0, modesxb = 41, nbrec = -1
!
! Radiative Vector Potential Diagnostic Parameters:
! ntar = number of time steps between radiative vector potential
!        diagnostic
! ndar = (0,1,2,3) = display (nothing,radiative vector potential,
!                             spectrum,both)
! modesxar = number of modes in x to keep for radiative vector
!            potential diagnostic
! narrec = current record number for radiative vector potential writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntar = 0, ndar = 1, modesxar = 41, narrec = -1
!
! Electromagnetic Power Spectrum Diagnostic Parameters:
! wrmin/wrmax = minimum/maximum high frequency used in power spectrum
! dwr = high frequency increment used in power spectrum
      real :: wrmin = 0.0, wrmax = 8.0, dwr = 0.01
!
! define namelist
      namelist /input1b/ omx, omy, omz, ntje, modesxje, njerec, nta,    &
     &nda, modesxa, narec, ntet, ndet, modesxet, netrec, ntb, modesxb,  &
     &nbrec, ntar, ndar, modesxar, narrec, wrmin, wrmax, dwr
!
! Darwin Namelist
! ndc = number of corrections in darwin iteration
      integer :: ndc = 2
!
! define namelist
      namelist /input1d/ ndc
!
! Ion Namelist
! Background Ion Parameters:
! npxi = number of background ions distributed in x direction
      integer :: npxi =  384
! qmi = charge on ion, in units of e
! rmass = ion/electron mass ratio
      real :: qmi = 1.0, rmass = 100.0
! rtempxi/rtempyi/rtempzi = electron/ion temperature ratio of background
! ions in x/y/z direction
      real :: rtempxi = 1.0, rtempyi = 1.0, rtempzi = 1.0
! vxi0/vyi0/vzi0 = drift velocity of ions in x/y/z direction
      real :: vxi0 = 0.0, vyi0 = 0.0, vzi0 = 0.0
!
! Beam Ion Parameters:
! npxbi = number of beam ions in x direction
      integer :: npxbi =   0
! vdxi/vdyi/vdzi = drift velocity of beam ions in x/y/z direction
      real :: vdxi = 0.0, vdyi = 0.0, vdzi = 0.0
! rtempdxi/rtempdyi/rtempdzi = electron/ion temperature ratio of beam
! ions in x/y/z direction
      real :: rtempdxi = 1.0, rtempdyi = 1.0, rtempdzi = 1.0
!
! Initial Ion Density Parameters:
! ndprofi = ion profile (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
      integer :: ndprofi = 0
! ampdxi = amplitude of ion density compared to uniform in x
! scaledxi = scale length for spatial ion coordinate in x
! shiftdxi = shift of spatial ion coordinate in x
      real :: ampdxi = 0.0, scaledxi = 0.0, shiftdxi = 0.0
!
! Ion Density Diagnostic Parameters:
! ntdi = number of time steps between ion density diagnostic
! nddi = (0,1,2,3) = display (nothing,ion density,spectrum,both)
! modesxdi = number of modes in x to keep for ion density diagnostic
! ndrec = current record number for ion density writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: ntdi = 0, nddi = 1, modesxdi = 41, ndirec = -1
!
! Ion Current Diagnostic Parameters:
! ntji = number of time steps between ion current diagnostic
! ndji = (0,1,2,3) = display (nothing,ion current,spectrum,both)
! modesxji = number of modes in x to keep for ion current diagnostic
! njirec = current record number for ion current writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntji = 0, ndji = 1, modesxji = 41, njirec = -1
!
! Ion Power Spectrum Diagnostic Parameters:
! wimin/wimax = minimum/maximum frequency used in ion power spectrum
! dwi = frequency increment used in ion power spectrum
      real :: wimin = 0.0, wimax = 0.1, dwi = 0.001
!
! define namelist
      namelist /ions1/ npxi, npxbi, qmi, rmass, rtempxi, rtempyi,       &
     &rtempzi, vxi0, vyi0, vzi0, vdxi, vdyi, vdzi, rtempdxi, rtempdyi,  &
     &rtempdzi, ndprofi, ampdxi, scaledxi, shiftdxi, ntdi, nddi,        &
     &modesxdi, ndirec, ntji, ndji, modesxji, njirec, wimin, wimax, dwi
!
! Output namelists
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
!
! Namelist output for electron density diagnostic
! fdename = file name for electron density diagnostic
      character(len=32) :: fdename = 'denek1.0'
! define namelist
      namelist /dene1d/ idrun, indx, ntde, modesxde, nderec, t0, tend,  &
     &dt, ceng, fdename
!
! Namelist output for potential diagnostic
! fpname = file name for potential diagnostic
      character(len=32) :: fpname = 'potk1.0'
! define namelist
      namelist /pot1d/ idrun, indx, ntp, modesxp, nprec, t0, tend, dt,  &
     &ceng, fpname
!
! Namelist output for longitudinal efield diagnostic
! felname = file name for longitudinal efield diagnostic
      character(len=32) :: felname = 'elk1.0'
! define namelist
      namelist /el1d/ idrun, indx, ntel, modesxel, nelrec, t0, tend, dt,&
     &ceng, felname
!
! Namelist output for electron current diagnostic
! fjename = file name for electron current diagnostic
      character(len=32) :: fjename = 'vcurek1.0'
! define namelist
      namelist /vcure1d/ idrun, indx, ntje, modesxje, ndim, omx, omy,   &
     &omz, ci, njerec, t0, tend, dt, ceng, fjename
!
! Namelist output for vector potential diagnostic
! faname = file name for vector potential diagnostic
      character(len=32) :: faname = 'vpotk1.0'
! define namelist
      namelist /vpot1d/ idrun, indx, nta, modesxa, ndim, omx, omy, omz, &
     &ci, narec, t0, tend, dt, ceng, faname
!
! Namelist output for transverse efield diagnostic
! fetname = file name for transverse efield diagnostic
      character(len=32) :: fetname = 'etk1.0'
! define namelist
      namelist /et1d/ idrun, indx, ntet, modesxet, netrec, t0, tend, dt,&
     &ceng, fetname
!
! Namelist output for magnetic field diagnostic
! fetname = file name for magnetic field diagnostic
      character(len=32) :: fbname = 'bk1.0'
! define namelist
      namelist /b1d/ idrun, indx, ntb, modesxb, nbrec, t0, tend, dt,    &
     &ceng, fbname
!
! Namelist output for radiative vector potential diagnostic
! farname = file name for vector potential diagnostic
      character(len=32) :: farname = 'vpotrk1.0'
! define namelist
      namelist /vpotr1d/ idrun, indx, ntar, modesxar, ndim, omx, omy,   &
     &omz, ci, narrec, t0, tend, dt, ceng, farname
!
! Namelist output for fluid moments diagnostic
! nprd = dimension of fluid moment arrays fmse and fmsi
      integer :: nprd = 0
! ffename/ffiname = file name for electron/ion fluid moments diagnostic
      character(len=32) :: ffename = 'fmer1.0', ffiname = 'fmir1.0'
! define namelist
      namelist /fm1d/ idrun, indx, ntfm, npro, ndim, nprd, nferec,      &
     &nfirec, t0, tend, dt, ffename, ffiname
!
! Namelist output for velocity-space diagnostic
! nfvd = dimension of velocity distribution arrays fv and fvi
! nfed = dimension of energy distribution arrays fe and fei
      integer :: nfvd = 0, nfed = 0
! fvename/fviname = file name for electron/ion velocity-space diagnostic
      character(len=32) :: fvename = 'fve1.0', fviname = 'fvi1.0'
! define namelist
      namelist /fv1d/ idrun, indx, ntv, nmv, nvft, ndim, nfvd, nfed,    &
     &omx, omy, omz, nverec, nvirec, t0, tend, dt, fvename, fviname
!
! Namelist output for trajectory diagnostic
! ndimp = size of phase space trajectories
      integer :: ndimp = 0
! ftname = file name for trajectory diagnostic
      character(len=32) :: ftname = 'tr1.0'
! define namelist
      namelist /tr1d/ idrun, indx, ntt, ndt, nst, nmv, ndim, ndimp,     &
    &nprobt, ntrec, t0, tend, dt, ftname
!
! Namelist output for phase space diagnostic
! nsxb = number of segments in x for global velocity distribution
      integer :: nsxb = 0
! fsename/fsiname = file name for electron/ion phase space diagnostic
      character(len=32) :: fsename = 'pse1.0', fsiname = 'psi1.0'
! define namelist
      namelist /ps1d/ idrun, indx, nts, nmv, ndim, nsxb, nserec, nsirec,&
     & t0, tend, dt, fsename, fsiname
!
! Namelist output for ion density diagnostic
! fdname = file name for ion density diagnostic
      character(len=32) :: fdiname = 'denik1.0'
! define namelist
      namelist /deni1d/ idrun, indx, ntdi, modesxdi, ndirec, t0, tend,  &
     &dt, ceng, fdiname 
!
! Namelist output for ion current diagnostic
! fjiname = file name for ion current diagnostic
      character(len=32) :: fjiname = 'vcurik1.0'
! define namelist
      namelist /vcuri1d/ idrun, indx, ntji, modesxji, ndim, omx, omy,   &
     &omz, ci, njirec, t0, tend, dt, ceng, fjiname
!
      contains
!
      subroutine readnml1(iuin)
! read namelist from unit iuin
      implicit none
      integer, intent(in) :: iuin
! local data
      integer :: ios
      open(unit=iuin,file='input1',form='formatted',status='old')
! read global input parameters
      read (iuin,input1)
! read electromagnetic input parameters
      if ((emf==1).or.(emf==2)) then
         read (iuin,input1b,iostat=ios)
         if (ios /= 0) then
            write (*,*) 'error in reading input1b namelist'
            rewind iuin
         endif
      endif
! read darwin input parameters
      if (emf==2) then
         read (iuin,input1d,iostat=ios)
         if (ios /= 0) then
            write (*,*) 'error in reading input1d namelist'
            rewind iuin
         endif
      endif
! read ion input parameters
      if (movion==1) then
         read (iuin,ions1,iostat=ios)
         if (ios /= 0) write (*,*) 'error in reading ions1 namelist'
      endif
      rewind iuin
      end subroutine
!
      subroutine writnml1(iudm)
! write final diagnostic metafile to unit iudm
      implicit none
      integer, intent(in) :: iudm
! local data
      character(len=10) :: cdrun
      character(len=32) :: fname
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag1.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
! write out global input parameters
      write (iudm,input1)
! write out electromagnetic input parameters
      if ((emf==1).or.(emf==2)) then
         write (iudm,input1b)
      endif
! write out darwin input parameters
      if (emf==2) then
         write (iudm,input1d)
      endif
! electron density diagnostic
      if (ntde > 0) then
         write (iudm,dene1d)
      endif
! potential diagnostic
      if (ntp > 0) then
         write (iudm,pot1d)
      endif
! longitudinal efield diagnostic
      if (ntel > 0) then
         write (iudm,el1d)
      endif
! electron current diagnostic
      if (ntje > 0) then
         ceng = 0.0
         write (iudm,vcure1d)
      endif
! vector potential diagnostic
      if (nta > 0) then
         write (iudm,vpot1d)
      endif
! transverse efield diagnostic
      if (ntet > 0) then
         write (iudm,et1d)
      endif
! magnetic field diagnostic
      if (ntb > 0) then
         write (iudm,b1d)
      endif
! radiative vector potential diagnostic
      if (ntar > 0) then
         write (iudm,vpotr1d)
      endif
! fluid moments diagnostic
      if (ntfm > 0) then
         write (iudm,fm1d)
      endif
! velocity-space diagnostic
      if (ntv > 0) then
         write (iudm,fv1d)
      endif
! trajectory diagnostic
      if (ntt > 0) then
         write (iudm,tr1d)
      endif
! phase space diagnostic
      if (nts > 0) then
         write (iudm,ps1d)
      endif
! ion parameters
      if (movion==1) then
! write out ion input parameters
         write (iudm,ions1)
! ion density diagnostic
         if (ntdi > 0) then
            ceng = 0.0
            write (iudm,deni1d)
         endif
! ion current diagnostic
         if (ntji > 0) then
            ceng = 0.0
            write (iudm,vcuri1d)
         endif
      endif
      close(unit=iudm)
      end subroutine
!
      subroutine closeff(iunit)
! close Fortran file with unit number iunit
      implicit none
      integer, intent(in) :: iunit
      close(unit=iunit)
      end subroutine
!
      end module
