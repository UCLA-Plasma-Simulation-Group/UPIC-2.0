!-----------------------------------------------------------------------
! Input for 2-1/2D MPI/OpenMP PIC codes
!
      module in2
!
! input2mod.f defines namelists containing input and output variables:
! readnml2 reads namelist from unit iuin
! writnml2 writes final diagnostic metafile to unit iudm
! written by viktor k. decyk, ucla
! copyright 2011, regents of the university of california
! update: March 12, 2018
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
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer :: indx =   9, indy =   9
! psolve = type of poisson solver = (1,2,3)
!     integer :: psolve = PERIODIC_2D
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: relativity = 0
! ci = reciprocal of velocity of light
      real :: ci = 0.1
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! ndim = number of velocity dimensions = 2 or 3
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
! npx/npy = number of background electrons distributed in x/y direction
      integer :: npx = 3072, npy = 3072
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
!
! Beam Electron Parameters:
! npxb/npyb = number of beam electrons distributed in x/y direction
      integer :: npxb = 0, npyb =  0
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
!
! Time Parameters:
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
      real :: tend = 10.000, dt = 0.1
!
! Numerical Parameters:
! inorder = interpolation order
! popt = particle optimization scheme
! dopt = charge deposit optimization scheme
! djopt = current deposit optimization scheme
!     integer :: inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD
!     integer :: djopt = STANDARD
! ax/ay = smoothed particle size in x/y direction
!     real :: ax = .816497, ay = .816497
!     real :: ax = .866025, ay = .866025
      real :: ax = .912871, ay = .912871
! mx/my = number of grids in x/y in sorting tiles
! should be less than or equal to 32
      integer :: mx = 16, my = 16
! nextrand = (0,N) = generate (default,Nth block) of random numbers
      integer :: nextrand = 0
!
! Initial Electron Density Parameters:
! density profile is of form n(x)*n(y), where
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
! for ndprof = 0, n(x) = n0
! for ndprof = 1, n(x) = n0*(1 + ampdx*(x/nx - 0.5))
! for ndprof = 2, n(x) = n0*(1 + ampdx*sin(x/scaledx - shift))
! for ndprof = (3,4,5), n(x) = n0*(1 + ampdx*f((x - shift)/scaledx))
! where f = (exp(-x**2/2),sech(x)**2,exp(x))
! n0 is determined by the requirement that the integral over density
! equals the total number of particles distributed by the function
! n(y) is the same function as n(x), but with different parameters.
      integer :: ndprof = 0
! ampdx/ampdy = amplitude of density compared to uniform in x/y
! scaledx/scaledy = scale length for spatial coordinate in x/y
! shiftdx/shiftdy = shift of spatial coordinate in x/y
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
      real :: ampdy = 0.0, scaledy = 0.0, shiftdy = 0.0
!
! Zero Force Parameter:
! mzf = (0,1) = (no,yes) set forces to zero
      integer :: mzf = 0
!
! External Traveling Wave Driver:
! Ext(x) = e1*sin(k0*x + freq*time) + e2*cos(k0*x - freq*time)
! amodex = wave number which determines k0 = 2*pi*amodex/NX
! freq = frequency of external wave driver
! trmp = ramp-up time for external wave driver
! toff = shut-off time for external wave driver
      real :: amodex = 0.0, freq = 0.0, trmp = 0.0, toff = 0.0
! el0/er0 = external pump amplitude for left-going/right-going wave
!     e1 = el0*(time/trmp), e2 = er0*(time/trmp), if time < trmp
!     e1 = el0,             e2 = er0,             if trmp < time < toff
!     e1 = 0,               e2 = 0,               if time > toff
      real :: el0 = 0.0, er0 = 0.0
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
! modesxde/modesyde = number of modes in x/y to keep for electron
!                     density diagnostic
! nderec = current record number for electron density writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntde = 0, modesxde = 41, modesyde = 41, nderec = -1
!
! Potential Diagnostic Parameters:
! ntp = number of time steps between potential diagnostic
! ndp = (0,1,2,3) = display (nothing,potential,spectrum,both)
! modesxp/modesyp = number of modes in x/y to keep for potential
!                   diagnostic
      integer :: ntp = 0, ndp = 1, modesxp = 41, modesyp = 41
! nprec = current record number for potential writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: nprec = -1
!
! Longitudinal Efield Diagnostic Parameters:
! ntel = number of time steps between longitudinal efield diagnostic
! modesxel/modesyel = number of modes in x/y to keep for longitudinal
!                     efield diagnostic
! nelrec = current record number for longitudinal efield writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntel = 0, modesxel = 41, modesyel = 41, nelrec = -1
!
! Power Spectrum Diagnostic Parameters:
! wmin/wmax = minimum/maximum frequency used in power spectrum
! dw = frequency increment used in power spectrum
      real :: wmin = 0.0, wmax = 2.0, dw = 0.01
!
! Fluid Moments Diagnostic Parameter:
! ntfm = number of time steps between fluid moments diagnostic
! ndfm = (0,1,2,3) = process (nothing,electrons,ions,both)
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
! ndv = (0,1,2,3) = process (nothing,electrons,ions,both)
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
! nst = type of test particle distribution:
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
      integer :: nts = 0, nds = 3
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
! Display Parameters:
! nplot = maximum number of plots per page
      integer :: nplot = 4
! ndstyle = (1,2,3) = display (color map,contour plot,both)
      integer :: ndstyle = 1
!
! Multi-tasking Parameter:
! nvp = number of distributed memory nodes
! nvpp = number of shared memory nodes (0=default)
      integer :: nvp = 1, nvpp = 0
! imbalance = load imbalance fraction repartition trigger
! (< 0.0  to suppress repartition)
      real :: imbalance = -1.0
!
! Error Processing Parameter
! monitor = (0,1,2) = (disable,normal,extended) error processing
      integer :: monitor = 0
!
! define namelist
      namelist /input2/ idrun, idcode, indx, indy, mx, my, npx, npy,    &
     &npxb, npyb, qme, vtx, vty, vtz, vx0, vy0, vz0, vdx, vdy, vdz,     &
     &vtdx, vtdy, vtdz, relativity, ci, xtras, ndim, nvdist, treverse,  &
     &tend, dt, ax, ay, nextrand, mzf, ndprof, ampdx, scaledx, shiftdx, &
     &ampdy, scaledy, shiftdy, amodex, freq, trmp, toff, el0, er0, ntw, &
     &ndw, ntde, modesxde, modesyde, nderec, ntp, ndp, modesxp, modesyp,&
     &nprec, ntel, modesxel, modesyel, nelrec, wmin, wmax, dw, ntfm,    &
     &ndfm, npro, nferec, nfirec, ntv, ndv, nmv, nvft, nverec, nvirec,  &
     &ntt, ndt, nst, nprobt, vtsx, dvtx, ntrec, nts, nds, mvx, mvy,     &
     &nserec, nsirec, movion, emf, nustrt, ntr, idrun0, nplot, ndstyle, &
     &nvpp, imbalance, monitor
!
! equivalence data to simplify MPI broadcast
      integer, parameter :: lnin2 = 100
      double precision, dimension(lnin2) :: ddin2
      private :: lnin2, ddin2
      equivalence (ddin2(1),idrun), (ddin2(2),idcode), (ddin2(3),indx)
      equivalence (ddin2(4),indy), (ddin2(5),mx), (ddin2(6),my)
      equivalence (ddin2(7),npx), (ddin2(8),npy), (ddin2(9),npxb)
      equivalence (ddin2(10),npyb), (ddin2(11),qme), (ddin2(12),vtx)
      equivalence (ddin2(13),vty), (ddin2(14),vtz), (ddin2(15),vx0)
      equivalence (ddin2(16),vy0), (ddin2(17),vz0), (ddin2(18),vdx)
      equivalence (ddin2(19),vdy), (ddin2(20),vdz), (ddin2(21),vtdx)
      equivalence (ddin2(22),vtdy), (ddin2(23),vtdz)
      equivalence (ddin2(24),relativity), (ddin2(25),ci)
      equivalence (ddin2(26),xtras), (ddin2(27),ndim)
      equivalence (ddin2(28),nvdist), (ddin2(29),treverse)
      equivalence (ddin2(30),tend), (ddin2(31),dt), (ddin2(32),ax)
      equivalence (ddin2(33),ay), (ddin2(34),nextrand), (ddin2(35),mzf)
      equivalence (ddin2(36),ndprof), (ddin2(37),ampdx)
      equivalence (ddin2(38),scaledx), (ddin2(39),shiftdx)
      equivalence (ddin2(40),ampdy), (ddin2(41),scaledy)
      equivalence (ddin2(42),shiftdy), (ddin2(43),amodex)
      equivalence (ddin2(44),freq), (ddin2(45),trmp), (ddin2(46),toff)
      equivalence (ddin2(47),el0), (ddin2(48),er0), (ddin2(49),ntw)
      equivalence (ddin2(50),ndw), (ddin2(51),ntde)
      equivalence (ddin2(52),modesxde), (ddin2(53),modesyde)
      equivalence (ddin2(54),nderec), (ddin2(55),ntp), (ddin2(56),ndp)
      equivalence (ddin2(57),modesxp), (ddin2(58),modesyp)
      equivalence (ddin2(59),nprec), (ddin2(60),ntel)
      equivalence (ddin2(61),modesxel), (ddin2(62),modesyel)
      equivalence (ddin2(63),nelrec), (ddin2(64),wmin), (ddin2(65),wmax)
      equivalence (ddin2(66),dw), (ddin2(67),ntfm), (ddin2(68),ndfm)
      equivalence (ddin2(69),npro), (ddin2(70),nferec)
      equivalence (ddin2(71),nfirec), (ddin2(72),ntv), (ddin2(73),ndv)
      equivalence (ddin2(74),nmv), (ddin2(75),nvft), (ddin2(76),nverec)
      equivalence (ddin2(77),nvirec), (ddin2(78),ntt), (ddin2(79),ndt)
      equivalence (ddin2(80),nst), (ddin2(81),nprobt), (ddin2(82),vtsx)
      equivalence (ddin2(83),dvtx), (ddin2(84),ntrec), (ddin2(85),nts)
      equivalence (ddin2(86),nds), (ddin2(87),mvx), (ddin2(88),mvy)
      equivalence (ddin2(89),nserec), (ddin2(90),nsirec)
      equivalence (ddin2(91),movion), (ddin2(92),emf)
      equivalence (ddin2(93),nustrt), (ddin2(94),ntr)
      equivalence (ddin2(95),idrun0), (ddin2(96),nplot)
      equivalence (ddin2(97),ndstyle), (ddin2(98),nvpp)
      equivalence (ddin2(99),imbalance), (ddin2(100),monitor)
!
! Electromagnetic Namelist
! External Magnetic Field Parameters:
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
!
! Electron Current Diagnostic Parameters:
! ntje = number of time steps between electron current diagnostic
! modesxje/modesyje = number of modes in x/y to keep for electron
!                     current diagnostic
! njerec = current record number for electron current writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: ntje = 0, modesxje = 41, modesyje = 41, njerec = -1
!
! Vector Potential Diagnostic Parameters:
! nta = number of time steps between vector potential diagnostic
! nda = (0,1,2,3) = display (nothing,vector potential,spectrum,both)
! modesxa/modesya = number of modes in x/y to keep for vector potential
!                   diagnostic
      integer :: nta = 0, nda = 1, modesxa = 41, modesya = 41
! narec = current record number for vector potential writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: narec = -1
!
! Transverse Efield Diagnostic Parameters:
! ntet = number of time steps between transverse efield diagnostic
! ndet = (0,1,2,3) = display (nothing,transverse efield,spectrum,both)
! modesxet/modesyet = number of modes in x/y to keep for transverse
!                     efield diagnostic
      integer :: ntet = 0, ndet = 1, modesxet = 41, modesyet = 41
! netrec = current record number for ltransverse efield writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: netrec = -1
!
! Magnetic Field Diagnostic Parameters:
! ntb = number of time steps between magnetic field diagnostic
! modesxb/modesyb = number of modes in x/y to keep for magnetic field
!                   diagnostic
! nbrec = current record number for magnetic field writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: ntb = 0, modesxb = 41, modesyb = 41, nbrec = -1
!
! Radiative Vector Potential Diagnostic Parameters:
! ntar = number of time steps between radiative vector potential
!        diagnostic
! ndar = (0,1,2,3) = display (nothing,radiative vector potential,
!                             spectrum,both)
! modesxar/modesyar = number of modes in x/y to keep for radiative
!                     vector potential diagnostic
      integer :: ntar = 0, ndar = 1, modesxar = 41, modesyar = 41
! narrec = current record number for radiative vector potential writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: narrec = -1
!
! Electromagnetic Power Spectrum Diagnostic Parameters:
! wrmin/wrmax = minimum/maximum high frequency used in power spectrum
! dwr = high frequency increment used in power spectrum
      real :: wrmin = 0.0, wrmax = 8.0, dwr = 0.01
!
! define namelist
      namelist /input2b/ omx, omy, omz, ntje, modesxje, modesyje,       &
     &njerec, nta, nda, modesxa, modesya, narec, ntet, ndet, modesxet,  &
     &modesyet, netrec, ntb, modesxb, modesyb, nbrec, ntar, ndar,       &
     &modesxar, modesyar, narrec, wrmin, wrmax, dwr
!
! equivalence data to simplify MPI broadcast
      integer, parameter :: lnin2b = 29
      double precision, dimension(lnin2b) :: ddin2b
      private :: lnin2b, ddin2b
      equivalence (ddin2b(1),omx), (ddin2b(2),omy), (ddin2b(3),omz)
      equivalence (ddin2b(4),ntje), (ddin2b(5),modesxje)
      equivalence (ddin2b(6),modesyje), (ddin2b(7),njerec)
      equivalence (ddin2b(8),nta), (ddin2b(9),nda), (ddin2b(10),modesxa)
      equivalence (ddin2b(11),modesya), (ddin2b(12),narec)
      equivalence (ddin2b(13),ntet), (ddin2b(14),ndet)
      equivalence (ddin2b(15),modesxet), (ddin2b(16),modesyet)
      equivalence (ddin2b(17),netrec), (ddin2b(18),ntb)
      equivalence (ddin2b(19),modesxb), (ddin2b(20),modesyb)
      equivalence (ddin2b(21),nbrec), (ddin2b(22),ntar)
      equivalence (ddin2b(23),ndar), (ddin2b(24),modesxar)
      equivalence (ddin2b(25),modesyar), (ddin2b(26),narrec)
      equivalence (ddin2b(27),wrmin), (ddin2b(28),wrmax)
      equivalence (ddin2b(29),dwr)
!
! Darwin Namelist
! ndc = number of corrections in darwin iteration
      integer :: ndc = 2
!
! define namelist
      namelist /input2d/ ndc
!
! equivalence data to simplify MPI broadcast
      integer, parameter :: lnin2d = 1
      double precision, dimension(lnin2d) :: ddin2d
      private :: lnin2d, ddin2d
      equivalence (ddin2d(1),ndc)
!
! Ion Namelist
! Background Ion Parameters:
! npxi/npyi = number of background ions distributed in x/y direction
      integer :: npxi = 3072, npyi = 3072
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
! npxbi/npybi = number of beam ions distributed in x/y direction
      integer :: npxbi = 0, npybi = 0
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
! ampdxi/ampdyi = amplitude of ion density compared to uniform in x/y
! scaledxi/scaledyi = scale length for spatial ion coordinate in x/y
! shiftdxi/shiftdyi = shift of spatial ion coordinate in x/y
      real :: ampdxi = 0.0, scaledxi = 0.0, shiftdxi = 0.0
      real :: ampdyi = 0.0, scaledyi = 0.0, shiftdyi = 0.0
!
! Ion Density Diagnostic Parameters:
! ntdi = number of time steps between ion density diagnostic
! nddi = (0,1,2,3) = display (nothing,ion density,spectrum,both)
! modesxdi/modesydi = number of modes in x/y to keep for ion density
!                     diagnostic
      integer :: ntdi = 0, nddi = 1, modesxdi = 41, modesydi = 41
! ndrec = current record number for ion density writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: ndirec = -1
!
! Ion Current Diagnostic Parameters:
! ntji = number of time steps between ion current diagnostic
! ndji = (0,1,2,3) = display (nothing,ion current,spectrum,both)
! modesxji/modesyji = number of modes in x/y to keep for ion current
!                     diagnostic
      integer :: ntji = 0, ndji = 1, modesxji = 41, modesyji = 41
! njirec = current record number for ion current writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: njirec = -1
!
! Ion Power Spectrum Diagnostic Parameters:
! wimin/wimax = minimum/maximum frequency used in ion power spectrum
! dwi = frequency increment used in ion power spectrum
      real :: wimin = 0.0, wimax = 0.1, dwi = 0.001
!
! define namelist
      namelist /ions2/ npxi, npyi, npxbi, npybi, qmi, rmass, rtempxi,   &
     &rtempyi, rtempzi, vxi0, vyi0, vzi0, vdxi, vdyi, vdzi, rtempdxi,   &
     &rtempdyi, rtempdzi, ndprofi, ampdxi, scaledxi, shiftdxi, ampdyi,  &
     &scaledyi, shiftdyi, ntdi, nddi, modesxdi, modesydi, ndirec, ntji, &
     &ndji, modesxji, modesyji, njirec, wimin, wimax, dwi
!
! equivalence data to simplify MPI broadcast
      integer, parameter :: lnion2 = 38
      double precision, dimension(lnion2) :: ddion2
      private :: lnion2, ddion2
      equivalence (ddion2(1),npxi), (ddion2(2),npyi), (ddion2(3),npxbi)
      equivalence (ddion2(4),npybi), (ddion2(5),qmi), (ddion2(6),rmass)
      equivalence (ddion2(7),rtempxi), (ddion2(8),rtempyi)
      equivalence (ddion2(9),rtempzi), (ddion2(10),vxi0)
      equivalence (ddion2(11),vyi0), (ddion2(12),vzi0)
      equivalence (ddion2(13),vdxi), (ddion2(14),vdyi)
      equivalence (ddion2(15),vdzi), (ddion2(16),rtempdxi)
      equivalence (ddion2(17),rtempdyi), (ddion2(18),rtempdzi)
      equivalence (ddion2(19),ndprofi), (ddion2(20),ampdxi)
      equivalence (ddion2(21),scaledxi), (ddion2(22),shiftdxi)
      equivalence (ddion2(23),ampdyi), (ddion2(24),scaledyi)
      equivalence (ddion2(25),shiftdyi), (ddion2(26),ntdi)
      equivalence (ddion2(27),nddi), (ddion2(28),modesxdi)
      equivalence (ddion2(29),modesydi), (ddion2(30),ndirec)
      equivalence (ddion2(31),ntji), (ddion2(32),ndji)
      equivalence (ddion2(33),modesxji), (ddion2(34),modesyji)
      equivalence (ddion2(35),njirec), (ddion2(36),wimin)
      equivalence (ddion2(37),wimax), (ddion2(38),dwi)
!
! Output namelists
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
!
! Namelist output for electron density diagnostic
! fdename = file name for electron density diagnostic
      character(len=32) :: fdename = 'denek2.0'
! define namelist
      namelist /dene2d/ idrun, indx, indy, ntde, modesxde, modesyde,    &
     &nderec, nvp, t0, tend, dt, ceng, fdename
!
! Namelist output for potential diagnostic
! fpname = file name for potential diagnostic
      character(len=32) :: fpname = 'potk2.0'
! define namelist
      namelist /pot2d/ idrun, indx, indy, ntp, modesxp, modesyp, nprec, &
     &nvp, t0, tend, dt, ceng, fpname
!
! Namelist output for longitudinal efield diagnostic
! felname = file name for longitudinal efield diagnostic
      character(len=32) :: felname = 'elk2.0'
! define namelist
      namelist /el2d/ idrun, indx, indy, ntel, modesxel, modesyel, ndim,&
     &nelrec, nvp, t0, tend, dt, ceng, felname
!
! Namelist output for electron current diagnostic
! fjename = file name for electron current diagnostic
      character(len=32) :: fjename = 'vcurek2.0'
! define namelist
      namelist /vcure2d/ idrun, indx, indy, ntje, modesxje, modesyje,   &
     &ndim, omx, omy, omz, ci, njerec, nvp, t0, tend, dt, ceng, fjename
!
! Namelist output for vector potential diagnostic
! faname = file name for vector potential diagnostic
      character(len=32) :: faname = 'vpotk2.0'
! define namelist
      namelist /vpot2d/ idrun, indx, indy, nta, modesxa, modesya, ndim, &
     &omx, omy, omz, ci, narec, nvp, t0, tend, dt, ceng, faname
!
! Namelist output for transverse efield diagnostic
! fetname = file name for transverse efield diagnostic
      character(len=32) :: fetname = 'etk2.0'
! define namelist
      namelist /et2d/ idrun, indx, indy, ntet, modesxet, modesyet, ndim,&
     &netrec, nvp, t0, tend, dt, ceng, fetname
!
! Namelist output for magnetic field diagnostic
! fetname = file name for magnetic field diagnostic
      character(len=32) :: fbname = 'bk2.0'
! define namelist
      namelist /b2d/ idrun, indx, indy, ntb, modesxb, modesyb, ndim,    &
     &nbrec, nvp, t0, tend, dt, ceng, fbname
!
! Namelist output for radiative vector potential diagnostic
! farname = file name for vector potential diagnostic
      character(len=32) :: farname = 'vpotrk2.0'
! define namelist
      namelist /vpotr2d/ idrun, indx, indy, ntar, modesxar, modesyar,   &
     &ndim, omx, omy, omz, ci, narrec, nvp, t0, tend, dt, ceng, farname
!
! Namelist output for fluid moments diagnostic
! nprd = dimension of fluid moment arrays fmse and fmsi
      integer :: nprd = 0
! ffename/ffiname = file name for electron/ion fluid moments diagnostic
      character(len=32) :: ffename = 'fmer2.0', ffiname = 'fmir2.0'
! define namelist
      namelist /fm2d/ idrun, indx, indy, ntfm, npro, ndim, nprd, nferec,&
     &nfirec, nvp, t0, tend, dt, ceng, ffename, ffiname
!
! Namelist output for velocity-space diagnostic
! nfvd = dimension of velocity distribution arrays fv and fvi
! nfed = dimension of energy distribution arrays fe and fei
      integer :: nfvd = 0, nfed = 0
! fvename/fviname = file name for electron/ion velocity-space diagnostic
      character(len=32) :: fvename = 'fve2.0', fviname = 'fvi2.0'
! define namelist
      namelist /fv2d/ idrun, indx, indy, ntv, nmv, nvft, ndim, nfvd,    &
     &nfed, omx, omy, omz, nverec, nvirec, nvp, t0, tend, dt, fvename,  &
     &fviname
!
! Namelist output for trajectory diagnostic
! ndimp = size of phase space trajectories
      integer :: ndimp = 0
! ftname = file name for trajectory diagnostic
      character(len=32) :: ftname = 'tr2.0'
! define namelist
      namelist /tr2d/ idrun, indx, indy, ntt, ndt, nst, nmv, ndim,      &
     &ndimp, nprobt, ntrec, nvp, t0, tend, dt, ftname
!
! Namelist output for phase space diagnostic
! nsxb/nsyb = number of segments in x/y for global velocity distribution
      integer :: nsxb = 0, nsyb= 0
! fsename/fsiname = file name for electron/ion phase space diagnostic
      character(len=32) :: fsename = 'pse2.0', fsiname = 'psi2.0'
! define namelist
      namelist /ps2d/ idrun, indx, indy, nts, nmv, ndim, nsxb, nsyb,    &
     &nserec, nsirec, nvp, t0, tend, dt, fsename, fsiname
!
! Namelist output for ion density diagnostic
! fdname = file name for ion density diagnostic
      character(len=32) :: fdiname = 'denik2.0'
! define namelist
      namelist /deni2d/ idrun, indx, indy, ntdi, modesxdi, modesydi,    &
     &ndirec, nvp, t0, tend, dt, ceng, fdiname
!
! Namelist output for ion current diagnostic
! fjiname = file name for ion current diagnostic
      character(len=32) :: fjiname = 'vcurik2.0'
! define namelist
      namelist /vcuri2d/ idrun, indx, indy, ntji, modesxji, modesyji,   &
     &ndim, omx, omy, omz, ci, njirec, nvp, t0, tend, dt, ceng, fjiname
!
      contains
!
      subroutine readnml2(iuin)
! read namelist from unit iuin
      implicit none
      integer, intent(in) :: iuin
! local data
      integer :: ios
      open(unit=iuin,file='input2',form='formatted',status='old')
! read global input parameters
      read (iuin,input2)
! read electromagnetic input parameters
      if ((emf==1).or.(emf==2)) then
         read (iuin,input2b,iostat=ios)
         if (ios /= 0) then
            write (*,*) 'error in reading input2b namelist'
            rewind iuin
         endif
      endif
! read darwin input parameters
      if (emf==2) then
         read (iuin,input2d,iostat=ios)
         if (ios /= 0) then
            write (*,*) 'error in reading input2d namelist'
            rewind iuin
         endif
      endif
! read ion input parameters
      if (movion==1) then
         read (iuin,ions2,iostat=ios)
         if (ios /= 0) write (*,*) 'error in reading ions2 namelist'
      endif
      rewind iuin
      end subroutine
!
      subroutine writnml2(iudm)
! write final diagnostic metafile to unit iudm
      implicit none
      integer, intent(in) :: iudm
! local data
      character(len=10) :: cdrun
      character(len=32) :: fname
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
! write out global input parameters
      write (iudm,input2)
! write out electromagnetic input parameters
      if ((emf==1).or.(emf==2)) then
         write (iudm,input2b)
      endif
! write out darwin input parameters
      if (emf==2) then
         write (iudm,input2d)
      endif
! electron density diagnostic
      if (ntde > 0) then
         write (iudm,dene2d)
      endif
! potential diagnostic
      if (ntp > 0) then
         write (iudm,pot2d)
      endif
! longitudinal efield diagnostic
      if (ntel > 0) then
         write (iudm,el2d)
      endif
! electron current diagnostic
      if (ntje > 0) then
         ceng = 0.0
         write (iudm,vcure2d)
      endif
! vector potential diagnostic
      if (nta > 0) then
         write (iudm,vpot2d)
      endif
! transverse efield diagnostic
      if (ntet > 0) then
         write (iudm,et2d)
      endif
! magnetic field diagnostic
      if (ntb > 0) then
         write (iudm,b2d)
      endif
! radiative vector potential diagnostic
      if (ntar > 0) then
         write (iudm,vpotr2d)
      endif
! fluid moments diagnostic
      if (ntfm > 0) then
         write (iudm,fm2d)
      endif
! velocity-space diagnostic
      if (ntv > 0) then
         write (iudm,fv2d)
      endif
! trajectory diagnostic
      if (ntt > 0) then
         write (iudm,tr2d)
      endif
! phase space diagnostic
      if (nts > 0) then
         write (iudm,ps2d)
      endif
! ion parameters
      if (movion==1) then
! write out ion input parameters
         write (iudm,ions2)
! ion density diagnostic
         if (ntdi > 0) then
            ceng = 0.0
            write (iudm,deni2d)
         endif
! ion current diagnostic
         if (ntji > 0) then
            ceng = 0.0
            write (iudm,vcuri2d)
         endif
      endif
      end subroutine
!
      subroutine sendnmls2()
! this subroutine broadcasts namelists to other nodes
! namelists are equivalenced to a double precision array
      implicit none
! broadcase namelist input2
      call PPBDCAST(ddin2,lnin2)
! broadcast namelist input2b
      call PPBDCAST(ddin2b,lnin2b)
! broadcast namelist input2d
      call PPBDCAST(ddin2d,lnin2d)
! broadcast namelist ions2
      call PPBDCAST(ddion2,lnion2)
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
