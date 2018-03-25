!-----------------------------------------------------------------------
! Input for 3D MPI/OpenMP PIC codes
!
      module in3
!
! input3mod.f defines namelists containing input and output variables:
! readnml3 reads namelist from unit iuin
! writnml3 writes final diagnostic metafile to unit iudm
! written by viktor k. decyk, ucla
! copyright 2011, regents of the university of california
! update: March 14, 2018
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
! indx/indy/indz = exponent which determines grid points in x/y/z
! direction: nx = 2**indx, ny = 2**indy, nz = 2**indz.
      integer :: indx =   7, indy =   7, indz =   7
! psolve = type of poisson solver = (1,2,3)
!     integer :: psolve = PERIODIC_2D
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: relativity = 0
! ci = reciprocal of velocity of light
      real :: ci = 0.1
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! ndim = number of velocity dimensions = 3
      integer :: ndim = 3
! nvdist = velocity distribution type
! nvdist = (1,2) = (maxwellian/juttner,waterbag) distribution
! for nvdist=2, maximum velocity in x/y/z is (vtx/vty/vtz)*sqrt(3)
      integer :: nvdist = 1
! treverse = (0,1) = (no,yes) reverse simulation at end back to start
      integer :: treverse = 0
!
! Background Electron Parameters:
! npx/npy/npz = number of background electrons distributed in x/y/z
! direction
      integer :: npx = 384, npy = 384, npz = 384
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
!
! Beam Electron Parameters:
! npxb/npyb/npzb = number of beam electrons distributed in x/y/z
! direction
      integer :: npxb = 0, npyb = 0, npzb = 0
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
! ax/ay/az = smoothed particle size in x/y/z direction
!     real :: ax = .816497, ay = .816497, az = .816497
!     real :: ax = .866025, ay = .866025, az = .866025
      real :: ax = .912871, ay = .912871, az = .912871
! mx/my/mz = number of grids in x/y/z in sorting tiles
! should be less than or equal to 16
      integer :: mx = 8, my = 8, mz = 8
! nextrand = (0,N) = generate (default,Nth block) of random numbers
      integer :: nextrand = 0
!
! Initial Electron Density Parameters:
! density profile is of form n(x)*n(y)*n(z), where
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
! for ndprof = 0, n(x) = n0
! for ndprof = 1, n(x) = n0*(1 + ampdx*(x/nx - 0.5))
! for ndprof = 2, n(x) = n0*(1 + ampdx*sin(x/scaledx - shift))
! for ndprof = (3,4,5), n(x) = n0*(1 + ampdx*f((x - shift)/scaledx))
! where f = (exp(-x**2/2),sech(x)**2,exp(x))
! n0 is determined by the requirement that the integral over density
! equals the total number of particles distributed by the function
! n(y) and n(z) are the same function as n(x), but with different
! parameters.
      integer :: ndprof = 0
! ampdx/ampdy/ampdz = amplitude of density compared to uniform in x/y/z
! scaledx/scaledy/scaledz = scale length for spatial coordinate in x/y/z
! shiftdx/shiftdy/shiftdz = shift of spatial coordinate in x/y/z
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
      real :: ampdy = 0.0, scaledy = 0.0, shiftdy = 0.0
      real :: ampdz = 0.0, scaledz = 0.0, shiftdz = 0.0
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
! modesxde/modesyde/modeszde = number of modes in x/y/z to keep for
!                              electron density diagnostic
      integer :: ntde = 0, modesxde = 41, modesyde = 41, modeszde = 41
! nderec = current record number for electron density writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: nderec = -1
!
! Potential Diagnostic Parameters:
! ntp = number of time steps between potential diagnostic
! ndp = (0,1,2,3) = display (nothing,potential,spectrum,both)
      integer :: ntp = 0, ndp = 1
! modesxp/modesyp/modeszp = number of modes in x/y/z to keep for
!                           potential diagnostic
      integer :: modesxp = 41, modesyp = 41, modeszp = 41
! nprec = current record number for potential writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: nprec = -1
!
! Longitudinal Efield Diagnostic Parameters:
! ntel = number of time steps between longitudinal efield diagnostic
! modesxel/modesyel/modeszel = number of modes in x/y/z to keep for
!                              longitudinal efield diagnostic
      integer :: ntel = 0, modesxel = 41, modesyel = 41, modeszel = 41
! nelrec = current record number for longitudinal efield writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: nelrec = -1
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
!           (0 for beginnning of file, -1 to disable writes)
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
      integer :: nts = 0, nds = 3
! mvx/mvy/mvz = number of grids in x/y/z for phase space aggregation
      integer :: mvx = 3, mvy = 3, mvz = 3
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
! Multi-tasking Parameter:
! nvpy/nvpz = number of distributed memory nodes in x/y
! nvpp = number of shared memory nodes (0=default)
      integer :: nvpy = 1, nvpz = 1, nvpp = 0
! imbalance = load imbalance fraction repartition trigger
! (< 0.0  to suppress repartition)
      real :: imbalance = -1.0
!
! Error Processing Parameter
! monitor = (0,1,2) = (disable,normal,extended) error processing
      integer :: monitor = 0
!
! define namelist
      namelist /input3/ idrun, idcode, indx, indy, indz, mx, my, mz,    &
     &npx, npy, npz, npxb, npyb, npzb, qme, vtx, vty, vtz, vx0, vy0,    &
     &vz0, vdx, vdy, vdz, vtdx, vtdy, vtdz, relativity, ci, xtras, ndim,&
     &nvdist, treverse, tend, dt, ax, ay, az, nextrand, mzf, ndprof,    &
     &ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy, ampdz, scaledz, &
     &shiftdz, amodex, freq, trmp, toff, el0, er0, ntw, ndw, ntde,      &
     &modesxde, modesyde, modeszde, nderec, ntp, ndp, modesxp, modesyp, &
     &modeszp, nprec, ntel, modesxel, modesyel, modeszel, nelrec, wmin, &
     &wmax, dw, ntfm, ndfm, npro, nferec, nfirec, ntv, ndv, nmv, nvft,  &
     &nverec, nvirec, ntt, ndt, nst, nprobt, vtsx, dvtx, ntrec, nts,    &
     &nds, mvx, mvy, mvz, nserec, nsirec, movion, emf, nustrt, ntr,     &
     &idrun0, nvpp, imbalance, monitor
!
! equivalence data to simplify MPI broadcast
      integer, parameter :: lnin3 = 110
      double precision, dimension(lnin3) :: ddin3
      private :: lnin3, ddin3
      equivalence (ddin3(1),idrun), (ddin3(2),idcode), (ddin3(3),indx)
      equivalence (ddin3(4),indy), (ddin3(5),indz), (ddin3(6),mx)
      equivalence (ddin3(7),my), (ddin3(8),mz), (ddin3(9),npx)
      equivalence (ddin3(10),npy), (ddin3(11),npz), (ddin3(12),npxb)
      equivalence (ddin3(13),npyb), (ddin3(14),npzb), (ddin3(15),qme)
      equivalence (ddin3(16),vtx), (ddin3(17),vty), (ddin3(18),vtz)
      equivalence (ddin3(19),vx0), (ddin3(20),vy0), (ddin3(21),vz0)
      equivalence (ddin3(22),vdx), (ddin3(23),vdy), (ddin3(24),vdz)
      equivalence (ddin3(25),vtdx), (ddin3(26),vtdy), (ddin3(27),vtdz)
      equivalence (ddin3(28),relativity), (ddin3(29),ci)
      equivalence (ddin3(30),xtras), (ddin3(31),ndim)
      equivalence (ddin3(32),nvdist), (ddin3(33),treverse)
      equivalence (ddin3(34),tend), (ddin3(35),dt), (ddin3(36),ax)
      equivalence (ddin3(37),ay), (ddin3(38),az), (ddin3(39),nextrand)
      equivalence (ddin3(40),mzf), (ddin3(41),ndprof), (ddin3(42),ampdx)
      equivalence (ddin3(43),scaledx), (ddin3(44),shiftdx)
      equivalence (ddin3(45),ampdy), (ddin3(46),scaledy)
      equivalence (ddin3(47),shiftdy), (ddin3(48),ampdz)
      equivalence (ddin3(49),scaledz), (ddin3(50),shiftdz)
      equivalence (ddin3(51),amodex), (ddin3(52),freq), (ddin3(53),trmp)
      equivalence (ddin3(54),toff), (ddin3(55),el0), (ddin3(56),er0)
      equivalence (ddin3(57),ntw), (ddin3(58),ndw), (ddin3(59),ntde)
      equivalence (ddin3(60),modesxde), (ddin3(61),modesyde)
      equivalence (ddin3(62),modeszde), (ddin3(63),nderec)
      equivalence (ddin3(64),ntp), (ddin3(65),ndp), (ddin3(66),modesxp)
      equivalence (ddin3(67),modesyp), (ddin3(68),modeszp)
      equivalence (ddin3(69),nprec), (ddin3(70),ntel)
      equivalence (ddin3(71),modesxel), (ddin3(72),modesyel)
      equivalence (ddin3(73),modeszel), (ddin3(74),nelrec)
      equivalence (ddin3(75),wmin), (ddin3(76),wmax), (ddin3(77),dw)
      equivalence (ddin3(78),ntfm), (ddin3(79),ndfm), (ddin3(80),npro)
      equivalence (ddin3(81),nferec), (ddin3(82),nfirec)
      equivalence (ddin3(83),ntv), (ddin3(84),ndv), (ddin3(85),nmv)
      equivalence (ddin3(86),nvft), (ddin3(87),nverec)
      equivalence (ddin3(88),nvirec), (ddin3(89),ntt), (ddin3(90),ndt)
      equivalence (ddin3(91),nst), (ddin3(92),nprobt), (ddin3(93),vtsx)
      equivalence (ddin3(94),dvtx), (ddin3(95),ntrec), (ddin3(96),nts)
      equivalence (ddin3(97),nds), (ddin3(98),mvx), (ddin3(99),mvy)
      equivalence (ddin3(100),mvz), (ddin3(101),nserec)
      equivalence (ddin3(102),nsirec), (ddin3(103),movion)
      equivalence (ddin3(104),emf), (ddin3(105),nustrt)
      equivalence (ddin3(106),ntr), (ddin3(107),idrun0)
      equivalence (ddin3(108),nvpp), (ddin3(109),imbalance)
      equivalence (ddin3(110),monitor)
!
! Electromagnetic Namelist
! External Magnetic Field Parameters:
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
!
! Electron Current Diagnostic Parameters:
! ntje = number of time steps between electron current diagnostic
! modesxje/modesyje/modeszje = number of modes in x/y/z to keep for
!                              electron current diagnostic
      integer :: ntje = 0, modesxje = 41, modesyje = 41, modeszje = 41
! njerec = current record number for electron current writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: njerec = -1
!
! Vector Potential Diagnostic Parameters:
! nta = number of time steps between vector potential diagnostic
! nda = (0,1,2,3) = display (nothing,vector potential,spectrum,both)
      integer :: nta = 0, nda = 1
! modesxa/modesya/modesza = number of modes in x/y/z to keep for vector
!                           potential diagnostic
      integer :: modesxa = 41, modesya = 41, modesza = 41
! narec = current record number for vector potential writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: narec = -1
!
! Transverse Efield Diagnostic Parameters:
! ntet = number of time steps between transverse efield diagnostic
! ndet = (0,1,2,3) = display (nothing,transverse efield,spectrum,both)
      integer :: ntet = 0, ndet = 1
! modesxet/modesyet/modeszet = number of modes in x/y/z to keep for
!                              transverse efield diagnostic
      integer :: modesxet = 41, modesyet = 41, modeszet = 41
! netrec = current record number for ltransverse efield writes
!          (0 for beginnning of file, -1 to disable writes)
      integer :: netrec = -1
!
! Magnetic Field Diagnostic Parameters:
! ntb = number of time steps between magnetic field diagnostic
! modesxb/modesyb/modeszb = number of modes in x/y/z to keep for
!                           magnetic field diagnostic
      integer :: ntb = 0, modesxb = 41, modesyb = 41, modeszb = 41
! nbrec = current record number for magnetic field writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: nbrec = -1
!
! Radiative Vector Potential Diagnostic Parameters:
! ntar = number of time steps between radiative vector potential
!        diagnostic
! ndar = (0,1,2,3) = display (nothing,radiative vector potential,
!                             spectrum,both)
      integer :: ntar = 0, ndar = 1
! modesxar/modesyar/modeszar = number of modes in x/y/z to keep for
!                              radiative vector potential diagnostic
      integer :: modesxar = 41, modesyar = 41, modeszar = 41
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
      namelist /input3b/ omx, omy, omz, ntje, modesxje, modesyje,       &
     &modeszje, njerec, nta, nda, modesxa, modesya, modesza, narec,     &
     &ntet, ndet, modesxet, modesyet, modeszet, netrec, ntb, modesxb,   &
     &modesyb, modeszb, nbrec, ntar, ndar, modesxar, modesyar, modeszar,&
     &narrec, wrmin, wrmax, dwr
!
! equivalence data to simplify MPI broadcast
      integer, parameter :: lnin3b = 34
      double precision, dimension(lnin3b) :: ddin3b
      private :: lnin3b, ddin3b
      equivalence (ddin3b(1),omx), (ddin3b(2),omy), (ddin3b(3),omz)
      equivalence (ddin3b(4),ntje), (ddin3b(5),modesxje)
      equivalence (ddin3b(6),modesyje), (ddin3b(7),modeszje)
      equivalence (ddin3b(8),njerec), (ddin3b(9),nta), (ddin3b(10),nda)
      equivalence (ddin3b(11),modesxa), (ddin3b(12),modesya)
      equivalence (ddin3b(13),modesza), (ddin3b(14),narec)
      equivalence (ddin3b(15),ntet), (ddin3b(16),ndet)
      equivalence (ddin3b(17),modesxet), (ddin3b(18),modesyet)
      equivalence (ddin3b(19),modeszet), (ddin3b(20),netrec)
      equivalence (ddin3b(21),ntb), (ddin3b(22),modesxb)
      equivalence (ddin3b(23),modesyb), (ddin3b(24),modeszb)
      equivalence (ddin3b(25),nbrec), (ddin3b(26),ntar)
      equivalence (ddin3b(27),ndar), (ddin3b(28),modesxar)
      equivalence (ddin3b(29),modesyar), (ddin3b(30),modeszar)
      equivalence (ddin3b(31),narrec), (ddin3b(32),wrmin)
      equivalence (ddin3b(33),wrmax), (ddin3b(34),dwr)
!
! Darwin Namelist
! ndc = number of corrections in darwin iteration
      integer :: ndc = 2
!
! define namelist
      namelist /input3d/ ndc
!
! equivalence data to simplify MPI broadcast
      integer, parameter :: lnin3d = 1
      double precision, dimension(lnin3d) :: ddin3d
      private :: lnin3d, ddin3d
      equivalence (ddin3d(1),ndc)
!
! Ion Namelist
! Background Ion Parameters:
! npxi/npyi/npzi = number of background ions distributed in x/y/z
! direction
      integer :: npxi = 384, npyi = 384, npzi = 384
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
! npxbi/npybi/npzbi = number of beam ions distributed in x/y/z direction
      integer :: npxbi = 0, npybi = 0, npzbi = 0
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
! ampdxi/ampdyi/ampdzi = amplitude of ion density compared to uniform
!                        in x/y/z
! scaledxi/scaledyi/scaledzi = scale length for spatial ion coordinate!
!                              in x/y/z
! shiftdxi/shiftdyi/shiftdzi = shift of spatial ion coordinate in x/y/z
      real :: ampdxi = 0.0, scaledxi = 0.0, shiftdxi = 0.0
      real :: ampdyi = 0.0, scaledyi = 0.0, shiftdyi = 0.0
      real :: ampdzi = 0.0, scaledzi = 0.0, shiftdzi = 0.0
!
! Ion Density Diagnostic Parameters:
! ntdi = number of time steps between ion density diagnostic
! nddi = (0,1,2,3) = display (nothing,ion density,spectrum,both)
      integer :: ntdi = 0, nddi = 1
! modesxdi/modesydi/modeszdi = number of modes in x/y/z to keep for ion
!                              density diagnostic
      integer :: modesxdi = 41, modesydi = 41, modeszdi = 41
! ndrec = current record number for ion density writes
!         (0 for beginnning of file, -1 to disable writes)
      integer :: ndirec = -1
!
! Ion Current Diagnostic Parameters:
! ntji = number of time steps between ion current diagnostic
! ndji = (0,1,2,3) = display (nothing,ion current,spectrum,both)
      integer :: ntji = 0, ndji = 1
! modesxji/modesyji/modeszji = number of modes in x/y/z to keep for ion
!                              current diagnostic
      integer :: modesxji = 41, modesyji = 41, modeszji = 41
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
      namelist /ions3/ npxi, npyi, npzi, npxbi, npybi, npzbi, qmi,      &
     &rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0, vdxi, vdyi,   &
     &vdzi, rtempdxi, rtempdyi, rtempdzi, ndprofi, ampdxi, scaledxi,    &
     &shiftdxi, ampdyi, scaledyi, shiftdyi, ampdzi, scaledzi, shiftdzi, &
     &ntdi, nddi, modesxdi, modesydi, modeszdi, ndirec, ntji, ndji,     &
     &modesxji, modesyji, modeszji, njirec, wimin, wimax, dwi
!
! equivalence data to simplify MPI broadcast
      integer, parameter :: lnion3 = 45
      double precision, dimension(lnion3) :: ddion3
      private :: lnion3, ddion3
      equivalence (ddion3(1),npxi), (ddion3(2),npyi), (ddion3(3),npzi)
      equivalence (ddion3(4),npxbi), (ddion3(5),npybi)
      equivalence (ddion3(6),npzbi), (ddion3(7),qmi), (ddion3(8),rmass)
      equivalence (ddion3(9),rtempxi), (ddion3(10),rtempyi)
      equivalence (ddion3(11),rtempzi), (ddion3(12),vxi0)
      equivalence (ddion3(13),vyi0), (ddion3(14),vzi0)
      equivalence (ddion3(15),vdxi), (ddion3(16),vdyi)
      equivalence (ddion3(17),vdzi), (ddion3(18),rtempdxi)
      equivalence (ddion3(19),rtempdyi), (ddion3(20),rtempdzi)
      equivalence (ddion3(21),ndprofi), (ddion3(22),ampdxi)
      equivalence (ddion3(23),scaledxi), (ddion3(24),shiftdxi)
      equivalence (ddion3(25),ampdyi), (ddion3(26),scaledyi)
      equivalence (ddion3(27),shiftdyi), (ddion3(28),ampdzi)
      equivalence (ddion3(29),scaledzi), (ddion3(30),shiftdzi)
      equivalence (ddion3(31),ntdi), (ddion3(32),nddi)
      equivalence (ddion3(33),modesxdi), (ddion3(34),modesydi)
      equivalence (ddion3(35),modeszdi), (ddion3(36),ndirec)
      equivalence (ddion3(37),ntji), (ddion3(38),ndji)
      equivalence (ddion3(39),modesxji), (ddion3(40),modesyji)
      equivalence (ddion3(41),modeszji), (ddion3(42),njirec)
      equivalence (ddion3(43),wimin), (ddion3(44),wimax)
      equivalence (ddion3(45),dwi)
!
! Output namelists
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
!
! Namelist output for electron density diagnostic
! fdename = file name for electron density diagnostic
      character(len=32) :: fdename = 'denek3.0'
! define namelist
      namelist /dene3d/ idrun, indx, indy, indz, ntde, modesxde,        &
     &modesyde, modeszde, nderec, nvpy, nvpz, t0, tend, dt, ceng,       &
     &fdename
!
! Namelist output for potential diagnostic
! fpname = file name for potential diagnostic
      character(len=32) :: fpname = 'potk3.0'
! define namelist
      namelist /pot3d/ idrun, indx, indy, indz, ntp, modesxp, modesyp,  &
     &modeszp, nprec, nvpy, nvpz, t0, tend, dt, ceng, fpname
!
! Namelist output for longitudinal efield diagnostic
! felname = file name for longitudinal efield diagnostic
      character(len=32) :: felname = 'elk3.0'
! define namelist
      namelist /el3d/ idrun, indx, indy, indz, ntel, modesxel, modesyel,&
     &modeszel, ndim, nelrec, nvpy, nvpz, t0, tend, dt, ceng, felname
!
! Namelist output for electron current diagnostic
! fjename = file name for electron current diagnostic
      character(len=32) :: fjename = 'vcurek3.0'
! define namelist
      namelist /vcure3d/ idrun, indx, indy, indz, ntje, modesxje,       &
     &modesyje, modeszje, ndim, omx, omy, omz, ci, njerec, nvpy, nvpz,  &
     &t0, tend, dt, ceng, fjename
!
! Namelist output for vector potential diagnostic
! faname = file name for vector potential diagnostic
      character(len=32) :: faname = 'vpotk3.0'
! define namelist
      namelist /vpot3d/ idrun, indx, indy, indz, nta, modesxa, modesya, &
     &modesza, ndim, omx, omy, omz, ci, narec, nvpy, nvpz, t0, tend, dt,&
     &ceng, faname
!
! Namelist output for transverse efield diagnostic
! fetname = file name for transverse efield diagnostic
      character(len=32) :: fetname = 'etk3.0'
! define namelist
      namelist /et3d/ idrun, indx, indy, indz, ntet, modesxet, modesyet,&
     &modeszet, ndim, netrec, nvpy, nvpz, t0, tend, dt, ceng, fetname
!
! Namelist output for magnetic field diagnostic
! fetname = file name for magnetic field diagnostic
      character(len=32) :: fbname = 'bk3.0'
! define namelist
      namelist /b3d/ idrun, indx, indy, indz, ntb, modesxb, modesyb,    &
     &modeszb, ndim, nbrec, nvpy, nvpz, t0, tend, dt, ceng, fbname
!
! Namelist output for radiative vector potential diagnostic
! farname = file name for vector potential diagnostic
      character(len=32) :: farname = 'vpotrk3.0'
! define namelist
      namelist /vpotr3d/ idrun, indx, indy, indz, ntar, modesxar,       &
     &modesyar, modeszar, ndim, omx, omy, omz, ci, narrec, nvpy, nvpz,  &
     &t0, tend, dt, ceng, farname
!
! Namelist output for fluid moments diagnostic
! nprd = dimension of fluid moment arrays fmse and fmsi
      integer :: nprd = 0
! ffename/ffiname = file name for electron/ion fluid moments diagnostic
      character(len=32) :: ffename = 'fmer3.0', ffiname = 'fmir3.0'
! define namelist
      namelist /fm3d/ idrun, indx, indy, indz, ntfm, npro, ndim, nprd,  &
     &nferec, nfirec, nvpy, nvpz, t0, tend, dt, ceng, ffename, ffiname
!
! Namelist output for velocity-space diagnostic
! nfvd = dimension of velocity distribution arrays fv and fvi
! nfed = dimension of energy distribution arrays fe and fei
      integer :: nfvd = 0, nfed = 0
! fvename/fviname = file name for electron/ion velocity-space diagnostic
      character(len=32) :: fvename = 'fve3.0', fviname = 'fvi3.0'
! define namelist
      namelist /fv3d/ idrun, indx, indy, indz, ntv, nmv, nvft, ndim,    &
     &nfvd, nfed, omx, omy, omz, nverec, nvirec, nvpy, nvpz, t0, tend,  &
     &dt, fvename, fviname
!
! Namelist output for trajectory diagnostic
! ndimp = size of phase space trajectories
      integer :: ndimp = 0
! ftname = file name for trajectory diagnostic
      character(len=32) :: ftname = 'tr3.0'
! define namelist
      namelist /tr3d/ idrun, indx, indy, indz, ntt, ndt, nst, nmv, ndim,&
     &ndimp, nprobt, ntrec, nvpy, nvpz, t0, tend, dt, ftname
!
! Namelist output for phase space diagnostic
! nsxb/nsyb/nszb = number of segments in x/y/z for global velocity distribution
      integer :: nsxb = 0, nsyb= 0, nszb= 0
! fsename/fsiname = file name for electron/ion phase space diagnostic
      character(len=32) :: fsename = 'pse3.0', fsiname = 'psi3.0'
! define namelist
      namelist /ps3d/ idrun, indx, indy, indz, nts, nmv, ndim, nsxb,    &
     &nsyb, nszb, nserec, nsirec, nvpy, nvpz, t0, tend, dt, fsename,    &
     &fsiname
!
! Namelist output for ion density diagnostic
! fdname = file name for ion density diagnostic
      character(len=32) :: fdiname = 'denik3.0'
! define namelist
      namelist /deni3d/ idrun, indx, indy, indz, ntdi, modesxdi,        &
     &modesydi, modeszdi, ndirec, nvpy, nvpz, t0, tend, dt, ceng,       &
     &fdiname 
!
! Namelist output for ion current diagnostic
! fjiname = file name for ion current diagnostic
      character(len=32) :: fjiname = 'vcurik3.0'
! define namelist
      namelist /vcuri3d/ idrun, indx, indy, indz, ntji, modesxji,       &
     &modesyji, modeszji, ndim, omx, omy, omz, ci, njirec, nvpy, nvpz,  &
     &t0, tend, dt, ceng, fjiname
!
      contains
!
      subroutine readnml3(iuin)
! read namelist from unit iuin
      implicit none
      integer, intent(in) :: iuin
! local data
      integer :: ios
      open(unit=iuin,file='input3',form='formatted',status='old')
! read global input parameters
      read (iuin,input3)
! read electromagnetic input parameters
      if ((emf==1).or.(emf==2)) then
         read (iuin,input3b,iostat=ios)
         if (ios /= 0) then
            write (*,*) 'error in reading input3b namelist'
            rewind iuin
         endif
      endif
! read darwin input parameters
      if (emf==2) then
         read (iuin,input3d,iostat=ios)
         if (ios /= 0) then
            write (*,*) 'error in reading input3d namelist'
            rewind iuin
         endif
      endif
! read ion input parameters
      if (movion==1) then
         read (iuin,ions3,iostat=ios)
         if (ios /= 0) write (*,*) 'error in reading ions2 namelist'
      endif
      rewind iuin
      end subroutine
!
      subroutine writnml3(iudm)
! write final diagnostic metafile to unit iudm
      implicit none
      integer, intent(in) :: iudm
! local data
      character(len=10) :: cdrun
      character(len=32) :: fname
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag3.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
! write out global input parameters
      write (iudm,input3)
! write out electromagnetic input parameters
      if ((emf==1).or.(emf==2)) then
         write (iudm,input3b)
      endif
! write out darwin input parameters
      if (emf==2) then
         write (iudm,input3d)
      endif
! electron density diagnostic
      if (ntde > 0) then
         write (iudm,dene3d)
      endif
! potential diagnostic
      if (ntp > 0) then
         write (iudm,pot3d)
      endif
! longitudinal efield diagnostic
      if (ntel > 0) then
         write (iudm,el3d)
      endif
! electron current diagnostic
      if (ntje > 0) then
         ceng = 0.0
         write (iudm,vcure3d)
      endif
! vector potential diagnostic
      if (nta > 0) then
         write (iudm,vpot3d)
      endif
! transverse efield diagnostic
      if (ntet > 0) then
         write (iudm,et3d)
      endif
! magnetic field diagnostic
      if (ntb > 0) then
         write (iudm,b3d)
      endif
! radiative vector potential diagnostic
      if (ntar > 0) then
         write (iudm,vpotr3d)
      endif
! fluid moments diagnostic
      if (ntfm > 0) then
         write (iudm,fm3d)
      endif
! velocity-space diagnostic
      if (ntv > 0) then
         write (iudm,fv3d)
      endif
! trajectory diagnostic
      if (ntt > 0) then
         write (iudm,tr3d)
      endif
! phase space diagnostic
      if (nts > 0) then
         write (iudm,ps3d)
      endif
! ion parameters
      if (movion==1) then
! write out ion input parameters
         write (iudm,ions3)
! ion density diagnostic
         if (ntdi > 0) then
            ceng = 0.0
            write (iudm,deni3d)
         endif
! ion current diagnostic
         if (ntji > 0) then
            ceng = 0.0
            write (iudm,vcuri3d)
         endif
      endif
      end subroutine
!
      subroutine sendnmls3()
! this subroutine broadcasts namelists to other nodes
! namelists are equivalenced to a double precision array
      implicit none
! broadcase namelist input3
      call PPBDCAST(ddin3,lnin3)
! broadcast namelist input3b
      call PPBDCAST(ddin3b,lnin3b)
! broadcast namelist input3d
      call PPBDCAST(ddin3d,lnin3d)
! broadcast namelist ions3
      call PPBDCAST(ddion3,lnion3)
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
