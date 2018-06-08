!-----------------------------------------------------------------------
! 3D Electromagnetic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpbbeps3
      use in3
      use modmpinit3
      use modmppush3
      use modmpbpush3
      use modmpcurd3
      use modmpfield3
      use mpdiag3
      use mppmod3
      use omplib
      use ompplib3
      use mpsimul3
      use mpbsimul3
      implicit none
!
! idimp = number of particle coordinates = 6
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 6, ipbc = 1
! idps = number of partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      integer :: idps = 4, idds =    2
! wke/wki/we = particle kinetic/electrostatic field energy
! wf/wb = transverse electric field/magnetic field
      real :: wke = 0.0, wki = 0.0, we = 0.0, wf = 0.0, wb = 0.0
      real :: zero = 0.0
! plist = (true,false) = list of particles leaving tiles found in push
      logical :: plist = .true.
!
! declare scalars for standard code
      integer :: n
      integer :: nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh
      integer :: nxyzh, nxhyz, mx1, ntime, nloop, isign, ierr
      integer :: ntime0 = 0
      real :: qbme, affp, dth, omt, ws
      real :: qbmi, vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      double precision :: npxyz, npxyzb, np, npxyzi, npxyzbi
      double precision :: npi = 0.0d0
!
! declare scalars for MPI code
      integer :: nvp, idproc, kstrt, npmax, kyp, kzp
      integer :: kxyp, kyzp, nypmx, nzpmx, nypmn, nzpmn
      integer :: npp, nps, myp1, mzp1, mxyzp1, mxzyp1, ntmax
      integer :: npimax, nppi, maxnp
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nppmx1, nbmaxp, ntmaxp, npbmx
      integer :: irc = 0
      integer, dimension(2) :: irc2 = 0
!
! declare scalars for diagnostics
      integer :: js, ks, kyps, kzps, kxyps, kyzps, nmv21, numtp, it
      integer :: nyb, nybmx, nzb, nzbmx, mtw, mtv, mtt, itw, itv, itt
      real :: wef, eci, wk
! default Fortran unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19
      integer :: iude = 10, iup = 11, iuel = 12
      integer :: iur = 17, iur0 = 27, iscr = 99
      integer :: iua = 13, iuet = 14, iub = 15, iuar = 16
      integer :: iuje = 21, iufe = 23, iuve = 25, iut = 28, iuse = 29
      integer :: iudi = 20, iuji = 22, iufi = 24, iuvi = 26, iusi = 30
! dimensions for fourier data
      integer :: modesxpd, modesz2de, modesz2p, modesz2el, modesz2di
      integer :: modesz2ar, modesz2a, modesz2et, modesz2b, modesz2je
      integer :: modesz2ji
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe/qi = electron/ion charge density with guard cells
      real, dimension(:,:,:), allocatable :: qe, qi
! cue/cui = electron/ion current density with guard cells
      real, dimension(:,:,:,:), allocatable :: cue, cui
! fxyze/bxyze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:,:,:), allocatable :: fxyze, bxyze
! qt = scalar charge density field arrays in fourier space
      complex, dimension(:,:,:), allocatable :: qt
! cut = vector current density field arrays in fourier space
      complex, dimension(:,:,:,:), allocatable :: cut
! fxyzt = vector electric field arrays in fourier space
! bxyzt = vector magnetic field arrays in fourier space
      complex, dimension(:,:,:,:), allocatable :: fxyzt, bxyzt
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:,:), allocatable :: exyz, bxyz
! ffc = form factor array for poisson solver
      complex, dimension(:,:,:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
      integer, dimension(2) :: itot
      double precision, dimension(7) :: wtot
      integer, dimension(2) :: mterf
! arrays required for darwin initial fields
      real, dimension(:,:,:,:), allocatable :: amu
      complex, dimension(:,:,:,:), allocatable :: dcut, amut
!
! declare arrays for MPI code:
! edges(1:2) = lower:upper y boundaries of particle partition
! edges(3:4) = back:front z boundaries of particle partition
      real, dimension(:), allocatable  :: edges
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1:2) = lowermost global gridpoint in y/z
      integer, dimension(:), allocatable :: nyzp, noff
! iholep = location of hole left in linear particle arrays
      integer, dimension(:,:), allocatable :: iholep
!
! declare arrays for OpenMP code
! ppart/pparti = tiled electron/ion particle arrays
      real, dimension(:,:,:), allocatable :: ppart, pparti
! kpic/kipic = number of electrons/ions in each tile
      integer, dimension(:), allocatable :: kpic, kipic
! ncl = number of particles departing tile in each direction
      integer, dimension(:,:), allocatable :: ncl
! ihole = location/destination of each particle departing tile
      integer, dimension(:,:,:), allocatable :: ihole
!
! diagnostic arrays
! wt = energy time history array
      real, dimension(:,:), allocatable :: wt
! s = scratch array for energies
      double precision, dimension(:), allocatable :: s
! scratch arrays for scalar field
      complex, dimension(:,:,:), allocatable :: sfieldc
      real, dimension(:,:,:), allocatable :: sfield
! scratch arrays for vector field
      complex, dimension(:,:,:,:), allocatable :: vfieldc
      real, dimension(:,:,:,:), allocatable :: vfield
! denet/denit = store selected fourier modes for electron/ion density
      complex, dimension(:,:,:), allocatable :: denet, denit
! pott = store selected fourier modes for potential
      complex, dimension(:,:,:), allocatable :: pott
! elt = store selected fourier modes for longitudinal efield
      complex, dimension(:,:,:,:), allocatable :: elt
! curet = store selected fourier modes for electron current density
      complex, dimension(:,:,:,:), allocatable :: curet
! curit = store selected fourier modes for ion current density
      complex, dimension(:,:,:,:), allocatable :: curit
! oldcut = previous current density
      complex, dimension(:,:,:,:), allocatable :: oldcut
! vpotr = store selected fourier modes for radiative vector potential
      complex, dimension(:,:,:,:), allocatable :: vpotr
! vpott = store selected fourier modes for vector potential
      complex, dimension(:,:,:,:), allocatable :: vpott
! ett = store selected fourier modes for transverse efield
      complex, dimension(:,:,:,:), allocatable :: ett
! bt = store selected fourier modes for magnetic field
      complex, dimension(:,:,:,:), allocatable :: bt
! fmse/fmsi = electron/ion fluid moments
      real, dimension(:,:,:,:), allocatable :: fmse, fmsi
! fv/fvi = global electron/ion velocity distribution functions
      real, dimension(:,:), allocatable :: fv, fvi
! fe/fei = global electron/ion energy distribution functions
      real, dimension(:,:), allocatable :: fe, fei
! sfv = electron/ion velocity distribution functions in tile
      real, dimension(:,:,:), allocatable :: sfv
! fvm/fvmi = electron/ion vdrift, vth, entropy for global distribution
      real, dimension(:,:), allocatable :: fvm, fvmi
! fvtm/fvtmi = time history of electron/ion vdrift, vth, and entropy
!              for global distribution
      real, dimension(:,:,:), allocatable :: fvtm, fvtmi
! tedges(1:2) = back:front z boundaries of particle tags
! tedges(3:4) = lower:upper y boundaries of particle tags
      real, dimension(:), allocatable  :: tedges
! iprobt = scratch array 
      integer, dimension(:), allocatable :: iprobt
! partt = particle trajectories tracked
      real, dimension(:,:), allocatable :: partt
! fvtp = velocity distribution function for test particles
! fvmtp = vdrift, vth, and entropy for test particles
      real, dimension(:,:), allocatable :: fvtp, fvmtp, fetp
! partd = trajectory time history array
      real, dimension(:,:,:), allocatable :: partd
! fvs/fvsi = global electron/ion phase space distribution functions
      real, dimension(:,:,:,:,:), allocatable :: fvs, fvsi
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tguard = 0.0, tfield = 0.0, tdjpost = 0.0
      real :: tpush = 0.0, tsort = 0.0, tmov = 0.0, tfmov = 0.0
      real :: tdiag = 0.0
      real, dimension(2) :: tfft = 0.0
      double precision :: dtime
!
! start timing initialization
      call dtimer(dtime,itime,-1)
!   
! nvp = number of MPI ranks
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
!
      irc = 0
! nvpp = number of shared memory nodes (0=default)
      nvpp = 0
!     if (kstrt==1) then
!        write (*,*) 'enter number of nodes:'
!        read (5,*) nvpp
!     endif
! initialize for shared memory parallel processing
      call INIT_OMP(nvpp)
!   
! read namelists
      if (kstrt==1) then
! override default input data
         emf = 1
         relativity = 1
         call readnml3(iuin)
! override input data
         idcode = 2
! create string from idrun
         write (cdrun,'(i10)') idrun
         cdrun = adjustl(cdrun)
! text output file
         fname = 'output3.'//cdrun
         open(unit=iuot,file=trim(fname),form='formatted',              &
     &status='replace')
      endif
!
! broadcast namelists to other nodes
      call sendnmls3()
!
! increase number of coordinates for particle tag
      if (ntt > 0) idimp = idimp + 1
!
! initialize scalars for standard code
! np = total number of particles in simulation
      npxyz =  dble(npx)*dble(npy)*dble(npz)
      npxyzb =  dble(npxb)*dble(npyb)*dble(npzb)
      np = npxyz + npxyzb
! npi = total number of ions in simulation
      if (movion > 0) then
         npxyzi = dble(npxi)*dble(npyi)*dble(npzi)
         npxyzbi = dble(npxbi)*dble(npybi)*dble(npzbi)
         npi = npxyzi + npxyzbi
      endif
! nx/ny/nz = number of grid points in x/y direction
      nx = 2**indx; ny = 2**indy; nz = 2**indz
      nxh = nx/2; nyh = max(1,ny/2); nzh = max(1,nz/2)
      nxe = nx + 2; nye = ny + 2; nze = nz + 2
      nxeh = nxe/2
      nxyzh = max(nx,ny,nz)/2; nxhyz = max(nxh,ny,nz)
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = dble(nx)*dble(ny)*dble(nz)/np
      if (movion==1) then
         qbmi = qmi/rmass
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
         vtdyi = vtdy/sqrt(rmass*rtempdyi)
         vtdzi = vtdz/sqrt(rmass*rtempdzi)
      endif
      dth = 0.0
!
! obtain 2D partition (nvpy,nvpz) from nvp:
! nvpy/nvpz = number of processors in y/z
      call FCOMP32(nvp,nx,ny,nz,nvpy,nvpz,ierr)
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'FCOMP32 error: nvp,nvpy,nvpz=', nvp, nvpy, nvpz
         endif
         call PPEXIT(); stop
      endif
!
! initialize data for MPI code
      allocate(edges(idps),nyzp(idds),noff(idds))
! calculate partition variables:
! edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn
! edges(1:2) = lower:upper boundary of particle partition in y
! edges(3:4) = back:front boundary of particle partition in z
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1:2) = lowermost global gridpoint in y/z in particle partition
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! nypmn = minimum value of nyzp(1)
! nzpmn = minimum value of nyzp(2)
! find new 2d partition for uniform density distribution
!     call mpdcomp3(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,&
!    &nvpy,nvpz)
! find new 2d partition from initial analytic distribution function
      call mpfedges3(edges,nyzp,noff,ampdy,scaledy,shiftdy,ampdz,scaledz&
     &,shiftdz,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy,nvpz,ipbc,ndprof&
     &,ierr)
      if (ierr /= 0) then
         call PPEXIT(); stop
      endif
!
! check for unimplemented features
      if ((plist).and.(ipbc.ne.1)) then
         if (kstrt==1) then
            write (*,*) 'ipbc /= 1 and plist = .true. not yet supported'
            write (*,*) 'plist reset to .false.'
            plist = .false.
         endif
      endif
!
! initialize additional scalars for MPI code
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvpy + 1
! kzp = number of complex grids in each field partition in z direction
      kzp = (nz - 1)/nvpz + 1
! kxyp = number of complex grids in each field partition in x direction
! in transposed data
      kxyp = (nxh - 1)/nvpy + 1
! kyzp = number of complex grids in each field partition in y direction,
! in transposed data
      kyzp = (ny - 1)/nvpz + 1
! npmax/npimax = maximum number of electrons/ions in each partition
      npmax = (np/nvp)*1.25; npimax = (npi/nvp)*1.25
      maxnp = max(npmax,npimax)
! myp1/mzp1 = number of tiles in y/z direction
      myp1 = (nyzp(1) - 1)/my + 1; mzp1 = (nyzp(2) - 1)/mz + 1
! mxzyp1 = mx1*max(max(mzp1),max(myp1))
      mxzyp1 = mx1*max((nzpmx-2)/mz+1,(nypmx-2)/my+1)
      mxyzp1 = mx1*myp1*mzp1
! mterf = number of shifts required by field manager in y/z (0=search)
      mterf = 0
! ntmax = size of iholep buffer for particles leaving node
      ntmax = 0.2*npmax
!
! allocate data for standard code
      allocate(part(idimp,maxnp))
      allocate(qe(nxe,nypmx,nzpmx),qi(nxe,nypmx,nzpmx))
      allocate(cue(ndim,nxe,nypmx,nzpmx))
      allocate(fxyze(ndim,nxe,nypmx,nzpmx),bxyze(ndim,nxe,nypmx,nzpmx))
      allocate(qt(nze,kxyp,kyzp),cut(ndim,nze,kxyp,kyzp))
      allocate(fxyzt(ndim,nze,kxyp,kyzp),bxyzt(ndim,nze,kxyp,kyzp))
      allocate(exyz(ndim,nze,kxyp,kyzp),bxyz(ndim,nze,kxyp,kyzp))
      allocate(ffc(nzh,kxyp,kyzp),mixup(nxhyz),sct(nxyzh))
      allocate(kpic(mxyzp1),ncl(26,mxyzp1),iholep(ntmax+1,2))
      if (movion==1) allocate(kipic(mxyzp1),cui(ndim,nxe,nypmx,nzpmx))
!
! prepare fft tables
      call mpfft3_init(mixup,sct,indx,indy,indz)
! calculate form factors
      call mppois3_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,nvpy,nvpz)
! initialize different ensemble of random numbers
      if (nextrand > 0) call mnextran3(nextrand,ndim,npmax+npimax)
!
! open restart files
      if (kstrt==1) call open_restart3(iur,iur0,cdrun)
!
! new start
      if (nustrt==1) then
! initialize electrons
         nps = 1
         npp = 0
! background electrons
         if (npxyz > 0.0d0) then
! calculates initial electron co-ordinates with uniform density and
! velocities or momenta
!           call wmpdistr3(part,edges,npp,vtx,vty,vtz,vx0,vy0,vz0,ci,npx&
!    &,npy,npz,nx,ny,nz,kstrt,ipbc,relativity,ierr)
! calculates initial electron co-ordinates with various density profiles
            call mpfdistr3(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,&
     &shiftdy,ampdz,scaledz,shiftdz,npx,npy,npz,nx,ny,nz,kstrt,nvpy,nvpz&
     &,ipbc,ndprof,ierr)
! initialize electron velocities or momenta
            if (ierr==0) then
               call wmpvdistr3(part,nps,npp,vtx,vty,vtz,vx0,vy0,vz0,ci, &
     &npx,npy,npz,kstrt,nvpy,nvpz,relativity,ierr)
            endif
! check for background electron initialization error
            if (ierr /= 0) then
               call PPEXIT(); stop
            endif
         endif
! beam electrons
         if (npxyzb > 0.0d0) then
            nps = npp + 1
! calculates initial electron co-ordinates with uniform density and
! velocities or momenta
!           call wmpdistr3(part,edges,npp,vtdx,vtdy,vtdz,vdx,vdy,vdz,ci,&
!    &npxb,npyb,npzb,nx,ny,nz,kstrt,ipbc,relativity,ierr)
! calculates initial electron co-ordinates with various density profiles
            call mpfdistr3(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,&
     &shiftdy,ampdz,scaledz,shiftdz,npxb,npyb,npzb,nx,ny,nz,kstrt,nvpy, &
     &nvpz,ipbc,ndprof,ierr)
! initialize electron velocities or momenta
            if (ierr==0) then
               call wmpvdistr3(part,nps,npp,vtdx,vtdy,vtdz,vdx,vdy,vdz, &
     &ci,npxb,npyb,npzb,kstrt,nvpy,nvpz,relativity,ierr)
            endif
! check for beam electron initialization error
            if (ierr /= 0) then
               call PPEXIT(); stop
            endif
         endif
!
! check if any electrons are in the wrong node
         call mpfholes3(part,edges,npp,iholep,1)
! iholep overflow
         if (iholep(1,1) < 0) then
            ntmax = -iholep(1,1)
            ntmax = 1.5*ntmax
            deallocate(iholep)
            if (kstrt==1) then
               write (*,*) 'reallocating electron iholep: ntmax=', ntmax
            endif
            allocate(iholep(ntmax+1,2))
            call mpfholes3(part,edges,npp,iholep,1)
            if (iholep(1,1) < 0) then
               if (kstrt==1) write (*,*) 'iholep overflow: ntmax=',ntmax
               call PPEXIT(); stop
            endif
         endif
! more electrons to correct node
         call ipmove3(part,edges,npp,iholep,ny,nz,tmov,kstrt,nvpy,nvpz,1&
     &,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
!
! find number of electrons in each of mx, my, mz tiles:
! updates kpic, nppmx
         call mpdblkp3(part,kpic,npp,noff,nppmx,mx,my,mz,mx1,myp1,irc)
!
! allocate vector electron data
         nppmx0 = (1.0 + xtras)*nppmx
         ntmaxp = xtras*nppmx
         npbmx = xtras*nppmx
         nbmaxp = 0.125*mxzyp1*npbmx
         allocate(ppart(idimp,nppmx0,mxyzp1))
         allocate(ihole(2,ntmaxp+1,mxyzp1))
! copy ordered electron data for OpenMP: updates ppart and kpic
         call mpmovin3(part,ppart,kpic,npp,noff,mx,my,mz,mx1,myp1,irc)
!
! sanity check for electrons
         call mpcheck3(ppart,kpic,noff,nyzp,nx,mx,my,mz,mx1,myp1,irc)
!
! initialize background charge density: updates qi
         if (movion==0) then
            qi = 0.0
            call mppost3(ppart,qi,kpic,noff,-qme,tdpost,mx,my,mz,mx1,   &
     &myp1)
            call wmpaguard3(qi,nyzp,tguard,nx,kstrt,nvpy,nvpz)
         endif
!
! initialize ions
         if (movion==1) then
            cui = 0.0
            nps = 1
            nppi = 0
! background ions
            if (npxyzi > 0.0d0) then
! calculates initial ion co-ordinates with uniform density and
! velocities or momenta
!              call wmpdistr3(part,edges,nppi,vtxi,vtyi,vtzi,vxi0,vyi0, &
!    &vzi0,ci,npxi,npyi,npzi,nx,ny,nz,kstrt,ipbc,relativity,ierr)
! calculates initial ion co-ordinates with various density profiles
               call mpfdistr3(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,&
     &scaledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxi,npyi,npzi,nx,ny,nz&
     &,kstrt,nvpy,nvpz,ipbc,ndprofi,ierr)
! initialize ion velocities or momenta
               if (ierr==0) then
                  call wmpvdistr3(part,nps,nppi,vtxi,vtyi,vtzi,vxi0,vyi0&
     &,vzi0,ci,npxi,npyi,npzi,kstrt,nvpy,nvpz,relativity,ierr)
               endif
! check for background ion initialization error
               if (ierr /= 0) then
                  call PPEXIT(); stop
               endif
            endif
! beam ions
            if (npxyzbi > 0.0d0) then
               nps = nppi + 1
! calculates initial ion co-ordinates with uniform density and
! velocities or momenta
!              call wmpdistr3(part,edges,nppi,vtdxi,vtdyi,vtdzi,vdxi,   &
!    &vdyi,vdzi,ci,npxbi,npybi,npzbi,nx,ny,nz,kstrt,ipbc,relativity,ierr&
!    &)
! calculates initial ion co-ordinates with various density profiles
               call mpfdistr3(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,&
     &scaledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxbi,npybi,npzbi,nx,ny&
     &,nz,kstrt,nvpy,nvpz,ipbc,ndprofi,ierr)
! initialize ion velocities or momenta
               if (ierr==0) then
                  call wmpvdistr3(part,nps,nppi,vtdxi,vtdyi,vtdzi,vdxi, &
     &vdyi,vdzi,ci,npxbi,npybi,npzbi,kstrt,nvpy,nvpz,relativity,ierr)
               endif
! check for background ion initialization error
               if (ierr /= 0) then
                  call PPEXIT(); stop
               endif
            endif
!
! check if any ions are in the wrong node
            call mpfholes3(part,edges,nppi,iholep,1)
! iholep overflow
            if (iholep(1,1) < 0) then
               ntmax = -iholep(1,1)
               ntmax = 1.5*ntmax
               deallocate(iholep)
               if (kstrt==1) then
                  write (*,*) 'reallocating ion iholep: ntmax=', ntmax
               endif
               allocate(iholep(ntmax+1,2))
               call mpfholes3(part,edges,nppi,iholep,1)
               if (iholep(1,1) < 0) then
                  if (kstrt==1) then
                     write (*,*) 'iholep overflow: ntmax=', ntmax
                  endif
                  call PPEXIT(); stop
               endif
            endif
! more ions to correct node
            call ipmove3(part,edges,nppi,iholep,ny,nz,tmov,kstrt,nvpy,  &
     &nvpz,1,ierr)
            if (ierr /= 0) then
               call PPEXIT(); stop
            endif
!
! find number of ions in each of mx, my, mz tiles:
! updates kipic, nppmx
            call mpdblkp3(part,kipic,nppi,noff,nppmx,mx,my,mz,mx1,myp1, &
     &irc)
!
! allocate vector ion data
            nppmx1 = (1.0 + xtras)*nppmx
            allocate(pparti(idimp,nppmx1,mxyzp1))
            if (.not.allocated(ihole)) then
               ntmaxp = xtras*nppmx
               npbmx = xtras*nppmx
               nbmaxp = 0.25*mxzyp1*npbmx
               allocate(ihole(2,ntmaxp+1,mxyzp1))
            endif
! copy ordered ion data for OpenMP: updates pparti and kipic
            call mpmovin3(part,pparti,kipic,nppi,noff,mx,my,mz,mx1,myp1,&
     &irc)
!
! sanity check for ions
            call mpcheck3(pparti,kipic,noff,nyzp,nx,mx,my,mz,mx1,myp1,  &
     &irc)
         endif
!
! initialize transverse electromagnetic fields
         exyz = cmplx(0.0,0.0)
         bxyz = cmplx(0.0,0.0)
         cut = cmplx(0.0,0.0)
!
! restart to continue a run which was interrupted
      else if (nustrt==2) then
         write (*,*) 'nustrt = 2 not yet supported'
         call PPEXIT(); stop
!        if ((ntime+ntime0)> 0) dth = 0.5*dt
!        nstart = ntime + 1
! start a new run with data from a previous run
      else if (nustrt==0) then
! read in basic restart file for electrostatic code
         if (movion==1) cui = 0.0
! read first part of data:
! updates ntime, ntime0, part, npp, kpic, nppmx, ierr
         call bread_restart3a(part,kpic,noff,tdiag,kstrt,iur0,iscr,ntime&
     &,ntime0,npp,nppmx,mx1,myp1,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
! allocate vector electron data
         nppmx0 = (1.0 + xtras)*nppmx
         ntmaxp = xtras*nppmx
         npbmx = xtras*nppmx
         nbmaxp = 0.25*mxzyp1*npbmx
         allocate(ppart(idimp,nppmx0,mxyzp1))
         allocate(ihole(2,ntmaxp+1,mxyzp1))
! read second part of data:
! updates ppart, kpic, part, nppi, kipic, nppmx, ierr
         call bread_restart3b(part,ppart,kpic,kipic,noff,nyzp,tdiag,    &
     &kstrt,iur0,iscr,npp,nppi,nppmx,nx,mx1,myp1,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
! allocate vector ion data
         if (movion==1) then
            nppmx1 = (1.0 + xtras)*nppmx
            allocate(pparti(idimp,nppmx1,mxyzp1))
            if (.not.allocated(ihole)) then
               ntmaxp = xtras*nppmx
               npbmx = xtras*nppmx
               nbmaxp = 0.25*mxzyp1*npbmx
               allocate(ihole(2,ntmaxp+1,mxyzp1))
            endif
         endif
! read third part of data: updates pparti, kipic, qi, ierr
         call bread_restart3c(part,pparti,kipic,noff,nyzp,qi,tdiag,kstrt&
     &,iur0,ntime,ntime0,nppi,nx,mx1,myp1,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
! read in basic restart file for electromagnetic code:
! updates exyz, bxyz
         call bread_restart3(exyz,bxyz,tdiag,kstrt,iur0,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
         if ((ntime+ntime0)> 0) dth = 0.5*dt
      endif
!
! set magnitude of external magnetic field
      omt = sqrt(omx*omx + omy*omy + omz*omz)
!
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
! kyps/kzps = actual grids used in field partitions in y/z
      kyps = min(kyp,max(0,ny-kyp*js))
      kzps = min(kzp,max(0,nz-kzp*ks))
! kxyps/kyzps = actual grids in each transposed field partition in x/y
      kxyps = min(kxyp,max(0,nxh-kxyp*js))
      kyzps = min(kyzp,max(0,ny-kyzp*ks))
!
! allocate diagnostic arrays
! reverse simulation at end back to start
      if (treverse==1) nloop = 2*nloop
!
! energy time history
      if (ntw > 0) then
         mtw = (nloop - 1)/ntw + 1; itw = 0
         allocate(wt(mtw,7),s(7))
         wt = 0.0; s = 0.0d0
      endif
!
! allocate scratch arrays for scalar fields
      if ((ntde > 0).or.(ntp > 0).or.(ntdi > 0)) then
         allocate(sfieldc(nze,kxyp,kyzp),sfield(nxe,nypmx,nzpmx))
      endif
!
! allocate scratch arrays for vector fields
      if ((ntel>0).or.(ntje>0).or.(nta>0).or.(ntet>0).or.(ntb>0).or.    &
     &(ntar>0).or.(ntji>0)) then
         allocate(vfieldc(ndim,nze,kxyp,kyzp))
         allocate(vfield(ndim,nxe,nypmx,nzpmx))
      endif
!
! initialize electron density diagnostic
      if (ntde > 0) then
         modesxde = min(modesxde,nxh+1)
         modesyde = min(modesyde,nyh+1)
         modeszde = min(modeszde,nzh+1)
         modesz2de = min(2*modeszde-1,nz)
         modesxpd = min(modesxde,kxyp)
         if (modesxde==(nxh+1)) modesxpd = modesxpd + 1
         allocate(denet(modesz2de,modesxpd,kyzp))
! open file for real data: updates nderec and possibly iude
         if (kstrt==1) then
            fdename = 'dener3.'//cdrun
            if (nderec==0) then
               call dafopen3(sfield,nx,kyp,kzp,iude,nderec,trim(fdename)&
     &)
            endif
         endif
      endif
!
! initialize ion density diagnostic
      if (movion==1) then
         if (ntdi > 0) then
            modesxdi = min(modesxdi,nxh+1)
            modesydi = min(modesydi,nyh+1)
            modeszdi = min(modeszdi,nzh+1)
            modesz2di = min(2*modeszdi-1,nz)
            modesxpd = min(modesxdi,kxyp)
            if (modesxdi==(nxh+1)) modesxpd = modesxpd + 1
            allocate(denit(modesz2di,modesxpd,kyzp))
! open file for real data: updates ndirec and possibly iudi
            if (kstrt==1) then
               fdiname = 'denir3.'//cdrun
               if (ndirec==0) then
                  call dafopen3(sfield,nx,kyp,kzp,iudi,ndirec,          &
     &trim(fdiname))
               endif
            endif
         endif
      endif
!
! initialize potential diagnostic
      if (ntp > 0) then
         modesxp = min(modesxp,nxh+1)
         modesyp = min(modesyp,nyh+1)
         modeszp = min(modeszp,nzh+1)
         modesz2p = min(2*modeszp-1,nz)
         modesxpd = min(modesxp,kxyp)
         if (modesxp==(nxh+1)) modesxpd = modesxpd + 1
         allocate(pott(modesz2p,modesxpd,kyzp))
! open file for real data: updates nprec and possibly iup
         if (kstrt==1) then
            fpname = 'potr3.'//cdrun
            if (nprec==0) then
               call dafopen3(sfield,nx,kyp,kzp,iup,nprec,trim(fpname))
            endif
         endif
      endif
!
! initialize longitudinal efield diagnostic
      if (ntel > 0) then
         modesxel = min(modesxel,nxh+1)
         modesyel = min(modesyel,nyh+1)
         modeszel = min(modeszel,nzh+1)
         modesz2el = min(2*modeszel-1,nz)
         modesxpd = min(modesxel,kxyp)
         if (modesxel==(nxh+1)) modesxpd = modesxpd + 1
         allocate(elt(ndim,modesz2el,modesxpd,kyzp))
! open file for real data: updates nelrec and possibly iuel
         if (kstrt==1) then
            felname = 'elr3.'//cdrun
            if (nelrec==0) then
               call dafopenv3(vfield,nx,kyp,kzp,iuel,nelrec,            &
     &trim(felname))
            endif
         endif
      endif
!
! initialize electron current density diagnostic
      if (ntje > 0) then
         modesxje = min(modesxje,nxh+1)
         modesyje = min(modesyje,nyh+1)
         modeszje = min(modeszje,nzh+1)
         modesz2je = min(2*modeszje-1,nz)
         modesxpd = min(modesxje,kxyp)
         if (modesxje==(nxh+1)) modesxpd = modesxpd + 1
         allocate(curet(ndim,modesz2je,modesxpd,kyzp))
! open file for real data: updates njerec and possibly iuje
         if (kstrt==1) then
            fjename = 'curer3.'//cdrun
            if (njerec==0) then
               call dafopenv3(vfield,nx,kyp,kzp,iuje,njerec,            &
     &trim(fjename))
            endif
         endif
      endif
!
! initialize ion current density diagnostic
      if (movion==1) then
         if (ntji > 0) then
            modesxji = min(modesxji,nxh+1)
            modesyji = min(modesyji,nyh+1)
            modeszji = min(modeszji,nzh+1)
            modesz2ji = min(2*modeszji-1,nz)
            modesxpd = min(modesxji,kxyp)
            if (modesxji==(nxh+1)) modesxpd = modesxpd + 1
            allocate(curit(ndim,modesz2ji,modesxpd,kyzp))
! open file for real data: updates njirec and possibly iuji
            if (kstrt==1) then
               fjiname = 'curir3.'//cdrun
               if (njirec==0) then
                  call dafopenv3(vfield,nx,kyp,kzp,iuji,njirec,         &
     &trim(fjiname))
               endif
            endif
         endif
      endif
!
! initialize radiative vector potential diagnostic
      if (ntar > 0) then
         modesxar = min(modesxar,nxh+1)
         modesyar = min(modesyar,nyh+1)
         modeszar = min(modeszar,nzh+1)
         modesz2ar = min(2*modeszar-1,nz)
         modesxpd = min(modesxar,kxyp)
         if (modesxar==(nxh+1)) modesxpd = modesxpd + 1
         allocate(vpotr(ndim,modesz2ar,modesxpd,kyzp))
         allocate(oldcut(ndim,nze,kxyp,kyzp))
! open file for real data: updates narrec and possibly iuar
         if (kstrt==1) then
            farname = 'vpotrr3.'//cdrun
            if (narrec==0) then
               call dafopenv3(vfield,nx,kyp,kzp,iuar,narrec,            &
     &trim(farname))
            endif
         endif
      endif
!
! initialize vector potential diagnostic
      if (nta > 0) then
         modesxa = min(modesxa,nxh+1)
         modesya = min(modesya,nyh+1)
         modesza = min(modesza,nzh+1)
         modesz2a = min(2*modesza-1,nz)
         modesxpd = min(modesxa,kxyp)
         if (modesxa==(nxh+1)) modesxpd = modesxpd + 1
         allocate(vpott(ndim,modesz2a,modesxpd,kyzp))
! open file for real data: updates narec and possibly iua
         if (kstrt==1) then
            faname = 'vpotr3.'//cdrun
            if (narec==0) then
               call dafopenv3(vfield,nx,kyp,kzp,iua,narec,trim(faname))
            endif
         endif
      endif
!
! initialize transverse efield diagnostic
      if (ntet > 0) then
         modesxet = min(modesxet,nxh+1)
         modesyet = min(modesyet,nyh+1)
         modeszet = min(modeszet,nzh+1)
         modesz2et = min(2*modeszet-1,nz)
         modesxpd = min(modesxet,kxyp)
         if (modesxet==(nxh+1)) modesxpd = modesxpd + 1
         allocate(ett(ndim,modesz2et,modesxpd,kyzp))
! open file for real data: updates naetrec and possibly iuet
         if (kstrt==1) then
            fetname = 'etr3.'//cdrun
            if (netrec==0) then
               call dafopenv3(vfield,nx,kyp,kzp,iuet,netrec,            &
     &trim(fetname))
            endif
         endif
      endif
!
! initialize magnetic field diagnostic
      if (ntb > 0) then
         modesxb = min(modesxb,nxh+1)
         modesyb = min(modesyb,nyh+1)
         modeszb = min(modeszb,nzh+1)
         modesz2b = min(2*modeszb-1,nz)
         modesxpd = min(modesxb,kxyp)
         if (modesxb==(nxh+1)) modesxpd = modesxpd + 1
         allocate(bt(ndim,modesz2b,modesxpd,kyzp))
! open file for real data: updates nbrec and possibly iub
         if (kstrt==1) then
            fbname = 'br3.'//cdrun
            if (nbrec==0) then
               call dafopenv3(vfield,nx,kyp,kzp,iub,nbrec,trim(fbname))
            endif
         endif
      endif
!
! initialize fluid moments diagnostic
      if (ntfm > 0) then
         nprd = 0
         if (npro==1) then
            nprd = 1
         else if (npro==2) then
            nprd = 4
         else if (npro==3) then
            nprd = 10
         else if (npro==4) then
            nprd = 14
         endif
! electron moments
         if ((ndfm==1).or.(ndfm==3)) then
            allocate(fmse(nprd,nxe,nypmx,nzpmx))
! open file for real data: updates nferec and possibly iufe
            if (kstrt==1) then
               ffename = 'fmer3.'//cdrun
               if (nferec==0) then
                  call dafopenv3(fmse,nx,kyp,kzp,iufe,nferec,           &
     &trim(ffename))
               endif
            endif
         endif
! ion moments
         if (movion==1) then
            if ((ndfm==2).or.(ndfm==3)) then
               allocate(fmsi(nprd,nxe,nypmx,nzpmx))
! open file for real data: updates nfirec and possibly iufi
               if (kstrt==1) then
                  ffiname = 'fmir3.'//cdrun
                  if (nfirec==0) then
                     call dafopenv3(fmsi,nx,kyp,kzp,iufi,nfirec,           &
     &trim(ffiname))
                  endif
               endif
            endif
         endif
      endif
!
! initialize velocity-space diagnostic
      if (ntv > 0) then
         nfvd = 0; nfed = 0
         if ((nvft==1).or.(nvft==3)) then
            nfvd = ndim
         else if ((nvft==4).or.(nvft==5)) then
            nfvd = 2
         endif
         if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
            nfed = 1
         endif
         nmv21 = 2*nmv + 1
         mtv = (nloop - 1)/ntv + 1; itv = 0
         eci = ci; if (relativity==0) eci = 0.0
         wk = 0.0
! electron velocity diagnostic
         if ((ndv==1).or.(ndv==3)) then
! estimate maximum velocity or momentum
            ws = 0.0
            if (npxyz.gt.0.0d0) then
               ws = 4.0*vtx+abs(vx0)
               ws = max(ws,4.0*vty+abs(vy0))
               ws = max(ws,4.0*vtz+abs(vz0))
            endif
            if (npxyzb.gt.0.0d0) then
               ws = max(ws,4.0*vtdx+abs(vdx))
               ws = max(ws,4.0*vtdy+abs(vdy))
               ws = max(ws,4.0*vtdz+abs(vdz))
            endif
            allocate(fvm(ndim,3),sfv(nmv21+1,ndim,mxyzp1))
            allocate(fv(nmv21+1,nfvd),fe(nmv21+1,nfed))
            fvm = 0.0
! open file for electron velocity data: updates nverec and possibly iuve
            if (kstrt==1) then
               fvename = 'fve3.'//cdrun
               if (nverec==0) then
                  call dafopenfv3(fvm,fv,fe,wk,iuve,nverec,             &
     &trim(fvename))
               endif
            endif
! cartesian distribution
            if ((nvft==1).or.(nvft==3)) then
               allocate(fvtm(mtv,ndim,3))
               fvtm = 0.0
! set velocity or momentum scale
               fv(nmv21+1,:) = 2.0*ws
            endif
! cylindrical distribution
            if ((nvft==4).or.(nvft==5)) then
! set velocity or momentum scale
               fv(nmv21+1,:) = 2.0*ws
            endif
! energy distribution
            if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
! set energy scale for electrons
               ws = ws*ws
               fe(nmv21+1,:) = ws/(1.0 + sqrt(1.0 + ws*eci*eci))
            endif
         endif
! ion velocity diagnostic
         if (movion==1) then
            if ((ndv==2).or.(ndv==3)) then
               allocate(fvmi(ndim,3))
               allocate(fvi(nmv21+1,nfvd),fei(nmv21+1,nfed))
! estimate maximum ion velocity or momentum
               ws = 0.0
               if (npxyzi.gt.0.0d0) then
                   ws = 4.0*vtxi+abs(vxi0)
                   ws = max(ws,4.0*vtyi+abs(vyi0))
                   ws = max(ws,4.0*vtzi+abs(vzi0))
               endif
               if (npxyzbi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
                  ws = max(ws,4.0*vtdzi+abs(vdzi))
               endif
! open file for ion velocity data: updates nvirec and possibly iuvi
               if (kstrt==1) then
                  fviname = 'fvi3.'//cdrun
                  if (nvirec==0) then
                     call dafopenfv3(fvmi,fvi,fei,wk,iuvi,nvirec,       &
     &trim(fviname))
                  endif
               endif
! cartesian distribution
               if ((nvft==1).or.(nvft==3)) then
                  allocate(fvtmi(mtv,ndim,3))
                  fvtmi = 0.0
! set velocity or momentum scale
                  fvi(nmv21+1,:) = 2.0*ws
               endif
! cylindrical distribution
               if ((nvft==4).or.(nvft==5)) then
! set velocity or momentum scale
                  fvi(nmv21+1,:) = 2.0*ws
               endif
! energy distribution
               if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
! set energy scale for ions
                  ws = ws*ws
                  fei(nmv21+1,1) = ws/(1.0 + sqrt(1.0 + ws*eci*eci))
               endif
            endif
         endif
      endif
!
! initialize trajectory diagnostic
      if (ntt > 0) then
         if ((ndt==2).and.(movion==0)) ndt = 0
         if ((ndt==1).or.(ndt==2)) then
            allocate(iprobt(nprobt),tedges(idps))
         endif
! electron trajectories
         if (ndt==1) then
! set particle tags: updates nprobt, tedges, ppart and possibly iprobt
            call psetptraj3(ppart,tedges,kpic,iprobt,kstrt,nst,nvpy,nvpz&
     &,vtx,vtsx,dvtx,np,nprobt,irc)
! estimate maximum electron velocity or momentum
            if (nst==3) then
               ws = 0.0
               if (npxyz.gt.0.0d0) then
                  ws = 4.0*vtx+abs(vx0)
                  ws = max(ws,4.0*vty+abs(vy0))
                  ws = max(ws,4.0*vtz+abs(vz0))
               endif
               if (npxyzb.gt.0.0d0) then
                  ws = max(ws,4.0*vtdx+abs(vdx))
                  ws = max(ws,4.0*vtdy+abs(vdy))
                  ws = max(ws,4.0*vtdz+abs(vdz))
               endif
            endif
! ion trajectories
         else if (ndt==2) then
! set particle tags: updates nprobt, tedges, pparti and possibly iprobt
            call psetptraj3(pparti,tedges,kipic,iprobt,kstrt,nst,nvpy,  &
     &nvpz,vtxi,vtsx,dvtx,npi,nprobt,irc)
! estimate maximum ion velocity or momentum
            if (nst==3) then
               ws = 0.0
               if (npxyzi.gt.0.0d0) then
                  ws = 4.0*vtxi+abs(vxi0)
                  ws = max(ws,4.0*vtyi+abs(vyi0))
                  ws = max(ws,4.0*vtzi+abs(vzi0))
               endif
               if (npxyzbi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
                  ws = max(ws,4.0*vtzi+abs(vzi0))
               endif
            endif
         endif
! electron or ion trajectories
         if ((ndt==1).or.(ndt==2)) then
            if (nprobt.gt.16777215) then
               write(*,*) 'nprobt overflow = ', nprobt
               call PPEXIT(); stop
            endif
            ndimp = idimp
            allocate(partt(idimp,nprobt))
            ftname = 'tr3.'//cdrun
            if ((nst==1).or.(nst==2)) then
               mtt = (nloop - 1)/ntt + 1
               itt = 0
               allocate(partd(mtt,idimp,nprobt))
               partd = 0.0
! open file for trajectory data: updates ntrec and possibly iut
               if (kstrt==1) then
                  if (ntrec==0) then
                     call dafopentr3(partt,iut,ntrec,trim(ftname))
                  endif
               endif
            else if (nst==3) then
               allocate(fvtp(2*nmv+2,ndim),fvmtp(ndim,3))
               allocate(fetp(ndim,0))
               fvtp(2*nmv+2,:) = 2.0*ws
! open file for test particle diagnostic: updates ntrec and possibly iut
               if (kstrt==1) then
                  if (ntrec==0) then
                     ws = 0.0
                     call dafopenfv3(fvmtp,fvtp,fetp,ws,iut,ntrec,      &
     &trim(ftname))
                  endif
               endif
            endif
         endif
      endif
!
! initialize phase space diagnostic
      if (nts > 0) then
         nmv21 = 2*nmv + 1
         mvx = min(mvx,nx); mvy = min(mvy,nypmn); mvz = min(mvz,nzpmn)
         nsxb = (nx - 1)/mvx + 1; nsyb = (ny - 1)/mvy + 1
         nszb = (nz - 1)/mvz + 1
         nyb = (noff(1) + nyzp(1) - 1)/mvy - (noff(1) - 1)/mvy
         nzb = (noff(2) + nyzp(2) - 1)/mvz - (noff(2) - 1)/mvz
         if (js==0) nyb = nyb + 1
         if (ks==0) nzb = nzb + 1
         itot(1) = nyb; itot(2) = nzb
         call mpimax(itot(:2),tdiag)
         nybmx = itot(1); nzbmx = itot(2)
! electron phase space diagnostic
         if ((nds==1).or.(nds==3)) then
! estimate maximum electron velocity or momentum
            ws = 0.0
            if (npxyz.gt.0.0d0) then
               ws = 4.0*vtx+abs(vx0)
               ws = max(ws,4.0*vty+abs(vy0))
               ws = max(ws,4.0*vtz+abs(vz0))
            endif
            if (npxyzb.gt.0.0d0) then
               ws = max(ws,4.0*vtdx+abs(vdx))
               ws = max(ws,4.0*vtdy+abs(vdy))
               ws = max(ws,4.0*vtdz+abs(vdz))
            endif
            allocate(fvs(nmv21+1,ndim,nsxb,nybmx+1,nzb+1))
            fvs(nmv21+1,:,1,1,1) = 1.25*ws
! open file for electron phase space data:
! updates nserec and possibly iuse
! opens a new fortran unformatted stream file
            if (nserec==0) then
               if (kstrt==1) then
                  fsename = 'pse3.'//cdrun
                  iuse =  get_funit(iuse)
                  call fnopens3(iuse,trim(fsename))
               endif
! writes distributed non-uniform partition information
               call mpwrncomp3(nyb,nzb,nvpy,nvpz,iuse)
               nserec = 1
            endif
         endif
! ion phase space diagnostic
         if (movion==1) then
            if ((nds==2).or.(nds==3)) then
! estimate maximum ion velocity or momentum
               ws = 0.0
               if (npxyzi.gt.0.0d0) then
                   ws = 4.0*vtxi+abs(vxi0)
                   ws = max(ws,4.0*vtyi+abs(vyi0))
                   ws = max(ws,4.0*vtzi+abs(vzi0))
               endif
               if (npxyzbi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
                  ws = max(ws,4.0*vtdzi+abs(vdzi))
               endif
               allocate(fvsi(nmv21+1,ndim,nsxb,nybmx+1,nzb+1))
               fvsi(nmv21+1,:,1,1,1) = 1.25*ws
! open file for ion phase space data:
! updates nsirec and possibly iusi
! opens a new fortran unformatted stream file
               if (nsirec==0) then
                  if (kstrt==1) then
                     fsiname = 'psi3.'//cdrun
                     iusi =  get_funit(iusi)
                     call fnopens3(iusi,trim(fsiname))
                  endif
! writes distributed non-uniform partition information
                  call mpwrncomp3(nyb,nzb,nvpy,nvpz,iusi)
                  nsirec = 1
               endif
            endif
         endif
      endif
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (dt > 0.37*ci) then
         if (kstrt==1) then
            write (*,*) 'Warning: Courant condition may be exceeded!'
         endif
      endif
!
      if (kstrt==1) write (iuot,*) 'program mpbbeps3'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
      if (kstrt==1) write (iuot,*) 'ntime = ', ntime
!
! save previous current in fourier space for radiative vector potential
      if (ntar > 0) then
         it = ntime/ntar
         if (ntime==ntar*it) oldcut = cut
      endif
!
! fluid moments diagnostic
      if (ntfm > 0) then
         it = ntime/ntfm
         if (ntime==ntfm*it) then
! calculate electron fluid moments
            if ((ndfm==1).or.(ndfm==3)) then
               call dtimer(dtime,itime,-1)
               fmse = 0.0
               call dtimer(dtime,itime,1)
               tdiag = tdiag + real(dtime)
               call wmprofx3(ppart,fmse,kpic,noff,ci,tdiag,npro,mx,my,mz&
     &,mx1,myp1,relativity)
! add guard cells with OpenMP: updates fmse
               call wmpnacguard3(fmse,nyzp,tdiag,nx,kstrt,nvpy,nvpz)
! moves vector grid fmse from non-uniform to uniform partition
               isign = -1
               call mpfnmove3(fmse,noff,nyzp,isign,tdiag,kyp,kzp,ny,nz, &
     &kstrt,nvpy,nvpz,mterf,irc)
               if (irc /= 0) then
                  call PPEXIT(); stop
               endif
! calculates fluid quantities from fluid moments
               call mpfluidqs3(fmse,tdiag,npro,nx,ny,nz,kstrt,nvpy,kyp, &
     &kzp)
! write real space diagnostic output: updates nferec
               call mpvwrite3(fmse,tdiag,nx,ny,nz,kyp,kzp,nvpy,iufe,    &
     &nferec)
!              call mpvread3(fmse,tdiag,nx,ny,nz,kyp,kzp,nvpy,iufe,     &
!    &nferec,irc)
            endif
! calculate ion fluid moments
            if (movion==1) then
               if ((ndfm==2).or.(ndfm==3)) then
                  call dtimer(dtime,itime,-1)
                  fmsi = 0.0
                  call dtimer(dtime,itime,1)
                  tdiag = tdiag + real(dtime)
                  call wmprofx3(pparti,fmsi,kipic,noff,ci,tdiag,npro,mx,&
     &my,mz,mx1,myp1,relativity)
                  fmsi = rmass*fmsi
! add guard cells with OpenMP: updates fmsi
                  call wmpnacguard3(fmsi,nyzp,tdiag,nx,kstrt,nvpy,nvpz)
! moves vector grid fmsi from non-uniform to uniform partition
                  isign = -1
                  call mpfnmove3(fmsi,noff,nyzp,isign,tdiag,kyp,kzp,ny, &
     &nz,kstrt,nvpy,nvpz,mterf,irc)
                  if (irc /= 0) then
                     call PPEXIT(); stop
                  endif
! calculates fluid quantities from fluid moments: updates fmsi
                  call mpfluidqs3(fmsi,tdiag,npro,nx,ny,nz,kstrt,nvpy,  &
     &kyp,kzp)
! write real space diagnostic output: updates nfirec
                  call mpvwrite3(fmsi,tdiag,nx,ny,nz,kyp,kzp,nvpy,iufi, &
     &nfirec)
!                 call mpvread3(fmsi,tdiag,nx,ny,nz,kyp,kzp,nvpy,iufi,     &
!    &nfirec,irc)
               endif
            endif
         endif
      endif
!
! deposit electron current with OpenMP:
! updates ppart and cue, and possibly ncl, ihole, irc
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      call wmpdjpost3(ppart,cue,kpic,ncl,ihole,noff,nyzp,qme,dth,ci,    &
     &tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
! add guard cells with OpenMP: updates cue
      call wmpnacguard3(cue,nyzp,tguard,nx,kstrt,nvpy,nvpz)
!
! reorder particles by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
      if (irc==0) then
         call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,xtras,tsort,tmov, &
     &kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,     &
     &mxzyp1,plist,irc2)
      else
         irc2(1) = 1; irc2(2) = irc; irc = 0
      endif
!
      do while (irc2(1) /= 0)
! ihole overflow
         if (irc2(1)==1) then
            ntmaxp = (1.0 + xtras)*irc2(2)
            deallocate(ihole)
            allocate(ihole(2,ntmaxp+1,mxyzp1))
            irc2 = 0
            call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,xtras,tsort,   &
     &tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,&
     &mxzyp1,.false.,irc2)
! ppart overflow
         else if (irc2(1)==4) then
! restores electron coordinates from ppbuff: updates ppart, ncl
            call mprstor3(ppart,ppbuff,ncl,ihole,tsort)
! copy ordered electrons to linear array: updates part
            call mpcopyout3(part,ppart,kpic,it,irc)
            deallocate(ppart)
            nppmx0 = (1.0 + xtras)*irc2(2)
            allocate(ppart(idimp,nppmx0,mxyzp1))
! copies unordered electrons to ordered array: updates ppart
            call mpcopyin3(part,ppart,kpic,irc)
            call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,xtras,tsort,   &
     &tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,&
     &mxzyp1,plist,irc2)
         endif
      enddo
!
! sanity check for electrons
      if (monitor > 0) then
         call mpcheck3(ppart,kpic,noff,nyzp,nx,mx,my,mz,mx1,myp1,irc)
      endif
!
! electron current density diagnostic
      if (ntje > 0) then
         it = ntime/ntje
         if (ntime==ntje*it) then
            vfield = cue
! transform electron current density to fourier space: updates cut
! moves data to uniform partition
            isign = -1
            call wmpfft3rn(vfield,cut,noff,nyzp,isign,mixup,sct,tfft,   &
     &tfmov,indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf&
     &,ierr)
! calculate smoothed electron current in fourier space: updates vfieldc
            call mpsmooth33(cut,vfieldc,ffc,tfield,nx,ny,nz,kstrt,nvpy, &
     &nvpz)
! store selected fourier modes: updates curet
            call mprdvmodes3(vfieldc,curet,tfield,nx,ny,nz,modesxje,    &
     &modesyje,modeszje,kstrt,nvpy,nvpz)
! transform smoothed electron current to real space: updates vfield
            isign = 1
            call mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates njerec
            call mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuje,     &
     &njerec)
!           call mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuje,njerec&
!    &,irc)
         endif
      endif
!
! deposit ion current with OpenMP:
      if (movion==1) then
! updates pparti and cui, and possibly ncl, ihole, irc
         call dtimer(dtime,itime,-1)
         cui = 0.0
         call dtimer(dtime,itime,1)
         tdjpost = tdjpost + real(dtime)
         call wmpdjpost3(pparti,cui,kipic,ncl,ihole,noff,nyzp,qmi,dth,ci&
     &,tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
! add guard cells with OpenMP: updates cui
         call wmpnacguard3(cui,nyzp,tguard,nx,kstrt,nvpy,nvpz)
!
! reorder particles by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
         if (irc==0) then
            call ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,xtras,tsort, &
     &tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,&
     &mxzyp1,plist,irc2)
         else
            irc2(1) = 1; irc2(2) = irc; irc = 0
         endif
!
         do while (irc2(1) /= 0)
! ihole overflow
            if (irc2(1)==1) then
               ntmaxp = (1.0 + xtras)*irc2(2)
               deallocate(ihole)
               allocate(ihole(2,ntmaxp+1,mxyzp1))
               irc2 = 0
               call ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,xtras,    &
     &tsort,tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1&
     &,mzp1,mxzyp1,.false.,irc2)
! pparti overflow
            else if (irc2(1)==4) then
! restores ion coordinates from ppbuff: updates pparti, ncl
               call mprstor3(pparti,ppbuff,ncl,ihole,tsort)
! copy ordered electrons to linear array: updates part
               call mpcopyout3(part,pparti,kipic,it,irc)
               deallocate(pparti)
               nppmx1 = (1.0 + xtras)*irc2(2)
               allocate(pparti(idimp,nppmx1,mxyzp1))
! copies unordered ions to ordered array: updates pparti
               call mpcopyin3(part,pparti,kipic,irc)
               call ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,xtras,    &
     &tsort,tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1&
     &,mzp1,mxzyp1,plist,irc2)
            endif
         enddo
!
! sanity check for ions
         if (monitor > 0) then
            call mpcheck3(pparti,kipic,noff,nyzp,nx,mx,my,mz,mx1,myp1,  &
     &irc)
         endif
      endif
!
! ion current density diagnostic
      if (movion==1) then
         if (ntji > 0) then
            it = ntime/ntji
            if (ntime==ntji*it) then
               vfield = cui
! transform ion current density to fourier space: updates cut
! moves data to uniform partition
               isign = -1
               call wmpfft3rn(vfield,cut,noff,nyzp,isign,mixup,sct,tfft,&
     &tfmov,indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf&
     &,ierr)
! calculate smoothed ion current in fourier space: updates vfieldc
               call mpsmooth33(cut,vfieldc,ffc,tfield,nx,ny,nz,kstrt,   &
     &nvpy,nvpz)
! store selected fourier modes: updates curit
               call mprdvmodes3(vfieldc,curit,tfield,nx,ny,nz,modesxji, &
     &modesyji,modeszji,kstrt,nvpy,nvpz)
! transform smoothed ion current to real space: updates vfield
               isign = 1
               call mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,  &
     &indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates njirec
               call mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuji,  &
     &njirec)
!              call mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuji,   &
!    &njirec,irc)
            endif
         endif
      endif
!
! deposit electron charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      call mppost3(ppart,qe,kpic,noff,qme,tdpost,mx,my,mz,mx1,myp1)
! add guard cells with OpenMP: updates qe
      call wmpaguard3(qe,nyzp,tguard,nx,kstrt,nvpy,nvpz)
!
! electron density diagnostic
      if (ntde > 0) then
         it = ntime/ntde
         if (ntime==ntde*it) then
            sfield = -qe
! transform electron density to fourier space: updates qt
! moves data to uniform partition
            isign = -1
            call wmpfft3r(sfield,qt,noff,nyzp,isign,mixup,sct,tfft,tfmov&
     &,indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr&
     &)
! calculate smoothed electron density in fourier space: updates sfieldc
            call mpsmooth3(qt,sfieldc,ffc,tfield,nx,ny,nz,kstrt,nvpy,   &
     &nvpz)
! store selected fourier modes: updates denet
            call mprdmodes3(sfieldc,denet,tfield,nx,ny,nz,modesxde,     &
     &modesyde,modeszde,kstrt,nvpy,nvpz)
! transform smoothed electron density to real space: updates sfield
            isign = 1
            call mpfft3r(sfield,sfieldc,isign,mixup,sct,tfft,indx,indy, &
     &indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates nderec
            call mpwrite3(sfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iude,nderec&
     &)
!           call mpread3(sfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iude,nderec,&
!    &irc)
         endif
      endif
!
! deposit ion charge with OpenMP: updates qi
      if (movion==1) then
         call dtimer(dtime,itime,-1)
         qi = 0.0
         call dtimer(dtime,itime,1)
         tdpost = tdpost + real(dtime)
         call mppost3(pparti,qi,kipic,noff,qmi,tdpost,mx,my,mz,mx1,myp1)
! add guard cells with OpenMP: updates qi
         call wmpaguard3(qi,nyzp,tguard,nx,kstrt,nvpy,nvpz)
      endif
!
! ion density diagnostic
      if (movion==1) then
         if (ntdi > 0) then
            it = ntime/ntdi
            if (ntime==ntdi*it) then
               sfield = qi
! transform ion density to fourier space: updates qt
! moves data to uniform partition
               isign = -1
               call wmpfft3r(sfield,qt,noff,nyzp,isign,mixup,sct,tfft,  &
     &tfmov,indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf&
     &,ierr)
! calculate smoothed ion density in fourier space: updates sfieldc
               call mpsmooth3(qt,sfieldc,ffc,tfield,nx,ny,nz,kstrt,nvpy,&
     &nvpz)
! store selected fourier modes: updates denit
               call mprdmodes3(sfieldc,denit,tfield,nx,ny,nz,modesxdi,  &
     &modesydi,modeszdi,kstrt,nvpy,nvpz)
! transform smoothed ion density to real space: updates sfield
               isign = 1
               call mpfft3r(sfield,sfieldc,isign,mixup,sct,tfft,indx,   &
     &indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates ndirec
               call mpwrite3(sfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iudi,   &
     &ndirec)
!              call mpread3(sfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iudi,    &
!    &ndirec,irc)
            endif
         endif
      endif
!
! add electron and ion densities: updates qe
      call mpaddqei3(qe,qi,nyzp,tfield,nx)
!
! add electron and ion current densities: updates cue
      if (movion==1) call mpaddcuei3(cue,cui,nyzp,tfield,nx)
!
! transform charge to fourier space with OpenMP:
! moves data to uniform partition
! updates qt, mterf, and ierr, modifies qe
      isign = -1
      call wmpfft3r(qe,qt,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx,indy&
     &,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! transform current to fourier space with OpenMP:
! moves data to uniform partition
! updates cut, mterf, and ierr, modifies cue
      isign = -1
      call wmpfft3rn(cue,cut,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! take transverse part of current with OpenMP: updates cut
      call mpcuperp3(cut,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! radiative vector potential diagnostic
      if (ntar > 0) then
         it = ntime/ntar
         if (ntime==ntar*it) then
! average current: updates vfieldc = 0.5*(cut + oldcut)
            call mcuave3(vfieldc,cut,oldcut,tfield,nz)
! calculate radiative vector potential in fourier space: updates vfieldc
! vfieldc should contain averaged current on entry
            call mpavrpot3(vfieldc,bxyz,ffc,affp,ci,tfield,nx,ny,nz,    &
     &kstrt,nvpy,nvpz)
! store selected fourier modes: updates vpotr
            call mprdvmodes3(vfieldc,vpotr,tfield,nx,ny,nz,modesxar,    &
     &modesyar,modeszar,kstrt,nvpy,nvpz)
! transform radiative vector potential to real space: updates vfield
            isign = 1
            call mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates narrec
            call mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuar,     &
     &narrec)
!           call mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuar,narrec&
!    &,irc)
         endif
      endif
!
! calculate electromagnetic fields in fourier space with OpenMP:
! updates exyz, bxyz, wf, wb
      if ((ntime+ntime0)==0) then
! initialize electromagnetic fields from darwin fields
! calculate initial darwin magnetic field
         call mpibpois3(cut,bxyz,ffc,ci,wb,tfield,nx,ny,nz,kstrt,nvpy,  &
     &nvpz)
         wf = 0.0
! calculate initial darwin electric field
         allocate(amu(6,nxe,nypmx,nzpmx),amut(6,nze,kxyp,kyzp))
         allocate(dcut(ndim,nze,kxyp,kyzp))
         amu = 0.0
         call wmpgmjpost3(ppart,amu,kpic,noff,qme,ci,tdjpost,mx,my,mz,  &
     &mx1,myp1,relativity)
         call wmpnacguard3(amu,nyzp,tguard,nx,kstrt,nvpy,nvpz)
         isign = -1
         call wmpfft3rn(amu,amut,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
         deallocate(amu)
         call mpdcuperp3(dcut,amut,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
         deallocate(amut)
         call mpetfield3(dcut,exyz,ffc,affp,ci,wf,tfield,nx,ny,nz,kstrt,&
     &nvpy,nvpz)
         deallocate(dcut)
         dth = 0.5*dt
! update electromagnetic fields
      else
         call mpmaxwel3(exyz,bxyz,cut,ffc,affp,ci,dt,wf,wb,tfield,nx,ny,&
     &nz,kstrt,nvpy,nvpz)
      endif
!
! calculate longitudinal force/charge in fourier space with OpenMP:
! updates fxyzt, we
      call mppois3(qt,fxyzt,ffc,we,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! add longitudinal and transverse electric fields with OpenMP:
! updates fxyzt
      isign = 1
      call mpemfield3(fxyzt,exyz,ffc,isign,tfield,nx,ny,nz,kstrt,nvpy,  &
     &nvpz)
! copy magnetic field with standard procedure: updates bxyzt
      isign = -1
      call mpemfield3(bxyzt,bxyz,ffc,isign,tfield,nx,ny,nz,kstrt,nvpy,  &
     &nvpz)
!
! transform electric force to real space with OpenMP:
! moves data to non-uniform partition
! updates fxyze, mterf, and ierr, modifies fxyzt
      isign = 1
      call wmpfft3rn(fxyze,fxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! transform magnetic force to real space with OpenMP:
! moves data to non-uniform partition
! updates bxyze, modifies bxyzt
      isign = 1
      call wmpfft3rn(bxyze,bxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! add constant to magnetic field with OpenMP: updates bxyze
      if (omt > 0.0) call mpbaddext3(bxyze,nyzp,tfield,omx,omy,omz,nx)
!
! copy guard cells with OpenMP: updates fxyze, bxyze
      call wmpncguard3(fxyze,nyzp,tguard,nx,kstrt,nvpy,nvpz)
      call wmpncguard3(bxyze,nyzp,tguard,nx,kstrt,nvpy,nvpz)
!
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! calculate potential in fourier space: updates sfieldc
            call mppot3(qt,sfieldc,ffc,ws,tfield,nx,ny,nz,kstrt,nvpy,   &
     &nvpz)
! store selected fourier modes: updates pott
            call mprdmodes3(sfieldc,pott,tfield,nx,ny,nz,modesxp,modesyp&
     &,modeszp,kstrt,nvpy,nvpz)
! transform potential to real space: updates sfield
            isign = 1
            call mpfft3r(sfield,sfieldc,isign,mixup,sct,tfft,indx,indy, &
     &indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates nprec
            call mpwrite3(sfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iup,nprec)
!           call mpread3(sfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iup,nprec,  &
!    &irc)
         endif
      endif
!
! longitudinal efield diagnostic
      if (ntel > 0) then
         it = ntime/ntel
         if (ntime==ntel*it) then
! calculate longitudinal efield in fourier space: updates vfieldc
            call mpelfield3(qt,vfieldc,ffc,ws,tfield,nx,ny,nz,kstrt,nvpy&
     &,nvpz)
! store selected fourier modes: updates elt
            call mprdvmodes3(vfieldc,elt,tfield,nx,ny,nz,modesxel,      &
     &modesyel,modeszel,kstrt,nvpy,nvpz)
! transform longitudinal efield to real space: updates vfield
            isign = 1
            call mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates nelrec
            call mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuel,     &
     &nelrec)
!           call mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuel,nelrec&
!    &,irc)
         endif
      endif
!
! vector potential diagnostic
      if (nta > 0) then
         it = ntime/nta
         if (ntime==nta*it) then
! calculate vector potential in fourier space: updates vfieldc
            call mpavpot3(bxyz,vfieldc,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! store selected fourier modes: updates vpott
            call mprdvmodes3(vfieldc,vpott,tfield,nx,ny,nz,modesxa,     &
     &modesya,modesza,kstrt,nvpy,nvpz)
! transform vector potential to real space: updates vfield
            isign = 1
            call mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates narec
            call mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iua,narec)
!           call mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iua,narec, &
!    &irc)
         endif
      endif
!
! transverse efield diagnostic
      if (ntet > 0) then
         it = ntime/ntet
         if (ntime==ntet*it) then
            vfieldc = exyz
! store selected fourier modes: updates ett
            call mprdvmodes3(vfieldc,ett,tfield,nx,ny,nz,modesxet,      &
     &modesyet,modeszet,kstrt,nvpy,nvpz)
! transform transverse efield to real space: updates vfield
            isign = 1
            call mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates netrec
            call mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuet,     &
     &netrec)
!           call mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iuet,netrec&
!    &,irc)
         endif
      endif
!
! magnetic field diagnostic
      if (ntb > 0) then
         it = ntime/ntb
         if (ntime==ntb*it) then
            vfieldc = bxyz
! store selected fourier modes: updates bt
            call mprdvmodes3(vfieldc,bt,tfield,nx,ny,nz,modesxb,modesyb,&
     &modeszb,kstrt,nvpy,nvpz)
! transform magnetic field to real space: updates vfield
            isign = 1
            call mpfft3rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! write real space diagnostic output: updates nbrec
            call mpvwrite3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iub,nbrec)
!           call mpvread3(vfield,tdiag,nx,ny,nz,kyp,kzp,nvpy,iub,nbrec, &
!    &irc)
         endif
      endif
!
! velocity-space diagnostic
      if (ntv > 0) then
         it = ntime/ntv
         if (ntime==ntv*it) then
            if ((ndv==1).or.(ndv==3)) then
! calculate electron cartesian distribution function and moments
               if ((nvft==1).or.(nvft==3)) then
                  call mpvpdist3(ppart,kpic,fv,sfv,fvm,tdiag,nvp,nmv)
! store time history electron vdrift, vth, and entropy
                  itv = itv + 1
                  fvtm(itv,:,:) = fvm
               endif
! calculate electron cylindrical distribution function and moments
               if ((nvft==4).or.(nvft==5)) then
                  call mpvbpdist3(ppart,kpic,fv,sfv,fvm,omx,omy,omz,    &
     &tdiag,nmv)
               endif
! electron energy distribution
               if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
                  call mperpdist3(ppart,kpic,fe,sfv,eci,wk,tdiag,nmv)
               endif
! write electron velocity-space diagnostic output: updates nverec
               if (kstrt==1) then
                  call dafwritefv3(fvm,fv,fe,wk,tdiag,iuve,nverec)
               endif
            endif
! ion distribution functions
            if (movion==1) then
               if ((ndv==2).or.(ndv==3)) then
! calculate ion cartesian distribution function and moments
                  if ((nvft==1).or.(nvft==3)) then
                     call mpvpdist3(pparti,kipic,fvi,sfv,fvmi,tdiag,nvp,&
     &nmv)
! store time history of ion vdrift, vth, and entropy
                     fvtmi(itv,:,:) = fvmi
                  endif
! calculate ion cylindrical distribution function and moments
                  if ((nvft==4).or.(nvft==5)) then
                     call mpvbpdist3(pparti,kipic,fvi,sfv,fvmi,omx,omy, &
     &omz,tdiag,nmv)
                  endif
! ion energy distribution
                  if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
                     call mperpdist3(pparti,kipic,fei,sfv,eci,wk,tdiag, &
     &nmv)
                     wk = rmass*wk
                  endif
! write ion velocity-space diagnostic output: updates nvirec
                  if (kstrt==1) then
                     call dafwritefv3(fvmi,fvi,fei,wk,tdiag,iuvi,nvirec)
                  endif
               endif
            endif
         endif
      endif
!
! trajectory diagnostic
      if (ntt > 0) then
         it = ntime/ntt
         if (ntime==ntt*it) then
            ierr = 0
! copies tagged electrons in ppart to array partt: updates partt, numtp
            if (ndt==1) then
               call mptraj3(ppart,kpic,partt,tdiag,numtp,ierr)
! copies tagged ions in ppart to array partt: updates partt, numtp
            else if (ndt==2) then
               if (movion==1) then
                  call mptraj3(pparti,kipic,partt,tdiag,numtp,ierr)
               endif
            endif
            if (ierr /= 0) then
! partt overflow
               if (ierr <  0) then
                  deallocate(partt)
                  nprobt = numtp
                  allocate(partt(idimp,nprobt))
                  ierr = 0
! copies tagged electrons in ppart to array partt: updates partt, numtp
                  if (ndt==1) then
                     call mptraj3(ppart,kpic,partt,tdiag,numtp,ierr)
! copies tagged ions in ppart to array partt: updates partt, numtp
                  else if (ndt==2) then
                     if (movion==1) then
                        call mptraj3(pparti,kipic,partt,tdiag,numtp,ierr&
    &)
                     endif
                  endif
               endif
               if (ierr /= 0) then
                  call PPEXIT(); stop
               endif
            endif
! electron or ion trajectories
            if ((ndt==1).or.(ndt==2)) then
! reorder tagged particles
               if ((nst==1).or.(nst==2)) then
! determines list of tagged particles leaving this node: updates iholep
                  call mpfholes3(partt,tedges,numtp,iholep,2)
! iholep overflow
                  if (iholep(1,1) < 0) then
                     ntmax = -iholep(1,1)
                     ntmax = 1.5*ntmax
                     deallocate(iholep)
                     if (kstrt==1) then
                        write (*,*) 'info:reallocating iholep:ntmax=',  &
     &ntmax
                     endif
                     allocate(iholep(ntmax+1,2))
                     call mpfholes3(partt,tedges,numtp,iholep,2)
                     if (iholep(1,1) < 0) then
                        if (kstrt==1) then
                           write (*,*) 'iholep overflow: ntmax=', ntmax
                           call PPEXIT(); stop
                        endif
                     endif
                  endif
! copies tagged particles: updates part
                  call mpcpytraj3(partt,part,tdiag,numtp)
! moves tagged electrons into original spatial region:
! updates part, numtp
                  call ipmove3(part,tedges,numtp,iholep,ny,nz,tmov,kstrt&
     &,nvpy,nvpz,2,ierr)
                  if (ierr /= 0) then
                     call PPEXIT(); stop
                  endif
! reorders tagged particles: updates partt
                  call mpordtraj3(part,partt,tedges,tdiag,numtp,ierr)
! collects distributed test particle data onto node 0
                  call mppartt3(partt,tdiag,numtp,nvpy,nvpz,ierr)
                  if (ierr /= 0) then
                     call PPEXIT(); stop
                  endif
! write trajectory diagnostic output: updates ntrec
                  if (kstrt==1) then
                     call dafwritetr3(partt,tdiag,iut,ntrec)
                     itt = itt + 1
                     partd(itt,:,:) = partt
                  endif
               else if (nst==3) then
! calculate test particle distribution function and moments
                  call mpvdist3(partt,fvtp,fvmtp,tdiag,numtp,nvp,nmv)
! write test particle diagnostic output: updates ntrec
                  if (kstrt==1) then
                     ws = 0.0
                     call dafwritefv3(fvmtp,fvtp,fetp,ws,tdiag,iut,ntrec)
                  endif
               endif
            endif
         endif
      endif
!
! phase space diagnostic
      if (nts > 0) then
         it = ntime/nts
         if (ntime==nts*it) then
! electron phase space diagnostic
            if ((nds==1).or.(nds==3)) then
! calculates velocity distribution in different regions of space:
! updates fvs
               call mpvspdist3(ppart,kpic,fvs,noff,tdiag,nmv,mvx,mvy,mvz&
     &,nyb)
! adjusts 3d velocity distribution in different regions of space:
! updates fvs
               call mpadjfvs3(fvs,tdiag,noff,nyzp,nmv,mvy,mvz,nyb,nvpy, &
     &nzbmx)
! write phase space diagnostic output: updates nserec
               if (nserec > 0) then
                  call mpwrfvsdata3(fvs,tdiag,nyb,nzb,nzbmx,iuse)
                  nserec = nserec + 1
               endif
            endif
! ion phase space
            if (movion==1) then
               if ((nds==2).or.(nds==3)) then
! calculates velocity distribution in different regions of space:
! updates fvsi
                  call mpvspdist3(pparti,kipic,fvsi,noff,tdiag,nmv,mvx, &
     &mvy,mvz,nyb)
! adjusts 3d velocity distribution in different regions of space:
! updates fvsi
                  call mpadjfvs3(fvsi,tdiag,noff,nyzp,nmv,mvy,mvz,nyb,  &
     &nvpy,nzbmx)
! write phase space diagnostic output: updates nsirec
                  if (nsirec > 0) then
                     call mpwrfvsdata3(fvsi,tdiag,nyb,nzb,nzbmx,iusi)
                     nsirec = nsirec + 1
                  endif
               endif
            endif
         endif
      endif
!
! push electrons with OpenMP:
! updates ppart and wke, and possibly ncl, ihole, irc
      wke = 0.0
      call wmpbpush3(ppart,fxyze,bxyze,kpic,ncl,ihole,noff,nyzp,qbme,dt,&
     &dth,ci,wke,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,&
     &irc)
!
! reorder electrons by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
      if (irc==0) then
         call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,xtras,tsort,tmov, &
     &kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,     &
     &mxzyp1,plist,irc2)
      else
         irc2(1) = 1; irc2(2) = irc; irc = 0
      endif
!
      do while (irc2(1) /= 0)
! ihole overflow
         if (irc2(1)==1) then
            ntmaxp = (1.0 + xtras)*irc2(2)
            deallocate(ihole)
            allocate(ihole(2,ntmaxp+1,mxyzp1))
            irc2 = 0
            call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,xtras,tsort,   &
     &tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,&
     &mxzyp1,.false.,irc2)
! ppart overflow
         else if (irc2(1)==4) then
! restores electron coordinates from ppbuff: updates ppart, ncl
            call mprstor3(ppart,ppbuff,ncl,ihole,tsort)
! copy ordered electrons to linear array: updates part
            call mpcopyout3(part,ppart,kpic,it,irc)
            deallocate(ppart)
            nppmx0 = (1.0 + xtras)*irc2(2)
            allocate(ppart(idimp,nppmx0,mxyzp1))
! copies unordered electrons to ordered array: updates ppart
            call mpcopyin3(part,ppart,kpic,irc)
            call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,xtras,tsort,   &
     &tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,&
     &mxzyp1,plist,irc2)
         endif
      enddo
!
! sanity check for electrons
      if (monitor > 0) then
         call mpcheck3(ppart,kpic,noff,nyzp,nx,mx,my,mz,mx1,myp1,irc)
      endif
!
! push ions with OpenMP:
      if (movion==1) then
! updates pparti and wki, and possibly ncl, ihole, irc
         wki = 0.0
         call wmpbpush3(pparti,fxyze,bxyze,kipic,ncl,ihole,noff,nyzp,   &
     &qbmi,dt,dth,ci,wki,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,         &
     &relativity,plist,irc)
         wki = wki*rmass
!
! reorder ions by tile with OpenMP and MPI
! updates: pparti, kipic, and irc and possibly ncl and ihole
         if (irc==0) then
            call ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,xtras,tsort, &
     &tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,&
     &mxzyp1,plist,irc2)
         else
            irc2(1) = 1; irc2(2) = irc; irc = 0
         endif
!
         do while (irc2(1) /= 0)
! ihole overflow
            if (irc2(1)==1) then
               ntmaxp = (1.0 + xtras)*irc2(2)
               deallocate(ihole)
               allocate(ihole(2,ntmaxp+1,mxyzp1))
               irc2 = 0
               call ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,xtras,    &
     &tsort,tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1&
     &,mzp1,mxzyp1,.false.,irc2)
! pparti overflow
            else if (irc2(1)==4) then
! restores ion coordinates from ppbuff: updates pparti, ncl
               call mprstor3(pparti,ppbuff,ncl,ihole,tsort)
! copy ordered electrons to linear array: updates part
               call mpcopyout3(part,pparti,kipic,it,irc)
               deallocate(pparti)
               nppmx1 = (1.0 + xtras)*irc2(2)
               allocate(pparti(idimp,nppmx1,mxyzp1))
! copies unordered ions to ordered array: updates pparti
               call mpcopyin3(part,pparti,kipic,irc)
               call ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,xtras,    &
     &tsort,tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1&
     &,mzp1,mxzyp1,plist,irc2)
            endif
         enddo
!
! sanity check for ions
         if (monitor > 0) then
            call mpcheck3(pparti,kipic,noff,nyzp,nx,mx,my,mz,mx1,myp1,  &
     &irc)
         endif
      endif
!
! start running simulation backwards:
! need to advance maxwell field solver one step ahead
      if (treverse==1) then
         if (((ntime+1)==(nloop/2)).or.((ntime+1)==nloop)) then
! deposit electron current: updates cue
            cue = 0.0
            call wmpdjpost3(ppart,cue,kpic,ncl,ihole,noff,nyzp,qme,zero,&
     &ci,tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
            call wmpnacguard3(cue,nyzp,tguard,nx,kstrt,nvpy,nvpz)
! deposit ion current: updates cui
            if (movion==1) then
               cui = 0.0
               call wmpdjpost3(pparti,cui,kipic,ncl,ihole,noff,nyzp,qmi,&
     &zero,ci,tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist, &
     &irc)
               call wmpnacguard3(cui,nyzp,tguard,nx,kstrt,nvpy,nvpz)
               call mpaddcuei3(cue,cui,nyzp,tfield,nx)
            endif
            isign = -1
            call wmpfft3rn(cue,cut,noff,nyzp,isign,mixup,sct,tfft,tfmov,&
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
            call mpcuperp3(cut,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! updates exyz, bxyz, wf, wb
            call mpmaxwel3(exyz,bxyz,cut,ffc,affp,ci,dt,wf,wb,tfield,nx,&
     &ny,nz,kstrt,nvpy,nvpz)
! reverse time
            dt = -dt; dth = -dth
         endif
      endif
!
! energy diagnostic
      if (ntw > 0) then
         it = ntime/ntw
         if (ntime==ntw*it) then
            wef = we + wf + wb
            wtot(1) = wef
            wtot(2) = wke
            wtot(3) = wki
            wtot(4) = wef + wke
            wtot(5) = we
            wtot(6) = wf
            wtot(7) = wb
            call mpdsum(wtot,tdiag)
            wke = wtot(2)
            wki = wtot(3)
            we = wtot(5)
            wf = wtot(6)
            wb = wtot(7)
            wef = we + wf + wb
            ws = wef + wke + wki
            if (ntime==0) s(6) = ws
            if (kstrt==1) then
               write (iuot,*) 'Total Field, Kinetic and Total Energies:'
               if (movion==0) then
                  write (iuot,'(3e14.7)') wef, wke, ws
               else
                  write (iuot,'(4e14.7)') wef, wke, wki, ws
               endif
               write (iuot,*) 'Electrostatic, Transverse Electric and Ma&
     &gnetic Field Energies:'
               write (iuot,'(3e14.7)') we, wf, wb
            endif
            itw = itw + 1
! store energies in time history array
            wt(itw,:) = (/wef,wke,wki,ws,we,wf,wb/)
            s(1) = s(1) + we
            s(2) = s(2) + wke
            s(3) = s(3) + wf
            s(4) = s(4) + wb
            s(5) = s(5) + wki
            s(6) = min(s(6),dble(ws))
            s(7) = max(s(7),dble(ws))
         endif
      endif
!
! restart file
      if (ntr > 0) then
         it = n/ntr
         if (n==ntr*it) then
            call dtimer(dtime,itime,-1)
! write out basic restart file for electrostatic code
            call bwrite_restart3(part,ppart,pparti,qi,kpic,kipic,tdiag, &
     &kstrt,iur,iscr,n,ntime0,irc)
! write out basic restart file for electromagnetic code
            call bwrite_brestart3(exyz,bxyz,tdiag,kstrt,iur)
            call dtimer(dtime,itime,1)
            tfield = tfield + real(dtime)
         endif
      endif
!
      enddo
      ntime = ntime + 1
!
! loop time
      call dtimer(dtime,ltime,1)
      tloop = tloop + real(dtime)
!
! * * * end main iteration loop * * *
!
      if (kstrt==1) then
         write (iuot,*)
         write (iuot,*) 'ntime, relativity = ', ntime, relativity
         if (treverse==1) write (iuot,*) 'treverse = ', treverse
         write (iuot,*) 'MPI nodes nvpy, nvpz = ', nvpy, nvpz
!
! energy diagnostic
         if (ntw > 0) then
            s(6) = (s(7) - s(6))/wt(1,4)
            write (iuot,*) 'Energy Conservation = ', real(s(6))
            s(1) = s(1)/real(itw)
            write (iuot,*) 'Average Field Energy <WE> = ', real(s(1))
            s(2) = s(2)/real(itw)
            write (iuot,*) 'Average Electron Kinetic Energy <WKE> = ',  &
     &real(s(2))
            write (iuot,*) 'Ratio <WE>/<WKE>= ', real(s(1)/s(2))
            s(3) = s(3)/real(itw)
            write (iuot,*) 'Average Transverse EField Energy <WF> = ',  &
     &real(s(3))
            write (iuot,*) 'Ratio <WF>/<WKE>= ', real(s(3)/s(2))
            s(4) = s(4)/real(itw)
            write (iuot,*) 'Average Magnetic Field Energy <WB> = ',     &
     &real(s(4))
            write (iuot,*) 'Ratio <WB>/<WKE>= ', real(s(4)/s(2))
         endif
!
         write (iuot,*)
         write (iuot,*) 'initialization time = ', tinit
         write (iuot,*) 'deposit time = ', tdpost
         write (iuot,*) 'current deposit time = ', tdjpost
         tdpost = tdpost + tdjpost
         write (iuot,*) 'total deposit time = ', tdpost
         write (iuot,*) 'guard time = ', tguard
         write (iuot,*) 'solver time = ', tfield
         write (iuot,*) 'field move time = ', tfmov
         write (iuot,*) 'fft and transpose time = ', tfft(1), tfft(2)
         write (iuot,*) 'push time = ', tpush
         write (iuot,*) 'particle move time = ', tmov
         write (iuot,*) 'sort time = ', tsort
         tfield = tfield + tguard + tfft(1) + tfmov
         write (iuot,*) 'total solver time = ', tfield
         tsort = tsort + tmov
         time = tdpost + tpush + tsort
         write (iuot,*) 'total particle time = ', time
         write (iuot,*) 'total diagnostic time = ', tdiag
         ws = time + tfield + tdiag
         tloop = tloop - ws
         write (iuot,*) 'total and additional time = ', ws, tloop
         write (iuot,*)
!
         ws = 1.0e+09/(real(nloop)*real(np+npi))
         write (iuot,*) 'Push Time (nsec) = ', tpush*ws
         write (iuot,*) 'Deposit Time (nsec) = ', tdpost*ws
         write (iuot,*) 'Sort Time (nsec) = ', tsort*ws
         write (iuot,*) 'Total Particle Time (nsec) = ', time*ws
!
! reset parameters for final diagnostic metafile
! electron density diagnostic
         if (ntde > 0) then
            nderec = nderec - 1
         endif
! potential diagnostic
         if (ntp > 0) then
            nprec = nprec - 1; ceng = affp
         endif
! longitudinal efield diagnostic
         if (ntel > 0) then
            nelrec = nelrec - 1; ceng = affp
         endif
! electron current diagnostic
         if (ntje > 0) then
            njerec = njerec - 1
         endif
! radiative vector potential diagnostic
         if (ntar > 0) then
            narrec = narrec - 1; ceng = affp
         endif
! vector potential diagnostic
         if (nta > 0) then
            narec = narec - 1; ceng = affp
         endif
! transverse efield diagnostic
         if (ntet > 0) then
            netrec = netrec - 1; ceng = affp
         endif
! magnetic field diagnostic
         if (ntb > 0) then
            nbrec = nbrec - 1; ceng = affp
         endif
! fluid moments diagnostic
         if (ntfm > 0) then
            if ((ndfm==1).or.(ndfm==3)) nferec = nferec - 1
            if (movion==1) then
               if ((ndfm==2).or.(ndfm==3)) nfirec = nfirec - 1
            endif
            ceng = affp
         endif
! velocity-space diagnostic
         if (ntv > 0) then
            if ((ndv==1).or.(ndv==3)) nverec = nverec - 1
            if (movion==1) then
               if ((ndv==2).or.(ndv==3)) nvirec = nvirec - 1
            endif
         endif
! trajectory diagnostic
         if (ntt > 0) then
            ntrec = ntrec - 1
         endif
! phase space diagnostic
         if (nts > 0) then
            if ((nds==1).or.(nds==3)) nserec = nserec - 1
            if (movion==1) then
               if ((nds==2).or.(nds==3)) nsirec = nsirec - 1
            endif
         endif
! ion diagnostics
         if (movion==1) then
! ion density diagnostic
            if (ntdi > 0) then
               ndirec = ndirec - 1
            endif
! ion current diagnostic
            if (ntji > 0) then
               njirec = njirec - 1
            endif
         endif
! write final diagnostic metafile
         call writnml3(iudm)
         close(unit=iudm)
! close restart files
         call close_restart3(iur,iur0)
! close output file
         write (iuot,*) ' * * * q.e.d. * * *'
         close(unit=iuot)
      endif
!
      call PPEXIT()
      end program
