!-----------------------------------------------------------------------
! 2-1/2D Electromagnetic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpbbeps2
      use in2
      use minit2
      use mpush2
      use mbpush2
      use mcurd2
      use mfield2
      use mdiag2
      use mppmod2
      use pgraf2
      use pgraf1
      use omplib
      use ompplib2
      use f2
      use fb2
      implicit none
!
! idimp = dimension of phase space = 5
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 5, ipbc = 1
! idps = number of partition boundaries
      integer :: idps = 2
! wke/wki/we = particle kinetic/electrostatic field energy
! wf/wb = transverse electric field/magnetic field
      real :: wke = 0.0, wki = 0.0, we = 0.0, wf = 0.0, wb = 0.0
      real :: zero = 0.0
! plist = (true,false) = list of particles leaving tiles found in push
      logical :: plist = .true.
!
! declare scalars for standard code
      integer :: n
      integer :: nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy
      integer :: mx1, ntime, nloop, isign, ierr
      integer :: ntime0 = 0
      real :: qbme, affp, dth, omt, ws
      real :: qbmi, vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      real :: xmin, xmax, ymin, ymax
      double precision :: npxy, npxyb, np, npxyi, npxybi
      double precision :: npi = 0.0d0
!
! declare scalars for MPI code
      integer :: idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn
      integer :: nyp, noff, npp, nps, myp1, mxyp1, ntmax
      integer :: npimax, nppi, maxnp
      integer :: nterf
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nppmx1, nbmaxp, ntmaxp, npbmx
      integer :: irc = 0
      integer, dimension(2) :: irc2 = 0
!
! declare scalars for diagnostics
      integer :: kxps, kyps, nmv21, numtp, nyb, nybmx, it, mtw, mtv, mtt
      integer :: itw, itv, itt
      real :: wef, ts, eci, wk
! default Fortran unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19
      integer :: iur = 17, iur0 = 27, iscr = 99
      integer :: iude = 10, iup = 11, iuel = 12
      integer :: iua = 13, iuet = 14, iub = 15, iuar = 16
      integer :: iuje = 21, iufe = 23, iuve = 25, iut = 28, iuse = 29
      integer :: iudi = 20, iuji = 22, iufi = 24, iuvi = 26, iusi = 30
! dimensions for fourier data
      integer :: modesxpd, modesy2de, modesy2p, modesy2el, modesy2di
      integer :: modesy2ar, modesy2a, modesy2et, modesy2b, modesy2je
      integer :: modesy2ji
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! declare arrays for standard code
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe/qi = electron/ion charge density with guard cells
      real, dimension(:,:), allocatable :: qe, qi
! cue/cui = electron/ion current density with guard cells
      real, dimension(:,:,:), allocatable :: cue, cui
! fxyze/bxyze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:,:), allocatable :: fxyze, bxyze
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:), allocatable :: exyz, bxyz
! qt = scalar charge density field array in fourier space
      complex, dimension(:,:), allocatable :: qt
! cut = vector current density field array in fourier space
! fxyt/bxyt = vector electric/magnetic field in fourier space
      complex, dimension(:,:,:), allocatable :: cut, fxyt, bxyt
! ffc = form factor array for poisson solver
      complex, dimension(:,:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
      integer, dimension(1) :: itot
      double precision, dimension(7) :: wtot
! arrays required for darwin initial fields
      real, dimension(:,:,:), allocatable :: amu
      complex, dimension(:,:,:), allocatable :: dcut, amut
!
! declare arrays for MPI code
! edges(1:2) = lower:upper y boundaries of particle partition
      real, dimension(:), allocatable  :: edges
! iholep = location of hole left in linear particle arrays
      integer, dimension(:), allocatable :: iholep
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
      complex, dimension(:,:), allocatable :: sfieldc
      real, dimension(:,:), allocatable :: sfield
! scratch arrays for vector field
      complex, dimension(:,:,:), allocatable :: vfieldc
      real, dimension(:,:,:), allocatable :: vfield
! denet/denit = store selected fourier modes for electron/ion density
      complex, dimension(:,:), allocatable :: denet, denit
! pott = store selected fourier modes for potential
      complex, dimension(:,:), allocatable :: pott
! elt = store selected fourier modes for longitudinal efield
      complex, dimension(:,:,:), allocatable :: elt
! curet = store selected fourier modes for electron current density
      complex, dimension(:,:,:), allocatable :: curet
! curit = store selected fourier modes for ion current density
      complex, dimension(:,:,:), allocatable :: curit
! oldcut = previous current density
      complex, dimension(:,:,:), allocatable :: oldcut
! vpotr = store selected fourier modes for radiative vector potential
      complex, dimension(:,:,:), allocatable :: vpotr
! vpott = store selected fourier modes for vector potential
      complex, dimension(:,:,:), allocatable :: vpott
! ett = store selected fourier modes for transverse efield
      complex, dimension(:,:,:), allocatable :: ett
! bt = store selected fourier modes for magnetic field
      complex, dimension(:,:,:), allocatable :: bt
! fmse/fmsi = electron/ion fluid moments
      real, dimension(:,:,:), allocatable :: fmse, fmsi
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
! tedges(1:2) = lower:upper boundary of particle tags
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
      real, dimension(:,:,:,:), allocatable :: fvs, fvsi
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tpush = 0.0, tsort = 0.0, tmov = 0.0
      real :: tfmov = 0.0, tdiag = 0.0
      real, dimension(2) :: tfft = 0.0
      double precision :: dtime
!
! start timing initialization
      call dtimer(dtime,itime,-1)
!      
! nvp = number of distributed memory nodes
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
! read namelist
         call readnml2(iuin)
! override input data
         idcode = 2
         ndim = 3
! create string from idrun
         write (cdrun,'(i10)') idrun
         cdrun = adjustl(cdrun)
! text output file
         fname = 'output2.'//cdrun
         open(unit=iuot,file=trim(fname),form='formatted',              &
     &status='replace')
      endif
!
! broadcast namelists to other nodes
      call sendnmls2()
!
! open graphics device
      call IPLTCOMM(nplot)
      if (kstrt==1) then
         irc = open_pgraphs(nplot)
! set palette to color wheel
         call STPALIT(2)
      endif
!
! increase number of coordinates for particle tag
      if (ntt > 0) idimp = idimp + 1
!
! initialize scalars for standard code
! np = total number of particles in simulation
      npxy = dble(npx)*dble(npy); npxyb = dble(npxb)*dble(npyb)
      np = npxy + npxyb
! npi = total number of ions in simulation
      if (movion > 0) then
         npxyi = dble(npxi)*dble(npyi); npxybi = dble(npxbi)*dble(npybi)
         npi = npxyi + npxybi
      endif
! nx/ny = number of grid points in x/y direction
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 2; nxeh = nxe/2; nnxe = ndim*nxe
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = dble(nx)*dble(ny)/np
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
! check if too many processors
      if (nvp > ny) then
         if (kstrt==1) then
         write (*,*) 'Too many processors requested: ny, nvp=', ny, nvp
         endif
         call PPEXIT(); stop
      endif
!
! initialize data for MPI code
      allocate(edges(idps))
! calculate partition variables: edges, nyp, noff, nypmx
! edges(1:2) = lower:upper boundary of particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! noff = lowermost global gridpoint in particle partition
! nypmx = maximum size of particle partition, including guard cells
! nypmn = minimum value of nyp
! find new 1d partition for uniform density distribution
!     call mpdcomp2(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp)
! find new 1d partition from initial analytic distribution function
!     call mpfedges2(edges,nyp,noff,ampdy,scaledy,shiftdy,nypmx,nypmn,ny&
!    &,kstrt,nvp,ipbc,ndprof,ierr)
! find new 1d partition from initial analytic distribution function
! where initial plasma density is confined to some subregion
      xmin = dxmin*real(nx); xmax = dxmax*real(nx)
      ymin = dymin*real(ny); ymax = dymax*real(ny)
      call mpgfedges2(edges,nyp,noff,ampdy,scaledy,shiftdy,ymin,ymax,   &
     &nypmx,nypmn,ny,kstrt,nvp,ndprof,ierr)
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
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! npmax/npimax = maximum number of electrons/ions in each partition
      npmax = (np/nvp)*1.25; npimax = (npi/nvp)*1.25
      maxnp = max(npmax,npimax)
! myp1 = number of tiles in y direction
      myp1 = (nyp - 1)/my + 1; mxyp1 = mx1*myp1
! nterf = number of shifts required by field manager (0=search)
      nterf = 0
! ntmax = size of iholep buffer for particles leaving node
      ntmax = 0.2*npmax
!
! allocate and initialize data for standard code
      allocate(part(idimp,maxnp))
      allocate(qe(nxe,nypmx),qi(nxe,nypmx))
      allocate(cue(ndim,nxe,nypmx))
      allocate(fxyze(ndim,nxe,nypmx),bxyze(ndim,nxe,nypmx))
      allocate(exyz(ndim,nye,kxp),bxyz(ndim,nye,kxp))
      allocate(qt(nye,kxp),fxyt(ndim,nye,kxp))
      allocate(cut(ndim,nye,kxp),bxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),mixup(nxhy),sct(nxyh))
      allocate(kpic(mxyp1),ncl(8,mxyp1),iholep(ntmax+1))
      if (movion==1) allocate(kipic(mxyp1),cui(ndim,nxe,nypmx))
!
! prepare fft tables
      call mpfft2_init(mixup,sct,indx,indy)
! calculate form factors
      call mppois2_init(ffc,ax,ay,affp,nx,ny,kstrt)
! initialize different ensemble of random numbers
      if (nextrand > 0) call mnextran2(nextrand,ndim,npmax+npimax)
!
! open restart files
      if (kstrt==1) call open_restart2(iur,iur0,cdrun)
!
! new start
      if (nustrt==1) then
! initialize electrons
         nps = 1
         npp = 0
! background electrons
         if (npxy > 0.0d0) then
! calculates initial electron co-ordinates with various density profiles
!           call mpfdistr2(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,&
!    &shiftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,ierr)
! where initial plasma density is confined to some subregion
            call mpgfdistr2(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy&
     &,shiftdy,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,ndprof,ierr)
! initialize electron velocities or momenta
            if (ierr==0) then
! special cases
               if (nvdist==3) then
                  call mpvbdistr2h(part,nps,npp,vtx,vtz,vx0,vz0,omx,omy,&
     &omz,npx,npy,kstrt,nvp,ierr)
! general cases
               else
                  call wmpvdistr2h(part,nps,npp,vtx,vty,vtz,vx0,vy0,vz0,&
     &ci,npx,npy,kstrt,nvp,relativity,ierr)
               endif
            endif
! check for background electron initialization error
            if (ierr /= 0) then
               call PPEXIT(); stop
            endif
         endif
! beam electrons
         if (npxyb > 0.0d0) then
            nps = npp + 1
! calculates initial electron co-ordinates with various density profiles
!           call mpfdistr2(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,&
!    &shiftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,ierr)
! where initial plasma density is confined to some subregion
            call mpgfdistr2(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy&
     &,shiftdy,xmin,xmax,ymin,ymax,npxb,npyb,nx,ny,kstrt,nvp,ndprof,ierr&
     &)
! initialize electron velocities or momenta
            if (ierr==0) then
! special cases
               if (nvdist==3) then
                  call mpvbdistr2h(part,nps,npp,vtdx,vtdz,vdx,vdz,omx,  &
     &omy,omz,npxb,npyb,kstrt,nvp,ierr)
! general cases
               else
                  call wmpvdistr2h(part,nps,npp,vtdx,vtdy,vtdz,vdx,vdy, &
     &vdz,ci,npxb,npyb,kstrt,nvp,relativity,ierr)
               endif
            endif
! check for beam electron initialization error
            if (ierr /= 0) then
               call PPEXIT(); stop
            endif
         endif
!
! check if any electrons are in the wrong node
         call mpfholes2(part,edges,npp,iholep,ndim,1)
! iholep overflow
         if (iholep(1) < 0) then
            ntmax = -iholep(1)
            ntmax = 1.5*ntmax
            deallocate(iholep)
            if (kstrt==1) then
               write (*,*) 'reallocating electron iholep: ntmax=', ntmax
            endif
            allocate(iholep(ntmax+1))
            call mpfholes2(part,edges,npp,iholep,ndim,1)
            if (iholep(1) < 0) then
               if (kstrt==1) write (*,*) 'iholep overflow: ntmax=',ntmax
               call PPEXIT(); stop
            endif
         endif
! move electrons to correct node
         call ipmove2(part,edges,npp,iholep,ny,tinit,kstrt,nvp,ndim,1,  &
     &ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
!
! find number of electrons in each of mx, my tiles: updates kpic, nppmx
         call mpdblkp2(part,kpic,npp,noff,nppmx,mx,my,mx1,irc)
!
! allocate vector electron data
         nppmx0 = (1.0 + xtras)*nppmx
         ntmaxp = xtras*nppmx
         npbmx = xtras*nppmx
         nbmaxp = 0.25*mx1*npbmx
         allocate(ppart(idimp,nppmx0,mxyp1))
         allocate(ihole(2,ntmaxp+1,mxyp1))
! copy ordered electron data for OpenMP
         call mpmovin2p(part,ppart,kpic,npp,noff,mx,my,mx1,irc)
!
! sanity check for electrons
         call mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,irc)
!
! initialize background charge density: updates qi
         if (movion==0) then
            call mpset_pszero2(qi,tinit,mx,my,mx1,myp1)
! set background to negative of initial electron density
            if (ionbkg==1) then
               call mppost2(ppart,qi,kpic,noff,-qme,tinit,mx,my,mx1,popt&
     &)
               call wmpaguard2(qi,nyp,tinit,nx,kstrt,nvp)
            endif
         endif
!
! initialize ions
         if (movion==1) then
            call mpset_pvzero2(cui,tinit,mx,my,mx1,myp1)
            nps = 1
            nppi = 0
! background ions
            if (npxyi > 0.0d0) then
! calculates initial ion co-ordinates with various density profiles
!              call mpfdistr2(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,&
!    &scaledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,ierr)
! where initial plasma density is confined to some subregion
               call mpgfdistr2(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi&
     &,scaledyi,shiftdyi,xmin,xmax,ymin,ymax,npxi,npyi,nx,ny,kstrt,nvp, &
     &ndprofi,ierr)
! initialize ion velocities or momenta
               if (ierr==0) then
! special cases
                  if (nvdist==3) then
                     call mpvbdistr2h(part,nps,nppi,vtxi,vtzi,vxi0,vzi0,&
     &omx,omy,omz,npxi,npyi,kstrt,nvp,ierr)
! general cases
                  else
                     call wmpvdistr2h(part,nps,nppi,vtxi,vtyi,vtzi,vxi0,&
     &vyi0,vzi0,ci,npxi,npyi,kstrt,nvp,relativity,ierr)
                  endif
               endif
! check for background ion initialization error
               if (ierr /= 0) then
                  call PPEXIT(); stop
               endif
            endif
! beam ions
            if (npxybi > 0.0d0) then
               nps = nppi + 1
! calculates initial ion co-ordinates with various density profiles
!              call mpfdistr2(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,&
!    &scaledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,ierr)
! where initial plasma density is confined to some subregion
               call mpgfdistr2(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi&
     &,scaledyi,shiftdyi,xmin,xmax,ymin,ymax,npxbi,npybi,nx,ny,kstrt,nvp&
     &,ndprofi,ierr)
! initialize ion velocities or momenta
               if (ierr==0) then
! special cases
                  if (nvdist==3) then
                     call mpvbdistr2h(part,nps,nppi,vtdxi,vtdzi,vdxi,   &
     &vdzi,omx,omy,omz,npxbi,npybi,kstrt,nvp,ierr)
! general cases
                  else
                     call wmpvdistr2h(part,nps,nppi,vtdxi,vtdyi,vtdzi,  &
     &vdxi,vdyi,vdzi,ci,npxbi,npybi,kstrt,nvp,relativity,ierr)
                  endif
               endif
! check for beam ion initialization error
               if (ierr /= 0) then
                  call PPEXIT(); stop
               endif
            endif
!
! check if any ions are in the wrong node
            call mpfholes2(part,edges,nppi,iholep,ndim,1)
! iholep overflow
            if (iholep(1) < 0) then
               ntmax = -iholep(1)
               ntmax = 1.5*ntmax
               deallocate(iholep)
               if (kstrt==1) then
                  write (*,*) 'reallocating ion iholep: ntmax=', ntmax
               endif
               allocate(iholep(ntmax+1))
               call mpfholes2(part,edges,nppi,iholep,ndim,1)
               if (iholep(1) < 0) then
                  if (kstrt==1) then
                     write (*,*) 'iholep overflow: ntmax=', ntmax
                  endif
                  call PPEXIT(); stop
               endif
            endif
! move ions to correct node
            call ipmove2(part,edges,nppi,iholep,ny,tinit,kstrt,nvp,ndim,&
     &1,ierr)
            if (ierr /= 0) then
               call PPEXIT(); stop
            endif
!
! find number of ions in each of mx, my tiles: updates kipic, nppmx
            call mpdblkp2(part,kipic,nppi,noff,nppmx,mx,my,mx1,irc)
!
! allocate vector ion data
            nppmx1 = (1.0 + xtras)*nppmx
            allocate(pparti(idimp,nppmx1,mxyp1))
            if (.not.allocated(ihole)) then
               ntmaxp = xtras*nppmx
               npbmx = xtras*nppmx
               nbmaxp = 0.25*mx1*npbmx
               allocate(ihole(2,ntmaxp+1,mxyp1))
            endif
! copy ordered ion data for OpenMP
            call mpmovin2p(part,pparti,kipic,nppi,noff,mx,my,mx1,irc)
!
! sanity check for ions
            call mpcheck2(pparti,kipic,noff,nyp,nx,mx,my,mx1,irc)
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
         if (movion==1) then
            call mpset_pvzero2(cui,tinit,mx,my,mx1,myp1)
         endif
! read first part of data:
! updates ntime, ntime0, part, npp, kpic, nppmx, ierr
         call bread_restart2a(part,kpic,tinit,kstrt,iur0,iscr,ntime,    &
     &ntime0,npp,nppmx,noff,mx1,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
! allocate vector electron data
         nppmx0 = (1.0 + xtras)*nppmx
         ntmaxp = xtras*nppmx
         npbmx = xtras*nppmx
         nbmaxp = 0.25*mx1*npbmx
         allocate(ppart(idimp,nppmx0,mxyp1))
         allocate(ihole(2,ntmaxp+1,mxyp1))
! read second part of data:
! updates ppart, kpic, part, nppi, kipic, nppmx, ierr
         call bread_restart2b(part,ppart,kpic,kipic,tinit,kstrt,iur0,   &
     &iscr,npp,nppi,nppmx,noff,nyp,nx,mx1,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
! allocate vector ion data
         if (movion==1) then
            nppmx1 = (1.0 + xtras)*nppmx
            allocate(pparti(idimp,nppmx1,mxyp1))
            if (.not.allocated(ihole)) then
               ntmaxp = xtras*nppmx
               npbmx = xtras*nppmx
               nbmaxp = 0.25*mx1*npbmx
               allocate(ihole(2,ntmaxp+1,mxyp1))
            endif
         endif
! read third part of data: updates pparti, kipic, qi, ierr
         call bread_restart2c(part,pparti,kipic,qi,tinit,kstrt,iur0,    &
     &ntime,ntime0,nppi,noff,nyp,nx,mx1,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
! read in basic restart file for electromagnetic code:
! updates exyz, bxyz
         call bread_restart23(exyz,bxyz,tinit,kstrt,iur0,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
         if ((ntime+ntime0)> 0) dth = 0.5*dt
      endif
!
! set magnitude of external magnetic field
      omt = sqrt(omx*omx + omy*omy + omz*omz)
!
! kxps/kyps = actual grids used in field partitions in x/y direction
      kxps = min(kxp,max(0,nxh-kxp*(kstrt-1)))
      kyps = min(kyp,max(0,ny-kyp*(kstrt-1)))
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
         allocate(sfieldc(nye,kxp),sfield(nxe,nypmx))
      endif
!
! allocate scratch arrays for vector fields
      if ((ntel>0).or.(ntje>0).or.(nta>0).or.(ntet>0).or.(ntb>0).or.    &
     &(ntar>0).or.(ntji>0)) then
         allocate(vfieldc(ndim,nye,kxp),vfield(ndim,nxe,nypmx))
      endif
!
! initialize electron density diagnostic
      if (ntde > 0) then
         modesxde = min(modesxde,nxh)
         modesyde = min(modesyde,nyh)
         modesy2de = 2*modesyde - 1
         modesxpd = min(modesxde,kxp)
         allocate(denet(modesy2de,modesxpd))
! open file
         if (kstrt==1) then
! open file for complex data: updates nderec and possibly iude
!           fdename = 'denek2.'//cdrun
!           if (nderec==0) then
!              call dafopenc2(denet,iude,nderec,trim(fdename))
!           endif
! open file for real data: updates nderec and possibly iude
            fdename = 'dener2.'//cdrun
            if (nderec==0) then
               call dafopen2(sfield,nx,kyp,iude,nderec,trim(fdename))
            endif
         endif
      endif
!
! initialize ion density diagnostic
      if (movion==1) then
         if (ntdi > 0) then
            modesxdi = min(modesxdi,nxh)
            modesydi = min(modesydi,nyh)
            modesy2di = 2*modesydi - 1
            modesxpd = min(modesxdi,kxp)
            allocate(denit(modesy2di,modesxpd))
! open file
            if (kstrt==1) then
! open file for complex data: updates ndirec and possibly iudi
!              fdiname = 'denik2.'//cdrun
!              if (ndirec==0) then
!                 call dafopenc2(denit,iudi,ndirec,trim(fdiname))
!              endif
! open file for real data: updates ndirec and possibly iudi
               fdiname = 'denir2.'//cdrun
               if (ndirec==0) then
                  call dafopen2(sfield,nx,kyp,iudi,ndirec,trim(fdiname))
               endif
            endif
         endif
      endif
!
! initialize potential diagnostic
      if (ntp > 0) then
         modesxp = min(modesxp,nxh)
         modesyp = min(modesyp,nyh)
         modesy2p = 2*modesyp - 1
         modesxpd = min(modesxp,kxp)
         allocate(pott(modesy2p,modesxpd))
! open file
         if (kstrt==1) then
! open file for complex data: updates nprec and possibly iup
!           fpname = 'potk2.'//cdrun
!           if (nprec==0) call dafopenc2(pott,iup,nprec,trim(fpname))
! open file for real data: updates nprec and possibly iup
            fpname = 'potr2.'//cdrun
            if (nprec==0) then
               call dafopen2(sfield,nx,kyp,iup,nprec,trim(fpname))
            endif
         endif
      endif
!
! initialize longitudinal efield diagnostic
      if (ntel > 0) then
         modesxel = min(modesxel,nxh)
         modesyel = min(modesyel,nyh)
         modesy2el = 2*modesyel - 1
         modesxpd = min(modesxel,kxp)
         allocate(elt(ndim,modesy2el,modesxpd))
! open file
         if (kstrt==1) then
! open file for complex data: updates nelrec and possibly iuel
!           felname = 'elk2.'//cdrun
!           if (nelrec==0) then
!              call dafopenvc2(elt,iuel,nelrec,trim(felname))
!           endif
! open file for real data: updates nelrec and possibly iuel
            felname = 'elr2.'//cdrun
            if (nelrec==0) then
               call dafopenv2(vfield,nx,kyp,iuel,nelrec,trim(felname))
            endif
         endif
      endif
!
! initialize electron current density diagnostic
      if (ntje > 0) then
         modesxje = min(modesxje,nxh)
         modesyje = min(modesyje,nyh)
         modesy2je = 2*modesyje - 1
         modesxpd = min(modesxje,kxp)
         allocate(curet(ndim,modesy2je,modesxpd))
! open file
         if (kstrt==1) then
! open file for complex data: updates njerec and possibly iuje
!           fjename = 'curek2.'//cdrun
!           if (njerec==0) then
!              call dafopenvc2(curet,iuje,njerec,trim(fjename))
!            endif
! open file for real data: updates njerec and possibly iuje
            fjename = 'curer2.'//cdrun
            if (njerec==0) then
               call dafopenv2(vfield,nx,kyp,iuje,njerec,trim(fjename))
            endif
         endif
      endif
!
! initialize ion current density diagnostic
      if (movion==1) then
         if (ntji > 0) then
            modesxji = min(modesxji,nxh)
            modesyji = min(modesyji,nyh)
            modesy2ji = 2*modesyji - 1
            modesxpd = min(modesxji,kxp)
            allocate(curit(ndim,modesy2ji,modesxpd))
! open file
            if (kstrt==1) then
! open file for complex data: updates njirec and possibly iuji
!              fjiname = 'curik2.'//cdrun
!              if (njirec==0) then
!                 call dafopenvc2(curit,iuji,njirec,trim(fjiname))
!              endif
! open file for real data: updates njirec and possibly iuji
               fjiname = 'curir2.'//cdrun
               if (njirec==0) then
                  call dafopenv2(vfield,nx,kyp,iuji,njirec,trim(fjiname)&
     &)
               endif
            endif
         endif
      endif
!
! initialize radiative vector potential diagnostic
      if (ntar > 0) then
         modesxar = min(modesxar,nxh)
         modesyar = min(modesyar,nyh)
         modesy2ar= 2*modesyar - 1
         modesxpd = min(modesxar,kxp)
         allocate(vpotr(ndim,modesy2ar,modesxpd))
         allocate(oldcut(ndim,nye,kxp))
! open file
         if (kstrt==1) then
! open file for complex data: updates narrec and possibly iuar
!           farname = 'vpotrk2.'//cdrun
!           if (narrec==0) then
!              call dafopenvc2(vpotr,iuar,narrec,trim(farname))
!           endif
! open file for real data: updates narrec and possibly iuar
            farname = 'vpotrr2.'//cdrun
            if (narrec==0) then
               call dafopenv2(vfield,nx,kyp,iuar,narrec,trim(farname))
            endif
         endif
      endif
!
! initialize vector potential diagnostic
      if (nta > 0) then
         modesxa = min(modesxa,nxh)
         modesya = min(modesya,ny)
         modesy2a = 2*modesya - 1
         modesxpd = min(modesxa,kxp)
         allocate(vpott(ndim,modesy2a,modesxpd))
! open file
         if (kstrt==1) then
! open file for complex data: updates narec and possibly iua
!           faname = 'vpotk2.'//cdrun
!           if (narec==0) call dafopenvc2(vpott,iua,narec,trim(faname))
! open file for real data: updates narec and possibly iua
            faname = 'vpotr2.'//cdrun
            if (narec==0) then
               call dafopenv2(vfield,nx,kyp,iua,narec,trim(faname))
            endif
         endif
      endif
!
! initialize transverse efield diagnostic
      if (ntet > 0) then
         modesxet = min(modesxet,nxh)
         modesyet = min(modesyet,nyh)
         modesy2et = 2*modesyet - 1
         modesxpd = min(modesxet,kxp)
         allocate(ett(ndim,modesy2et,modesxpd))
! open file
         if (kstrt==1) then
! open file for complex data: updates netrec and possibly iuet
!           fetname = 'etk2.'//cdrun
!           if (netrec==0) then
!              call dafopenvc2(ett,iuet,netrec,trim(fetname))
!           endif
! open file for real data: updates netrec and possibly iuet
            fetname = 'etr2.'//cdrun
            if (netrec==0) then
               call dafopenv2(vfield,nx,kyp,iuet,netrec,trim(fetname))
            endif
         endif
      endif
!
! initialize magnetic field diagnostic
      if (ntb > 0) then
         modesxb = min(modesxb,nxh)
         modesyb = min(modesyb,nyh)
         modesy2b = 2*modesyb - 1
         modesxpd = min(modesxb,kxp)
         allocate(bt(ndim,modesy2b,modesxpd))
! open file
         if (kstrt==1) then
! open file for complex data: updates nbrec and possibly iub
!           fbname = 'bk2.'//cdrun
!           if (nbrec==0) call dafopenvc2(bt,iub,nbrec,trim(fbname))
! open file for real data: updates nbrec and possibly iub
            fbname = 'br2.'//cdrun
            if (nbrec==0) then
               call dafopenv2(vfield,nx,kyp,iub,nbrec,trim(fbname))
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
            allocate(fmse(nprd,nxe,nypmx))
! open file for real data: updates nferec and possibly iufe
            if (kstrt==1) then
               ffename = 'fmer2.'//cdrun
               if (nferec==0) then
                  call dafopenv2(fmse,nx,kyp,iufe,nferec,trim(ffename))
               endif
            endif
         endif
! ion moments
         if (movion==1) then
            if ((ndfm==2).or.(ndfm==3)) then
               allocate(fmsi(nprd,nxe,nypmx))
! open file for real data: updates nfirec and possibly iufi
               if (kstrt==1) then
                  ffiname = 'fmir2.'//cdrun
                  if (nfirec==0) then
                     call dafopenv2(fmsi,nx,kyp,iufi,nfirec,            &
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
            if (npxy.gt.0.0d0) then
               ws = 4.0*vtx+abs(vx0)
               ws = max(ws,4.0*vty+abs(vy0))
               ws = max(ws,4.0*vtz+abs(vz0))
            endif
            if (npxyb.gt.0.0d0) then
               ws = max(ws,4.0*vtdx+abs(vdx))
               ws = max(ws,4.0*vtdy+abs(vdy))
               ws = max(ws,4.0*vtdz+abs(vdz))
            endif
            allocate(fvm(ndim,3),sfv(nmv21+1,ndim,mxyp1))
            allocate(fv(nmv21+1,nfvd),fe(nmv21+1,nfed))
            fvm = 0.0
! open file for electron velocity data: updates nverec and possibly iuve
            if (kstrt==1) then
               fvename = 'fve2.'//cdrun
               if (nverec==0) then
                  call dafopenfv2(fvm,fv,fe,wk,iuve,nverec,             &
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
               if (npxyi.gt.0.0d0) then
                   ws = 4.0*vtxi+abs(vxi0)
                   ws = max(ws,4.0*vtyi+abs(vyi0))
                   ws = max(ws,4.0*vtzi+abs(vzi0))
               endif
               if (npxybi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
                  ws = max(ws,4.0*vtdzi+abs(vdzi))
               endif
! open file for ion velocity data: updates nvirec and possibly iuvi
               if (kstrt==1) then
                  fviname = 'fvi2.'//cdrun
                  if (nvirec==0) then
                     call dafopenfv2(fvmi,fvi,fei,wk,iuvi,nvirec,       &
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
! set electron tags: updates nprobt, tedges, ppart and possibly iprobt
            call psetptraj2(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,vtsx,&
     &dvtx,np,nprobt,irc)
! estimate maximum electron velocity or momentum
            if (nst==3) then
               ws = 0.0
               if (npxy.gt.0.0d0) then
                  ws = 4.0*vtx+abs(vx0)
                  ws = max(ws,4.0*vty+abs(vy0))
                  ws = max(ws,4.0*vtz+abs(vz0))
               endif
               if (npxyb.gt.0.0d0) then
                  ws = max(ws,4.0*vtdx+abs(vdx))
                  ws = max(ws,4.0*vtdy+abs(vdy))
                  ws = max(ws,4.0*vtdz+abs(vdz))
               endif
            endif
! ion trajectories
         else if (ndt==2) then
! set ion tags: updates nprobt, tedges, pparti and possibly iprobt
            call psetptraj2(pparti,tedges,kipic,iprobt,kstrt,nst,vtxi,  &
     &vtsx,dvtx,npi,nprobt,irc)
! estimate maximum ion velocity or momentum
            if (nst==3) then
               ws = 0.0
               if (npxyi.gt.0.0d0) then
                  ws = 4.0*vtxi+abs(vxi0)
                  ws = max(ws,4.0*vtyi+abs(vyi0))
                  ws = max(ws,4.0*vtzi+abs(vzi0))
               endif
               if (npxybi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
                  ws = max(ws,4.0*vtdzi+abs(vdzi))
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
            ftname = 'tr2.'//cdrun
            if ((nst==1).or.(nst==2)) then
               mtt = (nloop - 1)/ntt + 1
               itt = 0
               allocate(partd(mtt,idimp,nprobt))
               partd = 0.0
! open file for trajectory data: updates ntrec and possibly iut
               if (kstrt==1) then
                  if (ntrec==0) then
                     call dafopentr2(partt,iut,ntrec,trim(ftname))
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
                     call dafopenfv2(fvmtp,fvtp,fetp,ws,iut,ntrec,      &
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
         mvx = min(mvx,nx); mvy = min(mvy,nypmn)
         nsxb = (nx - 1)/mvx + 1; nsyb = (ny - 1)/mvy + 1
         nyb = (noff + nyp - 1)/mvy - (noff - 1)/mvy
         if (kstrt==1) nyb = nyb + 1
         itot(1) = nyb
         call mpimax(itot(:1),tinit)
         nybmx = itot(1)
! electron phase space diagnostic
         if ((nds==1).or.(nds==3)) then
! estimate maximum electron velocity or momentum
            ws = 0.0
            if (npxy.gt.0.0d0) then
               ws = 4.0*vtx+abs(vx0)
               ws = max(ws,4.0*vty+abs(vy0))
               ws = max(ws,4.0*vtz+abs(vz0))
            endif
            if (npxyb.gt.0.0d0) then
               ws = max(ws,4.0*vtdx+abs(vdx))
               ws = max(ws,4.0*vtdy+abs(vdy))
               ws = max(ws,4.0*vtdz+abs(vdz))
            endif
            allocate(fvs(nmv21+1,ndim,nsxb,nyb+1))
            fvs(nmv21+1,:,1,1) = 1.25*ws
! open file for electron phase space data:
! updates nserec and possibly iuse
! opens a new fortran unformatted stream file
            if (nserec==0) then
               if (kstrt==1) then
                  fsename = 'pse2.'//cdrun
                  iuse =  get_funit(iuse)
                  call fnopens2(iuse,trim(fsename))
               endif
               nserec = 1
            endif
         endif
! ion phase space diagnostic
         if (movion==1) then
            if ((nds==2).or.(nds==3)) then
! estimate maximum ion velocity or momentum
               ws = 0.0
               if (npxyi.gt.0.0d0) then
                   ws = 4.0*vtxi+abs(vxi0)
                   ws = max(ws,4.0*vtyi+abs(vyi0))
                   ws = max(ws,4.0*vtzi+abs(vzi0))
               endif
               if (npxybi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
                  ws = max(ws,4.0*vtdzi+abs(vdzi))
               endif
               allocate(fvsi(nmv21+1,ndim,nsxb,nyb+1))
               fvsi(nmv21+1,:,1,1) = 1.25*ws
! open file for ion phase space data:
! updates nsirec and possibly iusi
! opens a new fortran unformatted stream file
               if (nsirec==0) then
                  if (kstrt==1) then
                     fsiname = 'psi2.'//cdrun
                     iusi =  get_funit(iusi)
                     call fnopens2(iusi,trim(fsiname))
                  endif
                  nsirec = 1
               endif
            endif
         endif
      endif
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = 0.0
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (dt > 0.45*ci) then
         if (kstrt==1) then
            write (*,*) 'Warning: Courant condition may be exceeded!'
         endif
      endif
!
      if (kstrt==1) write (iuot,*) 'program mpbbeps2'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
!
! print time step
      if (kstrt==1) then
         if (ntw > 0) then
            it = ntime/ntw
            if (ntime==ntw*it) write (iuot,*) 'ntime = ', ntime
         endif
      endif
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
               call mpset_pvzero2(fmse,tdiag,mx,my,mx1,myp1)
               call wmprofx23(ppart,fmse,kpic,noff,ci,tdiag,npro,mx,my, &
     &mx1,relativity)
! add guard cells with OpenMP: updates fmse
               call wmpnacguard2(fmse,nyp,tdiag,nx,kstrt,nvp)
! moves vector grid fmse from non-uniform to uniform partition
               isign = -1
               call mpfnmove2(fmse,noff,nyp,isign,tdiag,kyp,ny,kstrt,nvp&
     &,nterf,irc)
               if (irc /= 0) then
                  call PPEXIT(); stop
               endif
! calculates fluid quantities from fluid moments: updates fmsi
               call mpfluidqs23(fmse,tdiag,npro,nx,ny,kstrt,kyp)
! write real space diagnostic output: updates nferec
               call mpvwrite2(fmse,tdiag,nx,ny,kyp,iufe,nferec)
!              call mpvread2(fmse,tdiag,nx,ny,kyp,iufe,nferec,irc)
            endif
! calculate ion fluid moments
            if (movion==1) then
               if ((ndfm==2).or.(ndfm==3)) then
                  call mpset_pvzero2(fmsi,tdiag,mx,my,mx1,myp1)
                  call wmprofx23(pparti,fmsi,kipic,noff,ci,tdiag,npro,mx&
     &,my,mx1,relativity)
                  fmsi = rmass*fmsi
! add guard cells with OpenMP: updates fmsi
                  call wmpnacguard2(fmsi,nyp,tdiag,nx,kstrt,nvp)
! moves vector grid fmsi from non-uniform to uniform partition
                  isign = -1
                  call mpfnmove2(fmsi,noff,nyp,isign,tdiag,kyp,ny,kstrt,&
     &nvp,nterf,irc)
                  if (irc /= 0) then
                     call PPEXIT(); stop
                  endif
! calculates fluid quantities from fluid moments: updates fmsi
                  call mpfluidqs23(fmsi,tdiag,npro,nx,ny,kstrt,kyp)
! write real space diagnostic output: updates nfirec
                  call mpvwrite2(fmsi,tdiag,nx,ny,kyp,iufi,nfirec)
!                 call mpvread2(fmsi,tdiag,nx,ny,kyp,iufi,nfirec,irc)
               endif
            endif
         endif
      endif
!
! deposit electron current with OpenMP:
! updates ppart and cue, and possibly ncl, ihole, irc
      call mpset_pvzero2(cue,tdjpost,mx,my,mx1,myp1)
      call wmpdjpost2(ppart,cue,kpic,ncl,ihole,noff,nyp,qme,dth,ci,     &
     &tdjpost,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
!
! add guard cells with OpenMP: updates cue
      call wmpnacguard2(cue,nyp,tguard,nx,kstrt,nvp)
!
! reorder electrons by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
      if (irc==0) then
         call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov,  &
     &kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,plist,irc2)
      else
         irc2(1) = 1; irc2(2) = irc; irc = 0
      endif
!
      do while (irc2(1) /= 0)
! ihole overflow
         if (irc2(1)==1) then
            ntmaxp = (1.0 + xtras)*irc2(2)
            deallocate(ihole)
            allocate(ihole(2,ntmaxp+1,mxyp1))
            irc2 = 0
            call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov&
     &,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,.false.,irc2)
! ppart overflow
         else if (irc2(1)==4) then
! restores electron coordinates from ppbuff: updates ppart, ncl
            call mprstor2(ppart,ppbuff,ncl,ihole,tsort)
! copy ordered electrons to linear array: updates part
            call mpcopyout2(part,ppart,kpic,it,irc)
! part overflow
            if (irc > 0) then
               deallocate(part)
               npmax = irc
               maxnp = max(npmax,npimax)
               allocate(part(idimp,maxnp))
               irc = 0
               call mpcopyout2(part,ppart,kpic,it,irc)
            endif
            deallocate(ppart)
            nppmx0 = (1.0 + xtras)*irc2(2)
            allocate(ppart(idimp,nppmx0,mxyp1))
! copies unordered electrons to ordered array: updates ppart
            call mpcopyin2(part,ppart,kpic,irc)
            call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov&
     &,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,plist,irc2)
         endif
      enddo
!
! sanity check for electrons
      if (monitor > 0) then
         call mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,irc)
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
            call wmpfft2rn(vfield,cut,noff,nyp,isign,mixup,sct,tfft,    &
     &tfmov,indx,indy,kstrt,nvp,kyp,ny,nterf,ierr)
! calculate smoothed electron current in fourier space: updates vfieldc
            call mpsmooth23(cut,vfieldc,ffc,tfield,nx,ny,kstrt)
! store selected fourier modes: updates curet
            call mprdvmodes2(vfieldc,curet,tfield,nx,ny,modesxje,       &
     &modesyje,kstrt)
! write fourier space diagnostic output: updates njerec
!           call mpvcwrite2(curit,tdiag,modesxje,modesy2je,kxp,iuje,    &
!    &njerec)
!           call mpvcread2(curet,tdiag,modesxje,modesy2je,kxp,iuje,     &
!    &njerec,irc)
! transform smoothed electron current to real space: updates vfield
            isign = 1
            call mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &kstrt,nvp,kyp)
! write real space diagnostic output: updates njerec
            call mpvwrite2(vfield,tdiag,nx,ny,kyp,iuje,njerec)
!           call mpvread2(vfield,tdiag,nx,ny,kyp,iuje,njerec,irc)
! move vector data to non-uniform partition, updates: vfield, ierr
            isign = 1
            call mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,&
     &nterf,ierr)
! display smoothed electron current
            call wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
            call pdvector2(vfield,nyp,nvp,' ELECTRON CURRENT',ntime,999,&
     &0,ndstyle,1,nx,ny,ierr)
            if (ierr==1) then
                  call PPEXIT(); stop
            endif
            ierr = 0
         endif
      endif
!
! deposit ion current with OpenMP:
      if (movion==1) then
! updates pparti and cui, and possibly ncl, ihole, irc
         call mpset_pvzero2(cui,tdjpost,mx,my,mx1,myp1)
         call wmpdjpost2(pparti,cui,kipic,ncl,ihole,noff,nyp,qmi,dth,ci,&
     &tdjpost,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
! add guard cells with OpenMP: updates cui
         call wmpnacguard2(cui,nyp,tguard,nx,kstrt,nvp)
!
! reorder ions by tile with OpenMP and MPI
! updates: pparti, kipic, and irc and possibly ncl and ihole
         if (irc==0) then
            call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort,  &
     &tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,plist,irc2)
         else
            irc2(1) = 1; irc2(2) = irc; irc = 0
         endif
!
         do while (irc2(1) /= 0)
! ihole overflow
            if (irc2(1)==1) then
               ntmaxp = (1.0 + xtras)*irc2(2)
               deallocate(ihole)
               allocate(ihole(2,ntmaxp+1,mxyp1))
               irc2 = 0
               call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort&
     &,tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,.false.,irc2)
! pparti overflow
            else if (irc2(1)==4) then
! restores ion coordinates from ppbuff: updates pparti, ncl
               call mprstor2(pparti,ppbuff,ncl,ihole,tsort)
! copy ordered ions to linear array: updates part
               call mpcopyout2(part,pparti,kipic,it,irc)
! part overflow
               if (irc > 0) then
                  deallocate(part)
                  npimax = irc
                  maxnp = max(npmax,npimax)
                  allocate(part(idimp,maxnp))
                  irc = 0
                  call mpcopyout2(part,pparti,kipic,it,irc)
               endif
               deallocate(pparti)
               nppmx1 = (1.0 + xtras)*irc2(2)
               allocate(pparti(idimp,nppmx1,mxyp1))
! copies unordered ions to ordered array: updates pparti
               call mpcopyin2(part,pparti,kipic,irc)
               call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort&
     &,tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,plist,irc2)
            endif
         enddo
!
! sanity check for ions
         if (monitor > 0) then
            call mpcheck2(pparti,kipic,noff,nyp,nx,mx,my,mx1,irc)
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
               call wmpfft2rn(vfield,cut,noff,nyp,isign,mixup,sct,tfft, &
     &tfmov,indx,indy,kstrt,nvp,kyp,ny,nterf,ierr)
! calculate smoothed ion current in fourier space: updates vfieldc
               call mpsmooth23(cut,vfieldc,ffc,tfield,nx,ny,kstrt)
! store selected fourier modes: updates curit
               call mprdvmodes2(vfieldc,curit,tfield,nx,ny,modesxji,    &
     &modesyji,kstrt)
! write fourier space diagnostic output: updates njirec
!              call mpvcwrite2(curit,tdiag,modesxji,modesy2ji,kxp,iuji, &
!    &njirec)
!              call mpvcread2(curit,tdiag,modesxji,modesy2ji,kxp,iuji,  &
!    &njirec,irc)
! transform smoothed ion current to real space: updates vfield
               isign = 1
               call mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,  &
     &indy,kstrt,nvp,kyp)
! write real space diagnostic output: updates njirec
               call mpvwrite2(vfield,tdiag,nx,ny,kyp,iuji,njirec)
!              call mpvread2(vfield,tdiag,nx,ny,kyp,iuji,njirec,irc)
! move vector data to non-uniform partition, updates: vfield, ierr
               isign = 1
               call mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt, &
     &nvp,nterf,ierr)
! display smoothed ion current
               call wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
               call pdvector2(vfield,nyp,nvp,'ION CURRENT',ntime,999,0, &
     &ndstyle,1,nx,ny,ierr)
               if (ierr==1) then
                  call PPEXIT(); stop
               endif
               ierr = 0
            endif
         endif
      endif
!
! deposit electron charge with OpenMP: updates qe
      call mpset_pszero2(qe,tdpost,mx,my,mx1,myp1)
      call mppost2(ppart,qe,kpic,noff,qme,tdpost,mx,my,mx1,popt)
! add guard cells with OpenMP: updates qe
      call wmpaguard2(qe,nyp,tguard,nx,kstrt,nvp)
!
! electron density diagnostic
      if (ntde > 0) then
         it = ntime/ntde
         if (ntime==ntde*it) then
            sfield = -qe
! transform electron density to fourier space: updates qt
! moves data to uniform partition
            isign = -1
            call wmpfft2r(sfield,qt,noff,nyp,isign,mixup,sct,tfft,tfmov,&
     &indx,indy,kstrt,nvp,kyp,ny,nterf,ierr)
! calculate smoothed electron density in fourier space: updates sfieldc
            call mpsmooth2(qt,sfieldc,ffc,tfield,nx,ny,kstrt)
! store selected fourier modes: updates denet
            call mprdmodes2(sfieldc,denet,tfield,nx,ny,modesxde,modesyde&
     &,kstrt)
! write fourier space diagnostic output: updates nderec
!           call mpcwrite2(denet,tdiag,modesxde,modesy2de,kxp,iude,     &
!    &nderec)
!           call mpcread2(denet,tdiag,modesxde,modesy2de,kxp,iude,nderec&
!    &,irc)
! transform smoothed electron density to real space: updates sfield
            isign = 1
            call mpfft2r(sfield,sfieldc,isign,mixup,sct,tfft,indx,indy, &
     &kstrt,nvp,kyp)
! write real space diagnostic output: updates nderec
            call mpwrite2(sfield,tdiag,nx,ny,kyp,iude,nderec)
!           call mpread2(sfield,tdiag,nx,ny,kyp,iude,nderec,irc)
! move data to non-uniform partition, updates: sfield, ierr
            isign = 1
            call mpfmove2(sfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp, &
     &nterf,ierr)
! display smoothed electron density
            call wmpcguard2(sfield,nyp,tguard,nx,kstrt,nvp)
            call pdscaler2(sfield,nyp,nvp,'EDENSITY',ntime,999,1,ndstyle&
     &,nx,ny,ierr)
            if (ierr==1) then
               call PPEXIT(); stop
            endif
            ierr = 0
         endif
      endif
!
! deposit ion charge with OpenMP: updates qi
      if (movion==1) then
         call mpset_pszero2(qi,tdpost,mx,my,mx1,myp1)
         call mppost2(pparti,qi,kipic,noff,qmi,tdpost,mx,my,mx1,popt)
! add guard cells with OpenMP: updates qi
         call wmpaguard2(qi,nyp,tguard,nx,kstrt,nvp)
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
               call wmpfft2r(sfield,qt,noff,nyp,isign,mixup,sct,tfft,   &
     &tfmov,indx,indy,kstrt,nvp,kyp,ny,nterf,ierr)
! calculate smoothed ion density in fourier space: updates sfieldc
               call mpsmooth2(qt,sfieldc,ffc,tfield,nx,ny,kstrt)
! store selected fourier modes: updates denit
               call mprdmodes2(sfieldc,denit,tfield,nx,ny,modesxdi,     &
     &modesydi,kstrt)
! write fourier space diagnostic output: updates ndirec
!              call mpcwrite2(denit,tdiag,modesxdi,modesy2di,kxp,iudi,  &
!    &ndirec)
!              call mpcread2(denit,tdiag,modesxdi,modesy2di,kxp,iudi,   &
!    &ndirec,irc)
! transform smoothed ion density to real space: updates sfield
               isign = 1
               call mpfft2r(sfield,sfieldc,isign,mixup,sct,tfft,indx,   &
     &indy,kstrt,nvp,kyp)
! write real space diagnostic output: updates ndirec
               call mpwrite2(sfield,tdiag,nx,ny,kyp,iudi,ndirec)
!              call mpread2(sfield,tdiag,nx,ny,kyp,iudi,ndirec,irc)
! move data to non-uniform partition, updates: sfield, ierr
               isign = 1
               call mpfmove2(sfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,  &
     &nvp,nterf,ierr)
! display smoothed ion density
               call wmpcguard2(sfield,nyp,tguard,nx,kstrt,nvp)
               call pdscaler2(sfield,nyp,nvp,'ION DENSITY',ntime,999,1, &
     &ndstyle,nx,ny,ierr)
               if (ierr==1) then
                  call PPEXIT(); stop
               endif
               ierr = 0
            endif
         endif
      endif
!
! add electron and ion densities: updates qe
      call mpaddqei2(qe,qi,nyp,tfield,nx)
!
! add electron and ion current densities: updates cue
      if (movion==1) call mpaddcuei2(cue,cui,nyp,tfield,nx)
!
! transform charge to fourier space with OpenMP:
! moves data to uniform partition
! updates qt, nterf, and ierr, modifies qe
      isign = -1
      call wmpfft2r(qe,qt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,indy,&
     &kstrt,nvp,kyp,ny,nterf,ierr)
!
! transform current to fourier space with OpenMP:
! moves data to uniform partition
! updates cut, nterf, and ierr, modifies cue
      isign = -1
      call wmpfft2rn(cue,cut,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,  &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! take transverse part of current with OpenMP: updates cut
      call mpcuperp2(cut,tfield,nx,ny,kstrt)
!
! radiative vector potential diagnostic
      if (ntar > 0) then
         it = ntime/ntar
         if (ntime==ntar*it) then
! average current: updates vfieldc = 0.5*(cut + oldcut)
            call mcuave2(vfieldc,cut,oldcut,tfield,ny)
! calculate radiative vector potential in fourier space: updates vfieldc
! vfieldc should contain averaged current on entry
            call mpavrpot2(vfieldc,bxyz,ffc,affp,ci,tfield,nx,ny,kstrt)
! store selected fourier modes: updates vpotr
            call mprdvmodes2(vfieldc,vpotr,tfield,nx,ny,modesxar,       &
     &modesyar,kstrt)
! write fourier space diagnostic output: updates narrec
!           call mpvcwrite2(vpotr,tdiag,modesxar,modesy2ar,kxp,iuar,    &
!    &narrec)
!           call mpvcread2(vpotr,tdiag,modesxar,modesy2ar,kxp,iuar,     &
!    &narrec,irc)
! transform radiative vector potential to real space: updates vfield
            isign = 1
            call mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &kstrt,nvp,kyp)
! write real space diagnostic output: updates narrec
            call mpvwrite2(vfield,tdiag,nx,ny,kyp,iuar,narrec)
!           call mpvread2(vfield,tdiag,nx,ny,kyp,iuar,narrec,irc)
! move vector data to non-uniform partition, updates: vfield, ierr
            isign = 1
            call mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,&
     &nterf,ierr)
! display radiative vector potential
            call wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
            call pdvector2(vfield,nyp,nvp,'RADIATIVE VPOTENTIAL',ntime, &
     &999,0,ndstyle,1,nx,ny,ierr)
            if (ierr==1) then
               call PPEXIT(); stop
            endif
            ierr = 0
         endif
      endif
!
! calculate electromagnetic fields in fourier space with OpenMP:
! updates exyz, bxyz, wf, wb
      if ((ntime+ntime0)==0) then
! initialize electromagnetic fields from darwin fields
! calculate initial darwin magnetic field
         call mpibpois2(cut,bxyz,ffc,ci,wb,tfield,nx,ny,kstrt)
         wf = 0.0
! calculate initial darwin electric field
         allocate(amu(4,nxe,nypmx),amut(4,nye,kxp),dcut(ndim,nye,kxp))
         call mpset_pvzero2(amu,tdjpost,mx,my,mx1,myp1)
         call wmpgmjpost2(ppart,amu,kpic,noff,qme,ci,tdjpost,mx,my,mx1, &
     &popt,relativity)
         call wmpnacguard2(amu,nyp,tguard,nx,kstrt,nvp)
         isign = -1
         call wmpfft2rn(amu,amut,noff,nyp,isign,mixup,sct,tfft,tfmov,   &
     &indx,indy,kstrt,nvp,kyp,ny,nterf,ierr)
         deallocate(amu)
         call mpdcuperp2(dcut,amut,tfield,nx,ny,kstrt)
         deallocate(amut)
         call mpetfield2(dcut,exyz,ffc,affp,ci,wf,tfield,nx,ny,kstrt)
         deallocate(dcut)
         dth = 0.5*dt
! update electromagnetic fields
      else
         call mpmaxwel2(exyz,bxyz,cut,ffc,affp,ci,dt,wf,wb,tfield,nx,ny,&
     &kstrt)
      endif
!
! calculate longitudinal force/charge in fourier space with OpenMP:
! updates fxyt, we
      call mppois2(qt,fxyt,ffc,we,tfield,nx,ny,kstrt)
!
! add longitudinal and transverse electric fields with OpenMP:
! updates fxyt
      isign = 1
      call mpemfield2(fxyt,exyz,ffc,isign,tfield,nx,ny,kstrt)
! copy magnetic field with OpenMP: updates bxyt
      isign = -1
      call mpemfield2(bxyt,bxyz,ffc,isign,tfield,nx,ny,kstrt)
!
! transform force to real space with OpenMP:
! moves data to non-uniform partition
! updates fxyze, nterf, and ierr, modifies fxyt
      isign = 1
      call wmpfft2rn(fxyze,fxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! transform magnetic field to real space with OpenMP:
! moves data to non-uniform partition
! updates bxyze, nterf, and ierr, modifies bxyt
      isign = 1
      call wmpfft2rn(bxyze,bxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! add constant to magnetic field with OpenMP: updates bxyze
      if (omt > 0.0) call mpbaddext2(bxyze,nyp,tfield,omx,omy,omz,nx)
!
! copy guard cells with OpenMP: updates fxyze, bxyze
      call wmpncguard2(fxyze,nyp,tguard,nx,kstrt,nvp)
      call wmpncguard2(bxyze,nyp,tguard,nx,kstrt,nvp)
!
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! calculate potential in fourier space: updates sfieldc
            call mppot2(qt,sfieldc,ffc,ws,tfield,nx,ny,kstrt)
! store selected fourier modes: updates pott
            call mprdmodes2(sfieldc,pott,tfield,nx,ny,modesxp,modesyp,  &
     &kstrt)
! write fourier space diagnostic output: updates nprec
!           call mpcwrite2(pott,tdiag,modesxp,modesy2p,kxp,iup,nprec)
!           call mpcread2(pott,tdiag,modesxp,modesy2p,kxp,iup,nprec,irc)
! transform potential to real space: updates sfield
            isign = 1
            call mpfft2r(sfield,sfieldc,isign,mixup,sct,tfft,indx,indy, &
     &kstrt,nvp,kyp)
! write real space diagnostic output: updates nprec
            call mpwrite2(sfield,tdiag,nx,ny,kyp,iup,nprec)
!           call mpread2(sfield,tdiag,nx,ny,kyp,iup,nprec,irc)
! move data to non-uniform partition, updates: sfield, ierr
            isign = 1
            call mpfmove2(sfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp, &
     &nterf,ierr)
! display potential
            call wmpcguard2(sfield,nyp,tguard,nx,kstrt,nvp)
            call pdscaler2(sfield,nyp,nvp,'POTENTIAL',ntime,999,0,      &
     &ndstyle,nx,ny,ierr)
            if (ierr==1) then
               call PPEXIT(); stop
            endif
            ierr = 0
         endif
      endif
!
! longitudinal efield diagnostic
      if (ntel > 0) then
         it = ntime/ntel
         if (ntime==ntel*it) then
! calculate longitudinal efield in fourier space: updates vfieldc
            call mpelfield2(qt,vfieldc,ffc,ws,tfield,nx,ny,kstrt)
! store selected fourier modes: updates elt
            call mprdvmodes2(vfieldc,elt,tfield,nx,ny,modesxel,modesyel,&
     &kstrt)
! write fourier space diagnostic output: updates nelrec
!           call mpvcwrite2(elt,tdiag,modesxel,modesy2el,kxp,iuel,nelrec&
!    &)
!           call mpvcread2(elt,tdiag,modesxel,modesy2el,kxp,iuel,nelrec,&
!    &irc)
! transform longitudinal efield to real space: updates vfield
            isign = 1
            call mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &kstrt,nvp,kyp)
! write real space diagnostic output: updates nelrec
            call mpvwrite2(vfield,tdiag,nx,ny,kyp,iuel,nelrec)
!           call mpvread2(vfield,tdiag,nx,ny,kyp,iuel,nelrec,irc)
! move data to non-uniform partition, updates: vfield, ierr
            isign = 1
            call mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,&
     &nterf,ierr)
! display longitudinal efield 
            call wmpncguard2(vfield,nyp,tguard,nx,kstrt,nvp)
            call pdvector2(vfield,nyp,nvp,'ELFIELD',ntime,999,0,ndstyle,&
     &1,nx,ny,ierr)
            if (ierr==1) then
               call PPEXIT(); stop
            endif
            ierr = 0
         endif
      endif
!
! vector potential diagnostic
      if (nta > 0) then
         it = ntime/nta
         if (ntime==nta*it) then
! calculate vector potential in fourier space: updates vfieldc
            call mpavpot2(bxyz,vfieldc,tfield,nx,ny,kstrt)
! store selected fourier modes: updates vpott
            call mprdvmodes2(vfieldc,vpott,tfield,nx,ny,modesxa,modesya,&
     &kstrt)
! write fourier space diagnostic output: updates narec
!           call mpvcwrite2(vpott,tdiag,modesxa,modesy2a,kxp,iua,narec)
!           call mpvcread2(vpott,tdiag,modesxa,modesy2a,kxp,iua,narec,  &
!    &irc)
! transform vector potential to real space: updates vfield
            isign = 1
            call mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &kstrt,nvp,kyp)
! write real space diagnostic output: updates narec
            call mpvwrite2(vfield,tdiag,nx,ny,kyp,iua,narec)
!           call mpvread2(vfield,tdiag,nx,ny,kyp,iua,narec,irc)
! move vector data to non-uniform partition, updates: vfield, ierr
            isign = 1
            call mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,&
     &nterf,ierr)
! display vector potential
            call wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
            call pdvector2(vfield,nyp,nvp,'VECTOR POTENTIAL',ntime,999,0&
     &,ndstyle,1,nx,ny,ierr)
            if (ierr==1) then
               call PPEXIT(); stop
            endif
            ierr = 0
         endif
      endif
!
! transverse efield diagnostic
      if (ntet > 0) then
         it = ntime/ntet
         if (ntime==ntet*it) then
            vfieldc = exyz
! store selected fourier modes: updates ett
            call mprdvmodes2(vfieldc,ett,tfield,nx,ny,modesxet,modesyet,&
     &kstrt)
! write fourier space diagnostic output: updates netrec
!           call mpvcwrite2(ett,tdiag,modesxet,modesy2et,kxp,iuet,netrec&
!    &)
!           call mpvcread2(ett,tdiag,modesxet,modesy2et,kxp,iuet,netrec,&
!    &irc)
! transform transverse efield to real space: updates vfield
            isign = 1
            call mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &kstrt,nvp,kyp)
! write real space diagnostic output: updates netrec
            call mpvwrite2(vfield,tdiag,nx,ny,kyp,iuet,netrec)
!           call mpvread2(vfield,tdiag,nx,ny,kyp,iuet,netrec,irc)
! move data to non-uniform partition, updates: vfield, ierr
            isign = 1
            call mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,&
     &nterf,ierr)
! display transverse efield 
            call wmpncguard2(vfield,nyp,tguard,nx,kstrt,nvp)
            call pdvector2(vfield,nyp,nvp,'TRANSVERSE EFIELD',ntime,999,&
     &0,ndstyle,1,nx,ny,ierr)
            if (ierr==1) then
               call PPEXIT(); stop
            endif
            ierr = 0
         endif
      endif
!
! magnetic field diagnostic
      if (ntb > 0) then
         it = ntime/ntb
         if (ntime==ntb*it) then
            vfieldc = bxyz
! store selected fourier modes: updates bt
            call mprdvmodes2(vfieldc,bt,tfield,nx,ny,modesxb,modesyb,   &
     &kstrt)
! write fourier space diagnostic output: updates nbrec
!           call mpvcwrite2(bt,tdiag,modesxb,modesy2b,kxp,iub,nbrec)
!           call mpvcread2(bt,tdiag,modesxb,modesy2b,kxp,iub,nbrec,irc)
! transform magnetic field to real space: updates vfield
            isign = 1
            call mpfft2rn(vfield,vfieldc,isign,mixup,sct,tfft,indx,indy,&
     &kstrt,nvp,kyp)
! write real space diagnostic output: updates nbrec
            call mpvwrite2(vfield,tdiag,nx,ny,kyp,iub,nbrec)
!           call mpvread2(vfield,tdiag,nx,ny,kyp,iub,nbrec,irc)
! move vector data to non-uniform partition, updates: vfield, ierr
            isign = 1
            call mpfnmove2(vfield,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,&
     &nterf,ierr)
! display magnetic field
            call wmpncguard2(vfield,nyp,tfield,nx,kstrt,nvp)
            call pdvector2(vfield,nyp,nvp,'MAGNETIC FIELD',ntime,999,0, &
     &ndstyle,1,nx,ny,ierr)
            if (ierr==1) then
               call PPEXIT(); stop
            endif
            ierr = 0
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
                  call mpvpdist2(ppart,kpic,fv,sfv,fvm,tdiag,nvp,nmv)
! store time history electron vdrift, vth, and entropy
                  itv = itv + 1
                  fvtm(itv,:,:) = fvm
! display electron velocity distributions
                  call pdisplayfv1(fv,fvm,' ELECTRON',kstrt,ntime,nmv,2,&
     &ierr)
                  if (ierr==1) then
                     call PPEXIT(); stop
                  endif
                  ierr = 0
               endif
! calculate electron cylindrical distribution function and moments
               if ((nvft==4).or.(nvft==5)) then
                  call mpvbpdist2(ppart,kpic,fv,sfv,fvm,omx,omy,omz,    &
     &tdiag,nmv)
! display electron velocity distributions in cylindrical co-ordinates
                  call pdisplayfvb1(fv,fvm,' ELECTRON',kstrt,ntime,nmv, &
     &2,ierr)
                  if (ierr==1) then
                     call PPEXIT(); stop
                  endif
                  ierr = 0
               endif
! electron energy distribution
               if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
                  call mperpdist2(ppart,kpic,fe,sfv,eci,wk,tdiag,ndim,  &
     &nmv)
! display electron energy distribution
                  call pdisplayfe1(fe,wk,' ELECTRON',kstrt,ntime,nmv,   &
     &ierr)
                  if (ierr==1) then
                     call PPEXIT(); stop
                  endif
                  ierr = 0
               endif
! write electron velocity-space diagnostic output: updates nverec
               if (kstrt==1) then
                  call dafwritefv2(fvm,fv,fe,wk,tdiag,iuve,nverec)
               endif
            endif
! ion distribution functions
            if (movion==1) then
               if ((ndv==2).or.(ndv==3)) then
! calculate ion cartesian distribution function and moments
                  if ((nvft==1).or.(nvft==3)) then
                     call mpvpdist2(pparti,kipic,fvi,sfv,fvmi,tdiag,nvp,&
     &nmv)
! store time history of ion vdrift, vth, and entropy
                     fvtmi(itv,:,:) = fvmi
! display ion velocity distributions
                     call pdisplayfv1(fvi,fvmi,' ION',kstrt,ntime,nmv,2,&
     &ierr)
                     if (ierr==1) then
                        call PPEXIT(); stop
                     endif
                     ierr = 0
                  endif
! calculate ion cylindrical distribution function and moments
                  if ((nvft==4).or.(nvft==5)) then
                     call mpvbpdist2(pparti,kipic,fvi,sfv,fvmi,omx,omy, &
     &omz,tdiag,nmv)
! display ion velocity distributions in cylindrical co-ordinates
                     call pdisplayfvb1(fvi,fvmi,' ION',kstrt,ntime,nmv,2&
     &,ierr)
                     if (ierr==1) then
                        call PPEXIT(); stop
                     endif
                     ierr = 0
                  endif
! ion energy distribution
                  if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
                     call mperpdist2(pparti,kipic,fei,sfv,eci,wk,tdiag, &
     &ndim,nmv)
                     wk = rmass*wk
! display ion energy distribution
                     ts = fei(nmv21+1,1)
                     fei(nmv21+1,1) = rmass*fei(nmv21+1,1)
                     call pdisplayfe1(fei,wk,' ION',kstrt,ntime,nmv,ierr)
                     fei(nmv21+1,1) = ts
                     if (ierr==1) then
                        call PPEXIT(); stop
                     endif
                     ierr = 0
                  endif
! write ion velocity-space diagnostic output: updates nvirec
                  if (kstrt==1) then
                     call dafwritefv2(fvmi,fvi,fei,wk,tdiag,iuvi,nvirec)
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
               call mptraj2(ppart,kpic,partt,tdiag,numtp,ierr)
! copies tagged ions in ppart to array partt: updates partt, numtp
            else if (ndt==2) then
               if (movion==1) then
                  call mptraj2(pparti,kipic,partt,tdiag,numtp,ierr)
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
                     call mptraj2(ppart,kpic,partt,tdiag,numtp,ierr)
! copies tagged ions in ppart to array partt: updates partt, numtp
                  else if (ndt==2) then
                     if (movion==1) then
                        call mptraj2(pparti,kipic,partt,tdiag,numtp,ierr&
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
! determines list of tagged pareticles leaving this node: updates iholep
                  call mpfholes2(partt,tedges,numtp,iholep,ndim,2)
! iholep overflow
                  if (iholep(1) < 0) then
                     ntmax = -iholep(1)
                     ntmax = 1.5*ntmax
                     deallocate(iholep)
                     if (kstrt==1) then
                        write (*,*) 'info:reallocating iholep:ntmax=',  &
     &ntmax
                     endif
                     allocate(iholep(ntmax+1))
                     call mpfholes2(partt,tedges,numtp,iholep,ndim,2)
                     if (iholep(1) < 0) then
                        if (kstrt==1) then
                           write (*,*) 'iholep overflow: ntmax=', ntmax
                           call PPEXIT(); stop
                        endif
                     endif
                  endif
! copies tagged particles: updates part
                  call mpcpytraj2(partt,part,tdiag,numtp)
! moves tagged electrons into original spatial region:
! updates part, numtp
                  call ipmove2(part,tedges,numtp,iholep,ny,tmov,kstrt,  &
     &nvp,ndim,2,ierr)
                  if (ierr /= 0) then
                     call PPEXIT(); stop
                  endif
! reorders tagged particles: updates partt
                  call mpordtraj2(part,partt,tedges,tdiag,numtp,ierr)
! collects distributed test particle data onto node 0
                  call mppartt2(partt,tdiag,numtp,ierr)
                  if (ierr /= 0) then
                     call PPEXIT(); stop
                  endif
! write trajectory diagnostic output: updates ntrec
                  if (kstrt==1) then
                     call dafwritetr2(partt,tdiag,iut,ntrec)
                     itt = itt + 1
                     partd(itt,:,:) = partt
                  endif
               else if (nst==3) then
! calculate test particle distribution function and moments
                  call mpvdist2(partt,fvtp,fvmtp,tdiag,numtp,nvp,nmv)
! write test particle diagnostic output: updates ntrec
                  if (kstrt==1) then
                     ws = 0.0
                     call dafwritefv2(fvmtp,fvtp,fetp,ws,tdiag,iut,ntrec&
     &)
                  endif
! display test particle velocity distributions
                  if (ndt==1) then
                     call pdisplayfv1(fvtp,fvmtp,' ELECTRON',kstrt,ntime&
     &,nmv,1,ierr)
                  else if (ndt==2) then
                     call pdisplayfv1(fvtp,fvmtp,' ION',kstrt,ntime,nmv,&
     &1,ierr)
                  endif
                  if (ierr==1) then
                     call PPEXIT(); stop
                  endif
                  ierr = 0
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
               call mpvspdist2(ppart,kpic,fvs,tdiag,noff,nmv,mvx,mvy)
! adjusts 3d velocity distribution in different regions of space:
! updates fvs
               call mpadjfvs2(fvs,tdiag,noff,nyp,nmv,mvy)
! write phase space diagnostic output: updates nserec
               if (nserec > 0) then
                  call mpwrfvsdata2(fvs,tdiag,nyb,nybmx,iuse)
                  nserec = nserec + 1
               endif
            endif
! ion phase space
            if (movion==1) then
               if ((nds==2).or.(nds==3)) then
! calculates velocity distribution in different regions of space:
! updates fvsi
                  call mpvspdist2(pparti,kipic,fvsi,tdiag,noff,nmv,mvx, &
     &mvy)
! adjusts 3d velocity distribution in different regions of space:
! updates fvsi
                  call mpadjfvs2(fvsi,tdiag,noff,nyp,nmv,mvy)
! write phase space diagnostic output: updates nsirec
                  if (nsirec > 0) then
                     call mpwrfvsdata2(fvsi,tdiag,nyb,nybmx,iusi)
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
      call wmpbpush2(ppart,fxyze,bxyze,kpic,ncl,ihole,noff,nyp,qbme,dt, &
     &dth,ci,wke,tpush,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
!
! reorder electrons by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
      if (irc==0) then
         call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov,  &
     &kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,plist,irc2)
      else
         irc2(1) = 1; irc2(2) = irc; irc = 0
      endif
!
      do while (irc2(1) /= 0)
! ihole overflow
         if (irc2(1)==1) then
            ntmaxp = (1.0 + xtras)*irc2(2)
            deallocate(ihole)
            allocate(ihole(2,ntmaxp+1,mxyp1))
            irc2 = 0
            call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov&
     &,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,.false.,irc2)
! ppart overflow
         else if (irc2(1)==4) then
! restores electron coordinates from ppbuff: updates ppart, ncl
            call mprstor2(ppart,ppbuff,ncl,ihole,tsort)
! copy ordered electrons to linear array: updates part
            call mpcopyout2(part,ppart,kpic,it,irc)
! part overflow
            if (irc > 0) then
               deallocate(part)
               npmax = irc
               maxnp = max(npmax,npimax)
               allocate(part(idimp,maxnp))
               irc = 0
               call mpcopyout2(part,ppart,kpic,it,irc)
            endif
            deallocate(ppart)
            nppmx0 = (1.0 + xtras)*irc2(2)
            allocate(ppart(idimp,nppmx0,mxyp1))
! copies unordered electrons to ordered array: updates ppart
            call mpcopyin2(part,ppart,kpic,irc)
            call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov&
     &,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,plist,irc2)
         endif
      enddo
!
! sanity check for electrons
      if (monitor > 0) then
         call mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,irc)
      endif
!
! push ions with OpenMP:
      if (movion==1) then
! updates pparti and wki, and possibly ncl, ihole, irc
         wki = 0.0
         call wmpbpush2(pparti,fxyze,bxyze,kipic,ncl,ihole,noff,nyp,qbmi&
     &,dt,dth,ci,wki,tpush,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,  &
     &irc)
         wki = wki*rmass
!
! reorder ions by tile with OpenMP and MPI
! updates: pparti, kipic, and irc and possibly ncl and ihole
         if (irc==0) then
            call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort,  &
     &tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,plist,irc2)
         else
            irc2(1) = 1; irc2(2) = irc; irc = 0
         endif
!
         do while (irc2(1) /= 0)
! ihole overflow
            if (irc2(1)==1) then
               ntmaxp = (1.0 + xtras)*irc2(2)
               deallocate(ihole)
               allocate(ihole(2,ntmaxp+1,mxyp1))
               irc2 = 0
               call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort&
     &,tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,.false.,irc2)
! pparti overflow
            else if (irc2(1)==4) then
! restores ion coordinates from ppbuff: updates pparti, ncl
               call mprstor2(pparti,ppbuff,ncl,ihole,tsort)
! copy ordered ions to linear array: updates part
               call mpcopyout2(part,pparti,kipic,it,irc)
! part overflow
               if (irc > 0) then
                  deallocate(part)
                  npimax = irc
                  maxnp = max(npmax,npimax)
                  allocate(part(idimp,maxnp))
                  irc = 0
                  call mpcopyout2(part,pparti,kipic,it,irc)
               endif
               deallocate(pparti)
               nppmx1 = (1.0 + xtras)*irc2(2)
               allocate(pparti(idimp,nppmx1,mxyp1))
! copies unordered ions to ordered array: updates pparti
               call mpcopyin2(part,pparti,kipic,irc)
               call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort&
     &,tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,popt,plist,irc2)
            endif
         enddo
!
! sanity check for ions
         if (monitor > 0) then
            call mpcheck2(pparti,kipic,noff,nyp,nx,mx,my,mx1,irc)
         endif
      endif
!
! start running simulation backwards:
! need to advance maxwell field solver one step ahead
      if (treverse==1) then
         if (((ntime+1)==(nloop/2)).or.((ntime+1)==nloop)) then
! deposit electron current: updates cue
            call mpset_pvzero2(cue,tdjpost,mx,my,mx1,myp1)
            call wmpdjpost2(ppart,cue,kpic,ncl,ihole,noff,nyp,qme,zero, &
     &ci,tdjpost,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
            call wmpnacguard2(cue,nyp,tguard,nx,kstrt,nvp)
! deposit ion current: updates cui
            if (movion==1) then
               call mpset_pvzero2(cui,tdjpost,mx,my,mx1,myp1)
               call wmpdjpost2(pparti,cui,kipic,ncl,ihole,noff,nyp,qmi, &
     &zero,ci,tdjpost,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
               call wmpnacguard2(cui,nyp,tguard,nx,kstrt,nvp)
               call mpaddcuei2(cue,cui,nyp,tfield,nx)
            endif
            isign = -1
            call wmpfft2rn(cue,cut,noff,nyp,isign,mixup,sct,tfft,tfmov, &
     &indx,indy,kstrt,nvp,kyp,ny,nterf,ierr)
            call mpcuperp2(cut,tfield,nx,ny,kstrt)
! updates exyz, bxyz, wf, wb
            call mpmaxwel2(exyz,bxyz,cut,ffc,affp,ci,dt,wf,wb,tfield,nx,&
     &ny,kstrt)
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
! write out basic restart file for electrostatic code
            call bwrite_restart2(part,ppart,pparti,qi,kpic,kipic,tdiag, &
     &kstrt,iur,iscr,n,ntime0,irc)
! write out basic restart file for electromagnetic code
            call bwrite_restart23(exyz,bxyz,tdiag,kstrt,iur)
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
! reset graphs
      if ((ntw > 0).or.(ntt > 0).or.(ntv > 0)) then
         if (nplot > 0) call reset_pgraphs(kstrt,irc)
      endif
!
      if (kstrt.eq.1) then
         write (iuot,*)
         write (iuot,*) 'ntime, relativity = ', ntime, relativity
         if (treverse==1) write (iuot,*) 'treverse = ', treverse
         write (iuot,*) 'MPI nodes nvp = ', nvp
      endif
!
! trajectory diagnostic
      if (ntt > 0) then
         if ((ndt==1).or.(ndt==2)) then
            if ((nst==1).or.(nst==2)) then
               ts = t0 + dt*real(ntt)
! displays time history of trajectories on node 0
               if (kstrt==1) then
                  if (nplot > 0) call reset_nplot(1,irc)
                  call pdisplaytr1(partd,ts,dt*real(ntt),kstrt,itt,3,999&
     &,ierr)
                  if (ierr==1) then
                     call PPEXIT(); stop
                  endif
                  ierr = 0
                  if (nplot > 0) call reset_nplot(nplot,irc)
               endif
            endif
         endif
      endif
!
! energy diagnostic
      if (ntw > 0) then
         ts = t0 + dt*real(ntw)
         call pdisplayw1(wt,ts,dt*real(ntw),kstrt,itw,ierr)
         if (kstrt==1) then
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
      endif
!
! velocity-space diagnostic
      if (ntv > 0) then
         if ((nvft==1).or.(nvft==3)) then
            ts = t0 + dt*real(ntv)
            call pdisplayfvt1(fvtm,' ELECT',ts,dt*real(ntv),kstrt,itv,  &
     &ierr)
! ions
            if (movion==1) then
               call pdisplayfvt1(fvtmi,' ION',ts,dt*real(ntv),kstrt,itv,&
     &ierr)
            endif
         endif
      endif
!
      if (kstrt==1) then
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
         call writnml2(iudm)
         close(unit=iudm)
! close restart files
         call close_restart2(iur,iur0)
! close output file
         write (iuot,*) ' * * * q.e.d. * * *'
         close(unit=iuot)
! close graphics device
         call close_pgraphs
      endif
!
      call PPEXIT()
      end program
