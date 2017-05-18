!-----------------------------------------------------------------------
! 3D Electrostatic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpbeps3
      use in3
      use modmpinit3
      use modmppush3
      use modmpfield3
      use mpdiag3
      use mppmod3
      use omplib
      use ompplib3
      implicit none
!
! idimp = number of particle coordinates = 6
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 6, ipbc = 1
! idps = number of partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      integer :: idps = 4, idds =    2
! wke/wki/we = particle kinetic/electric field
      real :: wke = 0.0, wki = 0.0, we = 0.0
! plist = (true,false) = list of particles leaving tiles found in push
      logical :: plist = .true.
!
! declare scalars for standard code
      integer :: n
      integer :: nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh
      integer :: nxyzh, nxhyz, mx1, ntime, nloop, isign, ierr
      real :: qbme, affp, ws
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
      integer :: it, mtw
      integer :: itw
! default Fortran unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19
      integer :: iude = 10, iup = 11, iuel = 12
      integer :: iudi = 20
! dimensions for fourier data
      integer :: modesxpd, modesz2de, modesz2p, modesz2el, modesz2di
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe/qi = electron/ion charge density with guard cells
      real, dimension(:,:,:), allocatable :: qe, qi
! fxyze = smoothed electric field with guard cells
      real, dimension(:,:,:,:), allocatable :: fxyze
! qt = scalar charge density field array in fourier space
      complex, dimension(:,:,:), allocatable :: qt
! fxyzt = vector electric field array in fourier space
      complex, dimension(:,:,:,:), allocatable :: fxyzt
! ffc = form factor array for poisson solver
      complex, dimension(:,:,:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
      double precision, dimension(4) :: wtot, work
      integer, dimension(2) :: mterf
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
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0, tmov = 0.0, tfmov = 0.0
      real :: tdiag = 0.0
      real, dimension(2) :: tfft = 0.0
      double precision :: dtime
!
! start timing initialization
      call dtimer(dtime,itime,-1)
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
! nvp = number of MPI ranks
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
!
! read namelists
      if (kstrt==1) then
         call readnml3(iuin)
! override input data
         idcode = 1
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
! nx/ny/nz = number of grid points in x/y/z direction
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
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
         vtdzi = vtz/sqrt(rmass*rtempdzi)
      endif
!
! obtain 2D partition (nvpy,nvpz) from nvp:
! nvpy/nvpz = number of processors in y/z
      call FCOMP32(nvp,nx,ny,nz,nvpy,nvpz,ierr)
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'FCOMP32 error: nvp,nvpy,nvpz=', nvp, nvpy, nvpz
         endif
         call PPEXIT()
         stop
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
      call mpdcomp3(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,&
     &nvpy,nvpz)
! find new 2d partition from initial analytic distribution function
      call mpfedges3(edges,nyzp,noff,ampdy,scaledy,shiftdy,ampdz,scaledz&
     &,shiftdz,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy,nvpz,ipbc,ndprof&
     &,ierr)
      if (ierr /= 0) then
         call PPEXIT()
         stop
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
      allocate(fxyze(ndim,nxe,nypmx,nzpmx))
      allocate(qt(nze,kxyp,kyzp))
      allocate(fxyzt(ndim,nze,kxyp,kyzp))
      allocate(ffc(nzh,kxyp,kyzp),mixup(nxhyz),sct(nxyzh))
      allocate(kpic(mxyzp1),iholep(ntmax+1,2))
!
! prepare fft tables
      call mpfft3_init(mixup,sct,indx,indy,indz)
! calculate form factors
      call mppois3_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,nvpy,nvpz)
!
! initialize electrons
      nps = 1
      npp = 0
! background electrons
      if (npxyz > 0.0d0) then
! calculates initial electron co-ordinates and velocities with
! uniform density and maxwellian velocity with drift
!        call mpdistr3(part,edges,npp,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,  &
!    &npz,nx,ny,nz,kstrt,ipbc,ierr)
! calculates initial electron co-ordinates with various density profiles
         call mpfdistr3(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,   &
     &shiftdy,ampdz,scaledz,shiftdz,npx,npy,npz,nx,ny,nz,kstrt,nvpy,nvpz&
     &,ipbc,ndprof,ierr)
! initialize electron velocities or momenta
         if (ierr==0) then
            call wmpvdistr3(part,nps,npp,vtx,vty,vtz,vx0,vy0,vz0,ci,npx,&
     &npy,npz,kstrt,nvpy,nvpz,relativity,ierr)
         endif
! check for background electron initialization error
         if (ierr /= 0) then
            call PPEXIT()
            stop
         endif
      endif
! beam electrons
      if (npxyzb > 0.0d0) then
         nps = npp + 1
! calculates initial electron co-ordinates and velocities with
! uniform density and maxwellian velocity with drift
!        call mpdistr3(part,edges,npp,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,  &
!    &npyb,npzb,nx,ny,nz,kstrt,ipbc,ierr)
! calculates initial electron co-ordinates with various density profiles
         call mpfdistr3(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,   &
     &shiftdy,ampdz,scaledz,shiftdz,npxb,npyb,npzb,nx,ny,nz,kstrt,nvpy, &
     &nvpz,ipbc,ndprof,ierr)
! initialize electron velocities or momenta
         if (ierr==0) then
            call wmpvdistr3(part,nps,npp,vtdx,vtdy,vtdz,vdx,vdy,vdz,ci, &
     &npxb,npyb,npzb,kstrt,nvpy,nvpz,relativity,ierr)
         endif
! check for beam electron initialization error
         if (ierr /= 0) then
            call PPEXIT()
            stop
         endif
      endif
!
! check if any electrons are in the wrong node
      call mpfholes3(part,edges,npp,iholep)
! iholep overflow
      if (iholep(1,1) < 0) then
         ntmax = -iholep(1,1)
         ntmax = 1.5*ntmax
         deallocate(iholep)
         if (kstrt==1) then
            write (*,*) 'reallocating electron iholep: ntmax=', ntmax
         endif
         allocate(iholep(ntmax+1,2))
         call mpfholes3(part,edges,npp,iholep)
         if (iholep(1,1) < 0) then
            if (kstrt==1) write (*,*) 'iholep overflow: ntmax=', ntmax
            call PPEXIT()
            stop
         endif
      endif
! more electrons to correct node
      call ipmove3(part,edges,npp,iholep,ny,nz,tmov,kstrt,nvpy,nvpz,ierr&
     &)
      if (ierr /= 0) then
         call PPEXIT()
         stop
      endif
!
! delete ipmove3 buffers, no longer needed
      if (movion==0) then
         deallocate(iholep)
         call mppdelszbuf()
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
      nbmaxp = 0.25*mxzyp1*npbmx
      allocate(ppart(idimp,nppmx0,mxyzp1))
      allocate(ncl(26,mxyzp1),ihole(2,ntmaxp+1,mxyzp1))
! copy ordered electron data for OpenMP: updates ppart and kpic
      call mpmovin3(part,ppart,kpic,npp,noff,mx,my,mz,mx1,myp1,irc)
!
! sanity check for electrons
      call mpcheck3(ppart,kpic,noff,nyzp,nx,mx,my,mz,mx1,myp1,irc)
!
! initialize background charge density: updates qi
      if (movion==0) then
         qi = 0.0
         call mppost3(ppart,qi,kpic,noff,-qme,tdpost,mx,my,mz,mx1,myp1)
         call wmpaguard3(qi,nyzp,tguard,nx,kstrt,nvpy,nvpz)
      endif
!
! initialize ions
      if (movion==1) then
         allocate(kipic(mxyzp1))
         nps = 1
         nppi = 0
! background ions
         if (npxyzi > 0.0d0) then
! calculates initial ion co-ordinates and velocities with
! uniform density and maxwellian velocity with drift
!           call mpdistr3(part,edges,nppi,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,&
!    &npxi,npyi,npzi,nx,ny,nz,kstrt,ipbc,ierr)
! calculates initial ion co-ordinates with various density profiles
            call mpfdistr3(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,   &
     &scaledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxi,npyi,npzi,nx,ny,nz&
     &,kstrt,nvpy,nvpz,ipbc,ndprofi,ierr)
! initialize ion velocities or momenta
            if (ierr==0) then
               call wmpvdistr3(part,nps,nppi,vtxi,vtyi,vtzi,vxi0,vyi0,  &
     &vzi0,ci,npxi,npyi,npzi,kstrt,nvpy,nvpz,relativity,ierr)
            endif
! check for background ion initialization error
            if (ierr /= 0) then
               call PPEXIT()
               stop
            endif
         endif
! beam ions
         if (npxyzbi > 0.0d0) then
! calculates initial ion co-ordinates and velocities with
! uniform density and maxwellian velocity with drift
!           call mpdistr3(part,edges,nppi,vtdxi,vtdyi,vtdzi,vdxi,vdyi,  &
!    &vdzi,npxbi,npybi,npzbi,nx,ny,nz,kstrt,ipbc,ierr)
! calculates initial ion co-ordinates with various density profiles
            call mpfdistr3(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,   &
     &scaledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxbi,npybi,npzbi,nx,ny&
     &,nz,kstrt,nvpy,nvpz,ipbc,ndprofi,ierr)
! initialize ion velocities or momenta
            if (ierr==0) then
               call wmpvdistr3(part,nps,nppi,vtdxi,vtdyi,vtdzi,vdxi,vdyi&
     &,vdzi,ci,npxbi,npybi,npzbi,kstrt,nvpy,nvpz,relativity,ierr)
            endif
! check for background ion initialization error
            if (ierr /= 0) then
               call PPEXIT()
               stop
            endif
         endif
!
! check if any ions are in the wrong node
         call mpfholes3(part,edges,nppi,iholep)
! iholep overflow
         if (iholep(1,1) < 0) then
            ntmax = -iholep(1,1)
            ntmax = 1.5*ntmax
            deallocate(iholep)
            if (kstrt==1) then
               write (*,*) 'reallocating ion iholep: ntmax=', ntmax
            endif
            allocate(iholep(ntmax+1,2))
            call mpfholes3(part,edges,nppi,iholep)
            if (iholep(1,1) < 0) then
               if (kstrt==1) write (*,*) 'iholep overflow: ntmax=',ntmax
               call PPEXIT()
               stop
            endif
         endif
! more ions to correct node
         call ipmove3(part,edges,nppi,iholep,ny,nz,tmov,kstrt,nvpy,nvpz,&
     &ierr)
         if (ierr /= 0) then
            call PPEXIT()
            stop
         endif
!
! delete ipmove3 buffers, no longer needed
         deallocate(iholep)
         call mppdelszbuf()
!
! find number of ions in each of mx, my, mz tiles:
! updates kipic, nppmx
         call mpdblkp3(part,kipic,nppi,noff,nppmx,mx,my,mz,mx1,myp1,irc)
!
! allocate vector ion data
         nppmx1 = (1.0 + xtras)*nppmx
         allocate(pparti(idimp,nppmx1,mxyzp1))
! copy ordered ion data for OpenMP: updates pparti and kipic
         call mpmovin3(part,pparti,kipic,nppi,noff,mx,my,mz,mx1,myp1,irc&
     &)
!
! sanity check for ions
         call mpcheck3(pparti,kipic,noff,nyzp,nx,mx,my,mz,mx1,myp1,irc)
      endif
!
! allocate diagnostic arrays
! reverse simulation at end back to start
      if (treverse==1) nloop = 2*nloop
!
! energy time history
      if (ntw > 0) then
         mtw = (nloop - 1)/ntw + 1; itw = 0
         allocate(wt(mtw,4),s(4))
         wt = 0.0; s = 0.0d0
      endif
!
! allocate scratch arrays for scalar fields
      if ((ntde > 0).or.(ntp > 0).or.(ntdi > 0)) then
         allocate(sfieldc(nze,kxyp,kyzp),sfield(nxe,nypmx,nzpmx))
      endif
!
! allocate scratch arrays for vector fields
      if (ntel > 0) then
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
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (kstrt==1) write (iuot,*) 'program mpbeps3'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
      if (kstrt==1) write (iuot,*) 'ntime = ', ntime
!
! deposit charge with OpenMP: updates qe
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
! transform charge to fourier space with OpenMP:
! moves data to uniform partition
! updates qt, mterf, and ierr, modifies qe
      isign = -1
      call wmpfft3r(qe,qt,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx,indy&
     &,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! calculate force/charge in fourier space with OpenMP: updates fxyzt, we
      call mppois3(qt,fxyzt,ffc,we,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! transform force to real space with OpenMP:
! moves data to non-uniform partition
! updates fxyze, mterf, and ierr, modifies fxyzt
      isign = 1
      call wmpfft3rn(fxyze,fxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! copy guard cells with OpenMP: updates fxyze
      call wmpncguard3(fxyze,nyzp,tguard,nx,kstrt,nvpy,nvpz)
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
! push electrons with OpenMP:
! updates ppart and wke, and possibly ncl, ihole, irc
      wke = 0.0
      call wmppush3(ppart,fxyze,kpic,ncl,ihole,noff,nyzp,qbme,dt,ci,wke,&
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
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
         call wmppush3(pparti,fxyze,kipic,ncl,ihole,noff,nyzp,qbmi,dt,ci&
     &,wki,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
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
! need to reverse time lag in leap-frog integration scheme
      if (treverse==1) then
         if (((ntime+1)==(nloop/2)).or.((ntime+1)==nloop)) then
            dt = -dt
            ws = 0.0
            call wmppush3zf(ppart,kpic,ncl,ihole,noff,nyzp,dt,ci,ws,    &
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
            call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,xtras,tsort,   &
     &tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,&
     &mxzyp1,plist,irc2)
            if (movion==1) then
               call wmppush3zf(pparti,kipic,ncl,ihole,noff,nyzp,dt,ci,ws&
     &,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
               call ompmove3(pparti,kipic,ncl,ihole,noff,nyzp,xtras,    &
     &tsort,tmov,kstrt,nvpy,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1&
     &,mzp1,mxzyp1,plist,irc2)
            endif
         endif
      endif
!
! energy diagnostic
      if (ntw > 0) then
         it = ntime/ntw
         if (ntime==ntw*it) then
            wtot(1) = we
            wtot(2) = wke
            wtot(3) = wki
            wtot(4) = we + wke
            call PPDSUM(wtot,work,4)
            we = wtot(1)
            wke = wtot(2)
            wki = wtot(3)
            ws = we + wke + wki
            if (ntime==0) s(3) = ws
            if (kstrt==1) then
               write (iuot,*) 'Field, Kinetic and Total Energies:'
               if (movion==0) then
                  write (iuot,'(3e14.7)') we, wke, ws
               else
                  write (iuot,'(4e14.7)') we, wke, wki, ws
               endif
            endif
            itw = itw + 1
! store energies in time history array
            wt(itw,:) = (/we,wke,wki,ws/)
            s(1) = s(1) + we
            s(2) = s(2) + wke
            s(3) = min(s(3),dble(ws))
            s(4) = max(s(4),dble(ws))
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
!           ts = t0 + dt*real(ntw)
!           call displayw1(wt,ts,dt*real(ntw),itw,irc)
            s(3) = (s(4) - s(3))/wt(1,4)
            write (iuot,*) 'Energy Conservation = ', real(s(3))
            s(1) = s(1)/real(itw)
            write (iuot,*) 'Average Field Energy <WE> = ', real(s(1))
            s(2) = s(2)/real(itw)
            write (iuot,*) 'Average Electron Kinetic Energy <WKE> = ',  &
     &real(s(2))
            write (iuot,*) 'Ratio <WE>/<WKE>= ', real(s(1)/s(2))
         endif
!
         write (iuot,*)
         write (iuot,*) 'initialization time = ', tinit
         write (iuot,*) 'deposit time = ', tdpost
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
         if (movion==1) then
! ion density diagnostic
            if (ntdi > 0) then
               ndirec = ndirec - 1
            endif
         endif
! write final diagnostic metafile
         call writnml3(iudm)
         write (iuot,*) ' * * * q.e.d. * * *'
         close(unit=iudm)
         close(unit=iuot)
      endif
!
      call PPEXIT()
      end program
