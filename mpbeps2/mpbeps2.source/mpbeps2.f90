!-----------------------------------------------------------------------
! 2D Electrostatic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpbeps2
      use in2
      use modmpinit2
      use modmppush2
      use modmpfield2
      use mpdiag2
      use mppmod2
      use pgraf2
      use pgraf1
      use omplib
      use ompplib2
      use mpsimul2
      implicit none
!
! idimp = dimension of phase space = 4
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 4, ipbc = 1
! idps = number of partition boundaries
      integer :: idps = 2
! wke/wki/we = particle kinetic/electric field
      real :: wke = 0.0, wki = 0.0, we = 0.0
! plist = (true,false) = list of particles leaving tiles found in push
      logical :: plist = .true.
!
! declare scalars for standard code
      integer :: n
      integer :: nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy
      integer :: mx1, ntime, nloop, isign, ierr
      integer :: ntime0 = 0
      real :: qbme, affp, ws
      real :: qbmi, vtxi, vtyi, vtdxi, vtdyi
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
      real :: ts, eci, wk
! default Fortran unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19
      integer :: iur = 17, iur0 = 27, iscr = 99
      integer :: iude = 10, iup = 11, iuel = 12
      integer :: iufe = 23, iuve = 25, iut = 28, iuse = 29
      integer :: iudi = 20, iufi = 24, iuvi = 26, iusi = 30
! dimensions for fourier data
      integer :: modesxpd, modesy2de, modesy2p, modesy2el, modesy2di
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! declare arrays for standard code
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe/qi = electron/ion charge density with guard cells
      real, dimension(:,:), allocatable :: qe, qi
! fxye = smoothed electric field with guard cells
      real, dimension(:,:,:), allocatable :: fxye
! qt = scalar charge density field array in fourier space
      complex, dimension(:,:), allocatable :: qt
! fxyt = vector electric field array in fourier space
      complex, dimension(:,:,:), allocatable :: fxyt
! ffc = form factor array for poisson solver
      complex, dimension(:,:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
      integer, dimension(1) :: itot
      double precision, dimension(4) :: wtot
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
      real :: tdpost = 0.0, tguard = 0.0, tfield = 0.0
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
! nvp = number of distributed memory nodes
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
!
! read namelists
      if (kstrt==1) then
         call readnml2(iuin)
! override input data
         idcode = 1
         ndim = 2
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
! np = total number of electrons in simulation
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
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
         vtdyi = vtdy/sqrt(rmass*rtempdyi)
      endif
!
! check if too many node
      if (nvp > ny) then
         if (kstrt==1) then
         write (*,*) 'Too many nodes requested: ny, nvp=', ny, nvp
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
      call mpfedges2(edges,nyp,noff,ampdy,scaledy,shiftdy,nypmx,nypmn,ny&
     &,kstrt,nvp,ipbc,ndprof,ierr)
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
! allocate data for standard code
      allocate(part(idimp,maxnp))
      allocate(qe(nxe,nypmx),qi(nxe,nypmx),fxye(ndim,nxe,nypmx))
      allocate(qt(nye,kxp),fxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),mixup(nxhy),sct(nxyh))
      allocate(kpic(mxyp1),ncl(8,mxyp1),iholep(ntmax+1))
      if (movion==1) allocate(kipic(mxyp1))
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
! calculates initial electron co-ordinates with uniform density
!           call mpudistr2(part,edges,npp,npx,npy,nx,ny,kstrt,ipbc,ierr)
! calculates initial electron co-ordinates with various density profiles
            call mpfdistr2(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,&
     &shiftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,ierr)
! initialize electron velocities or momenta
            if (ierr==0) then
               call wmpvdistr2(part,nps,npp,vtx,vty,vx0,vy0,ci,npx,npy, &
     &kstrt,nvp,relativity,ierr)
            endif
! check for background electron initialization error
            if (ierr /= 0) then
               call PPEXIT(); stop
            endif
         endif
! beam electrons
         if (npxyb > 0.0d0) then
            nps = npp + 1
! calculates initial electron co-ordinates with uniform density
!           call mpudistr2(part,edges,npp,npxb,npyb,nx,ny,kstrt,ipbc,   &
!    &ierr)
! calculates initial electron co-ordinates with various density profiles
            call mpfdistr2(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,&
     &shiftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,ierr)
! initialize electron velocities or momenta
            if (ierr==0) then
               call wmpvdistr2(part,nps,npp,vtdx,vtdy,vdx,vdy,ci,npxb,  &
     &npyb,kstrt,nvp,relativity,ierr)
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
         call ipmove2(part,edges,npp,iholep,ny,tmov,kstrt,nvp,ndim,1,   &
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
!
! copy ordered electron data for OpenMP
         call mpmovin2(part,ppart,kpic,npp,noff,mx,my,mx1,irc)
!
! sanity check for electrons
         call mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,irc)
!
! initialize background charge density: updates qi
         if (movion==0) then
            qi = 0.0
            call mppost2(ppart,qi,kpic,noff,-qme,tdpost,mx,my,mx1)
            call wmpaguard2(qi,nyp,tguard,nx,kstrt,nvp)
         endif
!
! initialize ions
         if (movion==1) then
            nps = 1
            nppi = 0
! background ions
            if (npxyi > 0.0d0) then
! calculates initial ion co-ordinates with uniform density
!              call mpudistr2(part,edges,nppi,npxi,npyi,nx,ny,kstrt,ipbc&
!    &,ierr)
! calculates initial ion co-ordinates with various density profiles
               call mpfdistr2(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,&
     &scaledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,ierr)
! initialize ion velocities or momenta
               if (ierr==0) then
                  call wmpvdistr2(part,nps,nppi,vtxi,vtyi,vxi0,vyi0,ci, &
     &npxi,npyi,kstrt,nvp,relativity,ierr)
               endif
! check for background ion initialization error
               if (ierr /= 0) then
                  call PPEXIT(); stop
               endif
            endif
! beam ions
            if (npxybi > 0.0d0) then
               nps = nppi + 1
! calculates initial ion co-ordinates with uniform density
!              call mpudistr2(part,edges,nppi,npxbi,npybi,nx,ny,kstrt,  &
!    &ipbc,ierr)
! calculates initial ion co-ordinates with various density profiles
               call mpfdistr2(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,&
     &scaledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,ierr)
! initialize ion velocities or momenta
               if (ierr==0) then
                  call wmpvdistr2(part,nps,nppi,vtdxi,vtdyi,vdxi,vdyi,ci&
     &,npxbi,npybi,kstrt,nvp,relativity,ierr)
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
            call ipmove2(part,edges,nppi,iholep,ny,tmov,kstrt,nvp,ndim,1&
     &,ierr)
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
!
! copy ordered ion data for OpenMP
            call mpmovin2(part,pparti,kipic,nppi,noff,mx,my,mx1,irc)
!
! sanity check for ions
            call mpcheck2(pparti,kipic,noff,nyp,nx,mx,my,mx1,irc)
         endif
!
! restart to continue a run which was interrupted
      else if (nustrt==2) then
         write (*,*) 'nustrt = 2 not yet supported'
         call PPEXIT(); stop
!        nstart = ntime + 1
! start a new run with data from a previous run
      else if (nustrt==0) then
! read in basic restart file for electrostatic code
! read first part of data:
! updates ntime, ntime0, part, npp, kpic, nppmx, ierr
         call bread_restart2a(part,kpic,tdiag,kstrt,iur0,iscr,ntime,    &
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
         call bread_restart2b(part,ppart,kpic,kipic,tdiag,kstrt,iur0,   &
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
         call bread_restart2c(part,pparti,kipic,qi,tdiag,kstrt,iur0,    &
     &ntime,ntime0,nppi,noff,nyp,nx,mx1,ierr)
         if (ierr /= 0) then
            call PPEXIT(); stop
         endif
      endif
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
         allocate(wt(mtw,4),s(4))
         wt = 0.0; s = 0.0d0
      endif
!
! allocate scratch arrays for scalar fields
      if ((ntde > 0).or.(ntp > 0).or.(ntdi > 0)) then
         allocate(sfieldc(nye,kxp),sfield(nxe,nypmx))
      endif
!
! allocate scratch arrays for vector fields
      if (ntel > 0) then
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
! initialize fluid moments diagnostic
      if (ntfm > 0) then
         nprd = 0
         if (npro==1) then
            nprd = 1
         else if (npro==2) then
            nprd = 3
         else if (npro==3) then
            nprd = 6
         else if (npro==4) then
            nprd = 9
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
         endif
         if ((nvft==2).or.(nvft==3)) then
            nfed = 1
         endif
         nmv21 = 2*nmv + 1
         mtv = (nloop - 1)/ntv + 1; itv = 0
         eci = ci; if (relativity==0) eci = 0.0
         wk = 0.0
! electron velocity diagnostic
         if ((ndv==1).or.(ndv==3)) then
! estimate maximum electron velocity or momentum
            ws = 0.0
            if (npxy.gt.0.0d0) then
               ws = 4.0*vtx+abs(vx0)
               ws = max(ws,4.0*vty+abs(vy0))
            endif
            if (npxyb.gt.0.0d0) then
               ws = max(ws,4.0*vtdx+abs(vdx))
               ws = max(ws,4.0*vtdy+abs(vdy))
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
! energy distribution
            if ((nvft==2).or.(nvft==3)) then
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
               endif
               if (npxybi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
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
! energy distribution
               if ((nvft==2).or.(nvft==3)) then
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
               endif
               if (npxyb.gt.0.0d0) then
                  ws = max(ws,4.0*vtdx+abs(vdx))
                  ws = max(ws,4.0*vtdy+abs(vdy))
               endif
            endif
! ion trajectories
         else if (ndt==2) then
! set ione tags: updates nprobt, tedges, pparti and possibly iprobt
            call psetptraj2(pparti,tedges,kipic,iprobt,kstrt,nst,vtxi,  &
     &vtsx,dvtx,npi,nprobt,irc)
! estimate maximum ion velocity or momentum
            if (nst==3) then
               ws = 0.0
               if (npxyi.gt.0.0d0) then
                  ws = 4.0*vtxi+abs(vxi0)
                  ws = max(ws,4.0*vtyi+abs(vyi0))
               endif
               if (npxybi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
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
         call mpimax(itot(:1),tdiag)
         nybmx = itot(1)
! electron phase space diagnostic
         if ((nds==1).or.(nds==3)) then
! estimate maximum electron velocity or momentum
            ws = 0.0
            if (npxy.gt.0.0d0) then
               ws = 4.0*vtx+abs(vx0)
               ws = max(ws,4.0*vty+abs(vy0))
            endif
            if (npxyb.gt.0.0d0) then
               ws = max(ws,4.0*vtdx+abs(vdx))
               ws = max(ws,4.0*vtdy+abs(vdy))
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
               endif
               if (npxybi.gt.0.0d0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
                  ws = max(ws,4.0*vtdyi+abs(vdyi))
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
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (kstrt==1) write (iuot,*) 'program mpbeps2'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
      if (kstrt==1) write (iuot,*) 'ntime = ', ntime
!
! deposit electron charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      call mppost2(ppart,qe,kpic,noff,qme,tdpost,mx,my,mx1)
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
         call dtimer(dtime,itime,-1)
         qi = 0.0
         call dtimer(dtime,itime,1)
         tdpost = tdpost + real(dtime)
         call mppost2(pparti,qi,kipic,noff,qmi,tdpost,mx,my,mx1)
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
! transform charge to fourier space with OpenMP:
! moves data to uniform partition
! updates qt, nterf, and ierr, modifies qe
      isign = -1
      call wmpfft2r(qe,qt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,indy,&
     &kstrt,nvp,kyp,ny,nterf,ierr)
!
! calculate force/charge in fourier space with OpenMP: updates fxyt, we
      call mppois2(qt,fxyt,ffc,we,tfield,nx,ny,kstrt)
!
! transform force to real space with OpenMP:
! moves data to non-uniform partition
! updates fxye, nterf, and ierr, modifies fxyt
      isign = 1
      call wmpfft2rn(fxye,fxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! copy guard cells with OpenMP: updates fxye
      call wmpncguard2(fxye,nyp,tguard,nx,kstrt,nvp)
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
               call wmgprofx2(ppart,fxye,fmse,kpic,noff,nyp,qbme,dt,ci, &
     &tdiag,npro,nx,mx,my,mx1,relativity)
! add guard cells with OpenMP: updates fmse
               call wmpnacguard2(fmse,nyp,tdiag,nx,kstrt,nvp)
! moves vector grid fmse from non-uniform to uniform partition
               isign = -1
               call mpfnmove2(fmse,noff,nyp,isign,tdiag,kyp,ny,kstrt,nvp&
     &,nterf,irc)
               if (irc /= 0) then
                  call PPEXIT(); stop
               endif
! calculates fluid quantities from fluid moments
               call mpfluidqs22(fmse,tdiag,npro,nx,ny,kstrt,kyp)
! write real space diagnostic output: updates nferec
               call mpvwrite2(fmse,tdiag,nx,ny,kyp,iufe,nferec)
!              call mpvread2(fmse,tdiag,nx,ny,kyp,iufe,nferec,irc)
            endif
! calculate ion fluid moments
            if (movion==1) then
               if ((ndfm==2).or.(ndfm==3)) then
                  call dtimer(dtime,itime,-1)
                  fmsi = 0.0
                  call dtimer(dtime,itime,1)
                  tdiag = tdiag + real(dtime)
                  call wmgprofx2(pparti,fxye,fmsi,kipic,noff,nyp,qbmi,dt&
     &,ci,tdiag,npro,nx,mx,my,mx1,relativity)
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
! calculates fluid quantities from fluid moments
                  call mpfluidqs22(fmsi,tdiag,npro,nx,ny,kstrt,kyp)
! write real space diagnostic output: updates nfirec
                  call mpvwrite2(fmsi,tdiag,nx,ny,kyp,iufi,nfirec)
!                 call mpvread2(fmsi,tdiag,nx,ny,kyp,iufi,nfirec,irc)
               endif
            endif
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
     &irc)
                  if (irc==1) exit; irc = 0
               endif
! electron energy distribution
               if ((nvft==2).or.(nvft==3)) then
                  call mperpdist2(ppart,kpic,fe,sfv,eci,wk,tdiag,ndim,  &
     &nmv)
! display electron energy distribution
                  call pdisplayfe1(fe,wk,' ELECTRON',kstrt,ntime,nmv,   &
     &irc)
                  if (irc==1) exit; irc = 0
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
     &irc)
                     if (irc==1) exit; irc = 0
                  endif
! ion energy distribution
                  if ((nvft==2).or.(nvft==3)) then
                     call mperpdist2(pparti,kipic,fei,sfv,eci,wk,tdiag, &
     &ndim,nmv)
                     wk = rmass*wk
! display ion energy distribution
                     ts = fei(nmv21+1,1)
                     fei(nmv21+1,1) = rmass*fei(nmv21+1,1)
                     call pdisplayfe1(fei,wk,' ION',kstrt,ntime,nmv,irc)
                     fei(nmv21+1,1) = ts
                     if (irc==1) exit; irc = 0
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
! determines list of tagged particles leaving this node: updates iholep
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
                     call dafwritefv2(fvmtp,fvtp,fetp,ws,tdiag,iut,ntrec)
                  endif
! display test particle velocity distributions
                  call pdisplayfv1(fvtp,fvmtp,' ELECTRON',kstrt,ntime,  &
     &nmv,1,irc)
                  if (irc==1) exit; irc = 0
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
      call wmppush2(ppart,fxye,kpic,ncl,ihole,noff,nyp,qbme,dt,ci,wke,  &
     &tpush,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
!
! reorder electrons by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
      if (irc==0) then
         call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov,  &
     &kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,plist,irc2)
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
     &,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,.false.,irc2)
! ppart overflow
         else if (irc2(1)==4) then
! restores electron coordinates from ppbuff: updates ppart, ncl
            call mprstor2(ppart,ppbuff,ncl,ihole,tsort)
! copy ordered electrons to linear array: updates part
            call mpcopyout2(part,ppart,kpic,it,irc)
            deallocate(ppart)
            nppmx0 = (1.0 + xtras)*irc2(2)
            allocate(ppart(idimp,nppmx0,mxyp1))
! copies unordered electrons to ordered array: updates ppart
            call mpcopyin2(part,ppart,kpic,irc)
            call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov&
     &,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,plist,irc2)
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
         call wmppush2(pparti,fxye,kipic,ncl,ihole,noff,nyp,qbmi,dt,ci, &
     &wki,tpush,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
         wki = wki*rmass
!
! reorder ions by tile with OpenMP and MPI
! updates: pparti, kipic, and irc and possibly ncl and ihole
         if (irc==0) then
            call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort,  &
     &tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,plist,irc2)
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
     &,tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,.false.,irc2)
! pparti overflow
            else if (irc2(1)==4) then
! restores ion coordinates from ppbuff: updates pparti, ncl
               call mprstor2(pparti,ppbuff,ncl,ihole,tsort)
! copy ordered ions to linear array: updates part
               call mpcopyout2(part,pparti,kipic,it,irc)
               deallocate(pparti)
               nppmx1 = (1.0 + xtras)*irc2(2)
               allocate(pparti(idimp,nppmx1,mxyp1))
! copies unordered ions to ordered array: updates pparti
               call mpcopyin2(part,pparti,kipic,irc)
               call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort&
     &,tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,plist,irc2)
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
! need to reverse time lag in leap-frog integration scheme
      if (treverse==1) then
         if (((ntime+1)==(nloop/2)).or.((ntime+1)==nloop)) then
            dt = -dt
            ws = 0.0
            call wmppush2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ws,tpush&
     &,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
            call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov&
     &,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,plist,irc2)
            if (movion==1) then
               call wmppush2zf(pparti,kipic,ncl,ihole,noff,nyp,dt,ci,ws,&
     &tpush,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
               call ompmove2(pparti,kipic,ncl,ihole,noff,nyp,xtras,tsort&
     &,tmov,kstrt,nvp,nx,ny,mx,my,npbmx,nbmaxp,mx1,plist,irc2)
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
            call mpdsum(wtot,tdiag)
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
! restart file
      if (ntr > 0) then
         it = n/ntr
         if (n==ntr*it) then
            call dtimer(dtime,itime,-1)
! write out basic restart file for electrostatic code
            call bwrite_restart2(part,ppart,pparti,qi,kpic,kipic,tdiag, &
     &kstrt,iur,iscr,n,ntime0,irc)
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
! reset graphs
      if ((ntw > 0).or.(ntt > 0).or.(ntv > 0)) then
         if (nplot > 0) call reset_pgraphs(kstrt,irc)
      endif
!
      if (kstrt==1) then
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
     &,irc)
                  if (irc==1) stop
                  if (nplot > 0) call reset_nplot(nplot,irc)
               endif
            endif
         endif
      endif
!
! energy diagnostic
      if (ntw > 0) then
         ts = t0 + dt*real(ntw)
         call pdisplayw1(wt,ts,dt*real(ntw),kstrt,itw,irc)
         if (kstrt==1) then
            s(3) = (s(4) - s(3))/wt(1,4)
            write (iuot,*) 'Energy Conservation = ', real(s(3))
            s(1) = s(1)/real(itw)
            write (iuot,*) 'Average Field Energy <WE> = ', real(s(1))
            s(2) = s(2)/real(itw)
            write (iuot,*) 'Average Electron Kinetic Energy <WKE> = ',  &
     &real(s(2))
            write (iuot,*) 'Ratio <WE>/<WKE>= ', real(s(1)/s(2))
         endif
      endif
!
! velocity diagnostic
      if (ntv > 0) then
         ts = t0 + dt*real(ntv)
         call pdisplayfvt1(fvtm,' ELECT',ts,dt*real(ntv),kstrt,itv,irc)
! ions
         if (movion==1) then
            call pdisplayfvt1(fvtmi,' ION',ts,dt*real(ntv),kstrt,itv,   &
     &irc)
         endif
      endif
!
      if (kstrt==1) then
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
