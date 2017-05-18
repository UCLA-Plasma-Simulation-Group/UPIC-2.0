!-----------------------------------------------------------------------
! 2-1/2D Darwin MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpdbeps2
      use in2
      use modmpinit2
      use modmppush2
      use modmpbpush2
      use modmpcurd2
      use modmpdpush2
      use modmpfield2
      use mpdiag2
      use mppmod2
      use pgraf2
      use omplib
      use ompplib2
      use pgraf2
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
      integer :: n, k
      integer :: nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy
      integer :: mdim, mx1, ntime, nloop, isign, ierr
      real :: qbme, affp, q2m0, wpm, wpmax, wpmin, omt, ws
      real :: qbmi, vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
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
      integer :: it, mtw
      integer :: itw
! default Fortran unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19
      integer :: iude = 10, iup = 11, iuel = 12
      integer :: iua = 13, iuet = 14, iub = 15
      integer :: iudi = 20, iuji = 21
! dimensions for fourier data
      integer :: modesxpd, modesy2de, modesy2p, modesy2el, modesy2di
      integer :: modesy2a, modesy2et, modesy2b, modesy2ji
      real :: wef
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe/qi = electron/ion charge density with guard cells
      real, dimension(:,:), allocatable :: qe, qi
! cue/cui = electron/ion current density with guard cells
      real, dimension(:,:,:), allocatable :: cue, cui
! dcu/dcui = electron/ion acceleration density with guard cells
      real, dimension(:,:,:), allocatable :: dcu, dcui
! amu/amui = electron/ion momentum flux with guard cells
      real, dimension(:,:,:), allocatable :: amu, amui
! cus = smoothed transverse electric field with guard cells
      real, dimension(:,:,:), allocatable :: cus
! exyze = smoothed total electric field with guard cells
! fxyze = smoothed longitudinal electric field with guard cells
! bxyze = smoothed magnetic field with guard cells
      real, dimension(:,:,:), allocatable :: fxyze, exyze, bxyze
! qt = scalar charge density field array in fourier space
      complex, dimension(:,:), allocatable :: qt
! cut = vector current density in fourier space
! dcut = vector acceleration density in fourier space
! amut = tensor momentum flux in fourier space
      complex, dimension(:,:,:), allocatable :: cut, dcut, amut
! exyt = vector transverse electric field in fourier space
! fxyt = vector longitudinal electric field in fourier space
! bxyt = vector magnetic field in fourier space
      complex, dimension(:,:,:), allocatable :: exyt, fxyt, bxyt
! ffc, ffe = form factor arrays for poisson solvers
      complex, dimension(:,:), allocatable :: ffc, ffe
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
      double precision, dimension(7) :: wtot, work
!
! declare arrays for MPI code:
! edges(1:2) = lower:upper y boundaries of particle partition
      real, dimension(:), allocatable  :: edges
! iholep = location of hole left in linear particle arrays
      integer, dimension(:), allocatable :: iholep
!
! declare arrays for OpenMP code:
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
! curit = store selected fourier modes for ion current density
      complex, dimension(:,:,:), allocatable :: curit
! vpott = store selected fourier modes for vector potential
      complex, dimension(:,:,:), allocatable :: vpott
! ett = store selected fourier modes for transverse efield
      complex, dimension(:,:,:), allocatable :: ett
! bt = store selected fourier modes for magnetic field
      complex, dimension(:,:,:), allocatable :: bt
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0
      real :: tmov = 0.0, tfmov = 0.0, tdiag = 0.0
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
! override default input data
         emf = 2
! read namelist
         call readnml2(iuin)
! override input data
         idcode = 3
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
! mdim = dimension of amu array
      mdim = 2*ndim - 2
      qbme = qme
      affp = dble(nx)*dble(ny)/np
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
! check if too many processors
      if (nvp > ny) then
         if (kstrt==1) then
         write (*,*) 'Too many processors requested: ny, nvp=', ny, nvp
         endif
         call PPEXIT()
         stop
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
      allocate(qe(nxe,nypmx),qi(nxe,nypmx),fxyze(ndim,nxe,nypmx))
      allocate(exyze(ndim,nxe,nypmx),cue(ndim,nxe,nypmx))
      allocate(dcu(ndim,nxe,nypmx),cus(ndim,nxe,nypmx))
      allocate(amu(mdim,nxe,nypmx),bxyze(ndim,nxe,nypmx))
      allocate(qt(nye,kxp),fxyt(ndim,nye,kxp),exyt(ndim,nye,kxp))
      allocate(cut(ndim,nye,kxp),dcut(ndim,nye,kxp))
      allocate(amut(mdim,nye,kxp),bxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),ffe(nyh,kxp),mixup(nxhy),sct(nxyh))
      allocate(kpic(mxyp1),iholep(ntmax+1))
!
! prepare fft tables
      call mpfft2_init(mixup,sct,indx,indy)
! calculate form factor: ffc
      call mppois2_init(ffc,ax,ay,affp,nx,ny,kstrt)
!
! initialize electrons
      nps = 1
      npp = 0
! background electrons
      if (npxy > 0.0d0) then
! calculates initial electron co-ordinates with uniform density
!        call mpudistr2(part,edges,npp,npx,npy,nx,ny,kstrt,ipbc,ierr)
! calculates initial electron co-ordinates with various density profiles
         call mpfdistr2(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,   &
     &shiftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,ierr)
! initialize electron velocities or momenta
         if (ierr==0) then
            call wmpvdistr2h(part,nps,npp,vtx,vty,vtz,vx0,vy0,vz0,ci,   &
     &npx,npy,kstrt,nvp,relativity,ierr)
         endif
! check for background electron initialization error
         if (ierr /= 0) then
            call PPEXIT()
            stop
         endif
      endif
! beam electrons
      if (npxyb > 0.0d0) then
         nps = npp + 1
! calculates initial electron co-ordinates with uniform density
!        call mpudistr2(part,edges,npp,npxb,npyb,nx,ny,kstrt,ipbc,ierr)
! calculates initial electron co-ordinates with various density profiles
         call mpfdistr2(part,npp,ampdx,scaledx,shiftdx,ampdy,scaledy,   &
     &shiftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,ierr)
! initialize electron velocities or momenta
         if (ierr==0) then
            call wmpvdistr2h(part,nps,npp,vtdx,vtdy,vtdz,vdx,vdy,vdz,ci,&
     &npxb,npyb,kstrt,nvp,relativity,ierr)
         endif
! check for beam electron initialization error
         if (ierr /= 0) then
            call PPEXIT()
            stop
         endif
      endif
!
! check if any electrons are in the wrong node
      call mpfholes2(part,edges,npp,iholep)
! iholep overflow
      if (iholep(1) < 0) then
         ntmax = -iholep(1)
         ntmax = 1.5*ntmax
         deallocate(iholep)
         if (kstrt==1) then
            write (*,*) 'reallocating electron iholep: ntmax=', ntmax
         endif
         allocate(iholep(ntmax+1))
         call mpfholes2(part,edges,npp,iholep)
         if (iholep(1) < 0) then
            if (kstrt==1) write (*,*) 'iholep overflow: ntmax=', ntmax
            call PPEXIT()
            stop
         endif
      endif
! more electrons to correct node
      call ipmove2(part,edges,npp,iholep,ny,tmov,kstrt,nvp,ierr)
      if (ierr /= 0) then
         call PPEXIT()
         stop
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
      allocate(ncl(8,mxyp1))
      allocate(ihole(2,ntmaxp+1,mxyp1))
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
         allocate(kipic(mxyp1),cui(ndim,nxe,nypmx))
         allocate(dcui(ndim,nxe,nypmx),amui(mdim,nxe,nypmx))
         cui = 0.0
         nps = 1
         nppi = 0
! background ions
         if (npxyi > 0.0d0) then
! calculates initial ion co-ordinates with uniform density
!           call mpudistr2(part,edges,nppi,npxi,npyi,nx,ny,kstrt,ipbc,  &
!    &ierr)
! calculates initial ion co-ordinates with various density profiles
            call mpfdistr2(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,   &
     &scaledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,ierr)
! initialize ion velocities or momenta
            if (ierr==0) then
               call wmpvdistr2h(part,nps,nppi,vtxi,vtyi,vtzi,vxi0,vyi0, &
     &vzi0,ci,npxi,npyi,kstrt,nvp,relativity,ierr)
            endif
! check for background ion initialization error
            if (ierr /= 0) then
               call PPEXIT()
               stop
            endif
         endif
! beam ions
         if (npxybi > 0.0d0) then
            nps = nppi + 1
! calculates initial ion co-ordinates with uniform density
!           call mpudistr2(part,edges,nppi,npxbi,npybi,nx,ny,kstrt,ipbc,&
!    &ierr)
! calculates initial ion co-ordinates with various density profiles
            call mpfdistr2(part,nppi,ampdxi,scaledxi,shiftdxi,ampdyi,   &
     &scaledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,ierr)
! initialize ion velocities or momenta
            if (ierr==0) then
               call wmpvdistr2h(part,nps,nppi,vtdxi,vtdyi,vtdzi,vdxi,   &
     &vdyi,vdzi,ci,npxbi,npybi,kstrt,nvp,relativity,ierr)
            endif
! check for beam ion initialization error
            if (ierr /= 0) then
               call PPEXIT()
               stop
            endif
         endif
!
! check if any ions are in the wrong node
         call mpfholes2(part,edges,nppi,iholep)
! iholep overflow
      if (iholep(1) < 0) then
         ntmax = -iholep(1)
         ntmax = 1.5*ntmax
         deallocate(iholep)
         if (kstrt==1) then
            write (*,*) 'reallocating ion iholep: ntmax=', ntmax
         endif
         allocate(iholep(ntmax+1))
         call mpfholes2(part,edges,nppi,iholep)
         if (iholep(1) < 0) then
            if (kstrt==1) write (*,*) 'iholep overflow: ntmax=', ntmax
            call PPEXIT()
            stop
         endif
      endif
! more ions to correct node
         call ipmove2(part,edges,nppi,iholep,ny,tmov,kstrt,nvp,ierr)
         if (ierr /= 0) then
            call PPEXIT()
            stop
         endif
!
! find number of ions in each of mx, my tiles: updates kipic, nppmx
         call mpdblkp2(part,kipic,nppi,noff,nppmx,mx,my,mx1,irc)
!
! allocate vector ion data
         nppmx1 = (1.0 + xtras)*nppmx
         allocate(pparti(idimp,nppmx1,mxyp1))
! copy ordered ion data for OpenMP
         call mpmovin2(part,pparti,kipic,nppi,noff,mx,my,mx1,irc)
!
! sanity check for ions
         call mpcheck2(pparti,kipic,noff,nyp,nx,mx,my,mx1,irc)
      endif
!
! find maximum and minimum initial electron density
      qe = 0.0
      call mppost2(ppart,qe,kpic,noff,qme,time,mx,my,mx1)
      call wmpaguard2(qe,nyp,tguard,nx,kstrt,nvp)
      if (movion==1) then
         qi = 0.0
         call mppost2(pparti,qi,kipic,noff,qmi,time,mx,my,mx1)
         call wmpaguard2(qi,nyp,tguard,nx,kstrt,nvp)
         call mpfwptminx2(qe,qi,nyp,qbme,qbmi,wpmax,wpmin,nx)
      else
         call mpfwpminx2(qe,nyp,qbme,wpmax,wpmin,nx)
      endif
      wtot(1) = wpmax
      wtot(2) = -wpmin
      call PPDMAX(wtot,work,2)
      wpmax = wtot(1)
      wpmin = -wtot(2)
      wpm = 0.5*(wpmax + wpmin)*affp
! accelerate convergence: update wpm
      if (wpm <= 10.0) wpm = 0.75*wpm
      if (kstrt==1) write (iuot,*) 'wpm=',wpm
      q2m0 = wpm/affp
!
! calculate form factor: ffe
      call mpepois2_init(ffe,ax,ay,affp,wpm,ci,nx,ny,kstrt)
!
! initialize electric fields
      cus = 0.0; fxyze = 0.0
! set magnitude of external transverse magnetic field
      omt = sqrt(omx*omx + omy*omy + omz*omz)
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
      if ((ntel>0).or.(nta>0).or.(ntet>0).or.(ntb>0).or.(ntji>0)) then
         allocate(vfieldc(ndim,nye,kxp),vfield(ndim,nxe,nypmx))
      endif
!
! initialize electron density diagnostic
      if (ntde > 0) then
         modesxde = min(modesxde,nxh+1)
         modesyde = min(modesyde,nyh+1)
         modesy2de = min(2*modesyde-1,ny)
         modesxpd = min(modesxde,kxp)
         if (modesxde==(nxh+1)) modesxpd = modesxpd + 1
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
            modesxdi = min(modesxdi,nxh+1)
            modesydi = min(modesydi,nyh+1)
            modesy2di = min(2*modesydi-1,ny)
            modesxpd = min(modesxdi,kxp)
            if (modesxdi==(nxh+1)) modesxpd = modesxpd + 1
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
         modesxp = min(modesxp,nxh+1)
         modesyp = min(modesyp,nyh+1)
         modesy2p = min(2*modesyp-1,ny)
         modesxpd = min(modesxp,kxp)
         if (modesxp==(nxh+1)) modesxpd = modesxpd + 1
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
         modesxel = min(modesxel,nxh+1)
         modesyel = min(modesyel,nyh+1)
         modesy2el = min(2*modesyel-1,ny)
         modesxpd = min(modesxel,kxp)
         if (modesxel==(nxh+1)) modesxpd = modesxpd + 1
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
! initialize ion current density diagnostic
      if (movion==1) then
         if (ntji > 0) then
            modesxji = min(modesxji,nxh+1)
            modesyji = min(modesyji,nyh+1)
            modesy2ji = min(2*modesyji-1,ny)
            modesxpd = min(modesxji,kxp)
            if (modesxji==(nxh+1)) modesxpd = modesxpd + 1
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
! initialize vector potential diagnostic
      if (nta > 0) then
         modesxa = min(modesxa,nxh+1)
         modesya = min(modesya,nyh+1)
         modesy2a = min(2*modesya-1,ny)
         modesxpd = min(modesxa,kxp)
         if (modesxa==(nxh+1)) modesxpd = modesxpd + 1
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
         modesxet = min(modesxet,nxh+1)
         modesyet = min(modesyet,nyh+1)
         modesy2et = min(2*modesyet-1,ny)
         modesxpd = min(modesxet,kxp)
         if (modesxet==(nxh+1)) modesxpd = modesxpd + 1
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
         modesxb = min(modesxb,nxh+1)
         modesyb = min(modesyb,nyh+1)
         modesy2b = min(2*modesyb-1,ny)
         modesxpd = min(modesxb,kxp)
         if (modesxb==(nxh+1)) modesxpd = modesxpd + 1
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
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (kstrt==1) write (iuot,*) 'program mpdbeps2'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
      if (kstrt==1) write (iuot,*) 'ntime = ', ntime
!
! deposit electron current with OpenMP: updates cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      call wmpdjpost2(ppart,cue,kpic,ncl,ihole,noff,nyp,qme,zero,ci,    &
     &tdjpost,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
! add guard cells with OpenMP: updates cue
      call wmpnacguard2(cue,nyp,tguard,nx,kstrt,nvp)
!
! deposit ion current with OpenMP: updates cui
      if (movion==1) then
         call dtimer(dtime,itime,-1)
         cui = 0.0
         call dtimer(dtime,itime,1)
         tdjpost = tdjpost + real(dtime)
         call wmpdjpost2(pparti,cui,kipic,ncl,ihole,noff,nyp,qmi,zero,ci&
     &,tdjpost,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
! add guard cells with OpenMP: updates cui
         call wmpnacguard2(cui,nyp,tguard,nx,kstrt,nvp)
      endif
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
               call PPEXIT()
               stop
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
                  call PPEXIT()
                  stop
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
! calculate longitudinal force/charge in fourier space with OpenMP:
! updates fxyt, we
      call mppois2(qt,fxyt,ffc,we,tfield,nx,ny,kstrt)
!
! transform force to real space with OpenMP:
! moves data to non-uniform partition
! updates fxyze, nterf, and ierr, modifies fxyt
      isign = 1
      call wmpfft2rn(fxyze,fxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
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
! calculate magnetic field in fourier space with OpenMP:
! updates bxyt, wm
      call mpbbpois2(cut,bxyt,ffc,ci,wb,tfield,nx,ny,kstrt)
!
! transform magnetic field to real space with OpenMP:
! moves data to non-uniform partition
! updates bxyze, nterf, and ierr, modifies bxyt
      isign = 1
      call wmpfft2rn(bxyze,bxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! add constant to magnetic field with OpenMP: updates bxyze
      call mpbaddext2(bxyze,nyp,tfield,omx,omy,omz,nx)
!
! copy guard cells with OpenMP: updates fxyze, bxyze
      call wmpncguard2(fxyze,nyp,tguard,nx,kstrt,nvp)
      call wmpncguard2(bxyze,nyp,tguard,nx,kstrt,nvp)
!
! add longitudinal and old transverse electric fields with OpenMP:
! updates exyze
      call mpaddvrfield2(exyze,cus,fxyze,tfield)
!
! deposit electron acceleration density and momentum flux with OpenMP:
! updates dcu, amu
      call dtimer(dtime,itime,-1)
      dcu = 0.0; amu = 0.0
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      call wmpgdjpost2(ppart,exyze,bxyze,dcu,amu,kpic,noff,nyp,qme,qbme,&
     &dt,ci,tdcjpost,nx,mx,my,mx1,relativity)
! add guard cells with OpenMP: updates dcu, amu
      call wmpnacguard2(dcu,nyp,tguard,nx,kstrt,nvp)
      call wmpnacguard2(amu,nyp,tguard,nx,kstrt,nvp)
! deposit ion acceleration density and momentum flux with OpenMP:
! updates dcui, amui
      if (movion==1) then
         call dtimer(dtime,itime,-1)
         dcui = 0.0; amui = 0.0
         call dtimer(dtime,itime,1)
         tdcjpost = tdcjpost + real(dtime)
         call wmpgdjpost2(pparti,exyze,bxyze,dcui,amui,kipic,noff,nyp,  &
     &qmi,qbmi,dt,ci,tdcjpost,nx,mx,my,mx1,relativity)
! add guard cells with OpenMP: updates dcui, amui
         call wmpnacguard2(dcui,nyp,tguard,nx,kstrt,nvp)
         call wmpnacguard2(amui,nyp,tguard,nx,kstrt,nvp)
! add electron and ion momentum flux densities: updates dcu, amu
         call mpaddcuei2(dcu,dcui,nyp,tfield,nx)
         call mpaddamui2(amu,amui,nyp,tfield,nx)
      endif
!
! add old scaled electric field with OpenMP: updates dcu
      call mpascfguard2(dcu,cus,nyp,q2m0,tdcjpost,nx)
!
! transform acceleration density and momentum flux to fourier space
! with OpenMP and moves data to uniform partition
! updates dcut, amut, nterf, and ierr, modifies dcu, amu
      isign = -1
      call wmpfft2rn(dcu,dcut,noff,nyp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
      call wmpfft2rn(amu,amut,noff,nyp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! take transverse part of time derivative of current with OpenMP:
! updates dcut
      call mpadcuperp2(dcut,amut,tfield,nx,ny,kstrt)
!
! calculate transverse electric field with OpenMP: updates exyt, wf
      call mpepois2(dcut,exyt,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
!
! transform transverse electric field to real space with OpenMP:
! moves data to non-uniform partition
! updates cus, nterf, and ierr, modifies cur
      isign = 1
      call wmpfft2rn(cus,exyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! copy guard cells with OpenMP: updates cus
      call wmpncguard2(cus,nyp,tguard,nx,kstrt,nvp)
!
! add longitudinal and transverse electric fields with OpenMP:
! exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call mpaddvrfield2(exyze,cus,fxyze,tfield)
!
! inner iteration loop
      do k = 1, ndc
!
! deposit electron current and acceleration density and momentum flux
! with OpenMP: updates cue, dcu, amu
      call dtimer(dtime,itime,-1)
      cue = 0.0; dcu = 0.0; amu = 0.0
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      call wmpgdcjpost2(ppart,exyze,bxyze,cue,dcu,amu,kpic,noff,nyp,qme,&
     &qbme,dt,ci,tdcjpost,nx,mx,my,mx1,relativity)
! add guard cells for current, acceleration density, and momentum flux
! with OpenMP: updates cue, dcu, amu
      call wmpnacguard2(cue,nyp,tguard,nx,kstrt,nvp)
      call wmpnacguard2(dcu,nyp,tguard,nx,kstrt,nvp)
      call wmpnacguard2(amu,nyp,tguard,nx,kstrt,nvp)
!
! deposit ion current and acceleration density and momentum flux
! with OpenMP: updates cui, dcui, amui
      if (movion==1) then
         call dtimer(dtime,itime,-1)
         cui = 0.0; dcui = 0.0; amui = 0.0
         call dtimer(dtime,itime,1)
         tdcjpost = tdcjpost + real(dtime)
         call wmpgdcjpost2(pparti,exyze,bxyze,cui,dcui,amui,kipic,noff, &
     &nyp,qmi,qbmi,dt,ci,tdcjpost,nx,mx,my,mx1,relativity)
! add guard cells for current, acceleration density, and momentum flux
! with OpenMP: updates cui, dcui, amui
         call wmpnacguard2(cui,nyp,tguard,nx,kstrt,nvp)
         call wmpnacguard2(dcui,nyp,tguard,nx,kstrt,nvp)
         call wmpnacguard2(amui,nyp,tguard,nx,kstrt,nvp)
! add electron and ion densities: updates cue, dcu, amu
         call mpaddcuei2(cue,cui,nyp,tfield,nx)
         call mpaddcuei2(dcu,dcui,nyp,tfield,nx)
         call mpaddamui2(amu,amui,nyp,tfield,nx)
      endif
!
! add scaled electric field with OpenMP: updates dcu
      call mpascfguard2(dcu,cus,nyp,q2m0,tdcjpost,nx)
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
! calculate magnetic field in fourier space with OpenMP:
! updates bxyt, wm
      call mpbbpois2(cut,bxyt,ffc,ci,wb,tfield,nx,ny,kstrt)
!
! transform magnetic field to real space with OpenMP:
! moves data to non-uniform partition
! updates bxyze, nterf, and ierr, modifies bxyt
      isign = 1
      call wmpfft2rn(bxyze,bxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! add constant to magnetic field with OpenMP: updates bxyze
      call mpbaddext2(bxyze,nyp,tfield,omx,omy,omz,nx)
!
! transform acceleration density and momentum flux to fourier space
! with OpenMP and moves data to uniform partition
! updates dcut, amut, nterf, and ierr, modifies dcu, amu
      isign = -1
      call wmpfft2rn(dcu,dcut,noff,nyp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
      call wmpfft2rn(amu,amut,noff,nyp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! take transverse part of time derivative of current with OpenMP:
! updates dcut
      call mpadcuperp2(dcut,amut,tfield,nx,ny,kstrt)
!
! calculate transverse electric field with OpenMP: updates exyt, wf
      call mpepois2(dcut,exyt,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
!
! transform transverse electric field to real space with OpenMP:
! moves data to non-uniform partition
! updates cus, nterf, and ierr, modifies cur
      isign = 1
      call wmpfft2rn(cus,exyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! copy guard cells with OpenMP: updates bxyze, cus
      call wmpncguard2(bxyze,nyp,tguard,nx,kstrt,nvp)
      call wmpncguard2(cus,nyp,tguard,nx,kstrt,nvp)
!
! add longitudinal and transverse electric fields with OpenMP:
! exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call mpaddvrfield2(exyze,cus,fxyze,tfield)
!
      enddo
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
               call PPEXIT()
               stop
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
               call PPEXIT()
               stop
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
            call mpapot2(cut,vfieldc,ffc,ci,ws,tfield,nx,ny,kstrt)
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
               call PPEXIT()
               stop
            endif
            ierr = 0
         endif
      endif
!
! transverse efield diagnostic
      if (ntet > 0) then
         it = ntime/ntet
         if (ntime==ntet*it) then
! calculate transverse efield in fourier space: updates vfieldc
            call mpetfield2(dcut,vfieldc,ffe,affp,ci,ws,tfield,nx,ny,   &
     &kstrt)
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
               call PPEXIT()
               stop
            endif
            ierr = 0
         endif
      endif
!
! magnetic field diagnostic
      if (ntb > 0) then
         it = ntime/ntb
         if (ntime==ntb*it) then
! calculate magnetic field in fourier space: updates vfieldc
            call mpibpois2(cut,vfieldc,ffc,ci,ws,tfield,nx,ny,kstrt)
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
               call PPEXIT()
               stop
            endif
            ierr = 0
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
                  call PPEXIT()
                  stop
               endif
               ierr = 0
            endif
         endif
      endif
!
! push electrons with OpenMP:
! updates ppart and wke, and possibly ncl, ihole, irc
      wke = 0.0
      call wmpbpush2(ppart,exyze,bxyze,kpic,ncl,ihole,noff,nyp,qbme,dt, &
     &dt,ci,wke,tpush,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
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
         call wmpbpush2(pparti,exyze,bxyze,kipic,ncl,ihole,noff,nyp,qbmi&
     &,dt,dt,ci,wki,tpush,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
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
            wef = we + wb
            wtot(1) = wef
            wtot(2) = wke
            wtot(3) = wki
            wtot(4) = wef + wke
            wtot(5) = we
            wtot(6) = wf
            wtot(7) = wb
            call PPDSUM(wtot,work,7)
            wke = wtot(2)
            wki = wtot(3)
            we = wtot(5)
            wf = wtot(6)
            wb = wtot(7)
            wef = we + wb
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
      enddo
      ntime = ntime + 1
!
! loop time
      call dtimer(dtime,ltime,1)
      tloop = tloop + real(dtime)
!
! * * * end main iteration loop * * *
!
      if (kstrt.eq.1) then
         write (iuot,*)
         write (iuot,*) 'ntime, relativity, ndc = ', ntime, relativity, &
     &ndc
         if (treverse==1) write (iuot,*) 'treverse = ', treverse
         write (iuot,*) 'MPI nodes nvp = ', nvp
!
! energy diagnostic
         if (ntw > 0) then
!           ts = t0 + dt*real(ntw)
!           call displayw1(wt,ts,dt*real(ntw),itw,irc)
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
         write (iuot,*) 'current derivative deposit time = ', tdcjpost
         tdpost = tdpost + tdjpost + tdcjpost
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
! ion diagnostics
         if (movion==1) then
! ion density diagnostic
            if (ntdi > 0) then
               ndirec = ndirec - 1
            endif
! ion current diagnostic
            if (ntdi > 0) then
               njirec = njirec - 1
            endif
         endif
! write final diagnostic metafile
         call writnml2(iudm)
         write (iuot,*) ' * * * q.e.d. * * *'
         close(unit=iudm)
         close(unit=iuot)
      endif
!
      call PPEXIT()
      end program
