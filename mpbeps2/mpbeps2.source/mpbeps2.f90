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
      use omplib
      use ompplib2
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
!     logical :: plist = .true.
      logical :: plist = .false.
!
! declare scalars for standard code
      integer :: n
      integer :: nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy
      integer :: mx1, ntime, nloop, isign, ierr
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
      integer :: it, mtw
      integer :: itw
! default Fortran unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19
      integer :: iude = 10, iup = 11, iuel = 12
      integer :: iudi = 20
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
      double precision, dimension(4) :: wtot, work
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
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
      endif
!
! check if too many node
      if (nvp > ny) then
         if (kstrt==1) then
         write (*,*) 'Too many nodes requested: ny, nvp=', ny, nvp
         endif
         call PPEXIT()
         stop
      endif
!
! open graphics device
      call IPLTCOMM(nplot)
      if (kstrt==1) then
         irc = open_pgraphs(nplot)
! set palette to color wheel
         call STPALIT(2)
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
! allocate data for standard code
      allocate(part(idimp,maxnp))
      allocate(qe(nxe,nypmx),qi(nxe,nypmx),fxye(ndim,nxe,nypmx))
      allocate(qt(nye,kxp),fxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),mixup(nxhy),sct(nxyh))
      allocate(kpic(mxyp1),iholep(ntmax+1))
!
! prepare fft tables
      call mpfft2_init(mixup,sct,indx,indy)
! calculate form factors
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
            call wmpvdistr2(part,nps,npp,vtx,vty,vx0,vy0,ci,npx,npy,    &
     &kstrt,nvp,relativity,ierr)
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
            call wmpvdistr2(part,nps,npp,vtdx,vtdy,vdx,vdy,ci,npxb,npyb,&
     &kstrt,nvp,relativity,ierr)
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
      allocate(ncl(8,mxyp1),ihole(2,ntmaxp+1,mxyp1))
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
         allocate(kipic(mxyp1))
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
               call wmpvdistr2(part,nps,nppi,vtxi,vtyi,vxi0,vyi0,ci,npxi&
     &,npyi,kstrt,nvp,relativity,ierr)
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
               call wmpvdistr2(part,nps,nppi,vtdxi,vtdyi,vdxi,vdyi,ci,  &
     &npxbi,npybi,kstrt,nvp,relativity,ierr)
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
!
! copy ordered ion data for OpenMP
         call mpmovin2(part,pparti,kipic,nppi,noff,mx,my,mx1,irc)
!
! sanity check for ions
         call mpcheck2(pparti,kipic,noff,nyp,nx,mx,my,mx1,irc)
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
         write (iuot,*) 'MPI nodes nvp = ', nvp
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
         call writnml2(iudm)
         write (iuot,*) ' * * * q.e.d. * * *'
         close(unit=iudm)
         close(unit=iuot)
      endif
!
      call PPEXIT()
      end program
