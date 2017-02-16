!-----------------------------------------------------------------------
! 3D Electrostatic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpbeps3
      use modmpinit3
      use modmppush3
      use modmpfield3
      use mppmod3
      use omplib
      use ompplib3
      implicit none
! indx/indy/indz = exponent which determines grid points in x/y/z
! direction: nx = 2**indx, ny = 2**indy, nz = 2**indz.
      integer :: indx =   7, indy =   7, indz =   7
! npx/npy/npz = number of background electrons distributed in x/y/z
! direction
      integer :: npx = 384, npy = 384, npz = 384
! npxb/npyb/npzb = number of beam electrons distributed in x/y/z
! direction
      integer :: npxb = 0, npyb = 0, npzb = 0
! ndim = number of velocity coordinates = 3
      integer :: ndim = 3
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
! qme = charge on electron, in units of e.
      real :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx/vty/vtz = thermal velocity of background electrons in x/y/z
! direction
      real :: vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of background electrons in x/y/z
! direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
! ax/ay/az = smoothed particle size in x/y/z direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, az = .912871, ci = 0.1
! idimp = number of particle coordinates = 6
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 6, ipbc = 1, relativity = 0
! idps = number of partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      integer :: idps = 4, idds =    2
! wke/wki/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, wki = 0.0, we = 0.0, wt = 0.0
! mx/my/mz = number of grids in x/y/z in sorting tiles
! sorting tiles, should be less than or equal to 16
      integer :: mx = 8, my = 8, mz = 8
! fraction of extra particles needed for particle management
      real :: xtras = 0.2
! plist = (true,false) = list of particles leaving tiles found in push
      logical :: plist = .true.
!
! movion = (0,1) = (no,yes) move the ions
      integer :: movion = 1
! npxi/npyi/npzi = number of background ions distributed in x/y/z
! direction
      integer :: npxi = 384, npyi = 384, npzi = 384
! npxbi/npybi/npzbi = number of beam ions distributed in x/y/z direction
      integer :: npxbi = 0, npybi = 0, npzbi = 0
! qmi = charge on ion, in units of e
! rmass = ion/electron mass ratio
      real :: qmi = 1.0, rmass = 100.0
! rtempxi/rtempyi/rtempzi = electron/ion temperature ratio of background
! ions in x/y/z direction
      real :: rtempxi = 1.0, rtempyi = 1.0, rtempzi = 1.0
! vxi0/vyi0/vzi0 = drift velocity of ions in x/y/z direction
      real :: vxi0 = 0.0, vyi0 = 0.0, vzi0 = 0.0
! rtempdxi/rtempdyi/rtempdzi = electron/ion temperature ratio of beam
! ions in x/y/z direction
      real :: rtempdxi = 1.0, rtempdyi = 1.0, rtempdzi = 1.0
! vdxi/vdyi/vdzi = drift velocity of beam ions in x/y/z direction
      real :: vdxi = 0.0, vdyi = 0.0, vdzi = 0.0
!
! monitor = (0,1,2) = (disable,normal,extended) error processing
      integer :: monitor = 0
!
! declare scalars for standard code
      integer :: n
      integer :: nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh, nnxe
      integer :: nxyzh, nxhyz, mx1, ntime, nloop, isign, it, ierr
      real :: qbme, affp
      real :: qbmi, vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      double precision :: npxyz, npxyzb, np, npxyzi, npxyzbi
      double precision :: npi = 0.0d0
!
! declare scalars for MPI code
      integer :: nvpy, nvpz, nvp, idproc, kstrt, npmax, kyp, kzp
      integer :: kxyp, kyzp, nypmx, nzpmx, nypmn, nzpmn
      integer :: npp, nps, myp1, mzp1, mxyzp1, mxzyp1
      integer :: npimax, nppi
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nppmx1, nbmaxp, ntmaxp, npbmx, nvpp
      integer :: irc = 0
      integer, dimension(2) :: irc2 = 0
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
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0, tmov = 0.0, tfmov = 0.0
      real, dimension(2) :: tfft = 0.0
      double precision :: dtime
!
! start timing initialization
      call dtimer(dtime,itime,-1)
!
      irc = 0
! nvpp = number of shared memory nodes (0=default)
      nvpp = 0
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvpp
! initialize for shared memory parallel processing
      call INIT_OMP(nvpp)
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
      nxeh = nxe/2; nnxe = ndim*nxe
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
! nvp = number of MPI ranks
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
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
      call mpdcomp3(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,&
     &nvpy,nvpz)
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
! myp1/mzp1 = number of tiles in y/z direction
      myp1 = (nyzp(1) - 1)/my + 1; mzp1 = (nyzp(2) - 1)/mz + 1
! mxzyp1 = mx1*max(max(mzp1),max(myp1))
      mxzyp1 = mx1*max((nzpmx-2)/mz+1,(nypmx-2)/my+1)
      mxyzp1 = mx1*myp1*mzp1
! mterf = number of shifts required by field manager in y/z (0=search)
      mterf = 0
!
! allocate data for standard code
      allocate(part(idimp,max(npmax,npimax)))
      allocate(qe(nxe,nypmx,nzpmx),qi(nxe,nypmx,nzpmx))
      allocate(fxyze(ndim,nxe,nypmx,nzpmx))
      allocate(qt(nze,kxyp,kyzp))
      allocate(fxyzt(ndim,nze,kxyp,kyzp))
      allocate(ffc(nzh,kxyp,kyzp),mixup(nxhyz),sct(nxyzh))
      allocate(kpic(mxyzp1))
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
         call mpdistr3(part,edges,npp,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,  &
        &npz,nx,ny,nz,kstrt,ipbc,ierr)
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
         call mpdistr3(part,edges,npp,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,  &
     &npyb,npzb,nx,ny,nz,kstrt,ipbc,ierr)
! check for beam electron initialization error
         if (ierr /= 0) then
            call PPEXIT()
            stop
         endif
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
            call mpdistr3(part,edges,nppi,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,&
     &npxi,npyi,npzi,nx,ny,nz,kstrt,ipbc,ierr)
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
            call mpdistr3(part,edges,nppi,vtdxi,vtdyi,vtdzi,vdxi,vdyi,  &
     &vdzi,npxbi,npybi,npzbi,nx,ny,nz,kstrt,ipbc,ierr)
! check for background ion initialization error
            if (ierr /= 0) then
               call PPEXIT()
               stop
            endif
         endif
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
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (kstrt==1) write (*,*) 'program mpbeps3'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
!     if (kstrt==1) write (*,*) 'ntime = ', ntime
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
! add electron and ion densities: updates qe
      call mpaddqei3(qe,qi,nyzp,tfield,nx)
!
! transform charge to fourier space with OpenMP:
! updates qt, mterf, and ierr, modifies qe
      isign = -1
      call wmpfft3r(qe,qt,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx,indy&
     &,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! calculate force/charge in fourier space with OpenMP: updates fxyzt, we
      call mppois3(qt,fxyzt,ffc,we,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! transform force to real space with OpenMP:
! updates fxyze, mterf, and ierr, modifies fxyzt
      isign = 1
      call wmpfft3rn(fxyze,fxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! copy guard cells with OpenMP: updates fxyze
      call wmpncguard3(fxyze,nyzp,tguard,nx,kstrt,nvpy,nvpz)
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
! energy diagnostic
      wtot(1) = we
      wtot(2) = wke
      wtot(3) = wki
      wtot(4) = we + wke
      call PPDSUM(wtot,work,4)
      we = wtot(1)
      wke = wtot(2)
      wki = wtot(3)
      if (ntime==0) then
         if (kstrt==1) then
            write (*,*) 'Initial Field, Kinetic and Total Energies:'
            if (movion==0) then
               write (*,'(3e14.7)') we, wke, wke + we
            else
               write (*,'(4e14.7)') we, wke, wki, wke + wki + we
            endif
         endif
      endif
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
         write (*,*)
         write (*,*) 'ntime, relativity = ', ntime, relativity
         write (*,*) 'MPI nodes nvpy, nvpz = ', nvpy, nvpz
         write (*,*) 'Final Field, Kinetic and Total Energies:'
         if (movion==0) then
            write (*,'(3e14.7)') we, wke, wke + we
         else
            write (*,'(4e14.7)') we, wke, wki, wke + wki + we
         endif
!
         write (*,*)
         write (*,*) 'initialization time = ', tinit
         write (*,*) 'deposit time = ', tdpost
         write (*,*) 'guard time = ', tguard
         write (*,*) 'solver time = ', tfield
         write (*,*) 'field move time = ', tfmov
         write (*,*) 'fft and transpose time = ', tfft(1), tfft(2)
         write (*,*) 'push time = ', tpush
         write (*,*) 'particle move time = ', tmov
         write (*,*) 'sort time = ', tsort
         tfield = tfield + tguard + tfft(1) + tfmov
         write (*,*) 'total solver time = ', tfield
         tsort = tsort + tmov
         time = tdpost + tpush + tsort
         write (*,*) 'total particle time = ', time
         wt = time + tfield
         tloop = tloop - wt
         write (*,*) 'total and additional time = ', wt, tloop
         write (*,*)
!
         wt = 1.0e+09/(real(nloop)*real(np+npi))
         write (*,*) 'Push Time (nsec) = ', tpush*wt
         write (*,*) 'Deposit Time (nsec) = ', tdpost*wt
         write (*,*) 'Sort Time (nsec) = ', tsort*wt
         write (*,*) 'Total Particle Time (nsec) = ', time*wt
      endif
!
      call PPEXIT()
      end program
