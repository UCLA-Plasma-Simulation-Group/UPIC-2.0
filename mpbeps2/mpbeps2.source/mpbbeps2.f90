!-----------------------------------------------------------------------
! 2-1/2D Electromagnetic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpbbeps2
      use modmpinit2
      use modmppush2
      use modmpbpush2
      use modmpcurd2
      use modmpfield2
      use mppmod2
      use omplib
      use ompplib2
      implicit none
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer :: indx =   9, indy =   9
! npx/npy = number of background electrons distributed in x/y direction
! npxb/npyb = number of beam electrons distributed in x/y direction
      integer :: npx =  3072, npy =   3072, npxb = 0, npyb =  0
! ndim = number of velocity coordinates = 3
      integer :: ndim = 3
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations
! qme = charge on electron, in units of e.
      real :: tend = 10.0, dt = 0.04, qme = -1.0
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
! ax/ay = smoothed particle size in x/y direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, ci = 0.1
! idimp = dimension of phase space = 5
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 5, ipbc = 1, relativity = 1
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! idps = number of partition boundaries
      integer :: idps = 2
! wke/wki/we = particle kinetic/electrostatic field energy
! wf/wm = transverse electric field/magnetic field
      real :: wke = 0.0, wki = 0.0, we = 0.0, wf = 0.0, wm = 0.0
! wt = total energy
      real :: wt = 0.0
! sorting tiles, should be less than or equal to 32
      integer :: mx = 16, my = 16
! fraction of extra particles needed for particle management
      real :: xtras = 0.2
! plist = (true,false) = list of particles leaving tiles found in push
      logical :: plist = .true.
!
! movion = (0,1) = (no,yes) move the ions
      integer :: movion = 1
! npxi/npyi = number of background ions distributed in x/y direction
! npxib/npyib = number of beam ions distributed in x/y direction
      integer :: npxi = 3072, npyi = 3072, npxbi = 0, npybi = 0
! qmi = charge on ion, in units of e
! rmass = ion/electron mass ratio
      real :: qmi = 1.0, rmass = 100.0
! rtempxi/rtempyi/rtempzi = electron/ion temperature ratio of background
! ions in x/y/z direction
      real :: rtempxi = 1.0, rtempyi = 1.0, rtempzi = 1.0
! vxi0/vyi0/vzi0 = drift velocity of background ions in x/y/z direction
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
      integer :: nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy
      integer :: mx1, ntime, nloop, isign, it, ierr
      real :: qbme, affp, dth, omt
      real :: qbmi, vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      double precision :: npxy, npxyb, np, npxyi, npxybi
      double precision :: npi = 0.0d0
!
! declare scalars for MPI code
      integer :: nvp, idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn
      integer :: nyp, noff, npp, nps, myp1, mxyp1
      integer :: npimax, nppi
      integer :: nterf
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nppmx1, nbmaxp, ntmaxp, npbmx, nvpp
      integer :: irc = 0
      integer, dimension(2) :: irc2 = 0
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
      double precision, dimension(7) :: wtot, work
! arrays required for darwin initial fields
      real, dimension(:,:,:), allocatable :: amu
      complex, dimension(:,:,:), allocatable :: dcut, amut
!
! declare arrays for MPI code
! edges(1:2) = lower:upper y boundaries of particle partition
      real, dimension(:), allocatable  :: edges
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
      real :: tdjpost = 0.0, tpush = 0.0, tsort = 0.0, tmov = 0.0
      real :: tfmov = 0.0
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
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
         vtdzi = vtz/sqrt(rmass*rtempdzi)
      endif
      dth = 0.0
!      
! nvp = number of distributed memory nodes
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
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
      call mpdcomp2(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp)
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
! myp1 = number of tiles in y direction
      myp1 = (nyp - 1)/my + 1; mxyp1 = mx1*myp1
! nterf = number of shifts required by field manager (0=search)
      nterf = 0
!
! allocate and initialize data for standard code
      allocate(part(idimp,max(npmax,npimax)))
      allocate(qe(nxe,nypmx),qi(nxe,nypmx))
      allocate(cue(ndim,nxe,nypmx))
      allocate(fxyze(ndim,nxe,nypmx),bxyze(ndim,nxe,nypmx))
      allocate(exyz(ndim,nye,kxp),bxyz(ndim,nye,kxp))
      allocate(qt(nye,kxp),fxyt(ndim,nye,kxp))
      allocate(cut(ndim,nye,kxp),bxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),mixup(nxhy),sct(nxyh))
      allocate(kpic(mxyp1))
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
         call mudistr2(part,edges,npp,npx,npy,nx,ny,kstrt,ipbc,irc)
! initialize electron velocities or momenta
         call wmvdistr2h(part,nps,npp,vtx,vty,vtz,vx0,vy0,vz0,ci,npx,npy&
     &,kstrt,relativity,irc)
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
         call mudistr2(part,edges,npp,npxb,npyb,nx,ny,kstrt,ipbc,irc)
! initialize electron velocities or momenta
         call wmvdistr2h(part,nps,npp,vtdx,vtdy,vtdz,vdx,vdy,vdz,ci,npxb&
     &,npyb,kstrt,relativity,irc)
! check for beam electron initialization error
         if (ierr /= 0) then
            call PPEXIT()
            stop
         endif
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
         cui = 0.0
         nps = 1
         nppi = 0
! background ions
         if (npxyi > 0.0d0) then
! calculates initial ion co-ordinates with uniform density
            call mudistr2(part,edges,nppi,npxi,npyi,nx,ny,kstrt,ipbc,irc&
     &)
! initialize ion velocities or momenta
            call wmvdistr2h(part,nps,nppi,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,&
     &ci,npxi,npyi,kstrt,relativity,irc)
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
            call mudistr2(part,edges,nppi,npxbi,npybi,nx,ny,kstrt,ipbc, &
     &irc)
! initialize ion velocities or momenta
            call wmvdistr2h(part,nps,nppi,vtdxi,vtdyi,vtdzi,vdxi,vdyi,  &
     &vdzi,ci,npxbi,npybi,kstrt,relativity,irc)
! check for beam ion initialization error
            if (ierr /= 0) then
               call PPEXIT()
               stop
            endif
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
! initialize transverse electromagnetic fields
      exyz = cmplx(0.0,0.0)
      bxyz = cmplx(0.0,0.0)
! set magnitude of external magnetic field
      omt = sqrt(omx*omx + omy*omy + omz*omz)
!
      if (dt > 0.45*ci) then
         if (kstrt==1) then
            write (*,*) 'Warning: Courant condition may be exceeded!'
         endif
      endif
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (kstrt==1) write (*,*) 'program mpbbeps2'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
!     if (kstrt==1) write (*,*) 'ntime = ', ntime
!
! deposit electron current with OpenMP:
! updates ppart and cue, and possibly ncl, ihole, irc
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      call wmpdjpost2(ppart,cue,kpic,ncl,ihole,noff,nyp,qme,dth,ci,     &
     &tdjpost,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
! add guard cells with OpenMP: updates cue
      call wmpnacguard2(cue,nyp,tguard,nx,kstrt,nvp)
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
! deposit ion current with OpenMP:
      if (movion==1) then
! updates pparti and cui, and possibly ncl, ihole, irc
         call dtimer(dtime,itime,-1)
         cui = 0.0
         call dtimer(dtime,itime,1)
         tdjpost = tdjpost + real(dtime)
         call wmpdjpost2(pparti,cui,kipic,ncl,ihole,noff,nyp,qmi,dth,ci,&
     &tdjpost,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
! add guard cells with OpenMP: updates cui
         call wmpnacguard2(cui,nyp,tguard,nx,kstrt,nvp)
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
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      call mppost2(ppart,qe,kpic,noff,qme,tdpost,mx,my,mx1)
! add guard cells with OpenMP: updates qe
      call wmpaguard2(qe,nyp,tguard,nx,kstrt,nvp)
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
! add electron and ion densities: updates qe
      call mpaddqei2(qe,qi,nyp,tfield,nx)
!
! add electron and ion current densities: updates cue
      if (movion==1) call mpaddcuei2(cue,cui,nyp,tfield,nx)
!
! transform charge to fourier space with OpenMP:
! updates qt, nterf, and ierr, modifies qe
      isign = -1
      call wmpfft2r(qe,qt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,indy,&
     &kstrt,nvp,kyp,ny,nterf,ierr)
!
! transform current to fourier space with OpenMP:
! updates cut, nterf, and ierr, modifies cue
      isign = -1
      call wmpfft2rn(cue,cut,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,  &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! take transverse part of current with OpenMP: updates cut
      call mpcuperp2(cut,tfield,nx,ny,kstrt)
!
! calculate electromagnetic fields in fourier space with OpenMP:
! updates exyz, bxyz, wf, wm
      if (ntime==0) then
! initialize electromagnetic fields from darwin fields
! calculate initial darwin magnetic field
         call mpibpois2(cut,bxyz,ffc,ci,wm,tfield,nx,ny,kstrt)
         wf = 0.0
! calculate initial darwin electric field
         allocate(amu(4,nxe,nypmx),amut(4,nye,kxp),dcut(ndim,nye,kxp))
         amu = 0.0
         call wmpgmjpost2(ppart,amu,kpic,noff,qme,ci,tdjpost,mx,my,mx1, &
     &relativity)
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
         call mpmaxwel2(exyz,bxyz,cut,ffc,affp,ci,dt,wf,wm,tfield,nx,ny,&
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
! updates fxyze, nterf, and ierr, modifies fxyt
      isign = 1
      call wmpfft2rn(fxyze,fxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! transform magnetic field to real space with OpenMP:
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
! push electrons with OpenMP:
! updates ppart and wke, and possibly ncl, ihole, irc
      wke = 0.0
      call wmpbpush2(ppart,fxyze,bxyze,kpic,ncl,ihole,noff,nyp,qbme,dt, &
     &dth,ci,wke,tpush,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
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
         call wmpbpush2(pparti,fxyze,bxyze,kipic,ncl,ihole,noff,nyp,qbmi&
     &,dt,dth,ci,wki,tpush,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
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
! energy diagnostic
      wt = we + wf + wm
      wtot(1) = wt
      wtot(2) = wke
      wtot(3) = wki
      wtot(4) = wke + wt
      wtot(5) = we
      wtot(6) = wf
      wtot(7) = wm
      call PPDSUM(wtot,work,7)
      wke = wtot(2)
      wki = wtot(3)
      we = wtot(5)
      wf = wtot(6)
      wm = wtot(7)
      if ((ntime==0).and.(kstrt==1)) then
         wt = we + wf + wm
         write (*,*) 'Initial Total Field, Kinetic and Total Energies:'
         if (movion==0) then
            write (*,'(3e14.7)') wt, wke, wke + wt
         else
            write (*,'(4e14.7)') wt, wke, wki, wke + wki + wt
         endif
         write (*,*) 'Initial Electrostatic, Transverse Electric and Mag&
     &netic Field Energies:'
         write (*,'(3e14.7)') we, wf, wm
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
      if (kstrt.eq.1) then
         write (*,*)
         write (*,*) 'ntime, relativity = ', ntime, relativity
         write (*,*) 'MPI nodes nvp = ', nvp
         wt = we + wf + wm
         write (*,*) 'Final Total Field, Kinetic and Total Energies:'
         if (movion==0) then
            write (*,'(3e14.7)') wt, wke, wke + wt
         else
            write (*,'(4e14.7)') wt, wke, wki, wke + wki + wt
         endif
         write (*,*) 'Final Electrostatic, Transverse Electric and Magne&
     &tic Field Energies:'
         write (*,'(3e14.7)') we, wf, wm
!
         write (*,*)
         write (*,*) 'initialization time = ', tinit
         write (*,*) 'deposit time = ', tdpost
         write (*,*) 'current deposit time = ', tdjpost
         tdpost = tdpost + tdjpost
         write (*,*) 'total deposit time = ', tdpost
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
