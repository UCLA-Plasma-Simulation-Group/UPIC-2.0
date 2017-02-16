!-----------------------------------------------------------------------
! 3D Electromagnetic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpbbeps3
      use modmpinit3
      use modmppush3
      use modmpbpush3
      use modmpcurd3
      use modmpfield3
      use mppmod3
      use omplib
      use ompplib3
      implicit none
! indx/indy/indz = exponent which determines grid points in x/y/z
! direction: nx = 2**indx, ny = 2**indy, nz = 2**indz.
      integer, parameter :: indx =   7, indy =   7, indz =   7
! npx/npy/npz = number of electrons distributed in x/y/z direction.
      integer, parameter :: npx =  384, npy =   384, npz =   384
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 3
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.035, qme = -1.0
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real, parameter :: vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real, parameter :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! ax/ay/az = smoothed particle size in x/y/z direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, az = .912871, ci = 0.1
! idimp = number of particle coordinates = 6
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 6, ipbc = 1, relativity = 1
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! idps = number of partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      integer :: idps = 4, idds =    2
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! mx/my/mz = number of grids in x/y/z in sorting tiles
! sorting tiles, should be less than or equal to 16
      integer :: mx = 8, my = 8, mz = 8
! fraction of extra particles needed for particle management
      real :: xtras = 0.2
! list = (true,false) = list of particles leaving tiles found in push
      logical :: list = .true.
! declare scalars for standard code
      integer :: n
      integer :: nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh, nnxe
      integer :: nxyzh, nxhyz, mx1, ntime, nloop, isign, ierr
      real :: qbme, affp, dth, omt
      double precision :: np
!
! declare scalars for MPI code
      integer :: nvpy, nvpz, nvp, idproc, kstrt, npmax, kyp, kzp
      integer :: kxyp, kyzp, nypmx, nzpmx, nypmn, nzpmn
      integer :: npp, nps, myp1, mzp1, mxyzp1, mxzyp1
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nbmaxp, ntmaxp, npbmx, irc
      integer :: nvpp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe = electron charge density with guard cells
      real, dimension(:,:,:), allocatable :: qe
! cue = electron current density with guard cells
! fxyze/bxyze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:,:,:), allocatable :: cue, fxyze, bxyze
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
      double precision, dimension(7) :: wtot, work
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
!
! declare arrays for OpenMP code
! ppart = tiled particle array
      real, dimension(:,:,:), allocatable :: ppart
! kpic = number of particles in each tile
      integer, dimension(:), allocatable :: kpic
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
      np =  dble(npx)*dble(npy)*dble(npz)
! nx/ny/nz = number of grid points in x/y direction
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
      dth = 0.0
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
      if ((list).and.(ipbc.ne.1)) then
         if (kstrt==1) then
            write (*,*) 'ipbc /= 1 and list = .true. not yet supported'
            write (*,*) 'list reset to .false.'
            list = .false.
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
! npmax = maximum number of electrons in each partition
      npmax = (np/nvp)*1.25
! myp1/mzp1 = number of tiles in y/z direction
      myp1 = (nyzp(1) - 1)/my + 1; mzp1 = (nyzp(2) - 1)/mz + 1
! mxzyp1 = mx1*max(max(mzp1),max(myp1))
      mxzyp1 = mx1*max((nzpmx-2)/mz+1,(nypmx-2)/my+1)
      mxyzp1 = mx1*myp1*mzp1
! mterf = number of shifts required by field manager in y/z (0=search)
      mterf = 0
!
! allocate data for standard code
      allocate(part(idimp,npmax))
      allocate(qe(nxe,nypmx,nzpmx),fxyze(ndim,nxe,nypmx,nzpmx))
      allocate(cue(ndim,nxe,nypmx,nzpmx),bxyze(ndim,nxe,nypmx,nzpmx))
      allocate(qt(nze,kxyp,kyzp),cut(ndim,nze,kxyp,kyzp))
      allocate(fxyzt(ndim,nze,kxyp,kyzp),bxyzt(ndim,nze,kxyp,kyzp))
      allocate(exyz(ndim,nze,kxyp,kyzp),bxyz(ndim,nze,kxyp,kyzp))
      allocate(ffc(nzh,kxyp,kyzp),mixup(nxhyz),sct(nxyzh))
      allocate(kpic(mxyzp1))
!
! prepare fft tables
      call mpfft3_init(mixup,sct,indx,indy,indz)
! calculate form factors
      call mppois3_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,nvpy,nvpz)
! initialize electrons
      nps = 1
      npp = 0
      call mpdistr3(part,edges,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy, &
     &npz,nx,ny,nz,ipbc,ierr)
! check for particle initialization error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'particle initialization error: ierr=', ierr
         endif
         call PPEXIT()
         stop
      endif
!
! initialize transverse electromagnetic fields
      exyz = cmplx(0.0,0.0)
      bxyz = cmplx(0.0,0.0)
! set magnitude of external magnetic field
      omt = sqrt(omx*omx + omy*omy + omz*omz)
!
! find number of particles in each of mx, my, mz tiles:
! updates kpic, nppmx
      call mpdblkp3(part,kpic,npp,noff,nppmx,mx,my,mz,mx1,myp1,irc)
!
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmaxp = xtras*nppmx
      npbmx = xtras*nppmx
      nbmaxp = 0.125*mxzyp1*npbmx
      allocate(ppart(idimp,nppmx0,mxyzp1))
      allocate(ncl(26,mxyzp1),ihole(2,ntmaxp+1,mxyzp1))
!
! copy ordered particle data for OpenMP: updates ppart and kpic
      call mpmovin3(part,ppart,kpic,npp,noff,mx,my,mz,mx1,myp1,irc)
! sanity check
      call mpcheck3(ppart,kpic,noff,nyzp,nx,mx,my,mz,mx1,myp1,irc)
!
      if (dt > 0.37*ci) then
         write (*,*) 'Warning: Courant condition may be exceeded!'
      endif
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (kstrt==1) write (*,*) 'program mpbbeps3'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
!     if (kstrt==1) write (*,*) 'ntime = ', !
! deposit current with OpenMP:
! updates ppart and cue, and possibly ncl, ihole, irc
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      call wmpdjpost3(ppart,cue,kpic,ncl,ihole,noff,nyzp,qme,dth,ci,    &
     &tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,list,irc)
!
! reorder particles by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
      call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,tsort,tmov,kstrt,nvpy&
     &,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,mxzyp1,list,irc&
     &)
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      call mppost3(ppart,qe,kpic,noff,qme,tdpost,mx,my,mz,mx1,myp1)
!
! add guard cells with OpenMP: updates cue, qe
      call wmpnacguard3(cue,nyzp,tguard,nx,kstrt,nvpy,nvpz)
      call wmpaguard3(qe,nyzp,tguard,nx,kstrt,nvpy,nvpz)
!
! transform charge to fourier space with OpenMP:
! updates qt, mterf, and ierr, modifies qe
      isign = -1
      call wmpfft3r(qe,qt,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx,indy&
     &,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! transform current to fourier space with OpenMP:
! updates cut, mterf, and ierr, modifies cue
      isign = -1
      call wmpfft3rn(cue,cut,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! take transverse part of current with OpenMP: updates cut
      call mpcuperp3(cut,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! calculate electromagnetic fields in fourier space with OpenMP:
! updates exyz, bxyz, wf, wm
      if (ntime==0) then
! initialize electromagnetic fields from darwin fields
! calculate initial darwin magnetic field
         call mpibpois3(cut,bxyz,ffc,ci,wm,tfield,nx,ny,nz,kstrt,nvpy,  &
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
         call mpmaxwel3(exyz,bxyz,cut,ffc,affp,ci,dt,wf,wm,tfield,nx,ny,&
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
! transform electric force to real space with OpenMP:!
! updates fxyze, mterf, and ierr, modifies fxyzt
      isign = 1
      call wmpfft3rn(fxyze,fxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! transform magnetic force to real space with OpenMP: updates bxyze,
! modifies bxyzt
      isign = 1
      call wmpfft3rn(bxyze,bxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,ierr)
!
! add constant to magnetic field with OpenMP: updates bxyze
!     if (omt > 0.0) call mpbaddext2(bxyze,nyp,tfield,omx,omy,omz,nx)
!
! copy guard cells with OpenMP: updates fxyze, bxyze
      call wmpncguard3(fxyze,nyzp,tguard,nx,kstrt,nvpy,nvpz)
      call wmpncguard3(bxyze,nyzp,tguard,nx,kstrt,nvpy,nvpz)
!
! push particles with OpenMP:
! updates ppart and wke, and possibly ncl, ihole, irc
      wke = 0.0
      call wmpbpush3(ppart,fxyze,bxyze,kpic,ncl,ihole,noff,nyzp,qbme,dt,&
     &dth,ci,wke,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,list, &
     &irc)
!
! reorder particles by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
      call ompmove3(ppart,kpic,ncl,ihole,noff,nyzp,tsort,tmov,kstrt,nvpy&
     &,nvpz,nx,ny,nz,mx,my,mz,npbmx,nbmaxp,mx1,myp1,mzp1,mxzyp1,list,irc&
     &)
!
! energy diagnostic
      wt = we + wf + wm
      wtot(1) = wt
      wtot(2) = wke
      wtot(3) = 0.0
      wtot(4) = wke + wt
      wtot(5) = we
      wtot(6) = wf
      wtot(7) = wm
      call PPDSUM(wtot,work,7)
      wke = wtot(2)
      we = wtot(5)
      wf = wtot(6)
      wm = wtot(7)
      if ((ntime==0).and.(kstrt==1)) then
         wt = we + wf + wm
         write (*,*) 'Initial Total Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') wt, wke, wke + wt
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
      if (kstrt==1) then
         write (*,*)
         write (*,*) 'ntime, relativity = ', ntime, relativity
         write (*,*) 'MPI nodes nvpy, nvpz = ', nvpy, nvpz
         wt = we + wf + wm
         write (*,*) 'Final Total Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') wt, wke, wke + wt
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
         wt = 1.0e+09/(real(nloop)*real(np))
         write (*,*) 'Push Time (nsec) = ', tpush*wt
         write (*,*) 'Deposit Time (nsec) = ', tdpost*wt
         write (*,*) 'Sort Time (nsec) = ', tsort*wt
         write (*,*) 'Total Particle Time (nsec) = ', time*wt
      endif
!
      call PPEXIT()
      end program
