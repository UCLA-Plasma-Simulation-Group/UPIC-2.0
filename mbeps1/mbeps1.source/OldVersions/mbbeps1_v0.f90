!-----------------------------------------------------------------------
! 1-2/2D Electromagnetic OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mbbeps1
      use minit1
      use mpush1
      use msort1
      use mgard1
      use mfft1
      use mfield1
      use mbpush1
      use mcurd1
      use omplib
      implicit none
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx =  18432
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.05, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! vtx/vz0 = thermal/drift velocity of electrons in z direction
      real, parameter :: vtz = 1.0, vz0 = 0.0
! ax = smoothed particle size in x direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ci = 0.1
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 4, ipbc = 1, relativity = 1
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! mx = number of grids in x in sorting tiles
      integer :: mx = 32
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! list = (true,false) = list of particles leaving tiles found in push
      logical :: list = .true.
! declare scalars for standard code
      integer :: n
      integer :: np, nx, nxh, nxe, nxeh
      integer :: mx1, ntime, nloop, isign
      real :: qbme, affp, dth, omt
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, ntmax, npbmx, irc
      integer :: nvp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe = electron charge density with guard cells
! fxe = smoothed longitudinal electric field with guard cells
      real, dimension(:), allocatable :: qe, fxe
! cue = electron current density with guard cells
! fxyze/byze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:), allocatable :: cue, fxyze, byze
! eyz/byz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:), allocatable :: eyz, byz
! ffc = form factor array for poisson solver
      complex, dimension(:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
! arrays required for darwin initial fields
      real, dimension(:,:), allocatable :: amu, dcu
!
! declare arrays for OpenMP (tiled) code:
! ppart = tiled particle array
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), allocatable :: ppart, ppbuff
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
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tpush = 0.0, tsort = 0.0
      double precision :: dtime
!
! start timing initialization
      call dtimer(dtime,itime,-1)
!
      irc = 0
! nvp = number of shared memory nodes (0=default)
      nvp = 0
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvp
! initialize for shared memory parallel processing
      call INIT_OMP(nvp)
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx = number of grid points in x direction
      np = npx; nx = 2**indx; nxh = nx/2
      nxe = nx + 2; nxeh = nxe/2
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx)/real(np)
      dth = 0.0
!
! check for unimplemented features
      if ((list).and.(ipbc.ne.1)) then
         write (*,*) 'ipbc /= 1 and list = .true. not yet supported'
         write (*,*) 'list reset to .false.'
         list = .false.
      endif
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(qe(nxe),fxe(nxe))
      allocate(fxyze(3,nxe),cue(2,nxe),byze(2,nxe))
      allocate(eyz(2,nxeh),byz(2,nxeh))
      allocate(ffc(nxh),mixup(nxh),sct(nxh))
      allocate(kpic(mx1))
!
! prepare fft tables
      call mfft1_init(mixup,sct,indx)
! calculate form factors
      call mpois1_init(ffc,ax,affp,nx)
! initialize electrons
      call mdistr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,nx,ipbc)
!
! initialize transverse electromagnetic fields
      eyz = cmplx(0.0,0.0)
      byz = cmplx(0.0,0.0)
! set magnitude of external transverse magnetic field
      omt = sqrt(omy*omy + omz*omz)
!
! find number of particles in each of mx, tiles: updates kpic, nppmx
      call mdblkp2(part,kpic,nppmx,mx,irc)
!
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmax = xtras*nppmx
      npbmx = xtras*nppmx
      allocate(ppart(idimp,nppmx0,mx1))
      allocate(ppbuff(idimp,npbmx,mx1))
      allocate(ncl(2,mx1))
      allocate(ihole(2,ntmax+1,mx1))
! copy ordered particle data for OpenMP: updates ppart and kpic
      call mpmovin1(part,ppart,kpic,mx,irc)
!
! sanity check
      call mcheck1(ppart,kpic,nx,mx,irc)
!
      if (dt > 0.64*ci) then
         write (*,*) 'Warning: Courant condition may be exceeded!'
      endif
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      write (*,*) 'program mbbeps1'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
!     write (*,*) 'ntime = ', ntime
!
! deposit current with OpenMP:
! updates ppart and cue, and possibly ncl, ihole, irc
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      call wmdjpost1(ppart,cue,kpic,ncl,ihole,qme,dth,ci,tdjpost,nx,mx, &
     &ipbc,relativity,list,irc)
!
! reorder particles by tile with OpenMP:
! updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
      call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      call mpost1(ppart,qe,kpic,qme,tdpost,mx)
!
! add guard cells with standard procedure: updates cue, qe
      call macguard1(cue,tguard,nx)
      call maguard1(qe,tguard,nx)
!
! transform charge to fourier space with standard procedure:
! updates qe
      isign = -1
      call mfft1r(qe,isign,mixup,sct,tfft,indx)
!
! transform current to fourier space with standard procedure:
! updates cue
      isign = -1
      call mfft1rn(cue,isign,mixup,sct,tfft,indx)
!
! calculate electromagnetic fields in fourier space with standard
! procedure: updates eyz, byz
      if (ntime==0) then
! initialize electromagnetic fields from darwin fields
! calculate initial darwin magnetic field
         call mibpois1(cue,byz,ffc,ci,wm,tfield,nx)
         wf = 0.0
! calculate initial darwin electric field
         allocate(amu(2,nxe),dcu(2,nxe))
         amu = 0.0
         call wmgmjpost1(ppart,amu,kpic,qme,ci,tdjpost,mx,relativity)
         call macguard1(amu,tguard,nx)
         isign = -1
         call mfft1rn(amu,isign,mixup,sct,tfft,indx)
         call mdcuperp1(dcu,amu,tfield,nx)
         deallocate(amu)
         call metfield1(dcu,eyz,ffc,ci,wf,tfield,nx)
         deallocate(dcu)
         dth = 0.5*dt
      else
         call mmaxwel1(eyz,byz,cue,ffc,ci,dt,wf,wm,tfield,nx)
      endif
!
! calculate longitudinal force/charge in fourier space with standard
! procedure: updates fxe, we
      call mpois1(qe,fxe,ffc,we,tfield,nx)
!
! add longitudinal and transverse electric fields with standard
! procedure: updates fxyze
      call memfield1(fxyze,fxe,eyz,ffc,tfield,nx)
! copy magnetic field with standard procedure: updates byze
      call mbmfield1(byze,byz,ffc,tfield,nx)
!
! transform electric force to real space with standard procedure:
! updates fxyze
      isign = 1
      call mfft1rn(fxyze,isign,mixup,sct,tfft,indx)
!
! transform magnetic force to real space with standard procedure:
! updates byze
      isign = 1
      call mfft1rn(byze,isign,mixup,sct,tfft,indx)
!
! add constant to magnetic field with OpenMP: updates bxyze
      if (omt > 0.0) call mbaddext1(byze,tfield,omy,omz,nx)
!
! copy guard cells with standard procedure: updates fxyze, byze
      call mcguard1(fxyze,tguard,nx)
      call mcguard1(byze,tguard,nx)
!
! push particles with OpenMP:
! updates ppart and wke, and possibly ncl, ihole, irc
      wke = 0.0
      call wmbpush1(ppart,fxyze,byze,kpic,ncl,ihole,omx,qbme,dt,dth,ci, &
     &wke,tpush,nx,mx,ipbc,relativity,list,irc)
!
! reorder particles by tile with OpenMP:
! updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
      call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)
!
      if (ntime==0) then
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
      write (*,*) 'ntime, relativity = ', ntime, relativity
      wt = we + wf + wm
      write (*,*) 'Final Total Field, Kinetic and Total Energies:'
      write (*,'(3e14.7)') wt, wke, wke + wt
      write (*,*) 'Final Electrostatic, Transverse Electric and Magnetic&
     & Field Energies:'
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
      write (*,*) 'fft time = ', tfft
      write (*,*) 'push time = ', tpush
      write (*,*) 'sort time = ', tsort
      tfield = tfield + tguard + tfft
      write (*,*) 'total solver time = ', tfield
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
!
      stop
      end program
