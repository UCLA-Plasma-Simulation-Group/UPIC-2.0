!-----------------------------------------------------------------------
! 1D Electrostatic OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mbeps1
      use minit1
      use mpush1
      use msort1
      use mgard1
      use mfft1
      use mfield1
      use omplib
      implicit none
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx = 18432
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx = thermal velocity of electrons in x direction
! vx0 = drift velocity of electrons in x direction.
      real, parameter :: vtx = 1.0, vx0 = 0.0
! ax = smoothed particle size in x direction
! ci = reciprocal of velocity of light
      real :: ax = .912871, ci = 0.1
! idimp = number of particle coordinates = 2
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 2, ipbc = 1, relativity = 0
! wke/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, we = 0.0, wt = 0.0
! mx = number of grids in x in sorting tiles
      integer :: mx = 32
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! list = (true,false) = list of particles leaving tiles found in push
      logical :: list = .true.
! declare scalars for standard code
      integer :: n
      integer :: np, nx, nxh, nxe
      integer :: mx1, ntime, nloop, isign
      real :: qbme, affp
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, ntmax, npbmx, irc
      integer :: nvp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe = electron charge density with guard cells
      real, dimension(:), allocatable :: qe
! fxe = smoothed electric field with guard cells
      real, dimension(:), allocatable :: fxe
! ffc = form factor array for poisson solver
      complex, dimension(:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
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
      real :: tpush = 0.0, tsort = 0.0
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
      nxe = nx + 2
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx)/real(np)
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
      allocate(ffc(nxh),mixup(nxh),sct(nxh))
      allocate(kpic(mx1))
!
! prepare fft tables
      call mfft1_init(mixup,sct,indx)
! calculate form factors
      call mpois1_init(ffc,ax,affp,nx)
! initialize electrons
      call mdistr2(part,vtx,vx0,npx,nx,ipbc)
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
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      write (*,*) 'program mbeps1'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
!     write (*,*) 'ntime = ', ntime
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      call mpost1(ppart,qe,kpic,qme,tdpost,mx)
!
! add guard cells with standard procedure: updates qe
      call maguard1(qe,tguard,nx)
!
! transform charge to fourier space with standard procedure:
! updates qe
      isign = -1
      call mfft1r(qe,isign,mixup,sct,tfft,indx)
!
! calculate force/charge in fourier space with standard procedure:
! updates fxe, we
      isign = -1
      call mpois1(qe,fxe,ffc,we,tfield,nx)
!
! transform force to real space with standard procedure: updates fxe
      isign = 1
      call mfft1r(fxe,isign,mixup,sct,tfft,indx)
!
! copy guard cells with standard procedure: updates fxe
      call mdguard1(fxe,tguard,nx)
!
! push particles with OpenMP:
      wke = 0.0
! updates part, wke and possibly ncl, ihole, and irc
      call wmpush1(ppart,fxe,kpic,ncl,ihole,qbme,dt,ci,wke,tpush,nx,mx, &
     &ipbc,relativity,list,irc)
!
! reorder particles by tile with OpenMP:
! updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
      call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)
!
      if (ntime==0) then
         write (*,*) 'Initial Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') we, wke, wke + we
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
      write (*,*) 'Final Field, Kinetic and Total Energies:'
      write (*,'(3e14.7)') we, wke, wke + we
!
      write (*,*)
      write (*,*) 'initialization time = ', tinit
      write (*,*) 'deposit time = ', tdpost
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
      write (*,*)
!
      stop
      end program
