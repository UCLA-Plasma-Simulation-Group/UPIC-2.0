!-----------------------------------------------------------------------
! 1-2/2D Darwin OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mdbeps1
      use minit1
      use mpush1
      use msort1
      use mgard1
      use mfft1
      use mfield1
      use mbpush1
      use mcurd1
      use mdpush1
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
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
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
      integer :: idimp = 4, ipbc = 1, relativity = 0
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.4, omy = 0.0, omz = 0.0
! ndc = number of corrections in darwin iteration
      integer :: ndc = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
      real :: zero = 0.0
! mx = number of grids in x in sorting tiles
      integer :: mx = 32
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! list = (true,false) = list of particles leaving tiles found in push
      logical :: list = .true.
! declare scalars for standard code
      integer :: n, k
      integer :: np, nx, nxh, nxe, nxeh
      integer :: mx1, ntime, nloop, isign
      real :: qbme, affp, q2m0, wpm, wpmax, wpmin
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
! dcu = acceleration density with guard cells
! cus = transverse electric field
! amu = momentum flux with guard cells
      real, dimension(:,:), allocatable :: cue, dcu, cus, amu
! exyze/byze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:), allocatable :: exyze, byze
! ffc, ffe = form factor arrays for poisson solvers
      complex, dimension(:), allocatable :: ffc, ffe
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
      real :: tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0
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
      allocate(cue(2,nxe),dcu(2,nxe),cus(2,nxe),amu(2,nxe))
      allocate(exyze(3,nxe),byze(2,nxe))
      allocate(ffc(nxh),ffe(nxh),mixup(nxh),sct(nxh))
      allocate(kpic(mx1))
!
! prepare fft tables
      call mfft1_init(mixup,sct,indx)
! calculate form factor: ffc
      call mpois1_init(ffc,ax,affp,nx)
! initialize electrons
      call mdistr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,nx,ipbc)
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
! find maximum and minimum initial electron density
      qe = 0.0
      call mpost1(ppart,qe,kpic,qme,tdpost,mx)
      call maguard1(qe,tguard,nx)
      call mfwpminx1(qe,qbme,wpmax,wpmin,nx)
      wpm = 0.5*(wpmax + wpmin)*affp
! accelerate convergence: update wpm
      if (wpm <= 10.0) wpm = 0.75*wpm
      write (*,*) 'wpm=',wpm
      q2m0 = wpm/affp
! calculate form factor: ffe
      call mepois1_init(ffe,ax,affp,wpm,ci,nx)
!
! initialize electric fields
      cus = 0.0; fxe = 0.0
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      write (*,*) 'program mdbeps1'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
!     write (*,*) 'ntime = ', ntime
!
! deposit current with OpenMP: updates cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      call wmdjpost1(ppart,cue,kpic,ncl,ihole,qme,zero,ci,tdjpost,nx,mx,&
     &ipbc,relativity,.false.,irc)
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      call mpost1(ppart,qe,kpic,qme,tdpost,mx)
!
! add guard cells with standard procedure: updates qe, cue
      call maguard1(qe,tguard,nx)
      call macguard1(cue,tguard,nx)
!
! transform charge to fourier space with standard procedure:
! updates qe
      isign = -1
      call mfft1r(qe,isign,mixup,sct,tfft,indx)
!
! calculate longitudinal force/charge in fourier space with standard
! procedure: updates fxe, we
      call mpois1(qe,fxe,ffc,we,tfield,nx)
!
! transform longitudinal electric force to real space with standard
! procedure: updates fxe
      isign = 1
      call mfft1r(fxe,isign,mixup,sct,tfft,indx)
!
! transform current to fourier space with standard procedure:
! updates cue
      isign = -1
      call mfft1rn(cue,isign,mixup,sct,tfft,indx)
!
! calculate magnetic field in fourier space with standard procedure:
! updates byze, wm
      call mbbpois1(cue,byze,ffc,ci,wm,tfield,nx)
!
! transform magnetic force to real space with standard procedure:
! updates byze
      isign = 1
      call mfft1rn(byze,isign,mixup,sct,tfft,indx)
!
! add constant to magnetic field with standard procedure: updates byze
      call mbaddext1(byze,tfield,omy,omz,nx)
!
! copy guard cells with standard procedure: updates fxe, byze
      call mdguard1(fxe,tguard,nx)
      call mcguard1(byze,tguard,nx)
!
! add longitudinal and old transverse electric fields with standard
! procedure: updates exyze
      call maddvrfield1(exyze,cus,fxe,tfield)
!
! deposit electron acceleration density and momentum flux with OpenMP:
! updates dcu, amu
      call dtimer(dtime,itime,-1)
      dcu = 0.0; amu = 0.0
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      call wmgdjpost1(ppart,exyze,byze,dcu,amu,kpic,omx,qme,qbme,dt,ci, &
     &tdcjpost,nx,mx,relativity)
! add old scaled electric field with standard procedure: updates dcu
      call mascfguard1(dcu,cus,q2m0,tdcjpost,nx)
!
! add guard cells with standard procedure: updates dcu, amu
      call macguard1(dcu,tguard,nx)
      call macguard1(amu,tguard,nx)
!
! transform acceleration density and momentum flux to fourier space
! with standard procedure: updates dcu, amu
      isign = -1
      call mfft1rn(dcu,isign,mixup,sct,tfft,indx)
      call mfft1rn(amu,isign,mixup,sct,tfft,indx)
!
! take transverse part of time derivative of current with standard
! procedure: updates dcu
      call madcuperp1(dcu,amu,tfield,nx)
!
! calculate transverse electric field with standard procedure:
! updates cus, wf
      call mepois1(dcu,cus,ffe,affp,ci,wf,tfield,nx)
!
! transform transverse electric field to real space with standard
! procedure: updates cus
      isign = 1
      call mfft1rn(cus,isign,mixup,sct,tfft,indx)
!
! copy guard cells with standard procedure: updates cus
      call mcguard1(cus,tguard,nx)
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxe, updates exyze
! cus needs to be retained for next time step
      call maddvrfield1(exyze,cus,fxe,tfield)
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
      call wmgdcjpost1(ppart,exyze,byze,cue,dcu,amu,kpic,omx,qme,qbme,dt&
     &,ci,tdcjpost,nx,mx,relativity)
! add scaled electric field with standard procedure: updates dcu
      call mascfguard1(dcu,cus,q2m0,tdcjpost,nx)
!
! add guard cells for current, acceleration density, and momentum flux
! with standard procedure: updates cue, dcu, amu
      call macguard1(cue,tguard,nx)
      call macguard1(dcu,tguard,nx)
      call macguard1(amu,tguard,nx)
!
! transform current to fourier space with standard procedure:
! update cue
      isign = -1
      call mfft1rn(cue,isign,mixup,sct,tfft,indx)
!
! calculate magnetic field in fourier space with standard procedure:
! updates byze, wm
      call mbbpois1(cue,byze,ffc,ci,wm,tfield,nx)
!
! transform magnetic force to real space with standard procedure:
! updates byze
      isign = 1
      call mfft1rn(byze,isign,mixup,sct,tfft,indx)
!
! add constant to magnetic field with standard procedure: updates bzye
      call mbaddext1(byze,tfield,omy,omz,nx)
!
! transform acceleration density and momentum flux to fourier space
! with standard procedure: updates dcu, amu
      isign = -1
      call mfft1rn(dcu,isign,mixup,sct,tfft,indx)
      call mfft1rn(amu,isign,mixup,sct,tfft,indx)
!
! take transverse part of time derivative of current with standard
! procedure: updates dcu
      call madcuperp1(dcu,amu,tfield,nx)
!
! calculate transverse electric field with standard procedure:
! updates cus, wf
      call mepois1(dcu,cus,ffe,affp,ci,wf,tfield,nx)
!
! transform transverse electric field to real space with standard
! procedure: updates cus
      isign = 1
      call mfft1rn(cus,isign,mixup,sct,tfft,indx)
!
! copy guard cells with standard procedure: updates byze, cus
      call mcguard1(byze,tguard,nx)
      call mcguard1(cus,tguard,nx)
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxe, updates exyze
! cus needs to be retained for next time step
      call maddvrfield1(exyze,cus,fxe,tfield)
!
      enddo
!
! push particles with OpenMP:
! updates ppart and wke, and possibly ncl, ihole, irc
      wke = 0.0
      call wmbpush1(ppart,exyze,byze,kpic,ncl,ihole,omx,qbme,dt,dt,ci,  &
     &wke,tpush,nx,mx,ipbc,relativity,list,irc)
!
! reorder particles by tile with OpenMP:
! updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
      call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)
!
      if (ntime==0) then
         wt = we + wm
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
      write (*,*) 'ntime, relativity, ndc = ', ntime, relativity, ndc
      wt = we + wm
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
      write (*,*) 'current derivative deposit time = ', tdcjpost
      tdpost = tdpost + tdjpost + tdcjpost
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
