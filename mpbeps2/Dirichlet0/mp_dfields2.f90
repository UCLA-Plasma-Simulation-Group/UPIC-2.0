!-----------------------------------------------------------------------
! 2-1/2D Dirichlet Darwin Field MPI/OpenMP PIC test code
! Solves for Potential, Longitudinal Electric Field, Vector Potential,
! and Magnetic Field without smoothing, using spectral methods with
! Dirichlet (conducting) boundary conditions.
! The solution for the Electric Field, Vector Potential, and Current
! has the following form:
! Fx(x,y) = sum(Fx(kn,km)*cos(kn*x)*sin(km*y)),
! Fy(x,y) = sum(Fy(kn,km)*sin(kn*x)*cos(km*y)),
! Fz(x,y) = sum(Fz(kn,km)*sin(kn*x)*sin(km*y)),
! where kn = n*pi/Nx, and km = m*pi/Ny
! The solution for the Magnetic Field has the following form:
! Bx(x,y) = sum(Bx(kn,km)*sin(kn*x)*cos(km*y)),
! By(x,y) = sum(By(kn,km)*cos(kn*x)*sin(km*y)),
! Bz(x,y) = sum(Bz(kn,km)*cos(kn*x)*cos(km*y)),
! The solution for the Scalar Potential is of the form:
! P(x,y) = sum(P(kn,km)*sin(kn*x)*sin(km*y))
! Two methods of solution are illustrated.
! The first method doubles the grid in each dimension and creates
! symmetric and anti-symmetric image sources.  A normal FFT can then be
! used to solve the equations.  However, the cost is 4 times larger than
! the cost of a periodic solver.
! The second method uses fast sine-cosine transforms directly and
! requires new solvers.  The speed is comparable to that of a periodic
! solver
! The parameter ksolve selects between the two solvers.
! written by Viktor K. Decyk, UCLA
      program mp_dfields2
      use modmpfield2
      use modmpdfield2
      use mppmod2
      use omplib
      use omppflib2
      use mpdplib2
      use modmpfsct2
      use pgraf2
      implicit none
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer, parameter :: indx =   9, indy =   9
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 3
! ax/ay = smoothed particle size in x/y direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, ci = 0.1
! idps = number of partition boundaries
      integer, parameter :: idps = 2
! we/wm = electrostatic field energy/magnetic field
      real :: we = 0.0, wm = 0.0
! idt = (1,2,3) = display (color map,contour plot,both)
! nplt = number of plots per page
      integer :: idt = 2, nplt = 1
! ksolve = (1,2) use method (1,2) to solve the equations
      integer :: ksolve = 1
!
! declare scalars for standard code
      integer :: k
      integer :: nx, ny, nxh, nyh, nxe, nye, nxyh, nxhy, isign, irc
      real :: affp
!
! dirichlet boundary conditions
      integer :: indx1, indy1, nxy, nx2, ny2, nx2e, kxp2, kyp2
!
! declare scalars for MPI code
      integer :: nvp, idproc, kstrt, kxp, kyp, nypmx, nyp, noff, nterf
      integer, dimension(1) :: msg
!
! declare scalars for OpenMP code
      integer :: nvpp
!
! declare arrays for standard code:
! qe = charge density with guard cells
! pot = potential with guard cells
      real, dimension(:,:), allocatable :: qe, pot
! cue = current density with guard cells
      real, dimension(:,:,:), allocatable :: cue
! axyze = vector potential with guard cells
! fxyze = longitudinal electric field with guard cells
! bxyze = magnetic field with guard cells
      real, dimension(:,:,:), allocatable :: fxyze, axyze, bxyze
! qt = charge density in fourier space
! pott = potential in fourier space
      real, dimension(:,:), allocatable :: qt, pott
! cut = current density in fourier space
      real, dimension(:,:,:), allocatable :: cut
! axyt = vector potential in fourier space
! fxyt = longitudinal electric field in fourier space
! bxyt = magnetic field in fourier space
      real, dimension(:,:,:), allocatable :: fxyt, axyt, bxyt
! ffd = form factor array for poisson solvers
      complex, dimension(:,:), allocatable :: ffd
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sctd = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sctd
!
! dirichlet boundary conditions with method 1
      real, dimension(:,:), allocatable :: qe2, pot2
      real, dimension(:,:,:), allocatable :: cue2
      real, dimension(:,:,:), allocatable :: fxye2, axye2, bxye2
      complex, dimension(:,:), allocatable :: qt2, pott2
      complex, dimension(:,:,:), allocatable :: cut2
      complex, dimension(:,:,:), allocatable :: fxyt2, axyt2, bxyt2
      integer, dimension(:), allocatable :: mixup2
      complex, dimension(:), allocatable :: sct2
!
! declare arrays for MPI code:
! edges(1:2) = lower:upper y boundaries of particle partition
      real, dimension(:), allocatable  :: edges
!
! declare and initialize timing data
      real :: tfield = 0.0, tfmov = 0.0
      real, dimension(2) :: tfft = 0.0
!
! nvpp = number of shared memory nodes (0=default)
      nvpp = 0
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvpp
! initialize for shared memory parallel processing
      call INIT_OMP(nvpp)
!
! initialize scalars for standard code
! nx/ny = number of grid points in x/y direction
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 2
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny); nxy = max(nx,ny)
! affp = normalization constant for field solvers
      affp = 1.0
!
! dirichlet boundary conditions with method 1
      indx1 = indx + 1; indy1 = indy + 1
      nx2e = 2*nxe; ny2 = 2*ny; nx2 = 2*nx
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
! open graphics device
      call IPLTCOMM(nplt)
      if (kstrt==1) then
         irc = open_pgraphs(nplt)
! set palette to color wheel
         call STPALIT(2)
      endif
!
! initialize partition data for MPI code
      allocate(edges(idps))
! create real space partitions, boundaries may not be uniform
! edges(1:2) = lower:upper boundary of partition
      edges(1) = real(int(real(ny/nvp)*real(kstrt-1)))
      if (kstrt==nvp) then
         edges(2) = real(ny)
      else
         edges(2) = real(int(real(ny/nvp)*real(kstrt)))
      endif
! parameters describing non-uniform partitions
! noff = lowermost global gridpoint in partition
! nyp = number of primary gridpoints in partition
      noff = int(edges(1))
      nyp = int(edges(2)) - noff
! nypmx = maximum value of nyp + 1
      msg = nyp
      call mpimax(msg,tfmov)
      nypmx = msg(1) + 1
!
! initialize additional scalars for MPI code
! parameters describing uniform partitions
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! nterf = number of shifts required by field manager (0=search)
      nterf = 0
!
! dirichlet boundary conditions
      kxp2 = (nx - 1)/nvp + 1; kyp2 = 2*kyp
!
! allocate and initialize data for standard code
      allocate(qe(nxe,nypmx),pot(nxe,nypmx))
      allocate(cue(ndim,nxe,nypmx))
      allocate(fxyze(ndim,nxe,nypmx),axyze(ndim,nxe,nypmx))
      allocate(bxyze(ndim,nxe,nypmx))
      allocate(qt(nye,kxp2+1),pott(nye,kxp2+1))
      allocate(cut(ndim,nye,kxp2+1))
      allocate(fxyt(ndim,nye,kxp2+1),axyt(ndim,nye,kxp2+1))
      allocate(bxyt(ndim,nye,kxp2+1))
      allocate(ffd(ny,kxp2))
      allocate(mixup(nxhy),sctd(nxy))
!
! dirichlet boundary conditions with method 1
      allocate(qe2(nx2e,kyp2),pot2(nx2e,kyp2))
      allocate(cue2(ndim,nx2e,kyp2))
      allocate(fxye2(ndim,nx2e,kyp2),axye2(ndim,nx2e,kyp2))
      allocate(bxye2(ndim,nx2e,kyp2))
      allocate(qt2(ny2,kxp2),pott2(ny2,kxp2))
      allocate(cut2(ndim,ny2,kxp2))
      allocate(fxyt2(ndim,ny2,kxp2),axyt2(ndim,ny2,kxp2))
      allocate(bxyt2(ndim,ny2,kxp2))
      allocate(mixup2(2*nxhy),sct2(2*nxyh))
!
! prepare fft tables
      call mpfsct2_init(mixup,sctd,indx,indy)
! calculate form factor: ffd
      call mppoisd2_init(ffd,ax,ay,affp,nx,ny,kstrt)
!
! dirichlet boundary conditions with method 1
      call mpfft2_init(mixup2,sct2,indx1,indy1)
!     call mppois2_init(ffd,ax,ay,affp,nx2,ny2,kstrt)
!
! deposit point charge (-0.5) at global location x = nx/4, y = ny/4
      qe = 0.0
      k = ny/4 - noff
! deposit to correct local node
      if ((k >= 1).and.(k <= nyp)) then
         qe(nx/4,k) = -0.5
      endif
!
! moves scalar grid qe from non-uniform to uniform partition
      isign = -1
      call mpfmove2(qe,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,nterf,irc)
      if (irc /= 0) then
         call PPEXIT()
         stop
      endif
!
! Fourier Transform charge from Real Space with Method 1:
      if (ksolve==1) then
! create image charges for dirichlet boundary conditions: updates qe2
         call mpdblsin2d(qe,qe2,tfmov,nx,ny,kstrt,nvp,kyp)
!
! fourier transform charge to fourier space: updates qt2, modifies qe2
         isign = -1
         call mpfft2r(qe2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt, &
     &nvp,kyp2)
!
! Fourier Transform charge from Real Space with Method 2:
      else if (ksolve==2) then
! transform charge with fast sine transform: updates qt, modifies qe
         isign = -1
         call mpfsst2r(qe,qt,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp, &
     &kxp2,kyp)
      endif
!
! Solve Potential with Method 1:
      if (ksolve==1) then
! solves 2d poisson's equation for potential in fourier space:
! updates pott2, we
         call mppot2(qt2,pott2,ffd,we,tfield,nx2,ny2,kstrt)
!
! Solve Potential with Method 2:
      else if (ksolve==2) then
! solves 2d poisson's equation for potential with fast sine-cosine
! transforms: updates pott, we
         call mppotd2(qt,pott,ffd,we,tfield,nx,ny,kstrt)
      endif
!
! Fourier Transform to Real Space with Method 1:
      if (ksolve==1) then
! transform potential to real space: updates pot2, modifies pott2
         isign = 1
         call mpfft2r(pot2,pott2,isign,mixup2,sct2,tfft,indx1,indy1,    &
     &kstrt,nvp,kyp2)
!
! copy solution to primary array pot: updates pot, pot2
         call mphafdbl2d(pot,pot2,tfmov,nx,ny,kstrt,nvp,kyp)
!
! Fourier Transform to Real Space with Method 2:
      else if (ksolve==2) then
! transform potential to real space with fast sine transform:
! updates pot, modifies pott
         isign = 1
         call mpfsst2r(pot,pott,isign,mixup,sctd,tfft,indx,indy,kstrt,  &
     &nvp,kxp2,kyp)
      endif
!
! moves scalar grid pot from uniform to non-uniform partition
      isign = 1
      call mpfmove2(pot,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,nterf,irc)
!
! display potential
      call mplcguard2(pot,nyp,tfield,kstrt,nvp)
      call pdscaler2(pot,nyp,nvp,'POTENTIAL',0,999,0,idt,nx,ny,irc)
!
! Solve Electric Field with Method 1:
      if (ksolve==1) then
! solves 2d poisson's equation for longitudinal electric field
! in fourier space: updates fxyt2, we
         call mpelfield2(qt2,fxyt2,ffd,we,tfield,nx2,ny2,kstrt)
!
! Solve Electric Field with Method 2:
      else if (ksolve==2) then
! solves 2d poisson's equation for longitudinal electric field with
! fast sine-cosine transforms: updates pott, we
        call mpelfieldd2(qt,fxyt,ffd,we,tfield,nx,ny,kstrt)
      endif
!
! Fourier Transform to Real Space with Method 1:
      if (ksolve==1) then
! transform longitudinal electric field to real space:
! updates fxye2, modifies fxyt2
         isign = 1
         call mpfft2rn(fxye2,fxyt2,isign,mixup2,sct2,tfft,indx1,indy1,  &
     &kstrt,nvp,kyp2)
!
! copy solution to primary array fxyze: updates fxyze, fxye2
         call mphafdbl2c(fxyze,fxye2,tfmov,nx,ny,kstrt,nvp,kyp)
!
! Fourier Transform to Real Space with Method 2:
      else if (ksolve==2) then
! transform longitudinal electric field to real space with fast
! sine/cosine transforms: updates fxyze, modifies fxyt
         isign = 1
         call mpfcst2rn(fxyze,fxyt,isign,mixup,sctd,tfft,indx,indy,kstrt&
     &,nvp,kxp2,kyp)
      endif
!
! moves vector grid fxyze from uniform to non-uniform partition
      isign = 1
      call mpfnmove2(fxyze,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,nterf, &
     &irc)
!
! display longitudinal electric field
      call mpnlcguard2(fxyze,nyp,tfield,kstrt,nvp)
      call pdvector2(fxyze,nyp,nvp,'EL FIELD',0,999,0,idt,1,nx,ny,irc)
!
! deposit point current (1,-1,0.5) at global location x = nx/4, y = ny/4
      cue = 0.0
      k = ny/4 - noff
! deposit to correct local node
      if ((k >= 1).and.(k <= nyp)) then
         cue(1,nx/4,k) = 1.0
         cue(2,nx/4,k) = -1.0
         cue(3,nx/4,k) = 0.5
      endif
!
! moves vector grid cue from non-uniform to uniform partition
      isign = -1
      call mpfnmove2(cue,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,nterf,irc&
     &)
      if (irc /= 0) then
         call PPEXIT()
         stop
      endif
!
! Fourier Transform current from Real Space with Method 1:
      if (ksolve==1) then
! create image currents for dirichlet boundary conditions: updates cue2
         call mpdblsin2c(cue,cue2,tfmov,nx,ny,kstrt,nvp,kyp)
!
! fourier transform current to fourier space:
! updates cut2, modifies cue2
         isign = -1
         call mpfft2rn(cue2,cut2,isign,mixup2,sct2,tfft,indx1,indy1,    &
     &kstrt,nvp,kyp2)
!
! Fourier Transform current from Real Space with Method 2:
      else if (ksolve==2) then
! transform current with fast sine/cosine transforms:
! updates cut, modifies cue
         isign = -1
         call mpfcst2rn(cue,cut,isign,mixup,sctd,tfft,indx,indy,kstrt,  &
     &nvp,kxp2,kyp)
      endif
!
! Take transverse part of current with Method 1
! using fourier series: updates cut2
      if (ksolve==1) then
         call mpcuperp2(cut2,tfield,nx2,ny2,kstrt)
!
! Take transverse part of current with Method 2
      else if (ksolve==2) then
! using sine/cosine transforms: updates cut
         call mpcuperpd2(cut,tfield,nx,ny,kstrt)
      endif
!
! Solve Vector Potential with Method 1:
      if (ksolve==1) then
! solves 2d poisson's equation for darwin vector potential
! in fourier space: updates axyt2, wm
         call mpapot2(cut2,axyt2,ffd,ci,wm,tfield,nx2,ny2,kstrt)
!
! Solve Vector Potential with Method 2:
      else if (ksolve==2) then
! solves 2d dirichlet poisson's equation for darwin vector potential
! with fast sine-cosine transforms: updates axyt, wm
         call mpapotd2(cut,axyt,ffd,ci,wm,tfield,nx,ny,kstrt)
      endif
!
! Fourier Transform to Real Space with Method 1:
      if (ksolve==1) then
! transform vector potential to real space:
! updates axye2, modifies axyt2
         isign = 1
         call mpfft2rn(axye2,axyt2,isign,mixup2,sct2,tfft,indx1,indy1,  &
     &kstrt,nvp,kyp2)
!
! copy solution to primary array axyze: updates axyze, axye2
         call mphafdbl2c(axyze,axye2,tfmov,nx,ny,kstrt,nvp,kyp)
!
! Fourier Transform to Real Space with Method 2:
      else if (ksolve==2) then
! transform vector potential to real space with fast sine/cosine
! transforms: updates axyze, modifies axyt
         isign = 1
         call mpfcst2rn(axyze,axyt,isign,mixup,sctd,tfft,indx,indy,kstrt&
     &,nvp,kxp2,kyp)
      endif
!
! moves vector grid axyze from uniform to non-uniform partition
      isign = 1
      call mpfnmove2(axyze,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,nterf, &
     &irc)
!
! display darwin vector potential
      call mpnlcguard2(axyze,nyp,tfield,kstrt,nvp)
      call pdvector2(axyze,nyp,nvp,'VPOTENTIAL',0,999,0,idt,1,nx,ny,irc)
!
! Solve Magnetic Field with Method 1:
      if (ksolve==1) then
! solves 2d poisson's equation for darwin magnetic field
! in fourier space: updates bxyt2, wm
         call mpibpois2(cut2,bxyt2,ffd,ci,wm,tfield,nx2,ny2,kstrt)
!
! Solve Magnetic Field with Method 2:
      else if (ksolve==2) then
! solves 2d dirichlet poisson's equation for darwin magnetic field
! with fast sine-cosine transforms: updates bxyt, wm
         call mpibpoisd2(cut,bxyt,ffd,ci,wm,tfield,nx,ny,kstrt)
      endif
!
! Fourier Transform to Real Space with Method 1:
      if (ksolve==1) then
! transform magnetic field to real space:
! updates bxye2, modifies bxyt2
         isign = 1
         call mpfft2rn(bxye2,bxyt2,isign,mixup2,sct2,tfft,indx1,indy1,  &
     &kstrt,nvp,kyp2)
!
! copy solution to primary array bxyze: updates bxyze, bxye2
         call mphafdbl2c(bxyze,bxye2,tfmov,nx,ny,kstrt,nvp,kyp)
!
! Fourier Transform to Real Space with Method 2:
      else if (ksolve==2) then
! transform transform magnetic to real space with fast sine/cosine
! transforms: updates axyze, modifies axyt
         isign = 1
         call mpfsct2rn(bxyze,bxyt,isign,mixup,sctd,tfft,indx,indy,kstrt&
     &,nvp,kxp2,kyp)
      endif
!
! moves vector grid bxyze from uniform to non-uniform partition
      isign = 1
      call mpfnmove2(bxyze,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,nterf, &
     &irc)
!
! display darwin darwin magnetic field
      call mpnlcguard2(bxyze,nyp,tfield,kstrt,nvp)
      call pdvector2(bxyze,nyp,nvp,'BFIELD',0,999,0,idt,1,nx,ny,irc)
!
      if (kstrt==1) then
         write (*,*)
         write (*,*) 'MPI nodes nvp = ', nvp
!
         write (*,*)
         write (*,*) 'solver time = ', tfield
         write (*,*) 'field move time = ', tfmov
         write (*,*) 'fft and transpose time = ', tfft(1), tfft(2)
         tfield = tfield + tfft(1) + tfmov
         write (*,*) 'total solver time = ', tfield
         write (*,*)
!
      endif
!
! close graphics device
      if (kstrt==1) call close_pgraphs
      call PPEXIT()
      end program
