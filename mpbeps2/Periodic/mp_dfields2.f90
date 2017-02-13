!-----------------------------------------------------------------------
! 2-1/2D Periodic Darwin Field MPI/OpenMP PIC test code
! Solves for Potential, Longitudinal Electric Field, Vector Potential,
! and Magnetic Field without smoothing
! written by Viktor K. Decyk, UCLA
      program mp_dfields2
      use modmpfield2
      use mppmod2
      use omplib
      use omppflib2
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
!
! declare scalars for standard code
      integer :: k
      integer :: nx, ny, nxh, nyh, nxe, nye, nxyh, nxhy, isign, irc
      real :: affp
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
      real, dimension(:,:,:), allocatable :: axyze, fxyze, bxyze
! qt = charge density in fourier space
! pott = potential in fourier space
      complex, dimension(:,:), allocatable :: qt, pott
! cut = current density in fourier space
      complex, dimension(:,:,:), allocatable :: cut
! axyt = vector potential in fourier space
! fxyt = longitudinal electric field in fourier space
! bxyt = magnetic field in fourier space
      complex, dimension(:,:,:), allocatable :: axyt, fxyt, bxyt
! ffc = form factor array for poisson solver
      complex, dimension(:,:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
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
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
! affp = normalization constant for field solvers
      affp = 1.0
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
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
!
! allocate and initialize data for standard code
      allocate(qe(nxe,nypmx),pot(nxe,nypmx),fxyze(ndim,nxe,nypmx))
      allocate(cue(ndim,nxe,nypmx))
      allocate(axyze(ndim,nxe,nypmx),bxyze(ndim,nxe,nypmx))
      allocate(qt(nye,kxp),pott(nye,kxp),fxyt(ndim,nye,kxp))
      allocate(cut(ndim,nye,kxp))
      allocate(axyt(ndim,nye,kxp),bxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),mixup(nxhy),sct(nxyh))
!
! prepare fft tables
      call mpfft2_init(mixup,sct,indx,indy)
! calculate form factor: ffc
      call mppois2_init(ffc,ax,ay,affp,nx,ny,kstrt)
!
! deposit point charge (-0.5) at global location x = nx/4, y = ny/4
      qe = 0.0
      k = ny/4 - noff
! deposit to correct local node
      if ((k.ge.1).and.(k.le.nyp)) then
         qe(nx/4,k) = -0.5
      endif
!
! transform charge to fourier space:
! updates qt, nterf, and irc, modifies qe
! moves data in non-uniform partition to uniform partition
      isign = -1
      call wmpfft2r(qe,qt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,indy,&
     &kstrt,nvp,kyp,ny,nterf,irc)
!
! solves 2d poisson's equation for potential in fourier space:
! updates pott, we
      call mppot2(qt,pott,ffc,we,tfield,nx,ny,kstrt)
!
! transform potential to real space:
! updates pot, nterf, and irc, modifies pott
! moves data in uniform partition to non-uniform partition
      isign = 1
      call wmpfft2r(pot,pott,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,  &
     &indy,kstrt,nvp,kyp,ny,nterf,irc)
!
! display potential
      call wmpcguard2(pot,nyp,tfield,nx,kstrt,nvp)
      call pdscaler2(pot,nyp,nvp,'POTENTIAL',0,999,0,idt,nx,ny,irc)
!
! solves 2d poisson's equation for longitudinal electric field
! in fourier space: updates fxyt, we
      call mpelfield2(qt,fxyt,ffc,we,tfield,nx,ny,kstrt)
!
! transform longitudinal electric field to real space:
! updates fxyze, nterf, and irc, modifies fxyt
! moves data in uniform partition to non-uniform partition
      isign = 1
      call wmpfft2rn(fxyze,fxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,irc)
!
! display longitudinal electric field
      call wmpncguard2(fxyze,nyp,tfield,nx,kstrt,nvp)
      call pdvector2(fxyze,nyp,nvp,'EL FIELD',0,999,0,idt,1,nx,ny,irc)
!
! deposit point current (1,-1,0.5) at global location x = nx/4, y = ny/4
      cue = 0.0
      k = ny/4 - noff
      if ((k.ge.1).and.(k.le.nyp)) then
         cue(1,nx/4,k) = 1.0
         cue(2,nx/4,k) = -1.0
         cue(3,nx/4,k) = 0.5
      endif
!
! transform current to fourier space with OpenMP:
! updates cut, nterf, and irc, modifies cue
! moves data in non-uniform partition to uniform partition
      isign = -1
      call wmpfft2rn(cue,cut,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,  &
     &indy,kstrt,nvp,kyp,ny,nterf,irc)
!
! take transverse part of current with OpenMP:
! updates cut
      call mpcuperp2(cut,tfield,nx,ny,kstrt)
!
! solves 2d poisson's equation for darwin vector potential
! in fourier space: updates axyt, wm
      call mpapot2(cut,axyt,ffc,ci,wm,tfield,nx,ny,kstrt)
!
! transform vector potential to real space:
! updates axyze, nterf, and irc, modifies axyt
! moves data in uniform partition to non-uniform partition
      isign = 1
      call wmpfft2rn(axyze,axyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,irc)
!
! display darwin vector potential
      call wmpncguard2(axyze,nyp,tfield,nx,kstrt,nvp)
      call pdvector2(axyze,nyp,nvp,'VPOTENTIAL',0,999,0,idt,1,nx,ny,irc)
!
! solves 2d poisson's equation for darwin magnetic field
! in fourier space: updates bxyt, wm
      call mpibpois2(cut,bxyt,ffc,ci,wm,tfield,nx,ny,kstrt)
!
! transform magnetic field to real space:
! updates bxyze, nterf, and irc, modifies bxyt
! moves data in uniform partition to non-uniform partition
      isign = 1
      call wmpfft2rn(bxyze,bxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,irc)
!
! display darwin darwin magnetic field
      call wmpncguard2(bxyze,nyp,tfield,nx,kstrt,nvp)
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
