!-----------------------------------------------------------------------
! 3D Periodic Darwin Field MPI/OpenMP PIC test code
! Solves for Potential, Longitudinal Electric Field, Vector Potential,
! and Magnetic Field without smoothing
! written by Viktor K. Decyk, UCLA
      program mp_dfields3
      use modmpfield3
      use mppmod3
      use omplib
      use omppflib3
      use libmpinit3_h
      implicit none
! indx/indy/indz = exponent which determines grid points in x/y/z
! direction: nx = 2**indx, ny = 2**indy, nz = 2**indz.
      integer :: indx =   7, indy =   7, indz =   7
! ndim = number of velocity coordinates = 3
      integer :: ndim = 3
! ax/ay/az = smoothed particle size in x/y/z direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, az = .912871, ci = 0.1
! idps = number of partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      integer :: idps = 4, idds =    2
! we/wm = electrostatic field energy/magnetic field
      real :: we = 0.0, wm = 0.0
!
! declare scalars for standard code
      integer :: j, k, nyps, nzps
      integer :: nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxyzh, nxhyz
      integer :: isign, irc
      real :: affp
!
! declare scalars for MPI code
      integer :: nvpy, nvpz, nvp, idproc, kstrt, kyp, kzp
      integer :: kxyp, kyzp, nypmx, nzpmx
      integer, dimension(2) :: mterf, msg
!
! declare scalars for OpenMP code
      integer :: nvpp
!
! declare arrays for standard code:
! qe = charge density with guard cells
! pot = potential with guard cells
      real, dimension(:,:,:), allocatable :: qe, pot
! cue = ecurrent density with guard cells
      real, dimension(:,:,:,:), allocatable :: cue
! axyze = vector potential with guard cells
! fxyze = longitudinal electric field with guard cells
! bxyze = magnetic field with guard cells
      real, dimension(:,:,:,:), allocatable :: axyze, fxyze, bxyze
! qt = scalar charge density field array in fourier space
! pott = potential in fourier space
      complex, dimension(:,:,:), allocatable :: qt, pott
! cut = scalar charge density field arrays in fourier space
      complex, dimension(:,:,:,:), allocatable :: cut
! axyzt = vector potential in fourier space
! fxyzt = longitudinal electric field in fourier space
! bxyzt = magnetic field in fourier space
      complex, dimension(:,:,:,:), allocatable :: axyzt, fxyzt, bxyzt
! ffc = form factor arrays for poisson solver
      complex, dimension(:,:,:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
!
! declare arrays for MPI code:
! edges(1:2) = lower:upper y boundaries of particle partition
! edges(3:4) = back:front z boundaries of particle partition
      real, dimension(:), allocatable  :: edges
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1:2) = lowermost global gridpoint in y/z
      integer, dimension(:), allocatable :: nyzp, noff
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
! nx/ny/nz = number of grid points in x/y/z direction
      nx = 2**indx; ny = 2**indy; nz = 2**indz
      nxh = nx/2; nyh = max(1,ny/2); nzh = max(1,nz/2)
      nxe = nx + 2; nye = ny + 2; nze = nz + 2
      nxyzh = max(nx,ny,nz)/2; nxhyz = max(nxh,ny,nz)
! affp = normalization constant for field solvers
      affp = 1.0
!      
! nvp = number of MPI ranks
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
! obtain 2D partition (nvpy,nvpz) from nvp:
! nvpy/nvpz = number of processors in y/z
      call FCOMP32(nvp,nx,ny,nz,nvpy,nvpz,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'FCOMP32 error: nvp,nvpy,nvpz=', nvp, nvpy, nvpz
         endif
         call PPEXIT()
         stop
      endif
!
! initialize partition data for MPI code
      allocate(edges(idps),nyzp(idds),noff(idds))
! create real space partitions, boundaries may not be uniform
! edges(1:2) = lower:upper boundary of particle partition in y
! edges(3:4) = back:front boundary of particle partition in z
      k = (kstrt - 1)/nvpy
      j = kstrt - nvpy*k
      k = k + 1
      edges(1) = real(int(real(ny/nvpy)*real(j-1)))
      edges(3) = real(int(real(nz/nvpz)*real(k-1)))
      edges(2) = real(int(real(ny/nvpy)*real(j)))
      edges(4) = real(int(real(nz/nvpz)*real(k)))
      if (j==nvpy) edges(2) = real(ny)
      if (k==nvpz) edges(4) = real(nz)
! noff(1:2) = lowermost global gridpoint in y/z in particle partition
! nyzp(1:2) = number of primary gridpoints in y/z
      noff(1) = int(edges(1))
      noff(2) = int(edges(3))
      nyzp(1) = int(edges(2)) - noff(1)
      nyzp(2) = int(edges(4)) - noff(2)
! nypmx = maximum value of nyp + 1
! nzpmx = maximum value of nzp + 1
      msg = nyzp
      call mpimax(msg,tfmov)
      nypmx = msg(1) + 1
      nzpmx = msg(2) + 1
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
! mterf = number of shifts required by field manager in y/z (0=search)
      mterf = 0
!
! allocate data for standard code
      allocate(qe(nxe,nypmx,nzpmx),pot(nxe,nypmx,nzpmx))
      allocate(fxyze(ndim,nxe,nypmx,nzpmx))
      allocate(cue(ndim,nxe,nypmx,nzpmx))
      allocate(axyze(ndim,nxe,nypmx,nzpmx),bxyze(ndim,nxe,nypmx,nzpmx))
      allocate(qt(nze,kxyp,kyzp),pott(nze,kxyp,kyzp))
      allocate(fxyzt(ndim,nze,kxyp,kyzp))
      allocate(cut(ndim,nze,kxyp,kyzp))
      allocate(axyzt(ndim,nze,kxyp,kyzp),bxyzt(ndim,nze,kxyp,kyzp))
      allocate(ffc(nzh,kxyp,kyzp))
      allocate(mixup(nxhyz),sct(nxyzh))
!
! prepare fft tables
      call mpfft3_init(mixup,sct,indx,indy,indz)
! calculate form factor: ffc
      call mppois3_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,nvpy,nvpz)
!
! deposit point charge (-0.5) at global location
! x = nx/4, y = ny/4, z = nz/4
      qe = 0.0
      j = ny/4 - noff(1)
      k = nz/4 - noff(2)
! deposit to correct local node
      if ((j.ge.1).and.(j.le.nyzp(1)).and.(k.ge.1).and.(k.le.nyzp(2)))  &
     &then
         qe(nx/4,j,k) = -0.5
      endif
!
! transform charge to fourier space:
! updates qt, mterf, and irc, modifies qe
      isign = -1
      call wmpfft3r(qe,qt,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx,indy&
     &,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,irc)
!
! solves 3d poisson's equation for potential in fourier space:
! updates pott, we
      call mppot3(qt,pott,ffc,we,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! transform potential to real space:
! updates pot, mterf, and irc, modifies pott
! moves data in uniform partition to non-uniform partition
      isign = 1
      call wmpfft3r(pot,pott,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,irc)
!
! solves 3d poisson's equation for longitudinal electric field
! updates fxyzt, we
      call mpelfield3(qt,fxyzt,ffc,we,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! transform longitudinal electric force to real space:
! updates fxyze, mterf, and irc, modifies fxyzt
      isign = 1
      call wmpfft3rn(fxyze,fxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,irc)
!
! deposit point current (1,-1,0.5) at global location
! x = nx/4, y = ny/4, z = nz/4
      cue = 0.0
      j = ny/4 - noff(1)
      k = nz/4 - noff(2)
! deposit to correct local node
      if ((j.ge.1).and.(j.le.nyzp(1)).and.(k.ge.1).and.(k.le.nyzp(2)))  &
     &then
         cue(1,nx/4,j,k) = 1.0
         cue(2,nx/4,j,k) = -1.0
         cue(3,nx/4,j,k) = 0.5
      endif
!
! transform current to fourier space:
! updates cut, mterf, and irc, modifies cue
      isign = -1
      call wmpfft3rn(cue,cut,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,irc)
!
! take transverse part of current: updates cut
      call mpcuperp3(cut,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! solves 3d poisson's equation for darwin vector potential
! in fourier space: updates axyzt, wm
      call mpapot3(cut,axyzt,ffc,ci,wm,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
!
! transform vector potential to real space:
! updates axyze, mterf, and irc, modifies axyzt
! moves data in uniform partition to non-uniform partition
      isign = 1
      call wmpfft3rn(axyze,axyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,irc)
!
! solves 3d poisson's equation for darwin magnetic field
! updates bxyzt, wm
      call mpibpois3(cut,bxyzt,ffc,ci,wm,tfield,nx,ny,nz,kstrt,nvpy,nvpz&
     &)
!
! transform magnetic force to real space: updates bxyze, modifies bxyzt
      isign = 1
      call wmpfft3rn(bxyze,bxyzt,noff,nyzp,isign,mixup,sct,tfft,tfmov,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mterf,irc)
!
      if (kstrt==1) then
         write (*,*)
         write (*,*) 'MPI nodes nvpy, nvpz = ', nvpy, nvpz
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
      call PPEXIT()
      end program
