!-----------------------------------------------------------------------
! Fortran Library for processing local guard cells
! 2D MPI/OpenMP PIC Codes:
! written by Viktor K. Decyk, UCLA
! PPDGUARD2XL replicate local periodic scalar field with linear
!             interpolation
! PPCGUARD2XL replicate local periodic vector field with linear
!             interpolation
! PPAGUARD2XL accumulate local periodic scalar field with linear
!             interpolation
! PPACGUARD2XL accumulate local periodic vector field with linear
!              interpolation
! copyright 2016, regents of the university of california
! update: january 30, 2016
!-----------------------------------------------------------------------
      subroutine PPDGUARD2XL(q,nyp,nx,nxe,nypmx)
! replicate extended periodic scalar field in x direction
! linear interpolation, for distributed data
! nyp = number of primary (complete) gridpoints in particle partition
! nx = system length in x direction
! nxe = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells
      implicit none
      integer nyp, nx, nxe, nypmx
      real q
      dimension q(nxe,nypmx)
! local data
      integer k, myp1
! replicate edges of extended field
      myp1 = nyp + 1
      do 10 k = 1, myp1
      q(nx+1,k) = q(1,k)
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PPCGUARD2XL(fxy,nyp,nx,ndim,nxe,nypmx)
! replicate extended periodic vector field in x direction
! linear interpolation, for distributed data
! nyp = number of primary (complete) gridpoints in particle partition
! nx = system length in x direction
! ndim = leading dimension of array fxy
! nxe = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells
      implicit none
      integer nyp, nx, ndim, nxe, nypmx
      real fxy
      dimension fxy(ndim,nxe,nypmx)
! local data
      integer i, k, myp1
! replicate edges of extended field
      myp1 = nyp + 1
      do 20 k = 1, myp1
      do 10 i = 1, ndim
      fxy(i,nx+1,k) = fxy(i,1,k)
   10 continue
   20 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PPAGUARD2XL(q,nyp,nx,nxe,nypmx)
! accumulate extended periodic scalar field in x direction
! linear interpolation, for distributed data
! nyp = number of primary (complete) gridpoints in particle partition
! nx = system length in x direction
! nxe = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells
      implicit none
      integer nyp, nx, nxe, nypmx
      real q
      dimension q(nxe,nypmx)
! local data
      integer k, myp1
! accumulate edges of extended field
      myp1 = nyp + 1
      do 10 k = 1, myp1
      q(1,k) = q(1,k) + q(nx+1,k)
      q(nx+1,k) = 0.0
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PPACGUARD2XL(cu,nyp,nx,ndim,nxe,nypmx)
! accumulate extended periodic vector field in x direction
! linear interpolation, for distributed data
! nyp = number of primary (complete) gridpoints in particle partition
! nx = system length in x direction
! ndim = leading dimension of array fxy
! nxe = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells
      implicit none
      real cu
      integer nyp, nx, ndim, nxe, nypmx
      dimension cu(ndim,nxe,nypmx)
! local data
      integer i, k, myp1
! accumulate edges of extended field
      myp1 = nyp + 1
      do 20 k = 1, myp1
      do 10 i = 1, ndim
      cu(i,1,k) = cu(i,1,k) + cu(i,nx+1,k)
      cu(i,nx+1,k) = 0.0
   10 continue
   20 continue
      return
      end
