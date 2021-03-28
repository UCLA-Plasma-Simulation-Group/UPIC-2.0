!-----------------------------------------------------------------------
! Fortran Library for processing guard cells
! 1D OpenMP PIC Codes:
! written by Viktor K. Decyk, UCLA
! DGUARD1L replicate periodic scalar field with linear interpolation
! CGUARD1L replicate periodic 2 component vector field with linear
!          interpolation
! BGUARD1L replicate periodic 3 component vector field with linear
!          interpolation
! AGUARD1L accumulate periodic scalar field with linear interpolation
! ACGUARD1L accumulate periodic vector field with linear interpolation
! AMCGUARD1L accumulate periodic tensor field with linear interpolation
! copyright 2016, regents of the university of california
! update: december 6, 2017
!-----------------------------------------------------------------------
      subroutine DGUARD1L(fx,nx,nxe)
! replicate extended periodic scalar field fx
! linear interpolation
! nx = system length in x direction
! nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real fx
      integer nx, nxe
      dimension fx(nxe)
      fx(nx+1) = fx(1)
      return
      end
!-----------------------------------------------------------------------
      subroutine CGUARD1L(byz,nx,nxe)
! replicate extended periodic vector field byz
! linear interpolation
! nx = system length in x direction
! nxe = second dimension of field arrays, must be >= nx+1
      implicit none
      real byz
      integer nx, nxe
      dimension byz(2,nxe)
! local data
      integer i
! copy edges of extended field
      do 10 i = 1, 2
      byz(i,nx+1) = byz(i,1)
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine BGUARD1L(fxyz,nx,nxe)
! replicate extended periodic vector field fxyz
! linear interpolation
! nx = system length in x direction
! nxe = second dimension of field arrays, must be >= nx+1
      implicit none
      real fxyz
      integer nx, nxe
      dimension fxyz(3,nxe)
! local data
      integer i
! copy edges of extended field
      do 10 i = 1, 3
      fxyz(i,nx+1) = fxyz(i,1)
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine AGUARD1L(q,nx,nxe)
! accumulate extended periodic field
! linear interpolation
! nx = system length in x direction
! nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real q
      integer nx, nxe
      dimension q(nxe)
! accumulate edge of extended field
      q(1) = q(1) + q(nx+1)
      q(nx+1) = 0.0
      return
      end
!-----------------------------------------------------------------------
      subroutine ACGUARD1L(cu,nx,nxe)
! accumulate extended periodic vector field cu
! linear interpolation
! nx = system length in x direction
! nxe = second dimension of field arrays, must be >= nx+1
      implicit none
      real cu
      integer nx, nxe
      dimension cu(2,nxe)
! local data
      integer i
! accumulate edges of extended field
      do 10 i = 1, 2
      cu(i,1) = cu(i,1) + cu(i,nx+1)
      cu(i,nx+1) = 0.0
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine AMCGUARD1L(amu,nx,ndim,nxe)
! accumulate extended periodic tensor field amu
! linear interpolation
! nx = system length in x direction
! ndim = first dimension of field arrays
! nxe = second dimension of field arrays, must be >= nx+1
      implicit none
      real amu
      integer nx, ndim, nxe
      dimension amu(ndim,nxe)
! local data
      integer i
! accumulate edges of extended field
      do 10 i = 1, ndim
      amu(i,1) = amu(i,1) + amu(i,nx+1)
      amu(i,nx+1) = 0.0
   10 continue
      return
      end
