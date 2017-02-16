!-----------------------------------------------------------------------
! Fortran Library for processing local guard cells
! 3D MPI/OpenMP PIC Codes:
! written by Viktor K. Decyk, UCLA
! PPDGUARD32XL replicate local periodic scalar field with linear
!              interpolation
! PPCGUARD32XL replicate local periodic vector field with linear
!              interpolation
! PPAGUARD32XL accumulate local periodic scalar field with linear
!              interpolation
! PPACGUARD32XL accumulate local periodic vector field with linear
!               interpolation
! copyright 2016, regents of the university of california
! update: february 15, 2016
!-----------------------------------------------------------------------
      subroutine PPDGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds)
! replicate extended periodic scalar field in x direction
! linear interpolation, for distributed data with 2D decomposition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      integer nyzp
      real q
      dimension q(nxe,nypmx,nzpmx), nyzp(idds)
! local data
      integer k, l, myp1, mzp1
! replicate edges of extended field
      myp1 = nyzp(1) + 1
      mzp1 = nyzp(2) + 1
!$OMP PARALLEL DO PRIVATE(k,l)
      do 20 l = 1, mzp1
      do 10 k = 1, myp1
      q(nx+1,k,l) = q(1,k,l)
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPCGUARD32XL(fxyz,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
! replicate extended periodic vector field in x direction
! linear interpolation, for distributed data with 2D decomposition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! nx = system length in x direction
! ndim = leading dimension of field array fxyz
! nxe = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, ndim, nxe, nypmx, nzpmx, idds
      integer nyzp
      real fxyz
      dimension fxyz(ndim,nxe,nypmx,nzpmx), nyzp(idds)
! local data
      integer i, k, l, myp1, mzp1
! replicate edges of extended field
      myp1 = nyzp(1) + 1
      mzp1 = nyzp(2) + 1
!$OMP PARALLEL DO PRIVATE(i,k,l)
      do 30 l = 1, mzp1
      do 20 k = 1, myp1
      do 10 i = 1, ndim
      fxyz(i,nx+1,k,l) = fxyz(i,1,k,l)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPAGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds)
! accumulate extended periodic scalar field in x direction
! linear interpolation, for distributed data with 2D decomposition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      integer nyzp
      real q
      dimension q(nxe,nypmx,nzpmx), nyzp(idds)
      integer k, l, myp1, mzp1
! accumulate edges of extended field
      myp1 = nyzp(1) + 1
      mzp1 = nyzp(2) + 1
!$OMP PARALLEL DO PRIVATE(k,l)
      do 20 l = 1, mzp1
      do 10 k = 1, myp1
      q(1,k,l) = q(1,k,l) + q(nx+1,k,l)
      q(nx+1,k,l) = 0.0
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPACGUARD32XL(cu,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
! accumulate extended periodic vector field in x direction
! linear interpolation, for distributed data with 2D decomposition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! nx = system length in x direction
! ndim = leading dimension of array cu
! nxe = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, ndim, nxe, nypmx, nzpmx, idds
      integer nyzp
      real cu
      dimension cu(ndim,nxe,nypmx,nzpmx), nyzp(idds)
      integer i, k, l, myp1, mzp1
! accumulate edges of extended field
      myp1 = nyzp(1) + 1
      mzp1 = nyzp(2) + 1
!$OMP PARALLEL DO PRIVATE(i,k,l)
      do 30 l = 1, mzp1
      do 20 k = 1, myp1
      do 10 i = 1, ndim
      cu(i,1,k,l) = cu(i,1,k,l) + cu(i,nx+1,k,l)
      cu(i,nx+1,k,l) = 0.0
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end