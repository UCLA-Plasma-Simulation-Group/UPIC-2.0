!-----------------------------------------------------------------------
!
      module modmpgard2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpgard2.f
! mpdguard2x replicates local periodic scalar field
!            calls PPDGUARD2XL
! mpcguard2x replicates local periodic vector field
!            calls PPCGUARD2XL
! mpaguard2x accumulates local periodic scalar field
!            calls PPAGUARD2XL
! mpacguard2x accumulates local periodic vector field
!             calls PPACGUARD2XL
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: january 30, 2016
!
      use libmpgard2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpdguard2x(q,nyp,tguard,nx)
! replicates local periodic scalar field
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: q
! local data
      integer :: nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(q,1); nypmx = size(q,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPDGUARD2XL(q,nyp,nx,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcguard2x(fxy,nyp,tguard,nx)
! replicates local periodic vector field
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: fxy
! local data
      integer :: ndim, nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(fxy,1); nxe = size(fxy,2); nypmx = size(fxy,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPCGUARD2XL(fxy,nyp,nx,ndim,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaguard2x(q,nyp,tguard,nx)
! accumulates local periodic scalar field
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: q
! local data
      integer :: nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(q,1); nypmx = size(q,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPAGUARD2XL(q,nyp,nx,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpacguard2x(cu,nyp,tguard,nx)
! accumulates local periodic vector field
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: cu
! local data
      integer :: ndim, nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nxe = size(cu,2); nypmx = size(cu,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPACGUARD2XL(cu,nyp,nx,ndim,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
      end module
