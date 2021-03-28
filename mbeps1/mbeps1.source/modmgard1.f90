!-----------------------------------------------------------------------
!
      module mgard1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmgard1.f
! mdguard1 replicates local periodic scalar field
!          calls DGUARD1L
! mcguard1 replicates local periodic vector field
!          calls CGUARD1L or BGUARD1L
! maguard1 accumulates local periodic scalar field
!          calls AGUARD1L
! macguard1 accumulates local periodic vector field
!           calls ACGUARD1L
! mamcguard1 accumulates local periodic tensor field
!            calls AMCGUARD1L
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: december 6, 2017
!
      use libmgard1_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mdguard1(fx,tguard,nx)
! replicates local periodic scalar field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:), intent(inout) :: fx
! local data
      integer :: nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(fx,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call DGUARD1L(fx,nx,nxe)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcguard1(fxy,tguard,nx)
! replicates local periodic vector field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: fxy
! local data
      integer :: ndim, nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(fxy,1); nxe = size(fxy,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call CGUARD1L(fxy,nx,nxe)
      case (3)
         call BGUARD1L(fxy,nx,nxe)
      end select
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine maguard1(q,tguard,nx)
! accumulates local periodic scalar field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:), intent(inout) :: q
! local data
      integer :: nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(q,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call AGUARD1L(q,nx,nxe)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine macguard1(cu,tguard,nx)
! accumulates local periodic vector field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: cu
! local data
      integer :: ndim, nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nxe = size(cu,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call ACGUARD1L(cu,nx,nxe)
      end select
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mamcguard1(amu,tguard,nx)
! accumulates local periodic tensor field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: amu
! local data
      integer :: ndim, nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(amu,1); nxe = size(amu,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call AMCGUARD1L(amu,nx,ndim,nxe)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
      end module
