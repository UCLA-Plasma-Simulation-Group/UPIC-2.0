!-----------------------------------------------------------------------
!
      module mgard3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpgard3.f
! mpdguard3x replicates local periodic scalar field
!            calls PPDGUARD32XL
! mpcguard3x replicates local periodic vector field
!            calls PPCGUARD32XL
! mpaguard3x accumulates local periodic scalar field
!            calls PPAGUARD32XL
! mpacguard3x accumulates local periodic vector field
!             calls PPACGUARD32XL
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 16, 2016
!
      use libmpgard3_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpdguard3x(q,nyzp,tguard,nx)
! replicates local periodic scalar field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: q
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPDGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcguard3x(fxyz,nyzp,tguard,nx)
! replicates local periodic vector field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:,:,:,:), intent(inout) :: fxyz
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: ndim, nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(fxyz,1); nxe = size(fxyz,2)
      nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPCGUARD32XL(fxyz,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaguard3x(q,nyzp,tguard,nx)
! accumulates local periodic scalar field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: q
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPAGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpacguard3x(cu,nyzp,tguard,nx)
! accumulates local periodic vector field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tguard
      real, dimension(:,:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: ndim, nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nxe = size(cu,2)
      nypmx = size(cu,3); nzpmx = size(cu,4)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPACGUARD32XL(cu,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
      end module
