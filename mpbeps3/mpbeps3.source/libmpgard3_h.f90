!-----------------------------------------------------------------------
! Interface file for libmpgard3.f
      module libmpgard3_h
      implicit none
!
      interface
         subroutine PPDGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, dimension(nxe,nypmx,nzpmx), intent(inout) :: q
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPCGUARD32XL(fxyz,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, ndim, nxe, nypmx, nzpmx, idds
         real, dimension(ndim,nxe,nypmx,nzpmx), intent(inout) :: fxyz
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPAGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, dimension(nxe,nypmx,nzpmx), intent(inout) :: q
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPACGUARD32XL(cu,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, ndim, nxe, nypmx, nzpmx, idds
         real, dimension(ndim,nxe,nypmx,nzpmx), intent(inout) :: cu
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      end module