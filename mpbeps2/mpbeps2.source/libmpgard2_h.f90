!-----------------------------------------------------------------------
! Interface file for libmpgard2.f
      module libmpgard2_h
      implicit none
!
      interface
         subroutine PPDGUARD2XL(q,nyp,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, dimension(nxe,nypmx), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine PPCGUARD2XL(fxy,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx), intent(inout) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine PPAGUARD2XL(q,nyp,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, dimension(nxe,nypmx), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine PPACGUARD2XL(cu,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx), intent(inout) :: cu
         end subroutine
      end interface
!
      end module