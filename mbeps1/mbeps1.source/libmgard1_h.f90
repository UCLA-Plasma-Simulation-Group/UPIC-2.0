!-----------------------------------------------------------------------
! Interface file for libmgard1.f
      module libmgard1_h
      implicit none
!
      interface
         subroutine DGUARD1L(fx,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(nxe), intent(inout) :: fx
         end subroutine
      end interface
!
      interface
         subroutine CGUARD1L(byz,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(2,nxe), intent(inout) :: byz
         end subroutine
      end interface
!
      interface
         subroutine BGUARD1L(fxyz,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(3,nxe), intent(inout) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine AGUARD1L(q,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(nxe), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD1L(cu,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(2,nxe), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine AMCGUARD1L(amu,nx,ndim,nxe)
         implicit none
         integer, intent(in) :: nx, ndim, nxe
         real, dimension(ndim,nxe), intent(inout) :: amu
         end subroutine
      end interface
!
      end module
