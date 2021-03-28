!-----------------------------------------------------------------------
! Interface file for libmdpush1.f
      module libmdpush1_h
      implicit none
!
      interface
         subroutine FWPMINMX1(qe,qbme,wpmax,wpmin,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe), intent(in) :: qe
         end subroutine
      end interface
!
      interface
         subroutine FWPTMINMX1(qe,qi,qbme,qbmi,wpmax,wpmin,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: qbme, qbmi
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe), intent(in) :: qe, qi
         end subroutine
      end interface
!
      interface
         subroutine GDJPPOST1L(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt&
     &,idimp,nppmx,nx,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1
         real, intent(in) :: omx, qm, qbm, dt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         real, dimension(2,nxv), intent(inout) :: dcu, amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GDCJPPOST1L(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,  &
     &qbm,dt,idimp,nppmx,nx,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1
         real, intent(in) :: omx, qm, qbm, dt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         real, dimension(2,nxv), intent(inout) :: cu, dcu, amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRDJPPOST1L(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,ci,  &
     &qbm,dt,idimp,nppmx,nx,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1
         real, intent(in) :: omx, qm, qbm, dt, ci
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         real, dimension(2,nxv), intent(inout) :: dcu, amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRDCJPPOST1L(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm, &
     &ci,qbm,dt,idimp,nppmx,nx,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1
         real, intent(in) :: omx, qm, qbm, dt, ci
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         real, dimension(2,nxv), intent(inout) :: cu, dcu, amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: q2m0
         real, dimension(2,nxe), intent(inout) :: dcu
         real, dimension(2,nxe), intent(in) :: cus
         end subroutine
      end interface
!
      end module
