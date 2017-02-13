!-----------------------------------------------------------------------
! Interface file for libmpdpush2.f
      module libmpdpush2_h
      implicit none
!
      interface
         subroutine PPFWPMINMX2(qe,nyp,qbme,wpmax,wpmin,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe,nypmx), intent(in) :: qe
         end subroutine
      end interface
!
      interface
         subroutine PPFWPTMINMX2(qe,qi,nyp,qbme,qbmi,wpmax,wpmin,nx,nxe,&
     &nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, intent(in) :: qbme, qbmi
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe,nypmx), intent(in) :: qe, qi
         end subroutine
      end interface
!
      interface
         subroutine PPGDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,&
     &qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1
         real, intent(in) ::  qm, qbm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(3,nxv,nypmx), intent(inout) :: dcu
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp&
     &,qm,qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1
         real, intent(in) ::  qm, qbm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(3,nxv,nypmx), intent(inout) :: cu, dcu
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm&
     &,qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1
         real, intent(in) ::  qm, qbm, dt, ci
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(3,nxv,nypmx), intent(inout) :: dcu
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,  &
     &nyp,qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1
         real, intent(in) ::  qm, qbm, dt, ci
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(3,nxv,nypmx), intent(inout) :: cu, dcu
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPASCFGUARD2L(dcu,cus,nyp,q2m0,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, intent(in) :: q2m0
         real, dimension(3,nxe,nypmx), intent(inout) :: dcu
         real, dimension(3,nxe,nypmx), intent(in) :: cus
         end subroutine
      end interface
!
      end module