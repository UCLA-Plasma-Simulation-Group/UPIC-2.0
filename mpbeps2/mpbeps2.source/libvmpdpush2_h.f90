!-----------------------------------------------------------------------
! Interface file for libvmpdpush2.f
      module libvmpdpush2_h
      implicit none
!
      interface
         subroutine VPPGDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm&
     &,qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
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
         subroutine VPPGDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,  &
     &nyp,qm,qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
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
         subroutine VPPGRDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp, &
     &qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
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
         subroutine VPPGRDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff, &
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
      end module