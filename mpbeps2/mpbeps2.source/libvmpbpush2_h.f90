!-----------------------------------------------------------------------
! Interface file for libvmpbpush2.f
      module libvmpbpush2_h
      implicit none
!
      interface
         subroutine VPPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc&
     &,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VPPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp&
     &,qbm,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax, &
     &irc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout)  :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout)  :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout)  :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VPPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,  &
     &dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VPPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,  &
     &nyp,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1, &
     &ntmax,irc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         end subroutine
      end interface
!
      end module