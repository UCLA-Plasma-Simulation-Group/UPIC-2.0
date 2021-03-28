!-----------------------------------------------------------------------
! Interface file for libmbpush1.f
      module libmbpush1_h
      implicit none
!
      interface
         subroutine GBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,   &
     &idimp,nppmx,nx,mx,nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt&
     &,dtc,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GABPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,  &
     &idimp,nppmx,nx,mx,nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GABPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm, &
     &dt,dtc,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek&
     &,idimp,nppmx,nx,mx,nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm, &
     &dt,dtc,ci,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GARBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci, &
     &ek,idimp,nppmx,nx,mx,nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GARBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,&
     &dt,dtc,ci,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GEARBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,&
     &ek,idimp,nppmx,nx,mx,nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GEARBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm&
     &,dt,dtc,ci,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      end module
