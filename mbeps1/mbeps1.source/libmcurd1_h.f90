!-----------------------------------------------------------------------
! Interface file for libmcurd1.f
      module libmcurd1_h
      implicit none
!
      interface
         subroutine GJPPOST1L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,mx,nxv,&
     &mx1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(2,nxv), intent(inout) :: cu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GJPPOSTF1L(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,idimp&
     &,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(2,nxv), intent(inout) :: cu
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRJPPOST1L(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,mx,&
     &nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(2,nxv), intent(inout) :: cu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRJPPOSTF1L(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx, &
     &idimp,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(2,nxv), intent(inout) :: cu
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GMJPPOST1L(ppart,amu,kpic,qm,nppmx,idimp,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, nxv, mx1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(2,nxv), intent(inout) :: amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRMJPPOST1L(ppart,amu,kpic,qm,ci,nppmx,idimp,mx,nxv,&
     &mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, nxv, mx1
         real, intent(in) :: qm, ci
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(2,nxv), intent(inout) :: amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      end module
