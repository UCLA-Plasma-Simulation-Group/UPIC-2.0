!-----------------------------------------------------------------------
! Interface file for libmpush1.f
      module libmpush1_h
      implicit none
!
      interface
         subroutine PPMOVIN1L(part,ppart,kpic,nppmx,idimp,nop,mx,mx1,irc&
     &)
         implicit none
         integer, intent(in) :: nppmx, idimp, nop, mx, mx1
         integer, intent(inout) :: irc
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCOPYOUT1(part,ppart,kpic,np,nop,nppmx,idimp,mx1,  &
     &irc)
         implicit none
         integer, intent(in) :: nop, nppmx, idimp, mx1
         integer, intent(inout) :: np, irc
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCOPYIN1(part,ppart,kpic,nop,nppmx,idimp,mx1,irc)
         implicit none
         integer, intent(in) :: nop, nppmx, idimp, mx1
         integer, intent(inout) :: irc
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK1L(ppart,kpic,idimp,nppmx,nx,mx,mx1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSH1L(ppart,fx,kpic,qbm,dt,ek,idimp,nppmx,nx,mx, &
     &nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(nxv), intent(in) :: fx
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSHF1L(ppart,fx,kpic,ncl,ihole,qbm,dt,ek,idimp,  &
     &nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(nxv), intent(in) :: fx
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRPPUSH1L(ppart,fx,kpic,qbm,dt,ci,ek,idimp,nppmx,nx,&
     &mx,nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: qbm, dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(nxv), intent(in) :: fx
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRPPUSHF1L(ppart,fx,kpic,ncl,ihole,qbm,dt,ci,ek,    &
     &idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(nxv), intent(in) :: fx
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPUSH1ZF(ppart,kpic,dt,ek,idimp,nppmx,nx,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx1, ipbc
         real, intent(in) :: dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPUSHF1ZF(ppart,kpic,ncl,ihole,dt,ek,idimp,nppmx,nx,&
     &mx,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1,ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine RPPUSH1ZF(ppart,kpic,dt,ci,ek,idimp,nppmx,nx,mx1,   &
     &ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx1, ipbc
         real, intent(in) :: ci, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine RPPUSHF1ZF(ppart,kpic,ncl,ihole,dt,ci,ek,idimp,nppmx&
     &,nx,mx,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: ci, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST1L(ppart,q,kpic,qm,nppmx,idimp,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, nxv, mx1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nxv), intent(inout) :: q
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      end module
