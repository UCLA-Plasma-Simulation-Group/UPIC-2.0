!-----------------------------------------------------------------------
! Interface file for libmppush2.f
      module libmppush2_h
      implicit none
!
      interface
         subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,    &
     &npmax,mx,my,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, npmax, mx, my, mx1, mxyp1
         integer, intent(in) :: npp, noff
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVIN2LP(part,ppart,kpic,kp,npp,noff,nppmx,idimp,&
     &npmax,mx,my,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, npmax, mx, my, mx1, mxyp1
         integer, intent(in) :: npp, noff
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(inout) :: kpic
         integer, dimension(nppmx,mxyp1), intent(inout) :: kp
         end subroutine
      end interface
!
      interface
         subroutine PPPCOPYOUT2(part,ppart,kpic,npp,npmax,nppmx,idimp,  &
     &mxyp1,irc)
         implicit none
         integer, intent(in) :: npmax, nppmx, idimp, mxyp1
         integer, intent(inout) :: npp, irc
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPCOPYIN2(part,ppart,kpic,npmax,nppmx,idimp,mxyp1, &
     &irc)
         implicit none
         integer, intent(in) :: npmax, nppmx, idimp, mxyp1
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,&
     &mx1,myp1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mx1, myp1
         integer, intent(in) :: noff, nyp
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(in) :: ppart
         integer, dimension(mx1*myp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,nx,ny, &
     &mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: nx, ny, mx, my, idimp, nppmx, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1, ipbc
         integer, intent(in):: noff, nyp
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(2,nxv,nypmx), intent(in) :: fxy
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
     &,ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer, intent(in) :: nx, ny, mx, my, idimp, nppmx, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1, ntmax
         integer, intent(in) :: noff, nyp
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(2,nxv,nypmx), intent(in) :: fxy
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGRPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,nx,&
     &ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: nx, ny, mx, my, idimp, nppmx, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1, ipbc
         integer, intent(in):: noff, nyp
         real, intent(in) :: qbm, dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(2,nxv,nypmx), intent(in) :: fxy
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm, &
     &dt,ci,ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer, intent(in) :: nx, ny, mx, my, idimp, nppmx, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1, ntmax
         integer, intent(in) :: noff, nyp
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(2,nxv,nypmx), intent(in) :: fxy
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGPPUSH2ZF(ppart,kpic,dt,ek,nx,ny,idimp,nppmx,mxyp1&
     &,ipbc)
         implicit none
         integer, intent(in) :: nx, ny, idimp, nppmx, mxyp1, ipbc
         real, intent(in) :: dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,nx,&
     &ny,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer, intent(in) :: nx, ny, mx, my, idimp, nppmx, mx1, mxyp1
         integer, intent(in) :: ntmax, noff, nyp
         integer, intent(inout) :: irc
         real, intent(in) :: dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGRPPUSH2ZF(ppart,kpic,dt,ci,ek,nx,ny,idimp,nppmx, &
     &mxyp1,ipbc)
         implicit none
         integer, intent(in) :: nx, ny, idimp, nppmx, mxyp1, ipbc
         real, intent(in) :: dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek&
     &,nx,ny,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer, intent(in) :: nx, ny, mx, my, idimp, nppmx, mx1, mxyp1
         integer, intent(in) :: ntmax, noff, nyp
         integer, intent(inout) :: irc
         real, intent(in) :: dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,  &
     &nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx, my, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1
         integer, intent(in) :: noff
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nxv,nypmx), intent(inout) :: q
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine SET_PSZERO2(q,mx,my,nxv,nypmx,mx1,myp1)
         implicit none
         integer, intent(in) :: mx, my, nxv, nypmx, mx1, myp1
         real, dimension(nxv,nypmx), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine SET_PVZERO2(cu,mx,my,ndim,nxv,nypmx,mx1,myp1)
         implicit none
         integer, intent(in) :: mx, my, ndim, nxv, nypmx, mx1, myp1
         real, dimension(ndim,nxv,nypmx), intent(inout) :: cu
         end subroutine
      end interface
!
      end module