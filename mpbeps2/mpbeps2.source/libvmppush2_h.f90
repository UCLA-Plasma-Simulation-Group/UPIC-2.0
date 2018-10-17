!-----------------------------------------------------------------------
! Interface file for libvmppush2.f
      module libvmppush2_h
      implicit none
!
      interface
         subroutine VPPGPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,nx,ny,&
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
         subroutine VPPGPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm, &
     &dt,ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
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
         subroutine VPPGRPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,  &
     &nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
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
         subroutine VPPGRPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,&
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
         subroutine VPPGPPUSH2ZF(ppart,kpic,dt,ek,nx,ny,idimp,nppmx,    &
     &mxyp1,ipbc)
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
         subroutine VPPGPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,  &
     &nx,ny,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
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
         subroutine VPPGRPPUSH2ZF(ppart,kpic,dt,ci,ek,nx,ny,idimp,nppmx,&
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
         subroutine VPPGRPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ci, &
     &ek,nx,ny,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
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
         subroutine VPPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my, &
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
      end module