!-----------------------------------------------------------------------
! Interface file for libmpcurd2.f
      module libmpcurd2_h
      implicit none
!
      interface
         subroutine PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,&
     &ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(inout) :: cu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,&
     &nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer, intent(in) :: noff, nyp, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(inout) :: cu
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp&
     &,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(inout) :: cu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt&
     &,ci,nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer, intent(in) :: noff, nyp, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(inout) :: cu
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGMJPPOST2L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx, &
     &my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, mx, my, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1
         real, intent(in) ::  qm
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRMJPPOST2L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,&
     &mx,my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, mx, my, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1
         real, intent(in) ::  qm, ci
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      end module