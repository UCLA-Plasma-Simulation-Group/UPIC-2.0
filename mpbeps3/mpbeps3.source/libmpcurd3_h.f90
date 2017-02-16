!-----------------------------------------------------------------------
! Interface file for libmpcurd3.f
      module libmpcurd3_h
      implicit none
!
      interface
         subroutine PPGJPPOST32L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx&
     &,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: idds, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(inout) :: cu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PPGJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm, &
     &dt,nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1, &
     &ntmax,idds,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(inout) :: cu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(26,mxyzp1), intent(inout)  :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(inout)  :: ihole
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPGRJPPOST32L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,    &
     &idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: idds, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(inout) :: cu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PPGRJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,&
     &dt,ci,nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,     &
     &mxyzp1,ntmax,idds,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(inout) :: cu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(26,mxyzp1), intent(inout)  :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(inout)  :: ihole
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPGMJPPOST32L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,&
     &my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, mz, nxv
         integer, intent(in) :: nypmx, nzpmx, mx1, myp1, mxyzp1, idds
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(6,nxv,nypmx,nzpmx), intent(inout) :: amu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PPGRMJPPOST32L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp&
     &,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, mz, nxv
         integer, intent(in) :: nypmx, nzpmx, mx1, myp1, mxyzp1, idds
         real, intent(in) :: qm, ci
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(6,nxv,nypmx,nzpmx), intent(inout) :: amu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      end module