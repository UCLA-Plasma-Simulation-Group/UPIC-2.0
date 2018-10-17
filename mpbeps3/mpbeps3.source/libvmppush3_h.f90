!-----------------------------------------------------------------------
! Interface file for libvmppush3.f
      module libvmppush3_h
      implicit none
!
      interface
         subroutine VPPGPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,   &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds&
     &,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: idds, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,  &
     &qbm,dt,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1, &
     &mxyzp1,ntmax,idds,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(26,mxyzp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(inout) :: ihole
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGRPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ci,ek&
     &,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,   &
     &idds,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: idds, ipbc
         real, intent(in) :: qbm, ci, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGRPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp, &
     &qbm,dt,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,   &
     &myp1,mxyzp1,ntmax,idds,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, ci, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(26,mxyzp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(inout) :: ihole
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGPPUSH32ZF(ppart,kpic,dt,ek,idimp,nppmx,nx,ny,nz,&
     &mxyzp1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mxyzp1, ipbc
         real, intent(in) :: dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         integer, dimension(mxyzp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VPPGPPUSHF32ZF(ppart,kpic,ncl,ihole,noff,nyzp,dt,ek,&
     &idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mxyzp1,ntmax,idds,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: mx1, myp1, mxyzp1, ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: dt
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(26,mxyzp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(inout) :: ihole
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGRPPUSH32ZF(ppart,kpic,dt,ci,ek,idimp,nppmx,nx,ny&
     &,nz,mxyzp1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mxyzp1, ipbc
         real, intent(in) :: dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         integer, dimension(mxyzp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VPPGRPPUSHF32ZF(ppart,kpic,ncl,ihole,noff,nyzp,dt,ci&
     &,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mxyzp1,ntmax,idds,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: mx1, myp1, mxyzp1, ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(26,mxyzp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(inout) :: ihole
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGPPOST32L(ppart,q,kpic,noff,qm,nppmx,idimp,mx,my,&
     &mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, mz, nxv
         integer, intent(in) :: nypmx, nzpmx, mx1, myp1, mxyzp1, idds
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(nxv,nypmx,nzpmx), intent(inout) :: q
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      end module