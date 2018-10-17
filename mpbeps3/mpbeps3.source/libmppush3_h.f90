!-----------------------------------------------------------------------
! Interface file for libmppush3.f
      module libmppush3_h
      implicit none
!
      interface
         subroutine PPPMOVIN3L(part,ppart,kpic,npp,noff,nppmx,idimp,    &
     &npmax,mx,my,mz,mx1,myp1,mxyzp1,idds,irc)
         implicit none
         integer, intent(in) :: npp, nppmx, idimp, npmax, mx, my, mz
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         integer, dimension(mxyzp1), intent(inout) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVIN3LP(part,ppart,kpic,kp,npp,noff,nppmx,idimp,&
     &npmax,mx,my,mz,mx1,myp1,mxyzp1,idds,irc)
         implicit none
         integer, intent(in) :: npp, nppmx, idimp, npmax, mx, my, mz
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         integer, dimension(mxyzp1), intent(inout) :: kpic
         integer, dimension(nppmx,mxyzp1), intent(inout) :: kp
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PPPCOPYOUT3(part,ppart,kpic,npp,npmax,nppmx,idimp,  &
     &mxyzp1,irc)
         implicit none
         integer, intent(in) :: npmax, nppmx, idimp, mxyzp1
         integer, intent(inout) :: npp, irc
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         integer, dimension(mxyzp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPCOPYIN3(part,ppart,kpic,npmax,nppmx,idimp,mxyzp1,&
     &irc)
         implicit none
         integer, intent(in) :: npmax, nppmx, idimp,mxyzp1
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         integer, dimension(mxyzp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPCHECK3L(ppart,kpic,noff,nyzp,idimp,nppmx,nx,mx,my&
     &,mz,mx1,myp1,mzp1,idds,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mz
         integer, intent(in) :: mx1, myp1, mzp1, idds
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1*mzp1), intent(in) :: ppart
         integer, dimension(mx1*myp1*mzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPGPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,    &
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
         subroutine PPGPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm&
     &,dt,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,    &
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
         subroutine PPGRPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ci,ek,&
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds&
     &,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: idds, ipbc
         real, intent(in) :: qbm, dt, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPGRPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,  &
     &qbm,dt,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,   &
     &myp1,mxyzp1,ntmax,idds,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, ci
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
         subroutine PPGPPUSH32ZF(ppart,kpic,dt,ek,idimp,nppmx,nx,ny,nz, &
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
         subroutine PPGPPUSHF32ZF(ppart,kpic,ncl,ihole,noff,nyzp,dt,ek, &
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
         subroutine PPGRPPUSH32ZF(ppart,kpic,dt,ci,ek,idimp,nppmx,nx,ny,&
     &nz,mxyzp1,ipbc)
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
         subroutine PPGRPPUSHF32ZF(ppart,kpic,ncl,ihole,noff,nyzp,dt,ci,&
     &ek,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mxyzp1,ntmax,idds,irc)
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
         subroutine PPGPPOST32L(ppart,q,kpic,noff,qm,nppmx,idimp,mx,my, &
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
      interface
         subroutine SET_PSZERO3(q,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mzp1&
     &)
         implicit none
         integer, intent(in) :: mx, my, mz, nxv, nypmx, nzpmx
         integer, intent(in) :: mx1, myp1, mzp1
         real, dimension(nxv,nypmx,nzpmx), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine SET_PVZERO3(cu,mx,my,mz,ndim,nxv,nypmx,nzpmx,mx1,   &
     &myp1,mzp1)
         implicit none
         integer, intent(in) :: mx, my, mz, ndim, nxv, nypmx, nzpmx
         integer, intent(in) :: mx1, myp1, mzp1
         real, dimension(ndim,nxv,nypmx,nzpmx), intent(inout) :: cu
         end subroutine
      end interface
!
      end module