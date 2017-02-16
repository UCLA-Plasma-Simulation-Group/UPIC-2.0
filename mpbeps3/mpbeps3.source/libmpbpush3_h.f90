!-----------------------------------------------------------------------
! Interface file for libmpbpush3.f
      module libmpbpush3_h
      implicit none
!
      interface
         subroutine PPGBPPUSH32L(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt, &
     &dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,    &
     &mxyzp1,idds,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: idds, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPGBPPUSHF32L(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,  &
     &nyzp,qbm,dt,dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx, &
     &mx1,myp1,mxyzp1,ntmax,idds,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout)  :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(26,mxyzp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(inout) :: ihole
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPGRBPPUSH32L(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,&
     &dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1, &  
     &mxyzp1,idds,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: idds, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPGRBPPUSHF32L(ppart,fxyz,bxyz,kpic,ncl,ihole,noff, &
     &nyzp,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,    &
     &nzpmx,mx1,myp1,mxyzp1,ntmax,idds,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nypmx, nzpmx, mx1, myp1, mxyzp1
         integer, intent(in) :: ntmax, idds
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(26,mxyzp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(inout) :: ihole
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      end module