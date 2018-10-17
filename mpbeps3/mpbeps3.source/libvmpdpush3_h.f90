!-----------------------------------------------------------------------
! Interface file for libvmpdpush3.f
      module libvmpdpush3_h
      implicit none
!
      interface
         subroutine VPPGDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,  &
     &amu,qm,qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,   &
     &mxyzp1,idds)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mz, nxv
         integer, intent(in) :: nypmx, nzpmx, mx1, myp1, mxyzp1, idds
         real, intent(in) :: qm, qbm, dt
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         real, dimension(3,nxv,nypmx,nzpmx), intent(inout) :: dcu
         real, dimension(6,nxv,nypmx,nzpmx), intent(inout) :: amu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,  &
     &dcu,amu,qm,qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1&
     &,mxyzp1,idds)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mz, nxv
         integer, intent(in) :: nypmx, nzpmx, mx1, myp1, mxyzp1, idds
         real, intent(in) :: qm, qbm, dt
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         real, dimension(3,nxv,nypmx,nzpmx), intent(inout) :: cu, dcu
         real, dimension(6,nxv,nypmx,nzpmx), intent(inout) :: amu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGRDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu, &
     &amu,qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,&
     &mxyzp1,idds)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mz, nxv
         integer, intent(in) :: nypmx, nzpmx, mx1, myp1, mxyzp1, idds
         real, intent(in) :: qm, qbm, dt, ci
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         real, dimension(3,nxv,nypmx,nzpmx), intent(inout) :: dcu
         real, dimension(6,nxv,nypmx,nzpmx), intent(inout) :: amu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine VPPGRDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu, &
     &dcu,amu,qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1, &
     &myp1,mxyzp1,idds)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mz, nxv
         integer, intent(in) :: nypmx, nzpmx, mx1, myp1, mxyzp1, idds
         real, intent(in) :: qm, qbm, dt, ci
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         real, dimension(3,nxv,nypmx,nzpmx), intent(inout) :: cu, dcu
         real, dimension(6,nxv,nypmx,nzpmx), intent(inout) :: amu
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      end module