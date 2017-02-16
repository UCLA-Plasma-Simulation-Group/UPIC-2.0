!-----------------------------------------------------------------------
! Interface file for libmpdpush3.f
      module libmpdpush3_h
      implicit none
!
      interface
         subroutine PPFWPMINMX32(qe,nyzp,qbme,wpmax,wpmin,nx,nxe,nypmx, &  
     &nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe,nypmx,nzpmx), intent(in) :: qe
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPFWPTMINMX32(qe,qi,nyzp,qbme,qbmi,wpmax,wpmin,nx,  &
     &nxe,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, intent(in) :: qbme, qbmi
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe,nypmx,nzpmx), intent(in) :: qe, qi
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPGDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu&
     &,qm,qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1&
     &,idds)
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
         subroutine PPGDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu&
     &,amu,qm,qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,  &
     &mxyzp1,idds)
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
         subroutine PPGRDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,  &
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
         subroutine PPGRDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,  &
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
      interface
         subroutine MPPASCFGUARD32L(dcu,cus,nyzp,q2m0,nx,nxe,nypmx,nzpmx&
     &,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, intent(in) :: q2m0
         real, dimension(3,nxe,nypmx,nzpmx), intent(inout) :: dcu
         real, dimension(3,nxe,nypmx,nzpmx), intent(in) :: cus
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      end module