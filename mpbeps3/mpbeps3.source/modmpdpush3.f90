!-----------------------------------------------------------------------
!
      module modmpdpush3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpdpush3.f
! mpfwpminx3 calculates maximum and minimum plasma frequency
!            calls PPFWPMINMX32
! mpfwptminx3 calculates maximum and minimum total plasma frequency
!             calls PPFWPTMINMX32
! mpgdjpost3 calculates particle momentum flux and acceleration density
!            calls PPGDCJPPOST2L
! mpgdcjpost3 calculates particle momentum flux, acceleration density
!             and current density
!             calls PPGDCJPPOST32L
! mpgrdjpost3 calculates particle momentum flux and acceleration density
!             calls PPGRDJPPOST32L
! mpgrdcjpost3 calculates particle momentum flux, acceleration density
!              and current density
!              calls PPGRDCJPPOST32L
! mpascfguard3 add scaled vector field to extended periodic field
!              calls MPPASCFGUARD32L
! wmpgdjpost3 generic procedure to calculate particle momentum flux and
!             acceleration density
!             calls mpgrdjpost3 or mpgdjpost3
! wmpgdcjpost3 generic procedure to calculate particle momentum flux,
!              acceleration density, and current
!              calls mpgrdcjpost3 or mpgdcjpost3
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 2, 2017
!
      use libmpdpush3_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpfwpminx3(qe,nyzp,qbme,wpmax,wpmin,nx)
! calculates maximum and minimum plasma frequency
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: qbme
      real, intent(inout) :: wpmax, wpmin
      real, dimension(:,:,:), intent(in) :: qe
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
! extract dimensions
      nxe = size(qe,1); nypmx = size(qe,2); nzpmx = size(qe,3)
      idds = size(nyzp,1)
! call low level procedure
      call PPFWPMINMX32(qe,nyzp,qbme,wpmax,wpmin,nx,nxe,nypmx,nzpmx,idds&
     &)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfwptminx3(qe,qi,nyzp,qbme,qbmi,wpmax,wpmin,nx)
! calculates maximum and minimum total plasma frequency
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: qbme, qbmi
      real, intent(inout) :: wpmax, wpmin
      real, dimension(:,:,:), intent(in) :: qe, qi
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
! extract dimensions
      nxe = size(qe,1); nypmx = size(qe,2); nzpmx = size(qe,3)
      idds = size(nyzp,1)
! call low level procedure
      call PPFWPTMINMX32(qe,qi,nyzp,qbme,qbmi,wpmax,wpmin,nx,nxe,nypmx, &
     &nzpmx,idds)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgdjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,qm,  &
     &qbm,dt,tdcjpost,nx,mx,my,mz,mx1,myp1)
! deposit time derivative of current
      implicit none
      integer, intent(in) :: nx, mx, my, mz, mx1, myp1
      real, intent(in) :: qm, qbm, dt
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,qm,qbm, &
     &dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgdcjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,amu, &
     &qm,qbm,dt,tdcjpost,nx,mx,my,mz,mx1,myp1)
! deposit current and time derivative of current
      implicit none
      integer, intent(in) :: nx, mx, my, mz, mx1, myp1
      real, intent(in) :: qm, qbm, dt
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: cu, dcu
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,amu,qm, &
     &qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,   &
     &idds)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrdjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,qm, &
     &qbm,dt,ci,tdcjpost,nx,mx,my,mz,mx1,myp1)
! deposit relativistic time derivative of current
      implicit none
      integer, intent(in) :: nx, mx, my, mz, mx1, myp1
      real, intent(in) :: qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,qm,qbm,&
     &dt,ci,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds&
     &)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrdcjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,amu,&
     &qm,qbm,dt,ci,tdcjpost,nx,mx,my,mz,mx1,myp1)
! deposit relativistic current and time derivative of current
      implicit none
      integer, intent(in) :: nx, mx, my, mz, mx1, myp1
      real, intent(in) :: qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: cu, dcu
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,amu,qm,&
     &qbm,dt,ci,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &idds)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpascfguard3(dcu,cus,nyzp,q2m0,tdcjpost,nx)
! add scaled vector field to extended periodic field
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: q2m0
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:,:), intent(in) :: cus
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(dcu,2); nypmx = size(dcu,3); nzpmx = size(dcu,4)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPASCFGUARD32L(dcu,cus,nyzp,q2m0,nx,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpgdjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,qm, &
     &qbm,dt,ci,tdcjpost,nx,mx,my,mz,mx1,myp1,relativity)
! generic procedure to calculate particle momentum flux and acceleration
! density
      implicit none
      integer, intent(in) :: nx, mx, my, mz, mx1, myp1, relativity
      real, intent(in) :: qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! updates dcu, amu
      if (relativity==1) then
         call mpgrdjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,qm,qbm,&
     &dt,ci,tdcjpost,nx,mx,my,mz,mx1,myp1)
      else
         call mpgdjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,qm,qbm, &
     &dt,tdcjpost,nx,mx,my,mz,mx1,myp1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpgdcjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,amu,&
     &qm,qbm,dt,ci,tdcjpost,nx,mx,my,mz,mx1,myp1,relativity)
! generic procedure to calculate particle momentum flux, acceleration
! density and current
      implicit none
      integer, intent(in) :: nx, mx, my, mz, mx1, myp1, relativity
      real, intent(in) :: qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: cu, dcu
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! updates cue, dcu, amu
      if (relativity==1) then
         call mpgrdcjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,amu,qm,&
     &qbm,dt,ci,tdcjpost,nx,mx,my,mz,mx1,myp1)
      else
         call mpgdcjpost3(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,amu,qm, &
     &qbm,dt,tdcjpost,nx,mx,my,mz,mx1,myp1)
      endif
      end subroutine
!
      end module
