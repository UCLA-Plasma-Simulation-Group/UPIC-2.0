!-----------------------------------------------------------------------
!
      module modmpdpush2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpdpush2.f
! mpfwpminx2 calculates maximum and minimum plasma frequency
!            calls PPFWPMINMX2
! mpfwptminx2 calculates maximum and minimum total plasma frequency
!             calls PPFWPTMINMX2
! mpgdjpost2 calculates particle momentum flux and acceleration density
!            calls PPGDJPPOST2L
! mpgdcjpost2 calculates particle momentum flux, acceleration density
!             and current density
!             calls PPGDCJPPOST2L
! mpgrdjpost2 calculates particle momentum flux and acceleration density
!             with relativistic particles
!             calls PPGRDJPPOST2L
! mpgrdcjpost2 calculates particle momentum flux, acceleration density
!              and current density with relativistic particles
!              calls PPGRDCJPPOST2L
! mpascfguard2 add scaled vector field to extended periodic field
!              calls PPASCFGUARD2L
! wmpgdjpost2 generic procedure to calculate particle momentum flux and
!             acceleration density
!             calls mpgrdjpost2 or mpgdjpost2
! wmpgdcjpost2 generic procedure to calculate particle momentum flux,
!              acceleration density, and current
!              calls mpgrdcjpost2 or mpgdcjpost2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 2, 2017
!
      use libmpdpush2_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpfwpminx2(qe,nyp,qbme,wpmax,wpmin,nx)
! calculates maximum and minimum plasma frequency
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(in) :: qbme
      real, intent(inout) :: wpmax, wpmin
      real, dimension(:,:), intent(in) :: qe
! local data
      integer :: nxe, nypmx
! extract dimensions
      nxe = size(qe,1); nypmx = size(qe,2)
! call low level procedure
      call PPFWPMINMX2(qe,nyp,qbme,wpmax,wpmin,nx,nxe,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfwptminx2(qe,qi,nyp,qbme,qbmi,wpmax,wpmin,nx)
! calculates maximum and minimum total plasma frequency
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(in) :: qbme, qbmi
      real, intent(inout) :: wpmax, wpmin
      real, dimension(:,:), intent(in) :: qe, qi
! local data
      integer :: nxe, nypmx
! extract dimensions
      nxe = size(qe,1); nypmx = size(qe,2)
! call low level procedure
      call PPFWPTMINMX2(qe,qi,nyp,qbme,qbmi,wpmax,wpmin,nx,nxe,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgdjpost2(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,qbm, &
     &dt,tdcjpost,nx,mx,my,mx1)
! deposit time derivative of current
      implicit none
      integer, intent(in) :: nx, mx, my, mx1
      integer, intent(in) :: noff, nyp
      real, intent(in) ::  qm, qbm, dt
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,qbm,dt,  &
     &idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgdcjpost2(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,qm, &
     &qbm,dt,tdcjpost,nx,mx,my,mx1)
! deposit current and time derivative of current
      implicit none
      integer, intent(in) :: nx, mx, my, mx1
      integer, intent(in) :: noff, nyp
      real, intent(in) ::  qm, qbm, dt
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: cu, dcu
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,qm,qbm, &
     &dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrdjpost2(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,qbm,&
     &dt,ci,tdcjpost,nx,mx,my,mx1)
! deposit relativistic time derivative of current with relativistic
! particles
      implicit none
      integer, intent(in) :: nx, mx, my, mx1
      integer, intent(in) :: noff, nyp
      real, intent(in) ::  qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,qbm,dt, &
     &ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrdcjpost2(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,qm,&
     &qbm,dt,ci,tdcjpost,nx,mx,my,mx1)
! deposit relativistic current and time derivative of current with
! relativistic particles
      implicit none
      integer, intent(in) :: nx, mx, my, mx1
      integer, intent(in) :: noff, nyp
      real, intent(in) ::  qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: cu, dcu
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,qm,qbm,&
     &dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpascfguard2(dcu,cus,nyp,q2m0,tdcjpost,nx)
! add scaled vector field to extended periodic field
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(in) :: q2m0
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:), intent(in) :: cus
! local data
      integer :: nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(dcu,2); nypmx = size(dcu,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPASCFGUARD2L(dcu,cus,nyp,q2m0,nx,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpgdjpost2(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,qbm,&
     &dt,ci,tdcjpost,nx,mx,my,mx1,relativity)
! generic procedure to calculate particle momentum flux and acceleration
! density
      implicit none
      integer, intent(in) :: noff, nyp, nx, mx, my, mx1, relativity
      real, intent(in) ::  qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! updates dcu, amu
      if (relativity==1) then
         call mpgrdjpost2(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,qbm,dt,&
     &ci,tdcjpost,nx,mx,my,mx1)
      else
         call mpgdjpost2(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,qbm,dt, &
     &tdcjpost,nx,mx,my,mx1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpgdcjpost2(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,qm,&
     &qbm,dt,ci,tdcjpost,nx,mx,my,mx1,relativity)
! generic procedure to calculate particle momentum flux, acceleration
! density and current
      implicit none
      integer, intent(in) :: noff, nyp, nx, mx, my, mx1, relativity
      real, intent(in) ::  qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: cu, dcu
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! updates cue, dcu, amu
      if (relativity==1) then
         call mpgrdcjpost2(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,qm,qbm&
     &,dt,ci,tdcjpost,nx,mx,my,mx1)
      else
         call mpgdcjpost2(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,qm,qbm,&
     &dt,tdcjpost,nx,mx,my,mx1)
      endif
      end subroutine
!
      end module
