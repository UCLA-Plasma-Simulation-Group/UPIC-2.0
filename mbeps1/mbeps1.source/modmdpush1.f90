!-----------------------------------------------------------------------
!
      module mdpush1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmdpush1.f
! mfwpminx1 calculates maximum and minimum plasma frequency
!           calls FWPMINMX1
! mfwptminx1 calculates maximum and minimum total plasma frequency
!            calls FWPTMINMX1
! mgdjpost1 calculates particle momentum flux and acceleration density
!           calls GDJPPOST1L
! mgdcjpost1 calculates particle momentum flux, acceleration density
!            and current density
!            calls GDCJPPOST1L
! mgrdjpost1 calculates particle momentum flux and acceleration density
!            with relativistic particles
!            calls GRDJPPOST1L
! mgrdcjpost1 calculates particle momentum flux, acceleration density
!             and current density with relativistic particles
!             calls GRDCJPPOST1L
! mascfguard1 add scaled vector field to extended periodic field
!             calls ASCFGUARD1L
! wmgdjpost1 generic procedure to calculate particle momentum flux and
!            acceleration density
!            calls mgrdjpost1 or mgdjpost1
! wmgdcjpost1 generic procedure to calculate particle momentum flux,
!             acceleration density, and current
!             calls mpgrdcjpost2 or mpgdcjpost2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 1, 2017
!
      use libmdpush1_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mfwpminx1(qe,qbme,wpmax,wpmin,nx)
! calculates maximum and minimum plasma frequency
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: qbme
      real, intent(inout) :: wpmax, wpmin
      real, dimension(:), intent(in) :: qe
! local data
      integer :: nxe
! extract dimensions
      nxe = size(qe,1)
! call low level procedure
      call FWPMINMX1(qe,qbme,wpmax,wpmin,nx,nxe)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfwptminx1(qe,qi,qbme,qbmi,wpmax,wpmin,nx)
! calculates maximum and minimum total plasma frequency
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: qbme, qbmi
      real, intent(inout) :: wpmax, wpmin
      real, dimension(:), intent(in) :: qe, qi
! local data
      integer :: nxe
! extract dimensions
      nxe = size(qe,1)
! call low level procedure
      call FWPTMINMX1(qe,qi,qbme,qbmi,wpmax,wpmin,nx,nxe)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgdjpost1(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt,   &
     &tdcjpost,nx,mx)
! deposit time derivative of current
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(in) ::  omx, qm, qbm, dt
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: dcu
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GDJPPOST1L(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt,idimp,  &
     &nppmx,nx,mx,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgdcjpost1(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,qbm,dt&
     &,tdcjpost,nx,mx)
! deposit current and time derivative of current
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(in) ::  omx, qm, qbm, dt
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: cu, dcu
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GDCJPPOST1L(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,qbm,dt,    &
     &idimp,nppmx,nx,mx,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgrdjpost1(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt,ci&
     &,tdcjpost,nx,mx)
! deposit relativistic time derivative of current with relativistic
! particles
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(in) ::  omx, qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: dcu
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRDJPPOST1L(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,ci,qbm,dt,    &
     &idimp,nppmx,nx,mx,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgrdcjpost1(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,qbm, &
     &dt,ci,tdcjpost,nx,mx)
! deposit relativistic current and time derivative of current with
! relativistic particles
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(in) ::  omx, qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: cu, dcu
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRDCJPPOST1L(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,ci,qbm,dt,&
     &idimp,nppmx,nx,mx,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mascfguard1(dcu,cus,q2m0,tdcjpost,nx)
! add scaled vector field to extended periodic field
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: q2m0
      real, intent(inout) :: tdcjpost
      real, dimension(:,:), intent(inout) :: dcu
      real, dimension(:,:), intent(in) :: cus
! local data
      integer :: nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(dcu,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgdjpost1(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt,ci&
     &,tdcjpost,nx,mx,relativity)
! generic procedure to calculate particle momentum flux and acceleration
! density
      implicit none
      integer, intent(in) :: nx, mx, relativity
      real, intent(in) :: omx, qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: dcu
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! updates dcu, amu
      if (relativity==1) then
         call mgrdjpost1(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt,ci,  &
     &tdcjpost,nx,mx)
      else
         call mgdjpost1(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt,      &
     &tdcjpost,nx,mx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgdcjpost1(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,qbm, &
     &dt,ci,tdcjpost,nx,mx,relativity)
! generic procedure to calculate particle momentum flux, acceleration
! density and current
      implicit none
      integer, intent(in) :: nx, mx, relativity
      real, intent(in) :: omx, qm, qbm, dt, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: cu, dcu
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! updates cue, dcu, amu
      if (relativity==1) then
         call mgrdcjpost1(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,qbm,dt, &
     &ci,tdcjpost,nx,mx)
      else
         call mgdcjpost1(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,qbm,dt,  &
     &tdcjpost,nx,mx)
      endif
      end subroutine
!
      end module
