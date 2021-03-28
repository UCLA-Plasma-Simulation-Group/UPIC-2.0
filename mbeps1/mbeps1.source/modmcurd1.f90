!-----------------------------------------------------------------------
!
      module mcurd1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmcurd1.f
! mdjpost1 deposit current and update particle positions
!          calls GJPPOST1L
! mdjpostf1 deposit current and update particle positions
!           determine which particles are leaving tile
!           calls GJPPOSTF1L
! mrdjpost1 deposit current and update relativistic particle positions
!           calls GRJPPOST1L
! mrdjpostf1 deposit current and update relativistic particle positions
!            determine which particles are leaving tile
!            calls GRJPPOSTF1L
! mgmjpost1 deposits momentum flux
!           calls GMJPPOST1L
! mgrmjpost1 deposits relativistic momentum flux
!            calls GRMJPPOST1L
! wmdjpost1 generic procedure to deposit current and update particle
!           positions
!           calls mrdjpostf1, mdjpostf1, mrdjpost1, or mdjpost1
! wmgmjpost1 generic procedure to deposit momentum flux
!            calls mgrmjpost1, or mgmjpost1
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: december 19, 2016
!
      use libmcurd1_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mdjpost1(ppart,cu,kpic,qm,dt,tdjpost,nx,mx,ipbc)
! deposit current and update particle positions
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: qm, dt
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GJPPOST1L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,mx,nxv,mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mdjpostf1(ppart,cu,kpic,ncl,ihole,qm,dt,tdjpost,nx,mx, &
     &irc)
! deposit current and update particle positions
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: qm, dt
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GJPPOSTF1L(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,idimp,nx,mx,  &
     &nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrdjpost1(ppart,cu,kpic,qm,dt,ci,tdjpost,nx,mx,ipbc)
! deposit current and update relativistic particle positions
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRJPPOST1L(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,mx,nxv,mx1, &
     &ipbc)
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrdjpostf1(ppart,cu,kpic,ncl,ihole,qm,dt,ci,tdjpost,nx,&
     &mx,irc)
! deposit current and update relativistic particle positions
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRJPPOSTF1L(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx,idimp,nx, &
     &mx,nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgmjpost1(ppart,amu,kpic,qm,tdcjpost,mx)
! deposit momentum flux
      implicit none
      integer, intent(in) :: mx
      real, intent(in) ::  qm
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(amu,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GMJPPOST1L(ppart,amu,kpic,qm,nppmx,idimp,mx,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgrmjpost1(ppart,amu,kpic,qm,ci,tdcjpost,mx)
! deposit relativistic momentum flux
      implicit none
      integer, intent(in) :: mx
      real, intent(in) ::  qm, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(amu,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRMJPPOST1L(ppart,amu,kpic,qm,ci,nppmx,idimp,mx,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmdjpost1(ppart,cu,kpic,ncl,ihole,qm,dt,ci,tdjpost,nx, &
     &mx,ipbc,relativity,list,irc)
! generic procedure to deposit current and update particle positions
! list = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, mx, ipbc, relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: list
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! also calculate list of particles leaving tile
      if (list) then
! updates ppart, cue, ncl, ihole, irc
         if (relativity==1) then
            call mrdjpostf1(ppart,cu,kpic,ncl,ihole,qm,dt,ci,tdjpost,nx,&
     &mx,irc)
         else
            call mdjpostf1(ppart,cu,kpic,ncl,ihole,qm,dt,tdjpost,nx,mx, &
     &irc)
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart, cue
         if (relativity==1) then
            call mrdjpost1(ppart,cu,kpic,qm,dt,ci,tdjpost,nx,mx,ipbc)
         else
            call mdjpost1(ppart,cu,kpic,qm,dt,tdjpost,nx,mx,ipbc)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgmjpost1(ppart,amu,kpic,qm,ci,tdcjpost,mx,relativity)
! generic procedure to deposit momentum flux
      implicit none
      integer, intent(in) :: mx, relativity
      real, intent(in) ::  qm, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! updates amu
      if (relativity==1) then
         call mgrmjpost1(ppart,amu,kpic,qm,ci,tdcjpost,mx)
      else
         call mgmjpost1(ppart,amu,kpic,qm,tdcjpost,mx)
      endif
      end subroutine
!
      end module
