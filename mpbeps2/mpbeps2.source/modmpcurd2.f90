!-----------------------------------------------------------------------
!
      module mcurd2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpcurd2.f
! mpdjpost2 deposit current and update particle positions
!           calls PPGJPPOST2L or VPPGJPPOST2L
! mpdjpostf2 deposit current and update particle positions
!            determine which particles are leaving tile
!            calls PPGJPPOSTF2L or VPPGJPPOSTF2L
! mprdjpost2 deposit current and update relativistic particle positions
!            calls PPGRJPPOST2L or VPPGRJPPOST2L
! mprdjpostf2 deposit current and update relativistic particle positions
!             determine which particles are leaving tile
!             calls PPGRJPPOSTF2L or VPPGRJPPOSTF2L
! mpgmjpost2 deposits momentum flux
!            calls PPGMJPPOST2L or VPPGMJPPOST2L
! mpgrmjpost2 deposits relativistic momentum flux
!             calls PPGRMJPPOST2L or VPPGRMJPPOST2L
! wmpdjpost2 generic procedure to deposit current and update particle
!            positions
!            calls mprdjpostf2, mpdjpostf2, mprdjpost2, or mpdjpost2
! wmpgmjpost2 generic procedure to deposit momentum flux
!             calls mpgrmjpost2, or mpgmjpost2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: august 1, 2018
!
      use libmpcurd2_h
      use libvmpcurd2_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpdjpost2(ppart,cu,kpic,noff,qm,dt,tdjpost,nx,ny,mx,my,&
     &mx1,ipbc,popt)
! deposit current and update particle positions
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt
      integer, intent(in) :: noff
      real, intent(in) :: qm, dt
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector current deposit
      case (2)
         call VPPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny,mx&
     &,my,nxv,nypmx,mx1,mxyp1,ipbc)
! standard current deposit
      case default
         call PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny,mx,&
     &my,nxv,nypmx,mx1,mxyp1,ipbc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,     &
     &tdjpost,nx,ny,mx,my,mx1,popt,irc)
! deposit current and update particle positions
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, popt
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qm, dt
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector current deposit
      case (2)
         call VPPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,nppmx&
     &,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! standard current deposit
      case default
         call PPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,nppmx,&
     &idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdjpost2(ppart,cu,kpic,noff,qm,dt,ci,tdjpost,nx,ny,mx&
     &,my,mx1,ipbc,popt)
! deposit current and update relativistic particle positions
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt
      integer, intent(in) :: noff
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector current deposit
      case (2)
         call VPPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,nx, &
     &ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! standard current deposit
      case default
         call PPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,nx,ny&
     &,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci, &
     &tdjpost,nx,ny,mx,my,mx1,popt,irc)
! deposit current and update relativistic particle positions
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, popt
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector current deposit
      case (2)
         call VPPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci, &
     &nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! standard current deposit
      case default
         call PPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci,  &
     &nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgmjpost2(ppart,amu,kpic,noff,qm,tdcjpost,mx,my,mx1,  &
     &popt)
! deposit momentum flux
      implicit none
      integer, intent(in) :: mx, my, mx1, popt
      integer, intent(in) :: noff
      real, intent(in) ::  qm
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mdim, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mdim = size(amu,1); nxv = size(amu,2); nypmx = size(amu,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(mdim)
      case (4)
         select case(popt)
! vector momentum flux deposit
         case (2)
            call VPPGMJPPOST2L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my,&
     &nxv,nypmx,mx1,mxyp1)
! standard momentum flux deposit
         case default
            call PPGMJPPOST2L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my, &
     &nxv,nypmx,mx1,mxyp1)
         end select
      case default
         write (*,*) 'mpgmjpost2: unsupported dimension mdim = ', mdim
      end select
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,  &
     &mx1,popt)
! deposit relativistic momentum flux
      implicit none
      integer, intent(in) :: mx, my, mx1, popt
      integer, intent(in) :: noff
      real, intent(in) ::  qm, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mdim, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mdim = size(amu,1); nxv = size(amu,2); nypmx = size(amu,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(mdim)
      case (4)
         select case(popt)
! vector momentum flux deposit
         case (2)
            call VPPGRMJPPOST2L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,mx&
     &,my,nxv,nypmx,mx1,mxyp1)
! standard momentum flux deposit
         case default
            call PPGRMJPPOST2L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,mx,&
     &my,nxv,nypmx,mx1,mxyp1)
         end select
      case default
         write (*,*) 'mpgrmjpost2: unsupported dimension mdim = ', mdim
      end select
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpdjpost2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci,  &
     &tdjpost,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
! generic procedure to deposit current and update particle positions
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt, relativity
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, cue, ncl, ihole, irc
         if (relativity==1) then
            call mprdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci, &
     &tdjpost,nx,ny,mx,my,mx1,popt,irc)
         else
            call mpdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,     &
     &tdjpost,nx,ny,mx,my,mx1,popt,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmpdjpost2 ihole overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart, cue
         if (relativity==1) then
            call mprdjpost2(ppart,cu,kpic,noff,qm,dt,ci,tdjpost,nx,ny,mx&
     &,my,mx1,ipbc,popt)
         else
            call mpdjpost2(ppart,cu,kpic,noff,qm,dt,tdjpost,nx,ny,mx,my,&
     &mx1,ipbc,popt)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpgmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,  &
     &mx1,popt,relativity)
! generic procedure to deposit momentum flux
      implicit none
      integer, intent(in) :: noff, mx, my, mx1, popt, relativity
      real, intent(in) ::  qm, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
! updates amu
      if (relativity==1) then
         call mpgrmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,mx1, &
     &popt)
      else
         call mpgmjpost2(ppart,amu,kpic,noff,qm,tdcjpost,mx,my,mx1,popt)
      endif
      end subroutine
!
      end module
