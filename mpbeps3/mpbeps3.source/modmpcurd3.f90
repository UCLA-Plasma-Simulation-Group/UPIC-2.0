!-----------------------------------------------------------------------
!
      module mcurd3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpcurd3.f
! mpdjpost3 deposit current and update particle positions
!           calls PPGJPPOST32L or VPPGJPPOST32L
! mpdjpostf3 deposit current and update particle positions
!            determine which particles are leaving tile
!            calls PPGJPPOSTF32L or VPPGJPPOSTF32L
! mprdjpost3 deposit current and update relativistic particle positions
!            calls PPGRJPPOSTF32L or VPPGRJPPOSTF32L
! mprdjpostf3 deposit current and update relativistic particle positions
!             determine which particles are leaving tile
!             calls PPGRJPPOSTF32L or VPPGRJPPOSTF32L
! mpgmjpost3 deposits momentum flux
!            calls PPGMJPPOST32L or VPPGMJPPOST32L
! mpgrmjpost3 deposits relativistic momentum flux
!             calls PPGRMJPPOST32L or VPPGRMJPPOST32L
! wmpdjpost3 generic procedure to deposit current and update particle
!            positions
!            calls mprdjpostf3, mpdjpostf3, mprdjpost3, or mpdjpost3
! wmpgmjpost3 generic procedure to deposit momentum flux
!             calls mpgrmjpost3 or mpgmjpost3
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: May 16, 2018
!
      use libmpcurd3_h
      use libvmpcurd3_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpdjpost3(ppart,cu,kpic,noff,qm,dt,tdjpost,nx,ny,nz,mx,&
     &my,mz,mx1,myp1,ipbc,popt)
! deposit current and update particle positions
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      integer, intent(in) :: popt
      real, intent(in) :: qm, dt
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector current deposit
      case (2)
         call VPPGJPPOST32L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny, &
     &nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! standard current deposit
      case default
         call PPGJPPOST32L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny,nz&
     &,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdjpostf3(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,    &
     &tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,popt,irc)
! deposit current and update particle positions
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, popt
      integer, intent(inout) :: irc
      real, intent(in) :: qm, dt
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, ntmax, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
      ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector current deposit
      case (2)
         call VPPGJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,   &
     &nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,    &
     &ntmax,idds,irc)
! standard current deposit
      case default
         call PPGJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,    &
     &nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,    &
     &ntmax,idds,irc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdjpost3(ppart,cu,kpic,noff,qm,dt,ci,tdjpost,nx,ny,nz&
     &,mx,my,mz,mx1,myp1,ipbc,popt)
! deposit current and update relativistic particle positions
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      integer, intent(in) :: popt
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector current deposit
      case (2)
         call VPPGRJPPOST32L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,nx,&
     &ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! standard current deposit
      case default
         call PPGRJPPOST32L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,nx, &
     &ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdjpostf3(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,ci,&
     &tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,popt,irc)
! deposit current and update relativistic particle positions
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, popt
      integer, intent(inout) :: irc
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, ntmax, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
      ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector current deposit
      case (2)
         call VPPGRJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,ci&
     &,nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,   &
     &ntmax,idds,irc)
! standard current deposit
      case default
         call PPGRJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,ci,&
     &nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,    &
     &ntmax,idds,irc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgmjpost3(ppart,amu,kpic,noff,qm,tdcjpost,mx,my,mz,mx1&
     &,myp1,popt)
! deposit momentum flux
      implicit none
      integer, intent(in) :: mx, my, mz, mx1, myp1, popt
      real, intent(in) :: qm
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, nppmx, mdim, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mdim = size(amu,1); nxv = size(amu,2)
      nypmx = size(amu,3); nzpmx = size(amu,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(mdim)
      case (6)
! vector momentum flux deposit
         select case(popt)
         case (2)
            call VPPGMJPPOST32L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my&
     &,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! standard momentum flux deposit
         case default
            call PPGMJPPOST32L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my,&
     &mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
         end select
      case default
         write (*,*) 'mpgmjpost3: unsupported dimension mdim = ', mdim
      end select
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrmjpost3(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,mz&
     &,mx1,myp1,popt)
! deposit relativistic momentum flux
      implicit none
      integer, intent(in) :: mx, my, mz, mx1, myp1, popt
      real, intent(in) :: qm, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, nppmx, mdim, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mdim = size(amu,1); nxv = size(amu,2)
      nypmx = size(amu,3); nzpmx = size(amu,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(mdim)
      case (6)
! vector momentum flux deposit
         select case(popt)
         case (2)
            call VPPGRMJPPOST32L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp, &
     &mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! standard momentum flux deposit
         case default
            call PPGRMJPPOST32L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,mx&
     &,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
         end select
      case default
         write (*,*) 'mpgrmjpost3: unsupported dimension mdim = ', mdim
      end select
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpdjpost3(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,ci, &
     &tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,popt,relativity,plist,irc)
! generic procedure to deposit current and update particle positions
! list = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      integer, intent(in) :: popt, relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, cue, ncl, ihole, irc
         if (relativity==1) then
            call mprdjpostf3(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,ci,&
     &tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,popt,irc)
         else
            call mpdjpostf3(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,    &
     &tdjpost,nx,ny,nz,mx,my,mz,mx1,myp1,popt,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmpdjpost3 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart, cue
         if (relativity==1) then
            call mprdjpost3(ppart,cu,kpic,noff,qm,dt,ci,tdjpost,nx,ny,nz&
     &,mx,my,mz,mx1,myp1,ipbc,popt)
         else
            call mpdjpost3(ppart,cu,kpic,noff,qm,dt,tdjpost,nx,ny,nz,mx,&
     &my,mz,mx1,myp1,ipbc,popt)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpgmjpost3(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,mz&
     &,mx1,myp1,popt,relativity)
! generic procedure to deposit momentum flux
      implicit none
      integer, intent(in) :: mx, my, mz, mx1, myp1, relativity, popt
      real, intent(in) ::  qm, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! updates amu
      if (relativity==1) then
         call mpgrmjpost3(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,mz,  &
     &mx1,myp1,popt)
      else
         call mpgmjpost3(ppart,amu,kpic,noff,qm,tdcjpost,mx,my,mz,mx1,  &
     &myp1,popt)
      endif
      end subroutine
!
      end module
