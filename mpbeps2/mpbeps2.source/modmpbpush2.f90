!-----------------------------------------------------------------------
!
      module mbpush2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpbpush2.f
! mpbpush2 push magnetized particles
!          calls PPGBPPUSH23L or VPPGBPPUSH23L
! mpbpushf2 push magnetized particles and determine which particles are
!           leaving tile
!           calls PPGBPPUSHF23L or VPPGBPPUSHF23L
! mprbpush2 push relativistic, magnetized particles
!           calls PPGRBPPUSH23L or VPPGRBPPUSH23L
! mprbpushf2 push relativistic, magnetized particles and determine which
!            particles are leaving tile
!            calls PPGRBPPUSHF23L or VPPGRBPPUSHF23L
! wmpbpush2 generic procedure to push magnetized particles
!           calls mprbpushf2, mpbpushf2, mprbpush2, or mpbpush2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: august 1, 2018
!
      use libmpbpush2_h
      use libvmpbpush2_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,    &
     &tpush,nx,ny,mx,my,mx1,ipbc,popt)
! push magnetized particles with 2d electromagnetic fields
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt
      integer, intent(in) :: noff, nyp
      real, intent(in) :: qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
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
      select case(popt)
! vector push
      case (2)
         call VPPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,  &
     &idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! standard push
      case default
         call PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,   &
     &idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
     &dtc,ek,tpush,nx,ny,mx,my,mx1,popt,irc)
! push magnetized particles with 2d electromagnetic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, popt
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector push
      case (2)
         call VPPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm, &
     &dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! standard push
      case default
         call PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
     &,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek,&
     &tpush,nx,ny,mx,my,mx1,ipbc,popt)
! push relativistic, magnetized particles with 2d electromagnetic fields
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt
      integer, intent(in) :: noff, nyp
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
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
      select case(popt)
! vector push
      case (2)
         call VPPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci, &
     &ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! standard push
      case default
         call PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek&
     &,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
     &,dtc,ci,ek,tpush,nx,ny,mx,my,mx1,popt,irc)
! push relativistic, magnetized particles with 2d electromagnetic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, popt
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector push
      case (2)
         call VPPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,&
     &dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc&
     &)
! standard push
      case default
         call PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm, &
     &dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc&
     &)
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpbpush2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
     &dtc,ci,ek,tpush,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
! generic procedure to push magnetized particles
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt, relativity
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, wke, ncl, ihole, irc
         if (relativity==1) then
            call mprbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
     &,dtc,ci,ek,tpush,nx,ny,mx,my,mx1,popt,irc)
         else
            call mpbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
     &dtc,ek,tpush,nx,ny,mx,my,mx1,popt,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmpbpush2 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and wke
         if (relativity==1) then
            call mprbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek,&
     &tpush,nx,ny,mx,my,mx1,ipbc,popt)
         else
            call mpbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,    &
     &tpush,nx,ny,mx,my,mx1,ipbc,popt)
         endif
      endif
      end subroutine
!
      end module
