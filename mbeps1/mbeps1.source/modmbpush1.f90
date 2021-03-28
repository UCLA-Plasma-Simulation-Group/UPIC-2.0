!-----------------------------------------------------------------------
!
      module mbpush1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmbpush1.f
! mbpush13 push magnetized particles with Boris pusher
!          calls GBPPUSH13L
! mbpushf13 push magnetized particles and determine which particles are
!           leaving tile
!           calls GBPPUSHF13L
! mabpush13 push magnetized particles with analytic Boris pusher
!           calls GABPPUSH13L
! mabpushf13 push magnetized particles with analytic Boris pusher and
!            determine which particles are leaving tile
!            calls GABPPUSHF13L
! mrbpush13 push relativistic, magnetized particles with Boris pusher
!           calls GRBPPUSH13L
! mrbpushf13 push relativistic, magnetized particles and determine which
!            particles are leaving tile
!            calls GRBPPUSHF13L
! marbpush13 push relativistic, magnetized particles with analytic Boris
!            pusher
!            calls GARBPPUSH13L
! marbpushf13 push relativistic, magnetized particles with analytic
!             Boris pusher and determine which particles are leaving
!             tile
!             calls GARBPPUSHF13L
! mearbpush13 push relativistic, magnetized particles with exact
!             analytic pusher
!             calls GEARBPPUSH13L
! mearbpush13 push relativistic, magnetized particles with exact
!             analytic pusher and determine which particles are leaving
!             tile
!             calls GEARBPPUSHF13L
! wmbpush1 generic procedure to push magnetized particles with Boris
!          pusher
!          calls mrbpushf13, mbpushf13, mrbpush13, or mbpush13
! wmabpush1 generic procedure to push magnetized particles with analytic
!           Boris pusher
!           calls marbpushf13, mabpushf13, marbpush13, or mabpush13
! wmeabpush1 generic procedure to push magnetized particles with exact 
!            analytic pusher
!            calls mearbpushf13 or mearbpush13
! written by  viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: december 16, 2020
!
      use libmbpush1_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,tpush,nx&
     &,mx,ipbc)
! push magnetized particles with 1d electromagnetic fields
! with Boris pusher
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: omx, qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
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
      call GBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,idimp,nppmx,&
     &nx,mx,nxv,mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc,&
     &ek,tpush,nx,mx,irc)
! push magnetized particles with 1d electromagnetic fields
! with Boris pusher and determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: omx, qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc,ek, &
     &idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mabpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,tpush, &
     &nx,mx,ipbc)
! push magnetized particles with 1d electromagnetic fields with
! analytic Boris pusher
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: omx, qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
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
      call GABPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,idimp,nppmx&
     &,nx,mx,nxv,mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mabpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc&
     &,ek,tpush,nx,mx,irc)
! push magnetized particles with 1d electromagnetic fields with
! analytic Boris pusher and determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: omx, qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GABPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc,ek,&
     &idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,    &
     &tpush,nx,mx,ipbc)
! push relativistic, magnetized particles with 1d electromagnetic fields
! with Boris pusher
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
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
      call GRBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,idimp,  &
     &nppmx,nx,mx,nxv,mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc&
     &,ci,ek,tpush,nx,mx,irc)
! push relativistic, magnetized particles with 1d electromagnetic fields
! with Boris pusher and determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc,ci,&
     &ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine marbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,   &
     &tpush,nx,mx,ipbc)
! push relativistic, magnetized particles with 1d electromagnetic fields
! with analytic Boris pusher
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
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
      call GARBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,idimp, &
     &nppmx,nx,mx,nxv,mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine marbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,  &
     &dtc,ci,ek,tpush,nx,mx,irc)
! push relativistic, magnetized particles with 1d electromagnetic fields
! with analytic Boris pusher, determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GARBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc,ci&
     &,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mearbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,  &
     &tpush,nx,mx,ipbc)
! push relativistic, magnetized particles with 1d electromagnetic fields
! with exact analytic pusher
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
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
      call GEARBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,idimp,&
     &nppmx,nx,mx,nxv,mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mearbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt, &
     &dtc,ci,ek,tpush,nx,mx,irc)
! push relativistic, magnetized particles with 1d electromagnetic fields
! with exact analytic pusher, determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GEARBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc, &
     &ci,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmbpush1(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc, &
     &ci,ek,tpush,nx,mx,ipbc,relativity,plist,irc)
! generic procedure to push magnetized particles with Boris pusher
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, mx, ipbc, relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, wke, ncl, ihole, irc
         if (relativity==1) then
            call mrbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc&
     &,ci,ek,tpush,nx,mx,irc)
         else
            call mbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc,&
     &ek,tpush,nx,mx,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmbpush1 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and wke
         if (relativity==1) then
            call mrbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,    &
     &tpush,nx,mx,ipbc)
         else
            call mbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,tpush,nx&
     &,mx,ipbc)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmabpush1(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc,&
     &ci,ek,tpush,nx,mx,ipbc,relativity,plist,irc)
! generic procedure to push magnetized particles with analytic Boris
! pusher
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, mx, ipbc, relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, wke, ncl, ihole, irc
         if (relativity==1) then
            call marbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,  &
     &dtc,ci,ek,tpush,nx,mx,irc)
         else
            call mabpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc&
     &,ek,tpush,nx,mx,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmabpush1 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and wke
         if (relativity==1) then
            call marbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,   &
     &tpush,nx,mx,ipbc)
         else
            call mabpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,tpush, &
     &nx,mx,ipbc)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmeabpush1(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,dtc&
     &,ci,ek,tpush,nx,mx,ipbc,relativity,plist,irc)
! generic procedure to push magnetized particles with exact analytic
! pusher
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, mx, ipbc, relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: omx, qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, wke, ncl, ihole, irc
         if (relativity==1) then
            call mearbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt, &
     &dtc,ci,ek,tpush,nx,mx,irc)
         else
            call mearbpushf13(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt, &
     &dtc,0.0,ek,tpush,nx,mx,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmeabpush1 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and wke
         if (relativity==1) then
            call mearbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,  &
     &tpush,nx,mx,ipbc)
         else
            call mearbpush13(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,0.0,ek, &
     &tpush,nx,mx,ipbc)
         endif
      endif
      end subroutine
!
      end module
