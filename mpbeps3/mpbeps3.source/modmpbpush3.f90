!-----------------------------------------------------------------------
!
      module modmpbpush3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpbpush3.f
! mpbpush3 push magnetized particles
!          calls PPGBPPUSH32L
! mpbpushf3 push magnetized particles and determine which particles are
!           leaving tile
!           calls PPGBPPUSHF32L
! mprbpush3 push relativistic, magnetized particles
!           calls PPGRBPPUSH32L
! mprbpushf2 push relativistic, magnetized particles and determine which
!            particles are leaving tile
!            calls PPGRBPPUSHF32L
! wmpbpush3 generic procedure to push magnetized particles
!           calls mprbpushf3, mpbpushf3, mprbpush3, or mpbpush3
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: january 28, 2017
!
      use libmpbpush3_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpbpush3(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,dtc,ek, &
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc)
! push magnetized particles with 3d electromagnetic fields
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      real, intent(in) :: qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
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
      call PPGBPPUSH32L(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,dtc,ek,   &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds&
     &,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpbpushf3(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,qbm,&
     &dt,dtc,ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
! push magnetized particles with 3d electromagnetic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
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
      nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
      ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGBPPUSHF32L(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt&
     &,dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,   &
     &mxyzp1,ntmax,idds,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprbpush3(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,dtc,ci,&
     &ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc)
! push relativistic, magnetized particles with 3d electromagnetic fields
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
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
      call PPGRBPPUSH32L(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,dtc,ci,ek&
     &,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,   &
     &idds,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprbpushf3(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,qbm&
     &,dt,dtc,ci,ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
! push relativistic, magnetized particles with 3d electromagnetic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
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
      nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
      ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRBPPUSHF32L(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,qbm, &
     &dt,dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,   &
     &myp1,mxyzp1,ntmax,idds,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpbpush3(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,qbm,&
     &dt,dtc,ci,ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,    &
     &plist,irc)
! generic procedure to push magnetized particles
! list = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      integer, intent(in) :: relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, wke, ncl, ihole, irc
         if (relativity==1) then
            call mprbpushf3(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,qbm&
     &,dt,dtc,ci,ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
         else
            call mpbpushf3(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,qbm,&
     &dt,dtc,ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
         endif
! do not also calculate list of particles leaving tile
         if (irc /= 0) then
            write (*,*) 'info:wmpbpush3 overflow: irc=', irc
         endif
      else
! updates ppart and wke
         if (relativity==1) then
            call mprbpush3(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,dtc,ci,&
     &ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc)
         else
            call mpbpush3(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,dtc,ek, &
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc)
         endif
      endif
      end subroutine
!
      end module
