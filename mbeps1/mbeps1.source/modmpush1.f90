!-----------------------------------------------------------------------
!
      module mpush1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmpush1.f
! mpmovin1 reorder and copy particles to ordered array
!          calls PPMOVIN1L
! mpcopyout1 copy ordered particles to linear array
!            calls PPCOPYOUT1
! mpcopyin1 copy unordered particles from linear array to ordered array
!           calls PPCOPYIN1
! mcheck1 verify particles are all in correct tiles
!         calls PPCHECK1L
! mpush11 push particles
!         calls GPPUSH1L
! mpushf11 push particles and determine which particles are leaving tile
!          calls GPPUSHF1L
! mrpush11 push relativistic particles with 1d electrostatic fields
!          calls GRPPUSH1L
! mrpushf11 push relativistic particles and determine which particles
!           are leaving tile
!           calls GRPPUSHF1L
! mpush1zf push particles in 1d with fixed velocities
!          calls PPUSH1ZF
! mpushf1zf push particles in 1d with fixed velocities and determines
!           which particles are leaving tile
!           calls PPUSHF1ZF
! mrpush1zf push relativistic particles in 1d with fixed momenta
!           calls RPPUSH1ZF
! mrpushf1zf push relativistic particles in 1d with fixed momenta and 
!            determines which particles are leaving tile
!            calls RPPUSHF1ZF
! mpost1 deposits charge density
!        calls GPPOST1L
! wmpush1 generic procedure to push particles
!         calls mrpushf11, mpushf11, mrpush11, or mpush11
! wmpush1zf generic procedure to push particles with fixed velocity
!           calls mpush1zf, mpushf1zf, mrpush1zf, or mrpushf1zf
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 16, 2017
!
      use libmpush1_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpmovin1(part,ppart,kpic,mx,irc)
! order particles in part by tiles and copy to ppart
! store number of particles in each tile in kpic
      implicit none
      integer, intent(in) :: mx
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(inout) :: kpic
! local data
      integer :: idimp, nop, nppmx, mx1
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! call low level procedure
      call PPMOVIN1L(part,ppart,kpic,nppmx,idimp,nop,mx,mx1,irc)
! check for errors
      if (irc /= 0) write (*,*) 'mpmovin1 overflow error, irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcopyout1(part,ppart,kpic,np,irc)
! copy ordered particles in array ppart to unordered array part
      implicit none
      integer, intent(inout) :: np, irc
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nop, nppmx, mx1
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! call low level procedure
      call PPCOPYOUT1(part,ppart,kpic,np,nop,nppmx,idimp,mx1,irc)
! check for errors
      if (irc /= 0) write (*,*) 'mpcopyout1 overflow error, irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcopyin1(part,ppart,kpic,irc)
! copy unordered particles in array part to ordered array ppart
      implicit none
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nop, nppmx, mx1
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! call low level procedure
      call PPCOPYIN1(part,ppart,kpic,nop,nppmx,idimp,mx1,irc)
! check for errors
      if (irc /= 0) write (*,*) 'mpcopyin1 overflow error, irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcheck1(ppart,kpic,nx,mx,irc)
! perform a sanity check to make sure particles ordered by tiles are all
! within bounds.
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mx1
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! call low level procedure
      call PPCHECK1L(ppart,kpic,idimp,nppmx,nx,mx,mx1,irc)
! check error
      if (irc /= 0) write (*,*) 'mcheck1 error: irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpush11(ppart,fx,kpic,qbm,dt,ek,tpush,nx,mx,ipbc)
! push particles with 1d electrostatic fields
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: qbm, dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:), intent(in) :: fx
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fx,1)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GPPUSH1L(ppart,fx,kpic,qbm,dt,ek,idimp,nppmx,nx,mx,nxv,mx1,  &
     &ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpushf11(ppart,fx,kpic,ncl,ihole,qbm,dt,ek,tpush,nx,mx,&
     &irc)
! push particles with 1d electrostatic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:), intent(in) :: fx
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fx,1)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GPPUSHF1L(ppart,fx,kpic,ncl,ihole,qbm,dt,ek,idimp,nppmx,nx,mx&
     &,nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrpush11(ppart,fx,kpic,qbm,dt,ci,ek,tpush,nx,mx,ipbc)
! push relativistic particles with 1d electrostatic fields
      implicit none
      integer, intent(in) :: nx, mx, ipbc
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:), intent(in) :: fx
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fx,1)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRPPUSH1L(ppart,fx,kpic,qbm,dt,ci,ek,idimp,nppmx,nx,mx,nxv,  &
     &mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrpushf11(ppart,fx,kpic,ncl,ihole,qbm,dt,ci,ek,tpush,nx&
     &,mx,irc)
! push relativistic particles with 1d electrostatic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:), intent(in) :: fx
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, nxv, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fx,1)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRPPUSHF1L(ppart,fx,kpic,ncl,ihole,qbm,dt,ci,ek,idimp,nppmx, &
     &nx,mx,nxv,mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpush1zf(ppart,kpic,dt,ek,tpush,nx,ipbc)
! push particles in 1d with fixed velocities
      implicit none
      integer, intent(in) :: nx, ipbc
      real, intent(in) :: dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPUSH1ZF(ppart,kpic,dt,ek,idimp,nppmx,nx,mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpushf1zf(ppart,kpic,ncl,ihole,dt,ek,tpush,nx,mx,irc)
! push particles in 1d with fixed velocities
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPUSHF1ZF(ppart,kpic,ncl,ihole,dt,ek,idimp,nppmx,nx,mx,mx1,  &
     &ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrpush1zf(ppart,kpic,dt,ci,ek,tpush,nx,ipbc)
! push relativistic particles in 1d with fixed momenta
      implicit none
      integer, intent(in) :: nx, ipbc
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call RPPUSH1ZF(ppart,kpic,dt,ci,ek,idimp,nppmx,nx,mx1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrpushf1zf(ppart,kpic,ncl,ihole,dt,ci,ek,tpush,nx,mx,  &
     &irc)
! push relativistic particles in 1d with fixed momenta
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      integer, intent(inout) :: irc
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, mx1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call RPPUSHF1ZF(ppart,kpic,ncl,ihole,dt,ci,ek,idimp,nppmx,nx,mx,  &
     &mx1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpost1(ppart,q,kpic,qm,tdpost,mx)
! deposit charge
      implicit none
      integer, intent(in) :: mx
      real, intent(in) :: qm
      real, intent(inout) :: tdpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:), intent(inout) :: q
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(q,1)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GPPOST1L(ppart,q,kpic,qm,nppmx,idimp,mx,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpush1(ppart,fx,kpic,ncl,ihole,qbm,dt,ci,ek,tpush,nx, &
     &mx,ipbc,relativity,plist,irc)
! generic procedure to push particles
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, mx, ipbc, relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:), intent(in) :: fx
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, ek, ncl, ihole, irc
         if (relativity==1) then
            call mrpushf11(ppart,fx,kpic,ncl,ihole,qbm,dt,ci,ek,tpush,nx&
     &,mx,irc)
         else
            call mpushf11(ppart,fx,kpic,ncl,ihole,qbm,dt,ek,tpush,nx,mx,&
     &irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmpush1 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and ek
         if (relativity==1) then
            call mrpush11(ppart,fx,kpic,qbm,dt,ci,ek,tpush,nx,mx,ipbc)
         else
            call mpush11(ppart,fx,kpic,qbm,dt,ek,tpush,nx,mx,ipbc)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpush1zf(ppart,kpic,ncl,ihole,dt,ci,ek,tpush,nx,mx,   &
     &ipbc,relativity,plist,irc)
! generic procedure to push particles with fixed velocity
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, mx, ipbc, relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, ek, ncl, ihole, irc
         if (relativity==1) then
            call mrpushf1zf(ppart,kpic,ncl,ihole,dt,ci,ek,tpush,nx,mx,  &
     &irc)
         else
            call mpushf1zf(ppart,kpic,ncl,ihole,dt,ek,tpush,nx,mx,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmpush1zf overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and ek
         if (relativity==1) then
            call mrpush1zf(ppart,kpic,dt,ci,ek,tpush,nx,ipbc)
         else
            call mpush1zf(ppart,kpic,dt,ek,tpush,nx,ipbc)
         endif
      endif
      end subroutine
!
      end module
