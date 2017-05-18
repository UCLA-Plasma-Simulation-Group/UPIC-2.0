!-----------------------------------------------------------------------
!
      module modmppush2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmppush2.f
! mpmovin2 reorder and copy particles to ordered array
!          calls PPPMOVIN2L
! mpcopyout2 copy ordered particles to unordered array
!            calls PPPCOPYOUT2
! mpcopyin2 copy unordered particles from linear array to ordered array
!           calls PPPCOPYIN2
! mpcheck2 verify particles are all in correct tiles
!          calls PPPCHECK2L
! mppush2 push particles
!         calls PPGPPUSH2L
! mppushf2 push particles and determine which particles are leaving tile
!          calls PPGPPUSHF2L
! mprpush2 push relativistic particles with 2d electrostatic fields
!          calls PPGRPPUSH2L
! mprpushf2 push relativistic particles and determine which particles
!           are leaving tile
!           calls PPGRPPUSHF2L
! mppush2zf push particles in 2d with fixed velocities
!           calls PPGPPUSH2ZF
! mppushf2zf push particles in 2d with fixed velocities and determines
!            which particles are leaving tile
!            calls PPGPPUSHF2ZF
! mprpush2zf push relativistic particles in 2d with fixed momenta
!            calls PPGRPPUSH2ZF
! mprpushf2zf push relativistic particles in 2d with fixed momenta and 
!             determines which particles are leaving tile
!             calls PPGRPPUSHF2ZF
! mppost2 deposits charge density
!         calls PPGPPOST2L
! wmppush2 generic procedure to push particles
!          calls mprpushf2, mppushf2, mprpush2, or mppush2
! wmppush2zf generic procedure to push particles with fixed velocity
!            calls mppush2zf, mppushf2zf, mprpush2zf, or mprpushf2zf
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 16, 2017
!
      use libmppush2_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpmovin2(part,ppart,kpic,npp,noff,mx,my,mx1,irc)
! order particles in part by tiles and copy to ppart
! store number of particles in each tile in kpic
      implicit none
      integer, intent(in) :: mx, my, mx1
      integer, intent(in) :: npp, noff
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(inout) :: kpic
! local data
      integer :: idimp, npmax, nppmx, mxyp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyp1 = size(kpic,1)
! call low level procedure
      call PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax,mx,my, &
     &mx1,mxyp1,irc)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpmovin2 overflow error, irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcopyout2(part,ppart,kpic,npp,irc)
! copy ordered particles in array ppart to unordered array part
! store total particle number in npp
      implicit none
      integer, intent(inout) :: npp, irc
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, npmax, nppmx, mxyp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyp1 = size(kpic,1)
! call low level procedure
      call PPPCOPYOUT2(part,ppart,kpic,npp,npmax,nppmx,idimp,mxyp1,irc)
! check for errors
      if (irc /= 0) write (*,*) 'mpcopyout2 overflow error, irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcopyin2(part,ppart,kpic,irc)
! copy unordered particles in array part to ordered array ppart
      implicit none
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, npmax, nppmx, mxyp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyp1 = size(kpic,1)
! call low level procedure
      call PPPCOPYIN2(part,ppart,kpic,npmax,nppmx,idimp,mxyp1,irc)
! check for errors
      if (irc /= 0) write (*,*) 'mpcopyin2 overflow error, irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,irc)
! perform a sanity check to make sure particles ordered by tiles are all
! within bounds.
      implicit none
      integer, intent(in) :: nx, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, myp1
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      myp1 = size(kpic,1)/mx1
! call low level procedure
      call PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,mx1,myp1,&
     &irc)
! check error
      if (irc /= 0) then
         write (*,*) 'mpcheck2 error: irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppush2(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,tpush,nx,ny, &
     &mx,my,mx1,ipbc)
! push particles with 2d electrostatic fields
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc
      integer, intent(in):: noff, nyp
      real, intent(in) :: qbm, dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, ndim, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      ndim = size(fxy,1); nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call PPGPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,nx,ny,mx,my, &
     &idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
      case default
         write (*,*) 'mppush2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppushf2(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ek,  &
     &tpush,nx,ny,mx,my,mx1,irc)
! push particles with 2d electrostatic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, ndim, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      ndim = size(fxy,1); nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call PPGPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ek,nx&
     &,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
      case default
         write (*,*) 'mppushf2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpush2(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,tpush,nx,&
     &ny,mx,my,mx1,ipbc)
! push relativistic particles with 2d electrostatic fields
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc
      integer, intent(in):: noff, nyp
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, ndim, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      ndim = size(fxy,1); nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call PPGRPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,nx,ny,mx,&
     &my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
      case default
         write (*,*) 'mprpush2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpushf2(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ci, &
     &ek,tpush,nx,ny,mx,my,mx1,irc)
! push relativistic particles with 2d electrostatic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, ndim, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      ndim = size(fxy,1); nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call PPGRPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ci, &
     &ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
      case default
         write (*,*) 'mprpushf2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppush2zf(ppart,kpic,dt,ek,tpush,nx,ny,ipbc)
! push particles in 2d with fixed velocities
      implicit none
      integer, intent(in) :: nx, ny, ipbc
      real, intent(in) :: dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGPPUSH2ZF(ppart,kpic,dt,ek,nx,ny,idimp,nppmx,mxyp1,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppushf2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,tpush,nx&
     &,ny,mx,my,mx1,irc)
! push particles in 2d with fixed velocities
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
      call PPGPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,nx,ny,mx,my,&
     &idimp,nppmx,mx1,mxyp1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpush2zf(ppart,kpic,dt,ci,ek,tpush,nx,ny,ipbc)
! push relativistic particles in 2d with fixed momenta
      implicit none
      integer, intent(in) :: nx, ny, ipbc
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRPPUSH2ZF(ppart,kpic,dt,ci,ek,nx,ny,idimp,nppmx,mxyp1,ipbc&
     &)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpushf2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,    &
     &tpush,nx,ny,mx,my,mx1,irc)
! push relativistic particles in 2d with fixed momenta
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,nx,ny,mx&
     &,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppost2(ppart,q,kpic,noff,qm,tdpost,mx,my,mx1)
! deposit charge
      implicit none
      integer, intent(in) :: mx, my, mx1
      integer, intent(in) :: noff
      real, intent(in) :: qm
      real, intent(inout) :: tdpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: q
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(q,1); nypmx = size(q,2)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv,nypmx, &
     &mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmppush2(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ci,ek&
     &,tpush,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
! generic procedure to push particles
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, relativity
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, ek, ncl, ihole, irc
         if (relativity==1) then
            call mprpushf2(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ci, &
     &ek,tpush,nx,ny,mx,my,mx1,irc)
         else
            call mppushf2(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ek,  &
     &tpush,nx,ny,mx,my,mx1,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmppush2 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and ek
         if (relativity==1) then
            call mprpush2(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,tpush,nx,&
     &ny,mx,my,mx1,ipbc)
         else
            call mppush2(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,tpush,nx,ny, &
     &mx,my,mx1,ipbc)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmppush2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,tpush&
     &,nx,ny,mx,my,mx1,ipbc,relativity,plist,irc)
! generic procedure to push particles with fixed velocity
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, relativity
      integer, intent(in) :: noff, nyp
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
            call mprpushf2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,    &
     &tpush,nx,ny,mx,my,mx1,irc)
         else
            call mppushf2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,tpush,nx&
     &,ny,mx,my,mx1,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmppush2zf overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and ek
         if (relativity==1) then
            call mprpush2zf(ppart,kpic,dt,ci,ek,tpush,nx,ny,ipbc)
         else
            call mppush2zf(ppart,kpic,dt,ek,tpush,nx,ny,ipbc)
         endif
      endif
      end subroutine
!
      end module
