!-----------------------------------------------------------------------
!
      module mpush2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmppush2.f
! mpmovin2 reorder and copy particles to ordered array
!          calls PPPMOVIN2L
! mpmovin2p reorder and copy particles to ordered array designed for
!           NUMA first touch architectures
!           calls PPPMOVIN2LP
! mpcopyout2 copy ordered particles to unordered array
!            calls PPPCOPYOUT2
! mpcopyin2 copy unordered particles from linear array to ordered array
!           calls PPPCOPYIN2
! mpcheck2 verify particles are all in correct tiles
!          calls PPPCHECK2L
! mppush2 push particles
!         calls PPGPPUSH2L or VPPGPPUSH2L
! mppushf2 push particles and determine which particles are leaving tile
!          calls PPGPPUSHF2L or VPPGPPUSHF2L
! mprpush2 push relativistic particles with 2d electrostatic fields
!          calls PPGRPPUSH2L or VPPGRPPUSH2L
! mprpushf2 push relativistic particles and determine which particles
!           are leaving tile
!           calls PPGRPPUSHF2L or VPPGRPPUSHF2L
! mppush2zf push particles in 2d with fixed velocities
!           calls PPGPPUSH2ZF or VPPGPPUSH2ZF
! mppushf2zf push particles in 2d with fixed velocities and determines
!            which particles are leaving tile
!            calls PPGPPUSHF2ZF or VPPGPPUSHF2ZF
! mprpush2zf push relativistic particles in 2d with fixed momenta
!            calls PPGRPPUSH2ZF or VPPGRPPUSH2ZF
! mprpushf2zf push relativistic particles in 2d with fixed momenta and 
!             determines which particles are leaving tile
!             calls PPGRPPUSHF2ZF
! mppost2 deposits charge density
!         calls PPGPPOST2L or VPPGPPOST2L
! wmppush2 generic procedure to push particles
!          calls mprpushf2, mppushf2, mprpush2, or mppush2
! wmppush2zf generic procedure to push particles with fixed velocity
!            calls mppush2zf, mppushf2zf, mprpush2zf, or mprpushf2zf
! mpset_pszero2 zeros out charge density array.
!               calls SET_PSZERO2
! mpset_pvzero2 zeros out current density array.
!               calls SET_PVZERO2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: august 6, 2018
!
      use libmppush2_h
      use libvmppush2_h
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
      subroutine mpmovin2p(part,ppart,kpic,npp,noff,mx,my,mx1,irc)
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
      integer, dimension(:,:), allocatable :: kp
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyp1 = size(kpic,1)
      allocate(kp(nppmx,mxyp1))
! call low level procedure
      call PPPMOVIN2LP(part,ppart,kpic,kp,npp,noff,nppmx,idimp,npmax,mx,&
     &my,mx1,mxyp1,irc)
      deallocate(kp)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpmovin2p overflow error, irc=', irc
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
     &mx,my,mx1,ipbc,popt)
! push particles with 2d electrostatic fields
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt
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
         select case(popt)
! vector push
         case (2)
            call VPPGPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,nx,ny,mx,&
     &my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
! standard push
         case default
            call PPGPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,nx,ny,mx, &
     &my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
         end select
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
     &tpush,nx,ny,mx,my,mx1,popt,irc)
! push particles with 2d electrostatic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, popt
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
         select case(popt)
! vector push
         case (2)
            call VPPGPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt, &
     &ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! standard push
         case default
            call PPGPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ek&
     &,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         end select
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
     &ny,mx,my,mx1,ipbc,popt)
! push relativistic particles with 2d electrostatic fields
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt
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
         select case(popt)
! vector push
         case (2)
            call VPPGRPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,nx,ny&
     &,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
! standard push
         case default
            call PPGRPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,nx,ny,&
     &mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
         end select
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
     &ek,tpush,nx,ny,mx,my,mx1,popt,irc)
! push relativistic particles with 2d electrostatic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, popt
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
         select case(popt)
! vector push
         case (2)
            call VPPGRPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
     &ci,ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! standard push
         case default
            call PPGRPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt, &
     &ci,ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         end select
      case default
         write (*,*) 'mprpushf2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppush2zf(ppart,kpic,dt,ek,tpush,nx,ny,ipbc,popt)
! push particles in 2d with fixed velocities
      implicit none
      integer, intent(in) :: nx, ny, ipbc, popt
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
      select case(popt)
! vector push
      case (2)
         call VPPGPPUSH2ZF(ppart,kpic,dt,ek,nx,ny,idimp,nppmx,mxyp1,ipbc&
     &)
! standard push
      case default
         call PPGPPUSH2ZF(ppart,kpic,dt,ek,nx,ny,idimp,nppmx,mxyp1,ipbc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppushf2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,tpush,nx&
     &,ny,mx,my,mx1,popt,irc)
! push particles in 2d with fixed velocities
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, popt
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
      select case(popt)
! vector push
      case (2)
         call VPPGPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,nx,ny,mx&
     &,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
! standard push
      case default
         call PPGPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,nx,ny,mx,&
     &my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpush2zf(ppart,kpic,dt,ci,ek,tpush,nx,ny,ipbc,popt)
! push relativistic particles in 2d with fixed momenta
      implicit none
      integer, intent(in) :: nx, ny, ipbc, popt
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
      select case(popt)
! vector push
      case (2)
         call VPPGRPPUSH2ZF(ppart,kpic,dt,ci,ek,nx,ny,idimp,nppmx,mxyp1,&
     &ipbc)
! standard push
      case default
         call PPGRPPUSH2ZF(ppart,kpic,dt,ci,ek,nx,ny,idimp,nppmx,mxyp1, &
     &ipbc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpushf2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,    &
     &tpush,nx,ny,mx,my,mx1,popt,irc)
! push relativistic particles in 2d with fixed momenta
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, popt, noff, nyp
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
      select case(popt)
! vector push
      case (2)
         call VPPGRPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,nx, &
     &ny,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
! standard push
      case default
         call PPGRPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,nx,ny&
     &,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
      end select
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppost2(ppart,q,kpic,noff,qm,tdpost,mx,my,mx1,popt)
! deposit charge
      implicit none
      integer, intent(in) :: mx, my, mx1, popt
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
      select case(popt)
! vector charge deposit
      case (2)
         call VPPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv,   &
     &nypmx,mx1,mxyp1)
! standard charge deposit
      case default
         call PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv,    &
     &nypmx,mx1,mxyp1)
      end select
! record time
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmppush2(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ci,ek&
     &,tpush,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
! generic procedure to push particles
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt, relativity
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
     &ek,tpush,nx,ny,mx,my,mx1,popt,irc)
         else
            call mppushf2(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ek,  &
     &tpush,nx,ny,mx,my,mx1,popt,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmppush2 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and ek
         if (relativity==1) then
            call mprpush2(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,tpush,nx,&
     &ny,mx,my,mx1,ipbc,popt)
         else
            call mppush2(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,tpush,nx,ny, &
     &mx,my,mx1,ipbc,popt)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmppush2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,tpush&
     &,nx,ny,mx,my,mx1,ipbc,popt,relativity,plist,irc)
! generic procedure to push particles with fixed velocity
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, popt, relativity
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
     &tpush,nx,ny,mx,my,mx1,popt,irc)
         else
            call mppushf2zf(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,tpush,nx&
     &,ny,mx,my,mx1,popt,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmppush2zf overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and ek
         if (relativity==1) then
            call mprpush2zf(ppart,kpic,dt,ci,ek,tpush,nx,ny,ipbc,popt)
         else
            call mppush2zf(ppart,kpic,dt,ek,tpush,nx,ny,ipbc,popt)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpset_pszero2(q,tdpost,mx,my,mx1,myp1)
! zeros out charge density array.
      implicit none
      integer, intent(in) :: mx, my, mx1, myp1
      real, intent(inout) :: tdpost
      real, dimension(:,:), intent(inout) :: q
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(q,1); nypmx = size(q,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call SET_PSZERO2(q,mx,my,nxv,nypmx,mx1,myp1)
! record time
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpset_pvzero2(cu,tdjpost,mx,my,mx1,myp1)
! zeros out current density array.
      implicit none
      integer, intent(in) :: mx, my, mx1, myp1
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: cu
! local data
      integer :: ndim, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nxv = size(cu,2); nypmx = size(cu,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call SET_PVZERO2(cu,mx,my,ndim,nxv,nypmx,mx1,myp1)
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
      end module
