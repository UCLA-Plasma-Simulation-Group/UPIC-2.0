!-----------------------------------------------------------------------
!
      module modmppush3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmppush3.f
! mpmovin3 reorder and copy particles to ordered array
!          calls PPPMOVIN3L
! mpcopyout3 copy ordered particles to unordered array
!            calls PPPCOPYOUT3
! mpcopyin3 copy unordered particles from linear array to ordered array
!           calls PPPCOPYIN3
! mpcheck3 verify particles are all in correct tiles
!          calls PPPCHECK3L
! mppush3 push particles
!         calls PPGPPUSH32L
! mppushf3 push particles and determine which particles are leaving tile
!          calls PPGPPUSHF32L
! mprpush3 push relativistic particles with 3d electrostatic fields
!          calls PPGRPPUSH32L
! mprpushf3 push relativistic particles and determine which particles
!           are leaving tile
!           calls PPGRPPUSHF32L
! mppush3zf push particles in 3d with fixed velocities
!           calls PPGPPUSH32ZF
! mppushf3zf push particles in 3d with fixed velocities and determines
!            which particles are leaving tile
!            calls PPGPPUSHF32ZF
! mprpush3zf push relativistic particles in 3d with fixed momenta
!            calls PPGRPPUSH32ZF
! mprpushf3zf push relativistic particles in 3d with fixed momenta and 
!             determines which particles are leaving tile
!             calls PPGRPPUSHF32ZF
! mppost3 deposits charge density
!         calls PPGPPOST32L
! wmppush3 generic procedure to push particles
!          calls mprpushf3, mppushf3, mprpush3, or mppush3
! wmppush3zf generic procedure to push particles with fixed velocity
!            calls mppush3zf, mppushf3zf, mprpush3zf, or mprpushf3zf
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 16, 2017
!
      use libmppush3_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine mpmovin3(part,ppart,kpic,npp,noff,mx,my,mz,mx1,myp1,irc&
     &)
! order particles in part by tiles and copy to ppart
! store number of particles in each tile in kpic
      implicit none
      integer, intent(in) :: npp, mx, my, mz, mx1, myp1
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, npmax, nppmx, mxyzp1, idds
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! call low level procedure
      call PPPMOVIN3L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax,mx,my, &
     &mz,mx1,myp1,mxyzp1,idds,irc)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpmovin3 overflow error, irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcopyout3(part,ppart,kpic,npp,irc)
! copy ordered particles in array ppart to unordered array part
! store total particle number in npp
      implicit none
      integer, intent(inout) :: npp, irc
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, npmax, nppmx, mxyzp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1)
! call low level procedure
      call PPPCOPYOUT3(part,ppart,kpic,npp,npmax,nppmx,idimp,mxyzp1,irc)
! check for errors
      if (irc /= 0) write (*,*) 'mpcopyout3 overflow error, irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcopyin3(part,ppart,kpic,irc)
! copy unordered particles in array part to ordered array ppart
      implicit none
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, npmax, nppmx, mxyzp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1)
! call low level procedure
      call PPPCOPYIN3(part,ppart,kpic,npmax,nppmx,idimp,mxyzp1,irc)
! check for errors
      if (irc /= 0) write (*,*) 'mpcopyin3 overflow error, irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcheck3(ppart,kpic,noff,nyzp,nx,mx,my,mz,mx1,myp1,irc)
! perform a sanity check to make sure particles ordered by tiles are all
! within bounds.
      implicit none
      integer, intent(in) :: nx, mx, my, mz, mx1, myp1
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, mzp1, idds
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mzp1 = size(kpic,1)/(mx1*myp1); idds = size(noff,1)
! call low level procedure
      call PPPCHECK3L(ppart,kpic,noff,nyzp,idimp,nppmx,nx,mx,my,mz,mx1, &
     &myp1,mzp1,idds,irc)
! check error
      if (irc /= 0) then
         write (*,*) 'mpcheck3 error: irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppush3(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,tpush,nx,ny&
     &,nz,mx,my,mz,mx1,myp1,ipbc)
! push particles with 3d electrostatic fields
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      real, intent(in) :: qbm, dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz
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
      call PPGPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,idimp,nppmx, &
     &nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppushf3(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt,ek,&
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
! push particles with 3d electrostatic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, ntmax, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      mxyzp1 = size(kpic,1); ntmax = size(ihole,2) - 1
      idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt,ek,  &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,    &
     &ntmax,idds,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpush3(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ci,ek,tpush, &
     &nx,ny,nz,mx,my,mz,mx1,myp1,ipbc)
! push relativistic particles with 3d electrostatic fields
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz
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
      call PPGRPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ci,ek,idimp,   &
     &nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpushf3(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt,ci&
     &,ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
! push relativistic particles with 3d electrostatic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, ntmax, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
      mxyzp1 = size(kpic,1); ntmax = size(ihole,2) - 1
      idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt,ci, &
     &ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1, &
     &ntmax,idds,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppush3zf(ppart,kpic,dt,ek,tpush,nx,ny,nz,ipbc)
! push particles in 3d with fixed velocities
      implicit none
      integer, intent(in) :: nx, ny, nz, ipbc
      real, intent(in) :: dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mxyzp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGPPUSH32ZF(ppart,kpic,dt,ek,idimp,nppmx,nx,ny,nz,mxyzp1,   &
     &ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppushf3zf(ppart,kpic,ncl,ihole,noff,nyzp,dt,ek,tpush, &
     &nx,ny,nz,mx,my,mz,mx1,myp1,irc)
! push particles in 3d with fixed velocities
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1
      integer, intent(inout) :: irc
      real, intent(in) :: dt
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, mxyzp1, ntmax, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1); ntmax = size(ihole,2) - 1
      idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PPGPPUSHF32ZF(ppart,kpic,ncl,ihole,noff,nyzp,dt,ek,idimp,    &
     &nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mxyzp1,ntmax,idds,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpush3zf(ppart,kpic,dt,ci,ek,tpush,nx,ny,nz,ipbc)
! push relativistic particles in 3d with fixed momenta
      implicit none
      integer, intent(in) :: nx, ny, nz, ipbc
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mxyzp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRPPUSH32ZF(ppart,kpic,dt,ci,ek,idimp,nppmx,nx,ny,nz,mxyzp1&
     &,ipbc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprpushf3zf(ppart,kpic,ncl,ihole,noff,nyzp,dt,ci,ek,   &
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
! push relativistic particles in 3d with fixed momenta
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1
      integer, intent(inout) :: irc
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, mxyzp1, ntmax, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1); ntmax = size(ihole,2) - 1
      idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGRPPUSHF32ZF(ppart,kpic,ncl,ihole,noff,nyzp,dt,ci,ek,idimp,&
     &nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mxyzp1,ntmax,idds,irc)
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppost3(ppart,q,kpic,noff,qm,tdpost,mx,my,mz,mx1,myp1)
! deposit charge
      implicit none
      integer, intent(in) :: mx, my, mz, mx1, myp1
      real, intent(in) :: qm
      real, intent(inout) :: tdpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: q
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, nppmx, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPGPPOST32L(ppart,q,kpic,noff,qm,nppmx,idimp,mx,my,mz,nxv,   &
     &nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! record time
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmppush3(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt,ci,&
     &ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
! generic procedure to push particles
! list = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      integer, intent(in) :: relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, wke, ncl, ihole, irc
         if (relativity==1) then
            call mprpushf3(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt,ci&
     &,ek,tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
         else
            call mppushf3(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt,ek,&
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmppush3 overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and wke
         if (relativity==1) then
            call mprpush3(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ci,ek,tpush, &
     &nx,ny,nz,mx,my,mz,mx1,myp1,ipbc)
         else
            call mppush3(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,tpush,nx,ny&
     &,nz,mx,my,mz,mx1,myp1,ipbc)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmppush3zf(ppart,kpic,ncl,ihole,noff,nyzp,dt,ci,ek,    &
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,ipbc,relativity,plist,irc)
! generic procedure to push particles with fixed velocity
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, nz, mx, my, mz, mx1, myp1, ipbc
      integer, intent(in) :: relativity
      integer, intent(inout) :: irc
      logical, intent(in) :: plist
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:), intent(in) :: noff, nyzp
! also calculate list of particles leaving tile
      if (plist) then
! updates ppart, ek, ncl, ihole, irc
         if (relativity==1) then
            call mprpushf3zf(ppart,kpic,ncl,ihole,noff,nyzp,dt,ci,ek,   &
     &tpush,nx,ny,nz,mx,my,mz,mx1,myp1,irc)
         else
            call mppushf3zf(ppart,kpic,ncl,ihole,noff,nyzp,dt,ek,tpush, &
     &nx,ny,nz,mx,my,mz,mx1,myp1,irc)
         endif
         if (irc /= 0) then
            write (*,*) 'info:wmppush3zf overflow: irc=', irc
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and ek
         if (relativity==1) then
            call mprpush3zf(ppart,kpic,dt,ci,ek,tpush,nx,ny,nz,ipbc)
         else
            call mppush3zf(ppart,kpic,dt,ek,tpush,nx,ny,nz,ipbc)
         endif
      endif
      end subroutine
!
      end module
