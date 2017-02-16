!-----------------------------------------------------------------------
!
      module modmpinit3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpinit3.f
! mpdcomp3 determines spatial decomposition for uniform distribution
!          calls PDNICOMP32L
! mpdistr3 calculates initial particle co-ordinates and velocities
!          with uniform density and maxwellian velocity with drift
!          for 3d code
!          calls PDISTR32
! mprdistr3 calculates initial particle co-ordinates and momenta with
!           uniform density and maxwell-juttner momentum with drift
!           for 3d code with relativistic particles
!           calls PRDISTR32
! mpudistr3 calculates initial particle co-ordinates with uniform
!           density for 3d code
!           calls PUDISTR32
! mpvdistr3 calculates initial particle velocities with maxwellian
!           velocity with drift for 3d code
!           calls PVDISTR32
! mpvrdistr3 calculates initial particle momenta with maxwell-juttner
!            distribution with drift for 3d code
!            calls PVRDISTR32
! mpdblkp3 finds the maximum number of particles in each tile
!          calls PPDBLKP3L
! wmpdistr3 generic procedure to initialize particle co-ordinates and
!           velocities/momenta with uniform density and
!           maxwellian/juttner distributions with drift for 3d code
!           calls mpdistr3 or mprdistr3
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 15, 2017
!
      use libmpinit3_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpdcomp3(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny,nz,&
     &kstrt,nvpy,nvpz)
! determines spatial decomposition for uniform density distribution
      implicit none
      integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz
      integer, intent(inout) :: nypmx, nzpmx, nypmn, nzpmn
      integer, dimension(:), intent(inout) :: nyzp, noff
      real, dimension(:), intent(inout) :: edges
! local data
      integer :: idps, idds
! extract dimensions
      idps = size(edges,1); idds = size(noff,1)
! call low level procedure
      call PDNICOMP32L(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny,nz,   &
     &kstrt,nvpy,nvpz,idps,idds)
! check for errors
      if (kstrt==1) then
         if (nypmn < 1) then
            write (*,*) 'combination not supported nvpy, ny =',nvpy,ny
         endif
         if (nzpmn < 1) then
            write (*,*) 'combination not supported nvpz, nz =',nvpz,nz
         endif
      endif
      if ((nypmn < 1).or.(nzpmn < 1)) then
         call PPEXIT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdistr3(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy&
     &,npz,nx,ny,nz,kstrt,ipbc,irc)
! calculates initial particle co-ordinates and velocities in 3d
! with uniform density and maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: npx, npy, npz, nx, ny, nz, kstrt, ipbc
      integer, intent(inout) :: npp, irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
! local data
      integer :: idimp, npmax, idps
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1)
! call low level procedure
      call PDISTR32(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz, &
     &nx,ny,nz,idimp,npmax,idps,ipbc,irc)
      if (irc > 0) then
         if (kstrt==1) write (*,*) 'mpdistr3:buffer overflow, irc=', irc
      else if (irc < 0) then
         if (kstrt==1) write (*,*) 'mpdistr3:particle number error'
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdistr3(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx&
     &,npy,npz,nx,ny,nz,kstrt,ipbc,irc)
! calculates initial particle co-ordinates and momenta with uniform
! density and maxwell-juttner momentum with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: npx, npy, npz, nx, ny, nz, kstrt, ipbc
      integer, intent(inout) :: npp, irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
! local data
      integer :: idimp, npmax, idps
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1)
! call low level procedure
      call PRDISTR32(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,npy, &
     &npz,nx,ny,nz,idimp,npmax,idps,ipbc,irc)
      if (irc > 0) then
         if (kstrt==1) write (*,*) 'mprdistr3:buffer overflow, irc=',irc
      else if (irc < 0) then
         if (kstrt==1) write (*,*) 'mprdistr3:particle number error'
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpudistr3(part,edges,npp,npx,npy,npz,nx,ny,nz,kstrt,   &
     &ipbc,irc)
! calculates initial particle co-ordinates in 3d with uniform density
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: npx, npy, npz, nx, ny, nz, kstrt, ipbc
      integer, intent(inout) :: npp, irc
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
! local data
      integer :: idimp, npmax, idps
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1)
! call low level procedure
      call PUDISTR32(part,edges,npp,npx,npy,npz,nx,ny,nz,idimp,npmax,   &
     &idps,ipbc,irc)
      if (irc > 0) then
         if (kstrt==1) write (*,*) 'mupdistr3:buffer overflow, irc=',irc
      else if (irc < 0) then
         if (kstrt==1) write (*,*) 'mpudistr3:particle number error'
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvdistr3(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,&
     &npz,kstrt,irc)
! calculates initial particle velocities in 3d
! with maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, npz, kstrt
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,  &
     &idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvdistr3:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvrdistr3(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,&
     &npy,npz,kstrt,irc)
! calculates initial particle momenta in 3d
! with maxwell-juttner distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, npz, kstrt
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVRDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,npy,  &
     &npz,idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvrdistr3:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdblkp3(part,kpic,npp,noff,nppmx,mx,my,mz,mx1,myp1,irc&
     &)
! finds the maximum number of particles in each tile
      implicit none
      integer, intent(in) :: npp, mx, my, mz, mx1, myp1
      integer, intent(inout) :: nppmx, irc
      real, dimension(:,:), intent(in) :: part
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, npmax, mxyzp1, idds
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! call low level procedure
      call PPDBLKP3L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mz,mx1, &
     &myp1,mxyzp1,idds,irc)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpdblkp3 error, irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpdistr3(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx&
     &,npy,npz,nx,ny,nz,kstrt,ipbc,relativity,irc)
! calculates initial particle co-ordinates and velocities/momenta with
! uniform density and maxwellian/maxwell-juttner momentum with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: npx, npy, npz, nx, ny, nz, kstrt, ipbc
      integer, intent(in) :: relativity
      integer, intent(inout) :: npp, irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
! maxwell-juttner distribution
      if (relativity==1) then
         call mprdistr3(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,  &
     &npy,npz,nx,ny,nz,kstrt,ipbc,irc)
! maxwellian distribution
      else
         call mpdistr3(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,  &
     &npz,nx,ny,nz,kstrt,ipbc,irc)
      endif
      end subroutine
!
      end module
