!-----------------------------------------------------------------------
!
      module modmpinit2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpinit2.f
! mpdcomp2 determines spatial decomposition for uniform distribution
!          calls PDNICOMP2L
! mpudistr2 calculates initial particle co-ordinates with uniform
!           density
!           for 2d or 2-1/2d code
!           calls PUDISTR2
! mpvdistr2 calculates initial particle velocities with maxwellian
!           velocity with drift for 2d code
!           calls PVDISTR2
! mpvdistr2h calculates initial particle velocities with maxwellian
!            velocity with drift for 2-1/2d code
!            calls PVDISTR2H
! mpvrdistr2 calculates initial particle momentum with maxwell-juttner
!            distribution with drift for 2d code
!            calls PVRDISTR2
! mpvrdistr2h calculates initial particle momenta with maxwell-juttner
!             distribution with drift for 2-1/2d code
!             calls PVRDISTR2H
! mpdblkp2 finds the maximum number of particles in each tile
!          calls PPDBLKP2L
! wmpvdistr2 generic procedure to initialize particle velocities in 2d
!            calls mpvdistr2 or mpvrdistr2
! wmpvdistr2h generic procedure to initialize particle velocities in
!             2-1/2d
!             calls mpvdistr2h or mpvrdistr2h
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 14, 2017
!
      use libmpinit2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpdcomp2(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp)
! determines spatial decomposition for uniform density distribution
      implicit none
      integer, intent(in) :: ny, kstrt, nvp
      integer, intent(inout) :: nyp, noff, nypmx, nypmn
      real, dimension(:), intent(inout) :: edges
! local data
      integer :: idps
! extract dimensions
      idps = size(edges,1)
! call low level procedure
      call PDNICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps)
! check for errors
      if (nypmn < 1) then
         if (kstrt==1) then
            write (*,*) 'combination not supported nvp, ny =',nvp,ny
         endif
         call PPEXIT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdistr2(part,edges,npp,vtx,vty,vdx,vdy,npx,npy,nx,ny, &
     &ipbc,ierr)
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
      implicit none
      integer, intent(in) :: npx, npy, nx, ny, ipbc
      integer, intent(inout) :: npp, ierr
      real, intent(in) :: vtx, vty, vdx, vdy
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
! local data
      integer :: idimp, npmax, idps
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1)
! call low level procedure
      call PDISTR2(part,edges,npp,vtx,vty,vdx,vdy,npx,npy,nx,ny,idimp,  &
     &npmax,idps,ipbc,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdistr2h(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,  &
     &npy,nx,ny,ipbc,ierr)
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
      implicit none
      integer, intent(in) :: npx, npy, nx, ny, ipbc
      integer, intent(inout) :: npp, ierr
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
! local data
      integer :: idimp, npmax, idps
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1)
! call low level procedure
      call PDISTR2H(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,nx,ny&
     &,idimp,npmax,idps,ipbc,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpudistr2(part,edges,npp,npx,npy,nx,ny,kstrt,ipbc,irc)
! calculates initial particle co-ordinates in 2d or 2-1/2d
! with uniform density
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: npx, npy, nx, ny, kstrt, ipbc
      integer, intent(inout) :: npp, irc
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
! local data
      integer :: idimp, npmax, idps
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1)
! call low level procedure
      call PUDISTR2(part,edges,npp,npx,npy,nx,ny,idimp,npmax,idps,ipbc, &
     &irc)
      if (irc > 0) then
         if (kstrt==1) write (*,*) 'mupdistr2:buffer overlow, irc=', irc
      else if (irc < 0) then
         if (kstrt==1) write (*,*) 'mpudistr2:particle number error'
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvdistr2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,kstrt,  &
     &irc)
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vdx, vdy
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVDISTR2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,idimp,npmax,   &
     &irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvdistr2:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy&
     &,kstrt,irc)
! calculates initial particle velocities in 2-1/2d
! with maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,&
     &npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvdistr2h:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvrdistr2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,    &
     &kstrt,irc)
! calculates initial particle momenta in 2d
! with maxwell-juttner distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vdx, vdy, ci
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVRDISTR2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,idimp,npmax&
     &,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvrdistr2:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvrdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx&
     &,npy,kstrt,irc)
! calculates initial particle momenta in 2-1/2d
! with maxwell-juttner distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVRDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,npy,  &
     &idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvrdistr2h:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdblkp2(part,kpic,npp,noff,nppmx,mx,my,mx1,irc)
! finds the maximum number of particles in each tile
      implicit none
      integer, intent(in) :: mx, my, mx1, npp
      integer, intent(in) :: noff
      integer, intent(inout) :: nppmx, irc
      real, dimension(:,:), intent(in) :: part
      integer, dimension(:), intent(inout) :: kpic
! local data
      integer :: idimp, npmax, mxyp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      mxyp1 = size(kpic,1)
! call low level procedure
      call PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mx1,    &
     &mxyp1,irc)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpdblkp2 error, irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpvdistr2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,    &
     &kstrt,relativity,irc)
! generic procedure to initialize particle velocities in 2d
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, relativity
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vdx, vdy, ci
      real, dimension(:,:), intent(inout) :: part
! maxwell-juttner distribution
      if (relativity==1) then
         call mpvrdistr2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,kstrt,  &
     &irc)
! maxwellian distribution
      else
         call mpvdistr2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,kstrt,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpvdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx&
     &,npy,kstrt,relativity,irc)
! generic procedure to initialize particle velocities in 2-1/2d
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, relativity
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! maxwell-juttner distribution
      if (relativity==1) then
         call mpvrdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,  &
     &npy,kstrt,irc)
! maxwellian distribution
      else
         call mpvdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,   &
     &kstrt,irc)
      endif
      end subroutine
!
      end module
