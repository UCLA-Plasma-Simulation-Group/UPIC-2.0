!-----------------------------------------------------------------------
!
      module msort1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmsort1.f
! morder1 creates list of particles which are leaving tile, buffers
!         outgoing particles, and copies buffers into particle array
!         calls PPORDER1L
! morderf1 buffers outgoing particles and copies buffers into particle
!          array
!          calls PPORDERF1L
! wmporder1 generic procedure to perform particle reordering into tiles
!           calls mporder1 or morderf1
! mprsncl1 restores initial values of address offset ncl
!          calls PPRSNCL1L
! mprstor1 restores particle coordinates from ppbuff
!          calls PPRSTOR1L
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: january 30, 2017
!
      use libmsort1_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,irc2)
! performs particle reordering into tiles
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, nppmx, npbmx, ntmax, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,mx,mx1, &
     &npbmx,ntmax,irc2)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine morderf1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,irc2)
! performs particle reordering into tiles,
! does not create list of particles which are leaving tile
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, nppmx, npbmx, ntmax, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,mx,mx1,&
     &npbmx,ntmax,irc2)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,plist&
     &,irc2)
! generic procedure to perform particle reordering into tiles
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(inout) :: tsort
      logical, intent(in) :: plist
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(2), intent(inout) :: irc2
! do not calculate list of particles leaving tile
      if (plist) then
         call morderf1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,irc2)
! calculate list of particles leaving tile
      else
         call mporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,irc2)
      endif
      if (irc2(1) /= 0) then
         write (*,*) 'info:wmporder1 overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprsncl1(ncl,tsort)
! restores initial values of address offset ncl
      implicit none
      real, intent(inout) :: tsort
      integer, dimension(:,:), intent(inout) :: ncl
! local data
      integer :: mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      mx1 = size(ncl,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRSNCL1L(ncl,mx1)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprstor1(ppart,ppbuff,ncl,ihole,tsort)
! restores particle coordinates from ppbuff
      implicit none
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: ppbuff
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
! local data
      integer :: idimp, nppmx, npbmx, ntmax, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2); mx1 = size(ppart,3)
      npbmx = size(ppbuff,2); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRSTOR1L(ppart,ppbuff,ncl,ihole,idimp,nppmx,mx1,npbmx,ntmax)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
      end module
