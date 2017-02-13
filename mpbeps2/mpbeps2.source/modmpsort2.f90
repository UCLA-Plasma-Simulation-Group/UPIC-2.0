!-----------------------------------------------------------------------
!
      module modmpsort2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpsort2.f
! mporder2a creates list of particles which are leaving tile, and
!           buffers outgoing particles
!           calls PPPORDER2LA
! mporderf2a buffers outgoing particles
!            calls PPPORDERF2LA
! mporder2b copies buffers into particle array
!           calls PPPORDER2LB
! mprsncl2 restores initial values of address offset ncl
!          calls PPPRSNCL2L
! mporderf2af performs final section of first part of particle sort
!             PPPORDERF2LAF
! mprstor2 restores particle coordinates from ppbuff
!          calls PPPRSTOR2L
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: january 30, 2017
!
      use libmpsort2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mporder2a(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,&
     &nclr,noff,nyp,tsort,kstrt,nx,ny,mx,my,irc2)
! performs first part of particle reordering into tiles
      implicit none
      integer, intent(in) :: kstrt, nx, ny, mx, my
      integer, intent(in) :: noff, nyp
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      real, dimension(:,:), intent(inout) :: sbufl, sbufr
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:,:), intent(inout) :: ncll, nclr
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, nppmx, npbmx, nbmax, ntmax, mx1, myp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(sbufl,2)
      ntmax = size(ihole,2) - 1
      mx1 = size(ncll,2); myp1 = size(kpic,1)/mx1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,nclr&
     &,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,nbmax,irc2)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      if (irc2(1) /= 0) then
         write (*,*) kstrt, 'info:mporder2a overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mporderf2a(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr&
     &,tsort,kstrt,irc2)
! performs first part of particle reordering into tiles,
! does not create list of particles which are leaving tile
      implicit none
      integer, intent(in) :: kstrt
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      real, dimension(:,:), intent(inout) :: sbufl, sbufr
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
      integer, dimension(:,:), intent(inout) :: ncll, nclr
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, nppmx, npbmx, nbmax, mx1, ntmax, myp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(sbufl,2)
      mx1 = size(ncll,2)
      ntmax = size(ihole,2) - 1; myp1 = size(ihole,3)/mx1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr,   &
     &idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc2)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      if (irc2(1) /= 0) then
         write (*,*) kstrt, 'info:mporderf2a overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mporder2b(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,&
     &mclr,tsort,kstrt,nx,ny,irc2)
! performs second part of particle reordering into tiles
      implicit none
      integer, intent(in) :: kstrt, nx, ny
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: ppbuff
      real, dimension(:,:), intent(in) :: rbufl, rbufr
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(in) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
      integer, dimension(:,:), intent(in) :: mcll, mclr
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, nppmx, npbmx, nbmax, mx1, ntmax, myp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(rbufl,2)
      mx1 = size(mcll,2)
      ntmax = size(ihole,2) - 1; myp1 = size(kpic,1)/mx1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr&
     &,idimp,nppmx,nx,ny,mx1,myp1,npbmx,ntmax,nbmax,irc2)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      if (irc2(1) /= 0) then
         write (*,*) kstrt, 'info:mporder2b overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprsncl2(ncl,tsort)
! restores initial values of address offset ncl
      implicit none
      real, intent(inout) :: tsort
      integer, dimension(:,:), intent(inout) :: ncl
! local data
      integer :: mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      mxyp1 = size(ncl,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPRSNCL2L(ncl,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mporderf2af(ppbuff,sbufl,sbufr,ncl,ncll,nclr,tsort,    &
     &kstrt,irc2)
! performs final section of first part of particle reordering into tiles
! buffers outgoing particles
      implicit none
      integer, intent(in) :: kstrt
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppbuff
      real, dimension(:,:), intent(inout) :: sbufl, sbufr
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:), intent(inout) :: ncll, nclr
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, npbmx, nbmax, mx1, myp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppbuff,1); npbmx = size(ppbuff,2)
      nbmax = size(sbufl,2); mx1 = size(ncll,2); myp1 = size(ncl,2)/mx1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPORDERF2LAF(ppbuff,sbufl,sbufr,ncl,ncll,nclr,idimp,mx1,myp1&
     &,npbmx,nbmax,irc2)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      if (irc2(1) /= 0) then
         write (*,*) kstrt, 'info:mporderf2af overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprstor2(ppart,ppbuff,ncl,ihole,tsort)
! restores particle coordinates from ppbuff
      implicit none
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: ppbuff
      integer, dimension(:,:), intent(in) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
! local data
      integer :: idimp, nppmx, mxyp1, npbmx, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(ppart,3)
      npbmx = size(ppbuff,2); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPRSTOR2L(ppart,ppbuff,ncl,ihole,idimp,nppmx,mxyp1,npbmx,   &
     &ntmax)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
      end module
