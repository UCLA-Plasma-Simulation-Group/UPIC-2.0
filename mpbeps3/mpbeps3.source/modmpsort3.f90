!-----------------------------------------------------------------------
!
      module msort3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpsort3.f
! mporder3a creates list of particles which are leaving tile, and
!           buffers outgoing particles
!           calls PPPORDER32LA or VPPPORDER32LA
! mporderf3a buffers outgoing particles
!            calls PPPORDERF32LA or VPPPORDERF32LA
! mporder3b copies buffers into particle array
!           calls PPPORDER32LB or VPPPORDER32LB
! mprsncl3 restores initial values of address offset ncl
!          calls PPPRSNCL3L
! mporderf3af performs final section of first part of particle sort
!             PPPORDERF32LAF
! mprstor3 restores particle coordinates from ppbuff
!          calls PPPRSTOR3L
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: july 5, 2018
!
      use libmpsort3_h
      use libvmpsort3_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mporder3a(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,&
     &nclr,noff,nyzp,tsort,kstrt,nx,ny,nz,mx,my,mz,mx1,myp1,popt,irc2)
! performs first part of particle reordering into tiles
      implicit none
      integer, intent(in) :: kstrt, nx, ny, nz, mx, my, mz, mx1, myp1
      integer, intent(in) :: popt
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      real, dimension(:,:,:), intent(inout) :: sbufl, sbufr
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:,:,:,:), intent(inout) :: ncll, nclr
      integer, dimension(:), intent(in) :: noff, nyzp
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, nppmx, npbmx, nbmax, ntmax, mzp1, mxzyp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(sbufl,2)
      ntmax = size(ihole,2) - 1
      mzp1 = size(kpic,1)/(mx1*myp1)
      mxzyp1 = size(ncll,2)
      idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector sort
      case (2)
         call VPPPORDER32LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll&
     &,nclr,noff,nyzp,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mzp1,mxzyp1&
     &,npbmx,ntmax,nbmax,idds,irc2)
! standard sort
      case default
         call PPPORDER32LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,&
     &nclr,noff,nyzp,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mzp1,mxzyp1,&
     &npbmx,ntmax,nbmax,idds,irc2)
      end select
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      if (irc2(1) /= 0) then
         write (*,*) kstrt, 'info:mporder3a overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mporderf3a(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr&
     &,tsort,kstrt,mx1,myp1,popt,irc2)
! performs first part of particle reordering into tiles,
! does not create list of particles which are leaving tile
      implicit none
      integer, intent(in) :: kstrt, mx1, myp1, popt
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      real, dimension(:,:,:), intent(inout) :: sbufl, sbufr
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
      integer, dimension(:,:,:,:), intent(inout) :: ncll, nclr
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, nppmx, npbmx, nbmax, ntmax, mzp1, mxzyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(sbufl,2)
      ntmax = size(ihole,2) - 1; mzp1 = size(ihole,3)/(mx1*myp1)
      mxzyp1 = size(ncll,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector sort
      case (2)
         call VPPPORDERF32LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,   &
     &nclr,idimp,nppmx,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax,nbmax,irc2)
! standard sort
      case default
         call PPPORDERF32LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr&
     &,idimp,nppmx,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax,nbmax,irc2)
      end select
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      if (irc2(1) /= 0) then
         write (*,*) kstrt, 'info:mporderf3a overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mporder3b(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,&
     &mclr,mcls,tsort,kstrt,nx,ny,nz,myp1,popt,irc2)
! performs second part of particle reordering into tiles
      implicit none
      integer, intent(in) :: kstrt, nx, ny, nz, myp1, popt
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: ppbuff
      real, dimension(:,:,:), intent(in) :: rbufl, rbufr
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(in) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
      integer, dimension(:,:,:,:), intent(in) :: mcll, mclr
      integer, dimension(:,:,:), intent(in) :: mcls
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, nppmx, npbmx, nbmax, ntmax, mx1, mzp1, mxzyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(rbufl,2)
      ntmax = size(ihole,2) - 1
      mx1 = size(mcls,2) - 1; mzp1 = size(kpic,1)/(mx1*myp1)
      mxzyp1 = size(mcll,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(popt)
! vector sort
      case (2)
         call VPPPORDER32LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll&
     &,mclr,mcls,idimp,nppmx,nx,ny,nz,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax, &
     &nbmax,irc2)
! standard sort
      case default
         call PPPORDER32LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,&
     &mclr,mcls,idimp,nppmx,nx,ny,nz,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax,  &
     &nbmax,irc2)
      end select
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      if (irc2(1) /= 0) then
         write (*,*) kstrt, 'info:mporder3b overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprsncl3(ncl,tsort)
! restores initial values of address offset ncl
      implicit none
      real, intent(inout) :: tsort
      integer, dimension(:,:), intent(inout) :: ncl
! local data
      integer :: mxzyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      mxzyp1 = size(ncl,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPRSNCL3L(ncl,mxzyp1)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mporderf3af(ppbuff,sbufl,sbufr,ncl,ncll,nclr,tsort,    &
     &kstrt,mx1,myp1,mzp1,irc2)
! performs final section of first part of particle sort
      implicit none
      integer, intent(in) :: kstrt, mx1, myp1, mzp1
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppbuff
      real, dimension(:,:,:), intent(inout) :: sbufl, sbufr
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:,:), intent(inout) :: ncll, nclr
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: idimp, npbmx, nbmax, mxzyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppbuff,1); npbmx = size(ppbuff,2); 
      nbmax = size(sbufl,2); mxzyp1 = size(ncll,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPORDERF32LAF(ppbuff,sbufl,sbufr,ncl,ncll,nclr,idimp,mx1,   &
     &myp1,mzp1,mxzyp1,npbmx,nbmax,irc2)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      if (irc2(1) /= 0) then
         write (*,*) kstrt, 'info:mporderf3af overflow: irc2=', irc2
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprstor3(ppart,ppbuff,ncl,ihole,tsort)
! restores particle coordinates from ppbuff
      implicit none
      real, intent(inout) :: tsort
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: ppbuff
      integer, dimension(:,:), intent(in) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
! local data
      integer :: idimp, nppmx, mxyzp1, npbmx, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyzp1 = size(ppart,3)
      npbmx = size(ppbuff,2); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPRSTOR3L(ppart,ppbuff,ncl,ihole,idimp,nppmx,mxyzp1,npbmx,  &
     &ntmax)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
      end module
