!-----------------------------------------------------------------------
!
      module minit2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpinit2.f
! mnextran2 skips over nextrand groups of random numbers in order to
!           initialize different random ensembles
!           calls NEXTRAN2
! mpdcomp2 determines spatial decomposition for uniform distribution
!          calls PDNICOMP2L
! mpudistr2 calculates initial particle co-ordinates with uniform
!           density for 2d code
!           calls PUDISTR2
! mfdistr2 calculates initial particle co-ordinates with various density
!          profiles for 2d code
!          calls PLDISTR2 or PFDISTR2
! mgfdistr2 calculates initial particle co-ordinates with general
!           distribution in space with limited x and y coordinate range
!           for 2d code
!           calls PGFDISTR2
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
! mpvbdistr2h calculates initial particle velocities in 2-1/2d for
!             magnetized plasma with maxwellian velocity with drift in
!             Bparallel direction and ring distribution in Bperp
!             calls PVBDISTR2H
! mpdblkp2 finds the maximum number of particles in each tile
!          calls PPDBLKP2L
! mpfedges2 finds new 1d partitions from initial analytic distribution
!           function
!           calls PFEDGES2
! mpgfedges2 finds new 1d partitions from initial analytic distribution
!            function with limited y coordinate range
!            calls PGFEDGES2
! mpfholes2 determines list of particles leaving this processor
!           calls PFHOLES2
! wmpvdistr2 generic procedure to initialize particle velocities in 2d
!            calls mpvdistr2 or mpvrdistr2
! wmpvdistr2h generic procedure to initialize particle velocities in
!             2-1/2d
!             calls mpvdistr2h or mpvrdistr2h
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: august 1, 2018
!
      use libmpinit2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mnextran2(nextrand,ndim,np)
! for 2d code, this subroutine skips over nextrand groups of random
! numbers in order to initialize different random ensembles
      implicit none
      integer, intent(in) :: nextrand, ndim, np
! call low level procedure
      call NEXTRAN2(nextrand,ndim,np)
      end subroutine
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
         if (kstrt==1) write (*,*) 'mupdistr2:buffer overflow, irc=',irc
      else if (irc < 0) then
         if (kstrt==1) write (*,*) 'mpudistr2:particle number error'
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfdistr2(part,npp,ampx,scalex,shiftx,ampy,scaley,     &
     &shifty,npx,npy,nx,ny,kstrt,nvp,ipbc,ndpro,irc)
! calculates initial particle co-ordinates in 2d
! with various density profiles
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
! ampdx/ampdy = amplitude of density compared to uniform in x/y
! scaledx/scaledy = scale length for spatial coordinate in x/y
! shiftdx/shiftdy = shift of spatial coordinate in x/y
! irc = (0,N) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: npx, npy, nx, ny, kstrt, nvp, ipbc, ndpro
      integer, intent(inout) :: npp, irc
      real, intent(in) :: ampx, scalex, shiftx, ampy, scaley, shifty
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax, ierr
      real :: zero
      double precision :: dmpx, sxi, dshiftx, dmpy, syi, dshifty
      double precision, external :: DSDISTR1, DGDISTR1, DHDISTR1
      double precision, external :: DEDISTR1
      idimp = size(part,1); npmax = size(part,2)
      irc = 0
      ierr = 0
      dmpx = dble(ampx)
      sxi = 0.0d0
      if (scalex /= 0.0) sxi = 1.0d0/dble(scalex)
      dshiftx = dble(shiftx)
      dmpy = dble(ampy)
      syi = 0.0d0
      if (scaley /= 0.0) syi = 1.0d0/dble(scaley)
      dshifty = dble(shifty)
      zero = 0.0
! call low level procedure
      select case(ndpro)
! uniform density
      case (0)
         call PLDISTR2(part,npp,zero,zero,npx,npy,nx,ny,kstrt,nvp,idimp,&
     &npmax,ipbc,ierr)
! linear density
      case (1)
         if ((ampx.le.2.0).and.(ampy.le.2.0)) then
            call PLDISTR2(part,npp,ampx,ampy,npx,npy,nx,ny,kstrt,nvp,   &
     &idimp,npmax,ipbc,ierr)
         else
            ierr = -4
         endif
! sinusoidal density
      case (2)
         if ((ampx.le.1.0).and.(ampy.le.1.0)) then
            call PFDISTR2(part,npp,DSDISTR1,dmpx,sxi,dshiftx,DSDISTR1,  &
     &dmpy,syi,dshifty,npx,npy,nx,ny,kstrt,nvp,idimp,npmax,ipbc,ierr)
         else
            ierr = -5
         endif
! gaussian density
      case (3)
         call PFDISTR2(part,npp,DGDISTR1,dmpx,sxi,dshiftx,DGDISTR1,dmpy,&
     &syi,dshifty,npx,npy,nx,ny,kstrt,nvp,idimp,npmax,ipbc,ierr)
! hyperbolic secant squared density
      case (4)
         call PFDISTR2(part,npp,DHDISTR1,dmpx,sxi,dshiftx,DHDISTR1,dmpy,&
     &syi,dshifty,npx,npy,nx,ny,kstrt,nvp,idimp,npmax,ipbc,ierr)
! exponential density
      case (5)
         call PFDISTR2(part,npp,DEDISTR1,dmpx,sxi,dshiftx,DEDISTR1,dmpy,&
     &syi,dshifty,npx,npy,nx,ny,kstrt,nvp,idimp,npmax,ipbc,ierr)
      case default
         ierr = -3
      end select
! check error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpfdistr2 error: ndpro, ierr=', ndpro, ierr
         endif
         irc = ierr
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgfdistr2(part,npp,ampx,scalex,shiftx,ampy,scaley,    &
     &shifty,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,ndpro,irc)
! calculates initial particle co-ordinates in 2d with general
! distribution in space with limited x and y coordinate range
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
! ampdx/ampdy = amplitude of density compared to uniform in x/y
! scaledx/scaledy = scale length for spatial coordinate in x/y
! shiftdx/shiftdy = shift of spatial coordinate in x/y
! irc = (0,N) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: npx, npy, nx, ny, kstrt, nvp, ndpro
      integer, intent(inout) :: npp, irc
      real, intent(in) :: ampx, scalex, shiftx, ampy, scaley, shifty
      real, intent(in) :: xmin, xmax, ymin, ymax
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax, ierr
      double precision :: zero
      double precision :: dmpx, sxi, dshiftx, dmpy, syi, dshifty
      double precision, external :: DLDISTR1, DSDISTR1, DGDISTR1
      double precision, external :: DHDISTR1, DEDISTR1
      idimp = size(part,1); npmax = size(part,2)
      irc = 0
      ierr = 0
      dmpx = dble(ampx)
      sxi = 0.0d0
      if (scalex /= 0.0) sxi = 1.0d0/dble(scalex)
      dshiftx = dble(shiftx)
      dmpy = dble(ampy)
      syi = 0.0d0
      if (scaley /= 0.0) syi = 1.0d0/dble(scaley)
      dshifty = dble(shifty)
      zero = 0.0d0
! call low level procedure
      select case(ndpro)
! uniform density
      case (0)
         call PGFDISTR2(part,npp,DLDISTR1,zero,zero,zero,DLDISTR1,zero, &
     &zero,zero,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,idimp,npmax,&
     &ierr)
! linear density
      case (1)
         if ((ampx.le.2.0).and.(ampy.le.2.0)) then
            sxi = 1.0d0/dble(nx); dshiftx = 0.5d0
            syi = 1.0d0/dble(ny); dshifty = 0.5d0
            call PGFDISTR2(part,npp,DLDISTR1,dmpx,sxi,dshiftx,DLDISTR1, &
     &dmpy,syi,dshifty,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,idimp&
     &,npmax,ierr)
         else
            ierr = -4
         endif
! sinusoidal density
      case (2)
         if ((ampx.le.1.0).and.(ampy.le.1.0)) then
            call PGFDISTR2(part,npp,DSDISTR1,dmpx,sxi,dshiftx,DSDISTR1, &
     &dmpy,syi,dshifty,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,idimp&
     &,npmax,ierr)
         else
            ierr = -5
         endif
! gaussian density
      case (3)
         call PGFDISTR2(part,npp,DGDISTR1,dmpx,sxi,dshiftx,DGDISTR1,dmpy&
     &,syi,dshifty,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,idimp,   &
     &npmax,ierr)
! hyperbolic secant squared density
      case (4)
         call PGFDISTR2(part,npp,DHDISTR1,dmpx,sxi,dshiftx,DHDISTR1,dmpy&
     &,syi,dshifty,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,idimp,   &
     &npmax,ierr)
! exponential density
      case (5)
         call PGFDISTR2(part,npp,DEDISTR1,dmpx,sxi,dshiftx,DEDISTR1,dmpy&
     &,syi,dshifty,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,idimp,   &
     &npmax,ierr)
      case default
         ierr = -3
      end select
! check error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpgfdistr2 error: ndpro, ierr=', ndpro, ierr
         endif
         irc = ierr
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvdistr2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,kstrt,  &
     &nvp,irc)
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vdx, vdy
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVDISTR2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,kstrt,nvp,idimp&
     &,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvdistr2:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy&
     &,kstrt,nvp,irc)
! calculates initial particle velocities in 2-1/2d
! with maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,kstrt,&
     &nvp,idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvdistr2h:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvrdistr2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,    &
     &kstrt,nvp,irc)
! calculates initial particle momenta in 2d
! with maxwell-juttner distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vdx, vdy, ci
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVRDISTR2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,kstrt,nvp, &
     &idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvrdistr2:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvrdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx&
     &,npy,kstrt,nvp,irc)
! calculates initial particle momenta in 2-1/2d
! with maxwell-juttner distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVRDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,npy,  &
     &kstrt,nvp,idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvrdistr2h:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvbdistr2h(part,nps,npp,vtr,vtz,vdr,vdz,omx,omy,omz,  &
     &npx,npy,kstrt,nvp,irc)
! calculates initial particle velocities in 2-1/2d
! for magnetized plasma with maxwellian velocity with drift in Bparallel
! direction and ring distribution in Bperp
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
      integer, intent(inout) :: irc
      real, intent(in) :: vtr, vtz, vdr, vdz, omx, omy, omz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
      irc = 0
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVBDISTR2H(part,nps,npp,vtr,vtz,vdr,vdz,omx,omy,omz,npx,npy, &
     &kstrt,nvp,idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvbdistr2h:particle number error, irc=', irc
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
      subroutine mpfedges2(edges,nyp,noff,ampy,scaley,shifty,nypmx,nypmn&
     &,ny,kstrt,nvp,ipbc,ndpro,ierr)
! finds new 1d partitions from initial analytic distribution function
      implicit none
      integer, intent(in) :: ny, kstrt, nvp, ipbc, ndpro
      integer, intent(inout) :: nyp, noff, nypmx, nypmn, ierr
      real, intent(in) :: ampy, scaley, shifty
      real, dimension(:), intent(inout) :: edges
! local data
      integer :: idps
      double precision :: dmpy, syi, dshifty, zero
      double precision, external :: DLDISTR1, DSDISTR1, DGDISTR1
      double precision, external :: DHDISTR1, DEDISTR1
      idps = size(edges,1)
      ierr = 0
      dmpy = dble(ampy)
      syi = 0.0d0
      if (scaley /= 0.0) syi = 1.0d0/dble(scaley)
      dshifty = dble(shifty)
      zero = 0.0d0
! call low level procedure
      select case(ndpro)
! uniform density
      case (0)
         call PFEDGES2(edges,nyp,noff,DLDISTR1,zero,zero,zero,nypmx,    &
     &nypmn,ny,kstrt,nvp,idps,ipbc)
! linear density
      case (1)
         if (ampy.le.2.0) then
            call PFEDGES2(edges,nyp,noff,DLDISTR1,dmpy,syi,dshifty,nypmx&
     &,nypmn,ny,kstrt,nvp,idps,ipbc)
         else
            ierr = -4
         endif
! sinusoidal density
      case (2)
         if (ampy.le.1.0) then
            syi = 1.0d0/dble(ny); dshifty = 0.5d0
            call PFEDGES2(edges,nyp,noff,DSDISTR1,dmpy,syi,dshifty,nypmx&
     &,nypmn,ny,kstrt,nvp,idps,ipbc)
         else
            ierr = -5
         endif
! gaussian density
      case (3)
         call PFEDGES2(edges,nyp,noff,DGDISTR1,dmpy,syi,dshifty,nypmx,  &
     &nypmn,ny,kstrt,nvp,idps,ipbc)
! hyperbolic secant squared density
      case (4)
         call PFEDGES2(edges,nyp,noff,DHDISTR1,dmpy,syi,dshifty,nypmx,  &
     &nypmn,ny,kstrt,nvp,idps,ipbc)
! exponential density
      case (5)
         call PFEDGES2(edges,nyp,noff,DEDISTR1,dmpy,syi,dshifty,nypmx,  &
     &nypmn,ny,kstrt,nvp,idps,ipbc)
      case default
         ierr = -3
      end select
! check error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpfedges2 error: ndpro, ierr=', ndpro, ierr
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgfedges2(edges,nyp,noff,ampy,scaley,shifty,ymin,ymax,&
     &nypmx,nypmn,ny,kstrt,nvp,ndpro,ierr)
! finds new 1d partitions from initial analytic distribution function
      implicit none
      integer, intent(in) :: ny, kstrt, nvp, ndpro
      integer, intent(inout) :: nyp, noff, nypmx, nypmn, ierr
      real, intent(in) :: ampy, scaley, shifty
      real, intent(in) :: ymin, ymax
      real, dimension(:), intent(inout) :: edges
! local data
      integer :: idps
      double precision :: dmpy, syi, dshifty, zero
      double precision, external :: DLDISTR1, DSDISTR1, DGDISTR1
      double precision, external :: DHDISTR1, DEDISTR1
      idps = size(edges,1)
      ierr = 0
      dmpy = dble(ampy)
      syi = 0.0d0
      if (scaley /= 0.0) syi = 1.0d0/dble(scaley)
      dshifty = dble(shifty)
      zero = 0.0d0
! call low level procedure
      select case(ndpro)
! uniform density
      case (0)
         call PGFEDGES2(edges,nyp,noff,DLDISTR1,zero,zero,zero,ymin,ymax&
     &,nypmx,nypmn,ny,kstrt,nvp,idps)
! linear density
      case (1)
         if (ampy.le.2.0) then
            syi = 1.0d0/dble(ny); dshifty = 0.5d0
            call PGFEDGES2(edges,nyp,noff,DLDISTR1,dmpy,syi,dshifty,ymin&
     &,ymax,nypmx,nypmn,ny,kstrt,nvp,idps)
         else
            ierr = -4
         endif
! sinusoidal density
      case (2)
         if (ampy.le.1.0) then
            call PGFEDGES2(edges,nyp,noff,DSDISTR1,dmpy,syi,dshifty,ymin&
     &,ymax,nypmx,nypmn,ny,kstrt,nvp,idps)
         else
            ierr = -5
         endif
! gaussian density
      case (3)
         call PGFEDGES2(edges,nyp,noff,DGDISTR1,dmpy,syi,dshifty,ymin,  &
     &ymax,nypmx,nypmn,ny,kstrt,nvp,idps)
! hyperbolic secant squared density
      case (4)
         call PGFEDGES2(edges,nyp,noff,DHDISTR1,dmpy,syi,dshifty,ymin,  &
     &ymax,nypmx,nypmn,ny,kstrt,nvp,idps)
! exponential density
      case (5)
         call PGFEDGES2(edges,nyp,noff,DEDISTR1,dmpy,syi,dshifty,ymin,  &
     &ymax,nypmx,nypmn,ny,kstrt,nvp,idps)
      case default
         ierr = -3
      end select
! check error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpgfedges2 error: ndpro, ierr=', ndpro, ierr
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfholes2(part,edges,npp,iholep,ndim,nc)
! determines list of particles which are leaving this node
      implicit none
      integer, intent(in) :: npp, ndim, nc
      real, dimension(:,:), intent(in) :: part
      real, dimension(:), intent(in) :: edges
      integer, dimension(:), intent(inout) :: iholep
! local data
      integer :: idimp, npmax, idps, ntmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1); ntmax = size(iholep,1) - 1
! call low level procedure
      if ((nc==1).or.((nc==2).and.(idimp > (ndim+2)))) then
         call PFHOLES2(part,edges,npp,iholep,ndim,nc,idimp,npmax,idps,  &
     &ntmax)
      else
         write (*,*) 'mpfholes2 error: nc, idimp=', nc, idimp
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpvdistr2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,    &
     &kstrt,nvp,relativity,irc)
! generic procedure to initialize particle velocities in 2d
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp, relativity
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vdx, vdy, ci
      real, dimension(:,:), intent(inout) :: part
! maxwell-juttner distribution
      if (relativity==1) then
         call mpvrdistr2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,kstrt, &
     &nvp,irc)
! maxwellian distribution
      else
         call mpvdistr2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,kstrt,nvp, &
     &irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpvdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx&
     &,npy,kstrt,nvp,relativity,irc)
! generic procedure to initialize particle velocities in 2-1/2d
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp, relativity
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! maxwell-juttner distribution
      if (relativity==1) then
         call mpvrdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,  &
     &npy,kstrt,nvp,irc)
! maxwellian distribution
      else
         call mpvdistr2h(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,  &
     &kstrt,nvp,irc)
      endif
      end subroutine
!
      end module
