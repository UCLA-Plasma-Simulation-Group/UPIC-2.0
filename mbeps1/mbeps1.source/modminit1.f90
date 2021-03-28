!-----------------------------------------------------------------------
!
      module minit1
!
! Fortran90 wrappers to 1d OpenMP PIC library libminit1.f
! mnextran1 skips over nextrand groups of random numbers in order to
!           initialize different random ensembles
!           calls NEXTRAN1
! mudistr1 calculates initial particle co-ordinates with uniform density
!          for 1d or 1-2/2d code
!          calls UDISTR1
! mfdistr1 calculates initial particle co-ordinates with various density
!          profiles for 1d or 1-2/2d code
!          calls UDISTR1, LDISTR1, or FDISTR1
! mgfdistr1 calculates initial particle co-ordinates with various
!           density profiles for 1d or 1-2/2d code where co-ordinates
!           are restricted to xmin <= x < xmax
!           calls UDISTR1, LDISTR1, or FDISTR1
! mvdistr1 calculates initial particle velocities with maxwellian
!          velocity with drift for 1d code
!          calls VDISTR1
! mvdistr1h calculates initial particle velocities with maxwellian
!           velocity with drift for 1-2/2d code
!           calls VDISTR1H
! mvrdistr1 calculates initial particle momentum with maxwell-juttner
!           distribution with drift for 1d code
!           calls VRDISTR1
! mvrdistr1h calculates initial particle momenta with maxwell-juttner
!            distribution with drift for 1-2/2d code
!            calls VRDISTR1H
! mwdistr1 calculates initial particle velocities with waterbag velocity
!          distribution with drift for 1d code
!          calls WDISTR1
! mwdistr1h calculates initial particle velocities with waterbag 
!           velocity distribution with drift for 1-2/2d code
!           calls WDISTR1H
! mvbdistr1h calculates initial particle velocities in 1-2/2d for
!            magnetized plasma with maxwellian velocity with drift in
!            Bparallel direction and ring distribution in Bperp
!            calls VBDISTR1H
! mdblkp1 finds the maximum number of particles in each tile
!         calls PPDBLKP2L
! wmvdistr1 generic procedure to initialize particle velocities in 1d
!           calls mvdistr1, mvrdistr1  or mwdistr1
! fprecision determines if default reals are actually doubles
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: december 1, 2017
!
      use libminit1_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mnextran1(nextrand,ndim,np)
! for 1d code, this subroutine skips over nextrand groups of random
! numbers in order to initialize different random ensembles
      implicit none
      integer, intent(in) :: nextrand, ndim, np
! call low level procedure
      call NEXTRAN1(nextrand,ndim,np)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mudistr1(part,nstart,npx,nx,ipbc,irc)
! calculates initial particle co-ordinates in 1d or 1-2/2d
! with uniform density
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx, nx, ipbc
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt
      irc = 0
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nt = nstart + npx - 1
! call low level procedure
      if (nt <= nop) then
         call UDISTR1(part,nstart,npx,idimp,nop,nx,ipbc)
      else
         write (*,*) 'mudistr1 overflow: nt, nop =', nt, nop
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfdistr1(part,ampx,scalex,shiftx,nstart,npx,nx,ipbc,   &
     &ndpro,irc)
! calculates initial particle co-ordinates in 1d
! with various density profiles
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
! ampdx = amplitude of density compared to uniform in x
! scaledx = scale length for spatial coordinate in x
! shiftdx = shift of spatial coordinate in x
! irc = (0,N) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx, nx, ipbc, ndpro
      integer, intent(inout) :: irc
      real, intent(in) :: ampx, scalex, shiftx
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt, ierr
      double precision :: dmpx, sxi, dshiftx
      double precision, external :: DSDISTR1, DGDISTR1, DHDISTR1
      double precision, external :: DEDISTR1
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      irc = 0
      ierr = 0
      nt = nstart + npx - 1
      dmpx = dble(ampx)
      sxi = 0.0d0
      if (scalex /= 0.0) sxi = 1.0d0/dble(scalex)
      dshiftx = dble(shiftx)
      if (nt > nop) then
         write (*,*) 'mfdistr1 overflow: nt, nop =', nt, nop
         irc = 1
      endif
! call low level procedure
      select case(ndpro)
! uniform density
      case (0)
         call UDISTR1(part,nstart,npx,idimp,nop,nx,ipbc)
! linear density
      case (1)
         call LDISTR1(part,ampx,nstart,npx,idimp,nop,nx,ipbc)
! sinusoidal density
      case (2)
         call FDISTR1(part,DSDISTR1,dmpx,sxi,dshiftx,nstart,npx,idimp,  &
     &nop,nx,ipbc,ierr)
! gaussian density
      case (3)
         call FDISTR1(part,DGDISTR1,dmpx,sxi,dshiftx,nstart,npx,idimp,  &
     &nop,nx,ipbc,ierr)
! hyperbolic secant squared density
      case (4)
         call FDISTR1(part,DHDISTR1,dmpx,sxi,dshiftx,nstart,npx,idimp,  &
     &nop,nx,ipbc,ierr)
! exponential density
      case (5)
         call FDISTR1(part,DEDISTR1,dmpx,sxi,dshiftx,nstart,npx,idimp,  &
     &nop,nx,ipbc,ierr)
      case default
         ierr = -3
      end select
! check error
      if (ierr /= 0) then
         write (*,*) 'mfdistr1 error: ndpro, ierr=', ndpro, ierr
         irc = ierr
      endif
      end subroutine 
!
!-----------------------------------------------------------------------
      subroutine mgfdistr1(part,ampx,scalex,shiftx,xmin,xmax,nstart,npx,&
     &nx,ndpro,irc)
! calculates initial particle co-ordinates in 1d
! with various density profiles
! where co-ordinates are restricted to xmin <= x < xmax
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
! ampdx = amplitude of density compared to uniform in x
! scaledx = scale length for spatial coordinate in x
! shiftdx = shift of spatial coordinate in x
! irc = (0,N) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx, nx, ndpro
      integer, intent(inout) :: irc
      real, intent(in) :: ampx, scalex, shiftx, xmin, xmax
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt, ierr
      double precision zero, dmpx, sxi, half, dshiftx
      double precision, external :: DLDISTR1, DSDISTR1, DGDISTR1
      double precision, external :: DHDISTR1, DEDISTR1
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      irc = 0
      ierr = 0
      nt = nstart + npx - 1
      dmpx = dble(ampx)
      sxi = 0.0d0
      if (scalex /= 0.0) sxi = 1.0d0/dble(scalex)
      dshiftx = dble(shiftx)
      if (nt > nop) then
         write (*,*) 'mgfdistr1 overflow: nt, nop =', nt, nop
         irc = 1
      endif
! call low level procedure
      select case(ndpro)
! uniform density
      case (0)
         zero = 0.0d0
         call GFDISTR1(part,DLDISTR1,zero,zero,zero,xmin,xmax,nstart,npx&
     &,idimp,nop,nx,ierr)
! linear density
      case (1)
         sxi = 1.0d0/dble(nx); half = 0.5d0
         call GFDISTR1(part,DLDISTR1,dmpx,sxi,half,xmin,xmax,nstart,npx,&
     &idimp,nop,nx,ierr)
! sinusoidal density
      case (2)
         call GFDISTR1(part,DSDISTR1,dmpx,sxi,dshiftx,xmin,xmax,nstart, &
     &npx,idimp,nop,nx,ierr)
! gaussian density
      case (3)
         call GFDISTR1(part,DGDISTR1,dmpx,sxi,dshiftx,xmin,xmax,nstart, &
     &npx,idimp,nop,nx,ierr)
! hyperbolic secant squared density
      case (4)
         call GFDISTR1(part,DHDISTR1,dmpx,sxi,dshiftx,xmin,xmax,nstart, &
     &npx,idimp,nop,nx,ierr)
! exponential density
      case (5)
         call GFDISTR1(part,DEDISTR1,dmpx,sxi,dshiftx,xmin,xmax,nstart, &
     &npx,idimp,nop,nx,ierr)
      case default
         ierr = -3
      end select
! check error
      if (ierr /= 0) then
         write (*,*) 'mgfdistr1 error: ndpro, ierr=', ndpro, ierr
         irc = ierr
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvdistr1(part,nstart,vtx,vdx,npx,irc)
! calculates initial particle velocities in 1d
! with maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vdx
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt
      irc = 0
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nt = nstart + npx - 1
! call low level procedure
      if (nt <= nop) then
         call VDISTR1(part,vtx,vdx,nstart,npx,idimp,nop)
      else
         write (*,*) 'mvdistr1 overflow: nt, nop =', nt, nop
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvdistr1h(part,nstart,vtx,vty,vtz,vdx,vdy,vdz,npx,irc)
! calculates initial particle velocities in 1-2/2d
! with maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt
      irc = 0
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nt = nstart + npx - 1
! call low level procedure
      if (nt <= nop) then
         call VDISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,nstart,npx,idimp,nop&
     &)
      else
         write (*,*) 'mvdistr1h overflow: nt, nop =', nt, nop
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvrdistr1(part,nstart,vtx,vdx,ci,npx,irc)
! calculates initial particle momenta in 1d
! with maxwell-juttner distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vdx, ci
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt
      irc = 0
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nt = nstart + npx - 1
! call low level procedure
      if (nt <= nop) then
         call VRDISTR1(part,vtx,vdx,ci,nstart,npx,idimp,nop)
      else
         write (*,*) 'mvrdistr1 overflow: nt, nop =', nt, nop
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvrdistr1h(part,nstart,vtx,vty,vtz,vdx,vdy,vdz,ci,npx, &
     &irc)
! calculates initial particle momenta in 1-2/2d
! with maxwell-juttner distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt
      irc = 0
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nt = nstart + npx - 1
! call low level procedure
      if (nt <= nop) then
         call VRDISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,ci,nstart,npx,idimp&
     &,nop)
      else
         write (*,*) 'mvrdistr1h overflow: nt, nop =', nt, nop
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwdistr1(part,nstart,vtx,vdx,npx,irc)
! calculates initial particle velocities in 1d
! with waterbag velocity distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vdx
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt
      irc = 0
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nt = nstart + npx - 1
! call low level procedure
      if (nt <= nop) then
         call WDISTR1(part,vtx,vdx,nstart,npx,idimp,nop)
      else
         write (*,*) 'mwdistr1 overflow: nt, nop =', nt, nop
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwdistr1h(part,nstart,vtx,vty,vtz,vdx,vdy,vdz,npx,irc)
! calculates initial particle velocities in 1-2/2d
! with waterbag velocity distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt
      irc = 0
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nt = nstart + npx - 1
! call low level procedure
      if (nt <= nop) then
         call WDISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,nstart,npx,idimp,nop&
     &)
      else
         write (*,*) 'mwdistr1h overflow: nt, nop =', nt, nop
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvbdistr1h(part,nstart,vtr,vtz,vdr,vdz,omx,omy,omz,npx,&
     &irc)
! calculates initial particle velocities in 1-2/2d
! for magnetized plasma with maxwellian velocity with drift in Bparallel
! direction and ring distribution in Bperp
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nstart, npx
      integer, intent(inout) :: irc
      real, intent(in) :: vtr, vtz, vdr, vdz, omx, omy, omz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop, nt
      irc = 0
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nt = nstart + npx - 1
! call low level procedure
      if (nt <= nop) then
         call VBDISTR1H(part,vtr,vtz,vdr,vdz,omx,omy,omz,nstart,npx,    &
     &idimp,nop)
      else
         write (*,*) 'mvbdistr1h overflow: nt, nop =', nt, nop
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mdblkp2(part,kpic,nppmx,np,mx,irc)
! finds the maximum number of particles in each tile
      implicit none
      integer, intent(in) :: np, mx
      integer, intent(inout) :: nppmx, irc
      real, dimension(:,:), intent(in) :: part
      integer, dimension(:), intent(inout) :: kpic
! local data
      integer :: idimp, nop, mx1
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      mx1 = size(kpic,1)
! call low level procedure
      call DBLKP1L(part,kpic,nppmx,idimp,np,nop,mx,mx1,irc)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mdblkp1 error, irc=', irc
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmvdistr1(part,nstart,vtx,vdx,ci,npx,nvdist,relativity,&
     &irc)
! generic procedure to initialize particle velocities in 1d
      implicit none
      integer, intent(in) :: nstart, npx, nvdist, relativity
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vdx, ci
      real, dimension(:,:), intent(inout) :: part
! equilibrium distribution
      if (nvdist==1) then
! maxwell-juttner distribution
         if (relativity==1) then
            call mvrdistr1(part,nstart,vtx,vdx,ci,npx,irc)
! maxwellian distribution
         else
            call mvdistr1(part,nstart,vtx,vdx,npx,irc)
         endif
! waterbag distribution
      else if (nvdist==2) then
          call mwdistr1(part,nstart,vtx,vdx,npx,irc)
      else
         write (*,*) 'wmvdistr1: invalid nvdist = ', nvdist
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmvdistr1h(part,nstart,vtx,vty,vtz,vdx,vdy,vdz,ci,npx, &
     &nvdist,relativity,irc)
! generic procedure to initialize particle velocities in 1-2/2d
      implicit none
      integer, intent(in) :: nstart, npx, nvdist, relativity
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! equilibrium distribution
      if (nvdist==1) then
! maxwell-juttner distribution
         if (relativity==1) then
            call mvrdistr1h(part,nstart,vtx,vty,vtz,vdx,vdy,vdz,ci,npx, &
     &irc)
! maxwellian distribution
         else
            call mvdistr1h(part,nstart,vtx,vty,vtz,vdx,vdy,vdz,npx,irc)
         endif
! waterbag distribution
      else if (nvdist==2) then
          call mwdistr1h(part,nstart,vtx,vty,vtz,vdx,vdy,vdz,npx,irc)
      else
         write (*,*) 'wmvdistr1h: invalid nvdist = ', nvdist
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      integer function fprecision()
! function determines if default reals are actually doubles
      implicit none
      real :: prec
! ndprec = (0,1) = (no,yes) = (normal,autodouble) precision used
      if (digits(prec) > 24) then
         fprecision = 1
      else
         fprecision = 0
      endif
      end function fprecision
!
      end module
