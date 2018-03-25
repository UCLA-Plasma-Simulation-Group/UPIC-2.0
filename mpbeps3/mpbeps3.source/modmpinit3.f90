!-----------------------------------------------------------------------
!
      module modmpinit3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpinit3.f
! mnextran3 skips over nextrand groups of random numbers in order to
!           initialize different random ensembles
!           calls NEXTRAN3
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
! mfdistr3 calculates initial particle co-ordinates with various density
!          profiles for 3d code
!          calls PLDISTR32 or PFDISTR32
! mpvdistr3 calculates initial particle velocities with maxwellian
!           velocity with drift for 3d code
!           calls PVDISTR32
! mpvrdistr3 calculates initial particle momenta with maxwell-juttner
!            distribution with drift for 3d code
!            calls PVRDISTR32
! mpvbdistr3 calculates initial particle velocities in 2-1/2d for
!            magnetized plasma with maxwellian velocity with drift in
!            Bparallel direction and ring distribution in Bperp
!            calls PVBDISTR32
! mpdblkp3 finds the maximum number of particles in each tile
!          calls PPDBLKP3L
! mpfedges3 finds new 2d partitions from initial analytic distribution
!           function
!           calls PFEDGES32
! mpfholes3 determines list of particles leaving this processor
!           calls PFHOLES32
! wmpdistr3 generic procedure to initialize particle co-ordinates and
!           velocities/momenta with uniform density and
!           maxwellian/juttner distributions with drift for 3d code
!           calls mpdistr3 or mprdistr3
! wmpvdistr3 generic procedure to initialize particle velocities in 3d
!            calls mpvdistr3 or mpvrdistr3
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: march 23, 2018
!
      use libmpinit3_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mnextran3(nextrand,ndim,np)
! for 3d code, this subroutine skips over nextrand groups of random
! numbers in order to initialize different random ensembles
      implicit none
      integer, intent(in) :: nextrand, ndim, np
! call low level procedure
      call NEXTRAN3(nextrand,ndim,np)
      end subroutine
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
      subroutine mpfdistr3(part,npp,ampx,scalex,shiftx,ampy,scaley,     &
     &shifty,ampz,scalez,shiftz,npx,npy,npz,nx,ny,nz,kstrt,nvpy,nvpz,   &
     &ipbc,ndpro,irc)
! calculates initial particle co-ordinates in 3d
! with various density profiles
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
! ampdx/ampdy/ampdz = amplitude of density compared to uniform in x/y/z
! scaledx/scaledy/scaledz = scale length for spatial coordinate in x/y/z
! shiftdx/shiftdy/shiftdz = shift of spatial coordinate in x/y/z
! irc = (0,N) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: npx, npy, npz, nx, ny, nz, kstrt
      integer, intent(in) :: nvpy, nvpz, ipbc, ndpro
      integer, intent(inout) :: npp, irc
      real, intent(in) :: ampx, scalex, shiftx, ampy, scaley, shifty
      real, intent(in) :: ampz, scalez, shiftz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: nvp, idimp, npmax, ierr
      real :: zero
      double precision :: dmpx, sxi, dshiftx, dmpy, syi, dshifty
      double precision :: dmpz, szi, dshiftz
      double precision, external :: DSDISTR1, DGDISTR1, DHDISTR1
      double precision, external :: DEDISTR1
      idimp = size(part,1); npmax = size(part,2)
      nvp = nvpy*nvpz
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
      dmpz = dble(ampz)
      szi = 0.0d0
      if (scalez /= 0.0) szi = 1.0d0/dble(scalez)
      dshiftz = dble(shiftz)
      zero = 0.0
! call low level procedure
      select case(ndpro)
! uniform density
      case (0)
         call PLDISTR32(part,npp,zero,zero,zero,npx,npy,npz,nx,ny,nz,   &
     &kstrt,nvpy,nvpz,idimp,npmax,ipbc,ierr)
! linear density
      case (1)
         if ((ampx.le.2.0).and.(ampy.le.2.0).and.(ampz.le.2.0)) then
            call PLDISTR32(part,npp,ampx,ampy,ampz,npx,npy,npz,nx,ny,nz,&
     &kstrt,nvpy,nvpz,idimp,npmax,ipbc,ierr)
        endif
! sinusoidal density
      case (2)
         if ((ampx.le.1.0).and.(ampy.le.1.0).and.(ampz.le.1.0)) then
            call PFDISTR32(part,npp,DSDISTR1,dmpx,sxi,dshiftx,DSDISTR1, &
     &dmpy,syi,dshifty,DSDISTR1,dmpz,szi,dshiftz,npx,npy,npz,nx,ny,nz,  &
     &kstrt,nvpy,nvpz,idimp,npmax,ipbc,ierr)
         endif
! gaussian density
      case (3)
         call PFDISTR32(part,npp,DGDISTR1,dmpx,sxi,dshiftx,DGDISTR1,dmpy&
     &,syi,dshifty,DGDISTR1,dmpz,szi,dshiftz,npx,npy,npz,nx,ny,nz,kstrt,&
     &nvpy,nvpz,idimp,npmax,ipbc,ierr)
! hyperbolic secant squared density
      case (4)
         call PFDISTR32(part,npp,DHDISTR1,dmpx,sxi,dshiftx,DHDISTR1,dmpy&
     &,syi,dshifty,DHDISTR1,dmpz,szi,dshiftz,npx,npy,npz,nx,ny,nz,kstrt,&
     &nvpy,nvpz,idimp,npmax,ipbc,ierr)
! exponential density
      case (5)
         call PFDISTR32(part,npp,DEDISTR1,dmpx,sxi,dshiftx,DEDISTR1,dmpy&
     &,syi,dshifty,DEDISTR1,dmpz,szi,dshiftz,npx,npy,npz,nx,ny,nz,kstrt,&
     &nvpy,nvpz,idimp,npmax,ipbc,ierr)
      case default
         ierr = -3
      end select
! check error
      if (ierr /= 0) then
         write (*,*) 'mpfdistr3 error: ndpro, ierr=', ndpro, ierr
         irc = ierr
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvdistr3(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,&
     &npz,kstrt,nvpy,nvpz,irc)
! calculates initial particle velocities in 3d
! with maxwellian velocity with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, npz, kstrt, nvpy, nvpz
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,  &
     &kstrt,nvpy,nvpz,idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvdistr3:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvrdistr3(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,&
     &npy,npz,kstrt,nvpy,nvpz,irc)
! calculates initial particle momenta in 3d
! with maxwell-juttner distribution with drift
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, npz, kstrt, nvpy, nvpz
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVRDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,npy,  &
     &npz,kstrt,nvpy,nvpz,idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvrdistr3:particle number error, irc=', irc
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvbdistr3(part,nps,npp,vtr,vtz,vdr,vdz,omx,omy,omz,npx&
     &,npy,npz,kstrt,nvpy,nvpz,irc)
! calculates initial particle velocities in 3d
! for magnetized plasma with maxwellian velocity with drift in Bparallel
! direction and ring distribution in Bperp
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, npz, kstrt, nvpy, nvpz
      integer, intent(inout) :: irc
      real, intent(in) :: vtr, vtz, vdr, vdz, omx, omy, omz
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
      irc = 0
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! call low level procedure
      call PVBDISTR32(part,nps,npp,vtr,vtz,vdr,vdz,omx,omy,omz,npx,npy, &
     &npz,kstrt,nvpy,nvpz,idimp,npmax,irc)
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpvbdistr3:particle number error, irc=', irc
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
      subroutine mpfedges3(edges,nyzp,noff,ampy,scaley,shifty,ampz,     &
     &scalez,shiftz,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy,nvpz,ipbc, &
     &ndpro,ierr)
! finds new 2d partitions from initial analytic distribution function
      implicit none
      integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, ipbc, ndpro
      integer, intent(inout) :: nypmx, nypmn, nzpmx, nzpmn, ierr
      real, intent(in) :: ampy, scaley, shifty, ampz, scalez, shiftz
      real, dimension(:), intent(inout) :: edges
      integer, dimension(:), intent(inout) :: nyzp, noff
! local data
      integer :: idps, idds
      double precision :: dmpy, syi, dshifty, dmpz, szi, dshiftz, zero
      double precision, external :: DLDISTR1, DSDISTR1, DGDISTR1
      double precision, external :: DHDISTR1, DEDISTR1
      idps = size(edges,1); idds = size(noff)
      ierr = 0
      dmpy = dble(ampy)
      syi = 0.0d0
      if (scaley /= 0.0) syi = 1.0d0/dble(scaley)
      dshifty = dble(shifty)
      dmpz = dble(ampz)
      szi = 0.0d0
      if (scalez /= 0.0) szi = 1.0d0/dble(scalez)
      dshiftz = dble(shiftz)
      zero = 0.0d0
! call low level procedure
      select case(ndpro)
! uniform density
      case (0)
         call PFEDGES32(edges,nyzp,noff,DLDISTR1,zero,zero,zero,DLDISTR1&
     &,zero,zero,zero,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy,nvpz,idps&
     &,idds,ipbc)
! linear density
      case (1)
         if ((ampy.le.2.0).and.(ampy.le.2.0)) then
            call PFEDGES32(edges,nyzp,noff,DLDISTR1,dmpy,syi,dshifty,   &
     &DLDISTR1,dmpz,szi,dshiftz,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy&
     &,nvpz,idps,idds,ipbc)
         else
            ierr = -4
         endif
! sinusoidal density
      case (2)
         if ((ampy.le.1.0).and.(ampy.le.1.0)) then
            call PFEDGES32(edges,nyzp,noff,DSDISTR1,dmpy,syi,dshifty,   &
     &DSDISTR1,dmpz,szi,dshiftz,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy&
     &,nvpz,idps,idds,ipbc)
         else
            ierr = -5
         endif
! gaussian density
      case (3)
         call PFEDGES32(edges,nyzp,noff,DGDISTR1,dmpy,syi,dshifty,      &
     &DGDISTR1,dmpz,szi,dshiftz,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy&
     &,nvpz,idps,idds,ipbc)
! hyperbolic secant squared density
      case (4)
         call PFEDGES32(edges,nyzp,noff,DHDISTR1,dmpy,syi,dshifty,      &
     &DHDISTR1,dmpz,szi,dshiftz,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy&
     &,nvpz,idps,idds,ipbc)
! exponential density
      case (5)
         call PFEDGES32(edges,nyzp,noff,DEDISTR1,dmpy,syi,dshifty,      &
     &DEDISTR1,dmpz,szi,dshiftz,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy&
     &,nvpz,idps,idds,ipbc)
      case default
         ierr = -3
      end select
! check error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpfedges3 error: ndpro, ierr=', ndpro, ierr
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfholes3(part,edges,npp,iholep,nc)
! determines list of particles which are leaving this node
      implicit none
      integer, intent(in) :: npp, nc
      real, dimension(:,:), intent(in) :: part
      real, dimension(:), intent(in) :: edges
      integer, dimension(:,:), intent(inout) :: iholep
! local data
      integer :: idimp, npmax, idps, ntmax
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1); ntmax = size(iholep,1) - 1
! call low level procedure
      if ((nc==1).or.((nc==2).and.(idimp > 6))) then
         call PFHOLES32(part,edges,npp,iholep,nc,idimp,npmax,idps,ntmax)
      else
         write (*,*) 'mpfholes3 error: nc, idimp=', nc, idimp
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
!-----------------------------------------------------------------------
      subroutine wmpvdistr3(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,&
     &npy,npz,kstrt,nvpy,nvpz,relativity,irc)
! generic procedure to initialize particle velocities in 3d
      implicit none
      integer, intent(in) :: nps, npp, npx, npy, npz, kstrt, nvpy, nvpz
      integer, intent(in) :: relativity
      integer, intent(inout) :: irc
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
      real, dimension(:,:), intent(inout) :: part
! maxwell-juttner distribution
      if (relativity==1) then
         call mpvrdistr3(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,npy&
     &,npz,kstrt,nvpy,nvpz,irc)
! maxwellian distribution
      else
         call mpvdistr3(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz&
     &,kstrt,nvpy,nvpz,irc)
      endif
      end subroutine
!
      end module
