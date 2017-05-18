!-----------------------------------------------------------------------
!
      module mppmod2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library mpplib2.f90
! these wrappers encapulate MPI communications
! mpcguard2 copies scalar data to guard cells
!           calls PPNCGUARD2L
! mpncguard2 copies vector data to guard cells
!            calls PPNCGUARD2L
! mpnaguard2 adds scalar data from guard cells
!            calls PPNAGUARD2L
! mpnacguard2 adds vector data from guard cells
!             calls PPNACGUARD2L
! mpfmove2 moves scalar grids into appropriate spatial regions
!          calls PPFMOVE2
! mpfnmove2 moves vector grids into appropriate spatial regions
!           calls PPFMOVE2
! ipmove2 moves particles to appropriate node
!         calls PPMOVE2
! mpmove2 moves particles into appropriate tiles
!         calls PPPMOVE2
! mpimax finds parallel maximum for each element of an integer vector
!        calls PPIMAX
! mpwrite2 collects a subset of a distributed real 2d scalar array and
!          writes it to a direct access binary file
!          calls PPWRITE2
! mpread2 reads a real 2d scalar array from a direct access binary file
!         and distributes it
!         calls PPREAD2
! mpvwrite2 collects a subset of a distributed real 2d scalar array and
!           writes it to a direct access binary file
!           calls PPVWRITE2
! mpvread2 reads a real 2d vector array from a direct access binary file
!          and distributes it
!          calls PPVREAD2
! mpcwrite2 collects a subset of a distributed complex 2d scalar array
!           and writes it to a direct access binary file
!           calls PPCWRITE2
! mpcread2 reads a complex 2d scalar array from a direct access binary
!          file and distributes it
!          calls PPCREAD2
! mpvcwrite2 collects a subset of a distributed complex 2d scalar array
!            and writes it to a direct access binary file
!            calls PPVCWRITE2
! mpvcread2 reads a complex 2d vector array from a direct access binary
!           file and distributes it
!           calls PPVCREAD2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 10, 2017
!
      use mpplib2
      implicit none
!
! scr = guard cell buffer received from nearby nodes
      real, dimension(:), allocatable  :: scr
      integer :: szscr = -1
! fbuf = scratch memory for field manager
      real, dimension(:,:), allocatable :: fbuf
      integer :: szfbuf = -1
! ig = scratch buffer for integer reduction
      integer, dimension(:), allocatable :: ig
      integer :: szig = -1
! fg = scratch buffer for real scalar file output
      real, dimension(:,:), allocatable :: fg
      integer :: szfg = -1
! fvg = scratch buffer for real vector file output
      real, dimension(:,:,:), allocatable :: fvg
      integer :: szfvg = -1
! fgc = scratch buffer for complex scalar file output
      complex, dimension(:,:), allocatable :: fgc
      integer :: szfgc = -1
! fvgc = scratch buffer for complex vector file output
      complex, dimension(:,:,:), allocatable :: fvgc
      integer :: szfvgc = -1
! buffer data for particle managers
      real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer :: szbuf = -1
      save
!
      private :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      private :: scr, fbuf, ig, fg, fvg, fgc, fvgc
      private :: szscr, szfbuf, szig, szfg, szfvg, szfgc, szfvgc, szbuf
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpcguard2(f,nyp,tguard,kstrt,nvp)
! copies scalar data to guard cells in non-uniform partitions
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpncguard2(f,nyp,tguard,kstrt,nvp)
! copies vector data to guard cells in non-uniform partitions
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1)*size(f,2); nypmx = size(f,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpnaguard2(f,nyp,tguard,nx,kstrt,nvp)
! adds scalar data from guard cells in non-uniform partitions
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2)
! check if required size of buffer has increased
      if (szscr < nxv) then
         if (szscr > 0) deallocate(scr)
! allocate new buffer
         allocate(scr(nxv))
         szscr = nxv
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpnacguard2(f,nyp,tguard,nx,kstrt,nvp)
! adds vector data from guard cells in non-uniform partitions
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxv = size(f,2); nypmx = size(f,3)
! check if required size of buffer has increased
      if (szscr < ndim*nxv) then
         if (szscr > 0) deallocate(scr)
! allocate new buffer
         allocate(scr(ndim*nxv))
         szscr = ndim*nxv
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter, &
     &ierr)
! field manager: moves scalar grids into appropriate spatial regions
      implicit none
      integer, intent(in) :: noff, nyp, isign, kyp, ny, kstrt, nvp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:), intent(inout) :: f
! local data
      integer :: noffs, nyps, noffd, nypd, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
      if (isign==0) return
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2)
! check if required size of buffer has increased
      if (szfbuf < nxv*nypmx) then
         if (szfbuf > 0) deallocate(fbuf)
! allocate new buffer
         allocate(fbuf(nxv,nypmx))
         szfbuf = nxv*nypmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPFMOVE2(f,fbuf,noff,nyp,noffs,nyps,noffd,nypd,isign,kyp,ny, &
     &kstrt,nvp,nxv,nypmx,mter,ierr)
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
! check for error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpfmove2 repartition error: ierr=', ierr
         endif
         call PPEXIT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,&
     &ierr)
! field manager: moves vector grids into appropriate spatial regions
      implicit none
      integer, intent(in) :: noff, nyp, isign, kyp, ny, kstrt, nvp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: noffs, nyps, noffd, nypd, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
      if (isign==0) return
! extract dimensions
      nxv = size(f,1)*size(f,2); nypmx = size(f,3)
! check if required size of buffer has increased
      if (szfbuf < nxv*nypmx) then
         if (szfbuf > 0) deallocate(fbuf)
! allocate new buffer
         allocate(fbuf(nxv,nypmx))
         szfbuf = nxv*nypmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPFMOVE2(f,fbuf,noff,nyp,noffs,nyps,noffd,nypd,isign,kyp,ny, &
     &kstrt,nvp,nxv,nypmx,mter,ierr)
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
! check for error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpfnmove2 repartition error: ierr=', ierr
         endif
         call PPEXIT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ipmove2(part,edges,npp,iholep,ny,tmov,kstrt,nvp,ierr)
! particle manager: moves particles into appropriate spatial regions
      implicit none
      integer, intent(in) :: ny, kstrt, nvp
      real, intent(inout) :: tmov
      integer, intent(inout) :: npp, ierr
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
      integer, dimension(:), intent(inout) :: iholep
! local data
      integer :: idimp, npmax, idps, nbmax, ntmax
      integer, dimension(5) :: info
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1); ntmax = size(iholep,1) - 1
      nbmax = ntmax/2
! check if size of buffers has changed
      if (szbuf < idimp*nbmax) then
         if (szbuf > 0) deallocate(sbufl,sbufr,rbufl,rbufr)
! allocate buffers
         allocate(sbufl(idimp,nbmax),sbufr(idimp,nbmax))
         allocate(rbufl(idimp,nbmax),rbufr(idimp,nbmax))
         szbuf = idimp*nbmax
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,iholep,ny,    &
     &kstrt,nvp,idimp,npmax,idps,nbmax,ntmax,info)
      ierr = info(1)
! record time
      call dtimer(dtime,itime,1)
      tmov = tmov + real(dtime)
! check for errors
      if (ierr /= 0) then
         if (kstrt==1) then
            if (ierr > 0) then
               write (*,*) 'ipmove2 particle overflow error: ierr=',ierr
            else if (ierr < 0) then
               write (*,*) 'ipmove2 Iteration overflow error: ierr=',   &
     &ierr
            endif
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,   &
     &tmov,kstrt,nvp)
! tiled particle manager: moves particles into appropriate tiles
      implicit none
      integer, intent(in) :: kstrt, nvp
      real, intent(inout) :: tmov
      real, dimension(:,:), intent(in) :: sbufl, sbufr
      real, dimension(:,:), intent(inout) :: rbufl, rbufr
      integer, dimension(:,:), intent(in) :: ncll, nclr
      integer, dimension(:,:), intent(inout) :: mcll, mclr
! local data
      integer :: idimp, nbmax, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(sbufl,1); nbmax = size(sbufl,2)
      mx1 = size(ncll,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,  &
     &nvp,idimp,nbmax,mx1)
! record time
      call dtimer(dtime,itime,1)
      tmov = tmov + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpimax(if,tmov)
! finds parallel maximum for each element of an integer vector
      implicit none
      real, intent(inout) :: tmov
      integer, dimension(:), intent(inout) :: if
! local data
      integer :: nxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxp = size(if,1)
! check if required size of buffer has increased
      if (szig < nxp) then
         if (szig > 0) deallocate(ig)
! allocate new buffer
         allocate(ig(nxp))
         szig = nxp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPIMAX(if,ig,nxp)
! record time
      call dtimer(dtime,itime,1)
      tmov = tmov + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrite2(f,tdiag,nx,ny,kyp,iunit,nrec)
! collects a subset of a distributed real 2d scalar array and writes it
! to a direct access binary file, for uniform partitions
      implicit none
      integer, intent(in) :: nx, ny, kyp, iunit
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: f
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2)
! check if required size of buffer has increased
      if (szfg < nxv*nypmx) then
         if (szfg > 0) deallocate(fg)
! allocate new buffer
         allocate(fg(nxv,nypmx))
         szfg = nxv*nypmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRITE2(f,fg,nx,ny,kyp,nxv,nypmx,iunit,nrec)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpread2(f,tdiag,nx,ny,kyp,iunit,nrec,irc)
! reads a real 2d scalar array from a direct access binary file and
! distributes it, for uniform partitions
      implicit none
      integer, intent(in) :: nx, ny, kyp, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2)
! check if required size of buffer has increased
      if (szfg < nxv*nypmx) then
         if (szfg > 0) deallocate(fg)
! allocate new buffer
         allocate(fg(nxv,nypmx))
         szfg = nxv*nypmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPREAD2(f,fg,nx,ny,kyp,nxv,nypmx,iunit,nrec,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpread2 file read error: irc=', irc
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvwrite2(f,tdiag,nx,ny,kyp,iunit,nrec)
! collects a subset of a distributed real 2d vector array and writes it
! to a direct access binary file, for uniform partitions
      implicit none
      integer, intent(in) :: nx, ny, kyp, iunit
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: f
! local data
      integer :: ndim, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxv = size(f,2); nypmx = size(f,3)
! check if required size of buffer has increased
      if (szfvg < ndim*nxv*nypmx) then
         if (szfvg > 0) deallocate(fvg)
! allocate new buffer
         allocate(fvg(ndim,nxv,nypmx))
         szfvg = ndim*nxv*nypmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPVWRITE2(f,fvg,nx,ny,kyp,ndim,nxv,nypmx,iunit,nrec)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvread2(f,tdiag,nx,ny,kyp,iunit,nrec,irc)
! reads a real 2d vector array from a direct access binary file and
! distributes it, for uniform partitions
      implicit none
      integer, intent(in) :: nx, ny, kyp, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxv = size(f,2); nypmx = size(f,3)
! check if required size of buffer has increased
      if (szfvg < ndim*nxv*nypmx) then
         if (szfvg > 0) deallocate(fvg)
! allocate new buffer
         allocate(fvg(ndim,nxv,nypmx))
         szfvg = ndim*nxv*nypmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPVREAD2(f,fvg,nx,ny,kyp,ndim,nxv,nypmx,iunit,nrec,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpvread2 file read error: irc=', irc
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcwrite2(f,tdiag,nx,ny,kxp,iunit,nrec)
! collects a subset of a distributed complex 2d scalar array and writes
! it to a direct access binary file, for uniform partitions
      implicit none
      integer, intent(in) :: nx, ny, kxp, iunit
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      complex, dimension(:,:), intent(in) :: f
! local data
      integer :: nyv, kxpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(f,1); kxpd = size(f,2)
! check if required size of buffer has increased
      if (szfgc < nyv*kxpd) then
         if (szfgc > 0) deallocate(fgc)
! allocate new buffer
         allocate(fgc(nyv,kxpd))
         szfgc = nyv*kxpd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPCWRITE2(f,fgc,nx,ny,kxp,nyv,kxpd,iunit,nrec)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcread2(f,tdiag,nx,ny,kxp,iunit,nrec,irc)
! reads a complex 2d scalar array from a direct access binary file and
! distributes it, for uniform partitions
      implicit none
      integer, intent(in) :: nx, ny, kxp, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout) :: tdiag
      complex, dimension(:,:), intent(inout) :: f
! local data
      integer :: nyv, kxpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(f,1); kxpd = size(f,2)
! check if required size of buffer has increased
      if (szfgc < nyv*kxpd) then
         if (szfgc > 0) deallocate(fgc)
! allocate new buffer
         allocate(fgc(nyv,kxpd))
         szfgc = nyv*kxpd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPCREAD2(f,fgc,nx,ny,kxp,nyv,kxpd,iunit,nrec,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpcread2 file read error: irc=', irc
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvcwrite2(f,tdiag,nx,ny,kxp,iunit,nrec)
! collects a subset of a distributed complex 2d vector array and writes
! it to a direct access binary file, for uniform partitions
      implicit none
      integer, intent(in) :: nx, ny, kxp, iunit
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      complex, dimension(:,:,:), intent(in) :: f
! local data
      integer :: ndim, nyv, kxpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nyv = size(f,2); kxpd = size(f,3)
! check if required size of buffer has increased
      if (szfvgc < ndim*nyv*kxpd) then
         if (szfvgc > 0) deallocate(fvgc)
! allocate new buffer
         allocate(fvgc(ndim,nyv,kxpd))
         szfvgc = ndim*nyv*kxpd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPVCWRITE2(f,fvgc,nx,ny,kxp,ndim,nyv,kxpd,iunit,nrec)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvcread2(f,tdiag,nx,ny,kxp,iunit,nrec,irc)
! reads a complex 2d vector array from a direct access binary file and
! distributes it, for uniform partitions
      implicit none
      integer, intent(in) :: nx, ny, kxp, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout) :: tdiag
      complex, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nyv, kxpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nyv = size(f,2); kxpd = size(f,3)
! check if required size of buffer has increased
      if (szfvgc < ndim*nyv*kxpd) then
         if (szfvgc > 0) deallocate(fvgc)
! allocate new buffer
         allocate(fvgc(ndim,nyv,kxpd))
         szfvgc = ndim*nyv*kxpd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPVCREAD2(f,fvgc,nx,ny,kxp,ndim,nyv,kxpd,iunit,nrec,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpvcread2 file read error: irc=', irc
      endif
      end subroutine
!
      end module
