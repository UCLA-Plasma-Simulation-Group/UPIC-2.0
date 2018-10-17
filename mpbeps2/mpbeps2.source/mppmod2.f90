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
! mpsum finds parallel sum for each element of a real 1d vector
!       calls PPSUM
! mpdsum finds parallel sum for each element of a double precision 1d
!        vector
!        calls PPDSUM
! mpdmax finds parallel maximum for each element of a double precision
!        1d vector
!        calls PPDMAX
! mpsum2 finds parallel sum for each element of a real 2d vector
!        calls PPSUM
! mpimax finds parallel maximum for each element of an integer vector
!        calls PPIMAX
! mpwrite2 collects a subset of a distributed real 2d scalar array and
!          writes it to a direct access binary file
!          calls PPWRITE2
! mpread2 reads a real 2d scalar array from a direct access binary file
!         and distributes it
!         calls PPREAD2
! mpvwrite2 collects a subset of a distributed real 2d vector array and
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
! mpvcwrite2 collects a subset of a distributed complex 2d vector array
!            and writes it to a direct access binary file
!            calls PPVCWRITE2
! mpvcread2 reads a complex 2d vector array from a direct access binary
!           file and distributes it
!           calls PPVCREAD2
! mpwrpart2 collects distributed particle data and writes to a fortran
!           unformatted file
!           calls PPWRPART2
! mprdpart2 reads particle data from a fortran unformatted file and
!           distributes it
!           calls PPRDPART2
! mpwrdata2 collects distributed real scalar field data and writes to a
!           fortran unformatted file
!           calls PPWRDATA2
! mprddata2 reads real scalar field data from a fortran unformatted file
!           and distributes it
!           calls PPRDDATA2
! mpwrvdata2 collects distributed real vector field data and writes
!             to a fortran unformatted file
!             calls PPWRVDATA2
! mprdvdata2 reads real vector field data from a fortran unformatted
!             file and distributes it
!             calls PPRDVDATA2
! mpwrvcdata2 collects distributed complex vector field data and writes
!             to a fortran unformatted file
!             calls PPWRVCDATA2
! mprdvcdata2 reads complex vector field data from a fortran unformatted
!             file and distributes it
!             calls PPRDVCDATA2
! mppartt2 collects distributed test particle data
!          calls PPARTT2
! mpadjfvs2 adjusts 3d velocity distribution in different regions of
!           space, so that partial regions have equal grid points
!           calls PPADJFVS2
! mpwrfvsdata2 collects distributed 4d real vector non-uniform data and
!              writes to a fortran unformatted file
!              calls PPWRVNDATA2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: august 10, 2018
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
! g = scratch buffer for real reduction
      real, dimension(:), allocatable :: g
      integer :: szg = -1
! dg = scratch buffer for real reduction
      double precision, dimension(:), allocatable :: dg
      integer :: szdg = -1
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
! gvs = scratch buffer for phase space calculation
      real, dimension(:,:,:), allocatable :: gvs
      integer :: szgvs = -1
      save
!
      private :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      private :: scr, fbuf, g, dg, ig, fg, fvg, fgc, fvgc, gvs
      private :: szscr, szfbuf, szg, szdg, szig, szfg, szfvg
      private :: szfgc, szfvgc, szbuf, szgvs
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
      subroutine ipmove2(part,edges,npp,iholep,ny,tmov,kstrt,nvp,ndim,nc&
     &,ierr)
! particle manager: moves particles into appropriate spatial regions
      implicit none
      integer, intent(in) :: ny, kstrt, nvp, ndim, nc
      real, intent(inout) :: tmov
      integer, intent(inout) :: npp, ierr
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
      integer, dimension(:), intent(in) :: iholep
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
      if ((nc==1).or.((nc==2).and.(idimp > (ndim+2)))) then
         call PPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,iholep,ny, &
     &kstrt,nvp,ndim,nc,idimp,npmax,idps,nbmax,ntmax,info)
      ierr = info(1)
      else
         write (*,*) 'ipmove2 error: nc, idimp=', nc, idimp
      endif
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
      subroutine mpsum(f,tdiag)
! finds parallel sum for each element of a real 1d vector
      implicit none
      real, intent(inout) :: tdiag
      real, dimension(:), intent(inout) :: f
! local data
      integer :: nxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxp = size(f)
! check if required size of buffer has increased
      if (szg < nxp) then
         if (szg > 0) deallocate(g)
! allocate new buffer
         allocate(g(nxp))
         szg = nxp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPSUM(f,g,nxp)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdsum(f,tdiag)
! finds parallel sum for each element of a double precision 1d vector
      implicit none
      real, intent(inout) :: tdiag
      double precision, dimension(:), intent(inout) :: f
! local data
      integer :: nxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxp = size(f)
! check if required size of buffer has increased
      if (szdg < nxp) then
         if (szdg > 0) deallocate(dg)
! allocate new buffer
         allocate(dg(nxp))
         szdg = nxp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPDSUM(f,dg,nxp)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdmax(f,tdiag)
! finds parallel maximum for each element of a double precision 1d
! vector
      implicit none
      real, intent(inout) :: tdiag
      double precision, dimension(:), intent(inout) :: f
! local data
      integer :: nxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxp = size(f)
! check if required size of buffer has increased
      if (szdg < nxp) then
         if (szdg > 0) deallocate(dg)
! allocate new buffer
         allocate(dg(nxp))
         szdg = nxp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPDMAX(f,dg,nxp)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpsum2(f,tdiag)
! finds parallel sum for each element of a real 2d vector
      implicit none
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(inout) :: f
! local data
      integer :: nxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxp = size(f)
! check if required size of buffer has increased
      if (szg < nxp) then
         if (szg > 0) deallocate(g)
! allocate new buffer
         allocate(g(nxp))
         szg = nxp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPSUM(f,g,nxp)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
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
      if (nrec < 0) return
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
      if (nrec < 0) return
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
      if (nrec < 0) return
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
      if (nrec < 0) return
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
      if (nrec < 0) return
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
      if (nrec < 0) return
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
      if (nrec < 0) return
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
      if (nrec < 0) return
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
!-----------------------------------------------------------------------
      subroutine mpwrpart2(part,tdiag,npp,iunit,iscr)
! collects distributed particle data and writes to a fortran unformatted
! file
      implicit none
      integer, intent(in) :: npp, iunit, iscr
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRPART2(part,npp,idimp,npmax,iunit,iscr)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdpart2(part,tdiag,npp,iunit,iscr,irc)
! reads particle data from a fortran unformatted file and distributes it
      implicit none
      integer, intent(in) :: iunit, iscr
      integer, intent(inout) :: npp, irc
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, npmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRDPART2(part,npp,idimp,npmax,iunit,iscr,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrdata2(f,tdiag,iunit)
! collects distributed real scalar field data and writes to a fortran
! unformatted file
      implicit none
      integer, intent(in) :: iunit
      real, intent(inout) :: tdiag
      real, intent(in), dimension(:,:) :: f
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
      call PPWRDATA2(f,fg,nxv,nypmx,iunit)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprddata2(f,tdiag,iunit,irc)
! reads real scalar field data from a fortran unformatted file and
! distributes it
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: irc
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
      call PPRDDATA2(f,fg,nxv,nypmx,iunit,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrvdata2(f,tdiag,iunit)
! collects distributed real vector field data and writes to a fortran
! unformatted file
      implicit none
      integer, intent(in) :: iunit
      real, intent(inout) :: tdiag
      real, intent(in), dimension(:,:,:) :: f
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
      call PPWRVDATA2(f,fvg,ndim,nxv,nypmx,iunit)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdvdata2(f,tdiag,iunit,irc)
! reads real vector field data from a fortran unformatted file and
! distributes it
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: irc
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
      call PPRDVDATA2(f,fvg,ndim,nxv,nypmx,iunit,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrvcdata2(f,tdiag,iunit)
! collects distributed complex vector field data and writes to a fortran
! unformatted file
      implicit none
      integer, intent(in) :: iunit
      real, intent(inout) :: tdiag
      complex, intent(in), dimension(:,:,:) :: f
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
      call PPWRVCDATA2(f,fvgc,ndim,nyv,kxpd,iunit)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdvcdata2(f,tdiag,iunit,irc)
! reads complex vector field data from a fortran unformatted file and
! distributes it
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: irc
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
      call PPRDVCDATA2(f,fvgc,ndim,nyv,kxpd,iunit,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppartt2(partt,tdiag,numtp,irc)
! collects distributed test particle data
      implicit none
      integer, intent(in) :: numtp
      integer, intent(inout) :: irc
      real, intent(inout) :: tdiag
      real, intent(inout), dimension(:,:) :: partt
! local data
      integer :: idimp, nprobt
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(partt,1); nprobt = size(partt,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPARTT2(partt,numtp,idimp,nprobt,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! check for errors
      if (irc /= 0) write (*,*) 'mppartt2: irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpadjfvs2(fvs,tdiag,noff,nyp,nmv,mvy)
! adjusts 3d velocity distribution in different regions of space
      implicit none
      integer, intent(in) :: noff, nyp, nmv, mvy
      real, intent(inout) :: tdiag
      real, dimension(:,:,:,:), intent(inout) :: fvs
! local data
      integer :: nmvf, ndim, nxb, nyb
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nmvf = size(fvs,1); ndim = size(fvs,2)
      nxb = size(fvs,3); nyb = size(fvs,4) - 1
! check if required size of buffer has increased
      if (szgvs < nmvf*ndim*nxb) then
         if (szgvs > 0) deallocate(gvs)
! allocate new buffer
         allocate(gvs(nmvf,ndim,nxb))
         szgvs = nmvf*ndim*nxb
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPADJFVS2(fvs,gvs,noff,nyp,ndim,nmv,mvy,nxb,nyb,nmvf)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrfvsdata2(fvs,tdiag,nyp,nypmx,iunit)
! collects distributed 4d real vector non-uniform data and writes to a
! fortran unformatted file
      implicit none
      integer, intent(in) :: iunit, nyp, nypmx
      real, intent(inout) :: tdiag
      real, intent(in), dimension(:,:,:,:) :: fvs
! local data
      integer :: nndim, nxv
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nndim = size(fvs,1)*size(fvs,2); nxv = size(fvs,3)
! check if required size of buffer has increased
      if (szfvg < nndim*nxv*nypmx) then
         if (szfvg > 0) deallocate(fvg)
! allocate new buffer
         allocate(fvg(nndim,nxv,nypmx))
         szfvg = nndim*nxv*nypmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRVNDATA2(fvs,fvg,nndim,nxv,nyp,nypmx,iunit)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
      end module
