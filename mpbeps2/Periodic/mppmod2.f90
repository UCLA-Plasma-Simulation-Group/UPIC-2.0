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
! mpmove2 moves particles into appropriate spatial regions
!         calls PPPMOVE2
! mpimax finds parallel maximum for each element of an integer vector
!        calls PPIMAX
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 4, 2017
!
      use mpplib2
      implicit none
!
! scr = guard cell buffer received from nearby processors
      real, dimension(:), allocatable  :: scr
      integer :: szscr = 0
! fbuf = scratch memory for field manager
      real, dimension(:,:), allocatable :: fbuf
      integer :: szfbuf = 0
! ig = scratch buffer for integer reduction
      integer, dimension(:), allocatable :: ig
      integer :: szig = 0
      save
!
      private :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      private :: fbuf, scr, szscr, szfbuf
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
         if (szscr /= 0) deallocate(scr)
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
         if (szscr /= 0) deallocate(scr)
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
         if (szfbuf /= 0) deallocate(fbuf)
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
      if (ierr > 0) then
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
         if (szfbuf /= 0) deallocate(fbuf)
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
      if (ierr > 0) then
         if (kstrt==1) then
            write (*,*) 'mpfnmove2 repartition error: ierr=', ierr
         endif
         call PPEXIT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,   &
     &tmov,kstrt,nvp)
! particle manager: moves particles into appropriate spatial regions
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
         if (szig /= 0) deallocate(ig)
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
      end module
