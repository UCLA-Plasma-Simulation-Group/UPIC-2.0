!-----------------------------------------------------------------------
!
      module mppmod3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library mpplib3.f90
! these wrappers encapulate MPI communications
! mpcguard3 copies scalar data to guard cells in y and z
!           calls PPNCGUARD32L
! mpncguard3 copies vector data to guard cells in y and z
!            calls PPNCGUARD32L
! mpnaguard3 adds scalar data from guard cells in y and z
!            calls PPNAGUARD32L
! mpnacguard3 adds vector data from guard cells in y and z
!             calls PPNACGUARD32L
! mpfmove3 moves scalar grids into appropriate spatial regions in y/z
!          calls PPFYMOVE32 and PPFZMOVE32
! mpfnmove3 moves vector grids into appropriate spatial regions in y/z
!           calls PPFYMOVE32 and PPFZMOVE32
! mpmove3 moves particles into appropriate spatial regions in y and z
!         calls PPPMOVE32
! mpimax finds parallel maximum for each element of an integer vector
!        calls PPIMAX
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 16, 2017
!
      use mpplib3
      implicit none
!
! scr/scs = guard cell buffers received from nearby processors
      real, dimension(:,:), allocatable  :: scr
      real, dimension(:,:,:), allocatable  :: scs
      integer :: szscr = 0, szscs = 0
! gbuf/hbuf = scratch memorys for field manager
      real, dimension(:,:), allocatable :: gbuf, hbuf
      integer :: szghbuf = 0
! ig = scratch buffer for integer reduction
      integer, dimension(:), allocatable :: ig
      integer :: szig = 0
      save
!
      private :: gbuf, hbuf, scr, scs, szscr, szghbuf
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpcguard3(f,nyzp,tguard,kstrt,nvpy,nvpz)
! copies scalar data to guard cells in non-uniform partitions in y and z
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxv, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
      idds = size(nyzp,1)
! check if required size of buffer has increased
      if (szscs < nxv*nzpmx*2) then
         if (szscs /= 0) deallocate(scs)
! allocate new buffer
         allocate(scs(nxv,nzpmx,2))
         szscs = nxv*nzpmx*2
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpncguard3(f,nyzp,tguard,kstrt,nvpy,nvpz)
! copies vector data to guard cells in non-uniform partitions in y and z
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz
      real, intent(inout) :: tguard
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxv, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1)*size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
      idds = size(nyzp,1)
! check if required size of buffer has increased
      if (szscs < nxv*nzpmx*2) then
         if (szscs /= 0) deallocate(scs)
! allocate new buffer
         allocate(scs(nxv,nzpmx,2))
         szscs = nxv*nzpmx*2
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpnaguard3(f,nyzp,tguard,nx,kstrt,nvpy,nvpz)
! adds scalar data from guard cells in non-uniform partitions in y and z
      implicit none
      integer, intent(in) :: nx, kstrt, nvpy, nvpz
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxv, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
      idds = size(nyzp,1)
! check if required size of buffers has increased
      if (szscr < nxv*nypmx) then
         if (szscr /= 0) deallocate(scr)
! allocate new buffer
         allocate(scr(nxv,nypmx))
         szscr = nxv*nypmx
      endif
      if (szscs < nxv*nzpmx*2) then
         if (szscs /= 0) deallocate(scs)
! allocate new buffer
         allocate(scs(nxv,nzpmx,2))
         szscs = nxv*nzpmx*2
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx,    &
     &nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpnacguard3(f,nyzp,tguard,nx,kstrt,nvpy,nvpz)
! adds vector data from guard cells in non-uniform partitions in y and z
      implicit none
      integer, intent(in) :: nx, kstrt, nvpy, nvpz
      real, intent(inout) :: tguard
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: ndim, nxv, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1)
      nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
      idds = size(nyzp,1)
! check if required size of buffers has increased
      if (szscr < ndim*nxv*nypmx) then
         if (szscr /= 0) deallocate(scr)
! allocate new buffer
         allocate(scr(ndim*nxv,nypmx))
         szscr = ndim*nxv*nypmx
      endif
      if (szscs < ndim*nxv*nzpmx*2) then
         if (szscs /= 0) deallocate(scs)
! allocate new buffer
         allocate(scs(ndim*nxv,nzpmx,2))
         szscs = ndim*nxv*nzpmx*2
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNACGUARD32L(f,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxv,    &
     &nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfmove3(f,noff,nyzp,isign,tfmov,kyp,kzp,ny,nz,kstrt,  &
     &nvpy,nvpz,mter,ierr)
! field manager: moves scalar grids into appropriate spatial regions
! in y and z
      implicit none
      integer, intent(in) :: isign, kyp, kzp, ny, nz, kstrt, nvpy, nvpz
      integer, intent(inout) :: ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: noff, nyzp
      integer, dimension(2), intent(inout) :: mter
! local data
      integer :: nxv, nypmx, nzpmx, idds
      integer, dimension(size(noff,1)) :: noffs, nyzps, noffd, nyzpd
      integer, dimension(4) :: itime
      double precision :: dtime
      if (isign==0) return
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
      idds = size(noff,1)
! check if required size of buffer has increased
      if (szghbuf < nxv*nypmx*nzpmx) then
         if (szghbuf /= 0) deallocate(gbuf,hbuf)
! allocate new buffer
         allocate(gbuf(nxv,nypmx*nzpmx),hbuf(nxv,nypmx*nzpmx))
         szghbuf = nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedures
! move from non-uniform to uniform field
      if (isign < 0) then
         call PPFZMOVE32(f,gbuf,hbuf,noff,nyzp,noffs,nyzps,noffd,nyzpd, &
     &isign,kyp,kzp,ny,nz,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds,mter(2), &
     &ierr)
         call PPFYMOVE32(f,gbuf,hbuf,noff,nyzp,noffs,nyzps,noffd,nyzpd, &
     &isign,kyp,kzp,ny,nz,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds,mter(1), &
     &ierr)
! move from uniform to non-uniform field
      else if (isign > 0) then
         call PPFYMOVE32(f,gbuf,hbuf,noff,nyzp,noffs,nyzps,noffd,nyzpd, &
     &isign,kyp,kzp,ny,nz,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds,mter(1), &
     &ierr)
         call PPFZMOVE32(f,gbuf,hbuf,noff,nyzp,noffs,nyzps,noffd,nyzpd, &
     &isign,kyp,kzp,ny,nz,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds,mter(2), &
     &ierr)
! move from non-uniform to non-uniform field
      else
         if (kstrt==1) then
            write (*,*) 'isign = 0 requires a different function'
            call PPEXIT()
            stop
         endif
      endif
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
! check for error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpfmove3 repartition error: ierr=', ierr
         endif
         call PPEXIT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfnmove3(f,noff,nyzp,isign,tfmov,kyp,kzp,ny,nz,kstrt, &
     &nvpy,nvpz,mter,ierr)
! field manager: moves vector grids into appropriate spatial regions
! in y and z
      implicit none
      integer, intent(in) :: isign, kyp, kzp, ny, nz, kstrt, nvpy, nvpz
      integer, intent(inout) :: ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: noff, nyzp
      integer, dimension(2), intent(inout) :: mter
! local data
      integer :: nxv, nypmx, nzpmx, idds
      integer, dimension(size(noff,1)) :: noffs, nyzps, noffd, nyzpd
      integer, dimension(4) :: itime
      double precision :: dtime
      if (isign==0) return
! extract dimensions
      nxv = size(f,1)*size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
      idds = size(noff,1)
! check if required size of buffer has increased
      if (szghbuf < nxv*nypmx*nzpmx) then
         if (szghbuf /= 0) deallocate(gbuf,hbuf)
! allocate new buffer
         allocate(gbuf(nxv,nypmx*nzpmx),hbuf(nxv,nypmx*nzpmx))
         szghbuf = nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedures
! move from non-uniform to uniform field in z
      if (isign < 0) then
         call PPFZMOVE32(f,gbuf,hbuf,noff,nyzp,noffs,nyzps,noffd,nyzpd, &
     &isign,kyp,kzp,ny,nz,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds,mter(2), &
     &ierr)
         call PPFYMOVE32(f,gbuf,hbuf,noff,nyzp,noffs,nyzps,noffd,nyzpd, &
     &isign,kyp,kzp,ny,nz,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds,mter(1), &
     &ierr)
! move from uniform in z to non-uniform fields
      else if (isign > 0) then
         call PPFYMOVE32(f,gbuf,hbuf,noff,nyzp,noffs,nyzps,noffd,nyzpd, &
     &isign,kyp,kzp,ny,nz,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds,mter(1), &
     &ierr)
         call PPFZMOVE32(f,gbuf,hbuf,noff,nyzp,noffs,nyzps,noffd,nyzpd, &
     &isign,kyp,kzp,ny,nz,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds,mter(2), &
     &ierr)
      endif
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
! check for error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'mpfnmove3 repartition error: ierr=', ierr
         endif
         call PPEXIT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpmove3(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,   &
     &mcls,tmov,kstrt,nvpy,nvpz,myp1,mzp1,irc)
! particle manager: moves particles into appropriate spatial regions
! in y and z
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz, myp1, mzp1
      integer, intent(inout) :: irc
      real, intent(inout) :: tmov
      real, dimension(:,:,:), intent(in) :: sbufr, sbufl
      real, dimension(:,:,:), intent(inout) :: rbufr, rbufl
      integer, dimension(:,:,:,:), intent(inout) :: ncll, nclr
      integer, dimension(:,:,:,:), intent(inout) :: mcll, mclr
      integer, dimension(:,:,:), intent(inout) :: mcls
! local data
      integer :: idimp, nbmax, mxzyp1, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(sbufl,1); nbmax = size(sbufl,2)
      mxzyp1 = size(ncll,2); mx1 = size(mcls,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPMOVE32(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,mcls,  &
     &kstrt,nvpy,nvpz,idimp,nbmax,mx1,myp1,mzp1,mxzyp1,irc)
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
