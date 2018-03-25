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
! ipmove3 moves particles to appropriate node
!         calls PPMOVE32
! mpmove3 moves particles into appropriate spatial regions in y and z
!         calls PPPMOVE32
! mpsum finds parallel sum for each element of a real 1d vector
!       calls PPSUM
! mpdsum finds parallel sum for each element of a double precision 1d
!        vector
!        calls PPDSUM
! mpsum2 finds parallel sum for each element of a real 2d vector
!        calls PPSUM
! mpimax finds parallel maximum for each element of an integer vector
!        calls PPIMAX
! mpwrite3 collects a subset of a distributed real 3d scalar array and
!          writes it to a direct access binary file
!          calls PPWRITE32
! mpread3 reads a real 3d scalar array from a direct access binary file
!         and distributes it
!         calls PPREAD32
! mpvwrite3 collects a subset of a distributed real 3d scalar array and
!           writes it to a direct access binary file
!           calls PPVWRITE32
! mpvread3 reads a real 3d vector array from a direct access binary file
!          and distributes it
!          calls PPVREAD32
! mpwrpart3 collects distributed particle data and writes to a fortran
!           unformatted file
!           calls PPWRPART3
! mprdpart3 reads particle data from a fortran unformatted file and
!           distributes it
!           calls PPRDPART3
! mpwrdata3 collects distributed real scalar field data and writes to a
!           fortran unformatted file
!           calls PPWRDATA3
! mprddata3 reads real scalar field data from a fortran unformatted file
!           and distributes it
!           calls PPRDDATA3
! mpwrvdata3 collects distributed real vector field data and writes to a
!            fortran unformatted file
!            calls PPWRVDATA3
! mprdvdata3 reads real vector field data from a fortran unformatted
!            file and distributes it
!            calls PPRDVDATA3
! mpwrvcdata3 collects distributed complex vector field data and writes
!             to a fortran unformatted file
!             calls PPWRVCDATA3
! mprdvcdata3 reads complex vector field data from a fortran unformatted
!             file and distributes it
!             calls PPRDVCDATA3
! mppartt3 collects distributed test particle data
!          calls PPARTT3
! mpadjfvs3 adjusts 3d velocity distribution in different regions of
!           space, so that partial regions have equal grid points
!           calls PPADJFVS3
! mpwrncomp3 collects distributed non-uniform partition information
!            writes to a fortran unformatted file
!            calls PPWRNCOMP3
! mpwrfvsdata3 collects distributed 5d real vector non-uniform data and
!              writes to a fortran unformatted file
!              calls PPWRVNDATA3
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: march 21, 2018
!
      use mpplib3
      implicit none
!
! scr/scs = guard cell buffers received from nearby processors
      real, dimension(:,:), allocatable  :: scr
      real, dimension(:,:,:), allocatable  :: scs
      integer :: szscr = 0, szscs = -1
! gbuf/hbuf = scratch memorys for field manager
      real, dimension(:,:), allocatable :: gbuf, hbuf
      integer :: szghbuf = -1
! g = scratch buffer for real reduction
      real, dimension(:), allocatable :: g
      integer :: szg = -1
! dg = scratch buffer for real reduction
      double precision, dimension(:), allocatable :: dg
      integer :: szdg = -1
! ig = scratch buffer for integer reduction
      integer, dimension(:), allocatable :: ig
      integer :: szig = -1
! fg = scratch buffer for scalar file output
      real, dimension(:,:,:), allocatable :: fg
      integer :: szfg = -1
! fvg = scratch buffer for vector file output
      real, dimension(:,:,:,:), allocatable :: fvg
      integer :: szfvg = -1
! fvgc = scratch buffer for complex vector file output
      complex, dimension(:,:,:,:), allocatable :: fvgc
      integer :: szfvgc = -1
! buffer data for particle managers
      real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer :: szbuf = -1
! gvs = scratch buffer for phase space calculation
      real, dimension(:,:,:,:,:), allocatable :: gvs
      integer :: szgvs = -1
! hvs = scratch buffer for phase space calculation
      real, dimension(:,:,:,:), allocatable :: hvs
      integer :: szhvs = -1
      save
!
      private :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      private :: gbuf, hbuf, scr, scs, g, dg, ig, fg, fvg, fvgc, gvs
      private :: hvs
      private :: szscr, szghbuf, szg, szdg, szig, szfg, szfvg, szfvgc
      private :: szbuf, szgvs, szhvs
!
      contains
!
      subroutine mppdelszbuf()
! allocate szbuf buffers used by ipmove3
      implicit none
      deallocate(sbufl,sbufr,rbufl,rbufr)
      szbuf = 0
      end subroutine
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
         if (szscs > 0) deallocate(scs)
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
         if (szscs > 0) deallocate(scs)
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
         if (szscr > 0) deallocate(scr)
! allocate new buffer
         allocate(scr(nxv,nypmx))
         szscr = nxv*nypmx
      endif
      if (szscs < nxv*nzpmx*2) then
         if (szscs > 0) deallocate(scs)
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
         if (szscr > 0) deallocate(scr)
! allocate new buffer
         allocate(scr(ndim*nxv,nypmx))
         szscr = ndim*nxv*nypmx
      endif
      if (szscs < ndim*nxv*nzpmx*2) then
         if (szscs > 0) deallocate(scs)
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
         if (szghbuf > 0) deallocate(gbuf,hbuf)
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
         if (szghbuf > 0) deallocate(gbuf,hbuf)
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
      subroutine ipmove3(part,edges,npp,iholep,ny,nz,tmov,kstrt,nvpy,   &
     &nvpz,nc,ierr)
! particle manager: moves particles into appropriate spatial regions
      implicit none
      integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, nc
      real, intent(inout) :: tmov
      integer, intent(inout) :: npp, ierr
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
      integer, dimension(:,:), intent(inout) :: iholep
! local data
      integer :: idimp, npmax, idps, nbmax, ntmax
      integer, dimension(7) :: info
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
      if ((nc==1).or.((nc==2).and.(idimp > 6))) then
         call PPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,iholep,ny,&
     &nz,kstrt,nvpy,nvpz,nc,idimp,npmax,idps,nbmax,ntmax,info)
         ierr = info(1)
      else
         write (*,*) 'ipmove3 error: nc, idimp=', nc, idimp
      endif
! record time
      call dtimer(dtime,itime,1)
      tmov = tmov + real(dtime)
! check for errors
      if (ierr /= 0) then
         if (kstrt==1) then
            if (ierr > 0) then
               write (*,*) 'ipmove3 particle overflow error: ierr=',ierr
            else if (ierr < 0) then
               write (*,*) 'ipmove3 iholep overflow error: ierr=',ierr
            endif
         endif
      endif
      end subroutine
!
!!-----------------------------------------------------------------------
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
      subroutine mpwrite3(f,tdiag,nx,ny,nz,kyp,kzp,nvpy,iunit,nrec)
! collects a subset of a distributed real 3d scalar array and writes it
! to a direct access binary file, for uniform 2d partitions
      implicit none
      integer, intent(in) :: nx, ny, nz, kyp, kzp, nvpy, iunit
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: f
! local data
      integer :: nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
      if (nrec < 0) return
! check if required size of buffer has increased
      if (szfg < nxv*nypmx*nzpmx) then
         if (szfg > 0) deallocate(fg)
! allocate new buffer
         allocate(fg(nxv,nypmx,nzpmx))
         szfg = nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRITE32(f,fg,nx,ny,nz,kyp,kzp,nvpy,nxv,nypmx,nzpmx,iunit,  &
     &nrec)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpread3(f,tdiag,nx,ny,nz,kyp,kzp,nvpy,iunit,nrec,irc)
! reads a real 3d scalar array from a direct access binary file and
! distributes it, for uniform 2d partitions
      implicit none
      integer, intent(in) :: nx, ny, nz, kyp, kzp, nvpy, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
      if (nrec < 0) return
! check if required size of buffer has increased
      if (szfg < nxv*nypmx*nzpmx) then
         if (szfg > 0) deallocate(fg)
! allocate new buffer
         allocate(fg(nxv,nypmx,nzpmx))
         szfg = nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure-
      call PPREAD32(f,fg,nx,ny,nz,kyp,kzp,nvpy,nxv,nypmx,nzpmx,iunit,   &
     &nrec,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpread3 file read error: irc=', irc
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvwrite3(f,tdiag,nx,ny,nz,kyp,kzp,nvpy,iunit,nrec)
! collects a subset of a distributed real 3d vector array and writes it
! to a direct access binary file, for uniform 2d partitions
      implicit none
      integer, intent(in) :: nx, ny, nz, kyp, kzp, nvpy, iunit
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      real, dimension(:,:,:,:), intent(in) :: f
! local data
      integer :: ndim, nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxv = size(f,2)
      nypmx = size(f,3); nzpmx = size(f,4)
      if (nrec < 0) return
! check if required size of buffer has increased
      if (szfvg < ndim*nxv*nypmx*nzpmx) then
         if (szfvg > 0) deallocate(fvg)
! allocate new buffer
         allocate(fvg(ndim,nxv,nypmx,nzpmx))
         szfvg = ndim*nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPVWRITE32(f,fvg,nx,ny,nz,kyp,kzp,nvpy,ndim,nxv,nypmx,nzpmx, &
     &iunit,nrec)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvread3(f,tdiag,nx,ny,nz,kyp,kzp,nvpy,iunit,nrec,irc)
! reads a real 3d vector array from a direct access binary file and
! distributes it, for uniform 2d partitions
      implicit none
      integer, intent(in) :: nx, ny, nz, kyp, kzp, nvpy, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout) :: tdiag
      real, dimension(:,:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxv = size(f,2)
      nypmx = size(f,3); nzpmx = size(f,4)
      if (nrec < 0) return
! check if required size of buffer has increased
      if (szfvg < ndim*nxv*nypmx*nzpmx) then
         if (szfvg > 0) deallocate(fvg)
! allocate new buffer
         allocate(fvg(ndim,nxv,nypmx,nzpmx))
         szfvg = ndim*nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPVREAD32(f,fvg,nx,ny,nz,kyp,kzp,nvpy,ndim,nxv,nypmx,nzpmx,  &
     &iunit,nrec,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpvread3 file read error: irc=', irc
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrpart3(part,tdiag,npp,iunit,iscr)
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
      call PPWRPART3(part,npp,idimp,npmax,iunit,iscr)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdpart3(part,tdiag,npp,iunit,iscr,irc)
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
      call PPRDPART3(part,npp,idimp,npmax,iunit,iscr,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrdata3(f,tdiag,iunit)
! collects distributed real scalar field data and writes to a fortran
! unformatted file
      implicit none
      integer, intent(in) :: iunit
      real, intent(inout) :: tdiag
      real, intent(in), dimension(:,:,:) :: f
! local data
      integer :: nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
! check if required size of buffer has increased
      if (szfg < nxv*nypmx*nzpmx) then
         if (szfg > 0) deallocate(fg)
! allocate new buffer
         allocate(fg(nxv,nypmx,nzpmx))
         szfg = nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRDATA3(f,fg,nxv,nypmx,nzpmx,iunit)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprddata3(f,tdiag,iunit,irc)
! reads real scalar field data from a fortran unformatted file and
! distributes it
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: irc
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
! check if required size of buffer has increased
      if (szfg < nxv*nypmx*nzpmx) then
         if (szfg > 0) deallocate(fg)
! allocate new buffer
         allocate(fg(nxv,nypmx,nzpmx))
         szfg = nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRDDATA3(f,fg,nxv,nypmx,nzpmx,iunit,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrvdata3(f,tdiag,iunit)
! collects distributed real scalar field data and writes to a fortran
! unformatted file
      implicit none
      integer, intent(in) :: iunit
      real, intent(inout) :: tdiag
      real, intent(in), dimension(:,:,:,:) :: f
! local data
      integer :: ndim, nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxv = size(f,2)
      nypmx = size(f,3); nzpmx = size(f,4)
! check if required size of buffer has increased
      if (szfvg < ndim*nxv*nypmx*nzpmx) then
         if (szfvg > 0) deallocate(fvg)
! allocate new buffer
         allocate(fvg(ndim,nxv,nypmx,nzpmx))
         szfvg = ndim*nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRVDATA3(f,fvg,ndim,nxv,nypmx,nzpmx,iunit)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdvdata3(f,tdiag,iunit,irc)
! reads real scalar field data from a fortran unformatted file and
! distributes it
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: irc
      real, intent(inout) :: tdiag
      real, dimension(:,:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxv = size(f,2)
      nypmx = size(f,3); nzpmx = size(f,4)
! check if required size of buffer has increased
      if (szfvg < ndim*nxv*nypmx*nzpmx) then
         if (szfvg > 0) deallocate(fvg)
! allocate new buffer
         allocate(fvg(ndim,nxv,nypmx,nzpmx))
         szfvg = ndim*nxv*nypmx*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRDVDATA3(f,fvg,ndim,nxv,nypmx,nzpmx,iunit,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrvcdata3(f,tdiag,iunit)
! collects distributed complex vector field data and writes to a fortran
! unformatted file
      implicit none
      integer, intent(in) :: iunit
      real, intent(inout) :: tdiag
      complex, intent(in), dimension(:,:,:,:) :: f
! local data
      integer :: ndim, nzv, kxypd, kyzpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nzv = size(f,2)
      kxypd = size(f,3); kyzpd = size(f,4)
! check if required size of buffer has increased
      if (szfvgc < ndim*nzv*kxypd*kyzpd) then
         if (szfvgc > 0) deallocate(fvgc)
! allocate new buffer
         allocate(fvgc(ndim,nzv,kxypd,kyzpd))
         szfvgc = ndim*nzv*kxypd*kyzpd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRVCDATA3(f,fvgc,ndim,nzv,kxypd,kyzpd,iunit)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdvcdata3(f,tdiag,iunit,irc)
! reads complex vector field data from a fortran unformatted file and
! distributes it
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: irc
      real, intent(inout) :: tdiag
      complex, dimension(:,:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nzv, kxypd, kyzpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nzv = size(f,2)
      kxypd = size(f,3); kyzpd = size(f,4)
! check if required size of buffer has increased
      if (szfvgc < ndim*nzv*kxypd*kyzpd) then
         if (szfvgc > 0) deallocate(fvgc)
! allocate new buffer
         allocate(fvgc(ndim,nzv,kxypd,kyzpd))
         szfvgc = ndim*nzv*kxypd*kyzpd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRDVCDATA3(f,fvgc,ndim,nzv,kxypd,kyzpd,iunit,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppartt3(partt,tdiag,numtp,nvpy,nvpz,irc)
! collects distributed test particle data
      implicit none
      integer, intent(in) :: numtp, nvpy, nvpz
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
      call PPARTT3(partt,numtp,nvpy,nvpz,idimp,nprobt,irc)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! check for errors
      if (irc /= 0) write (*,*) 'mppartt3: irc=', irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpadjfvs3(fvs,tdiag,noff,nyzp,nmv,mvy,mvz,nyb,nvpy,    &
     &nzbmx)
! adjusts 3d velocity distribution in different regions of space
      implicit none
      integer, intent(in) :: nmv, mvy, mvz, nyb, nvpy, nzbmx
      real, intent(inout) :: tdiag
      real, dimension(:,:,:,:,:), intent(inout) :: fvs
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: nmvf, nxb, nybmx, nzb, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nmvf = size(fvs,1)
      nxb = size(fvs,3); nybmx = size(fvs,4) - 1; nzb = size(fvs,5) - 1
      idds = size(noff,1)
! check if required size of buffers has increased
      if (szgvs < 6*nmvf*nxb*(nzbmx+1)) then
         if (szgvs > 0) deallocate(gvs)
! allocate new buffer
         allocate(gvs(nmvf,3,nxb,nzbmx+1,2))
         szgvs = 6*nmvf*nxb*(nzbmx+1)
      endif
      if (szhvs < 3*nmvf*nxb*nyb) then
         if (szhvs > 0) deallocate(hvs)
! allocate new buffer
         allocate(hvs(nmvf,3,nxb,nyb))
         szhvs = 3*nmvf*nxb*nyb
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPADJFVS3(fvs,gvs,hvs,noff,nyzp,nmv,mvy,mvz,nvpy,nxb,nyb,nzb,&
     &nybmx,nzbmx,nmvf,idds)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrncomp3(nyp,nzp,nvpy,nvpz,iunit)
! this subroutine collects distributed non-uniform partition information
! and writes to a fortran unformatted file
      implicit none
      integer, intent(in) :: nyp, nzp, nvpy, nvpz, iunit
! call low level procedure
      call PPWRNCOMP3(nyp,nzp,nvpy,nvpz,iunit)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrfvsdata3(fvs,tdiag,nyp,nzp,nzpmx,iunit)
! collects distributed 5d real vector non-uniform data and writes to a
! fortran unformatted file
      implicit none
      integer, intent(in) :: iunit, nyp, nzp, nzpmx
      real, intent(inout) :: tdiag
      real, intent(in), dimension(:,:,:,:,:) :: fvs
! local data
      integer :: nndim, nxv, nypmx1, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nndim = size(fvs,1)*size(fvs,2); nxv = size(fvs,3)
      nypmx1 = size(fvs,4); nypmx = nypmx1 - 1
! check if required size of buffer has increased
      if (szfvg < nndim*nxv*nypmx1*nzpmx) then
         if (szfvg > 0) deallocate(fvg)
! allocate new buffer
         allocate(fvg(nndim,nxv,nypmx1,nzpmx))
         szfvg = nndim*nxv*nypmx1*nzpmx
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRVNDATA3(fvs,fvg,nndim,nxv,nyp,nzp,nypmx,nzpmx,iunit)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
      end module
