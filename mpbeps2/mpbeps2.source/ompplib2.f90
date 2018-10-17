!-----------------------------------------------------------------------
!
      module ompplib2
!
! ompmove2 reorder particles by tile with OpenMP and MPI
!          calls mporder2a, mporderf2a and mporder2b
!          as well as mpmove2, mprsncl2
! wmpfft2r performs 2d real/complex FFT for scalar data,
!          moving data between uniform and non-uniform partitions
!          calls mpfmove2 and mpfft2r
! wmpfft2rn performs 2d real/complex FFT for n component vector data,
!           moving data between uniform and non-uniform partitions
!           calls mpfnmove2 and mpfft2rn
! wmpcguard2 copy scalar guard cells to local and remote partitions
! wmpncguard2 copy vector guard cells to local and remote partitions
!             calls mpncguard2, mpcguard2x
! wmpaguard2 add scalar guard cells from local and remote partitions
!            calls mpaguard2x, mpnaguard2
! wmpnacguard2 add vector guard cells from local and remote partitions
!              calls mpacguard2x, mpnacguard2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: august 1, 2018
!
      use msort2
      use mfft2
      use mgard2
      use mppmod2, only: mpmove2, mpfmove2, mpfnmove2, mpcguard2,       &
     &mpncguard2, mpnaguard2, mpnacguard2, mpimax
      implicit none
!
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), allocatable :: ppbuff
      integer :: szpbuf = -1
! sbufl/sbufr = particle buffers sent to nearby processors
! rbufl/rbufr = particle buffers received from nearby processors
      real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer :: szbufs = -1
! ncll/nclr/mcll/mclr = number offsets send/received from processors
      integer, dimension(:,:), allocatable :: ncll, nclr, mcll, mclr
      integer :: sznbufs = -1
! msg = mpi message buffer
      integer, dimension(1) :: msg
      save
!
!     private :: ppbuff, szpbuf
      private :: szpbuf
      private :: sbufl, sbufr, rbufl, rbufr, szbufs
      private :: ncll, nclr, mcll, mclr, sznbufs
!
      contains
!
!-----------------------------------------------------------------------
      subroutine ompmove2(ppart,kpic,ncl,ihole,noff,nyp,xtras,tsort,tmov&
     &,kstrt,nvp,nx,ny,mx,my,npbmx,nbmax,mx1,popt,plist,irc2)
! reorder particles by tile with OpenMP and MPI
! plist = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: kstrt, nvp, nx, ny, mx, my, popt
      integer, intent(inout) :: npbmx, nbmax
      integer, intent(in) :: mx1
      integer, intent(in) :: noff, nyp
      real, intent(inout) :: tsort, tmov
      real, intent(in) :: xtras
      logical, intent(in) :: plist
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: nter = 10
      integer :: iter
      integer :: idimp, nppmx, mxyp1, ntmax
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(ppart,3)
      ntmax = size(ihole,2) - 1
! check if required size of buffer has increased
      if (szpbuf < idimp*npbmx*mxyp1) then
         if (szpbuf > 0) deallocate(ppbuff)
! allocate new buffer
         allocate(ppbuff(idimp,npbmx,mxyp1))
         szpbuf = idimp*npbmx*mxyp1
      endif
! check if required size of buffers has increased
      if (szbufs < idimp*nbmax) then
         if (szbufs > 0) deallocate(sbufl,sbufr,rbufl,rbufr)
! allocate new buffers
         allocate(sbufl(idimp,nbmax),sbufr(idimp,nbmax))
         allocate(rbufl(idimp,nbmax),rbufr(idimp,nbmax))
         szbufs = idimp*nbmax
      endif
! check if required size of buffers has increased
      if (sznbufs < 3*mx1) then
         if (sznbufs > 0) deallocate(ncll,nclr,mcll,mclr)
! allocate new buffers
         allocate(ncll(3,mx1),nclr(3,mx1),mcll(3,mx1),mclr(3,mx1))
         sznbufs = 3*mx1
      endif
!
! ppart overflow
! second part of particle reorder on x and y cell with mx, my tiles:
      if (irc2(1)==4) then
         irc2 = 0
         call mporder2b(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,   &
     &mclr,tsort,kstrt,nx,ny,popt,irc2)
         return
! first part of particle reorder on x and y cell with mx, my tiles:
! list of particles leaving tile already calculated by push
      else if (plist) then
! updates: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
         call mporderf2a(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr,  &
     &tsort,kstrt,popt,irc2)
! calculate list of particles leaving tile
      else
! updates ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
         call mporder2a(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,   &
     &nclr,noff,nyp,tsort,kstrt,nx,ny,mx,my,popt,irc2)
      endif
!
      iter = 0
      do while (irc2(1) /= 0)
         iter = iter + 1
         if (iter > nter) then
            write (*,*) kstrt, 'ompmove2: iteration exceeded'
            call PPABORT()
            stop
         endif
! ihole overflow
         if (irc2(1)==1) then
            return
! ppbuff overflow
         else if (irc2(1)==2) then
            npbmx = (1.0 + xtras)*irc2(2)
            deallocate(ppbuff)
! allocate new buffer
            allocate(ppbuff(idimp,npbmx,mxyp1))
            szpbuf = idimp*npbmx*mxyp1
! restores initial values of ncl
            call mprsncl2(ncl,tsort)
            irc2 = 0
            call mporderf2a(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr&
     &,tsort,kstrt,popt,irc2)
! sbufr/sbufl overflow
         else if (irc2(1)==3) then
            nbmax = (1.0 + xtras)*irc2(2)
            deallocate(sbufl,sbufr)
! allocate new buffer
            allocate(sbufl(idimp,nbmax),sbufr(idimp,nbmax))
            szbufs = idimp*nbmax
            irc2 = 0
            call mporderf2af(ppbuff,sbufl,sbufr,ncl,ncll,nclr,tsort,    &
     &kstrt,irc2)
         endif
      enddo
!
! verify that mpi send/receive buffers are large enough
      msg = nbmax
      call mpimax(msg,tmov)
! rbufr/rbufl overflow
      if ((msg(1) > size(rbufl,2)).or.(msg(1) > size(rbufr,2))) then
         deallocate(rbufl,rbufr)
         allocate(rbufl(idimp,msg(1)),rbufr(idimp,msg(1)))
! reallocate sbufl/sbufr
         if (msg(1) > nbmax) then
            rbufl(:,1:nbmax) = sbufl
            rbufr(:,1:nbmax) = sbufr
            nbmax = msg(1)
            deallocate(sbufl,sbufr)
            allocate(sbufl(idimp,nbmax),sbufr(idimp,nbmax))
            szbufs = idimp*nbmax
            sbufl = rbufl
            sbufr = rbufr
         endif
      endif
!
! move particles into appropriate spatial regions with MPI:
! updates rbufr, rbufl, mcll, mclr
      call mpmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,tmov,    &
     &kstrt,nvp)
!
! second part of particle reorder on x and y cell with mx, my tiles:
! updates ppart, kpic
      call mporder2b(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr, &
     &tsort,kstrt,nx,ny,popt,irc2)
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpfft2r(f,g,noff,nyp,isign,mixup,sct,tfft,tfmov,indx, &
     &indy,kstrt,nvp,kyp,ny,mter,ierr)
! performs 2d real/complex FFT for scalar data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kyp, ny
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:), intent(inout) :: f
      complex, dimension(:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
! inverse fourier transform: from real to complex
      if (isign < 0) then
! moves scalar grids from non-uniform to uniform partition
         call mpfmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,ierr&
     &)
! wrapper function for scalar 2d real/complex FFT
         call mpfft2r(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,kyp)
! forward fourier transform: from complex to real
      else if (isign > 0) then
! wrapper function for scalar 2d real/complex FFT
         call mpfft2r(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,kyp)
! moves scalar grids from uniform to non-uniform partition
         call mpfmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,ierr&
     &)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpfft2rn(f,g,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
     &indy,kstrt,nvp,kyp,ny,mter,ierr)
! performs 2d real/complex FFT for n component vector data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kyp, ny
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:), intent(inout) :: f
      complex, dimension(:,:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
! inverse fourier transform
      if (isign < 0) then
! moves vector grids from non-uniform to uniform partition
         call mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,   &
     &ierr)
! wrapper function for n component vector 2d real/complex FFT
         call mpfft2rn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,kyp)
! forward fourier transform
      else if (isign > 0) then
! wrapper function for n component vector 2d real/complex FFT
         call mpfft2rn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,kyp)
! moves vector grids from uniform to non-uniform partition
         call mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,   &
     &ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpcguard2(f,nyp,tguard,nx,kstrt,nvp)
! copy scalar guard cells to local and remote partitions
      implicit none
      integer, intent(in) :: nyp, nx, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: f
! copies data to guard cells in non-uniform partitions
      call mpcguard2(f,nyp,tguard,kstrt,nvp)
! replicates local periodic scalar field
      call mpdguard2x(f,nyp,tguard,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpncguard2(f,nyp,tguard,nx,kstrt,nvp)
! copy vector guard cells to local and remote partitions
      implicit none
      integer, intent(in) :: nyp, nx, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
! copies data to guard cells in non-uniform partitions
      call mpncguard2(f,nyp,tguard,kstrt,nvp)
! replicates local periodic vector field
      call mpcguard2x(f,nyp,tguard,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpaguard2(f,nyp,tguard,nx,kstrt,nvp)
! add scalar guard cells from local and remote partitions
      implicit none
      integer, intent(in) :: nyp, nx, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: f
! accumulates local periodic scalar field
      call mpaguard2x(f,nyp,tguard,nx)
! adds scalar data from guard cells in non-uniform partitions
      call mpnaguard2(f,nyp,tguard,nx,kstrt,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpnacguard2(f,nyp,tguard,nx,kstrt,nvp)
! add vector guard cells from local and remote partitions
      implicit none
      integer, intent(in) :: nyp, nx, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
! accumulates local periodic vector field
      call mpacguard2x(f,nyp,tguard,nx)
! adds vector data from guard cells in non-uniform partitions
      call mpnacguard2(f,nyp,tguard,nx,kstrt,nvp)
      end subroutine
!
      end module
