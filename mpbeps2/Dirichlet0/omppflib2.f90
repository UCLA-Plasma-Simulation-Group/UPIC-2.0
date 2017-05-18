!-----------------------------------------------------------------------
!
      module omppflib2
!
! wmpfft2r performs 2d real/complex FFT for scalar data,
!          moving data between uniform and non-uniform partitions
!          calls mpfmove2 and mpfft2r
! wmpfft2rn performs 2d real/complex FFT for n component vector data,
!           moving data between uniform and non-uniform partitions
!           calls mpfnmove2 and mpfft2rn
! wmpcguard2 copy scalar guard cells to local and remote partitions
! wmpncguard2 copy vector guard cells to local and remote partitions
!             calls mpncguard2, mpcguard2x
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 8, 2017
!
      use modmpfft2
      use modmpgard2
      use mppmod2, only: mpfmove2, mpfnmove2, mpcguard2, mpncguard2
      implicit none
!
      contains
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
      end module
