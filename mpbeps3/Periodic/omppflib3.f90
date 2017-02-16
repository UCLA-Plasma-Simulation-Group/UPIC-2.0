!-----------------------------------------------------------------------
!
      module omppflib3
!
! wmpfft3r performs 3d real/complex FFT for scalar data,
!          moving data between uniform and non-uniform partitions
!          calls mpfmove3 and mpfft3r
! wmpfft3rn performs 3d real/complex FFT for n component vector data,
!           moving data between uniform and non-uniform partitions
!           calls mpfnmove3 and mpfft3rn
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 15, 2017
!
      use modmpfft3
      use mppmod3, only: mpfmove3, mpfnmove3
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine wmpfft3r(f,h,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx,&
     &indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mter,ierr)
! performs 3d real/complex FFT for scalar data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy, nvpz
      integer, intent(in) :: kxyp, kyp, kyzp, kzp, ny, nz
      integer, intent(inout) :: ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:), intent(inout) :: f
      complex, dimension(:,:,:), intent(inout) :: h
      integer, dimension(:), intent(in) :: noff, nyzp
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
      integer, dimension(2), intent(inout) :: mter
! inverse fourier transform: from real to complex
      if (isign < 0) then
! moves scalar grids from non-uniform to uniform partition
         call mpfmove3(f,noff,nyzp,isign,tfmov,kyp,kzp,ny,nz,kstrt,nvpy,&
     &nvpz,mter,ierr)
! wrapper function for scalar 3d real/complex FFT
         call mpfft3r(f,h,isign,mixup,sct,tfft,indx,indy,indz,kstrt,nvpy&
     &,nvpz,kxyp,kyp,kyzp,kzp)
! forward fourier transform: from complex to real
      else if (isign > 0) then
! wrapper function for scalar 3d real/complex FFT
         call mpfft3r(f,h,isign,mixup,sct,tfft,indx,indy,indz,kstrt,nvpy&
     &,nvpz,kxyp,kyp,kyzp,kzp)
! moves scalar grids from uniform to non-uniform partition
         call mpfmove3(f,noff,nyzp,isign,tfmov,kyp,kzp,ny,nz,kstrt,nvpy,&
     &nvpz,mter,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpfft3rn(f,h,noff,nyzp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,indz,kstrt,nvpy,nvpz,kxyp,kyp,kyzp,kzp,ny,nz,mter,ierr)
! performs 3d real/complex FFT for n component vector data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy, nvpz
      integer, intent(in) :: kxyp, kyp, kyzp, kzp, ny, nz
      integer, intent(inout) :: ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:,:), intent(inout) :: f
      complex, dimension(:,:,:,:), intent(inout) :: h
      integer, dimension(:), intent(in) :: noff, nyzp
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
      integer, dimension(2), intent(inout) :: mter
! inverse fourier transform
      if (isign < 0) then
! moves vector grids from non-uniform to uniform partition
         call mpfnmove3(f,noff,nyzp,isign,tfmov,kyp,kzp,ny,nz,kstrt,nvpy&
     &,nvpz,mter,ierr)
! wrapper function for n component vector 3d real/complex FFT
         call mpfft3rn(f,h,isign,mixup,sct,tfft,indx,indy,indz,kstrt,   &
     &nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! forward fourier transform
      else if (isign > 0) then
! wrapper function for n component vector 3d real/complex FFT
         call mpfft3rn(f,h,isign,mixup,sct,tfft,indx,indy,indz,kstrt,   &
     &nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! moves vector grids from uniform to non-uniform partition
         call mpfnmove3(f,noff,nyzp,isign,tfmov,kyp,kzp,ny,nz,kstrt,nvpy&
     &,nvpz,mter,ierr)
      endif
      end subroutine
!
      end module
