!-----------------------------------------------------------------------
!
      module omppdflib2
!
! wmpfsst2r performs 2d real sine-sine transform for scalar data,
!           moving data between uniform and non-uniform partitions
!           calls mpfmove2 and mpfsst2r
! wmpfcst2rn performs 2d real cosine-sine transforms for n component
!            vector data, moving data between uniform and non-uniform
!            partitions
!            calls mpfnmove2 and mpfcst2rn
! wmpfsct2rn performs 2d real sine-cosine transforms for n component
!            vector data, moving data between uniform and non-uniform
!            partitions
!            calls mpfnmove2 and mpfsct2rn
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: april 28, 2017
!
      use modmpfsct2
      use mppmod2, only: mpfmove2, mpfnmove2
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine wmpfsst2r(f,g,noff,nyp,isign,mixup,sctd,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kxp2,kyp,ny,mter,ierr)
! performs 2d real sine-sine transform for scalar data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      integer, intent(in) :: ny
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:), intent(inout) :: f, g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! inverse fast sine-sine transform transform
      if (isign < 0) then
! moves scalar grids from non-uniform to uniform partition
         call mpfmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,ierr&
     &)
! wrapper function for scalar 2d real sine-sine transform
         call mpfsst2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,   &
     &kxp2,kyp)
! forward fast sine-sine transform transform
      else if (isign > 0) then
! wrapper function for scalar 2d real sine-sine transform
         call mpfsst2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,   &
     &kxp2,kyp)
! moves scalar grids from uniform to non-uniform partition
         call mpfmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,ierr&
     &)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpfcst2rn(f,g,noff,nyp,isign,mixup,sctd,tfft,tfmov,   &
     &indx,indy,kstrt,nvp,kxp2,kyp,ny,mter,ierr)
! performs 2d real cosine-sine transforms for n component vector data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      integer, intent(in) :: ny
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:), intent(inout) :: f, g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! inverse fourier transform
      if (isign < 0) then
! moves vector grids from non-uniform to uniform partition
         call mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,   &
     &ierr)
! wrapper function for n component vector 2d real cosine-sine transforms
         call mpfcst2rn(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,  &
     &kxp2,kyp)
! forward fourier transform
      else if (isign > 0) then
! wrapper function for n component vector 2d real cosine-sine transforms
         call mpfcst2rn(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,  &
     &kxp2,kyp)
! moves vector grids from uniform to non-uniform partition
         call mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,   &
     &ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpfsct2rn(f,g,noff,nyp,isign,mixup,sctd,tfft,tfmov,   &
     &indx,indy,kstrt,nvp,kxp2,kyp,ny,mter,ierr)
! performs 2d real sine-cosine transforms for n component vector data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      integer, intent(in) :: ny
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:), intent(inout) :: f, g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! inverse fourier transform
      if (isign < 0) then
! moves vector grids from non-uniform to uniform partition
         call mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,   &
     &ierr)
! wrapper function for n component vector 2d real sine-cosine transforms
         call mpfsct2rn(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,  &
     &kxp2,kyp)
! forward fourier transform
      else if (isign > 0) then
! wrapper function for n component vector 2d real sine-cosine transforms
         call mpfsct2rn(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,  &
     &kxp2,kyp)
! moves vector grids from uniform to non-uniform partition
         call mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,   &
     &ierr)
      endif
      end subroutine
!
      end module
