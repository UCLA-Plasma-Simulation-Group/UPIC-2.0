!-----------------------------------------------------------------------
!
      module modmpfft2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpfft2.f
! mpfft2_init calculates tables needed by 2d FFTs
!             calls WPFFT2RINIT
! mpfft2r wrapper function for scalar 2d real/complex FFT
!         calls WPPFFT2RM
! mpfft2rn wrapper function for vector 2d real/complex FFT
!          calls WPPFFT2RM2, WPPFFT2RM3, or WPPFFT2RMN
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 11, 2016
!
      use libmpfft2_h
      implicit none
!
! ntpose = (0,1) = (no,yes) input, output data are transposed
      integer, parameter :: ntpose = 1
!
! bs/br = complex send/receive buffers for data transpose
      complex, dimension(:,:), allocatable :: bs, br
      integer :: szbuf = 0
! ss = scratch array for WPPFFT2RMN
      complex, dimension(:,:), allocatable :: ss
      integer :: szss = 0
      save
!
      private :: ntpose, bs, br, ss, szbuf, szss
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpfft2_init(mixup,sct,indx,indy)
! calculates tables needed by 2d FFTs
      implicit none
      integer, intent(in) :: indx, indy
      integer, dimension(:), intent(inout) :: mixup
      complex, dimension(:), intent(inout) :: sct
! local data
      integer :: nxhyd, nxyhd
! extract dimensions
      nxhyd = size(mixup,1); nxyhd = size(sct,1)
! call low level procedure
      call WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfft2r(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,  &
     &kyp)
! wrapper function for scalar 2d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kyp
      real, dimension(:,:), intent(inout) :: f
      complex, dimension(:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: nxvh, kypd, nyv, kxp, nxhyd, nxyhd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      nxvh = size(f,1)/2; kypd = size(f,2)
      nyv = size(g,1); kxp = size(g,2)
      nxhyd = size(mixup,1); nxyhd = size(sct,1)
! check if required size of buffers has increased
      if (szbuf < kxp*kyp) then
         if (szbuf /= 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(kxp,kyp),br(kxp,kyp))
         szbuf = kxp*kyp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call WPPFFT2RM(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,    &
     &kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfft2rn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp, &
     &kyp)
! wrapper function for n component vector 2d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kyp
      real, dimension(:,:,:), intent(inout) :: f
      complex, dimension(:,:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: ndim, nxvh, kypd, nyv, kxp, nxhyd, nxyhd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
      nyv = size(g,2); kxp = size(g,3)
      nxhyd = size(mixup,1); nxyhd = size(sct,1)
! check if required size of buffers has increased
      if (szbuf < ndim*kxp*kyp) then
         if (szbuf /= 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(ndim*kxp,kyp),br(ndim*kxp,kyp))
         szbuf = ndim*kxp*kyp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call WPPFFT2RM2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,&
     &kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
      case (3)
         call WPPFFT2RM3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,&
     &kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
      case default
! check if required size of buffer has increased
         if (szss < ndim*nxvh*kypd) then
            if (szss /= 0) deallocate(ss)
! allocate new buffer
            allocate(ss(ndim*nxvh,kypd))
            szss = ndim*nxvh*kypd
         endif
         call WPPFFT2RMN(f,g,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx,  &
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,ndim,nxhyd,nxyhd)
      end select
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
      end module
