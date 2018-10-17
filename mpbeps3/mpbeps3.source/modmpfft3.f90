!-----------------------------------------------------------------------
!
      module mfft3
!
! Fortran90 wrappers to 3d MPI/OpenMP/Vector PIC library libvmpfft3.f
! mpfft3_init calculates tables needed by 3d FFTs
!             calls WPFFT32RINIT
! mpfft3r wrapper function for scalar 3d real/complex FFT
!         calls WPPFFT32RVVM
! mpfft3rn wrapper function for vector 3d real/complex FFT
!          calls WPPFFT32RVM3 or WPPFFT32RVMN
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 16, 2018
!
      use libvmpfft3_h
      implicit none
!
! ntpose = (0,1) = (no,yes) input, output data are transposed
      integer, parameter :: ntpose = 1
!
! gs = complex scratch array for intermediate fft in y
      complex, dimension(:,:,:,:), allocatable :: gs
      integer :: szgs = -1
! bs/br = complex send/receive buffers for data transpose
      complex, dimension(:,:,:), allocatable :: bs, br
      integer :: szbuf = -1
! ss = scratch array for WPPFFT2RMN
      complex, dimension(:,:), allocatable :: ss
      integer :: szss = -1
      save
!
      private :: ntpose, gs, bs, br, ss, szgs, szbuf, szss
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpfft3_init(mixup,sct,indx,indy,indz)
! calculates tables needed by 3d FFTs
      implicit none
      integer, intent(in) :: indx, indy, indz
      integer, dimension(:), intent(inout) :: mixup
      complex, dimension(:), intent(inout) :: sct
! local data
      integer :: nxhyzd, nxyzhd
! extract dimensions
      nxhyzd = size(mixup,1); nxyzhd = size(sct,1)
! call low level procedure
      call WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfft3r(f,h,isign,mixup,sct,tfft,indx,indy,indz,kstrt, &
     &nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! wrapper function for scalar 3d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy, nvpz
      integer, intent(in) :: kxyp, kyp, kyzp, kzp
      real, dimension(:,:,:), intent(inout) :: f
      complex, dimension(:,:,:), intent(inout) :: h
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: nxvh, kypd, kzpd, nyv, kxypd, nzv, kyzpd, kzyp
      integer :: nxhyzd, nxyzhd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      nxvh = size(f,1)/2; kypd = size(f,2); kzpd = size(f,3)
      nzv = size(h,1); kxypd = size(h,2); kyzpd = size(h,3)
      nxhyzd = size(mixup,1); nxyzhd = size(sct,1)
! calculate sizes
      nyv = 2**indy; kzyp = max(kyzp,kyp)
! check if required size of buffers has increased
      if (szgs < nyv*kxypd*kzpd) then
         if (szgs > 0) deallocate(gs)
! allocate new buffer
         allocate(gs(1,nyv,kxypd,kzpd))
         szgs = nyv*kxypd*kzpd
      endif
      if (szbuf < kxyp*kzyp*kzp) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(1,kxyp*kzyp,kzp),br(1,kxyp*kzyp,kzp))
         szbuf = kxyp*kzyp*kzp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call WPPFFT32RVM(f,gs,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy&
     &,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,  &
     &kyzpd,kzpd,kzyp,nxhyzd,nxyzhd)
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfft3rn(f,h,isign,mixup,sct,tfft,indx,indy,indz,kstrt,&
     &nvpy,nvpz,kxyp,kyp,kyzp,kzp)
! wrapper function for n component vector 3d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy, nvpz
      integer, intent(in) :: kxyp, kyp, kyzp, kzp
      real, dimension(:,:,:,:), intent(inout) :: f
      complex, dimension(:,:,:,:), intent(inout) :: h
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: ndim, nxvh, kypd, kzpd, nyv, kxypd, nzv, kyzpd, kzyp
      integer :: nxhyzd, nxyzhd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      ndim = size(f,1); nxvh = size(f,2)/2
      kypd = size(f,3); kzpd = size(f,4)
      nzv = size(h,2); kxypd = size(h,3); kyzpd = size(h,4)
      nxhyzd = size(mixup,1); nxyzhd = size(sct,1)
! calculate sizes
      nyv = 2**indy; kzyp = max(kyzp,kyp)
! check if required size of buffers has increased
      if (szgs < ndim*nyv*kxypd*kzpd) then
         if (szgs > 0) deallocate(gs)
! allocate new buffer
         allocate(gs(ndim,nyv,kxypd,kzpd))
         szgs = ndim*nyv*kxypd*kzpd
      endif
      if (szbuf < ndim*kxyp*kzyp*kzp) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(ndim,kxyp*kzyp,kzp),br(ndim,kxyp*kzyp,kzp))
         szbuf = ndim*kxyp*kzyp*kzp
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call WPPFFT32RVM3(f,gs,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,   &
     &kypd,kyzpd,kzpd,kzyp,nxhyzd,nxyzhd)
      case default
! check if required size of buffer has increased
         if (szss < ndim*nxvh*kzpd) then
            if (szss > 0) deallocate(ss)
! allocate new buffer
            allocate(ss(ndim*nxvh,kzpd))
            szss = ndim*nxvh*kzpd
         endif
         call WPPFFT32RVMN(f,gs,h,bs,br,ss,isign,ntpose,mixup,sct,ttp,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,    &
     &kxypd,kypd,kyzpd,kzpd,kzyp,ndim,nxhyzd,nxyzhd)
      end select
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
      end module
