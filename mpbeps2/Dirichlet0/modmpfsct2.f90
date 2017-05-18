!-----------------------------------------------------------------------
!
      module modmpfsct2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpfsct2.f
! mpfsct2_init calculates tables needed by 2d Fast Sine/Cosine
!              Transforms
!              calls WPFST2RINIT
! mpfsst2r wrapper function for 2d scalar real sine-sine transform
!          calls WPPFSST2RM
! mpfsct2r wrapper function for 2d scalar real sine-cosine transform
!          calls WPPFSCT2RM
! mpfcst2r wrapper function for 2d scalar real cosine-sine transform
!          calls WPPFCST2RM
! mpfcct2r wrapper function for 2d scalar real cosine-cosine transform
!          calls WPPFCCT2RM
! mpfcst2rn wrapper function for 2d vector real cosine-sine
!           transforms
!           calls WPPFCST2RM2 or WPPFCST2RM3
! mpfsct2rn wrapper function for 2d vector real cosine-sine
!           transforms
!           calls WPPFSCT2RM2 or WPPFSCT2RM3
! mpfsct2rm wrapper function for 2d tensor real sine-cosine transforms
!           calls WPPFSCT2RM4 or WPPFSCT2RM22 or WPPFSST2RM23
! written by viktor k. decyk, ucla
! copyright 2017, regents of the university of california
! update: april 28, 2017
!
      use libmpfft2_h
      implicit none
!
! ntpose = (0,1) = (no,yes) input, output data are transposed
      integer, parameter :: ntpose = 1
!
! bs/br = complex send/receive buffers for data transpose
      complex, dimension(:,:), allocatable :: bs, br
      integer :: szbuf = -1
      save
!
      private :: ntpose, bs, br, szbuf
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpfsct2_init(mixup,sctd,indx,indy)
! initialize 2d real sine-cosine transforms
      implicit none
      integer, intent(in) :: indx, indy
      integer, dimension(:), intent(inout) :: mixup
      complex, dimension(:), intent(inout) :: sctd
! local data
      integer :: nxhyd, nxyd
      nxhyd = size(mixup); nxyd = size(sctd)
      call WPFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfsst2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,&
     &kxp2,kyp)
! wrapper function for 2d scalar real sine-sine transform
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      real, dimension(:,:), intent(inout) :: f
      real, dimension(:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: nxvh, kypd, nyv, kxp2d, nxhyd, nxyd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      nxvh = size(f,1)/2; kypd = size(f,2)
      nyv = size(g,1); kxp2d = size(g,2)
      nxhyd = size(mixup,1); nxyd = size(sctd,1)
! check if required size of buffers has increased
      if (szbuf < (kxp2+1)*(kyp+1)) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(kxp2+1,kyp+1),br(kxp2+1,kyp+1))
         szbuf = (kxp2+1)*(kyp+1)
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call WPPFSST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,indy,  &
     &kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfsct2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,&
     &kxp2,kyp)
! wrapper function for 2d scalar real sine-cosine transform
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      real, dimension(:,:), intent(inout) :: f
      real, dimension(:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: nxvh, kypd, nyv, kxp2d, nxhyd, nxyd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      nxvh = size(f,1)/2; kypd = size(f,2)
      nyv = size(g,1); kxp2d = size(g,2)
      nxhyd = size(mixup,1); nxyd = size(sctd,1)
! check if required size of buffers has increased
      if (szbuf < (kxp2+1)*(kyp+1)) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(kxp2+1,kyp+1),br(kxp2+1,kyp+1))
         szbuf = (kxp2+1)*(kyp+1)
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call WPPFSCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,indy,  &
     &kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfcst2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,&
     &kxp2,kyp)
! wrapper function for 2d scalar real cosine-sine transform
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      real, dimension(:,:), intent(inout) :: f
      real, dimension(:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: nxvh, kypd, nyv, kxp2d, nxhyd, nxyd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      nxvh = size(f,1)/2; kypd = size(f,2)
      nyv = size(g,1); kxp2d = size(g,2)
      nxhyd = size(mixup,1); nxyd = size(sctd,1)
! check if required size of buffers has increased
      if (szbuf < (kxp2+1)*(kyp+1)) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(kxp2+1,kyp+1),br(kxp2+1,kyp+1))
         szbuf = (kxp2+1)*(kyp+1)
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call WPPFCST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,indy,  &
     &kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfcct2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp,&
     &kxp2,kyp)
! wrapper function for 2d scalar real cosine-cosine transform
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      real, dimension(:,:), intent(inout) :: f
      real, dimension(:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: nxvh, kypd, nyv, kxp2d, nxhyd, nxyd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      nxvh = size(f,1)/2; kypd = size(f,2)
      nyv = size(g,1); kxp2d = size(g,2)
      nxhyd = size(mixup,1); nxyd = size(sctd,1)
! check if required size of buffers has increased
      if (szbuf < (kxp2+1)*(kyp+1)) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(kxp2+1,kyp+1),br(kxp2+1,kyp+1))
         szbuf = (kxp2+1)*(kyp+1)
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call WPPFCCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,indy,  &
     &kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfcst2rn(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp&
     &,kxp2,kyp)
! wrapper function for 2d vector real cosine-sine transforms
! for the electric field with dirichlet or magnetic field with neumann
! boundary conditions
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      real, dimension(:,:,:), intent(inout) :: f
      real, dimension(:,:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: ndim, nxvh, kypd, nyv, kxp2d, nxhyd, nxyd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
      nyv = size(g,2); kxp2d = size(g,3)
      nxhyd = size(mixup,1); nxyd = size(sctd,1)
! check if required size of buffers has increased
      if (szbuf < ndim*(kxp2+1)*(kyp+1)) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(ndim,(kxp2+1)*(kyp+1)),br(ndim,(kxp2+1)*(kyp+1)))
         szbuf = ndim*(kxp2+1)*(kyp+1)
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (ndim==2) then
         call WPPFCST2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,   &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
      else if (ndim==3) then
         call WPPFCST2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,   &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
      else
         if (kstrt==1) write (*,*) 'mpfcst2rn error, ndim=',ndim
      endif
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfsct2rn(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp&
     &,kxp2,kyp)
! wrapper function for 2d vector real sine-cosine transforms
! for the magnetic field with dirichlet or electric field with neumann
! boundary conditions
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      real, dimension(:,:,:), intent(inout) :: f
      real, dimension(:,:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: ndim, nxvh, kypd, nyv, kxp2d, nxhyd, nxyd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
      nyv = size(g,2); kxp2d = size(g,3)
      nxhyd = size(mixup,1); nxyd = size(sctd,1)
! check if required size of buffers has increased
      if (szbuf < ndim*(kxp2+1)*(kyp+1)) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(ndim,(kxp2+1)*(kyp+1)),br(ndim,(kxp2+1)*(kyp+1)))
         szbuf = ndim*(kxp2+1)*(kyp+1)
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (ndim==2) then
         call WPPFSCT2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,   &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
      else if (ndim==3) then
         call WPPFSCT2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,   &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
      else
         if (kstrt==1) write (*,*) 'mpfsct2rn error, ndim=',ndim
      endif
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfsct2rm(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,nvp&
     &,kxp2,kyp)
! wrapper function for 2d tensor real cosine-sine transforms
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxp2, kyp
      real, dimension(:,:,:), intent(inout) :: f
      real, dimension(:,:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sctd
      real, dimension(2), intent(inout) :: tfft
! local data
      integer :: mdim, nxvh, kypd, nyv, kxp2d, nxhyd, nxyd
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: ttp
! extract dimensions
      mdim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
      nyv = size(g,2); kxp2d = size(g,3)
      nxhyd = size(mixup,1); nxyd = size(sctd,1)
! check if required size of buffers has increased
      if (szbuf < mdim*(kxp2+1)*(kyp+1)) then
         if (szbuf > 0) deallocate(bs,br)
! allocate new buffers
         allocate(bs(mdim,(kxp2+1)*(kyp+1)),br(mdim,(kxp2+1)*(kyp+1)))
         szbuf = mdim*(kxp2+1)*(kyp+1)
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (mdim==2) then
         call WPPFSCT2RM22(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,  &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
      else if (mdim==3) then
         call WPPFSST2RM23(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,  &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
      else if (mdim==4) then
         call WPPFSCT2RM4(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,   &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
      else
         if (kstrt==1) write (*,*) 'mpfsct2rm error, mdim=',mdim
      endif
! record time
      call dtimer(dtime,itime,1)
      tfft(1) = tfft(1) + real(dtime)
      tfft(2) = tfft(2) + ttp
      end subroutine
!
      end module
