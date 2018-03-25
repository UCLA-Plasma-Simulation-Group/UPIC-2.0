!-----------------------------------------------------------------------
!
      module pgraf2
!
! Fortran90 interface to 2d PIC Fortran77 library plibgks2.f
! open_pgraphs open graphics device
!              calls GROPEN and SETNPLT
! close_pgraphs close graphics device
!               calls GRCLOSE
! set_ppalit selects one from three available palettes
!           calls STPALIT
! pdscaler2 displays 2d parallel scalar field in real space
!           calls PCARPET or PCONTUR
! pdvector2 displays 2d parallel vector field in real space
!           calls PCARPET or PCONTUR
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 9, 2017
!
      use plibgks2, only: IPLTCOMM, PGRCLOSE, PCARPET, PCONTUR, PGRASP23
      implicit none
!
! pfs = scratch array for scalar displays
      real, dimension(:,:), allocatable :: pfs
      integer :: szpfs = -1
! plf = scratch array for scalar displays
      integer, dimension(:,:), allocatable :: plf
      integer :: szplf = -1
! pfvs = scratch array for vector displays
      real, dimension(:,:), allocatable :: pfvs
      integer :: szpfvs = -1
      save
!
      private :: pfs, szpfs, plf, szplf
!
      contains
!
!-----------------------------------------------------------------------
      function open_pgraphs(nplot) result(irc)
! open graphics device
      integer, intent(in) :: nplot
      integer :: irc
      call GROPEN
      call SETNPLT(nplot,irc)
      end function
!
!-----------------------------------------------------------------------
      subroutine close_pgraphs
! close graphics device
      call PGRCLOSE
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine set_ppalit(idpal)
! selects one from three available palettes
! idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
      implicit none
      integer, intent(in) :: idpal
      call STPALIT(idpal)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pdscaler2(f,nyp,nvp,label,itime,isc,ist,idt,nx,ny,irc)
! displays 2d parallel scalar field in real space
! f = 2d parallel scalar field in real space
! nyp = number of primary gridpoints in field partition
! nvp = number of real or virtual processors requested
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! the range of values of f are given by fmax and fmin.
! if ist = 0, then fmax = 2**isc and fmin = -2**isc.
! if ist = 1, then fmax = 2**isc and fmin = 0.
! if ist = -1, then fmax = 0 and fmin = -2**isc.
! if ist = 2, then fmax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! idt = (1,2,3) = display (color map,contour plot,both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nyp, nvp, itime, isc, ist, idt, nx, ny
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: f
! local data
! ntc = number of valid colors, should be power of 2, <= 256
! nc = number of contour lines
      integer :: ntc = 16, nc = 16
      integer :: nxv, nypmx, lx
      character(len=12) :: lbl
   91 format(' T = ',i7)
      nxv = size(f,1); nypmx = size(f,2)
      lx = nx
! plot guard cells if present
      if ((lx+1) <= nxv) lx = lx + 1
! check if required size of buffer has increased
      if (szpfs < nxv*nypmx) then
         if (szpfs > 0) deallocate(pfs)
! allocate new buffer
         allocate(pfs(nxv,nypmx))
         szpfs = nxv*nypmx
      endif
      write (lbl,91) itime
! color map plot for all values
      if (idt /= 2) then
         call PCARPET(f,pfs,nyp,nvp,label,isc,ist,lx,ny,nxv,nypmx,lbl,  &
     &ntc,irc)
      endif
! contour map for all values
      if (idt /= 1) then
! check if required size of buffer has increased
         if (szplf < nxv*(nypmx+1)) then
            if (szplf > 0) deallocate(plf)
! allocate new buffer
            allocate(plf(nxv,nypmx+1))
            szplf = nxv*(nypmx+1)
         endif
         call PCONTUR(f,pfs,plf,nyp,nvp,label,isc,ist,lx,ny,nxv,nypmx,  &
     &lbl,nc,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pdvector2(fv,nyp,nvp,label,itime,isc,ist,idt,idm,nx,ny,&
     &irc)
! displays 2d parallel vector field in real space
! fv = 2d parallel vector field in real space
! nyp = number of primary gridpoints in field partition
! nvp = number of real or virtual processors requested
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of fv
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! ist = flag for choosing positive and/or negative values
! the range of values of f are given by fmax and fmin.
! if ist = 0, then fmax = 2**isc and fmin = -2**isc.
! if ist = 1, then fmax = 2**isc and fmin = 0.
! if ist = -1, then fmax = 0 and fmin = -2**isc.
! if ist = 2, then fmax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! idt = (1,2,3) = display (color map,contour plot,both)
! idm = (1,2,3) = display (components,sum(abs),both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nyp, nvp, itime, isc, ist, idt, idm, nx, ny
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: fv
! local data
      integer :: i, j, k, lx, nxv, nypmx
      real :: sum1
      character(len=2) :: c
! calculate array size with guard cells
      nxv = size(fv,2); nypmx = size(fv,3)
      lx = nx
! plot guard cells if present
      if ((lx+1) <= nxv) lx = lx + 1
! check if required size of buffer has increased
      if (szpfvs < nxv*nypmx) then
         if (szpfvs > 0) deallocate(pfvs)
! allocate new buffer
         allocate(pfvs(nxv,nypmx))
         szpfvs = nxv*nypmx
      endif
! display components
      if (idm /= 2) then
         do i = 1, size(fv,1)
            do k = 1, nypmx
            do j = 1, lx
               pfvs(j,k) = fv(i,j,k)
            enddo
            enddo
            if (i==1) then
               c = ':X'
            else if (i==2) then
               c = ':Y'
            else if (i==3) then
               c = ':Z'
            else
               write (c,'(":",i1)') i
            endif
! display i component
            call pdscaler2(pfvs,nyp,nvp,label//c,itime,isc,ist,idt,nx,ny&
     &,irc)
            if (irc /= 0) exit
         enddo
      endif
! display sum of absolute values
      if (idm /= 1) then
         do k = 1, nypmx
         do j = 1, lx
            sum1 = 0.0
            do i = 1, size(fv,1)
            sum1 = sum1 + fv(i,j,k)**2
            enddo
            pfvs(j,k) = sqrt(sum1)
         enddo
         enddo
! display amplitude
         call pdscaler2(pfvs,nyp,nvp,label,itime,isc,ist,idt,nx,ny,irc)
      endif
      end subroutine
!
      end module
