!-----------------------------------------------------------------------
!
      module graf2
!
! Fortran90 interface to 2d PIC Fortran77 library libgks2.f
! pdscaler2 displays 2d scalar field in real space
!           calls CARPET or CONTUR
! pdvector2 displays 2d vector field in real space
!           calls CARPET or CONTUR
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: may 11, 2017
!
      implicit none
!
! lf = scratch array for scalar displays
      integer, dimension(:,:), allocatable :: lf
      integer :: szlf = -1
! fvs = scratch array for vector displays
      real, dimension(:,:), allocatable :: fvs
      integer :: szfvs = -1
      save
!
      private :: lf, szlf, fvs, szfvs
!
      contains
!
!-----------------------------------------------------------------------
      function open_graphs(nplot) result(irc)
! open graphics device
      integer, intent(in) :: nplot
      integer :: irc
      call GROPEN
      call SETNPLT(nplot,irc)
      end function
!
!-----------------------------------------------------------------------
      subroutine close_graphs
! close graphics device
      call GRCLOSE
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dscaler2(f,label,itime,isc,ist,idt,nx,ny,irc)
! displays 2d scalar field in real space
! f = 2d scalar field in real space
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
      integer, intent(in) :: itime, isc, ist, idt, nx, ny
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: f
! local data
! ntc = number of valid colors, should be power of 2, <= 256
! nc = number of contour lines
      integer :: ntc = 16, nc = 16
      integer :: nxv, nyv, lx
      character(len=12) :: lbl
   91 format(' T = ',i7)
      nxv = size(f,1); nyv = size(f,2)
      lx = nx
! plot guard cells if present
      if ((lx+1) <= nxv) lx = lx + 1
      write (lbl,91) itime
! color map plot for all values
      if (idt /= 2) then
         call CARPET(f,label,isc,ist,lx,ny,nxv,lbl,ntc,irc)
      endif
! contour map for all values
      if (idt /= 1) then
! check if required size of buffer has increased
         if (szlf < nxv*(nyv+1)) then
            if (szlf > 0) deallocate(lf)
! allocate new buffer
            allocate(lf(nxv,nyv+1))
            szlf = nxv*(nyv+1)
         endif
         call CONTUR(f,lf,label,isc,ist,lx,ny,nxv,lbl,nc,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dvector2(fv,label,itime,isc,ist,idt,idm,nx,ny,irc)
! displays 2d pvector field in real space
! fv = 2d vector field in real space
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
      integer, intent(in) :: itime, isc, ist, idt, idm, nx, ny
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: fv
! local data
      integer :: i, j, k, lx, nxv, nyv
      real :: sum1
      character(len=2) :: c
! calculate array size with guard cells
      nxv = size(fv,2); nyv = size(fv,3)
      lx = nx
! plot guard cells if present
      if ((lx+1) <= nxv) lx = lx + 1
! check if required size of buffer has increased
      if (szfvs < nxv*nyv) then
         if (szfvs > 0) deallocate(fvs)
! allocate new buffer
         allocate(fvs(nxv,nyv))
         szfvs = nxv*nyv
      endif
! display components
      if (idm /= 2) then
         do i = 1, size(fv,1)
            do k = 1, nyv
            do j = 1, lx
               fvs(j,k) = fv(i,j,k)
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
            call dscaler2(fvs,label//c,itime,isc,ist,idt,nx,ny,irc)
            if (irc /= 0) exit
         enddo
      endif
! display sum of absolute values
      if (idm /= 1) then
         do k = 1, nyv
         do j = 1, lx
            sum1 = 0.0
            do i = 1, size(fv,1)
            sum1 = sum1 + fv(i,j,k)**2
            enddo
            fvs(j,k) = sqrt(sum1)
         enddo
         enddo
! display amplitude
         call dscaler2(fvs,label,itime,isc,ist,idt,nx,ny,irc)
      endif
      end subroutine
!
      end module
