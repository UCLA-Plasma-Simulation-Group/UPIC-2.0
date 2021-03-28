!-----------------------------------------------------------------------
!
      module graf2
!
! Fortran90 interface to 2d PIC Fortran77 library libgks2.f
! open_graphs2 open graphics device
!              calls GROPEN and SETNPLT
! close_graphs2 close graphics device
!               calls GRCLOSE
! set_palit selects one from three available palettes
!           calls STPALIT
! dscaler2 displays 2d scalar field in real space.
!          calls CARPET or CONTUR
! dscalerl2 displays color map of 2d scalar field in real space
!           calls CARPETL
! dvector2 displays 2d vector field in real space.
!          calls CARPET or CONTUR
! dvectorl2 displays color map of 2d vector field in real space.
!           calls CARPETL
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 6, 2017
!
      use libgraf1_h
      use libgraf2_h
      implicit none
!
! fs = scratch array for vector displays
      real, dimension(:,:), allocatable :: fs
      integer :: szfs = -1
      save
!
      private :: fs, szfs
!
      contains
!
!-----------------------------------------------------------------------
      function open_graphs2(nplot) result(irc)
      implicit none
! open graphics device
      integer, intent(in) :: nplot
      integer :: irc
      call GROPEN
      call SETNPLT(nplot,irc)
      end function
!
!-----------------------------------------------------------------------
      subroutine close_graphs2
! close graphics device
      call GRCLOSE
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine set_palit(idpal)
! selects one from three available palettes
! idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
      implicit none
      integer, intent(in) :: idpal
      call STPALIT(idpal)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dscaler2(f,label,itime,isc,ist,idt,nx,ny,irc)
! displays 2d scalar field in real space
! f = 2d scalar field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
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
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, ist, idt, nx, ny
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: f
! local data
      integer :: nxv, lx, ly
! ntc = number of valid colors, should be power of 2, <= 256
! nc = number of contour lines
      integer :: ntc = 16, nc = 16
      character(len=12) :: lbl
      integer, dimension(size(f,1),size(f,2)) :: lf
   91 format(' T = ',i7)
      nxv = size(f,1)
      lx = nx; ly = ny
! plot guard cells if present
      if ((lx+1) <= nxv) lx = lx + 1
      if ((ly+1) <= size(f,2)) ly = ly + 1
      write (lbl,91) itime
! color map plot for all values
      if (idt /= 2) then
         call CARPET(f,label,isc,ist,lx,ly,nxv,lbl,ntc,irc)
      endif
! contour map for all values
      if (idt /= 1) then
         call CONTUR(f,lf,label,isc,ist,lx,ly,nxv,lbl,nc,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dscalerl2(f,label,xmin,xmax,ymin,ymax,itime,isc,ist,nx,&
     &ny,irc)
! displays color map of 2d scalar field in real space
! f = 2d scalar field in real space
! label = field label
! xmin/xmax = numerical labels for x axis
! ymin/ymax = numerical labels for y axis
! itime = current time step
! isc = power of 2 scale of range of values of f
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
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, ist, nx, ny
      integer, intent(inout) :: irc
      real, intent(in) :: xmin, xmax, ymin, ymax
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: f
! local data
      integer :: nxv, lx, ly
! ntc = number of valid colors, should be power of 2, <= 256
      integer :: ntc = 16
      character(len=12) :: lbl
   91 format(' T = ',i7)
      nxv = size(f,1)
      lx = nx; ly = ny
! plot guard cells if present
      if ((lx+1) <= nxv) lx = lx + 1
      if ((ly+1) <= size(f,2)) ly = ly + 1
      write (lbl,91) itime
! color map plot for all values
      call CARPETL(f,label,xmin,xmax,ymin,ymax,isc,ist,lx,ly,nxv,lbl,ntc&
     &,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dvector2(fv,label,itime,isc,ist,idt,idm,nx,ny,irc)
! displays 2d vector field in real space
! fv = 2d vector field in real space
! fs = 2d scratch field array
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
      integer :: i, j, k, lx, ly
      real :: sum1
      character(len=2) :: c
! calculate array size with guard cells
      lx = nx + 1; ly = ny + 1
      lx = min(lx,size(fv,2))
      ly = min(ly,size(fv,3))
! check if required size of buffer has increased
      if (szfs < lx*ly) then
         if (szfs > 0) deallocate(fs)
! allocate new buffer
         allocate(fs(lx,ly))
         szfs = lx*ly
      endif
! display components
      if (idm /= 2) then
         do i = 1, size(fv,1)
            do k = 1, ly
            do j = 1, lx
               fs(j,k) = fv(i,j,k)
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
            call dscaler2(fs,label//c,itime,isc,ist,idt,nx,ny,irc)
         enddo
      endif
! display sum of absolute values
      if (idm /= 1) then
         do k = 1,ly
         do j = 1, lx
            sum1 = 0.0
            do i = 1, size(fv,1)
            sum1 = sum1 + fv(i,j,k)**2
            enddo
            fs(j,k) = sqrt(sum1)
         enddo
         enddo
! display amplitude
         call dscaler2(fs,label,itime,isc,ist,idt,nx,ny,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dvectorl2(fv,label,xmin,xmax,ymin,ymax,itime,isc,ist,  &
     &idm,nx,ny,irc)
! displays color map of 2d vector field in real space
! fv = 2d vector field in real space
! fs = 2d scratch field array
! label = field label
! xmin/xmax = numerical labels for x axis
! ymin/ymax = numerical labels for y axis
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
! idm = (1,2,3) = display (components,sum(abs),both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, ist, idm, nx, ny
      integer, intent(inout) :: irc
      real, intent(in) :: xmin, xmax, ymin, ymax
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: fv
! local data
      integer :: i, j, k, lx, ly
      real :: sum1
      character(len=2) :: c
! calculate array size with guard cells
      lx = nx + 1; ly = ny + 1
      lx = min(lx,size(fv,2))
      ly = min(ly,size(fv,3))
! check if required size of buffer has increased
      if (szfs < lx*ly) then
         if (szfs > 0) deallocate(fs)
! allocate new buffer
         allocate(fs(lx,ly))
         szfs = lx*ly
      endif
! display components
      if (idm /= 2) then
         do i = 1, size(fv,1)
            do k = 1, ly
            do j = 1, lx
               fs(j,k) = fv(i,j,k)
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
            call dscalerl2(fs,label//c,xmin,xmax,ymin,ymax,itime,isc,ist&
     &,nx,ny,irc)
         enddo
      endif
! display sum of absolute values
      if (idm /= 1) then
         do k = 1,ly
         do j = 1, lx
            sum1 = 0.0
            do i = 1, size(fv,1)
            sum1 = sum1 + fv(i,j,k)**2
            enddo
            fs(j,k) = sqrt(sum1)
         enddo
         enddo
! display amplitude
         call dscalerl2(fs,label,xmin,xmax,ymin,ymax,itime,isc,ist,nx,ny&
     &,irc)
      endif
      end subroutine
!
      end module
