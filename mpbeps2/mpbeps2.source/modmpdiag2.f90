!-----------------------------------------------------------------------
!
      module mpdiag2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpdiag2.f
! get_funit returns an unconnected fortran unit number
! fnrecl find record length of direct access file
! dafopen2 opens new binary file for real 2d scalar data.
! dafopenv2 opens new binary file for real 2d vector data.
! dafopenc2 opens new binary file for complex 2d scalar data.
! dafopenvc2 opens new binary file for complex 2d vector data.
! pgwrite2 copies distributed scalar real data collected from segments
!          into a 2d global scalar array
! pgvwrite2 copies distributed vector real data collected from segments
!           into a 2d global vector array
! pgcwrite2 copies distributed scalar complex data collected from
!           segments into a 2d global scalar array
! pgvcwrite2 copies distributed vector complexdata collected from
!            segments into a 2d global vector array
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 1, 2017
!
!     use libmpdiag2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      function get_funit(start) result(funit)
! this function returns an unconnected fortran unit number,
! starting with unit = start.  returns -1 if none found
      integer, intent(in) :: start
      integer :: funit
! local data
      integer :: i
      logical :: connected
      funit = -1
! check connection status
      do i = start, 99
         inquire(unit=i,opened=connected)
         if (.not.connected) then
            funit = i
            exit
         endif
      enddo
      end function
!
!-----------------------------------------------------------------------
      function fnrecl(fname) result(it)
! find record length of direct access file
      character(len=*), intent(in) :: fname
      integer :: it, ios
      inquire(file=fname,recl=it,iostat=ios)
      if (ios /= 0) it = 0
      end function
!
!-----------------------------------------------------------------------
      subroutine dafopen2(f,nx,kyp,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 2d real scalar data array to be written in each record
! nx/kyp = number of data elements per record to be written in x/y 
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx, kyp
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1); lrec = lrec*nx*kyp
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenv2(f,nx,kyp,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 2d real vector data array to be written in each record
! nx/kyp = number of data elements per record to be written in x/y 
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx, kyp
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1,1); lrec = lrec*size(f,1)*nx*kyp
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenc2(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = 2d complex scalar data array to be written in each record
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      complex, dimension(:,:), intent(in) :: fc
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fc(1,1); lrec = lrec*size(fc)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenvc2(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = 2d complex vector data array to be written in each record
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      complex, dimension(:,:,:), intent(in) :: fc
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fc(1,1,1); lrec = lrec*size(fc)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pgwrite2(f,g,nx,ny,nvp)
! this subroutine copies into a 2d global scalar real array, data which
! has been written in blocks with spatial decomposition by the procedure
! PPWRITE2. data must have been written with a uniform partition
! f = 2d real scalar data blocked input array
! g = 2d real scalar data global output array
! nx/ny = system length in x/y direction
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) :: nx, ny, nvp
      real, dimension(nx,(ny-1)/nvp+1,nvp), intent(in) :: f
      real, dimension(nx,ny), intent(inout) :: g
! local data
      integer :: j, k, m, kyp, kyb, koff, kyps
! kyp = number of real grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! kyb = minimum number of processors in distributed array
      kyb = (ny - 1)/kyp + 1
      do m = 1, kyb
      koff = kyp*(m - 1)
      kyps = min(kyp,max(0,ny-koff))
      do k = 1, kyps
      do j = 1, nx
      g(j,k+koff) = f(j,k,m)
      enddo
      enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pgvwrite2(f,g,nx,ny,ndim,nvp)
! this subroutine copies into a 2d global vector real array, data which
! has been written in blocks with spatial decomposition by the procedure
! PPVWRITE2. data must have been written with a uniform partition
! f = 2d real vector data blocked input array
! g = 2d real vector data global output array
! nx/ny = system length in x/y direction
! ndim = first dimension of data array f and g
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) :: nx, ny, ndim, nvp
      real, dimension(ndim,nx,(ny-1)/nvp+1,nvp), intent(in) :: f
      real, dimension(ndim,nx,ny), intent(inout) :: g
! local data
      integer :: i, j, k, m, kyp, kyb, koff, kyps
! kyp = number of real grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! kyb = minimum number of processors in distributed array
      kyb = (ny - 1)/kyp + 1
      do m = 1, nvp
      koff = kyp*(m - 1)
      kyps = min(kyp,max(0,ny-koff))
      do k = 1, kyps
      do j = 1, nx
      do i = 1, ndim
      g(i,j,k+koff) = f(i,j,k,m)
      enddo
      enddo
      enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pgcwrite2(f,g,nx,ny,kxp,nvp)
! this subroutine copies and transposes into a 2d global scalar complex
! array, data which has been written in blocks with spatial
! decomposition by the procedure PPCWRITE2.
! data must have been written with a uniform partition
! f = 2d complex scalar data blocked input array
! g = 2d complex scalar data global output array
! nx/ny = system length in x/y direction
! kxp = number of complex grids in each field partition in x direction
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) :: nx, ny, kxp, nvp
      complex, dimension(ny,kxp,nvp), intent(in) :: f
      complex, dimension(nx,ny), intent(inout) :: g
! local data
      integer :: j, k, m, joff, kxpp
      do m = 1, nvp
      joff = kxp*(m - 1)
      kxpp = min(kxp,max(0,nx-joff))
      do j = 1, kxpp
      do k = 1, ny
      g(j+joff,k) = f(k,j,m)
      enddo
      enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pgcread2(f,g,nx,ny,kxp,nvp)
! this subroutine copies and transposes a 2d global scalar complex array
! into data which will written in blocks with spatial decomposition
! by the procedure PPCREAD2.
! data is written with a uniform partition
! f = 2d complex scalar data blocked output array
! g = 2d complex scalar data global input array
! nx/ny = system length in x/y direction
! kxp = number of complex grids in each field partition in x direction
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) :: nx, ny, kxp, nvp
      complex, dimension(ny,kxp,nvp), intent(inout) :: f
      complex, dimension(nx,ny), intent(in) :: g
! local data
      integer :: j, k, m, joff, kxpp
      do m = 1, nvp
      joff = kxp*(m - 1)
      kxpp = min(kxp,max(0,nx-joff))
      do j = 1, kxpp
      do k = 1, ny
      f(k,j,m) = g(j+joff,k)
      enddo
      enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pgvcwrite2(f,g,nx,ny,kxp,ndim,nvp)
! this subroutine copies and transposes into a 2d global vector complex
! array, data which has been written in blocks with spatial
! decomposition by the procedure PPVCWRITE2.
! data must have been written with a uniform partition
! f = 2d complex vector data blocked input array
! g = 2d complex vector data global output array
! nx/ny = system length in x/y direction
! kxp = number of complex grids in each field partition in x direction
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) :: nx, ny, kxp, ndim, nvp
      complex, dimension(ndim,ny,kxp,nvp), intent(in) :: f
      complex, dimension(ndim,nx,ny), intent(inout) :: g
! local data
      integer :: i, j, k, m, joff, kxpp
      do m = 1, nvp
      joff = kxp*(m - 1)
      kxpp = min(kxp,max(0,nx-joff))
      do j = 1, kxpp
      do k = 1, ny
      do i = 1, ndim
      g(i,j+joff,k) = f(i,k,j,m)
      enddo
      enddo
      enddo
      enddo
      end subroutine
!
      end module
