!-----------------------------------------------------------------------
!
      module mpdiag3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpdiag3.f
! get_funit returns an unconnected fortran unit number
! fnrecl find record length of direct access file
! dafopen3 opens new binary file for real 1d scalar data.
! dafopenv3 opens new binary file for real 1d vector data.
! dafopenc3 opens new binary file for complex 1d scalar data.
! dafopenvc3 opens new binary file for complex 1d vector data.
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 1, 2017
!
!     use libmpdiag3_h
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
      subroutine dafopen3(f,nx,kyp,kzp,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 3d real scalar data array to be written in each record
! nx/kyp/kzp = number of data elements per record to be written in x/y/z
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx, kyp, kzp
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1,1); lrec = lrec*nx*kyp*kzp
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenv3(f,nx,kyp,kzp,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 3d real vector data array to be written in each record
! nx/kyp/kzp = number of data elements per record to be written in x/y/z
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx, kyp, kzp
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:,:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1,1,1)
      lrec = lrec*size(f,1)*nx*kyp*kzp
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenc3(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = 3d complex scalar data array to be written in each record
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
      subroutine dafopenvc3(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = 3d complex vector data array to be written in each record
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      complex, dimension(:,:,:,:), intent(in) :: fc
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fc(1,1,1,1); lrec = lrec*size(fc)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pgwrite3(f,g,nx,ny,nz,nvpy,nvpz)
! this subroutine copies into a 3d global real array, data which has
! been written in blocks with 2D spatial decomposition by the procedure
! PPWRITE32. data must have been written with a uniform partition
! f = 3d real scalar data blocked input array
! g = 3d real scalar data global output array
! nx/ny/nz = system length in x/y/z direction
! nvpy/nvpz = number of real or virtual processors in y/z
      implicit none
      integer, intent(in) :: nx, ny, nz, nvpy, nvpz
      real, dimension(nx,(ny-1)/nvpy+1,(nz-1)/nvpz+1,nvpy,nvpz),        &
     &intent(in) :: f
      real, dimension(nx,ny,nz), intent(inout) :: g
! local data
      integer :: j, k, l, my, mz, kyp, kzp, koff, loff, kyps, kzps
! kyp = number of complex grids in each field partition in y direction
! kzp = number of complex grids in each field partition in z direction
      kyp = (ny - 1)/nvpy + 1
      kzp = (nz - 1)/nvpz + 1
      do mz = 1, nvpz
      loff = kzp*(mz - 1)
      kzps = min(kzp,max(0,nz-loff))
      do my = 1, nvpy
      koff = kyp*(my - 1)
      kyps = min(kyp,max(0,ny-koff))
      do l = 1, kzps
      do k = 1, kyps
      do j = 1, nx
      g(j,k+koff,l+loff) = f(j,k,l,my,mz)
      enddo
      enddo
      enddo
      enddo
      enddo
      end subroutine
!
      end module
