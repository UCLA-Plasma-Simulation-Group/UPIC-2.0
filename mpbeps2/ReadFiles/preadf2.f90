!-----------------------------------------------------------------------
! This program reads real periodic 2d scalar data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program preadf2
      use in2, only: idrun, indx, indy, nvp, tend, dt, ci,              &
     &pot2d, ntp, fpname, modesxp, modesyp, nprec,                      &
     &dene2d, ntde, fdename, modesxde, modesyde, nderec,                &
     &deni2d, ntdi, fdiname, modesxdi, modesydi, ndirec
      use graf2
!
      implicit none
      integer, parameter :: ns = 3
      integer :: iudm = 19, ius = 11
      integer :: i, n, m, nx, ny, kyp, kyb, nyv, lrec, nrec, ios, ierr
      integer :: nplot = 1
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:), allocatable :: sfield
      character(len=16), dimension(ns) :: dname = (/'potential       ', &
     &'electron density','ion density     '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
      write (*,*) 'enter idrun:'
      read (5,*) idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      open(unit=iudm,file=fname,form='formatted',status='old')
!
! determine which scalar diagnostics are available
      do n = 1, ns
      select case(n)
! load metadata for potential data
      case (1)
         read (iudm,pot2d,iostat=ios)
! load metadata for electron density data
      case (2)
         read (iudm,dene2d,iostat=ios)
! load metadata for ion density data
      case (3)
         read (iudm,deni2d,iostat=ios)
      end select
      if (ios==0) nscalars(n) = 1
      rewind iudm
      enddo
!
! select diagnostic
      m = sum(nscalars)
      if (m > 1) then
         n = -1
         do while (n < 0)
            do i = 1, ns
               if (nscalars(i)==1) then
                  write (*,*) 'enter ', i, 'for ', trim(dname(i))
               endif
            enddo
            read (5,*) n
            if (n==0) stop
            if ((n >= 1).and.(n <= ns)) then
               if (nscalars(n)==0) n = -1
            else
               n = -1
            endif
            if (n < 0) then
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            endif
         enddo
      else if (m==1) then
         do i = 1, ns
            if (nscalars(i)==1) then
               n = i
               exit
            endif
         enddo
      else
         write (*,*) 'no scalar diagnostic files found'
         stop
      endif
!
      write (*,*) trim(dname(n)), ' diagnostic selected'
!
      select case(n)
      case (1)
         read (iudm,pot2d,iostat=ios)
         fname = fpname; nrec = nprec
      case (2)
         read (iudm,dene2d,iostat=ios)
         fname = fdename; nrec = nderec
      case (3)
         read (iudm,deni2d,iostat=ios)
         fname = fdiname; nrec = ndirec
      end select
!
! nx/ny = number of global grid points in x/y direction
      nx = 2**indx; ny = 2**indy
! kyp = number of real grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! kyb = minimum number of processors in distributed array
      kyb = (ny - 1)/kyp + 1
! nyv = second dimension of scalar field array, >= ny
      nyv = kyp*kyb
!
! allocate scalar array
      allocate(sfield(nx,nyv))
! open direct access file for scalar field
      inquire(iolength=lrec) sfield(1,1); lrec = lrec*nx*nyv
      open(unit=ius,file=fname,form='unformatted',access='direct',      &
     &recl=lrec,status='old')
!
! nrec = number of complete records
      nrec = nrec/kyb
      write (*,*) 'records found: nrec = ', nrec
!
! open graphics device
      ierr = open_graphs(nplot)
! set palette to color wheel
      call STPALIT(2)
!
! read and display data
      do i = 1, nrec
         read (unit=ius,rec=i) sfield
         call dscaler2(sfield,trim(dname(n)),i,999,0,1,nx,ny,ierr)
         if (ierr==1) exit
      enddo
!
      call close_graphs
!
      end program
!
! unneeded function in input2mod.f90
      subroutine PPBDCAST(f,nxp)
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f
      end subroutine
!
