!-----------------------------------------------------------------------
! This program reads real periodic 2d scalar data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program preadf2
      use in2, only: idrun, indx, indy, nvp, tend, dt
      use cmfield2
      use graf2
!
      implicit none
      integer, parameter :: ns = 3
      integer :: iudm = 19, ius = 11
      integer :: i, n, m, ii, it, nx, ny, kyp, kyb, nyv, nrec, ierr
      integer :: modesx, modesy
      integer :: nplot = 1, nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:), allocatable :: sfield
      character(len=16), dimension(ns) :: dname = (/'POTENTIAL       ', &
     &'ELECTRON DENSITY','ION DENSITY     '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
      write (*,*) 'enter idrun:'
      read (5,*) idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which scalar diagnostics are available
      call readsdiags2(iudm,nscalars)
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
! return parameters for selected scalar diagnostic
      call sdiagparams2(iudm,n,nts,modesx,modesy,nrec,fname)
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
      dt = dt*real(nts)
!
! open stream file for scalar field
      call fsopen2(ius,fname)
!
! nrec = number of complete records
      nrec = nrec/kyb
      write (*,*) 'records found: nrec = ', nrec
!
! open graphics device
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! read and display data
      do ii = 1, nrec
! read real scalar field
         call fread2(ius,sfield,nx,nyv)
         it = nts*(ii - 1)
         time = dt*real(ii - 1)
! show time
         write (*,*) 'it,time=',it,time
! display real space data
         call dscaler2(sfield,trim(dname(n)),it,999,0,1,nx,ny,ierr)
         if (ierr==1) exit
      enddo
!
      call closeff2(iudm)
      call closeff2(ius)
      call close_graphs2()
!
      end program
