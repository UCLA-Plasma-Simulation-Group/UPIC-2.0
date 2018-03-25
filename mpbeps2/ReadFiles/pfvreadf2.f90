!-----------------------------------------------------------------------
! This program reads real 2d velocity data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program pflreadf2
      use in2, only: idrun, nvp, ndim, tend, dt
      use cmfield2
      use graf2
      use graf1
!
      implicit none
      integer, parameter :: ns = 2
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, ii, it, nmv, nfvd, nfed, nrec, nmv21, nmvf
      integer :: ierr
      integer :: nplot = 1, nts = 0
      real :: time, wk
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:), allocatable :: fvm, fv, fe
      character(len=20), dimension(ns) :: dname = (/                    &
     &'electron velocity   ','ion velocity        '/)
      character(len=8), dimension(ns) :: sname = (/'ELECTRON','ION     '&
     &/)
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
! determine which velocity diagnostics are available
      call readfvdiags2(iudm,nscalars)
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
         write (*,*) 'no velocity-space diagnostic files found'
         stop
      endif
!
      write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected velocity diagnostic
      call fvdiagparams2(iudm,n,nts,nmv,nfvd,nfed,nrec,fname)
      nplot = 4
!
! nmv21 = number of velocity bins
      nmv21 = 2*nmv + 1; nmvf = nmv21 + 1
!
! allocate velocity distribution arrays
      allocate(fvm(ndim,3),fv(nmvf,nfvd),fe(nmvf,nfed))
      dt = dt*real(nts)
!
! open stream file for distribution field
      call fsopen2(iuv,fname)
!
! nrec = number of complete records
      write (*,*) 'records found: nrec = ', nrec
!
! open graphics device
      ierr = open_graphs2(nplot)
!
! read and display data
      do ii = 1, nrec
! read real velocity fields
         call freadfv2(iuv,fvm,fv,fe,wk,ndim,nmvf,nfvd,nfed)
         it = nts*(ii - 1)
         time = dt*real(ii - 1)
         if (nfvd > 0) then
! display cylindrical distribution
            if ((ndim==3).and.(nfvd < ndim)) then
               write (*,*) sname(n),' cylindrical:it,time=',it,time
               call displayfvb1(fv,fvm,trim(sname(n)),it,nmv,2,ierr)
               if (ierr==1) exit
! display cartesian distribution
            else
               write (*,*) sname(n),' cartesian:it,time=',it,time
               call displayfv1(fv,fvm,trim(sname(n)),it,nmv,2,ierr)
               if (ierr==1) exit
            endif
         endif
! display energy distribution
         if (nfed > 0) then
            write (*,*) sname(n),' energy:it,time=',it,time
            call displayfe1(fe,wk,trim(sname(n)),it,nmv,ierr)
            if (ierr==1) exit
         endif
      enddo
!
      call closeff2(iudm)
      call closeff2(iuv)
      call close_graphs2()
!
      end program
