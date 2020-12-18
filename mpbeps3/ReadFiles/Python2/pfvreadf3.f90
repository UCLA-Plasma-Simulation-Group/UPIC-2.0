!-----------------------------------------------------------------------
! This program reads real periodic 3d velocity data
! written for 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program pfvreadf3
      use in3, only: idrun, nvpy, nvpz, ndim, tend, dt
      use cmfield3
!
      implicit none
      integer, parameter :: ns = 2
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, ii, it, nmv, nfvd, nfed, nrec, nmv21, nmvf
      integer :: nts = 0
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
      fname = 'diag3.'//cdrun
      call ffopen3(iudm,fname)
!
! determine which velocity diagnostics are available
      call readfvdiags3(iudm,nscalars)
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
! return parameters for selected velocity diagnostic
      call fvdiagparams3(iudm,n,nts,nmv,nfvd,nfed,nrec,fname)
!
! nmv21 = number of velocity bins
      nmv21 = 2*nmv + 1; nmvf = nmv21 + 1
!
! allocate velocity distribution arrays
      allocate(fvm(ndim,3),fv(nmvf,nfvd),fe(nmvf,nfed))
      dt = dt*real(nts)
!
! open stream file for velocity field
      call fsopen3(iuv,fname)
!
! nrec = number of complete records
      write (*,*) 'records found: nrec = ', nrec
!
! read and display data
      do ii = 1, nrec
! read real velocity fields
         call freadfv3(iuv,fvm,fv,fe,wk,ndim,nmvf,nfvd,nfed)
         it = nts*(ii - 1)
         time = dt*real(ii - 1)
! show time
         if (nfvd > 0) then
! display cylindrical distribution
            if (nfvd < ndim) then
               write (*,*) sname(n),' cylindrical:it,time=',it,time
! display cartesian distribution
            else
               write (*,*) sname(n),' cartesian:it,time=',it,time
            endif
         endif
! display energy distribution
         if (nfed > 0) then
            write (*,*) sname(n),' energy:it,time=',it,time
         endif
      enddo
!
      call closeff3(iudm)
      call closeff3(iuv)
!
      end program
