!-----------------------------------------------------------------------
! This program reads 1d data written by 1D OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program diag1
      use cmfield1
      use graf1
      use pread1
      implicit none
      integer, parameter :: ms = 6, ns = 6
      integer :: idrun = 0, iudm = 19
      integer :: i, m, n
! mscalars = table of available diagnostic types
      integer, dimension(ms) :: mscalars = 0
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      character(len=22), dimension(ms) ::                               &
     &dname = (/'COMPLEX SCALAR FIELDS ','COMPLEX VECTOR FIELDS ',      &
     &'REAL FLUID DATA       ','VELOCITY DATA         ',                &
     &'TRAJECTORY DATA       ','PHASE SPACE DATA      '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! create string from idrun
      write (*,*) 'enter idrun:'
      read (5,*) idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag1.'//cdrun
      call ffopen1(iudm,fname)
!
! determine which scalar diagnostics are available
      call readsdiags1(iudm,nscalars)
      mscalars(1) = sum(nscalars)
! determine which vector diagnostics are available
      nscalars = 0
      call readvdiags1(iudm,nscalars)
      mscalars(2) = sum(nscalars)
! determine which fluid diagnostics are available
      nscalars = 0
      call readfldiags1(iudm,nscalars)
      mscalars(3) = sum(nscalars)
! determine which velocity diagnostics are available
      nscalars = 0
      call readfvdiags1(iudm,nscalars)
      mscalars(4) = sum(nscalars)
! determine which trajectory diagnostics are available
      nscalars = 0
      call readtrdiags1(iudm,nscalars)
      mscalars(5) = sum(nscalars)
! determine which phase space diagnostics are available
      nscalars = 0
      call readpsdiags1(iudm,nscalars)
      mscalars(6) = sum(nscalars)
!
! select diagnostic
      m = sum(mscalars)
      do
         if (m > 0) then
            n = -1
            do while (n < 0)
               do i = 1, ms
                  if (mscalars(i) > 0) then
                     write (*,*) 'enter ', i, 'for ', trim(dname(i))
                  endif
               enddo
               write (*,*) 'enter ', 0, 'for EXIT'
               read (5,*) n
               if (n==0) then
                  call closeff1(iudm)
                  stop
               endif
               if ((n >= 1).and.(n <= ns)) then
                  if (mscalars(n)==0) n = -1
               else
                  n = -1
               endif
               if (n > 0) exit
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            enddo
         else
            write (*,*) 'no diagnostic data found'
            stop
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
         if (n==1) then
            call pcreadf1(idrun)
         else if (n==2) then
            call pvcreadf1(idrun)
         else if (n==3) then
            call pflreadf1(idrun)
         else if (n==4) then
            call pfvreadf1(idrun)
         else if (n==5) then
            call ptreadf1(idrun)   
         else if (n==6) then
            call psreadf1(idrun)
         endif
      enddo
!
      end program
