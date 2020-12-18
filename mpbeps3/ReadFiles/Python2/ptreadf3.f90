!-----------------------------------------------------------------------
! This program reads real 3d trajectory data
! written for 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ptreadf3
      use in3, only: idrun, nvpy, nvpz, ndim, tend, dt
      use cmfield3
!
      implicit none
      integer, parameter :: ns = 1
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, ii, it, ix, ndt, nst, nmv, ndimp, nprobt, nrec
      integer :: nmvf, ierr
      integer :: nts = 0
      real :: time, wk
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:), allocatable :: partt
      real, dimension(:,:,:), allocatable :: partd
      real, dimension(:,:), allocatable :: fvmtp, fvtp, fetp
      character(len=8), dimension(2) :: sname = (/'ELECTRON','ION     '/&
     &)
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
! determine which trajectory diagnostics are available
      call readtrdiags3(iudm,nscalars)
!
! select diagnostic
      m = sum(nscalars)
      if (m==1) then
         do i = 1, ns
            if (nscalars(i)==1) then
               n = i
               exit
            endif
         enddo
      else
         write (*,*) 'no trajectory diagnostic files found'
         stop
      endif
!
      write (*,*) ' trajectory diagnostic selected'
!
! return parameters for selected trajectory diagnostic
      call trdiagparams3(iudm,n,nts,ndt,nst,nmv,ndimp,nprobt,nrec,fname)
!
      write (*,*) trim(sname(ndt)), ' trajectory available'
!
! nmvf = size of velocity distribution array
      nmvf = 2*nmv + 2
!
! allocate trajectory arrays
      if ((nst==1).or.(nst==2)) then
         allocate(partt(ndimp,nprobt))
         allocate(partd(nrec,ndimp,nprobt))
! allocate trajectory distribution arrays
      else if (nst==3) then
         allocate(fvmtp(ndim,3),fvtp(nmvf,ndim),fetp(nmvf,0))
      endif
      dt = dt*real(nts)
!
! open stream file for trajectory field
      call fsopen3(iuv,fname)
!
! nrec = number of complete records
      write (*,*) 'records found: nrec = ', nrec
!
! read and display data
      do ii = 1, nrec
         it = nts*(ii - 1)
         time = dt*real(ii - 1)
! read real trajectory fields
         if ((nst==1).or.(nst==2)) then
            call freadtr3(iuv,partt,ndimp,nprobt)
            partd(ii,:,:) = partt
            write (*,*) sname(ndt),' trajectory:it,time=',it,time
! display cartesian distribution
         else if (nst==3) then
            call freadfv3(iuv,fvmtp,fvtp,fetp,wk,ndim,nmvf,ndim,0)
            write (*,*) sname(ndt),' distribution:it,time=',it,time
         endif
      enddo
!
      call closeff3(iudm)
      call closeff3(iuv)
!
      end program
