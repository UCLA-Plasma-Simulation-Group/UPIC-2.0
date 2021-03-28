!-----------------------------------------------------------------------
! This program reads compressed complex periodic 1d scalar data
! written for 1D OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program pcreadf1
      use in1, only: idrun, indx, dt
      use cmfield1
      use graf1
      use graf2
!
      implicit none
      integer, parameter :: ns = 3
      integer :: iudm = 19, ius = 11
      integer :: i, n, m, ii, nx, nrec, norm
      integer :: ierr
      integer :: nxh, modesx
      integer :: nplot = 1, iw = 200, nts = 0
      real :: wmin = 0.0, wmax = 2.0
      real :: time, dnx, akmin, akmax, pmin
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      complex, dimension(:), allocatable :: sfieldc
! data for decompression
      real, dimension(:), allocatable :: sfield
      integer, dimension(:), allocatable :: mixup
      complex, dimension(:), allocatable :: sct
! data for frequency analysis
      real, dimension(:), allocatable :: wm
      real, dimension(:,:,:), allocatable :: pkw
      double precision, dimension(:,:,:), allocatable :: pks
      real, dimension(:,:), allocatable :: wk
      character(len=16), dimension(ns) :: dname = (/'POTENTIAL       ', &
     &'ELECTRON DENSITY','ION DENSITY     '/)
      character(len=10), dimension(2) :: cwk = (/'   W > 0  ',          &
     &                                           '   W < 0  ' /)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
      write (*,*) 'enter idrun:'
      read (5,*) idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag1.'//cdrun
      call ffopen1(iudm,fname)
!
! determine which scalar diagnostics are available
      call readsdiags1(iudm,nscalars)
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
      call sdiagparams1(iudm,n,nts,modesx,nrec,norm,fname)
!
! nx = number of grid points in x direction
      nx = 2**indx; nxh = nx/2
!
! allocate complex scalar array
      allocate(sfieldc(modesx))
! allocate scalar array
      allocate(sfield(nx))
!
! open stream file for scalar field
      call fsopen1(ius,fname)
!
! nrec = number of complete records
      write (*,*) 'records found: nrec = ', nrec
!
! allocate and initialize data for frequency analysis
      allocate(wm(iw),pkw(modesx,iw,2))
      allocate(pks(4,modesx,iw),wk(modesx,2))
      do m = 1, iw
         wm(m) = ((wmax-wmin)/real(iw))*real(m-1)
      enddo
      dt = dt*real(nts)
      dnx = 6.28318530717959/real(nx)
      akmin = 0.0; akmax = dnx*real(modesx - 1)
      pks = 0.0d0
!
! open graphics device
      ierr = open_graphs(nplot)
! set palette to color wheel
      call set_palit(2)
!
! prepare fft tables for decompression
      allocate(mixup(nxh),sct(nxh))
      call mfft1_init(mixup,sct,indx)
!
! read complex scalar data and display
      do ii = 1, nrec
         call freadc1(ius,sfieldc,modesx)
         time = dt*real(ii - 1)
! perform incremental frequency analysis
         call micspect1(sfieldc,wm,pkw,pks,time,0.0,nrec,iw,modesx,nx,1)
! decompress field data
         call mwrmodes1(sfield,sfieldc,nx,modesx)
! fourier transform to real space
         call mfft1r(sfield,1,mixup,sct,indx)
! display real space data
         call dscaler1(sfield,trim(dname(n)),ii,999,0,nx,ierr)
         if (ierr==1) exit
      enddo
!
      call reset_graphs()
      call reset_nplot(2,ierr)
!
! find the frequency with the maximum power for each mode
!     wk(:,1) = wm(maxloc(pkw(:,:,1),dim=2))
!     wk(:,2) = wm(maxloc(pkw(:,:,2),dim=2))
! display frequencies as a function of mode number
!     call dmscaler1(wk,trim(dname(n)),nrec,999,1,modesx,cwk,ierr)
!
      pmin = minval(pkw,pkw.gt.0.0)
      where (pkw > 0.0)
         pkw = alog(pkw)
      else where
         pkw = alog(pmin)
      end where
!
! display positive frequencies as a function of mode number
      fname = dname(n)//cwk(1)
      call dscalerl2(pkw(:,:,1),trim(fname),akmin,akmax,wmin,wmax,nrec, &
     &999,2,modesx,iw,ierr)
! display negative frequencies as a function of mode number
      fname = dname(n)//cwk(2)
      call dscalerl2(pkw(:,:,2),trim(fname),akmin,akmax,wmin,wmax,nrec, &
     &999,2,modesx,iw,ierr)
!
      call closeff(iudm)
      call closeff(ius)
      call close_graphs()
!
      end program
