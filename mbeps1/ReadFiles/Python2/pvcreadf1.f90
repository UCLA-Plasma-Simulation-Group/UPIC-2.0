!-----------------------------------------------------------------------
! This program reads compressed complex periodic 1d vector data
! written for 1D OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program pvcreadf1
      use in1, only: idrun, indx, ndim, dt
      use cmfield1
      use graf1
      use graf2
!
      implicit none
      integer, parameter :: ns = 7
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, ii, nx, mdim, nrec, norm
      integer :: ierr
      integer :: nxh, modesx
      integer :: nplot = 1, iw = 800, nts = 0
      real :: wmin = 0.0, wmax = 8.0
      real :: time, dnx, akmin, akmax, pmin
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      complex, dimension(:,:), allocatable :: vfieldc
! data for decompression
      real, dimension(:,:), allocatable :: vfield
      integer, dimension(:), allocatable :: mixup
      complex, dimension(:), allocatable :: sct
! data for frequency analysis
      real, dimension(:), allocatable :: wm
      real, dimension(:,:,:,:), allocatable :: vpkw
      double precision, dimension(:,:,:,:), allocatable :: vpks
      real, dimension(:,:,:), allocatable :: vwk
      character(len=20), dimension(ns) :: dname = (/                    &
     &'LONGITUDINAL EFIELD ','ELEC CURRENT DENSITY',                    &
     &'VECTOR POTENTIAL    ','TRANSVERSE EFIELD   ',                    &
     &'MAGNETIC FIELD      ','RADIATIVE VPOTENTIAL',                    &
     &'ION CURRENT DENSITY '/)
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
! determine which vector diagnostics are available
      call readvdiags1(iudm,nscalars)
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
! return parameters for selected vector diagnostic
      call vdiagparams1(iudm,n,nts,modesx,nrec,norm,fname)
      mdim = ndim - 1
      nplot = mdim
!
! nx = number of global grid points in x direction
      nx = 2**indx; nxh = nx/2
!
! allocate complex vector array
      allocate(vfieldc(mdim,modesx))
! allocate vector array
      allocate(vfield(mdim,nx))
!
! open stream file for vector field
      call fsopen1(iuv,fname)
!
! nrec = number of complete records
      write (*,*) 'records found: nrec = ', nrec
!
! allocate and initialize data for frequency analysis
      allocate(wm(iw),vpkw(mdim,modesx,iw,2))
      allocate(vpks(mdim,4,modesx,iw))
      allocate(vwk(mdim,modesx,2))
      do m = 1, iw
         wm(m) = ((wmax-wmin)/real(iw))*real(m-1)
      enddo
      dt = dt*real(nts)
      dnx = 6.28318530717959/real(nx)
      akmin = 0.0; akmax = dnx*real(modesx - 1)
      vpks = 0.0d0
!
! open graphics device
      ierr = open_graphs(nplot)
! set palette to color wheel
      call STPALIT(2)
!
! prepare fft tables for decompression
      allocate(mixup(nxh),sct(nxh))
      call mfft1_init(mixup,sct,indx)
!
! read and transpose complex vector data and display
      do ii = 1, nrec
         call freadvc1(iuv,vfieldc,mdim,modesx)
         time = dt*real(ii - 1)
! perform incremental frequency analysis
         call mivcspect1(vfieldc,wm,vpkw,vpks,time,0.0,nrec,iw,modesx,nx&
     &,1)
! decompress field data
         call mwrvmodes1(vfield,vfieldc,nx,modesx)
! fourier transform to real space
         call mfft1rn(vfield,1,mixup,sct,indx)
! display real space data
         call dvector1(vfield,trim(dname(n)),ii,999,0,1,nx,ierr)
         if (ierr==1) exit
      enddo
!
      call reset_graphs()
      call reset_nplot(1,ierr)
      call reset_nplot(4,ierr)
!
! find the frequency with the maximum power for each mode
!     vwk(1,:,1) = wm(maxloc(vpkw(1,:,:,1),dim=2))
!     vwk(2,:,1) = wm(maxloc(vpkw(2,:,:,1),dim=2))
!     vwk(1,:,2) = wm(maxloc(vpkw(1,:,:,2),dim=2))
!     vwk(2,:,2) = wm(maxloc(vpkw(2,:,:,2),dim=2))
! display frequencies as a function of mode number
!     call dmvector1(vwk,trim(dname(n)),nrec,999,2,2,modesx,cwk,ierr)
!
      pmin = minval(vpkw,vpkw.gt.0.0)
      where (vpkw > 0.0)
         vpkw = alog(vpkw)
      else where
         vpkw = alog(pmin)
      end where
!
! display positive frequencies as a function of mode number
      fname = dname(n)//'('//cwk(1)//')'
! display positive frequencies as a function of mode number
      call dvectorl2(vpkw(:,:,:,1),trim(fname),akmin,akmax,wmin,wmax,   &
     &nrec,999,2,1,modesx,iw,ierr)
! display negative frequencies as a function of mode number
      fname = dname(n)//'('//cwk(2)//')'
      call dvectorl2(vpkw(:,:,:,2),trim(fname),akmin,akmax,wmin,wmax,   &
     &nrec,999,2,1,modesx,iw,ierr)
!
      call closeff(iudm)
      call closeff(iuv)
      call close_graphs()
!
      end program

