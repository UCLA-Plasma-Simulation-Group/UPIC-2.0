!-----------------------------------------------------------------------
! This program reads compressed complex periodic 2d scalar data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program pcreadf2
      use in2, only: idrun, indx, indy, nvp, tend, dt, ci,              &
     &pot2d, ntp, fpname, modesxp, modesyp, nprec,                      &
     &dene2d, ntde, fdename, modesxde, modesyde, nderec,                &
     &deni2d, ntdi, fdiname, modesxdi, modesydi, ndirec
      use cmfield2
      use graf2
!
      implicit none
      integer, parameter :: ns = 3
      integer :: iudm = 19, ius = 11
      integer :: i, j, k, n, m, it, nx, ny, kxp, kxb, mxvh, lrec, nrec
      integer :: ios, ierr
      integer :: nxh, nyh, modesx, modesy, modesy2
      integer :: nyv, nxyh, nxhy
      integer :: nplot = 1, iw = 100, nts = 0
      real :: time
      integer, dimension(ns) :: nscalars = 0
      complex, dimension(:,:), allocatable :: sfieldc
! data for decompression
      real, dimension(:,:), allocatable :: sfield
      integer, dimension(:), allocatable :: mixup
      complex, dimension(:), allocatable :: sct
! data for frequency analysis
      real, dimension(:), allocatable :: wm
      real, dimension(:,:,:,:), allocatable :: pkw
      double precision, dimension(:,:,:,:), allocatable :: pks
      real, dimension(:,:,:), allocatable :: wk
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
         nts = ntp
         modesx = modesxp; modesy = modesyp
         fname = fpname; nrec = nprec
      case (2)
         read (iudm,dene2d,iostat=ios)
         nts = ntde
         modesx = modesxde; modesy = modesyde
         fname = fdename; nrec = nderec
      case (3)
         read (iudm,deni2d,iostat=ios)
         nts = ntdi
         modesx = modesxdi; modesy = modesydi
         fname = fdiname; nrec = ndirec
      end select
!
! nx/ny = number of grid points in x/y direction
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kxb = minimum number of processors in distributed complex array
      kxb = (modesx - 1)/kxp + 1
! modesy2 = total number of modes in y
      modesy2 = 2*modesy - 1
! nyv = second dimension of scalar field array, >= ny
      nyv = ny + 1
!
! allocate complex scalar array
      allocate(sfieldc(modesx,modesy2))
! allocate scalar array
      allocate(sfield(nx,nyv))
!
! open direct access file for scalar field
      mxvh = min(modesx,kxp)*kxb
      inquire(iolength=lrec) sfieldc(1,1); lrec = lrec*mxvh*modesy2
      open(unit=ius,file=fname,form='unformatted',access='direct',      &
     &recl=lrec,status='old')
!
! nrec = number of complete records
      nrec = nrec/kxb
      write (*,*) 'records found: nrec = ', nrec
!
! allocate and initialize data for frequency analysis
      allocate(wm(iw),pkw(modesx,modesy2,iw,2))
      allocate(pks(4,modesx,modesy2,iw),wk(modesx,modesy2,2))
      do m = 1, iw
         wm(m) = (2.0/real(iw))*real(m-1)
      enddo
      dt = dt*real(nts)
      pks = 0.0d0
!
! open graphics device
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! prepare fft tables for decompression
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
      allocate(mixup(nxhy),sct(nxyh))
      call mfft2_init(mixup,sct,indx,indy)
!
! read and transpose complex scalar data and display
      do i = 1, nrec
         read (unit=ius,rec=i) ((sfieldc(j,k),k=1,modesy2),j=1,modesx)
         it = nts*(ii - 1)
! display fourier space data
!        call dcscaler2(sfieldc,trim(dname(n)),it,999,1,1,2,modesx,     &
!    &modesy2,ierr)
!        if (ierr==1) exit
         time = dt*real(i - 1)
! perform incremental frequency analysis
         call micspect2(sfieldc,wm,pkw,pks,time,0.0,nrec,iw,modesx,     &
     &modesy)
! decompress field data
         call mwrmodes2(sfield,sfieldc,nx,ny,modesx,modesy)
! fourier transform to real space
         call mfft2r(sfield,1,mixup,sct,indx,indy)
! display real space data
         call dscaler2(sfield,trim(dname(n)),it,999,0,1,nx,ny,ierr)
         if (ierr==1) exit
      enddo
!
! find the frequency with the maximum power for each mode, omega > 0
      wk(:,:,1) = reshape(wm(reshape(maxloc(pkw(:,:,:,1),dim=3),        &
     &(/modesx*modesy2/))),(/modesx,modesy2/))
! display positive frequencies as a function of mode number
      call dscaler2(wk(:,:,1),trim(dname(n)),1,999,1,1,modesx,modesy2,  &
     &ierr)
! find the frequency with the maximum power for each mode, omega < 0
      wk(:,:,2) = reshape(wm(reshape(maxloc(pkw(:,:,:,2),dim=3),        &
     &(/modesx*modesy2/))),(/modesx,modesy2/))
! display negative frequencies as a function of mode number
      call dscaler2(wk(:,:,2),trim(dname(n)),-1,999,1,1,modesx,modesy2, &
     &ierr)
!
      call closeff2(iudm)
      call closeff2(ius)
      call close_graphs2()
!
      end program
