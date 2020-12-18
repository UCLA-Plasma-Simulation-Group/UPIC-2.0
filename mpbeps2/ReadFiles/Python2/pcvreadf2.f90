!-----------------------------------------------------------------------
! This program reads compressed complex periodic 2d vector data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program pcvreadf2
      use in2, only: idrun, indx, indy, nvp, ndim, tend, dt,            &
     &ci, omx, omy, omz,                                                &
     &el2d, ntel, felname, modesxel, modesyel, nelrec,                  &
     &vpot2d, nta, faname, modesxa, modesya, narec,                     &
     &et2d, ntet, fetname, modesxet, modesyet, netrec,                  &
     &b2d, ntb, fbname, modesxb, modesyb, nbrec,                        &
     &vpotr2d, ntar, farname, modesxar, modesyar, narrec,               &
     &vcuri2d, ntji, fjiname, modesxji, modesyji, njirec
      use cmfield2
      use graf2
!
      implicit none
      integer, parameter :: ns = 6
      integer :: iudm = 19, iuv = 12
      integer :: i, j, k, n, m, ii, it, nx, ny, kxp, kxb, mxvh, lrec
      integer :: nrec, ios, ierr
      integer :: nxh, nyh, modesx, modesy, modesy2
      integer :: nyv, nxyh, nxhy
      integer :: nplot = 1, iw = 100, nts = 0
      real :: time
      integer, dimension(ns) :: nscalars = 0
      complex, dimension(:,:,:), allocatable :: vfieldc
! data for decompression
      real, dimension(:,:,:), allocatable :: vfield
      integer, dimension(:), allocatable :: mixup
      complex, dimension(:), allocatable :: sct
! data for frequency analysis
      real, dimension(:), allocatable :: wm
      real, dimension(:,:,:,:,:), allocatable :: vpkw
      double precision, dimension(:,:,:,:,:), allocatable :: vpks
      real, dimension(:,:,:,:), allocatable :: vwk
      character(len=20), dimension(ns) :: dname = (/                    &
     &'longitudinal efield ','vector potential    ',                    &
     &'transverse efield   ','magnetic field      ',                    &
     &'radiative vpotential','ion current density '/)
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
! determine which vector diagnostics are available
      do n = 1, ns
      select case(n)
! load metadata for longitudinal efield data
      case (1)
         read (iudm,el2d,iostat=ios)
! load metadata for vector potential data
      case (2)
         read (iudm,vpot2d,iostat=ios)
! load metadata for transverse efield data
      case (3)
         read (iudm,et2d,iostat=ios)
! load metadata for magnetic field data
      case (4)
         read (iudm,b2d,iostat=ios)
! load metadata for radiative vector potential data
      case (5)
         read (iudm,vpotr2d,iostat=ios)
! load metadata for ion current density data
      case (6)
         read (iudm,vcuri2d,iostat=ios)
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
         read (iudm,el2d,iostat=ios)
         nts = ntel
         modesx = modesxel; modesy = modesyel
         fname = felname; nrec = nelrec
      case (2)
         read (iudm,vpot2d,iostat=ios)
         nts = nta
         modesx = modesxa; modesy = modesya
         fname = faname; nrec = narec
      case (3)
         read (iudm,et2d,iostat=ios)
         nts = ntet
         modesx = modesxet; modesy = modesyet
         fname = fetname; nrec = netrec
      case (4)
         read (iudm,b2d,iostat=ios)
         nts = ntb
         modesx = modesxb; modesy = modesyb
         fname = fbname; nrec = nbrec
      case (5)
         read (iudm,vpotr2d,iostat=ios)
         nts = ntar
         modesx = modesxar; modesy = modesyar
         fname = farname; nrec = narrec
      case (6)
         read (iudm,vcuri2d,iostat=ios)
         nts = ntji
         modesx = modesxji; modesy = modesyji
         fname = fjiname; nrec = njirec
      end select
      rewind iudm
      nplot = ndim
!
! nx/ny = number of global grid points in x/y direction
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
      allocate(vfieldc(ndim,modesx,modesy2))
! allocate scalar array
      allocate(vfield(ndim,nx,nyv))
!
! open direct access file for vector field
      mxvh = min(modesx,kxp)*kxb
      inquire(iolength=lrec) vfieldc(1,1,1)
      lrec = lrec*ndim*mxvh*modesy2
      open(unit=iuv,file=fname,form='unformatted',access='direct',      &
     &recl=lrec,status='old')
!
! nrec = number of complete records
      nrec = nrec/kxb
      write (*,*) 'records found: nrec = ', nrec
!
! allocate and initialize data for frequency analysis
      allocate(wm(iw),vpkw(ndim,modesx,modesy2,iw,2))
      allocate(vpks(ndim,4,modesx,modesy2,iw))
      allocate(vwk(ndim,modesx,modesy2,2))
      do m = 1, iw
         wm(m) = (2.0/real(iw))*real(m-1)
      enddo
      dt = dt*real(nts)
      vpks = 0.0d0
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
! read and transpose complex vector data and display
      do ii = 1, nrec
         read (unit=iuv,rec=ii) (((vfieldc(i,j,k),i=1,ndim),k=1,modesy2)&
     &,j=1,modesx)
         it = nts*(ii - 1)
! display fourier space data
!        call dcvector2(vfieldc,trim(dname(n)),it,999,1,1,2,mxvh,modesy2&
!    &,ierr)
!        if (ierr==1) exit
         time = dt*real(ii - 1)
! perform incremental frequency analysis
         call mivcspect2(vfieldc,wm,vpkw,vpks,time,0.0,nrec,iw,modesx,  &
     &modesy)
! decompress field data
         call mwrvmodes2(vfield,vfieldc,nx,ny,modesx,modesy)
! fourier transform to real space
         call mfft2rn(vfield,1,mixup,sct,indx,indy)
! display real space data
         call dvector2(vfield,trim(dname(n)),it,999,0,1,1,nx,ny,ierr)
         if (ierr==1) exit
      enddo
!
! find the frequency with the maximum power for each mode, omega > 0
      vwk(:,:,:,1) = reshape(wm(reshape(maxloc(vpkw(:,:,:,:,1),dim=4),  &
     &(/ndim*modesx*modesy2/))),(/ndim,modesx,modesy2/))
! display positive frequencies as a function of mode number
      call dvector2(vwk(:,:,:,1),trim(dname(n)),1,999,1,1,1,modesx,     &
     &modesy2,ierr)
! find the frequency with the maximum power for each mode, omega < 0
      vwk(:,:,:,2) = reshape(wm(reshape(maxloc(vpkw(:,:,:,:,2),dim=4),  &
     &(/ndim*modesx*modesy2/))),(/ndim,modesx,modesy2/))
! display negative frequencies as a function of mode number
      call dvector2(vwk(:,:,:,2),trim(dname(n)),-1,999,1,1,1,modesx,    &
     &modesy2,ierr)
!
      call closeff2(iudm)
      call closeff2(iuv)
      call close_graphs2()
