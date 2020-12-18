!-----------------------------------------------------------------------
! Module for processing 2d data written by 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      module pread2
!
      contains
!
!-----------------------------------------------------------------------
      subroutine pcreadf2(idrun)
! reads compressed complex periodic 2d scalar data
! idrun = run identifier for current run, -1 if unknown
      use in2, only: indx, indy, nvp, tend, dt
      use cmfield2
      use graf2
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 3
      integer :: nplot = 1, iudm = 19, ius = 11
      integer :: i, j, k, n, m, id, ii, it, nx, ny, kxp, kxb, nrec, ierr
      integer :: nxh, nyh, modesx, modesy, modesy2, nyv, nxyh, nxhy
      integer :: iw = 100, nts = 0
      real :: time
! nscalars = table of available diagnostics
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
      character(len=16), dimension(ns) :: dname = (/'POTENTIAL       ', &
     &'ELECTRON DENSITY','ION DENSITY     '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! create string from idrun
      id = idrun
      if (idrun < 0) then
         write (*,*) 'enter idrun:'
         read (5,*) id
      endif
      write (cdrun,'(i10)') id
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which complex scalar diagnostics are available
      call readsdiags2(iudm,nscalars)
!
! open graphics device
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! select diagnostic
      m = sum(nscalars)
      do
         if (m > 0) then
            n = -1
            do while (n < 0)
               do i = 1, ns
                  if (nscalars(i)==1) then
                     write (*,*) 'enter ', i, 'for ', trim(dname(i))
                  endif
               enddo
               write (*,*) 'enter ', 0, 'for EXIT'
               read (5,*) n
               if (n==0) exit
               if ((n >= 1).and.(n <= ns)) then
                  if (nscalars(n)==0) n = -1
               else
                  n = -1
               endif
               if (n > 0) exit
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            enddo
         else
            write (*,*) 'no complex scalar diagnostic files found'
            n = 0
         endif
! exit procedure
         if (n==0) then
            if (allocated(sfieldc)) deallocate(sfieldc)
            if (allocated(sfield)) deallocate(sfield)
            if (allocated(wm)) deallocate(wm)
            if (allocated(pkw)) deallocate(pkw)
            if (allocated(pks)) deallocate(pks)
            if (allocated(wk)) deallocate(wk)
            if (allocated(mixup)) deallocate(mixup)
            if (allocated(sct)) deallocate(sct)
            call closeff2(iudm)
            call close_graphs2()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected scalar diagnostic:
! nts, modesx, modesy, nrec, fname
         call sdiagparams2(iudm,n,nts,modesx,modesy,nrec,fname)
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
         if (.not.allocated(sfieldc)) allocate(sfieldc(modesx,modesy2))
! allocate scalar array
         if (.not.allocated(sfield)) allocate(sfield(nx,nyv))
!
! open stream file for scalar field
         call fsopen2(ius,fname)
!
! nrec = number of complete records
         nrec = nrec/kxb
         write (*,*) 'records found: nrec = ', nrec
!
! allocate and initialize data for frequency analysis
         if (.not.allocated(wm)) allocate(wm(iw))
         if (.not.allocated(pkw)) allocate(pkw(modesx,modesy2,iw,2))
         if (.not.allocated(pks)) allocate(pks(4,modesx,modesy2,iw))
         if (.not.allocated(wk)) allocate(wk(modesx,modesy2,2))
         do m = 1, iw
            wm(m) = (2.0/real(iw))*real(m-1)
         enddo
         dt = dt*real(nts)
         pks = 0.0d0
!
! prepare fft tables for decompression
         nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
         if (.not.allocated(mixup)) allocate(mixup(nxhy))
         if (.not.allocated(sct)) allocate(sct(nxyh))
         call mfft2_init(mixup,sct,indx,indy)
!
! read and transpose complex scalar data and display
         do ii = 1, nrec
            call freadc2(ius,sfieldc,modesx,modesy2)
            it = nts*(ii - 1)
! display fourier space data
!           call dcscaler2(sfieldc,trim(dname(n)),it,999,1,1,2,modesx,  &
!    &modesy2,ierr)
!           if (ierr==1) exit
            time = dt*real(ii - 1)
! perform incremental frequency analysis
            call micspect2(sfieldc,wm,pkw,pks,time,0.0,nrec,iw,modesx,  &
     &modesy)
! decompress field data
            call mwrmodes2(sfield,sfieldc,nx,ny,modesx,modesy)
! fourier transform to real space
            call mfft2r(sfield,1,mixup,sct,indx,indy)
! show time
            write (*,*) 'it,time=',it,time
! display real space data
            call dscaler2(sfield,trim(dname(n)),it,999,0,1,nx,ny,ierr)
            if (ierr==1) exit
         enddo
!
! find the frequency with the maximum power for each mode, omega > 0
         wk(:,:,1) = reshape(wm(reshape(maxloc(pkw(:,:,:,1),dim=3),     &
     &(/modesx*modesy2/))),(/modesx,modesy2/))
! display positive frequencies as a function of mode number
         call dscaler2(wk(:,:,1),trim(dname(n)),1,999,1,1,modesx,modesy2&
     &,ierr)
! find the frequency with the maximum power for each mode, omega < 0
         wk(:,:,2) = reshape(wm(reshape(maxloc(pkw(:,:,:,2),dim=3),     &
     &(/modesx*modesy2/))),(/modesx,modesy2/))
! display negative frequencies as a function of mode number
         call dscaler2(wk(:,:,2),trim(dname(n)),-1,999,1,1,modesx,modesy2&
     &,ierr)
         call closeff2(ius)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pcvreadf2(idrun)
! reads compressed complex periodic 2d vector data
! idrun = run identifier for current run, -1 if unknown
      use in2, only: indx, indy, nvp, ndim, tend, dt
      use cmfield2
      use graf2
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 7
      integer :: nplot = 1, iudm = 19, iuv = 12
      integer :: i, j, k, n, m, id, ii, it, nx, ny, kxp, kxb, nrec, ierr
      integer :: nxh, nyh, modesx, modesy, modesy2, nyv, nxyh, nxhy
      integer :: iw = 100, nts = 0
      real :: time
! nscalars = table of available diagnostics
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
     &'LONGITUDINAL EFIELD ','ELEC CURRENT DENSITY',                    &
     &'VECTOR POTENTIAL    ','TRANSVERSE EFIELD   ',                    &
     &'MAGNETIC FIELD      ','RADIATIVE VPOTENTIAL',                    &
     &'ION CURRENT DENSITY '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! create string from idrun
      id = idrun
      if (idrun < 0) then
         write (*,*) 'enter idrun:'
         read (5,*) id
      endif
      write (cdrun,'(i10)') id
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which complex vector diagnostics are available
      call readvdiags2(iudm,nscalars)
!
! open graphics device
      nplot = ndim
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! select diagnostic
      m = sum(nscalars)
      do
         if (m > 0) then
            n = -1
            do while (n < 0)
               do i = 1, ns
                  if (nscalars(i)==1) then
                     write (*,*) 'enter ', i, 'for ', trim(dname(i))
                  endif
               enddo
               write (*,*) 'enter ', 0, 'for EXIT'
               read (5,*) n
               if (n==0) exit
               if ((n >= 1).and.(n <= ns)) then
                  if (nscalars(n)==0) n = -1
               else
                  n = -1
               endif
               if (n > 0) exit
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            enddo
         else
            write (*,*) 'no complex vector diagnostic files found'
            n = 0
         endif
! exit procedure
         if (n==0) then
            if (allocated(vfieldc)) deallocate(vfieldc)
            if (allocated(vfield)) deallocate(vfield)
            if (allocated(wm)) deallocate(wm)
            if (allocated(vpkw)) deallocate(vpkw)
            if (allocated(vpks)) deallocate(vpks)
            if (allocated(vwk)) deallocate(vwk)
            if (allocated(mixup)) deallocate(mixup)
            if (allocated(sct)) deallocate(sct)
            call closeff2(iudm)
            call close_graphs2()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected vector diagnostic:
! nts, modesx, modesy, nrec, fname
         call vdiagparams2(iudm,n,nts,modesx,modesy,nrec,fname)
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
! allocate complex vector array
         if (.not.allocated(vfieldc)) then
            allocate(vfieldc(ndim,modesx,modesy2))
         endif
! allocate vector array
         if (.not.allocated(vfield)) allocate(vfield(ndim,nx,nyv))
!
! open stream file for vector field
         call fsopen2(iuv,fname)
!
! nrec = number of complete records
         nrec = nrec/kxb
         write (*,*) 'records found: nrec = ', nrec
!
! allocate and initialize data for frequency analysis
         if (.not.allocated(wm)) allocate(wm(iw))
         if (.not.allocated(vpkw)) then
            allocate(vpkw(ndim,modesx,modesy2,iw,2))
         endif
         if (.not.allocated(vpks)) then
            allocate(vpks(ndim,4,modesx,modesy2,iw))
         endif
         if (.not.allocated(vwk)) allocate(vwk(ndim,modesx,modesy2,2))
         do m = 1, iw
            wm(m) = (2.0/real(iw))*real(m-1)
         enddo
         dt = dt*real(nts)
         vpks = 0.0d0
!
! prepare fft tables for decompression
         nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
         if (.not.allocated(mixup)) allocate(mixup(nxhy))
         if (.not.allocated(sct)) allocate(sct(nxyh))
         call mfft2_init(mixup,sct,indx,indy)
!
! read and transpose complex vector data and display
         do ii = 1, nrec
            call freadvc2(iuv,vfieldc,ndim,modesx,modesy2)
            it = nts*(ii - 1)
! display fourier space data
!           call dcvector2(vfieldc,trim(dname(n)),it,999,1,1,2,mxvh,    &
!    &modesy2,ierr)
!           if (ierr==1) exit
            time = dt*real(ii - 1)
! perform incremental frequency analysis
            call mivcspect2(vfieldc,wm,vpkw,vpks,time,0.0,nrec,iw,modesx&
     &,modesy)
! decompress field data
            call mwrvmodes2(vfield,vfieldc,nx,ny,modesx,modesy)
! fourier transform to real space
            call mfft2rn(vfield,1,mixup,sct,indx,indy)
! show time
            write (*,*) 'it,time=',it,time
! display real space data
            call dvector2(vfield,trim(dname(n)),it,999,0,1,1,nx,ny,ierr)
            if (ierr==1) exit
         enddo
!
! find the frequency with the maximum power for each mode, omega > 0
          vwk(:,:,:,1) = reshape(wm(reshape(maxloc(vpkw(:,:,:,:,1),     &
     &dim=4),(/ndim*modesx*modesy2/))),(/ndim,modesx,modesy2/))
! display positive frequencies as a function of mode number
         call dvector2(vwk(:,:,:,1),trim(dname(n)),1,999,1,1,1,modesx,  &
     &modesy2,ierr)
! find the frequency with the maximum power for each mode, omega < 0
         vwk(:,:,:,2) = reshape(wm(reshape(maxloc(vpkw(:,:,:,:,2),      &
     &dim=4),(/ndim*modesx*modesy2/))),(/ndim,modesx,modesy2/))
! display negative frequencies as a function of mode number
         call dvector2(vwk(:,:,:,2),trim(dname(n)),-1,999,1,1,1,modesx, &
     &modesy2,ierr)
         call closeff2(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine preadf2(idrun)
! reads real periodic 2d scalar data
! idrun = run identifier for current run, -1 if unknown
      use in2, only: indx, indy, nvp, tend, dt
      use cmfield2
      use graf2
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 3
      integer :: nplot = 1, iudm = 19, ius = 11
      integer :: i, n, m, id, ii, it, nx, ny, kyp, kyb, nyv, nrec, ierr
      integer :: modesx, modesy
      integer :: nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:), allocatable :: sfield
      character(len=16), dimension(ns) :: dname = (/'POTENTIAL       ', &
     &'ELECTRON DENSITY','ION DENSITY     '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! create string from idrun
      id = idrun
      if (idrun < 0) then
         write (*,*) 'enter idrun:'
         read (5,*) id
      endif
      write (cdrun,'(i10)') id
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which scalar diagnostics are available
      call readsdiags2(iudm,nscalars)
!
! open graphics device
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! select diagnostic
      m = sum(nscalars)
      do
         if (m > 0) then
            n = -1
            do while (n < 0)
               do i = 1, ns
                  if (nscalars(i)==1) then
                     write (*,*) 'enter ', i, 'for ', trim(dname(i))
                  endif
               enddo
               write (*,*) 'enter ', 0, 'for EXIT'
               read (5,*) n
               if (n==0) exit
               if ((n >= 1).and.(n <= ns)) then
                  if (nscalars(n)==0) n = -1
               else
                  n = -1
               endif
               if (n > 0) exit
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            enddo
         else
            write (*,*) 'no scalar diagnostic files found'
            n = 0
         endif
! exit procedure
         if (n==0) then
            if (allocated(sfield)) deallocate(sfield)
            call closeff2(iudm)
            call close_graphs2()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected scalar diagnostic:
! nts, modesx, modesy, nrec, fname
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
         if (.not.allocated(sfield)) allocate(sfield(nx,nyv))
         dt = dt*real(nts)
!
! open stream file for scalar field
         call fsopen2(ius,fname)
!
! nrec = number of complete records
         nrec = nrec/kyb
         write (*,*) 'records found: nrec = ', nrec
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
         call closeff2(ius)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pvreadf2(idrun)
! reads real periodic 2d vector data
! idrun = run identifier for current run, -1 if unknown
      use in2, only: indx, indy, nvp, ndim, tend, dt
      use cmfield2
      use graf2
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 7
      integer :: nplot = 1, iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nx, ny, kyp, kyb, nyv, nrec, ierr
      integer :: modesx, modesy
      integer :: nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:), allocatable :: vfield
      character(len=20), dimension(ns) :: dname = (/                    &
     &'LONGITUDINAL EFIELD ','ELEC CURRENT DENSITY',                    &
     &'VECTOR POTENTIAL    ','TRANSVERSE EFIELD   ',                    &
     &'MAGNETIC FIELD      ','RADIATIVE VPOTENTIAL',                    &
     &'ION CURRENT DENSITY '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! create string from idrun
      id = idrun
      if (idrun < 0) then
         write (*,*) 'enter idrun:'
         read (5,*) id
      endif
      write (cdrun,'(i10)') id
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which vector diagnostics are available
      call readvdiags2(iudm,nscalars)
!
! open graphics device
      nplot = ndim
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! select diagnostic
      m = sum(nscalars)
      do
         if (m > 0) then
            n = -1
            do while (n < 0)
               do i = 1, ns
                  if (nscalars(i)==1) then
                     write (*,*) 'enter ', i, 'for ', trim(dname(i))
                  endif
               enddo
               write (*,*) 'enter ', 0, 'for EXIT'
               read (5,*) n
               if (n==0) exit
               if ((n >= 1).and.(n <= ns)) then
                  if (nscalars(n)==0) n = -1
               else
                  n = -1
               endif
               if (n > 0) exit
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            enddo
         else
            write (*,*) 'no vector diagnostic files found'
            n = 0
         endif
! exit procedure
         if (n==0) then
            if (allocated(vfield)) deallocate(vfield)
            call closeff2(iudm)
            call close_graphs2()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected vector diagnostic:
! nts, modesx, modesy, nrec, fname
         call vdiagparams2(iudm,n,nts,modesx,modesy,nrec,fname)
!
! nx/ny = number of global grid points in x/y direction
         nx = 2**indx; ny = 2**indy
! kyp = number of real grids in each field partition in y direction
         kyp = (ny - 1)/nvp + 1
! kyb = minimum number of processors in distributed array
         kyb = (ny - 1)/kyp + 1
! nyv = third dimension of vector field array, >= ny
         nyv = kyp*kyb
!
! allocate vector array
         if (.not.allocated(vfield)) allocate(vfield(ndim,nx,nyv))
         dt = dt*real(nts)
!
! open stream file for vector field
         call fsopen2(iuv,fname)
!
! nrec = number of complete records
         nrec = nrec/kyb
         write (*,*) 'records found: nrec = ', nrec
!
! read and display data
         do ii = 1, nrec
! read real vector field
            call freadv2(iuv,vfield,ndim,nx,nyv)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
! show time
            write (*,*) 'it,time=',it,time
! display real space data
            call dvector2(vfield,trim(dname(n)),it,999,0,1,1,nx,ny,ierr)
            if (ierr==1) exit
         enddo
         call closeff2(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pflreadf2(idrun)
! reads real periodic 2d fluid data
! idrun = run identifier for current run, -1 if unknown
      use in2, only: indx, indy, nvp, ndim, tend, dt
      use cmfield2
      use graf2
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 2
      integer :: nplot = 1, iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nx, ny, kyp, kyb, nyv, npro, nprd
      integer :: nrec, ierr
      integer :: nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:), allocatable :: fms, vfield, ufield
      real, dimension(:,:), allocatable :: sfield
      character(len=20), dimension(ns) :: dname = (/                    &
     &'ELECT FLUID MOMENTS ','ION FLUID MOMENTS   '/)
      character(len=8), dimension(ns) :: sname = (/'ELECTRON','ION     '&
     &/)
      character(len=16), dimension(5) :: ename = (/' DENSITY        ',  &
     &' VELOCITY FIELD ',' PRESSURE TENSOR',' ENERGY         ',         &
     &' HEAT FLUX      '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! create string from idrun
      id = idrun
      if (idrun < 0) then
         write (*,*) 'enter idrun:'
         read (5,*) id
      endif
      write (cdrun,'(i10)') id
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which fluid diagnostics are available
      call readfldiags2(iudm,nscalars)
!
! open graphics device
      nplot = 4
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! select diagnostic
      m = sum(nscalars)
      do
         if (m > 0) then
            n = -1
            do while (n < 0)
               do i = 1, ns
                  if (nscalars(i)==1) then
                     write (*,*) 'enter ', i, 'for ', trim(dname(i))
                  endif
               enddo
               write (*,*) 'enter ', 0, 'for EXIT'
               read (5,*) n
               if (n==0) exit
               if ((n >= 1).and.(n <= ns)) then
                  if (nscalars(n)==0) n = -1
               else
                  n = -1
               endif
               if (n > 0) exit
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            enddo
         else
            write (*,*) 'no fluid moments diagnostic files found'
            n = 0
         endif
! exit procedure
         if (n==0) then
            if (allocated(sfield)) deallocate(sfield)
            if (allocated(fms)) deallocate(fms)
            if (allocated(vfield)) deallocate(vfield)
            if (allocated(ufield)) deallocate(ufield)
            call closeff2(iudm)
            call close_graphs2()
            return
         endif

!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected fluid diagnostic:
! nts, npro, nprd, nrec, fname
         call fldiagparams2(iudm,n,nts,npro,nprd,nrec,fname)
!
! nx/ny = number of global grid points in x/y direction
         nx = 2**indx; ny = 2**indy
! kyp = number of real grids in each field partition in y direction
         kyp = (ny - 1)/nvp + 1
! kyb = minimum number of processors in distributed array
         kyb = (ny - 1)/kyp + 1
! nyv = third dimension of vector field array, >= ny
         nyv = kyp*kyb
!
! allocate vector array
         if (.not.allocated(sfield)) allocate(sfield(nx,nyv))
         if (.not.allocated(fms)) allocate(fms(nprd,nx,nyv))
         if (ndim==3) then
            if (.not.allocated(vfield)) allocate(vfield(ndim,nx,nyv))
            if (npro > 2) then
               if (.not.allocated(ufield)) then
                  allocate(ufield(2*ndim,nx,nyv))
               endif
            endif
         endif
         dt = dt*real(nts)
!
! open stream file for vector field
         call fsopen2(iuv,fname)
!
! nrec = number of complete records
         nrec = nrec/kyb
         write (*,*) 'records found: nrec = ', nrec
!
! read and display data
         do ii = 1, nrec
! read real vector field
            call freadv2(iuv,fms,nprd,nx,nyv)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
! show time
            write (*,*) sname(n),' it,time=',it,time
! electrostatic case
            if (ndim==2) then
! display density in real space
               if (npro > 0) then
                  sfield = fms(1,:,:)
                  call dscaler2(sfield,trim(sname(n))//trim(ename(1)),it&
     &,999,1,1,nx,ny,ierr)
                  if (ierr==1) exit
               endif
! display velocity field in real space
               if (npro > 1) then
                  sfield = fms(2,:,:)
                  call dscaler2(sfield,trim(sname(n))//trim(ename(2)),it&
     &,999,0,1,nx,ny,ierr)
                  if (ierr==1) exit
               endif
! display pressure tensor in real space
               if (npro > 2) then
                  sfield= fms(3,:,:)
                  call dscaler2(sfield,trim(sname(n))//trim(ename(3)),it&
     &,999,0,1,nx,ny,ierr)
                  if (ierr==1) exit
               endif
! display electron heat flux in real space
               if (npro==4) then
                  sfield = fms(5,:,:)
                  call dscaler2(sfield,trim(sname(n))//trim(ename(5)),it&
     &,999,0,1,nx,ny,ierr)
                  if (ierr==1) exit
               endif
! electromagnetic case
            else if (ndim==3) then
! display density in real space
               if (npro > 0) then
                  sfield = fms(1,:,:)
                  call dscaler2(sfield,trim(sname(n))//trim(ename(1)),it&
     &,999,1,1,nx,ny,ierr)
                  if (ierr==1) exit
               endif
! display velocity field in real space
               if (npro > 1) then
                  vfield  = fms(2:4,:,:)
                  call dvector2(vfield,trim(sname(n))//trim(ename(2)),it&
     &,999,0,1,1,nx,ny,ierr)
                  if (ierr==1) exit
               endif
! display pressure tensor in real space
               if (npro > 2) then
                  ufield= fms(5:10,:,:)
                  call dvector2(ufield,trim(sname(n))//trim(ename(3)),it&
     &,999,0,1,1,nx,ny,ierr)
                  if (ierr==1) exit
               endif
! display electron heat flux in real space
               if (npro==4) then
                  vfield = fms(12:14,:,:)
                  call dvector2(vfield,trim(sname(n))//trim(ename(5)),it&
     &,999,0,1,1,nx,ny,ierr)
                  if (ierr==1) exit
               endif
            endif
         enddo
         call closeff2(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
     subroutine pfvreadf2(idrun)
! reads real periodic 2d velocity data
! idrun = run identifier for current run, -1 if unknown
      use in2, only: nvp, ndim, tend, dt
      use cmfield2
      use graf2
      use graf1
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 2
      integer :: nplot = 1, iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nmv, nfvd, nfed, nrec, nmv21, nmvf
      integer :: ierr
      integer :: nts = 0
      real :: time, wk
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:), allocatable :: fvm, fv, fe
      character(len=20), dimension(ns) :: dname = (/                    &
     &'ELECTRON VELOCITY   ','ION VELOCITY        '/)
      character(len=8), dimension(ns) :: sname = (/'ELECTRON','ION     '&
     &/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! create string from idrun
      id = idrun
      if (idrun < 0) then
         write (*,*) 'enter idrun:'
         read (5,*) id
      endif
      write (cdrun,'(i10)') id
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which velocity diagnostics are available
      call readfvdiags2(iudm,nscalars)
!
! open graphics device
      nplot = 4
      ierr = open_graphs2(nplot)
!
! select diagnostic
      m = sum(nscalars)
      do
         if (m > 0) then
            n = -1
            do while (n < 0)
               do i = 1, ns
                  if (nscalars(i)==1) then
                     write (*,*) 'enter ', i, 'for ', trim(dname(i))
                  endif
               enddo
               write (*,*) 'enter ', 0, 'for EXIT'
               read (5,*) n
               if (n==0) exit
               if ((n >= 1).and.(n <= ns)) then
                  if (nscalars(n)==0) n = -1
               else
                  n = -1
               endif
               if (n > 0) exit
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            enddo
         else
            write (*,*) 'no velocity-space diagnostic files found'
            n = 0
         endif
! exit procedure
         if (n==0) then
            if (allocated(fvm)) deallocate(fvm)
            if (allocated(fv)) deallocate(fv)
            if (allocated(fe)) deallocate(fe)
            call closeff2(iudm)
            call close_graphs2()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected velocity diagnostic:
! nts, nmv, nfvd, nfed, nrec, fname
         call fvdiagparams2(iudm,n,nts,nmv,nfvd,nfed,nrec,fname)
!
! nmv21 = number of velocity bins
         nmv21 = 2*nmv + 1; nmvf = nmv21 + 1
!
! allocate velocity distribution arrays
         if (.not.allocated(fvm)) allocate(fvm(ndim,3))
         if (.not.allocated(fv)) allocate(fv(nmvf,nfvd))
         if (.not.allocated(fe)) allocate(fe(nmvf,nfed))
         dt = dt*real(nts)
!
! open stream file for distribution field
         call fsopen2(iuv,fname)
!
! nrec = number of complete records
         write (*,*) 'records found: nrec = ', nrec
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
         call closeff2(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ptreadf2(idrun)
! reads real periodic 2d velocity data
! idrun = run identifier for current run, -1 if unknown
      use in2, only: nvp, ndim, tend, dt
      use cmfield2
      use graf2
      use graf1
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 1
      integer :: nplot = 1, iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, ix, ndt, nst, nmv, ndimp, nprobt
      integer :: nrec, nmvf, ierr
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
! create string from idrun
      id = idrun
      if (idrun < 0) then
         write (*,*) 'enter idrun:'
         read (5,*) id
      endif
      write (cdrun,'(i10)') id
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which trajectory diagnostics are available
      call readtrdiags2(iudm,nscalars)
!
! open graphics device
      ierr = open_graphs2(nplot)
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
         n = 0
      endif
! exit procedure
      if (n==0) then
         call closeff2(iudm)
         return
      endif
!
      write (*,*) ' trajectory diagnostic selected'
!
! return parameters for selected trajectory diagnostic
      call trdiagparams2(iudm,n,nts,ndt,nst,nmv,ndimp,nprobt,nrec,fname)
      nplot = 1
!
      write (*,*) trim(sname(ndt)), ' trajectory available'
!
! nmvf = size of velocity distribution array
      nmvf = 2*nmv + 2
!
! allocate trajectory arrays
      if ((nst==1).or.(nst==2)) then
         if (.not.allocated(partt)) allocate(partt(ndimp,nprobt))
         if (.not.allocated(partd)) allocate(partd(nrec,ndimp,nprobt))
! allocate trajectory distribution arrays
      else if (nst==3) then
         if (.not.allocated(fvmtp)) allocate(fvmtp(ndim,3))
         if (.not.allocated(fvtp)) allocate(fvtp(nmvf,ndim))
         if (.not.allocated(fetp)) allocate(fetp(nmvf,0))
      endif
      dt = dt*real(nts)
!
! open stream file for trajectory field
      call fsopen2(iuv,fname)
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
            call freadtr2(iuv,partt,ndimp,nprobt)
            partd(ii,:,:) = partt
            write (*,*) sname(ndt),' trajectory:it,time=',it,time
! display cartesian distribution
         else if (nst==3) then
            call freadfv2(iuv,fvmtp,fvtp,fetp,wk,ndim,nmvf,ndim,0)
            write (*,*) sname(ndt),' distribution:it,time=',it,time
            call displayfv1(fvtp,fvmtp,trim(sname(ndt)),it,nmv,2,ierr)
            if (ierr==1) exit
         endif
      enddo
!
! display time history of trajectories
      if ((nst==1).or.(nst==2)) then
! display vx
         ix = 3
         call displaytr1(partd,0.0,dt,nrec,ix,999,ierr)
      endif
      if (allocated(partt)) deallocate(partt)
      if (allocated(partd)) deallocate(partd)
      if (allocated(fvmtp)) deallocate(fvmtp)
      if (allocated(fvtp)) deallocate(fvtp)
      if (allocated(fetp)) deallocate(fetp)
      call closeff2(iuv)
      call close_graphs2()
      write (*,*)
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine psreadf2(idrun)
! reads real periodic 2d phase space data
! idrun = run identifier for current run, -1 if unknown
      use in2, only: indx, indy, nvp, ndim, tend, dt
      use cmfield2
      use graf2
      use graf1
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 2
      integer :: nplot = 1, iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nmv, nsxb, nsyb, nrec, nx, ny
      integer :: nmv21, nmvf, ierr
      integer :: nts = 0
      real :: time, xmax, ymin, ymax
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:,:), allocatable :: fvs
      real, dimension(:,:,:), allocatable :: pvs
      real, dimension(:,:), allocatable :: ps, psx, psy
      character(len=20), dimension(ns) :: dname = (/                    &
     &'ELECTRON PHASE SPACE','ION PHASE SPACE     '/)
      character(len=8), dimension(ns) :: sname = (/'ELECTRON','ION     '&
     &/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! create string from idrun
      id = idrun
      if (idrun < 0) then
         write (*,*) 'enter idrun:'
         read (5,*) id
      endif
      write (cdrun,'(i10)') id
      cdrun = adjustl(cdrun)
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which phase space diagnostics are available
      call readpsdiags2(iudm,nscalars)
!
! open graphics device
      nplot = 4
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! select diagnostic
      m = sum(nscalars)
      do
         if (m > 0) then
            n = -1
            do while (n < 0)
               do i = 1, ns
                  if (nscalars(i)==1) then
                     write (*,*) 'enter ', i, 'for ', trim(dname(i))
                  endif
               enddo
               write (*,*) 'enter ', 0, 'for EXIT'
               read (5,*) n
               if (n==0) exit
               if ((n >= 1).and.(n <= ns)) then
                  if (nscalars(n)==0) n = -1
               else
               n = -1
               endif
               if (n > 0) exit
               write (*,*) 'invalid entry, try again or enter 0 to quit'
            enddo
         else
            write (*,*) 'no phase space diagnostic files found'
            n = 0
         endif
! exit procedure
         if (n==0) then
            if (allocated(fvs)) deallocate(fvs)
            if (allocated(pvs)) deallocate(pvs)
            if (allocated(ps)) deallocate(ps)
            if (allocated(psx)) deallocate(psx)
            if (allocated(psy)) deallocate(psy)
            call closeff2(iudm)
            call close_graphs2()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected phase space diagnostic:
! nts, nmv, nsxb, nsyb, nrec, fname
         call psdiagparams2(iudm,n,nts,nmv,nsxb,nsyb,nrec,fname)
!
! nx/ny = number of global grid points in x/y direction
         nx = 2**indx; ny = 2**indy
! nmv21 = number of velocity bins
         nmv21 = 2*nmv + 1; nmvf = nmv21 + 1
!
! allocate velocity distribution arrays
         if (.not.allocated(fvs)) allocate(fvs(nmvf,ndim,nsxb,nsyb))
         if (.not.allocated(pvs)) allocate(pvs(nmvf,ndim,max(nsxb,nsyb)))
         if (.not.allocated(ps)) allocate(ps(nmvf,max(nsxb,nsyb)))
         if (.not.allocated(psx)) allocate(psx(nsxb,nmvf))
         if (.not.allocated(psy)) allocate(psy(nsyb,nmvf))
         dt = dt*real(nts)
!
! open stream file for phase space field
         call fsopen2(iuv,fname)
!
! nrec = number of complete records
         write (*,*) 'records found: nrec = ', nrec
!
! read and display data
         do ii = 1, nrec
! read real velocity distribution fields
            call freadps2(iuv,fvs,ndim,nmvf,nsxb,nsyb)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
            write (*,*) sname(n),' it,time=',it,time
! sum over y
            pvs = sum(fvs,dim=4)
! select x-vx
            ps = pvs(:,1,:nsxb)
            psx = transpose(ps)
            xmax = real(nx)
            ymax = fvs(nmv21+1,1,1,1); ymin = -ymax
! display real space data
            fname = trim(sname(n))//' X-VX'
            call dscalerl2(psx,fname,0.0,xmax,ymin,ymax,it,999,1,nsxb,nmv21&
     &,ierr)
            if (ierr==1) exit
! select x-vy
            ps = pvs(:,2,:nsxb)
            psx = transpose(ps)
            ymax = fvs(nmv21+1,2,1,1); ymin = -ymax
! display real space data
            fname = trim(sname(n))//' X-VY'
            call dscalerl2(psx,fname,0.0,xmax,ymin,ymax,it,999,1,nsxb,nmv21&
     &,ierr)
            if (ierr==1) exit
! sum over x
            pvs = sum(fvs,dim=3)
! select y-vx
            ps = pvs(:,1,:nsyb)
            psy = transpose(ps)
            xmax = real(ny)
            ymax = fvs(nmv21+1,1,1,1); ymin = -ymax
! display real space data
            fname = trim(sname(n))//' Y-VX'
            call dscalerl2(psy,fname,0.0,xmax,ymin,ymax,it,999,1,nsxb,nmv21&
     &,ierr)
            if (ierr==1) exit
! select y-vy
            ps = pvs(:,2,:nsyb)
            psy = transpose(ps)
            ymax = fvs(nmv21+1,2,1,1); ymin = -ymax
! display real space data
            fname = trim(sname(n))//' Y-VY'
            call dscalerl2(psy,fname,0.0,xmax,ymin,ymax,it,999,1,nsxb,nmv21&
     &,ierr)
            if (ierr==1) exit
         enddo
         call closeff2(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
      end module
