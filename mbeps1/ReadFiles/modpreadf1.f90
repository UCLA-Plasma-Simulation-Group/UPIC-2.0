!-----------------------------------------------------------------------
! Module for processing 1d data written by 1D OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      module pread1
!
      contains
!
!-----------------------------------------------------------------------
      subroutine pcreadf1(idrun)
! reads compressed complex periodic 1d scalar data
! idrun = run identifier for current run, -1 if unknown
      use in1, only: indx, tend,  dt
      use cmfield1
      use graf1
      use graf2
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 4
      integer :: nplot = 1, iudm = 19, ius = 11
      integer :: i, n, m, id, ii, it, nx, nrec, norm, ierr
      integer :: nxh, modesx
      integer :: iw = 200, nts = 0
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
      character(len=20), dimension(ns) :: dname = (/                    &
     &'POTENTIAL           ','LONGITUDINAL EFIELD ',                    &
     &'ELECTRON DENSITY    ','ION DENSITY         '/)
      character(len=10), dimension(2) :: cwk = (/'   W > 0  ',          &
     &                                           '   W < 0  ' /)
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
      fname = 'diag1.'//cdrun
      call ffopen1(iudm,fname)
!
! determine which scalar diagnostics are available
      call readsdiags1(iudm,nscalars)
!
! open graphics device
      ierr = open_graphs(nplot)
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
            call closeff1(iudm)
            call close_graphs()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected scalar diagnostic:
! nts, modesx, nrec, norm, fname
         call sdiagparams1(iudm,n,nts,modesx,nrec,norm,fname)
!
! nx = number of grid points in x direction
         nx = 2**indx; nxh = nx/2
!
! allocate complex scalar array
         if (.not.allocated(sfieldc)) allocate(sfieldc(modesx))
! allocate scalar array
         if (.not.allocated(sfield)) allocate(sfield(nx))
!
! open stream file for scalar field
         call fsopen1(ius,fname)
!
! nrec = number of complete records
         write (*,*) 'records found: nrec = ', nrec
!
! allocate and initialize data for frequency analysis
         if (.not.allocated(wm)) allocate(wm(iw))
         if (.not.allocated(pkw)) allocate(pkw(modesx,iw,2))
         if (.not.allocated(pks)) allocate(pks(4,modesx,iw))
         if (.not.allocated(wk)) allocate(wk(modesx,2))
         do m = 1, iw
            wm(m) = ((wmax-wmin)/real(iw))*real(m-1)
         enddo
         dt = dt*real(nts)
         dnx = 6.28318530717959/real(nx)
         akmin = 0.0; akmax = dnx*real(modesx - 1)
         pks = 0.0d0
!
! prepare fft tables for decompression
         if (.not.allocated(mixup)) allocate(mixup(nxh))
         if (.not.allocated(sct)) allocate(sct(nxh))
         call mfft1_init(mixup,sct,indx)
!
! read complex scalar data and display
         do ii = 1, nrec
            call freadc1(ius,sfieldc,modesx)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
! perform incremental frequency analysis
            call micspect1(sfieldc,wm,pkw,pks,time,0.0,nrec,iw,modesx,nx,1)
! decompress field data
            call mwrmodes1(sfield,sfieldc,nx,modesx)
! fourier transform to real space
            call mfft1r(sfield,1,mixup,sct,indx)
! display real space data
            call dscaler1(sfield,trim(dname(n)),it,999,0,nx,ierr)
            if (ierr==1) exit
         enddo
!
         call reset_graphs()
         call reset_nplot(2,ierr)
!
! find the frequency with the maximum power for each mode
!        wk(:,1) = wm(maxloc(pkw(:,:,1),dim=2))
!        wk(:,2) = wm(maxloc(pkw(:,:,2),dim=2))
! display frequencies as a function of mode number
!        call dmscaler1(wk,trim(dname(n)),nrec,999,1,modesx,cwk,ierr)
!
         pmin = minval(pkw,pkw.gt.0.0)
         where (pkw > 0.0)
            pkw = alog(pkw)
         else where
            pkw = alog(pmin)
         end where
!
! display positive frequencies as a function of mode number
         it = nts*nrec
         fname = dname(n)//cwk(1)
         call dscalerl2(pkw(:,:,1),trim(fname),akmin,akmax,wmin,wmax,   &
     &it,999,2,modesx,iw,ierr)
! display negative frequencies as a function of mode number
         fname = dname(n)//cwk(2)
         call dscalerl2(pkw(:,:,2),trim(fname),akmin,akmax,wmin,wmax,   &
     &it,999,2,modesx,iw,ierr)
!
         call reset_nplot(1,ierr)
         call closeff1(ius)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pvcreadf1(idrun)
! reads compressed complex periodic 1d vector data
! idrun = run identifier for current run, -1 if unknown
      use in1, only: indx, ndim, tend, dt
      use cmfield1
      use graf1
      use graf2
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 6
      integer :: nplot = 1, iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nx, mdim, nrec, norm, ierr
      integer :: nxh, modesx
      integer :: iw = 800, nts = 0
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
     &'ELEC CURRENT DENSITY','VECTOR POTENTIAL    ',                    &
     &'TRANSVERSE EFIELD   ','MAGNETIC FIELD      ',                    &
     &'RADIATIVE VPOTENTIAL','ION CURRENT DENSITY '/)
      character(len=10), dimension(2) :: cwk = (/'   W > 0  ',          &
     &                                           '   W < 0  ' /)
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
      fname = 'diag1.'//cdrun
      call ffopen1(iudm,fname)
!
! determine which vector diagnostics are available
      call readvdiags1(iudm,nscalars)
!
! open graphics device
      mdim = ndim - 1
      nplot = mdim
      ierr = open_graphs(nplot)
! set palette to color wheel
      call STPALIT(2)
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
            call closeff1(iudm)
            call close_graphs()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected vector diagnostic:
! nts, modesx, nrec, norm, fname
         call vdiagparams1(iudm,n,nts,modesx,nrec,norm,fname)
!
! nx = number of global grid points in x direction
         nx = 2**indx; nxh = nx/2
!
! allocate complex vector array
         if (.not.allocated(vfieldc))allocate(vfieldc(mdim,modesx))
! allocate vector array
         if (.not.allocated(vfield)) allocate(vfield(mdim,nx))
!
! open stream file for vector field
         call fsopen1(iuv,fname)
!
! nrec = number of complete records
         write (*,*) 'records found: nrec = ', nrec
!
! allocate and initialize data for frequency analysis
         if (.not.allocated(wm)) allocate(wm(iw))
         if (.not.allocated(vpkw)) allocate(vpkw(mdim,modesx,iw,2))
         if (.not.allocated(vpks)) allocate(vpks(mdim,4,modesx,iw))
         if (.not.allocated(vwk)) allocate(vwk(mdim,modesx,2))
         do m = 1, iw
            wm(m) = ((wmax-wmin)/real(iw))*real(m-1)
         enddo
         dt = dt*real(nts)
         dnx = 6.28318530717959/real(nx)
         akmin = 0.0; akmax = dnx*real(modesx - 1)
         vpks = 0.0d0
!
! prepare fft tables for decompression
         if (.not.allocated(mixup)) allocate(mixup(nxh),)
         if (.not.allocated(sct)) allocate(sct(nxh))
         call mfft1_init(mixup,sct,indx)
!
! read and transpose complex vector data and display
         do ii = 1, nrec
            call freadvc1(iuv,vfieldc,mdim,modesx)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
! perform incremental frequency analysis
            call mivcspect1(vfieldc,wm,vpkw,vpks,time,0.0,nrec,iw,modesx&
     &,nx,1)
! decompress field data
            call mwrvmodes1(vfield,vfieldc,nx,modesx)
! fourier transform to real space
            call mfft1rn(vfield,1,mixup,sct,indx)
! display real space data
            call dvector1(vfield,trim(dname(n)),it,999,0,1,nx,ierr)
            if (ierr==1) exit
         enddo
!
         call reset_graphs()
         call reset_nplot(1,ierr)
         call reset_nplot(nplot,ierr)
!
! find the frequency with the maximum power for each mode
!        vwk(1,:,1) = wm(maxloc(vpkw(1,:,:,1),dim=2))
!        vwk(2,:,1) = wm(maxloc(vpkw(2,:,:,1),dim=2))
!        vwk(1,:,2) = wm(maxloc(vpkw(1,:,:,2),dim=2))
!        vwk(2,:,2) = wm(maxloc(vpkw(2,:,:,2),dim=2))
! display frequencies as a function of mode number
!        call dmvector1(vwk,trim(dname(n)),nrec,999,2,2,modesx,cwk,ierr)
!
         pmin = minval(vpkw,vpkw.gt.0.0)
         where (vpkw > 0.0)
            vpkw = alog(vpkw)
         else where
            vpkw = alog(pmin)
         end where
!
! display positive frequencies as a function of mode number
         it = nts*nrec
         fname = dname(n)//'('//cwk(1)//')'
         call dvectorl2(vpkw(:,:,:,1),trim(fname),akmin,akmax,wmin,wmax,&
     &it,999,2,1,modesx,iw,ierr)
! display negative frequencies as a function of mode number
         fname = dname(n)//'('//cwk(2)//')'
         call dvectorl2(vpkw(:,:,:,2),trim(fname),akmin,akmax,wmin,wmax,&
     &it,999,2,1,modesx,iw,ierr)
!
         call reset_nplot(1,ierr)
         call reset_nplot(nplot,ierr)
         call closeff1(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pflreadf1(idrun)
! reads real periodic 1d fluid data
! idrun = run identifier for current run, -1 if unknown
      use in1, only: indx, ndim, dt
      use cmfield1
      use graf1
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 2
      integer :: nplot = 4, iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nx, npro, nprd, nrec, ierr
      integer :: nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:), allocatable :: fms, vfield, ufield
      real, dimension(:), allocatable :: sfield
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
      fname = 'diag1.'//cdrun
      call ffopen1(iudm,fname)
!
! determine which fluid diagnostics are available
      call readfldiags1(iudm,nscalars)
!
! open graphics device
      ierr = open_graphs(nplot)
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
            call closeff1(iudm)
            call close_graphs()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected fluid diagnostic:
! nts, npro, nprd, nrec, fname
         call fldiagparams1(iudm,n,nts,npro,nprd,nrec,fname)
!
! nx = number of global grid points in x direction
         nx = 2**indx
!
! allocate vector arrays
         if (.not.allocated(sfield)) allocate(sfield(nx))
         if (.not.allocated(fms)) allocate(fms(nprd,nx))
         if (ndim==3) then
            if (.not.allocated(vfield)) allocate(vfield(ndim,nx))
            if (npro > 2) then
               if (.not.allocated(ufield)) allocate(ufield(2*ndim,nx))
            endif
         endif
         dt = dt*real(nts)
!
! open stream file for vector field
         call fsopen1(iuv,fname)
!
! nrec = number of complete records
         write (*,*) 'records found: nrec = ', nrec
!
! read and display data
         do ii = 1, nrec
! read real vector field
            call freadv1(iuv,fms,nprd,nx)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
! electrostatic case
            if (ndim==1) then
! display density in real space
               if (npro > 0) then
                  sfield = fms(1,:)
                  call dscaler1(sfield,trim(sname(n))//trim(ename(1)),it&
     &,999,1,nx,ierr)
                  if (ierr==1) exit
               endif
! display velocity field in real space
               if (npro > 1) then
                  sfield = fms(2,:)
                  call dscaler1(sfield,trim(sname(n))//trim(ename(2)),it&
     &,999,0,nx,ierr)
                  if (ierr==1) exit
               endif
! display pressure tensor in real space
               if (npro > 2) then
                  sfield= fms(3,:)
                  call dscaler1(sfield,trim(sname(n))//trim(ename(3)),it&
     &,999,0,nx,ierr)
                  if (ierr==1) exit
               endif
! display electron heat flux in real space
               if (npro==4) then
                  sfield = fms(5,:)
                  call dscaler1(sfield,trim(sname(n))//trim(ename(5)),it&
     &,999,0,nx,ierr)
                  if (ierr==1) exit
               endif
! electromagnetic case
            else if (ndim==3) then
! display density in real space
               if (npro > 0) then
                  sfield = fms(1,:)
                  call dscaler1(sfield,trim(sname(n))//trim(ename(1)),it&
     &,999,1,nx,ierr)
                  if (ierr==1) exit
               endif
! display velocity field in real space
               if (npro > 1) then
                  vfield  = fms(2:4,:)
                  call dvector1(vfield,trim(sname(n))//trim(ename(2)),it&
     &,999,0,2,nx,ierr)
                  if (ierr==1) exit
               endif
! display pressure tensor in real space
               if (npro > 2) then
                  ufield= fms(5:10,:)
                  call dvector1(ufield,trim(sname(n))//trim(ename(3)),it&
     &,999,0,2,nx,ierr)
                  if (ierr==1) exit
               endif
! display electron heat flux in real space
               if (npro==4) then
                  vfield = fms(12:14,:)
                  call dvector1(vfield,trim(sname(n))//trim(ename(5)),it&
     &,999,0,2,nx,ierr)
                  if (ierr==1) exit
               endif
            endif
         enddo
!
         call reset_graphs()
         call reset_nplot(4,ierr)
         call closeff1(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
     subroutine pfvreadf1(idrun)
! reads real periodic 1d velocity data
! idrun = run identifier for current run, -1 if unknown
      use in1, only: ndim, tend, dt
      use cmfield1
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
      fname = 'diag1.'//cdrun
      call ffopen1(iudm,fname)
!
! determine which velocity diagnostics are available
      call readfvdiags1(iudm,nscalars)
!
! open graphics device
      nplot = 4
      ierr = open_graphs(nplot)
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
            call closeff1(iudm)
            call close_graphs()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected velocity diagnostic:
! nts, nmv, nfvd, nfed, nrec, fname
         call fvdiagparams1(iudm,n,nts,nmv,nfvd,nfed,nrec,fname)
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
         call fsopen1(iuv,fname)
!
! nrec = number of complete records
         write (*,*) 'records found: nrec = ', nrec
!
! read and display data
         do ii = 1, nrec
! read real velocity fields
            call freadfv1(iuv,fvm,fv,fe,wk,ndim,nmvf,nfvd,nfed)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
            if (nfvd > 0) then
! display cylindrical distribution
               if ((ndim==3).and.(nfvd < ndim)) then
                  call displayfvb1(fv,fvm,trim(sname(n)),it,nmv,2,ierr)
                  if (ierr==1) exit
! display cartesian distribution
               else
                  call displayfv1(fv,fvm,trim(sname(n)),it,nmv,2,ierr)
                  if (ierr==1) exit
               endif
            endif
! display energy distribution
            if (nfed > 0) then
               call displayfe1(fe,wk,trim(sname(n)),it,nmv,ierr)
               if (ierr==1) exit
            endif
         enddo
!
         call reset_graphs()
         call reset_nplot(4,ierr)
         call closeff1(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ptreadf1(idrun)
! reads real periodic 1d velocity data
! idrun = run identifier for current run, -1 if unknown
      use in1, only: ndim, tend, dt
      use cmfield1
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
      fname = 'diag1.'//cdrun
      call ffopen1(iudm,fname)
!
! determine which trajectory diagnostics are available
      call readtrdiags1(iudm,nscalars)
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
         call closeff1(iudm)
         return
      endif
!
      write (*,*) ' trajectory diagnostic selected'
!
! open graphics device
      ierr = open_graphs(nplot)
!
! return parameters for selected trajectory diagnostic:
! nts, ndt, nst, nmv, ndimp, nprobt, nrec, fname
      call trdiagparams1(iudm,n,nts,ndt,nst,nmv,ndimp,nprobt,nrec,fname)
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
      call fsopen1(iuv,fname)
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
            call freadtr1(iuv,partt,ndimp,nprobt)
            partd(ii,:,:) = partt
! read and display cartesian distribution
         else if (nst==3) then
            call freadfv1(iuv,fvmtp,fvtp,fetp,wk,ndim,nmvf,ndim,0)
            call displayfv1(fvtp,fvmtp,trim(sname(ndt)),it,nmv,2,ierr)
            if (ierr==1) exit
         endif
      enddo
!
! display time history of trajectories
      if ((nst==1).or.(nst==2)) then
! display vx
         ix = 2
         call displaytr1(partd,0.0,dt,nrec,ix,999,ierr)
      endif
!
      if (allocated(partt)) deallocate(partt)
      if (allocated(partd)) deallocate(partd)
      if (allocated(fvmtp)) deallocate(fvmtp)
      if (allocated(fvtp)) deallocate(fvtp)
      if (allocated(fetp)) deallocate(fetp)
!
      call reset_graphs()
      call closeff1(iudm)
      call closeff1(iuv)
      call close_graphs()
      write (*,*)
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine psreadf1(idrun)
! reads real periodic 1d phase space data
! idrun = run identifier for current run, -1 if unknown
      use in1, only: indx, ndim, tend, dt
      use cmfield1
      use graf2
      use graf1
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 2
      integer :: nplot = 1, iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nmv, nsxb, nrec, nx
      integer :: nmv21, nmvf, ierr
      integer :: nts = 0
      real :: time, xmax, ymin, ymax
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:), allocatable :: fvs
      real, dimension(:,:,:), allocatable :: pvs
      real, dimension(:,:), allocatable :: ps, psx
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
      fname = 'diag1.'//cdrun
      call ffopen1(iudm,fname)
!
! determine which phase space diagnostics are available
      call readpsdiags1(iudm,nscalars)
!
! open graphics device
      nplot = ndim
      ierr = open_graphs(nplot)
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
            call closeff1(iudm)
            call close_graphs()
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected phase space diagnostic:
! nts, nmv, nsxb, nrec, fname
         call psdiagparams1(iudm,n,nts,nmv,nsxb,nrec,fname)
!
! nx = number of global grid points in x direction
         nx = 2**indx
! nmv21 = number of velocity bins
         nmv21 = 2*nmv + 1; nmvf = nmv21 + 1
!
! allocate velocity distribution arrays
         if (.not.allocated(fvs)) allocate(fvs(nmvf,ndim,nsxb))
         if (.not.allocated(pvs)) allocate(pvs(nmvf,ndim,nsxb))
         if (.not.allocated(ps)) allocate(ps(nmvf,nsxb))
         if (.not.allocated(psx)) allocate(psx(nsxb,nmvf))
         dt = dt*real(nts)
!
! open stream file for phase space field
         call fsopen1(iuv,fname)
!
! nrec = number of complete records
         write (*,*) 'records found: nrec = ', nrec
!
! read and display data
         do ii = 1, nrec
! read real velocity distribution fields
            call freadps1(iuv,fvs,ndim,nmvf,nsxb)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
            pvs = fvs
            time = dt*real(ii - 1)
! electrostatic case
            if (ndim==1) then
! select x-vx
               ps = pvs(:,1,:nsxb)
               psx = transpose(ps)
               xmax = real(nx)
               ymax = fvs(nmv21+1,1,1); ymin = -ymax
! display real space data
               fname = trim(sname(n))//' X-VX'
               call dscalerl2(psx,fname,0.0,xmax,ymin,ymax,it,999,1,nsxb&
     &,nmv21,ierr)
               if (ierr==1) exit
! electromagnetic case
            else if (ndim==3) then
! select x-vx
               ps = pvs(:,1,:nsxb)
               psx = transpose(ps)
               xmax = real(nx)
               ymax = fvs(nmv21+1,1,1); ymin = -ymax
! display real space data
               fname = trim(sname(n))//' X-VX'
               call dscalerl2(psx,fname,0.0,xmax,ymin,ymax,it,999,1,nsxb&
     &,nmv21,ierr)
               if (ierr==1) exit
! select x-vy
               ps = pvs(:,2,:nsxb)
               psx = transpose(ps)
               ymax = fvs(nmv21+1,2,1); ymin = -ymax
! display real space data
               fname = trim(sname(n))//' X-VY'
               call dscalerl2(psx,fname,0.0,xmax,ymin,ymax,it,999,1,nsxb&
     &,nmv21,ierr)
               if (ierr==1) exit
! select x-vz
               ps = pvs(:,3,:nsxb)
               psx = transpose(ps)
               ymax = fvs(nmv21+1,2,1); ymin = -ymax
! display real space data
               fname = trim(sname(n))//' X-VZ'
               call dscalerl2(psx,fname,0.0,xmax,ymin,ymax,it,999,1,nsxb&
     &,nmv21,ierr)
               if (ierr==1) exit
            endif
         enddo
!
         call reset_graphs()
         call closeff1(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
      end module

