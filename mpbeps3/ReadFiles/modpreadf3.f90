!-----------------------------------------------------------------------
! Module for processing 3d data written by 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      module pread3
!
      contains
!
!-----------------------------------------------------------------------
      subroutine preadf3(idrun)
! reads real periodic 3d scalar data
! idrun = run identifier for current run, -1 if unknown
      use in3, only: indx, indy, indz, nvpy, nvpz, tend, dt
      use cmfield3
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 3
      integer :: iudm = 19, ius = 11
      integer :: i, m, n, id, ii, it, nx, ny, nz, kyp, kzp, kyb, kzb
      integer :: nyv, nzv, nrec
      integer :: modesx, modesy, modesz
      integer :: nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:), allocatable :: sfield
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
      fname = 'diag3.'//cdrun
      call ffopen3(iudm,fname)
!
! determine which scalar diagnostics are available
      call readsdiags3(iudm,nscalars)
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
            call closeff3(iudm)
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected scalar diagnostic:
! nts, modesx, modesy, modesz, nrec, fname
         call sdiagparams3(iudm,n,nts,modesx,modesy,modesz,nrec,fname)
!
! nx/ny/nz = number of global grid points in x/y/z direction
         nx = 2**indx; ny = 2**indy; nz = 2**indz
! kyp/kzp = number of real grids in each field partition in y/z
         kyp = (ny - 1)/nvpy + 1; kzp = (nz - 1)/nvpz + 1
! kyb/kzb = minimum number of processors in distributed array in y/z
         kyb = (ny - 1)/kyp + 1; kzb = (nz - 1)/kzp + 1
! nyv = second dimension of scalar field array, >= ny
! nzv = third dimension of scalar field array, >= nz
         nyv = kyp*kyb; nzv = kzp*kzb
!
! allocate scalar array
         if (.not.allocated(sfield)) allocate(sfield(nx,nyv,nzv))
         dt = dt*real(nts)
!
! open stream file for scalar field
         call fsopen3(ius,fname)
!
! nrec = number of complete records
         nrec = nrec/(kyb*kzb)
         write (*,*) 'records found: nrec = ', nrec
!
! read and transpose scalar data
         do ii = 1, nrec
! read real scalar field
            call fread3(ius,sfield,nx,kyp,kyb,kzp,kzb)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
! show time
            write (*,*) 'it,time=',it,time
         enddo
         call closeff3(ius)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pvreadf3(idrun)
! reads real periodic 3d vector data
! idrun = run identifier for current run, -1 if unknown
      use in3, only: indx, indy, indz, nvpy, nvpz, ndim, tend, dt
      use cmfield3
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 7
      integer :: iudm = 19, iuv = 12
      integer :: i, m, n, id, ii, it, nx, ny, nz, kyp, kzp, kyb, kzb
      integer :: nyv, nzv, nrec
      integer :: modesx, modesy, modesz
      integer :: nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:,:), allocatable :: vfield
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
      fname = 'diag3.'//cdrun
      call ffopen3(iudm,fname)
!
! determine which vector diagnostics are available
      call readvdiags3(iudm,nscalars)
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
            call closeff3(iudm)
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected vector diagnostic:
! nts, modesx, modesy, modesz, nrec, fname
         call vdiagparams3(iudm,n,nts,modesx,modesy,modesz,nrec,fname)
!
! nx/ny/nz = number of global grid points in x/y/z direction
         nx = 2**indx; ny = 2**indy; nz = 2**indz
! kyp/kzp = number of real grids in each field partition in y/z
         kyp = (ny - 1)/nvpy + 1; kzp = (nz - 1)/nvpz + 1
! kyb/kzb = minimum number of processors in distributed array in y/z
         kyb = (ny - 1)/kyp + 1; kzb = (nz - 1)/kzp + 1
! nyv = second dimension of scalar field array, >= ny
! nzv = third dimension of scalar field array, >= nz
         nyv = kyp*kyb; nzv = kzp*kzb
!
! allocate vector array
         if (.not.allocated(vfield)) allocate(vfield(ndim,nx,nyv,nzv))
         dt = dt*real(nts)
!
! open stream file for vector field
         call fsopen3(iuv,fname)
!
! nrec = number of complete records
         nrec = nrec/(kyb*kzb)
         write (*,*) 'records found: nrec = ', nrec
!
! read and transpose vector data
         do ii = 1, nrec
! read real vector field
            call freadv3(iuv,vfield,ndim,nx,kyp,kyb,kzp,kzb)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
! show time
            write (*,*) 'it,time=',it,time
         enddo
         call closeff3(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pflreadf3(idrun)
! reads real periodic 3d fluid data
! idrun = run identifier for current run, -1 if unknown
      use in3, only: indx, indy, indz, nvpy, nvpz, ndim, tend, dt
      use cmfield3
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 2
      integer :: iudm = 19, iuv = 12
      integer :: i, m, n, id, ii, it, nx, ny, nz, kyp, kzp, kyb, kzb
      integer :: nyv, nzv, npro, nprd, nrec
      integer :: nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:,:), allocatable :: fms, vfield, ufield
      real, dimension(:,:,:), allocatable :: sfield
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
      fname = 'diag3.'//cdrun
      call ffopen3(iudm,fname)
!
! determine which fluid diagnostics are available:
      call readfldiags3(iudm,nscalars)
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
            call closeff3(iudm)
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected fluid diagnostic:
! nts, npro, nprd, nrec, fname
         call fldiagparams3(iudm,n,nts,npro,nprd,nrec,fname)
!
! nx/ny/nz = number of global grid points in x/y/z direction
         nx = 2**indx; ny = 2**indy; nz = 2**indz
! kyp/kzp = number of real grids in each field partition in y/z
         kyp = (ny - 1)/nvpy + 1; kzp = (nz - 1)/nvpz + 1
! kyb/kzb = minimum number of processors in distributed array in y/z
         kyb = (ny - 1)/kyp + 1; kzb = (nz - 1)/kzp + 1
! nyv = second dimension of scalar field array, >= ny
! nzv = third dimension of scalar field array, >= nz
         nyv = kyp*kyb; nzv = kzp*kzb
!
! allocate vector array
         if (.not.allocated(sfield)) allocate(sfield(nx,nyv,nzv))
         if (.not.allocated(fms)) allocate(fms(nprd,nx,nyv,nzv))
         if (ndim==3) then
            if (.not.allocated(vfield)) then
               allocate(vfield(ndim,nx,nyv,nzv))
            endif
            if (npro > 2) then
               if (.not.allocated(ufield)) then
                  allocate(ufield(2*ndim,nx,nyv,nzv))
               endif
            endif
         endif
         dt = dt*real(nts)
!
! open stream file for vector field
         call fsopen3(iuv,fname)
!
! nrec = number of complete records
         nrec = nrec/(kyb*kzb)
         write (*,*) 'records found: nrec = ', nrec
!
! read and transpose vector data
         do ii = 1, nrec
! read real vector field
            call freadv3(iuv,fms,nprd,nx,kyp,kyb,kzp,kzb)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
! show time
            write (*,*) sname(n),' it,time=',it,time
         enddo
         call closeff3(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pfvreadf3(idrun)
! reads real periodic 3d velocity data
! idrun = run identifier for current run, -1 if unknown
      use in3, only: nvpy, nvpz, ndim, tend, dt
      use cmfield3
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 2
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nmv, nfvd, nfed, nrec, nmv21, nmvf
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
      fname = 'diag3.'//cdrun
      call ffopen3(iudm,fname)
!
! determine which velocity diagnostics are available
      call readfvdiags3(iudm,nscalars)
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
            call closeff3(iudm)
            return
         endif
!
! return parameters for selected velocity diagnostic:
! nts, nmv, nfvd, nfed, nrec, fname
         call fvdiagparams3(iudm,n,nts,nmv,nfvd,nfed,nrec,fname)
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
         call closeff3(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ptreadf3(idrun)
! reads real periodic 3d velocity data
! idrun = run identifier for current run, -1 if unknown
      use in3, only: nvpy, nvpz, ndim, tend, dt
      use cmfield3
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 1
      integer :: iudm = 19, iuv = 12
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
         n = 0
      endif
! exit procedure
      if (n==0) then
         call closeff3(iudm)
         return
      endif
!
      write (*,*) ' trajectory diagnostic selected'
!
! return parameters for selected trajectory diagnostic:
! nts, ndt, nst, nmv, ndimp, nprobt, nrec, fname
      call trdiagparams3(iudm,n,nts,ndt,nst,nmv,ndimp,nprobt,nrec,fname)
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
      if (allocated(partt)) deallocate(partt)
      if (allocated(partd)) deallocate(partd)
      if (allocated(fvmtp)) deallocate(fvmtp)
      if (allocated(fvtp)) deallocate(fvtp)
      if (allocated(fetp)) deallocate(fetp)
      call closeff3(iuv)
      write (*,*)
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine psreadf3(idrun)
! reads real periodic 3d phase space data
! idrun = run identifier for current run, -1 if unknown
      use in3, only: indx, indy, indz, nvpy, nvpz, ndim, tend, dt
      use cmfield3
      implicit none
      integer, intent(in) :: idrun
! local data
      integer, parameter :: ns = 2
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, id, ii, it, nmv, nsxb, nsyb, nszb, nrec
      integer :: nx, ny, nz, nmv21, nmvf, ierr
      integer :: nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      integer, dimension(:,:,:), allocatable :: noff, nyzp
      real, dimension(:,:,:,:,:), allocatable :: fvs, gvs
      real, dimension(:,:,:), allocatable :: pvs
      real, dimension(:,:), allocatable :: ps, psx, psy, psz
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
      fname = 'diag3.'//cdrun
      call ffopen3(iudm,fname)
!
! determine which phase space diagnostics are available
      call readpsdiags3(iudm,nscalars)
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
            if (allocated(noff)) deallocate(noff)
            if (allocated(nyzp)) deallocate(nyzp)
            if (allocated(fvs)) deallocate(fvs)
            if (allocated(pvs)) deallocate(pvs)
            if (allocated(ps)) deallocate(ps)
            if (allocated(psx)) deallocate(psx)
            if (allocated(psy)) deallocate(psy)
            if (allocated(psz)) deallocate(psz)
            call closeff3(iudm)
            return
         endif
!
         write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected phase space diagnostic:
! nts, nmv, nsxb, nsyb, nszb, nrec, fname
         call psdiagparams3(iudm,n,nts,nmv,nsxb,nsyb,nszb,nrec,fname)
!
! nx/ny/nz = number of global grid points in x/y/z direction
         nx = 2**indx; ny = 2**indy; nz = 2**indz
! nmv21 = number of velocity bins
         nmv21 = 2*nmv + 1; nmvf = nmv21 + 1
!
! allocate velocity distribution arrays
         if (.not.allocated(noff)) allocate(noff(2,nvpy,nvpz))
         if (.not.allocated(nyzp)) allocate(nyzp(2,nvpy,nvpz))
         if (.not.allocated(fvs)) allocate(fvs(nmvf,ndim,nsxb,nsyb,nszb))
         if (.not.allocated(pvs)) allocate(pvs(nmvf,ndim,max(nsxb,nsyb,nszb)))
         if (.not.allocated(ps)) allocate(ps(nmvf,max(nsxb,nsyb,nszb)))
         if (.not.allocated(psx)) allocate(psx(nsxb,nmvf))
         if (.not.allocated(psy)) allocate(psy(nsyb,nmvf))
         if (.not.allocated(psz)) allocate(psz(nszb,nmvf))
         dt = dt*real(nts)
!
! open stream file for phase space field
         call fsopen3(iuv,fname)
!
! reads distributed non-uniform partition information: noff, nyzp, ierr
         call freadncomp3(iuv,noff,nyzp,ierr)
         if (ierr /= 0) return
!
! nrec = number of complete records
         write (*,*) 'records found: nrec = ', nrec
!
! read and display data
         do ii = 1, nrec
! read real velocity distribution fields
            call freadps3(iuv,fvs,noff,nyzp,ndim,nmvf,nsxb,nvpy,nvpz)
            it = nts*(ii - 1)
            time = dt*real(ii - 1)
            write (*,*) sname(n),' it,time=',it,time
! sum over y and z
            call summit2(fvs,pvs,4,5)
! select x-vx
            ps = pvs(:,1,:nsxb)
            psx = transpose(ps)
! select x-vy
            ps = pvs(:,2,:nsxb)
            psx = transpose(ps)
! select x-vz
            ps = pvs(:,3,:nsxb)
            psx = transpose(ps)
! sum over x and z
            call summit2(fvs,pvs,3,5)
! select y-vx
            ps = pvs(:,1,:nsyb)
            psy = transpose(ps)
! select y-vy
            ps = pvs(:,2,:nsyb)
            psy = transpose(ps)
! select y-vz
            ps = pvs(:,3,:nsyb)
            psy = transpose(ps)
! sum over x and y
            call summit2(fvs,pvs,3,4)
! select z-vx
            ps = pvs(:,1,:nszb)
            psz = transpose(ps)
! select z-vy
            ps = pvs(:,2,:nszb)
            psz = transpose(ps)
! select z-vz
            ps = pvs(:,3,:nszb)
            psz = transpose(ps)
         enddo
         call closeff3(iuv)
         write (*,*)
      enddo
!
      end subroutine
!
      end module
