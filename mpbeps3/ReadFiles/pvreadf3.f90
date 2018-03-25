!-----------------------------------------------------------------------
! This program reads real periodic 3d vector data
! written by 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program preadf3
      use in3, only: idrun, indx, indy, indz, nvpy, nvpz, ndim, tend, dt
      use cmfield3
!
      implicit none
      integer, parameter :: ns = 7
      integer :: iudm = 19, iuv = 12
      integer :: i, m, n, ii, it, nx, ny, nz, kyp, kzp, kyb, kzb
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
      write (*,*) 'enter idrun:'
      read (5,*) idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag3.'//cdrun
      call ffopen3(iudm,fname)
!
! determine which vector diagnostics are available
      call readvdiags3(iudm,nscalars)
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
         write (*,*) 'no vector diagnostic files found'
         stop
      endif
!
      write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected vector diagnostic
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
      allocate(vfield(ndim,nx,nyv,nzv))
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
!
      call closeff3(iudm)
      call closeff3(iuv)
!
      end program
