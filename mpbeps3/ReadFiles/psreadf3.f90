!-----------------------------------------------------------------------
! This program reads real 3d phase space data
! written for 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program psreadf3
      use in3, only: idrun, indx, indy, indz, nvpy, nvpz, ndim, tend, dt
      use cmfield3
!
      implicit none
      integer, parameter :: ns = 2
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, ii, it, nmv, nsxb, nsyb, nszb, nrec
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
     &'electron phase space','ion phase space     '/)
      character(len=8), dimension(ns) :: sname = (/'ELECTRON','ION     '&
     &/)
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
! determine which phase space diagnostics are available
      call readpsdiags3(iudm,nscalars)
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
         write (*,*) 'no phase space diagnostic files found'
         stop
      endif
!
      write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected phase space diagnostic
      call psdiagparams3(iudm,n,nts,nmv,nsxb,nsyb,nszb,nrec,fname)
!
! nx/ny/nz = number of global grid points in x/y/z direction
      nx = 2**indx; ny = 2**indy; nz = 2**indz
! nmv21 = number of velocity bins
      nmv21 = 2*nmv + 1; nmvf = nmv21 + 1
!
! allocate velocity distribution arrays
      allocate(noff(2,nvpy,nvpz),nyzp(2,nvpy,nvpz))
      allocate(fvs(nmvf,ndim,nsxb,nsyb,nszb))
      allocate(pvs(nmvf,ndim,max(nsxb,nsyb,nszb)))
      allocate(ps(nmvf,max(nsxb,nsyb,nszb)))
      allocate(psx(nsxb,nmvf),psy(nsyb,nmvf),psz(nszb,nmvf))
      dt = dt*real(nts)
!
! open stream file for phase space field
      call fsopen3(iuv,fname)
!
! reads distributed non-uniform partition information
      call freadncomp3(iuv,noff,nyzp,ierr)
      if (ierr /= 0) stop
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
!
      call closeff3(iudm)
      call closeff3(iuv)
!
      end program
