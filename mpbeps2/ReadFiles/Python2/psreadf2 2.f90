!-----------------------------------------------------------------------
! This program reads real 2d phase space data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program psreadf2
      use in2, only: idrun, indx, indy, nvp, ndim, tend, dt
      use cmfield2
      use graf2
      use graf1
!
      implicit none
      integer, parameter :: ns = 2
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, ii, it, nmv, nsxb, nsyb, nrec, nx, ny
      integer :: nmv21, nmvf, ierr
      integer :: nplot = 1, nts = 0
      real :: time, xmax, ymin, ymax
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:,:), allocatable :: fvs
      real, dimension(:,:,:), allocatable :: pvs
      real, dimension(:,:), allocatable :: ps, psx, psy
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
      fname = 'diag2.'//cdrun
      call ffopen2(iudm,fname)
!
! determine which phase space diagnostics are available
      call readpsdiags2(iudm,nscalars)
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
      call psdiagparams2(iudm,n,nts,nmv,nsxb,nsyb,nrec,fname)
      nplot = 4
!
! nx/ny = number of global grid points in x/y direction
      nx = 2**indx; ny = 2**indy
! nmv21 = number of velocity bins
      nmv21 = 2*nmv + 1; nmvf = nmv21 + 1
!
! allocate velocity distribution arrays
      allocate(fvs(nmvf,ndim,nsxb,nsyb))
      allocate(pvs(nmvf,ndim,max(nsxb,nsyb)),ps(nmvf,max(nsxb,nsyb)))
      allocate(psx(nsxb,nmvf),psy(nsyb,nmvf))
      dt = dt*real(nts)
!
! open stream file for phase space field
      call fsopen2(iuv,fname)
!
! nrec = number of complete records
      write (*,*) 'records found: nrec = ', nrec
!
! open graphics device
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
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
!
      call closeff2(iudm)
      call closeff2(iuv)
      call close_graphs2()
!
      end program
