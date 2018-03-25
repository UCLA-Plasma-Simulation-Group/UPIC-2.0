!-----------------------------------------------------------------------
! This program reads real periodic 2d fluid data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program pflreadf2
      use in2, only: idrun, indx, indy, nvp, ndim, tend, dt
      use cmfield2
      use graf2
!
      implicit none
      integer, parameter :: ns = 2
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, ii, it, nx, ny, kyp, kyb, nyv, npro, nprd
      integer :: nrec, ierr
      integer :: nplot = 1, nts = 0
      real :: time
! nscalars = table of available diagnostics
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:), allocatable :: fms, vfield, ufield
      real, dimension(:,:), allocatable :: sfield
      character(len=20), dimension(ns) :: dname = (/                    &
     &'elect fluid moments ','ion fluid moments   '/)
      character(len=8), dimension(ns) :: sname = (/'ELECTRON','ION     '&
     &/)
      character(len=16), dimension(5) :: ename = (/' DENSITY        ',  &
     &' VELOCITY FIELD ',' PRESSURE TENSOR',' ENERGY         ',         &
     &' HEAT FLUX      '/)
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
! determine which fluid diagnostics are available
      call readfldiags2(iudm,nscalars)
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
         write (*,*) 'no fluid moments diagnostic files found'
         stop
      endif
!
      write (*,*) trim(dname(n)), ' diagnostic selected'
!
! return parameters for selected fluid diagnostic
      call fldiagparams2(iudm,n,nts,npro,nprd,nrec,fname)
      nplot = npro
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
      allocate(fms(nprd,nx,nyv),sfield(nx,nyv))
      if (ndim==3) then
         allocate(vfield(ndim,nx,nyv))
         if (npro > 2) allocate(ufield(2*ndim,nx,nyv))
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
! open graphics device
      ierr = open_graphs2(nplot)
! set palette to color wheel
      call set_palit(2)
!
! read and display data
      do ii = 1, nrec
! read real vector field
         call freadv2(iuv,fms,nprd,nx,nyv)
         it = nts*(ii - 1)
         time = dt*real(ii - 1)
! show time
         write (*,*) sname(n),'it,time=',it,time
! electrostatic case
         if (ndim==1) then
! display density in real space
            if (npro > 0) then
               sfield = fms(1,:,:)
               call dscaler2(sfield,trim(sname(n))//trim(ename(1)),it,  &
     &999,1,1,nx,ny,ierr)
               if (ierr==1) exit
            endif
! display velocity field in real space
            if (npro > 1) then
               sfield = fms(2,:,:)
               call dscaler2(sfield,trim(sname(n))//trim(ename(2)),it,  &
     &999,0,1,nx,ny,ierr)
               if (ierr==1) exit
            endif
! display pressure tensor in real space
            if (npro > 2) then
               sfield= fms(3,:,:)
               call dscaler2(sfield,trim(sname(n))//trim(ename(3)),it,  &
     &999,0,1,nx,ny,ierr)
               if (ierr==1) exit
            endif
! display electron heat flux in real space
            if (npro==4) then
               sfield = fms(5,:,:)
               call dscaler2(sfield,trim(sname(n))//trim(ename(5)),it,  &
     &999,0,1,nx,ny,ierr)
               if (ierr==1) exit
            endif
! electromagnetic case
         else if (ndim==3) then
! display density in real space
            if (npro > 0) then
               sfield = fms(1,:,:)
               call dscaler2(sfield,trim(sname(n))//trim(ename(1)),it,  &
     &999,1,1,nx,ny,ierr)
               if (ierr==1) exit
            endif
! display velocity field in real space
            if (npro > 1) then
               vfield  = fms(2:4,:,:)
               call dvector2(vfield,trim(sname(n))//trim(ename(2)),it,  &
     &999,0,1,1,nx,ny,ierr)
               if (ierr==1) exit
            endif
! display pressure tensor in real space
            if (npro > 2) then
               ufield= fms(5:10,:,:)
               call dvector2(ufield,trim(sname(n))//trim(ename(3)),it,  &
     &999,0,1,1,nx,ny,ierr)
               if (ierr==1) exit
            endif
! display electron heat flux in real space
            if (npro==4) then
               vfield = fms(12:14,:,:)
               call dvector2(vfield,trim(sname(n))//trim(ename(5)),it,  &
     &999,0,1,1,nx,ny,ierr)
               if (ierr==1) exit
            endif
         endif
      enddo
!
      call closeff2(iudm)
      call closeff2(iuv)
      call close_graphs2()
!
      end program
