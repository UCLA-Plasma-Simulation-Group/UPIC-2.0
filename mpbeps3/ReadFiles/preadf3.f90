!-----------------------------------------------------------------------
! This program reads real periodic 3d scalar data
! written by 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program preadf3
      use in3, only: idrun, indx, indy, indz, nvpy, nvpz, tend, dt, ci, &
     &pot3d, ntp, fpname, modesxp, modesyp, modeszp, nprec,             &
     &dene3d, ntde, fdename, modesxde, modesyde, modeszde, nderec,      &
     &deni3d, ntdi, fdiname, modesxdi, modesydi, modeszdi, ndirec

!
      implicit none
      integer, parameter :: ns = 3
      integer :: iudm = 19, ius = 11
      integer :: i, j, k, l, m, n, ii, nx, ny, nz, kyp, kzp, kyb, kzb
      integer :: nyv, nzv, lrec, nrec, ios
      integer :: nplot = 1
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:), allocatable :: sfield
      character(len=16), dimension(ns) :: dname = (/'potential       ', &
     &'electron density','ion density     '/)
      character(len=10) :: cdrun
      character(len=32) :: fname
!
      write (*,*) 'enter idrun:'
      read (5,*) idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag3.'//cdrun
      open(unit=iudm,file=fname,form='formatted',status='old')
!
! determine which scalar diagnostics are available
      do n = 1, ns
      select case(n)
! load metadata for potential data
      case (1)
         read (iudm,pot3d,iostat=ios)
! load metadata for electron density data
      case (2)
         read (iudm,dene3d,iostat=ios)
! load metadata for ion density data
      case (3)
         read (iudm,deni3d,iostat=ios)
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
         read (iudm,pot3d,iostat=ios)
         fname = fpname; nrec = nprec
      case (2)
         read (iudm,dene3d,iostat=ios)
         fname = fdename; nrec = nderec
      case (3)
         read (iudm,deni3d,iostat=ios)
         fname = fdiname; nrec = ndirec
      end select
!
! nx
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
      allocate(sfield(nx,nyv,nzv))
! open direct access file for scalar field
      inquire(iolength=lrec) sfield(1,1,1); lrec = lrec*nx*nyv*nzv
      open(unit=ius,file=fname,form='unformatted',access='direct',      &
     &recl=lrec,status='old')
!
! nrec = number of complete records
      nrec = nrec/(kyb*kzb)
      write (*,*) 'records found: nrec = ', nrec
!
! read and transpose scalar data
      do ii = 1, nrec
         read (unit=ius,rec=ii) (((((sfield(j,k+kyp*(n-1),l+kzp*(m-1)), &
     &j=1,nx),k=1,kyp),l=1,kzp),n=1,kyb),m=1,kzb)
      enddo
!
      end program
!
! unneeded function in input3mod.f90
      subroutine PPBDCAST(f,nxp)
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f
      end subroutine
!
