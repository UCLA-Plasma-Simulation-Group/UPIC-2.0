!-----------------------------------------------------------------------
! This program reads real periodic 3d vector data
! written by 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program preadf3
      use in3, only: idrun, indx, indy, indz, nvpy, nvpz, ndim, tend,   &
     &dt, ci, omx, omy, omz,                                            &
     &el3d, ntel, felname, modesxel, modesyel, modeszel, nelrec,        &
     &vpot3d, nta, faname, modesxa, modesya, modesza, narec,            &
     &et3d, ntet, fetname, modesxet, modesyet, modeszet, netrec,        &
     &b3d, ntb, fbname, modesxb, modesyb, modeszb, nbrec,               &
     &vpotr3d, ntar, farname, modesxar, modesyar, modeszar, narrec,     &
     &vcuri3d, ntji, fjiname, modesxji, modesyji, modeszji, njirec
!
      implicit none
      integer, parameter :: ns = 6
      integer :: iudm = 19, iuv = 12
      integer :: i, j, k, l, m, n, ii, nx, ny, nz, kyp, kzp, kyb, kzb
      integer :: nyv, nzv, lrec, nrec, ios
      integer :: nplot = 1
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:,:), allocatable :: vfield
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
      fname = 'diag3.'//cdrun
      open(unit=iudm,file=fname,form='formatted',status='old')
!
! determine which vector diagnostics are available
      do n = 1, ns
      select case(n)
! load metadata for longitudinal efield data
      case (1)
         read (iudm,el3d,iostat=ios)
! load metadata for vector potential data
      case (2)
         read (iudm,vpot3d,iostat=ios)
! load metadata for transverse efield data
      case (3)
         read (iudm,et3d,iostat=ios)
! load metadata for magnetic field data
      case (4)
         read (iudm,b3d,iostat=ios)
! load metadata for radiative vector potential data
      case (5)
         read (iudm,vpotr3d,iostat=ios)
! load metadata for ion current density data
      case (6)
         read (iudm,vcuri3d,iostat=ios)
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
         read (iudm,el3d,iostat=ios)
         fname = felname; nrec = nelrec
      case (2)
         read (iudm,vpot3d,iostat=ios)
         fname = faname; nrec = narec
      case (3)
         read (iudm,et3d,iostat=ios)
         fname = fetname; nrec = netrec
      case (4)
         read (iudm,b3d,iostat=ios)
         fname = fbname; nrec = nbrec
      case (5)
         read (iudm,vpotr3d,iostat=ios)
         fname = farname; nrec = narrec
      case (6)
         read (iudm,vcuri3d,iostat=ios)
         fname = fjiname; nrec = njirec
      end select
      rewind iudm
      nplot = ndim
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
! open direct access file for vector field
      inquire(iolength=lrec) vfield(1,1,1,1)
      lrec = lrec*ndim*nx*nyv*nzv
      open(unit=iuv,file=fname,form='unformatted',access='direct',      &
     &recl=lrec,status='old')
!
! nrec = number of complete records
      nrec = nrec/(kyb*kzb)
      write (*,*) 'records found: nrec = ', nrec
!
! read and transpose vector data
      do ii = 1, nrec
         read (unit=iuv,rec=ii) ((((((vfield(i,j,k+kyp*(n-1),           &
     &l+kzp*(m-1)),i=1,ndim),j=1,nx),k=1,kyp),l=1,kzp),n=1,kyb),m=1,kzb)
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
