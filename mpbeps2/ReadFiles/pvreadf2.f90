!-----------------------------------------------------------------------
! This program reads real periodic 2d vector data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program preadf2
      use in2, only: idrun, indx, indy, nvp, ndim, tend, dt,            &
     &ci, omx, omy, omz,                                                &
     &el2d, ntel, felname, modesxel, modesyel, nelrec,                  &
     &vpot2d, nta, faname, modesxa, modesya, narec,                     &
     &et2d, ntet, fetname, modesxet, modesyet, netrec,                  &
     &b2d, ntb, fbname, modesxb, modesyb, nbrec,                        &
     &vpotr2d, ntar, farname, modesxar, modesyar, narrec,               &
     &vcuri2d, ntji, fjiname, modesxji, modesyji, njirec
      use graf2
!
      implicit none
      integer, parameter :: ns = 6
      integer :: iudm = 19, iuv = 12
      integer :: i, n, m, nx, ny, kyp, kyb, nyv, lrec, nrec, ios, ierr
      integer :: nplot = 1
      integer, dimension(ns) :: nscalars = 0
      real, dimension(:,:,:), allocatable :: vfield
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
         fname = felname; nrec = nelrec
      case (2)
         read (iudm,vpot2d,iostat=ios)
         fname = faname; nrec = narec
      case (3)
         read (iudm,et2d,iostat=ios)
         fname = fetname; nrec = netrec
      case (4)
         read (iudm,b2d,iostat=ios)
         fname = fbname; nrec = nbrec
      case (5)
         read (iudm,vpotr2d,iostat=ios)
         fname = farname; nrec = narrec
      case (6)
         read (iudm,vcuri2d,iostat=ios)
         fname = fjiname; nrec = njirec
      end select
      rewind iudm
      nplot = ndim
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
! allocate vector array
      if (.not.allocated(vfield)) allocate(vfield(ndim,nx,nyv))
! open direct access file for vector field
      inquire(iolength=lrec) vfield(1,1,1); lrec = lrec*ndim*nx*nyv
      open(unit=iuv,file=fname,form='unformatted',access='direct',      &
     &recl=lrec,status='old')
!
! nrec = number of complete records
      nrec = nrec/kyb
      write (*,*) 'records found: nrec = ', nrec
!
! open graphics device
      ierr = open_graphs(nplot)
! set palette to color wheel
      call STPALIT(2)
!
! read and display data
      do i = 1, nrec
         read (unit=iuv,rec=i) vfield
         call dvector2(vfield,trim(dname(n)),i,999,0,1,1,nx,ny,ierr)
         if (ierr==1) exit
      enddo
!
      call close_graphs
!
      end program
!
! unneeded function in input2mod.f90
      subroutine PPBDCAST(f,nxp)
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f
      end subroutine
!
