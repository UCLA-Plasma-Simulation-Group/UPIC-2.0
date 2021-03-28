!-----------------------------------------------------------------------
! This program reads real periodic 1d velocity distribution data
! written for 1D OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program pfvreadf1
      use graf1
!
      implicit none
      integer, parameter :: nmv = 80, ndim = 3, ius = 64
      integer :: ii, it, is, ntime, nmv21, ierr
      real :: ci = 0.0
      real :: ci2, fmin, de, e, gam
      double precision :: dsum
      real, dimension(:,:), allocatable :: fv, dfv
      real, dimension(:), allocatable :: g
      character(len=6) :: chr
      character(len=10) :: cdrun
      character(len=32) :: fname
      character(len=10), dimension(2) :: chs = (/'INITIAL   ',          &
     &'FINAL     '/)
!     character(len=10), dimension(3) :: chs = (/'INITIAL   ',          &
!    &'FINAL     ','POST FINAL'/)
!
      nmv21 = 2*nmv + 1
      allocate(fv(nmv21+1,ndim),dfv(nmv21,2),g(nmv21))
!     allocate(fv(nmv21+1,ndim),dfv(nmv21,3),g(nmv21))
      ci2 = ci*ci
!
      write (cdrun,'(i10)') ius
      fname = 'fort.'//adjustl(cdrun)
      open(unit=ius,file=fname,form='formatted',status='old')
!
      ierr = open_graphs(1)

!     fname = 'LOG(GAMMA-1): WATERBAG/JUTTNER'
!     fname = 'LOG(GAMMA-1): MAXWELLIAN/JUTTNER'
      fname = 'LOG(GAMMA-1): JUTTNER/JUTTNER'
      do ii = 1, 1000
!     do ii = 1, 4000
!     do ii = 1, 1500
      read (ius,*) chr, ntime
      dsum = 0.0d0
      do it = 1, size(fv,1)
      read (ius,*) is,fv(it,1)
      if (it.gt.1) dsum = dsum + dble(fv(it,1))
      enddo
      if (ii.eq.1) write (*,*) 'sum=',real(dsum)
!
      if (ii==1) then
         de = fv(1,1)/real(nmv)
         do it = 1, nmv21
         e = de*(real(it)-0.5)
         gam = 1.0 + e*ci2
!        g(it) = gam/sqrt((gam+1.0)*e)
         g(it) = gam*sqrt((gam+1.0)*e)
         enddo
      endif
!     fmin = minval(fv(2:,1),fv(2:,1).gt.0.0)
      fmin = 1
      where (fv(2:,1) > 0.0)
         fv(2:,1) = alog(fv(2:,1)/g)
      else where
         fv(2:,1) = alog(fmin)
      end where
      if (ii==1) dfv(:,1) = fv(2:,1)
      if (ii==1000) dfv(:,2) = fv(2:,1)
!     if (ii==4000) dfv(:,2) = fv(2:,1)
!     if (ii==1500) dfv(:,3) = fv(3:,1)
      call dscaler1(fv(2:,1),fname,ntime,999,1,nmv21,ierr)
      if (ierr==1) exit
      enddo
!
      call dmscaler1(dfv,fname,ntime,999,1,nmv21,chs,ierr)
!
      write (50,*) 'Initial'
      do it = 1, nmv21
      write (50,*) it,dfv(it,1),exp(dfv(it,1)),0.5*g(it)**2
      enddo
!
      write (50,*) 'Final'
      do it = 1, nmv21
      write (50,*) it,dfv(it,2),exp(dfv(it,2)),0.5*g(it)**2
      enddo
!
!     write (50,*) 'Post Final'
!     do it = 1, nmv21
!     write (50,*) it,dfv(it,3),exp(dfv(it,3)),0.5*g(it)**2
!     enddo
!
      read (5,*) chr
!
      call close_graphs()
!
      end program
