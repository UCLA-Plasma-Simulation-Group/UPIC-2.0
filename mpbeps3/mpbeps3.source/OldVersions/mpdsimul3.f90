!-----------------------------------------------------------------------
! High Level library for 3D Darwin MPI/OpenMP PIC code
!
! subroutines defined:
!
! bwrite_drestart3: write out basic restart file for darwin code
! bread_drestart3: read in basic restart file for darwin code
!
! written by Viktor K. Decyk, UCLA
! copyright 1999-2016, regents of the university of california
! update: june 13, 2018
      module fd3
      use f3
      use fb3
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine bwrite_drestart3(cus,wpm,q2m0,tdiag,kstrt,iur)
! write out basic restart file for darwin code
      implicit none
! cus = smoothed transverse electric field
      real, dimension(:,:,:,:), intent(in) :: cus
! wpm = normalized total plasma frequency squared
! q2m0 = shift constant in darwin iteration
      real, intent(in) :: wpm, q2m0
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptors
      integer, intent(in) :: kstrt, iur
! local data
      integer :: ndim, nxv, nypmx, nzpmx
      ndim = size(cus,1); nxv = size(cus,2)
      nypmx = size(cus,3); nzpmx = size(cus,4)
! write out shift constants for iteration
      if (kstrt==1) write (iur) wpm, q2m0
! write out darwin electric field
      if (kstrt==1) write (iur) ndim, nxv, nypmx, nzpmx
      call mpwrvdata3(cus,tdiag,iur)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_drestart3(cus,wpm,q2m0,tdiag,kstrt,iur,irc)
! read in basic restart file for darwin code
      implicit none
! cus = smoothed transverse electric field
      real, dimension(:,:,:,:), intent(inout) :: cus
! wpm = normalized total plasma frequency squared
! q2m0 = shift constant in darwin iteration
      real, intent(inout) :: wpm, q2m0
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptors
      integer, intent(in) :: kstrt, iur
      integer, intent(inout) :: irc
!
! local data
      integer :: ios, it, is, ir, iq
      integer :: ndim, nxv, nypmx, nzpmx
      integer, dimension(1) :: kval = 0
      double precision, dimension(2) :: dval = 0.0d0
      irc = 0
!
      ndim = size(cus,1); nxv = size(cus,2)
      nypmx = size(cus,3); nzpmx = size(cus,4)
! read in shift constant for iteration
      if (kstrt==1) then
         read (iur,iostat=ios) wpm, q2m0
         if (ios /= 0) then
            write (*,*) 'wpm, q2m0 restart error, ios = ', ios
            kval(1) = 1
         endif
      endif
! broadcast error condition
      call PPBICAST(kval,1)
      irc = kval(1)
      if (irc /= 0) return
! broadcase shift constants
      dval(1) = wpm
      dval(2) = q2m0
      call PPBDCAST(dval,2)
      wpm = dval(1)
      q2m0 = dval(2)
! read in darwin electric field field
      if (kstrt==1) then
         read (iur,iostat=ios) it, is, ir, iq
         if (ios /= 0) then
            write (*,*) 'cus size restart error, ios = ', ios
            kval(1) = 1
         endif
         if (it /= ndim) then
            write (*,*) 'cus restart error, size(cus,1)=', it, ndim
            kval(1) = 1
         else if (is /= nxv) then
            write (*,*) 'cus restart error, size(cus,2)=', is, nxv
            kval(1) = 1
         else if (ir /= nypmx) then
            write (*,*) 'cus restart error, size(cus,3)=', ir, nypmx
            kval(1) = 1
         else if (iq /= nzpmx) then
            write (*,*) 'cus restart error, size(cus,4)=', iq, nzpmx
            kval(1) = 1
         endif
      endif
! broadcast error condition
      call PPBICAST(kval,1)
      irc = kval(1)
      if (irc /= 0) return
! read in field array
      call mprdvdata3(cus,tdiag,iur,irc)
      end subroutine
!
      end module
