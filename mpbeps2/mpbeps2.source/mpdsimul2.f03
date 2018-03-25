!-----------------------------------------------------------------------
! High Level library for 2-2/2D Darwin MPI/OpenMP PIC code
!
! subroutines defined:
!
! bwrite_drestart13: write out basic restart file for darwin code
! bread_drestart13: read in basic restart file for darwin code
! dwrite_drestart13: write out restart diagnostic file for darwin code
! dread_drestart13: read in restart diagnostic file for darwin code
!
! written by Viktor K. Decyk, UCLA
! copyright 1999-2016, regents of the university of california
! update: october 30, 2017
      module mpdsimul2
      use mpsimul2
      use mpbsimul2
!     use mdpush1
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine bwrite_drestart23(cus,wpm,q2m0,tdiag,kstrt,iur)
! write out basic restart file for darwin code
      implicit none
! cus = smoothed transverse electric field
      real, dimension(:,:,:), intent(in) :: cus
! wpm = normalized total plasma frequency squared
! q2m0 = shift constant in darwin iteration
      real, intent(in) :: wpm, q2m0
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptors
      integer, intent(in) :: kstrt, iur
! local data
      integer :: ndim, nxv, nypmx
      ndim = size(cus,1); nxv = size(cus,2); nypmx = size(cus,3)
! write out shift constants for iteration
      if (kstrt==1) write (iur) wpm, q2m0
! write out darwin electric field
      if (kstrt==1) write (iur) ndim, nxv, nypmx
      call mpwrvdata2(cus,tdiag,iur)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_drestart23(cus,wpm,q2m0,tdiag,kstrt,iur,irc)
! read in basic restart file for darwin code
      implicit none
! cus = smoothed transverse electric field
      real, dimension(:,:,:), intent(inout) :: cus
! q2m0 = shift constant in darwin iteration
      real, intent(inout) :: wpm, q2m0
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptors
      integer, intent(in) :: kstrt, iur
      integer, intent(inout) :: irc
!
! local data
      integer :: ios, it, is, ir
      integer :: ndim, nxv, nypmx
      integer, dimension(1) :: kval = 0
      double precision, dimension(2) :: dval = 0.0d0
      irc = 0
!
      ndim = size(cus,1); nxv = size(cus,2); nypmx = size(cus,3)
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
         read (iur,iostat=ios) it, is, ir
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
         endif
      endif
! broadcast error condition
      call PPBICAST(kval,1)
      irc = kval(1)
      if (irc /= 0) return
! read in field array
      call mprdvdata2(cus,tdiag,iur,irc)
      end subroutine
!
      end module
