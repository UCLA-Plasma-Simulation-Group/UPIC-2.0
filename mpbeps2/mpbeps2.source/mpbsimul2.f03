!-----------------------------------------------------------------------
! High Level library for 2-1/2D Electromagnetic MPI/OpenMP PIC code
!
! subroutines defined:
!
! bwrite_restart23: write out basic restart file for electromagnetic
!                   code
! bread_restart23:  read in basic restart file for electromagnetic code
! dwrite_restart23: write out restart diagnostic file for
!                   electromagnetic code
! dread_restart23: read in restart diagnostic file for electromagnetic
!                  code
!
! written by Viktor K. Decyk, UCLA
! copyright 1999-2017, regents of the university of california
! update: august 6, 2018
      module fb2
      use f2
!     use mbpush1
!     use mcurd1
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine bwrite_restart23(exyz,bxyz,tdiag,kstrt,iur)
! write out basic restart file for electromagnetic code
      implicit none
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:), intent(in) :: exyz, bxyz
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptors
      integer, intent(in) :: kstrt, iur
! local data
      integer :: ndim, nyv, kxpd
      ndim = size(exyz,1); nyv = size(exyz,2); kxpd = size(exyz,3)
! write out electromagnetic fields
      if (kstrt==1) write (iur) ndim, nyv, kxpd
      call mpwrvcdata2(exyz,tdiag,iur)
      call mpwrvcdata2(bxyz,tdiag,iur)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_restart23(exyz,bxyz,tdiag,kstrt,iur,irc)
! read in basic restart file for electromagnetic code
      implicit none
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:), intent(inout) :: exyz, bxyz
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptors
      integer, intent(in) :: kstrt, iur
      integer, intent(inout) :: irc
!
! local data
      integer :: ios, it, is, ir
      integer :: ndim, nyv, kxpd
      integer, dimension(1) :: kval = 0
      irc = 0
!
      ndim = size(exyz,1); nyv = size(exyz,2); kxpd = size(exyz,3)
! read in electromagnetic fields
      if (kstrt==1) then
         read (iur,iostat=ios) it, is, ir
         if (ios /= 0) then
            write (*,*) 'exyz/bxyz size restart error, ios = ', ios
            kval(1) = 1
         endif
         if (it /= ndim) then
            write (*,*) 'exyz/bxyz restart error, size(exyz,1)=',it,ndim
            kval(1) = 1
         else if (is /= nyv) then
            write (*,*) 'exyz/bxyz restart error, size(exyz,2)=',is,nyv
            kval(1) = 1
         else if (ir /= kxpd) then
            write (*,*) 'exyz/bxyz restart error, size(exyz,3)=',ir,kxpd
            kval(1) = 1
         endif
      endif
! broadcast error condition
      call PPBICAST(kval,1)
      irc = kval(1)
      if (irc /= 0) return
! read in field arrays
      call mprdvcdata2(exyz,tdiag,iur,irc)
      call mprdvcdata2(bxyz,tdiag,iur,irc)
      end subroutine
!
      end module
