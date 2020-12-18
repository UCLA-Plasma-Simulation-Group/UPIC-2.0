!-----------------------------------------------------------------------
! High Level library for 3D Electromagnetic MPI/OpenMP PIC code
!
! subroutines defined:
!
! bwrite_restart3: write out basic restart file for electromagnetic code
! bread_restart3:  read in basic restart file for electromagnetic code
! dwrite_restart3: write out restart diagnostic file for
!                  electromagnetic code
! dread_restart3: read in restart diagnostic file for electromagnetic
!                 code
!
! written by Viktor K. Decyk, UCLA
! copyright 1999-2017, regents of the university of california
! update: may 16, 2018
      module fb3
      use f3
!     use mbpush3
!     use mcurd3
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine bwrite_brestart3(exyz,bxyz,tdiag,kstrt,iur)
! write out basic restart file for electromagnetic code
      implicit none
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:,:), intent(in) :: exyz, bxyz
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptors
      integer, intent(in) :: kstrt, iur
! local data
      integer :: ndim, nzv, kxypd, kyzpd
      ndim = size(exyz,1); nzv = size(exyz,2); kxypd = size(exyz,3)
      kyzpd = size(exyz,4)
! write out electromagnetic fields
      if (kstrt==1) write (iur) ndim, nzv, kxypd, kyzpd 
      call mpwrvcdata3(exyz,tdiag,iur)
      call mpwrvcdata3(bxyz,tdiag,iur)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_restart3(exyz,bxyz,tdiag,kstrt,iur,irc)
! read in basic restart file for electromagnetic code
      implicit none
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:,:), intent(inout) :: exyz, bxyz
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptors
      integer, intent(in) :: kstrt, iur
      integer, intent(inout) :: irc
!
! local data
      integer :: ios, it, is, ir, iq
      integer :: ndim, nzv, kxypd, kyzpd
      integer, dimension(1) :: kval = 0
      irc = 0
!
      ndim = size(exyz,1); nzv = size(exyz,2); kxypd = size(exyz,3)
      kyzpd = size(exyz,4)
! read in electromagnetic fields
      if (kstrt==1) then
         read (iur,iostat=ios) it, is, ir, iq
         if (ios /= 0) then
            write (*,*) 'exyz/bxyz size restart error, ios = ', ios
            kval(1) = 1
         endif
         if (it /= ndim) then
            write (*,*) 'exyz/bxyz restart error, size(exyz,1)=',it,ndim
            kval(1) = 1
         else if (is /= nzv) then
            write (*,*) 'exyz/bxyz restart error, size(exyz,2)=',is,nzv
            kval(1) = 1
         else if (ir /= kxypd) then
            write (*,*) 'exyz/bxyz restart error, size(exyz,3)=',ir,    &
     &kxypd
            kval(1) = 1
         else if (iq /= kyzpd) then
            write (*,*) 'exyz/bxyz restart error, size(exyz,4)=',iq,    &
     &kyzpd
            kval(1) = 1
         endif
      endif
! broadcast error condition
      call PPBICAST(kval,1)
      irc = kval(1)
      if (irc /= 0) return
! read in field arrays
      call mprdvcdata3(exyz,tdiag,iur,irc)
      call mprdvcdata3(bxyz,tdiag,iur,irc)
      end subroutine
!
      end module
