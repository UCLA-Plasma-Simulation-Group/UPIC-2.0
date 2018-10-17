!-----------------------------------------------------------------------
! High Level library for 2D Electrostatic MPI/OpenMP PIC code
!
! subroutines defined:
!
! open_restart2: open restart files
! bwrite_restart2: write out basic restart file for electrostatic code
! bwrite_restart2x: write out basic restart file for electrostatic code
!                   with two ion species
! bread_restart2a: read in basic restart file for electrostatic code,
!                  headers and electron data
! bread_restart2b: read in basic restart file for electrostatic code,
!                  copy electron data and read in ion data if needed
! bread_restart2x: read in basic restart file for electrostatic code,
!                  copy ion data and read in second ion data if needed
! bread_restart2c: read in basic restart file for electrostatic code,
!                  copy ion data or read in fixed ion density
! dwrite_restart2: write out restart diagnostic file for electrostatic
!                  code
! dread_restart2: read in restart diagnostic file for electrostatic code
! close_restart2: close restart files
!
! written by Viktor K. Decyk, UCLA
! copyright 1999-2017, regents of the university of california
! update: august 15, 2018
      module f2
      use in2
      use minit2
      use mpush2
!     use msort1
!     use mgard1
!     use mfft1
!     use mfield1
      use mdiag2
      use mppmod2
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine open_restart2(iur,iur0,cdrun)
! open restart files
      implicit none
      integer, intent(in) :: iur0
      integer, intent(inout) :: iur
      character(len=10), intent(in) :: cdrun
! iur, iur0 = restart, old restart file descriptors
! cdrun = idrun string
! local data
      character(len=10) :: cdrun0
      character(len=32) :: fname
!     character(len=10) :: a = 'sequential'
      character(len=6) :: a = 'stream'
      character(len=11) :: f = 'unformatted'
!
! reset file
!     fname = 'reset2'
!     iurr = get_funit(iurr)
!     open(iurr,file=trim(fname),access=a,form=f,status='unknown')
! restart file
      fname = 'rstrt2.'//cdrun
      iur = get_funit(iur)
! start a new run from random numbers
      if (nustrt==1) then
         if (ntr > 0) then
            open(iur,file=trim(fname),access=a,form=f,status= 'unknown')
         endif
! continue a run which was interrupted
      else if (nustrt==2) then
         open(iur,file=trim(fname),access=a,form=f,status='old')
! start a new run with data from a previous run
      else if (nustrt==0) then
         if (ntr > 0) then
            open(iur,file=trim(fname),access=a,form=f,status= 'unknown')
         endif
         if (idrun /= idrun0) then
            write (cdrun0,'(i10)') idrun0
            cdrun0 = adjustl(cdrun0)
            fname = 'rstrt2.'//cdrun0
            open(iur0,file=trim(fname),access=a,form=f,status='old')
         else
            write (*,*) 'restart warning: old, new idruns identical'
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bwrite_restart2(part,ppart,pparti,qi,kpic,kipic,tdiag, &
     &kstrt,iur,iscr,ntime,ntime0,irc)
! write out basic restart file for electrostatic code
      implicit none
! part = particle array
      real, dimension(:,:), intent(inout) :: part
! ppart/pparti = tiled electron/ion particle arrays
      real, dimension(:,:,:), intent(in) :: ppart, pparti
! qi = ion charge density with guard cells
      real, dimension(:,:), intent(in) :: qi
! kpic/kipic = number of electrons/ions in each tile
      integer, dimension(:), intent(in) :: kpic, kipic
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur/iscr = restart/scratch file descriptors
! ntime/ntime0 = current/initial time step
      integer, intent(in) :: kstrt, iur, iscr, ntime, ntime0
      integer, intent(inout) :: irc
! local data
      integer :: idimp, npp, nppi, nxv, nypmx
      idimp = size(part,1)
      nxv = size(qi,1); nypmx = size(qi,2)
      if (kstrt==1) rewind iur
! write out current and initial time
      if (kstrt==1) write (iur) ntime, ntime0
! write out number of processors
      if (kstrt==1) write (iur) nvp
! copy ordered particles to linear array: updates part, npp, irc
      call mpcopyout2(part,ppart,kpic,npp,irc)
! write out size of electron array
      if (kstrt==1) write (iur) idimp
! write out electrons: updates part
      call mpwrpart2(part,tdiag,npp,iur,iscr)
! write out if ions are moving
      if (kstrt==1) write (iur) movion
      if (movion > 0) then
! copy ordered particles to linear array: updates part, nppi, irc
         call mpcopyout2(part,pparti,kipic,nppi,irc)
! write out size of ion array
         if (kstrt==1) write (iur) idimp
! write out ions: updates part
         call mpwrpart2(part,tdiag,nppi,iur,iscr)
! write out ion density, if ions are not moving
      else
         nxv = size(qi,1); nypmx = size(qi,2)
         if (kstrt==1) write (iur) nxv, nypmx
         call mpwrdata2(qi,tdiag,iur)
      endif
! write out electric field parameter
      if (kstrt==1) write (iur) emf
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bwrite_restart2x(part,ppart,pparti,pparti2,qi,kpic,    &
     &kipic,kipic2,tdiag,kstrt,iur,iscr,ntime,ntime0,irc)
! write out basic restart file for electrostatic code
! with two ions species
      implicit none
! part = particle array
      real, dimension(:,:), intent(inout) :: part
! ppart/pparti = tiled electron/ion particle arrays
      real, dimension(:,:,:), intent(in) :: ppart, pparti, pparti2
! qi = ion charge density with guard cells
      real, dimension(:,:), intent(in) :: qi
! kpic/kipic = number of electrons/ions in each tile
      integer, dimension(:), intent(in) :: kpic, kipic, kipic2
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur/iscr = restart/scratch file descriptors
! ntime/ntime0 = current/initial time step
      integer, intent(in) :: kstrt, iur, iscr, ntime, ntime0
      integer, intent(inout) :: irc
! local data
      integer :: idimp, npp, nppi, nppi2, nxv, nypmx
      idimp = size(part,1)
      nxv = size(qi,1); nypmx = size(qi,2)
      if (kstrt==1) rewind iur
! write out current and initial time
      if (kstrt==1) write (iur) ntime, ntime0
! write out number of processors
      if (kstrt==1) write (iur) nvp
! copy ordered particles to linear array: updates part, npp, irc
      call mpcopyout2(part,ppart,kpic,npp,irc)
! write out size of electron array
      if (kstrt==1) write (iur) idimp
! write out electrons: updates part
      call mpwrpart2(part,tdiag,npp,iur,iscr)
! write out if ions are moving
      if (kstrt==1) write (iur) movion
      if (movion > 0) then
! copy ordered particles to linear array: updates part, nppi, irc
         call mpcopyout2(part,pparti,kipic,nppi,irc)
! write out size of ion array
         if (kstrt==1) write (iur) idimp
! write out ions: updates part
         call mpwrpart2(part,tdiag,nppi,iur,iscr)
! write out second ions species
         if (movion > 1) then
! copy ordered particles to linear array: updates part, nppi2, irc
            call mpcopyout2(part,pparti2,kipic2,nppi2,irc)
! write out size of ion array
            if (kstrt==1) write (iur) idimp
! write out ions: updates part
            call mpwrpart2(part,tdiag,nppi2,iur,iscr)
         endif
! write out ion density, if ions are not moving
      else
         nxv = size(qi,1); nypmx = size(qi,2)
         if (kstrt==1) write (iur) nxv, nypmx
         call mpwrdata2(qi,tdiag,iur)
      endif
! write out electric field parameter
      if (kstrt==1) write (iur) emf
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_restart2a(part,kpic,tdiag,kstrt,iur,iscr,ntime,  &
     &ntime0,npp,nppmx,noff,mx1,irc)
! read in basic restart file for electrostatic code, first part:
! read in headers and electron data
      implicit none
! part = particle array
      real, dimension(:,:), intent(inout) :: part
! kpic = number of electrons in each tile
      integer, dimension(:), intent(inout) :: kpic
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur/iscr = restart/scratch file descriptors
! noff = lowermost global gridpoint in particle partition
! mx1 = (system length in x direction - 1)/mx + 1, where
! mx = number of grids in sorting cell in x
      integer, intent(in) :: kstrt, iur, iscr, noff, mx1
! ntime/ntime0 = current/initial time step
! npp = number of particles in partition
! nppmx = maximum number of particles in tile
      integer, intent(inout) :: ntime, ntime0, npp, nppmx, irc
!
! local data
      integer :: ndimp, ios, it
      integer, dimension(3) :: kval = 0
      irc = 0
!
! read in current and initial time
      if (kstrt==1) then
         rewind iur
         read (iur,iostat=ios) ntime, ntime0
         if (ios /= 0) then
            write (*,*) 'ntime restart error, ios = ', ios
            kval(1) = 1
         else
            kval(2) = ntime; kval(3) = ntime0
            write (*,*) 'restarting from ntime, idrun0=', ntime, idrun0
         endif
      endif
! read in number of processors
      if (kstrt==1) then
         read (iur,iostat=ios) it
         if (ios /= 0) then
            write (*,*) 'nvp restart error, ios = ', ios
            kval(1) = 1
         endif
         if (it /= nvp) then
            write (*,*) 'restart error, nvp=', it, nvp
            kval(1) = 1
         endif
      endif
! read in size of electron array
      if (kstrt==1) then
         read (iur,iostat=ios) ndimp
         if (ios /= 0) then
            write (*,*) 'idimp restart error, ios = ', ios
            kval(1) = 1
         endif
         if (ndimp /= size(part,1)) then
            write (*,*) 'restart error, idimp=', ndimp, size(part,1)
            kval(1) = 1
         endif
      endif
! broadcast time at restart and error condition
      call PPBICAST(kval,3)
      irc = kval(1); ntime = kval(2); ntime0 = kval(3); irc = kval(1)
      if (irc /= 0) return
! read in electrons: updates part, npp, irc
      call mprdpart2(part,tdiag,npp,iur,iscr,irc)
      if (irc /= 0) return
! find number of electrons in each of mx, my tiles: updates kpic, nppmx
      call mpdblkp2(part,kpic,npp,noff,nppmx,mx,my,mx1,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_restart2b(part,ppart,kpic,kipic,tdiag,kstrt,iur, &
     &iscr,npp,nppi,nppmx,noff,nyp,nx,mx1,irc)
! read in basic restart file for electrostatic code, second part:!
! copy electron data to ppart and read in ion data if needed
      implicit none
! part = particle array
      real, dimension(:,:), intent(inout) :: part
! ppart = tiled electron particle arrays
      real, dimension(:,:,:), intent(inout) :: ppart
! kpic/kipic = number of electrons/ions in each tile
      integer, dimension(:), intent(inout) :: kpic, kipic
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur/iscr = restart/scratch file descriptors
! npp = number of electrons in partition
! noff = lowermost global gridpoint in particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! nx = system length in x direction
! mx1 = (system length in x direction - 1)/mx + 1, where
! mx = number of grids in sorting cell in x
      integer, intent(in) :: kstrt, iur, iscr, npp, noff, nyp, nx, mx1
! nppi = number of ions in partition
! nppmx = maximum number of particles in tile
      integer, intent(inout) :: nppi, nppmx, irc
!
! local data
      integer :: ndimp, ios, it
      integer, dimension(1) :: kval = 0
      irc = 0
!
! copy ordered electron data for OpenMP: updates ppart and kpic
      call mpmovin2p(part,ppart,kpic,npp,noff,mx,my,mx1,irc)
! sanity check for electrons
      call mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,irc)
! read in movion to determine if ions are moving
      if (kstrt==1) then
         read (iur,iostat=ios) it
         if (ios /= 0) then
            write (*,*) 'movion restart error, ios = ', ios
            kval(1) = 1
         endif
         if (it /= movion) then
            write (*,*) 'movion restart error, movion = ', it, movion
            kval(1) = 1
         endif
      endif
! broadcast error condition
      call PPBICAST(kval,1)
      irc = kval(1)
      if (irc /= 0) return
! ions are moving
      if (movion > 0) then
! read in size of ion array
         if (kstrt==1) then
            read (iur,iostat=ios) ndimp
            if (ios /= 0) then
               write (*,*) 'ion idimp restart error, ios = ', ios
               kval(1) = 1
            endif
            if (ndimp /= size(part,1)) then
               write (*,*) 'ion restart error,idimp=',ndimp,size(part,1)
               kval(1) = 1
            endif
         endif
! broadcast error condition
         call PPBICAST(kval,1)
         irc = kval(1)
         if (irc /= 0) return
! read in ions: updates part, nppi, irc
         call mprdpart2(part,tdiag,nppi,iur,iscr,irc)
! find number of ions in each of mx, my tiles: updates kipic, nppmx
         call mpdblkp2(part,kipic,nppi,noff,nppmx,mx,my,mx1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_restart2x(part,pparti,kipic,kipic2,tdiag,kstrt,  &
     &iur,iscr,nppi,nppi2,nppmx,noff,nyp,nx,mx1,irc)
! read in basic restart file for electrostatic code, extra part:
! copy ion data to pparti and read in additional ion data if needed
      implicit none
! part = particle array
      real, dimension(:,:), intent(inout) :: part
! ppart = tiled electron particle arrays
      real, dimension(:,:,:), intent(inout) :: pparti
! kipic/kipic2 = number of ions/and additional ions in each tile
      integer, dimension(:), intent(inout) :: kipic, kipic2
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur/iscr = restart/scratch file descriptors
! nppi = number of ions in partition
! noff = lowermost global gridpoint in particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! nx = system length in x direction
! mx1 = (system length in x direction - 1)/mx + 1, where
! mx = number of grids in sorting cell in x
      integer, intent(in) :: kstrt, iur, iscr, nppi, noff, nyp, nx, mx1
! nppi = number of ions in partition
! nppmx = maximum number of particles in tile
      integer, intent(inout) :: nppi2, nppmx, irc
!
! local data
      integer :: ndimp, ios
      integer, dimension(1) :: kval = 0
      irc = 0
!
! ions are not moving
      if (movion==0) return
! copy ordered ion data for OpenMP: updates pparti and kipic
      call mpmovin2p(part,pparti,kipic,nppi,noff,mx,my,mx1,irc)
! sanity check for ions
      call mpcheck2(pparti,kipic,noff,nyp,nx,mx,my,mx1,irc)
! additional ions are present
      nppi2 = 0
      if (movion > 1) then
! read in size of ion array
         if (kstrt==1) then
            read (iur,iostat=ios) ndimp
            if (ios /= 0) then
               write (*,*) 'ion idimp restart error, ios = ', ios
               kval(1) = 1
            endif
            if (ndimp /= size(part,1)) then
               write (*,*) 'ion restart error,idimp=',ndimp,size(part,1)
               kval(1) = 1
            endif
         endif
! broadcast error condition
         call PPBICAST(kval,1)
         irc = kval(1)
         if (irc /= 0) return
! read in additional ions: updates part, nppi2, irc
         call mprdpart2(part,tdiag,nppi2,iur,iscr,irc)
! find number of ions in each of mx, my tiles: updates kipic2, nppmx
         call mpdblkp2(part,kipic2,nppi2,noff,nppmx,mx,my,mx1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_restart2c(part,pparti,kipic,qi,tdiag,kstrt,iur,  &
     &ntime,ntime0,nppi,noff,nyp,nx,mx1,irc)
! read in basic restart file for electrostatic code, third part:
! copy ion data to pparti or read in fixed ion density
      implicit none
! part = particle array
      real, dimension(:,:), intent(in) :: part
! pparti = tiled ion particle arrays
      real, dimension(:,:,:), intent(inout) :: pparti
! kipic = number of ions in each tile
      integer, dimension(:), intent(inout) :: kipic
! qi = ion charge density with guard cells
      real, dimension(:,:), intent(inout) :: qi
      real, intent(inout) :: tdiag
! kstrt = starting data block number
! iur = restart file descriptor
! nppi = number of ions in partition
! noff = lowermost global gridpoint in particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! nx = system length in x direction
! mx1 = (system length in x direction - 1)/mx + 1, where
! mx = number of grids in sorting cell in x
      integer, intent(in) :: kstrt, iur, nppi, noff, nyp, nx, mx1
! ntime/ntime0 = current/initial time step
      integer, intent(inout) :: ntime, ntime0, irc
!
! local data
      integer :: ios, it, is
      integer, dimension(1) :: kval = 0
      irc = 0
!
! ions are moving
      if (movion > 0) then
! copy ordered ion data for OpenMP: updates pparti and kipic
         call mpmovin2p(part,pparti,kipic,nppi,noff,mx,my,mx1,irc)
! sanity check for ions
         call mpcheck2(pparti,kipic,noff,nyp,nx,mx,my,mx1,irc)
! ions are not moving, read in ion density qi
      else
         if (kstrt==1) then
            read (iur,iostat=ios) it, is
            if (ios /= 0) then
               write (*,*) 'qi size restart error, ios = ', ios
               kval(1) = 1
            endif
            if (it /= size(qi,1)) then
               write (*,*) 'qi restart error, size(qi,1)=',it,size(qi,1)
               kval(1) = 1
            else if (is /= size(qi,2)) then
               write (*,*) 'qi restart error, size(qi,2)=',is,size(qi,2)
               kval(1) = 1
            endif
         endif
         call mprddata2(qi,tdiag,iur,irc)
      endif
! read in electric field parameter
      if (kstrt==1) then
         read (iur,iostat=ios) it
         if (ios /= 0) then
            write (*,*) 'emf restart error, ios = ', ios
            kval(1) = 1
         endif
         if (it /= emf) then
            write (*,*) 'warning: emf values differ, emf=',it,emf
         endif
      endif
! broadcast error condition
      call PPBICAST(kval,1)
      irc = kval(1)
      ntime0 = ntime0 + ntime
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine close_restart2(iur,iur0)
! close restart files
      implicit none
      integer, intent(in) :: iur, iur0
! iur, iur0 = restart, old restart file descriptors
!     close(unit=iurr)
      if (nustrt==1) then
         if (ntr > 0) close(unit=iur)
      else if (nustrt==2) then
         close(unit=iur)
      else if (nustrt==0) then
         if (ntr > 0) close(unit=iur)
         if (idrun /= idrun0) close(unit=iur0)
      endif
      end subroutine
!
      end module
