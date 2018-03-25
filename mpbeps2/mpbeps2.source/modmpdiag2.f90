!-----------------------------------------------------------------------
!
      module mpdiag2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpdiag2.f
! get_funit returns an unconnected fortran unit number
! fnrecl find record length of direct access file
! dafopen2 opens new binary file for real 2d scalar data.
! dafopenv2 opens new binary file for real 2d vector data.
! dafopenc2 opens new binary file for complex 2d scalar data.
! dafopenvc2 opens new binary file for complex 2d vector data.
! dafopenfv2 opens new binary file for real velocity-space data.
! dafopentr2 opens new binary file for real trajectory data.
! mpvpdist2 calculates 2 or 3 component velocity distribution, velocity
!           moments, and entropy with segmented particle array
!           calls PPVDIST2 or PPVDIST23
! mperpdist2 calculates 2d or 2-1/2d energy distribution for
!            relativistic particles with segmented particle array
!            calls PERPDIST2 or PERPDIST23
! mpvbpdist2 calculates 3d velocity distribution and velocity moments
!            for magnetized plasma with segmented particle array
!            calls PPVBDIST23
! mpvspdist2 calculates 2d velocity distribution in different regions of
!            space
!            calls PPVSDIST2 or PPVSDIST23
! mpvdist2 calculates 2 or 3 component velocity distribution, velocity
!          moments, and entropy with standard particle array
!          calls PVDIST2 or PVDIST23
! mpvbdist2 calculates 3d velocity distribution and velocity moments
!          for magnetized plasma with standard particle array
! mpprofx23 calculates fluid moments from particle quantities
!           assumes particle positions and velocities at same time level
!           for 2-1/2d code
!           calls PPROFX23L
! mprprofx23 calculates fluid moments from relativistic particle
!            quantities.  assumes particle positions and momenta at same
!            time level, for 2-1/2d code
!            calls PRPROFX23L
! mpprofx22 calculates fluid moments from particle quantities
!           assumes particle positions and velocities at same time level
!           for 2d code
!           calls PPROFX22L
! mprprofx22 calculates fluid moments from relativistic particle
!            quantities.  assumes particle positions and momenta at same
!            time level, for 2d code
!            calls PRPROFX22L
! mpgprofx2 calculates fluid moments from particle quantities
!           and electrostatic fields assumes particle positions and
!           velocities not at same time levels
!           calls PGPROFX2L
! mpgrprofx2 calculates fluid moments from relativistic particle
!            quantities and electrostatic fields.  assumes particle
!            positions and velocities not at same time levels
!            calls PGRPROFX2L
! mpgbprofx2 calculates fluid moments from particle quantities and 
!            electromagnetic fields assumes particle positions and
!            velocities not at same time levels
!            calls PGBPROFX23L
! mpgrbprofx2 calculates fluid moments from relativistic particle
!             quantities. assumes particle positions and velocities
!             not at same time levels
!             calls PGRBPROFX23L
! mpfluidqs23 calculates fluid quantities from fluid moments
!             for 2-1/2d code
!             calls FLUIDQS23
! mpfluidqs22 calculates fluid quantities from fluid moments for 2d code
!             calls FLUIDQS22
! wmprofx23 generic procedure to calculate fluid moments
!           calls mpprofx23 or mprprofx23
! wmprofx2 generic procedure to calculate fluid moments
!          calls mpprofx22 or mprprofx22
! wmgprofx2 generic procedure to calculate fluid moments with
!           electrostatic fields
!           calls mpgprofx2 or mpgrprofx2
! wmgbprofx2 generic procedure to calculate fluid moments with
!            electromagnetic fields
!            calls mpgbprofx2 or mpgrbprofx2
! psettraj2 sets test charge distribution by setting a particle id in
!           particle location 5 or 6
!           calls PSTPTRAJ2 or PSTPTRAJ23
! mptraj2 copies tagged particles in ppart to array partt
!         calls PPTRAJ2 or PPTRAJ23
! mpordtraj2 reorders tagged particles in partt to array spartt
!            calls PORDTRAJ2 or PORDTRAJ23
! mpcpytraj2 copies tagged particles in partt to array part
!            calls PCPYTRAJ2
! dafwritefv2 writes velocity-space record in direct access binary file
! dafwritetr2  writes trajectory record in direct access binary file
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: march 9, 2018
!
      use mppmod2, only: mpsum2
      use libmpdiag2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      function get_funit(start) result(funit)
! this function returns an unconnected fortran unit number,
! starting with unit = start.  returns -1 if none found
      implicit none
      integer, intent(in) :: start
      integer :: funit
! local data
      integer :: i
      logical :: connected
      funit = -1
! check connection status
      do i = start, 99
         inquire(unit=i,opened=connected)
         if (.not.connected) then
            funit = i
            exit
         endif
      enddo
      end function
!
!-----------------------------------------------------------------------
      function fnrecl(fname) result(it)
! find record length of direct access file
      implicit none
      character(len=*), intent(in) :: fname
      integer :: it, ios
      inquire(file=fname,recl=it,iostat=ios)
      if (ios /= 0) it = 0
      end function
!
!-----------------------------------------------------------------------
      subroutine dafopen2(f,nx,kyp,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 2d real scalar data array to be written in each record
! nx/kyp = number of data elements per record to be written in x/y 
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx, kyp
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1); lrec = lrec*nx*kyp
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenv2(f,nx,kyp,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 2d real vector data array to be written in each record
! nx/kyp = number of data elements per record to be written in x/y 
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx, kyp
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1,1); lrec = lrec*size(f,1)*nx*kyp
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenc2(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = 2d complex scalar data array to be written in each record
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      complex, dimension(:,:), intent(in) :: fc
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fc(1,1); lrec = lrec*size(fc)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenvc2(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = 2d complex vector data array to be written in each record
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      complex, dimension(:,:,:), intent(in) :: fc
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fc(1,1,1); lrec = lrec*size(fc)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenfv2(fvm,fv,fe,wk,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fvm,fv,fe = real velocity data to be written in each record
! wk = total energy contained in distribution
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      real, intent(in) :: wk
      real, dimension(:,:), intent(inout) :: fv, fvm, fe
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec, lwrec
      inquire(iolength=lrec) fvm(1,1)
      inquire(iolength=lwrec) wk
      lrec = lrec*(size(fvm) + size(fv) + size(fe)) + lwrec
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopentr2(partt,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! partt = real trajectory data to be written in each record 
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:), intent(in) :: partt
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) partt(1,1); lrec = lrec*size(partt)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenps2(fvs,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fvs = spatially resolved distribution function
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:,:,:), intent(inout) :: fvs
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fvs(1,1,1,1)
      lrec = lrec*size(fvs)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvpdist2(ppart,kpic,fv,sfv,fvm,tdiag,nvp,nmv)
! calculates 2d velocity distribution, velocity moments, and entropy
! with segmented particle array
      implicit none
      integer, intent(in) :: nvp, nmv
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fv, fvm
      real, dimension(:,:,:), intent(inout) :: sfv
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nmvf, idimv, mxyp1
      integer, dimension(4) :: itime
      real, dimension(3) :: scale
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(fv,1); idimv = size(fv,2)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimv==2) then
         call PPVDIST2(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,mxyp1,nmv, &
     &nmvf)
      else if (idimv==3) then
         call PPVDIST23(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,mxyp1,nmv,&
     &nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! sum distribution over processors, but save velocity scale
      scale(1:idimv) = fv(2*nmv+2,:)
      call mpsum2(fv,tdiag)
      fv(2*nmv+2,:) = scale(1:idimv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mperpdist2(ppart,kpic,fv,sfv,ci,wk,tdiag,ndim,nmv)
! calculates 2d or 2-1/2d energy distribution for relativistic particles
! with segmented particle array
      implicit none
      integer, intent(in) :: ndim, nmv
      real, intent(in) :: ci
      real, intent(inout) :: wk, tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fv
      real, dimension(:,:,:), intent(inout) :: sfv
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nmvf, idimv, mxyp1
      integer, dimension(4) :: itime
      real, dimension(3) :: scale
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(fv,1); idimv = size(fv,2)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (ndim==2) then
         call PERPDIST2(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyp1,nmv,  &
     &nmvf)
      else if (ndim==3) then
         call PERPDIST23(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyp1,nmv, &
     &nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! sum distribution and energy over processors, but save energy scale
      scale(1:idimv) = fv(2*nmv+2,:)
      call mpsum2(fv,tdiag)
      fv(2*nmv+2,:) = scale(1:idimv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvbpdist2(ppart,kpic,fv,sfv,fvm,omx,omy,omz,tdiag,nmv)
! calculates 3d velocity distribution and velocity moments
! for magnetized plasma with segmented particle array
      implicit none
      integer, intent(in) :: nmv
      real, intent(in) :: omx, omy, omz
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fv, fvm
      real, dimension(:,:,:), intent(inout) :: sfv
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nmvf,  mxyp1
      integer, dimension(4) :: itime
      real, dimension(2) :: scale
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(fv,1)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPVBDIST23(ppart,kpic,fv,sfv,fvm,omx,omy,omz,idimp,nppmx,    &
     &mxyp1,nmv,nmvf)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! sum distribution over processors, but save velocity scale
      scale = fv(2*nmv+2,:)
      call mpsum2(fv,tdiag)
      fv(2*nmv+2,:) = scale
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvspdist2(ppart,kpic,fvs,tdiag,noff,nmv,mvx,mvy)
! calculates 2d velocity distribution in different regions of space
      implicit none
      integer, intent(in) :: noff, nmv, mvx, mvy
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: fvs
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nmvf, idimv, nxb, nyb, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(fvs,1); idimv = size(fvs,2)
      nxb = size(fvs,3); nyb = size(fvs,4) - 1
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimv==2) then
         call PPVSDIST2(ppart,kpic,fvs,noff,nmv,mvx,mvy,nxb,nyb,idimp,  &
     &nppmx,mxyp1,nmvf)
      else if (idimv==3) then
         call PPVSDIST23(ppart,kpic,fvs,noff,nmv,mvx,mvy,nxb,nyb,idimp, &
     &nppmx,mxyp1,nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvdist2(part,fv,fvm,tdiag,npp,nvp,nmv)
! calculates 2d velocity distribution, velocity moments, and entropy
! with standard particle array
      implicit none
      integer, intent(in) :: npp, nvp, nmv
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:), intent(inout) :: fv, fvm
! local data
      integer :: idimp, npmax, nmvf, idimv
      integer, dimension(4) :: itime
      real, dimension(3) :: scale
      double precision :: dtime
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nmvf = size(fv,1); idimv = size(fv,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimv==2) then
         call PVDIST2(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
      else if (idimv==3) then
         call PVDIST23(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! sum distribution over processors, but save velocity scale
      scale(1:idimv) = fv(2*nmv+2,:)
      call mpsum2(fv,tdiag)
      fv(2*nmv+2,:) = scale(1:idimv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvbdist2(part,fv,fvm,omx,omy,omz,tdiag,npp,nmv)
! calculates 3d velocity distribution and velocity moments
! for magnetized plasma with standard particle array
      implicit none
      integer, intent(in) :: npp, nmv
      real, intent(in) :: omx, omy, omz
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:), intent(inout) :: fv, fvm
! local data
      integer :: idimp, npmax, nmvf
      integer, dimension(4) :: itime
      real, dimension(2) :: scale
      double precision :: dtime
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nmvf = size(fv,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PVBDIST23(part,fv,fvm,omx,omy,omz,npp,idimp,npmax,nmv,nmvf)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! sum distribution over processors, but save velocity scale
      scale = fv(2*nmv+2,:)
      call mpsum2(fv,tdiag)
      fv(2*nmv+2,:) = scale
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpprofx23(ppart,fms,kpic,noff,tdiag,npro,mx,my,mx1)
! calculates fluid moments from particle quantities
! assumes particle positions and velocities at same time level
! for 2-1/2d code
      implicit none
      integer, intent(in) :: noff, npro, mx, my, mx1
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPROFX23L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,my,nprd,nxv&
     &,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprprofx23(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mx1)
! calculates fluid moments from relativistic particle quantities
! assumes particle positions and momenta at same time level
! for 2-1/2d code
      implicit none
      integer, intent(in) :: noff, npro, mx, my, mx1
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PRPROFX23L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro,mx,my,nprd&
     &,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpprofx22(ppart,fms,kpic,noff,tdiag,npro,mx,my,mx1)
! calculates fluid moments from particle quantities
! assumes particle positions and velocities at same time level
! for 2d code
      implicit none
      integer, intent(in) :: noff, npro, mx, my, mx1
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPROFX22L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,my,nprd,nxv&
     &,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprprofx22(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mx1)
! calculates fluid moments from relativistic particle quantities
! assumes particle positions and momenta at same time level
! for 2d code
      implicit none
      integer, intent(in) :: noff, npro, mx, my, mx1
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PRPROFX22L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro,mx,my,nprd&
     &,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgprofx2(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,tdiag,npro&
     &,nx,mx,my,mx1)
! calculates fluid moments from particle quantities
! and electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: noff, nyp, npro, nx, mx, my, mx1
      real, intent(in) :: qbm, dt
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PGPROFX2L(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,idimp,nppmx,npro&
     &,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrprofx2(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,ci,tdiag,&
     &npro,nx,mx,my,mx1)
! calculates fluid moments from relativistic particle quantities
! and electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: noff, nyp, npro, nx, mx, my, mx1
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PGRPROFX2L(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,ci,idimp,nppmx,&
     &npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgbprofx2(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,tdiag&
     &,npro,nx,mx,my,mx1)
! calculates fluid moments from particle quantities
! and electromagnetic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: noff, nyp, npro, nx, mx, my, mx1
      real, intent(in) :: qbm, dt
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PGBPROFX23L(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,idimp,    &
     &nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrbprofx2(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,ci, &
     &tdiag,npro,nx,mx,my,mx1)
! calculates fluid moments from relativistic particle quantities
! and electromagnetic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: noff, nyp, npro, nx, mx, my, mx1
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PGRBPROFX23L(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,ci,idimp,&
     &nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfluidqs23(fms,tdiag,npro,nx,ny,kstrt,kyp)
! calculates fluid quantities from fluid moments for 2-1/2d code
      implicit none
      integer, intent(in) :: npro, nx, ny, kstrt, kyp
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(inout) :: fms
! local data
      integer :: nprd, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call FLUIDQS23(fms,npro,nx,ny,kstrt,kyp,nprd,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfluidqs22(fms,tdiag,npro,nx,ny,kstrt,kyp)
! calculates fluid quantities from fluid moments for 2d code
      implicit none
      integer, intent(in) :: npro, nx, ny, kstrt, kyp
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(inout) :: fms
! local data
      integer :: nprd, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nprd = size(fms,1); nxv = size(fms,2); nypmx = size(fms,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call FLUIDQS22(fms,npro,nx,ny,kstrt,kyp,nprd,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmprofx23(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mx1, &
     &relativity)
! generic procedure to calculate fluid moments
! assumes particle positions and velocities at same time levels
! for 2-1/2d code
      implicit none
      integer, intent(in) :: noff, npro, mx, my, mx1
      integer, intent(in) :: relativity
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! updates fms
      if (relativity==1) then
         call mprprofx23(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mx1)
      else
         call mpprofx23(ppart,fms,kpic,noff,tdiag,npro,mx,my,mx1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmprofx2(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mx1,  &
     &relativity)
! generic procedure to calculate fluid moments
! assumes particle positions and velocities at same time levels
! for 2 code
      implicit none
      integer, intent(in) :: noff, npro, mx, my, mx1
      integer, intent(in) :: relativity
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! updates fms
      if (relativity==1) then
         call mprprofx22(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mx1)
      else
         call mpprofx22(ppart,fms,kpic,noff,tdiag,npro,mx,my,mx1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgprofx2(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,ci,tdiag, &
     &npro,nx,mx,my,mx1,relativity)
! generic procedure to calculate fluid moments with electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: noff, nyp, npro, nx, mx, my, mx1
      integer, intent(in) :: relativity
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! updates fms
      if (relativity==1) then
         call mpgrprofx2(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,ci,tdiag,   &
     &npro,nx,mx,my,mx1)
      else
         call mpgprofx2(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,tdiag,npro,  &
     &nx,mx,my,mx1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgbprofx2(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,ci,  &
     &tdiag,npro,nx,mx,my,mx1,relativity)
! generic procedure to calculate fluid moments with electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: noff, nyp, npro, nx, mx, my, mx1
      integer, intent(in) :: relativity
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      real, dimension(:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! updates fms
      if (relativity==1) then
         call mpgrbprofx2(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,ci,    &
     &tdiag,npro,nx,mx,my,mx1)
      else
         call mpgbprofx2(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,tdiag,  &
     &npro,nx,mx,my,mx1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine psetptraj2(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,vtsx,&
     &dvtx,np,nprobt,irc)
! sets test charge distribution by setting a particle id in particle
! location 6
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! nprobt = number of test charges whose trajectories will be stored.
! irc = (0,1) = (no,yes) error condition exists
      integer, intent(in) :: kstrt, nst
      integer, intent(inout) :: nprobt, irc
      real, intent(in) :: vtx, vtsx, dvtx
      double precision, intent(in) :: np
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:), intent(inout) :: tedges
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(inout) :: iprobt
! local data
      integer :: idimp, nppmx, mxyp1, idps
      irc = 0
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(kpic,1); idps = size(tedges,1)
! call low level procedure
      if (idimp > 5) then
         call PSTPTRAJ23(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,vtsx,   &
     &dvtx,idimp,nppmx,mxyp1,idps,np,nprobt)
      else if (idimp > 4) then
         call PSTPTRAJ2(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,vtsx,dvtx&
     &,idimp,nppmx,mxyp1,idps,np,nprobt)
      else
         if (kstrt==1) then
            write (*,*) 'psetptraj2 error: idimp=', idimp
            irc = 1
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mptraj2(ppart,kpic,partt,tdiag,numtp,irc)
! this copies tagged particles in ppart to array partt
! numtp = number of test particles found on this node
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(inout) :: numtp, irc
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:,:), intent(inout) :: partt
! local data
      integer :: idimp, nppmx, mxyp1, nprobt
      integer, dimension(4) :: itime
      double precision :: dtime
      irc = 0
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(kpic,1); nprobt = size(partt,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimp > 5) then
         call PPTRAJ23(ppart,kpic,partt,numtp,idimp,nppmx,mxyp1,nprobt)
         if (numtp > nprobt) irc = -1
      else if (idimp > 4) then
         call PPTRAJ2(ppart,kpic,partt,numtp,idimp,nppmx,mxyp1,nprobt)
         if (numtp > nprobt) irc = -2
      else
         write (*,*) 'mptraj2 error: idimp=', idimp
         irc = 1
      endif
      if (irc < 0) then
         write (*,*) 'info:mptraj2 overflow: numpt,nprobt=',numtp,nprobt
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpordtraj2(partt,spartt,tedges,tdiag,numtp,irc)
! this procedure reorders tagged particles in partt to array spartt
      implicit none
      integer, intent(in) :: numtp
      integer, intent(inout) :: irc
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: partt
      real, dimension(:,:), intent(inout) :: spartt
      real, dimension(:), intent(in) :: tedges
! local data
      integer :: idimp, nprobt, idps
      integer, dimension(4) :: itime
      double precision :: dtime
      irc = 0
! extract dimensions
      idimp = size(partt,1); nprobt = size(partt,2)
      idps = size(tedges,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimp > 5) then
         call PORDTRAJ23(partt,spartt,tedges,numtp,idimp,idps,nprobt)
      else if (idimp > 4) then
         call PORDTRAJ2(partt,spartt,tedges,numtp,idimp,idps,nprobt)
      else
         write (*,*) 'mpordtraj2 error: idimp=', idimp
         irc = 1
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcpytraj2(partt,part,tdiag,numtp)
! this procedure copies tagged particles in partt to array part
      implicit none
      integer, intent(in) :: numtp
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: partt
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nprobt
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(partt,1); nprobt = size(partt,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PCPYTRAJ2(partt,part,numtp,idimp,nprobt)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwritefv2(fvm,fv,fe,wk,tdiag,iunit,nrec)
! this subroutine writes real velocity record in direct access binary
! file
! fvm,fv,fe = real velocity data to be written in each record
! wk = total energy contained in distribution
! iunit = fortran unit number to be used 
! nrec = record number for write (then updated to next record)
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: nrec
      real, intent(in) :: wk
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: fvm, fv, fe
! local data
      integer :: j, k, ndim, nmvf, nfvd, nfed
      integer, dimension(4) :: itime
      double precision :: dtime
      ndim = size(fvm,1)
      nmvf = size(fv,1); nfvd = size(fv,2); nfed = size(fe,2)
      if (nrec < 1) return
      call dtimer(dtime,itime,-1)
      write (unit=iunit,rec=nrec) ((fvm(j,k),j=1,ndim),k=1,3),          &
     &((fv(j,k),j=1,nmvf),k=1,nfvd), ((fe(j,k),j=1,nmvf),k=1,nfed), wk
      nrec = nrec + 1
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwritetr2(partt,tdiag,iunit,nrec)
! this subroutine writes real trajetory record in direct access binary
! file
! partt = real trajectory data to be written in each record 
! iunit = fortran unit number to be used 
! nrec = record number for write (then updated to next record)
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: partt
! local data
      integer :: j, k, idimp, nprobt
      integer, dimension(4) :: itime
      double precision :: dtime
      idimp = size(partt,1); nprobt = size(partt,2)
      if (nrec < 1) return
      call dtimer(dtime,itime,-1)
      write (unit=iunit,rec=nrec) ((partt(j,k),j=1,idimp),k=1,nprobt)
      nrec = nrec + 1
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
      end module
