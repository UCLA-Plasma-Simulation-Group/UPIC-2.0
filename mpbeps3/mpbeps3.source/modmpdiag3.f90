!-----------------------------------------------------------------------
!
      module mpdiag3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpdiag3.f
! get_funit returns an unconnected fortran unit number
! fnrecl find record length of direct access file
! dafopen3 opens new binary file for real 1d scalar data.
! dafopenv3 opens new binary file for real 1d vector data.
! dafopenc3 opens new binary file for complex 1d scalar data.
! dafopenvc3 opens new binary file for complex 1d vector data.
! dafopenfv3 opens new binary file for real velocity-space data.
! dafopentr3 opens new binary file for real trajectory data.
! mpvpdist3 calculates 3 component velocity distribution, velocity
!           moments, and entropy with segmented particle array
!           PPVDIST32
! mperpdist3 calculates 3d energy distribution for relativistic
!            particles with segmented particle array
!            calls PERPDIST32
! mpvbpdist3 calculates 3d velocity distribution and velocity moments
!            for magnetized plasma with segmented particle array
!            calls PPVBDIST32
! mpvspdist3 calculates 3d velocity distribution in different regions of
!            space
!            calls PPVSDIST32
! mpvdist3 calculates 3 component velocity distribution, velocity
!          moments, and entropy with standard particle array
!          calls PVDIST32
! mvbdist3 calculates 3d velocity distribution and velocity moments
!          for magnetized plasma with standard particle array
! mpprofx3 calculates fluid moments from particle quantities
!          assumes particle positions and velocities at same time level
!          calls PPROFX32L
! mprprofx3 calculates fluid moments from relativistic particle
!           quantities.  assumes particle positions and momenta at same
!           time level
!           calls PRPROFX32L
! mpgprofx3 calculates fluid moments from particle quantities
!           and electrostatic fields assumes particle positions and
!           velocities not at same time levels
!           calls PGPROFX32L
! mpgrprofx3 calculates fluid moments from relativistic particle
!            quantities and electrostatic fields.  assumes particle
!            positions and velocities not at same time levels
!            calls PGRPROFX32L
! mpgbprofx3 calculates fluid moments from particle quantities and 
!            electromagnetic fields assumes particle positions and
!            velocities not at same time levels
!            calls PGBPROFX32L
! mpgrbprofx3 calculates fluid moments from relativistic particle
!             quantities. assumes particle positions and velocities
!             not at same time levels
!             calls PGRBPROFX32L
! mpfluidqs3 calculates fluid quantities from fluid moments
!            calls FLUIDQS3
! wmprofx3 generic procedure to calculate fluid moments
!          calls mpprofx3 or mprprofx3
! wmgprofx3 generic procedure to calculate fluid moments with
!           electrostatic fields
!           calls mpgprofx3 or mpgrprofx3
! wmgbprofx3 generic procedure to calculate fluid moments with
!            electromagnetic fields
!            calls mpgbprofx3 or mpgrbprofx3
! psettraj3 sets test charge distribution by setting a particle id in
!           particle location 7
!           calls PSTPTRAJ3
! mptraj3 copies tagged particles in ppart to array partt
!         calls PPTRAJ3
! mpordtraj3 reorders tagged particles in partt to array spartt
!            calls PORDTRA23
! mpcpytraj3 copies tagged particles in partt to array part
!            calls PCPYTRAJ32
! dafwritefv3 writes velocity-space record in direct access binary file
! dafwritetr3  writes trajectory record in direct access binary file
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: march 21, 2018
!
      use mppmod3, only: mpsum2
      use libmpdiag3_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      function get_funit(start) result(funit)
! this function returns an unconnected fortran unit number,
! starting with unit = start.  returns -1 if none found
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
      character(len=*), intent(in) :: fname
      integer :: it, ios
      inquire(file=fname,recl=it,iostat=ios)
      if (ios /= 0) it = 0
      end function
!
!-----------------------------------------------------------------------
      subroutine dafopen3(f,nx,kyp,kzp,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 3d real scalar data array to be written in each record
! nx/kyp/kzp = number of data elements per record to be written in x/y/z
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx, kyp, kzp
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1,1); lrec = lrec*nx*kyp*kzp
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenv3(f,nx,kyp,kzp,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 3d real vector data array to be written in each record
! nx/kyp/kzp = number of data elements per record to be written in x/y/z
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx, kyp, kzp
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:,:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1,1,1)
      lrec = lrec*size(f,1)*nx*kyp*kzp
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenc3(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = 3d complex scalar data array to be written in each record
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
      subroutine dafopenvc3(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = 3d complex vector data array to be written in each record
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      complex, dimension(:,:,:,:), intent(in) :: fc
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fc(1,1,1,1); lrec = lrec*size(fc)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenfv3(fvm,fv,fe,wk,iunit,nrec,fname)
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
      subroutine dafopentr3(partt,iunit,nrec,fname)
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
      subroutine mpvpdist3(ppart,kpic,fv,sfv,fvm,tdiag,nvp,nmv)
! calculates 3d velocity distribution, velocity moments, and entropy
! with segmented particle array
      integer, intent(in) :: nvp, nmv
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:,:), intent(inout) :: fv, fvm
      real, dimension(:,:,:), intent(inout) :: sfv
! local data
      integer :: idimp, nppmx, nmvf, idimv, mxyzp1
      integer, dimension(4) :: itime
      real, dimension(3) :: scale
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(fv,1); idimv = size(fv,2)
      mxyzp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PPVDIST32(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,mxyzp1,nmv,  &
     &nmvf)
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
      subroutine mperpdist3(ppart,kpic,fv,sfv,ci,wk,tdiag,nmv)
! calculates 3d energy distribution for relativistic particles
! with segmented particle array
      integer, intent(in) :: nmv
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
      call PERPDIST32(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyp1,nmv,nmvf&
     &)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
! sum distribution over processors, but save energy scale
      scale(1:idimv) = fv(2*nmv+2,:)
      call mpsum2(fv,tdiag)
      fv(2*nmv+2,:) = scale(1:idimv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvbpdist3(ppart,kpic,fv,sfv,fvm,omx,omy,omz,tdiag,nmv)
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
      integer :: idimp, nppmx, nmvf,  mxyzp1
      integer, dimension(4) :: itime
      real, dimension(2) :: scale
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(fv,1)
      mxyzp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PPVBDIST32(ppart,kpic,fv,sfv,fvm,omx,omy,omz,idimp,nppmx,    &
     &mxyzp1,nmv,nmvf)
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
      subroutine mpvspdist3(ppart,kpic,fvs,noff,tdiag,nmv,mvx,mvy,mvz,  &
     &nyb)
! calculates 3d velocity distribution in different regions of space
      implicit none
      integer, intent(in) :: nmv, mvx, mvy, mvz, nyb
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:,:), intent(inout) :: fvs
      integer, dimension(:), intent(in) :: kpic, noff
! local data
      integer :: idimp, nppmx, nmvf, idimv, nxb, nybmx, nzb, mxyzp1
      integer :: idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(fvs,1); idimv = size(fvs,2)
      nxb = size(fvs,3); nybmx = size(fvs,4) - 1; nzb = size(fvs,5) - 1
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimv==3) then
         call PPVSDIST32(ppart,kpic,fvs,noff,nmv,mvx,mvy,mvz,nxb,nyb,nzb&
     &,idimp,nppmx,mxyzp1,nmvf,nybmx,idds)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpvdist3(part,fv,fvm,tdiag,npp,nvp,nmv)
! calculates 3d velocity distribution, velocity moments, and entropy
! with standard particle array
      integer, intent(in) :: npp, nvp, nmv
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:), intent(inout) :: fv, fvm
! local data
      integer :: idimp, npmax, nmvf
      integer, dimension(4) :: itime
      real, dimension(3) :: scale
      double precision :: dtime
! extract dimensions
      idimp = size(part,1)
      nmvf = size(fv,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PVDIST32(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
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
      subroutine mvbdist3(part,fv,fvm,omx,omy,omz,tdiag,npp,nmv)
! calculates 3d velocity distribution and velocity moments
! for magnetized plasma with standard particle array
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
      call PVBDIST32(part,fv,fvm,omx,omy,omz,npp,idimp,npmax,nmv,nmvf)
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
      subroutine mpprofx3(ppart,fms,kpic,noff,tdiag,npro,mx,my,mz,mx1,  &
     &myp1)
! calculates fluid moments from particle quantities
! assumes particle positions and velocities at same time level
      implicit none
      integer, intent(in) :: npro, mx, my, mz, mx1, myp1
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      nypmx = size(fms,3); nzpmx = size(fms,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PPROFX32L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,my,mz,nprd,&
     &nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprprofx3(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mz,  &
     &mx1,myp1)
! calculates fluid moments from relativistic particle quantities
! assumes particle positions and momenta at same time level
      integer, intent(in) :: npro, mx, my, mz, mx1, myp1
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      nypmx = size(fms,3); nzpmx = size(fms,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PRPROFX32L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro,mx,my,mz, &
     &nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgprofx3(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,tdiag,  &
     &npro,nx,mx,my,mz,mx1,myp1)
! calculates fluid moments from particle quantities
! and electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx, my, mz, mx1, myp1
      real, intent(in) :: qbm, dt
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      nypmx = size(fms,3); nzpmx = size(fms,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PGPROFX32L(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,idimp,nppmx, &
     &npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrprofx3(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,ci,    &
     &tdiag,npro,nx,mx,my,mz,mx1,myp1)
! calculates fluid moments from relativistic particle quantities
! and electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx, my, mz, mx1, myp1
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      nypmx = size(fms,3); nzpmx = size(fms,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PGRPROFX32L(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,ci,idimp,   &
     &nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgbprofx3(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt,  &
     &tdiag,npro,nx,mx,my,mz,mx1,myp1)
! calculates fluid moments from particle quantities
! and electromagnetic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx, my, mz, mx1, myp1
      real, intent(in) :: qbm, dt
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      nypmx = size(fms,3); nzpmx = size(fms,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PGBPROFX32L(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt,idimp, &
     &nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgrbprofx3(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt, &
     &ci,tdiag,npro,nx,mx,my,mz,mx1,myp1)
! calculates fluid moments from relativistic particle quantities
! and electromagnetic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx, my, mz, mx1, myp1
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! local data
      integer :: idimp, nppmx, nprd, nxv, nypmx, nzpmx, mxyzp1, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      nypmx = size(fms,3); nzpmx = size(fms,4)
      mxyzp1 = size(kpic,1); idds = size(noff,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PGRBPROFX32L(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt,ci,   &
     &idimp,nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &idds)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpfluidqs3(fms,tdiag,npro,nx,ny,nz,kstrt,nvpy,kyp,kzp)
! calculates fluid quantities from fluid moments
      implicit none
      integer, intent(in) :: npro, nx, ny, nz, kstrt, nvpy, kyp, kzp
      real, intent(inout) :: tdiag
      real, dimension(:,:,:,:), intent(inout) :: fms
! local data
      integer :: nprd, nxv, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nprd = size(fms,1); nxv = size(fms,2)
      nypmx = size(fms,3); nzpmx = size(fms,4)
! initialize timer
      call dtimer(dtime,itime,-1)
      call FLUIDQS3(fms,npro,nx,ny,nz,kstrt,nvpy,kyp,kzp,nprd,nxv,nypmx,&
     &nzpmx)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmprofx3(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mz,mx1&
     &,myp1,relativity)
! generic procedure to calculate fluid moments
! assumes particle positions and velocities at same time levels
      integer, intent(in) :: npro, mx, my, mz, mx1, myp1
      integer, intent(in) :: relativity
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff
! updates fms
      if (relativity==1) then
         call mprprofx3(ppart,fms,kpic,noff,ci,tdiag,npro,mx,my,mz,mx1, &
     &myp1)
      else
         call mpprofx3(ppart,fms,kpic,noff,tdiag,npro,mx,my,mz,mx1,myp1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgprofx3(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,ci,tdiag&
     &,npro,nx,mx,my,mz,mx1,myp1,relativity)
! generic procedure to calculate fluid moments with electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx, my, mz, mx1, myp1
      integer, intent(in) :: relativity
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! updates fms
      if (relativity==1) then
         call mpgrprofx3(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,ci,tdiag, &
     &npro,nx,mx,my,mz,mx1,myp1)
      else
         call mpgprofx3(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,tdiag,npro,&
     &nx,mx,my,mz,mx1,myp1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgbprofx3(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt,ci&
     &,tdiag,npro,nx,mx,my,mz,mx1,myp1,relativity)
! generic procedure to calculate fluid moments with electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx, my, mz, mx1, myp1
      integer, intent(in) :: relativity
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:,:), intent(in) :: fxyz, bxyz
      real, dimension(:,:,:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(in) :: noff, nyzp
! updates fms
      if (relativity==1) then
         call mpgrbprofx3(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt,ci, &
     &tdiag,npro,nx,mx,my,mz,mx1,myp1)
      else
         call mpgbprofx3(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt,tdiag&
     &,npro,nx,mx,my,mz,mx1,myp1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine psetptraj3(ppart,tedges,kpic,iprobt,kstrt,nst,nvpy,nvpz&
     &,vtx,vtsx,dvtx,np,nprobt,irc)
! sets test charge distribution by setting a particle id in particle
! location 7
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! nprobt = number of test charges whose trajectories will be stored.
! irc = (0,1) = (no,yes) error condition exists
      integer, intent(in) :: kstrt, nst, nvpy, nvpz
      integer, intent(inout) :: nprobt, irc
      real, intent(in) :: vtx, vtsx, dvtx
      double precision, intent(in) :: np
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:), intent(inout) :: tedges
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(inout) :: iprobt
! local data
      integer :: idimp, nppmx, mxyzp1, idps
      irc = 0
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1); idps = size(tedges,1)
! call low level procedure
      if (idimp > 6) then
         call PSTPTRAJ3(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,vtsx,   &
     &dvtx,nvpy,nvpz,idimp,nppmx,mxyzp1,idps,np,nprobt)
      else
         if (kstrt==1) then
            write (*,*) 'psetptraj3 error: idimp=', idimp
            irc = 1
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mptraj3(ppart,kpic,partt,tdiag,numtp,irc)
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
      integer :: idimp, nppmx, mxyzp1, nprobt
      integer, dimension(4) :: itime
      double precision :: dtime
      irc = 0
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyzp1 = size(kpic,1); nprobt = size(partt,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimp > 6) then
         call PPTRAJ3(ppart,kpic,partt,numtp,idimp,nppmx,mxyzp1,nprobt)
         if (numtp > nprobt) irc = -1
      else
         write (*,*) 'mptraj3 error: idimp=', idimp
         irc = 1
      endif
      if (irc < 0) then
         write (*,*) 'info:mptraj3 overflow: numpt,nprobt=',numtp,nprobt
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpordtraj3(partt,spartt,tedges,tdiag,numtp,irc)
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
      if (idimp > 6) then
         call PORDTRAJ3(partt,spartt,tedges,numtp,idimp,idps,nprobt)
      else
         write (*,*) 'mpordtraj3 error: idimp=', idimp
         irc = 1
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcpytraj3(partt,part,tdiag,numtp)
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
      call PCPYTRAJ3(partt,part,numtp,idimp,nprobt)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwritefv3(fvm,fv,fe,wk,tdiag,iunit,nrec)
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
      subroutine dafwritetr3(partt,tdiag,iunit,nrec)
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
