!-----------------------------------------------------------------------
!
      module mdiag1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmdiag1.f
! get_funit returns an unconnected fortran unit number
! fnrecl find record length of direct access file
! dafopen1 opens new binary file for real 1d scalar data.
! dafopenv1 opens new binary file for real 1d vector data.
! dafopenc1 opens new binary file for complex 1d scalar data.
! dafopenvc1 opens new binary file for complex 1d vector data.
! dafwrite1 writes real scalar record in direct access binary file
! dafopenfv1 opens new binary file for real velocity-space data.
! dafopentr1 opens new binary file for real trajectory data.
! fnopens1: opens a new fortran unformatted stream file
! fnsets1: repositions a fortran unformatted stream file
! dafwritev1 writes real vector record in direct access binary file
! dafwritec1 writes complex scalar record in direct access binary file
! dafwritevc1 writes complex vector record in direct access binary file
! mcspect1 performs frequency analysis of complex time series
!          calls CSPECT1
! micspect1 performs incremental frequency analysis of complex scalar
!           time series for one time step
!           calls ICSPECT1
! mivcspect1 performs incremental frequency analysis of complex vector 
!            time series for one time step
!            calls IVCSPECT1
! mvpdist1 calculates 1 component velocity distribution, velocity
!          moments, and entropy with segmented particle array
!          calls VPDIST1 or VPDIST13
! mvbpdist1 calculates 3d velocity distribution and velocity moments
!           for magnetized plasma with segmented particle array
!           calls VBPDIST13
! merpdist1 calculates 1d or 1-2/2d energy distribution for relativistic
!           particles with segmented particle array
!           calls ERPDIST1 or ERPDIST13
! mvspdist1 calculates 1d velocity distribution in different regions of
!           space
!           calls PVSDIST1 or PVSDIST13
! mvdist1 calculates 1 or 3 component velocity distribution, velocity
!         moments, and entropy with standard particle array
!         calls VPDIST1 or VDIST13
! mvbdist1 calculates 3d velocity distribution and velocity moments
!          for magnetized plasma with standard particle array
!          calls VBDIST13
! mprofx13 calculates fluid moments from particle quantities
!          assumes particle positions and velocities at same time level
!          for 1-2/2d code
!          calls PROFX13L
! mrprofx13 calculates fluid moments from relativistic particle
!           quantities.  assumes particle positions and momenta at same
!           time level, for 1-2/2d code
!           calls RPROFX13L
! mprofx1 calculates fluid moments from particle quantities
!         assumes particle positions and velocities at same time level
!         for 1d code
!         calls PROFX1L
! mrprofx1 calculates fluid moments from relativistic particle
!          quantities.  assumes particle positions and momenta at same
!          time level, for 1d code
!          calls RPROFX1L
! mgprofx1 calculates fluid moments from particle quantities
!          and electrostatic fields assumes particle positions and
!          velocities not at same time levels
!          calls GPROFX1L
! mgrprofx1 calculates fluid moments from relativistic particle
!           quantities and electrostatic fields.  assumes particle
!           positions and velocities not at same time levels
!           calls GRPROFX1L
! mgbprofx1 calculates fluid moments from particle quantities and 
!           electromagnetic fields assumes particle positions and
!           velocities not at same time levels
!           calls GBPROFX13L
! mgrbprofx1 calculates fluid moments from relativistic particle
!            quantities. assumes particle positions and velocities
!            not at same time levels
!            calls GRBPROFX13L
! mfluidqs13 calculates fluid quantities from fluid moments
!             for 1-2/2d code
!             calls FLUIDQS13
! mfluidqs1 calculates fluid quantities from fluid moments for 1d code
!             calls FLUIDQS1
! wmprofx13 generic procedure to calculate fluid moments for 1-2/2d code
!           calls mprofx13 or mrprofx13
! wmprofx1 generic procedure to calculate fluid moments
!          calls mprofx1 or mrprofx1 for 1d code
! wmgprofx1 generic procedure to calculate fluid moments with
!           electrostatic fields
!           calls mgprofx1 or mgrprofx1
! wmgbprofx1 generic procedure to calculate fluid moments with
!            electromagnetic fields
!            calls mgbprofx1 or mgrbprofx1
! settraj1 sets test charge distribution by setting a particle id in
!          particle location 3 or 5
!          calls STPTRAJ1 or STPTRAJ13
! mfnptraj1 finds tagged particles in ppart
!           calls FNPTRAJ1 or FNPTRAJ13
! mptraj1 copies tagged particles in ppart to array partt
!         calls PTRAJ1 or PTRAJ13
! dafwritefv1 writes velocity-space record in direct access binary file
! dafwritetr1 writes trajectory record in direct access binary file
! setmbeam1 marks beam particles by setting a particle id in particle
!           location 3 or 5
!           calls STPBEAM1 or STPBEAM13
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 1, 2021
!
      use libmdiag1_h
      use in1, only: nustrt
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
      subroutine dafopen1(f,nx,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 1d real scalar data array to be written in each record
! nx = number of data elements per record to be written in x 
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx
      integer, intent(inout) :: iunit, nrec
      real, dimension(:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1); lrec = lrec*nx
      iunit = get_funit(iunit)
      if (nustrt==2) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old')
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace')
      endif
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenv1(f,nx,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! f = 1d real vector data array to be written in each record
! nx = number of data elements per record to be written in x 
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(in) :: nx
      integer, intent(inout) :: iunit, nrec
      real, dimension(:,:), intent(in) :: f
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1,1); lrec = lrec*size(f,1)*nx
      iunit = get_funit(iunit)
      if (nustrt==2) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old')
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace')
      endif
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenc1(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = scalar data array to be written in each record
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      complex, dimension(:), intent(in) :: fc
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fc(1); lrec = lrec*size(fc)
      iunit = get_funit(iunit)
      if (nustrt==2) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old')
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace')
      endif
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenvc1(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = vector data array to be written in each record
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
      if (nustrt==2) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old')
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace')
      endif
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopenfv1(fvm,fv,fe,wk,iunit,nrec,fname)
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
      if (nustrt==2) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old')
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace')
      endif
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafopentr1(partt,iunit,nrec,fname)
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
      if (nustrt==2) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old')
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace')
      endif
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fnopens1(iunit,fname)
! this subroutine opens a new fortran unformatted stream file
! this is a Fortran2003 feature
! iunit = fortran unit number to be used 
! fname = file name
      implicit none
      integer, intent(in) :: iunit
      character(len=*), intent(in) :: fname
      if (nustrt==2) then
         open(unit=iunit,file=fname,access='stream',form='unformatted', &
     &status='old')
      else
         open(unit=iunit,file=fname,access='stream',form='unformatted', &
     &status='replace')
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fnsets1(nrec,nsize,fname)
! this subroutine repositions a fortran unformatted stream file
! this is a Fortran2003 feature
! nrec = number of writes
! nsize = number of words in each write
! fname = file name
      implicit none
      integer, intent(in) :: nrec, nsize
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec, iunit
      integer(kind=selected_int_kind(8)) :: loc
      real :: f
      logical :: opn
      inquire(iolength=lrec) f
      loc = lrec*nrec*nsize + 1
      inquire (file=fname,opened=opn,number=iunit)
      if (opn) then
         write (unit=iunit,pos=loc)
      else
         write (*,*) 'file not repositioned (not open):', fname
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwrite1(f,tdiag,iunit,nrec,nx)
! this subroutine writes real scalar record in direct access binary file
! f = real scalar data array to be written
! iunit = fortran unit number to be used 
! nrec = record number for write (then updated to next record)
! nx = number of elements to be written in record
      implicit none
      integer, intent(in) :: iunit, nx
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      real, dimension(:), intent(in) :: f
! local data
      integer :: j
      integer, dimension(4) :: itime
      double precision :: dtime
      if (nrec < 1) return
      call dtimer(dtime,itime,-1)
      write (unit=iunit,rec=nrec) (f(j),j=1,nx)
      nrec = nrec + 1
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwritev1(f,tdiag,iunit,nrec,nx)
! this subroutine writes real vector record in direct access binary file
! f = real vector data array to be written
! iunit = fortran unit number to be used 
! nrec = record number for write (then updated to next record)
! nx = number of elements to be written in record
      implicit none
      integer, intent(in) :: iunit, nx
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: f
! local data
      integer :: j, k, ndim
      integer, dimension(4) :: itime
      double precision :: dtime
      ndim = size(f,1)
      if (nrec < 1) return
      call dtimer(dtime,itime,-1)
      write (unit=iunit,rec=nrec) ((f(j,k),j=1,ndim),k=1,nx)
      nrec = nrec + 1
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwritec1(fc,tdiag,iunit,nrec,nx)
! this subroutine writes scalar record in direct access binary file
! fc = complex scalar data array to be written
! iunit = fortran unit number to be used 
! nrec = record number for write (then updated to next record)
! nx = number of elements to be written in record
      implicit none
      integer, intent(in) :: iunit, nx
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      complex, dimension(:), intent(in) :: fc
! local data
      integer :: j
      integer, dimension(4) :: itime
      double precision :: dtime
      if (nrec < 1) return
      call dtimer(dtime,itime,-1)
      write (unit=iunit,rec=nrec) (fc(j),j=1,nx)
      nrec = nrec + 1
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwritevc1(fc,tdiag,iunit,nrec,nx)
! this subroutine writes vector record in direct access binary file
! fc = complex vector data array to be written
! iunit = fortran unit number to be used 
! nrec = record number for write (then updated to next record)
! nx = number of elements to be written in record
      implicit none
      integer, intent(in) :: iunit, nx
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      complex, dimension(:,:), intent(in) :: fc
! local data
      integer :: j, k, ndim
      integer, dimension(4) :: itime
      double precision :: dtime
      ndim = size(fc,1)
      if (nrec < 1) return
      call dtimer(dtime,itime,-1)
      write (unit=iunit,rec=nrec) ((fc(j,k),j=1,ndim),k=1,nx)
      nrec = nrec + 1
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcspect1(fc,wm,pkw,t0,dt,tdiag,nt,iw,modesx)
! performs incremental frequency analysis of complex time series
      integer, intent(in) :: nt, iw, modesx
      real, intent(in) :: t0, dt
      real, intent(inout) :: tdiag
      complex, dimension(:,:), intent(in) :: fc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:), intent(inout) :: pkw
! local data
      integer :: ntd, iwd, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ntd = size(fc,1); modesxd = size(fc,2)
      iwd = size(wm,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call CSPECT1(fc,wm,pkw,t0,dt,nt,iw,modesx,ntd,iwd,modesxd)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine micspect1(fc,wm,pkw,pks,time,t0,tdiag,nt,iw,modesx,nx, &
     &norm)
! performs incremental frequency analysis of complex scalar time series
! for one time step
! norm = (-1,0,1) = normalize with (inverse gradient,null,gradient) op
      integer, intent(in) :: nt, iw, modesx, nx, norm
      real, intent(in) :: time, t0
      real, intent(inout) :: tdiag
      complex, dimension(:), intent(in) :: fc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:), intent(inout) :: pkw
      double precision, dimension(:,:,:), intent(inout) :: pks
! local data
      integer :: iwd, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      modesxd = size(fc,1); iwd = size(wm,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ICSPECT1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm,iwd,     &
     &modesxd)
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mivcspect1(fvc,wm,vpkw,vpks,time,t0,tdiag,nt,iw,modesx,&
     &nx,norm)
! performs incremental frequency analysis of complex vector time series
! for one time step
! norm = (-1,0,1) = normalize with (inverse curl,null,curl) op
      integer, intent(in) :: nt, iw, modesx, nx, norm
      real, intent(in) :: time, t0
      real, intent(inout) :: tdiag
      complex, dimension(:,:), intent(in) :: fvc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:,:), intent(inout) :: vpkw
      double precision, dimension(:,:,:,:), intent(inout) :: vpks
! local data
      integer :: iwd, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      modesxd = size(fvc,2); iwd = size(wm,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call IVCSPECT1(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,nx,norm,iwd, &
     &modesxd)
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvpdist1(ppart,kpic,sfv,fvm,tdiag,np,nmv)
! calculates 1d velocity distribution, velocity moments, and entropy
! with segmented particle array
      integer, intent(in) :: np, nmv
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:,:,:), intent(inout) :: sfv
      real, dimension(:,:), intent(inout) :: fvm
! local data
      integer :: idimp, nppmx, nmvf, idimv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(sfv,1); idimv = size(sfv,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      if (idimv==1) then
         call VPDIST1(ppart,kpic,sfv,fvm,idimp,nppmx,mx1,np,nmv,nmvf)
      else if (idimv==3) then
         call VPDIST13(ppart,kpic,sfv,fvm,idimp,nppmx,mx1,np,nmv,nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvbpdist1(ppart,kpic,sfv,fvm,omx,omy,omz,tdiag,np,nmv)
! calculates 3d velocity distribution and velocity moments
! for magnetized plasma with segmented particle array
      integer, intent(in) :: np, nmv
      real, intent(in) :: omx, omy, omz
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:,:,:), intent(inout) :: sfv
      real, dimension(:,:), intent(inout) :: fvm
! local data
      integer :: idimp, nppmx, nmvf, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(sfv,1); mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call VBPDIST13(ppart,kpic,sfv,fvm,omx,omy,omz,idimp,nppmx,mx1,np, &
    &nmv,nmvf)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine merpdist1(ppart,kpic,sfv,ci,wk,tdiag,nmv)
! calculates 1d or 1-2/2d energy distribution for relativistic particles
! with segmented particle array
      integer, intent(in) :: nmv
      real, intent(in) :: ci
      real, intent(inout) :: wk, tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:,:,:), intent(inout) :: sfv
! local data
      integer :: idimp, nppmx, nmvf, idimv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(sfv,1); idimv = size(sfv,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      if (idimv==1) then
         call ERPDIST1(ppart,kpic,sfv,ci,wk,idimp,nppmx,mx1,nmv,nmvf)
      else if (idimv==3) then
         call ERPDIST13(ppart,kpic,sfv,ci,wk,idimp,nppmx,mx1,nmv,nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvspdist1(ppart,kpic,fvs,tdiag,nmv,mvx)
! calculates 1d velocity distribution in different regions of space
      implicit none
      integer, intent(in) :: nmv, mvx
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: fvs
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nmvf, idimv, nxb, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(fvs,1); idimv = size(fvs,2)
      nxb = size(fvs,3)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimv==1) then
         call PVSDIST1(ppart,kpic,fvs,nmv,mvx,nxb,idimp,nppmx,mx1,nmvf)
      else if (idimv==3) then
         call PVSDIST13(ppart,kpic,fvs,nmv,mvx,nxb,idimp,nppmx,mx1,nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvdist1(part,fv,fvm,tdiag,np,nmv)
! calculates 1d velocity distribution, velocity moments, and entropy
! with standard particle array
      integer, intent(in) :: np, nmv
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:), intent(inout) :: fv, fvm
! local data
      integer :: idimp, nmvf, idimv
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(part,1)
      nmvf = size(fv,1); idimv = size(fv,2)
! initialize timer
      call dtimer(dtime,itime,-1)
      if (idimv==1) then
         call VDIST1(part,fv,fvm,idimp,np,nmv,nmvf)
      else if (idimv==3) then
         call VDIST13(part,fv,fvm,idimp,np,nmv,nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvbdist1(part,fv,fvm,omx,omy,omz,tdiag,np,nmv)
! calculates 3d velocity distribution and velocity moments
! for magnetized plasma with standard particle array
      integer, intent(in) :: np, nmv
      real, intent(in) :: omx, omy, omz
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:), intent(inout) :: fv, fvm
! local data
      integer :: idimp, nmvf
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(part,1)
      nmvf = size(fv,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call VBDIST13(part,fv,fvm,omx,omy,omz,idimp,np,nmv,nmvf)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprofx13(ppart,fms,kpic,tdiag,npro,mx)
! calculates fluid moments from particle quantities
! assumes particle positions and velocities at same time level
! for 1-2/2d code
      implicit none
      integer, intent(in) :: npro, mx
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PROFX13L(ppart,fms,kpic,nppmx,idimp,npro,mx,nprd,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrprofx13(ppart,fms,kpic,ci,tdiag,npro,mx)
! calculates fluid moments from relativistic particle quantities
! assumes particle positions and momenta at same time level
! for 1-2/2d code
      integer, intent(in) :: npro, mx
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call RPROFX13L(ppart,fms,kpic,ci,nppmx,idimp,npro,mx,nprd,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprofx1(ppart,fms,kpic,tdiag,npro,mx)
! calculates fluid moments from particle quantities
! assumes particle positions and velocities at same time level
! for 1d code
      implicit none
      integer, intent(in) :: npro, mx
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call PROFX1L(ppart,fms,kpic,nppmx,idimp,npro,mx,nprd,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrprofx1(ppart,fms,kpic,ci,tdiag,npro,mx)
! calculates fluid moments from relativistic particle quantities
! assumes particle positions and momenta at same time level
! for 2d code
      integer, intent(in) :: npro, mx
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call RPROFX1L(ppart,fms,kpic,ci,nppmx,idimp,npro,mx,nprd,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgprofx1(ppart,fx,fms,kpic,qbm,dt,tdiag,npro,nx,mx)
! calculates fluid moments from particle quantities
! and electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx
      real, intent(in) :: qbm, dt
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:), intent(in) :: fx
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call GPROFX1L(ppart,fx,fms,kpic,qbm,dt,idimp,nppmx,npro,nx,mx,nprd&
     &,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgrprofx1(ppart,fx,fms,kpic,qbm,dt,ci,tdiag,npro,nx,mx)
! calculates fluid moments from relativistic particle quantities
! and electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:), intent(in) :: fx
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call GRPROFX1L(ppart,fx,fms,kpic,qbm,dt,ci,idimp,nppmx,npro,nx,mx,&
     &nprd,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgbprofx1(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,tdiag,npro&
     &,nx,mx)
! calculates fluid moments from particle quantities
! and electromagnetic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx
      real, intent(in) :: omx, qbm, dt
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call GBPROFX13L(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,idimp,nppmx,   &
     &npro,nx,mx,nprd,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgrbprofx1(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,ci,tdiag,&
     &npro,nx,mx)
! calculates fluid moments from relativistic particle quantities
! and electromagnetic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx
      real, intent(in) :: omx, qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, nprd, nxv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nprd = size(fms,1); nxv = size(fms,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      call GRBPROFX13L(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,ci,idimp,nppmx&
     &,npro,nx,mx,nprd,nxv,mx1)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfluidqs13(fms,tdiag,npro,nx)
! calculates fluid quantities from fluid moments for 1-2/2d code
      implicit none
      integer, intent(in) :: npro, nx
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(inout) :: fms
! local data
      integer :: nprd, nxv
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nprd = size(fms,1); nxv = size(fms,2)
! initialize timer
      call dtimer(dtime,itime,-1)
      call FLUIDQS13(fms,npro,nx,nprd,nxv)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfluidqs1(fms,tdiag,npro,nx)
! calculates fluid quantities from fluid moments for 1d code
      implicit none
      integer, intent(in) :: npro, nx
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(inout) :: fms
! local data
      integer :: nprd, nxv
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nprd = size(fms,1); nxv = size(fms,2)
! initialize timer
      call dtimer(dtime,itime,-1)
      call FLUIDQS1(fms,npro,nx,nprd,nxv)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmprofx13(ppart,fms,kpic,ci,tdiag,npro,mx,relativity)
! generic procedure to calculate fluid moments
! assumes particle positions and velocities at same time levels
! for 1-2/2d code
      integer, intent(in) :: npro, mx, relativity
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! updates fms
      if (relativity==1) then
         call mrprofx13(ppart,fms,kpic,ci,tdiag,npro,mx)
      else
         call mprofx13(ppart,fms,kpic,tdiag,npro,mx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmprofx1(ppart,fms,kpic,ci,tdiag,npro,mx,relativity)
! generic procedure to calculate fluid moments
! assumes particle positions and velocities at same time levels
! for 1d code
      integer, intent(in) :: npro, mx, relativity
      real, intent(in) :: ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! updates fms
      if (relativity==1) then
         call mrprofx1(ppart,fms,kpic,ci,tdiag,npro,mx)
      else
         call mprofx1(ppart,fms,kpic,tdiag,npro,mx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgprofx1(ppart,fx,fms,kpic,qbm,dt,ci,tdiag,npro,nx,mx,&
     &relativity)
! generic procedure to calculate fluid moments with electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx, relativity
      real, intent(in) :: qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:), intent(in) :: fx
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! updates fms
      if (relativity==1) then
         call mgrprofx1(ppart,fx,fms,kpic,qbm,dt,ci,tdiag,npro,nx,mx)
      else
         call mgprofx1(ppart,fx,fms,kpic,qbm,dt,tdiag,npro,nx,mx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmgbprofx1(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,ci,tdiag,&
     &npro,nx,mx,relativity)
! generic procedure to calculate fluid moments with electrostatic fields
! assumes particle positions and velocities not at same time levels
      implicit none
      integer, intent(in) :: npro, nx, mx, relativity
      real, intent(in) :: omx, qbm, dt, ci
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(in) :: fxyz, byz
      real, dimension(:,:), intent(inout) :: fms
      integer, dimension(:), intent(in) :: kpic
! updates fms
      if (relativity==1) then
         call mgrbprofx1(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,ci,tdiag,   &
     &npro,nx,mx)
      else
         call mgbprofx1(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,tdiag,npro,nx&
     &,mx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine setptraj1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,np,nprobt&
     &,irc)
! sets test charge distribution by setting a particle id in particle
! location 3 or 5
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! nprobt = number of test charges whose trajectories will be stored.
! irc = (0,1) = (no,yes) error condition exists
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(inout) :: iprobt
      integer, intent(in) :: nst, np
      integer, intent(inout) :: nprobt, irc
      real, intent(in) :: vtx, vtsx, dvtx
! local data
      integer :: idimp, mx1, nppmx
      irc = 0
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! call low level procedure
      if (idimp > 4) then
         call STPTRAJ13(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,nppmx,&
     &mx1,np,nprobt)
      else if (idimp > 2) then
         call STPTRAJ1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,nppmx, &
     &mx1,np,nprobt)
      else
         write (*,*) 'setptraj1 error: idimp=', idimp
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfnptraj1(ppart,kpic,nprobt,irc)
! this finds tagged particles in ppart
! nprobt = number of test charges whose trajectories will be stored.
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(inout) :: nprobt, irc
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mx1
      irc = 0
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! call low level procedure
      if (idimp > 4) then
         call FNPTRAJ13(ppart,kpic,idimp,nppmx,mx1,nprobt)
      else if (idimp > 2) then
         call FNPTRAJ1(ppart,kpic,idimp,nppmx,mx1,nprobt)
      else
         write (*,*) 'mfnptraj1 error: idimp=', idimp
         irc = 1
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mptraj1(ppart,kpic,partt,tdiag,irc)
! this copies tagged particles in ppart to array partt
! irc = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(inout) :: irc
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:,:), intent(inout) :: partt
! local data
      integer :: idimp, nppmx, mx1, nprobt
      integer, dimension(4) :: itime
      double precision :: dtime
      irc = 0
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1); nprobt = size(partt,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (idimp > 4) then
         call PTRAJ13(ppart,kpic,partt,idimp,nppmx,mx1,nprobt)
      else if (idimp > 2) then
         call PTRAJ1(ppart,kpic,partt,idimp,nppmx,mx1,nprobt)
      else
         write (*,*) 'mptraj1 error: idimp=', idimp
         irc = 1
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwritefv1(fvm,fv,fe,wk,tdiag,iunit,nrec)
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
      subroutine dafwritetr1(partt,tdiag,iunit,nrec)
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
!-----------------------------------------------------------------------
      subroutine mwrfvsdata1(fvs,tdiag,iunit)
! writes 3d real vector data to a fortran unformatted file
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: fvs
! local data
      integer :: i, j, k, nmvf, ndim, nsxb
      integer, dimension(4) :: itime
      double precision :: dtime
      nmvf = size(fvs,1); ndim = size(fvs,2); nsxb = size(fvs,3)
      call dtimer(dtime,itime,-1)
      write (unit=iunit) (((fvs(i,j,k),i=1,nmvf),j=1,ndim),k=1,nsxb)
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine setmbeam1(part,npx,irc)
! marks beam particles by setting a particle id in particle location
! 3 or 5
! irc = (0,1) = (no,yes) error condition exists
      integer, intent(in) :: npx
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop
      idimp = size(part,1); nop = size(part,2)
      irc = 0
! call low level procedure
      if (idimp > 4) then
         call STPBEAM13(part,npx,idimp,nop)
      else if (idimp > 2) then
         call STPBEAM1(part,npx,idimp,nop)
      else
         write (*,*) 'setbeamc1 error: idimp=', idimp
         irc = 1
      endif
      end subroutine
!
      end module
