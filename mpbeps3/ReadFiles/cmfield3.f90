!-----------------------------------------------------------------------
!
      module cmfield3
!
! Fortran90 wrappers to 3d OpenMP library cmfield3.f
! ffopen3 opens an old fortran formatted file
! fsopen3 opens an old fortran unformatted stream file
! readsdiags3 determine which scalar diagnostics are available
! readvdiags3 determine which vector diagnostics are available
! readfldiags3 determine which fluid diagnostics are available
! readfvdiags3 determine which velocity diagnostics are available
! readtrdiags3 determine which trajectory diagnostics are available
! readpsdiags3 determine which phase space diagnostics are available
! sdiagparams3 return parameters for selected scalar diagnostic
! vdiagparams3 return parameters for selected vector diagnostic
! fldiagparams3 return parameters for selected fluid diagnostic
! fvdiagparams3 return parameters for selected velocity diagnostic
! trdiagparams3 return parameters for selected trajectory diagnostic
! psdiagparams3 return parameters for selected phase space diagnostic
! freadc3 read complex scalar field
! freadvc3 read complex vector field
! fread3 read real scalar field
! freadv3 read real vector field
! freadfv3 read real velocity record
! freadtr3 read real trajectory record
! freadncomp3 read distributed non-uniform partition information
! freadps3 read real phase space record
! rewindff3 rewind Fortran file
! closeff3 close Fortran file
! mwrmodes3 reads and copies lowest order scalar modes to packed array,
!           writing zeroes to high order modes
!           calls WRMODES2
! mwrvmodes3 reads and copies lowest order vector modes to packed array,
!            writing zeroes to high order modes
!            calls WRVMODES2
! mcspect2 performs frequency analysis of complex time series
!          calls CSPECT2
! micspect2 performs incremental frequency analysis of complex scalar
!           time series for one time step
!           calls ICSPECT2
! mivcspect2 performs incremental frequency analysis of complex vector 
!            time series for one time step
!            calls IVCSPECT2
! mfft3_init calculates tables needed by 3d FFTs
!            calls WFFT3RINIT
! mfft3r wrapper function for scalar 3d real/complex FFT
!        calls  WFFT3RMX
! mfft3r3 wrapper function for vector 3d real/complex FFT
!         calls WFFT3RM3
! written by Viktor K. Decyk, UCLA
!
      use cmfield3_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine ffopen3(iunit,fname)
! this subroutine opens an old fortran formatted file
! iunit = fortran unit number to be used 
! fname = file name
      implicit none
      integer, intent(in) :: iunit
      character(len=*), intent(in) :: fname
! local data
      open(unit=iunit,file=fname,form='formatted',status='old')
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fsopen3(iunit,fname)
! this subroutine opens an old fortran unformatted stream file
! iunit = fortran unit number to be used 
! fname = file name
      implicit none
      integer, intent(in) :: iunit
      character(len=*), intent(in) :: fname
! local data
      open(unit=iunit,file=fname,access='stream',form='unformatted',    &
     &status='old')
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readsdiags3(iunit,nscalars)
! determine which scalar diagnostics are available
      use in3, only: pot3d, dene3d, deni3d
      implicit none
! iunit = fortran unit number to be used 
      integer, intent(in) :: iunit
! nscalars = table of available diagnositcs
      integer, dimension(:), intent(inout) :: nscalars
! local data
      integer :: n, ns, ios
      ns = size(nscalars)
! check each case
      do n = 1, ns
      select case(n)
! load metadata for potential data
      case (1)
         read (iunit,pot3d,iostat=ios)
! load metadata for electron density data
      case (2)
         read (iunit,dene3d,iostat=ios)
! load metadata for ion density data
      case (3)
         read (iunit,deni3d,iostat=ios)
      end select
      if (ios==0) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readvdiags3(iunit,nscalars)
! determine which vector diagnostics are available
      use in3, only: el3d, vcure3d, vpot3d, et3d, b3d, vpotr3d, vcuri3d
      implicit none
! iunit = fortran unit number to be used 
      integer, intent(in) :: iunit
! nscalars = table of available diagnositcs
      integer, dimension(:), intent(inout) :: nscalars
! local data
      integer :: n, ns, ios
      ns = size(nscalars)
! check each case
      do n = 1, ns
      select case(n)
! load metadata for longitudinal efield data
      case (1)
         read (iunit,el3d,iostat=ios)
! load metadata for electron current density data
      case (2)
         read (iunit,vcuri3d,iostat=ios)
! load metadata for vector potential data
      case (3)
         read (iunit,vpot3d,iostat=ios)
! load metadata for transverse efield data
      case (4)
         read (iunit,et3d,iostat=ios)
! load metadata for magnetic field data
      case (5)
         read (iunit,b3d,iostat=ios)
! load metadata for radiative vector potential data
      case (6)
         read (iunit,vpotr3d,iostat=ios)
! load metadata for ion current density data
      case (7)
         read (iunit,vcuri3d,iostat=ios)
      end select
      if (ios==0) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readfldiags3(iunit,nscalars)
! determine which fluid diagnostics are available
      use in3, only: fm3d, nferec, nfirec
      implicit none
! iunit = fortran unit number to be used 
      integer, intent(in) :: iunit
! nscalars = table of available diagnositcs
      integer, dimension(:), intent(inout) :: nscalars
! local data
      integer :: n, ns, ios, nrec
      ns = size(nscalars)
! check each case
      do n = 1, ns
      nrec = -1
      select case(n)
! load metadata for electron fluid moments data
      case (1)
         read (iunit,fm3d,iostat=ios)
         nrec = nferec
! load metadata for ion fluid moments data
      case (2)
         read (iunit,fm3d,iostat=ios)
         nrec = nfirec
      end select
      if ((ios==0).and.(nrec >= 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readfvdiags3(iunit,nscalars)
! determine which velocity diagnostics are available
      use in3, only: fv3d, nverec, nvirec
      implicit none
! iunit = fortran unit number to be used 
      integer, intent(in) :: iunit
! nscalars = table of available diagnositcs
      integer, dimension(:), intent(inout) :: nscalars
! local data
      integer :: n, ns, ios, nrec
      ns = size(nscalars)
! check each case
      do n = 1, ns
      nrec = -1
      select case(n)
! load metadata for electron velocity-space data
      case (1)
         read (iunit,fv3d,iostat=ios)
         nrec = nverec
! load metadata for ion velocity-space data
      case (2)
         read (iunit,fv3d,iostat=ios)
         nrec = nvirec
      end select
      if ((ios==0).and.(nrec >= 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readtrdiags3(iunit,nscalars)
! determine which trajectory diagnostics are available
      use in3, only: tr3d, ntrec
      implicit none
! iunit = fortran unit number to be used 
      integer, intent(in) :: iunit
! nscalars = table of available diagnositcs
      integer, dimension(:), intent(inout) :: nscalars
! local data
      integer :: n, ns, ios, nrec
      ns = size(nscalars)
! check each case
      do n = 1, ns
      nrec = -1
      select case(n)
! load metadata for trajectory data
      case (1)
         read (iunit,tr3d,iostat=ios)
         nrec = ntrec
      end select
      if ((ios==0).and.(nrec >= 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readpsdiags3(iunit,nscalars)
! determine which phase space diagnostics are available
      use in3, only: ps3d, nserec, nsirec
      implicit none
! iunit = fortran unit number to be used 
      integer, intent(in) :: iunit
! nscalars = table of available diagnositcs
      integer, dimension(:), intent(inout) :: nscalars
! local data
      integer :: n, ns, ios, nrec
      ns = size(nscalars)
! check each case
      do n = 1, ns
      nrec = -1
      select case(n)
! load metadata for electron phase-space data
      case (1)
         read (iunit,ps3d,iostat=ios)
         nrec = nserec
! load metadata for ion phase-space data
      case (2)
         read (iunit,ps3d,iostat=ios)
         nrec = nsirec
      end select
      if ((ios==0).and.(nrec >= 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sdiagparams3(iunit,n,mts,modesx,modesy,modesz,nrec,    &
     &fname)
! return parameters for selected scalar diagnostic
      use in3, only:                                                    &
     &pot3d, ntp, fpname, modesxp, modesyp, modeszp, nprec,             &
     &dene3d, ntde, fdename, modesxde, modesyde, modeszde, nderec,      &
     &deni3d, ntdi, fdiname, modesxdi, modesydi, modeszdi, ndirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2,3)
! mts = number of time steps between selected diagnostic
! modesx/modesy/modesz = number of modes in selected diagnostic in x/y/z
! nrec = number of records in selected diagnostic
! fname = file name for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, modesx, modesy, modesz, nrec
      character(len=*), intent(inout) :: fname
! local data
      select case(n)
! load metadata for potential data
      case (1)
         read (iunit,pot3d)
         mts = ntp
         modesx = modesxp; modesy = modesxp; modesz = modeszp
         fname = fpname; nrec = nprec
! load metadata for electron density data
      case (2)
         read (iunit,dene3d)
         mts = ntde
         modesx = modesxde; modesy = modesyde; modesz = modeszde
         fname = fdename; nrec = nderec
! load metadata for ion density data
      case (3)
         read (iunit,deni3d)
         mts = ntdi
         modesx = modesxdi; modesy = modesydi; modesz = modeszdi
         fname = fdiname; nrec = ndirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine vdiagparams3(iunit,n,mts,modesx,modesy,modesz,nrec,    &
     &fname)
! return parameters for selected vector diagnostic
      use in3, only:                                                    &
     &el3d, ntel, felname, modesxel, modesyel, modeszel, nelrec,        &
     &vcure3d, ntje, fjename, modesxje, modesyje, modeszje, njerec,     &
     &vpot3d, nta, faname, modesxa, modesya, modesza, narec,            &
     &et3d, ntet, fetname, modesxet, modesyet, modeszet, netrec,        &
     &b3d, ntb, fbname, modesxb, modesyb, modeszb, nbrec,               &
     &vpotr3d, ntar, farname, modesxar, modesyar, modeszar, narrec,     &
     &vcuri3d, ntji, fjiname, modesxji, modesyji, modeszji, njirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2,3,4,5,6,7)
! mts = numbe of time steps between selected diagnostic
! modesx/modesy/modesz = number of modes in selected diagnostic in x/y/z
! nrec = number of records in selected diagnostic
! fname = file name for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, modesx, modesy, modesz, nrec
      character(len=*), intent(inout) :: fname
! local data
      select case(n)
! load metadata for longitudinal efield data
      case (1)
         read (iunit,el3d)
         mts = ntel
         modesx = modesxel; modesy = modesyel; modesz = modeszel
         fname = felname; nrec = nelrec
! load metadata for electron current density data
      case (2)
         read (iunit,vcure3d)
         mts = ntje
         modesx = modesxje; modesy = modesyje; modesz = modeszje
         fname = fjename; nrec = njerec
! load metadata for vector potential data
      case (3)
         read (iunit,vpot3d)
         mts = nta
         modesx = modesxa; modesy = modesya; modesz = modesza
         fname = faname; nrec = narec
! load metadata for transverse efield data
      case (4)
         read (iunit,et3d)
         mts = ntet
         modesx = modesxet; modesy = modesyet; modesz = modeszet
         fname = fetname; nrec = netrec
! load metadata for magnetic field data
      case (5)
         read (iunit,b3d)
         mts = ntb
         modesx = modesxb; modesy =  modesyb; modesz =  modeszb
         fname = fbname; nrec = nbrec
! load metadata for radiative vector potential data
      case (6)
         read (iunit,vpotr3d)
         mts = ntar
         modesx = modesxar; modesy = modesyar; modesz = modeszar
         fname = farname; nrec = narrec
! load metadata for ion current density data
      case (7)
         read (iunit,vcuri3d)
         mts = ntji
         modesx = modesxji; modesy = modesyji; modesz = modeszji
         fname = fjiname; nrec = njirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fldiagparams3(iunit,n,mts,mpro,mprd,nrec,fname)
! return parameters for selected fluid diagnostic
      use in3, only: fm3d, ntfm, npro, nprd, ffename, ffiname, nferec,  &
     &nfirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2)
! mts = number of time steps between selected diagnostic
! mpro = number of moments calculated
! mprd = dimension of fluid moment arrays fms
! nrec = number of records in selected diagnostic
! fname = file names for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, mpro, mprd, nrec
      character(len=*), intent(inout) :: fname
! local data
      integer :: ios
      select case(n)
! load metadata for electron fluid moments data
      case (1)
         read (iunit,fm3d,iostat=ios)
         mts = ntfm
         mpro = npro; mprd = nprd
         fname = ffename; nrec = nferec
! load metadata for ion fluid moments data
      case (2)
         read (iunit,fm3d,iostat=ios)
         mts = ntfm
         mpro = npro; mprd = nprd
         fname = ffiname; nrec = nfirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fvdiagparams3(iunit,n,mts,mmv,mfvd,mfed,nrec,fname)
! return parameters for selected velocity diagnostic
      use in3, only: fv3d, ntv, nmv, nfvd, nfed, fvename, fviname,&
     &nverec, nvirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2)
! mts = number of time steps between selected diagnostic
! mmv = number of segments in v for velocity distribution
! mfvd = number of moments calculated
! mfed = dimension of fluid moment arrays fms
! nrec = number of records in selected diagnostic
! fname = file names for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, mmv, mfvd, mfed, nrec
      character(len=*), intent(inout) :: fname
! local data
      integer :: ios
      select case(n)
! load metadata for electron velocity-space data
      case (1)
         read (iunit,fv3d,iostat=ios)
         mts = ntv; mmv = nmv
         mfvd = nfvd; mfed = nfed
         fname = fvename; nrec = nverec
! load metadata for ion velocity-space data
      case (2)
         read (iunit,fv3d,iostat=ios)
         mts = ntv; mmv = nmv
         mfvd = nfvd; mfed = nfed
         fname = fviname; nrec = nvirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine trdiagparams3(iunit,n,mts,mdt,mst,mmv,mdimp,mprobt,nrec&
     &,fname)
! return parameters for selected trajectory diagnostic
      use in3, only: tr3d, ntt, ndt, nst, nmv, ndimp, nprobt, ftname,   &
     &ntrec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1)
! mts = number of time steps between selected diagnostic
! mdt = electron or ion selector
! mst = type of test particle distribution:
! mmv = number of segments in v for test particle distribution
! mdimp = dimension of trajectory array partt
! nrec = number of records in selected diagnostic
! fname = file names for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, mdt, mst, mmv, mdimp, mprobt, nrec
      character(len=*), intent(inout) :: fname
! local data
      select case(n)
! load metadata for trajectory data
      case (1)
         read (iunit,tr3d)
         mts = ntt; mdt = ndt
         mst = nst; mmv = nmv
         mdimp = ndimp; mprobt = nprobt
         fname = ftname; nrec = ntrec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine psdiagparams3(iunit,n,mts,mmv,msxb,msyb,mszb,nrec,fname&
     &)
! return parameters for selected phase space diagnostic
      use in3, only: ps3d, nts, nmv, nsxb, nsyb, nszb, fsename, fsiname,&
     &nserec, nsirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2)
! mts = number of time steps between selected diagnostic
! mmv = number of segments in v for velocity distribution
! mxb/msyb/mszb = number of segments in x/y/z for velocity distribution
! nrec = number of records in selected diagnostic
! fname = file names for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, mmv, msxb, msyb, mszb, nrec
      character(len=*), intent(inout) :: fname
! local data
      select case(n)
! load metadata for electron phase-space data
      case (1)
         read (iunit,ps3d)
         mts = nts; mmv = nmv
         msxb = nsxb; msyb = nsyb; mszb = nszb
         fname = fsename; nrec = nserec
! load metadata for ion phase-space data
      case (2)
         read (iunit,ps3d)
         mts = nts; mmv = nmv
         msxb = nsxb; msyb = nsyb; mszb = nszb
         fname = fsiname; nrec = nsirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fread3(iunit,f,nx,kyp,kyb,kzp,kzb)
! read and transpose scalar data
      implicit none
      integer, intent(in) :: iunit, nx, kyp, kyb, kzp, kzb
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: j, k, l, n, m
      read (unit=iunit) (((((f(j,k+kyp*(n-1),l+kzp*(m-1)),j=1,nx),      &
     &k=1,kyp),l=1,kzp),n=1,kyb),m=1,kzb)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadv3(iunit,fv,ndim,nx,kyp,kyb,kzp,kzb)
! read and transpose vector data
      implicit none
      integer, intent(in) :: iunit, ndim, nx, kyp, kyb, kzp, kzb
      real, dimension(:,:,:,:), intent(inout) :: fv
! local data
      integer :: i, j, k, l, n, m
      read (unit=iunit) ((((((fv(i,j,k+kyp*(n-1),l+kzp*(m-1)),i=1,ndim),&
     &j=1,nx),k=1,kyp),l=1,kzp),n=1,kyb),m=1,kzb)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadfv3(iunit,fvm,fv,fe,wk,ndim,nmvf,nfvd,nfed)
! read real velocity record from binary file
! fvm,fv,fe = real velocity data to be read in each record
! wk = total energy contained in distribution
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit, ndim, nmvf,nfvd,nfed
      real, intent(inout) :: wk
      real, dimension(:,:), intent(inout) :: fvm, fv, fe
! local data
      integer :: j, k
      read (unit=iunit) ((fvm(j,k),j=1,ndim),k=1,3),                    &
     &((fv(j,k),j=1,nmvf),k=1,nfvd), ((fe(j,k),j=1,nmvf),k=1,nfed), wk
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadtr3(iunit,partt,idimp,nprobt)
! read real trajectory record from binary file
! partt = real trajectory data to be wto be read in each record
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit, idimp, nprobt
      real, dimension(:,:), intent(inout) :: partt
! local data
      integer :: j, k
      read (unit=iunit) ((partt(j,k),j=1,idimp),k=1,nprobt)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadncomp3(iunit,noff,nyzp,ierr)
! reads distributed non-uniform partition information from binary file
      implicit none
      integer, intent(in) :: iunit
      integer, intent(inout) :: ierr
      integer, dimension(:,:,:), intent(inout) :: noff, nyzp
! local data
      integer :: j, k, nvpy ,nvpz, mvpy, mvpz, nps
      ierr = 0
      nvpy = size(noff,2); nvpz = size(noff,3)
      read (unit=iunit) mvpy, mvpz
! check for domain size error
      if ((mvpy /= nvpy).or.(mvpz /= nvpz).or.(size(noff,1) /= 2)) then
         if (mvpy /= nvpy) then
            write (*,*) 'partition error: nvpy, mvpy', nvpy, mvpy
            ierr = 1
         endif
         if (mvpz /= nvpz) then
            write (*,*) 'partition error: nvpz, mvpz', nvpz, mvpz
            ierr = ierr + 2
         endif
         if (size(noff,1) /= 2) then
            ierr = ierr + 4
         endif
         return
      endif
! calculate offsets in y
      do k = 1, nvpz
      nps = 0
      do j = 1, nvpy
      read (unit=iunit) nyzp(:,j,k)
      noff(1,j,k) = nps
      nps = nps + nyzp(1,j,k)
      enddo
      enddo
! calculate offsets in z
      do j = 1, nvpy
      nps = 0
      do k = 1, nvpz
      noff(2,j,k) = nps
      nps = nps + nyzp(2,j,k)
      enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadps3(iunit,fvs,noff,nyzp,ndim,nmvf,nsxb,nvpy,nvpz)
! read real phase space record from binary file
! fvs = real velocity data to be read in each record
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit, ndim, nmvf, nsxb, nvpy, nvpz
      real, dimension(:,:,:,:,:), intent(inout) :: fvs
      integer, dimension(2,nvpy,nvpz), intent(in) :: noff, nyzp
! local data
      integer :: jj, i, j, k, l, n, m
      read (unit=iunit) (((((((fvs(jj,i,j,k+noff(1,n,m),l+noff(2,n,m)), &
     &jj=1,nmvf),i=1,ndim),j=1,nsxb),k=1,nyzp(1,n,m)),l=1,nyzp(2,n,m)), &
     &n=1,nvpy),m=1,nvpz)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine rewindff3(iunit)
! rewind Fortran file with unit number iunit
      implicit none
      integer, intent(in) :: iunit
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine closeff3(iunit)
! close Fortran file with unit number iunit
      implicit none
      integer, intent(in) :: iunit
      close(unit=iunit)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine summit2(fvs,pvs,in1,in2)
! sum 5d array fvs over two indices, where 3 <= in1, in2 <= 5
      implicit none
      integer, intent(in) :: in1, in2
      real, dimension(:,:,:,:,:), intent(in) :: fvs
      real, dimension(:,:,:), intent(inout) :: pvs
! local data
      integer :: i, j, k, l, m
      double precision :: dsum
      if (((in1==3).and.(in2==4)).or.((in1==4).and.(in2==3))) then
         do m = 1, size(fvs,5)
         do j = 1, size(fvs,2)
         do i = 1, size(fvs,1)
         dsum = 0.0d0
         do l = 1, size(fvs,4)
         do k = 1, size(fvs,3)
         dsum = dsum + fvs(i,j,k,l,m)
         enddo
         enddo
         pvs(i,j,m) = dsum
         enddo
         enddo
         enddo
      else if (((in1==3).and.(in2==5)).or.((in1==5).and.(in2==3))) then
         do l = 1, size(fvs,4)
         do j = 1, size(fvs,2)
         do i = 1, size(fvs,1)
         dsum = 0.0d0
         do m = 1, size(fvs,5)
         do k = 1, size(fvs,3)
         dsum = dsum + fvs(i,j,k,l,m)
         enddo
         enddo
         pvs(i,j,l) = dsum
         enddo
         enddo
         enddo
      else if (((in1==4).and.(in2==5)).or.((in1==5).and.(in2==4))) then
         do k = 1, size(fvs,3)
         do j = 1, size(fvs,2)
         do i = 1, size(fvs,1)
         dsum = 0.0d0
         do m = 1, size(fvs,5)
         do l = 1, size(fvs,4)
         dsum = dsum + fvs(i,j,k,l,m)
         enddo
         enddo
         pvs(i,j,k) = dsum
         enddo
         enddo
         enddo
      else
         pvs = 0.0
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwrmodes3(pot,pott,nx,ny,nz,modesx,modesy,modesz)
! reads and copies lowest order scalar modes to packed array, writing
! zeroes to high order modes
      implicit none
      integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
      real, dimension(:,:,:), intent(inout) :: pot
      complex, dimension(:,:,:), intent(in) :: pott
! local data
      integer :: nxvh, nyv, nzv, modesxd, modesyd, modeszd
! extract dimensions
      nxvh = size(pot,1)/2; nyv = size(pot,2); nzv = size(pot,3)
      modesxd = size(pott,1); modesyd = size(pott,2)
! call low level procedure
      call WRMODES3(pot,pott,nx,ny,nz,modesx,modesy,modesz,nxvh,nyv,nzv,&
     &modesxd,modesyd,modeszd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwrvmodes3(vpot,vpott,nx,ny,nz,modesx,modesy,modesz)
! reads and copies lowest order vector modes to packed array, writing
! zeroes to high order modes
      implicit none
      integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
      real, dimension(:,:,:,:), intent(inout) :: vpot
      complex, dimension(:,:,:,:), intent(in) :: vpott
! local data
      integer :: ndim, nxvh, nyv, nzv, modesxd, modesyd, modeszd
! extract dimensions
      ndim = size(vpot,1); nxvh = size(vpot,2)/2
      nyv = size(vpot,3); nzv = size(vpot,4)
      modesxd = size(vpott,2); modesyd = size(vpott,3)
      modeszd = size(vpott,4)
! call low level procedure
      call WRVMODES3(vpot,vpott,nx,ny,nz,modesx,modesy,modesz,ndim,nxvh,&
     &nyv,nzv,modesxd,modesyd,modeszd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcspect2(fc,wm,pkw,t0,dt,nt,iw,modesx,modesy2)
! performs incremental frequency analysis of complex time series
      integer, intent(in) :: nt, iw, modesx, modesy2
      real, intent(in) :: t0, dt
      complex, dimension(:,:,:), intent(in) :: fc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:,:), intent(inout) :: pkw
! local data
      integer :: ntd, iwd, modesxd, modesyd
! extract dimensions
      ntd = size(fc,1); modesxd = size(fc,2); modesyd = size(fc,3)
      iwd = size(wm,1)
! call low level procedure
      call CSPECT2(fc,wm,pkw,t0,dt,nt,iw,modesx,modesy2,ntd,iwd,modesxd,&
     &modesyd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine micspect2(fc,wm,pkw,pks,time,t0,nt,iw,modesx,modesy)
! performs incremental frequency analysis of complex scalar time series
! for one time step
      integer, intent(in) :: nt, iw, modesx, modesy
      real, intent(in) :: time, t0
      complex, dimension(:,:), intent(in) :: fc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:,:), intent(inout) :: pkw
      double precision, dimension(:,:,:,:), intent(inout) :: pks
! local data
      integer :: iwd, modesxd, modesyd
! extract dimensions
      modesxd = size(fc,1); modesyd = size(fc,2); iwd = size(wm,1)
! call low level procedure
      call ICSPECT2(fc,wm,pkw,pks,time,t0,nt,iw,modesx,modesy,iwd,      &
     &modesxd,modesyd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mivcspect2(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,      &
     &modesy)
! performs incremental frequency analysis of complex vector time series
! for one time step
      integer, intent(in) :: nt, iw, modesx, modesy
      real, intent(in) :: time, t0
      complex, dimension(:,:,:), intent(in) :: fvc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:,:,:), intent(inout) :: vpkw
      double precision, dimension(:,:,:,:,:), intent(inout) :: vpks
! local data
      integer :: ndim, iwd, modesxd, modesyd
! extract dimensions
      ndim = size(fvc,1); modesxd = size(fvc,2); modesyd = size(fvc,3)
      iwd = size(wm,1)
! call low level procedure
      call IVCSPECT2(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,modesy,ndim, &
     &iwd,modesxd,modesyd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft3_init(mixup,sct,indx,indy,indz)
! calculates tables needed by 3d FFTs
      implicit none
      integer, intent(in) :: indx, indy, indz
      integer, dimension(:), intent(inout) :: mixup
      complex, dimension(:), intent(inout) :: sct
! local data
      integer :: nxhyzd, nxyzhd
! extract dimensions
      nxhyzd = size(mixup,1); nxyzhd = size(sct,1)
! call low level procedure
      call WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft3r(f,isign,mixup,sct,indx,indy,indz)
! wrapper function for scalar 3d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy, indz
      real, dimension(:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
! local data
      integer :: nxhd, nyd, nzd, nxhyzd, nxyzhd
! extract dimensions
      nxhd = size(f,1)/2; nyd = size(f,2); nzd = size(f,3)
      nxhyzd = size(mixup,1); nxyzhd = size(sct,1)
! call low level procedure
      call WFFT3RMX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd,nxhyzd&
     &,nxyzhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft3r3(f,isign,mixup,sct,indx,indy,indz)
! wrapper function for vector 3d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy, indz
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
! local data
      integer :: ndim, nxhd, nyd, nzd, nxhyzd, nxyzhd
! extract dimensions
      ndim = size(f,1); nxhd = size(f,2)/2
      nyd = size(f,3); nzd = size(f,4)
      nxhyzd = size(mixup,1); nxyzhd = size(sct,1)
! call low level procedure
      call WFFT3RM3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd,nxhyzd&
     &,nxyzhd)
      end subroutine
!
      end module
