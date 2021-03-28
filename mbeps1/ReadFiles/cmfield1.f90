!-----------------------------------------------------------------------
!
      module cmfield1
!
! Fortran90 wrappers to 1d OpenMP library cmfield1.f
! ffopen1 opens an old fortran formatted file
! fsopen1 opens an old fortran unformatted stream file
! readsdiags1 determine which scalar diagnostics are available
! readvdiags1 determine which vector diagnostics are available
! readfldiags1 determine which fluid diagnostics are available
! readfvdiags1 determine which velocity diagnostics are available
! readtrdiags1 determine which trajectory diagnostics are available
! readpsdiags1 determine which phase space diagnostics are available
! sdiagparams1 return parameters for selected scalar diagnostic
! vdiagparams1 return parameters for selected vector diagnostic
! fldiagparams1 return parameters for selected fluid diagnostic
! fvdiagparams1 return parameters for selected velocity diagnostic
! trdiagparams1 return parameters for selected trajectory diagnostic
! psdiagparams1 return parameters for selected phase space diagnostic
! freadc1 read complex scalar field
! freadvc1 read complex vector field
! fread1 read real scalar field
! freadv1 read real vector field
! freadfv1 read real velocity record
! freadtr1 read real trajectory record
! freadps1 read real phase space record
! rewindff1 rewind Fortran file
! closeff1 close Fortran file
! mwrmodes1 reads and copies lowest order scalar modes to packed array,
!           writing zeroes to high order modes
!           calls WRMODES1
! mwrvmodes1 reads and copies lowest order vector modes to packed array
!            writing zeroes to high order modes
!            calls WRVMODES1
! mcspect1 performs frequency analysis of complex time series
!          calls CSPECT1
! micspect1 performs incremental frequency analysis of complex scalar
!           time series for one time step
!           calls ICSPECT1
! mivcspect1 performs incremental frequency analysis of complex vector 
!            time series for one time step
!            calls IVCSPECT1
! mfft1_init calculates tables needed by 1d FFTs
!            calls WFFT1RINIT
! mfft1r wrapper function for in place scalar 1d real/complex FFT
!        calls FFT1RXX
! mfft1rn wrapper function for in place vector 1d real/complex FFT
!         calls FFT1R2X or FFT1R3X
! mdivf1 calculates the divergence in fourier space
!        calls DIVF1
! written by Viktor K. Decyk, UCLA
!
      use cmfield1_h
      implicit none
!
! t = scratch array for mfft1rn
      complex, dimension(:), allocatable :: t
      integer :: szt = 0
      save
!
      private :: t, szt
!
      contains
!
!-----------------------------------------------------------------------
      subroutine ffopen1(iunit,fname)
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
      subroutine fsopen1(iunit,fname)
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
      subroutine readsdiags1(iunit,nscalars)
! determine which scalar diagnostics are available
      use in1, only: pot1d, el1d, dene1d, deni1d
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
         read (iunit,pot1d,iostat=ios)
! load metadata for longitudinal efield data
      case (2)
         read (iunit,el1d,iostat=ios)
! load metadata for electron density data
      case (3)
         read (iunit,dene1d,iostat=ios)
! load metadata for ion density data
      case (4)
         read (iunit,deni1d,iostat=ios)
      end select
      if (ios==0) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readvdiags1(iunit,nscalars)
! determine which vector diagnostics are available
      use in1, only: vcure1d, vpot1d, et1d, b1d, vpotr1d, vcuri1d
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
! load metadata for electron current density data
      case (1)
         read (iunit,vcuri1d,iostat=ios)
! load metadata for vector potential data
      case (2)
         read (iunit,vpot1d,iostat=ios)
! load metadata for transverse efield data
      case (3)
         read (iunit,et1d,iostat=ios)
! load metadata for magnetic field data
      case (4)
         read (iunit,b1d,iostat=ios)
! load metadata for radiative vector potential data
      case (5)
         read (iunit,vpotr1d,iostat=ios)
! load metadata for ion current density data
      case (6)
         read (iunit,vcuri1d,iostat=ios)
      end select
      if (ios==0) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readfldiags1(iunit,nscalars)
! determine which fluid diagnostics are available
      use in1, only: fm1d, nferec, nfirec
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
         read (iunit,fm1d,iostat=ios)
         nrec = nferec
! load metadata for ion fluid moments data
      case (2)
         read (iunit,fm1d,iostat=ios)
         nrec = nfirec
      end select
      if ((ios==0).and.(nrec > 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readfvdiags1(iunit,nscalars)
! determine which velocity diagnostics are available
      use in1, only: fv1d, nverec, nvirec
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
         read (iunit,fv1d,iostat=ios)
         nrec = nverec
! load metadata for ion velocity-space data
      case (2)
         read (iunit,fv1d,iostat=ios)
         nrec = nvirec
      end select
      if ((ios==0).and.(nrec > 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readtrdiags1(iunit,nscalars)
! determine which trajectory diagnostics are available
      use in1, only: tr1d, ntrec
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
         read (iunit,tr1d,iostat=ios)
         nrec = ntrec
      end select
      if ((ios==0).and.(nrec > 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readpsdiags1(iunit,nscalars)
! determine which phase space diagnostics are available
      use in1, only: ps1d, nserec, nsirec
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
         read (iunit,ps1d,iostat=ios)
         nrec = nserec
! load metadata for ion phase-space data
      case (2)
         read (iunit,ps1d,iostat=ios)
         nrec = nsirec
      end select
      if ((ios==0).and.(nrec > 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sdiagparams1(iunit,n,mts,modesx,nrec,norm,fname)
! return parameters for selected scalar diagnostic
      use in1, only: pot1d, ntp, fpname, modesxp, nprec,                &
     &el1d, ntel, felname, modesxel, nelrec,                            &
     &dene1d, ntde, fdename, modesxde, nderec,                          &
     &deni1d, ntdi, fdiname, modesxdi, ndirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2,3)
! mts = number of time steps between selected diagnostic
! modesx = number of modes in selected diagnostic
! nrec = number of records in selected diagnostic
! norm = (-1,0,1) = normalize with (inverse gradient,null,gradient) op
! fname = file name for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, modesx, nrec, norm
      character(len=*), intent(inout) :: fname
! local data
      integer :: ios
      norm = 0
      select case(n)
! load metadata for potential data
      case (1)
         read (iunit,pot1d,iostat=ios)
         mts = ntp
         modesx = modesxp
         fname = fpname; nrec = nprec
         norm = 1
! load metadata for longitudinal efield data
      case (2)
         read (iunit,el1d,iostat=ios)
         mts = ntel
         modesx = modesxel
         fname = felname; nrec = nelrec
! load metadata for electron density data
      case (3)
         read (iunit,dene1d,iostat=ios)
         mts = ntde
         modesx = modesxde
         fname = fdename; nrec = nderec
! load metadata for ion density data
      case (4)
         read (iunit,deni1d,iostat=ios)
         mts = ntdi
         modesx = modesxdi
         fname = fdiname; nrec = ndirec
         norm = -1
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine vdiagparams1(iunit,n,mts,modesx,nrec,norm,fname)
! return parameters for selected vector diagnostic
      use in1, only: vcure1d, ntje, fjename, modesxje, njerec,          &
     &vpot1d, nta, faname, modesxa, narec,                              &
     &et1d, ntet, fetname, modesxet, netrec,                            &
     &b1d, ntb, fbname, modesxb, nbrec,                                 &
     &vpotr1d, ntar, farname, modesxar, narrec,                         &
     &vcuri1d, ntji, fjiname, modesxji, njirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2,3,4,5,6,7)
! mts = numbe of time steps between selected diagnostic
! modesx = number of modes in selected diagnostic
! nrec = number of records in selected diagnostic
! norm = (-1,0,1) = normalize with (inverse curl,null,curl) op
! fname = file name for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, modesx, nrec, norm
      character(len=*), intent(inout) :: fname
! local data
      integer :: ios
      norm = 0
      select case(n)
! load metadata for electron current density data
      case (1)
         read (iunit,vcure1d,iostat=ios)
         mts = ntje
         modesx = modesxje
         fname = fjename; nrec = njerec
         norm = -1
! load metadata for vector potential data
      case (2)
         read (iunit,vpot1d,iostat=ios)
         mts = nta
         modesx = modesxa
         fname = faname; nrec = narec
         norm = 1
! load metadata for transverse efield data
      case (3)
         read (iunit,et1d,iostat=ios)
         mts = ntet
         modesx = modesxet
         fname = fetname; nrec = netrec
! load metadata for magnetic field data
      case (4)
         read (iunit,b1d,iostat=ios)
         mts = ntb
         modesx = modesxb
         fname = fbname; nrec = nbrec
! load metadata for radiative vector potential data
      case (5)
         read (iunit,vpotr1d,iostat=ios)
         mts = ntar
         modesx = modesxar
         fname = farname; nrec = narrec
         norm = 1
! load metadata for ion current density data
      case (6)
         read (iunit,vcuri1d,iostat=ios)
         mts = ntji
         modesx = modesxji
         fname = fjiname; nrec = njirec
         norm = -1
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fldiagparams1(iunit,n,mts,mpro,mprd,nrec,fname)
! return parameters for selected fluid diagnostic
      use in1, only: fm1d, ntfm, npro, nprd, ffename, ffiname, nferec,  &
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
         read (iunit,fm1d,iostat=ios)
         mts = ntfm
         mpro = npro; mprd = nprd
         fname = ffename; nrec = nferec
! load metadata for ion fluid moments data
      case (2)
         read (iunit,fm1d,iostat=ios)
         mts = ntfm
         mpro = npro; mprd = nprd
         fname = ffiname; nrec = nfirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fvdiagparams1(iunit,n,mts,mmv,mfvd,mfed,nrec,fname)
! return parameters for selected velocity diagnostic
      use in1, only: fv1d, ntv, nmv, nfvd, nfed, fvename, fviname,&
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
      select case(n)
! load metadata for electron velocity-space data
      case (1)
         read (iunit,fv1d)
         mts = ntv; mmv = nmv
         mfvd = nfvd; mfed = nfed
         fname = fvename; nrec = nverec
! load metadata for ion velocity-space data
      case (2)
         read (iunit,fv1d)
         mts = ntv; mmv = nmv
         mfvd = nfvd; mfed = nfed
         fname = fviname; nrec = nvirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine trdiagparams1(iunit,n,mts,mdt,mst,mmv,mdimp,mprobt,nrec&
     &,fname)
! return parameters for selected trajectory diagnostic
      use in1, only: tr1d, ntt, ndt, nst, nmv, ndimp, nprobt, ftname,   &
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
         read (iunit,tr1d)
         mts = ntt; mdt = ndt
         mst = nst; mmv = nmv
         mdimp = ndimp; mprobt = nprobt
         fname = ftname; nrec = ntrec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine psdiagparams1(iunit,n,mts,mmv,msxb,nrec,fname)
! return parameters for selected phase space diagnostic
      use in1, only: ps1d, nts, nmv, nsxb, fsename, fsiname, nserec,    &
     &nsirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2)
! mts = number of time steps between selected diagnostic
! mmv = number of segments in v for velocity distribution
! mxb = number of segments in x for velocity distribution
! nrec = number of records in selected diagnostic
! fname = file names for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, mmv, msxb, nrec
      character(len=*), intent(inout) :: fname
! local data
      select case(n)
! load metadata for electron phase-space data
      case (1)
         read (iunit,ps1d)
         mts = nts; mmv = nmv
         msxb = nsxb;
         fname = fsename; nrec = nserec
! load metadata for ion phase-space data
      case (2)
         read (iunit,ps1d)
         mts = nts; mmv = nmv
         msxb = nsxb
         fname = fsiname; nrec = nsirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadc1(iunit,fc,nx)
! read complex scalar field
      implicit none
      integer, intent(in) :: iunit, nx
      complex, dimension(:), intent(inout) :: fc
! local data
      integer :: j, ios
      read (unit=iunit,iostat=ios) (fc(j),j=1,nx)
      if (ios /= 0) write (*,*) 'IO error in unit=',iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadvc1(iunit,fvc,ndim,nx)
! read complex vector field
      implicit none
      integer, intent(in) :: iunit, ndim, nx
      complex, dimension(:,:), intent(inout) :: fvc
! local data
      integer :: i, j, ios
      read (unit=iunit,iostat=ios) ((fvc(i,j),i=1,ndim),j=1,nx)
      if (ios /= 0) write (*,*) 'IO error in unit=',iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fread1(iunit,f,nx)
! read real scalar field
      implicit none
      integer, intent(in) :: iunit, nx
      real, dimension(:), intent(inout) :: f
! local data
      integer :: j, ios
      read (unit=iunit,iostat=ios) (f(j),j=1,nx)
      if (ios /= 0) write (*,*) 'IO error in unit=',iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadv1(iunit,fv,ndim,nx)
! read real vector field
      implicit none
      integer, intent(in) :: iunit, ndim, nx
      real, dimension(:,:), intent(inout) :: fv
! local data
      integer :: i, j, ios
      read (unit=iunit,iostat=ios) ((fv(i,j),i=1,ndim),j=1,nx)
      if (ios /= 0) write (*,*) 'IO error in unit=',iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadfv1(iunit,fvm,fv,fe,wk,ndim,nmvf,nfvd,nfed)
! read real velocity record from binary file
! fvm,fv,fe = real velocity data to be read in each record
! wk = total energy contained in distribution
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit, ndim, nmvf, nfvd, nfed
      real, intent(inout) :: wk
      real, dimension(:,:), intent(inout) :: fvm, fv, fe
! local data
      integer :: j, k, ios
      read (unit=iunit,iostat=ios) ((fvm(j,k),j=1,ndim),k=1,3),                    &
     &((fv(j,k),j=1,nmvf),k=1,nfvd), ((fe(j,k),j=1,nmvf),k=1,nfed), wk
      if (ios /= 0) write (*,*) 'IO error in unit=',iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadtr1(iunit,partt,idimp,nprobt)
! read real trajectory record from binary file
! partt = real trajectory data to be wto be read in each record
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit, idimp, nprobt
      real, dimension(:,:), intent(inout) :: partt
! local data
      integer :: j, k, ios
      read (unit=iunit,iostat=ios) ((partt(j,k),j=1,idimp),k=1,nprobt)
      if (ios /= 0) write (*,*) 'IO error in unit=',iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadps1(iunit,fvs,ndim,nmvf,nsxb)
! read real phase space record from binary file
! fvs = real velocity data to be read in each record
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit, ndim, nmvf, nsxb
      real, dimension(:,:,:), intent(inout) :: fvs
! local data
      integer :: i, j, k, ios
      read (unit=iunit,iostat=ios) (((fvs(i,j,k),i=1,nmvf),j=1,ndim),   &
     &k=1,nsxb)
      if (ios /= 0) write (*,*) 'IO error in unit=',iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine rewindff1(iunit)
! rewind Fortran file with unit number iunit
      implicit none
      integer, intent(in) :: iunit
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine closeff1(iunit)
! close Fortran file with unit number iunit
      implicit none
      integer, intent(in) :: iunit
      close(unit=iunit)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwrmodes1(pot,pott,nx,modesx)
! reads and copies lowest order scalar modes to packed array, writing
! zeroes to high order modes
      implicit none
      integer, intent(in) :: nx, modesx
      real, dimension(:), intent(inout) :: pot
      complex, dimension(:), intent(in) :: pott
! local data
      integer :: nxvh, modesxd
! extract dimensions
      nxvh = size(pot,1)/2
      modesxd = size(pott,1)
! call low level procedure
      call WRMODES1(pot,pott,nx,modesx,nxvh,modesxd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwrvmodes1(vpot,vpott,nx,modesx)
! extracts and stores lowest order vector modes to unpacked array
      implicit none
      integer, intent(in) :: nx, modesx
      real, dimension(:,:), intent(inout) :: vpot
      complex, dimension(:,:), intent(in) :: vpott
! local data
      integer :: ndim, nxvh, modesxd
! extract dimensions
      ndim = size(vpot,1); nxvh = size(vpot,2)/2
      modesxd = size(vpott,2)
! call low level procedure
      call WRVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcspect1(fc,wm,pkw,t0,dt,nt,iw,modesx)
! performs incremental frequency analysis of complex time series
      integer, intent(in) :: nt, iw, modesx
      real, intent(in) :: t0, dt
      complex, dimension(:,:), intent(in) :: fc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:), intent(inout) :: pkw
! local data
      integer :: ntd, iwd, modesxd
! extract dimensions
      ntd = size(fc,1); modesxd = size(fc,2)
      iwd = size(wm,1)
! call low level procedure
      call CSPECT1(fc,wm,pkw,t0,dt,nt,iw,modesx,ntd,iwd,modesxd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine micspect1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm)
! performs incremental frequency analysis of complex scalar time series
! for one time step
! norm = (-1,0,1) = normalize with (inverse gradient,null,gradient) op
      integer, intent(in) :: nt, iw, modesx, nx, norm
      real, intent(in) :: time, t0
      complex, dimension(:), intent(in) :: fc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:), intent(inout) :: pkw
      double precision, dimension(:,:,:), intent(inout) :: pks
! local data
      integer :: iwd, modesxd
! extract dimensions
      modesxd = size(fc,1); iwd = size(wm,1)
! call low level procedure
      call ICSPECT1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm,iwd,     &
     &modesxd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mivcspect1(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,nx,   &
     &norm)
! performs incremental frequency analysis of complex vector time series
! for one time step
! norm = (-1,0,1) = normalize with (inverse curl,null,curl) op
      integer, intent(in) :: nt, iw, modesx, nx, norm
      real, intent(in) :: time, t0
      complex, dimension(:,:), intent(in) :: fvc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:,:), intent(inout) :: vpkw
      double precision, dimension(:,:,:,:), intent(inout) :: vpks
! local data
      integer :: iwd, modesxd
! extract dimensions
      modesxd = size(fvc,2); iwd = size(wm,1)
! call low level procedure
      call IVCSPECT1(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,nx,norm,iwd, &
     &modesxd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft1_init(mixup,sct,indx)
! calculates tables needed by 1d FFTs
      implicit none
      integer, intent(in) :: indx
      integer, dimension(:), intent(inout) :: mixup
      complex, dimension(:), intent(inout) :: sct
! local data
      integer :: nxhd
! extract dimensions
      nxhd = size(mixup,1)
! call low level procedure
      call WFFT1RINIT(mixup,sct,indx,nxhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft1r(f,isign,mixup,sct,indx)
! wrapper function for scalar 1d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx
      real, dimension(:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
! local data
      integer :: nxd, nxhd
! extract dimensions
      nxd = size(f,1)
      nxhd = size(mixup,1)
! check if required size of buffer has increased
      if (szt < nxhd) then
         if (szt /= 0) deallocate(t)
! allocate new buffers
         allocate(t(nxhd))
         szt = nxhd
      endif
! call low level procedure
      call FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft1rn(f,isign,mixup,sct,indx)
! wrapper function for n component vector 1d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx
      real, dimension(:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
! local data
      integer :: ndim, nxd, nxhd
! extract dimensions
      ndim = size(f,1); nxd = size(f,2)
      nxhd = size(mixup,1)
! check if required size of buffer has increased
      if (szt < ndim*nxhd) then
         if (szt /= 0) deallocate(t)
! allocate new buffers
         allocate(t(ndim*nxhd))
         szt = ndim*nxhd
      endif
! call low level procedure
      select case(ndim)
      case (2)
         call FFT1R2X(f,t,isign,mixup,sct,indx,nxd,nxhd)
      case (3)
         call FFT1R3X(f,t,isign,mixup,sct,indx,nxd,nxhd)
      end select
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mdivf1(f,df,nx)
! calculates the divergence in fourier space
      implicit none
      integer, intent(in) :: nx
      real, dimension(:,:), intent(in) :: f
      real, dimension(:), intent(inout) :: df
! local data
      integer :: ndim, nxvh
! extract dimensions
      ndim = size(f,1); nxvh = size(f,2)/2
! call low level procedure
      call DIVF1(f,df,nx,ndim,nxvh)
      end subroutine
!
      end module
