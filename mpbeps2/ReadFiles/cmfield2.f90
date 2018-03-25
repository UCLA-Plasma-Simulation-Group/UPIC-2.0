!-----------------------------------------------------------------------
!
      module cmfield2
!
! Fortran90 wrappers to 2d OpenMP library cmfield2.f
! ffopen2 opens an old fortran formatted file
! fsopen2 opens an old fortran unformatted stream file
! readsdiags2 determine which scalar diagnostics are available
! readvdiags2 determine which vector diagnostics are available
! readfldiags2 determine which fluid diagnostics are available
! readfvdiags2 determine which velocity diagnostics are available
! readtrdiags2 determine which trajectory diagnostics are available
! readpsdiags2 determine which phase space diagnostics are available
! sdiagparams2 return parameters for selected scalar diagnostic
! vdiagparams2 return parameters for selected vector diagnostic
! fldiagparams2 return parameters for selected fluid diagnostic
! fvdiagparams2 return parameters for selected velocity diagnostic
! trdiagparams2 return parameters for selected trajectory diagnostic
! psdiagparams2 return parameters for selected phase space diagnostic
! freadc2 read complex scalar field
! freadvc2 read complex vector field
! fread2 read real scalar field
! freadv2 read real vector field
! freadfv2 read real velocity record
! freadtr3 read real trajectory record
! freadps2 read real phase space record
! rewindff2 rewind Fortran file
! closeff2 close Fortran file
! mwrmodes2 reads and copies lowest order scalar modes to packed array,
!           writing zeroes to high order modes
!           calls WRMODES2
! mwrvmodes2 reads and copies lowest order vector modes to packed array,
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
! mfft2_init calculates tables needed by 2d FFTs
!            calls WFFT2RINIT
! mfft2r wrapper function for scalar 2d real/complex FFT
!        calls  WFFT2RMX
! mfft2rn wrapper function for vector 2d real/complex FFT
!         calls WWFFT2RM2 or WFFT2RM3
! written by Viktor K. Decyk, UCLA
!
      use cmfield2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine ffopen2(iunit,fname)
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
      subroutine fsopen2(iunit,fname)
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
      subroutine readsdiags2(iunit,nscalars)
! determine which scalar diagnostics are available
      use in2, only: pot2d, dene2d, deni2d
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
         read (iunit,pot2d,iostat=ios)
! load metadata for electron density data
      case (2)
         read (iunit,dene2d,iostat=ios)
! load metadata for ion density data
      case (3)
         read (iunit,deni2d,iostat=ios)
      end select
      if (ios==0) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readvdiags2(iunit,nscalars)
! determine which vector diagnostics are available
      use in2, only: el2d, vcure2d, vpot2d, et2d, b2d, vpotr2d, vcuri2d
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
         read (iunit,el2d,iostat=ios)
! load metadata for electron current density data
      case (2)
         read (iunit,vcuri2d,iostat=ios)
! load metadata for vector potential data
      case (3)
         read (iunit,vpot2d,iostat=ios)
! load metadata for transverse efield data
      case (4)
         read (iunit,et2d,iostat=ios)
! load metadata for magnetic field data
      case (5)
         read (iunit,b2d,iostat=ios)
! load metadata for radiative vector potential data
      case (6)
         read (iunit,vpotr2d,iostat=ios)
! load metadata for ion current density data
      case (7)
         read (iunit,vcuri2d,iostat=ios)
      end select
      if (ios==0) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readfldiags2(iunit,nscalars)
! determine which fluid diagnostics are available
      use in2, only: fm2d, nferec, nfirec
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
         read (iunit,fm2d,iostat=ios)
         nrec = nferec
! load metadata for ion fluid moments data
      case (2)
         read (iunit,fm2d,iostat=ios)
         nrec = nfirec
      end select
      if ((ios==0).and.(nrec >= 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readfvdiags2(iunit,nscalars)
! determine which velocity diagnostics are available
      use in2, only: fv2d, nverec, nvirec
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
         read (iunit,fv2d,iostat=ios)
         nrec = nverec
! load metadata for ion velocity-space data
      case (2)
         read (iunit,fv2d,iostat=ios)
         nrec = nvirec
      end select
      if ((ios==0).and.(nrec >= 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readtrdiags2(iunit,nscalars)
! determine which trajectory diagnostics are available
      use in2, only: tr2d, ntrec
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
         read (iunit,tr2d,iostat=ios)
         nrec = ntrec
      end select
      if ((ios==0).and.(nrec >= 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine readpsdiags2(iunit,nscalars)
! determine which phase space diagnostics are available
      use in2, only: ps2d, nserec, nsirec
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
         read (iunit,ps2d,iostat=ios)
         nrec = nserec
! load metadata for ion phase-space data
      case (2)
         read (iunit,ps2d,iostat=ios)
         nrec = nsirec
      end select
      if ((ios==0).and.(nrec >= 0)) nscalars(n) = 1
      rewind iunit
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sdiagparams2(iunit,n,mts,modesx,modesy,nrec,fname)
! return parameters for selected scalar diagnostic
      use in2, only: pot2d, ntp, fpname, modesxp, modesyp, nprec,       &
     &dene2d, ntde, fdename, modesxde, modesyde, nderec,                &
     &deni2d, ntdi, fdiname, modesxdi, modesydi, ndirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2,3)
! mts = number of time steps between selected diagnostic
! modesx/modesy = number of modes in selected diagnostic in x/y
! nrec = number of records in selected diagnostic
! fname = file name for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, modesx, modesy, nrec
      character(len=*), intent(inout) :: fname
! local data
      select case(n)
! load metadata for potential data
      case (1)
         read (iunit,pot2d)
         mts = ntp
         modesx = modesxp; modesy = modesxp
         fname = fpname; nrec = nprec
! load metadata for electron density data
      case (2)
         read (iunit,dene2d)
         mts = ntde
         modesx = modesxde; modesy = modesyde
         fname = fdename; nrec = nderec
! load metadata for ion density data
      case (3)
         read (iunit,deni2d)
         mts = ntdi
         modesx = modesxdi; modesy = modesydi
         fname = fdiname; nrec = ndirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine vdiagparams2(iunit,n,mts,modesx,modesy,nrec,fname)
! return parameters for selected vector diagnostic
      use in2, only: el2d, ntel, felname, modesxel, modesyel, nelrec,   &
     &vcure2d, ntje, fjename, modesxje, modesyje, njerec,               &
     &vpot2d, nta, faname, modesxa, modesya, narec,                     &
     &et2d, ntet, fetname, modesxet, modesyet, netrec,                  &
     &b2d, ntb, fbname, modesxb, modesyb, nbrec,                        &
     &vpotr2d, ntar, farname, modesxar, modesyar, narrec,               &
     &vcuri2d, ntji, fjiname, modesxji, modesyji, njirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2,3,4,5,6,7)
! mts = numbe of time steps between selected diagnostic
! modesx/modesy = number of modes in selected diagnostic in x/y
! nrec = number of records in selected diagnostic
! fname = file name for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, modesx, modesy, nrec
      character(len=*), intent(inout) :: fname
! local data
      select case(n)
! load metadata for longitudinal efield data
      case (1)
         read (iunit,el2d)
         mts = ntel
         modesx = modesxel; modesy = modesyel
         fname = felname; nrec = nelrec
! load metadata for electron current density data
      case (2)
         read (iunit,vcure2d)
         mts = ntje
         modesx = modesxje; modesy = modesyje
         fname = fjename; nrec = njerec
! load metadata for vector potential data
      case (3)
         read (iunit,vpot2d)
         mts = nta
         modesx = modesxa; modesy = modesya
         fname = faname; nrec = narec
! load metadata for transverse efield data
      case (4)
         read (iunit,et2d)
         mts = ntet
         modesx = modesxet; modesy = modesyet
         fname = fetname; nrec = netrec
! load metadata for magnetic field data
      case (5)
         read (iunit,b2d)
         mts = ntb
         modesx = modesxb; modesy =  modesyb
         fname = fbname; nrec = nbrec
! load metadata for radiative vector potential data
      case (6)
         read (iunit,vpotr2d)
         mts = ntar
         modesx = modesxar; modesy = modesyar
         fname = farname; nrec = narrec
! load metadata for ion current density data
      case (7)
         read (iunit,vcuri2d)
         mts = ntji
         modesx = modesxji; modesy = modesyji
         fname = fjiname; nrec = njirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fldiagparams2(iunit,n,mts,mpro,mprd,nrec,fname)
! return parameters for selected fluid diagnostic
      use in2, only: fm2d, ntfm, npro, nprd, ffename, ffiname, nferec,  &
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
      select case(n)
! load metadata for electron fluid moments data
      case (1)
         read (iunit,fm2d)
         mts = ntfm
         mpro = npro; mprd = nprd
         fname = ffename; nrec = nferec
! load metadata for ion fluid moments data
      case (2)
         read (iunit,fm2d)
         mts = ntfm
         mpro = npro; mprd = nprd
         fname = ffiname; nrec = nfirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fvdiagparams2(iunit,n,mts,mmv,mfvd,mfed,nrec,fname)
! return parameters for selected velocity diagnostic
      use in2, only: fv2d, ntv, nmv, nfvd, nfed, fvename, fviname,&
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
         read (iunit,fv2d)
         mts = ntv; mmv = nmv
         mfvd = nfvd; mfed = nfed
         fname = fvename; nrec = nverec
! load metadata for ion velocity-space data
      case (2)
         read (iunit,fv2d)
         mts = ntv; mmv = nmv
         mfvd = nfvd; mfed = nfed
         fname = fviname; nrec = nvirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine trdiagparams2(iunit,n,mts,mdt,mst,mmv,mdimp,mprobt,nrec&
     &,fname)
! return parameters for selected trajectory diagnostic
      use in2, only: tr2d, ntt, ndt, nst, nmv, ndimp, nprobt, ftname,   &
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
         read (iunit,tr2d)
         mts = ntt; mdt = ndt
         mst = nst; mmv = nmv
         mdimp = ndimp; mprobt = nprobt
         fname = ftname; nrec = ntrec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine psdiagparams2(iunit,n,mts,mmv,msxb,msyb,nrec,fname)
! return parameters for selected phase space diagnostic
      use in2, only: ps2d, nts, nmv, nsxb, nsyb, fsename, fsiname,      &
     &nserec, nsirec
      implicit none
! iunit = fortran unit number to be used 
! n = diagnostic selected (1,2)
! mts = number of time steps between selected diagnostic
! mmv = number of segments in v for velocity distribution
! mxb/msyb = number of segments in x/y for velocity distribution
! nrec = number of records in selected diagnostic
! fname = file names for selected diagnostic
      integer, intent(in) :: iunit, n
      integer, intent(inout) :: mts, mmv, msxb, msyb, nrec
      character(len=*), intent(inout) :: fname
! local data
      select case(n)
! load metadata for electron phase-space data
      case (1)
         read (iunit,ps2d)
         mts = nts; mmv = nmv
         msxb = nsxb; msyb = nsyb
         fname = fsename; nrec = nserec
! load metadata for ion phase-space data
      case (2)
         read (iunit,ps2d)
         mts = nts; mmv = nmv
         msxb = nsxb; msyb = nsyb
         fname = fsiname; nrec = nsirec
      end select
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadc2(iunit,fc,nx,ny)
! read complex scalar field
      implicit none
      integer, intent(in) :: iunit, nx, ny
      complex, dimension(:,:), intent(inout) :: fc
! local data
      integer :: j, k
      read (unit=iunit) ((fc(j,k),j=1,nx),k=1,ny)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadvc2(iunit,fvc,ndim,nx,ny)
! read complex vector field
      implicit none
      integer, intent(in) :: iunit, ndim, nx, ny
      complex, dimension(:,:,:), intent(inout) :: fvc
! local data
      integer :: i, j, k
      read (unit=iunit) (((fvc(i,j,k),i=1,ndim),j=1,nx),k=1,ny)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fread2(iunit,f,nx,ny)
! read real scalar field
      implicit none
      integer, intent(in) :: iunit, nx, ny
      real, dimension(:,:), intent(inout) :: f
! local data
      integer :: j, k
      read (unit=iunit) ((f(j,k),j=1,nx),k=1,ny)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadv2(iunit,fv,ndim,nx,ny)
! read real vector field
      implicit none
      integer, intent(in) :: iunit, ndim, nx, ny
      real, dimension(:,:,:), intent(inout) :: fv
! local data
      integer :: i, j, k
      read (unit=iunit) (((fv(i,j,k),i=1,ndim),j=1,nx),k=1,ny)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadfv2(iunit,fvm,fv,fe,wk,ndim,nmvf,nfvd,nfed)
! read real velocity record from binary file
! fvm,fv,fe = real velocity data to be read in each record
! wk = total energy contained in distribution
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit, ndim, nmvf, nfvd, nfed
      real, intent(inout) :: wk
      real, dimension(:,:), intent(inout) :: fvm, fv, fe
! local data
      integer :: j, k
      read (unit=iunit) ((fvm(j,k),j=1,ndim),k=1,3),                    &
     &((fv(j,k),j=1,nmvf),k=1,nfvd), ((fe(j,k),j=1,nmvf),k=1,nfed), wk
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine freadtr2(iunit,partt,idimp,nprobt)
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
      subroutine freadps2(iunit,fvs,ndim,nmvf,nsxb,nsyb)
! read real velocity record from binary file
! fvs = real velocity data to be read in each record
! iunit = fortran unit number to be used 
      implicit none
      integer, intent(in) :: iunit, ndim, nmvf, nsxb, nsyb
      real, dimension(:,:,:,:), intent(inout) :: fvs
! local data
      integer :: i, j, k, l
      read (unit=iunit) ((((fvs(i,j,k,l),i=1,nmvf),j=1,ndim),k=1,nsxb), &
     &l=1,nsyb)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine rewindff2(iunit)
! rewind Fortran file with unit number iunit
      implicit none
      integer, intent(in) :: iunit
      rewind iunit
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine closeff2(iunit)
! close Fortran file with unit number iunit
      implicit none
      integer, intent(in) :: iunit
      close(unit=iunit)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwrmodes2(pot,pott,nx,ny,modesx,modesy)
! reads and copies lowest order scalar modes to packed array, writing
! zeroes to high order modes
      implicit none
      integer, intent(in) :: nx, ny, modesx, modesy
      real, dimension(:,:), intent(inout) :: pot
      complex, dimension(:,:), intent(in) :: pott
! local data
      integer :: nxvh, nyv, modesxd, modesyd
! extract dimensions
      nxvh = size(pot,1)/2; nyv = size(pot,2)
      modesxd = size(pott,1); modesyd = size(pott,2)
! call low level procedure
      call WRMODES2(pot,pott,nx,ny,modesx,modesy,nxvh,nyv,modesxd,      &
     &modesyd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwrvmodes2(vpot,vpott,nx,ny,modesx,modesy)
! reads and copies lowest order vector modes to packed array, writing
! zeroes to high order modes
      implicit none
      integer, intent(in) :: nx, ny, modesx, modesy
      real, dimension(:,:,:), intent(inout) :: vpot
      complex, dimension(:,:,:), intent(in) :: vpott
! local data
      integer :: ndim, nxvh, nyv, modesxd, modesyd
! extract dimensions
      ndim = size(vpot,1); nxvh = size(vpot,2)/2; nyv = size(vpot,3)
      modesxd = size(vpott,2); modesyd = size(vpott,3)
! call low level procedure
      call WRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,nxvh,nyv,      &
     &modesxd,modesyd)
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
      subroutine mfft2_init(mixup,sct,indx,indy)
! calculates tables needed by 2d FFTs
      implicit none
      integer, intent(in) :: indx, indy
      integer, dimension(:), intent(inout) :: mixup
      complex, dimension(:), intent(inout) :: sct
! local data
      integer :: nxhyd, nxyhd
! extract dimensions
      nxhyd = size(mixup,1); nxyhd = size(sct,1)
! call low level procedure
      call WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft2r(f,isign,mixup,sct,indx,indy)
! wrapper function for scalar 2d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy
      real, dimension(:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
! local data
      integer :: nxhd, nyd, nxhyd, nxyhd
! extract dimensions
      nxhd = size(f,1)/2; nyd = size(f,2)
      nxhyd = size(mixup,1); nxyhd = size(sct,1)
! call low level procedure
      call WFFT2RMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft2rn(f,isign,mixup,sct,indx,indy)
! wrapper function for vector 2d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy
      real, dimension(:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
! local data
      integer :: ndim, nxhd, nyd, nxhyd, nxyhd
! extract dimensions
      ndim = size(f,1); nxhd = size(f,2)/2; nyd = size(f,3)
      nxhyd = size(mixup,1); nxyhd = size(sct,1)
! call low level procedure
      select case(ndim)
      case (2)
         call WFFT2RM2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd)
      case (3)
         call WFFT2RM3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd)
      case default
         write (*,*) 'invalid dimension: ndim = ', ndim
      end select
      end subroutine
!
      end module
