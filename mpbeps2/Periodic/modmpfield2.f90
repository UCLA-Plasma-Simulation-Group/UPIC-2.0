!-----------------------------------------------------------------------
!
      module modmpfield2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpfield2.f
! mppois2_init calculates table needed by 2d poisson solver
!              calls MPPOIS22
! mppois2 solves 2d or 2-1/2d poisson's equation for smoothed electric
!         field
!         calls MPPOIS22 or MPPOIS23
! mpaddqei2 adds electron and ion densities
!           calls MPADDQEI2
! mpcuperp2 calculates the transverse current in fourier space
!           calls MPPCUPERP2
! mpibpois2 solves 2-1/2d poisson's equation for unsmoothed magnetic
!           field
!           calls MIPPBPOISP23
! mpmaxwel2 solves 2-1/2d maxwell's equation for unsmoothed transverse
!           electric and magnetic fields
!           calls MPPMAXWEL2
! mpemfield2 adds and smooths or copies and smooths complex vector
!            fields in fourier space
!            calls MPPEMFIELD2
! mpaddcuei2 adds electron and ion current densities
!            calls MPADDCUEI23
! mpaddamui2 adds electron and ion momentum flux densities
!            calls MPADDAMUI23
! mpbaddext2 adds constant to magnetic field for 2-1/2d code
!            calls PPBADDEXT2
! mpaddvrfield2 calculates a = b + c
!               calls PPADDVRFIELD2
! mpbbpois2 solves 2-1/2d poisson's equation in fourier space for
!           smoothed magnetic field
!           calls MPPBBPOISP23
! mpdcuperp2 calculates transverse part of the derivative of the current
!            density from the momentum flux
!            calls MPPDCUPERP23
! mpadcuperp2 calculates transverse part of the derivative of the
!             current density from the momentum flux and acceleration
!             density
!             calls MPPADCUPERP23
! mpepois2_init calculates table needed by 2-1/2d poisson solver for
!               transverse electric field
!               calls MPPEPOISP23
! mpepois2 solves 2-1/2d poisson's equation for smoothed transverse
!          electric field
!          calls MPPEPOISP23
! mppot2 solves 2d poisson's equation for potential
!        calls MPPOTP2
! mpelfield2 solves 2d or 2-1/2d poisson's equation for unsmoothed
!            longitudinal electric field
!            calls MPPELFIELD22 or MPPELFIELD23
! mpdivf2 calculates the divergence in fourier space
!         calls MPPDIVF2
! mpgradf2 calculates the gradient in fourier space
!          calls MPPGRADF2
! mpcurlf2 calculates the curl in fourier space
!          calls MPPCURLF2
! mpavpot2 calculates 2-1/2d vector potential from magnetic field
!          calls MPPAVPOT23
! mpavrpot2 solves 2-1/2d poisson's equation for the radiative part of
!           the vector potential
!           calls MPPAVRPOT23
! mpapot2 solves 2-1/2d poisson's equation for vector potential
!         calls MPPAPOTP23
! mpetfield2 solves 2-1/2d poisson's equation for unsmoothed transverse
!            electric field
!            calls MPPETFIELD23
! mpsmooth2 provides a 2d scalar smoothing function
!           calls MPPSMOOTH2
! mpsmooth23 provides a 2d vector smoothing function
!            calls MPPSMOOTH23
! mprdmodes2 extracts and stores lowest order scalar modes to unpacked
!            array
!            calls PPRDMODES2
! mpwrmodes2 reads and copies lowest order scalar modes to packed array,
!            writing zeroes to high order modes
!            calls PPWRMODES2
! mprdvmodes2 extracts and stores lowest order vector modes to unpacked
!             array
!             calls PPRDVMODES2
! mpwrvmodes2 reads and copies lowest order vector modes to packed array
!             writing zeroes to high order modes
!             calls PPWRVMODES2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: january 12, 2017
!
      use libmpfield2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mppois2_init(ffc,ax,ay,affp,nx,ny,kstrt)
! calculates table needed by 2d poisson solver
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ax, ay, affp
      complex, dimension(:,:), intent(inout) :: ffc
! local data
      integer :: isign = 0
      integer :: nyv, nyhd, kxp
      real :: we
      complex, dimension(1,1)  :: q
      complex, dimension(2,1,1) :: fxy
      nyv = size(q,1)
      nyhd = size(ffc,1); kxp = size(ffc,2)
      call MPPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp,  &
     &nyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppois2(q,fxy,ffc,we,tfield,nx,ny,kstrt)
! solves 2d or 2-1/2d poisson's equation for smoothed electric field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: we, tfield
      complex, dimension(:,:), intent(in)  :: q
      complex, dimension(:,:,:), intent(inout) :: fxy
      complex, dimension(:,:), intent(inout) :: ffc
! local data
      integer :: isign = -1
      integer :: nyv, kxp, ndim, nyhd
      real :: ax, ay, affp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(q,1); kxp = size(q,2)
      ndim = size(fxy,1); nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call MPPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp&
     &,nyhd)
      case (3)
         call MPPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp&
     &,nyhd)
      case default
         write (*,*) 'mppois2: unsupported dimension ndim = ' , ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaddqei2(qe,qi,nyp,tfield,nx)
! adds electron and ion densities
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(inout) :: qe
      real, dimension(:,:), intent(in) :: qi
! local data
      integer :: nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(qe,1); nypmx = size(qe,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPADDQEI2(qe,qi,nyp,nx,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcuperp2(cu,tfield,nx,ny,kstrt)
! calculates the transverse current in fourier space
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(inout) :: cu
! local data
      integer :: nyv, kxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(cu,2); kxp = size(cu,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpibpois2(cu,bxy,ffc,ci,wm,tfield,nx,ny,kstrt)
! solves 2-1/2d poisson's equation for unsmoothed magnetic field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      complex, dimension(:,:,:), intent(in) :: cu
      complex, dimension(:,:,:), intent(inout) :: bxy
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(cu,2); kxp = size(cu,3)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MIPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpmaxwel2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,tfield,nx,ny,&
     &kstrt)
! solves 2-1/2d maxwell's equation for unsmoothed transverse electric
! and magnetic fields
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: affp, ci, dt
      real, intent(inout) :: wf, wm, tfield
      complex, dimension(:,:,:), intent(inout) :: exy, bxy
      complex, dimension(:,:,:), intent(in)  :: cu
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(cu,2); kxp = size(cu,3)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt,nyv,  &
     &kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpemfield2(fxy,exy,ffc,isign,tfield,nx,ny,kstrt)
! adds and smooths or copies and smooths complex vector fields in
! fourier space
      implicit none
      integer, intent(in) :: isign, nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(inout) :: fxy
      complex, dimension(:,:,:), intent(in) :: exy
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(fxy,2); kxp = size(fxy,3)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaddcuei2(cue,cui,nyp,tfield,nx)
! adds electron and ion current densities
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(inout) :: cue
      real, dimension(:,:,:), intent(in) :: cui
! local data
      integer :: nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(cue,2); nypmx = size(cue,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPADDCUEI23(cue,cui,nyp,nx,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaddamui2(amu,amui,nyp,tfield,nx)
! adds electron and ion momentum flux densities
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(inout) :: amu
      real, dimension(:,:,:), intent(in) :: amui
! local data
      integer :: nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(amu,2); nypmx = size(amu,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPADDAMUI23(amu,amui,nyp,nx,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpbaddext2(bxy,nyp,tfield,omx,omy,omz,nx)
! adds constant to magnetic field for 2-1/2d code
      implicit none
      integer, intent(in) :: nyp, nx
      real, intent(inout) :: tfield
      real, intent(in) :: omx, omy, omz
      real, dimension(:,:,:), intent(inout) :: bxy
! local data
      integer :: nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(bxy,2); nypmx = size(bxy,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPBADDEXT2(bxy,nyp,omx,omy,omz,nx,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaddvrfield2(a,b,c,tfield)
! calculates a = b + c
      implicit none
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(inout) :: a
      real, dimension(:,:,:), intent(in) :: b, c
! local data
      integer :: ndim, nxe, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(a,1); nxe = size(a,2); nypmx = size(a,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPADDVRFIELD2(a,b,c,ndim,nxe,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpbbpois2(cu,bxy,ffc,ci,wm,tfield,nx,ny,kstrt)
! solves 2-1/2d poisson's equation in fourier space for smoothed
! magnetic field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      complex, dimension(:,:,:), intent(in) :: cu
      complex, dimension(:,:,:), intent(inout) :: bxy
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(cu,2); kxp = size(cu,3)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPBBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdcuperp2(dcu,amu,tfield,nx,ny,kstrt)
! calculates transverse part of the derivative of the current density
! from the momentum flux
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(inout) :: dcu
      complex, dimension(:,:,:), intent(in) :: amu
! local data
      integer :: nyv, kxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(dcu,2); kxp = size(dcu,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpadcuperp2(dcu,amu,tfield,nx,ny,kstrt)
! calculates transverse part of the derivative of the current density
! from the momentum flux and acceleration density
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(inout) :: dcu
      complex, dimension(:,:,:), intent(in) :: amu
! local data
      integer :: nyv, kxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(dcu,2); kxp = size(dcu,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPADCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpepois2_init(ffe,ax,ay,affp,wp0,ci,nx,ny,kstrt)
! calculates table needed by 2-1/2d poisson solver for transverse
! electric field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ax, ay, affp, wp0, ci
      complex, dimension(:,:), intent(inout) :: ffe
! local data
      integer :: isign = 0
      integer :: nyv, nyhd, kxp
      real :: wf
      complex, dimension(3,1,1)  :: dcu
      complex, dimension(3,1,1) :: exy
      nyv = size(dcu,2)
      nyhd = size(ffe,1); kxp = size(ffe,2)
      call MPPEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,    &
     &kstrt,nyv,kxp,nyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpepois2(dcu,exy,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
! solves 2-1/2d poisson's equation for smoothed transverse electric
! field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: affp, ci
      real, intent(inout) :: wf, tfield
      complex, dimension(:,:,:), intent(in) :: dcu
      complex, dimension(:,:,:), intent(inout) :: exy
      complex, dimension(:,:), intent(inout) :: ffe
! local data
      integer :: isign = -1
      integer :: nyv, kxp, nyhd
      real :: ax, ay, wp0
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(dcu,2); kxp = size(dcu,3)
      nyhd = size(ffe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,    &
     &kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppot2(q,pot,ffc,we,tfield,nx,ny,kstrt)
! solves 2d poisson's equation for potential
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: we, tfield
      complex, dimension(:,:), intent(in) :: q
      complex, dimension(:,:), intent(inout) :: pot
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(q,1); kxp = size(q,2)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPOTP2(q,pot,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpelfield2(q,fxy,ffc,we,tfield,nx,ny,kstrt)
! solves 2d or 2-1/2d poisson's equation for unsmoothed longitudinal
! electric field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: we, tfield
      complex, dimension(:,:), intent(in)  :: q
      complex, dimension(:,:,:), intent(inout) :: fxy
      complex, dimension(:,:), intent(inout) :: ffc
! local data
      integer :: nyv, kxp, ndim, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(q,1); kxp = size(q,2)
      ndim = size(fxy,1); nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call MPPELFIELD22(q,fxy,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
      case (3)
         call MPPELFIELD23(q,fxy,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
      case default
         write (*,*) 'mpelfield2: unsupported dimension ndim = ' , ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdivf2(f,df,tfield,nx,ny,kstrt)
! calculates the divergence in fourier space
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(in) :: f
      complex, dimension(:,:), intent(inout) :: df
! local data
      integer :: ndim, nyv, kxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nyv = size(f,2); kxp = size(f,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgradf2(df,f,tfield,nx,ny,kstrt)
! calculates the gradient in fourier space
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(in) :: df
      complex, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nyv, kxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(df,1); kxp = size(df,2)
      ndim = size(f,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcurlf2(f,g,tfield,nx,ny,kstrt)
! calculates the curl in fourier space
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(in) :: f
      complex, dimension(:,:,:), intent(inout) :: g
! local data
      integer :: nyv, kxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(f,2); kxp = size(f,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPCURLF2(f,g,nx,ny,kstrt,nyv,kxp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpavpot2(bxy,axy,tfield,nx,ny,kstrt)
! calculates 2-1/2d vector potential from magnetic field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(in) :: bxy
      complex, dimension(:,:,:), intent(inout) :: axy
! local data
      integer :: nyv, kxp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(bxy,2); kxp = size(bxy,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpavrpot2(axy,bxy,ffc,affp,ci,tfield,nx,ny,kstrt)
! solves 2-1/2d poisson's equation for the radiative part of the vector
! potential
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: affp, ci
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(inout) :: axy
      complex, dimension(:,:,:), intent(in) :: bxy
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(bxy,2); kxp = size(bxy,3)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpapot2(cu,axy,ffc,ci,wm,tfield,nx,ny,kstrt)
! solves 2-1/2d poisson's equation in fourier space for vector potential
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      complex, dimension(:,:,:), intent(in) :: cu
      complex, dimension(:,:,:), intent(inout) :: axy
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(cu,2); kxp = size(cu,3)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPAPOTP23(cu,axy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpetfield2(dcu,exy,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
! solves 2-1/2d poisson's equation for unsmoothed transverse electric
! field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: affp, ci
      real, intent(inout) :: wf, tfield
      complex, dimension(:,:,:), intent(in) :: dcu
      complex, dimension(:,:,:), intent(inout) :: exy
      complex, dimension(:,:), intent(in) :: ffe
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(dcu,2); kxp = size(dcu,3)
      nyhd = size(ffe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPETFIELD23(dcu,exy,ffe,affp,ci,wf,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpsmooth2(q,qs,ffc,tfield,nx,ny,kstrt)
! provides a 2d scalar smoothing function
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(in) :: q
      complex, dimension(:,:), intent(inout) :: qs
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(q,1); kxp = size(q,2)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPSMOOTH2(q,qs,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
      subroutine mpsmooth23(cu,cus,ffc,tfield,nx,ny,kstrt)
! provides a 2d vector smoothing function
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(in) :: cu
      complex, dimension(:,:,:), intent(inout) :: cus
      complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nyv, kxp, nyhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(cu,2); kxp = size(cu,3)
      nyhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPSMOOTH23(cu,cus,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdmodes2(pot,pott,tfield,nx,ny,modesx,modesy,kstrt)
! extracts and stores lowest order scalar modes to unpacked array
      implicit none
      integer, intent(in) :: nx, ny, modesx, modesy, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(in) :: pot
      complex, dimension(:,:), intent(inout) :: pott
! local data
      integer :: nyv, kxp, modesyd, modesxpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(pot,1); kxp = size(pot,2)
      modesyd = size(pott,1); modesxpd = size(pott,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRDMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,kxp,       &
     &modesxpd,modesyd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrmodes2(pot,pott,tfield,nx,ny,modesx,modesy,kstrt)
! reads and copies lowest order scalar modes to packed array, writing
! zeroes to high order modes
      implicit none
      integer, intent(in) :: nx, ny, modesx, modesy, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(inout) :: pot
      complex, dimension(:,:), intent(in) :: pott
! local data
      integer :: nyv, kxp, modesyd, modesxpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(pot,1); kxp = size(pot,2)
      modesyd = size(pott,1); modesxpd = size(pott,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,kxp,       &
     &modesxpd,modesyd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdvmodes2(vpot,vpott,tfield,nx,ny,modesx,modesy,kstrt&
     &)
! extracts and stores lowest order vector modes to unpacked array
      implicit none
      integer, intent(in) :: nx, ny, modesx, modesy, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(in) :: vpot
      complex, dimension(:,:,:), intent(inout) :: vpott
! local data
      integer :: ndim, nyv, kxp, modesyd, modesxpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(vpot,1); nyv = size(vpot,2); kxp = size(vpot,3)
      modesyd = size(vpott,2); modesxpd = size(vpott,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRDVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,kstrt,nyv,kxp&
     &,modesxpd,modesyd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrvmodes2(vpot,vpott,tfield,nx,ny,modesx,modesy,kstrt&
     &)
! extracts and stores lowest order vector modes to unpacked array
      implicit none
      integer, intent(in) :: nx, ny, modesx, modesy, kstrt
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(inout) :: vpot
      complex, dimension(:,:,:), intent(in) :: vpott
! local data
      integer :: ndim, nyv, kxp, modesyd, modesxpd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(vpot,1); nyv = size(vpot,2); kxp = size(vpot,3)
      modesyd = size(vpott,2); modesxpd = size(vpott,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,kstrt,nyv,kxp&
     &,modesxpd,modesyd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
      end module
