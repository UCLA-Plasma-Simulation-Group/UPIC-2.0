!-----------------------------------------------------------------------
!
      module modmpdfield2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpdfield2.f
! mppoisd2_init calculates table needed by 2d dirichlet poisson solver
!               calls MPPOISD22
! mppoisd2 solves 2d or 2-1/2d dirichlet poisson's equation for smoothed
!          electric field
!          calls MPPOISD22 or MPPOISD23
! mpcuperpd2 calculates the transverse current in fourier space with
!            dirichlet boundary conditions
!            calls MPPCUPERPD2
! mpibpoisd2 solves 2-1/2d dirichlet poisson's equation in fourier space
!            for unsmoothed magnetic field
!            calls MIPPBPOISD23
! mpmaxweld2 solves 2-1/2d dirichlet maxwell's equation for unsmoothed
!            transverse electric and magnetic fields
!            calls MPPMAXWELD2
! mpemfieldr2 adds and smooths or copies and smooths real vector
!             fields in fourier space
!             calls MPPEMFIELDR2
! mpbbpoisd2 solves 2-1/2d dirichlet poisson's equation in fourier space
!            for smoothed magnetic field
!            calls MPPBBPOISD23
! mpdcuperpd2 calculates transverse part of the derivative of the
!             current density in fourier space from the momentum flux
!             for dirichlet boundary conditions
!             calls MPPDCUPERPD23
! mpadcuperpd2 calculates transverse part of the derivative of the
!              current density in fourier space from the momentum flux
!              and acceleration density for dirichlet boundary conditions
!              calls MPPADCUPERPD23
! mpepoisd2_init calculates table needed by 2-1/2d dirichlet poisson
!                solver for transverse electric field
!                calls MPPEPOISD23
! mpepoisd2 solves 2-1/2d dirichlet poisson's equation for smoothed
!           transverse electric field
!           calls MPPEPOISD23
! mppotd2 solves 2d dirichlet poisson's equation for potential
!         calls MPPOTPD2
! mpelfieldd2 solves 2d or 2-1/2d dirichlet poisson's equation for
!             unsmoothed
!             longitudinal electric field
!             calls MPPELFIELDD22 or MPPELFIELDD23
! mpdivfd2 calculates the divergence of periodic 2d vector field
!          calls MPPDIVFD2
! mpgradfd2 calculates the gradient of periodic 2d scalar field
!           calls MPPGRADFD2
! mpcurlfd2 calculates the curl of periodic 2-1/2d vector field
!           calls MPPCURLFD2
! mpavpotd2 calculates 2-1/2d vector potential from magnetic field
!           calls MPPAVPOTD23
! mpapotd2 solves 2-1/2d dirichlet poisson's equation for vector
!          potential
!          calls MPPAPOTD23
! mpetfieldd2 solves 2-1/2d dirichlet poisson's equation for unsmoothed
!             transverse electric field
!             calls MPPETFIELDD23
! mpsmoothd2 provides a 2d scalar dirichlet smoothing function
!           calls MPPSMOOTHD2
! mpsmoothd23 provides a 2d vector dirichlet smoothing function
!             calls MPPSMOOTHD23
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: april 27, 2017
!
      use libmpdfield2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mppoisd2_init(ffd,ax,ay,affp,nx,ny,kstrt)
! calculates table needed by 2d dirichlet  poisson solver
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ax, ay, affp
      complex, dimension(:,:), intent(inout) :: ffd
! local data
      integer :: isign = 0
      integer :: nyv, nyd, kxp2
      real :: we
      real, dimension(1,1)  :: q
      real, dimension(2,1,1) :: fxy
      nyv = size(q,1)
      nyd = size(ffd,1); kxp2 = size(ffd,2)
      call MPPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp2,&
     &nyd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppoisd2(q,fxy,ffd,we,tfield,nx,ny,kstrt)
! solves 2d or 2-1/2d dirichlet poisson's equation for smoothed electric
! field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: we, tfield
      real, dimension(:,:), intent(in)  :: q
      real, dimension(:,:,:), intent(inout) :: fxy
      complex, dimension(:,:), intent(inout) :: ffd
! local data
      integer :: isign = -1
      integer :: nyv, kxp2, ndim, nyd
      real :: ax, ay, affp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(q,1); kxp2 = size(q,2) - 1
      ndim = size(fxy,1); nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call MPPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,  &
     &kxp2,nyd)
      case (3)
         call MPPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,  &
     &kxp2,nyd)
      case default
         write (*,*) 'mppoisd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcuperpd2(cu,tfield,nx,ny,kstrt)
! calculates the transverse current in fourier space for dirichlet
! boundary conditions
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(inout) :: cu
! local data
      integer :: ndim, nyv, kxp2
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nyv = size(cu,2); kxp2 = size(cu,3) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPCUPERPD2(cu,nx,ny,kstrt,nyv,kxp2)
      case default
         write (*,*) 'mpcuperpd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpibpoisd2(cu,bxy,ffd,ci,wm,tfield,nx,ny,kstrt)
! solves 2-1/2d dirichlet poisson's equation in fourier space for
! unsmoothed magnetic field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      real, dimension(:,:,:), intent(in) :: cu
      real, dimension(:,:,:), intent(inout) :: bxy
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: ndim, nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nyv = size(cu,2); kxp2 = size(cu,3) - 1
      nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MIPPBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,nyd)
      case default
         write (*,*) 'mpibpoisd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpmaxweld2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,tfield,nx,ny&
     &,kstrt)
! solves 2-1/2d dirichlet maxwell's equation for unsmoothed transverse
! electric and magnetic fields 
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: affp, ci, dt
      real, intent(inout) :: wf, wm, tfield
      real, dimension(:,:,:), intent(inout) :: exy, bxy
      real, dimension(:,:,:), intent(in)  :: cu
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: ndim, nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nyv = size(cu,2); kxp2 = size(cu,3) - 1
      nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPMAXWELD2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstrt,  &
     &nyv,kxp2,nyd)
      case default
         write (*,*) 'mpmaxweld2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpemfieldr2(fxy,exy,ffd,isign,tfield,nx,ny,kstrt)
! adds and smooths or copies and smooths real vector fields in
! fourier space
      implicit none
      integer, intent(in) :: isign, nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(inout) :: fxy
      real, dimension(:,:,:), intent(in) :: exy
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: ndim, nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(fxy,1); nyv = size(fxy,2); kxp2 = size(fxy,3) - 1
      nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPEMFIELDR2(fxy,exy,ffd,isign,nx,ny,kstrt,nyv,kxp2,nyd)
      case default
         write (*,*) 'mpemfieldr2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpbbpoisd2(cu,bxy,ffd,ci,wm,tfield,nx,ny,kstrt)
! solves 2-1/2d dirichlet poisson's equation in fourier space for
! smoothed magnetic field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      real, dimension(:,:,:), intent(in) :: cu
      real, dimension(:,:,:), intent(inout) :: bxy
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: ndim, nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nyv = size(cu,2); kxp2 = size(cu,3) - 1
      nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPBBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,nyd)
      case default
         write (*,*) 'mpbbpoisd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdcuperpd2(dcu,amu,tfield,nx,ny,kstrt)
! calculates transverse part of the derivative of the current density
! in fourier space from the momentum flux
! for dirichlet boundary conditions
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:), intent(in) :: amu
! local data
      integer :: ndim, nyv, kxp2
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(dcu,1); nyv = size(dcu,2); kxp2 = size(dcu,3) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPDCUPERPD23(dcu,amu,nx,ny,kstrt,nyv,kxp2)
      case default
         write (*,*) 'mpdcuperpd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpadcuperpd2(dcu,amu,tfield,nx,ny,kstrt)
! calculates transverse part of the derivative of the current density
! in fourier space from the momentum flux and acceleration density for
! dirichlet boundary conditions
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(inout) :: dcu
      real, dimension(:,:,:), intent(in) :: amu
! local data
      integer :: ndim, nyv, kxp2
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(dcu,1); nyv = size(dcu,2); kxp2 = size(dcu,3) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPADCUPERPD23(dcu,amu,nx,ny,kstrt,nyv,kxp2)
      case default
         write (*,*) 'mpadcuperpd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpepoisd2_init(ffe,ax,ay,affp,wp0,ci,nx,ny,kstrt)
! calculates table needed by 2-1/2d dirichlet poisson solver for
! transverse electric field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ax, ay, affp, wp0, ci
      complex, dimension(:,:), intent(inout) :: ffe
! local data
      integer :: isign = 0
      integer :: nyv, nyd, kxp2
      real :: wf
      real, dimension(3,1,1)  :: dcu, exy
      nyv = size(dcu,2)
      nyd = size(ffe,1); kxp2 = size(ffe,2)
      call MPPEPOISD23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,    &
     &kstrt,nyv,kxp2,nyd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpepoisd2(dcu,exy,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
! solves 2-1/2d dirichlet poisson's equation for smoothed transverse
! electric field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: affp, ci
      real, intent(inout) :: wf, tfield
      real, dimension(:,:,:), intent(in) :: dcu
      real, dimension(:,:,:), intent(inout) :: exy
      complex, dimension(:,:), intent(inout) :: ffe
! local data
      integer :: isign = -1
      integer :: ndim, nyv, kxp2, nyd
      real :: ax, ay, wp0
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(dcu,1); nyv = size(dcu,2); kxp2 = size(dcu,3) - 1
      nyd = size(ffe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPEPOISD23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny, &
     &kstrt,nyv,kxp2,nyd)
      case default
         write (*,*) 'mpepoisd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppotd2(q,pot,ffd,we,tfield,nx,ny,kstrt)
! solves 2d dirichlet poisson's equation for potential
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: we, tfield
      real, dimension(:,:), intent(in) :: q
      real, dimension(:,:), intent(inout) :: pot
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(q,1); kxp2 = size(q,2) - 1
      nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPOTPD2(q,pot,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpelfieldd2(q,fxy,ffd,we,tfield,nx,ny,kstrt)
! solves 2d or 2-1/2d dirichlet poisson's equation for unsmoothed longitudinal
! electric field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: we, tfield
      real, dimension(:,:), intent(in)  :: q
      real, dimension(:,:,:), intent(inout) :: fxy
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: nyv, kxp2, ndim, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(q,1); kxp2 = size(q,2) - 1
      ndim = size(fxy,1); nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call MPPELFIELDD22(q,fxy,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
      case (3)
         call MPPELFIELDD23(q,fxy,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
      case default
         write (*,*) 'mpelfieldd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdivfd2(f,df,tfield,nx,ny,kstrt)
! calculates the divergence of periodic 2d vector field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(in) :: f
      real, dimension(:,:), intent(inout) :: df
! local data
      integer :: ndim, nyv, kxp2
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nyv = size(f,2); kxp2 = size(f,3) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2:)
         call MPPDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2)
      case default
         write (*,*) 'mpgradfd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgradfd2(df,f,tfield,nx,ny,kstrt)
! calculates the gradient of periodic 2d scalar field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(in) :: df
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nyv, kxp2
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nyv = size(f,2); kxp2 = size(f,3) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2:)
         call MPPGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2)
      case default
         write (*,*) 'mpgradfd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcurlfd2(f,g,tfield,nx,ny,kstrt)
! calculates the curl of periodic 2-1/2d vector field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(in) :: f
      real, dimension(:,:,:), intent(inout) :: g
! local data
      integer :: ndim, nyv, kxp2
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nyv = size(f,2); kxp2 = size(f,3) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2)
      case default
         write (*,*) 'mpcurlfd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpavpotd2(bxy,axy,tfield,nx,ny,kstrt)
! calculates 2-1/2d vector potential from magnetic field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(in) :: bxy
      real, dimension(:,:,:), intent(inout) :: axy
! local data
      integer :: ndim, nyv, kxp2
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(bxy,1); nyv = size(bxy,2); kxp2 = size(bxy,3) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPAVPOTD23(bxy,axy,nx,ny,kstrt,nyv,kxp2)
      case default
         write (*,*) 'mpavpotd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpapotd2(cu,axy,ffd,ci,wm,tfield,nx,ny,kstrt)
! solves 2-1/2d dirichlet poisson's equation in fourier space for
! vector potential
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      real, dimension(:,:,:), intent(in) :: cu
      real, dimension(:,:,:), intent(inout) :: axy
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: ndim, nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nyv = size(cu,2); kxp2 = size(cu,3) - 1
      nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPAPOTD23(cu,axy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,nyd)
      case default
         write (*,*) 'mpapotd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpetfieldd2(dcu,exy,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
! solves 2-1/2d dirichlet poisson's equation for unsmoothed transverse
! electric field
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(in) :: affp, ci
      real, intent(inout) :: wf, tfield
      real, dimension(:,:,:), intent(in) :: dcu
      real, dimension(:,:,:), intent(inout) :: exy
      complex, dimension(:,:), intent(in) :: ffe
! local data
      integer :: ndim, nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(dcu,1); nyv = size(dcu,2); kxp2 = size(dcu,3) - 1
      nyd = size(ffe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPETFIELDD23(dcu,exy,ffe,affp,ci,wf,nx,ny,kstrt,nyv,kxp2,&
     &nyd)
      case default
         write (*,*) 'mpetfieldd2: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpsmoothd2(q,qs,ffd,tfield,nx,ny,kstrt)
! provides a 2d dirichlet scalar smoothing function
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(in) :: q
      real, dimension(:,:), intent(inout) :: qs
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nyv = size(q,1); kxp2 = size(q,2) - 1
      nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPSMOOTHD2(q,qs,ffd,nx,ny,kstrt,nyv,kxp2,nyd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
      subroutine mpsmoothd23(cu,cus,ffd,tfield,nx,ny,kstrt)
! provides a 2d vector smoothing function
      implicit none
      integer, intent(in) :: nx, ny, kstrt
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(in) :: cu
      real, dimension(:,:,:), intent(inout) :: cus
      complex, dimension(:,:), intent(in) :: ffd
! local data
      integer :: ndim, nyv, kxp2, nyd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nyv = size(cu,2); kxp2 = size(cu,3) - 1
      nyd = size(ffd,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call MPPSMOOTHD23(cu,cus,ffd,nx,ny,kstrt,nyv,kxp2,nyd)
      case default
         write (*,*) 'mpsmoothd23: unsupported dimension ndim = ', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
      end module
