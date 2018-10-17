!-----------------------------------------------------------------------
!
      module mfield3
!
! Fortran90 wrappers to 3d MPI/OpenMP PIC library libmpfield3.f
! mppois3_init calculates table needed by 3d poisson solver
!              calls VMPPOIS332
! mppois3 solves 3d poisson's equation for smoothed electric field
!         calls VMPPOIS332
! mpaddqei3 adds electron and ion densities
!           calls MPADDQEI32
! mpcuperp3 calculates the transverse current in fourier space
!           calls MPPCUPERP32
! mpibpois3 solves 3d poisson's equation for unsmoothed magnetic field
!           calls VMIPPBPOISP332
! mpmaxwel3 solves 3d maxwell's equation for unsmoothed transvers
!           electric and magnetic fields
!           calls VMPPMAXWEL32
! mpemfield3 adds and smooths or copies and smooths complex vector
!            fields in fourier space
!            calls MPPEMFIELD32
! mpaddcuei3 adds electron and ion current densities
!            calls MPADDCUEI32
! mpaddamui3 adds electron and ion momentum flux densities
!            calls MPADDAMUI32
! mpbaddext3 adds constant to magnetic field for 3d code
!            calls MPPBADDEXT32
! mpaddvrfield3 calculates a = b + c
!               calls MPPADDVRFIELD32
! mpbbpois3 solves 3d poisson's equation in fourier space for smoothed
!           magnetic field
!           calls VMPPBBPOISP332
! mpdcuperp3 calculates transverse part of the derivative of the current
!            density from the momentum flux
!            calls MPPDCUPERP32
! mpadcuperp3 calculates transverse part of the derivative of the
!             current density from the momentum flux and acceleration
!             density
!             calls MPPADCUPERP32
! mpepois3_init calculates table needed by 3d poisson solver for
!               transverse electric field
!               calls VMPPEPOISP332
! mpepois3 solves 3d poisson's equation for smoothed transverse electric
!          field
!          calls MPPEPOISP332
! mppot3 solves 3d poisson's equation for potential
!        calls MPPOTP32
! mpelfield3 solves 3d poisson's equation for unsmoothed longitudinal
!            electric field
!            calls MPPELFIELD32
! mpdivf3 calculates the divergence in fourier space
!         calls MPPDIVF32
! mpgradf3 calculates the gradient in fourier space
!          calls MPPGRADF32
! mpcurlf3 calculates the curl in fourier space
!          calls MPPCURLF32
! mpavpot3 calculates 3d vector potential from magnetic field
!          calls MPPAVPOT332
! mcuave3 averages current in fourier space for 3d code
!         calls MCUAVE33
! mpavrpot3 solves 3d poisson's equation for the radiative part of the
!           vector potential
!           calls MPPAVRPOT332
! mpapot3 solves 3d poisson's equation for vector potential
!         calls MPPAPOTP32
! mpetfield3 solves 3d poisson's equation for unsmoothed transverse
!            electric field
!            calls VMPPETFIELD332
! mpsmooth3 provides a 3d scalar smoothing function
!           calls MPPSMOOTH32
! mpsmooth33 provides a 3d vector smoothing function
!            calls MPPSMOOTH332
! mprdmodes3 extracts and stores lowest order scalar modes to unpacked
!            array
!            calls PPRDMODES32
! mpwrmodes3 reads and copies lowest order scalar modes to packed array,
!            writing zeroes to high order modes
!            calls PPWRMODES32
! mprdvmodes3 extracts and stores lowest order vector modes to unpacked
!             array
!             calls PPRDVMODES32
! mpwrvmodes3 reads and copies lowest order vector modes to packed array
!             writing zeroes to high order modes
!             calls PPWRVMODES32
! mpset_pcvzero3 zeros out transverse field array.
!                calls SET_PCVZERO3
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: may 16, 2018
!
      use libmpfield3_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mppois3_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,nvpy,nvpz&
     &)
! calculates table needed by 3d poisson solver
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: ax, ay, az, affp
      complex, dimension(:,:,:), intent(inout) :: ffc
! local data
      integer :: isign = 0
      integer :: nzv, nzhd, kxyp, kyzp
      real :: we
      complex, dimension(1,1,1)  :: q
      complex, dimension(3,1,1,1) :: fxyz
      nzv = size(q,1)
      nzhd = size(ffc,1); kxyp = size(ffc,2); kyzp = size(ffc,3)
      call VMPPOIS332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt, &
     &nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppois3(q,fxyz,ffc,we,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! solves 3d poisson's equation for smoothed electric field
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: we, tfield
      complex, dimension(:,:,:), intent(in)  :: q
      complex, dimension(:,:,:,:), intent(inout) :: fxyz
      complex, dimension(:,:,:), intent(inout) :: ffc
! local data
      integer :: isign = -1
      integer :: nzv, kxyp, kyzp, nzhd
      real :: ax, ay, az, affp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(q,1); kxyp = size(q,2); kyzp = size(q,3)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call VMPPOIS332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt, &
     &nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaddqei3(qe,qi,nyzp,tfield,nx)
! adds electron and ion densities
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:,:), intent(inout) :: qe
      real, dimension(:,:,:), intent(in) :: qi
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(qe,1); nypmx = size(qe,2); nzpmx = size(qe,3)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPADDQEI32(qe,qi,nyzp,nx,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcuperp3(cu,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! calculates the transverse current in fourier space
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(inout) :: cu
! local data
      integer :: nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(cu,2); kxyp = size(cu,3); kyzp = size(cu,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPCUPERP32(cu,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpibpois3(cu,bxyz,ffc,ci,wm,tfield,nx,ny,nz,kstrt,nvpy,&
     &nvpz)
! solves 3d poisson's equation for unsmoothed magnetic field
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      complex, dimension(:,:,:,:), intent(in) :: cu
      complex, dimension(:,:,:,:), intent(inout) :: bxyz
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(cu,2); kxyp = size(cu,3); kyzp = size(cu,4)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call VMIPPBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,nvpz,nzv&
     &,kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpmaxwel3(exyz,bxyz,cu,ffc,affp,ci,dt,wf,wm,tfield,nx, &
     &ny,nz,kstrt,nvpy,nvpz)
! solves 3d maxwell's equation for unsmoothed transverse electric and
! magnetic fields
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: affp, ci, dt
      real, intent(inout) :: wf, wm, tfield
      complex, dimension(:,:,:,:), intent(inout) :: exyz, bxyz
      complex, dimension(:,:,:,:), intent(in) :: cu
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(cu,2); kxyp = size(cu,3); kyzp = size(cu,4)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call VMPPMAXWEL32(exyz,bxyz,cu,ffc,affp,ci,dt,wf,wm,nx,ny,nz,kstrt&
     &,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpemfield3(fxyz,exyz,ffc,isign,tfield,nx,ny,nz,kstrt,  &
     &nvpy,nvpz)
! adds and smooths or copies and smooths complex vector fields in
! fourier space
      implicit none
      integer, intent(in) :: isign, nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(inout) :: fxyz
      complex, dimension(:,:,:,:), intent(in) :: exyz
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(fxyz,2); kxyp = size(fxyz,3); kyzp = size(fxyz,4)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPEMFIELD32(fxyz,exyz,ffc,isign,nx,ny,nz,kstrt,nvpy,nvpz,nzv&
     &,kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaddcuei3(cue,cui,nyzp,tfield,nx)
! adds electron and ion current densities
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:,:,:), intent(inout) :: cue
      real, dimension(:,:,:,:), intent(in) :: cui
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(cue,2); nypmx = size(cue,3); nzpmx = size(cue,4)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPADDCUEI32(cue,cui,nyzp,nx,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaddamui3(amu,amui,nyzp,tfield,nx)
! adds electron and ion momentum flux densities
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:,:,:), intent(inout) :: amu
      real, dimension(:,:,:,:), intent(in) :: amui
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(amu,2); nypmx = size(amu,3); nzpmx = size(amu,4)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPADDAMUI32(amu,amui,nyzp,nx,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpbaddext3(bxyz,nyzp,tfield,omx,omy,omz,nx)
! adds constant to magnetic field for 3d code
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, intent(in) :: omx, omy, omz
      real, dimension(:,:,:,:), intent(inout) :: bxyz
      integer, dimension(:), intent(in) :: nyzp
! local data
      integer :: nxe, nypmx, nzpmx, idds
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(bxyz,2); nypmx = size(bxyz,3); nzpmx = size(bxyz,4)
      idds = size(nyzp,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPBADDEXT32(bxyz,nyzp,omx,omy,omz,nx,nxe,nypmx,nzpmx,idds)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpaddvrfield3(a,b,c,tfield)
! calculates a = b + c
      implicit none
      real, intent(inout) :: tfield
      real, dimension(:,:,:,:), intent(inout) :: a
      real, dimension(:,:,:,:), intent(in) :: b, c
! local data
      integer :: ndim, nxe, nypmx, nzpmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(a,1); nxe = size(a,2)
      nypmx = size(a,3); nzpmx = size(a,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPADDVRFIELD32(a,b,c,ndim,nxe,nypmx,nzpmx)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpbbpois3(cu,bxyz,ffc,ci,wm,tfield,nx,ny,nz,kstrt,nvpy,&
     &nvpz)
! solves 3d poisson's equation in fourier space for smoothed magnetic
! field
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      complex, dimension(:,:,:,:), intent(in) :: cu
      complex, dimension(:,:,:,:), intent(inout) :: bxyz
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(cu,2); kxyp = size(cu,3); kyzp = size(cu,4)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call VMPPBBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,nvpz,nzv&
     &,kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdcuperp3(dcu,amu,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! calculates transverse part of the derivative of the current density
! from the momentum flux
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(inout) :: dcu
      complex, dimension(:,:,:,:), intent(in) :: amu
! local data
      integer :: nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(dcu,2); kxyp = size(dcu,3); kyzp = size(dcu,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPDCUPERP32(dcu,amu,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpadcuperp3(dcu,amu,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! calculates transverse part of the derivative of the current density
! from the momentum flux and acceleration density
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(inout) :: dcu
      complex, dimension(:,:,:,:), intent(in) :: amu
! local data
      integer :: nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(dcu,2); kxyp = size(dcu,3); kyzp = size(dcu,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPADCUPERP32(dcu,amu,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpepois3_init(ffe,ax,ay,az,affp,wp0,ci,nx,ny,nz,kstrt, &
     &nvpy,nvpz)
! calculates table needed by 3d poisson solver for transverse electric
! field
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: ax, ay, az, affp, wp0, ci
      complex, dimension(:,:,:), intent(inout) :: ffe
! local data
      integer :: isign = 0
      integer :: nzv, nzhd, kxyp, kyzp
      real :: wf
      complex, dimension(3,1,1)  :: dcu
      complex, dimension(3,1,1) :: exyz
      nzv = size(dcu,2)
      nzhd = size(ffe,1); kxyp = size(ffe,2); kyzp = size(ffe,3)
      call VMPPEPOISP332(dcu,exyz,isign,ffe,ax,ay,az,affp,wp0,ci,wf,nx, &
     &ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpepois3(dcu,exyz,ffe,affp,ci,wf,tfield,nx,ny,nz,kstrt,&
     &nvpy,nvpz)
! solves 3d poisson's equation for smoothed transverse electric field
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: affp, ci
      real, intent(inout) :: wf, tfield
      complex, dimension(:,:,:,:), intent(in) :: dcu
      complex, dimension(:,:,:,:), intent(inout) :: exyz
      complex, dimension(:,:,:), intent(inout) :: ffe
! local data
      integer :: isign = -1
      integer :: nzv, kxyp, kyzp, nzhd
      real :: ax, ay, az, wp0
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(dcu,2); kxyp = size(dcu,3); kyzp = size(dcu,4)
      nzhd = size(ffe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call VMPPEPOISP332(dcu,exyz,isign,ffe,ax,ay,az,affp,wp0,ci,wf,nx, &
     &ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mppot3(q,pot,ffc,we,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! solves 3d poisson's equation for potential
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: we, tfield
      complex, dimension(:,:,:), intent(in) :: q
      complex, dimension(:,:,:), intent(inout) :: pot
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(q,1); kxyp = size(q,2); kyzp = size(q,3)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPOTP32(q,pot,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,&
     &nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpelfield3(q,fxyz,ffc,we,tfield,nx,ny,nz,kstrt,nvpy,   &
     &nvpz)
! solves 3d poisson's equation for unsmoothed longitudinal electric
! field
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: we, tfield
      complex, dimension(:,:,:), intent(in)  :: q
      complex, dimension(:,:,:,:), intent(inout) :: fxyz
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(q,1); kxyp = size(q,2); kyzp = size(q,3)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPELFIELD32(q,fxyz,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,&
     &kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdivf3(f,df,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! calculates the divergence in fourier space
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(in) :: f
      complex, dimension(:,:,:), intent(inout) :: df
! local data
      integer :: nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(f,2); kxyp = size(f,3); kyzp = size(f,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPDIVF32(f,df,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpgradf3(df,f,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! calculates the gradient in fourier space
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(in) :: df
      complex, dimension(:,:,:,:), intent(inout) :: f
! local data
      integer :: nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(df,1); kxyp = size(df,2); kyzp = size(df,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPGRADF32(df,f,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcurlf3(f,g,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! calculates the curl in fourier space
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(in) :: f
      complex, dimension(:,:,:,:), intent(inout) :: g
! local data
      integer :: nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(f,2); kxyp = size(f,3); kyzp = size(f,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPCURLF32(f,g,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpavpot3(bxyz,axyz,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! calculates 3d vector potential from magnetic field
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(in) :: bxyz
      complex, dimension(:,:,:,:), intent(inout) :: axyz
! local data
      integer :: nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(bxyz,2); kxyp = size(bxyz,3); kyzp = size(bxyz,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPAVPOT332(bxyz,axyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcuave3(cuave,cunew,cuold,tfield,nz)
! averages current in fourier space for 3d code
      implicit none
      integer, intent(in) :: nz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(in) :: cunew, cuold
      complex, dimension(:,:,:,:), intent(inout) :: cuave
! local data
      integer :: nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(cuave,2); kxyp = size(cuave,3); kyzp = size(cuave,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MCUAVE33(cuave,cunew,cuold,nz,kxyp,kyzp,nzv)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpavrpot3(axyz,bxyz,ffc,affp,ci,tfield,nx,ny,nz,kstrt, &
     &nvpy,nvpz)
! solves 3d poisson's equation for the radiative part of the vector
! potential
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: affp, ci
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(inout) :: axyz
      complex, dimension(:,:,:,:), intent(in) :: bxyz
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(bxyz,2); kxyp = size(bxyz,3); kyzp = size(bxyz,4)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPAVRPOT332(axyz,bxyz,ffc,affp,ci,nx,ny,nz,kstrt,nvpy,nvpz, &
     &nzv,kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpapot3(cu,axyz,ffc,ci,wm,tfield,nx,ny,nz,kstrt,nvpy,  &
     &nvpz)
! solves 3d poisson's equation in fourier space for vector potential
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      complex, dimension(:,:,:,:), intent(in) :: cu
      complex, dimension(:,:,:,:), intent(inout) :: axyz
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(cu,2); kxyp = size(cu,3); kyzp = size(cu,4)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPAPOTP32(cu,axyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,nvpz,nzv,   &
     &kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpetfield3(dcu,exyz,ffe,affp,ci,wf,tfield,nx,ny,nz,    &
     &kstrt,nvpy,nvpz)
! solves 3d poisson's equation for unsmoothed transverse electric field
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(in) :: affp, ci
      real, intent(inout) :: wf, tfield
      complex, dimension(:,:,:,:), intent(in) :: dcu
      complex, dimension(:,:,:,:), intent(inout) :: exyz
      complex, dimension(:,:,:), intent(in) :: ffe
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(dcu,2); kxyp = size(dcu,3); kyzp = size(dcu,4)
      nzhd = size(ffe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call VMPPETFIELD332(dcu,exyz,ffe,affp,ci,wf,nx,ny,nz,kstrt,nvpy,  &
     &nvpz,nzv,kxyp,kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpsmooth3(q,qs,ffc,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! provides a 3d scalar smoothing function
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(in) :: q
      complex, dimension(:,:,:), intent(inout) :: qs
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(q,1); kxyp = size(q,2); kyzp = size(q,3)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPSMOOTH32(q,qs,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp, &
     &nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpsmooth33(cu,cus,ffc,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! provides a 3d vector smoothing function
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(in) :: cu
      complex, dimension(:,:,:,:), intent(inout) :: cus
      complex, dimension(:,:,:), intent(in) :: ffc
! local data
      integer :: nzv, kxyp, kyzp, nzhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(cu,2); kxyp = size(cu,3); kyzp = size(cu,4)
      nzhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MPPSMOOTH332(cu,cus,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,   &
     &kyzp,nzhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdmodes3(pot,pott,tfield,nx,ny,nz,modesx,modesy,     &
     &modesz,kstrt,nvpy,nvpz)
! extracts and stores lowest order scalar modes to unpacked array
      implicit none
      integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
      integer, intent(in) :: kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(in) :: pot
      complex, dimension(:,:,:), intent(inout) :: pott
! local data
      integer :: nzv, kxyp, kyzp, modeszd, modesxpd, modesypd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(pot,1); kxyp = size(pot,2); kyzp = size(pot,3)
      modeszd = size(pott,1); modesxpd = size(pott,2)
      modesypd = size(pott,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRDMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz,kstrt,nvpy&
     &,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrmodes3(pot,pott,tfield,nx,ny,nz,modesx,modesy,     &
     &modesz,kstrt,nvpy,nvpz)
! reads and copies lowest order scalar modes to packed array, writing
! zeroes to high order modes
      implicit none
      integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
      integer, intent(in) :: kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:), intent(inout) :: pot
      complex, dimension(:,:,:), intent(in) :: pott
! local data
      integer :: nzv, kxyp, kyzp, modeszd, modesxpd, modesypd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nzv = size(pot,1); kxyp = size(pot,2); kyzp = size(pot,3)
      modeszd = size(pott,1); modesxpd = size(pott,2)
      modesypd = size(pott,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz,kstrt,nvpy&
     &,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mprdvmodes3(vpot,vpott,tfield,nx,ny,nz,modesx,modesy,  &
     &modesz,kstrt,nvpy,nvpz)
! extracts and stores lowest order vector modes to unpacked array
      implicit none
      integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
      integer, intent(in) :: kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(in) :: vpot
      complex, dimension(:,:,:,:), intent(inout) :: vpott
! local data
      integer :: ndim, nzv, kxyp, kyzp, modeszd, modesxpd, modesypd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(vpot,1); nzv = size(vpot,2)
      kxyp = size(vpot,3); kyzp = size(vpot,4)
      modeszd = size(vpott,2); modesxpd = size(vpott,3)
      modesypd = size(vpott,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPRDVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,modesz,ndim,  &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpwrvmodes3(vpot,vpott,tfield,nx,ny,nz,modesx,modesy,  &
     &modesz,kstrt,nvpy,nvpz)
! extracts and stores lowest order vector modes to unpacked array
      implicit none
      integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
      integer, intent(in) :: kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(inout) :: vpot
      complex, dimension(:,:,:,:), intent(in) :: vpott
! local data
      integer :: ndim, nzv, kxyp, kyzp, modeszd, modesxpd, modesypd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(vpot,1); nzv = size(vpot,2)
      kxyp = size(vpot,3); kyzp = size(vpot,4)
      modeszd = size(vpott,2); modesxpd = size(vpott,3)
      modesypd = size(vpott,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPWRVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,modesz,ndim,  &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpset_pcvzero3(exyz,tfield,nx,ny,nz,kstrt,nvpy,nvpz)
! zeros out transverse field array.
      implicit none
      integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
      real, intent(inout) :: tfield
      complex, dimension(:,:,:,:), intent(inout) :: exyz
! local data
      integer :: ndim, nzv, kxyp, kyzp
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(exyz,1); nzv = size(exyz,2)
      kxyp = size(exyz,3); kyzp = size(exyz,4)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call SET_PCVZERO3(exyz,nx,ny,nz,kstrt,nvpy,nvpz,ndim,nzv,kxyp,kyzp&
     &)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
      end module
