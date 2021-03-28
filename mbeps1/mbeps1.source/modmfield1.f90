!-----------------------------------------------------------------------
!
      module mfield1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmfield1.f
! mpois1_init calculates table needed by 2d poisson solver
!             calls POIS1
! mpois1 solves 1d poisson's equation for smoothed electric field
!        calls POIS1
! maddqei1 adds electron and ion densities
!          calls ADDQEI1
! mibpois1 solves 1-2/2d poisson's equation for unsmoothed magnetic
!          field
!          calls IBPOIS13
! mmaxwel1 solves 1-2/2d maxwell's equation for unsmoothed transverse
!          electric and magnetic fields using verlet algorithm
!          calls MAXWEL1
! mamaxwel1 solves 1-2/2d maxwell's equation for unsmoothed transverse
!           electric and magnetic fields using analytic algorithm due to
!           irving haber
!           calls AMAXWEL1
! memfield1 adds and smooths complex vector fields in fourier space
!           calls EMFIELD1
! mbmfield1 copies and smooths complex vector fields in fourier space
!           calls BMFIELD1
! maddcuei1 adds electron and ion current densities
!           calls ADDCUEI13
! meaddext1 add external traveling wave field to electric field for
!           1d code
!           calls EADDEXT1
! meaddext13 add external traveling and external circularly polarized
!            wave fields to electric field for 1-2/2d code
!            calls EADDEXT13
! mbaddext1 adds constant to magnetic field for 1-2/2d code
!           calls BADDEXT1
! maddvrfield1 calculates a = b + c
!              calls ADDVRFIELD13
! mbbpois1 solves 1-2/2d poisson's equation in fourier space for
!           smoothed magnetic field
!           calls BBPOIS13
! mdcuperp1 calculates transverse part of the derivative of the current
!           density from the momentum flux
!           calls DCUPERP13
! madcuperp1 calculates transverse part of the derivative of the
!            current density from the momentum flux and acceleration
!            density
!            calls ADCUPERP13
! mepois1_init calculates table needed by 1-2/2d poisson solver for
!              transverse electric field
!              calls EPOIS13
! mepois1 solves 1-2/2d poisson's equation for smoothed transverse
!         electric field
!         calls EPOIS13
! mpot1 solves 1d poisson's equation for potential
!       calls POTP1
! melfield1 solves 1d or 1d poisson's equation for unsmoothed
!           longitudinal electric field
!           calls ELFIELD1
! mdivf1 calculates the divergence in fourier space
!        calls DIVF1
! mgradf1 calculates the gradient in fourier space
!         calls GRADF1
! mcurlf1 calculates the curl in fourier space
!         calls CURLF1
! mavpot1 calculates 1-2/2d vector potential from magnetic field
!         calls AVPOT13
! mcuave1 averages current in fourier space for 1-2/2d code
!         calls CUAVE13
! mavrpot1 solves 1-2/2d poisson's equation for the radiative part of
!          the vector potential
!          calls AVRPOT13
! mapot1 solves 1-2/2d poisson's equation for vector potential
!        calls APOTP13
! metfield1 solves 1-2/2d poisson's equation for unsmoothed transverse
!           electric field
!           calls ETFIELD13
! msmooth1 provides a 1d scalar smoothing function
!          calls SMOOTH1
! msmooth13 provides a 1d vector smoothing function
!           calls SMOOTH13
! mrdmodes1 extracts and stores lowest order scalar modes to unpacked
!           array
!           calls RDMODES1
! mwrmodes1 reads and copies lowest order scalar modes to packed array,
!           writing zeroes to high order modes
!           calls WRMODES1
! mrdvmodes1 extracts and stores lowest order vector modes to unpacked
!            array
!            calls RDVMODES1
! mwrvmodes1 reads and copies lowest order vector modes to packed array
!            writing zeroes to high order modes
!            calls WRVMODES1
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: january 13, 2018
!
      use libmfield1_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mpois1_init(ffc,ax,affp,nx)
! calculates table needed by 1d poisson solver
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ax, affp
      complex, dimension(:), intent(inout) :: ffc
! local data
      integer :: isign = 0
      real :: we
      real, dimension(1)  :: q, fx
      call POIS1(q,fx,isign,ffc,ax,affp,we,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpois1(q,fx,ffc,we,tfield,nx)
! solves 1d poisson's equation for smoothed electric field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: we, tfield
      real, dimension(:), intent(in)  :: q
      real, dimension(:), intent(inout) :: fx
      complex, dimension(:), intent(inout) :: ffc
! local data
      integer :: isign = -1
      real :: ax, affp
      integer, dimension(4) :: itime
      double precision :: dtime
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call POIS1(q,fx,isign,ffc,ax,affp,we,nx)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine maddqei1(qe,qi,tfield,nx)
! adds electron and ion densities
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:), intent(inout) :: qe
      real, dimension(:), intent(in) :: qi
! local data
      integer :: nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(qe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ADDQEI1(qe,qi,nx,nxe)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mibpois1(cu,byz,ffc,ci,wm,tfield,nx)
! solves 2-1/2d poisson's equation for unsmoothed magnetic field
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      real, dimension(:,:), intent(in) :: cu
      complex, dimension(:,:), intent(inout) :: byz
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(cu,2)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call IBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mmaxwel1(eyz,byz,cu,ffc,ci,dt,wf,wm,tfield,nx)
! solves 1-2/2d maxwell's equation for unsmoothed transverse electric
! and magnetic fields using verlet algorithm
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ci, dt
      real, intent(inout) :: wf, wm, tfield
      complex, dimension(:,:), intent(inout) :: eyz, byz
      real, dimension(:,:), intent(in)  :: cu
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(cu,2)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call MAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mamaxwel1(eyz,byz,cu,ffc,ci,dt,wf,wm,tfield,nx)
! solves 1-2/2d maxwell's equation for unsmoothed transverse electric
! and magnetic fields using analytic algorithm due to irving haber
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ci, dt
      real, intent(inout) :: wf, wm, tfield
      complex, dimension(:,:), intent(inout) :: eyz, byz
      real, dimension(:,:), intent(in)  :: cu
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(cu,2)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call AMAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine memfield1(fxyz,fx,eyz,ffc,tfield,nx)
! adds and smooths complex vector fields in fourier space
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(inout) :: fxyz
      real, dimension(:), intent(in) :: fx
      complex, dimension(:,:), intent(in) :: eyz
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(fxyz,2)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call EMFIELD1(fxyz,fx,eyz,ffc,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mbmfield1(fyz,eyz,ffc,tfield,nx)
! copies and smooths complex vector fields in fourier space
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(inout) :: fyz
      complex, dimension(:,:), intent(in) :: eyz
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(fyz,2)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call BMFIELD1(fyz,eyz,ffc,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine maddcuei1(cue,cui,tfield,nx)
! adds electron and ion current densities
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(inout) :: cue
      real, dimension(:,:), intent(in) :: cui
! local data
      integer :: nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(cue,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ADDCUEI13(cue,cui,nx,nxe)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine meaddext1(fxe,tfield,amodex,freq,time,trmp,toff,el0,er0&
     &,nx)
! add external traveling wave field to electric field for 1-2/2d code
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, intent(in) :: amodex, freq, time, trmp, toff, el0, er0
      real, dimension(:), intent(inout) :: fxe
! local data
      integer :: nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(fxe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call EADDEXT1(fxe,amodex,freq,time,trmp,toff,el0,er0,nx,nxe)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine meaddext13(fxyze,tfield,amodex,freq,time,trmp,toff,el0,&
     &er0,ey0,ez0,nx)
! add external traveling and external circularly polarized wave fields
! to electric field for 1-2/2d code
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, intent(in) :: amodex, freq, time, trmp, toff, el0, er0
      real, intent(in) :: ey0, ez0
      real, dimension(:,:), intent(inout) :: fxyze
! local data
      integer :: ndim, nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(fxyze,1); nxe = size(fxyze,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (3)
         call EADDEXT13(fxyze,amodex,freq,time,trmp,toff,el0,er0,ey0,ez0&
     &,nx,nxe)
      case default
         write (*,*) 'meaddext13 error: ndim=', ndim
      end select
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mbaddext1(byz,tfield,omy,omz,nx)
! adds constant to magnetic field for 1-2/2d code
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, intent(in) :: omy, omz
      real, dimension(:,:), intent(inout) :: byz
! local data
      integer :: nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(byz,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call BADDEXT1(byz,omy,omz,nx,nxe)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine maddvrfield1(a,b,c,tfield)
! calculates a = b + c
      implicit none
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(inout) :: a
      real, dimension(:,:), intent(in) :: b
      real, dimension(:), intent(in) :: c
! local data
      integer :: nxe
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxe = size(a,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ADDVRFIELD13(a,b,c,nxe)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mbbpois1(cu,byz,ffc,ci,wm,tfield,nx)
! solves 1-2/2d poisson's equation in fourier space for smoothed
! magnetic field
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      real, dimension(:,:), intent(in) :: cu
      real, dimension(:,:), intent(inout) :: byz
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(cu,2)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call BBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mdcuperp1(dcu,amu,tfield,nx)
! calculates transverse part of the derivative of the current density
! from the momentum flux
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(inout) :: dcu
      real, dimension(:,:), intent(in) :: amu
! local data
      integer :: nxvh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(dcu,2)/2
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call DCUPERP13(dcu,amu,nx,nxvh)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine madcuperp1(dcu,amu,tfield,nx)
! calculates transverse part of the derivative of the current density
! from the momentum flux and acceleration density
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(inout) :: dcu
      real, dimension(:,:), intent(in) :: amu
! local data
      integer :: nxvh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(dcu,2)/2
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ADCUPERP13(dcu,amu,nx,nxvh)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mepois1_init(ffe,ax,affp,wp0,ci,nx)
! calculates table needed by 1-2/2d poisson solver for transverse
! electric field
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ax, affp, wp0, ci
      complex, dimension(:), intent(inout) :: ffe
! local data
      integer :: isign = 0
      integer :: nxvh, nxhd
      real :: wf
      real, dimension(2,1,1)  :: dcu
      real, dimension(2,1,1) :: eyz
      nxvh = size(dcu,2)/2
      nxhd = size(ffe,1)
      call EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,nxhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mepois1(dcu,eyz,ffe,affp,ci,wf,tfield,nx)
! solves 1-2/2d poisson's equation for smoothed transverse electric
! field
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: affp, ci
      real, intent(inout) :: wf, tfield
      real, dimension(:,:), intent(in) :: dcu
      real, dimension(:,:), intent(inout) :: eyz
      complex, dimension(:), intent(inout) :: ffe
! local data
      integer :: isign = -1
      integer :: nxvh, nxhd
      real :: ax, wp0
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(dcu,2)/2
      nxhd = size(ffe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpot1(q,pot,ffc,we,tfield,nx)
! solves 1d poisson's equation for potential
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: we, tfield
      real, dimension(:), intent(in) :: q
      complex, dimension(:), intent(inout) :: pot
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(q,1)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call POTP1(q,pot,ffc,we,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine melfield1(q,fx,ffc,we,tfield,nx)
! solves 1d poisson's equation for unsmoothed longitudinal electric
! field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: we, tfield
      real, dimension(:), intent(in)  :: q
      complex, dimension(:), intent(inout) :: fx
      complex, dimension(:), intent(inout) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(q,1)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ELFIELD1(q,fx,ffc,we,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mdivf1(f,df,tfield,nx)
! calculates the divergence in fourier space
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(in) :: f
      complex, dimension(:), intent(inout) :: df
! local data
      integer :: ndim, nxvh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxvh = size(f,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call DIVF1(f,df,nx,ndim,nxvh)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mgradf1(df,f,tfield,nx)
! calculates the gradient in fourier space
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      complex, dimension(:), intent(in) :: df
      complex, dimension(:,:), intent(inout) :: f
! local data
      integer :: ndim, nxvh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(df,1)
      ndim = size(f,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call GRADF1(df,f,nx,ndim,nxvh)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcurlf1(f,g,tfield,nx)
! calculates the curl in fourier space
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(in) :: f
      complex, dimension(:,:), intent(inout) :: g
! local data
      integer :: nxvh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(f,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call CURLF1(f,g,nx,nxvh)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mavpot1(byz,ayz,tfield,nx)
! calculates 1-2/2d vector potential from magnetic field
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(in) :: byz
      complex, dimension(:,:), intent(inout) :: ayz
! local data
      integer :: nxvh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(byz,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call AVPOT13(byz,ayz,nx,nxvh)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcuave1(cuave,cunew,cuold,tfield,nx)
! averages current in fourier space for 1-2/2d code
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(in) :: cunew, cuold
      complex, dimension(:,:), intent(inout) :: cuave
! local data
      integer :: nxvh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(cuave,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call CUAVE13(cuave,cunew,cuold,nx,nxvh)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mavrpot1(ayz,byz,ffc,ci,tfield,nx)
! solves 1-2/2d poisson's equation for the radiative part of the vector
! potential
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ci
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(inout) :: ayz
      complex, dimension(:,:), intent(in) :: byz
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(byz,2)
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call AVRPOT13(ayz,byz,ffc,ci,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mapot1(cu,ayz,ffc,ci,wm,tfield,nx)
! solves 1-2/2d poisson's equation in fourier space for vector potential
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ci
      real, intent(inout) :: wm, tfield
      real, dimension(:,:), intent(in) :: cu
      complex, dimension(:,:), intent(inout) :: ayz
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(cu,2)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call APOTP13(cu,ayz,ffc,ci,wm,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine metfield1(dcu,eyz,ffe,ci,wf,tfield,nx)
! solves 1-2/2d poisson's equation for unsmoothed transverse electric
! field
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: ci
      real, intent(inout) :: wf, tfield
      real, dimension(:,:), intent(in) :: dcu
      complex, dimension(:,:), intent(inout) :: eyz
      complex, dimension(:), intent(in) :: ffe
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(dcu,2)/2
      nxhd = size(ffe,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ETFIELD13(dcu,eyz,ffe,ci,wf,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine msmooth1(q,qs,ffc,tfield,nx)
! provides a 1d scalar smoothing function
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:), intent(in) :: q
      complex, dimension(:), intent(inout) :: qs
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(q,1)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call SMOOTH1(q,qs,ffc,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
      subroutine msmooth13(cu,cus,ffc,tfield,nx)
! provides a 1d vector smoothing function
      implicit none
      integer, intent(in) :: nx
      real, intent(inout) :: tfield
      real, dimension(:,:), intent(in) :: cu
      complex, dimension(:,:), intent(inout) :: cus
      complex, dimension(:), intent(in) :: ffc
! local data
      integer :: nxvh, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(cu,2)/2
      nxhd = size(ffc,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call SMOOTH13(cu,cus,ffc,nx,nxvh,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrdmodes1(pot,pott,tfield,nx,modesx)
! extracts and stores lowest order scalar modes to unpacked array
      implicit none
      integer, intent(in) :: nx, modesx
      real, intent(inout) :: tfield
      complex, dimension(:), intent(in) :: pot
      complex, dimension(:), intent(inout) :: pott
! local data
      integer :: nxvh, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(pot,1)
      modesxd = size(pott,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call RDMODES1(pot,pott,nx,modesx,nxvh,modesxd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwrmodes1(pot,pott,tfield,nx,modesx)
! reads and copies lowest order scalar modes to packed array, writing
! zeroes to high order modes
      implicit none
      integer, intent(in) :: nx, modesx
      real, intent(inout) :: tfield
      complex, dimension(:), intent(inout) :: pot
      complex, dimension(:), intent(in) :: pott
! local data
      integer :: nxvh, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxvh = size(pot,1)
      modesxd = size(pott,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call WRMODES1(pot,pott,nx,modesx,nxvh,modesxd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mrdvmodes1(vpot,vpott,tfield,nx,modesx)
! extracts and stores lowest order vector modes to unpacked array
      implicit none
      integer, intent(in) :: nx, modesx
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(in) :: vpot
      complex, dimension(:,:), intent(inout) :: vpott
! local data
      integer :: ndim, nxvh, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(vpot,1); nxvh = size(vpot,2)
      modesxd = size(vpott,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call RDVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mwrvmodes1(vpot,vpott,tfield,nx,modesx)
! extracts and stores lowest order vector modes to unpacked array
      implicit none
      integer, intent(in) :: nx, modesx
      real, intent(inout) :: tfield
      complex, dimension(:,:), intent(inout) :: vpot
      complex, dimension(:,:), intent(in) :: vpott
! local data
      integer :: ndim, nxvh, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(vpot,1); nxvh = size(vpot,2)
      modesxd = size(vpott,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call WRVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
! record time
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
      end subroutine
!
      end module
