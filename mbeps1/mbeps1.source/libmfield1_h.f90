!-----------------------------------------------------------------------
! Interface file for libmfield1.f
      module libmfield1_h
      implicit none
!
      interface
         subroutine POIS1(q,fx,isign,ffc,ax,affp,we,nx)
         implicit none
         integer, intent(in) :: isign, nx
         real, intent(in) :: ax, affp
         real, intent(inout) :: we
         real, dimension(nx), intent(in) :: q
         real, dimension(nx), intent(inout) :: fx
         complex, dimension(nx/2), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ADDQEI1(qe,qi,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(nxe), intent(inout) :: qe
         real, dimension(nxe), intent(in) :: qi
         end subroutine
      end interface
!
      interface
         subroutine IBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(2,nxvh), intent(inout) :: byz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(2,nxvh), intent(inout) :: eyz, byz
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine AMAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(2,nxvh), intent(inout) :: eyz, byz
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EMFIELD1(fxyz,fx,eyz,ffc,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, dimension(3,2*nxvh), intent(inout) :: fxyz
         real, dimension(2*nxvh), intent(in) :: fx
         complex, dimension(3,nxvh), intent(in) :: eyz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BMFIELD1(fyz,eyz,ffc,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, dimension(2,2*nxvh), intent(inout) :: fyz
         complex, dimension(2,nxvh), intent(in) :: eyz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ADDCUEI13(cue,cui,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(2,nxe), intent(inout) :: cue
         real, dimension(2,nxe), intent(in) :: cui
         end subroutine
      end interface
!
      interface
         subroutine EADDEXT1(fxe,amodex,freq,time,trmp,toff,el0,er0,nx, &
     &nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: amodex, freq, time, trmp, toff, el0, er0
         real, dimension(nxe), intent(inout) :: fxe
         end subroutine
      end interface
!
      interface
         subroutine EADDEXT13(fxyze,amodex,freq,time,trmp,toff,el0,er0, &
     &ey0,ez0,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: amodex, freq, time, trmp, toff, el0, er0
         real, intent(in) :: ey0, ez0
         real, dimension(3,nxe), intent(inout) :: fxyze
         end subroutine
      end interface
!
      interface
         subroutine BADDEXT1(byz,omy,omz,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: omy, omz
         real, dimension(2,nxe), intent(inout) :: byz
         end subroutine
      end interface
!
      interface
         subroutine ADDVRFIELD13(fxyze,eyze,fxe,nxe)
         implicit none
         integer, intent(in) :: nxe
         real, dimension(3,nxe), intent(inout) :: fxyze
         real, dimension(2,nxe), intent(in) :: eyze
         real, dimension(nxe), intent(in) :: fxe
         end subroutine
      end interface
!
      interface
         subroutine BBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(2,2*nxvh), intent(in) :: cu
         real, dimension(2,2*nxvh), intent(inout) :: byz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine DCUPERP13(dcu,amu,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
         real, dimension(2,2*nxvh), intent(inout) :: dcu
         real, dimension(2,2*nxvh), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine ADCUPERP13(dcu,amu,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
         real, dimension(2,2*nxvh), intent(inout) :: dcu
         real, dimension(2,2*nxvh), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,&
     &nxhd)
         implicit none
         integer, intent(in) :: isign, nx, nxvh, nxhd
         real, intent(in) :: ax, affp, wp0, ci
         real, intent(inout) :: wf
         real, dimension(2,2*nxvh), intent(in) :: dcu
         real, dimension(2,2*nxvh), intent(inout) :: eyz
         complex, dimension(nxhd), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine POTP1(q,pot,ffc,we,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(inout) :: we
         real, dimension(2*nxvh), intent(in) :: q
         complex, dimension(nxvh), intent(inout) :: pot
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ELFIELD1(q,fx,ffc,we,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(inout) :: we
         real, dimension(2*nxvh), intent(in) :: q
         complex, dimension(nxvh), intent(inout) :: fx
         complex, dimension(nxhd), intent(inout) :: ffc
         end subroutine
      end interface
! 
      interface
         subroutine DIVF1(f,df,nx,ndim,nxvh)
         implicit none
         integer, intent(in) :: nx, ndim, nxvh
         complex, dimension(ndim,nxvh), intent(in) :: f
         complex, dimension(nxvh), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine GRADF1(df,f,nx,ndim,nxvh)
         implicit none
         integer, intent(in) :: nx, ndim, nxvh
         complex, dimension(nxvh), intent(in) :: df
         complex, dimension(ndim,nxvh), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine CURLF1(f,g,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
         complex, dimension(2,nxvh), intent(in) :: f
         complex, dimension(2,nxvh), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine AVPOT13(byz,ayz,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
         complex, dimension(2,nxvh), intent(in) :: byz
         complex, dimension(2,nxvh), intent(inout) :: ayz
         end subroutine
      end interface
!
      interface
         subroutine CUAVE13(cuave,cunew,cuold,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
         real, dimension(2,2*nxvh), intent(in) :: cunew, cuold
         complex, dimension(2,nxvh), intent(inout) :: cuave
         end subroutine
      end interface
!
      interface
         subroutine AVRPOT13(ayz,byz,ffc,ci,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         complex, dimension(2,nxvh), intent(in) :: byz
         complex, dimension(2,nxvh), intent(inout) :: ayz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine APOTP13(cu,ayz,ffc,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(2,nxvh), intent(inout) :: ayz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ETFIELD13(dcu,eyz,ffe,ci,wf,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wf
         real, dimension(2,2*nxvh), intent(in) :: dcu
         complex, dimension(2,nxvh), intent(inout) :: eyz
         complex, dimension(nxhd), intent(in) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine SMOOTH1(q,qs,ffc,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, dimension(2*nxvh), intent(in) :: q
         complex, dimension(nxvh), intent(inout) :: qs
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine SMOOTH13(cu,cus,ffc,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(2,nxvh), intent(inout) :: cus
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine RDMODES1(pot,pott,nx,modesx,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, nxvh, modesxd
         complex, dimension(nxvh), intent(in) :: pot
         complex, dimension(modesxd), intent(inout) :: pott
         end subroutine
      end interface
!
      interface
         subroutine WRMODES1(pot,pott,nx,modesx,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, nxvh, modesxd
         complex, dimension(nxvh), intent(inout) :: pot
         complex, dimension(modesxd), intent(in) :: pott
         end subroutine
      end interface
!
      interface
         subroutine RDVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, ndim, nxvh, modesxd
         complex, dimension(ndim,nxvh), intent(in) :: vpot
         complex, dimension(ndim,modesxd), intent(inout) :: vpott
         end subroutine
      end interface
!
      interface
         subroutine WRVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, ndim, nxvh, modesxd
         complex, dimension(ndim,nxvh), intent(inout) :: vpot
         complex, dimension(ndim,modesxd), intent(in) :: vpott
         end subroutine
      end interface
!
      end module
