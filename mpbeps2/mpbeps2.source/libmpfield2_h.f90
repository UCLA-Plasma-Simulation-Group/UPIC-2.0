!-----------------------------------------------------------------------
! Interface file for libmpfield2.f
      module libmpfield2_h
      implicit none
!
      interface
         subroutine VMPPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in)  :: q
         complex, dimension(2,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine VMPPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPADDQEI2(qe,qi,nyp,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, dimension(nxe,nypmx), intent(inout) :: qe
         real, dimension(nxe,nypmx), intent(in) :: qi
         end subroutine
      end interface
!
      interface
         subroutine MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine VMIPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp, &
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: bxy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine VMPPMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,  &
     &kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(3,nyv,kxp), intent(inout) :: exy, bxy
         complex, dimension(3,nyv,kxp), intent(in)  :: cu
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,  &
     &nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(3,nyv,kxp), intent(in) :: exy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPADDCUEI23(cue,cui,nyp,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, dimension(3,nxe,nypmx), intent(inout) :: cue
         real, dimension(3,nxe,nypmx), intent(in) :: cui
         end subroutine
      end interface
!
      interface
         subroutine MPPADDAMUI23(amu,amui,nyp,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, dimension(6,nxe,nypmx), intent(inout) :: amu
         real, dimension(6,nxe,nypmx), intent(in) :: amui
         end subroutine
      end interface
!
      interface
         subroutine PPBADDEXT2(bxy,nyp,omx,omy,omz,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, intent(in) :: omx, omy, omz
         real, dimension(3,nxe,nypmx), intent(inout) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine PPADDVRFIELD2(a,b,c,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx), intent(inout) :: a
         real, dimension(ndim,nxe,nypmx), intent(in) :: b, c
         end subroutine
      end interface
!
      interface
         subroutine VMPPBBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp, &
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: bxy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: dcu
         complex, dimension(4,nyv,kxp), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MPPADCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: dcu
         complex, dimension(4,nyv,kxp), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine VMPPEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,&
     &nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp, wp0, ci
         real, intent(inout) :: wf
         complex, dimension(3,nyv,kxp), intent(in) :: dcu
         complex, dimension(3,nyv,kxp), intent(inout) :: exy
         complex, dimension(nyhd,kxp), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine MPPOTP2(q,pot,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(nyv,kxp), intent(inout) :: pot
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPELFIELD22(q,fxy,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in)  :: q
         complex, dimension(2,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPELFIELD23(q,fxy,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(ndim,nyv,kxp), intent(in) :: f
         complex, dimension(nyv,kxp), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine MPPGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(nyv,kxp), intent(in) :: df
         complex, dimension(ndim,nyv,kxp), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine MPPCURLF2(f,g,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(in) :: f
         complex, dimension(3,nyv,kxp), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine MPPAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(in) :: bxy
         complex, dimension(3,nyv,kxp), intent(inout) :: axy
         end subroutine
      end interface
!
      interface
         subroutine MCUAVE23(cuave,cunew,cuold,ny,kxp,nyv)
         implicit none
         integer, intent(in) :: ny, kxp, nyv
         complex, dimension(3,nyv,kxp), intent(in) :: cunew, cuold
         complex, dimension(3,nyv,kxp), intent(inout) :: cuave
         end subroutine
      end interface
!
      interface
         subroutine MPPAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,&
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci
         complex, dimension(3,nyv,kxp), intent(inout) :: axy
         complex, dimension(3,nyv,kxp), intent(in) :: bxy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPAPOTP23(cu,axy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd&
     &)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: axy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPETFIELD23(dcu,exy,ffe,affp,ci,wf,nx,ny,kstrt,nyv,&
     &kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci
         real, intent(inout) :: wf
         complex, dimension(3,nyv,kxp), intent(in) :: dcu
         complex, dimension(3,nyv,kxp), intent(inout) :: exy
         complex, dimension(nyhd,kxp), intent(in) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH2(q,qs,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(nyv,kxp), intent(inout) :: qs
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH23(cu,cus,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: cus
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine PPRDMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,  &
     &kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, kstrt, nyv, kxp
         integer, intent(in) :: modesxpd, modesyd
         complex, dimension(nyv,kxp), intent(in) :: pot
         complex, dimension(modesyd,modesxpd), intent(inout) :: pott
         end subroutine
      end interface
!
      interface
         subroutine PPWRMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,  &
     &kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, kstrt, nyv, kxp
         integer, intent(in) :: modesxpd, modesyd
         complex, dimension(nyv,kxp), intent(inout) :: pot
         complex, dimension(modesyd,modesxpd), intent(in) :: pott
         end subroutine
      end interface
!
      interface
         subroutine PPRDVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,    &
     &kstrt,nyv,kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, kstrt
         integer, intent(in) :: nyv, kxp, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp), intent(in) :: vpot
         complex, dimension(ndim,modesyd,modesxpd), intent(inout) ::     &
     &   vpott
         end subroutine
      end interface
!
      interface
         subroutine PPWRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,    &
     &kstrt,nyv,kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, kstrt
         integer, intent(in) :: nyv, kxp, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp), intent(inout) :: vpot
         complex, dimension(ndim,modesyd,modesxpd), intent(in) :: vpott
         end subroutine
      end interface
!
      interface
         subroutine SET_PCVZERO2(exy,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(ndim,nyv,kxp), intent(inout) :: exy
         end subroutine
      end interface
!
      end module
