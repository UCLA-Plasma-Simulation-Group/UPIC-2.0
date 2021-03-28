!-----------------------------------------------------------------------
! Interface file for libmdiag1.f
      module libmdiag1_h
      implicit none
!
      interface
         subroutine CSPECT1(fc,wm,pkw,t0,dt,nt,iw,modesx,ntd,iwd,modesxd&
     &)
         implicit none
         integer, intent(in) :: nt, iw, modesx, ntd, iwd, modesxd
         real, intent(in) :: t0, dt
         complex, dimension(ntd,modesxd), intent(in) :: fc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(modesxd,iwd,2), intent(inout) :: pkw
         end subroutine
      end interface
!
      interface
         subroutine ICSPECT1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm,&
     &iwd,modesxd)
         implicit none
         integer, intent(in) :: nt, iw, modesx, nx, norm, iwd, modesxd
         real, intent(in) :: time, t0
         complex, dimension(modesxd), intent(in) :: fc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(modesxd,iwd,2), intent(inout) :: pkw
         double precision, dimension(4,modesxd,iwd), intent(inout) ::   &
     &pks
         end subroutine
      end interface
!
      interface
         subroutine IVCSPECT1(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,nx, &
     &norm,iwd,modesxd)
         implicit none
         integer, intent(in) :: nt, iw, modesx, nx, norm, iwd, modesxd
         real, intent(in) :: time, t0
         complex, dimension(2,modesxd), intent(in) :: fvc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(2,modesxd,iwd,2), intent(inout) :: vpkw
         double precision, dimension(2,4,modesxd,iwd), intent(inout) :: &
     &vpks
         end subroutine
      end interface
!
      interface
         subroutine VPDIST1(ppart,kpic,sfv,fvm,idimp,nppmx,mx1,np,nmv,  &
     &nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, np, nmv, nmvf
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         real, dimension(nmvf,mx1+1), intent(inout) :: sfv
         real, dimension(3), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine VPDIST13(ppart,kpic,sfv,fvm,idimp,nppmx,mx1,np,nmv, &
     &nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, np, nmv, nmvf
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         real, dimension(nmvf,3,mx1+1), intent(inout) :: sfv
         real, dimension(3,3), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine VBPDIST13(ppart,kpic,sfv,fvm,omx,omy,omz,idimp,nppmx&
     &,mx1,np,nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, np, nmv, nmvf
         real, intent(in) :: omx, omy, omz
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         real, dimension(nmvf,3,mx1+1), intent(inout) :: sfv
         real, dimension(3,2), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine ERPDIST1(ppart,kpic,sfv,ci,wk,idimp,nppmx,mx1,nmv,  &
     &nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, nmv, nmvf
         real, intent(in) :: ci
         real, intent(inout) :: wk
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         real, dimension(nmvf,mx1+1), intent(inout) :: sfv
         end subroutine
      end interface
!
      interface
         subroutine ERPDIST13(ppart,kpic,sfv,ci,wk,idimp,nppmx,mx1,nmv, &
     &nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, nmv, nmvf
         real, intent(in) :: ci
         real, intent(inout) :: wk
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         real, dimension(nmvf,mx1+1), intent(inout) :: sfv
         end subroutine
      end interface
!
      interface
         subroutine PVSDIST1(ppart,kpic,fvs,nmv,mvx,nxb,idimp,nppmx,mx1,&
     &nmvf)
         implicit none
         integer, intent(in) :: nmv, mvx, nxb
         integer, intent(in) :: idimp, nppmx, mx1, nmvf
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nmvf,nxb), intent(inout) :: fvs
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PVSDIST13(ppart,kpic,fvs,nmv,mvx,nxb,idimp,nppmx,mx1&
       &,nmvf)
         implicit none
         integer, intent(in) :: nmv, mvx, nxb
         integer, intent(in) :: idimp, nppmx, mx1, nmvf
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nmvf,3,nxb), intent(inout) :: fvs
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VDIST1(part,fv,fvm,idimp,np,nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, np, nmv, nmvf
         real, dimension(idimp,np), intent(in) :: part
         real, dimension(nmvf), intent(inout) :: fv
         real, dimension(3), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine VDIST13(part,fv,fvm,idimp,np,nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, np, nmv, nmvf
         real, dimension(idimp,np), intent(in) :: part
         real, dimension(nmvf,3), intent(inout) :: fv
         real, dimension(3,3), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine VBDIST13(part,fv,fvm,omx,omy,omz,idimp,np,nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, np, nmv, nmvf
         real, intent(in) :: omx, omy, omz
         real, dimension(idimp,np), intent(in) :: part
         real, dimension(nmvf,2), intent(inout) :: fv
         real, dimension(2,2), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine PROFX13L(ppart,fms,kpic,nppmx,idimp,npro,mx,nprd,nxv&
     &,mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, mx
         integer, intent(in) :: nprd, nxv, mx1
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nprd,nxv), intent(inout) :: fms
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine RPROFX13L(ppart,fms,kpic,ci,nppmx,idimp,npro,mx,nprd&
     &,nxv,mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, mx
         integer, intent(in) :: nprd, nxv, mx1
         real, intent(in) :: ci
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nprd,nxv), intent(inout) :: fms
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PROFX1L(ppart,fms,kpic,nppmx,idimp,npro,mx,nprd,nxv,&
     &mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, mx
         integer, intent(in) :: nprd, nxv, mx1
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nprd,nxv), intent(inout) :: fms
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine RPROFX1L(ppart,fms,kpic,ci,nppmx,idimp,npro,mx,nprd,&
     &nxv,mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, mx
         integer, intent(in) :: nprd, nxv, mx1
         real, intent(in) :: ci
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nprd,nxv), intent(inout) :: fms
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPROFX1L(ppart,fx,fms,kpic,qbm,dt,idimp,nppmx,npro, &
     &nx,mx,nprd,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, npro, nx, mx
         integer, intent(in) :: nprd, nxv, mx1
         real, intent(in) :: qbm, dt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nxv), intent(in) :: fx
         real, dimension(nprd,nxv), intent(inout) :: fms
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRPROFX1L(ppart,fx,fms,kpic,qbm,dt,ci,idimp,nppmx,  &
     &npro,nx,mx,nprd,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, npro, nx, mx
         integer, intent(in) :: nprd, nxv, mx1
         real, intent(in) :: qbm, dt, ci
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nxv), intent(in) :: fx
         real, dimension(nprd,nxv), intent(inout) :: fms
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPROFX13L(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,idimp,&
     &nppmx,npro,nx,mx,nprd,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, npro, nx, mx
         integer, intent(in) :: nprd, nxv, mx1
         real, intent(in) :: omx, qbm, dt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         real, dimension(nprd,nxv), intent(inout) :: fms
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRBPROFX13L(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,ci,  &
     &idimp,nppmx,npro,nx,mx,nprd,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, npro, nx, mx
         integer, intent(in) :: nprd, nxv, mx1
         real, intent(in) :: omx, qbm, dt, ci
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         real, dimension(nprd,nxv), intent(inout) :: fms
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine FLUIDQS13(fms,npro,nx,nprd,nxv)
         implicit none
         integer, intent(in) :: npro, nx, nprd, nxv
         real, dimension(nprd,nxv), intent(inout) :: fms
         end subroutine
      end interface
!
      interface
         subroutine FLUIDQS1(fms,npro,nx,nprd,nxv)
         implicit none
         integer, intent(in) :: npro, nx, nprd, nxv
         real, dimension(nprd,nxv), intent(inout) :: fms
         end subroutine
      end interface
!
      interface
         subroutine STPTRAJ1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp, &
     &nppmx,mx1,np,nprobt)
         implicit none
         integer, intent(in) :: nst, idimp, nppmx, mx1, np
         integer, intent(inout) :: nprobt
         real, intent(in) :: vtx, vtsx, dvtx
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(nprobt), intent(inout) :: iprobt
         end subroutine
      end interface
!
      interface
         subroutine STPTRAJ13(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,&
     &nppmx,mx1,np,nprobt)
         implicit none
         integer, intent(in) :: nst, idimp, nppmx, mx1, np
         integer, intent(inout) :: nprobt
         real, intent(in) :: vtx, vtsx, dvtx
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(nprobt), intent(inout) :: iprobt
         end subroutine
      end interface
!
      interface
         subroutine FNPTRAJ1(ppart,kpic,idimp,nppmx,mx1,nprobt)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1
         integer, intent(inout) :: nprobt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine FNPTRAJ13(ppart,kpic,idimp,nppmx,mx1,nprobt)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1
         integer, intent(inout) :: nprobt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PTRAJ1(ppart,kpic,partt,idimp,nppmx,mx1,nprobt)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, nprobt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         real, dimension(idimp,nprobt), intent(inout) ::partt
         end subroutine
      end interface
!
      interface
         subroutine PTRAJ13(ppart,kpic,partt,idimp,nppmx,mx1,nprobt)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, nprobt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         real, dimension(idimp,nprobt), intent(inout) ::partt
         end subroutine
      end interface
!
      interface
         subroutine STPBEAM1(part,npx,idimp,nop)
         implicit none
         integer, intent(in) :: npx, idimp, nop
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine STPBEAM13(part,npx,idimp,nop)
         implicit none
         integer, intent(in) :: npx, idimp, nop
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      end module
