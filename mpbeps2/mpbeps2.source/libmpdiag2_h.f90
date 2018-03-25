!-----------------------------------------------------------------------
! Interface file for libmpdiag2.f
      module libmpdiag2_h
      implicit none
!
      interface
         subroutine PPVDIST2(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,mxyp1&
     &,nmv,nmvf)
         implicit none
         integer, intent(in) :: nvp, idimp, nppmx, mxyp1, nmv, nmvf
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nmvf,2), intent(inout) :: fv
         real, dimension(nmvf,2,mxyp1), intent(inout) :: sfv
         real, dimension(2,3), intent(inout) :: fvm
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPVDIST23(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,    &
     &mxyp1,nmv,nmvf)
         implicit none
         integer, intent(in) :: nvp, idimp, nppmx, mxyp1, nmv, nmvf
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nmvf,3), intent(inout) :: fv
         real, dimension(nmvf,3,mxyp1), intent(inout) :: sfv
         real, dimension(3,3), intent(inout) :: fvm
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PERPDIST2(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyp1,&
     &nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mxyp1, nmv, nmvf
         real, intent(in) :: ci
         real, intent(inout) :: wk
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nmvf,1), intent(inout) :: fv
         real, dimension(nmvf,2,mxyp1), intent(inout) :: sfv
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PERPDIST23(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyp1&
     &,nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mxyp1, nmv, nmvf
         real, intent(in) :: ci
         real, intent(inout) :: wk
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nmvf,1), intent(inout) :: fv
         real, dimension(nmvf,3,mxyp1), intent(inout) :: sfv
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPVBDIST23(ppart,kpic,fv,sfv,fvm,omx,omy,omz,idimp, &
     &nppmx,mxyp1,nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mxyp1, nmv, nmvf
         real, intent(in) :: omx, omy, omz
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nmvf,2), intent(inout) :: fv
         real, dimension(nmvf,3,mxyp1), intent(inout) :: sfv
         real, dimension(3,3), intent(inout) :: fvm
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPVSDIST2(ppart,kpic,fvs,noff,nmv,mvx,mvy,nxb,nyb,  &
     &idimp,nppmx,mxyp1,nmvf)
         implicit none
         integer, intent(in) :: noff, nmv, mvx, mvy, nxb, nyb
         integer, intent(in) :: idimp, nppmx, mxyp1, nmvf
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nmvf,2,nxb,nyb+1), intent(inout) :: fvs
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPVSDIST23(ppart,kpic,fvs,noff,nmv,mvx,mvy,nxb,nyb, &
     &idimp,nppmx,mxyp1,nmvf)
         implicit none
         integer, intent(in) :: noff, nmv, mvx, mvy, nxb, nyb
         integer, intent(in) :: idimp, nppmx, mxyp1, nmvf
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nmvf,3,nxb,nyb+1), intent(inout) :: fvs
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PVDIST2(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
         implicit none
         integer, intent(in) :: npp, nvp, idimp, npmax, nmv, nmvf
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(nmvf,2), intent(inout) :: fv
         real, dimension(2,3), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine PVDIST23(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
         implicit none
         integer, intent(in) :: npp, nvp, idimp, npmax, nmv, nmvf
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(nmvf,3), intent(inout) :: fv
         real, dimension(3,3), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine PVBDIST23(part,fv,fvm,omx,omy,omz,npp,idimp,npmax,  &
     &nmv,nmvf)
         implicit none
         integer, intent(in) :: npp, idimp, npmax, nmv, nmvf
         real, intent(in) :: omx, omy, omz
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(nmvf,2), intent(inout) :: fv
         real, dimension(3,3), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine PPROFX23L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,  &
     &my,nprd,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, npro, mx, my
         integer, intent(in) :: nprd, nxv, nypmx, mx1, mxyp1
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PRPROFX23L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro, &
     &mx,my,nprd,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, npro, mx, my
         integer, intent(in) :: nprd, nxv, nypmx, mx1, mxyp1
         real, intent(in) :: ci
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPROFX22L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,  &
     &my,nprd,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, npro, mx, my
         integer, intent(in) :: nprd, nxv, nypmx, mx1, mxyp1
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PRPROFX22L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro, &
     &mx,my,nprd,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, npro, mx, my
         integer, intent(in) :: nprd, nxv, nypmx, mx1, mxyp1
         real, intent(in) :: ci
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PGPROFX2L(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,idimp, &
     &nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp
         integer, intent(in) :: idimp, nppmx, npro, nx, mx, my
         integer, intent(in) :: nprd, nxv, nypmx, mx1, mxyp1
         real, intent(in) :: qbm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PGRPROFX2L(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,ci,   &
     &idimp,nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp
         integer, intent(in) :: idimp, nppmx, npro, nx, mx, my
         integer, intent(in) :: nprd, nxv, nypmx, mx1, mxyp1
         real, intent(in) :: qbm, dt, ci
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PGBPROFX23L(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt, &
     &idimp,nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp
         integer, intent(in) :: idimp, nppmx, npro, nx, mx, my
         integer, intent(in) :: nprd, nxv, nypmx, mx1, mxyp1
         real, intent(in) :: qbm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PGRBPROFX23L(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,&
     &ci,idimp,nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp
         integer, intent(in) :: idimp, nppmx, npro, nx, mx, my
         integer, intent(in) :: nprd, nxv, nypmx, mx1, mxyp1
         real, intent(in) :: qbm, dt, ci
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine FLUIDQS23(fms,npro,nx,ny,kstrt,kyp,nprd,nxv,nypmx)
         implicit none
         integer, intent(in) :: npro, nx, ny, kstrt, kyp
         integer, intent(in) :: nprd, nxv, nypmx
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         end subroutine
      end interface
!
      interface
         subroutine FLUIDQS22(fms,npro,nx,ny,kstrt,kyp,nprd,nxv,nypmx)
         implicit none
         integer, intent(in) :: npro, nx, ny, kstrt, kyp
         integer, intent(in) :: nprd, nxv, nypmx
         real, dimension(nprd,nxv,nypmx), intent(inout) :: fms
         end subroutine
      end interface
!
      interface
         subroutine PSTPTRAJ2(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,   &
     &vtsx,dvtx,idimp,nppmx,mxyp1,idps,np,nprobt)
         implicit none
         integer, intent(in) :: kstrt, nst, idimp, nppmx, mxyp1, idps
         integer, intent(inout) :: nprobt
         real, intent(in) :: vtx, vtsx, dvtx
         double precision, intent(in) :: np
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(idps), intent(inout) :: tedges
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(nprobt), intent(inout) :: iprobt
         end subroutine
      end interface
!
      interface
         subroutine PSTPTRAJ23(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,  &
     &vtsx,dvtx,idimp,nppmx,mxyp1,idps,np,nprobt)
         implicit none
         integer, intent(in) :: kstrt, nst, idimp, nppmx, mxyp1, idps
         integer, intent(inout) :: nprobt
         real, intent(in) :: vtx, vtsx, dvtx
         double precision, intent(in) :: np
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(idps), intent(inout) :: tedges
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(nprobt), intent(inout) :: iprobt
         end subroutine
      end interface
!
      interface
         subroutine PPTRAJ2(ppart,kpic,partt,numtp,idimp,nppmx,mxyp1,   &
     &nprobt)
         implicit none
         integer, intent(in) :: idimp, nppmx, mxyp1, nprobt
         integer, intent(inout) :: numtp
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         real, dimension(idimp,nprobt), intent(inout) :: partt
         end subroutine
      end interface
!
      interface
         subroutine PPTRAJ23(ppart,kpic,partt,numtp,idimp,nppmx,mxyp1,  &
     &nprobt)
         implicit none
         integer, intent(in) :: idimp, nppmx, mxyp1, nprobt
         integer, intent(inout) :: numtp
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         real, dimension(idimp,nprobt), intent(inout) :: partt
         end subroutine
      end interface
!
      interface
         subroutine PORDTRAJ2(partt,spartt,tedges,numtp,idimp,idps,     &
     &nprobt)
         implicit none
         integer, intent(in) :: numtp, idimp, idps, nprobt
         real, dimension(idimp,nprobt), intent(in) :: partt
         real, dimension(idimp,nprobt), intent(inout) :: spartt
         real, dimension(idps), intent(in) :: tedges
         end subroutine
      end interface
!
      interface
         subroutine PORDTRAJ23(partt,spartt,tedges,numtp,idimp,idps,    &
     &nprobt)
         implicit none
         integer, intent(in) :: numtp, idimp, idps, nprobt
         real, dimension(idimp,nprobt), intent(in) :: partt
         real, dimension(idimp,nprobt), intent(inout) :: spartt
         real, dimension(idps), intent(in) :: tedges
         end subroutine
      end interface
!
      interface
         subroutine PCPYTRAJ2(partt,part,numtp,idimp,nprobt)
         implicit none
         integer, intent(in) :: numtp, idimp, nprobt
         real, dimension(idimp,nprobt), intent(in) :: partt
         real, dimension(idimp,nprobt), intent(inout) :: part
         end subroutine
      end interface
!
      end module
