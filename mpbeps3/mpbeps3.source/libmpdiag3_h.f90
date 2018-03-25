!-----------------------------------------------------------------------
! Interface file for libmpdiag3.f
      module libmpdiag3_h
      implicit none
!
      interface
         subroutine PPVDIST32(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,    &
     &mxyzp1,nmv,nmvf)
         implicit none
         integer, intent(in) :: nvp, idimp, nppmx, mxyzp1, nmv, nmvf
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(nmvf,3), intent(inout) :: fv
         real, dimension(nmvf,3,mxyzp1), intent(inout) :: sfv
         real, dimension(3,3), intent(inout) :: fvm
         integer, dimension(mxyzp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PERPDIST32(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyp1&
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
         subroutine PPVBDIST32(ppart,kpic,fv,sfv,fvm,omx,omy,omz,idimp, &
     &nppmx,mxyzp1,nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, nppmx, mxyzp1, nmv, nmvf
         real, intent(in) :: omx, omy, omz
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(nmvf,2), intent(inout) :: fv
         real, dimension(nmvf,3,mxyzp1), intent(inout) :: sfv
         real, dimension(3,3), intent(inout) :: fvm
         integer, dimension(mxyzp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPVSDIST32(ppart,kpic,fvs,noff,nmv,mvx,mvy,mvz,nxb, &
     &nyb,nzb,idimp,nppmx,mxyzp1,nmvf,nybmx,idds)
         implicit none
         integer, intent(in) :: nmv, mvx, mvy, mvz, nxb, nyb, nzb
         integer, intent(in) :: idimp, nppmx, mxyzp1, nmvf, nybmx, idds
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(nmvf,3,nxb,nybmx+1,nzb+1), intent(inout) :: fvs
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PVDIST32(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
         implicit none
         integer, intent(in) :: npp, nvp, idimp, npmax, nmv, nmvf
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(nmvf,3), intent(inout) :: fv
         real, dimension(3,3), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine PVBDIST32(part,fv,fvm,omx,omy,omz,npp,idimp,npmax,  &
     &nmv,nmvf)
         implicit none
         integer, intent(in) :: npp, idimp, npmax, nmv, nmvf
         real, intent(in) :: omx, omy, omz
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(nmvf,2), intent(inout) :: fv
         real, dimension(3,2), intent(inout) :: fvm
         end subroutine
      end interface
!
      interface
         subroutine PPROFX32L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,my&
     &,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, mx, my, mz
         integer, intent(in) :: nprd, nxv, nypmx, nzpmx
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(nprd,nxv,nypmx,nzpmx), intent(inout) :: fms
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PRPROFX32L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro, &
     &mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, mx, my, mz
         integer, intent(in) :: nprd, nxv, nypmx, nzpmx
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         real, intent(in) :: ci
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(nprd,nxv,nypmx,nzpmx), intent(inout) :: fms
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PGPROFX32L(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,    &
     &idimp,nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, nx, mx, my, mz
         integer, intent(in) :: nprd, nxv, nypmx, nzpmx
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         real, intent(in) :: qbm, dt
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz
         real, dimension(nprd,nxv,nypmx,nzpmx), intent(inout) :: fms
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
        subroutine PGRPROFX32L(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,ci, &
     &idimp,nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, nx, mx, my, mz
         integer, intent(in) :: nprd, nxv, nypmx, nzpmx
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         real, intent(in) :: qbm, dt, ci
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz
         real, dimension(nprd,nxv,nypmx,nzpmx), intent(inout) :: fms
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PGBPROFX32L(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm, &
     &dt,idimp,nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,    &
     &mxyzp1,idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, nx, mx, my, mz
         integer, intent(in) :: nprd, nxv, nypmx, nzpmx
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         real, intent(in) :: qbm, dt
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         real, dimension(nprd,nxv,nypmx,nzpmx), intent(inout) :: fms
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine PGRBPROFX32L(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,&
     &dt,ci,idimp,nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1, &
     &mxyzp1,idds)
         implicit none
         integer, intent(in) :: nppmx, idimp, npro, nx, mx, my, mz
         integer, intent(in) :: nprd, nxv, nypmx, nzpmx
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         real, intent(in) :: qbm, dt, ci
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx,nzpmx), intent(in) :: fxyz, bxyz
         real, dimension(nprd,nxv,nypmx,nzpmx), intent(inout) :: fms
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(idds), intent(in) :: noff, nyzp
         end subroutine
      end interface
!
      interface
         subroutine FLUIDQS3(fms,npro,nx,ny,nz,kstrt,nvpy,kyp,kzp,nprd, &
     &nxv,nypmx,nzpmx)
         implicit none
         integer, intent(in) :: npro, nx, ny, nz, kstrt, nvpy, kyp, kzp
         integer, intent(in) :: nprd, nxv, nypmx, nzpmx
         real, dimension(nprd,nxv,nypmx,nzpmx), intent(inout) :: fms
         end subroutine
      end interface
!
      interface
!-----------------------------------------------------------------------
         subroutine PSTPTRAJ3(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,   &
     &vtsx,dvtx,nvpy,nvpz,idimp,nppmx,mxyzp1,idps,np,nprobt)
         implicit none
         integer, intent(in) :: kstrt, nst, nvpy, nvpz, idimp, nppmx
         integer, intent(in) :: mxyzp1, idps
         integer, intent(inout) :: nprobt
         real, intent(in) :: vtx, vtsx, dvtx
         double precision, intent(in) :: np
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(idps), intent(inout) :: tedges
         integer, dimension(mxyzp1), intent(in) :: kpic
         integer, dimension(nprobt), intent(inout) :: iprobt
         end subroutine
      end interface
!
      interface
         subroutine PPTRAJ3(ppart,kpic,partt,numtp,idimp,nppmx,mxyzp1,  &
     &nprobt)
         implicit none
         integer, intent(in) :: idimp, nppmx, mxyzp1, nprobt
         integer, intent(inout) :: numtp
         real, dimension(idimp,nppmx,mxyzp1), intent(in) :: ppart
         integer, dimension(mxyzp1), intent(in) :: kpic
         real, dimension(idimp,nprobt), intent(inout) :: partt
         end subroutine
      end interface
!
      interface
         subroutine PORDTRAJ3(partt,spartt,tedges,numtp,idimp,idps,     &
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
         subroutine PCPYTRAJ3(partt,part,numtp,idimp,nprobt)
         implicit none
         integer, intent(in) :: numtp, idimp, nprobt
         real, dimension(idimp,nprobt), intent(in) :: partt
         real, dimension(idimp,nprobt), intent(inout) :: part
         end subroutine
      end interface
!
      end module
