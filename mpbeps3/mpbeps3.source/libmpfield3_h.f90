!-----------------------------------------------------------------------
! Interface file for libmpfield3.f
      module libmpfield3_h
      implicit none
!
      interface
         subroutine MPPOIS332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz&
     &,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         real, intent(in) :: ax, ay, az, affp
         real, intent(inout) :: we
         complex, dimension(nzv,kxyp,kyzp), intent(in)  :: q
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: fxyz
         complex, dimension(nzhd,kxyp,kyzp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPADDQEI32(qe,qi,nyzp,nx,nxe,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, dimension(nxe,nypmx,nzpmx), intent(inout) :: qe
         real, dimension(nxe,nypmx,nzpmx), intent(in) :: qi
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine MPPCUPERP32(cu,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,   &
     &kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine MIPPBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,&
     &nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: cu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: bxyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPMAXWEL32(exyz,bxyz,cu,ffc,affp,ci,dt,wf,wm,nx,ny,&
     &nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         real, intent(in) :: affp, ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: exyz
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: bxyz
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: cu
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPEMFIELD32(fxyz,exyz,ffc,isign,nx,ny,nz,kstrt,nvpy&
     &,nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: fxyz
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: exyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPADDCUEI32(cue,cui,nyzp,nx,nxe,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, dimension(3,nxe,nypmx,nzpmx), intent(inout) :: cue
         real, dimension(3,nxe,nypmx,nzpmx), intent(in) :: cui
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine MPADDAMUI32(amu,amui,nyzp,nx,nxe,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, dimension(6,nxe,nypmx,nzpmx), intent(inout) :: amu
         real, dimension(6,nxe,nypmx,nzpmx), intent(in) :: amui
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine MPPBADDEXT32(bxyz,nyzp,omx,omy,omz,nx,nxe,nypmx,    &
     &nzpmx,idds)
         implicit none
         integer, intent(in) :: nx, nxe, nypmx, nzpmx, idds
         real, intent(in) :: omx, omy, omz
         real, dimension(3,nxe,nypmx,nzpmx), intent(inout) :: bxyz
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine MPPADDVRFIELD32(a,b,c,ndim,nxe,nypmx,nzpmx)
         implicit none
         integer, intent(in) :: ndim, nxe, nypmx, nzpmx
         real, dimension(ndim,nxe,nypmx,nzpmx), intent(inout) :: a
         real, dimension(ndim,nxe,nypmx,nzpmx), intent(in) :: b, c
         end subroutine
      end interface
!
      interface
         subroutine MPPBBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,&
     &nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: cu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: bxyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPDCUPERP32(dcu,amu,nx,ny,nz,kstrt,nvpy,nvpz,nzv,  &
     &kxyp,kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: dcu
         complex, dimension(6,nzv,kxyp,kyzp), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MPPADCUPERP32(dcu,amu,nx,ny,nz,kstrt,nvpy,nvpz,nzv, &
     &kxyp,kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: dcu
         complex, dimension(6,nzv,kxyp,kyzp), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MPPEPOISP332(dcu,exyz,isign,ffe,ax,ay,az,affp,wp0,ci&
     &,wf,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         real, intent(in) :: ax, ay, az, affp, wp0, ci
         real, intent(inout) :: wf
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: dcu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: exyz
         complex, dimension(nzhd,kxyp,kyzp), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine MPPOTP32(q,pot,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv, &
     &kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         real, intent(inout) :: we
         complex, dimension(nzv,kxyp,kyzp), intent(in) :: q
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: pot
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPELFIELD32(q,fxyz,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,&
     &nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         real, intent(inout) :: we
         complex, dimension(nzv,kxyp,kyzp), intent(in)  :: q
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: fxyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPDIVF32(f,df,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,   &
     &kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: f
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine MPPGRADF32(df,f,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,  &
     &kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp
         complex, dimension(nzv,kxyp,kyzp), intent(in) :: df
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine MPPCURLF32(f,g,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,   &
     &kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp,kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: f
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine MPPAVPOT332(bxyz,axyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv, &
     &kxyp,kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: bxyz
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: axyz
         end subroutine
      end interface
!
      interface
         subroutine MPPAVRPOT332(axyz,bxyz,ffc,affp,ci,nx,ny,nz,kstrt,  &
     &nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         real, intent(in) :: affp, ci
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: axyz
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: bxyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPAPOTP32(cu,axyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,   &
     &nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: cu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: axyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPETFIELD332(dcu,exyz,ffe,affp,ci,wf,nx,ny,nz,kstrt&
     &,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         real, intent(in) :: affp, ci
         real, intent(inout) :: wf
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: dcu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: exyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH32(q,qs,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,  &
     &kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         complex, dimension(nzv,kxyp,kyzp), intent(in) :: q
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: qs
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH332(cu,cus,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv&
     &,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: cu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: cus
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine PPRDMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz, &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: kstrt, nvpy, nvpz, nzv, kxyp, kyzp
         integer, intent(in) :: modesxpd, modesypd, modeszd
         complex, dimension(nzv,kxyp,kyzp), intent(in) :: pot
         complex, dimension(modeszd,modesxpd,modesypd), intent(inout) ::&
     & pott
         end subroutine
      end interface
!
      interface
         subroutine PPWRMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz, &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: kstrt, nvpy, nvpz, nzv, kxyp, kyzp
         integer, intent(in) :: modesxpd, modesypd, modeszd
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: pot
         complex, dimension(modeszd,modesxpd,modesypd), intent(in) ::   &
     & pott
         end subroutine
      end interface
!
      interface
         subroutine PPRDVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,     &
     &modesz,ndim,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,      &
     &modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: ndim, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
         integer, intent(in) :: modesxpd, modesypd, modeszd
         complex, dimension(ndim,nzv,kxyp,kyzp), intent(in) :: vpot
         complex, dimension(ndim,modeszd,modesxpd,modesypd),            &
     &intent(inout) :: vpott
         end subroutine
      end interface
!
      interface
         subroutine PPWRVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,     &
     &modesz,ndim,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd       &
     &,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: ndim, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
         integer, intent(in) :: modesxpd, modesypd, modeszd
         complex, dimension(ndim,nzv,kxyp,kyzp), intent(inout) :: vpot
         complex, dimension(ndim,modeszd,modesxpd,modesypd), intent(in) &
     &:: vpott
         end subroutine
      end interface
!
      end module
