!-----------------------------------------------------------------------
! Interface file for libmpinit3.f
      module libmpinit3_h
      implicit none
!
      interface
         subroutine PDICOMP32L(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn, &
     &ny,nz,kstrt,nvpy,nvpz,idps,idds)
         implicit none
         integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, idps, idds
         integer, intent(inout) :: nypmx, nzpmx, nypmn, nzpmn
         integer, dimension(idds), intent(inout) :: nyzp, noff
         real, dimension(idps), intent(inout) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PDNICOMP32L(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,&
     &ny,nz,kstrt,nvpy,nvpz,idps,idds)
         implicit none
         integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, idps, idds
         integer, intent(inout) :: nypmx, nzpmx, nypmn, nzpmn
         integer, dimension(idds), intent(inout) :: nyzp, noff
         real, dimension(idps), intent(inout) :: edges
         end subroutine
      end interface
!
      interface
         subroutine FCOMP32(nvp,nx,ny,nz,nvpy,nvpz,ierr)
         implicit none
         integer, intent(in) :: nvp, nx, ny, nz
         integer, intent(inout) :: nvpy, nvpz, ierr
         end subroutine
      end interface
!
      interface
         subroutine PDISTR32(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,&
     &npy,npz,nx,ny,nz,idimp,npmax,idps,ipbc,ierr)
         implicit none
         integer, intent(in) :: npx, npy, npz, nx, ny, nz, idimp
         integer, intent(in) :: npmax, idps, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PRDISTR32(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,&
     &npx,npy,npz,nx,ny,nz,idimp,npmax,idps,ipbc,ierr)
         implicit none
         integer, intent(in) :: npx, npy, npz, nx, ny, nz, idimp
         integer, intent(in) :: npmax, idps, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PUDISTR32(part,edges,npp,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,idps,ipbc,ierr)
         implicit none
         integer, intent(in) :: npx, npy, npz, nx, ny, nz, idimp, npmax
         integer, intent(in) :: idps, ipbc
         integer, intent(inout) :: npp, ierr
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PLDISTR32(part,npp,anlx,anly,anlz,npx,npy,npz,nx,ny,&
     &nz,kstrt,nvpy,nvpz,idimp,npmax,ipbc,ierr)
         implicit none
         integer, intent(in) :: npx, npy, npz, nx, ny, nz, kstrt
         integer, intent(in) :: nvpy, nvpz, idimp, npmax, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: anlx, anly, anlz
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PFDISTR32(part,npp,fnx,argx1,argx2,argx3,fny,argy1, &
     &argy2,argy3,fnz,argz1,argz2,argz3,npx,npy,npz,nx,ny,nz,kstrt,nvpy,&
     &nvpz,idimp,npmax,ipbc,ierr)
         implicit none
         integer, intent(in) :: npx, npy, npz, nx, ny, nz, kstrt
         integer, intent(in) :: nvpy, nvpz, idimp, npmax, ipbc
         integer, intent(inout) :: npp, ierr
         double precision, intent(in) :: argx1, argx2, argx3
         double precision, intent(in) :: argy1, argy2, argy3
         double precision, intent(in) :: argz1, argz2, argz3
         real, dimension(idimp,npmax), intent(inout) :: part
         interface
            function fnx(argx1,argx2,argx3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argx1, argx2, argx3
            double precision :: fnx
            end function
         end interface
         interface
            function fny(argy1,argy2,argy3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argy1, argy2, argy3
            double precision :: fny
            end function
         end interface
         interface
            function fnz(argz1,argz2,argz3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argz1, argz2, argz3
            double precision :: fnz
            end function
         end interface
         end subroutine
      end interface
!
      interface
         subroutine PVDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx, &
     &npy,npz,kstrt,nvpy,nvpz,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, npz, kstrt
         integer, intent(in) :: nvpy, nvpz, idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVRDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci, &
     &npx,npy,npz,kstrt,nvpy,nvpz,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, npz, kstrt
         integer, intent(in) :: nvpy, nvpz, idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PPDBLKP3L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my&
     &,mz,mx1,myp1,mxyzp1,idds,irc)
         implicit none
         integer, intent(in) :: npp, idimp, npmax, mx, my, mz
         integer, intent(in) :: mx1, myp1, mxyzp1, idds
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,npmax), intent(in) :: part
         integer, dimension(mxyzp1), intent(inout) :: kpic
         integer, dimension(idds), intent(in) :: noff
         end subroutine
      end interface
!
      interface
         subroutine PFEDGES32(edges,nyzp,noff,fny,argy1,argy2,argy3,fnz,&
     &argz1,argz2,argz3,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy,nvpz,  &
     &idps,idds,ipbc)
         implicit none
         integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, idps, idds
         integer, intent(in) :: ipbc
         integer, intent(inout) :: nypmx, nypmn, nzpmx, nzpmn
         double precision, intent(in) :: argy1, argy2, argy3
         double precision, intent(in) :: argz1, argz2, argz3
         real, dimension(idps), intent(inout) :: edges
         integer, dimension(idds), intent(inout) :: nyzp, noff
         interface
            function fny(argy1,argy2,argy3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argy1, argy2, argy3
            double precision :: fny
            end function
         end interface
         interface
            function fnz(argz1,argz2,argz3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argz1, argz2, argz3
            double precision :: fnz
            end function
         end interface
         end subroutine
      end interface
!
      interface
         subroutine PFHOLES32(part,edges,npp,ihole,idimp,npmax,idps,    &
     &ntmax)
         implicit none
         integer, intent(in) :: npp, idimp, npmax, idps, ntmax
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idps), intent(in) :: edges
         integer, dimension(ntmax+1,2), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         function ranorm()
         implicit none
         double precision :: ranorm
         end function
      end interface
!
      interface
         function randum()
         implicit none
         double precision :: randum
         end function
      end interface
!
      interface
         function DLDISTR1(x,anlx,anxi,shift,intg)
         implicit none
         integer, intent(in) :: intg
         double precision, intent(in) :: x, anlx, anxi, shift
         double precision :: DLDISTR1
         end function
      end interface
!
      interface
         function DSDISTR1(x,ans,dkx,phase,intg)
         implicit none
         integer, intent(in) :: intg
         double precision, intent(in) :: x, ans, dkx, phase
         double precision :: DSDISTR1
         end function
      end interface
!
      interface
         function DGDISTR1(x,ang,wi,x0,intg)
         implicit none
         integer, intent(in) :: intg
         double precision, intent(in) :: x, ang, wi, x0
         double precision :: DGDISTR1
         end function
      end interface
!
      interface
         function DHDISTR1(x,anh,wi,x0,intg)
         implicit none
         integer, intent(in) :: intg
         double precision, intent(in) :: x, anh, wi, x0
         double precision :: DHDISTR1
         end function
      end interface
!
      interface
         function DEDISTR1(x,ane,wi,x0,intg)
         implicit none
         integer, intent(in) :: intg
         double precision, intent(in) :: x, ane, wi, x0
         double precision :: DEDISTR1
         end function
      end interface
!
      interface
         function DGDISTR0(x,ang,wi,x0,intg)
         implicit none
         integer, intent(in) :: intg
         double precision, intent(in) :: x, ang, wi, x0
         double precision :: DGDISTR0
         end function
      end interface
!
      end module