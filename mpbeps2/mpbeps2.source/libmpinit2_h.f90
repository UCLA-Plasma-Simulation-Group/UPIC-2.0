!-----------------------------------------------------------------------
! Interface file for libmpinit2.f
      module libmpinit2_h
      implicit none
!
      interface
         subroutine PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,  &
     &idps)
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idps
         integer, intent(inout) :: nyp, noff, nypmx, nypmn
         real, dimension(idps), intent(inout) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PDNICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp, &
     &idps)
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idps
         integer, intent(inout) :: nyp, noff, nypmx, nypmn
         real, dimension(idps), intent(inout) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PUDISTR2(part,edges,npp,npx,npy,nx,ny,idimp,npmax,  &
     &idps,ipbc,ierr)
         implicit none
         integer, intent(in) :: npx, npy, nx, ny, idimp, npmax, idps
         integer, intent(in) :: ipbc
         integer, intent(inout) :: npp, ierr
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PLDISTR2(part,npp,anlx,anly,npx,npy,nx,ny,kstrt,nvp,&
     &idimp,npmax,ipbc,ierr)
         implicit none
         integer, intent(in) :: npx, npy, nx, ny, kstrt, nvp
         integer, intent(in) :: idimp, npmax, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: anlx, anly
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PFDISTR2(part,npp,fnx,argx1,argx2,argx3,fny,argy1,  &
     &argy2,argy3,npx,npy,nx,ny,kstrt,nvp,idimp,npmax,ipbc,ierr)
         implicit none
         integer, intent(in) :: npx, npy, nx, ny, kstrt, nvp
         integer, intent(in) :: idimp, npmax, ipbc
         integer, intent(inout) :: npp, ierr
         double precision, intent(in) :: argx1, argx2, argx3
         double precision, intent(in) :: argy1, argy2, argy3
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
         end subroutine
      end interface
!
      interface
         subroutine PGFDISTR2(part,npp,fnx,argx1,argx2,argx3,fny,argy1, &
     &argy2,argy3,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,idimp,    &
     &npmax,ierr)
         implicit none
         integer, intent(in) :: npx, npy, nx, ny, kstrt, nvp
         integer, intent(in) :: idimp, npmax
         integer, intent(inout) :: npp, ierr
         double precision, intent(in) :: argx1, argx2, argx3
         double precision, intent(in) :: argy1, argy2, argy3
         real, intent(in) :: xmin, xmax, ymin, ymax
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
         end subroutine
      end interface
!
      interface
         subroutine PVDISTR2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,kstrt,&
     &nvp,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
         integer, intent(in) :: idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx, &
     &npy,kstrt,nvp,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
         integer, intent(in) :: idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVRDISTR2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,  &
     &kstrt,nvp,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
         integer, intent(in) :: idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vdx, vdy, ci
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVRDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci, &
     &npx,npy,kstrt,nvp,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
         integer, intent(in) :: idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVBDISTR2H(part,nps,npp,vtr,vtz,vdr,vdz,omx,omy,omz,&
     &npx,npy,kstrt,nvp,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, kstrt, nvp
         integer, intent(in) :: idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtr, vtz, vdr, vdz, omx, omy, omz
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my&
     &,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: idimp, npmax, mx, my, mx1, mxyp1, npp
         integer, intent(in) :: noff
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,npmax), intent(in) :: part
         integer, dimension(mxyp1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PFEDGES2(edges,nyp,noff,fny,argy1,argy2,argy3,nypmx,&
     &nypmn,ny,kstrt,nvp,idps,ipbc)
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idps, ipbc
         integer, intent(inout) :: nyp, noff, nypmx, nypmn
         double precision, intent(in) :: argy1, argy2, argy3
         real, dimension(idps), intent(inout) :: edges
         interface
            function fny(argy1,argy2,argy3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argy1, argy2, argy3
            double precision :: fny
            end function
         end interface
         end subroutine
      end interface
!
      interface
         subroutine PGFEDGES2(edges,nyp,noff,fny,argy1,argy2,argy3,ymin,&
     &ymax,nypmx,nypmn,ny,kstrt,nvp,idps)
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idps
         integer, intent(inout) :: nyp, noff, nypmx, nypmn
         real, intent(in) :: ymin, ymax
         double precision, intent(in) :: argy1, argy2, argy3
         real, dimension(idps), intent(inout) :: edges
         interface
            function fny(argy1,argy2,argy3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argy1, argy2, argy3
            double precision :: fny
            end function
         end interface
         end subroutine
      end interface
!
      interface
         subroutine PFHOLES2(part,edges,npp,ihole,ndim,nc,idimp,npmax,  &
     &idps,ntmax)
         implicit none
         integer, intent(in) :: npp, ndim, nc, idimp, npmax, idps, ntmax
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idps), intent(in) :: edges
         integer, dimension(ntmax+1), intent(inout) :: ihole
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