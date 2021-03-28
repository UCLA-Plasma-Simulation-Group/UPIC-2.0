!-----------------------------------------------------------------------
! Interface file for libminit1.f
      module libminit1_h
      implicit none
!
      interface
         subroutine NEXTRAN1(nextrand,ndim,np)
         implicit none
         integer, intent(in) :: nextrand, ndim, np
         end subroutine
      end interface
!
      interface
         subroutine UDISTR1(part,jstart,npx,idimp,nop,nx,ipbc)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop, nx, ipbc
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine LDISTR2(part,anlx,jstart,npx,idimp,nop,nx,ipbc)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop, nx, ipbc
         real, intent(in) :: anlx
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine FDISTR1(part,fnx,argx1,argx2,argx3,jstart,npx,idimp,&
     &nop,nx,ipbc,ierr)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop, nx, ipbc
         integer, intent(inout) :: ierr
         double precision, intent(in) :: argx1, argx2, argx3
         real, dimension(idimp,nop), intent(inout) :: part
         interface
            function fnx(argx1,argx2,argx3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argx1, argx2, argx3
            double precision :: fnx
            end function
         end interface
         end subroutine
      end interface
!
      interface
         subroutine GFDISTR1(part,fnx,argx1,argx2,argx3,xmin,xmax,jstart&
     &,npx,idimp,nop,nx,ierr)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop, nx
         integer, intent(inout) :: ierr
         real, intent(in) :: xmin, xmax
         double precision, intent(in) :: argx1, argx2, argx3
         real, dimension(idimp,nop), intent(inout) :: part
         interface
            function fnx(argx1,argx2,argx3,intg)
            implicit none
            integer, intent(in) :: intg
            double precision, intent(in) :: argx1, argx2, argx3
            double precision :: fnx
            end function
         end interface
         end subroutine
      end interface
!
      interface
         subroutine VDISTR1(part,vtx,vdx,jstart,npx,idimp,nop)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop
         real, intent(in) :: vtx, vdx
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine VDISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,jstart,npx,   &
     &idimp,nop)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine VRDISTR1(part,vtx,vdx,ci,jstart,npx,idimp,nop)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop
         real, intent(in) :: vtx, vdx, ci
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine VRDISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,ci,jstart,npx&
     &,idimp,nop)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine WDISTR1(part,vtx,vdx,jstart,npx,idimp,nop)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop
         real, intent(in) :: vtx, vdx
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine WDISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,jstart,npx,   &
     &idimp,nop)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine VBDISTR1H(part,vtr,vtz,vdr,vdz,omx,omy,omz,jstart,  &
     &npx,idimp,nop)
         implicit none
         integer, intent(in) :: jstart, npx, idimp, nop
         real, intent(in) :: vtr, vtz, vdr, vdz, omx, omy, omz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP1L(part,kpic,nppmx,idimp,np,nop,mx,mx1,irc)
         implicit none
         integer, intent(in) :: idimp, np, nop, mx, mx1
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,nop), intent(in) :: part
         integer, dimension(mx1), intent(inout) :: kpic
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
