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
         subroutine PVDISTR2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,idimp,&
     &npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx, &
     &npy,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVRDISTR2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,  &
     &idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vdx, vdy, ci
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVRDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci, &
     &npx,npy,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz, ci
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
      end module