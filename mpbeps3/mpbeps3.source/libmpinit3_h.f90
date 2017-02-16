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
         subroutine PVDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx, &
     &npy,npz,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, npz, idimp, npmax
         integer, intent(inout) :: ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PVRDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci, &
     &npx,npy,npz,idimp,npmax,ierr)
         implicit none
         integer, intent(in) :: nps, npp, npx, npy, npz, idimp, npmax
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