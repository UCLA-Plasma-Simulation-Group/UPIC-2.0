!-----------------------------------------------------------------------
! Interface file for libgks1.f
      module libgraf1_h
      implicit none
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GROPEN
         implicit none
         end subroutine
      end interface
      interface
         subroutine GRCLOSE
         implicit none
         end subroutine
      end interface
      interface
         subroutine SETNPLT(nplt,irc)
         implicit none
         integer, intent(in) :: nplt
         integer, intent(inout) :: irc
         end subroutine
      end interface
      interface
         subroutine RSTSCRN
         implicit none
         end subroutine
      end interface
      interface
         subroutine DISPR(f,label,xmin,xmax,isc,ist,mks,nx,nxv,ngs,chr, &
     &chrs,irc)
         implicit none
         integer, intent(in) :: isc, ist, mks, nx, nxv, ngs
         integer, intent(inout) :: irc
         real, intent(in) :: xmin, xmax
         character(len=*), intent(in) :: label, chr
         character(len=*), dimension(ngs), intent(in) :: chrs
!        real, dimension(*), intent(in) :: f
         real, intent(inout) :: f
         end subroutine
      end interface
      interface
         subroutine DISPS(f,label,xmin,xmax,isc,ist,nx,chr,irc)
         implicit none
         integer, intent(in) :: isc, ist, nx
         integer, intent(inout) :: irc
         real, intent(in) :: xmin, xmax
         character(len=*), intent(in) :: label, chr
!        real, dimension(*), intent(in) :: f
         real, intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine GRASP13(part,label,itime,isc,nx,iyp,ixp,idimp,npx,np&
     &,irc)
         implicit none
         integer, intent(in) :: itime, isc, nx, iyp, ixp, idimp, npx, np
         integer, intent(inout) :: irc
         character(len=*), intent(in) :: label
         real, dimension(idimp,np), intent(in) :: part
         end subroutine
      end interface
!
      interface
         subroutine PGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp&
     &,nppmx,mx1,irc)
         implicit none
         integer, intent(in) :: itime, isc, nx, iyp, ixp, idimp
         integer, intent(in) :: nppmx, mx1
         integer, intent(inout) :: irc
         character(len=*), intent(in) :: label
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PMGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,    &
     &idimp,nppmx,mx1,ltag,irc)
         implicit none
         integer, intent(in) :: itime, isc, nx, iyp, ixp, idimp
         integer, intent(in) :: nppmx, mx1, ltag
         integer, intent(inout) :: irc
         character(len=*), intent(in) :: label
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      end module
