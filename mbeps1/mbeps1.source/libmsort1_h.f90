!-----------------------------------------------------------------------
! Interface file for libmsort1.f
      module libmsort1_h
      implicit none
!
      interface
         subroutine PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx&
     &,mx,mx1,npbmx,ntmax,irc2)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1, npbmx, ntmax
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1), intent(inout) :: ppbuff
         integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         integer, dimension(2), intent(inout) :: irc2
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &nx,mx,mx1,npbmx,ntmax,irc2)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1, npbmx, ntmax
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1), intent(inout) :: ppbuff
         integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(in) :: ihole
         integer, dimension(2), intent(inout) :: irc2
         end subroutine
      end interface
!
      interface
         subroutine PPRSNCL1L(ncl,mx1)
         implicit none
         integer, intent(in) :: mx1
         integer, dimension(2,mx1), intent(inout) :: ncl
         end subroutine
      end interface
!
      interface
         subroutine PPRSTOR1L(ppart,ppbuff,ncl,ihole,idimp,nppmx,mx1,   &
     &npbmx,ntmax)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, npbmx, ntmax
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1), intent(in) :: ppbuff
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(in) :: ihole
         end subroutine
      end interface
!
      end module
