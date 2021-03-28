!-----------------------------------------------------------------------
! Interface file for libgks2.f
      module libgraf2_h
      implicit none
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine CARPET(f,label,isc,ist,nx,ny,nxv,chr,ntc,irc)
         implicit none
         integer, intent(in) :: isc, ist, nx, ny, nxv, ntc
         integer, intent(inout) :: irc
         character(len=*), intent(in) :: label, chr
         real, dimension(nxv,ny), intent(in) :: f
         end subroutine
      end interface
!
      interface
         subroutine CARPETL(f,label,xmin,xmax,ymin,ymax,isc,ist,nx,ny,  &
     &nxv,chr,ntc,irc)
         implicit none
         integer, intent(in) :: isc, ist, nx, ny, nxv, ntc
         integer, intent(inout) :: irc
         real, intent(in) :: xmin, xmax, ymin, ymax
         character(len=*), intent(in) :: label, chr
         real, dimension(nxv,ny), intent(in) :: f
         end subroutine
      end interface
!
      interface
         subroutine CONTUR(f,lf,label,isc,ist,nx,ny,nxv,chr,nc,irc)
         implicit none
         integer, intent(in):: isc, ist, nx, ny, nxv, nc
         integer, intent(inout) :: irc
         character(len=*), intent(in) :: label, chr
         real, dimension(nxv,ny), intent(in) :: f
         integer, dimension(nxv,ny), intent(inout) :: lf
         end subroutine
      end interface
!
      end module
