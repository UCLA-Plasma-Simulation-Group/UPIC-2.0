!-----------------------------------------------------------------------
! unneeded function in input2mod.f90
      subroutine PPBDCAST(f,nxp)
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f
      end subroutine
