!-----------------------------------------------------------------------
! Interface file for libmfft1.f
      module libmfft1_h
      implicit none
!
      interface
         subroutine WFFT1RINIT(mixup,sct,indx,nxhd)
         implicit none
         integer, intent(in) :: indx, nxhd
         integer, dimension(nxhd), intent(inout) :: mixup
         complex, dimension(nxhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
         real, dimension(nxd), intent(inout) :: f
         complex, dimension(nxhd), intent(inout) :: t
         integer, dimension(nxhd), intent(in) :: mixup
         complex, dimension(nxhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1R2X(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
         real, dimension(2,nxd), intent(inout) :: f
         complex, dimension(2,nxhd), intent(inout) :: t
         integer, dimension(nxhd), intent(in) :: mixup
         complex, dimension(nxhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1R3X(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer , intent(in):: isign, indx, nxd, nxhd
         real, dimension(3,nxd), intent(inout) :: f
         complex, dimension(3,nxhd), intent(inout) :: t
         integer, dimension(nxhd), intent(in) :: mixup
         complex, dimension(nxhd), intent(in) :: sct
         end subroutine
      end interface
!
      end module
