!-----------------------------------------------------------------------
! Interface file for libmpfsct2.f
      module libmpfsct2_h
      implicit none
!
      interface
         subroutine WPFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: indx, indy, nxhyd, nxyd
         integer, dimension(nxhyd), intent(inout) :: mixup
         complex, dimension(nxyd), intent(inout) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(nyv,kxp2d), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(nyv,kxp2d), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFCST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(nyv,kxp2d), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFCCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(nyv,kxp2d), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFCST2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         real, dimension(2,nyv,kxp2d), intent(inout) :: g
         real, dimension(2,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         real, dimension(2,nyv,kxp2d), intent(inout) :: g
         real, dimension(2,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCST2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,&
     &kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCT2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,&
     &kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCT2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,&
     &kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(2,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCST2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,&
     &kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(2,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFCST2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         real, dimension(3,nyv,kxp2d), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         real, dimension(3,nyv,kxp2d), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCSST2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi&
     &,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCCT2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi&
     &,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCST2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi&
     &,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(3,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCSCT2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi&
     &,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(3,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM4(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(4,2*nxvh,kypd), intent(inout) :: f
         real, dimension(4,nyv,kxp2d), intent(inout) :: g
         real, dimension(4,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCCST2RM4X(f,isign,mixup,sctd,indx,indy,kstrt,   &
     &kypi,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(4,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCSCT2RM4Y(g,isign,mixup,sctd,indx,indy,kstrt,   &
     &kxpi,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(4,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM22(f,g,bs,br,isign,ntpose,mixup,sctd,ttp, &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         real, dimension(2,nyv,kxp2d), intent(inout) :: g
         real, dimension(2,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCCST2RM22X(f,isign,mixup,sctd,indx,indy,kstrt,  &
     &kypi,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCSCT2RM22Y(g,isign,mixup,sctd,indx,indy,kstrt,  &
     &kxpi,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(2,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSST2RM23(f,g,bs,br,isign,ntpose,mixup,sctd,ttp, &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         real, dimension(3,nyv,kxp2d), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSSCT2RM32X(f,isign,mixup,sctd,indx,indy,kstrt,   &
     &kypi,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSSCT2RM23Y(g,isign,mixup,sctd,indx,indy,kstrt,   &
     &kxpi,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(3,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      end module