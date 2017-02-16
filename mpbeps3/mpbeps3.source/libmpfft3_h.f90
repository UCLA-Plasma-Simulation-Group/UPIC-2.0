!-----------------------------------------------------------------------
! Interface file for libmpfft3.f
      module libmpfft3_h
      implicit none
!
      interface
         subroutine WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: indx, indy, indz, nxhyzd, nxyzhd
         integer, dimension(nxhyzd), intent(inout) :: mixup
         complex, dimension(nxyzhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT32RM(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,  &
     &indx,indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,    &
     &kxypd,kypd,kyzpd,kzpd,kzyp,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, indz, kstrt
         integer, intent(in) :: nvpy, nvpz, nxvh, nyv, nzv, kxyp, kyp
         integer, intent(in) :: kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer, intent(in) :: kzyp, nxhyzd, nxyzhd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd,kzpd), intent(inout) :: f
         complex, dimension(nyv,kxypd,kzpd), intent(inout) :: g
         complex, dimension(nzv,kxypd,kyzpd), intent(inout) :: h
         complex, dimension(kxyp*kzyp,kzp), intent(inout) :: bs, br
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT32RM3(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp, &
     &indx,indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,    &
     &kxypd,kypd,kyzpd,kzpd,kzyp,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, indz, kstrt
         integer, intent(in) :: nvpy, nvpz, nxvh, nyv, nzv, kxyp, kyp
         integer, intent(in) :: kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer, intent(in) :: kzyp, nxhyzd, nxyzhd
         real, intent(inout) :: ttp
         real, dimension(3,2*nxvh,kypd,kzpd), intent(inout) :: f
         complex, dimension(3,nyv,kxypd,kzpd), intent(inout) :: g
         complex, dimension(3,nzv,kxypd,kyzpd), intent(inout) :: h
         complex, dimension(3,kxyp*kzyp,kzp), intent(inout) :: bs, br
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT32RMN(f,g,h,bs,br,ss,isign,ntpose,mixup,sct,  &
     &ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,&
     &kxypd,kypd,kyzpd,kzpd,kzyp,ndim,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, indz, kstrt
         integer, intent(in) :: nvpy, nvpz, nxvh, nyv, nzv, kxyp, kyp
         integer, intent(in) :: kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer, intent(in) :: kzyp, ndim, nxhyzd, nxyzhd
         real, intent(inout) :: ttp
         real, dimension(ndim,2*nxvh,kypd,kzpd), intent(inout) :: f
         complex, dimension(ndim,nyv,kxypd,kzpd), intent(inout) :: g
         complex, dimension(ndim,nzv,kxypd,kyzpd), intent(inout) :: h
         complex, dimension(ndim,kxyp*kzyp,kzp), intent(inout) :: bs, br
         complex, dimension(ndim,nxvh,kzpd) ::  ss
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RMXX(f,isign,mixup,sct,indx,indy,indz,kstrt, &
     &nvp,kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvp
         integer, intent(in) :: kypi, kypp, nxvh, kzpp, kypd, kzpd
         integer, intent(in) :: nxhyzd, nxyzhd
         real, dimension(2*nxvh,kypd,kzpd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RMXY(g,isign,mixup,sct,indx,indy,indz,kstrt, &
     &nvpy,nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy
         integer, intent(in) :: nvpz, kxypi, kxypp, nyv, kzpp, kxypd
         integer, intent(in) :: kzpd, nxhyzd, nxyzhd
         complex, dimension(nyv,kxypd,kzpd), intent(inout) :: g
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RMXZ(h,isign,mixup,sct,indx,indy,indz,kstrt, &
     &nvpy,nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy
         integer, intent(in) :: nvpz, kxypi, kxypp, nzv, kyzp, kxypd
         integer, intent(in) :: kyzpd, nxhyzd, nxyzhd
         complex, dimension(nzv,kxypd,kyzpd), intent(inout) :: h
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RM2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kypi, kypp
         integer, intent(in) :: nxvh, kypd, nxhyd, nxyhd
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RM3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,&
     &nvp,kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvp
         integer, intent(in) :: kypi, kypp, nxvh, kzpp, kypd, kzpd
         integer, intent(in) :: nxhyzd, nxyzhd
         real, dimension(3,2*nxvh,kypd,kzpd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RM3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,&
     &nvpy,nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy
         integer, intent(in) :: nvpz, kxypi, kxypp, nyv, kzpp, kxypd
         integer, intent(in) :: kzpd, nxhyzd, nxyzhd
         complex, dimension(3,nyv,kxypd,kzpd), intent(inout) :: g
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RM3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,&
     &nvpy,nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy
         integer, intent(in) :: nvpz, kxypi, kxypp, nzv, kyzp, kxypd
         integer, intent(in) :: kyzpd, nxhyzd, nxyzhd
         complex, dimension(3,nzv,kxypd,kyzpd), intent(inout) :: h
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RMNXX(f,ss,isign,mixup,sct,indx,indy,indz,   &
     &kstrt,nvp,kypi,kypp,nxvh,kzpp,kypd,kzpd,ndim,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvp
         integer, intent(in) :: kypi, kypp, nxvh, kzpp, kypd, kzpd
         integer, intent(in) :: ndim, nxhyzd, nxyzhd
         real, dimension(ndim,2*nxvh,kypd,kzpd), intent(inout) :: f
         complex, dimension(ndim,nxvh,kzpd) :: ss
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RMNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,&
     &nvpy,nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,ndim,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy
         integer, intent(in) :: nvpz, kxypi, kxypp, nyv, kzpp, kxypd
         integer, intent(in) :: kzpd, ndim, nxhyzd, nxyzhd
         complex, dimension(ndim,nyv,kxypd,kzpd), intent(inout) :: g
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT32RMNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,&
     &nvpy,nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,ndim,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, kstrt, nvpy
         integer, intent(in) :: nvpz, kxypi, kxypp, nzv, kyzp, kxypd
         integer, intent(in) :: kyzpd, ndim, nxhyzd, nxyzhd
         complex, dimension(ndim,nzv,kxypd,kyzpd), intent(inout) :: h
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      end module