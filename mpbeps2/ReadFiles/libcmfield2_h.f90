!-----------------------------------------------------------------------
! Interface file for cmfield2.f
      module cmfield2_h
      implicit none
!
      interface
         subroutine WRMODES2(pot,pott,nx,ny,modesx,modesy,nxvh,nyv,     &
     &modesxd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, nxvh, nyv
         integer, intent(in) :: modesxd, modesyd
         real, dimension(2*nxvh,nyv), intent(inout) :: pot
         complex, dimension(modesxd,modesyd), intent(in) :: pott
         end subroutine
      end interface
!
      interface
         subroutine WRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,nxvh, &
     &nyv,modesxd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, nxvh, nyv
         integer, intent(in) :: modesxd, modesyd
         real, dimension(ndim,2*nxvh,nyv), intent(inout) :: vpot
         complex, dimension(ndim,modesxd,modesyd), intent(in) :: vpott
         end subroutine
      end interface
!
      interface
         subroutine CSPECT2(fc,wm,pkw,t0,dt,nt,iw,modesx,modesy2,ntd,iwd&
     &,modesxd,modesyd)
         implicit none
         integer, intent(in) :: nt, iw, modesx, modesy2, ntd, iwd
         integer, intent(in) :: modesxd, modesyd
         real, intent(in) :: t0, dt
         complex, dimension(ntd,modesxd,modesyd), intent(in) :: fc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(modesxd,modesyd,iwd,2), intent(inout) :: pkw
         end subroutine
      end interface
!
      interface
         subroutine ICSPECT2(fc,wm,pkw,pks,time,t0,nt,iw,modesx,modesy, &
     &iwd,modesxd,modesyd)
         implicit none
         integer, intent(in) :: nt, iw, modesx, modesy, iwd
         integer, intent(in) :: modesxd, modesyd
         real, intent(in) :: time, t0
         complex, dimension(modesxd,modesyd), intent(in) :: fc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(modesxd,modesyd,iwd,2), intent(inout) :: pkw
         double precision, dimension(4,modesxd,modesyd,iwd),            &
     &intent(inout) :: pks
         end subroutine
      end interface
!
      interface
         subroutine IVCSPECT2(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,    &
     &modesy,ndim,iwd,modesxd,modesyd)
         implicit none
         integer, intent(in) :: nt, iw, modesx, modesy, ndim, iwd
         integer, intent(in) :: modesxd, modesyd
         real, intent(in) :: time, t0
         complex, dimension(ndim,modesxd,modesyd), intent(in) :: fvc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(ndim,modesxd,modesyd,iwd,2), intent(inout) ::  &
     &vpkw
         double precision, dimension(ndim,4,modesxd,modesyd,iwd),       &
     &intent(inout) :: vpks
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: indx, indy, nxhyd, nxyhd
         integer, dimension(nxhyd), intent(inout) :: mixup
         complex, dimension(nxyhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RM2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer , intent(in):: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RM3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer , intent(in):: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer , intent(in):: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer , intent(in):: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RM2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RM2Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      end module