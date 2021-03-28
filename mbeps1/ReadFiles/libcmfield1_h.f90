!-----------------------------------------------------------------------
! Interface file for cmfield1.f
      module cmfield1_h
      implicit none
!
      interface
         subroutine WRMODES1(pot,pott,nx,modesx,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, nxvh, modesxd
         real, dimension(3*nxvh), intent(inout) :: pot
         complex, dimension(modesxd), intent(in) :: pott
         end subroutine
      end interface
!
      interface
         subroutine WRVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, ndim, nxvh, modesxd
         real, dimension(ndim,2*nxvh), intent(inout) :: vpot
         complex, dimension(ndim,modesxd), intent(in) :: vpott
         end subroutine
      end interface
!
      interface
         subroutine CSPECT1(fc,wm,pkw,t0,dt,nt,iw,modesx,ntd,iwd,modesxd&
     &)
         implicit none
         integer, intent(in) :: nt, iw, modesx, ntd, iwd, modesxd
         real, intent(in) :: t0, dt
         complex, dimension(ntd,modesxd), intent(in) :: fc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(modesxd,iwd,2), intent(inout) :: pkw
         end subroutine
      end interface
!
      interface
         subroutine ICSPECT1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm,&
     &iwd,modesxd)
         implicit none
         integer, intent(in) :: nt, iw, modesx, nx, norm, iwd, modesxd
         real, intent(in) :: time, t0
         complex, dimension(modesxd), intent(in) :: fc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(modesxd,iwd,2), intent(inout) :: pkw
         double precision, dimension(4,modesxd,iwd), intent(inout) ::   &
     &pks
         end subroutine
      end interface
!
      interface
         subroutine IVCSPECT1(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,nx, &
     &norm,iwd,modesxd)
         implicit none
         integer, intent(in) :: nt, iw, modesx, nx, norm, iwd, modesxd
         real, intent(in) :: time, t0
         complex, dimension(2,modesxd), intent(in) :: fvc
         real, dimension(iwd), intent(in) :: wm
         real, dimension(2,modesxd,iwd,2), intent(inout) :: vpkw
         double precision, dimension(2,4,modesxd,iwd), intent(inout) :: &
     &vpks
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
      interface
         subroutine DIVF1(f,df,nx,ndim,nxvh)
         implicit none
         integer, intent(in) :: nx, ndim, nxvh
         real, dimension(ndim,nxvh), intent(in) :: f
         real, dimension(nxvh), intent(inout) :: df
         end subroutine
      end interface
!
      end module
