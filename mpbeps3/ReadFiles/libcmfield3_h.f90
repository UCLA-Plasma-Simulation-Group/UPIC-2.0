!-----------------------------------------------------------------------
! Interface file for cmfield3.f
      module cmfield3_h
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
         subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
         implicit none
         integer , intent(in):: indx, indy, indz, nxhyzd, nxyzhd
         integer, dimension(nxhyzd), intent(inout) :: mixup
         complex, dimension(nxyzhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RMX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd, &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RM3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd, &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer , intent(in):: nxhyzd, nxyzhd
         real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
          subroutine FFT3RMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,  &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer , intent(in):: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,  &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      end module