!-----------------------------------------------------------------------
! Interface file for libmpdfield2.f
      module libmpdfield2_h
      implicit none
!
      interface
         subroutine MPPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         real, dimension(nyv,kxp2+1), intent(in)  :: q
         real, dimension(2,nyv,kxp2+1), intent(inout) :: fxy
         complex, dimension(nyd,kxp2), intent(inout) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         real, dimension(nyv,kxp2+1), intent(in)  :: q
         real, dimension(3,nyv,kxp2+1), intent(inout) :: fxy
         complex, dimension(nyd,kxp2), intent(inout) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPCUPERPD2(cu,nx,ny,kstrt,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2
         real, dimension(3,nyv,kxp2+1), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine MIPPBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2, &
     &nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,nyv,kxp2+1), intent(in) :: cu
         real, dimension(3,nyv,kxp2+1), intent(inout) :: bxy
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPMAXWELD2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,  &
     &kstrt,nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(in) :: affp, ci, dt
         real, intent(inout) :: wf, wm
         real, dimension(3,nyv,kxp2+1), intent(inout) :: exy, bxy
         real, dimension(3,nyv,kxp2+1), intent(in)  :: cu
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPEMFIELDR2(fxy,exy,ffd,isign,nx,ny,kstrt,nyv,kxp2,&
     &nyd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, nyd
         real, dimension(3,nyv,kxp2+1), intent(inout) :: fxy
         real, dimension(3,nyv,kxp2+1), intent(in) :: exy
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPBBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2, &
     &nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,nyv,kxp2+1), intent(in) :: cu
         real, dimension(3,nyv,kxp2+1), intent(inout) :: bxy
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPDCUPERPD23(dcu,amu,nx,ny,kstrt,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2
         real, dimension(3,nyv,kxp2+1), intent(inout) :: dcu
         real, dimension(4,nyv,kxp2+1), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MPPADCUPERPD23(dcu,amu,nx,ny,kstrt,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2
         real, dimension(3,nyv,kxp2+1), intent(inout) :: dcu
         real, dimension(4,nyv,kxp2+1), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MPPEPOISD23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf, &
     &nx,ny,kstrt,nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(in) :: ax, ay, affp, wp0, ci
         real, intent(inout) :: wf
         real, dimension(3,nyv,kxp2+1), intent(in) :: dcu
         real, dimension(3,nyv,kxp2+1), intent(inout) :: exy
         complex, dimension(nyd,kxp2), intent(in) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine MPPOTPD2(q,pot,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(inout) :: we
         real, dimension(nyv,kxp2+1), intent(in)  :: q
         real, dimension(nyv,kxp2+1), intent(inout) :: pot
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPELFIELDD22(q,fxy,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(inout) :: we
         real, dimension(nyv,kxp2+1), intent(in)  :: q
         real, dimension(2,nyv,kxp2+1), intent(inout) :: fxy
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPELFIELDD23(q,fxy,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(inout) :: we
         real, dimension(nyv,kxp2+1), intent(in)  :: q
         real, dimension(3,nyv,kxp2+1), intent(inout) :: fxy
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp2
         real, dimension(3,nyv,kxp2+1), intent(in) :: f
         real, dimension(nyv,kxp2+1), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine MPPGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp2
         real, dimension(nyv,kxp2+1), intent(in) :: df
         real, dimension(3,nyv,kxp2+1), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine MPPCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2
         real, dimension(3,nyv,kxp2+1), intent(in) :: f
         real, dimension(3,nyv,kxp2+1), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine MPPAVPOTD23(bxy,axy,nx,ny,kstrt,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2
         real, dimension(3,nyv,kxp2+1), intent(in) :: bxy
         real, dimension(3,nyv,kxp2+1), intent(inout) :: axy
         end subroutine
      end interface
!
      interface
         subroutine MPPAPOTD23(cu,axy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,nyd&
     &)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,nyv,kxp2+1), intent(in) :: cu
         real, dimension(3,nyv,kxp2+1), intent(inout) :: axy
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPETFIELDD23(dcu,exy,ffe,affp,ci,wf,nx,ny,kstrt,nyv&
     &,kxp2,nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, intent(in) :: affp, ci
         real, intent(inout) :: wf
         real, dimension(3,nyv,kxp2+1), intent(in) :: dcu
         real, dimension(3,nyv,kxp2+1), intent(inout) :: exy
         complex, dimension(nyd,kxp2), intent(in) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTHD2(q,qs,ffd,nx,ny,kstrt,nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, dimension(nyv,kxp2+1), intent(in)  :: q
         real, dimension(nyv,kxp2+1), intent(inout) :: qs
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTHD23(cu,cus,ffd,nx,ny,kstrt,nyv,kxp2,nyd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, nyd
         real, dimension(3,nyv,kxp2+1), intent(in) :: cu
         real, dimension(3,nyv,kxp2+1), intent(inout) :: cus
         complex, dimension(nyd,kxp2), intent(in) :: ffd
         end subroutine
      end interface
!
      end module
