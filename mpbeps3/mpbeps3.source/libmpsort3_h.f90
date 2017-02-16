!-----------------------------------------------------------------------
! Interface file for libmpsort3.f
      module libmpsort3_h
      implicit none
!
      interface
         subroutine PPPORDER32LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole&
     &,ncll,nclr,noff,nyzp,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mzp1, &
     &mxzyp1,npbmx,ntmax,nbmax,idds,irc2)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: mx1, myp1, mzp1, mxzyp1, npbmx, ntmax
         integer, intent(in) :: nbmax, idds
         real, dimension(idimp,nppmx,mx1*myp1*mzp1), intent(inout) ::   &
     &ppart
         real, dimension(idimp,npbmx,mx1*myp1*mzp1), intent(inout) ::   &
     &ppbuff
         real, dimension(idimp,nbmax,2), intent(inout) :: sbufl, sbufr
         integer, dimension(mx1*myp1*mzp1), intent(in) :: kpic
         integer, dimension(26,mx1*myp1*mzp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1*mzp1), intent(inout) ::  &
     &ihole
         integer, dimension(3,mxzyp1,3,2), intent(inout) :: ncll, nclr
         integer, dimension(idds), intent(in) :: noff, nyzp
         integer, dimension(2), intent(inout) :: irc2
         end subroutine
      end interface
!
      interface
         subroutine PPPORDERF32LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,   &
     &ncll,nclr,idimp,nppmx,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax,nbmax,irc2)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, mzp1, mxzyp1
         integer, intent(in) :: npbmx, ntmax, nbmax
         real, dimension(idimp,nppmx,mx1*myp1*mzp1), intent(inout) ::   &
     &ppart
         real, dimension(idimp,npbmx,mx1*myp1*mzp1), intent(inout) ::   &
     &ppbuff
         real, dimension(idimp,nbmax,2), intent(inout) :: sbufl, sbufr
         integer, dimension(26,mx1*myp1*mzp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1*mzp1), intent(in) ::     &
     &ihole
         integer, dimension(3,mxzyp1,3,2), intent(inout) :: ncll, nclr
         integer, dimension(2), intent(inout) :: irc2
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER32LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole&
     &,mcll,mclr,mcls,idimp,nppmx,nx,ny,nz,mx1,myp1,mzp1,mxzyp1,npbmx,  &
     &ntmax,nbmax,irc2)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx1, myp1
         integer, intent(in) :: mzp1, mxzyp1, npbmx, ntmax, nbmax
         real, dimension(idimp,nppmx,mx1*myp1*mzp1), intent(inout) ::   &
     &ppart
         real, dimension(idimp,npbmx,mx1*myp1*mzp1), intent(in) ::      &
     &ppbuff
         real, dimension(idimp,nbmax,2), intent(in) :: rbufl, rbufr
         integer, dimension(mx1*myp1*mzp1), intent(inout) :: kpic
         integer, dimension(26,mx1*myp1*mzp1), intent(in) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1*mzp1), intent(in) ::     &
     &ihole
         integer, dimension(3,mxzyp1,3,2), intent(in) :: mcll, mclr
         integer, dimension(3,mx1+1,4), intent(in) :: mcls
         integer, dimension(2), intent(inout) :: irc2
         end subroutine
      end interface
!
      interface
         subroutine PPPRSNCL3L(ncl,mxzyp1)
         implicit none
         integer, intent(in) :: mxzyp1
         integer, dimension(26,mxzyp1), intent(inout) :: ncl
         end subroutine
      end interface
!
      interface
         subroutine PPPORDERF32LAF(ppbuff,sbufl,sbufr,ncl,ncll,nclr,    &
     &idimp,mx1,myp1,mzp1,mxzyp1,npbmx,nbmax,irc2)
         implicit none
         integer, intent(in) :: idimp, mx1, myp1, mzp1, mxzyp1
         integer, intent(in) :: npbmx, nbmax
         real, dimension(idimp,npbmx,mx1*myp1*mzp1), intent(inout) ::   &
     &ppbuff
         real, dimension(idimp,nbmax,2), intent(inout) :: sbufl, sbufr
         integer, dimension(26,mx1*myp1*mzp1), intent(inout) :: ncl
         integer, dimension(3,mxzyp1,3,2), intent(inout) :: ncll, nclr
         integer, dimension(2), intent(inout) :: irc2
         end subroutine
      end interface
!
      interface
         subroutine PPPRSTOR3L(ppart,ppbuff,ncl,ihole,idimp,nppmx,mxyzp1&
     &,npbmx,ntmax)
         implicit none
         integer, intent(in) :: idimp, nppmx, mxyzp1, npbmx, ntmax
         real, dimension(idimp,nppmx,mxyzp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mxyzp1), intent(in) :: ppbuff
         integer, dimension(26,mxyzp1), intent(in) :: ncl
         integer, dimension(2,ntmax+1,mxyzp1), intent(in) :: ihole
         end subroutine
      end interface
!
      end module