!-----------------------------------------------------------------------
!
      module pgraf2
!
! Fortran90 interface to 2d PIC Fortran77 library plibgks2.f
! open_pgraphs open graphics device
!               calls GROPEN and SETNPLT
! reset_pgraphs
!               calls reset_graphs and PDSYNC
! reset graphics device
! close_pgraphs close graphics device
!               calls GRCLOSE
! set_ppalit selects one from three available palettes
!            calls STPALIT
! pdscaler2 displays 2d parallel scalar field in real space
!           calls PCARPET or PCONTUR
! pdvector2 displays 2d parallel vector field in real space
!           calls PCARPET or PCONTUR
! pdgrasp2 displays phase space
!          calls PGRASP23
! pdbgrasp2 displays phase space for magnetized plasma
!           calls PBGRASP23
! ppdgrasp2 displays phase space with segmented particle array
!           calls PPGRASP23
! ppdbgrasp2 displays phase space for magnetized plasma with segmented
!            particle array
!            calls PPGRASP23 or PPBGRASP23
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: february 3, 2018
!
      use plibgks2, only: IPLTCOMM, PGRCLOSE, PDSYNC, PCARPET, PCONTUR, &
     &PGRASP23, PBGRASP23, PPGRASP23, PPBGRASP23
      implicit none
!
! pfs = scratch array for scalar displays
      real, dimension(:,:), allocatable :: pfs
      integer :: szpfs = -1
! plf = scratch array for scalar displays
      integer, dimension(:,:), allocatable :: plf
      integer :: szplf = -1
! pfvs = scratch array for vector displays
      real, dimension(:,:), allocatable :: pfvs
      integer :: szpfvs = -1
! pfp = scratch array for phase space displays
      real, dimension(:), allocatable :: pfp
      integer :: szpfp = -1
      save
!
      private :: pfs, szpfs, plf, szplf, pfvs, szpfvs, pfp, szpfp
!
      contains
!
!-----------------------------------------------------------------------
      function open_pgraphs(nplot) result(irc)
      implicit none
! open graphics device
      integer, intent(in) :: nplot
      integer :: irc
      call GROPEN
      call SETNPLT(nplot,irc)
      end function
!
!-----------------------------------------------------------------------
      subroutine reset_pgraphs(kstrt,irc)
      implicit none
      integer, intent(in) :: kstrt
      integer, intent(inout) :: irc
! reset graphics device
      if (kstrt==1) call RSTSCRN
! synchronize irc and iplot
      call PDSYNC(irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine reset_nplot(nplt,irc)
! resets the maximum number of plots per page
      implicit none
      integer, intent(in) :: nplt
      integer, intent(inout) :: irc
      call SETNPLT(nplt,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine close_pgraphs
      implicit none
! close graphics device
      call PGRCLOSE
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine set_ppalit(idpal)
! selects one from three available palettes
! idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
      implicit none
      integer, intent(in) :: idpal
      call STPALIT(idpal)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pdscaler2(f,nyp,nvp,label,itime,isc,ist,idt,nx,ny,irc)
! displays 2d parallel scalar field in real space
! f = 2d parallel scalar field in real space
! nyp = number of primary gridpoints in field partition
! nvp = number of real or virtual processors requested
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! the range of values of f are given by fmax and fmin.
! if ist = 0, then fmax = 2**isc and fmin = -2**isc.
! if ist = 1, then fmax = 2**isc and fmin = 0.
! if ist = -1, then fmax = 0 and fmin = -2**isc.
! if ist = 2, then fmax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! idt = (1,2,3) = display (color map,contour plot,both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nyp, nvp, itime, isc, ist, idt, nx, ny
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: f
! local data
! ntc = number of valid colors, should be power of 2, <= 256
! nc = number of contour lines
      integer :: ntc = 16, nc = 16
      integer :: nxv, nypmx, lx
      character(len=12) :: lbl
   91 format(' T = ',i7)
      nxv = size(f,1); nypmx = size(f,2)
      lx = nx
! plot guard cells if present
      if ((lx+1) <= nxv) lx = lx + 1
! check if required size of buffer has increased
      if (szpfs < nxv*nypmx) then
         if (szpfs > 0) deallocate(pfs)
! allocate new buffer
         allocate(pfs(nxv,nypmx))
         szpfs = nxv*nypmx
      endif
      write (lbl,91) itime
! color map plot for all values
      if (idt /= 2) then
         call PCARPET(f,pfs,nyp,nvp,label,isc,ist,lx,ny,nxv,nypmx,lbl,  &
     &ntc,irc)
      endif
! contour map for all values
      if (idt /= 1) then
! check if required size of buffer has increased
         if (szplf < nxv*(nypmx+1)) then
            if (szplf > 0) deallocate(plf)
! allocate new buffer
            allocate(plf(nxv,nypmx+1))
            szplf = nxv*(nypmx+1)
         endif
         call PCONTUR(f,pfs,plf,nyp,nvp,label,isc,ist,lx,ny,nxv,nypmx,  &
     &lbl,nc,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pdvector2(fv,nyp,nvp,label,itime,isc,ist,idt,idm,nx,ny,&
     &irc)
! displays 2d parallel vector field in real space
! fv = 2d parallel vector field in real space
! nyp = number of primary gridpoints in field partition
! nvp = number of real or virtual processors requested
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of fv
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! ist = flag for choosing positive and/or negative values
! the range of values of f are given by fmax and fmin.
! if ist = 0, then fmax = 2**isc and fmin = -2**isc.
! if ist = 1, then fmax = 2**isc and fmin = 0.
! if ist = -1, then fmax = 0 and fmin = -2**isc.
! if ist = 2, then fmax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! idt = (1,2,3) = display (color map,contour plot,both)
! idm = (1,2,3) = display (components,sum(abs),both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nyp, nvp, itime, isc, ist, idt, idm, nx, ny
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: fv
! local data
      integer :: i, j, k, lx, nxv, nypmx
      real :: sum1
      character(len=2) :: c
! calculate array size with guard cells
      nxv = size(fv,2); nypmx = size(fv,3)
      lx = nx
! plot guard cells if present
      if ((lx+1) <= nxv) lx = lx + 1
! check if required size of buffer has increased
      if (szpfvs < nxv*nypmx) then
         if (szpfvs > 0) deallocate(pfvs)
! allocate new buffer
         allocate(pfvs(nxv,nypmx))
         szpfvs = nxv*nypmx
      endif
! display components
      if (idm /= 2) then
         do i = 1, size(fv,1)
            do k = 1, nypmx
            do j = 1, lx
               pfvs(j,k) = fv(i,j,k)
            enddo
            enddo
            if (i==1) then
               c = ':X'
            else if (i==2) then
               c = ':Y'
            else if (i==3) then
               c = ':Z'
            else
               write (c,'(":",i1)') i
            endif
! display i component
            call pdscaler2(pfvs,nyp,nvp,label//c,itime,isc,ist,idt,nx,ny&
     &,irc)
            if (irc /= 0) exit
         enddo
      endif
! display sum of absolute values
      if (idm /= 1) then
         do k = 1, nypmx
         do j = 1, lx
            sum1 = 0.0
            do i = 1, size(fv,1)
            sum1 = sum1 + fv(i,j,k)**2
            enddo
            pfvs(j,k) = sqrt(sum1)
         enddo
         enddo
! display amplitude
         call pdscaler2(pfvs,nyp,nvp,label,itime,isc,ist,idt,nx,ny,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pdgrasp2(part,npp,label,itime,isc,nx,ny,iyp,ixp,ierr)
! displays phase space
      implicit none
      integer, intent(in) :: npp, itime, isc, nx, ny, iyp, ixp
      integer, intent(inout) :: ierr
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: part
! local data
      integer :: idimp, npmax
      idimp = size(part,1); npmax = size(part,2)
! check for errors
      ierr = 0
      if ((ixp.lt.1).or.(ixp.gt.idimp)) then
         ierr = 1
      endif
      if ((iyp.lt.1).or.(iyp.gt.idimp)) then
         ierr = ierr + 2
      endif
      if (ierr.gt.0) return
! check if required size of buffer has increased
      if (szpfp < 2*npmax) then
         if (szpfp > 0) deallocate(pfp)
! allocate new buffer
         allocate(pfp(2*npmax))
         szpfp = 2*npmax
      endif
      call PGRASP23(part,pfp,npp,label,itime,isc,nx,ny,iyp,ixp,idimp,   &
     &npmax,ierr)
      end subroutine pdgrasp2
!
!-----------------------------------------------------------------------
      subroutine pdbgrasp2(part,npp,label,itime,isc,omx,omy,omz,nx,ny,  &
     &iyp,ixp,ierr)
! displays phase space for magnetized plasma
      implicit none
      integer, intent(in) :: npp, itime, isc, nx, ny, iyp, ixp
      real, intent(in) :: omx, omy, omz
      integer, intent(inout) :: ierr
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: part
! local data
      integer :: idimp, npmax
      idimp = size(part,1); npmax = size(part,2)
! check for errors
      ierr = 0
      if ((ixp.lt.1).or.(ixp.gt.idimp)) then
         ierr = 1
      endif
      if ((iyp.lt.1).or.(iyp.gt.idimp)) then
         ierr = ierr + 2
      endif
      if (ierr.gt.0) return
! check if required size of buffer has increased
      if (szpfp < 2*npmax) then
         if (szpfp > 0) deallocate(pfp)
! allocate new buffer
         allocate(pfp(2*npmax))
         szpfp = 2*npmax
      endif
      call PBGRASP23(part,pfp,npp,label,itime,isc,omx,omy,omz,nx,ny,iyp,&
     &ixp,idimp,npmax,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ppdgrasp2(ppart,kpic,label,itime,isc,nx,ny,iyp,ixp,ierr&
     &)
! displays phase space with segmented particle array
      implicit none
      integer, intent(in) :: itime, isc, nx, ny, iyp, ixp
      integer, intent(inout) :: ierr
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mxyp1, npmax
      idimp = size(ppart,1); nppmx = size(ppart,2);
      mxyp1 = size(kpic,1)
      npmax = nppmx*mxyp1
! check for errors
      ierr = 0
      if ((ixp.lt.1).or.(ixp.gt.idimp)) then
         ierr = 1
      endif
      if ((iyp.lt.1).or.(iyp.gt.idimp)) then
         ierr = ierr + 2
      endif
      if (ierr.gt.0) return
! check if required size of buffer has increased
      if (szpfp < 2*npmax) then
         if (szpfp > 0) deallocate(pfp)
! allocate new buffer
         allocate(pfp(2*npmax))
         szpfp = 2*npmax
      endif
      call PPGRASP23(ppart,pfp,kpic,label,itime,isc,nx,ny,iyp,ixp,idimp,&
     &nppmx,mxyp1,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ppdbgrasp2(ppart,kpic,label,itime,isc,omx,omy,omz,nx,ny&
     &,iyp,ixp,ierr)
! displays phase space for magnetized plasma with segmented particle
! array
      implicit none
      integer, intent(in) :: itime, isc, nx, ny, iyp, ixp
      integer, intent(inout) :: ierr
      real, intent(in) :: omx, omy, omz
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mxyp1, npmax, nomt
      idimp = size(ppart,1); nppmx = size(ppart,2);
      mxyp1 = size(kpic,1)
      npmax = nppmx*mxyp1
! exit if co-ordinates out of range
      ierr = 0
      if ((ixp.lt.1).or.(ixp.gt.idimp)) then
         ierr = 1
      endif
      if ((iyp.lt.1).or.(iyp.gt.idimp)) then
         ierr = ierr + 2
      endif
      if (ierr.gt.0) return
! determine if magnetic fields points in one of the cartesian directions
      nomt = 0
      if (omx.ne.0.0) nomt = nomt + 1
      if (omy.ne.0.0) nomt = nomt + 1
      if (omz.ne.0.0) nomt = nomt + 1
! check if required size of buffer has increased
      if (szpfp < 2*npmax) then
         if (szpfp > 0) deallocate(pfp)
! allocate new buffer
         allocate(pfp(2*npmax))
         szpfp = 2*npmax
      endif
! magnetic field points in cartesian direction
      if (nomt < 2) then
         call PPGRASP23(ppart,pfp,kpic,label,itime,isc,nx,ny,iyp,ixp,   &
     &idimp,nppmx,mxyp1,ierr)
! magnetic field points in non-cartesian direction
      else
         call PPBGRASP23(ppart,pfp,kpic,label,itime,isc,omx,omy,omz,nx, &
     &ny,iyp,ixp,idimp,nppmx,mxyp1,ierr)
      endif
      end subroutine
!
      end module
