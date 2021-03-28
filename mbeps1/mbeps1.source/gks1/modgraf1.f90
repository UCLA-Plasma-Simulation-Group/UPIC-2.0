!-----------------------------------------------------------------------
!
      module graf1
!
! Fortran90 interface to 1d PIC Fortran77 library libgks1.f
! open_graphs open graphics device
!             calls GROPEN and SETNPLT
! close_graphs close graphics device
!              calls GRCLOSE
! reset_graphs reset graphics device
!              calls RSTSCRN
! reset_nplot resets the maximum number of plots per page
!             calls SETNPLT
! dscaler1 displays 1d scalar field in real space.
!          calls DISPS
! dmscaler1 displays multiple 1d scalar field in real space.
!           calls DISPR
! dvector1 displays amplitude or components of a 1d vector field in real
!          space
!          calls DISPS
! dmvector1 displays amplitude  or components of multiple 1d vector
!           field in real space
!           calls DISPR
! displayfv1 displays velocity distribution functions.
!            calls DISPR
! displayfe1 displays energy distribution functions
!            calls DISPR
! displayfvb1 displays velocity distribution functions in cylindrical
!             co-ordinates
!             calls DISPR
! displayfvt1 displays time history of vdrift, vth, and entropy 
!             calls DISPR
! displaytr1 displays time history of trajectories
!            calls DISPR
! displayw1 displays time history of electric field,
!           kinetic and total energies.
!           calls DISPR
! dgrasp1 displays particle in phase space.
!         calls GRASP13
! dbgrasp1 displays phase space for magnetized plasma
!          calls BGRASP13
! dpmgrasp1 displays phase space with marked or unmarked particles
!           calls PMGRASP1 or PGRASP13
! dpmbgrasp1 displays phase space with marked or unmarked particles
!            for magnetized plasma
!            calls PMGRASP13, PGRASP13, PMBGRASP13, or PBGRASP13
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 20, 2020
!
      use libgraf1_h
      implicit none
!
! npl = (0,1) = display is (off,on)
      integer, save :: npl = 1
!
      contains
!
!-----------------------------------------------------------------------
      function open_graphs(nplot) result(irc)
! open graphics device
      integer, intent(in) :: nplot
      integer :: irc
      call GROPEN
      call SETNPLT(nplot,irc)
      end function
!
!-----------------------------------------------------------------------
      subroutine close_graphs
! close graphics device
      call GRCLOSE
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine reset_graphs
! reset graphics device
      call RSTSCRN
      npl = 1
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
      subroutine dscaler1(f,label,itime,isc,ist,nx,irc)
! displays 1d scalar field in real space
! f = 1d scalar field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! ist = flag for choosing positive and/or negative values
! the plot has a scale in y given by ymax and ymin.
! if ist = 0, then ymax = 2**isc and ymin = -2**isc.
! if ist = 1, then ymax = 2**isc and ymin = 0.
! if ist = -1, then ymax = 0 and ymin = -2**isc.
! if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! nx = system length in x direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, ist, nx
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:), intent(inout) :: f
! local data
      integer :: nxv, lx
      real :: xmin, xmax
      character(len=12) :: lbl
   91 format(' T = ',i7)
      if (npl==0) return
      nxv = size(f)
      xmin = 0.0
      write (lbl,91) itime
      lx = min(nx+1,nxv); xmax = real(lx - 1)
      call DISPS(f(1),label,xmin,xmax,isc,ist,lx,lbl,irc)
      if (irc > 127) then
         npl = irc - 128
         if (npl==0) call CLRSCRN
         irc = 0
      endif
      end subroutine
!
!
!-----------------------------------------------------------------------
      subroutine dmscaler1(f,label,itime,isc,ist,nx,chrs,irc)
! displays multiple 1d scalar fields in real space
! f = 1d scalar field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! ist = flag for choosing positive and/or negative values
! the plot has a scale in y given by ymax and ymin.
! if ist = 0, then ymax = 2**isc and ymin = -2**isc.
! if ist = 1, then ymax = 2**isc and ymin = 0.
! if ist = -1, then ymax = 0 and ymin = -2**isc.
! if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! nx = system length in x direction
! chrs = array of short labels for each field
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, ist, nx
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(inout) :: f
      character(len=10), dimension(:), intent(in) :: chrs
! local data
      integer :: num, nxv, lx
      real :: xmin, xmax
      character(len=12) :: lbl
   91 format(' T = ',i7)
      if (npl==0) return
      nxv = size(f,1); num = size(f,2)
      xmin = 0.0
      write (lbl,91) itime
      lx = min(nx+1,nxv); xmax = real(lx - 1)
      call DISPR(f(1,1),label,xmin,xmax,isc,ist,0,nx,nxv,num,lbl,chrs,  &
     &irc)
      if (irc > 127) then
         npl = irc - 128
         if (npl==0) call CLRSCRN
         irc = 0
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dvector1(f,label,itime,isc,ist,idm,nx,irc)
! displays 1d vector field in real space
! f = 1d vector field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! ist = flag for choosing positive and/or negative values
! the plot has a scale in y given by ymax and ymin.
! if ist = 0, then ymax = 2**isc and ymin = -2**isc.
! if ist = 1, then ymax = 2**isc and ymin = 0.
! if ist = -1, then ymax = 0 and ymin = -2**isc.
! if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! idm = (1,2,3,4) = display (each component,composite,both,amplitude)
! nx = system length in x direction
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, ist, idm, nx
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: f
! local data
      character(len=10), dimension(size(f,1)) :: chrs
      real, dimension(size(f,2),size(f,1)) :: g
      integer :: i, j, ndim, nxv, lx
      real :: xmin, xmax, sum1
      character(len=12) :: lbl
      character(len=2) :: c
   91 format(' T = ',i7)
      if (npl==0) return
      write (lbl,91) itime
      ndim = size(f,1)
      nxv = size(f,2)
      xmin = 0.0
! transpose array and set labels
      do j = 1, nxv
         do i = 1, ndim
            g(j,i) = f(i,j)
         enddo
      enddo
! set labels
      do i = 1, ndim
         if (i==1) then
            chrs(1) = '    Y     '
         else if (i==2) then
            chrs(2) = '    Z     '
         else
            write (chrs(i),'(4x,i1)') i
         endif
      enddo
! display each component separately
      if (idm /= 2) then
         do i = 1, ndim
         if (i==1) then
            c = ':Y'
         else if (i==2) then
            c = ':Z'
         else
            write (c,'(":",i1)') i
         endif
         lx = min(nx+1,nxv); xmax = real(lx - 1)
         call DISPS(g(1,i),label//c,xmin,xmax,isc,ist,lx,lbl,irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
            return
         endif
         if (irc==1) return
         enddo
      endif
! display components on a common plot
      if (idm /= 1) then
         lx = min(nx+1,nxv); xmax = real(lx - 1)
         call DISPR(g(1,1),label,xmin,xmax,isc,ist,0,lx,nxv,ndim,lbl,   &
     &chrs,irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
         endif
      endif
! display amplitude
      if (idm==4)then
         do j = 1, lx
            sum1 = 0.0
            do i = 1, size(f,1)
            sum1 = sum1 + f(i,j)**2
            enddo
            g(j,1) = sqrt(sum1)
         enddo
         lx = min(nx+1,nxv); xmax = real(lx - 1)
         call DISPS(g(1,1),label,xmin,xmax,isc,ist,lx,lbl,irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dmvector1(f,label,itime,isc,ist,idm,nx,chrs,irc)
! displays multiple 1d vector field in real space
! f = 1d vector field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! ist = flag for choosing positive and/or negative values
! the plot has a scale in y given by ymax and ymin.
! if ist = 0, then ymax = 2**isc and ymin = -2**isc.
! if ist = 1, then ymax = 2**isc and ymin = 0.
! if ist = -1, then ymax = 0 and ymin = -2**isc.
! if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! idm = (1,2,3) = display (components,average,both)
! nx = system length in x direction
! chrs = array of short labels for each field
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, ist, idm, nx
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: f
      character(len=10), dimension(:), intent(in) :: chrs
! local data
      real, dimension(size(f,2),size(f,3)) :: g
      integer :: i, j, k, nv, num, nxv, nnxv, lx
      real :: anv, xmin, xmax
      character(len=12) :: lbl
      character(len=2) :: c
   91 format(' T = ',i7)
      if (npl==0) return
      write (lbl,91) itime
      nv = size(f,1)
      nxv = size(f,2); num = size(f,3)
      nnxv = num*nxv
      anv = 1.0/real(nv)
      xmin = 0.0
! display components
      if (idm /= 2) then
         do i = 1, nv
            g = f(i,:,:)
         if (i==1) then
            c = ':Y'
         else if (i==2) then
            c = ':Z'
         else
            write (c,'(":",i1)') i
         endif
         lx = min(nx+1,nxv); xmax = real(lx - 1)
         call DISPR(g(1,1),label//c,xmin,xmax,isc,ist,0,lx,nxv,num,lbl, &
     &chrs,irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
            return
         endif
         if (irc==1) return
         enddo
      endif
! display sum of absolute values
      if (idm /= 1) then
         lx = min(nx+1,nxv)
         do k = 1, size(f,3)
            do j = 1, lx
! display sum of absolutes
               g(j,k) = f(1,j,k)
               do i = 2, nv
                  g(j,k) = g(j,k) + f(i,j,k)
               enddo
               g(j,k) = g(j,k)*anv
            enddo
         enddo
         lx = min(nx+1,nxv); xmax = real(lx - 1)
         call DISPR(g(1,1),label,xmin,xmax,isc,ist,0,lx,nxv,num,lbl,chrs&
     &,irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine displayfv1(fv,fvm,label,itime,nmv,idt,irc)
! displays velocity distribution functions
! fv = velocity distribution
! fvm = velocity moments
! label = long character string label for plot
! itime = current time step
! nmv = number of velocity intervals
! idt = (1,2,3) = display (individual,composite,both) functions
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, nmv, idt
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(inout) :: fv
      real, dimension(:,:), intent(in) :: fvm
      character(len=*), intent(in) :: label
! isc = 999 = use default display scale
! ist = 1 = display positive values only
! mks = 0 = cycle through line styles
! local data
      integer :: isc = 999, ist = 1, mks = 0
      integer :: i, nmvf, nmv21, idimv
      real :: vmax, vmin
      character(len=12) :: c
      character(len=2) :: cs
      character(len=54) :: lbl
      character(len=45) :: chr
      character(len=10), dimension(3) :: chrs
   91 format(', T =',i7)
   92 format(' VD =',f9.6,' VTH =',f9.5)
   93 format(' VTX =',f9.5)
   94 format(' VTX =',f9.5,' VTY =',f9.5,' VTZ =',f9.5)
! chrs = short array of characters to label individual line samples
      data chrs /'    VX    ','    VY    ','    VZ    '/
      if (npl==0) return
      idimv = size(fv,2)
      if ((idimv /= 1) .and. (idimv /= 3)) return
      nmvf = size(fv,1)
      nmv21 = 2*nmv + 1
      write (c,91) itime
! each velocity distributions on its own plot
      if (idt /= 2) then
         do i = 1, idimv
         cs = trim(adjustl(chrs(i)))
         lbl = trim(label)//' VELOCITY DISTR VS '//cs//c
         write (chr,92) fvm(i,1), fvm(i,2)
         vmax = fv(nmv21+1,i)
         vmin = -vmax
         call DISPR(fv(1,i),lbl,vmin,vmax,isc,ist,mks,nmv21,nmvf,1,chr, &
     &chrs(i),irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
            return
         endif
         if (irc==1) return
         enddo
      endif
! all velocity distributions on common plot
      if (idt /= 1) then
         lbl = trim(label)//' VELOCITY DISTRS VS '//'V'//c
         if (idimv==1) then
            write (chr,93) fvm(1,2)
            vmax = fv(nmv21+1,1)
            vmin = -vmax
         else if (idimv==3) then
            write (chr,94) fvm(1,2), fvm(2,2), fvm(3,2)
            vmax = max(fv(nmv21+1,1),fv(nmv21+1,2),fv(nmv21+1,3))
            vmin = -vmax
         endif
         call DISPR(fv(1,1),lbl,vmin,vmax,isc,ist,mks,nmv21,nmvf,idimv, &
     &chr,chrs,irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine displayfe1(fe,wk,label,itime,nmv,irc)
! displays energy distribution functions
! fe = energy distribution
! wk = total energy contained in distribution
! label = long character string label for plot
! itime = current time step
! nmv = number of velocity intervals
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, nmv
      integer, intent(inout) :: irc
      real, intent(in) :: wk
      real, dimension(:,:), intent(inout) :: fe
      character(len=*), intent(in) :: label
! isc = 999 = use default display scale
! ist = 1 = display positive values only
! mks = 0 = cycle through line styles
! local data
      integer :: isc = 999, ist = 1, mks = 0
      integer :: nmvf, nmv21
      real :: emax, emin
      character(len=12) :: c
      character(len=54) :: lbl
      character(len=45) :: chr
      character(len=10), dimension(1) :: chrs
   91 format(', T =',i7)
   92 format(' TOTAL ENERGY =',e14.7)
! chrs = short array of characters to label individual line samples
      data chrs /'  ENERGY  '/
      if (npl==0) return
      nmvf = size(fe,1)
      nmv21 = 2*nmv + 1
      write (c,91) itime
      lbl = trim(label)//' ENERGY DISTR VS ENERGY'//c
      write (chr,92) wk
      emax = 2.0*fe(nmv21+1,1)
      emin = 0.0
      call DISPR(fe(1,1),lbl,emin,emax,isc,ist,mks,nmv21,nmvf,1,chr,    &
     &chrs(1),irc)
      if (irc > 127) then
         npl = irc - 128
         if (npl==0) call CLRSCRN
         irc = 0
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine displayfvb1(fv,fvm,label,itime,nmv,idt,irc)
! displays velocity distribution functions
! fv = velocity distribution
! fvm = velocity moments
! label = long character string label for plot
! itime = current time step
! nmv = number of velocity intervals
! idt = (1,2,3) = display (individual,composite,both) functions
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, nmv, idt
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(inout) :: fv
      real, dimension(:,:), intent(in) :: fvm
      character(len=*), intent(in) :: label
! isc = 999 = use default display scale
! ist = 1 = display positive values only
! mks = 0 = cycle through line styles
! local data
      integer :: isc = 999, ist = 1, mks = 0
      integer :: i, nmvf, nmv21, idimv
      real :: vmax, vmin
      character(len=12) :: c
      character(len=3) :: cs
      character(len=54) :: lbl
      character(len=45) :: chr
      character(len=5), dimension(2) :: chrs
   91 format(', T =',i7)
   92 format(' VD =',f9.6,' VTH =',f9.5)
   93 format(' VTR =',f9.5,' VTL =',f9.5)
! chrs = short array of characters to label individual line samples
      data chrs /' VPR ',' VPL '/
      if (npl==0) return
      idimv = 2
      nmvf = size(fv,1)
      nmv21 = 2*nmv + 1
      write (c,91) itime
! each velocity distributions on its own plot
      if (idt /= 2) then
         do i = 1, idimv
         cs = trim(adjustl(chrs(i)))
         lbl = trim(label)//' VELOCITY DISTR VS '//cs//c
         write (chr,92) fvm(i,1), fvm(i,2)
         vmax = fv(nmv21+1,i)
         vmin = -vmax
         call DISPR(fv(1,i),lbl,vmin,vmax,isc,ist,mks,nmv21,nmvf,1,chr, &
     &chrs(i),irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
            return
         endif
         if (irc==1) return
         enddo
      endif
! all velocity distributions on common plot
      if (idt /= 1) then
         lbl = trim(label)//' VELOCITY DISTRS VS '//'V'//c
         write (chr,93) fvm(1,2), fvm(2,2)
         vmax = max(fv(nmv21+1,1),fv(nmv21+1,2))
         vmin = -vmax
         call DISPR(fv(1,1),lbl,vmin,vmax,isc,ist,mks,nmv21,nmvf,idimv, &
     &chr,chrs,irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine displaymfv1(fv,fvm,label,itime,nmv,idt,irc)
! displays multiple velocity distribution functions
! fv = velocity distribution
! fvm = velocity moments
! label = long character string label for plot
! itime = current time step
! nmv = number of velocity intervals
! idt = (1,2,3) = display (individual,composite,both) functions
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, nmv, idt
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(inout) :: fv
      real, dimension(:,:,:), intent(in) :: fvm
      character(len=*), intent(in) :: label
! isc = 999 = use default display scale
! ist = 1 = display positive values only
! mks = 0 = cycle through line styles
! local data
      integer :: isc = 999, ist = 1, mks = 0
      integer :: i, nmvf, nmv2, idimv, ngs, ngsc
      real :: vmax, vmin
      character(len=12) :: c
      character(len=2) :: cs
      character(len=54) :: lbl
      character(len=45) :: chr
      character(len=10), dimension(3) :: chrs
   91 format(', T =',i7)
   92 format(' VD =',f9.6,' VTH =',f9.5)
   93 format(' VTX =',f9.5)
   94 format(' VTX =',f9.5,' VTY =',f9.5,' VTZ =',f9.5)
! chrs = short array of characters to label individual line samples
      data chrs /'    VX    ','    VY    ','    VZ    '/
      if (npl==0) return
      ngs = size(fv,3)
      idimv = size(fv,2)
      ngsc = idimv*ngs
      if ((idimv /= 1) .and. (idimv /= 3)) return
      nmvf = size(fv,1)*idimv
      nmv2 = 2*nmv + 1
      write (c,91) itime
! each velocity distributions on its own plot
      if (idt /= 2) then
         do i = 1, idimv
         cs = trim(adjustl(chrs(i)))
         lbl = trim(label)//' VELOCITY DISTR VS '//cs//c
         write (chr,92) fvm(1,i,1), fvm(2,i,1)
         vmax = fv(1,i,1)
         vmin = -vmax
         call DISPR(fv(2,i,1),lbl,vmin,vmax,isc,ist,mks,nmv2,nmvf,ngs,  &
     &chr,chrs(i),irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
            return
         endif
         if (irc==1) return
         enddo
      endif
! all velocity distributions on common plot
      if (idt /= 1) then
         lbl = trim(label)//' VELOCITY DISTRS VS '//'V'//c
         if (idimv==1) then
            write (chr,93) fvm(2,1,1)
            vmax = fv(1,1,1)
            vmin = -vmax
         else if (idimv==3) then
            write (chr,94) fvm(2,1,1), fvm(2,2,1), fvm(2,3,1)
            vmax = max(fv(1,1,1),fv(1,2,1),fv(1,3,1))
            vmin = -vmax
         endif
         call DISPR(fv(2,1,1),lbl,vmin,vmax,isc,ist,mks,nmv2,nmvf,ngsc, &
     &chr,chrs,irc)
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine displayfvt1(fvtm,label,t0,dtv,nt,irc)
! displays time history of vdrift, vth, and entropy 
! fvtm = time history array for velocity values
! t0 = initial velocity time
! dtv = time between velocity values
! nt = number of velocity values to be displayed
! irc = return code (0 = normal return)
      integer, intent(in) :: nt
      integer, intent(inout) :: irc
      real, intent(in) :: t0, dtv
      real, dimension(:,:,:), intent(inout) :: fvtm
      character(len=*), intent(in) :: label
! isc = 999 = use default display scale
! ist = 2 = display minimum range
! mks = 0 = cycle through line styles
! local data
      real, dimension(size(fvtm,1)) :: g
      integer :: isc = 999, ist = 2, mks = 0
      integer :: i, j, ntvd, ns, idimv
      real :: tmin, tmax
      character(len=36) :: lbl
      character(len=8), dimension(3) :: cs
      character(len=4), dimension(3) :: chrs
      data cs /'VDRIFT  ','VTHERMAL','ENTROPY '/
      data chrs /' VX ',' VY ',' VZ '/
      if (npl==0) return
      idimv = size(fvtm,2)
! quit if array is empty or incorrect
      if (nt <= 0) return
      ntvd = size(fvtm,1)
      ns = min(size(fvtm,3),3)
! tmin/tmax = range of time values in plot
      tmin = t0
      tmax = t0 + dtv*(nt - 1)
! display individual moments
      do i = 1, ns-1
! loop over components
      do j = 1, idimv
      lbl = trim(label)//chrs(j)//cs(i)//' VERSUS TIME'
      call DISPR(fvtm(1,j,i),lbl,tmin,tmax,isc,ist,mks,nt,ntvd,1,' ',   &
     &cs(i),irc)
      if (irc==1) return
      enddo
      enddo
! display sum of entropies
      g = fvtm(:,1,3)
      if (idimv==3) g = g + fvtm(:,2,3) + fvtm(:,3,3)
      lbl = trim(label)//' '//cs(3)//' VERSUS TIME'
      call DISPR(g(1),lbl,tmin,tmax,isc,ist,mks,nt,ntvd,1,' ',cs(3),irc)
      if (idimv==1) return
! all drifts on common plot
      lbl = trim(label)//' VDRIFTS VERSUS TIME'
      call DISPR(fvtm(1,1,1),lbl,tmin,tmax,isc,ist,mks,nt,ntvd,idimv,' '&
     &,chrs,irc)
! all thermal velocities on common plot
      lbl = trim(label)//' VTHERMALS VERSUS TIME'
      call DISPR(fvtm(1,1,2),lbl,tmin,tmax,isc,ist,mks,nt,ntvd,idimv,' '&
     &,chrs,irc)
! all thermal entropies on common plot
      lbl = trim(label)//' ENTROPIES VERSUS TIME'
      call DISPR(fvtm(1,1,3),lbl,tmin,tmax,isc,ist,mks,nt,ntvd,idimv,' '&
     &,chrs,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine displaytr1(pt,t0,dtt,nt,it,isc,irc)
! displays time history of trajectories
! pt = time history array for trajectories
! t0 = initial trajectory time
! dtt = time between trajectories values
! nt = number of trajectory values to be displayed
! it = (1,2,3,4) plot (x,vx,vy,vz)
! isc = power of 2 scale of range of values of velocity
! irc = return code (0 = normal return)
      integer, intent(in) :: nt, it, isc
      integer, intent(inout) :: irc
      real, intent(in) :: t0, dtt
      real, dimension(:,:,:), intent(inout) :: pt
! ist = 0 = display maximum range
! mks = 0 = cycle through line styles
!     local data
      integer :: ist = 0, mks = 0
      integer :: i, idimp, nttd, ns
      real :: tmin, tmax
      character(len=36) :: lbl
      character(len=4), dimension(4) :: cs 
      character(len=10), dimension(size(pt,3)) :: chrs
! chrs = short array of characters to label individual line samples
      data cs /' X  ',' VX ',' VY ',' VZ '/
      if (npl==0) return
! quit if array is empty or incorrect
      if (nt <= 0) return
      idimp = size(pt,2)
! exit if co-ordinate out of range
      if (it.ge.idimp) return
      nttd = size(pt,1)*idimp
      ns = size(pt,3)
! tmin/tmax = range of time values in plot
      tmin = t0
      tmax = t0 + dtt*(nt - 1)
! display individual trajectories
!     lbl = cs(it)//' VERSUS TIME'
      do i = 1, ns
      write (chrs(i),'(i10)') i
      chrs(i) = '    '//adjustl(chrs(i))
!     call DISPR(pt(1,it,1),lbl,tmin,tmax,isc,ist,mks,nt,nttd,1,' ',    &
!    &chrs(i),irc)
!     if (irc==1) return
      enddo
! all trajectories on common plot
      lbl = cs(it)//' VERSUS TIME'
      call DISPR(pt(1,it,1),lbl,tmin,tmax,isc,ist,mks,nt,nttd,ns,       &
     &' TEST PARTICLES',chrs,irc)
       end subroutine
!
!-----------------------------------------------------------------------
      subroutine displayw1(wt,t0,dtw,nt,irc)
! displays time history of electric field, kinetic, and
! total energies
! wt = time history array for energies
! t0 = initial energy time
! dtw = time between energy values
! nt = number of energy values to be displayed
! irc = return code (0 = normal return)
      integer, intent(in) :: nt
      integer, intent(inout) :: irc
      real, intent(in) :: t0, dtw
      real, dimension(:,:), intent(inout) :: wt
! isc = 999 = use default display scale
! ist = 2 = display minimum range
! mks = 0 = cycle through line styles
! local data
      integer :: isc = 999, ist = 2, mks = 0
      integer :: i, ntwd, ns
      real :: tmin, tmax
      character(len=36) :: lbl
      character(len=20), dimension(7) :: cs 
      character(len=10), dimension(7) :: chrs
! chrs = short array of characters to label individual line samples
      data cs /' TOTAL FIELD        ',' ELECTRON KINETIC   ',           &
     &' ION KINETIC     ',' TOTAL              ',' ES FIELD           ',&
     &' ET FIELD        ',' MAGNETIC FIELD     '/
      data chrs /'TOT FIELD ','ELECT ENRG',' ION ENRG ','TOTAL ENRG',   &
     &' EL FIELD ',' ET FIELD ',' B FIELD  '/
      if (npl==0) return
! quit if array is empty or incorrect
      if (nt <= 0) return
      ntwd = size(wt,1)
      ns = min(size(wt,2),7)
! tmin/tmax = range of time values in plot
      tmin = t0
      tmax = t0 + dtw*(nt - 1)
! display individual energies
      do i = 1, ns
      lbl = trim(cs(i))//' ENERGY VERSUS TIME'
      call DISPR(wt(1,i),lbl,tmin,tmax,isc,ist,mks,nt,ntwd,1,' ',chrs(i)&
     &,irc)
      if (irc==1) return
      enddo
! all energies on common plot
      lbl = ' ENERGIES VERSUS TIME'
      call DISPR(wt(1,1),lbl,tmin,tmax,isc,ist,mks,nt,ntwd,ns,' ',chrs, &
     &irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dgrasp1(part,np,label,itime,isc,nx,iyp,ixp,npx,irc)
! displays phase space
      implicit none
      integer, intent(in) :: np, itime, isc, nx, iyp, ixp, npx
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: part
! local data
      integer :: idimp
      idimp = size(part,1)
! exit if co-ordinates out of range
      if ((ixp.gt.idimp).or.(iyp.gt.idimp)) return
      call GRASP13(part,label,itime,isc,nx,iyp,ixp,idimp,npx,np,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dbgrasp1(part,np,label,itime,isc,omx,omy,omz,nx,iyp,ixp&
     &,npx,irc)
! displays phase space for magnetized plasma
      implicit none
      integer, intent(in) :: np, itime, isc, nx, iyp, ixp, npx
      integer, intent(inout) :: irc
      real, intent(in) :: omx, omy, omz
      character(len=*), intent(in) :: label
      real, dimension(:,:), intent(in) :: part
! local data
      integer :: idimp
      idimp = size(part,1)
! exit if co-ordinates out of range
      if ((ixp.gt.idimp).or.(iyp.gt.idimp)) return
      call BGRASP13(part,label,itime,isc,omx,omy,omz,nx,iyp,ixp,idimp,  &
     &npx,np,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dpmgrasp1(ppart,kpic,label,itime,isc,nx,iyp,ixp,ntsc,  &
     &irc)
! displays phase space for with marked or unmarked particles
! ntsc = (0,1) = (no,yes) color beam particles
      implicit none
      integer, intent(in) :: itime, isc, nx, iyp, ixp, ntsc
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mx1, ltag
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! exit if co-ordinates out of range
      if ((ixp.gt.idimp).or.(iyp.gt.idimp)) return
! plot marked particles with color
      if (ntsc > 0) then
         ltag = idimp
         call PMGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp,    &
     &nppmx,mx1,ltag,irc)
! plot unmarked particles with no colors
      else
         call PGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp,nppmx&
     &,mx1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dpmbgrasp1(ppart,kpic,label,itime,isc,omx,omy,omz,nx,  &
     &iyp,ixp,ntsc,irc)
! displays phase space for with marked or unmarked particles
! for magnetized plasma
! ntsc = (0,1) = (no,yes) color beam particles
      implicit none
      integer, intent(in) :: itime, isc, nx, iyp, ixp, ntsc
      integer, intent(inout) :: irc
      real, intent(in) :: omx, omy, omz
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, nppmx, mx1, nomt, ltag
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! exit if co-ordinates out of range
      if ((ixp.gt.idimp).or.(iyp.gt.idimp)) return
! determine if magnetic fields points in one of the cartesian directions
      nomt = 0
      if (omx.ne.0.0) nomt = nomt + 1
      if (omy.ne.0.0) nomt = nomt + 1
      if (omz.ne.0.0) nomt = nomt + 1
! magnetic field points in cartesian direction
      if (nomt < 2) then
! plot marked particles with color
         if (ntsc > 0) then
            ltag = idimp
            call PMGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp, &
     &nppmx,mx1,ltag,irc)
! plot unmarked particles with no colors
         else
            call PGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp,  &
     &nppmx,mx1,irc)
         endif
! magnetic field points in non-cartesian direction
      else
! plot marked particles with color
         if (ntsc > 0) then
            ltag = idimp
            call PMBGRASP13(ppart,kpic,label,itime,isc,omx,omy,omz,nx,  &
     &iyp,ixp,idimp,nppmx,mx1,ltag,irc)
! plot unmarked particles with no colors
         else
            call PBGRASP13(ppart,kpic,label,itime,isc,omx,omy,omz,nx,iyp&
     &,ixp,idimp,nppmx,mx1,irc)
         endif
      endif
      end subroutine
!
      end module
