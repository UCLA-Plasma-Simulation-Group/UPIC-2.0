!-----------------------------------------------------------------------
!
      module graf1
!
! Fortran90 interface to 2d PIC Fortran77 library libgks1.f
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
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 23, 2020
!
!     use libgraf1_h
      implicit none
!
! npl = (0,1) = display is (off,on)
      integer, save :: npl = 1
!
      contains
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
   93 format(' VTX =',f9.5,' VTY =',f9.5)
   94 format(' VTX =',f9.5,' VTY =',f9.5,' VTZ =',f9.5)
! chrs = short array of characters to label individual line samples
      data chrs /'    VX    ','    VY    ','    VZ    '/
      if (npl==0) return
      idimv = size(fv,2)
      if ((idimv /= 2) .and. (idimv /= 3)) return
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
         if (idimv==2) then
            write (chr,93) fvm(1,2), fvm(2,2)
            vmax = max(fv(nmv21+1,1),fv(nmv21+1,2))
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
! it = (1,2,3,4,5) plot (x,y,vx,vy,vz)
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
      character(len=4), dimension(5) :: cs 
      character(len=10), dimension(size(pt,3)) :: chrs
! chrs = short array of characters to label individual line samples
      data cs /' X  ',' Y  ',' VX ',' VY ',' VZ '/
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
      end module
