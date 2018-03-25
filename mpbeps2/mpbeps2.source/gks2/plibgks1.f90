! general 1d parallel gks graphics library
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: february 10, 2018
      module plibgks1
      use mpplib2
      implicit none
!
! data in common blocks
      integer idwk, ncols, iplot, nplot, iclr, iupd, idstr, idloc, nclsp
      integer ifrg, isx, isy, kprime
      real rx, ry
! idwk = workstation identifier
! ncols = number of foreground colors available for line plotting
! rx, ry = ndc coordinates of upper-right corner of workstation window
! iplot = plot location on page, 0 <= iplot < nplot
! nplot = number of plots per page
! iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
! iupd = (-1,0,1) = (no,default,yes) end plot
! (default=when iplot=nplot-1)
! idstr = string device number, 0 if no string device available
! idloc = locator device number, 0 if no locator device available
! nclsp = number of foreground colors supported on device
! ifrg = index of foreground color
! isx, isy = display width, height, in raster units
! kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc&
     &,nclsp,ifrg,isx,isy,kprime(8)
      save
!
      private :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
!
      contains
!
!-----------------------------------------------------------------------
      subroutine PDISPR(f,g,nvp,label,xmin,xmax,isc,ist,mks,nx,nxv,ngs, &
     &chr,chrs,irc)
! this subroutine plots ngs subarrays of the array f, on a common graph,
! each plot with nx points, versus a linear function in x,
! where xmin < x < xmax, for distributed data.
! depending on the number of colors in the display device, each subarray
! is plotted with a different color, given in order by:
! blue, red, yellow, cyan, magenta, green, and foreground.
! after all the colors are cycled through, then different line styles
! are cycled through if mks=0, in order: solid, dash, dot, dash-dot,
! or different marker types if mks>0: dot, plus, star, circle, cross.
! multiple plots per page can be displayed by dividing the screen into
! n x n subregions, where n*n is the next largest integer >= nplot
! the location (ix,iy) of a plot in the subregions is determined by
! the parameter iplot = ix + iy*n
! f = distributed field array to be plotted
! g = scratch array for receiving messages
! nvp = number of real or virtual processors requested
! label = long character string label for plot
! xmin/xmax = range of x values in plot
! isc = power of 2 scale of y coordinate for plot
! ist = flag for choosing positive and/or negative values
! the plots have a common scale in y given by ymax and ymin.
! if ist = 0, then ymax = 2**isc and ymin = -2**isc.
! if ist = 1, then ymax = 2**isc and ymin = 0.
! if ist = -1, then ymax = 0 and ymin = -2**isc.
! if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! which will contain the plots, determined by the absolute value of f
! mks = flag to determine whether lines or markers are used,
! mks=0 means cycle through lines styles, mks > 0 means cycle through
! marker styles, using the value of mks as the initial marker,
! mks < 0 means draw the first subarray with a line, then subsequent
! subarrays with markers, using the value of abs(mks) as the initial
! marker.
! nx = number of points plotted in each subarray
! nxv = first dimension of field array f, must be >= nx/nvp
! ngs = second dimension of array f, number of subarrays to be plotted
! chr = additional long character string comment for plot
! chrs = array of ngs short character labels used by subroutine tickd
! to label individual line or marker samples
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nvp, isc, ist, mks, nx, nxv, ngs
      integer, intent(inout) :: irc
      real, intent(in) :: xmin, xmax
      character(len=*), intent(in) :: label, chr
      real, dimension(nxv,ngs), intent(in) :: f
      real, dimension(nxv,ngs), intent(inout) :: g
      character(len=*), dimension(ngs), intent(in) :: chrs
! idwk = workstation identifier
! rx, ry = ndc coordinates of upper-right corner of workstation window
! iplot = plot location on page, 0 <= iplot < nplot
! nplot = number of plots per page
! iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
! iupd = (-1,0,1) = (no,default,yes) end plot
! (default=when iplot=nplot-1)
! nproc = number of real or virtual processors obtained
! mreal = default datatype for reals
! lworld = MPI_COMM_WORLD communicator
! local data
      integer msg, istatus
      integer nlts, nmks, ntx, nty, nxp, nxp1, i, j, k, is, it, ir
      integer npl1, npl, iy, ix, nrt, idproc, joff, id, lvp, nl, ierr
      real gmax, hmax
      real dv, smin, tmin, smax, tmax, csize, algdvi, fmax, fmin
      real rmax, rmin, ymin, ymax, dx, apl, aplx, aply, orx, ory
      real smn, smx, tmn, tmx, chh, xmn
      dimension msg(2), istatus(lstat)
      dimension gmax(2), hmax(2)
! dv = scale will be set in powers of this parameter
      data dv /2.0/
! smin/smax = range of x values of plotting window
! tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
! nlts = number of line types available
! nmks = number of markers available
! ntx/nty = number of ticks in grid in x/y direction
      data nlts,nmks,ntx,nty /4,5,11,9/
! csize = vertical size of characters
      data csize /0.034/
! set return code to normal
      irc = 0
! exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
      nxp = (nx - 1)/nvp + 1
      nxp1 = nxp + 1
! this segment is used for shared memory computers
!     idproc = 0
!     lvp = nvp
! this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
! find scales for plot
      is = isc
! nodes with data find range
      if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
         fmax = f(1,1)
         fmin = fmax
         do 20 k = 1, ngs
         do 10 j = 1, nxp
         fmax = amax1(fmax,f(j,k))
         fmin = amin1(fmin,f(j,k))
   10    continue
   20    continue
! this line is used for mpi computers
         hmax(1) = -fmin
         hmax(2) = fmax
         call PMAX(hmax,gmax,2,1)
         fmin = -hmax(1)
         fmax = hmax(2)
         if (fmax.eq.0.) fmax = 1.0e-35
         rmax = fmax - fmin
         if (rmax.eq.0.) rmax = 1.0e-35
         rmin = fmin
         ymax = abs(fmax)
         is = alog(ymax)*algdvi
         if (ymax.ge.1.) is = is + 1
         if (ymax.le.dv**(is-1)) is = is - 1
         ymin = abs(fmin)  
         if (ymin.gt.0.) then
            it = alog(ymin)*algdvi
            if (ymin.ge.1.) it = it + 1
            if (ymin.le.dv**(it-1)) it = it - 1
         endif
         if (fmax.gt.0.) then
            if (fmin.gt.0.) then
               fmin = dv**(it - 1)
            else if (fmin.lt.0.) then
               fmin = -dv**it
            endif
            fmax = dv**is
         else
            fmax = -dv**(is - 1)
            fmin = -dv**it
         endif
         if (ist.eq.0) then
            if (ymin.gt.ymax) then
               fmax = dv**it
            else
               fmin = -dv**is
            endif
         else if (ist.eq.2) then
            ir = alog(rmax)*algdvi
            if (rmax.ge.1.) ir = ir + 1
            if (rmax.le.dv**(ir-1)) ir = ir - 1
            fmin = rmin
            fmax = rmin + dv**ir
         endif
      else
         fmax = dv**is
         fmin = -fmax
      endif
! clip range if necessary
      if (ist.eq.1) then
         ymin = 0.
      else if (ist.eq.(-1)) then
         ymax = 0.  
      endif
! parameters for plots
      dx = xmax - xmin
      if (nx.gt.1) dx = dx/real(nx)
! find location for plot
      npl1 = sqrt(real(nplot-1)) + 0.0001
      npl = npl1 + 1
      apl = 1./real(npl)
      iy = iplot/npl
      ix = iplot - iy*npl
      aplx = apl*rx
      aply = apl*ry
      orx = aplx*real(ix)
      ory = aply*real(npl1 - iy)
      smn = orx + aplx*smin
      smx = orx + aplx*smax
      tmn = ory + aply*tmin
      tmx = ory + aply*tmax
      chh = aply*csize
! initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
! clear workstation, always
         if (idproc.eq.0) call gclrwk(idwk,1)
      endif
      if (idproc.eq.0) then
! draw grid and labels, call identity transformation
         call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,&
     &label,chr,chh)
! display sample lines or markers
         call dsmpln(orx,smn,tmn,tmx,ngs,nlts,nmks,mks,chrs,chh)
! define transformation number 2
         nrt = 1
! set window
         call gswn(nrt,xmin,xmax,ymin,ymax)
! set viewport
         call gsvp(nrt,smn,smx,tmn,tmx)
! select normalization transformation
         call gselnt(nrt)
! set linewidth scale factor, 1.0 = nominal
         call gslwsc(1.0)
! use markers
         if (mks.ne.0) then
! set marker size scale factor, 1.0 = nominal
            call gsmksc(1.0)
         endif
! set clipping indicator, 1 = on
         call gsclip(1)
      endif
! plot curves
! this segment is used for mpi computers
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      nl = nxv*(ngs - 1) + nxp1
      if (idproc.eq.0) then
! processor 0 plots his (or her) own data
         call mplotit(f,xmin,dx,nxv,ngs,nxp1,mks,nlts,nmks)
! then collects data from remaining nodes
         do 30 i = 2, nvp
         id = i - joff
         call MPI_RECV(g,nl,mreal,id,102,lworld,istatus,ierr)
! plot curves
         xmn = xmin + dx*real(nxp*(i - 1))
         call mplotit(g,xmn,dx,nxv,ngs,nxp1,mks,nlts,nmks)
   30    continue
! other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         call MPI_SEND(f,nl,mreal,0,102,lworld,ierr)
      endif
! update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
! update workstation, perform
         if (idproc.eq.0) call guwk(idwk,1)
! read code from input device, if present
         if (idproc.eq.0) call readrc(irc)
! this segment is used for mpi computers
         msg(1) = irc
         msg(2) = nplot
         call PBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
! reset defaults
      iclr = 0
      iupd = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PDISPS(f,g,nvp,label,xmin,xmax,isc,ist,nx,nxv,chr,irc)
! this subroutine plots an array f versus a linear function in x,
! where xmin < x < xmax.  It is plotted in solid line style, in blue
! if color is available, for distributed data.
! multiple plots per page can be displayed by dividing the screen into
! n x n subregions, where n*n is the next largest integer >= nplot
! the location (ix,iy) of a plot in the subregions is determined by
! the parameter iplot = ix + iy*n
! f = distributed field array to be plotted
! g = scratch array for receiving messages
! nvp = number of real or virtual processors requested
! label = long character string label for plot
! xmin/xmax = range of x values in plot
! isc = power of 2 scale of y coordinate for plot
! ist = flag for choosing positive and/or negative values
! the plot has a scale in y given by ymax and ymin.
! if ist = 0, then ymax = 2**isc and ymin = -2**isc.
! if ist = 1, then ymax = 2**isc and ymin = 0.
! if ist = -1, then ymax = 0 and ymin = -2**isc.
! if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! which will contain the plot, determined by the absolute value of f
! nx = number of points plotted
! nxv = first dimension of field array f, must be >= nx/nvp
! chr = additional long character string comment for plot
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nvp, isc, ist, nx, nxv
      integer, intent(inout) :: irc
      real, intent(in) :: xmin, xmax
      character(len=*), intent(in) :: label, chr
      real, dimension(nxv), intent(in) :: f
      real, dimension(nxv), intent(inout) :: g
      character*1 chs(1)
      data chs /' '/
      call PDISPR(f,g,nvp,label,xmin,xmax,isc,ist,0,nx,nxv,1,chr,chs,irc&
     &)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mplotit(f,xmin,dx,nxv,ngs,nxp1,mks,nlts,nmks)
! this subroutine plots curves, first cycle through colors, then line or
! marker types
! f = distributed field array to be plotted
! xmin/dx = range/increment of x values in plot
! ngs = second dimension of array f, number of subarrays to be plotted
! nxp1 = number of points per processor plotted in each subarray
! mks = flag to determine whether lines or markers are used,
! mks=0 means cycle through lines styles, mks > 0 means cycle through
! marker styles, using the value of mks as the initial marker,
! mks < 0 means draw the first subarray with a line, then subsequent
! subarrays with markers, using the value of abs(mks) as the initial
! marker.
! nlts = number of line types available
! nmks = number of markers available
      implicit none
      integer, intent(in) :: nxv, ngs, nxp1, mks, nlts, nmks
      real, intent(in) :: xmin, dx
      real, dimension(nxv,ngs), intent(in) :: f
! ncols = number of foreground colors available for line plotting
! kprime = table of color indices for prime colors
! local data
      integer i, j, k, nxs, mkr, icol, nxb, npts, nptl, ltype, mtype
      integer npt, js
! nxbs = length of scratch variable for plotting
      integer nxbs
      real x, y
      parameter(nxbs=65)
      dimension x(nxbs), y(nxbs)
      nxs = nxbs - 1
      mkr = abs(mks)
! plot curves, first cycle through colors, then line or marker types
      do 50 k = 1, ngs
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
! use line types
      if ((mks.eq.0).or.((mks.lt.0).and.(k.eq.1))) then
! blocking parameters for plots
         nxb = (nxp1 - 2)/nxs + 1
         npts = nxbs
! length of last block
         nptl = nxp1 - nxs*(nxb - 1)
! set polyline color index
! 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gsplci(icol)
         ltype = 1 + (k - 1)/ncols - nlts*((k - 1)/(nlts*ncols))
! set linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
         call gsln(ltype)
! plot curve
         npt = npts
! loop over number of blocks
         do 20 j = 1, nxb
         js = nxs*(j - 1)
         if (j.eq.nxb) npt = nptl
! calculate x,y axes for block
         do 10 i = 1, npt
         x(i) = xmin + dx*real(i + js - 1)
         y(i) = f(i+js,k)
   10    continue
! draw polyline
         call gpl(npt,x,y)
   20    continue
! use markers
      else
! blocking parameters for plots
         nxb = (nxp1 - 1)/nxs + 1
         npts = nxs
! length of last block
         nptl = nxp1 - nxs*(nxb - 1)
! set polymarker color index
! 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
         mtype = mkr + (k - 1)/ncols - nmks*((mkr - 1 + (k - 1)/ncols)/ &
     &nmks)
! set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mtype)
! plot polymarkers
         npt = npts
! loop over number of blocks
         do 40 j = 1, nxb
         js = nxs*(j - 1)
         if (j.eq.nxb) npt = nptl
! calculate x,y axes for block
         do 30 i = 1, npt
         x(i) = xmin + dx*real(i + js - 1)
         y(i) = f(i+js,k)
   30    continue
! dots
         if (mtype.eq.1) then
! treat dots by drawing a line to itself
            call spdots(x,y,npt,icol,nxbs)
         else
! draw polymarker
            call gpm(npt,x,y)
         endif
   40    continue
      endif
   50 continue
      end subroutine
!
      end module