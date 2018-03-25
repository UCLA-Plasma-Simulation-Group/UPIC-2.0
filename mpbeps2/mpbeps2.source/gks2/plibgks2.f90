!-----------------------------------------------------------------------
! general 2d parallel gks graphics library
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: february 2, 2018
      module plibgks2
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
      subroutine IPLTCOMM(nplt)
! initialize required elements of common block plotcm,
! to allow only one MPI node to display graphics.
      implicit none
      integer, intent(in) :: nplt
      nplot = nplt
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PGRCLOSE
! this subroutine deactivates workstation and closes gks
      implicit none
! idwk = workstation identifier
! lworld = MPI_COMM_WORLD communicator
! local data
      integer idproc, ierr, irc
      dimension irc(1)
! pause if plots are still pending
      if (((iplot.ne.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
! update workstation, perform
         if (idproc.eq.0) call guwk(idwk,1)
! this segment is used for mpi computers
         call MPI_COMM_RANK(lworld,idproc,ierr)
! read code from input device, if present
         if (idproc.eq.0) call readrc(irc(1))
! this segment is used for mpi computers
         call PPBICAST(irc,1)
      endif
! deactivate workstation
      call gdawk(idwk)
! close workstation
      call gclwk(idwk)
! close gks
      call gclks
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PDSYNC(irc)
! this subroutine synchronizes irc and iplot element in plotcm
      implicit none
      integer, intent(inout) :: irc
! local data
      integer, dimension(2) :: msg
      msg(1) = irc
      msg(2) = iplot
      call PPBICAST(msg,2)
      irc = msg(1)
      iplot = msg(2)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PCARPET(f,g,nyp,nvp,label,isc,ist,nx,ny,nxv,nypmx,chr, &
     &ntc,irc)
! this subroutine displays an array f as a color raster image, for
! distributed data
! a 256 color palette must have been defined prior to this call. 
! multiple plots per page can be displayed by dividing the screen into
! n x n subregions, where n*n is the next largest integer >= nplot
! the location (ix,iy) of a plot in the subregions is determined by
! the parameter iplot = ix + iy*n
! f = distributed field array to be plotted
! g = scratch array for receiving messages
! nyp = number of primary gridpoints in field partition
! nvp = number of real or virtual processors requested
! label = long character string label for plot
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! the range of values of f are given by fmax and fmin.
! if ist = 0, then fmax = 2**isc and fmin = -2**isc.
! if ist = 1, then fmax = 2**isc and fmin = 0.
! if ist = -1, then fmax = 0 and fmin = -2**isc.
! if ist = 2, then fmax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! which will contain the plots, determined by the absolute value of f
! nx/ny = length of field f in x/y direction
! nxv = first dimension of field array f, must be >= nx
! nypmx = maximum size of particle partition, including guard cells
! chr = additional long character string comment for plot
! ntc = number of valid colors, should be power of 2, <= 256
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nyp, nvp, isc, ist, nx, ny, nxv, nypmx, ntc
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label, chr
      real, dimension(nxv,nypmx), intent(in) :: f
      real, dimension(nxv,nypmx), intent(inout) ::  g
! idwk = workstation identifier
! rx, ry = ndc coordinates of upper-right corner of workstation window
! iplot = plot location on page, 0 <= iplot < nplot
! nplot = number of plots per page
! iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
! iupd = (-1,0,1) = (no,default,yes) end plot
! (default=when iplot=nplot-1)
! isx, isy = display width, height, in raster units
      integer npald, lxm, lym, lupt, ipal, img8
! npald = number of palette entries
      parameter(npald=256)
! lupt = (0,1) = (no,yes) pixel lookup table needed
! ipal = integer pixel lookup table
      dimension ipal(npald)
! lxm, lym = maximum number of pixels in x, y
      parameter(lxm=720,lym=540)
! img8 = integer image array
      dimension img8(lxm*lym)
      common /movicm/ lupt,ipal,img8
! nproc = number of real or virtual processors obtained
! mreal = default datatype for reals
! lworld = MPI_COMM_WORLD communicator
! local data
      integer msg, istatus
      integer ntx, nty, istyle, nyp1, is, j, k, it, ir, npl1, npl
      integer iy, ix, id, lxs, lys, lxp, lyp, ncv, ic, i1
      integer joff, idproc, i, lvp, nyps, ierr
      real gmax, hmax
      real dv, smin, tmin, smax, tmax, csize, algdvi, fmax, fmin
      real rmax, rmin, xmin, xmax, ymin, ymax, apl, sx, sy, orx, ory
      real smn, smx, tmn, tmx, chh, xmn, ymn, ac
      dimension msg(2), istatus(lstat)
      dimension gmax(2), hmax(2)
! dv = scale will be set in powers of this parameter
      data dv /2.0/
! smin/smax = range of x values of plotting window
! tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
! ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /11,11/
! csize = vertical size of characters
      data csize /0.034/
! istyle = (0,1) = color map (fills area,preserves aspect ratio)
      data istyle /1/
! set return code to normal
      irc = 0
! exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
      nyp1 = nyp + 1
! this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
! find scales for plot
      is = isc
! nodes with data find range
      if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
         fmax = f(1,1)
         fmin = fmax
         do 20 k = 1, nyp1
         do 10 j = 1, nx
         fmax = amax1(fmax,f(j,k))
         fmin = amin1(fmin,f(j,k))
   10    continue
   20    continue
! this line is used for mpi computers
         hmax(1) = -fmin
         hmax(2) = fmax
         call PPMAX(hmax,gmax,2)
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
         fmin = 0.
      else if (ist.eq.(-1)) then
         fmax = 0.  
      endif
! parameters for plots
      xmin = 0.
      xmax = real(nx - 1)
      ymin = 0.
      ymax = real(ny)
! find location for plot
      npl1 = sqrt(real(nplot-1)) + 0.0001
      npl = npl1 + 1
      apl = 1./real(npl)
      iy = iplot/npl
      ix = iplot - iy*npl
      sx = apl*rx
      sy = apl*ry
      orx = sx*real(ix)
      ory = sy*real(npl1 - iy)
      smn = orx + sx*smin
      smx = orx + sx*smax
      tmn = ory + sy*tmin
      tmx = ory + sy*tmax
      chh = sy*csize
! fill area
      xmn = smn
      ymn = tmn
! preserve aspect ratio
      if (istyle.eq.1) then
         if (nx.gt.ny) ymn = tmx - (tmx - tmn)*real(ny)/real(nx)
         if (ny.gt.nx) xmn = smx - (smx - smn)*real(nx)/real(ny)
      endif
! set size of raster image
      lxs = lxm
      lys = lym
      lxp = min(2,lxs)
      lyp = lys/nvp
! initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
! clear workstation, always
         if (idproc.eq.0) call gclrwk(idwk,1)
      endif
! create sample palette
      ac = real(ntc - 1)/real(lys - 1)
! rescale factor ncv
      ncv = 256/ntc
      do 40 k = 1, lys
      joff = lxp*(k - 1)
      ic = ac*real(lys - k) + 0.999999
! rescale index
      ic = ic*ncv
! lookup table required
      if (lupt.eq.1) ic = ipal(ic+1)
      do 30 j = 1, lxp
      img8(j+joff) = ic
   30 continue
   40 continue
      if (idproc.eq.0) then
! draw grid and labels, call identity transformation
         call tickz(xmin,xmax,ymin,ymax,orx,ory,xmn,smx,ymn,tmx,ntx,nty,&
     &label,chr,chh)
! display sample palette
         call dsmpal(img8,fmin,fmax,orx,smn,tmn,tmx,lxp,lys,chh)
      endif
! map parts of f to color raster image
      ic = lyp
! this segment is used for mpi computers
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      if (idproc.eq.0) then
! processor 0 maps his (or her) own data
         i1 = lxs*(lys - lyp) + 1
         call mraster(f,img8(i1),fmin,fmax,nx,nyp1,nxv,lxs,ic,ntc)
! then collects data from remaining nodes
         do 50 i = 2, nvp
         id = i - joff
         call MPI_RECV(g,nxv*nypmx,mreal,id,101,lworld,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nyps,ierr)
         nyps = nyps/nxv
         i1 = lxs*(lys - lyp*i) + 1
! then map the remote data
         if (i.eq.nvp) ic = lyp
         call mraster(g,img8(i1),fmin,fmax,nx,nyps,nxv,lxs,ic,ntc)
   50    continue
! copy with lookup table
         if (lupt.eq.1) then
            do 70 k = 1, lys
            joff = lxs*(k - 1)
            do 60 j = 1, lxs
            img8(j+joff) = ipal(img8(j+joff)+1)
   60       continue
   70       continue
         endif
! other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff
         call MPI_SEND(f,nxv*nyp1,mreal,0,101,lworld,ierr)
      endif
! cell array
      if (idproc.eq.0) then
! special case for rs/6000 with graPHIGS gks
!        call gca(xmn,tmx,smx,ymn,lys,lxs,1,1,lys,lxs,img8)
         call gca(xmn,tmx,smx,ymn,lxs,lys,1,1,lxs,lys,img8)
! update workstation, perform
         call guwk(idwk,1)
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
         call PPBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
! reset defaults
      iclr = 0
      iupd = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PCONTUR(f,g,lf,nyp,nvp,label,isc,ist,nx,ny,nxv,nypmx,  &
     &chr,nc,irc)
! this subroutine displays an array f as a contour plot.
! a maximum of ncols colors are used, used in order from lowest to
! highest contour: blue, green, cyan, foreground, yellow, magenta, red
! multiple plots per page can be displayed by dividing the screen into
! n x n subregions, where n*n is the next largest integer >= nplot
! the location (ix,iy) of a plot in the subregions is determined by
! the parameter iplot = ix + iy*n
! f = field array to be plotted
! lf = scratch field array
! g = scratch array for receiving messages
! nyp = number of primary gridpoints in field partition
! nvp = number of real or virtual processors requested
! label = long character string label for plot
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! the range of values of f are given by fmax and fmin.
! if ist = 0, then fmax = 2**isc and fmin = -2**isc.
! if ist = 1, then fmax = 2**isc and fmin = 0.
! if ist = -1, then fmax = 0 and fmin = -2**isc.
! if ist = 2, then fmax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! which will contain the plots, determined by the absolute value of f
! nx/ny = length of field f in x/y direction
! nxv = first dimension of field array f, must be >= nx
! nypmx = maximum size of particle partition, including guard cells
! chr = additional long character string comment for plot
! nc = number of contour lines
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nyp, nvp, isc, ist, nx, ny, nxv, nypmx, nc
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label, chr
      real, dimension(nxv,nypmx), intent(in) :: f
      real, dimension(nxv,nypmx), intent(inout) ::  g
      integer, dimension(nxv,nypmx+1), intent(inout) ::  lf
! idwk = workstation identifier
! ncols = number of foreground colors available for line plotting
! rx, ry = ndc coordinates of upper-right corner of workstation window
! iplot = plot location on page, 0 <= iplot < nplot
! nplot = number of plots per page
! iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
! iupd = (-1,0,1) = (no,default,yes) end plot
! (default=when iplot=nplot-1)
! kprime = table of color indices for prime colors
! nproc = number of real or virtual processors obtained
! mreal = default datatype for reals
! lworld = MPI_COMM_WORLD communicator
! local data
      integer icv, icolor, msg, istatus
      integer ntx, nty, istyle, nyp1, is, j, k, it, ir, npl1, npl
      integer iy, ix, id, ntc, joff, idproc, i, lvp, nypp, nyps, ierr
      real gmax, hmax
      real dv, smin, tmin, smax, tmax, csize, algdvi, fmax, fmin
      real rmax, rmin, xmin, xmax, ymin, ymax, apl, sx, sy, orx, ory
      real smn, smx, tmx, tmn, chh, xmn, ymn, dyp, ymnp, tmxp
! icolor = color index, used in order from lowest to highest contour
      dimension icv(8), icolor(8), msg(2), istatus(lstat)
      dimension gmax(2), hmax(2)
! dv = scale will be set in powers of this parameter
      data dv /2.0/
! smin/smax = range of x values of plotting window
! tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
! ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /11,11/
! csize = vertical size of characters
      data csize /0.034/
! icv = location in kprime array of color indices needed for icolor
      data icv /1,3,8,6,2,5,7,4/
! istyle = (0,1) = contour plot (fills area,preserves aspect ratio)
      data istyle /1/
! set return code to normal
      irc = 0
! exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
      nyp1 = nyp + 1
! this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
! find scales for plot
      is = isc
! nodes with data find range
      if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
         fmax = f(1,1)
         fmin = fmax
         do 20 k = 1, nyp1
         do 10 j = 1, nx
         fmax = amax1(fmax,f(j,k))
         fmin = amin1(fmin,f(j,k))
   10    continue
   20    continue
! this line is used for mpi computers
         hmax(1) = -fmin
         hmax(2) = fmax
         call PPMAX(hmax,gmax,2)
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
         fmin = 0.
      else if (ist.eq.(-1)) then
         fmax = 0.  
      endif
! parameters for plots
      xmin = 0.
      xmax = real(nx - 1)
      ymin = 0.
      ymax = real(ny)
! find location for plot
      npl1 = sqrt(real(nplot-1)) + 0.0001
      npl = npl1 + 1
      apl = 1./real(npl)
      iy = iplot/npl
      ix = iplot - iy*npl
      sx = apl*rx
      sy = apl*ry
      orx = sx*real(ix)
      ory = sy*real(npl1 - iy)
      smn = orx + sx*smin
      smx = orx + sx*smax
      tmn = ory + sy*tmin
      tmx = ory + sy*tmax
      chh = sy*csize
! fill area
      xmn = smn
      ymn = tmn
! preserve aspect ratio
      if (istyle.eq.1) then
         if (nx.gt.ny) ymn = tmx - (tmx - tmn)*real(ny)/real(nx)
         if (ny.gt.nx) xmn = smx - (smx - smn)*real(nx)/real(ny)
      endif
      dyp = (tmx - ymn)/real(ny+1)
! calculate color indices
      ntc = ncols + 1
      do 30 i = 1, ntc
      icolor(i) =  kprime(icv(i))
   30 continue
! initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
! clear workstation, always
         if (idproc.eq.0) call gclrwk(idwk,1)
      endif
      if (idproc.eq.0) then
! draw grid and labels, call identity transformation
         call tickz(xmin,xmax,ymin,ymax,orx,ory,xmn,smx,ymn,tmx,ntx,nty,&
     &label,chr,chh)
! display sample contours
         call dsmpcn(icolor,fmin,fmax,orx,smn,tmn,tmx,chh,ntc,nc)
      endif
! display parts of f as contour plot
! this segment is used for mpi computers
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      if (idproc.eq.0) then
! processor 0 maps his (or her) own data
         ymnp = ymn
         tmxp = ymn + dyp*real(nyp+1)
         nypp = nyp
! draw contour map
         call dcontr(f,lf,icolor,fmin,fmax,xmn,smx,ymnp,tmxp,nx,nyp1,nxv&
     &,ntc,nc)
! then collects data from remaining nodes
         do 40 i = 2, nvp
         id = i - joff
         call MPI_RECV(g,nxv*nypmx,mreal,id,101,lworld,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nyps,ierr)
         nyps = nyps/nxv
         ymnp = ymnp + dyp*real(nypp)
         tmxp = tmxp + dyp*real(nyps)
         nypp = nyps - 1
! then draw contour map of remote data
         call dcontr(g,lf,icolor,fmin,fmax,xmn,smx,ymnp,tmxp,nx,nyps,nxv&
     &,ntc,nc)
   40    continue
! other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff - 1
         call MPI_SEND(f,nxv*nyp1,mreal,0,101,lworld,ierr)
      endif
!
! update workstation, perform
      if (idproc.eq.0) call guwk(idwk,1)
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
         call PPBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
! reset defaults
      iclr = 0
      iupd = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PGRASP23(part,f,npp,label,itime,isc,nx,ny,iyp,ixp,idimp&
     &,npmax,irc)
! for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! part(5,n) = velocity vz of particle n in partition
! f = scratch array for receiving messages
! npp = number of particles in partition
! label = species label
! itime = current time step
! isc = power of 2 scale of range of values of velocity
! nx/ny = system length in x/y direction
! iyp/ixp = phase space coordinates to be displayed
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: npp
      integer, intent(in) :: itime, isc, nx, ny, idimp, npmax, iyp, ixp
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(idimp,npmax), intent(in) :: part
      real, dimension(2*npmax), intent(inout) :: f
! idwk = workstation identifier
! ncols = number of foreground colors available for line plotting
! rx, ry = ndc coordinates of upper-right corner of workstation window
! iplot = plot location on page, 0 <= iplot < nplot
! nplot = number of plots per page
! iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
! iupd = (-1,0,1) = (no,default,yes) end plot
! (default=when iplot=nplot-1)
! kprime = table of color indices for prime colors
! nproc = number of real or virtual processors obtained
! mint = default datatype for integers
! mreal = default datatype for reals
! lworld = MPI_COMM_WORLD communicator
! local data
      integer msg, istatus
      integer ntx, nty, is, npl1, npl, iy, ix, mks, ngs, nrt
      integer i, j, nd, icol, joff, id, idproc, nvp, lvp, ierr
      real gmax, hmax
      real dv, smin, tmin, smax, tmax, zero, csize, algdvi, apl
      real fmax, ymin, ymax, xmin, xmax, aplx, aply, orx, ory, smn, smx
      real tmn, tmx, chh, dd
      dimension msg(2), istatus(lstat)
      dimension gmax(1), hmax(1)
      character(len=26) lbl
      character(len=2) lblsp(5)
      character(len=10) chrs(2)
      save lblsp
   91 format (1x,a2,' VERSUS ',a2,', T = ',i7)
      data lblsp /' X',' Y','VX','VY','VZ'/
! dv = scale will be set in powers of this parameter
      data dv /2.0/
! smin/smax = range of x values of plotting window
! tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
! ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
! csize = vertical size of characters
      data zero,csize /0.,0.034/
! sample labels
      data chrs /'BACKGROUND','   BEAM   '/
! set return code to normal
      irc = 0
! exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
! this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
! find y scale for plot
      if (iyp.le.1) then
         ymin = zero
         ymax = real(nx)
      elseif (iyp.eq.2) then
         ymin = zero
         ymax = real(ny)
      elseif (iyp.gt.2) then
         is = isc
         if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
            fmax = abs(part(iyp,1))
            do 10 j = 1, npp
            fmax = amax1(fmax,abs(part(iyp,j)))
   10       continue
! this line is used for mpi computers
            gmax(1) = fmax
            call PPMAX(gmax,hmax,1)
            ymax = gmax(1)
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
         endif
         ymax = dv**is
         ymin = -ymax
      endif
! find x scale for plot
      if (ixp.le.1) then
         xmin = zero
         xmax = real(nx)
      elseif (ixp.eq.2) then
         xmin = zero
         xmax = real(ny)
      elseif (ixp.gt.2) then
         is = isc
         if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
            fmax = abs(part(ixp,1))
            do 20 j = 1, npp
            fmax = amax1(fmax,abs(part(ixp,j)))
   20       continue
! this segment is used for mpi computers
            gmax(1) = fmax
            call PPMAX(gmax,hmax,1)
            xmax = gmax(1)
            if (xmax.eq.0.) xmax = 1.0e-35
            is = alog(xmax)*algdvi
            if (xmax.ge.1.) is = is + 1
            if (xmax.le.dv**(is-1)) is = is - 1
         endif
         xmax = dv**is
         xmin = -xmax
      endif
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
! write labels
      write (lbl,91) lblsp(iyp), lblsp(ixp), itime
      if (idproc.eq.0) then
! select point as marker
         mks = 1
! draw grid and labels, call identity transformation
         call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,&
     &label,lbl,chh)
! display sample markers with dots of minimum visible size
         dd = (smax - smin)*real(npl)*2.0e-3
         ngs = 1
!        call dsmpln(orx,smn,tmn,tmx,ngs,1,1,mks,chrs,chh)
         call dsdmpln(orx,smn,tmn,tmx,dd,ngs,1,1,mks,chrs,chh)
! define transformation number 2
         nrt = 1
! set window
         call gswn(nrt,xmin,xmax,ymin,ymax)
! set viewport
         call gsvp(nrt,smn,smx,tmn,tmx)
! select normalization transformation
         call gselnt(nrt)
! set marker size scale factor, 1.0 = nominal
         call gsmksc(1.0)
! set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mks)
! set clipping indicator, 1 = on
         call gsclip(1)
      endif
! plot particles
! this segment is used for mpi computers
      nvp = nproc
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      if (idproc.eq.0) then
         nd = npp
         do 30 j = 1, nd
         f(j) = part(ixp,j)
         f(nd+j) = part(iyp,j)
   30    continue
         icol = 2
         icol = kprime(icol+1)
! set polymarker color index
         call gspmci(icol)
         do 40 i = 1, nvp
! processor 0 plots his (or her) own data
         if (i.gt.1) then
! collects data from remaining nodes
            id = i - joff
            call MPI_RECV(f,2*npmax,mreal,id,103,lworld,istatus,ierr)
! determine how many particles to plot
            call MPI_GET_COUNT(istatus,mreal,nd,ierr)
            nd = nd/2
         endif
         if (nd.eq.0) go to 40
! treat dots by drawing a line to itself
!        call spdots(f(1),f(nd+1),nd,icol,nd)
! treat dots by drawing a line to itself with non-zero width in x
         dd = (xmax - xmin)*real(npl)*2.0e-3
         call spddots(f(1),f(nd+1),dd,nd,icol,nd)
! update workstation, perform
         call guwk(idwk,1)
   40    continue
! other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff - 1
         nd = npp
         do 50 j = 1, nd
         f(j) = part(ixp,j)
         f(nd+j) = part(iyp,j)
   50    continue
         call MPI_SEND(f,2*nd,mreal,0,103,lworld,ierr)
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
         call PPBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
! reset defaults
      iclr = 0
      iupd = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PBGRASP23(part,f,npp,label,itime,isc,omx,omy,omz,nx,ny,&
     &iyp,ixp,idimp,npmax,irc)
! for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
! for magnetized plasma, rotating cartesian co-ordinates so that B
! points in the z direction.
! if iyp=2, plot vperp1, if iyp=3, plot vperp2, if iyp=4, plot vparallel
! if ixp=2, plot vperp1, if ixp=3, plot vperp2, if ixp=4, plot vparallel
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! part(5,n) = velocity vz of particle n in partition
! f = scratch array for receiving messages
! npp = number of particles in partition
! label = species label
! itime = current time step
! isc = power of 2 scale of range of values of velocity
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! nx/ny = system length in x/y direction
! iyp/ixp = phase space coordinates to be displayed
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: npp
      integer, intent(in) :: itime, isc, nx, ny, idimp, npmax, iyp, ixp
      integer, intent(inout) :: irc
      real, intent(in) :: omx, omy, omz
      character(len=*), intent(in) :: label
      real, dimension(idimp,npmax), intent(in) :: part
      real, dimension(2*npmax), intent(inout) :: f
! idwk = workstation identifier
! ncols = number of foreground colors available for line plotting
! rx, ry = ndc coordinates of upper-right corner of workstation window
! iplot = plot location on page, 0 <= iplot < nplot
! nplot = number of plots per page
! iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
! iupd = (-1,0,1) = (no,default,yes) end plot
! (default=when iplot=nplot-1)
! kprime = table of color indices for prime colors
! nproc = number of real or virtual processors obtained
! mint = default datatype for integers
! mreal = default datatype for reals
! lworld = MPI_COMM_WORLD communicator
! local data
      integer msg, istatus
      integer ntx, nty, is, it, npl1, npl, iy, ix, mks, ngs, nrt
      integer i, j, ndir, nd, icol, joff, id, idproc, nvp, lvp, ierr
      real gmax, hmax, vt
      real at1, at2, ox, oy, oz, px, py, pz, qx, qy, qz, vx, vy, vz
      real dv, smin, tmin, smax, tmax, zero, csize, algdvi, apl
      real ymin, ymax, xmin, xmax, aplx, aply, orx, ory, smn, smx
      real tmn, tmx, chh, dd
      dimension msg(2), istatus(lstat)
      dimension gmax(2), hmax(2), vt(3)
      character(len=30) lbl
      character(len=4) lblspb(5), lblspn(5), lblsp(5)
      character(len=10) chrs(2)
      save lblsp
   91 format (1x,a4,' VERSUS ',a4,', T = ',i7)
      data lblspb /' X  ',' Y  ','VPR1','VPR2','VPL '/
      data lblspn /' X  ',' Y  ',' VX ',' VY ',' VZ '/
! dv = scale will be set in powers of this parameter
      data dv /2.0/
! smin/smax = range of x values of plotting window
! tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
! ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
! csize = vertical size of characters
      data zero,csize /0.,0.034/
! sample labels
      data chrs /'BACKGROUND','   BEAM   '/
! set return code to normal
      irc = 0
! exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
! find rotation to convert to cylindrical co-ordinates
      at1 = sqrt(omx*omx + omy*omy + omz*omz)
! no rotation if zero B field
      if (at1.eq.0.0) then
         ox = 0.0
         oy = 0.0
         oz = 1.0
         lblsp(1) = lblspn(1)
         lblsp(2) = lblspn(2)
         lblsp(3) = lblspn(3)
         lblsp(4) = lblspn(4)
         lblsp(5) = lblspn(5)
! create rotation vectors
      else
! first create unit vector in B direction
         at1 = 1.0/at1
         ox = omx*at1
         oy = omy*at1
         oz = omz*at1
         lblsp(1) = lblspb(1)
         lblsp(2) = lblspb(2)
         lblsp(3) = lblspb(3)
         lblsp(4) = lblspb(4)
         lblsp(5) = lblspb(5)
      endif
! then create unit vector in first perpendicular direction
! find direction with smallest component of B
      ndir = 1
      at1 = abs(omx)
      at2 = abs(omy)
      if (at2.le.at1) then
         ndir = 2
         at1 = at2
      endif
      if (abs(omz).lt.at1) ndir = 3
! take the cross product of that direction with B
! vpr1 = x cross B
      if (ndir.eq.1) then
         at1 = 1.0/sqrt(oy*oy + oz*oz)
         px = 0.0
         py = -oz*at1
         pz = oy*at1
! vpr1 = y cross B
      else if (ndir.eq.2) then
         at1 = 1.0/sqrt(ox*ox + oz*oz)
         px = oz*at1
         py = 0.0
         pz = -ox*at1
! vpr1 = z cross B
      else if (ndir.eq.3) then
         at1 = 1.0/sqrt(ox*ox + oy*oy)
         px = -oy*at1
         py = ox*at1
         pz = 0.0
      endif
! finally create unit vector in second perpendicular direction
! vpr2 = B cross vpr1
      qx = oy*pz - oz*py
      qy = oz*px - ox*pz
      qz = ox*py - oy*px
! this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
! find y scale for plot
      if (iyp.le.1) then
         ymin = zero
         ymax = real(nx)
      elseif (iyp.eq.2) then
         ymin = zero
         ymax = real(ny)
      else
         ymax = 0.0
      endif
! find x scale for plot
      if (ixp.le.1) then
         xmin = zero
         xmax = real(nx)
      elseif (ixp.eq.2) then
         xmin = zero
         xmax = real(ny)
      else
         xmax = 0.0
      endif
! find velocity scales for plot
      if ((iyp.gt.2).or.(ixp.gt.2)) then
         is = isc
         if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
            do 10 j = 1, npp
            vx = part(3,j)
            vy = part(4,j)
            vz = part(5,j)
! vperp1 co-ordinate
            vt(1) = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
            vt(2) = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
            vt(3) = vx*ox + vy*oy + vz*oz
            if (iyp.gt.2) ymax = amax1(ymax,abs(vt(iyp-2)))
            if (ixp.gt.2) xmax = amax1(xmax,abs(vt(ixp-2)))
   10       continue
! this line is used for mpi computers
            gmax(1) = ymax
            gmax(2) = xmax
            call PPMAX(gmax,hmax,2)
            ymax = gmax(1)
            xmax = gmax(2)
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
            if (xmax.eq.0.) xmax = 1.0e-35
            it = alog(xmax)*algdvi
            if (xmax.ge.1.) it = it + 1
            if (xmax.le.dv**(it-1)) it = it - 1
         endif
         if (iyp.gt.2) then
            ymax = dv**is
            ymin = -ymax
         endif
         if (ixp.gt.2) then
            xmax = dv**it
            xmin = -xmax
         endif
      endif
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
! write labels
      write (lbl,91) lblsp(iyp), lblsp(ixp), itime
      if (idproc.eq.0) then
! select point as marker
         mks = 1
! draw grid and labels, call identity transformation
         call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,&
     &label,lbl,chh)
! display sample markers with dots of minimum visible size
         dd = (smax - smin)*real(npl)*2.0e-3
         ngs = 1
!        call dsmpln(orx,smn,tmn,tmx,ngs,1,1,mks,chrs,chh)
         call dsdmpln(orx,smn,tmn,tmx,dd,ngs,1,1,mks,chrs,chh)
! define transformation number 2
         nrt = 1
! set window
         call gswn(nrt,xmin,xmax,ymin,ymax)
! set viewport
         call gsvp(nrt,smn,smx,tmn,tmx)
! select normalization transformation
         call gselnt(nrt)
! set marker size scale factor, 1.0 = nominal
         call gsmksc(1.0)
! set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mks)
! set clipping indicator, 1 = on
         call gsclip(1)
      endif
! plot particles
! this segment is used for mpi computers
      nvp = nproc
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      if (idproc.eq.0) then
         nd = npp
         do 20 j = 1, nd
         if ((iyp.gt.2).or.(ixp.gt.2)) then
            vx = part(3,j)
            vy = part(4,j)
            vz = part(5,j)
! vperp1 co-ordinate
            vt(1) = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
            vt(2) = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
            vt(3) = vx*ox + vy*oy + vz*oz
         endif
         if (ixp.le.2) then
            vx = part(ixp,j)
         else
            vx = vt(ixp-2)
         endif
         if (iyp.le.2) then
            vy = part(iyp,j)
         else
            vy = vt(iyp-2)
         endif
         f(j) = vx
         f(nd+j) = vy
   20    continue
         icol = 2
         icol = kprime(icol+1)
! set polymarker color index
         call gspmci(icol)
         do 30 i = 1, nvp
! processor 0 plots his (or her) own data
         if (i.gt.1) then
! collects data from remaining nodes
            id = i - joff
            call MPI_RECV(f,2*npmax,mreal,id,103,lworld,istatus,ierr)
! determine how many particles to plot
            call MPI_GET_COUNT(istatus,mreal,nd,ierr)
            nd = nd/2
         endif
         if (nd.eq.0) go to 30
! treat dots by drawing a line to itself
!        call spdots(f(1),f(nd+1),nd,icol,nd)
! treat dots by drawing a line to itself with non-zero width in x
         dd = (xmax - xmin)*real(npl)*2.0e-3
         call spddots(f(1),f(nd+1),dd,nd,icol,nd)
! update workstation, perform
         call guwk(idwk,1)
   30    continue
! other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff - 1
         nd = npp
         do 40 j = 1, nd
         if ((iyp.gt.2).or.(ixp.gt.2)) then
            vx = part(3,j)
            vy = part(4,j)
            vz = part(5,j)
! vperp1 co-ordinate
            vt(1) = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
            vt(2) = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
            vt(3) = vx*ox + vy*oy + vz*oz
         endif
         if (ixp.le.2) then
            vx = part(ixp,j)
         else
            vx = vt(ixp-2)
         endif
         if (iyp.le.2) then
            vy = part(iyp,j)
         else
            vy = vt(iyp-2)
         endif
         f(j) = vx
         f(nd+j) = vy
   40    continue
         call MPI_SEND(f,2*nd,mreal,0,103,lworld,ierr)
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
         call PPBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
! reset defaults
      iclr = 0
      iupd = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPGRASP23(ppart,f,kpic,label,itime,isc,nx,ny,iyp,ixp,  &
     &idimp,nppmx,mxyp1,irc)
! for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
! ppart(1,n,m) = position x of particle n in tile m in partition
! ppart(2,n,m) = position y of particle n in tile m in partition
! ppart(3,n,m) = velocity vx of particle n in tile m in partition
! ppart(4,n,m) = velocity vy of particle n in tile m in partition
! ppart(5,n,m) = velocity vz of particle n in tile m in partition
! kpic = number of particles per tile
! f = scratch array for receiving messages
! label = species label
! itime = current time step
! isc = power of 2 scale of range of values of velocity
! nx/ny = system length in x/y direction
! iyp/ixp = phase space coordinates to be displayed
! idimp = size of phase space = 4 or 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! and where mx1 = (system length in x direction - 1)/mx + 1
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, nx, ny, iyp, ixp, idimp, nppmx
      integer, intent(in) :: mxyp1
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
      real, dimension(2*nppmx*mxyp1), intent(inout) :: f
      integer, dimension(mxyp1) :: kpic
! idwk = workstation identifier
! ncols = number of foreground colors available for line plotting
! rx, ry = ndc coordinates of upper-right corner of workstation window
! iplot = plot location on page, 0 <= iplot < nplot
! nplot = number of plots per page
! iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
! iupd = (-1,0,1) = (no,default,yes) end plot
! (default=when iplot=nplot-1)
! kprime = table of color indices for prime colors
! nproc = number of real or virtual processors obtained
! mint = default datatype for integers
! mreal = default datatype for reals
! lworld = MPI_COMM_WORLD communicator
! local data
! nxbs = length of scratch variable for plotting
      integer msg, istatus
      integer ntx, nty, is, npl1, npl, iy, ix, mks, ngs, nrt
      integer i, j, k, nd, icol, joff, id, idproc, nvp, lvp, ierr
      integer npmax, npp, noff
      real gmax, hmax
      real dv, smin, tmin, smax, tmax, zero, csize, algdvi, apl
      real fmax, ymin, ymax, xmin, xmax, aplx, aply, orx, ory, smn, smx
      real tmn, tmx, chh, dd
      real sfmax
      dimension msg(2), istatus(lstat)
      dimension gmax(1), hmax(1)
      character(len=26) lbl
      character(len=2) lblsp(5)
      character(len=10) chrs(2)
      save lblsp
   91 format (1x,a2,' VERSUS ',a2,', T = ',i7)
      data lblsp /' X',' Y','VX','VY','VZ'/
! dv = scale will be set in powers of this parameter
      data dv /2.0/
! smin/smax = range of x values of plotting window
! tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
! ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
! csize = vertical size of characters
      data zero,csize /0.,0.034/
! sample labels
      data chrs /'BACKGROUND','   BEAM   '/
! set return code to normal
      irc = 0
! exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
      npmax = nppmx*mxyp1
! this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
! find y scale for plot
      if (iyp.le.1) then
         ymin = zero
         ymax = real(nx)
      elseif (iyp.eq.2) then
         ymin = zero
         ymax = real(ny)
      elseif (iyp.gt.2) then
         is = isc
         if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
            fmax = abs(ppart(iyp,1,1))
!$OMP PARALLEL DO PRIVATE(j,k,npp,sfmax) REDUCTION(MAX:fmax)
            do 20 k = 1, mxyp1
            npp = kpic(k)
            sfmax = abs(ppart(iyp,1,k))
            do 10 j = 1, npp
            sfmax = amax1(sfmax,abs(ppart(iyp,j,k)))
   10       continue
            fmax = amax1(fmax,sfmax)
   20       continue
!$OMP END PARALLEL DO
! this line is used for mpi computers
            gmax(1) = fmax
            call PPMAX(gmax,hmax,1)
            ymax = gmax(1)
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
         endif
         ymax = dv**is
         ymin = -ymax
      endif
! find x scale for plot
      if (ixp.le.1) then
         xmin = zero
         xmax = real(nx)
      elseif (ixp.eq.2) then
         xmin = zero
         xmax = real(ny)
      elseif (ixp.gt.2) then
         is = isc
         if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
            fmax = abs(ppart(ixp,1,1))
!$OMP PARALLEL DO PRIVATE(j,k,npp,sfmax) REDUCTION(MAX:fmax)
            do 40 k = 1, mxyp1
            npp = kpic(k)
            sfmax = abs(ppart(ixp,1,k))
            do 30 j = 1, npp
            sfmax = amax1(sfmax,abs(ppart(ixp,j,k)))
   30       continue
            fmax = amax1(fmax,sfmax)
   40       continue
!$OMP END PARALLEL DO
! this segment is used for mpi computers
            gmax(1) = fmax
            call PPMAX(gmax,hmax,1)
            xmax = gmax(1)
            if (xmax.eq.0.) xmax = 1.0e-35
            is = alog(xmax)*algdvi
            if (xmax.ge.1.) is = is + 1
            if (xmax.le.dv**(is-1)) is = is - 1
         endif
         xmax = dv**is
         xmin = -xmax
      endif
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
! write labels
      write (lbl,91) lblsp(iyp), lblsp(ixp), itime
      if (idproc.eq.0) then
! select point as marker
         mks = 1
! draw grid and labels, call identity transformation
         call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,&
     &label,lbl,chh)
! display sample markers with dots of minimum visible size
         dd = (smax - smin)*real(npl)*2.0e-3
         ngs = 1
!        call dsmpln(orx,smn,tmn,tmx,ngs,1,1,mks,chrs,chh)
         call dsdmpln(orx,smn,tmn,tmx,dd,ngs,1,1,mks,chrs,chh)
! define transformation number 2
         nrt = 1
! set window
         call gswn(nrt,xmin,xmax,ymin,ymax)
! set viewport
         call gsvp(nrt,smn,smx,tmn,tmx)
! select normalization transformation
         call gselnt(nrt)
! set marker size scale factor, 1.0 = nominal
         call gsmksc(1.0)
! set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mks)
! set clipping indicator, 1 = on
         call gsclip(1)
         icol = 2
         icol = kprime(icol+1)
! set polymarker color index
         call gspmci(icol)
      endif
! plot particles
! this segment is used for mpi computers
      nvp = nproc
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
! find how many particles on this node
      nd = 0
      do 50 k = 1, mxyp1
      nd = nd + kpic(k)
   50 continue
      if (idproc.eq.0) then
! copy to contiguous buffer
         noff = 0
         do 70 k = 1, mxyp1
         npp = kpic(k)
         do 60 j = 1, npp
         f(j+noff) = ppart(ixp,j,k)
         f(nd+j+noff) = ppart(iyp,j,k)
   60    continue
         noff = noff + npp
   70    continue
         do 80 i = 1, nvp
! processor 0 plots his (or her) own data
         if (i.gt.1) then
! collects data from remaining nodes
            id = i - joff
            call MPI_RECV(f,2*npmax,mreal,id,103,lworld,istatus,ierr)
! determine how many particles to plot
            call MPI_GET_COUNT(istatus,mreal,nd,ierr)
            nd = nd/2
         endif
         if (nd.eq.0) go to 80
! treat dots by drawing a line to itself
!        call spdots(f(1),f(nd+1),nd,icol,nd)
! treat dots by drawing a line to itself with non-zero width in x
         dd = (xmax - xmin)*real(npl)*2.0e-3
         call spddots(f(1),f(nd+1),dd,nd,icol,nd)
! update workstation, perform
         call guwk(idwk,1)
   80    continue
! other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff - 1
! copy to contiguous buffer
         noff = 0
         do 100 k = 1, mxyp1
         npp = kpic(k)
         do 90 j = 1, npp
         f(j+noff) = ppart(ixp,j,k)
         f(nd+j+noff) = ppart(iyp,j,k)
   90    continue
         noff = noff + npp
  100    continue
         call MPI_SEND(f,2*nd,mreal,0,103,lworld,ierr)
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
         call PPBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
! reset defaults
      iclr = 0
      iupd = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBGRASP23(ppart,f,kpic,label,itime,isc,omx,omy,omz,nx,&
     &ny,iyp,ixp,idimp,nppmx,mxyp1,irc)
! for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
! for magnetized plasma, rotating cartesian co-ordinates so that B
! points in the z direction.
! if iyp=2, plot vperp1, if iyp=3, plot vperp2, if iyp=4, plot vparallel
! if ixp=2, plot vperp1, if ixp=3, plot vperp2, if ixp=4, plot vparallel
! ppart(1,n,m) = position x of particle n in tile m in partition
! ppart(2,n,m) = position y of particle n in tile m in partition
! ppart(3,n,m) = velocity vx of particle n in tile m in partition
! ppart(4,n,m) = velocity vy of particle n in tile m in partition
! ppart(5,n,m) = velocity vz of particle n in tile m in partition
! kpic = number of particles per tile
! f = scratch array for receiving messages
! label = species label
! itime = current time step
! isc = power of 2 scale of range of values of velocity
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! nx/ny = system length in x/y direction
! iyp/ixp = phase space coordinates to be displayed
! idimp = size of phase space = 4 or 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! and where mx1 = (system length in x direction - 1)/mx + 1
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, nx, ny, iyp, ixp, idimp, nppmx
      integer, intent(in) :: mxyp1
      integer, intent(inout) :: irc
      real, intent(in) :: omx, omy, omz
      character(len=*), intent(in) :: label
      real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
      real, dimension(2*nppmx*mxyp1), intent(inout) :: f
      integer, dimension(mxyp1) :: kpic
! idwk = workstation identifier
! ncols = number of foreground colors available for line plotting
! rx, ry = ndc coordinates of upper-right corner of workstation window
! iplot = plot location on page, 0 <= iplot < nplot
! nplot = number of plots per page
! iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
! iupd = (-1,0,1) = (no,default,yes) end plot
! (default=when iplot=nplot-1)
! kprime = table of color indices for prime colors
! nproc = number of real or virtual processors obtained
! mint = default datatype for integers
! mreal = default datatype for reals
! lworld = MPI_COMM_WORLD communicator
! local data
! nxbs = length of scratch variable for plotting
      integer msg, istatus
      integer ntx, nty, is, it, npl1, npl, iy, ix, mks, ngs, nrt
      integer i, j, k, ndir, nd, icol, joff, id, idproc, nvp, lvp, ierr
      integer npmax, npp, noff
      real gmax, hmax, vt
      real at1, at2, ox, oy, oz, px, py, pz, qx, qy, qz, vx, vy, vz
      real dv, smin, tmin, smax, tmax, zero, csize, algdvi, apl
      real ymin, ymax, xmin, xmax, aplx, aply, orx, ory, smn, smx
      real tmn, tmx, chh, dd
      real sxmax, symax
      dimension msg(2), istatus(lstat)
      dimension gmax(2), hmax(2), vt(3)
      character(len=30) lbl
      character(len=4) lblspb(5), lblspn(5), lblsp(5)
      character(len=10) chrs(2)
      save lblsp
   91 format (1x,a4,' VERSUS ',a4,', T = ',i7)
      data lblspb /' X  ',' Y  ','VPR1','VPR2','VPL '/
      data lblspn /' X  ',' Y  ',' VX ',' VY ',' VZ '/
! dv = scale will be set in powers of this parameter
      data dv /2.0/
! smin/smax = range of x values of plotting window
! tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
! ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
! csize = vertical size of characters
      data zero,csize /0.,0.034/
! sample labels
      data chrs /'BACKGROUND','   BEAM   '/
! set return code to normal
      irc = 0
! exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
      npmax = nppmx*mxyp1
! find rotation to convert to cylindrical co-ordinates
      at1 = sqrt(omx*omx + omy*omy + omz*omz)
! no rotation if zero B field
      if (at1.eq.0.0) then
         ox = 0.0
         oy = 0.0
         oz = 1.0
         lblsp(1) = lblspn(1)
         lblsp(2) = lblspn(2)
         lblsp(3) = lblspn(3)
         lblsp(4) = lblspn(4)
         lblsp(5) = lblspn(5)
! create rotation vectors
      else
! first create unit vector in B direction
         at1 = 1.0/at1
         ox = omx*at1
         oy = omy*at1
         oz = omz*at1
         lblsp(1) = lblspb(1)
         lblsp(2) = lblspb(2)
         lblsp(3) = lblspb(3)
         lblsp(4) = lblspb(4)
         lblsp(5) = lblspb(5)
      endif
! then create unit vector in first perpendicular direction
! find direction with smallest component of B
      ndir = 1
      at1 = abs(omx)
      at2 = abs(omy)
      if (at2.le.at1) then
         ndir = 2
         at1 = at2
      endif
      if (abs(omz).lt.at1) ndir = 3
! take the cross product of that direction with B
! vpr1 = x cross B
      if (ndir.eq.1) then
         at1 = 1.0/sqrt(oy*oy + oz*oz)
         px = 0.0
         py = -oz*at1
         pz = oy*at1
! vpr1 = y cross B
      else if (ndir.eq.2) then
         at1 = 1.0/sqrt(ox*ox + oz*oz)
         px = oz*at1
         py = 0.0
         pz = -ox*at1
! vpr1 = z cross B
      else if (ndir.eq.3) then
         at1 = 1.0/sqrt(ox*ox + oy*oy)
         px = -oy*at1
         py = ox*at1
         pz = 0.0
      endif
! finally create unit vector in second perpendicular direction
! vpr2 = B cross vpr1
      qx = oy*pz - oz*py
      qy = oz*px - ox*pz
      qz = ox*py - oy*px
! this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
! find y scale for plot
      if (iyp.le.1) then
         ymin = zero
         ymax = real(nx)
      elseif (iyp.eq.2) then
         ymin = zero
         ymax = real(ny)
      else
         ymax = 0.0
      endif
! find x scale for plot
      if (ixp.le.1) then
         xmin = zero
         xmax = real(nx)
      elseif (ixp.eq.2) then
         xmin = zero
         xmax = real(ny)
      else
         xmax = 0.0
      endif
! find velocity scales for plot
      if ((iyp.gt.2).or.(ixp.gt.2)) then
         is = isc
         if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
!$OMP PARALLEL DO PRIVATE(j,k,npp,vx,vy,vz,vt,symax,sxmax)              &
!$OMP& REDUCTION(MAX:ymax) REDUCTION(MAX:xmax)
            do 20 k = 1, mxyp1
            npp = kpic(k)
            symax = 0.0
            sxmax = 0.0
            do 10 j = 1, npp
            vx = ppart(3,j,k)
            vy = ppart(4,j,k)
            vz = ppart(5,j,k)
! vperp1 co-ordinate
            vt(1) = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
            vt(2) = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
            vt(3) = vx*ox + vy*oy + vz*oz
            if (iyp.gt.2) symax = amax1(symax,abs(vt(iyp-2)))
            if (ixp.gt.2) sxmax = amax1(sxmax,abs(vt(ixp-2)))
   10       continue
            ymax = amax1(ymax,symax)
            xmax = amax1(xmax,sxmax)
   20       continue
!$OMP END PARALLEL DO
! this line is used for mpi computers
            gmax(1) = ymax
            gmax(2) = xmax
            call PPMAX(gmax,hmax,2)
            ymax = gmax(1)
            xmax = gmax(2)
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
            if (xmax.eq.0.) xmax = 1.0e-35
            it = alog(xmax)*algdvi
            if (xmax.ge.1.) it = it + 1
            if (xmax.le.dv**(it-1)) it = it - 1
         endif
         if (iyp.gt.2) then
            ymax = dv**is
            ymin = -ymax
         endif
         if (ixp.gt.2) then
            xmax = dv**it
            xmin = -xmax
         endif
      endif
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
! write labels
      write (lbl,91) lblsp(iyp), lblsp(ixp), itime
      if (idproc.eq.0) then
! select point as marker
         mks = 1
! draw grid and labels, call identity transformation
         call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,&
     &label,lbl,chh)
! display sample markers with dots of minimum visible size
         dd = (smax - smin)*real(npl)*2.0e-3
         ngs = 1
!        call dsmpln(orx,smn,tmn,tmx,ngs,1,1,mks,chrs,chh)
         call dsdmpln(orx,smn,tmn,tmx,dd,ngs,1,1,mks,chrs,chh)
! define transformation number 2
         nrt = 1
! set window
         call gswn(nrt,xmin,xmax,ymin,ymax)
! set viewport
         call gsvp(nrt,smn,smx,tmn,tmx)
! select normalization transformation
         call gselnt(nrt)
! set marker size scale factor, 1.0 = nominal
         call gsmksc(1.0)
! set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mks)
! set clipping indicator, 1 = on
         call gsclip(1)
         icol = 2
         icol = kprime(icol+1)
! set polymarker color index
         call gspmci(icol)
      endif
! plot particles
! this segment is used for mpi computers
      nvp = nproc
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
! find how many particles on this node
      nd = 0
      do 30 k = 1, mxyp1
      nd = nd + kpic(k)
   30 continue
      if (idproc.eq.0) then
! copy to contiguous buffer
         noff = 0
         do 50 k = 1, mxyp1
         npp = kpic(k)
         do 40 j = 1, npp
         if ((iyp.gt.2).or.(ixp.gt.2)) then
            vx = ppart(3,j,k)
            vy = ppart(4,j,k)
            vz = ppart(5,j,k)
! vperp1 co-ordinate
            vt(1) = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
            vt(2) = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
            vt(3) = vx*ox + vy*oy + vz*oz
         endif
         if (ixp.le.2) then
            vx = ppart(ixp,j,k)
         else
            vx = vt(ixp-2)
         endif
         if (iyp.le.2) then
            vy = ppart(iyp,j,k)
         else
            vy = vt(iyp-2)
         endif
         f(j+noff) = vx
         f(nd+j+noff) = vy
   40    continue
         noff = noff + npp
   50    continue
         do 60 i = 1, nvp
! processor 0 plots his (or her) own data
         if (i.gt.1) then
! collects data from remaining nodes
            id = i - joff
            call MPI_RECV(f,2*npmax,mreal,id,103,lworld,istatus,ierr)
! determine how many particles to plot
            call MPI_GET_COUNT(istatus,mreal,nd,ierr)
            nd = nd/2
         endif
         if (nd.eq.0) go to 60
! treat dots by drawing a line to itself
!        call spdots(f(1),f(nd+1),nd,icol,nd)
! treat dots by drawing a line to itself with non-zero width in x
         dd = (xmax - xmin)*real(npl)*2.0e-3
         call spddots(f(1),f(nd+1),dd,nd,icol,nd)
! update workstation, perform
         call guwk(idwk,1)
   60    continue
! other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff - 1
! copy to contiguous buffer
         noff = 0
         do 80 k = 1, mxyp1
         npp = kpic(k)
         do 70 j = 1, npp
         if ((iyp.gt.2).or.(ixp.gt.2)) then
            vx = ppart(3,j,k)
            vy = ppart(4,j,k)
            vz = ppart(5,j,k)
! vperp1 co-ordinate
            vt(1) = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
            vt(2) = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
            vt(3) = vx*ox + vy*oy + vz*oz
         endif
         if (ixp.le.2) then
            vx = ppart(ixp,j,k)
         else
            vx = vt(ixp-2)
         endif
         if (iyp.le.2) then
            vy = ppart(iyp,j,k)
         else
            vy = vt(iyp-2)
         endif
         f(j+noff) = vx
         f(nd+j+noff) = vy
   70    continue
         noff = noff + npp
   80    continue
         call MPI_SEND(f,2*nd,mreal,0,103,lworld,ierr)
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
         call PPBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
! reset defaults
      iclr = 0
      iupd = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dsdmpln(orgx,smin,tmin,tmax,dd,ngs,nlts,nmks,mks,chrs, &
     &chh)
! this subroutine displays line or marker samples, with short character
! labels placed underneath
! dots are displayed with with non-zero width in x
! orgx = x origin of window
! smin = minimum x value of plotting window
! tmin/tmax = range of y values of plotting window
! dd = smallest visible size, such as (xmax - xmin)*4.0e-3
! ngs = number of subarrays being plotted
! nlts = number of line types available
! nmks = number of markers available
! mks = flag to determine whether lines or markers are used
! chrs = array of ngs short character labels for line or marker samples,
! each should be less than or equal to 10 characters in length
! chh = character height
      implicit none
      integer ngs, nlts, nmks, mks
      real orgx, smin, tmin, tmax, dd, chh
      character(len=*) chrs(ngs)
! rx, ry = ndc coordinates of upper-right corner of workstation window
! ncols = number of foreground colors available for line plotting
! kprime = table of color indices for prime colors
! local data
      integer k, mkr, icol, ltype, mtype
      real stx, ax, ay, dx, dy
      real x, y
! x,y = scratch variables for plotting
      dimension x(4), y(2)
! omit samples if there is only one curve
      if (ngs.le.1) return
      mkr = abs(mks)
! set marker size scale factor, 1.0 = nominal
      if (mkr.ne.0) call gsmksc(1.0)
! draw line or marker samples
      stx = .4*chh*(rx/ry)
      x(1) = orgx + 4.*stx
      x(2) = smin - 4.*stx
      dx = (x(2) - x(1))/3.
      x(3) = x(1) + dx
      x(4) = x(2) - dx
      ax = orgx + 2.*stx
      dy = (tmax - tmin - 4.*chh)/real(ngs + 1)
      ay = tmax - chh
! draw samples, first cycle through colors, then line or marker types
      do 30 k = 1, ngs
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
      y(1) = ay - dy*real(k)
      y(2) = y(1)
      if ((mks.eq.0).or.((mks.lt.0).and.(k.eq.1))) then
! set polyline color index
! 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gsplci(icol)
         ltype = 1 + (k - 1)/ncols - nlts*((k - 1)/(nlts*ncols))
! set linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
         call gsln(ltype)
! draw polyline
         call gpl(2,x,y)
      else
! set polymarker color index
! 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
         mtype = mkr + (k - 1)/ncols - nmks*((mkr - 1 + (k - 1)/ncols)/ &
     &nmks)
! set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mtype)
! dots
         if (mtype.eq.1) then
! treat dots by drawing a line to itself
!           call spdots(x(3),y,2,icol,3)
! treat dots by drawing a line to itself with non-zero width in x
            call spddots(x(3),y,dd,2,icol,3)
         else
! draw polymarker
            call gpm(2,x(3),y)
         endif
      endif
! draw label underneath sample line or marker
      y(1) = y(1) - 2.*chh
! draw text
      call gtx(ax,y(1),chrs(k))
   30 continue
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine spddots(x,y,dd,npt,icol,nxbs)
! this subroutine draws dot markers by drawing a line to itself
! with non-zero width in x
! x, y = arrays to be plotted
! dd = smallest visible size, such as (xmax - xmin)*4.0e-3
! npt = number of points to be plotted
! icol = color index
! nxbs = dimension of x, y arrays
      implicit none
      integer npt, icol, nxbs
      real x, y, dd
      dimension x(nxbs), y(nxbs)
! local data
      integer j
      real xs, ys
! xs, ys = scratch arrays for plotting
      dimension xs(2), ys(2)
! set polyline color index
      call gsplci(icol)
      do 10 j = 1, npt
      xs(1) = x(j)
      ys(1) = y(j)
      xs(2) = xs(1) + dd
      ys(2) = ys(1)
! draw polyline
      call gpl(2,xs,ys)
   10 continue
      end subroutine
!
      end module

