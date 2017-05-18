c general 2d gks graphics library
c written by viktor k. decyk, ucla
c copyright 1996, regents of the university of california
c update: february 7, 2017
      subroutine CARPET(f,label,isc,ist,nx,ny,nxv,chr,ntc,irc)
c this subroutine displays an array f as a color raster image.
c a 256 color palette must have been defined prior to this call. 
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = field array to be plotted
c label = long character string label for plot
c isc = power of 2 scale of range of values of f
c ist = flag for choosing positive and/or negative values
c the range of values of f are given by fmax and fmin.
c if ist = 0, then fmax = 2**isc and fmin = -2**isc.
c if ist = 1, then fmax = 2**isc and fmin = 0.
c if ist = -1, then fmax = 0 and fmin = -2**isc.
c if ist = 2, then fmax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c chr = additional long character string comment for plot
c ntc = number of valid colors, should be power of 2, <= 256
c irc = return code (0 = normal return)
c npald = number of palette entries
      parameter(npald=256)
c lxm, lym = maximum number of pixels in x, y
      parameter(lxm=360,lym=270)
      character*(*) label, chr
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      dimension f(nxv,ny)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img8 = integer image array
      dimension img8(lxm*lym)
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /11,11/
c csize = vertical size of characters
      data csize /0.034/
c istyle = (0,1) = contour plot (fills area,preserves aspect ratio)
      data istyle /1/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c find scales for plot
      is = isc
      if (abs(is).gt.116) then
         fmax = f(1,1)
         fmin = fmax
         do 20 k = 1, ny
         do 10 j = 1, nx
         fmax = amax1(fmax,f(j,k))
         fmin = amin1(fmin,f(j,k))
   10    continue
   20    continue
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
      if (ist.eq.1) then
         fmin = 0.
      else if (ist.eq.(-1)) then
         fmax = 0.  
      endif
c parameters for plots
      xmin = 0.
      xmax = real(nx - 1)
      ymin = 0.
      ymax = real(ny - 1)
c find location for plot
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
c fill area
      xmn = smn
      ymn = tmn
c preserve aspect ratio
      if (istyle.eq.1) then
         if (nx.gt.ny) ymn = tmx - (tmx - tmn)*real(ny)/real(nx)
         if (ny.gt.nx) xmn = smx - (smx - smn)*real(nx)/real(ny)
      endif
c set size of raster image
      lxs = lxm
      lys = lym
      lxp = min(2,lxs)
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c create sample palette
      ac = real(ntc - 1)/real(lys - 1)
c rescale factor ncv
      ncv = 256/ntc
      do 40 k = 1, lys
      joff = lxp*(k - 1)
      ic = ac*real(lys - k) + 0.999999
c rescale index
      ic = ic*ncv
c lookup table required
      if (lupt.eq.1) ic = ipal(ic+1)
      do 30 j = 1, lxp
      img8(j+joff) = ic
   30 continue
   40 continue
c draw grid and labels, call identity transformation
      call tickz(xmin,xmax,ymin,ymax,orx,ory,xmn,smx,ymn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample palette
      call dsmpal(img8,fmin,fmax,orx,smn,tmn,tmx,lxp,lys,chh)
c map f to color raster image
      call mraster(f,img8,fmin,fmax,nx,ny,nxv,lxs,lys,ntc)
c copy with lookup table
      if (lupt.eq.1) then
         do 60 k = 1, lys
         joff = lxs*(k - 1)
         do 50 j = 1, lxs
         img8(j+joff) = ipal(img8(j+joff)+1)
   50    continue
   60    continue
      endif
c cell array
c special case for rs/6000 with graPHIGS gks
c     call gca(xmn,tmx,smx,ymn,lys,lxs,1,1,lys,lxs,img8)
      call gca(xmn,tmx,smx,ymn,lxs,lys,1,1,lxs,lys,img8)
c update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         call readrc(irc)
      endif
c reset defaults
      iclr = 0
      iupd = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine CARPETL(f,label,xmin,xmax,ymin,ymax,isc,ist,nx,ny,nxv,c
     1hr,ntc,irc)
c this subroutine displays an array f as a color raster image.
c a 256 color palette must have been defined prior to this call. 
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = field array to be plotted
c label = long character string label for plot
c xmin/xmax = numerical labels for x axis
c ymin/ymax = numerical labels for y axis
c isc = power of 2 scale of range of values of f
c ist = flag for choosing positive and/or negative values
c the range of values of f are given by fmax and fmin.
c if ist = 0, then fmax = 2**isc and fmin = -2**isc.
c if ist = 1, then fmax = 2**isc and fmin = 0.
c if ist = -1, then fmax = 0 and fmin = -2**isc.
c if ist = 2, then fmax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c chr = additional long character string comment for plot
c ntc = number of valid colors, should be power of 2, <= 256
c irc = return code (0 = normal return)
c npald = number of palette entries
      parameter(npald=256)
c lxm, lym = maximum number of pixels in x, y
      parameter(lxm=360,lym=270)
      character*(*) label, chr
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      dimension f(nxv,ny)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img8 = integer image array
      dimension img8(lxm*lym)
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /11,11/
c csize = vertical size of characters
      data csize /0.034/
c istyle = (0,1) = contour plot (fills area,preserves aspect ratio)
      data istyle /0/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c find scales for plot
      is = isc
      if (abs(is).gt.116) then
         fmax = f(1,1)
         fmin = fmax
         do 20 k = 1, ny
         do 10 j = 1, nx
         fmax = amax1(fmax,f(j,k))
         fmin = amin1(fmin,f(j,k))
   10    continue
   20    continue
         if (fmax.eq.0.) fmax = 1.0e-35
         rmax = fmax - fmin
         if (rmax.eq.0.) rmax = 1.0e-35
         rmin = fmin
         gmax = abs(fmax)
         is = alog(gmax)*algdvi
         if (gmax.ge.1.) is = is + 1
         if (gmax.le.dv**(is-1)) is = is - 1
         gmin = abs(fmin)  
         if (gmin.gt.0.) then
            it = alog(gmin)*algdvi
            if (gmin.ge.1.) it = it + 1
            if (gmin.le.dv**(it-1)) it = it - 1
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
            if (gmin.gt.gmax) then
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
      if (ist.eq.1) then
         fmin = 0.
      else if (ist.eq.(-1)) then
         fmax = 0.  
      endif
c find location for plot
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
c fill area
      xmn = smn
      ymn = tmn
c preserve aspect ratio
      if (istyle.eq.1) then
         if (nx.gt.ny) ymn = tmx - (tmx - tmn)*real(ny)/real(nx)
         if (ny.gt.nx) xmn = smx - (smx - smn)*real(nx)/real(ny)
      endif
c set size of raster image
      lxs = lxm
      lys = lym
      lxp = min(2,lxs)
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c create sample palette
      ac = real(ntc - 1)/real(lys - 1)
c rescale factor ncv
      ncv = 256/ntc
      do 40 k = 1, lys
      joff = lxp*(k - 1)
      ic = ac*real(lys - k) + 0.999999
c rescale index
      ic = ic*ncv
c lookup table required
      if (lupt.eq.1) ic = ipal(ic+1)
      do 30 j = 1, lxp
      img8(j+joff) = ic
   30 continue
   40 continue
c draw grid and labels, call identity transformation
      call tickl(xmin,xmax,ymin,ymax,orx,ory,xmn,smx,ymn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample palette
      call dsmpal(img8,fmin,fmax,orx,smn,tmn,tmx,lxp,lys,chh)
c map f to color raster image
      call lraster(f,img8,fmin,fmax,nx,ny,nxv,lxs,lys,ntc)
c copy with lookup table
      if (lupt.eq.1) then
         do 60 k = 1, lys
         joff = lxs*(k - 1)
         do 50 j = 1, lxs
         img8(j+joff) = ipal(img8(j+joff)+1)
   50    continue
   60    continue
      endif
c cell array
c special case for rs/6000 with graPHIGS gks
c     call gca(xmn,tmx,smx,ymn,lys,lxs,1,1,lys,lxs,img8)
      call gca(xmn,tmx,smx,ymn,lxs,lys,1,1,lxs,lys,img8)
c update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         call readrc(irc)
      endif
c reset defaults
      iclr = 0
      iupd = 0
      return
      end
      subroutine dsmpal(img,fmin,fmax,orgx,smin,tmin,tmax,lx,ly,chh)
c this subroutine displays sample palette, labeled with fmin and fmax.
c img = color index array which contains sample palette
c fmin/fmax = range of f values in plot
c orgx = x origin of window
c smin = minimum x value of plotting window
c tmin/tmax = range of y values of plotting window
c lx,ly = number of columns and rows in the cell array
c chh = character height
      dimension img(lx,ly)
      character*12 lbl
   91 format (' f=',e9.3)
c label f values
      dx = .33*(smin - orgx)
      ax = orgx
      write (lbl,91) fmin
      ay = tmin + 2.*chh
c draw text
      call gtx(ax,ay,lbl)
      write (lbl,91) fmax
      ay = tmax - 3.*chh
c draw text
      call gtx(ax,ay,lbl)
c edges of sample palette
      ax = orgx + dx
      bx = ax + dx
      ay = tmin + 4.*chh
      by = tmax - 4.*chh
c cell array
c special case for rs/6000 with graPHIGS gks
c     call gca(ax,by,bx,ay,ly,lx,1,1,ly,lx,img)
      call gca(ax,by,bx,ay,lx,ly,1,1,lx,ly,img)
      return
      end
      subroutine STPALIT(idpal)
c this subroutine selects one from three available palettes
c idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
      parameter(npald=256,lpald=3*npald)
      character*1 pal(lpald)
c nbit = the number of colors, pixel depth
      data nbit /8/
c ntc = number of colors requested
      ntc = 2**nbit
c select palette
      call gitpal(pal,idpal,lpald,ntc)
c set palette
      call grspal(nbit,pal,npald,ntc,lpald)
      return
      end
      subroutine gitpal(pal,idpal,npald,ipmx)
c this subroutine retrieves various palettes
c pal is a character array, with rgb values in successive bytes
c idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
c cold/hot: blue,green,cyan,foreground,yellow,magenta,red
c color wheel: magenta, blue, cyan, green, yellow, red, magenta
c rainbow: violet, indigo, blue, green, yellow, orange, red
c npald = three times the number of palette entries
c ipmx = maximum color value in the palette
c maxpal = maximum number of palettes
c ncgr = number of segments in palette
      parameter (maxpal=3,ncgr=7)
      parameter (ncgd=maxpal*ncgr)
      character*1 pal(npald)
      dimension red(ncgd), green(ncgd), blue(ncgd), is(ncgd)
      dimension iss(ncgr), reds(ncgr), greens(ncgr), blues(ncgr)
      save is, red, green, blue
c index segments
      data is /1,43,85,127,170,212,254,1,43,85,127,170,212,254,1,28,85,1
     141,198,254,255/
c color segments
      data red /0.,0.,0.,1.,1.,1.,1.,1.,0.,0.,0.,1.,1.,1.,0.5,0.,0.,0.,1
     1.,1.,1./
      data green /0.,1.,1.,1.,1.,0.,0.,0.,0.,1.,1.,1.,0.,0.,0.,0.,1.,1.,
     11.,0.,1./
      data blue /1.,0.,1.,1.,0.,1.,0.,1.,1.,1.,0.,0.,0.,1.,1.,1.,1.,0.,0
     1.,0.,1./
c return if palette id number is out of range of
      if ((idpal.lt.1).or.(idpal.gt.maxpal)) return
c return if palette size is too small
      if (npald.lt.768) return
c normalize color amplitudes
      apmx = real(ipmx - 1)
c extract index and color segments for palette selected
      do 10 i = 1, ncgr
      ioff = (idpal - 1)*ncgr
      iss(i) = is(i+ioff)
      reds(i) = red(i+ioff)*apmx
      greens(i) = green(i+ioff)*apmx
      blues(i) = blue(i+ioff)*apmx
   10 continue
c pad with background
      it = 3*iss(1)
      do 20 i = 1, it
      pal(i) = char(0)
   20 continue
c first point
      pal(it+1) = char(int(reds(1)))
      pal(it+2) = char(int(greens(1)))
      pal(it+3) = char(int(blues(1)))
c outer loop over sextants
      do 40 i = 2, ncgr
      it = iss(i) - iss(i-1)
      at = 1./real(it)
c inner loop does linear interpolation between segments
      do 30 j = 1, it
      ic = 3*(j + iss(i-1)) + 1
      att = real(j)*at
      pal(ic) = char(int((reds(i) - reds(i-1))*att + reds(i-1)))
      pal(ic+1) = char(int((greens(i) - greens(i-1))*att + greens(i-1)))
      pal(ic+2) = char(int((blues(i) - blues(i-1))*att + blues(i-1)))
   30 continue
   40 continue
c pad with foreground
      it = 3*iss(ncgr) + 4
      do 50 i = it, 768
      pal(i) = char(int(apmx))
   50 continue
   60 return
      end
      subroutine mraster(f,icola,fmin,fmax,nx,ny,nxv,lx,ly,ntc)
c this subroutine maps a floating point array to a color raster image.
c the array f is mapped to a pixel array, where the values of f between
c fmin and fmax are linearly interpolated to color index values between
c between 1 and ntc - 1, rescaled to values between 1 and 255, and then
c stored in the cell array icola.
c values of f outside the range are mapped to index 0.
c f = field array to be converted
c icola = color index array
c fmin/fmax = minimum/maximum range in f to be included
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c lx,ly = number of columns and rows in the cell array
c ntc = number of valid colors, should be power of 2, <= 256
      dimension f(nxv,ny)
      dimension icola(lx,ly)
      lx1 = lx - 1
      ly1 = ly - 1
c interpolation scales
      dxs = real(nx - 1)/real(lx1)
      dys = real(ny - 1)/real(ly1)
c color scale
      af = real(ntc - 1)/(fmax - fmin)
c rescale factor ncv
      ncv = 256/ntc
c loop over y pixels
      do 20 k = 1, ly1
      k1 = ly + 1 - k
c yt = location in f array corresponding to pixel k
      yt = dys*real(k - 1)
      m = yt
c dy = linear interpolation weight in y
      dy = yt - real(m)
c m = address in f array below or at pixel k
      m = m + 1
      dyt = 1. - dy
c loop over x pixels
      do 10 j = 1, lx1
c xt = location in f array corresponding to pixel j
      xt = dxs*real(j - 1)
      n = xt
c dx = linear interpolation weight in x
      dx = xt - real(n)
c n = address in f array below or at pixel j
      n = n + 1
      dxt = 1. - dx
c fc = interpolated f value at pixel (j,k)
      fc = dyt*(dxt*f(n,m) + dx*f(n+1,m)) + dy*(dxt*f(n,m+1) + dx*f(n+1,
     1m+1))
c find index if fc within range
      if ((fc.ge.fmin).and.(fc.lt.fmax)) then
         ic = (fc - fmin)*af + 1.0
         if (ic.gt.ntc) ic = ntc
c set index to background if out of range
      else
         ic = 0
      endif
c rescale and store index
      icola(j,k1) = ic*ncv
   10 continue
c fc = interpolated f value at pixel (nx,k)
      fc = dyt*f(nx,m) + dy*f(nx,m+1)
c find index if fc within range
      if ((fc.ge.fmin).and.(fc.lt.fmax)) then
         ic = (fc - fmin)*af + 1.0
         if (ic.gt.ntc) ic = ntc
c set index to background if out of range
      else
         ic = 0
      endif
c rescale and store index
      icola(lx,k1) = ic*ncv
   20 continue
c loop over x pixels for top row
      do 30 j = 1, lx1
c xt = location in f array corresponding to pixel j
      xt = dxs*real(j - 1)
      n = xt
c dx = linear interpolation weight in x
      dx = xt - real(n)
c n = address in f array below or at pixel j
      n = n + 1
      dxt = 1. - dx
c fc = interpolated f value at pixel (j,ny)
      fc = dxt*f(n,ny) + dx*f(n+1,ny)
c find index if fc within range
      if ((fc.ge.fmin).and.(fc.lt.fmax)) then
         ic = (fc - fmin)*af + 1.0
         if (ic.gt.ntc) ic = ntc
c set index to background if out of range
      else
         ic = 0
      endif
c rescale and store index
      icola(j,1) = ic*ncv
   30 continue
c fc = interpolated f value at pixel (nx,ny)
      fc = f(nx,ny)
c find index if fc within range
      if ((fc.ge.fmin).and.(fc.lt.fmax)) then
         ic = (fc - fmin)*af + 1.0
         if (ic.gt.ntc) ic = ntc
c set index to background if out of range
      else
         ic = 0
      endif
c rescale and store index
      icola(lx,1) = ic*ncv
      return
      end
      subroutine lraster(f,icola,fmin,fmax,nx,ny,nxv,lx,ly,ntc)
c this subroutine maps a floating point array to a color raster image.
c the array f is mapped to a pixel array, where the values of f between
c fmin and fmax are linearly interpolated to color index values between
c between 1 and ntc - 1, rescaled to values between 1 and 255, and then
c stored in the cell array icola.
c values of f outside the range are mapped to index 0.
c f = field array to be converted
c icola = color index array
c fmin/fmax = minimum/maximum range in f to be included
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c lx,ly = number of columns and rows in the cell array
c ntc = number of valid colors, should be power of 2, <= 256
      dimension f(nxv,ny)
      dimension icola(lx,ly)
      lx1 = lx - 1
      ly1 = ly - 1
c interpolation scales
      dxs = real(nx - 1)/real(lx1)
      dys = real(ny - 1)/real(ly1)
c color scale
      af = real(ntc - 1)/(fmax - fmin)
c rescale factor ncv
      ncv = 256/ntc
c loop over y pixels
      do 20 k = 1, ly1
      k1 = ly + 1 - k
c yt = location in f array corresponding to pixel k
      yt = dys*real(k - 1)
      m = yt
c dy = linear interpolation weight in y
      dy = yt - real(m)
c m = address in f array below or at pixel k
      m = m + 1
      dyt = 1. - dy
c loop over x pixels
      do 10 j = 1, lx1
c xt = location in f array corresponding to pixel j
      xt = dxs*real(j - 1)
      n = xt
c dx = linear interpolation weight in x
      dx = xt - real(n)
c n = address in f array below or at pixel j
      n = n + 1
      dxt = 1. - dx
c fc = interpolated f value at pixel (j,k)
      fc = dyt*(dxt*f(n,m) + dx*f(n+1,m)) + dy*(dxt*f(n,m+1) + dx*f(n+1,
     1m+1))
c find index if fc within range
      if ((fc.ge.fmin).and.(fc.lt.fmax)) then
         ic = (fc - fmin)*af + 1.0
         if (ic.gt.ntc) ic = ntc
c set index to first if out of range
      else
         ic = 1
      endif
c rescale and store index
      icola(j,k1) = ic*ncv
   10 continue
c fc = interpolated f value at pixel (nx,k)
      fc = dyt*f(nx,m) + dy*f(nx,m+1)
c find index if fc within range
      if ((fc.ge.fmin).and.(fc.lt.fmax)) then
         ic = (fc - fmin)*af + 1.0
         if (ic.gt.ntc) ic = ntc
c set index to first if out of range
      else
         ic = 1
      endif
c rescale and store index
      icola(lx,k1) = ic*ncv
   20 continue
c loop over x pixels for top row
      do 30 j = 1, lx1
c xt = location in f array corresponding to pixel j
      xt = dxs*real(j - 1)
      n = xt
c dx = linear interpolation weight in x
      dx = xt - real(n)
c n = address in f array below or at pixel j
      n = n + 1
      dxt = 1. - dx
c fc = interpolated f value at pixel (j,ny)
      fc = dxt*f(n,ny) + dx*f(n+1,ny)
c find index if fc within range
      if ((fc.ge.fmin).and.(fc.lt.fmax)) then
         ic = (fc - fmin)*af + 1.0
         if (ic.gt.ntc) ic = ntc
c set index to first if out of range
      else
         ic = 1
      endif
c rescale and store index
      icola(j,1) = ic*ncv
   30 continue
c fc = interpolated f value at pixel (nx,ny)
      fc = f(nx,ny)
c find index if fc within range
      if ((fc.ge.fmin).and.(fc.lt.fmax)) then
         ic = (fc - fmin)*af + 1.0
         if (ic.gt.ntc) ic = ntc
c set index to first if out of range
      else
         ic = 1
      endif
c rescale and store index
      icola(lx,1) = ic*ncv
      return
      end
      subroutine CONTUR(f,lf,label,isc,ist,nx,ny,nxv,chr,nc,irc)
c this subroutine displays an array f as a contour plot.
c a maximum of ncols colors are used, used in order from lowest to
c highest contour: blue, green, cyan, foreground, yellow, magenta, red
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = field array to be plotted
c lf = scratch field array
c label = long character string label for plot
c isc = power of 2 scale of range of values of f
c ist = flag for choosing positive and/or negative values
c the range of values of f are given by fmax and fmin.
c if ist = 0, then fmax = 2**isc and fmin = -2**isc.
c if ist = 1, then fmax = 2**isc and fmin = 0.
c if ist = -1, then fmax = 0 and fmin = -2**isc.
c if ist = 2, then fmax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c chr = additional long character string comment for plot
c nc = number of contour lines
c irc = return code (0 = normal return)
      character*(*) label, chr
c idwk = workstation identifier
c ncols = number of foreground colors available for line plotting
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      dimension f(nxv,ny), lf(nxv,ny)
c icolor = color index, used in order from lowest to highest contour
      dimension icv(8), icolor(8)
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /11,11/
c csize = vertical size of characters
      data csize /0.034/
c icv = location in kprime array of color indices needed for icolor
      data icv /1,3,8,6,2,5,7,4/
c istyle = (0,1) = contour plot (fills area,preserves aspect ratio)
      data istyle /1/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c find scales for plot
      is = isc
      if (abs(is).gt.116) then
         fmax = f(1,1)
         fmin = fmax
         do 20 k = 1, ny
         do 10 j = 1, nx
         fmax = amax1(fmax,f(j,k))
         fmin = amin1(fmin,f(j,k))
   10    continue
   20    continue
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
      if (ist.eq.1) then
         fmin = 0.
      else if (ist.eq.(-1)) then
         fmax = 0.  
      endif
c parameters for plots
      xmin = 0.
      xmax = real(nx - 1)
      ymin = 0.
      ymax = real(ny - 1)
c find location for plot
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
c fill area
      xmn = smn
      ymn = tmn
c preserve aspect ratio
      if (istyle.eq.1) then
         if (nx.gt.ny) ymn = tmx - (tmx - tmn)*real(ny)/real(nx)
         if (ny.gt.nx) xmn = smx - (smx - smn)*real(nx)/real(ny)
      endif
c calculate color indices
      ntc = ncols + 1
      do 30 i = 1, ntc
      icolor(i) =  kprime(icv(i))
   30 continue
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c draw grid and labels, call identity transformation
      call tickz(xmin,xmax,ymin,ymax,orx,ory,xmn,smx,ymn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample contours
      call dsmpcn(icolor,fmin,fmax,orx,smn,tmn,tmx,chh,ntc,nc)
c draw contour map
      call dcontr(f,lf,icolor,fmin,fmax,xmn,smx,ymn,tmx,nx,ny,nxv,ntc,nc
     1)
c update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         call readrc(irc)
      endif
c reset defaults
      iclr = 0
      iupd = 0
      return
      end
      subroutine dsmpcn(icolor,fmin,fmax,orgx,smin,tmin,tmax,chh,ntc,nc)
c this subroutine displays sample contour, labeled with fmin and fmax.
c icolor = color index array, used in order from lowest to highest
c if colors are not available, line styles are used instead
c fmin/fmax = range of f values in plot
c orgx = x origin of window
c smin = minimum x value of plotting window
c tmin/tmax = range of y values of plotting window
c chh = character height
c ntc = number of valid colors, should be in the range 2 <= ntc <= 8
c nc = number of contour lines
      dimension icolor(ntc)
      character*12 lbl
c x,y = scratch variables for plotting
      dimension x(2), y(2)
      dimension lns(4)
c lns = line style table, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      data lns /3,3,1,1/
   91 format (' f=',e9.3)
c label f values
      dx = .33*(smin - orgx)
      ax = orgx
      write (lbl,91) fmin
      ay = tmin + 2.*chh
c draw text
      call gtx(ax,ay,lbl)
      write (lbl,91) fmax
      ay = tmax - 3.*chh
c draw text
      call gtx(ax,ay,lbl)
c edges of sample contour
      x(1) = orgx + dx
      x(2) = x(1) + dx
      ay = tmin + 4.*chh
      by = tmax - 4.*chh 
c increment between contours
      dy = (by - ay)/real(nc + 1)
c increment between colors
      if (ntc.gt.2) then
         ac = real(ntc - 1)/real(nc)
c increment between line styles
      else
         ac = 4./real(nc)
      endif
      do 10 l = 1, nc
c calculate location
      y(1) = ay + dy*real(l)
      y(2) = y(1)
      il = ac*(real(l) - .5)
c set color
      if (ntc.gt.2) then
         icol = icolor(il+2)
c set polyline color index
         call gsplci(icol)
c set line style
      else
c set linetype
         call gsln(lns(il+1))
      endif
c draw polyline
      call gpl(2,x,y)
   10 continue
      return
      end
      subroutine dcontr(f,lf,icolor,fmin,fmax,smin,smax,tmin,tmax,nx,ny,
     1nxv,ntc,nc)
c this subroutine draws contour map of a floating point array, for 
c values of f between fmin and fmax.
c f = field array to be mapped
c lf = scratch array for storing locations where contour crosses grid,
c lf = (0,1,2,3) = contour crossed grid in (no,x,y,x and y) direction
c icolor = color index array, used in order from lowest to highest
c if colors are not available, line styles are used instead
c fmin/fmax = minimum/maximum range in f to be included
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c ntc = number of valid colors, should be in the range 2 <= ntc <= 8
c nc = number of contour lines, equally spaced between fmax and fmin
c nxbs = length of scratch variable for plotting (must be >=2)
      parameter(nxbs=65)
      logical ld, li
      dimension f(nxv,ny), lf(nxv,ny+1)
      dimension icolor(ntc)
c x,y = scratch arrays for plotting
      dimension x(nxbs), y(nxbs)
      dimension lns(4)
c lns = line style table, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      data lns /3,3,1,1/
      nx1 = nx - 1
      ny1 = ny - 1
c increment between contours
      dc = (fmax - fmin)/real(nc)
c increment between colors
      if (ntc.gt.2) then
         ac = real(ntc - 1)/real(nc)
c increment between line styles
      else
         ac = 4./real(nc)
      endif
c first contour
      ci = fmin + .5*dc
c scale factors for plot
      dsx = (smax - smin)/real(nx1)
      dsy = (tmax - tmin)/real(ny1)
c loop over contour levels
      do 120 l = 1, nc
      clev = ci + dc*real(l - 1)
      il = ac*(real(l) - .5)
c set color
      if (ntc.gt.2) then
         icol = icolor(il+2)
c set polyline color index
         call gsplci(icol)
c set line style
      else
c set linetype
         call gsln(lns(il+1))
      endif
c find if contour crosses x grid
      do 20 k = 1, ny
      do 10 j = 1, nx1
      ld = f(j,k).gt.clev
      li = f(j+1,k).gt.clev
      if ((ld.and.(.not.li)).or.((.not.ld).and.li)) then
         lf(j,k) = 1
      else
         lf(j,k) = 0
      endif
   10 continue
      lf(nx,k) = 0
   20 continue
c find if contour crosses y grid
      do 40 k = 1, ny1
      do 30 j = 1, nx
      ld = f(j,k).gt.clev
      li = f(j,k+1).gt.clev
      if ((ld.and.(.not.li)).or.((.not.ld).and.li)) then
         lf(j,k) = lf(j,k) + 2
      endif
   30 continue
   40 continue
c start following grid crossings
c first follow contours which intersect the four edges (i=1,4)
c then follow interior contours  (i=5)
      do 90 i = 1, 5
c idir = (0,1,2,3) = (bottom,right,top,left) edge
c idir = 5 = interior
      idir = i - 1
c initial counter
      n = 0
c initial crossing direction is undefined
      ndir = -1
c set initial offset
      joff = 0
      koff = 0
c set initial loop extent
      nlx = nx1
      nly = ny1
c special cases of edges
c bottom or top edge
      if ((idir.eq.0).or.(idir.eq.2)) then
         nly = 1
         if (idir.eq.2) koff = ny1
c right or left edge
      elseif ((idir.eq.1).or.(idir.eq.3)) then
         nlx = 1
         if (idir.eq.1) joff = nx1
      endif
c begin scanning grids
c in y direction
      do 80 kl = 1, nly
      k = kl + koff
c in x direction
      do 70 jl = 1, nlx
      j = jl + joff
c find initial location where contour crossed some grid
c set initial coordinates, crossing direction and marker removal flag
c ndir = (-1,0,1,2,3) = contour path is (undefined,up,left,down,right)
c iflg = (0,1) = (no,yes) remove grid crossing location from lf array
      if (lf(j,k).ne.0) then
         jj = j
         kk = k
c interior contour 
         if (idir.eq.4) then
c contour crossed y grid, going right
            if (lf(j,k).eq.2) then
               ndir = 3
c contour crossed x grid, going up
            else
               ndir = 0
            endif
            iflg = 0
            go to 60
c contour crosses edge
         elseif (((nlx.eq.1).and.(lf(j,k).ge.2)).or.((nly.eq.1).and.(lf(
     1j,k).ne.2))) then
            ndir = idir
            iflg = 1
            go to 60
c contour parallel to edge, skip until later
         else
            go to 70
         endif
c no grid crossed, try again
      else
         go to 70
      endif
c continue following contour, search for adjacent grid crossing
c if direction is up
   50 if (ndir.eq.0) then
         ndir = -1
c y grid crossed going right
         if ((jj.lt.nx).and.(lf(jj+1,kk).ge.2)) then
            jj = jj + 1
            ndir = 3
c x grid crossed going up
         elseif ((kk.lt.ny).and.((lf(jj,kk+1).eq.1).or.(lf(jj,kk+1).eq.3
     1))) then
            kk = kk + 1
            ndir = 0
c y grid crossed going left
         elseif (lf(jj,kk).ge.2) then
            ndir = 1
         endif
c if direction is left
      elseif (ndir.eq.1) then
         ndir = -1
         if (jj.gt.1) then
c x grid crossed going up
            if ((kk.lt.ny).and.((lf(jj-1,kk+1).eq.1).or.(lf(jj-1,kk+1).e
     1q.3))) then
               jj = jj - 1
               kk = kk + 1
               ndir = 0
c y grid crossed going left
            elseif (lf(jj-1,kk).ge.2) then
               jj = jj - 1
               ndir = 1
c x grid crossed going down
            elseif ((lf(jj-1,kk).eq.1).or.(lf(jj-1,kk).eq.3)) then
               jj = jj - 1
               ndir = 2
            endif
         endif
c if direction is down
      elseif (ndir.eq.2) then
         ndir = -1
         if (kk.gt.1) then
c y grid crossed going left
            if (lf(jj,kk-1).ge.2) then
               kk = kk - 1
               ndir = 1
c x grid crossed going down
            elseif ((jj.lt.nx).and.((lf(jj,kk-1).eq.1).or.(lf(jj,kk-1).e
     1q.3))) then
               kk = kk - 1
               ndir = 2
c y grid crossed going right
            elseif ((jj.lt.nx).and.(lf(jj+1,kk-1).ge.2)) then
               jj = jj + 1
               kk = kk - 1
               ndir = 3
            endif
         endif
c if direction is right
      elseif (ndir.eq.3) then
         ndir = -1
c x grid crossed going down
         if ((lf(jj,kk).eq.1).or.(lf(jj,kk).eq.3)) then
            ndir = 2
c y grid crossed going right
         elseif ((jj.lt.nx).and.(lf(jj+1,kk).ge.2)) then
            jj = jj + 1
            ndir = 3
c x grid crossed going up
         elseif ((kk.lt.ny).and.((lf(jj,kk+1).eq.1).or.(lf(jj,kk+1).eq.3
     1))) then
            kk = kk + 1
            ndir = 0
         endif
      endif
c next grid crossing not found
   60 if (ndir.lt.0) then
c draw polyline
         if (n.gt.0) call gpl(n,x,y)
c initial counter
         n = 0
c initial crossing direction is undefined
         ndir = -1
         go to 70
c x grid crossed
      elseif ((ndir.eq.0).or.(ndir.eq.2)) then
c interpolate between grids
         dx = (clev - f(jj,kk))/(f(jj+1,kk) - f(jj,kk))
         dy = 0.
c remove marker for crossing x grid
         if (iflg.eq.1) lf(jj,kk) = lf(jj,kk) - 1
c y grid crossed
      else
c interpolate between grids
         dx = 0.
         dy = (clev - f(jj,kk))/(f(jj,kk+1) - f(jj,kk))
c remove marker for crossing y grid
         if (iflg.eq.1) lf(jj,kk) = lf(jj,kk) - 2
      endif
c increment counter
      n = n + 1
c draw scratch array, if buffer would overflow
      if (n.gt.nxbs) then
c draw polyline
         call gpl(nxbs,x,y)
c store last coordinate for next time
         x(1) = x(nxbs)
         y(1) = y(nxbs)
c reset counter
         n = 2
      endif
c store next coordinate in scratch array
      x(n) = smin + dsx*(real(jj - 1) + dx)
      y(n) = tmin + dsy*(real(kk - 1) + dy)
c reset marker removal flag
      iflg = 1
c find next grid crossing
      go to 50
   70 continue
   80 continue
   90 continue
c sanity check
      isum = 0
      do 110 k = 1, ny
      do 100 j = 1, nx
      isum = isum + iabs(lf(j,k))
  100 continue
  110 continue
      if (isum.ne.0) then
         write (6,*) 'some contours missed, isum, l = ',isum,l
      endif
  120 continue
      return
      end
      subroutine tickz(xmin,xmax,ymin,ymax,orgx,orgy,smin,smax,tmin,tmax
     1,ntx,nty,label,chr,chh)
c this subroutine draws grid and labels
c first, a box is drawn around the plots, and tick marks are added.
c then, the axes are labeled with the range of x and y coordinates.
c next, long character string labels are placed underneath the plots.
c xmin/xmax = range of x values in plot
c ymin/ymax = range of y values in plot
c orgx/orgy = origin of window
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
c ntx/nty = number of ticks in grid in x/y direction
c label = long character string label for plot
c chr = additional long character string comment for plot
c chh = character height, also determines size of tick marks
c for clarity, should have 4.2*chh < tmin - orgy
      character*(*) label, chr
      character*7 lbl
c rx, ry = ndc coordinates of upper-right corner of workstation window
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c x,y = scratch variables for plotting
      dimension x(5), y(5)
c  91 format (f7.1)
   91 format (i5)
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 2 = stroke precision
      call dfplps (chh,1,2)
c draw box around plots
      x(1) = smin
      x(2) = smax
      x(3) = smax
      x(4) = smin
      x(5) = smin
      y(1) = tmin
      y(2) = tmin
      y(3) = tmax
      y(4) = tmax
      y(5) = tmin
c draw polyline
      call gpl(5,x,y)
c draw ticks
      dx = 0.
      if (ntx.gt.1) dx = (smax - smin)/real(ntx - 1)
c sty = size of tick mark in y direction
      sty = .4*chh
c draw ticks along x axis
      y(1) = tmin
      y(2) = tmin - sty
      y(3) = tmax
      y(4) = tmax + sty
      do 10 j = 1, ntx
      x(1) = smin + dx*(j - 1)
      x(2) = x(1)
c draw polylines
      call gpl(2,x,y)
      call gpl(2,x,y(3))
   10 continue
c draw ticks along y axis
      dy = 0.
      if (nty.gt.1) dy = (tmax - tmin)/real(nty - 1)
c stx = size of tick marks in x direction
      stx = sty*(rx/ry)
      x(1) = smin
      x(2) = smin - stx
      x(3) = smax
      x(4) = smax + stx
      do 20 j = 1, nty
      y(1) = tmin + dy*(j - 1)
      y(2) = y(1)
c draw polylines
      call gpl(2,x,y)
      call gpl(2,x(3),y)
   20 continue
c label x axes
      ay = tmin - 1.6*chh
      write (lbl,91) int(xmin)
      call trimc(lbl,is,ls)
      ax = smin
c draw text
      call gtx(ax,ay,lbl(is:ls))
c set text alignment to (right,normal)
      call gstxal(3,0)
      write (lbl,91) int(xmax)
      call trimc(lbl,is,ls)
      ax = smax
c draw text
      call gtx(ax,ay,lbl(is:ls))
c label y axes
      ax = smin - 2.*stx
      write (lbl,91) int(ymin)
      call trimc(lbl,is,ls)
      ay = tmin
c draw text
      call gtx(ax,ay,lbl(is:ls))
      write (lbl,91) int(ymax)
      call trimc(lbl,is,ls)
      ay = tmax - chh
c draw text
      call gtx(ax,ay,lbl(is:ls))
c set text alignment to (normal,normal)
      call gstxal(0,0)
c write character string labels
      ax = orgx
      ay = orgy + 1.35*chh
c draw text
      call gtx(ax,ay,label)
      ay = orgy + .25*chh
c draw text
      call gtx(ax,ay,chr)
      return
      end
      subroutine tickl(xmin,xmax,ymin,ymax,orgx,orgy,smin,smax,tmin,tmax
     1,ntx,nty,label,chr,chh)
c this subroutine draws grid and labels
c first, a box is drawn around the plots, and tick marks are added.
c then, the axes are labeled with the range of x and y coordinates.
c next, long character string labels are placed underneath the plots.
c xmin/xmax = range of x values in plot
c ymin/ymax = range of y values in plot
c orgx/orgy = origin of window
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
c ntx/nty = number of ticks in grid in x/y direction
c label = long character string label for plot
c chr = additional long character string comment for plot
c chh = character height, also determines size of tick marks
c for clarity, should have 4.2*chh < tmin - orgy
      character*(*) label, chr
      character*7 lbl
c rx, ry = ndc coordinates of upper-right corner of workstation window
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c x,y = scratch variables for plotting
      dimension x(5), y(5)
   91 format (f7.1)
c  91 format (i5)
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 2 = stroke precision
      call dfplps (chh,1,2)
c draw box around plots
      x(1) = smin
      x(2) = smax
      x(3) = smax
      x(4) = smin
      x(5) = smin
      y(1) = tmin
      y(2) = tmin
      y(3) = tmax
      y(4) = tmax
      y(5) = tmin
c draw polyline
      call gpl(5,x,y)
c draw ticks
      dx = 0.
      if (ntx.gt.1) dx = (smax - smin)/real(ntx - 1)
c sty = size of tick mark in y direction
      sty = .4*chh
c draw ticks along x axis
      y(1) = tmin
      y(2) = tmin - sty
      y(3) = tmax
      y(4) = tmax + sty
      do 10 j = 1, ntx
      x(1) = smin + dx*(j - 1)
      x(2) = x(1)
c draw polylines
      call gpl(2,x,y)
      call gpl(2,x,y(3))
   10 continue
c draw ticks along y axis
      dy = 0.
      if (nty.gt.1) dy = (tmax - tmin)/real(nty - 1)
c stx = size of tick marks in x direction
      stx = sty*(rx/ry)
      x(1) = smin
      x(2) = smin - stx
      x(3) = smax
      x(4) = smax + stx
      do 20 j = 1, nty
      y(1) = tmin + dy*(j - 1)
      y(2) = y(1)
c draw polylines
      call gpl(2,x,y)
      call gpl(2,x(3),y)
   20 continue
c label x axes
      ay = tmin - 1.6*chh
      write (lbl,91) xmin
      call trimc(lbl,is,ls)
      ax = smin
c draw text
      call gtx(ax,ay,lbl(is:ls))
c set text alignment to (right,normal)
      call gstxal(3,0)
      write (lbl,91) xmax
      call trimc(lbl,is,ls)
      ax = smax
c draw text
      call gtx(ax,ay,lbl(is:ls))
c label y axes
      ax = smin - 2.*stx
      write (lbl,91) ymin
      call trimc(lbl,is,ls)
      ay = tmin
c draw text
      call gtx(ax,ay,lbl(is:ls))
      write (lbl,91) ymax
      call trimc(lbl,is,ls)
      ay = tmax - chh
c draw text
      call gtx(ax,ay,lbl(is:ls))
c set text alignment to (normal,normal)
      call gstxal(0,0)
c write character string labels
      ax = orgx
      ay = orgy + 1.35*chh
c draw text
      call gtx(ax,ay,label)
      ay = orgy + .25*chh
c draw text
      call gtx(ax,ay,chr)
      return
      end
      subroutine trimc(lbl,is,ls)
c this subroutine identifies of non-blank text.  after execution, the
c substring lbl(is:ls) has leading and trailing blanks removed.
c lbl = character to be trimmed
c is = location of first non-blank character
c ls = location of last non-blank character
      character*(*) lbl
c identify leading blanks
      is = 0
      ls = len(lbl)
   10 is = is + 1
      if ((lbl(is:is).eq.' ').and.(is.lt.ls)) go to 10
c identify trailing blanks
      ls = ls + 1
   20 ls = ls - 1
      if ((lbl(ls:ls).eq.' ').and.(ls.gt.is)) go to 20
      return
      end
