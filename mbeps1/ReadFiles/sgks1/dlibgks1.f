c 1-2/2d PIC library for diagnostics graphics
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: october 7, 2016
c-----------------------------------------------------------------------
      subroutine GRASP13(part,label,itime,isc,nx,iyp,ixp,idimp,npx,np,  
     1irc)
c for 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c plots background particles in blue and beam particles in red
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 4
c npx = number of background particles
c np = number of particles
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      parameter(nxbs=65)
      dimension part(idimp,np)
      character*(*) label
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
c x,y = scratch arrays for plotting
      dimension x(nxbs), y(nxbs)
      character*26 chr
      character*2 lblsp(4)
      character*10 chrs(2)
      save lblsp
   91 format (1x,a2,' VERSUS ',a2,', T = ',i7)
      data lblsp /' X','VX','VY','VZ'/
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
c csize = vertical size of characters
      data zero,csize /0.,0.034/
c sample labels
      data chrs /'BACKGROUND','   BEAM   '/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c find y scale for plot
      if (iyp.le.1) then
         if (iyp.lt.1) iyp = 1
         ymin = zero
         ymax = real(nx)
      elseif (iyp.gt.1) then
         if (iyp.gt.idimp) iyp = idimp
         is = isc
         if (abs(is).gt.116) then
            ymax = abs(part(iyp,1))
            do 10 j = 1, np
            ymax = amax1(ymax,abs(part(iyp,j)))
   10       continue
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
         endif
         ymax = dv**is
         ymin = -ymax
      endif
c find x scale for plot
      if (ixp.le.1) then
         if (ixp.lt.1) ixp = 1
         xmin = zero
         xmax = real(nx)
      elseif (ixp.gt.1) then
         if (ixp.gt.idimp) ixp = idimp
         is = isc
         if (abs(is).gt.116) then
            xmax = abs(part(ixp,1))
            do 20 j = 1, np
            xmax = amax1(xmax,abs(part(ixp,j)))
   20       continue
            if (xmax.eq.0.) xmax = 1.0e-35
            is = alog(xmax)*algdvi
            if (xmax.ge.1.) is = is + 1
            if (xmax.le.dv**(is-1)) is = is - 1
         endif
         xmax = dv**is
         xmin = -xmax
      endif
c find location for plot
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
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c write labels
      write (chr,91) lblsp(iyp), lblsp(ixp), itime
c select point as marker
      mks = 1
c draw grid and labels, call identity transformation
      call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample markers with dots of minimum visible size
      dd = (smax - smin)*real(npl)*2.0e-3
      ngs = 2
      if ((npx.eq.0).or.(np.eq.npx)) ngs = 1
      call dsdmpln(orx,smn,tmn,tmx,dd,ngs,1,1,mks,chrs,chh)
c define transformation number 2
      nrt = 1
c set window
      call gswn(nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp(nrt,smn,smx,tmn,tmx)
c select normalization transformation
      call gselnt(nrt)
c set marker size scale factor, 1.0 = nominal
      call gsmksc(1.0)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
      call gsmk(mks)
c set clipping indicator, 1 = on
      call gsclip(1)
c plot particles
      do 50 k = 1, 2
c determine how many particles to plot
      if (k.eq.1) then
         nd = npx
      elseif (k.eq.2) then
         nd = np - npx
      endif
      if (nd.eq.0) go to 50
      koff = npx*(k - 1)
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
      call gspmci(icol)
c parameters for plots
      nxs = nxbs - 1
c blocking parameters for plots
      nxb = (nd - 1)/nxs + 1
      npts = nxs
c length of last block
      nptl = nd - nxs*(nxb - 1)
c plot polymarkers
      npt = npts
c loop over number of blocks
      do 40 j = 1, nxb
      js = nxs*(j - 1) + koff
      if (j.eq.nxb) npt = nptl
c calculate x,y axes for block
      do 30 i = 1, npt
      x(i) = part(ixp,i+js)
      y(i) = part(iyp,i+js)
   30 continue
c treat dots by drawing a line to itself
c     call spdots(x,y,npt,icol,nxbs)
c treat dots by drawing a line to itself with non-zero width in x
      dd = (xmax - xmin)*real(npl)*2.0e-3
      call spddots(x,y,dd,npt,icol,nxbs)
   40 continue
   50 continue
c update workstation, perform
      call guwk(idwk,1)
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
      subroutine PGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp,  
     1nppmx,mx1,irc)
c for 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c kpic = number of particles per tile
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      parameter(nxbs=65)
      dimension ppart(idimp,nppmx,mx1), kpic(mx1)
      character*(*) label
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
c x,y = scratch arrays for plotting
      dimension x(nxbs), y(nxbs)
      character*26 chr
      character*2 lblsp(4)
      character*10 chrs(2)
      save lblsp
   91 format (1x,a2,' VERSUS ',a2,', T = ',i7)
      data lblsp /' X','VX','VY','VZ'/
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
c csize = vertical size of characters
      data zero,csize /0.,0.034/
c sample labels
      data chrs /'BACKGROUND','   BEAM   '/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c find y scale for plot
      if (iyp.le.1) then
         if (iyp.lt.1) iyp = 1
         ymin = zero
         ymax = real(nx)
      elseif (iyp.gt.1) then
         if (iyp.gt.idimp) iyp = idimp
         is = isc
         if (abs(is).gt.116) then
            ymax = abs(ppart(iyp,1,1))
!$OMP PARALLEL DO PRIVATE(j,k,npp,symax)
!$OMP& REDUCTION(MAX:ymax)
            do 20 k = 1, mx1
            npp = kpic(k)
            symax = abs(ppart(iyp,1,k))
            do 10 j = 1, npp
            symax = amax1(symax,abs(ppart(iyp,j,k)))
   10       continue
            ymax = amax1(ymax,symax)
   20       continue
!$OMP END PARALLEL DO
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
         endif
         ymax = dv**is
         ymin = -ymax
      endif
c find x scale for plot
      if (ixp.le.1) then
         if (ixp.lt.1) ixp = 1
         xmin = zero
         xmax = real(nx)
      elseif (ixp.gt.1) then
         if (ixp.gt.idimp) ixp = idimp
         is = isc
         if (abs(is).gt.116) then
            xmax = abs(ppart(ixp,1,1))
!$OMP PARALLEL DO PRIVATE(j,k,npp,sxmax)
!$OMP& REDUCTION(MAX:xmax)
            do 40 k = 1, mx1
            npp = kpic(k)
            sxmax = abs(ppart(ixp,1,k))
            do 30 j = 1, npp
            sxmax = amax1(sxmax,abs(ppart(ixp,j,k)))
   30       continue
            xmax = amax1(xmax,sxmax)
   40       continue
!$OMP END PARALLEL DO
            if (xmax.eq.0.) xmax = 1.0e-35
            is = alog(xmax)*algdvi
            if (xmax.ge.1.) is = is + 1
            if (xmax.le.dv**(is-1)) is = is - 1
         endif
         xmax = dv**is
         xmin = -xmax
      endif
c find location for plot
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
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c write labels
      write (chr,91) lblsp(iyp), lblsp(ixp), itime
c select point as marker
      mks = 1
c draw grid and labels, call identity transformation
      call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample markers with dots of minimum visible size
      dd = (smax - smin)*real(npl)*2.0e-3
      ngs = 1
      call dsdmpln(orx,smn,tmn,tmx,dd,ngs,1,1,mks,chrs,chh)
c define transformation number 2
      nrt = 1
c set window
      call gswn(nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp(nrt,smn,smx,tmn,tmx)
c select normalization transformation
      call gselnt(nrt)
c set marker size scale factor, 1.0 = nominal
      call gsmksc(1.0)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
      call gsmk(mks)
c set clipping indicator, 1 = on
      call gsclip(1)
      icol = 1
      icol = kprime(icol+1)
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
      call gspmci(icol)
c plot particles
      do 70 k = 1, mx1
c determine how many particles to plot
      nd = kpic(k)
c parameters for plots
      nxs = nxbs - 1
c blocking parameters for plots
      nxb = (nd - 1)/nxs + 1
      npts = nxs
c length of last block
      nptl = nd - nxs*(nxb - 1)
c plot polymarkers
      npt = npts
c loop over number of blocks
      do 60 j = 1, nxb
      js = nxs*(j - 1)
      if (j.eq.nxb) npt = nptl
c calculate x,y axes for block
      do 50 i = 1, npt
      x(i) = ppart(ixp,i+js,k)
      y(i) = ppart(iyp,i+js,k)
   50 continue
c treat dots by drawing a line to itself
c     call spdots(x,y,npt,icol,nxbs)
c treat dots by drawing a line to itself with non-zero width in x
      dd = (xmax - xmin)*real(npl)*2.0e-3
      call spddots(x,y,dd,npt,icol,nxbs)
   60 continue
   70 continue
c update workstation, perform
      call guwk(idwk,1)
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
      subroutine dsdmpln(orgx,smin,tmin,tmax,dd,ngs,nlts,nmks,mks,chrs,c
     1hh)
c this subroutine displays line or marker samples, with short character
c labels placed underneath
c dots are displayed with with non-zero width in x
c orgx = x origin of window
c smin = minimum x value of plotting window
c tmin/tmax = range of y values of plotting window
c dd = smallest visible size, such as (xmax - xmin)*4.0e-3
c ngs = number of subarrays being plotted
c nlts = number of line types available
c nmks = number of markers available
c mks = flag to determine whether lines or markers are used
c chrs = array of ngs short character labels for line or marker samples,
c each should be less than or equal to 10 characters in length
c chh = character height
      character*(*) chrs(ngs)
c rx, ry = ndc coordinates of upper-right corner of workstation window
c ncols = number of foreground colors available for line plotting
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c x,y = scratch variables for plotting
      dimension x(4), y(2)
c omit samples if there is only one curve
      if (ngs.le.1) return
      mkr = abs(mks)
c set marker size scale factor, 1.0 = nominal
      if (mkr.ne.0) call gsmksc(1.0)
c draw line or marker samples
      stx = .4*chh*(rx/ry)
      x(1) = orgx + 4.*stx
      x(2) = smin - 4.*stx
      dx = (x(2) - x(1))/3.
      x(3) = x(1) + dx
      x(4) = x(2) - dx
      ax = orgx + 2.*stx
      dy = (tmax - tmin - 4.*chh)/real(ngs + 1)
      ay = tmax - chh
c draw samples, first cycle through colors, then line or marker types
      do 30 k = 1, ngs
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
      y(1) = ay - dy*real(k)
      y(2) = y(1)
      if ((mks.eq.0).or.((mks.lt.0).and.(k.eq.1))) then
c set polyline color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gsplci(icol)
         ltype = 1 + (k - 1)/ncols - nlts*((k - 1)/(nlts*ncols))
c set linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
         call gsln(ltype)
c draw polyline
         call gpl(2,x,y)
      else
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
         mtype = mkr + (k - 1)/ncols - nmks*((mkr - 1 + (k - 1)/ncols)/n
     1mks)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mtype)
c dots
         if (mtype.eq.1) then
c treat dots by drawing a line to itself
c           call spdots(x(3),y,2,icol,3)
c treat dots by drawing a line to itself with non-zero width in x
            call spddots(x(3),y,dd,2,icol,3)
         else
c draw polymarker
            call gpm(2,x(3),y)
         endif
      endif
c draw label underneath sample line or marker
      y(1) = y(1) - 2.*chh
c draw text
      call gtx(ax,y(1),chrs(k))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine spddots(x,y,dd,npt,icol,nxbs)
c this subroutine draws dot markers by drawing a line to itself
c with non-zero width in x
c x, y = arrays to be plotted
c dd = smallest visible size, such as (xmax - xmin)*4.0e-3
c npt = number of points to be plotted
c icol = color index
c nxbs = dimension of x, y arrays
      dimension x(nxbs), y(nxbs)
c xs, ys = scratch arrays for plotting
      dimension xs(2), ys(2)
c set polyline color index
      call gsplci(icol)
      do 10 j = 1, npt
      xs(1) = x(j)
      ys(1) = y(j)
      xs(2) = xs(1) + dd
      ys(2) = ys(1)
c draw polyline
      call gpl(2,xs,ys)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp, 
     1nppmx,mx1,ltag,irc)
c for 1d or 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c with beam  particles marked for color
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c ppart(ltag,n,m) = particle id of tagged particle n in tile m
c kpic = number of particles per tile
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 3
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c ltag = location in particle array of tags
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      parameter(nxbs=1)
      dimension ppart(idimp,nppmx,mx1), kpic(mx1)
      character*(*) label
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
c x,y = scratch arrays for plotting
      dimension x(nxbs), y(nxbs)
      character*26 chr
      character*2 lblsp(4)
      character*10 chrs(2)
      save lblsp
   91 format (1x,a2,' VERSUS ',a2,', T = ',i7)
      data lblsp /' X','VX','VY','VZ'/
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
c csize = vertical size of characters
      data zero,csize /0.,0.034/
c sample labels
      data chrs /'BACKGROUND','   BEAM   '/
      irc = 1
c exit if ltag is in error
      if (ltag.gt.idimp) return
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c find y scale for plot
      if (iyp.le.1) then
         if (iyp.lt.1) iyp = 1
         ymin = zero
         ymax = real(nx)
      elseif (iyp.gt.1) then
         if (iyp.gt.idimp) iyp = idimp
         is = isc
         if (abs(is).gt.116) then
            ymax = abs(ppart(iyp,1,1))
!$OMP PARALLEL DO PRIVATE(j,k,npp,symax)
!$OMP& REDUCTION(MAX:ymax)
            do 20 k = 1, mx1
            npp = kpic(k)
            symax = abs(ppart(iyp,1,k))
            do 10 j = 1, npp
            symax = amax1(symax,abs(ppart(iyp,j,k)))
   10       continue
            ymax = amax1(ymax,symax)
   20       continue
!$OMP END PARALLEL DO
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
         endif
         ymax = dv**is
         ymin = -ymax
      endif
c find x scale for plot
      if (ixp.le.1) then
         if (ixp.lt.1) ixp = 1
         xmin = zero
         xmax = real(nx)
      elseif (ixp.gt.1) then
         if (ixp.gt.idimp) ixp = idimp
         is = isc
         if (abs(is).gt.116) then
            xmax = abs(ppart(ixp,1,1))
!$OMP PARALLEL DO PRIVATE(j,k,npp,sxmax)
!$OMP& REDUCTION(MAX:xmax)
            do 40 k = 1, mx1
            npp = kpic(k)
            sxmax = abs(ppart(ixp,1,k))
            do 30 j = 1, npp
            sxmax = amax1(sxmax,abs(ppart(ixp,j,k)))
   30       continue
            xmax = amax1(xmax,sxmax)
   40       continue
!$OMP END PARALLEL DO
            if (xmax.eq.0.) xmax = 1.0e-35
            is = alog(xmax)*algdvi
            if (xmax.ge.1.) is = is + 1
            if (xmax.le.dv**(is-1)) is = is - 1
         endif
         xmax = dv**is
         xmin = -xmax
      endif
c find location for plot
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
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c write labels
      write (chr,91) lblsp(iyp), lblsp(ixp), itime
c select point as marker
      mks = 1
c draw grid and labels, call identity transformation
      call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample markers with dots of minimum visible size
      dd = (smax - smin)*real(npl)*2.0e-3
      ngs = 1
      call dsdmpln(orx,smn,tmn,tmx,dd,ngs,1,1,mks,chrs,chh)
c define transformation number 2
      nrt = 1
c set window
      call gswn(nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp(nrt,smn,smx,tmn,tmx)
c select normalization transformation
      call gselnt(nrt)
c set marker size scale factor, 1.0 = nominal
      call gsmksc(1.0)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
      call gsmk(mks)
c set clipping indicator, 1 = on
      call gsclip(1)
      icol = 1
      icol = kprime(icol+1)
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
      call gspmci(icol)
c plot particles
      do 70 k = 1, mx1
c determine how many particles to plot
      nd = kpic(k)
c plot polymarkers
      npt = 1
c loop over number of blocks
      do 60 j = 1, nd
c calculate x,y axes for block
      do 50 i = 1, npt
      x(i) = ppart(ixp,j,k)
      y(i) = ppart(iyp,j,k)
   50 continue
      if (ppart(ltag,j,k).lt.0.0) then
         icol = kprime(4)
      else
         icol = kprime(3)
      endif
c treat dots by drawing a line to itself
c     call spdots(x,y,npt,icol,nxbs)
c treat dots by drawing a line to itself with non-zero width in x
      dd = (xmax - xmin)*real(npl)*2.0e-3
      call spddots(x,y,dd,npt,icol,nxbs)
   60 continue
   70 continue
c update workstation, perform
      call guwk(idwk,1)
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
