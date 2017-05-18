c plot10 package for tektronix 4014
c viktor k. decyk, ucla
c copyright 1990, regents of the university of california
c update: december 9, 1995
      subroutine movabs (i,j)
c absolute line drawing in screen coordinates, move (or dark vector)
c input: all
c i, j = x, y coordinate of the point to which move is desired
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c igs = ascii 29, left justified
      data igs /486539264/
c set terminal to graph mode, circuitry for dark vector
      call buffpk(igs,1)
c encode coordinates
      call waddrs(i,j)
c dark vector set
      igra = 2
      return
      end
      subroutine drwabs (i,j)
c absolute line drawing in screen coordinates, draw (or bright vector)
c input: all
c i, j = x, y coordinate of the endpoint of line (from current location)
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c ltype = current dash type specification
c izxis = current z-axis mode
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c reset line style to solid, if necessary
      if (ltype.ne.0) then
         ltype = 0
c set z-axis mode and line style
         if (iterm.eq.3) call czaxis(izxis)
c set line index instead of line style on color tektronix
         if (iterm.gt.3) call stlclr(1)
      endif
c encode coordinates
      call waddrs(i,j)
c graph mode set
      igra = 1
      return
      end
      subroutine pntabs (i,j)
c absolute line drawing in screen coordinates, move and draw point
c input: all
c i, j = x, y coordinate where point is displayed
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c ifs = ascii 28, left justified
      data ifs /469762048/
c use point plot mode if supported
      if (iterm.ge.3) then
c set point plot mode, if necessary
         if (igra.ne.3) then
c enter point plot mode
            call buffpk(ifs,1)
c point plot mode set
            igra = 3
         endif
c draw line to itself if point plot mode is not supported
      else
c move cursor
         call movabs(i,j)
      endif
c encode coordinates
      call waddrs(i,j)
      return
      end
      subroutine dshabs(i,j,l)
c absolute dash line drawing in screen coordinate
c input: all
c i, j = x, y coordinate of the line endpoint or move destination
c l = dash type specification
c (-1=move, 0=solid, 1=dot, 2=dash-dot, 3=short-dash, 4=long-dash)
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c ltype = current dash type specification
c izxis = current z-axis mode
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c mssf = screen scale factor (1=normal,4=extended)
      dimension lt(4,4)
      save ns,nl
c lt = software dash type patterns
      data lt /5,5,5,5,14,6,4,6,9,6,9,6,23,7,23,7/
c move requested
      if (l.eq.(-1)) then
c move cursor
         call movabs(i,j)
         ns = 0
         return
      endif
c use current line style if line style is unknown
      if ((l.lt.(-1)).or.(l.gt.7)) l = ltype
c reset line style, if necessary
      if (l.ne.ltype) then
         ltype = l
c set z-axis mode and line style
         if (iterm.eq.3) call czaxis(izxis)
c set line index instead of line style on color tektronix
         if (iterm.gt.3) then
            index = ltype + 1
            if (index.eq.8) index = 1
            call stlclr(index)
         endif
c reset software line style
         if ((iterm.lt.3).and.(ltype.ge.1).and.(ltype.le.4)) then
c ns = index to current pattern segment
c nl = length of current pattern segment
            ns = 0
            nl = mssf*lt(ns+1,ltype)
         endif
      endif
c hardware line style should be used
      if ((ltype.eq.0).or.(ltype.gt.4).or.(iterm.ge.3)) then
c encode coordinates
         call waddrs(i,j)
c graph mode set
         igra = 1
         return
      endif
c software dashed line
      cost = float(i - icx)
      sint = float(j - icy)
      alen = sqrt(cost*cost + sint*sint)
      len = alen + .5
c current pattern segment is longer than line length
      if (nl.ge.len) go to 20
c find starting location and direction cosines
      ix0 = icx
      iy0 = icy
      cost = cost/alen
      sint = sint/alen
c iterate pattern segments
   10 anl = float(nl)
c find end coordinate of next segment
      ix = float(ix0) + anl*cost + .5
      iy = float(iy0) + anl*sint + .5
c dark or bright vector flag
      it = ns - (ns/2)*2
c encode coordinates
      if (it.eq.0) call waddrs(ix,iy)
c move cursor
      if (it.eq.1) call movabs(ix,iy)
c increment pattern segment index
      ns = ns + 1
      if (ns.eq.4) ns = 0
c add length of next pattern segment
      nl = nl + mssf*lt(ns+1,ltype)
c do next segment
      if (nl.lt.len) go to 10
c finish up last segment, which may be incomplete
   20 it = ns - (ns/2)*2
c encode coordinates
      if (it.eq.0) call waddrs(i,j)
c move cursor
      if (it.eq.1) call movabs(i,j)
c adjust length of next pattern segment
      nl = nl - len
c graph mode set
      igra = 1
c if segment complete, reset to next pattern segment
      if (nl.eq.0) then
c increment pattern segment index
         ns = ns + 1
         if (ns.eq.4) ns = 0
c find length of next pattern segment
         nl = mssf*lt(ns+1,ltype)
      endif
      return
      end
      subroutine waddrs(i,j)
c internal subroutine for encoding screen coordinates
c input: all
c i, j = x, y coordinate of the screen coordinate to be encoded
c each (x,y) coordinate is first decomposed into two 5 bit pieces,
c Hi X, Lo X, Hi Y, Lo Y. a tag is then added before transmission,
c 32 for Hi X and Hi Y, 64 for Lo X, and 96 for Lo Y.
c except for 2 special cases, only those five bit pieces which
c changed from the last transmission need to be sent.  The exceptions
c are Lo X, which is always sent, and Hi X, which requires the
c transmission of Lo Y along with Hi X.
c if extended addressing is used, then a 4 bit Extra byte is used, which
c contains 2 extra low order bits of Y, followed by 2 extra bits of X,
c followed by the tag 96.  if the Extra byte is sent, the Lo Y must also
c be sent.
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c ixhp = highest order 5 bits of previous x screen coordinate + 32
c iyhp = highest order 5 bits of previous y screen coordinate + 32
c iylp = lowest order 5 bits of previous y screen coordinate + 96
c iextp = 4 extra bits of previous x, y screen coordinates + 96
c mssf = screen scale factor (1=normal,4=extended)
c decompose coordinates to 5 bit pieces and add tags
      ii = j/mssf
      iyh = ii/32
      iyl = ii - iyh*32 + 96
      iyh = iyh + 32
      iext = 4*(j - mssf*ii)
      ii = i/mssf
      ixh = ii/32
      ixl = ii - ixh*32  + 64
      ixh = ixh + 32
      iext = iext + (i - mssf*ii) + 96
c if not in graphics mode, send everything
      if (igra.eq.0) then
c encode everything
         if (mssf.eq.1) then
            line = ixl + 256*(ixh + 256*(iyl + 256*iyh))
            len = 4
c include Extra
         else
            line = ixh + 256*(iyl + 256*(iext + 256*iyh))
            len = 4
c transmit coordinates
            call buffpk(line,len)
            line = 16777216*ixl
            len = 1
            iextp = iext
         endif
c save previous coordinates for next time
         ixhp = ixh
         iyhp = iyh
         iylp = iyl
         go to 10
      endif
c Hi Y not changed
      if (iyh.eq.iyhp) then
c Hi X not changed
         if (ixh.eq.ixhp) then
c Extra not changed
            if (iext.eq.iextp) then
c Lo Y not changed
               if (iyl.eq.iylp) then
c encode Lo X
                  line = 16777216*ixl
                  len = 1
c Lo Y changed
               else
c encode Lo Y, Lo X
                  line = 65536*(ixl + 256*iyl)
                  len = 2
c save previous coordinate for next time
                  iylp = iyl
               endif
c Extra changed
            else
c encode Extra, Lo Y, Lo X
               line = 256*(ixl + 256*(iyl + 256*iext))
               len = 3
c save previous coordinates for next time
               iylp = iyl
               iextp = iext
            endif
c Hi X changed
         else
c Extra not changed
            if (iext.eq.iextp) then
c encode Lo Y, Hi X, Lo X
               line = 256*(ixl + 256*(ixh + 256*iyl))
               len = 3
c save previous coordinates for next time
               ixhp = ixh
               iylp = iyl
c Extra changed
            else
c encode Extra, Lo Y, Hi X, Lo X
               line = ixl + 256*(ixh + 256*(iyl + 256*iext))
               len = 4
c save previous coordinates for next time
               ixhp = ixh
               iylp = iyl
               iextp = iext
            endif
         endif
c Hi Y changed
      else
c Hi X not changed
         if (ixh.eq.ixhp) then
c Extra not changed
            if (iext.eq.iextp) then
c Lo Y not changed
               if (iyl.eq.iylp) then
c encode Hi Y, Lo X
                  line = 65536*(ixl + 256*iyh)
                  len = 2
c save previous coordinate for next time
                  iyhp = iyh
               else
c encode Hi Y, Lo Y, Lo X
                  line = 256*(ixl + 256*(iyl + 256*iyh))
                  len = 3
c save previous coordinates for next time
                  iyhp = iyh
                  iylp = iyl
               endif
c Extra changed
            else
c encode Hi Y, Extra, Lo Y, Lo X
               line = ixl + 256*(iyl + 256*(iext + 256*iyh))
               len = 4
c save previous coordinates for next time
               iyhp = iyh
               iylp = iyl
               iextp = iext
            endif
c Hi X changed
         else
c encode everything
c Extra not changed
            if (iext.eq.iextp) then
               line = ixl + 256*(ixh + 256*(iyl + 256*iyh))
               len = 4
c include Extra
            else
               line = ixh + 256*(iyl + 256*(iext + 256*iyh))
               len = 4
c transmit coordinates
               call buffpk(line,len)
               line = 16777216*ixl
               len = 1
               iextp = iext
            endif
c save previous coordinates for next time
            ixhp = ixh
            iyhp = iyh
            iylp = iyl
         endif
      endif
c transmit coordinates
   10 call buffpk(line,len)
c save coordinates
      icx = i
      icy = j
      return
      end
      subroutine anmode
c enter alphanumeric mode
      call alfmod
c dump the output buffer
      call tsend
      return
      end
      subroutine erase
c erase screen
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icps = transmission rate in characters per second
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c icx, icy = current x, y cursor location in screen coordinates
c ltype = current dash type specification
c izxis = current z-axis mode
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c mssf = screen scale factor (1=normal,4=extended)
c iescff = ascii 27 and 12, left justified, ascii 22 repeated
c isyn = ascii 22 repeated
      data iescff,isyn /453776918,370546198/
c reset z-axis to normal
      izxis = 0
c reset line style to solid, if necessary
      if (ltype.ne.0) then
         ltype = 0
c set z-axis mode and line style
         if (iterm.eq.3) call czaxis(izxis)
c set line index instead of line style on color tektronix
         if (iterm.gt.3) call stlclr(1)
      endif
c erase screen and select alpha mode
      call buffpk(iescff,4)
c reset cursor to home
      icx = 0
      icy = 767*mssf
      i = icps/4
      do 10 j = 1, i
c send sync characters to pause after erase 
      call buffpk(isyn,4)
   10 continue
c alphanumeric mode set
      igra = 0
      return
      end
      subroutine newlin
c alphanumeric character handling:
c generate a line feed and carriage return
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c icsx, icsy = current character width and height in raster units
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c mssf = screen scale factor (1=normal,4=extended)
c icrlf = ascii 10 and 13, left justified
      data icrlf /218759168/
c move cursor down one line, resets terminal from graph to alpha mode
      call buffpk(icrlf,2)
c reset cursor to left boundary
      icx = 0
      icy = icy - icsy
      if (icy.lt.0) icy = 767*mssf
c alphanumeric mode set
      igra = 0
      return
      end
      subroutine seeloc(i,j)
c locate the position of the graphic beam
c input: none
c i, j = screen x, y coordinate of the beam
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
      i = icx
      j = icy
      return
      end
      subroutine initt(i)
c initialization
c input: all
c i = transmission rate in characters per second
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c xmin = minimum horizontal user coordinate of window
c xrange = horizontal extent of user window rectangle
c ymin = minimum vertical user coordinate of window
c yrange = vertical extent of user window rectangle
c icps = transmission rate in characters per second
c minx = minimum horizontal screen coordinate of window
c lenx = horizontal extent of screen window rectangle
c miny = minimum vertical screen coordinate of window
c leny = vertical extent of screen window rectangle
c icsx, icsy = current character width and height in raster units
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c mssf = screen scale factor (1=normal,4=extended)
      save /tcom/
      icps = i
      minx = 0
      lenx = 1023
      miny = 0
      leny = 779
      xmin = float(minx)
      xrange = float(lenx)
      ymin = float(miny)
      yrange = float(leny)
      icsx = 14
      icsy = 22
      iterm = 1
      mssf = 1
c erase screen
      call erase
      return
      end
      subroutine finitt(i,j)
c termination
c input: all
c i, j = x, y coordinate of the point to which move is desired
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c iesc2 = ascii 27 and 50
      data iesc2 /456261632/
c move cursor
      call movabs(i,j)
c set alphanumeric mode
      call alfmod
c select vt-100 mode, turn on dialogue screen, and make it visible
      if (iterm.ge.4) call selcvt
c send escape 2
      call buffpk(iesc2,2)
c dump the output buffer
      call buffpk(i,0)
      return
      end
      subroutine swindo(mx,lx,my,ly)
c define the screen window
c input: all
c mx = the minimum horizontal screen coordinate
c lx = the horizontal extent of the rectangle
c my = the minimum vertical screen coordinate
c ly = the vertical extent of the rectangle
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c minx = minimum horizontal screen coordinate of window
c lenx = horizontal extent of screen window rectangle
c miny = minimum vertical screen coordinate of window
c leny = vertical extent of screen window rectangle
c mssf = screen scale factor (1=normal,4=extended)
      minx = max0(mx,0)
      lenx = min0(lx,1024*mssf-1)
      miny = max0(my,0)
      leny = min0(ly,780*mssf-1)
      return
      end
      subroutine vwindo(xm,xr,ym,yr)
c define the virtual window
c input: all
c xm = the minimum horizontal user coordinate
c xr = the horizontal extent of the rectangle
c ym = the minimum vertical user coordinate
c yr = the vertical extent of the rectangle
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c xmin = minimum horizontal user coordinate of window
c xrange = horizontal extent of user window rectangle
c ymin = minimum vertical user coordinate of window
c yrange = vertical extent of user window rectangle
      xmin = xm
      xrange = xr
      ymin = ym
      yrange = yr
      return
      end
      subroutine movea(x,y)
c line drawing in user (virtual) units, move (or dark vector)
c input: all
c x, y = x, y coordinate of the point to which move is desired
c convert to screen coordinates
      call caddrs(x,y,i,j,ierr)
c move cursor
      call movabs(i,j)
      return
      end
      subroutine drawa(x,y)
c line drawing in user (virtual) units, draw (or bright vector)
c input: all
c x, y = x, y coordinate of the endpoint of line (from current location)
c convert to screen coordinates
      call caddrs(x,y,i,j,ierr)
c draw line
      if (ierr.eq.0) call drwabs(i,j)
      return
      end
      subroutine pointa(x,y)
c line drawing in user (virtual) units, move and draw point
c input: all
c x, y = x, y coordinate where point is displayed
c convert to screen coordinates
      call caddrs(x,y,i,j,ierr)
c draw point
      if (ierr.eq.0) call pntabs(i,j)
      return
      end
      subroutine dasha(x,y,l)
c dash line drawing in user (virtual) units
c input: all
c x, y = x, y coordinate of the line endpoint or move destination
c l = dash type specification
c (-1=move, 0=solid, 1=dot, 2=dash-dot, 3=short-dash, 4=long-dash)
c convert to screen coordinates
      call caddrs(x,y,i,j,ierr)
c draw dashed line
      if (ierr.eq.0) call dshabs(i,j,l)
      return
      end
      subroutine caddrs(x,y,i,j,ierr)
c this is an internal subroutine which converts from user to screen
c coordinates
c input: x, y
c x, y = x, y user coordinate
c i, j = x, y screen coordinate
c ierr = error return code (0=normal,1=current & previous point clipped)
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c xmin = minimum horizontal user coordinate of window
c xrange = horizontal extent of user window rectangle
c ymin = minimum vertical user coordinate of window
c yrange = vertical extent of user window rectangle
c icx, icy = current x, y cursor location in screen coordinates
c minx = minimum horizontal screen coordinate of window
c lenx = horizontal extent of screen window rectangle
c miny = minimum vertical screen coordinate of window
c leny = vertical extent of screen window rectangle
      ierr = 0
c find x screen coordinate
      i = ((x - xmin)/xrange)*float(lenx) + .5
c clip in x direction
      if (i.lt.0) then
         i = 0
         if ((icx-minx).eq.0) ierr = 1
      elseif (i.gt.lenx) then
         i = lenx
         if ((icx-minx).eq.lenx) ierr = 1
      endif
      i = i + minx
c find y screen coordinate
      j = ((y - ymin)/yrange)*float(leny) + .5
c clip in y direction
      if (j.lt.0) then
         j = 0
         if ((icy-miny).eq.0) ierr = 1
      elseif (j.gt.leny) then
         j = leny
         if ((icy-miny).eq.leny) ierr = 1
      endif
      j = j + miny
      return
      end
      subroutine alfmod
c internal subroutine for entering alphanumeric mode
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c ius = ascii 31, left justified
      data ius /520093696/
c reset terminal from graph to alpha mode
      call buffpk(ius,1)
c alphanumeric mode set
      igra = 0
      return
      end
      subroutine tsend
c dump the output buffer
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icps = transmission rate in characters per second
c omit flush buffer if special code (icps<2) is set
      if (icps.gt.1) call buffpk(i,0)
      return
      end
      subroutine seetw(ixmin,ixmax,iymin,iymax)
c return the current values of the screen window
c input: none
c ixmin = the minimum horizontal screen coordinate
c ixmax = the maximum horizontal screen coordinate
c iymin = the minimum vertical screen coordinate
c iymax = the maximum vertical screen coordinate
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c minx = minimum horizontal screen coordinate of window
c lenx = horizontal extent of screen window rectangle
c miny = minimum vertical screen coordinate of window
c leny = vertical extent of screen window rectangle
      ixmin = minx
      ixmax = minx + lenx
      iymin = miny
      iymax = miny + leny
      return
      end
      subroutine seedw(axmin,axmax,aymin,aymax)
c return the current values of the virtual window limits
c input: none
c axmin = the minimum horizontal user coordinate
c axmax = the maximum horizontal user coordinate
c aymin = the minimum vertical user coordinate
c aymax = the maximum vertical user coordinate
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c xmin = minimum horizontal user coordinate of window
c xrange = horizontal extent of user window rectangle
c ymin = minimum vertical user coordinate of window
c yrange = vertical extent of user window rectangle
      axmin = xmin
      axmax = xmin + xrange
      aymin = ymin
      aymax = ymin + yrange
      return
      end
      subroutine bell
c output an audible tone
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c ibel = ascii 7, left justified
      data ibel /117440512/
c ring bell
      call buffpk(ibel,1)
c graph mode reset
      if (igra.eq.2) igra = 1
      return
      end
      subroutine baksp
c alphanumeric character handling:
c generate a backspace
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c icsx, icsy = current character width and height in raster units
c mssf = screen scale factor (1=normal,4=extended)
c ibs = ascii 8, left justified
      data ibs /134217728/
c backspace
      call buffpk(ibs,1)
      icx = icx - icsx
      if (icx.lt.0) icx = icx + 1024*mssf
      return
      end
      subroutine linef
c alphanumeric character handling:
c generate a line feed
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c icsx, icsy = current character width and height in raster units
c mssf = screen scale factor (1=normal,4=extended)
c ilf = ascii 10, left justified
      data ilf /167772160/
c move cursor down one line
      call buffpk(ilf,1)
      icy = icy - icsy
      if (icy.lt.0) icy = 767*mssf
      return
      end
      subroutine cartn
c alphanumeric character handling:
c generate a carriage return
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c icr = ascii 13, left justified
      data icr /218103808/
c resets terminal from graph to alpha mode
      call buffpk(icr,1)
c reset cursor to home
      icx = 0
c alphanumeric mode set
      igra = 0
      return
      end
      subroutine home
c alphanumeric character handling:
c moves the alphanumeric cursor to the upper left corner of the screen
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c mssf = screen scale factor (1=normal,4=extended)
      i = 0
      j = 767*mssf
c absolute line drawing in screen coordinates, move (or dark vector)
      call movabs(i,j)
      return
      end
      subroutine csize(ihorz,ivert)
c measure the size of a character
c input: none
c ihorz =  horizontal character dimension, including spacing
c ivert = vertical character dimension, including spacing
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icsx, icsy = current character width and height in raster units
      ihorz = icsx
      ivert = icsy
      return
      end
      subroutine chrsiz(ichar)
c change the character size on the 4014/15 terminal
c input: all
c ichar = character size code
c (1=large,2=medium-large,3=medium-small,4=small)
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icsx, icsy = current character width and height in raster units
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c mssf = screen scale factor (1=normal,4=extended)
c iescs = ascii 27 and 55, right justified
      data iescs /6967/
      if ((iterm.eq.1).or.(ichar.lt.1).or.(ichar.gt.4)) return
c large characters
      if (ichar.eq.1) then
         icsx = 14*mssf
         icsy = 22*mssf
c medium-large characters
      elseif (ichar.eq.2) then
         icsx = 13*mssf
         icsy = 21*mssf
c medium-small characters
      elseif (ichar.eq.3) then
         icsx = 9*mssf
         icsy = 13*mssf
c small characters
      elseif (ichar.eq.4) then
         icsx = 8*mssf
         icsy = 12*mssf
      endif
c select character size
      if (iterm.eq.5) then
         call stgtxs(icsx,icsy)
      else
         line = 65536*(iescs + ichar)
         call buffpk(line,2)
      endif
      return
      end
      subroutine czaxis(icode)
c modify the z-axis of the 4014/15 terminal
c icode = z-axis mode
c (0=normal,1=defocussed,2=enabled write-through,3=enabled non-store)
c modes 3 and 4 are not implemented
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c izxis = current z-axis mode
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c iescz = ascii 27 and 96, right justified
      data iescz /7008/
      if ((iterm.lt.3).or.(icode.lt.0).or.(icode.gt.3)) return
      izxis = icode
c replace enabled write-through with defocussed mode
      if (izxis.eq.3) izxis = 2
c left justify command
      line = 65536*(iescz + 8*izxis + ltype)
c select z-axis mode
      call buffpk(line,2)
      return
      end
      subroutine term(itrm,iscal)
c identify the 4014/15 terminal
c input: all
c itrm = terminal type (1=4014,2=4014,3=4014 w/EGM)
c iscal = number of addressable points (1024 or 4096)
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c xmin = minimum horizontal user coordinate of window
c xrange = horizontal extent of user window rectangle
c ymin = minimum vertical user coordinate of window
c yrange = vertical extent of user window rectangle
c minx = minimum horizontal screen coordinate of window
c lenx = horizontal extent of screen window rectangle
c miny = minimum vertical screen coordinate of window
c leny = vertical extent of screen window rectangle
c icsx, icsy = current character width and height in raster units
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c mssf = screen scale factor (1=normal,4=extended)
      if ((itrm.lt.1).or.(itrm.gt.5)) return
      iterm = itrm
      if ((itrm.ge.3).and.(iscal.eq.4096)) then
         mssf = 4
         minx = 0
         lenx = 4095
         miny = 0
         leny = 3119
         xmin = float(minx)
         xrange = float(lenx)
         ymin = float(miny)
         yrange = float(leny)
         icsx = 56
         icsy = 88
      endif
c set default character size to large
      if (iterm.ge.2) call chrsiz(1)
c set default z-axis mode to normal
      if (iterm.eq.3) call czaxis(0)
c color tektronix
      if (iterm.ge.4) then
c select tek mode, turn off dialogue screen, and make it invisible
         call selctk
c set surface color map to standard colors
         call dclrmp
c set line index to white on color tektronix
         call stlclr(1)
      endif
      return
      end
      subroutine seemod(line,izaxis,mode)
c check terminal modes
c input: none
c line = hardware line type in effect
c izaxis = hardware z-axis mode
c mode = software mode
c (0=alphanumeric,1=vector,2=point plot,3=incremental,4=dash)
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c igra = terminal mode (0=alpha,1=graph,2=dark vector,3=point plot)
c ltype = current dash type specification
c izxis = current z-axis mode
      line = ltype
      izaxis = izxis
      mode = igra
      if (ltype.gt.0) mode = 4
      if (mode.eq.3) mode = 2
      return
      end
      subroutine seetrm(ispeed,itrm,isize,maxsr)
c check terminal
c input: none
c ispeed = baud rate in characters per second
c itrm = terminal type
c isize = character size set in chrsiz
c maxsr = screen address range
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icps = transmission rate in characters per second
c icsx, icsy = current character width and height in raster units
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c mssf = screen scale factor (1=normal,4=extended)
      ispeed = icps
      itrm = iterm
      isize = 0
      if ((icsx.eq.(14*mssf)).and.(icsy.eq.(22*mssf))) isize = 1
      if ((icsx.eq.(13*mssf)).and.(icsy.eq.(21*mssf))) isize = 2
      if ((icsx.eq.(9*mssf)).and.(icsy.eq.(13*mssf))) isize = 3
      if ((icsx.eq.(8*mssf)).and.(icsy.eq.(12*mssf))) isize = 4
      maxsr = 1024*mssf
      return
      end
      subroutine aoutst(len,chr)
c outputs an array of characters
c input: all
c len = number of characters to be output
c chr = character string to be output
c len should be < 129
c lw = number of bytes per word
      parameter(lw=4)
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c icsx, icsy = current character width and height in raster units
c iterm = current terminal (1=4010,2=4014,3=4014 w/EGM,4=4105,5=4207)
c mssf = screen scale factor (1=normal,4=extended)
      common /march/ iebc, irvb, longi, mtype
c iebc = (0,1) = (no,yes) input characters are in ebcdic
      character*(*) chr
      dimension lout(32), ieta(256)
      save ieta
c ebcdic/ascii translation with conventions at ucla oac's ibm 3090vf.
c ascii codes for ebcdic 74,79,95,113,139,155 are non-standard
c ebcdic codes 34,53,106,161,192,208,224 are added for ibm compatibility
      data ieta /0,1,2,3,-1,9,-1,127,-1,-1,-1,11,12,13,14,15,16,17,18,19
     1,-1,-1,8,-1,24,25,-1,-1,28,29,30,31,-1,-1,28,-1,-1,10,23,27,-1,-1,
     2-1,-1,-1,5,6,7,-1,-1,22,-1,-1,30,-1,4,-1,-1,-1,-1,20,21,-1,26,32,-
     31,-1,-1,-1,-1,-1,-1,-1,-1,92,46,60,40,43,124,38,-1,-1,-1,-1,-1,-1,
     4-1,-1,-1,33,36,42,41,59,126,45,47,-1,-1,-1,-1,-1,-1,-1,-1,124,44,3
     57,95,62,63,-1,94,-1,-1,-1,-1,-1,-1,-1,96,58,35,64,39,61,34,-1,97,9
     68,99,100,101,102,103,104,105,-1,123,-1,-1,-1,-1,-1,106,107,108,109
     7,110,111,112,113,114,-1,125,-1,-1,-1,-1,-1,126,115,116,117,118,119
     8,120,121,122,-1,-1,-1,91,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
     9,-1,93,-1,-1,123,65,66,67,68,69,70,71,72,73,-1,-1,-1,-1,-1,-1,125,
     a74,75,76,77,78,79,80,81,82,-1,-1,-1,-1,-1,-1,92,-1,83,84,85,86,87,
     b88,89,90,-1,-1,-1,-1,-1,-1,48,49,50,51,52,53,54,55,56,57,-1,-1,-1,
     c-1,-1,-1/
c set alphanumeric mode
      if (iterm.lt.5) call alfmod
c graphic text mode
      if (iterm.eq.5) call gtxmod(len)
      l = (len - 1)/lw + 1
c zero output array
      do 10 i = 1, l
      lout(i) = 0
   10 continue
c convert characters to integer
      do 20 i = 1, len
      j = (i - 1)/lw + 1
c input is in ascii
      if (iebc.eq.0) then
         lout(j) = ichar(chr(i:i)) + 256*lout(j)
c input is in ebcdic
      else
         lout(j) = ieta(ichar(chr(i:i))+1) + 256*lout(j)
      endif
   20 continue
c left shift if necessary
      it1 = 256**(lw*l - len)
      lout(l) = it1*lout(l)
c write output
      call buffpk(lout,len)
      call tsend
c move cursor to end of string
      icx = icx + icsx*len
      if (icx.ge.(1024*mssf)) then
         icx = icx - 1024*mssf
         icy = icy - icsy
         if (icy.lt.0) icy = 767*mssf
      endif
      return
      end
      subroutine ainst(n,chr)
c accept array of characters from terminal
c input: n
c n = number of characters expected
c chr = input character string
c using fortran io, maximum number of characters expected = 80
      character*80 input
      character*(*) chr
   91 format (a80)
c make sure read does not overflow arrays
      l = len(chr)
      m = min0(n,l)
      if (m.gt.80) m = 80
c pad input with blanks
      do 10 i = 1, 80
      input(i:i) = ' '
   10 continue
c pad output with blanks
      do 20 i = 1, l
      chr(i:i) = ' '
   20 continue
c read input
c     call itget(input,n)
      read (5,91,end=30) input
   30 chr(1:m) = input(1:m)
      return
      end
      subroutine tinput(ichr)
c read one character from terminal and put into an ascii decimal integer
c input: none
c ichr = ascii code of input character
c written for the ibm rs/6000
c     call adein(len,ichr)
      character*1 c
c read array of characters from terminal
      call ainst(1,c)
c convert to integer
      ichr = ichar(c)
      return
      end
      subroutine tinstr(len,iarray)
c read input from terminal and put into an ascii decimal array
c input: len
c len = number of characters expected
c iarray = input array into which ascii code of characters is placed
c iarray should be at least len words long
c written for the ibm rs/6000
      dimension iarray(*)
      character*80 input
      n = min0(len,80)
c read array of characters from terminal
      call ainst(n,input)
c copy input into integer array
      do 10 i = 1, n
      iarray(i) = ichar(input(i:i))
   10 continue
      return
      end
      subroutine hdcopy
c generate hardcopy of screen contents
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icps = transmission rate in characters per second
c iscetb = ascii 27 and 23, left justified, ascii 22 repeated
c isyn = ascii 22 repeated
      data iscetb,isyn /454497814,370546198/
c creates make copy signal to hardcopy unit
      call buffpk(iscetb,4)
      i = icps + icps/2
      do 10 j = 1, i
c send sync characters to pause after hardcopy 
      call buffpk(isyn,4)
   10 continue
c dump the output buffer
      call tsend
      return
      end
      subroutine scursr(ichr,i,j)
c use the screen cursor
c input: none
c ichr = a keyboard character, 7 bit ascii right-adjusted
c i, j = screen x, y coordinate of the graphic cursor
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c mssf = screen scale factor (1=normal,4=extended)
c iebc = (0,1) = (no,yes) input characters are in ebcdic
      common /march/ iebc, irvb, longi, mtype
      character*8 input
      dimension ieta(256)
      save ieta
c iscsbu = ascii 27 and 26, left justified
      data iscsub /454688768/
c ebcdic/ascii translation with conventions at ucla oac's ibm 3090vf.
c ascii codes for ebcdic 74,79,95,113,139,155 are non-standard
c ebcdic codes 34,53,106,161,192,208,224 are added for ibm compatibility
      data ieta /0,1,2,3,-1,9,-1,127,-1,-1,-1,11,12,13,14,15,16,17,18,19
     1,-1,-1,8,-1,24,25,-1,-1,28,29,30,31,-1,-1,28,-1,-1,10,23,27,-1,-1,
     2-1,-1,-1,5,6,7,-1,-1,22,-1,-1,30,-1,4,-1,-1,-1,-1,20,21,-1,26,32,-
     31,-1,-1,-1,-1,-1,-1,-1,-1,92,46,60,40,43,124,38,-1,-1,-1,-1,-1,-1,
     4-1,-1,-1,33,36,42,41,59,126,45,47,-1,-1,-1,-1,-1,-1,-1,-1,124,44,3
     57,95,62,63,-1,94,-1,-1,-1,-1,-1,-1,-1,96,58,35,64,39,61,34,-1,97,9
     68,99,100,101,102,103,104,105,-1,123,-1,-1,-1,-1,-1,106,107,108,109
     7,110,111,112,113,114,-1,125,-1,-1,-1,-1,-1,126,115,116,117,118,119
     8,120,121,122,-1,-1,-1,91,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
     9,-1,93,-1,-1,123,65,66,67,68,69,70,71,72,73,-1,-1,-1,-1,-1,-1,125,
     a74,75,76,77,78,79,80,81,82,-1,-1,-1,-1,-1,-1,92,-1,83,84,85,86,87,
     b88,89,90,-1,-1,-1,-1,-1,-1,48,49,50,51,52,53,54,55,56,57,-1,-1,-1,
     c-1,-1,-1/
c put terminal in graphics mode
      call movabs(icx,icy)
c sets gin mode and starts crosshair cursor
      call buffpk(iscsub,2)
c dump the output buffer
      call tsend
c read array of characters from terminal
      call ainst(8,input)
c enter alphanumeric mode, and flush buffer
      call anmode
c translate from ebcdic
      if (iebc.eq.1) then
         do 10 i = 1, 5
         input(i:i) = char(ieta(ichar(input(i:i))+1))
   10    continue
      endif
c find keyboard character
      ichr = ichar(input(1:1))
c calculate x coordinate
      ixh = ichar(input(2:2)) - 32
      ixl = ichar(input(3:3)) - 32
      i = mssf*(32*ixh + ixl)
c calculate y coordinate
      iyh = ichar(input(4:4)) - 32
      iyl = ichar(input(5:5)) - 32
      j = mssf*(32*iyh + iyl)
      return
      end
      subroutine vcursr(ichr,x,y)
c use the virtual cursor
c input: none
c ichr = a keyboard character, 7 bit ascii right-adjusted
c x, y = virtual (user) x, y coordinate of the graphic cursor
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c xmin = minimum horizontal user coordinate of window
c xrange = horizontal extent of user window rectangle
c ymin = minimum vertical user coordinate of window
c yrange = vertical extent of user window rectangle
c minx = minimum horizontal screen coordinate of window
c lenx = horizontal extent of screen window rectangle
c miny = minimum vertical screen coordinate of window
c leny = vertical extent of screen window rectangle
c use the screen cursor
      call scursr(ichr,i,j)
c convert to user coordinates
      x = (float(i - minx)/float(lenx))*xrange + xmin
      y = (float(j - miny)/float(leny))*yrange + ymin
      return
      end
      subroutine buffpk(line,n)
c this program packs characters stored as integers into a buffer of
c characters stored as integers.  characters are assumed to be 7 bits.
c when full, buffer is written to graphics device
c input: all
c line = input characters stored as integers
c n = number of characters to be packed into buffer (0 = flush buffer)
c lw = number of bytes per word
c lmax = maximum number of characters in lout buffer
      parameter(lw=4,lmax=240)
      dimension line(*)
c lout = buffer of particles packed as integers
      dimension lout((lmax-1)/lw+2)
      save len,lout
c len = number of characters currently stored in lout buffer 
      data len /0/
      data lout(1) /0/
c nc = number of characters to be written to lout
      nc = n
c ncr = number of characters remaining to be written
      ncr = n
c flush buffer requested
      if (n.eq.0) then
         if (len.ge.0) go to 40
c buffer already flushed
         return
      endif
c calculate position pointers in character arrays
      ncr = (nc + len) - lmax
      if (ncr.ge.0) nc = lmax - len
c m = last full word written in lout buffer
      m = len/lw
c l = last character location written in last word of lout buffer
      l = len - lw*m
c lr = number of characters spaces available in last word of lout array
      lr = lw - l
c ls = right shift operator for packing input data
      ls = 256**l
c ll = left shift operator for packing input data
      ll = 1
      if (l.gt.0) ll = 256**lr
c initialize loop counters
      nct = nc
      i = 1
c main loop for packing character data into integer array
   10 lt = line(i)/ls
c pack lr characters from current input word into last word in lout
      lout(m+i) = lout(m+i) + lt
c check if done
      if (nct.lt.lr) go to 30
c pack remaining l characters from current input into next word in lout
      lout(m+i+1) = (line(i) - ls*lt)*ll
c check if done
   20 if (nct.le.lw) go to 30
c update loop counters
      i = i + 1
      nct = nct - lw
      go to 10
c update last position written in lout buffer
   30 len = len + nc
c no more characters left to write
      if (ncr.lt.0) return
c flush buffer if full
   40 if (len.gt.0) call tputc(lout,len)
c reset last position written in lout buffer
      len = 0
c force flush buffer
      if (n.eq.0) call tputc(lout,len)
c more characters remaining to be written
      if (ncr.gt.0) then
c update position pointers in character arrays
         nc = ncr
         ncr = ncr - lmax
         if (ncr.ge.0) nc = lmax
c copy last (overflow) word in lout buffer to beginning of buffer
         lout(1) = lout(m+i+1)
c update loop counters
         m = -i
         nct = nct + nc
         go to 20
c no more characters left to be written
      else
         lout(1) = 0
      endif
      return
      end
c addendum to plot10 package for tektronix 4105
c viktor k. decyk, ucla
c copyright 1989, regents of the university of california
      subroutine dclrmp
c set surface color map
c internal subroutine to define standard graphics color map for
c tektronix 4105. the colors defined are:
c 0=black, 1=white, 2=red, 3=green, 4=blue, 5=cyan, 6=magenta, 7=yellow
      dimension icmp(3,8)
      save iestg1, ifour, icmp
c iestg1 = escape TG1 = ascii 27 and 84 and 71 and 49
c ifour = ascii 52, left justified
      data iestg1,ifour /458508081,872415232/
c HLS values for 8 standard colors
      data icmp/0,0,0,0,100,0,120,50,100,240,50,100,0,50,100,300,50,100,
     160,50,100,180,50,100/
      do 20 j = 1, 8
c send set surface color map escape code
      call buffpk(iestg1,4)
      j1 = j - 1
c encode integer parameters in host syntax
      call codei(j1,line,n)
      line = ifour + line/256
      j1 = n + 1
c transmit index
      call buffpk(line,j1)
      do 10 i = 1, 3
c encode integer parameters in host syntax
      call codei(icmp(i,j),line,n)
c transmit HLS coordinates
      call buffpk(line,n)
   10 continue
   20 continue
      return
      end
      subroutine stlclr(l)
c internal subroutine to set line index
c this subroutine sets line color index for tektronix 4105
c colors are: black, white, red, green, blue, cyan, magenta, and yellow
c input: all
c l = line color index
c iescml = escape ML0 = ascii 27 and 77 and 76 and 48
      data iescml /458050608/
      line = iescml + l
c set index
      call buffpk(line,4)
      return
      end
      subroutine sttclr(l)
c internal subroutine to set text index
c this subroutine sets text color index for tektronix 4105
c colors are: black, white, red, green, blue, cyan, magenta, and yellow
c input: all
c l = text color index
c iescmt = escape MT0 = ascii 27 and 77 and 84 and 48
      data iescmt /458052656/
      line = iescmt + l
c set index
      call buffpk(line,4)
      return
      end
      subroutine stgtxs(ichx,ichy)
c set graphtext size
c internal subroutine to set graphtext size for tektronix 4105
c input: all
c ichx, ichy = character width and height in 4010 tekpoint units
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c mssf = screen scale factor (1=normal,4=extended)
c iescmc = escape MC = ascii 27 and 77 and 67, left justified
c i5 = ascii 53, left justified
      data iescmc,i5 /458048256,889192448/
      call buffpk(iescmc,3)
c convert width to terminal space units
      it1 = 3.2*float(ichx/mssf) + .5
c encode integer parameters in host syntax
      call codei(it1,line,n)
c transmit width parameter
      call buffpk(line,n)
c convert height to terminal space units
      it1 = 3.2*float(ichy/mssf) + .5
c encode integer parameters in host syntax
      call codei(it1,line,n)
c transmit height parameter
      call buffpk(line,n)
c transmit 5
      call buffpk(i5,1)
      return
      end
      subroutine gtxmod(len)
c graphic text mode
c internal subroutine to send write graphtext command for tektronix 4105
c input: all
c len = the number of characters to be written
c iescmc = escape LT = ascii 27 and 76 and 84, left justified
      data iesclt /457987072/
c transmit graphtext write command
      call buffpk(iesclt,3)
c encode integer parameters in host syntax
      call codei(len,line,n)
c transmit number of integers to be written
      call buffpk(line,n)
      return
      end
      subroutine selctk
c internal subroutine to select tektronix mode, turn off dialogue screen
c and make it invisible, for tektronix 4105
c iescff = ascii 27 and 12, left justified
      data iescff /453771264/
c iessc0 = escape %!0 = ascii 27 and 37 and 33 and 48
c ieska0 = escape KA0 = ascii 27 and 75 and 65 and 48
c ieslv0 = escape LV0 = ascii 27 and 76 and 86 and 48
      data iessc0,ieska0,ieslv0 /455418160,457916720,457987632/
c erase the screen
      call buffpk(iescff,4)
c select terminal code to tek mode
      call buffpk(iessc0,4)
c disable dialog area
      call buffpk(ieska0,4)
c set dialog area visibility to no
      call buffpk(ieslv0,4)
      return
      end
      subroutine selcvt
c internal subroutine to select vt-100 mode, turn on dialogue screen,
c and make it visible, for tektronix 4105
c iessc0 = escape %!2 = ascii 27 and 37 and 33 and 50
c ieska0 = escape KA1 = ascii 27 and 75 and 65 and 49
c ieslv0 = escape LV1 = ascii 27 and 76 and 86 and 49
      data iessc2,ieska1,ieslv1 /455418162,457916721,457987633/
c set dialog area visibility to yes
      call buffpk(ieslv1,4)
c enable dialog area
      call buffpk(ieska1,4)
c select terminal code to ansi edit mode
      call buffpk(iessc2,4)
      return
      end
      subroutine htab
c generate a tab 
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c icsx, icsy = current character width and height in raster units
c mssf = screen scale factor (1=normal,4=extended)
c ascii 9, left justified
      data iht /150994944/
c space one space to right
      call buffpk(iht,1)
      icx = icx + icsx
      if (icx.ge.(1024*mssf)) icx = icx - 1024*mssf
      return
      end
      subroutine vtab
c generate a vertical tab
      common /tcom/ xmin,xrange,ymin,yrange,
     1icps,icx,icy,igra,iyhp,iylp,ixhp,iextp,minx,lenx,miny,leny,
     2icsx,icsy,ltype,izxis,iterm,mssf
c icx, icy = current x, y cursor location in screen coordinates
c icsx, icsy = current character width and height in raster units
c mssf = screen scale factor (1=normal,4=extended)
c ascii 11, left justified
      data ivt /184549376/
c cause reverse line feed
      call buffpk(ivt,1)
      icy = icy + icsy
      if (icy.gt.(780*mssf)) icy = 0
      return
      end
      subroutine codei(i,line,n)
c encode integer parameters in host syntax
c internal subroutine to pack an integer parameter i, -65535 <i< 65535,
c in a format required by tektronix 4105 escape sequences.
c integers are encoded as a series of one, two, or three characters.
c the 4 lowest order bits plus the sign bit (=16,if positive) plus a tag
c (=32) encoded in the Lo-I character.  If the number is > 15, then the
c next 6 higher order bits plus a tag (=64) are encoded as the first
c Hi-I character.  Finally, if the number is > 1023, the highest order 6
c bits plus a tag (=64) are encoded in the second Hi-I character.  The
c numbers are then written in the following order: second Hi-I (if any),
c first Hi-I (if any), and Lo-I.
c input: i
c i = integer input
c line = left justified encryption
c n = the number of characters used in the encryption
c lw = number of bytes per word
      parameter(lw=4)
      il = i
      line = 0
c set positive sign bit and tag
      is = 48
c assume at least one character output
      n = 1
c negative input, change sign bit and sign of number
      if (il.lt.0) then
         is = 32
         il = -il
      endif
c two character output
      if (il.ge.16) then
         n = 2
c three character output
         if (il.ge.1024) then
            n = 3
            ih = il/1024
            il = il - ih*1024
c second Hi-I character
            line = 64 + ih
         endif
         ih = il/16
         il = il - ih*16
c add first Hi-I character
         line = 256*line + (64 + ih)
      endif
c add Lo-I character
      line = 256*line + (is + il)
c left justify
      line = line*256**(lw - n)
      return
      end
