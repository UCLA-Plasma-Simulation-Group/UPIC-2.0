c general 1d gks graphics library
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: february 26, 2013
      subroutine GROPEN
c this subroutine opens gks and activates standard workstation
c colors and maximum size of display surface are also set
c if a string (keyboard) device is available, it is initialized
c if a locator (mouse) device is available, it is initialized
      call GROPEN0(0)
      end
      subroutine GROPEN0(iecho)
c this subroutine opens gks and activates workstation
c colors and maximum size of display surface are also set
c if a string (keyboard) device is available, it is initialized
c if a locator (mouse) device is available, it is initialized
c iecho = (0,1) = echo area is in (lower right, lower left)
c idwk = workstation identifier
c ncols = number of foreground colors available for line plotting
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c idstr = string device number, 0 if no string device available
c idloc = locator device number, 0 if no locator device available
c nclsp = number of foreground colors supported on device
c ifrg = index of foreground color
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c ndi = increment between frames
      common /pextra/ ndi
c scratch strings and arrays needed for initializing input devices
      character*80 datar
      character*8 str
      dimension earea(4)
      save /plotcm/, /pextra/
c nerrfl = error file unit number, 6 for terminal
c meml = storage limit for one segment, in bytes
      data nerrfl,meml /6,0/
      data zero /0./
c set default plot values
      iplot = 0
      nplot = 1
      iclr = 0
      iupd = 0
      ndi = 1
c open gks
      call gopks(nerrfl,meml)
c get workstation parameters
c idcon = connection identifier, iwtype = workstation type
      call gitwks(idcon,iwtype)
      idwk = 1
c open workstation
      call gopwk(idwk,idcon,iwtype)
c activate workstation
      call gacwk(idwk)
c ncoli = number of colors available
c iscol = color availability indicator, 0 = monochrmome, 1 = color
c npci = number of predefined color indices
c inquire color facilities
      call gqcf(iwtype,ierr,ncoli,iscol,npci)
c if metafile, set maximum number of colors to 256
      if (ierr.gt.0) then
         ncoli = 256
         iscol = 1
      endif
c set number of foreground colors, 0 < ncols < 8
      if (iscol.eq.1) then
c special case for IBM graPHIGS gks and DEC gks
c        ncoli = npci
         nclsp = ncoli - 1
         ncols = min0(ncoli,8) - 1
      else
         nclsp = 1
         ncols = 1
      endif
c set colors
      call grcols
c idcun = device coordinate units, 0 = meters
c dcx, dcy = display width, height in device coordinate units
c isx, isy = display width, height, in raster units
c inquire maximum display surface size
      call gqdsp(iwtype,ierr,idcun,dcx,dcy,isx,isy)
c if metafile, set square display
      if (ierr.gt.0) then
         dcx = 1.0
         dcy = 1.0
         isx = 1
         isy = 1
      endif
      if (isx.gt.isy) then
         rx = 1.0
         ry = dcy/dcx
      else
         rx = dcx/dcy
         ry = 1.0
      endif
c set workstation window
      call gswkwn(idwk,zero,rx,zero,ry)
c set workstation viewport
      call gswkvp(idwk,zero,dcx,zero,dcy)
      idstr = 1
c inquire string device state
      call gqsts(idwk,idstr,1,ierr,mode,iesw,lstr,str,ipet,earea,lenb,ip
     1os,ldr,datar)
      if (ierr.eq.0) then
c echo area is in lower right or lower left hand corner
         ecx = .9*dcx
         if (iecho.eq.1) ecx = .05*dcx
         ecy = .05*dcy
c initialize string device
         call ginst(idwk,idstr,lstr,str,ipet,ecx,dcx,zero,ecy,lenb,ipos,
     1ldr,datar)
      else
         idstr = 0
      endif
      idloc = 1
c inquire locator device state
      call gqlcs(idwk,idloc,0,1,ierr,mode,iesw,nrt,px,py,ipet,earea,ldr,
     1datar)
      if (ierr.eq.0) then
c initialize locator device
         call ginlc(idwk,idloc,nrt,px,py,ipet,zero,dcx,zero,dcy,ldr,data
     1r)
      else
         idloc = 0
      endif
      return
      end
      subroutine grcols
c if possible, this subroutine sets the color indices as follows:
c 0 = background, 1 = foreground, 2 = blue, 3 = red
c 4 = yellow, 5 = cyan, 6 = magenta, 7 = green
c the rgb values of these colors are stored in the arrays
c reds, greens, and blues, respectively
c idwk = workstation identifier
c ncols = number of foreground colors available for line plotting
c ifrg = index of foreground color
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      dimension reds(8), greens(8), blues(8)
c background = black, foreground = white
c     data reds /0.,1.,0.,1.,1.,0.,1.,0./
c     data greens /0.,1.,0.,0.,1.,1.,0.,1./
c     data blues /0.,1.,1.,0.,0.,1.,1.,0./
c 2=cyan, 3=green, 4=yellow, 5=blue, 6=red, 7=magenta
      data reds /0.,1.,0.,0.,1.,1.,0.,1./
      data greens /0.,1.,1.,1.,1.,0.,0.,0./
      data blues /0.,1.,1.,0.,0.,1.,1.,0./
c background = white, foreground = black
c     data reds /1.,0.,0.,1.,1.,0.,1.,0./
c     data greens /1.,0.,0.,0.,1.,1.,0.,1./
c     data blues /1.,0.,1.,0.,0.,1.,1.,0./
c     data reds /1.,0.,0.,1.,0.,0.,1.,1./
c     data greens /1.,0.,0.,0.,1.,1.,0.,1./
c     data blues /1.,0.,1.,0.,0.,1.,1.,0./
c set background color index
      kprime(1) = 0
c set color representation
      call gscr(idwk,0,reds(1),greens(1),blues(1))
c set foreground color index
      ifrg = 1
c set color representation
      call gscr(idwk,1,reds(2),greens(2),blues(2))
      do 10 j = 2, 8
      kprime(j) = ifrg
   10 continue
c define allowed color indices
      do 20 j = 2, ncols
      icol = j + 1
c set color representation
      call gscr(idwk,j,reds(icol),greens(icol),blues(icol))
c store index
      kprime(icol) = j
   20 continue
      return
      end
      subroutine grspal(nbit,pal,npal,ipmx,lpald)
c this subroutine sets palette of color indices for nbit color.
c if npal = 0, a default palette is set, if possible, as follows:
c for nbit = 1, 0 = background, 1 = foreground
c for nbit = 2, 0 = background, 1 = cyan, 2 = red, 3 = foreground
c for nbit = 3, 0 = background, 1 = blue, 2 = green, 3 = cyan, 4 = red,
c 5 = magenta, 6 = yellow, 7 = foreground
c for nbit = 4, 0 = background, 1 = dark grey, 2 = blue, 3 = light blue,
c 4 = green, 5 = light green, 6 = cyan, 7 = light cyan, 8 = red, 
c 9 = light red, 10 = magenta, 11 = light magenta, 12 = yellow,
c 13 = light yellow, 14 = white, 15 = intense white
c the rgb values of these colors are stored in the arrays
c reds, greens, and blues, respectively
c if npal > 0, then a specified palette array pal is used, if possible.
c if enough colors are not available, then a default palette with enough
c colors is used, with an integer pixel lookup table to interpolate
c from given palette pal to 7 bit color.  scheme used is as follows:
c for a given palette entry i, ipal(i) = 4*ir + 2*ig + ib, where
c ir = 64*irh+8*irm+irl, ig = 64*igh+8*igm+igl, ib = 64*ibh+8*ibm+ibl
c irh = 0, when the red value for the ith entry is < .5*max(red)
c igh = 0, when the green value for the ith entry is < .5*max(green)
c ibh = 0, when the blue value for the ith entry is < .5*max(blue)
c otherwise, irh, igh, and ibh are 1.
c irm = 0, when the red value - .5*irh*max(red) < .25*max(red)
c igm = 0, when the green value - .5*igh*max(green) < .25*max(green)
c irm = 0, when the blue value - .5*ibh*max(blue) < .25*max(blue)
c otherwise, irm, igm, and ibm are 1.
c irl = 0, when the red value - (.5*irh+.25*irm)*max(red)
c                                               < .125*max(red)
c igl = 0, when the green value - (.5*igh+.25*igm)*max(green)
c                                               < .125*max(green)
c irl = 0, when the blue value - (.5*ibh+.25*ibm)*max(blue)
c                                               < .125*max(blue)
c otherwise, irl, igl, and ibl are 1.
c ipal is then divided by 512/ntc.
c nbit = the number of colors, pixel depth
c pal = palette array, with rgb values in successive bytes
c npal = (0,n) = (default,n) palette entries
c ipmx = maximum color value in the palette
c lpald = size of palette array
c npald = number of palette entries
      parameter(npald=256)
c lxm, lym = maximum number of pixels in x, y
      parameter(lxm=720,lym=540)
      character*1 pal(lpald)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img8 = integer image array
      dimension img8(lxm*lym)
c idwk = workstation identifier
c nclsp = number of foreground colors supported on device
c ifrg = index of foreground color
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      dimension redp(8), greenp(8), bluep(8)
      dimension reds(254), greens(254), blues(254)
      dimension redx(256), greenx(256), bluex(256)
      save /movicm/
c prime colors, background = black, foreground = white
c     data redp /0.,1.,0.,1.,1.,0.,1.,0./
c     data greenp /0.,1.,0.,0.,1.,1.,0.,1./
c     data bluep /0.,1.,1.,0.,0.,1.,1.,0./
c 2=cyan, 3=green, 4=yellow, 5=blue, 6=red, 7=magenta
      data redp /0.,1.,0.,0.,1.,1.,0.,1./
      data greenp /0.,1.,1.,1.,1.,0.,0.,0./
      data bluep /0.,1.,1.,0.,0.,1.,1.,0./
c prime colors, background = black, foreground = white
c     data redp /1.,0.,0.,1.,1.,0.,1.,0./
c     data greenp /1.,0.,0.,0.,1.,1.,0.,1./
c     data bluep /1.,0.,1.,0.,0.,1.,1.,0./
c     data redp /1.,0.,0.,1.,0.,0.,1.,1./
c     data greenp /1.,0.,0.,0.,1.,1.,0.,1./
c     data bluep /1.,0.,1.,0.,0.,1.,1.,0./
      data reds /0.,1.,0.,0.,1.,1.,0.,0.,0.,0.,1.,1.,1.,1.,0.,.33,0.,.33
     1,0.,.33,0.,.33,.67,1.,.67,1.,.67,1.,.67,1.,0.,0.,.33,.33,0.,0.,.33
     2,.33,0.,0.,.33,.33,0.,0.,.33,.33,.67,.67,1.,1.,.67,.67,1.,1.,.67,.
     367,1.,1.,.67,.67,1.,1.,0.,0.,0.,0.,.33,.33,.33,.33,0.,0.,0.,0.,.33
     4,.33,.33,.33,0.,0.,0.,0.,.33,.33,.33,.33,0.,0.,0.,0.,.33,.33,.33,.
     533,.67,.67,.67,.67,1.,1.,1.,1.,.67,.67,.67,.67,1.,1.,1.,1.,.67,.67
     6,.67,.67,1.,1.,1.,1.,.67,.67,.67,.67,1.,1.,1.,1.,0.,.14,0.,.14,0.,
     7.14,0.,.14,.29,.43,.29,.43,.29,.43,.29,.43,0.,.14,0.,.14,0.,.14,0.
     8,.14,.29,.43,.29,.43,.29,.43,.29,.43,0.,.14,0.,.14,0.,.14,0.,.14,.
     929,.43,.29,.43,.29,.43,.29,.43,0.,.14,0.,.14,0.,.14,0.,.14,.29,.43
     a,.29,.43,.29,.43,.29,.43,.57,.71,.57,.71,.57,.71,.57,.71,.86,1.,.8
     b6,1.,.86,1.,.86,1.,.57,.71,.57,.71,.57,.71,.57,.71,.86,1.,.86,1.,.
     c86,1.,.86,1.,.57,.71,.57,.71,.57,.71,.57,.71,.86,1.,.86,1.,.86,1.,
     d.86,1.,.57,.71,.57,.71,.57,.71,.57,.71,.86,1.,.86,1.,.86,1.,.86,1.
     e/
      data greens /0.,1.,0.,1.,0.,1.,0.,0.,1.,1.,0.,0.,1.,1.,0.,.33,0.,.
     133,.67,1.,.67,1.,0.,.33,0.,.33,.67,1.,.67,1.,0.,.33,0.,.33,0.,.33,
     20.,.33,.67,1.,.67,1.,.67,1.,.67,1.,0.,.33,0.,.33,0.,.33,0.,.33,.67
     3,1.,.67,1.,.67,1.,.67,1.,0.,0.,.33,.33,0.,0.,.33,.33,0.,0.,.33,.33
     4,0.,0.,.33,.33,.67,.67,1.,1.,.67,.67,1.,1.,.67,.67,1.,1.,.67,.67,1
     5.,1.,0.,0.,.33,.33,0.,0.,.33,.33,0.,0.,.33,.33,0.,0.,.33,.33,.67,.
     667,1.,1.,.67,.67,1.,1.,.67,.67,1.,1.,.67,.67,1.,1.,0.,.14,0.,.14,.
     729,.43,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.,.14,.29,.43
     8,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,.57,.71,.57,.71,.86,1.,.86,
     91.,.57,.71,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1.,.86,1.,.57
     a,.71,.57,.71,.86,1.,.86,1.,0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.
     b,.14,.29,.43,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.,.14,.
     c29,.43,.29,.43,.57,.71,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1
     d.,.86,1.,.57,.71,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1.,.86,
     e1./
      data blues /0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,.33,.67,1
     1.,0.,.33,.67,1.,0.,.33,.67,1.,0.,.33,.67,1.,0.,.33,0.,.33,.67,1.,.
     267,1.,0.,.33,0.,.33,.67,1.,.67,1.,0.,.33,0.,.33,.67,1.,.67,1.,0.,.
     333,0.,.33,.67,1.,.67,1.,0.,.33,0.,.33,0.,.33,0.,.33,.67,1.,.67,1.,
     4.67,1.,.67,1.,0.,.33,0.,.33,0.,.33,0.,.33,.67,1.,.67,1.,.67,1.,.67
     5,1.,0.,.33,0.,.33,0.,.33,0.,.33,.67,1.,.67,1.,.67,1.,.67,1.,0.,.33
     6,0.,.33,0.,.33,0.,.33,.67,1.,.67,1.,.67,1.,.67,1.,0.,.14,.29,.43,0
     7.,.14,.29,.43,0.,.14,.29,.43,0.,.14,.29,.43,.57,.71,.86,1.,.57,.71
     8,.86,1.,.57,.71,.86,1.,.57,.71,.86,1.,0.,.14,.29,.43,0.,.14,.29,.4
     93,0.,.14,.29,.43,0.,.14,.29,.43,.57,.71,.86,1.,.57,.71,.86,1.,.57,
     a.71,.86,1.,.57,.71,.86,1.,0.,.14,.29,.43,0.,.14,.29,.43,0.,.14,.29
     b,.43,0.,.14,.29,.43,.57,.71,.86,1.,.57,.71,.86,1.,.57,.71,.86,1.,.
     c57,.71,.86,1.,0.,.14,.29,.43,0.,.14,.29,.43,0.,.14,.29,.43,0.,.14,
     d.29,.43,.57,.71,.86,1.,.57,.71,.86,1.,.57,.71,.86,1.,.57,.71,.86,1
     e./
      data redx /0.,0.,.14,.14,0.,0.,.14,.14,0.,0.,.14,.14,0.,0.,.14,.14
     1,.29,.29,.43,.43,.29,.29,.43,.43,.29,.29,.43,.43,.29,.29,.43,.43,0
     2.,0.,.14,.14,0.,0.,.14,.14,0.,0.,.14,.14,0.,0.,.14,.14,.29,.29,.43
     3,.43,.29,.29,.43,.43,.29,.29,.43,.43,.29,.29,.43,.43,0.,0.,.14,.14
     4,0.,0.,.14,.14,0.,0.,.14,.14,0.,0.,.14,.14,.29,.29,.43,.43,.29,.29
     5,.43,.43,.29,.29,.43,.43,.29,.29,.43,.43,0.,0.,.14,.14,0.,0.,.14,.
     614,0.,0.,.14,.14,0.,0.,.14,.14,.29,.29,.43,.43,.29,.29,.43,.43,.29
     7,.29,.43,.43,.29,.29,.43,.43,.57,.57,.71,.71,.57,.57,.71,.71,.57,.
     857,.71,.71,.57,.57,.71,.71,.86,.86,1.,1.,.86,.86,1.,1.,.86,.86,1.,
     91.,.86,.86,1.,1.,.57,.57,.71,.71,.57,.57,.71,.71,.57,.57,.71,.71,.
     a57,.57,.71,.71,.86,.86,1.,1.,.86,.86,1.,1.,.86,.86,1.,1.,.86,.86,1
     b.,1.,.57,.57,.71,.71,.57,.57,.71,.71,.57,.57,.71,.71,.57,.57,.71,.
     c71,.86,.86,1.,1.,.86,.86,1.,1.,.86,.86,1.,1.,.86,.86,1.,1.,.57,.57
     d,.71,.71,.57,.57,.71,.71,.57,.57,.71,.71,.57,.57,.71,.71,.86,.86,1
     e.,1.,.86,.86,1.,1.,.86,.86,1.,1.,.86,.86,1.,1./
      data greenx /0.,.14,0.,.14,0.,.14,0.,.14,.29,.43,.29,.43,.29,.43,.
     129,.43,0.,.14,0.,.14,0.,.14,0.,.14,.29,.43,.29,.43,.29,.43,.29,.43
     2,0.,.14,0.,.14,0.,.14,0.,.14,.29,.43,.29,.43,.29,.43,.29,.43,0.,.1
     34,0.,.14,0.,.14,0.,.14,.29,.43,.29,.43,.29,.43,.29,.43,.57,.71,.57
     4,.71,.57,.71,.57,.71,.86,1.,.86,1.,.86,1.,.86,1.,.57,.71,.57,.71,.
     557,.71,.57,.71,.86,1.,.86,1.,.86,1.,.86,1.,.57,.71,.57,.71,.57,.71
     6,.57,.71,.86,1.,.86,1.,.86,1.,.86,1.,.57,.71,.57,.71,.57,.71,.57,.
     771,.86,1.,.86,1.,.86,1.,.86,1.,0.,.14,0.,.14,0.,.14,0.,.14,.29,.43
     8,.29,.43,.29,.43,.29,.43,0.,.14,0.,.14,0.,.14,0.,.14,.29,.43,.29,.
     943,.29,.43,.29,.43,0.,.14,0.,.14,0.,.14,0.,.14,.29,.43,.29,.43,.29
     a,.43,.29,.43,0.,.14,0.,.14,0.,.14,0.,.14,.29,.43,.29,.43,.29,.43,.
     b29,.43,.57,.71,.57,.71,.57,.71,.57,.71,.86,1.,.86,1.,.86,1.,.86,1.
     c,.57,.71,.57,.71,.57,.71,.57,.71,.86,1.,.86,1.,.86,1.,.86,1.,.57,.
     d71,.57,.71,.57,.71,.57,.71,.86,1.,.86,1.,.86,1.,.86,1.,.57,.71,.57
     e,.71,.57,.71,.57,.71,.86,1.,.86,1.,.86,1.,.86,1./
      data bluex /0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.,.14,.29,.43,.2
     19,.43,0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,
     2.57,.71,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1.,.86,1.,.57,.7
     31,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1.,.86,1.,0.,.14,0.,.1
     44,.29,.43,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.,.14,.29,
     5.43,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,.57,.71,.57,.71,.86,1.,.
     686,1.,.57,.71,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1.,.86,1.,
     7.57,.71,.57,.71,.86,1.,.86,1.,0.,.14,0.,.14,.29,.43,.29,.43,0.,.14
     8,0.,.14,.29,.43,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.,.1
     94,.29,.43,.29,.43,.57,.71,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.8
     a6,1.,.86,1.,.57,.71,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1.,.
     b86,1.,0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,
     c0.,.14,0.,.14,.29,.43,.29,.43,0.,.14,0.,.14,.29,.43,.29,.43,.57,.7
     d1,.57,.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1.,.86,1.,.57,.71,.57,
     e.71,.86,1.,.86,1.,.57,.71,.57,.71,.86,1.,.86,1./
c clear pixel lookup table flag
      lupt = 0
c ntc = number of colors requested
      ntc = 2**nbit
c for default palette, only 256 colors allowed
      if ((npal.eq.0).and.(ntc.gt.256)) ntc = 256
c do not use more colors than specified in palette
      if ((npal.gt.0).and.(ntc.gt.npal)) ntc = npal
c save number of colors requested
      ncr = ntc
c do not use more colors than workstation supports
      ncoli = nclsp + 1
      if (ntc.gt.ncoli) then
c set pixel lookup table flag
         lupt = 1
c determine how many colors are possible
         if (npal.gt.0) ntc = 128
         if (ntc.gt.ncoli) ntc = 64
         if (ntc.gt.ncoli) ntc = 32
         if (ntc.gt.ncoli) ntc = 16
         if (ntc.gt.ncoli) ntc = 8
         if (ntc.gt.ncoli) ntc = 4
         if (ntc.gt.ncoli) ntc = 2
         write (6,*) 'too many colors requested: ncr,ncoli = ',ncr,ncoli
         write (6,*) 'only ',ntc,' will be used'
      endif
c set default palette
      if ((npal.eq.0).or.(lupt.eq.1)) then
         do 10 j = 1, ntc
         jc = j - 1
c 128 or fewer colors
         if (ntc.le.128) then
            icol = j + ntc - 2
c set color representation
            call gscr(idwk,jc,reds(icol),greens(icol),blues(icol))
c 256 colors
         else
c set color representation
            call gscr(idwk,jc,redx(j),greenx(j),bluex(j))
         endif
   10    continue
c set specified palette
      else
         apmx = 1./real(ipmx - 1)
         do 20 j = 1, ntc
         jc = j - 1
         j1 = 3*(j - 1)
         red = real(ichar(pal(j1+1)))*apmx
         green = real(ichar(pal(j1+2)))*apmx
         blue = real(ichar(pal(j1+3)))*apmx
c set color representation
         call gscr(idwk,jc,red,green,blue)
   20    continue
      endif
c find which indices correspond to prime colors
      apmx = 1./real(ipmx - 1)
      do 40 i = 1, 8
      jmin = 0
      cdm = 3.0
      do 30 j = 1, ntc
      jc = j - 1
c default palette
      if ((npal.eq.0).or.(lupt.eq.1)) then
c 128 or fewer colors
         if (ntc.le.128) then
            icol = j + ntc - 2
            cd = (reds(icol) - redp(i))**2 + (greens(icol) - greenp(i))*
     1*2 + (blues(icol) - bluep(i))**2
c 256 colors
         else
            cd = (redx(j) - redp(i))**2 + (greenx(j) - greenp(i))**2  + 
     1(bluex(j) - bluep(i))**2
         endif
c specified palette
      else
         j1 = 3*(j - 1)
         cd = (real(ichar(pal(j1+1)))*apmx - redp(i))**2 + (real(ichar(p
     1al(j1+2)))*apmx - greenp(i))**2  + (real(ichar(pal(j1+3)))*apmx - 
     2bluep(i))**2
      endif
      if (cd.le.cdm) then
         jmin = jc
         cdm = cd
      endif
   30 continue
      kprime(i) = jmin
   40 continue
c swap colors for index 0 and background
      if (kprime(1).ne.0) then
            jc = kprime(1)
      if ((npal.eq.0).or.(lupt.eq.1)) then
c 128 or fewer colors
            if (ntc.le.128) then
               icol = ntc - 1
               r0 = reds(icol)
               g0 = greens(icol)
               b0 = blues(icol)
               icol = icol + jc
c set color representation
               call gscr(idwk,0,reds(icol),greens(icol),blues(icol))
               call gscr(idwk,jc,r0,g0,b0)
c 256 colors
            else
               r0 = redx(1)
               g0 = greenx(1)
               b0 = bluex(1)
               icol = jc + 1
c set color representation
               call gscr(idwk,0,redx(icol),greenx(icol),bluex(icol))
               call gscr(idwk,jc,r0,g0,b0)
            endif
c specified palette
         else
            apmx = 1./real(ipmx - 1)
            r0 = real(ichar(pal(1)))*apmx
            g0 = real(ichar(pal(2)))*apmx
            b0 = real(ichar(pal(3)))*apmx
            j1 = 3*jc
            red = real(ichar(pal(j1+1)))*apmx
            green = real(ichar(pal(j1+2)))*apmx
            blue = real(ichar(pal(j1+3)))*apmx
c set color representation
            call gscr(idwk,0,red,green,blue)
            call gscr(idwk,jc,r0,g0,b0)
         endif
         if (kprime(2).eq.0) kprime(2) = jc
c set foreground color index
      endif
      ifrg = kprime(2)
c     ifrg = ntc - 1
c pixel lookup table not needed
      if (lupt.eq.0) return
c create pixel lookup table
      itc = 512/ntc
c mapping from default 128 colors to fewer
      if (npal.eq.0) then
         aptq = 8.
         do 50 j = 1, ncr
         icol = j + ncr - 2
         ip = reds(icol)*aptq
         if (ip.gt.7) ip = 7
         it = 48*(ip/4) + 6*(ip/2) + ip
         ip = greens(icol)*aptq
         if (ip.gt.3) ip = 3
         it = 2*it + (48*(ip/4) + 6*(ip/2) + ip)
         ip = blues(icol)*aptq
         if (ip.gt.3) ip = 3
         it = 2*it + (48*(ip/4) + 6*(ip/2) + ip)
         ipal(j) = it/itc
   50    continue
c mapping from 256 specified colors to 128 or fewer colors
      else
         len = npal
         if (len.gt.npald) len = npald
         iptq = ipmx/8
         do 70 j = 1, len
         j1 = 3*(j - 1)
         it = 0
         do 60 i = 1, 3
         ip = ichar(pal(j1+i))/iptq
         it = 2*it + (48*(ip/4) + 6*(ip/2) + ip)
   60    continue
         ipal(j) = it/itc
   70    continue
      endif
      return
      end
      subroutine GRCLOSE
c this subroutine deactivates workstation and closes gks
c idwk = workstation identifier
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c pause if plots are still pending
      if (((iplot.ne.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         call readrc(irc)
      endif
c deactivate workstation
      call gdawk(idwk)
c close workstation
      call gclwk(idwk)
c close gks
      call gclks
      return
      end
      subroutine SETNPLT (nplt,irc)
c this subroutine resets the maximum number of plots per page
c if requested number is negative, it is set to default (=1)
c the current plot location is also reset to the next available location
c if next available location is iplot = 0 and the old location was not
c 0, then the workstation is updated and a return code can be generated.
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c set return code to normal
      irc = 0
c suppress plots (nplt = 0)
      if (nplt.eq.0) then
         iplot = 0
         nplot = 0
         return
      endif
c calculate old screen partition
      npl1 = sqrt(real(nplot-1)) + 0.0001
      npl = npl1 + 1
      apl = 1./real(npl)
c find location of previous plot
      iplt = iplot - 1
      if (iplt.lt.0) iplt = iplt + nplot
c find coordinates of previous plot
      iy = iplt/npl
      ix = iplt - iy*npl + 1
      iy = npl1 - iy
c reset maximum number of plots per page
      if (nplt.ge.0) then
         nplot = nplt
      else
         nplot = 1
      endif
c calculate new screen partition
      npl1 = sqrt(real(nplot-1)) + 0.0001
      npl = npl1 + 1
c find plot coordinates of previous plot in new screen partition
      px = apl*real(ix) - 0.0001
      py = apl*real(iy) + 0.0001
      ix = abs(px)*real(npl)
      iy = abs(1. - py)*real(npl)
c find new plot location
      iplt = 1 + (ix + iy*npl)
      if (iplt.ge.nplot) then
c pause if plots are still pending
         if (((iplot.ne.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
            call guwk(idwk,1)
c read code from input device, if present
            call readrc(irc)
         endif
         iplot = 0
      else
         iplot = iplt
      endif
      return
      end
      subroutine DISPR (f,label,xmin,xmax,isc,ist,mks,nx,nxv,ngs,chr,chr
     1s,irc)
c this subroutine plots ngs subarrays of the array f, on a common graph,
c each plot with nx points, versus a linear function in x,
c where xmin < x < xmax.
c depending on the number of colors in the display device, each subarray
c is plotted with a different color, given in order by:
c blue, red, yellow, cyan, magenta, green, and foreground.
c after all the colors are cycled through, then different line styles
c are cycled through if mks=0, in order: solid, dash, dot, dash-dot,
c or different marker types if mks>0: dot, plus, star, circle, cross.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = array containing subarrays to be plotted
c label = long character string label for plot
c xmin/xmax = range of x values in plot
c isc = power of 2 scale of y coordinate for plot
c ist = flag for choosing positive and/or negative values
c the plots have a common scale in y given by ymax and ymin.
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist = 1, then ymax = 2**isc and ymin = 0.
c if ist = -1, then ymax = 0 and ymin = -2**isc.
c if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c mks = flag to determine whether lines or markers are used,
c mks=0 means cycle through lines styles, mks > 0 means cycle through
c marker styles, using the value of mks as the initial marker,
c mks < 0 means draw the first subarray with a line, then subsequent
c subarrays with markers, using the value of abs(mks) as the initial
c marker.
c nx = number of points plotted in each subarray
c nxv = first dimension of array f, nxv >= nx
c ngs = second dimension of array f, number of subarrays to be plotted
c chr = additional long character string comment for plot
c chrs = array of ngs short character labels used by subroutine tickd
c to label individual line or marker samples
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      parameter(nxbs=65)
      character*(*) label, chr
      character*(*) chrs(ngs)
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
      dimension f(nxv,ngs)
c x,y = scratch arrays for plotting
      dimension x(nxbs), y(nxbs)
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c nlts = number of line types available
c nmks = number of markers available
c ntx/nty = number of ticks in grid in x/y direction
      data nlts,nmks,ntx,nty /4,5,11,9/
c csize = vertical size of characters
      data csize /0.034/
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
         do 20 k = 1, ngs
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
            if (fmax.eq.fmin) fmax = fmin + 1.0
         endif
      else
         fmax = dv**is
         fmin = -fmax
      endif
      ymax = fmax
      ymin = fmin
      if (ist.eq.1) then
         ymin = 0.
      else if (ist.eq.(-1)) then
         ymax = 0.  
      endif
c parameters for plots
      dx = xmax - xmin
      if (nx.gt.1) dx = dx/real(nx - 1)
      nxs = nxbs - 1
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
c draw grid and labels, call identity transformation
      call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample lines or markers
      call dsmpln(orx,smn,tmn,tmx,ngs,nlts,nmks,mks,chrs,chh)
c define transformation number 2
      nrt = 1
c set window
      call gswn(nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp(nrt,smn,smx,tmn,tmx)
c select normalization transformation
      call gselnt(nrt)
c set linewidth scale factor, 1.0 = nominal
      call gslwsc(1.0)
      mkr = abs(mks)
c use markers
      if (mkr.ne.0) then
c set marker size scale factor, 1.0 = nominal
         call gsmksc(1.0)
      endif
c set clipping indicator, 1 = on
      call gsclip(1)
c plot curves, first cycle through colors, then line or marker types
      do 70 k = 1, ngs
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
c use line types
      if ((mks.eq.0).or.((mks.lt.0).and.(k.eq.1))) then
c blocking parameters for plots
         nxb = (nx - 2)/nxs + 1
         npts = nxbs
c length of last block
         nptl = nx - nxs*(nxb - 1)
c set polyline color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gsplci(icol)
         ltype = 1 + (k - 1)/ncols - nlts*((k - 1)/(nlts*ncols))
c set linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
         call gsln(ltype)
c plot curve
         npt = npts
c loop over number of blocks
         do 40 j = 1, nxb
         js = nxs*(j - 1)
         if (j.eq.nxb) npt = nptl
c calculate x,y axes for block
         do 30 i = 1, npt
         x(i) = xmin + dx*real(i + js - 1)
         y(i) = f(i+js,k)
   30    continue
c draw polyline
         call gpl(npt,x,y)
   40    continue
c use markers
      else
c blocking parameters for plots
         nxb = (nx - 1)/nxs + 1
         npts = nxs
c length of last block
         nptl = nx - nxs*(nxb - 1)
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
         mtype = mkr + (k - 1)/ncols - nmks*((mkr - 1 + (k - 1)/ncols)/n
     1mks)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mtype)
c plot polymarkers
         npt = npts
c loop over number of blocks
         do 60 j = 1, nxb
         js = nxs*(j - 1)
         if (j.eq.nxb) npt = nptl
c calculate x,y axes for block
         do 50 i = 1, npt
         x(i) = xmin + dx*real(i + js - 1)
         y(i) = f(i+js,k)
   50    continue
c dots
         if (mtype.eq.1) then
c treat dots by drawing a line to itself
            call spdots(x,y,npt,icol,nxbs)
         else
c draw polymarker
            call gpm(npt,x,y)
         endif
   60    continue
      endif
   70 continue
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
      subroutine DISPC (f,g,label,zsc,zst,mks,nx,nxv,ngs,chr,chrs,irc)
c this subroutine plots ngs subarrays of the array f, on a common graph,
c each plot with nx points, versus the corresponding subarray of the
c array g.
c depending on the number of colors in the display device, each subarray
c is plotted with a different color, given in order by:
c blue, red, yellow, cyan, magenta, green, and foreground
c after all the colors are cycled through, then different lines styles
c are cycled through if mks=0, in order: solid, dash, dot, dash-dot,
c or different marker types if mks>0: dot, plus, star, circle, cross.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f, g = arrays containing subarrays to be plotted
c label = long character string label for plot
c real(zsc)/aimag(zsc) = power of 2 scale of x/y coordinate for plot
c real(zst)/aimag(zst) = flag for positive and/or negative x/y values.
c the plots have a common scale in y given by ymax and ymin,
c where isc = int(aimag(zsc)) and ist = int(aimag(zst)), as follows:
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist = 1, then ymax = 2**isc and ymin = 0.
c if ist = -1, then ymax = 0 and ymin = -2**isc.
c if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c the plots have a common scale in x given by xmax and xmin,
c where isc = int(real(zsc)) and ist = int(real(zst)), as follows:
c if ist = 0, then xmax = 2**isc and xmin = -2**isc.
c if ist = 1, then xmax = 2**isc and xmin = 0.
c if ist = -1, then xmax = 0 and xmin = -2**isc.
c if ist = 2, then xmin = gmin, xmax = gmin + 2**ir,
c where gmin/gmax are the function minimum/maximum, 
c and ir = power of 2 scale for (gmax - gmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of g
c mks = flag to determine whether lines or markers are used,
c mks=0 means cycle through lines styles, mks > 0 means cycle through
c marker styles, using the value of mks as the initial marker,
c mks < 0 means draw the first subarray with a line, then subsequent
c subarrays with markers, using the value of abs(mks) as the initial
c marker.
c nx = number of points plotted in each subarray
c nxv = first dimension of array f, nxv >= nx
c ngs = second dimension of array f, number of subarrays to be plotted
c chr = additional character string comment for plot
c chrs = array of ngs short character labels used by subroutine tickd
c to label individual line or marker samples
c irc = return code (0 = normal return)
      character*(*) label, chr
      character*(*) chrs(ngs)
      complex zsc, zst
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
      dimension f(nxv,ngs), g(nxv,ngs)
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c nlts = number of line types available
c nmks = number of markers available
c ntx/nty = number of ticks in grid in x/y direction
      data nlts,nmks,ntx,nty /4,5,9,9/
c csize = vertical size of characters
      data csize /0.034/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c find y scales for plot
      is = aimag(zsc)
      ist = aimag(zst)
      if (abs(is).gt.116) then
         fmax = f(1,1)
         fmin = fmax
         do 20 k = 1, ngs
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
      ymax = fmax
      ymin = fmin
      if (ist.eq.1) then
         ymin = 0.
      else if (ist.eq.(-1)) then
         ymax = 0.  
      endif
c find x scales for plot
      is = real(zsc)
      ist = real(zst)
      if (abs(is).gt.116) then
         gmax = g(1,1)
         gmin = gmax
         do 40 k = 1, ngs
         do 30 j = 1, nx
         gmax = amax1(gmax,g(j,k))
         gmin = amin1(gmin,g(j,k))
   30    continue
   40    continue
         if (gmax.eq.0.) gmax = 1.0e-35
         rmax = gmax - gmin
         rmin = gmin
         xmax = abs(gmax)
         is = alog(xmax)*algdvi
         if (xmax.ge.1.) is = is + 1
         if (xmax.le.dv**(is-1)) is = is - 1
         xmin = abs(gmin)  
         if (xmin.gt.0.) then
            it = alog(xmin)*algdvi
            if (xmin.ge.1.) it = it + 1
            if (xmin.le.dv**(it-1)) it = it - 1
         endif
         if (gmax.gt.0.) then
            if (gmin.gt.0.) then
               gmin = dv**(it - 1)
            else if (gmin.lt.0.) then
               gmin = -dv**it
            endif
            gmax = dv**is
         else
            gmax = -dv**(is - 1)
            gmin = -dv**it
         endif
         if (ist.eq.0) then
            if (xmin.gt.xmax) then
               gmax = dv**it
            else
               gmin = -dv**is
            endif
         else if (ist.eq.2) then
            ir = alog(rmax)*algdvi
            if (rmax.ge.1.) ir = ir + 1
            if (rmax.le.dv**(ir-1)) ir = ir - 1
            gmin = rmin
            gmax = rmin + dv**ir
         endif
      else
         gmax = dv**is
         gmin = -gmax
      endif
      xmax = gmax
      xmin = gmin
      if (ist.eq.1) then
         xmin = 0.
      else if (ist.eq.(-1)) then
         xmax = 0.  
      endif
      npts = nx - 1
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
c draw grid and labels, call identity transformation
      call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample lines or markers
      call dsmpln(orx,smn,tmn,tmx,ngs,nlts,nmks,mks,chrs,chh)
c define transformation number 2
      nrt = 1
c set window
      call gswn(nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp(nrt,smn,smx,tmn,tmx)
c select normalization transformation
      call gselnt(nrt)
c set linewidth scale factor, 1.0 = nominal
      call gslwsc(1.0)
      mkr = abs(mks)
c use markers
      if (mkr.gt.0) then
c set marker size scale factor, 1.0 = nominal
         call gsmksc(1.0)
      endif
c set clipping indicator, 1 = on
      call gsclip(1)
c plot curves, first cycle through colors, then line or marker types
      do 50 k = 1, ngs
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
c use line types
      if ((mks.eq.0).or.((mks.lt.0).and.(k.eq.1))) then
c set polyline color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gsplci(icol)
         ltype = 1 + (k - 1)/ncols - nlts*((k - 1)/(nlts*ncols))
c set linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
         call gsln(ltype)
c plot curve
c draw polyline
         call gpl(nx,g(1,k),f(1,k))
c use markers
      else
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
         mtype = mkr + (k - 1)/ncols - nmks*((mkr - 1 + (k - 1)/ncols)/n
     1mks)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mtype)
c plot polymarkers
         if (mtype.eq.1) then
c treat dots by drawing a line to itself
            call spdots(g(1,k),f(1,k),nx,icol,nxv)
         else
c draw polymarker
            call gpm(nx,g(1,k),f(1,k))
         endif
      endif
   50 continue
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
      subroutine DISPS (f,label,xmin,xmax,isc,ist,nx,chr,irc)
c this subroutine plots an array f versus a linear function in x,
c where xmin < x < xmax.  It is plotted in solid line style, in blue
c if color is available.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = array to be plotted
c label = long character string label for plot
c xmin/xmax = range of x values in plot
c isc = power of 2 scale of y coordinate for plot
c ist = flag for choosing positive and/or negative values
c the plot has a scale in y given by ymax and ymin.
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist = 1, then ymax = 2**isc and ymin = 0.
c if ist = -1, then ymax = 0 and ymin = -2**isc.
c if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plot, determined by the absolute value of f
c nx = dimension of array f, the number of points to be plotted.
c chr = additional long character string comment for plot
c irc = return code (0 = normal return)
      dimension f(nx)
      character*(*) label, chr
      character*1 chs(1)
      data chs /' '/
      call DISPR (f,label,xmin,xmax,isc,ist,0,nx,nx,1,chr,chs,irc)
      return
      end
      subroutine DISPD (f,g,label,zsc,zst,nx,chr,irc)
c this subroutine plots an array f versus an array g, in a solid line
c style, in blue if color is available.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f, g = arrays to be plotted
c label = long character string label for plot
c real(zsc)/aimag(zsc) = power of 2 scale of x/y coordinate for plot
c real(zst)/aimag(zst) = flag for positive and/or negative x/y values.
c the plot has a scale in y given by ymax and ymin,
c where isc = int(aimag(zsc)) and ist = int(aimag(zst)), as follows:
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist > 0, then ymax = 2**isc and ymin = 0.
c if ist < 0, then ymax = 0 and ymin = -2**isc.
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plot, determined by the absolute value of f
c the plot has a scale in x given by xmax and xmin,
c where isc = int(real(zsc)) and ist = int(real(zst)), as follows:
c if ist = 0, then xmax = 2**isc and xmin = -2**isc.
c if ist > 0, then xmax = 2**isc and xmin = 0.
c if ist < 0, then xmax = 0 and xmin = -2**isc.
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plot, determined by the absolute value of g
c nx = dimension of arrays f and g, the number of points to be plotted.
c chr = additional character string comment for plot
c irc = return code (0 = normal return)
      dimension f(nx), g(nx)
      character*(*) label, chr
      complex zsc, zst
      character*1 chs(1)
      data chs /' '/
      call DISPC (f,g,label,zsc,zst,0,nx,nx,1,chr,chs,irc)
      return
      end
      subroutine tickd(xmin,xmax,ymin,ymax,orgx,orgy,smin,smax,tmin,tmax
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
      character*10 lbl
c rx, ry = ndc coordinates of upper-right corner of workstation window
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c x,y = scratch variables for plotting
      dimension x(5), y(5)
   91 format (e10.3)
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
c identify leading blanks
      is = 0
      ls = len(lbl)
   30 is = is + 1
      if ((lbl(is:is).eq.' ').and.(is.lt.ls)) go to 30
      ax = smin
c draw text
      call gtx(ax,ay,lbl(is:ls))
c set text alignment to (right,normal)
      call gstxal(3,0)
      write (lbl,91) xmax
      ax = smax
c draw text
      call gtx(ax,ay,lbl)
c label y axes
      ax = smin - 2.*stx
      write (lbl,91) ymin
      ay = tmin
c draw text
      call gtx(ax,ay,lbl)
      write (lbl,91) ymax
      ay = tmax - chh
c draw text
      call gtx(ax,ay,lbl)
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
      subroutine dsmpln(orgx,smin,tmin,tmax,ngs,nlts,nmks,mks,chrs,chh)
c this subroutine displays line or marker samples, with short character
c labels placed underneath
c orgx = x origin of window
c smin = minimum x value of plotting window
c tmin/tmax = range of y values of plotting window
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
            call spdots(x(3),y,2,icol,3)
         else
c draw polymarker
            call gpm(2,x(3),y)
         endif
      endif
c draw label underneath sample line or marker
      y(1) = y(1) - 1.5*chh
c draw text
      call gtx(ax,y(1),chrs(k))
   30 continue
      return
      end
      subroutine idtran
c this subroutine performs the identify transformation number 1
c rx, ry = ndc coordinates of upper-right corner of workstation window
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      data zero /0./
c define transformation number 1
      nrt = 1
c set window
      call gswn(nrt,zero,rx,zero,ry)
c set viewport
      call gsvp(nrt,zero,rx,zero,ry)
c select normalization transformation
      call gselnt(nrt)
      return
      end
      subroutine dfplps (chh,nfont,iprec)
c this subroutine sets default plotting parameters
c chh = character height, in world coordinates
c nfont = character font number
c iprec = text precision (0=string,1=character,2=stroke)
c ifrg = index of foreground color
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c set linewidth scale factor, 1.0 = nominal
      call gslwsc(1.0)
c set polyline color index to foreground
      call gsplci(ifrg)
c set linetype, 1 = solid
      call gsln(1)
c set character height
      call gschh(chh)
c set text color index to foreground
      call gstxci(ifrg)
c set text font
      call gstxfp(nfont,iprec)
c set text alignment to (normal,normal)
      call gstxal(0,0)
      return
      end
      subroutine spdots(x,y,npt,icol,nxbs)
c this subroutine draws dot markers by drawing a line to itself
c x, y = arrays to be plotted
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
      xs(2) = xs(1)
      ys(2) = ys(1)
c draw polyline
      call gpl(2,xs,ys)
   10 continue
      return
      end
      subroutine readrc(irc)
c this subroutine reads return code from input device
c if mouse = 0, string device is used if present, otherwise locator
c if mouse = 1, locator device is used if present, otherwise string
c if neither device is present, the subroutine exits with irc = 0
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c mouse = (0,1) = (no,yes) mouse has higher priority than keyboard
      data mouse /0/
      save mouse
   10 if (mouse.eq.0) then
         if (idstr.gt.0) then
c read input code from keyboard, if present
            call readcc(irc)
         else if (idloc.gt.0) then
c read mouse location, if present
            call readmc(irc)
         endif
      else
         if (idloc.gt.0) then
c read mouse location, if present
            call readmc(irc)
         else if (idstr.gt.0) then
c read input code from keyboard, if present
            call readcc(irc)
         endif
      endif
c change input device, if requested
      if (irc.eq.5) then
         mouse = 1 - mouse
c if input device changes to mouse, read mouse code
         if (mouse.eq.1) go to 10
      endif
      return
      end
      subroutine readcc(irc)
c this subroutine reads characters and outputs an integer code irc
c characters which return codes are:
c q (irc=1), s (irc=2), m (irc=3), a (irc=4), u (irc=5), r (irc=6),
c p (irc=7), and n (irc=8)
c if a non-negative number num is entered, irc = num + 128. 
c a carriage return or any other characters gives irc=0.
c character '?' gives a help, then tries to read characters again.
c string device can be either in request (wait) or event (nowait) mode
c to go from request to event mode, enter character "a" (for animate)
c to go from event to request mode, enter carriage return
c if irc = -1 on entry, string device is forced into request mode
c and help is displayed
c if irc = -2 on entry, string device is forced into event mode
c if irc = -3 on entry, string device is forced into request mode
c and readcc returns immediately
c on first entry, string device is assumed to be in request mode
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c idstr = string device number, 0 if no string device available
c ifrg = index of foreground color
      parameter(nlines=3)
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      character*40 chr(nlines)
      character*8 str
      character*1 c
c mode = mode of operation of string device (0=request,1=sample,2=event)
c iesw = echo switch (0=noecho,1=echo)
      data mode,iesw /0,1/
c csize = height of characters
c tout = timeout period, in seconds
      data zero,csize,tout /0.,0.033,0./
c help information
      data chr /' q=quit, s=save, m=modify, a=animate    ',' u=use mouse
     1, r=reset, p=plotparms      ',' cr=continue, # = display frame num
     2ber #'/
      save mode,iesw,tout
c if irc = -1 on input, put string device into request mode
      if ((idstr.gt.0).and.(irc.eq.(-1))) then
         mode = 0
c set string mode, 0 = request
         call gsstm(idwk,idstr,mode,iesw)
c set string to display help
         str(1:1) = '?'
c if irc = -2 on input, put string device into event mode
      elseif ((idstr.gt.0).and.(irc.eq.(-2))) then
         mode = 2
c set string mode, 2 = event
         call gsstm(idwk,idstr,mode,iesw)
c if irc = -3 on input, put string device into request mode and quit
      elseif ((idstr.gt.0).and.(irc.eq.(-3))) then
         mode = 0
c set string mode, 0 = request
         call gsstm(idwk,idstr,mode,iesw)
         return
c clear return code
      else
         irc = 0
      endif
c define and select identity transformation
      call idtran
c set text color index to foreground
      call gstxci(ifrg)
      if (idstr.gt.0) then
c if string device is in event mode, check for events
         if (mode.eq.2) then
c await event
   10       call gwait(tout,idwks,icl,idnr)
c return if no event available
            if (icl.eq.0) return
c check again if event is not string
            if (icl.ne.6) go to 10
c check again if workstation or device id is not correct
            if ((idwks.ne.idwk).or.(idnr.ne.idstr)) go to 10
            str = '        '
            mode = 0
c set string mode, 0 = request
            call gsstm(idwk,idstr,mode,iesw)
c get string
            call ggtst(lens,str)
            irc = -1
         endif
c if string device is in request mode, wait for input
         if (mode.eq.0) then
            if (irc.lt.0) go to 30
   20       str = '        '
c request string
            call grqst(idwk,idstr,istat,lens,str)
   30       irc = 0
c clear plotparm return code
            irp = 0
            c = str(1:1)
            if (c.eq.'?') then
c clear workstation, always
               call gclrwk(idwk,1)
c set default plotting parameters
               chh = csize*ry
c 1 = default font, 0 = string precision
               call dfplps (chh,1,0)
               ay = ry - 2.*chh
c write out help
               do 40 i = 1, nlines
c draw text
               call gtx(zero,ay,chr(i))
               ay = ay - 2.*chh
   40          continue
c update workstation, perform
               call guwk(idwk,1)
               irc = -1
            else
               if ((c.eq.'Q').or.(c.eq.'q')) irc = 1
               if ((c.eq.'S').or.(c.eq.'s')) irc = 2
               if ((c.eq.'M').or.(c.eq.'m')) irc = 3
               if ((c.eq.'A').or.(c.eq.'a')) then
                  irc = 4
                  mode = 2
c set string mode, 2 = event
                  call gsstm(idwk,idstr,mode,iesw)
               endif
               if ((c.eq.'U').or.(c.eq.'u')) irc = 5
               if ((c.eq.'R').or.(c.eq.'r')) irc = 6
               if ((c.eq.'P').or.(c.eq.'p')) then
                  irc = 7
c set plotparm return code
                  irp = 1
c change plot parameters
                  call plparmc(zero,ry)
c set string to display help
                  str(1:1) = '?'
               endif
               if ((c.eq.'N').or.(c.eq.'n')) irc = 8
            endif
c go read characters again
            if (irc.lt.0) go to 20
c display help after changing plot parameters
            if (irp.eq.1) go to 30
c check if characters represent a number
            if (irc.eq.0) then
               i = 1
   50          iv = ichar(str(i:i)) - ichar('0')
c abort if not a number
               if ((iv.lt.0).or.(iv.gt.9)) go to 60
               irc = iv + 10*irc
               i = i + 1
               if (i.le.lens) go to 50
   60          if (i.gt.1) irc = 128 + irc
            endif
         endif
      endif
      return
      end
      subroutine plparmc(orgx,orgy)
c this subroutine allows one to make changes to certain of the plotting
c parameters stored in the common blocks plotcm and pextra.
c currently, only the variables nplot and ndi can be changed.
c variables are changed by typing variable name = new value.  If the
c new value is acceptable, type run.  to abort changes, type quit.
c on entry, string device is assumed to be in request mode
c orgx/orgy = coordinates of upper left hand corner of writing area,
c in the range, 0 < orgx < rx and 0 < orgy < ry
c nvars = number of variables to be modified
      parameter (nvars=2)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c nplot = number of plots per page
c idstr = string device number, 0 if no string device available
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c ndi = increment between frames
      common /pextra/ ndi
      character*32 input
      character*30 output
      character*40 prompt
      character*8 codev
      dimension nvalue(nvars)
      character*8 code(nvars)
      character*40 helpv(nvars)
      character*44 helpr(11)
      save prompt,code,helpv,helpr
   91 format (1x,a6,' = ',i4,', ',a4,' = ',i4,3x)
      data prompt /' enter variable=value, or run, quit or ?'/
      data code /'NPLOT   ','NDI     '/
      data helpv  /' NPLOT = number of plots per page       ',' NDI = in
     1crement between frames         '/
      data helpr /'    this program interactively modifies and ','update
     1s plot parameters.                    ','    to change parameters, 
     2type name=number, ','where name = the parameter to be changed and'
     3,'number is the new value, e.g., ndi=7        ','    to obtain the
     4 meaning of a parameter,   ','type its name, e.g., ndi            
     5        ','   to display current parameters, hit return','or enter
     6 key.                               ','to abort the modification p
     7rocess, type quit','when parameters are OK, type run.           '/
c csize = height of characters
      data csize /0.033/
c exit if no string device present
      if (idstr.eq.0) return
c first copy variables
      nvalue(1) = nplot
      nvalue(2) = ndi
c normalize character sizes
      chh = csize*ry
      chh2 = 2.*chh
c set default plotting parameters
      ax = orgx
c clear workstation, always
   10 call gclrwk(idwk,1)
c set default plotting parameters
      ay = orgy - chh2
c 1 = default font, 0 = string precision
      call dfplps (chh,1,0)
c display variable names and values
      write (output,91) code(1), nvalue(1), code(2), nvalue(2)
c draw text
      call gtx(ax,ay,output)
c display prompt
   20 ay = ay - chh2
      if (ay.lt.0.) go to 10
c draw text
      call gtx(ax,ay,prompt)
c update workstation, perform
      call guwk(idwk,1)
   30 input = '                                '
c request string
      call grqst(idwk,idstr,istat,lens,input)
c find if there is an equal sign
      i = 0
   40 i = i + 1
      if (i.gt.lens) go to 50
      if (input(i:i).ne.'=') go to 40
c nv = length of variable name
   50 nv = i - 1
      nv = min0(nv,len(codev))
c no variable name entered, update display
      if (nv.lt.1) go to 10
c convert name to upper case and left justify
      codev = '        '
      i = 0
      ic = 0
   60 i = i + 1
      if (i.gt.nv) go to 70
      iv = ichar(input(i:i))
c ignore blanks
      if (iv.eq.ichar(' ')) go to 60
c ic = length of name, with blanks stripped
      ic = ic + 1
c convert to upper case, if necessary
      if ((iv.ge.ichar('a')).and.(iv.le.ichar('z'))) iv = iv - ichar('a' 
     1)+ ichar('A')
      codev(ic:ic) = char(iv)
      go to 60
c check if command or '?' was entered
   70 if ((codev(1:ic).eq.'QUIT').or.(codev(1:1).eq.'Q')) return
      if (codev(1:1).eq.'?') then
c write help instructions
c clear workstation, always
         call gclrwk(idwk,1)
c set default plotting parameters
         ay = orgy
c 1 = default font, 0 = string precision
         call dfplps (chh,1,0)
         do 80 i = 1, 11
         ay = ay - chh2
c draw text
         call gtx(ax,ay,helpr(i))
   80    continue
c update workstation, perform
         call guwk(idwk,1)
         go to 30
      elseif (codev(1:ic).eq.'RUN') then
c update variables
         nplot = nvalue(1)
         ndi = nvalue(2)
         return
      endif
c check if name is known variable
      nvar = 0
      i = 0
   90 i = i + 1
      if (i.gt.nvars) go to 100
      if (codev.ne.code(i)) go to 90
c found name
      nvar = i
c if name not found, update display
  100 if (nvar.eq.0) go to 10
c check if number entered
      id = 0
c no value found, display variable help
      if ((nv+2).gt.lens) go to 130
      i = nv + 1
      isign = 1
      num = 0
  110 i = i + 1
      if (i.gt.lens) go to 120
      iv = ichar(input(i:i))
c ignore blanks and + sign
      if ((iv.eq.ichar(' ')).or.(iv.eq.ichar('+'))) go to 110
c reverse sign
      if (iv.eq.ichar('-')) then
         isign = -isign 
         go to 110
      endif
      iv = iv - ichar('0')
c check if valid digit
      if ((iv.ge.0).and.(iv.le.9)) then
         id = id + 1
         num = iv + 10*num
         go to 110
      endif
  120 num = isign*num
c if no value found, display variable help
  130 if (id.eq.0) then 
c display variable help
         ay = ay - chh2
         if (ay.lt.0.) go to 10
c draw text
         call gtx(ax,ay,helpv(nvar))
         go to 20
      else
c store value and read next input
         nvalue(nvar) = num
         go to 10
      endif
      end
      subroutine readmc(irc)
c this subroutine reads mouse position and outputs an integer code irc
c which gives its location on the screen, with one exception. 
c the screen is divided into n x n subregions, where n*n is the largest
c integer >= nplot, numbered successively left to right, top to bottom.
c the code returned is: irc = 16 + iplot, where 0 <= iplot <= n*n - 1
c if irc > 127, then irc will be set to 127.
c if the position chosen was inside a special character echo region,
c which is a small rectangle outlined in the lower right hand corner,
c then a help menu is displayed with various choices which return codes:
c quit (irc=1), save (irc=2), modify (irc=3), animate (irc=4), 
c use keys (irc=5), reset(irc=6), and plotparms(irc=7). 
c a non-negative number num can be entered, digit by digit, terminated
c by selecting OK (or CLEAR to try again), and returns irc = num + 128.
c if an area outside the menu box was chosen, then irc=0.
c locator device can be either in request (wait) or event (nowait) mode.
c to go from request to event mode, choose animate.
c to go from event to request mode, click anywhere.
c if no locator event occurred, irc = 0
c if irc = -1 on entry, locator device is forced into request mode
c and help is displayed
c if irc = -2 on entry, locator device is forced into event mode
c if irc = -3 on entry, locator device is forced into request mode
c and readmc returns immediately
c on first entry, locator is assumed to be in request mode
c idwk = workstation identifier
c ncols = number of foreground colors available for line plotting
c rx, ry = ndc coordinates of upper-right corner of workstation window
c nplot = number of plots per page
c idloc = locator device number, 0 if no locator device available
c ifrg = index of foreground color
c kprime = table of color indices for prime colors
      parameter (nchmax=8)
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      character*10 chr(nchmax)
c x,y = scratch variables for plotting
c window = window of normalization transformation returned by locator
c viewpt = viewport of normalization transformation returned by locator
      dimension x(5), y(5), window(4), viewpt(4)
c mode = mode of locator device (0=request,1=sample,2=event)
c iesw = echo switch (0=noecho,1=echo)
      data mode,iesw /0,1/
c ecx,one = minimum/maximum x value of special character echo area
c zero,ecy = minimum/maximum y value of special character echo area
c tmax = upper extent of menu box
c tout = timeout period, in seconds
      data zero,ecy,ecx,one,tmax,tout /0.,.05,.95,1.,.975,0./
c chr = array of character variables containing choices
      data chr /' quit     ',' save     ',' modify   ',' animate  ',' us
     1e keys ',' reset    ',' plotparms',' next     '/
      save mode,iesw,zero,ecy,ecx,one,tmax,tout,chr
c if irc = -1 on input, put locator into request mode
      if ((idloc.gt.0).and.(irc.eq.(-1))) then
         mode = 0
c set locator mode, 0 = request
         call gslcm(idwk,idloc,mode,iesw)
c if irc = -2 on input, put locator into event mode
      elseif ((idloc.gt.0).and.(irc.eq.(-2))) then
         mode = 2
c set locator mode, 2 = event
         call gslcm(idwk,idloc,mode,iesw)
c if irc = -3 on input, put locator into request mode
      elseif ((idloc.gt.0).and.(irc.eq.(-3))) then
         mode = 0
c set locator mode, 0 = request
         call gslcm(idwk,idloc,mode,iesw)
         return
c clear return code
      else
         irc = 0
      endif
c define and select identity transformation
      call idtran
c set linewidth scale factor, 1.0 = nominal
      call gslwsc(1.0)
      if (ncols.ge.3) then
c set fill area color index, yellow
         call gsfaci(kprime(5))
      else
c set fill area color index to foreground
         call gsfaci(ifrg)
      endif
c  set fill area style, 1 = solid
      call gsfais(1)
      if (idloc.gt.0) then
c draw box around special character echo area
         x(1) = ecx*rx
         x(2) = one*rx
         x(3) = one*rx
         x(4) = ecx*rx
         x(5) = ecx*rx
         y(1) = zero*ry
         y(2) = zero*ry
         y(3) = ecy*ry
         y(4) = ecy*ry
         y(5) = zero*ry
c fill area
         call gfa(5,x,y)
c update workstation, perform
         call guwk(idwk,1)
c calculate screen partition
         npl1 = sqrt(real(nplot-1)) + 0.0001
         npl = npl1 + 1
c if locator device is in event mode, check for events
         if (mode.eq.2) then
c await event
   10       call gwait(tout,idwks,icl,idnr)
c return if no event available
            if (icl.eq.0) return
c check again if event is not locator
            if (icl.ne.1) go to 10
c check again if workstation or device id is not correct
            if ((idwks.ne.idwk).or.(idnr.ne.idloc)) go to 10
            px = 0.
            py = 0.
            mode = 0
c set locator mode, 0 = request
            call gslcm(idwk,idloc,mode,iesw)
c get locator
            call ggtlc(nrt,px,py)
            irc = -2
         endif
c if locator device is in request mode, wait for input
         if (mode.eq.0) then
c check if help should be displayed
            if (irc.eq.(-1)) then
c set locator to display help
               px = .5*(ecx + one)
               py = .5*(ecy + zero)
               irc = 0
c find locator position
            else
c normal request
               if (irc.eq.0) then
                  px = 0.
                  py = 0.
c request locator
                  call grqlc(idwk,idloc,istat,nrt,px,py)
               endif
c inquire normalization transformation
               call gqnt(nrt,ierr,window,viewpt)
c find subregion location
               px = (viewpt(1) + (px - window(1))*(viewpt(2) - viewpt(1)
     1)/(window(2) - window(1)))/rx - 0.0001
               py = (viewpt(3) + (py - window(3))*(viewpt(4) - viewpt(3)
     1)/(window(4) - window(3)))/ry + 0.0001
               ix = abs(px)*real(npl)
               iy = abs(1. - py)*real(npl)
c set return code
               irc = 16 + (ix + iy*npl)
               if (irc.gt.127) irc = 127
            endif
c clear plotparm return code
   20       irp = 0
c check if special character echo area was chosen
            if ((px.ge.ecx).and.(px.le.one).and.(py.ge.zero).and.(py.le.
     1ecy)) then
               orgy = tmax*ry
               call amenu (chr,zero,orgy,nchmax,irc)
c clear workstation, always
               call gclrwk(idwk,1)
c if animate was chosen, set locator to event mode
               if (irc.eq.4) then
                  mode = 2
c set locator mode, 2 = event
                  call gslcm(idwk,idloc,mode,iesw)
c change plot parameters
               elseif (irc.eq.7) then
c change plotting parameters
                  call plparmm(zero,ry)
c set plotparm return code
                  irp = 1
               endif
            endif
c display help after changing plot parameters
            if (irp.eq.1) go to 20
         endif
      endif
      return
      end
      subroutine amenu (chr,orgx,orgy,nchmax,irc)
c this subroutine displays menu items for locator device to choose from.
c the input choices are in an array of nchmax character variables, each
c 10 characters long.  the number of choices actually presented is
c limited by the size of the display screen to about 14.  if the locator
c chooses a location inside the menu display, then a return code irc is
c returned corresponding to the item in the character choice list chosen
c if enough space is available, a numerical panel is also displayed.
c a non-negative number num can be then be entered, digit by digit,
c terminated by selecting OK (or CLEAR to try again), and a code
c irc = num + 128 is returned.
c choosing a location outside the menu display returns code irc = 0.
c locator must be in request mode before this subroutine is called.
c chr = array of character variables containing choices
c orgx/orgy = coordinates of upper left hand corner of menu boxes,
c in the range, 0 < orgx < rx and 0 < orgy < ry
c nchmax = dimension of chr, the number of choices possible
c irc = return code, 0 <= 1 <= nchmax
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c idloc = locator device number, 0 if no locator device available
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      character*10 chr(nchmax)
      character*8 cnum
c x,y = scratch variables for plotting
      dimension x(5), y(5)
c slen = width of menu box
c csize = height of characters, also determines height of menu box
      data slen,csize /.25,.033/
c exit if no locator device present
      if (idloc.eq.0) return
c normalize character sizes
      chh = csize*ry
      chh2 = 2.*chh
c calculate extent of first menu box
   10 smin = orgx
      smax = orgx + slen*rx
      tmax = orgy
      tmin = tmax - chh2
c find maximum possible number of menu items which can fit in display
      mmax = tmax/chh2
      nmax = min0(mmax,nchmax)
c calculate location for printing first choice
      dxa = .01*(smax - smin)
      ax = smin + dxa
      ay = tmin + .5*chh
c store coordinates of menu box
      x(1) = smin
      x(2) = smin
      x(3) = smax
      x(4) = smax
      x(5) = smin
      y(1) = tmax
      y(2) = tmin
      y(3) = tmin
      y(4) = tmax
      y(5) = tmax
c find bottom of last menu box
      tmin = tmax - real(nmax)*chh2
c clear workstation, always
      call gclrwk(idwk,1)
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps (chh,1,0)
c draw line across top of menu items
      call gpl(2,x(4),y(4))
c draw box around menu items and print possible choices
      do 30 j = 1, nmax
c draw polyline
      call gpl(4,x,y)
c draw text
      call gtx(ax,ay,chr(j))
c reset values for next box
      do 20 i = 1, 4
      y(i) = y(i) - chh2
   20 continue
      ay = ay - chh2
   30 continue
c if there is enough room, display number, control,and echo panels.
      if ((mmax-nmax).ge.3) then
c draw box for number panel
         nflg = 1
c draw polyline
         call gpl(4,x,y)
c draw zero
         is = ichar('0')
c draw text
         call gtx(ax,ay,char(is))
c draw remaining numbers
         dxn = .1*(smax - smin)
         do 40 i = 1, 9
         x(1) = smin + dxn*real(i)
         x(2) = x(1)
c draw polyline
         call gpl(2,x,y)
c draw text
         x(1) = x(1) + dxa
         call gtx(x(1),ay,char(i+is))
   40    continue
c draw box for control panel
         x(1) = smin
         x(2) = smin
         do 50 i = 1, 4
         y(i) = y(i) - chh2
   50    continue
         ay = ay - chh2
c draw polyline
         call gpl(4,x,y)
c draw left side
c draw text
         call gtx(ax,ay,'CLEAR')
c draw right side
         dxc = .5*(smax - smin)
         x(1) = smin + dxc
         x(2) = x(1)
c draw polyline
         call gpl(2,x,y)
c draw text
         x(1) = x(1) + dxa
         call gtx(x(1),ay,'  OK ')
c draw echo panel
         x(1) = smin
         x(2) = smin
         do 60 i = 1, 4
         y(i) = y(i) - chh2
   60    continue
         ay = ay - chh2
c draw polyline
         call gpl(4,x,y)
         tmin = tmin - 2.*chh2
c initial values for displaying numbers
         num = 0
         id = 0
         cnum = '        '
      else
         nflg = 0
      endif
c update workstation, perform
      call guwk(idwk,1)
c wait for input
   70 call fndloc(px,py)
c check if choice was inside menu region and if so get item number
      if ((px.ge.smin).and.(px.le.smax).and.(py.ge.tmin).and.(py.le.tmax
     1)) then
         iy = (py - tmin)/chh2
         irc = nmax - iy
c check if numbers are included
         if (nflg.gt.0) then
            iy = iy - 2
c menu item selected
            if (iy.ge.0) then
               irc = nmax - iy
               return
            endif
c number selected
            if (iy.eq.(-1)) then
               ix = (px - smin)/dxn
               if (ix.gt.9) ix = 9
               id = id + 1
               if (id.gt.7) go to 70
               num = ix + 10*num 
c echo number
               cnum(id:id) = char(ix+is)
c draw text
               call gtx(ax,ay,cnum(1:id))
c update workstation, perform
               call guwk(idwk,1)
c get next digit
               go to 70
c control panel selected
            elseif (iy.eq.(-2)) then
               ix = (px - smin)/dxc
c clear: reset initial values for displaying numbers
               if (ix.eq.0) then
                  num = 0
                  id = 0
                  cnum = '        '
                  go to 10
c number accepted, if entered:
               else
                  if (id.eq.0) then
                     irc = 0
                  else
                     irc = 128 + num
                  endif
                  return
               endif
            endif
         endif
      else
         irc = 0
      endif
      return
      end
      subroutine plparmm(orgx,orgy)
c this subroutine allows one to make changes to certain of the plotting
c parameters stored in the common blocks plotcm and pextra.
c currently, only the variables nplot and ndi can be changed.
c variables are changed by selecting variable name, then numbers,
c then selecting OK or CLEAR to accept or reject the new number.  If the
c new value is acceptable, select run.  to abort changes, select quit.
c on entry, locator device is assumed to be in request mode
c orgx/orgy = coordinates of upper left hand corner of writing area,
c in the range, 0 < orgx < rx and 0 < orgy < ry
c nvars = number of variables to be modified
      parameter (nvars=2)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c nplot = number of plots per page
c idloc = locator device number, 0 if no locator device available
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c ndi = increment between frames
      common /pextra/ ndi
      character*16 output
      character*44 prompt
      character*16 cnum
      dimension nvalue(nvars)
      character*8 code(nvars)
      character*40 helpv(nvars)
      character*44 helpr(10)
c x,y = scratch variables for plotting
      dimension x(5), y(5)
      save prompt,code,helpv,helpr
   91 format (1x,a6,' = ',i4,1x)
   92 format (1x,a6,' = ')
      data prompt /'select variable name and numbers or command '/
      data code /'NPLOT   ','NDI     '/
      data helpv  /' NPLOT = number of plots per page       ',' NDI = in
     1crement between frames         '/
      data helpr /'    this program interactively modifies and ','update
     1s plot parameters.                    ',' to change parameters, se
     2lect variable name','then numbers, then select OK or CLEAR to    '
     3,'accept or reject the new number.            ','    to obtain the
     4 meaning of a parameter,   ','select its name then OK, or select i
     5t twice.','to display current parameters, select CLEAR ','to abort
     6 modification process, select quit  ',' when parameters are OK, se
     7lect run.        '/
c slen = width of menu box
c csize = height of characters, also determines height of menu box
      data slen,csize /.375,.033/
c exit if no locator device present
      if (idloc.eq.0) return
c first copy variables
      nvalue(1) = nplot
      nvalue(2) = ndi
c normalize character sizes and box width
      chh = csize*ry
      chh2 = 2.*chh
      dxb = slen*rx
      dxn = dxb/15.
      dxa = .1*dxn
c initial values for display numbers
      nvar = 0
      num = 0
      id = 0
      cnum = '                '
c calculate extent of first menu box
   10 smin = orgx
      smax = orgx + dxb
      tmax = orgy
      tmin = tmax - chh2
c calculate location for printing first choice
      ax = orgx + dxa
      ay = tmin + .5*chh
c store coordinates of menu box
      x(1) = smin
      x(2) = smin
      x(3) = smax
      x(4) = smax
      x(5) = smin
      y(1) = tmax
      y(2) = tmin
      y(3) = tmin
      y(4) = tmax
      y(5) = tmax
c clear workstation, always
      call gclrwk(idwk,1)
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps (chh,1,0)
c display variable names and values
c draw line to left of menu items
      call gpl(2,x,y)
c draw box around menu items and print possible choices
      do 30 j = 1, nvars
c draw polyline
      call gpl(4,x(2),y(2))
c write out current values
      write (output,91) code(j), nvalue(j)
c draw text
      call gtx(ax,ay,output(1:15))
c reset values for next box
      do 20 i = 2, 5
      x(i) = x(i) + dxb
   20 continue
      ax = ax + dxb
   30 continue
c draw boxes for numbers
      x(1) = smin
      x(2) = smin
      x(3) = smax
      x(4) = smax
      do 40 i = 1, 4
      y(i) = y(i) - chh2
   40 continue
      ax = smin + dxa
      ay = ay - chh2
c draw polyline
      call gpl(4,x,y)
c draw zero
      is = ichar('0') - 1
c draw remaining numbers
      do 50 i = 1, 10
      x(1) = smin + dxn*real(i)
      x(2) = x(1)
c draw polyline
      call gpl(2,x,y)
c draw text
      call gtx(ax,ay,char(i+is))
      ax = ax + dxn
   50 continue
      is = is + 1
c draw text
      call gtx(ax,ay,'CLEAR')
c draw boxes for commands
      x(2) = smax
      x(3) = smax + dxb
      x(4) = x(3)
c draw polyline
      call gpl(3,x(2),y(2))
c draw OK command
      ax = smax + 5.*dxa
c draw text
      call gtx(ax,ay, 'OK')
c draw help command
      x(1) = smax + 3.*dxn
      x(2) = x(1)
c draw polyline
      call gpl(2,x,y)
c draw text
      x(1) = x(1) + 5.*dxa
      call gtx(x(1),ay,'?')
c draw run command
      dxc = 5.*dxn
      x(1) = smax + dxc
      x(2) = x(1)
c draw polyline
      call gpl(2,x,y)
c draw text
      call gtx(x(1),ay,' run ')
c draw quit command
      x(1) = smax + 2.*dxc
      x(2) = x(1)
c draw polyline
      call gpl(2,x,y)
c draw text
      x(1) = x(1) + 5.*dxa
      call gtx(x(1),ay,'quit')
c draw echo panel
      x(1) = smin
      x(2) = smin
      x(3) = smax
      x(4) = smax
      do 60 i = 1, 4
      y(i) = y(i) - chh2
   60 continue
      ax = smin + dxa
      ay = ay - chh2
c draw polyline
      call gpl(4,x,y)
c draw text
      if (nvar.gt.0) call gtx(ax,ay,cnum)
c find menu selection area
      smax = smin + 2.*dxb
      tmin = tmax - 2.*chh2
      ayy = ay
c display prompt
   70 ax = orgx + dxa
      ayy = ayy - chh2
      if (ayy.lt.0.) go to 10
c draw text
      call gtx(ax,ayy,prompt)
c update workstation, perform
      call guwk(idwk,1)
c wait for input
   80 call fndloc(px,py)
c check if choice was inside menu region
      if ((px.ge.smin).and.(px.le.smax).and.(py.ge.tmin).and.(py.le.tmax
     1)) then
         ix = (px - smin)/dxb
         iy = (py - tmin)/chh2
c number, help, or command selected
         if (iy.eq.0) then
c help or command selected
            if (ix.eq.1) then
               ix = (px - (smin + dxb))/dxc
c quit selected
               if (ix.eq.2) return
c help or OK selected
               if (ix.eq.0) then
                  ix = (px - (smin + dxb))/(3.*dxn)
c OK selected
                  if (ix.eq.0) then
                     if (nvar.eq.0) go to 10
c display variable help
                     if (id.eq.0) go to 100
c store value
                     nvalue(nvar) = num
c reset variable name
                     nvar = 0
                     num = 0
                     id = 0
                     cnum = '                '
                     go to 10
c write help instructions
                  elseif (ix.eq.1) then
c clear workstation, always
                     call gclrwk(idwk,1)
c set default plotting parameters
                     ay = orgy
c 1 = default font, 0 = string precision
                     call dfplps (chh,1,0)
                     do 90 i = 1, 10
                     ay = ay - chh2
c draw text
                     call gtx(ax,ay,helpr(i))
   90                continue
                     ay = ay - chh2
c draw text
                     if (nvar.gt.0) call gtx(ax,ay,helpv(nvar))
c update workstation, perform
                     call guwk(idwk,1)
c request locator
                     call grqlc(idwk,idloc,istat,nrt,px,py)
                     go to 10
                  endif
c run selected
               elseif (ix.eq.1) then
c update pending variable
                  if ((nvar.gt.0).and.(id.gt.0)) nvalue(nvar) = num
c update variables
                  nplot = nvalue(1)
                  ndi = nvalue(2)
                  return
               endif
c number or CLEAR selected
            elseif (ix.eq.0) then
               if (nvar.eq.0) go to 10
c number selected
               ix = (px - smin)/dxn
               if (ix.lt.10) then
                  id = id + 1
                  if (id.gt.4) go to 80
                  num = ix + 10*num 
c echo number
                  cnum(10+id:10+id) = char(ix+is)
c draw text
                  call gtx(ax,ay,cnum(1:10+id))
c update workstation, perform
                  call guwk(idwk,1)
c get next digit
                  go to 80
c clear selected
               else 
                  nvar = 0
                  num = 0
                  id = 0
                  cnum = '                '
                  go to 10
               endif
            endif
c variable selected
         elseif (iy.eq.1) then
            if ((ix+1).eq.nvar) go to 100
            nvar = ix + 1
            num = 0
            id = 0
            cnum = '                '
            write (cnum,92) code(nvar)
            go to 10
         endif
c if choice is not inside menu region, update display
      else
         if (nvar.eq.0) go to 10
         go to 100
      endif
c display variable help
  100 ayy = ayy - chh2
      if (ayy.lt.0.) go to 10
c draw text
      call gtx(ax,ayy,helpv(nvar))
      go to 70
      end
      subroutine fndloc(px,py)
c this subroutine finds position of locator
c px,py = locator position
c locator is assumed in request mode, with identity transformation
c idwk = workstation identifier
c idloc = locator device number, 0 if no locator device available
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c window = window of normalization transformation returned by locator
c viewpt = viewport of normalization transformation returned by locator
      dimension window(4), viewpt(4)
      px = 0.
      py = 0.
c request locator
      call grqlc(idwk,idloc,istat,nrt,px,py)
c inquire normalization transformation
      call gqnt(nrt,ierr,window,viewpt)
c find location
      px = viewpt(1) + (px - window(1))*(viewpt(2) - viewpt(1))/(window(
     12) - window(1))
      py = viewpt(3) + (py - window(3))*(viewpt(4) - viewpt(3))/(window(
     14) - window(3))
      return
      end
      subroutine verstc (ia,ic,nc,irc)
c this program decodes tektronix 4014 information and plots result
c written for the ibm rs/6000 - viktor k. decyk, ucla
c ia = input array of ascii characters
c ic = number of entries in ia array to decode, ic = 0 means last plot
c nc = dimension of ia array
c irc = return code (0 = normal return)
c idwk = workstation identifier
c ncols = number of foreground colors available for line plotting
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c ifrg = index of foreground color
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c iebc = (0,1) = (no,yes) output characters are in ebcdic
      common /march/ iebc, irvb, longi
c ndi = increment between frames
      common /pextra/ ndi
      character*1 lbl
      dimension ia(nc)
      dimension lt(8), iate(128)
c x,y = scratch variables for plotting
      dimension x(2), y(2)
      save xlen,ylen,chxp
      save ig,nf,nd,ip,ird
      save if,ls,im,ie,iyf,nrt,irl,at
      save icflg,lt,iate
      save cx,cy,hy,by,xmin,xmax,ymin,ymax,alx,chl
      save yh,yl,xh,xl,aix,aiy
c  91 format (1x,i6,17h frame(s) plotted)
c xlen,ylen = width and height of tektronix screen
c chxp = character expansion factor
      data xlen,ylen,chxp /1023.,780,0.8/
c ig = graphics mode flag, 0 = alpha mode, >0 = graphics mode:
c 1=draw,2=dark,3=point plot,4=special point plot,5=incremental plot 
c nf = current frame number being processed
c nd = current frame number being displayed 
c ip = printer flag = (0,1) = (no,yes) printer is attached
c ird = return entry point flag:
c 1 = new data block,2 = escape sequence,3 = special point plot
c 4 = incremental plot,5 = decode graphics address
      data ig,nf,nd,ip,ird /0,0,0,0,0/
c if  = focus flag = (0,1) = spot is (focused,defocused)
c ls = line style flag,0=solid,1=dot,2=dash-dot,3=short-dash,4=long-dash
c im = beam flag = (0,1) = beam in incremental plot mode is (off,on)
c iyf = Lo Y byte read flag, used in decoding graphics address
c nrt = transformation number
c irl = escape sequence return entry point
c at = intensity parameter, (0. < at < 1.)
c alx = x coordinate for numerical labels, 10 characters from right edge
      data if,ls,im,ie,iyf,nrt,alx,at /0,0,0,0,0,1,882.,1./
c special flags
c icflg = (0,1,2) = use (line style,color,both) for line style code
      data zero,icflg /0.,1/
c line style table to convert from tektronix to gks codes
      data lt /1,3,4,2,2,1,1,1/
c ascii to ebcdic conversion table
c ebcdic code for ascii 124 is non-standard
      data iate /0,1,2,3,55,45,46,47,22,5,37,11,12,13,14,15,16,17,18,19,
     160,61,50,38,24,25,63,39,34,29,53,31,64,90,127,123,91,108,80,125,77
     2,93,92,78,107,96,75,97,240,241,242,243,244,245,246,247,248,249,122
     3,94,76,126,110,111,124,193,194,195,196,197,198,199,200,201,209,210
     4,211,212,213,214,215,216,217,226,227,228,229,230,231,232,233,173,2
     524,189,95,109,121,129,130,131,132,133,134,135,136,137,145,146,147,
     6148,149,150,151,152,153,162,163,164,165,166,167,168,169,192,79,208
     7,161,7/
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
c check if last plot
      if (ic.lt.1) go to 360
      k = 1
c check if first plot
      if (ird.eq.0) go to 350
c go to return entry point, continue where left off 
      go to (10,70,150,190,270), ird
c new data block, read character
   10 ird = 1
      if (k.gt.ic) return
      lc = ia(k)
      k = k + 1
c control character
   20 if (lc.lt.32) go to 30
c decode graphics address
      if (ig.gt.0) go to 250
c alpha mode
      if (iebc.eq.1) lbl = char(iate(lc+1))
      if (iebc.eq.0) lbl = char(lc)
c plot character
      if (nf.eq.nd) call gtx (aix,aiy,lbl)
      go to 200
c control characters
c set graph mode (dark vector)
   30 if (lc.eq.29) go to 40
c carriage return
      if (lc.eq.13) go to 50
c set alpha mode
      if (lc.eq.31) go to 60
c escape sequence
      if (lc.eq.27) go to 70
c point plotting mode
      if (lc.eq.28) go to 180
c incremental plot mode
      if (lc.eq.30) go to 190
c tab
      if (lc.eq.9) go to 200
c line feed
      if (lc.eq.10) go to 210
c backspace
      if (lc.eq.8) go to 220
c vertical tab
      if (lc.eq.11) go to 230
c bell (set vector to draw)
      if (lc.eq.7) go to 240
c unknown control character, read next data block
      go to 10
c set graph mode (dark vector)
   40 ig = 2
c read next data block
      go to 10
c carriage return
   50 aix = zero
      ie = 0
c set alpha mode
   60 ig = 0
      ls = 0
      if = 0
c set linetype, 1 = solid
      if (icflg.ne.1) call gsln(1)
c set polyline color index, 1 = foreground
      if (icflg.ne.0) call gsplci(ifrg)
c set linewidth scale factor, 1.0 = nominal
      call gslwsc(1.0)
c read next data block
      go to 10
c escape sequence
   70 if (ird.ne.2) irl = ird
      ird = 2
      if (k.gt.ic) return
      lc = ia(k)
      k = k + 1
c new frame (erase)
      if (lc.eq.12) go to 80
c make hardcopy
      if (lc.eq.23) go to 90
c large characters
      if (lc.eq.56) go to 100
c medium-large characters
      if (lc.eq.57) go to 110
c medium-small characters
      if (lc.eq.58) go to 120
c small characters
      if (lc.eq.59) go to 130
c set line style and focus
      if ((lc.ge.96).and.(lc.le.119)) go to 140
c special point plot mode
      if (lc.eq.28) go to 150
c special code for del
      if (lc.eq.63) go to 160
c other control character
      if ((lc.eq.7).or.(lc.eq.8).or.(lc.eq.9).or.(lc.eq.11).or.(lc.eq.29
     1).or.(lc.eq.30).or.(lc.eq.31)) go to 30
c other escape sequence
      if ((lc.eq.0).or.(lc.eq.10).or.(lc.eq.13).or.(lc.eq.27).or.(lc.eq.
     1127)) go to 70
c unsupported escape sequence, read next data block
      go to 10
c new frame
   80 nf = nf + 1
c display frame
      if (nf.gt.nd) then
c add label
         call labelf(nf,alx,zero,chl)
         nd = nd + ndi
c update plot number
         iplot = iplot + 1
         if (iplot.eq.nplot) iplot = 0
         if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
            call guwk(idwk,1)
c read code from input device, if present
            call readrc(irc)
         endif
c check return codes
c quit
         if (irc.eq.1) then
c           write (6,91) nf
            return
         endif
c reset
         if (irc.eq.6) nd = 0
c new display number entered
         if (irc.ge.128) then
            nd = irc - 128
            if (nd.gt.0) nd = nd - 1
            irc = 0
         endif
c reset, if necessary
         if (nf.gt.nd) then
            nf = 0
            ird = 1
            irc = 6
         endif
      endif
c next frame will be displayed
      if (nf.eq.nd) then
c initiate plot
         if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
            call gclrwk(idwk,1)
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
         aplx = orx + aplx
         aply = ory + aply
c define transformation for tektronix screen
c set window
         call gswn (nrt,xmin,xmax,ymin,ymax)
c set viewport
         call gsvp (nrt,orx,aplx,ory,aply)
c select normalization transformation
         call gselnt(nrt)
c set default plotting parameters
         chh = .7*cy
         chl = cy
c 1 = default font, 2 = stroke precision
         call dfplps (chh,1,2)
      endif
c reset plotting parameters
      if (nf.le.nd) then
         ig = 0
         ls = 0
         if = 0
         aix = zero
         aiy = hy
         ie = 0
      endif
c reset if necessary
      if (irc.eq.6) return
c read next data block
      go to 10
c make hardcopy
   90 if (ip.gt.0) go to 170
      go to 170
c large characters
  100 cx = 14.
      cy = 22.
      call gschh(.7*cy)
      go to 170
c medium-large characters
  110 cx = 13.
      cy = 21.
      call gschh(.7*cy)
      go to 170
c medium-small characters
  120 cx = 9.
      cy = 13.
      call gschh(.7*cy)
      go to 170
c small characters
  130 cx = 8.
      cy = 12.
      call gschh(.7*cy)
      go to 170
c set line style and focus
  140 ls = lc - 96
      if = ls/8
      ls = ls - if*8
      it = lt(ls+1)
c set linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      if (icflg.ne.1) call gsln(it)
      if (icflg.ne.0) then
c        it = ls + 1
         if (ls.gt.3) it = ls + 1
c cycle through available colors
         it = it - ncols*((it - 1)/ncols)
         it = kprime(it+1)
c set polyline color index,
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gsplci(it)
      endif
c if = (0,1) = spot is (focused,defocused)
      alw = real(if + 1)
      if (if.eq.2) alw = 1.
c set linewidth scale factor, 1.0 = nominal
      call gslwsc(alw)
      go to 170
c special point plot mode
  150 ird = 3
      if (k.gt.ic) return
      lc = ia(k)
      k = k + 1
      ig = 4
      if ((lc.lt.32).or.(lc.gt.125)) go to 150
      if = 1 - lc/64
      if (lc.gt.65) lc = lc - 60
      if (lc.gt.56) lc = lc - 9
c intensity parameter, (0. < at < 1.), not implemented
      at = real(lc)
      at = .015*10.**(at/31.)
c if = (0,1) = spot is (focused,defocused)
      alw = real(if + 1)
c set linewidth scale factor, 1.0 = nominal
      call gslwsc(alw)
c read next data block
      go to 10
c special code for del
  160 if ((ig.eq.0).or.(ig.eq.5)) go to 170
      lc = 127
c Extra Byte
      if (iyf.eq.1) go to 280
c Lo Y
      go to 260
c go to return entry point, continue where left off 
  170 go to (10,70,150,190,270), irl
c point plotting mode
  180 ig = 3
c read next data block
      go to 10
c incremental plot mode
  190 ird = 4
      if (k.gt.ic) return
      lc = ia(k)
      k = k + 1
      ig = 5
c alpha mode
      if ((lc.eq.13).or.(lc.eq.31)) go to 30
c escape sequence
      if (lc.eq.27) go to 70
      if (lc.eq.126) lc = 127
      lc = lc - 32
      is = lc/16
      lc = lc - is*16
c beam off
      if (is.ne.2) im = 0
c beam on
      if (is.eq.3) im = 1
      dx = 0.
      dy = 0.
      it = lc/4
c increment north
      if (it.eq.1) dy = .25
c increment south
      if (it.eq.2) dy = -.25
      it = lc - 4*it
c increment east
      if (it.eq.1) dx = .25
c increment west
      if (it.eq.2) dx = -.25
      ai = aix + dx
      aj = aiy + dy
      it = is - (is/2)*2
      if ((is.eq.2).and.(im.eq.1)) it = 1
c draw increment
      if (it.eq.1) go to 320
c perform move
      go to 330
c tab
  200 aix = aix + cx
c read next data block
      if (aix.lt.xlen) go to 10
      aix = aix - xlen
c line feed
  210 aiy = aiy - cy
      if (aiy.lt.zero) aiy = hy
c read next data block
      go to 10
c backspace
  220 aix = aix - cx
c read next data block
      if (aix.ge.zero) go to 10
      aix = aix + xlen
c vertical tab
  230 aiy = aiy + cy
      if (aiy.gt.ylen) aiy = by
c read next data block
      go to 10
c bell (set vector to draw)
  240 if (ig.eq.2) ig = 1
c read next data block
      go to 10
c graph mode
c decode graphics address
c in the following description of the decoding scheme, bits 6-7
c represent tag bits.  the address comes in 5 lowest order bits of each
c of four words, some of which may be left out:
c Hi Y, Lo Y, Hi X, Lo X.  Sometimes an Extra Byte may be used, in which
c the lowest order 4 bits are used to give added precision in X and Y,
c the 2 lowest order bits for X and the next 2 higher order bits for Y.
c if tag=32, lc must be Hi Y address, unless iyf flag is set (when Lo Y 
c    was previously read), in which case lc must be Hi X.
c if tag=64, it must be Lo X and the last address.
c if tag=96, it must be Lo Y, unless iyf flag is set, in which case, the
c    preceeding Lo Y was actually an Extra Byte, and current lc is Lo Y.
c Hi Y
  250 if (lc.lt.64) go to 300
c Lo X
      if (lc.lt.96) go to 310
c Lo Y
  260 yl = real(lc)
      iyf = 1
c read additional address bytes
  270 ird = 5
      if (k.gt.ic) return
      lc = ia(k)
      k = k + 1
c carriage return
      if (lc.eq.13) go to 50
c unknown control character, try again
      if (lc.lt.27) go to 270
c escape sequence
      if (lc.eq.27) go to 70
c set new mode, read next data block
      if (lc.lt.32) go to 10
c Hi X
      if (lc.lt.64) go to 290
c Lo X
      if (lc.lt.96) go to 310
c Lo Y
      if (iyf.eq.0) go to 260
c Extra Byte
  280 ie = yl
      ie = ie - (ie/16)*16
      go to 260
c Hi X
  290 if (iyf.eq.0) go to 300
      xh = real(lc)
      go to 270
c Hi Y
  300 yh = real(lc)
      iyf = 0
      go to 270
c Lo X
  310 xl = real(lc)
c calculate address
      ai = (xh - 34.)*32. + xl
      aj = (yh - 35.)*32. + yl
      it = ie/4
      ai = ai + .25*real(ie - 4*it)
      aj = aj + .25*real(it)
      ai = ai
      aj = aj
c perform clipping
      if (ai.lt.zero) ai = zero
      if (ai.gt.xlen) ai = xlen
      if (aj.lt.zero) aj = zero
      if (aj.gt.ylen) aj = ylen
c perform move
      if (ig.ne.1) go to 330
c draw vector
  320 x(1) = aix
      x(2) = ai
      y(1) = aiy
      y(2) = aj
      if (nf.eq.nd) call gpl(2,x,y)
      aix = ai
      aiy = aj
c read next data block
      if (ig.lt.4) go to 10
c special point plot mode
      if (ig.eq.4) go to 150
c incremental plot mode
      if (ig.eq.5) go to 190
c perform move
  330 aix = ai
      aiy = aj
c set vector to draw
      if (ig.eq.2) go to 340
c incremental plot mode
      if (ig.eq.5) go to 190
c draw vector
      go to 320
c set vector to draw
  340 ig = 1
c read next data block
      go to 10
c set scales for first plot
c default character sizes
  350 cx = 14.
      cy = 22.
c home location
      hy = 767.
c bottom location
      by = zero
c default cursor location
      aix = zero
      aiy = hy
c screen size, enlarge by one character size
      xmin = zero - cx
      xmax = xlen + cx
      ymin = zero - cy
      ymax = ylen + cy
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
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
      aplx = orx + aplx
      aply = ory + aply
c define transformation for tektronix screen
c set window
      call gswn (nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp (nrt,orx,aplx,ory,aply)
c select normalization transformation
      call gselnt(nrt)
c for monochrome device, use line styles
      if (ncols.eq.1) icflg = 0
c set default plotting parameters
      chh = .7*cy
      chl = cy
c 1 = default font, 2 = stroke precision
      call dfplps (chh,1,2)
c set character expansion factor
      call gschxp(chxp)
c read next data block
      go to 10
c last plot
  360 nf = nf + 1
      if (nf.le.nd) then
c initiate plot
         if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
            call gclrwk(idwk,1)
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
         aplx = orx + aplx
         aply = ory + aply
c define transformation for tektronix screen
c set window
         call gswn (nrt,xmin,xmax,ymin,ymax)
c set viewport
         call gsvp (nrt,orx,aplx,ory,aply)
c select normalization transformation
         call gselnt(nrt)
c set default plotting parameters
         chh = .7*cy
c 1 = default font, 2 = stroke precision
         call dfplps (chh,1,2)
      endif
c add label
      call labelf(nf,alx,zero,chl)
      nd = nd + ndi
c reset plot number
      iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         call readrc(irc)
      endif
c check return codes
c force device into request mode
      if (irc.eq.0) then
         irc = -1
c read code from input device, if present
         call readrc(irc)
      endif
c animate
      if (irc.eq.4) irc = 6
c reset
      if (irc.eq.6) nd = 0
c new display number entered
      if (irc.ge.128) then
         nd = irc - 128
         if (nd.gt.0) nd = nd - 1
         irc = 0
      endif
c reset, if necessary
      if (nf.gt.nd) then
         nf = 0
         ird = 1
         irc = 6
      endif
c reset 
      if (irc.eq.6) then
c next frame will be displayed
         if (nf.eq.nd) then
c initiate plot
            if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
               call gclrwk(idwk,1)
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
            aplx = orx + aplx
            aply = ory + aply
c define transformation for tektronix screen
c set window
            call gswn (nrt,xmin,xmax,ymin,ymax)
c set viewport
            call gsvp (nrt,orx,aplx,ory,aply)
c select normalization transformation
            call gselnt(nrt)
c set default plotting parameters
            chh = .7*cy
c 1 = default font, 2 = stroke precision
            call dfplps (chh,1,2)
         endif
c reset plotting parameters
         if (nf.le.nd) then
            ig = 0
            ls = 0
            if = 0
            aix = zero
            aiy = hy
            ie = 0
         endif
c     else
c        write (6,91) nf
      endif
      return
      end
      subroutine dimage(image,lx,ly,lz,lenb,npix,nf,irc)
c this subroutine displays raster image stored in character array
c the data must first be copied to an array of 32 bit integers
c an identity transformation has been assumed
c image = uncompressed single image
c lx, ly = the size of the image, in pixels
c lz = width of picture, in bytes
c lenb = size of picture, in bytes
c npix = number of pixels per byte
c nf = current frame number being processed
c irc = return code (0 = normal return)
c npald = number of palette entries
      parameter(npald=256)
c lxm, lym = maximum number of pixels in x, y
      parameter(lxm=720,lym=540)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      character*1 image(lenb)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img8 = integer image array
      dimension img8(lxm*lym)
      save csize,alx
c csize = vertical size of characters
c alx = x coordinate for numerical labels
      data csize,alx /.02,.86/
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
c nbit = the number of colors, pixel depth
      nbit = 8/npix
c calculate maximum size of image
      lxs = lx
      if (lxs.gt.lxm) lxs = lxm
      lys = ly
      if (lys.gt.lym) lys = lym
c eight bit color
      if (nbit.eq.8) then
c copy to integer array without lookup table
         if (lupt.eq.0) then
            do 20 k = 1, lys
            ioff = lz*(k - 1)
            joff = lxs*(k - 1)
            do 10 j = 1, lxs
            img8(j+joff) = ichar(image(j+ioff))
   10       continue
   20       continue
c copy to integer array with lookup table
         else
            do 40 k = 1, lys
            ioff = lz*(k - 1)
            joff = lxs*(k - 1)
            do 30 j = 1, lxs
            img8(j+joff) = ipal(ichar(image(j+ioff))+1)
   30       continue
   40       continue
         endif
c nbits per pixel
      else
c maximum width
         lzs = lxs/npix
c convert from nbits per pixel to 8 bits per pixel
         ntc = 2**nbit
         npixm = npix - 1
c copy to integer array without lookup table
         if (lupt.eq.0) then
            do 70 k = 1, lys
c loop over bytes in image
            j1 = lxs*(k - 1) + npix + 1
            ioff = lz*(k - 1)
            do 60 j = 1, lzs
            itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
            do 50 i = 1, npixm
            it1 = itc/ntc
            img8(j1-i) = itc-it1*ntc
            itc = it1
   50       continue
            img8(j1-npix) = itc
            j1 = j1 + npix
   60       continue
   70       continue
c copy to integer array with lookup table
         else
            do 100 k = 1, lys
c loop over bytes in image
            j1 = lxs*(k - 1) + npix + 1
            ioff = lz*(k - 1)
            do 90 j = 1, lzs
            itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
            do 80 i = 1, npixm
            it1 = itc/ntc
            img8(j1-i) = ipal(itc-it1*ntc+1)
            itc = it1
   80       continue
            img8(j1-npix) = ipal(itc+1)
            j1 = j1 + npix
   90       continue
  100       continue
         endif
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
      aplx = orx + sx
      aply = ory + sy
c normalize characters and location
      afx = orx + alx*sx
      chl = csize*ry
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c cell array
c special case for rs/6000 with graPHIGS gks
c     call gca(orx,aply,aply,ory,lys,lxs,1,1,lys,lxs,img8)
      call gca(orx,aply,aplx,ory,lxs,lys,1,1,lxs,lys,img8)
c set text font, 1 = default font, 2 = stroke precision
      call gstxfp(1,2)
c add label
      call labelf(nf,afx,ory,chl)
c update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         call readrc(irc)
      endif
      return
      end
      subroutine fimage(image,lx,ly,lz,lenb,npix,nf)
c this subroutine displays raster image stored in character array
c the data must first be copied to an array of 32 bit integers
c an identity transformation has been assumed
c image = uncompressed single image
c lx, ly = the size of the image, in pixels
c lz = width of picture, in bytes
c lenb = size of picture, in bytes
c npix = number of pixels per byte
c nf = current frame number being processed
c npald = number of palette entries
      parameter(npald=256)
c lxm, lym = maximum number of pixels in x, y
      parameter(lxm=720,lym=540)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      character*1 image(lenb)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img8 = integer image array
      dimension img8(lxm*lym)
      save nfl,zero,csize,alx
c nfl = previous frame number
      data nfl /1/
c csize = vertical size of characters
c alx = x coordinate for numerical labels
      data zero,csize,alx /0.,.02,.86/
c nbit = the number of colors, pixel depth
      nbit = 8/npix
c calculate maximum size of image
      lxs = lx
      if (lxs.gt.lxm) lxs = lxm
      lys = ly
      if (lys.gt.lym) lys = lym
c eight bit color
      if (nbit.eq.8) then
c copy to integer array without lookup table
         if (lupt.eq.0) then
            do 20 k = 1, lys
            ioff = lz*(k - 1)
            joff = lxs*(k - 1)
            do 10 j = 1, lxs
            img8(j+joff) = ichar(image(j+ioff))
   10       continue
   20       continue
c copy to integer array with lookup table
         else
            do 40 k = 1, lys
            ioff = lz*(k - 1)
            joff = lxs*(k - 1)
            do 30 j = 1, lxs
            img8(j+joff) = ipal(ichar(image(j+ioff))+1)
   30       continue
   40       continue
         endif
c nbits per pixel
      else
c maximum width
         lzs = lxs/npix
c convert from nbits per pixel to 8 bits per pixel
         ntc = 2**nbit
         npixm = npix - 1
c copy to integer array without lookup table
         if (lupt.eq.0) then
            do 70 k = 1, lys
c loop over bytes in image
            j1 = lxs*(k - 1) + npix + 1
            ioff = lz*(k - 1)
            do 60 j = 1, lzs
            itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
            do 50 i = 1, npixm
            it1 = itc/ntc
            img8(j1-i) = itc-it1*ntc
            itc = it1
   50       continue
            img8(j1-npix) = itc
            j1 = j1 + npix
   60       continue
   70       continue
c copy to integer array with lookup table
         else
            do 100 k = 1, lys
c loop over bytes in image
            j1 = lxs*(k - 1) + npix + 1
            ioff = lz*(k - 1)
            do 90 j = 1, lzs
            itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
            do 80 i = 1, npixm
            it1 = itc/ntc
            img8(j1-i) = ipal(itc-it1*ntc+1)
            itc = it1
   80       continue
            img8(j1-npix) = ipal(itc+1)
            j1 = j1 + npix
   90       continue
  100       continue
         endif
      endif
c calculate range for picture
c x range
      if (isx.ge.lxs) then
         sx = rx*real(lxs)/real(isx)
      else
         sx = rx
      endif
c y range
      if (isy.ge.lys) then
         sy = ry*(1. - real(lys)/real(isy))
      else
         sy = zero
      endif
c normalize characters and location
      afx = alx*sx
      chl = csize*ry
c clear workstation, if frame number is not advanced
      if (nf.le.nfl) call gclrwk(idwk,1)
c cell array
c special case for rs/6000 with graPHIGS gks
c     call gca(zero,ry,sx,sy,lys,lxs,1,1,lys,lxs,img8)
      call gca(zero,ry,sx,sy,lxs,lys,1,1,lxs,lys,img8)
c set text font, 1 = default font, 2 = stroke precision
      call gstxfp(1,2)
c add label
      call labelf(nf,afx,sy,chl)
c save current frame number
      nfl = nf
c update workstation, perform
      call guwk(idwk,1)
      return
      end
      subroutine labelf(nf,ax,ay,chh)
c this subroutine puts a numerical label into graphs.
c the label consists of the # sign followed by the integer nf,
c left justified.  only the first 7 digits of nf are used. 
c nf = number to be printed
c ax,ay = print location of lower left hand corner of label
c chh = character height
c ifrg = index of foreground color
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      character*8 lbl
c first find how many digits in nf
      id = 0
      n = nf
   10 id = id + 1
      n = n/10
      if (n.gt.0) go to 10
c create label template
      lbl = '#       '
c create left justified label
      is = ichar('0')
      if (id.gt.7) id = 7
      ls = 10**(id - 1)
      nt = nf
      do 20 i = 1, id
      i1 = i + 1
      n = nt/ls
      lbl(i1:i1) = char(n+is)
      nt = nt - n*ls
      ls = ls/10
   20 continue
c set character height
      call gschh(chh)
c set text color index to foreground
      call gstxci(ifrg)
c draw text
      call gtx(ax,ay,lbl)
      return
      end
      subroutine cleanp(irc)
c this subroutine cleans up display by resetting plot location,
c updating output, forcing device into request mode and displaying help
c idwk = workstation identifier
c iplot = plot location on page, 0 <= iplot < nplot
c iupd = (-1,0,1) = (no,default,yes) end plot
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      irc = 0
c pause if plots are still pending
      if (((iplot.ne.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         call readrc(irc)
      endif
c reset plot number
      iplot = 0
c display help for normal return code
      if ((irc.eq.0).or.((irc.ge.16).and.(irc.le.127))) then
c force device into request mode
         irc = -1
c read code from input device, if present
         call readrc(irc)
      endif
      return
      end
      subroutine CLRSCRN
c clear screen
c idwk = workstation identifier
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c clear workstation, always
      call gclrwk(idwk,1)
c update workstation, perform
      call guwk(idwk,1)
      return
      end
      subroutine RSTSCRN
c clear screen and force input device into request mode for next time
c idwk = workstation identifier
c iplot = plot location on page, 0 <= iplot < nplot
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      iplot = 0
c clear workstation, always
      call gclrwk(idwk,1)
c update workstation, perform
      call guwk(idwk,1)
      irc = -3
      call readrc(irc)
      return
      end
      subroutine GTINPUT(label,prompt,input,irc)
c display label and prompt and request input
c irc = return code (0 = normal return)
      integer irc
      character*(*) label, prompt, input
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c idstr = string device number, 0 if no string device available
      integer mode, iesw
      real csize, chh, space, ax, ay
c mode = mode of operation of string device (0=request,1=sample,2=event)
c iesw = echo switch (0=noecho,1=echo)
      data mode,iesw /0,1/
      data csize /0.024/
      irc = 0
c exit if not input device
      if (idstr.eq.0) return
c set string mode, 0 = request
      call gsstm(idwk,idstr,mode,iesw)
      if (ry.lt..3) csize = .032
      chh = ry*csize
c space = vertical spacing between characters
      space = 1.5*chh
      ay = ry
c exit if character height too large
      if (ay.lt.space) return
      ax = 0.
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps(chh,1,0)
c display current parameters
c clear workstation, always
      call gclrwk(idwk,1)
c display label
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,label)
c display prompt
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,prompt)
c read input
c update workstation, perform
      call guwk(idwk,1)
      input = ' '
c request string
      call grqst(idwk,idstr,istat,lens,input)
c     read (5,'(a24)',end=35) input
c abort
      if ((input(1:1).eq.'q').or.(input(1:1).eq.'Q')) then
         if (input(2:2).eq.' ') then
            irc = 1
            return
         endif
      endif
      end
      subroutine PTOTPUT(label,prompt)
c display label and prompt
      character*(*) label, prompt
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
      real csize, chh, space, ax, ay
      data csize /0.024/
      irc = 0
      if (ry.lt..3) csize = .032
      chh = ry*csize
c space = vertical spacing between characters
      space = 1.5*chh
      ay = ry
c exit if character height too large
      if (ay.lt.space) return
      ax = 0.
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps(chh,1,0)
c display current parameters
c clear workstation, always
      call gclrwk(idwk,1)
c display prompt
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,label)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,prompt)
c update workstation, perform
      call guwk(idwk,1)
      end
      subroutine GTMINPUT(label,prompt,input,ndim,irc)
c display labels and prompt and request input
c ndim = number of labels
c irc = return code (0 = normal return)
      integer ndim, irc
      character*(*) label, prompt, input
      dimension label(ndim)
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c idstr = string device number, 0 if no string device available
      integer mode, iesw
      real csize, chh, space, ax, ay
c mode = mode of operation of string device (0=request,1=sample,2=event)
c iesw = echo switch (0=noecho,1=echo)
      data mode,iesw /0,1/
      data csize /0.024/
      irc = 0
c exit if not input device
      if (idstr.eq.0) return
c set string mode, 0 = request
      call gsstm(idwk,idstr,mode,iesw)
      if (ry.lt..3) csize = .032
      chh = ry*csize
c space = vertical spacing between characters
      space = 1.5*chh
      ay = ry
c exit if character height too large
      if (ay.lt.space) return
      ax = 0.
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps(chh,1,0)
c display current parameters
c clear workstation, always
      call gclrwk(idwk,1)
c display labels
      do 10 i = 1, ndim
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,label(i))
   10 continue
c display prompt
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,prompt)
c read input
c update workstation, perform
      call guwk(idwk,1)
      input = ' '
c request string
      call grqst(idwk,idstr,istat,lens,input)
c     read (5,'(a24)',end=35) input
c abort
      if ((input(1:1).eq.'q').or.(input(1:1).eq.'Q')) then
         if (input(2:2).eq.' ') then
            irc = 1
            return
         endif
      endif
      end
