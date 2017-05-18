c gks device driver for tektronix graphics using plot10
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c version for ibm rs/6000
c update: november 14, 1997
      block data
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      save /gkspl10/
      data kosv /0/
      end
      subroutine gqops(istat)
c inquire operating state value
c input arguments: none
c istat = operating state (0=gks closed,1=gks open,2=workstation open,
c 3=workstation active,4=segment open)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      istat = kosv
      return
      end
      subroutine gopks(nerrfl,meml)
c open gks
c input arguments: all
c nerrfl = error file unit number, 6 for terminal
c meml = storage limit for one segment, in bytes
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 1
      return
      end
      subroutine gopwk(idwk,idcon,iwtype)
c open workstation
c input arguments: all
c idwk = workstation identifier
c idcon = connection identifier
c iwtype = workstation type
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c iwt = current workstation type
c     call itputg(0,0,irc)
      kosv = 2
      iwk = idwk
      idc = idcon
      iwt = iwtype
      return
      end
      subroutine gqopwk(n,ierr,now,idwk)
c inquire set number of open workstations
c input arguments: n
c n = set member requested
c ierr = error indicator
c now = number of open workstations
c idwk = workstation identifier of the nth member
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      now = 0
      ierr = 0
      if (kosv.ge.2) now = 1
      if ((n.lt.0).or.(n.gt.1)) ierr = -12
      if (n.eq.1) idwk = iwk
      return
      end
      subroutine gqwkca(iwtype,ierr,iwkca)
c inquire workstation category
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator
c iwkca = workstation category
c (0 = output, 1 = input, 2 = outin, 3 = wiss, 4 = mo, 5 = mi)
      ierr = 0
c tektronix terminal
      if ((iwtype.eq.4010).or.(iwtype.eq.4014).or.(iwtype.eq.4016)) then
         iwkca = 2
c color tektronix
      elseif ((iwtype.eq.4105).or.(iwtype.eq.4207)) then
         iwkca = 2
c metafile
      elseif ((iwtype.eq.1).or.(iwtype.eq.2)) then
         iwkca = 0
      else
         ierr = 22
      endif
      return
      end
      subroutine gscnid(idcon,connam)
c set connection identifier
c input arguments: all
c idcon = connection identifier
c connam = connection name
      character*8 connam
c     open(unit=idcon,file=connam,form='formatted',status='unknown')
      return
      end
      subroutine gqwkc(idwk,ierr,idcon,iwtype)
c inquire workstation connection and type
c input arguments: idwk
c idwk = workstation identifier
c ierr = error indicator
c idcon = connection identifier
c iwtype = workstation type
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
c idc = current connection identifier
c iwt = current workstation type
      ierr = 0
      idcon = idc
      if (idwk.eq.iwk) then
         iwtype = iwt
      else
         ierr = 20
      endif
      return
      end
      subroutine gacwk(idwk)
c activate workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwt = current workstation type
      dimension itw(5)
c itw = workstation type conversion table
      data itw /4010,4014,4016,4105,4207/
c icps = transmission rate in characters per second
c iscal = number of addressable points
      data icps,iscal /30,1024/
c special case for metafiles (suppress flush by tsend)
c icps = (0,1) = (no,yes) encode tektronix output data
      if ((iwt.eq.1).or.(iwt.eq.2)) icps = 2 - iwt
c initialization
      call initt(icps)
c find terminal type
      iterm = 0
      do 10 i = 1, 5
      if ((iterm.eq.0).and.(iwt.eq.itw(i))) iterm = i
   10 continue
c default is 4016
      if (iterm.eq.0) iterm = 3
c use largest addressing possible
      if (iterm.eq.3) iscal = 4096
c identify the terminal type
      call term(iterm,iscal)
c read machine architecture description
      call starch
c create tektronix plot10 defaults
      call tkdflts
      kosv = 3
      return
      end
      subroutine gqacwk(n,ierr,naw,idwk)
c inquire set of active workstations
c input arguments: n
c n = set member requested
c ierr = error indicator
c naw = number of active workstations
c idwk = workstation identifier of the nth member
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      naw = 0
      ierr = 0
      if (kosv.ge.3) naw = 1
      if ((n.lt.0).or.(n.gt.1)) ierr = -12
      if (n.eq.1) idwk = iwk
      return
      end
      subroutine gqwks(idwk,ierr,istat)
c inquire workstation state
c input arguments: idwk
c idwk = workstation identifier
c ierr = error indicator
c istat = workstation state (0 = inactive, 1 = active)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      ierr = 0
      istat = 0
      if (kosv.lt.2) ierr = 25
      if (kosv.ge.3) istat = 1
      if (idwk.ne.iwk) then
         istat = 0
         ierr = 20
      endif
      return
      end
      subroutine gqcf(iwtype,ierr,ncoli,iscol,npci)
c inquire color facilities
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator (0=inquiry successful)
c ncoli = number of colors available
c iscol = color availability indicator, 0 = monochrome, 1 = color
c npci = number of predefined color indices
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      ierr = 0
      if ((iwtype.eq.4105).or.(iwtype.eq.4207)) then
         ncoli = 8
      else
         ncoli = 2
      endif
      iscol = 0
      if (ncoli.gt.2) iscol = 1
      npci = ncoli
      return
      end
      subroutine gscr(idwk,ic,cr,cg,cb)
c set color representation
c input arguments: all
c idwk = workstation identifier
c ic = color index
c cr/cg/cb = red/green/blue component (0 < cr,cg,cb < 1)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c icc = color index conversion table
      dimension reds(8), greens(8), blues(8)
      data reds /0.,1.,1.,0.,0.,0.,1.,1./
      data greens /0.,1.,0.,1.,0.,1.,0.,1./
      data blues /0.,1.,0.,0.,1.,1.,1.,0./
c find which tektronix index most closely matches gks index
c and set corresponding element of color index conversion table
      imin = 1
      cdm = cr*cr + cg*cg + cb*cb
      do 10 i = 1, 8
      cd = (cr - reds(i))**2 + (cg - greens(i))**2 + (cb - blues(i))**2
      if (cd.lt.cdm) then
         imin = i
         cdm = cd
      endif
   10 continue
      if ((ic.ge.0).and.(ic.le.7)) icc(ic+1) = imin - 1
      return
      end
      subroutine gqeci(idwk,n,ierr,ncoli,icol)
c inquire list element of color indices
c input arguments: idwk, n
c idwk = workstation identifier
c n = requested list element
c ierr = error indicator
c ncoli = number of indices currently defined
c icol = color index of requested element
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
c iwt = current workstation type
      ierr = 0
      if ((iwt.eq.4105).or.(iwt.eq.4207)) then
         ncoli = 8
      else
         ncoli = 2
      endif
      if (idwk.ne.iwk) ierr = -12
      if ((n.lt.0).or.(n.gt.ncoli)) ierr = -12
      if ((n.gt.0).and.(n.le.ncoli)) icol = n - 1
      return
      end
      subroutine gqdsp(iwtype,ierr,idcun,dcx,dcy,lx,ly)
c inquire maximum display surface size
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator (0=inquiry successful)
c idcun = device coordinate units, (0=meters,1=other)
c dcx/dcy = width/height in device coordinate units
c lx/ly = width/height in device (raster) units
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwt = current workstation type
      idcun = 1
c if workstation requested is open, query actual size
      if ((kosv.eq.3).and.(iwtype.eq.iwt)) then
c check terminal
         call seetrm(ispeed,iterm,isize,maxsr)
         if (maxsr.eq.4096) then
            lx = 4096
            ly = 3120
         else
            lx = 1024
            ly = 780
         endif
c return maximum possible size, if not actually opened
      else
         if ((iwtype.eq.4010).or.(iwtype.eq.4014)) then
            lx = 1024
            ly = 780
         else
            lx = 4096
            ly = 3120
         endif
      endif
      dcx = float(lx) - 1.
      dcy = float(ly) - 1.
      ierr = 0
      return
      end
      subroutine gswkwn(idwk,xmin,xmax,ymin,ymax)
c set workstation window
c input arguments: all
c idwk = workstation identifier
c xmin/xmax = window x coordinates in ndc
c ymin/ymax = window y coordinates in ndc
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c trans = normalization transformations
c clip to maximum window size
      wwnvp(1) = amax1(xmin,0.)
      wwnvp(2) = amin1(xmax,1.)
      wwnvp(3) = amax1(ymin,0.)
      wwnvp(4) = amin1(ymax,1.)
c redefine viewport of transformation 0
      trans(1,1) = wwnvp(1)
      trans(2,1) = wwnvp(2)
      trans(3,1) = wwnvp(3)
      trans(4,1) = wwnvp(4)
      return
      end
      subroutine gswkvp(idwk,xmin,xmax,ymin,ymax)
c set workstation viewport
c input arguments: all
c idwk = workstation identifier
c xmin/xmax = viewport x coordinates in device coordinates
c ymin/ymax = viewport y coordinates in device coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c check terminal
      call seetrm(ispeed,iterm,isize,maxsr)
c clip to maximum screen size
      wwnvp(5) = amax1(xmin,0.)
      wwnvp(7) = amax1(ymin,0.)
      if (maxsr.eq.4096) then
         wwnvp(6) = amin1(xmax,4095.)
         wwnvp(8) = amin1(ymax,3119.)
      else
         wwnvp(6) = amin1(xmax,1023.)
         wwnvp(8) = amin1(ymax,779.)
      endif
      return
      end
      subroutine gqsts(idwk,idstr,mldr,ierr,mode,iesw,lstr,str,ipet,eare
     1a,lenb,ipos,ldr,datar)
c inquire string device state
c input arguments: idwk, idstr, mldr
c idwk = workstation identifier
c idstr = string device number
c mldr = maximum dimension of data record array
c ierr = error indicator (0=inquiry successful)
c mode = operating mode (0=request,1=sample,2=event)
c iesw = echo switch (0=no echo,1=echo)
c lstr = number of characters in initial string
c str = initial string
c ipet = prompt/echo type (1=normal)
c earea(1)/earea(2) = echo area x coordinates in device coordinates
c earea(3)/earea(4) = echo area y coordinates in device coordinates
c lenb = input buffer size
c ipos = inital edit position
c ldr = length of data record array
c datar = data record array
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c cea = character echo area
      character*(*) str
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 1
      lstr = 1
      str(1:1) = ' '
      ipet = 1
      do 10 j = 1, 4
      earea(j) = cea(j)
   10 continue
      lenb = 80
      ipos = 1
      ldr = 1
c special case for metafiles
      if ((iwt.eq.1).or.(iwt.eq.2)) then
         ierr = 140
      else
         ierr = 0
      endif
      return
      end
      subroutine gqlcs(idwk,idloc,itype,mldr,ierr,mode,iesw,nrt,pxi,pyi,
     1ipet,earea,ldr,datar)
c inquire locator device state
c input arguments: idwk, idloc, itype, mldr
c idwk = workstation identifier
c idloc = locator device number
c itype = type of returned value (0=set,1=realized)
c mldr = maximum dimension of data record array
c ierr = error indicator (0=inquiry successful)
c mode = operating mode (0=request,1=sample,2=event)
c iesw = echo switch (0=no echo,1=echo)
c nrt = normalization transformation number for initial position
c pxi/pyi = initial position, in world coordinates
c ipet = prompt/echo type (1=normal)
c earea(1)/earea(2) = echo area x coordinates in device coordinates
c earea(3)/earea(4) = echo area y coordinates in device coordinates
c ldr = length of data record array
c datar = data record array
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c cea = character echo area
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 1
      nrt = 0
      pxi = .5
      pyi = .5
      ipet = 1
      do 10 j = 1, 4
      earea(j) = cea(j)
   10 continue
      ldr = 1
c special case for metafiles
      if ((iwt.eq.1).or.(iwt.eq.2)) then
         ierr = 140
      else
         ierr = 0
      endif
      return
      end
      subroutine ginst(idwk,idstr,lstr,str,ipet,xmin,xmax,ymin,ymax,lenb
     1,ipos,ldr,datar)
c initialize string device
c input arguments: all
c idwk = workstation identifier
c idstr = string device number
c lstr = number of characters in initial string
c str = initial string
c ipet = prompt/echo type (1=normal)
c xmin/xmax = echo area x coordinates in device coordinates
c ymin/ymax = echo area y coordinates in device coordinates
c lenb = input buffer size
c ipos = inital edit position
c ldr = length of data record array
c datar = data record array
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c cea = character echo area
      character*80 datar(ldr)
      character*(*) str
c save characer echo area
      cea(1) = xmin
      cea(2) = xmax
      cea(3) = ymin
      cea(4) = ymax
      return
      end
      subroutine ginlc(idwk,idloc,nrt,px,py,ipet,xmin,xmax,ymin,ymax,ldr
     1,datar)
c initialize locator device
c input arguments: all
c idwk = workstation identifier
c idloc = locator device number
c nrt = initial normalization transformation number
c px/py = initial locator position, in world coordinates
c ipet = prompt/echo type (1=normal)
c xmin/xmax = echo area x coordinates in device coordinates
c ymin/ymax = echo area y coordinates in device coordinates
c ldr = length of data record array
c datar = data record array
      character*80 datar(ldr)
      return
      end
      subroutine gdawk(idwk)
c deactivate workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c termination
      call finitt(0,0)
      kosv = 2
      return
      end
      subroutine gclwk(idwk)
c close workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c     call itputg(0,-1,irc)
      kosv = 1
      return
      end
      subroutine gclks
c close gks
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 0
      return
      end
      subroutine gclrwk(idwk,icofl)
c clear workstation
c input arguments: all
c idwk = workstation identifier
c icofl = control flag (0=conditionally,1=always)
c check terminal modes
      call seemod(line,izaxis,mode)
c check terminal
      call seetrm(ispeed,iterm,isize,maxsr)
c erase screen
      call erase
c modify the z axis of the 4014 terminal
      if (iterm.eq.3) call czaxis(izaxis)
      return
      end
      subroutine gqcntn(ierr,nrt)
c inquire current normalization transformation number
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c nrt = transformation number (0 <= nrt <= 25)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
      ierr = 0
      nrt = nrtn
      return
      end
      subroutine gqnt(nrt,ierr,window,viewpt)
c inquire normalization transformation
c input arguments: nrt
c nrt = transformation number (0 <= nrt <= 25)
c ierr = error indicator (0=inquiry successful)
c window(1)/window(2) = window x coordinates in world coordinates
c window(3)/window(4) = window y coordinates in world coordinates
c viewpt(1)/viewpt(2) = viewport x coordinates in ndc
c viewpt(3)/viewpt(4) = viewport y coordinates in ndc
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      dimension window(4), viewpt(4)
      ierr = 0
      if ((nrt.lt.0).or.(nrt.ge.maxt)) then
         nrt = 0
         ierr = 1
      endif
      n = nrt + 1
      do 10 j = 1, 4
      window(j) = trans(j+4,n)
      viewpt(j) = trans(j,n)
   10 continue
      return
      end
      subroutine gswn(nrt,xmin,xmax,ymin,ymax)
c set window
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
c xmin/xmax = window x coordinates in world coordinates
c ymin/ymax = window y coordinates in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      n = nrt + 1
c store transformation
      trans(5,n) = xmin
      trans(6,n) = xmax
      trans(7,n) = ymin
      trans(8,n) = ymax
      return
      end
      subroutine gsvp(nrt,xmin,xmax,ymin,ymax)
c set viewport
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
c xmin/xmax = viewport x coordinates in ndc
c ymin/ymax = viewport y coordinates in ndc
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      n = nrt + 1
c store transformation
      trans(1,n) = xmin
      trans(2,n) = xmax
      trans(3,n) = ymin
      trans(4,n) = ymax
      return
      end
      subroutine gselnt(nrt)
c select normalization transformation
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c wwnvp = workstation window and viewport
c trans = normalization transformations
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      nrtn = nrt
      n = nrt + 1
c clip to workstation window
      xmin = amax1(wwnvp(1),trans(1,n))
      xmax = amin1(wwnvp(2),trans(2,n))
      ymin = amax1(wwnvp(3),trans(3,n))
      ymax = amin1(wwnvp(4),trans(4,n))
c convert from viewport to screen coordinates
      scx = (wwnvp(6) - wwnvp(5))/(wwnvp(2) - wwnvp(1))
      minx = (xmin - wwnvp(1))*scx + wwnvp(5) + .5
      lenx = (xmax - wwnvp(1))*scx + wwnvp(5) + .5
      lenx = lenx - minx
      scy = (wwnvp(8) - wwnvp(7))/(wwnvp(4) - wwnvp(3))
      miny = (ymin - wwnvp(3))*scy + wwnvp(7) + .5
      leny = (ymax - wwnvp(3))*scy + wwnvp(7) + .5
      leny = leny - miny
c define the screen window
      call swindo(minx,lenx,miny,leny)
c convert window units
      xrange = trans(6,n) - trans(5,n)
      yrange = trans(8,n) - trans(7,n)
c define the virtual window
      call vwindo(trans(5,n),xrange,trans(7,n),yrange)
      return
      end
      subroutine gqln(ierr,ltype)
c inquire linetype
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c ltc = current line type
      ierr = 0
      ltype = ltc
c find the corresponding gks index
      if ((ltype.ge.1).and.(ltype.le.3)) then
         if (ltype.eq.3) ltype = 0
         ltype = ltype + 2 
      else
         ltype = ltype + 1
      endif
      if ((ltype.lt.1).or.(ltype.gt.8)) then
         ierr = 1
         ltype = 1
      endif
      return
      end
      subroutine gslwsc(alwsc)
c set linewidth scale factor
c input arguments: all
c alwsc = linewidth scale factor, (alwsc > 1.0)
      icode = 0
      if (alwsc.gt.1) icode = 1
c check terminal
      call seetrm(ispeed,iterm,isize,maxsr)
c modify the z axis of the 4014 terminal
      if (iterm.eq.3) call czaxis(icode)
      return
      end
      subroutine gsplci(icol)
c set polyline color index
c input arguments: all
c icol = color index
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcl = current line color
c icc = color index conversion table
c convert color index to tektronix
      if ((icol.ge.0).and.(icol.le.7)) then
         lcl = icc(icol+1)
      else
         lcl = icc(2)
      endif
      return
      end
      subroutine gsln(ltype)
c set linetype
c input arguments: all
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c ltc = current line type
c convert from gks to tektronix index
      if ((ltype.ge.2).and.(ltype.le.4)) then
         ltc = ltype - 2
         if (ltc.eq.0) ltc = 3
      else
         ltc = ltype - 1
      endif
      return
      end
      subroutine gpl(n,px,py)
c draw polyline
c input arguments: all
c n = number of points to be connected by a polyline
c px/py = x/y coordinates of points in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c kcf = clipping flag
c lcl = current line color
c ltc = current line type
c wwnvp = workstation window and viewport
      dimension px(n), py(n)
c on color tektronix use color instead of line style codes
      if ((iwt.eq.4105).or.(iwt.eq.4207)) then
         l = lcl - 1
      else
         l = ltc
      endif
c clipping is on
      if (kcf.eq.1) then
c move cursor
         call movea (px(1),py(1))
         do 10 j = 2, n
c draw dashed line
         call dasha (px(j),py(j),l)
   10    continue
c clipping is off
      elseif (kcf.eq.0) then
c get current values of screen window
         call seetw(minx,maxx,miny,maxy)
c get current values of the virtual window limits
         call seedw(xmin,xmax,ymin,ymax)
c calculate transformation factors
         scx = float(maxx - minx)/(xmax - xmin)
         scy = float(maxy - miny)/(ymax - ymin)
         aminx = float(minx) - xmin*scx
         aminy = float(miny) - ymin*scy
c convert to screen coordinates
         ax = px(1)*scx + aminx
         ay = py(1)*scy + aminy
c clip to workstation viewport
         ix = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
         iy = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
c move cursor
         call movabs (ix,iy)
         do 20 j = 2, n
c convert to screen coordinates
         ax = px(j)*scx + aminx
         ay = py(j)*scy + aminy
c clip to workstation viewport
         ix = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
         iy = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
c draw dashed line
         call dshabs (ix,iy,l)
   20    continue
      endif
      return
      end
      subroutine gqmk(ierr,mtype)
c inquire marker type
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c mtype = marker type
c 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross, 6 = square,
c 7 = square with cross, 8 = diamond, 9 = diamond with cross,
c 10 = filled circle
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mtc = current marker type
      ierr = 0
      mtype = mtc
      if ((mtype.lt.1).or.(mtype.gt.10)) then
         ierr = 1
         mtype = 3
      endif
      return
      end
      subroutine gsmksc(amksc)
c set marker size scale factor
c input arguments: all
c amksc = linewidth scale factor, (amksc > 1.0)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c ams = current marker size scale factor
      ams = amksc
      return
      end
      subroutine gspmci(imcol)
c set polymarker color index
c input arguments: all
c imcol = polymarker color index
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcm = current marker color
c icc = color index conversion table
c convert color index to tektronix
      if ((imcol.ge.0).and.(imcol.le.7)) then
         lcm = icc(imcol+1)
      else
         lcm = icc(2)
      endif
      return
      end
      subroutine gsmk(mtype)
c set marker type
c input arguments: all
c mtype = marker type
c 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross, 6 = square,
c 7 = square with cross, 8 = diamond, 9 = diamond with cross,
c 10 = filled circle
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mtc = current marker type
      if ((mtype.ge.1).and.(mtype.le.10)) then
         mtc = mtype
      else
c set marker symbol to star
         mtc = 3
      endif
      return
      end
      subroutine gpm(n,px,py)
c draw polymarker
c input arguments: all
c n = number of points to be connected by a polymarker
c px/py = x/y coordinates of points in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c kcf = clipping flag
c lcm = current marker color
c mtc = current marker type
c ams = current marker size scale factor
c wwnvp = workstation window and viewport
      dimension px(n), py(n)
      character*1 amks(10)
      save amks
      data amks /'.','+','*','o','x','#','H','V','W','@'/
c on color tektronix use color instead of line style codes
      if ((iwt.eq.4105).or.(iwt.eq.4207)) then
         l = lcm - 1
      else
         l = 0
      endif
c get current values of screen window
      call seetw(minx,maxx,miny,maxy)
c get current values of the virtual window limits
      call seedw(xmin,xmax,ymin,ymax)
c calculate transformation factors
      scx = float(maxx - minx)/(xmax - xmin)
      scy = float(maxy - miny)/(ymax - ymin)
c special case of dot markers
      if (mtc.eq.1) then
c clipping is on
         if (kcf.eq.1) then
c on color tektronix use color by setting line style code
            if ((iwt.eq.4105).or.(iwt.eq.4207)) then
c move cursor
               call movea (px(1),py(1))
c draw dashed line
               call dasha (px(1),py(1),l)
            endif
            do 10 j = 1, n
c draw point
            call pointa (px(j),py(j))
   10       continue
c clipping is off
         elseif (kcf.eq.0) then
            aminx = float(minx) - xmin*scx
            aminy = float(miny) - ymin*scy
c on color tektronix use color by setting line style code
            if ((iwt.eq.4105).or.(iwt.eq.4207)) then
c convert to screen coordinates
               ax = px(1)*scx + aminx
               ay = py(1)*scy + aminy
c clip to workstation viewport
               ix = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
               iy = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
c move cursor
               call movabs (ix,iy)
c draw dashed line
               call dshabs (ix,iy,l)
            endif
            do 20 j = 1, n
c convert to screen coordinates
            ax = px(j)*scx + aminx
            ay = py(j)*scy + aminy
c clip to workstation viewport
            ix = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
            iy = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
c draw point
            call pntabs (ix,iy)
   20       continue
         endif
c other markers
      else
c check terminal
         call seetrm(ispeed,iterm,isize,maxsr)
c calculate scaling factor
         sy = .75*ams
         if (maxsr.eq.4096) sy = 4.*sy
c calculate scaling factor
         sx = sy
c set direction cosines of character up angle
         upx = 0.0
         upy = 1.0
c find center of character
         dx = 5.5*sy
         dy = 6.5*sy
c calculate offset
         aminx = float(minx) - (xmin*scx + dx) + .5
         aminy = float(miny) - (ymin*scy + dy) + .5
c clipping is off
         if (kcf.eq.0) then
c reset clipping region to workstation viewport
            minx = wwnvp(5) + .5
            maxx = wwnvp(6) + .5
            miny = wwnvp(7) + .5
            maxy = wwnvp(8) + .5
         endif
c loop over markers
         do 30 j = 1, n
c convert to screen coordinates
         ix = px(j)*scx + aminx
         iy = py(j)*scy + aminy
c draw string
         call csdraw(amks(mtc),ix,iy,sx,sy,upx,upy,minx,maxx,miny,maxy,l
     1)
   30    continue
      endif
      return
      end
      subroutine gqtxal(ierr,itxalh,itxalv)
c inquire text alignment
c input arguments: none
c ierr = error indicator
c itxalh = horizontal component
c 0 = normal, 1 = left, 2 = center, 3 = right
c itxalv = vertical component:
c = 0 normal, 1 = top, 2 = cap, 3 = half, 4 = base, 5 = bottom
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itx = current horizontal text alignment
c ity = current vertical text alignment
      ierr = 0
      itxalh = itx
      itxalv = ity
      return
      end
      subroutine gqchh(ierr,chh)
c inquire character height
c input arguments: none
c ierr = error indicator
c chh = character height, in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kpc = current character precision
c cch = current character height
c stroke precision
      if (kpc.eq.2) then
         chh = cch
c string or character precision
      else
c measure size of character
         call csize(ihorz,ivert)
c get current values of screen window
         call seetw(minx,maxx,miny,maxy)
c get current values of the virtual window limits
         call seedw(xmin,xmax,ymin,ymax)
c ignore white space
         ivt = .7*float(ivert)
c convert to user coordinates
         chh = (float(ivt)/float(maxy - miny))*(ymax - ymin)
      endif
      ierr = 0
      return
      end
      subroutine gqtxp(ierr,itxp)
c inquire text path
c input arguments: none
c ierr = error indicator
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itxtp = current text path
      itxp = itxtp
      ierr = 0
      return
      end
      subroutine gqchup(ierr,chux,chuy)
c inquire character up vector
c input arguments: none
c ierr = error indicator
c chux/chuy = up vector, in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c chuv = character up vector
      chux = chuv(1)
      chuy = chuv(2)
      ierr = 0
      return
      end
      subroutine gstxal(itxalh,itxalv)
c set text alignment
c input arguments: all
c itxalh = horizontal component:
c 0 = normal, 1 = left, 2 = center, 3 = right
c itxalv = vertical component:
c = 0 normal, 1 = top, 2 = cap, 3 = half, 4 = base, 5 = bottom
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itx = current horizontal text alignment
c ity = current vertical text alignment
      if ((itxalh.ge.0).and.(itxalh.le.3)) itx = itxalh
      if ((itxalv.ge.0).and.(itxalh.le.3)) ity = itxalv
      return
      end
      subroutine gstxfp(nfont,iprec)
c set text font
c input arguments: all
c nfont = character font number
c iprec = text precision (0=string,1=character,2=stroke)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kpc = current character precision
      kpc = iprec
      return
      end
      subroutine gstxp(itxp)
c set text path
c input arguments: all
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itxtp = current text path
      itxtp = itxp
      return
      end
      subroutine gstxci(itcol)
c set text color index
c input arguments: all
c itcol = text color index
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lct = current text color
c icc = color index conversion table
c convert color index to tektronix
      if ((itcol.ge.0).and.(itcol.le.7)) then
         lct = icc(itcol+1)
      else
         lct = icc(2)
      endif
      return
      end
      subroutine gschh(chh)
c set character height
c input arguments: all
c chh = character height, in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c cch = current character height
      cch = chh
c get current values of screen window
      call seetw(minx,maxx,miny,maxy)
c get current values of the virtual window limits
      call seedw(xmin,xmax,ymin,ymax)
c convert to screen coordinates
      scy = float(maxy - miny)/(ymax - ymin)
      ivt = (chh/.7)*scy + .5
c check terminal
      call seetrm(ispeed,iterm,isize,maxsr)
c correct character size if extended addressing is used
      if (maxsr.eq.4096) ivt = ivt/4
c calculate character size code
      isz = 4 - min0(3,max0(0,(ivt-7)/5))
c change the character size, if necessary
      if (isz.ne.isize) call chrsiz(isz)
      return
      end
      subroutine gschup(chux,chuy)
c set character up vector
c input arguments: all
c chux/chuy = up vector, in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c chuv = character up vector
      chuv(1) = chux
      chuv(2) = chuy
      return
      end
      subroutine gschxp(chxp)
c set character expansion factor
c input arguments: all
c chxp = character expansion factor (>0.)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c chx = character width expansion factor
      if (chxp.gt.0.) then
         chx = chxp
      else
         chx = 1.0
      endif
      return
      end
      subroutine gtx(px,py,chars)
c display text
c input arguments: all
c px/py = starting x/y position of text, in world coordinates
c chars = test string to be displayed
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itx = current horizontal text alignment
c ity = current vertical text alignment
c kcf = clipping flag
c lct = current text color
c kpc = current character precision
c cch = current character height
c wwnvp = workstation window and viewport
c chuv = character up vector
      character*(*) chars
      n = len(chars)
c on color tektronix use color instead of line style codes
      if ((iwt.eq.4105).or.(iwt.eq.4207)) then
         l = lct - 1
      else
         l = 0
      endif
c get current values of screen window
      call seetw(minx,maxx,miny,maxy)
c get current values of the virtual window limits
      call seedw(xmin,xmax,ymin,ymax)
c calculate transformation factors
      scx = float(maxx - minx)/(xmax - xmin)
      scy = float(maxy - miny)/(ymax - ymin)
c string precision
      if (kpc.eq.2) then
c calculate character box size, in user coordinates
         cbx = .93*cch*chx*(scy/scx)
         cby = cch
c calculate scaling factor
         sy = cch*scy/15.
         sx = sy*chx
c set direction cosines of character up angle
         achu = 1./sqrt(chuv(1)*chuv(1) + chuv(2)*chuv(2))
         upx = chuv(1)*achu
         upy = chuv(2)*achu
c stroke or character precision
      else
c measure size of character
         call csize(ihorz,ivert)
c convert box size to user coordinates, ignore white space in height
         cbx = float(ihorz)/scx
         cby = .7*float(ivert)/scy
      endif
c determine horizontal offset and termination
      cx = px
      dx = float(n)*cbx
      if (itx.eq.2) then
         cx = cx - .5*dx
      elseif (itx.eq.3) then
         cx = cx - dx
      endif
      dx = dx + cx
c determine vertical offset
      cy = py
      if (ity.eq.3) then
         cy = cy - .5*cby
      elseif ((ity.eq.1).or.(ity.eq.2)) then
         cy = cy - cby
      endif
c clipping is on
      if (kcf.eq.1) then
c string precision
         if (kpc.eq.2) then
c convert to screen coordinates
            ix = (cx - xmin)*scx + float(minx) + .5
            iy = (cy - ymin)*scy + float(miny) + .5
c draw string
            call csdraw(chars,ix,iy,sx,sy,upx,upy,minx,maxx,miny,maxy,l)
c stroke or character precision
         else
c determine if characters should be clipped
            ns = (amax1(cx,xmin) - cx)/cbx
            nr = (dx - amin1(dx,xmax))/cbx
            nc = n - (ns + nr)
c move cursor
            call movea (cx,cy)
c output characters
            if (nc.gt.0) call aoutst(nc,chars(1+ns:n-nr))
         endif
c clipping is off
      elseif (kcf.eq.0) then
c convert to screen coordinates
         ix = (cx - xmin)*scx + float(minx) + .5
         iy = (cy - ymin)*scy + float(miny) + .5
c set clipping region to workstation viewport
         mnx = wwnvp(5) + .5
         mxx = wwnvp(6) + .5
         mny = wwnvp(7) + .5
         mxy = wwnvp(8) + .5
c string precision
         if (kpc.eq.2) then
c draw string
            call csdraw(chars,ix,iy,sx,sy,upx,upy,mnx,mxx,mny,mxy,l)
c stroke or character precision
         else
c determine if characters should be clipped
            jx = (dx - xmin)*scx + float(minx) + .5
            ns = float(max0(ix,mnx) - ix)/float(ihorz)
            nr = float(jx - min0(jx,mxx))/float(ihorz)
            nc = n - (ns + nr)
c clip to workstation viewport
            ix = min0(max0(ix,mnx),mxx)
            iy = min0(max0(iy,mny),mxy)
c move cursor
            call movabs (ix,iy)
c output characters
            if (nc.gt.0) call aoutst(nc,chars(1+ns:n-nr))
         endif
      endif
      return
      end
      subroutine gsfaci(ifcol)
c set fill area color index
c input arguments: all
c ifcol = color index
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcf = current fill color
c icc = color index conversion table
      if ((ifcol.ge.0).and.(ifcol.le.7)) then
         lcf = icc(ifcol+1)
      else
         lcf = icc(2)
      endif
      return
      end
      subroutine gsfais(ints)
c set fill area interior style
c input arguments: all
c ints = desired interior style:
c 0 = hollow (default), 1 = solid, 2 = pattern, 3 = hatch
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c infs = interior fill style
      infs = ints
      return
      end
      subroutine gfa(n,px,py)
c fill area
c input arguments: all
c n = number of points in fill area
c px,py = arrays of points, in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcf = current fill color
c kcf = clipping flag
c infs = interior fill style
c wwnvp = workstation window and viewport
      dimension px(n), py(n)
c on color tektronix use color instead of line style codes
      if ((iwt.eq.4105).or.(iwt.eq.4207)) then
         l = lcf - 1
      else
         l = 0
      endif
c clipping is on
      if (kcf.eq.1) then
c move cursor
         call movea (px(1),py(1))
         do 10 j = 2, n
c draw dashed line
         call dasha (px(j),py(j),l)
   10    continue
c draw dashed line
         call dasha (px(1),py(1),l)
c clipping is off
      elseif (kcf.eq.0) then
c get current values of screen window
         call seetw(minx,maxx,miny,maxy)
c get current values of the virtual window limits
         call seedw(xmin,xmax,ymin,ymax)
c calculate transformation factors
         scx = float(maxx - minx)/(xmax - xmin)
         scy = float(maxy - miny)/(ymax - ymin)
         aminx = float(minx) - xmin*scx
         aminy = float(miny) - ymin*scy
c convert to screen coordinates
         ax = px(1)*scx + aminx
         ay = py(1)*scy + aminy
c clip to workstation viewport
         ix0 = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
         iy0 = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
c move cursor
         call movabs (ix0,iy0)
         do 20 j = 2, n
c convert to screen coordinates
         ax = px(j)*scx + aminx
         ay = py(j)*scy + aminy
c clip to workstation viewport
         ix = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
         iy = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
c draw dashed line
         call dshabs (ix,iy,l)
   20    continue
c draw dashed line
         call dshabs (ix0,iy0,l)
      endif
      return
      end
      subroutine gca(px,py,qx,qy,icxd,icyd,ncs,nrs,idx,idy,icola)
c cell array
c input arguments: all
c px,py = lower-left cell corner, in world coordinates
c qx,qy = upper-right cell corner, in world coordinates
c icxd,icyd = color index array dimensions
c ncs,nrs = starting column and row in the color index array
c idx,idy = number of columns and rows in the cell array
c icola = color index array
c restriction: only the 3 lowest order bits of icola are displayed
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c icc = color index conversion table
c wwnvp = workstation window and viewport
      dimension icola(icxd,icyd)
c lxm, lym = maximum address of pixels in x, y
      lxm = wwnvp(6)
      lym = wwnvp(8)
c find location of upper left and lower right hand corner of image
      xu = px
      xl = qx
      yu = amax1(py,qy)
      yl = amin1(py,qy)
c get current values of screen window
      call seetw(minx,maxx,miny,maxy)
c get current values of the virtual window limits
      call seedw(xmin,xmax,ymin,ymax)
c calculate transformation factors
      scx = float(maxx - minx)/(xmax - xmin)
      scy = float(maxy - miny)/(ymax - ymin)
c convert to screen coordinates
      ax = (xu - xmin)*scx + float(minx) 
      ay = (yu - ymin)*scy + float(miny)
      bx = (xl - xmin)*scx + float(minx) 
      by = (yl - ymin)*scy + float(miny)
c clip to workstation viewport
      ix0 = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
      iy0 = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
      ix1 = amin1(amax1(bx,wwnvp(5)),wwnvp(6)) + .5
      iy1 = amin1(amax1(by,wwnvp(7)),wwnvp(8)) + .5
c move cursor
      call movabs (ix0,iy0)
c calculate initial offsets
      ipx = ix0 - 1
      ipy = iy0 + 1
      apx = float(ix0) + .5
      apy = 1.5
c calculate maximum index
      lxs = idx
c     if ((ipx+idx).gt.lxm) lxs = lxm - ipx
      lys = idy
c     if ((ipy-idy).lt.0) lys = ipy
c calculate scalings for pixels (dx = dy = 1., for no rescaling) 
      dx = float(ix1 - ix0)/float(lxs - 1) 
      dy = float(lys - 1)/float(iy0 - iy1)
c invf = (0,1) = (no,yes) image should be inverted vertically
      if (py.ge.qy) then
         invf = 0
         koff = nrs - 1
      else
         invf = 1
         koff = lys - nrs + 2
      endif
      km = iy0 - iy1 + 1
      joff = ncs - 1
c outer loop over rows
      do 20 kk = 1, km
      k = dy*float(kk - 1) + apy
c normal image
      if (invf.eq.0) then
         k1 = k + koff
c inverted image
      else
         k1 = koff - k
      endif
      iy = ipy - kk
c reset previous color code to background
      idr = 0
c next loop over bytes in row of color plane
      do 10 j = 1, lxs
      ix = dx*float(j - 1) + apx
      itc = icola(j+joff,k1)
c clip to 8 colors
      itc = itc - (itc/8)*8
c no change in color
      if (itc.eq.idr) then
c draw line in current color if at end of picture
         if ((ix.eq.ix1).and.(itc.gt.0)) call dshabs(ix,iy,icc(idr+1)-1)
c color change
      else
c previous color was not background, draw to previous point
         if (idr.gt.0) call dshabs(ix-1,iy,icc(idr+1)-1)
c reset previous color code to current color
         idr = itc
c if current color is not background, move to current point
         if (itc.gt.0) then
            call dshabs(ix,iy,-1)
c draw line in current color if at end of picture
            if ((ix.eq.ix1).and.(itc.gt.0)) call dshabs(ix,iy,icc(idr+1)
     1-1)
         endif
      endif
   10 continue
   20 continue
      return
      end
      subroutine gqclip(ierr,indcl,clrect)
c inquire clipping indicator
c input arguments: none
c ierr = error indicator
c indcl = clipping indicator (0=no clip, 1=clip)
c clrect = clipping rectangle, in ndc
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c kcf = clipping flag
c trans = normalization transformations
      dimension clrect(4)
      indcl = kcf
c find clipping rectangle
      n = nrtn + 1
      do 10 j = 1, 4
      clrect(j) = trans(j,n)
   10 continue
      ierr = 0
      return
      end
      subroutine gsclip(iclsw)
c set clipping indicator
c input arguments: all
c iclsw = clipping switch (0=no clip,1=clip)
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kcf = clipping flag
      kcf = iclsw
      return
      end
      subroutine guwk(idwk,iregfl)
c update workstation
c input arguments: all
c idwk = workstation identifier
c iregfl = regeneration flag (0=postponed,1=perform)
c dump the output buffer
      call tsend
      return
      end
      subroutine grqst(idwk,idstr,istat,lostr,str)
c request string
c input arguments: idwk, idstr
c idwk = workstation identifier
c idstr = string device number
c istat = return status (0=none,1=ok)
c lostr = number of characters in string
c str = returned string
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c cea = character echo area
      character*(*) str
c measure size of character
      call csize(ihorz,ivert)
c move cursor to echo area
      ax = cea(1)
      ay = cea(4) - .7*float(ivert)
c clip to workstation viewport
      ix = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
      iy = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
c move cursor
      call movabs(ix,iy)
c enter alphanumeric mode, and flush buffer
      call anmode
c accept characters from terminal
      call ainst(80,str)
c determine how many characters in input
      n = len(str)
      lostr = n
      i = 0
   10 i = i + 1
      if (i.gt.n) go to 20
      its = ichar(str(i:i))
c try again if printable character found
      if ((its.ge.32).and.(its.ne.ichar(' '))) go to 10
c control character found
      if (its.lt.32) then
         lostr = i
      elseif (i.lt.n) then
c two blank characters in a row will be considered end of record
         if (str(i+1:i+1).eq.' ') then
            lostr = i - 1
c if only one blank is found, try again
         else
            go to 10
         endif
      endif
   20 istat = 1
      return
      end
      subroutine grqlc(idwk,idloc,istat,nrt,px,py)
c request locator
c input arguments: idwk, idloc
c idwk = workstation identifier
c idloc = locator device number
c istat = return status (0=none,1=ok)
c nrt = normalization transformation number
c px/py = position, in world coordinates
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c trans = normalization transformations
c use the screen cursor
      call scursr(ichr,ix,iy)
c calculate transformation factors
      scx = (trans(6,1) - trans(5,1))/(wwnvp(6) - wwnvp(5))
      scy = (trans(8,1) - trans(7,1))/(wwnvp(8) - wwnvp(7))
c convert to user coordinates
      px = (float(ix) - wwnvp(5))*scx + trans(5,1)
      py = (float(iy) - wwnvp(7))*scy + trans(7,1)
      nrt = 0
      istat = 1
      return
      end
      subroutine gsstm(idwk,idstr,mode,iesw)
c set string mode
c input arguments: all
c idwk = workstation identifier
c idstr = string device number
c mode  =  mode of operation (0 = request,1 = sample,2 = event)
c iesw  =  echo switch (0 = no echo,1 = echo)
      return
      end
      subroutine gslcm(idwk,idloc,mode,iesw)
c set locator mode
c input arguments: all
c idwk = workstation identifier
c idloc = locator device number
c mode  =  mode of operation (0 = request,1 = sample,2 = event)
c iesw  =  echo switch (0 = no echo,1 = echo)
      return
      end
      subroutine gwait(tout,idwk,icl,idnr)
c await event
c input arguments: tout
c tout = timeout period, seconds
c idwk = workstation identifier
c icl = input class:
c (0=no class,1=locator,2=stroke,3=valuator,4=choice,5=pick,6=string)
c idnr = device number
      idwk = 0
      icl = 0
      idnr = 0
      return
      end
      subroutine ggtst(lostr,str)
c get string
c input arguments: none
c lostr = number of characters in string
c str = input string
      character*(*) str
      lostr = 0
      return
      end
      subroutine ggtlc(nrt,px,py)
c get locator
c input arguments: none
c nrt = normalization transformation number
c px/py = position, in world coordinates
      nrt = 0
      return
      end
      subroutine gesc(idfct,ldi,datai,mldr,ldr,datar)
c escape function
c input arguments: idfct, ldi, datai, mldr
c idfct = escape function identifier
c ldi = length of input data record
c datai = input data record
c mldr = maximum length of output data record
c ldr = length of output data record
c datar = output data record
      character*80 datai(ldi), datar(mldr)
      ldr = 0
c escape functions not supported
      ierr = 180
      return
      end
      subroutine tkdflts
c this subroutines creates default tables for plot10 driver for gks
      parameter(maxt=3)
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c iwt = current workstation type
c nrtn = current normalization transformation number
c itx = current horizontal text alignment
c ity = current vertical text alignment
c kcf = clipping flag
c lcl = current line color
c lcm = current marker color
c lct = current text color
c lcf = current fill color
c ltc = current line type
c mtc = current marker type
c infs = interior fill style
c itxtp = current text path
c kpc = current character precision
c icc = color index conversion table
c ams = current marker size scale factor
c chx = character width expansion factor
c cch = current character height
c wwnvp = workstation window and viewport
c cea = character echo area
c trans = normalization transformations
c chuv = character up vector
c create color index conversion table
c maps gks indices to tektronix indices: tekt_index = icc(gks_index+1)
      icc(1) = 0
      icc(2) = 1
      icc(3) = 4
      icc(4) = 2
      icc(5) = 7
      icc(6) = 5
      icc(7) = 6
      icc(8) = 3
c set default workstation window to square
      wwnvp(1) = 0.
      wwnvp(2) = 1.
      wwnvp(3) = 0.
      wwnvp(4) = 1.
c set default workstation viewport to square
c check terminal
      call seetrm(ispeed,iterm,isize,maxsr)
c maxsr = maximum screen address
      if (maxsr.eq.4096) then
         wwnvp(5) = 488.
         wwnvp(6) = 3607.
         wwnvp(7) = 0.
         wwnvp(8) = 3119.
      else
         wwnvp(5) = 122.
         wwnvp(6) = 901.
         wwnvp(7) = 0.
         wwnvp(8) = 779.
      endif
c make default character echo area equal to workstation viewport
      do 10 j = 1, 4
      cea(j) = wwnvp(j+4)
   10 continue
c default transformation
      nrtn = 0
      trans(1,1) = 0.
      trans(2,1) = 1.
      trans(3,1) = 0.
      trans(4,1) = 1.
      trans(5,1) = 0.
      trans(6,1) = 1.
      trans(7,1) = 0.
      trans(8,1) = 1.
      do 30 k = 2, maxt
      do 20 j = 1, 8
      trans(j,k) = trans(j,1)
   20 continue
   30 continue
      n = nrtn + 1
c convert from viewport to screen coordinates
      minx = wwnvp(5)
      lenx = wwnvp(6) - wwnvp(5)
      miny = wwnvp(7)
      lenx = wwnvp(8) - wwnvp(7)
c define the screen window
      call swindo(minx,lenx,miny,leny)
c convert window units
      xrange = trans(6,n) - trans(5,n)
      yrange = trans(8,n) - trans(7,n)
c define the virtual window
      call vwindo(trans(5,n),xrange,trans(7,n),yrange)
c set default clipping state to on
      kcf = 1
c default colors
      lcl = 1
      lcm = 1
      lct = 1
      lcf = 1
c set default line type to solid
      ltc = 0
c set default marker symbol to star
      mtc = 3
c default text alignment
      itx = 0
      ity = 0
c default marker size scale factor
      ams = 1.0
c set current text path to right
      itxtp = 0
c set current character precision to string
      kpc = 0
c default character width expansion factor
      chx = 1.0
c default character height
      cch = 0.01
c set default character up vector
      chuv(1) = 0.
      chuv(2) = 1.
c set default character size to medium-small
      call chrsiz(3)
c set default fill area interior style to hollow
      infs = 0
      return
      end
      subroutine csdraw(chr,ix,iy,sx,sy,upx,upy,minx,maxx,miny,maxy,lwt)
c this subroutine draws a string at location ix,iy using a scalable,
c rotatable tektronix "stick" font
c input: all
c chr = input string to be drawn
c ix, iy = x, y coordinate of lower left hand corner of character
c sx, sy = scaling factor in x, y relative to a 11 x 13 raster
c upx, upy = x, y direction cosines of character up angle
c minx = the minimum horizontal screen coordinate
c maxx = the maximum horizontal screen coordinate
c miny = the minimum vertical screen coordinate
c maxy = the maximum vertical screen coordinate
c lwt = line type used in drawing characters
c the font is defined in an 11 x 13 raster by draws and moves with 4 bit
c addresses.  the initial location ix,iy is the lower left hand corner
c of the character.  addresses are stored as (x,y) pairs (0<x<10,0<y<12)
c and are preceded by a move flag (=15) if the following coordinate
c represents a move rather than a draw.  the 11 x 13 raster is enlarged
c or reduced by the scaling factors (sx,sy).  the character rotation is
c determined by the up vector direction cosines upx, upy, so that the
c actual address used is given by:
c xa = sy*y*upx + sx*x*upy, ya = sy*y*upy - sx*x*upx.
c if any part of the character address is outside the range
c minx < ix+xa < maxx and miny < iy+ya < maxy, it is clipped.
c the 4 bit addresses are stored in packed integers in the array icfont.
c the location of the addresses in icfont for a particular character are
c stored in the array icfloc, and the number of addresses for that
c character are stored in the array icflen.
c lw = number of 4 bit font coordinate addresses per integer word
      parameter(lw=8)
c iebc = (0,1) = (no,yes) input characters are in ebcdic
      common /march/ iebc, irvb, longi, mtype
      character*(*) chr
c lb = scratch array of 4 bit font coordinate addresses
      dimension lb(8)
      dimension icflen(94), icfloc(94)
c icfont = character font coordinate address array
      dimension icfon1(114), icfon2(114), icfont(228)
      dimension ieta(256)
      equivalence (icfon1(1), icfont(1)), (icfon2(1), icfont(115))
      save nw,icflen,icfloc,icfont,ieta
c nw = -2**(4*lw-1) + 1, used for adding and subtracting high bit
      data nw /-2147483647/
c icflen = array of numbers of coordinate addresses in font array
      data icflen /10,10,20,35,43,25,5,9,9,15,10,7,5,5,5,32,12,21,28,9,1
     19,23,9,35,23,10,12,7,10,7,18,47,13,23,25,16,14,9,32,12,15,16,9,7,8
     2,6,27,12,32,14,25,10,13,11,24,7,16,9,9,5,9,7,2,5,24,19,17,22,21,18
     3,28,13,10,11,9,11,24,13,19,22,22,11,21,18,16,11,24,7,22,9,23,10,23
     4,13/
c icfloc = array of locations of characters in font coordinate array
      data icfloc /1,3,5,8,13,19,23,24,26,28,30,32,33,34,35,36,40,42,45,
     149,51,54,57,59,64,67,69,71,72,74,75,78,84,86,89,93,95,97,99,103,10
     25,107,109,111,112,113,114,118,120,124,126,130,132,134,136,139,140,
     3142,144,146,147,149,150,151,152,155,158,161,164,167,170,174,176,17
     48,180,182,184,187,189,192,195,198,200,203,206,208,210,213,214,217,
     5219,222,224,227/
c first half of icfont array
      data icfon1 /-184213676,1543503872,-225259652,2030043136,-21785408
     14,1895448655,145227776,-266266598,977822304,-2137868103,-135079911
     27,62914560,-266682549,1008470793,406341963,-164148918,974718726,37
     32244480,-100622159,-1010518880,1612842506,1073741824,-174735360,-1
     467562346,-1073741824,-200853612,-1073741824,-266682532,1358565552,
     5-261460132,1342177280,-217766608,-244752384,-184217600,-266686464,
     6-264207692,-959923574,1209402370,273613227,-224017137,545259520,-2
     757767222,-1433902496,1074397184,-257767222,-1433902225,-2036030848
     8,537001984,-133644182,1610612736,-88031128,1783244801,1048576,-896
     903392,-1608382454,709386336,-255146715,0,-121454432,-2107086262,67
     a1219744,1114139274,-1463812096,-94216064,-1563899222,671219744,-17
     b4747820,1392508928,-174747820,1395720192,-134059840,-259354716,671
     c08864,-234331456,-257767222,-1433055407,1342177280,-145258860,-182
     d0109754,1196968266,1519957700,-1028616142,335939600,123512736,-261
     e464064,210545320,-2046363542,1244135424,-90655036,-1028620222,3359
     f39610,805306368,209492648,-1533906944,212660327,1872756736,2126621
     g12,1610612736,-90655036,-1028620222,335939610,974103395,217082479,
     h-1398800384,-217641136,1559480256,-266205688,684682412,217759850,0
     i,217710592,207137952,211856384,-264207692/
c second half of icfont array
      data icfon2 /-959923574,1209402370,272629760,210545320,-2046427136
     1,-264207692,-959923574,1209402370,273638560,210545320,-2046386176,
     2-267319286,709386848,-2136815158,-1342177280,-184168692,-140928614
     34,-255839736,171622400,-255830774,1522532352,-255839741,86335314,1
     4887478700,-1393505792,-255814299,257337772,-1408435712,0,-18869452
     54,0,-256241664,-154482170,0,-265979360,-1610612736,-171483136,-244
     6152182,1779409956,35684513,217084552,-1972754430,1048576,-92700032
     7,1612843274,268435456,-99954777,-2010642942,545300736,-263566744,-
     82105515998,134811648,1089602227,-976830715,1157627904,-89603392,-1
     9604171702,1605018240,268500992,217084552,-1972764672,-184184997,15
     a43503872,-137943807,1048576,217743434,0,-238894556,83886080,149975
     b683,-2056974506,2022221472,149975688,-1972764672,-266313080,-19727
     c54430,2097152,-268382453,747416230,-2078014208,-99954773,-19432709
     d06,612672768,149963895,-1990197248,-250476534,675430498,-198741606
     e4,-205318906,135868168,1744830464,-91615327,-2145385976,-260025078
     f,1518338048,-260034045,86327122,1887478696,-1460631040,-87414783,2
     g034694,612672768,-259358710,0,-137968204,-1804377004,890636032,-17
     h1602092,1342177280,-205208138,-1770559914,890503936,-261979257,121
     i2833792/
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
      n = len(chr)
c clip to viewport
      ix0 = min0(max0(ix,minx),maxx)
      iy0 = min0(max0(iy,miny),maxy)
c move cursor
      call movabs(ix0,iy0)
c main loop over characters in string
      do 50 k = 1, n
c ic = ascii code of character to be drawn (33 <= ic <= 126)
c input is in ascii
      if (iebc.eq.0) then
         ic = ichar(chr(k:k))
c input is in ebcdic
      else
         ic = ieta(ichar(chr(k:k))+1)
      endif
c is = offset index into character font coordinate address array
      is = ic - 32
c blanks
      if (is.eq.0) go to 40
c skip if characters are not printable
      if ((is.lt.0).or.(is.gt.94)) go to 50
      id = 0
c set vertical offset for characters with descenders (g,j,p,q,y)
      if ((is.eq.71).or.(is.eq.74).or.(is.eq.80).or.(is.eq.81).or.(is.eq
     1.89)) id = 4
c mw = high-bit operator
      mw = nw - 1
      lwm = lw - 1
      lwp = lw + 1
c lena = number of 4 bit font coordinates addresses for character ic 
      lena = icflen(is)
c lenw = number of integer words for addresses needed by character ic
      lenw = (lena - 1)/lw + 1
c ioff = offset into character font coordinate array for character ic
      ioff = icfloc(is) - 1
c il = number of addresses in packed integer word
      il = lw
c my = (0,1) = processing (x,y) coordinate address
      my = 0
c im = (0,1) = (move,draw) to new coordinate
      im = 0
c clear clipping flag
      ierr = 0
c main loop over packed integers
      do 30 i = 1, lenw
      it = icfont(i+ioff)
c unpack 4 bit addresses from packed integer
      lb(1) = it
c remove high-bit if present 
      if (it.lt.0) lb(1) = lb(1) - mw
      do 10 j = 1, lwm
      lb(j+1) = lb(j)/16
      lb(j) = lb(j) - lb(j+1)*16
   10 continue
c restore high-bit if necessary
      if (it.lt.0) lb(lw) = lb(lw) + 8
c adjust number of addreses in last packed integer
      if (i.eq.lenw) il = lena - (i - 1)*lw
c draw character
      do 20 j = 1, il
      it = lb(lwp-j)
c x coordinate or move flag
      if (my.eq.0) then
c set move flag
         if (it.eq.15) then
            im = 1
c x coordinate
         else
            at = sx*float(it)
            ax = upy*at + float(ix) + .5
            ay = -upx*at + float(iy) + .5
c set y coordinate flag
            my = 1
         endif
c y coordinate
      else
         at = sy*float(it - id)
         jx = ax + upx*at
         jy = ay + upy*at
c clip x coordinate to viewport
         if (jx.lt.minx) then
            jx = minx
c set clipping flag
            if (ix0.eq.minx) ierr = 1
         elseif (jx.gt.maxx) then
            jx = maxx
c set clipping flag
            if (ix0.eq.maxx) ierr = 1
         endif
         ix0 = jx
c clip y coordinate to viewport
         if (jy.lt.miny) then
            jy = miny
c set clipping flag
            if (iy0.eq.miny) ierr = 1
         elseif (jy.gt.maxy) then
            jy = maxy
c set clipping flag
            if (iy0.eq.maxy) ierr = 1
         endif
         iy0 = jy
c move or draw if visible
         if (ierr.eq.0) then
c draw dashed line
            if (im.eq.0) then
               call dshabs(jx,jy,lwt)
c move cursor
            else
               call movabs(jx,jy)
            endif
c move if invisible
         else
            call movabs(jx,jy)
         endif
c reset (x,y) coordinate flag
         my = 0
c reset (move,draw) flag
         im = 0
c clear clipping flag
         ierr = 0
      endif
   20 continue
   30 continue
c find location of next character
   40 at = sx*float(14)
      ix = upy*at + float(ix) + .5
      iy = -upx*at + float(iy) + .5
c clip to viewport
      ix0 = min0(max0(ix,minx),maxx)
      iy0 = min0(max0(iy,miny),maxy)
c move cursor
      call movabs(ix0,iy0)
   50 continue
      return
      end
      subroutine dimagx(image,lx,ly,lz,lenb,npix,nf)
c this subroutine displays raster image stored in character array
c an identity transformation has been assumed
c image = uncompressed single image
c lx, ly = the size of the image, in pixels
c lz = width of picture, in bytes
c lenb = size of picture, in bytes
c npix = number of pixels per byte
c nf = current frame number being processed
c optimized version for tektronix graphics using plot10
c npald = number of palette entries
      parameter(npald=256)
      parameter(maxt=3)
c ifrg = index of foreground color
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c icc = color index conversion table
      common /gkspl10/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,kpc,nscr,icc,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension icc(8), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      character*1 image(lenb)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img8 = scratch integer image array
      dimension img8(8)
c lbl = character variable for numerical label
      character*8 lbl
      save alx
c alx = x coordinate for numerical labels
      data alx /.86/
c nbit = the number of colors, pixel depth
      nbit = 8/npix
c calculate maximum size of image
      lxs = lx
      if (lxs.gt.isx) lxs = isx
      lys = ly
      if (lys.gt.isy) lys = isy
      ix0 = 0
      iy0 = 0
      ix1 = ix0 + lxs - 1
      iy1 = iy0 + lys - 1
      ipx = ix0 - 1
      ipy = iy1 + 1
c normalize label location
      lafx = alx*float(isx)
      lafy = 0
c erase screen
      call erase
c move cursor
      call movabs(ix0,iy1)
c eight bit color
      if (nbit.eq.8) then
c outer loop over rows
         do 20 k = 1, lys
         ioff = lz*(k - 1)
         iy = ipy - k
c reset previous color code to background
         idr = 0
c loop over bytes in row of color plane
         do 10 j = 1, lxs
         ix = j + ipx
c do not use lookup table
         if (lupt.eq.0) then
            itc = ichar(image(j+ioff))
c use lookup table
         else
            itc = ipal(ichar(image(j+ioff))+1)
         endif
c clip to 8 colors
         itc = itc - (itc/8)*8
c no change in color
         if (itc.eq.idr) then
c draw line in current color if at end of picture
            if ((ix.eq.ix1).and.(itc.gt.0)) call dshabs(ix,iy,icc(idr+1)
     1-1)
c color change
         else
c previous color was not background, draw to previous point
            if (idr.gt.0) call dshabs(ix-1,iy,icc(idr+1)-1)
c reset previous color code to current color
            idr = itc
c if current color is not background, move to current point
            if (itc.gt.0) then
               call dshabs(ix,iy,-1)
c draw line in current color if at end of picture
               if ((ix.eq.ix1).and.(itc.gt.0)) call dshabs(ix,iy,icc(idr
     1+1)-1)
            endif
         endif
   10    continue
   20    continue
c nbits per pixel
      else
c maximum width
         lzs = lxs/npix
c convert from nbits per pixel to 8 bits per pixel
         ntc = 2**nbit
         npixm = npix - 1
         npixp = npix + 1
c outer loop over rows
         do 60 k = 1, lys
         ioff = lz*(k - 1)
         iy = ipy - k
c reset previous color code to background
         idr = 0
c loop over bytes in row of color plane
         do 50 j = 1, lzs
         joff = (j - 1)*npix + ipx
         itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
         do 30 i = 1, npixm
         it1 = itc/ntc
         img8(npixp-i) = itc-it1*ntc
         itc = it1
   30    continue
         img8(1) = itc
         do 40 i = 1, npix
         ix = i + joff
c do not use lookup table
         if (lupt.eq.0) then
            itc = img8(i)
c use lookup table
         else
            itc = ipal(img8(i)+1)
         endif
c clip to 8 colors
         itc = itc - (itc/8)*8
c no change in color
         if (itc.eq.idr) then
c draw line in current color if at end of picture
            if ((ix.eq.ix1).and.(itc.gt.0)) call dshabs(ix,iy,icc(idr+1)
     1-1)
c color change
         else
c previous color was not background, draw to previous point
            if (idr.gt.0) call dshabs(ix-1,iy,icc(idr+1)-1)
c reset previous color code to current color
            idr = itc
c if current color is not background, move to current point
            if (itc.gt.0) then
               call dshabs(ix,iy,-1)
c draw line in current color if at end of picture
               if ((ix.eq.ix1).and.(itc.gt.0)) call dshabs(ix,iy,icc(idr
     1+1)-1)
            endif
         endif
   40    continue
   50    continue
   60    continue
      endif
c add label
c first find how many digits in nf
      id = 0
      n = nf
   70 id = id + 1
      n = n/10
      if (n.gt.0) go to 70
c create label template
      lbl = '#       '
c create left justified label
      is = ichar('0')
      if (id.gt.7) id = 7
      ls = 10**(id - 1)
      nt = nf
      do 80 i = 1, id
      i1 = i + 1
      n = nt/ls
      lbl(i1:i1) = char(n+is)
      nt = nt - n*ls
      ls = ls/10
   80 continue
c move cursor
      call movabs(lafx,lafy)
c set color to foreground
      call dshabs(lafx,lafy,icc(ifrg+1)-1)
c set length of string
      n = id + 1
c output characters
      call aoutst(n,lbl)
      return
      end
      integer function kywait()
c special function to request keystroke
c returns: ascii code for keystroke
      character*1 str
c set cursor to lower left corner
      call movabs(0,0)
c enter alphanumeric mode, and flush buffer
      call anmode
c accept characters from terminal
      call ainst(1,str)
      kywait = ichar(str)
      return
      end
