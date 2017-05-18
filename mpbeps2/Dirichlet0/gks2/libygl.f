      subroutine gitwks(idcon,iwtype)
c this is a site-dependent subroutine which returns
c connection identifier and workstation type
c version for Ygl graphics
c iwtype = workstation type
c 1 = gl console
      iwtype = 1
c idcon = connection identifier, 1 seems to work
      idcon = 1
      return
      end
c gks device driver for Ygl 4.0 graphics
c written by viktor k. decyk, ucla
c copyright 1997, regents of the university of california
c version for ibm rs/6000
c update: july 12, 2011
      block data
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      save /gksgl/
      data kosv /0/
      end
      subroutine gqops(istat)
c inquire operating state value
c input arguments: none
c istat = operating state (0=gks closed,1=gks open,2=workstation open,
c 3=workstation active,4=segment open)
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c minx/miny = location of lower left hand corner of window
      minx = 40
      miny = 40
c llx/lly = current window width/height
      llx = 720
      lly = 540
      maxx = minx + llx - 1
      maxy = miny + lly - 1
c constrain location and size of window
      call prefpo(minx,maxx,miny,maxy)
      kosv = 2
      iwk = idwk
      idc = idcon
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      if (iwtype.eq.1) then
         iwkca = 2
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
c idc = current connection identifier
      ierr = 0
      idcon = idc
      if (idwk.eq.iwk) then
         iwtype = 1
      else
         ierr = 20
      endif
      return
      end
      subroutine gacwk(idwk)
c activate workstation
c input arguments: all
c idwk = workstation identifier
      integer winope
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c special case for sgi
c prevent graphical process from being put in background
c     call foregr
c create a window
      idwind = winope(' gkslib ',8)
c create gl defaults
      call gldflts
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      ierr = 0
c     ncoli = 4096
      ncoli = 128
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
c convert to integer
      icr = 255.*cr
      icg = 255.*cg
      icb = 255.*cb
c change a color map entry to a specified RGB value
      call mapcol(ic,icr,icg,icb)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
      ierr = 0
c     ncoli = 4096
      ncoli = 128
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
      idcun = 1
c return size of a window
      call getsiz(llx,lly)
c llx/lly = current window width/height
      lx = llx
      ly = lly
      dcx = real(lx - 1)
      dcy = real(ly - 1)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c return size of a window
      call getsiz(llx,lly)
c llx/lly = current window width/height
      dcx = real(llx - 1)
      dcy = real(lly - 1)
c clip to maximum screen size
      wwnvp(5) = amax1(xmin,0.)
      wwnvp(7) = amax1(ymin,0.)
      wwnvp(6) = amin1(xmax,dcx)
      wwnvp(8) = amin1(ymax,dcy)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      ierr = 0
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
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 1
      nrt = 0
      pxi = .5
      pyi = .5
      ipet = 1
c return size of a window
      call getsiz(llx,lly)
c llx/lly = current window width/height
      earea(1) = 0.
      earea(2) = real(llx - 1)
      earea(3) = 0.
      earea(4) = real(lly - 1)
      ldr = 1
      ierr = 0
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c keydv = keyboard device number
c idst = current string device number
c cea = character echo area
      character*80 datar(ldr)
      character*(*) str
c save string device number
      idst = idstr
c save characer echo area
      cea(1) = xmin
      cea(2) = xmax
      cea(3) = ymin
      cea(4) = ymax
c keyboard device number
      KEYBD = 513
      keydv = KEYBD
c enable input device for event queuing
      call qdevic(keydv)
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
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mousdv = botton device number
c idlc = current locator device number
c mosx/mosy = x/y valuator device numbers
      character*80 datar(ldr)
c save locator device number
      idlc = idloc
c right mouse
      MOUSE1 = 101
c middle mouse
      MOUSE2 = 102
c left mouse
      MOUSE3 = 103
c button device number
      mousdv = MOUSE3
c enable input device for event queuing
      call qdevic(mousdv)
c x/y mouse location valuators
      MOUSEX = 266
      MOUSEY = 267
c valuator device numbers
      mosx = MOUSEX
      mosy = MOUSEY
c enable input device for event queuing
      call qdevic(mosx)
      call qdevic(mosy)
c tie two valuators to a button
      call tie(mousdv,mosx,mosy)
      return
      end
      subroutine gdawk(idwk)
c deactivate workstation
c input arguments: all
c idwk = workstation identifier
      integer winget
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c keydv = keyboard device number
c mousdv = botton device number
c idst = current string device number
c idlc = current locator device number
c mosx/mosy = x/y valuator device numbers
c disable mouse if enabled
      if (idlc.gt.0) then
c disable an input device for event queuing
         call unqdev(mosx)
         call unqdev(mosy)
c disable an input device for event queuing
         call unqdev(mousdv)
      endif
c disable an input device for event queuing
      if (idst.gt.0) call unqdev(keydv)
c return identifier of current window
      idwind = winget()
c close the identified window
      call winclo(idwind)
      kosv = 2
      return
      end
      subroutine gclwk(idwk)
c close workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c terminate the graphics program
      call gexit
      kosv = 1
      return
      end
      subroutine gclks
c close gks
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 0
      return
      end
      subroutine gclrwk(idwk,icofl)
c clear workstation
c input arguments: all
c idwk = workstation identifier
c icofl = control flag (0=conditionally,1=always)
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c idvt = device transformation
c wwnvp = workstation window and viewport
c set current color in color map mode to background
      call color(0)
c clipping is on
      if (kcf.eq.1) then
c reset viewport to maximum screen size
         minx = wwnvp(5)
         maxx = wwnvp(6)
         miny = wwnvp(7)
         maxy = wwnvp(8)
c set the area of window used for drawing
         call viewpo(minx,maxx,miny,maxy)
c clear to the screenmask
         call clear
c set the area of window used for drawing
         call viewpo(idvt(1),idvt(2),idvt(3),idvt(4))
c clipping is off
      else
c clear to the screenmask
         call clear
      endif
      return
      end
      subroutine gqcntn(ierr,nrt)
c inquire current normalization transformation number
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c nrt = transformation number (0 <= nrt <= 25)
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      real*4 xmin, xmax, ymin, ymax
c nrtn = current normalization transformation number
c kcf = clipping flag
c idvt = device transformation
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
      idvt(1) = (xmin - wwnvp(1))*scx + wwnvp(5) + .5
      idvt(2) = (xmax - wwnvp(1))*scx + wwnvp(5) + .5
      scy = (wwnvp(8) - wwnvp(7))/(wwnvp(4) - wwnvp(3))
      idvt(3) = (ymin - wwnvp(3))*scy + wwnvp(7) + .5
      idvt(4) = (ymax - wwnvp(3))*scy + wwnvp(7) + .5
c clipping is on
      if (kcf.eq.1) then
c set the area of window used for drawing
         call viewpo(idvt(1),idvt(2),idvt(3),idvt(4))
c define orthographic transformation
         xmin = trans(5,n)
         xmax = trans(6,n)
         ymin = trans(7,n)
         ymax = trans(8,n)
         call ortho2(xmin,xmax,ymin,ymax)
c clipping is off
      else
c reset viewport to maximum screen size
         minx = wwnvp(5)
         maxx = wwnvp(6)
         miny = wwnvp(7)
         maxy = wwnvp(8)
c set the area of window used for drawing
         call viewpo(minx,maxx,miny,maxy)
c calculate new window to preserve scaling
         scx = (trans(6,n) - trans(5,n))/real(idvt(2) - idvt(1))
         scy = (trans(8,n) - trans(7,n))/real(idvt(4) - idvt(3))
         xmin = trans(5,n) + (wwnvp(5) - real(idvt(1)))*scx
         xmax = trans(5,n) + (wwnvp(6) - real(idvt(1)))*scx
         ymin = trans(7,n) + (wwnvp(7) - real(idvt(3)))*scy
         ymax = trans(7,n) + (wwnvp(8) - real(idvt(3)))*scy
c define orthographic transformation
         call ortho2(xmin,xmax,ymin,ymax)
      endif
      return
      end
      subroutine gqln(ierr,ltype)
c inquire linetype
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      integer getlst
c returns the current linestyle
      ltc = getlst()
      ltype = ltc + 1
      ierr = 0
      return
      end
      subroutine gslwsc(alwsc)
c set linewidth scale factor
c input arguments: all
c alwsc = linewidth scale factor, (alwsc > 1.0)
      lws = alwsc
      if (lws.lt.1) lws = 1
      if (lws.gt.3) lws = 3
c specify the linewidth
      call linewi(lws)
      return
      end
      subroutine gsplci(icol)
c set polyline color index
c input arguments: all
c icol = color index
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcl = current line color
      lcl = icol
      return
      end
      subroutine gsln(ltype)
c set linetype
c input arguments: all
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      ltc = 0
      if ((ltype.ge.1).and.(ltype.le.5)) ltc = ltype - 1
c select linestyle 
      call setlin(ltc)
      return
      end
      subroutine gpl(n,px,py)
c draw polyline
c input arguments: all
c n = number of points to be connected by a polyline
c px/py = x/y coordinates of points in world coordinates
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      real*4 ax, ay
c nrtn = current normalization transformation number
c lcl = current line color
c idvt = device transformation
c trans = normalization transformations
      dimension px(n), py(n)
c set current color in color map mode
      call color(lcl)
c move cursor
      ax = px(1)
      ay = py(1)
c special case of a line to itself
      if (n.gt.1) then
         if ((px(2).eq.ax).and.(py(2).eq.ay)) then
c calculate transformation factors
            m = nrtn + 1
            eps = (trans(6,m) - trans(5,m))/real(idvt(2) - idvt(1))
            ax = ax - eps
            eps = (trans(8,m) - trans(7,m))/real(idvt(4) - idvt(3))
            ay = ay - eps
         endif
      endif
c move the current graphics position to a specified point
      call move2(ax,ay)
      do 10 j = 2, n
c draw a line
      ax = px(j)
      ay = py(j)
      call draw2(ax,ay)
   10 continue
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mtc = current marker type
      mtype = mtc
      ierr = 0
      return
      end
      subroutine gsmksc(amksc)
c set marker size scale factor
c input arguments: all
c amksc = linewidth scale factor, (amksc > 1.0)
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c ams = current marker size scale factor
      ams = amksc
      return
      end
      subroutine gspmci(imcol)
c set polymarker color index
c input arguments: all
c imcol = polymarker color index
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcm = current marker color
      lcm = imcol
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mtc = current marker type
      if ((mtype.ge.1).and.(mtype.le.10)) then
c set marker symbol
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
      integer gethei, strwid
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      real*4 ax, ay
c nrtn = current normalization transformation number
c lcm = current marker color
c idvt = device transformation
c trans = normalization transformations
      dimension px(n), py(n)
      character*1 amks(10)
      save amks
      data amks /'.','+','*','o','x','#','H','V','W','@'/
c set current color in color map mode
      call color(lcm)
c special case of dot markers
      if (mtc.eq.1) then
         do 10 j = 1, n
c draw a point in modeling coordinates
         ax = px(j)
         ay = py(j)
         call pnt2(ax,ay)
   10    continue
c other markers
      else
c return maximum character height in current raster font
         ich = gethei()
c return width of the specified string
         icw = strwid(amks(mtc),1)
c convert to world coordinates
         m = nrtn + 1
         scx = (trans(6,m) - trans(5,m))/real(idvt(2) - idvt(1))
         scy = (trans(8,m) - trans(7,m))/real(idvt(4) - idvt(3))
         dx = .5*real(icw)*scx
         dy = .5*real(ich)*scy
         do 20 j = 1, n
c shift origin of character
         ax = px(j) - dx
         ay = py(j) - dy
c update the current text character position
         call cmov2(ax,ay)
c draw string of raster characters on screen
         call charst(amks(mtc),1)
   20    continue
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      integer gethei
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c idvt = device transformation
c trans = normalization transformations
c return maximum character height in current raster font
      ich = gethei()
c convert to world coordinates
      m = nrtn + 1
      scy = (trans(8,m) - trans(7,m))/real(idvt(4) - idvt(3))
      chh = real(ich)*scy
      ierr = 0
      return
      end
      subroutine gqtxp(ierr,itxp)
c inquire text path
c input arguments: none
c ierr = error indicator
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      chux = 0.
      chuy = 1.
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      nft = 0
c select a raster font for drawing text strings
      call font(nft)
      return
      end
      subroutine gstxp(itxp)
c set text path
c input arguments: all
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itxtp = current text path
      itxtp = itxp
      return
      end
      subroutine gstxci(itcol)
c set text color index
c input arguments: all
c itcol = text color index
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lct = current text color
      lct = itcol
      return
      end
      subroutine gschh(chh)
c set character height
c input arguments: all
c chh = character height, in world coordinates
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c cch = current character height
      cch = chh
      return
      end
      subroutine gschup(chux,chuy)
c set character up vector
c input arguments: all
c chux/chuy = up vector, in world coordinates
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      integer gethei, strwid
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      real*4 cx, cy
c nrtn = current normalization transformation number
c itx = current horizontal text alignment
c ity = current vertical text alignment
c lct = current text color
c idvt = device transformation
c trans = normalization transformations
      character*(*) chars
c set current color in color map mode
      call color(lct)
      n = len(chars)
c return maximum character height in current raster font
      ich = gethei()
c return width of the specified string
      icw = strwid(chars,n)
c convert to world coordinates
      m = nrtn + 1
      scx = (trans(6,m) - trans(5,m))/real(idvt(2) - idvt(1))
      scy = (trans(8,m) - trans(7,m))/real(idvt(4) - idvt(3))
      cbx = real(icw)*scx
      cby = real(ich)*scy
c determine horizontal offset
      cx = px
      if (itx.eq.2) then
         cx = cx - .5*cbx
      elseif (itx.eq.3) then
         cx = cx - cbx
      endif
c determine vertical offset
      cy = py
      if (ity.eq.3) then
         cy = cy - .5*cby
      elseif ((ity.eq.1).or.(ity.eq.2)) then
         cy = cy - cby
      endif
c update the current text character position
      call cmov2(cx,cy)
c draw string of raster characters on screen
      call charst(chars,n)
      return
      end
      subroutine gsfaci(ifcol)
c set fill area color index
c input arguments: all
c ifcol = color index
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcf = current fill color
      lcf = ifcol
      return
      end
      subroutine gsfais(ints)
c set fill area interior style
c input arguments: all
c ints = desired interior style:
c 0 = hollow (default), 1 = solid, 2 = pattern, 3 = hatch
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      real*4 ax, ay
c lcf = current fill color
c infs = interior fill style
      dimension px(n), py(n)
c set current color in color map mode
      call color(lcf)
c hollow fill
      if (infs.eq.0) then
c move the current graphics position to a specified point
         ax = px(1)
         ay = py(1)
         call move2(ax,ay)
         do 10 j = 2, n
c draw a line
         ax = px(j)
         ay = py(j)
         call draw2(ax,ay)
   10    continue
c draw a line
         ax = px(1)
         ay = py(1)
         call draw2(ax,ay)
c solid fill
      else
c select pattern for filling polygons and rectangles
         call setpat(0)
c move to starting point for a filled polygon
         ax = px(1)
         ay = py(1)
         call pmv2(ax,ay)
         do 20 j = 2, n
c specify next point in a filled polygon
         ax = px(j)
         ay = py(j)
         call pdr2(ax,ay)
   20    continue
c close a filled polygon
         call pclos
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
      parameter(lxm=1024,lym=1024)
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c idvt = device transformation
c wwnvp = workstation window and viewport
c trans = normalization transformations
      dimension icola(icxd,icyd)
c 4096 color image
      integer*2 image(lxm*lym)
c lxn, lyn = maximum address of pixels in x, y
      lxn = wwnvp(6)
      lyn = wwnvp(8)
c find location of upper left and lower right hand corner of image
      xu = px
      xl = qx
      yu = amax1(py,qy)
      yl = amin1(py,qy)
c calculate transformation factors
      m = nrtn + 1
      scx = real(idvt(2) - idvt(1))/(trans(6,m) - trans(5,m))
      aminx = real(idvt(1)) - trans(5,m)*scx
      scy = real(idvt(4) - idvt(3))/(trans(8,m) - trans(7,m))
      aminy = real(idvt(3)) - trans(7,m)*scy
c convert to screen coordinates
      ax = xu*scx + aminx
      ay = yu*scy + aminy
      bx = xl*scx + aminx
      by = yl*scy + aminy
c clip to workstation viewport
      ix0 = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
      iy0 = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
      ix1 = amin1(amax1(bx,wwnvp(5)),wwnvp(6)) + .5
      iy1 = amin1(amax1(by,wwnvp(7)),wwnvp(8)) + .5
c calculate size of rescaled raster image
      kn = ix1 - ix0 + 1
      km = iy0 - iy1 + 1
      if (kn.gt.lxm) kn = lxm
      if (km.gt.lym) km = lym
c calculate initial offsets
      ipx = ix0 - 1
      apx = real(ix0) + .5
      apy = 1.5
c calculate maximum index
      lxs = idx
c     if (idx.gt.lxn) lxs = lxn
      lys = idy
c     if (idy.gt.lyn) lys = lyn
c calculate scalings for pixels (dx = dy = 1., for no rescaling)
      dx = real(ix1 - ix0)/real(lxs - 1)
      dy = real(lys - 1)/real(iy0 - iy1)
c invf = (0,1) = (no,yes) image should be inverted vertically
      if (py.ge.qy) then
         invf = 1
         koff = lys - nrs + 2
      else
         invf = 0
         koff = nrs - 1
      endif
      joff = ncs - 1
c outer loop over rows
      do 50 kk = 1, km
      k = dy*real(kk - 1) + apy
c normal image
      if (invf.eq.0) then
         k1 = k + koff
c inverted image
      else
         k1 = koff - k
      endif
      ioff = kn*(kk - 1) - ipx
c reset previous color code to background
      idr = 0
c move cursor
      ixp = ix0
c next loop over bytes in row of color plane
      do 40 j = 1, lxs
      ix = dx*real(j - 1) + apx
      itc = icola(j+joff,k1)
c no change in color
      if (itc.eq.idr) then
c draw line in current color if at end of picture
         if (ix.eq.ix1) then
            do 10 i = ixp, ix
            image(i+ioff) = idr
   10       continue
         endif
c color change
      else
c draw to previous point
         if (ix.gt.ix0) then
            do 20 i = ixp, ix - 1
            image(i+ioff) = idr
   20       continue
         endif
c reset previous color code to current color
         idr = itc
c move to current point
         ixp = ix
c draw line in current color if at end of picture
         if (ix.eq.ix1) then
            do 30 i = ixp, ix
            image(i+ioff) = idr
   30       continue
         endif
      endif
   40 continue
   50 continue
c draw a rectangular array of pixels into the frame buffer
      call rectwr(ix0,iy1,ix1,iy0,image)
      return
      end
      subroutine gqclip(ierr,indcl,clrect)
c inquire clipping indicator
c input arguments: none
c ierr = error indicator
c indcl = clipping indicator (0=no clip, 1=clip)
c clrect = clipping rectangle, in ndc
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
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
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      real*4 xmin, xmax, ymin, ymax
c nrtn = current normalization transformation number
c kcf = clipping flag
c idvt = device transformation
c wwnvp = workstation window and viewport
c trans = normalization transformations
      m = nrtn + 1
c clipping is on, but was off
      if ((kcf.eq.0).and.(iclsw.eq.1)) then
c set the area of window used for drawing
         call viewpo(idvt(1),idvt(2),idvt(3),idvt(4))
c define orthographic transformation
         xmin = trans(5,m)
         xmax = trans(6,m)
         ymin = trans(7,m)
         ymax = trans(8,m)
         call ortho2(xmin,xmax,ymin,ymax)
c clipping is off, but was on
      elseif ((kcf.eq.1).and.(iclsw.eq.0)) then
c reset viewport to maximum screen size
         minx = wwnvp(5)
         maxx = wwnvp(6)
         miny = wwnvp(7)
         maxy = wwnvp(8)
c set the area of window used for drawing
         call viewpo(minx,maxx,miny,maxy)
c calculate new window to preserve scaling
         scx = (trans(6,m) - trans(5,m))/real(idvt(2) - idvt(1))
         scy = (trans(8,m) - trans(7,m))/real(idvt(4) - idvt(3))
         xmin = trans(5,m) + (wwnvp(5) - real(idvt(1)))*scx
         xmax = trans(5,m) + (wwnvp(6) - real(idvt(1)))*scx
         ymin = trans(7,m) + (wwnvp(7) - real(idvt(3)))*scy
         ymax = trans(7,m) + (wwnvp(8) - real(idvt(3)))*scy
c define orthographic transformation
         call ortho2(xmin,xmax,ymin,ymax)
      endif
      kcf = iclsw
      return
      end
      subroutine guwk(idwk,iregfl)
c update workstation
c input arguments: all
c idwk = workstation identifier
c iregfl = regeneration flag (0=postponed,1=perform)
      if (iregfl.eq.1) call sleep(0)
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
      integer getcol, qread
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      real*4 cx, cy
c nrtn = current normalization transformation number
c keydv = keyboard device number
c lct = current text color
c cea = character echo area
c idvt = device transformation
c trans = normalization transformations
      character*(*) str
      integer*2 cs
      lostr = 0
      n = len(str)
c convert to world coordinates
      m = nrtn + 1
      scx = (trans(6,m) - trans(5,m))/real(idvt(2) - idvt(1))
      scy = (trans(8,m) - trans(7,m))/real(idvt(4) - idvt(3))
      cx = trans(5,m) + (cea(1) - real(idvt(1)))*scx
      cy = trans(7,m) + (cea(3) - real(idvt(3)))*scy
c returns the current color
      iccl = getcol()
c set current color in color map mode to current text color
      call color(lct)
c read first entry in event queue
   10 idev = qread(cs)
c try again if not keyboard
      if (idev.ne.keydv) go to 10
      is = cs
c quit if control character or too many characters
      if ((is.lt.32).or.(lostr.ge.n)) go to 20
      lostr = lostr + 1
      str(lostr:lostr) = char(is)
c echo character
c update the current text character position
      call cmov2(cx,cy)
c draw string of raster characters on screen
      call charst(str(1:lostr),lostr)
c look for more characters
      go to 10
c restore previous color
c set current color in color map mode
   20 call color(iccl)
      istat = 1
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
      integer qread
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mousdv = botton device number
c mosx/mosy = x/y valuator device numbers
c wwnvp = workstation window and viewport
c trans = normalization transformations
      integer*2 cs
c return position of a window
      call getori(ix,iy)
c read first entry in event queue
   10 idev = qread(cs)
c try again if not mouse
      if (idev.ne.mousdv) go to 10
c try again if mouse event is up
      if (cs.eq.0) go to 10
c read first entry in event queue
   20 jdev = qread(cs)
c try again if no valuators returned
      if ((jdev.ne.mosx).and.(jdev.ne.mosy)) go to 20
c x valuator read first
      if (jdev.eq.mosx) then
         qx = real(cs - ix)
c read first entry in event queue
   30    kdev = qread(cs)
c try again if not y valuator
         if (kdev.ne.mosy) go to 30
         qy = real(cs - iy)
c y valuator read first
      else
         qy = real(cs - iy)
c read first entry in event queue
   40    kdev = qread(cs)
c try again if not y valuator
         if (kdev.ne.mosx) go to 40
         qx = real(cs - ix)
      endif
c calculate transformation factors
      scx = (trans(6,1) - trans(5,1))/(wwnvp(6) - wwnvp(5))
      scy = (trans(8,1) - trans(7,1))/(wwnvp(8) - wwnvp(7))
c convert from device to world coordinates
      px = (qx - wwnvp(5))*scx + trans(5,1)
      py = (qy - wwnvp(7))*scy + trans(7,1)
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
      integer qtest, qread
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
c keydv = keyboard device number
c mousdv = botton device number
c idst = current string device number
c idlc = current locator device number
      integer*2 cs
      idwk = iwk
c checks the contents of the event queue
   10 idev = qtest()
c event exists
      if (idev.ne.0) then
c keyboard event
         if (idev.eq.keydv) then
            icl = 6
            idnr = idst
c mouse event
         elseif (idev.eq.mousdv) then
c read first entry in event queue
            jdev = qread(cs)
c check if mouse is down
            if ((jdev.eq.mousdv).and.(cs.eq.1)) then
               icl = 1
               idnr = idlc
c reject mouse up events
            else
               icl = -1
               idnr = 0
            endif
c reject any other event
         else
c read first entry in event queue
            jdev = qread(cs)
            icl = -1
            idnr = 0
         endif
c no event
      else
         icl = 0
         idnr = 0
      endif
c check for other events
      if (icl.lt.0) go to 10
      return
      end
      subroutine ggtst(lostr,str)
c get string
c input arguments: none
c lostr = number of characters in string
c str = input string
      integer qread
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c keydv = keyboard device number
      character*(*) str
      integer*2 cs
      lostr = 0
c read first entry in event queue
   10 idev = qread(cs)
c try again if not keyboard
      if (idev.ne.keydv) go to 10
      is = cs
      lostr = 1
      str(1:1) = char(is)
      return
      end
      subroutine ggtlc(nrt,px,py)
c get locator
c input arguments: none
c nrt = normalization transformation number
c px/py = position, in world coordinates
      integer qread
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mousdv = botton device number
c mosx/mosy = x/y valuator device numbers
c wwnvp = workstation window and viewport
c trans = normalization transformations
      integer*2 cs
c return position of a window
      call getori(ix,iy)
c read first entry in event queue
   10 idev = qread(cs)
c try again if no valuators returned
      if ((idev.ne.mosx).and.(idev.ne.mosy)) go to 10
      if (idev.eq.mosx) then
         qx = real(cs - ix)
c read first entry in event queue
   20    jdev = qread(cs)
c try again if not y valuator
         if (jdev.ne.mosy) go to 20
         qy = real(cs - iy)
      elseif (idev.eq.mosy) then
         qy = real(cs - iy)
c read first entry in event queue
   40    jdev = qread(cs)
c try again if not y valuator
         if (jdev.ne.mosx) go to 40
         qx = real(cs - ix)
      endif
c calculate transformation factors
      scx = (trans(6,1) - trans(5,1))/(wwnvp(6) - wwnvp(5))
      scy = (trans(8,1) - trans(7,1))/(wwnvp(8) - wwnvp(7))
c convert from device to world coordinates
      px = (qx - wwnvp(5))*scx + trans(5,1)
      py = (qy - wwnvp(7))*scy + trans(7,1)
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
      subroutine gldflts
c this subroutines creates default tables for gl driver for gks
      parameter(maxt=3)
      common /gksgl/ kosv,iwk,idc,keydv,mousdv,nrtn,itx,ity,kcf,lcl,lcm,
     1lct,lcf,mtc,infs,itxtp,idst,idlc,mosx,mosy,idvt,ams,chx,cch,wwnvp,
     2cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
      real*4 xmin, xmax, ymin, ymax
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c keydv = keyboard device number
c mousdv = botton device number
c nrtn = current normalization transformation number
c itx = current horizontal text alignment
c ity = current vertical text alignment
c kcf = clipping flag
c lcl = current line color
c lcm = current marker color
c lct = current text color
c lcf = current fill color
c mtc = current marker type
c infs = interior fill style
c itxtp = current text path
c idst = current string device number
c idlc = current locator device number
c mosx/mosy = x/y valuator device numbers
c idvt = device transformation
c ams = current marker size scale factor
c chx = character width expansion factor
c cch = current character height
c wwnvp = workstation window and viewport
c cea = character echo area
c trans = normalization transformations
c chuv = character up vector
      save /gksgl/
c set default workstation window to square
      wwnvp(1) = 0.
      wwnvp(2) = 1.
      wwnvp(3) = 0.
      wwnvp(4) = 1.
c set default workstation viewport to square
c return size of a window
      call getsiz(lx,ly)
c lx/ly = current window width/height
      if (lx.gt.ly) then
         wwnvp(5) = .5*real(lx - ly)
         wwnvp(8) = real(ly - 1)
         wwnvp(6) = wwnvp(5) + wwnvp(8)
         wwnvp(7) = 0.
      else
         wwnvp(6) = real(lx - 1)
         wwnvp(7) = .5*real(ly - lx)
         wwnvp(5) = 0.
         wwnvp(8) = wwnvp(6) + wwnvp(7)
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
c convert from viewport to screen coordinates
      idvt(1) = wwnvp(5) + .5
      idvt(2) = wwnvp(6) + .5
      idvt(3) = wwnvp(7) + .5
      idvt(4) = wwnvp(8) + .5
c set default clipping state to on
      kcf = 1
c default colors
c change a color map entry to a specified RGB value, 0 = black
      call mapcol(0,0,0,0)
c change a color map entry to a specified RGB value, 255 = white
      call mapcol(1,255,255,255)
      lcl = 1
      lcm = lcl
      lct = lcl
      lcf = lcl
c set current color in color map mode to background
      call color(0)
c clear to the screenmask
      call clear
c set the area of window used for drawing
      call viewpo(idvt(1),idvt(2),idvt(3),idvt(4))
c define orthographic transformation
      xmin = trans(5,1)
      xmax = trans(6,1)
      ymin = trans(7,1)
      ymax = trans(8,1)
      call ortho2(xmin,xmax,ymin,ymax)
c default line styles
c define line style as short-dash
      call deflin(1,7967)
c define line style as dot
      call deflin(2,13107)
c define line style as dash-dot
      call deflin(3,6399)
c define line style as long-dash
      call deflin(4,4095)
c select linestyle to solid
      call setlin(0)
c default marker
c set marker symbol to star
      mtc = 3
c default text alignment
      itx = 0
      ity = 0
c default line width
c specify the linewidth
      call linewi(1)
c default marker size scale factor
      ams = 1.0
c itxtp = current text path
      itxtp = 0
c default character width expansion factor
      chx = 1.0
c default character height 
      cch = 0.01
c set default character up vector
      chuv(1) = 0.
      chuv(2) = 1.
c set default fill area interior style to hollow
      infs = 0
c set default string and locator device numbers
      idst = 0
      idlc = 0
      return
      end
      subroutine dimagx(image,lx,ly,lz,lenb,npix,nf)
c this subroutine displays raster image stored in character array
c if necessary, data is copied to a character array
c an identity transformation has been assumed
c image = uncompressed single image
c lx, ly = the size of the image, in pixels
c lz = width of picture, in bytes
c lenb = size of picture, in bytes
c npix = number of pixels per byte
c nf = current frame number being processed
c optimized for gl graphics
c npald = number of palette entries
      parameter(npald=256)
c lxm, lym = maximum number of pixels in x, y
      parameter(lxm=720,lym=540)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c ifrg = index of foreground color
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      character*1 image(lenb)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img16 = integer*2 image array
      integer*2 img16(lxm*lym)
      character*8 lbl
      real*4 afx, zero
      save nfl,n,zero,csize,alx,lbl,img16
c nfl = previous frame number
      data nfl /1/
c csize = vertical size of characters
c alx = x coordinate for numerical labels
      data zero,csize,alx /0.,.02,.86/
c lbl = character variable for numerical label
      data lbl /'#       '/
c nbit = the number of colors, pixel depth
      nbit = 8/npix
c calculate maximum size of image
      lxs = lx
      if (lxs.gt.isx) lxs = isx
      lys = ly
      if (lys.gt.isy) lys = isy
      if (lxs.gt.lxm) lxs = lxm
      if (lys.gt.lym) lys = lym
c normalize characters and location
      afx = alx*rx
      chl = csize*ry
c eight bit color
      if (nbit.eq.8) then
c display raster without modification
         if (lupt.eq.0) then
            do 20 k = 1, lys
            ioff = lz*(lys + 1 - k)
            joff = lxs*(k - 1)
            do 10 j = 1, lxs
            img16(j+joff) = ichar(image(j+ioff))
   10       continue
   20       continue
c copy to character array with lookup table
         else
            do 40 k = 1, lys
            ioff = lz*(lys + 1 - k)
            joff = lxs*(k - 1)
            do 30 j = 1, lxs
            img16(j+joff) = ipal(ichar(image(j+ioff))+1)
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
c copy to character array without lookup table
         if (lupt.eq.0) then
            do 70 k = 1, lys
c loop over bytes in image
            j1 = lxs*(k - 1) + npix + 1
            ioff = lz*(lys + 1 - k)
            do 60 j = 1, lzs
            itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
            do 50 i = 1, npixm
            it1 = itc/ntc
            img16(j1-i) = itc-it1*ntc
            itc = it1
   50       continue
            img16(j1-npix) = itc
            j1 = j1 + npix
   60       continue
   70       continue
c copy to character array with lookup table
         else
            do 100 k = 1, lys
c loop over bytes in image
            j1 = lxs*(k - 1) + npix + 1
            ioff = lz*(lys + 1 - k)
            do 90 j = 1, lzs
            itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
            do 80 i = 1, npixm
            it1 = itc/ntc
            img16(j1-i) = ipal(itc-it1*ntc+1)
            itc = it1
   80       continue
            img16(j1-npix) = ipal(itc+1)
            j1 = j1 + npix
   90       continue
  100       continue
         endif
      endif
c clear workstation, if frame number is not advanced
      if (nf.le.nfl) then
         call gclrwk(idwk,1)
c otherwise, erase label
      else
c set current color in color map mode
         call color(0)
c update the current text character position
         call cmov2(afx,zero)
c draw string of raster characters on screen
         call charst(lbl,8)
      endif
c draw a rectangular array of pixels into the frame buffer
      call rectwr(0,0,lxs-1,lys-1,img16)
c first find how many digits in nf
      id = 0
      n = nf
  110 id = id + 1
      n = n/10
      if (n.gt.0) go to 110
c create label template
      lbl = '#       '
c create left justfied label
      is = ichar('0')
      if (id.gt.7) id = 7
      ls = 10**(id - 1)
      nt = nf
      do 120 i = 1, id
      i1 = i + 1
      n = nt/ls
      lbl(i1:i1) = char(n+is)
      nt = nt - n*ls
      ls = ls/10
  120 continue
c save current frame number
      nfl = nf
c set current color in color map mode
      call color(ifrg)
c set length of string
      n = id + 1
c update the current text character position
      call cmov2(afx,zero)
c draw string of raster characters on screen
      call charst(lbl,n)
      return
      end
      integer function kywait()
c special function to request keystroke
c returns: ascii code for keystroke
      logical isqueu
      integer qread
      external isqueu, qread
      integer*2 cs
c istrt = (0,1) = (first,subsequent) entry to this subroutine
c keydv = keyboard device number
      save istrt,keydv
      data istrt,keydv /0,513/
c check if keyboard device is enabled first time through
      if (istrt.eq.0) then
c indicate whether a specified device is enabled for queuing
c if not, enable input device for event queuing
         if (.not.isqueu(keydv)) call qdevic(keydv)
      endif
c read first entry in event queue
   10 idev = qread(cs)
c try again if not keyboard
      if (idev.ne.keydv) go to 10
      kywait = cs
      istrt = 1
      return
      end
c gks null device driver for missing subroutines in Ygl graphics
c written by viktor k. decyk, ucla
c copyright 1997, regents of the university of california
c update: september 20, 1997
      subroutine setpat(index)
c select pattern for filling polygons and rectangles
c input argument: index
      integer index
      return
      end
c     logical function isqueu(dev)
c indicates wehter a specific device is queued
c input arguments: dev
c dev = device identifier
c     integer dev
c     isqueu = .false.
c     return
c     end
