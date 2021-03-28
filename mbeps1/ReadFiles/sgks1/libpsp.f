      subroutine gitwks(idcon,iwtype)
c this is a site-dependent subroutine which returns
c connection identifier and workstation type
c version for postscript printer library
c iwtype = workstation type
c 1=psp, 2=psl, 3=gpsp, 4=gpsl, 5=cpsp, 6=cpsl, 7=psp2, 8=psl2
      character*1 c
   91 format (' enter format: (1=psp,2=psl,3=gpsp,4=gpsl,5=cpsp,6=cpsl,q
     1=quit,?=help)')
   92 format (a1)
c write prompt
   10 write (6,91)
      read (5,92,end=20) c
c help requested
      if (c.ne.'?') go to 20
      write (6,*) ' Postscript printer file'
      write (6,*) ' 1 = psp = 720x540 pixels, b&w postscript, portrait'
      write (6,*) ' 2 = psl = 720x540 pixels, b&w postscript, landscape'
      write (6,*) ' 3 = gpsp = 720x540 pixels, gry postscript, portrait'
      write (6,*) ' 4 = gpsl = 720x540 pixels, gry postscript,landscape'
      write (6,*) ' 5 = cpsp = 720x540 pixels, clr postscript, portrait'
      write (6,*) ' 6 = cpsl = 720x540 pixels, clr postscript,landscape'
      write (6,*) ' 7 = psp2 = 720x540 pixels, L2 postscript, portrait'
      write (6,*) ' 8 = psl2 = 720x540 pixels, L2 postscript, landscape'
      go to 10
c convert to integer
   20 id = ichar(c) - ichar('0')
c id = 0 means abort
      if ((c.eq.'q').or.(c.eq.'Q').or.(id.eq.0)) stop 1
c request again if format type is invalid
      if ((id.lt.0).or.(id.gt.8)) go to 10
c return workstation type
      iwtype = id
c idcon = connection identifier, 21 seems to work
      idcon = 21
c open file for postscript images
      open(unit=21,file='pgraph.ps',form='formatted',status='unknown')
c special case for vax
c     open(unit=21,file='pgraph.ps',form='formatted',status='unknown',ca
c    1rriagecontrol='none')
      return
      end
c gks device driver for postscript printer library
c written by viktor k. decyk, ucla
c copyright 1997, regents of the university of california
c version for ibm rs/6000
c update: may 12, 2017
      block data
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c kosv = operating state value
      save /gkspsp/
      data kosv /0/
      end
      subroutine gqops(istat)
c inquire operating state value
c input arguments: none
c istat = operating state (0=gks closed,1=gks open,2=workstation open,
c 3=workstation active,4=segment open)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c iwt = current workstation type
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      iwkca = 0
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c kosv = operating state value
c iwt = current workstation type
   91 format (i3,1x,i3,' translate')
   92 format (f8.6,1x,f8.6,' scale')
c initialize postscript printer device
      call gpinit()
c create postscript defaults
      call psdflts
c save context
      write (21,*) 'save'
c set scale and origin in landscape mode
      if (iwt.eq.(2*(iwt/2))) then
         lsm = 1
         scl = 1.0
         igx = 576
         igy = 36
c set scale and origin in portrait mode
      else
         lsm = 0
         scl = 11./15.
         igx = 36
         igy = 360
      endif
c write out origin
      write (21,91) igx, igy
c rotate if landscape
      if (lsm.eq.1) write (21,*) '90 rotate'
c write out scale
      write (21,92) scl, scl
c save graphics context
      write (21,*) 'gsave'
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
c black and white postscript format
      if ((iwtype.eq.1).or.(iwtype.eq.2)) then
         ncoli = 2
         iscol = 0
         npci = 2
c grayscale or color postscript format
      else
         ncoli = 256
         iscol = 1
         npci = 256
      endif
      return
      end
      subroutine gscr(idwk,ic,cr,cg,cb)
c set color representation
c input arguments: all
c idwk = workstation identifier
c ic = color index
c cr/cg/cb = red/green/blue component (0 < cr,cg,cb < 1)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
      save krvc
c krvc = (0,1) = (no,yes) reverse black and white
      data krvc /1/
c return if index is out of range
      if ((ic.lt.0).or.(ic.gt.255)) return
      jc = ic + 1
c describe palette entry
      vpal(1,jc) = cr
      vpal(2,jc) = cg
      vpal(3,jc) = cb
c reverse black and white
      if (krvc.eq.1) then
         if (((cr.eq.0.).and.(cg.eq.0.).and.(cb.eq.0.)).or.((cr.eq.1.).a
     1nd.(cg.eq.1.).and.(cb.eq.1.))) then
         vpal(1,jc) = 1. - cr
         vpal(2,jc) = 1. - cg
         vpal(3,jc) = 1. - cb
         endif
      endif
c describe grayscale entry
      if ((iwt.ge.1).and.(iwt.le.4)) then
         vpal(4,jc) = .3*vpal(1,jc) + .59*vpal(2,jc) + .11*vpal(3,jc)
      elseif ((iwt.ge.5).and.(iwt.le.8)) then
c describe color entry
         vpal(4,jc) = real(int(255.*vpal(3,jc)) + 256*(int(255.*vpal(2,j
     1c)) + 256*int(255.*vpal(1,jc))))
      endif
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
c iwk = current workstation identifier
c black and white postscript format
      if ((iwt.eq.1).or.(iwt.eq.2)) then
         ncoli = 2
c grayscale or color postscript format
      else
         ncoli = 256
      endif
      ierr = 0
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
c lxm, lym = maximum size of image array
      parameter(lxm=720,lym=540)
      idcun = 1
      lx = lxm
      ly = lym
      dcx = real(lx) - 1.
      dcy = real(ly) - 1.
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
c lxm, lym = maximum size of image array
      parameter(lxm=720,lym=540)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c wwnvp = workstation window and viewport
      dcx = real(lxm) - 1.
      dcy = real(lym) - 1.
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
      character*(*) str
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 1
      lstr = 1
      ipet = 1
      earea(1) = 0.
      earea(2) = 1.
      earea(3) = 0.
      earea(4) = 1.
      lenb = 1
      ipos = 1
      ldr = 1
c no input supported
      ierr = 140
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
      iesw = 0
      nrt = 0
      pxi = .5
      pyi = .5
      ipet = 1
      earea(1) = 0.
      earea(2) = 1.
      earea(3) = 0.
      earea(4) = 1.
      ldr = 1
      ierr = 140
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
      character*80 datar(ldr)
      character*(*) str
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c kosv = operating state value
c write out footer for postscript files
   91 format (a9)
c restore context
      write (21,*) 'grestore restore'
c required comments for minimally conforming applications
      write (21,91) '%%Trailer'
      kosv = 2
      return
      end
      subroutine gclwk(idwk)
c close workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c kosv = operating state value
      kosv = 1
      return
      end
      subroutine gclks
c close gks
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
      save nf
   91 format ('%%Page: ? ',i4)
   92 format (i3,1x,i3,' translate')
   93 format (f8.6,1x,f8.6,' scale')
   94 format ('%%BeginDocument')
   95 format ('%!PS-Adobe-3.0 EPSF-3.0')
   96 format ('%!PS-Adobe-2.1 EPSF-2.0')
   97 format ('%%BoundingBox: 0 0 612 792')
   98 format ('%%Extensions: CMYK')
   99 format ('%%LanguageLevel: 2')
c nf = frame number
      data nf /1/
c restore and resave context
      write (21,*) 'grestore restore save'
c write out new page comments
      write (21,91) nf
c begin encapsulated document
      write (21,94)
      if ((iwt.eq.7).or.(iwt.eq.8)) then
         write (21,95)
      else
         write (21,96)
      endif
      write (21,97)
c flag colorimage command as a level 1 extension
      if ((iwt.eq.5).or.(iwt.eq.6)) write (21,98)
      if ((iwt.eq.7).or.(iwt.eq.8)) write (21,99)
c set scale and origin in landscape mode
      if (iwt.eq.(2*(iwt/2))) then
         lsm = 1
         scl = 1.0
         igx = 576
         igy = 36
c set scale and origin in portrait mode
      else
         lsm = 0
         scl = 11./15.
         igx = 36
         igy = 360
      endif
c write out origin
      write (21,92) igx, igy
c rotate if landscape
      if (lsm.eq.1) write (21,*) '90 rotate'
c write out scale
      write (21,93) scl, scl
c write out current state
      call gpsets()
c save graphics context
      write (21,*) 'gsave'
c increment frame counter
      nf = nf + 1
      return
      end
      subroutine gqcntn(ierr,nrt)
c inquire current normalization transformation number
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c nrt = transformation number (0 <= nrt <= 25)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c nrtn = current normalization transformation number
c wwnvp = workstation window and viewport
c dvt = device transformation
c trans = normalization transformations
   91 format(1x,2(i6,1x),'moveto')
   92 format(1x,2(i6,1x),'lineto')
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
      dvt(1) = (xmin - wwnvp(1))*scx + wwnvp(5)
      dvt(2) = (xmax - wwnvp(1))*scx + wwnvp(5)
      scy = (wwnvp(8) - wwnvp(7))/(wwnvp(4) - wwnvp(3))
      dvt(3) = (ymin - wwnvp(3))*scy + wwnvp(7)
      dvt(4) = (ymax - wwnvp(3))*scy + wwnvp(7)
c clipping is on
      if (kcf.eq.1) then
         minx = dvt(1) + .5
         maxx = dvt(2) + .5
         miny = dvt(3) + .5
         maxy = dvt(4) + .5
c clipping is off
      else
         minx = wwnvp(5) + .5
         maxx = wwnvp(6) + .5
         miny = wwnvp(7) + .5
         maxy = wwnvp(8) + .5
      endif
c establish new clippath
      write (21,*) 'grestore gsave'
      write (21,91) minx, miny
      write (21,92) maxx, miny
      write (21,92) maxx, maxy
      write (21,92) minx, maxy
      write (21,92) minx, miny
      write (21,*) 'closepath clip newpath'
c write out current state
      call gpsets()
      return
      end
      subroutine gqln(ierr,ltype)
c inquire linetype
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c ltc = current line type
      ierr = 0
      ltype = ltc
      if ((ltype.lt.1).or.(ltype.gt.4)) then
         ierr = 1
         ltype = 1
      endif
      return
      end
      subroutine gslwsc(alwsc)
c set linewidth scale factor
c input arguments: all
c alwsc = linewidth scale factor, (alwsc > 1.0)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c alws = current linewidth scale factor
   91 format(f5.2,' setlinewidth')
      alws = alwsc
      if (alws.lt.1.) alws = 1.0
      alws = .25*alws
      write (21,91) alws
      return
      end
      subroutine gsplci(icol)
c set polyline color index
c input arguments: all
c icol = color index
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c lcl = current line color
      lcl = icol
      return
      end
      subroutine gsln(ltype)
c set linetype
c input arguments: all
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c ltc = current line type
c lt = software dash type patterns
c 1 = solid, 2 = short-dash, 3 = dot, 4 = dash-dot, 5 = long-dash
      dimension lt(4,4)
      save lt
   91 format (' [] 0 setdash')
   92 format (' [',4(f4.1,1x),'] 0 setdash')
      data lt /9,6,9,6,5,5,5,5,14,6,4,6,23,7,23,7/
c no change in status
      if (ltype.eq.ltc) return
      if ((ltype.ge.1).and.(ltype.le.4)) then
         ltc = ltype
      else
         ltc = 1
      endif
c set dash pattern
      l = ltc - 1
      if (l.eq.0) then
         write (21,91)
      else
         write (21,92) .5*lt(1,l), .5*lt(2,l), .5*lt(3,l), .5*lt(4,l)
      endif
      return
      end
      subroutine gpl(n,px,py)
c draw polyline
c input arguments: all
c n = number of points to be connected by a polyline
c px/py = x/y coordinates of points in world coordinates
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
c nrtn = current normalization transformation number
c lcl = current line color
c cgrs = current grayscale
c dvt = device transformation
c trans = normalization transformations
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
      dimension px(n), py(n)
   91 format(1x,2(f9.2,1x),'moveto')
   92 format(1x,2(i6,1x),'moveto')
   93 format(1x,2(i6,1x),'lineto')
   94 format(1x,f4.2,' setgray')
   95 format(1x,3(f4.2,1x),'setrgbcolor')
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx + .5
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy + .5
c convert to screen coordinates
      ix0 = px(1)*scx + aminx
      iy0 = py(1)*scy + aminy
c special case of a line to itself
      if (n.gt.1) then
c convert to screen coordinates
         ix = px(2)*scx + aminx
         iy = py(2)*scy + aminy
         if ((ix.eq.ix0).and.(iy.eq.iy0)) then
            ax0 = ix0 - .25
            ay0 = iy0 - .25
c move point
            write (21,91) ax0, ay0
         else
c move point
            write (21,92) ix0, iy0
         endif
      else
c move point
         write (21,92) ix0, iy0
      endif
      do 10 j = 2, n
c convert to screen coordinates
      ix = px(j)*scx + aminx
      iy = py(j)*scy + aminy
c draw line
      write (21,93) ix, iy
   10 continue
c set color
      it = lcl + 1
      gs = vpal(4,it)
c check if color has changed.
      if (gs.ne.cgrs) then
         cgrs = gs
c grayscale postscript format
         if ((iwt.eq.3).or.(iwt.eq.4)) then
            write (21,94) cgrs
c color postscript format
         elseif ((iwt.ge.5).and.(iwt.le.8)) then
            write (21,95) vpal(1,it), vpal(2,it), vpal(3,it)
         endif
      endif
      write (21,*) 'stroke'
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c ams = current marker size scale factor
      ams = amksc
      return
      end
      subroutine gspmci(imcol)
c set polymarker color index
c input arguments: all
c imcol = polymarker color index
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
c nrtn = current normalization transformation number
c mtc = current marker type
c lcm = current marker color
c nft = current font number
c nfts = current font size
c ams = current marker size scale factor
c cgrs = current grayscale
c dvt = device transformation
c trans = normalization transformations
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
      dimension px(n), py(n)
      character*1 amks(10)
      save amks
   91 format(1x,2(f9.2,1x),'moveto')
   92 format(1x,2(i6,1x),'lineto')
   93 format(1x,i6,' hw sub ',i6,' moveto')
   94 format(1x,f4.2,' setgray')
   95 format(1x,i3,' scalefont setfont')
   96 format(1x,3(f4.2,1x),'setrgbcolor')
      data amks /'.','+','*','o','x','#','H','V','W','@'/
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx + .5
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy + .5
c find character height
      is = 10.*ams
      if (is.ne.nfts) then
c load current font
         if (nft.eq.1) then
            write (21,*) '/Courier findfont'
         elseif (nft.eq.2) then
            write (21,*) '/Times-Roman findfont'
         elseif (nft.eq.3) then
            write (21,*) '/Symbol findfont'
         endif
c scale and select font
         write (21,95) is
         nfts = is
      endif
c find character width
      write (21,*) '('//amks(mtc)//') stringwidth pop .5 mul'
      write (21,*) '/hw exch def'
c set color
      it = lcm + 1
      gs = vpal(4,it)
c check if color has changed.
      if (gs.ne.cgrs) then
         cgrs = gs
c grayscale postscript format
         if ((iwt.eq.3).or.(iwt.eq.4)) then
            write (21,94) cgrs
c color postscript format
         elseif ((iwt.ge.5).and.(iwt.le.8)) then
            write (21,96) vpal(1,it), vpal(2,it), vpal(3,it)
         endif
      endif
      do 10 j = 1, n
c convert to screen coordinates
      ix = px(j)*scx + aminx
      iy = py(j)*scy + aminy
c special case of dot markers
      if (mtc.eq.1) then
c draw dot
         ax = ix - .25
         ay = iy - .25
         write (21,91) ax, ay
         write (21,92) ix, iy
c other markers
      else
         iy = iy - .5*real(is)
c draw character string
         write (21,93) ix, iy
         write (21,*) '('//amks(mtc)//') show'
      endif
   10 continue
c make dots visible
      if (mtc.eq.1) write (21,*) 'stroke'
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c cch = current character height
      chh = cch
      ierr = 0
      return
      end
      subroutine gqtxp(ierr,itxp)
c inquire text path
c input arguments: none
c ierr = error indicator
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c nft = current font number
c nfts = current font size
   91 format(1x,i3,' scalefont setfont')
      nft = nfont
c Font 1 = Courier
      if (nft.eq.1) then
         write (21,*) '/Courier findfont'
c Font 2 = Times-Roman
      elseif (nft.eq.2) then
         write (21,*) '/Times-Roman findfont'
c Font 3 = Symbol
      elseif (nft.eq.3) then
         write (21,*) '/Symbol findfont'
      endif
      write (21,91) nfts
      return
      end
      subroutine gstxp(itxp)
c set text path
c input arguments: all
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c itxtp = current text path
      itxtp = itxp
      return
      end
      subroutine gstxci(itcol)
c set text color index
c input arguments: all
c itcol = text color index
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c lct = current text color
      lct = itcol
      return
      end
      subroutine gschh(chh)
c set character height
c input arguments: all
c chh = character height, in world coordinates
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c cch = current character height
      cch = chh
      return
      end
      subroutine gschup(chux,chuy)
c set character up vector
c input arguments: all
c chux/chuy = up vector, in world coordinates
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
c nrtn = current normalization transformation number
c itx = current horizontal text alignment
c ity = current vertical text alignment
c lct = current text color
c nft = current font number
c nfts = current font size
c cch = current character height
c cgrs = current grayscale
c dvt = device transformation
c trans = normalization transformations
c chuv = character up vector
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
      character*(*) chars
   91 format(1x,i3,' scalefont setfont')
   92 format(1x,f4.2,' setgray')
   93 format(1x,3(f4.2,1x),'setrgbcolor')
   94 format(1x,i6,' exch sub ',i6,' moveto')
   95 format(1x,2(i6,1x),'moveto')
      n = len(chars)
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx + .5
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy + .5
c calculate font scaling
      is = 1.2*cch*scy
c convert to screen coordinates
      ix = px*scx + aminx
      iy = py*scy + aminy
c determine vertical offset
      if (ity.eq.3) then
         iy = iy - .5*is
      elseif ((ity.eq.1).or.(ity.eq.2)) then
         iy = iy - is
      endif
      if (is.ne.nfts) then
c load current font
         if (nft.eq.1) then
            write (21,*) '/Courier findfont'
         elseif (nft.eq.2) then
            write (21,*) '/Times-Roman findfont'
         elseif (nft.eq.3) then
            write (21,*) '/Symbol findfont'
         endif
c scale and select font
         write (21,91) is
         nfts = is
      endif
c draw character string
      if (itx.eq.2) then
         write (21,*) '('//chars//') dup stringwidth pop .5 mul'
         write (21,94) ix, iy
      elseif (itx.eq.3) then
         write (21,*) '('//chars//') dup stringwidth pop'
         write (21,94) ix, iy
      else
         write (21,*) '('//chars//')'
         write (21,95) ix, iy
      endif
c set color
      it = lct + 1
      gs = vpal(4,it)
c check if color has changed.
      if (gs.ne.cgrs) then
         cgrs = gs
c grayscale postscript format
         if ((iwt.eq.3).or.(iwt.eq.4)) then
            write (21,92) cgrs
c color postscript format
         elseif ((iwt.ge.5).and.(iwt.le.8)) then
            write (21,93) vpal(1,it), vpal(2,it), vpal(3,it)
         endif
      endif
c calculate character rotation angle
      if ((chuv(1).eq.0.).and.(chuv(2).ge.0.0)) then
         write (21,*) 'show'
      else
         nangle = (360./6.28318530717959)*atan2(-chuv(1),chuv(2)) + .5
         write (21,*) 'gsave', nangle, ' rotate show grestore'
      endif
      return
      end
      subroutine gsfaci(ifcol)
c set fill area color index
c input arguments: all
c ifcol = color index
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
c nrtn = current normalization transformation number
c lcf = current fill color
c infs = interior fill style
c cgrs = current grayscale
c dvt = device transformation
c trans = normalization transformations
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
      dimension px(n), py(n)
   92 format(1x,2(i6,1x),'moveto')
   93 format(1x,2(i6,1x),'lineto')
   94 format(1x,f4.2,' setgray')
   95 format(1x,3(f4.2,1x),'setrgbcolor')
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx + .5
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy + .5
c convert to screen coordinates
      ix0 = px(1)*scx + aminx
      iy0 = py(1)*scy + aminy
c move point
      write (21,92) ix0, iy0
      do 10 j = 2, n
c convert to screen coordinates
      ix = px(j)*scx + aminx
      iy = py(j)*scy + aminy
c draw line
      write (21,93) ix, iy
   10 continue
c draw line
      write (21,*) ix0, iy0, ' lineto closepath'
c set color
      it = lcf + 1
      gs = vpal(4,it)
c check if color has changed.
      if (gs.ne.cgrs) then
         cgrs = gs
c grayscale postscript format
         if ((iwt.eq.3).or.(iwt.eq.4)) then
            write (21,94) cgrs
c color postscript format
         elseif ((iwt.ge.5).and.(iwt.le.8)) then
            write (21,95) vpal(1,it), vpal(2,it), vpal(3,it)
         endif
      endif
c hollow fill
      if (infs.eq.0) then
         write (21,*) 'stroke'
c solid fill
      else
         write (21,*) 'fill'
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
c lxm, lym = maximum size of image array
      parameter(lxm=720,lym=540)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
c nrtn = current normalization transformation number
c dvt = device transformation
c trans = normalization transformations
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
      dimension icola(icxd,icyd)
c lzsm = group size for printing hex characters
      parameter(lzsm=40)
c gs = temporary array for hex output
      character*1 gs(2*lzsm)
   91 format (i3,1x,i3,' translate')
   92 format (i3,1x,i3,' scale')
   93 format (i3,1x,i3,1x,i1,' [',i4,' 0 0 ',i3,1x,i3,' 0]')
   94 format ('{currentfile ',i3,' string readhexstring pop} image')
   95 format ('{currentfile ',i3,' string readhexstring pop} false 3 col
     1orimage')
   96 format (80a1)
c find location of upper left and lower right hand corner of image
      xu = amax1(px,qx)
      xl = amin1(px,qx)
      yu = amax1(py,qy)
      yl = amin1(py,qy)
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx + .5
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy + .5
c convert to screen coordinates
      ixu = xu*scx + aminx
      iyu = yu*scy + aminy
      ixl = xl*scx + aminx
      iyl = yl*scy + aminy
c find image width and height
      il = ixu - ixl + 1
      ih = iyu - iyl + 1
c save context
      write (21,*) 'gsave'
c write out origin
      write (21,91) ixl, iyl
c write out scale
      write (21,92) il, ih
c nbit = number of bits per pixel, nbit < 9
      nbit = 8
c black and white postscript format is monochrome 
      if ((iwt.eq.1).or.(iwt.eq.2)) nbit = 1
c nbc = number of bytes per color
      nbc = 1
c color postscript colorimage command does not use a palette
      if ((iwt.eq.5).or.(iwt.eq.6)) nbc = 3
c npix = number of pixels per character
      npix = 8/nbit
c lz = number of bytes in scan line
      lz = (idx - 1)/npix + 1
c ndx = location of next-to-last byte
      ndx = npix*(lz - 1)
c npl = number of pixels in last byte
      npl = idx - ndx
c nls = left shift operator for last byte
      nls = 2**(npix - npl)
c write out image operator command
      if ((iwt.eq.7).or.(iwt.eq.8)) then
         call indimg(idx,idy)
      else
c write out arguments for image operator
         write (21,93) idx, idy, nbit, idx, idy, 0
c grey-scale image
         if (nbc.eq.1) then
            write (21,94) lz
c color image
         else
            write (21,95) lz
         endif
      endif
c process lines in groups of lzs
      lzs = lzsm/nbc
      if (lzs.gt.lz) lzs = lz
c lg = number of full or partial groups of lzs
      lg = (lz - 1)/lzs + 1
c lgr = size of last group
      lgr = lz - lzs*(lg - 1)
c location of base characters in hex conversion
      i0 = ichar('0')
      ia = ichar('A') - 10
      idy1 = idy + 1
c outer loop over number of lines
      do 50 k = 1, idy
      k1 = idy1 - k
      lzb = lzs
c loop over number of lzs blocks in row
      do 40 j = 1, lg
      ioff = lzs*(j - 1)
      if (j.eq.lg) lzb = lgr
      npm = npix
c loop over lzs items in row
      do 30 i = 1, lzb
      iioff = nbc*(i - 1)
      do 20 ii = 1, nbc
c grey-scale or color image
      if (nbit.eq.8) then
         it = icola(i+ioff,k1)
         if (nbc.eq.1) then
c map to grey-scale
            if ((iwt.eq.3).or.(iwt.eq.4)) then
               il = 255.*vpal(4,it+1)
c map to indexed color image
            else
               il = it
            endif
c map to color image
         else
            il = 255.*vpal(ii,it+1)
         endif
c monochrome image
      else
         noff = npix*(i + ioff - 1)
c last byte may not be full
         if (ioff.eq.ndx) npm = npl
c pack pixels into bytes
         it = 0
         do 10 n = 1, npm
         it = 2*it + icola(n+noff,k1)
   10    continue
c left shift last byte
         if (ioff.eq.ndx) it = nls*it
c invert image
         il = 255 - it
      endif
c ih = 4 high bits
      ih = il/16
c il = 4 low bits
      il = il - 16*ih
c find ascii code for high bits
      if (ih.gt.9) then
         ih = ih + ia
      else
         ih = ih + i0
      endif
c find ascii code for low bits
      if (il.gt.9) then
         il = il + ia
      else
         il = il + i0
      endif
c convert to printable character
      gs(2*(ii+iioff)-1) = char(ih)
      gs(2*(ii+iioff)) = char(il)
   20 continue
   30 continue
c write characters in hex
      write (21,96) (gs(i),i=1,2*nbc*lzb)
   40 continue
   50 continue
c add EOF code for ASCIIHexDecode filter
      if ((iwt.eq.7).or.(iwt.eq.8)) write (21,*) '>'
c restore context
      write (21,*) 'grestore'
      return
      end
      subroutine gqclip(ierr,indcl,clrect)
c inquire clipping indicator
c input arguments: none
c ierr = error indicator
c indcl = clipping indicator (0=no clip, 1=clip)
c clrect = clipping rectangle, in ndc
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c kcf = clipping flag
c wwnvp = workstation window and viewport
c dvt = device transformation
   91 format(1x,2(i6,1x),'moveto')
   92 format(1x,2(i6,1x),'lineto')
c no change in status
      if (iclsw.eq.kcf) return
      kcf = iclsw
c clipping is on
      if (kcf.eq.1) then
         minx = dvt(1) + .5
         maxx = dvt(2) + .5
         miny = dvt(3) + .5
         maxy = dvt(4) + .5
c clipping is off
      else
         minx = wwnvp(5) + .5
         maxx = wwnvp(6) + .5
         miny = wwnvp(7) + .5
         maxy = wwnvp(8) + .5
      endif
c establish new clippath
      write (21,*) 'grestore gsave'
      write (21,91) minx, miny
      write (21,92) maxx, miny
      write (21,92) maxx, maxy
      write (21,92) minx, maxy
      write (21,92) minx, miny
      write (21,*) 'closepath clip newpath'
c write out current state
      call gpsets()
      return
      end
      subroutine guwk(idwk,iregfl)
c update workstation
c input arguments: all
c idwk = workstation identifier
c iregfl = regeneration flag (0=postponed,1=perform)
   91 format ('%%EndDocument')
c end endcapsulated document
      write (21,91)
c write image
      write (21,*) 'showpage'
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
      character*(*) str
      istat = 0
      lostr = 0
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
      istat = 0
      nrt = 0
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
      subroutine psdflts
c this subroutine creates default tables for postscript driver for gks
c lxm, lym = maximum size of image array
      parameter(lxm=720,lym=540)
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
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
c nft = current font number
c nfts = current font size
c alws = current linewidth scale factor
c ams = current marker size scale factor
c chx = character width expansion factor
c cch = current character height
c cgrs = current grayscale
c wwnvp = workstation window and viewport
c dvt = device transformation
c trans = normalization transformations
c chuv = character up vector
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
c set default workstation window to square
      wwnvp(1) = 0.
      wwnvp(2) = 1.
      wwnvp(3) = 0.
      wwnvp(4) = 1.
c set default workstation viewport to square
c set screen size
      lx = lxm
      ly = lym
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
      do 20 k = 2, maxt
      do 10 j = 1, 8
      trans(j,k) = trans(j,1)
   10 continue
   20 continue
c convert from viewport to screen coordinates
      dvt(1) = wwnvp(5)
      dvt(2) = wwnvp(6)
      dvt(3) = wwnvp(7)
      dvt(4) = wwnvp(8)
c set default clipping state to on
      kcf = 1
c default colors
      lcl = 1
      lcm = 1
      lct = 1
      lcf = 1
c cgrs = current grayscale
      cgrs = 0.0
c set default index 0 = white background
      vpal(1,1) = 1.0
      vpal(2,1) = 1.0
      vpal(3,1) = 1.0
      vpal(4,1) = 1.0
c set remaining indices to black
      do 30 j = 2, 256
      vpal(1,j) = 0.0
      vpal(2,j) = 0.0
      vpal(3,j) = 0.0
      vpal(4,j) = 0.0
   30 continue  
c set default line type to solid
      ltc = 1
c set default marker symbol to star
      mtc = 3
c default text alignment
      itx = 0
      ity = 0
c nft = current font number
      nft = 1
c nfts = current font size
      nfts = 10
c default linewidth scale factor
      alws = .25
c default marker size scale factor
      ams = 1.0
c set current text path to right
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
      return
      end
c internal library for postscript printer images
      subroutine gpinit()
c this subroutine initializes postscript printer device
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
   91 format ('%!PS-Adobe-2.1')
   92 format ('%%BoundingBox: 0 0 612 792')
   93 format ('%%Creator: GLIB')
   94 format ('%%DocumentFonts: Courier Times-Roman Symbol')
   95 format ('%%Extensions: CMYK')
   96 format ('%%EndComments')
   97 format ('%%EndProlog')
   98 format ('%!PS-Adobe-3.0')
   99 format ('%%LanguageLevel: 2')
c write postscript header
      if ((iwt.eq.7).or.(iwt.eq.8)) then
         write (21,98)
      else
         write (21,91)
      endif
      write (21,92)
      write (21,93)
      write (21,94)
c flag colorimage command as a level 1 extension
      if ((iwt.eq.5).or.(iwt.eq.6)) write (21,95)
      if ((iwt.eq.7).or.(iwt.eq.8)) write (21,99)
      write (21,96)
      write (21,97)
c read machine architecture description
c     call starch
      return
      end
      subroutine gpsets()
c write out current font, linestyle, linewidth, and grayscale
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c iwt = current workstation type
c ltc = current line type
c nft = current font number
c nfts = current font size
c alws = current linewidth scale factor
c cgrs = current grayscale
c lt = software dash type patterns
      dimension lt(4,4)
      save lt
   91 format (' [] 0 setdash')
   92 format (' [',4(f4.1,1x),'] 0 setdash')
   93 format(f5.2,' setlinewidth')
   94 format(1x,f4.2,' setgray')
   95 format(1x,i3,' scalefont setfont')
   96 format(1x,3(f4.2,1x),'setrgbcolor')
      data lt /9,6,9,6,5,5,5,5,14,6,4,6,23,7,23,7/
c set current font
      if (nft.eq.1) then
         write (21,*) '/Courier findfont'
      elseif (nft.eq.2) then
         write (21,*) '/Times-Roman findfont'
      elseif (nft.eq.3) then
         write (21,*) '/Symbol findfont'
      endif
c scale and select font
      write (21,95) nfts
c set dash pattern
      l = ltc - 1
      if (l.eq.0) then
         write (21,91)
      else
         write (21,92) .5*lt(1,l), .5*lt(2,l), .5*lt(3,l), .5*lt(4,l)
      endif
c set current line width
      write (21,93) alws
c grayscale postscript format
      if ((iwt.eq.3).or.(iwt.eq.4)) then
         write (21,94) cgrs
c color postscript format
      elseif ((iwt.ge.5).and.(iwt.le.8)) then
         it = cgrs
         is = it/256
         blue = real(it - 256*is)/255.
         it = is/256
         green = real(is - 256*it)/255.
         is = it/256
         red = real(it - 256*is)/255.
         write (21,96) red, green, blue
      endif
      return
      end
      subroutine indimg(lx,ly)
c this subroutine stores palette and image operator for level 2
c color images
c lx , ly = width and height of image
      parameter(maxt=3)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
c lzsm = group size for printing hex characters
      parameter(lzsm=24)
c gs = temporary array for hex output
      character*1 gs(2*lzsm)
   91 format (48a1)
   92 format (' /Width ',i4,' /Height ',i4)
   93 format (' /ImageMatrix [',i4,' 0 0 ',i4,' 0 0]')
c define color space
      write (21,*) '[/Indexed /DeviceRGB 255 <'
c location of base characters in hex conversion
      i0 = ichar('0')
      ia = ichar('A') - 10
c write palette
      do 30 j = 1, 32
      ioff = 8*(j - 1)
      do 20 i = 1, 8
      iioff = 3*(i - 1)
      do 10 ii = 1, 3
      il = 255.*vpal(ii,i+ioff)
c ih = 4 high bits
      ih = il/16
c il = 4 low bits
      il = il - 16*ih
c find ascii code for high bits
      if (ih.gt.9) then
         ih = ih + ia
      else
         ih = ih + i0
      endif
c find ascii code for low bits
      if (il.gt.9) then
         il = il + ia
      else
         il = il + i0
      endif
c convert to printable character
      gs(2*(ii+iioff)-1) = char(ih)
      gs(2*(ii+iioff)) = char(il)
   10 continue
   20 continue
c write characters in hex
      write (21,91) (gs(i),i=1,48)
   30 continue
      write (21,*) '> ] setcolorspace'
c define image dictionary
      write (21,*) '<<'
      write (21,*) ' /ImageType 1'
      write (21,92) lx, ly
      write (21,*) ' /BitsPerComponent 8'
      write (21,*) ' /Decode [0 255]'
      write (21,93) lx, ly
      write (21,*) ' /DataSource currentfile /ASCIIHexDecode filter'
      write (21,*) '>>'
      write (21,*) 'image'
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
c optimized for postscript printer library
c npald = number of palette entries
      parameter(npald=256)
c lzsm =  group size for printing hex characters
      parameter(lzsm=40)
      parameter(maxt=3)
c ifrg = index of foreground color
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c iwt = current workstation type
c cgrs = current grayscale
c vpal(i,j+1) = ith component of jth palette entry (i=4 is grayscale)
      common /gkspsp/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,
     1ltc,mtc,infs,itxtp,nft,nfts,alws,ams,chx,cch,cgrs,wwnvp,dvt,trans,
     2chuv,vpal
      dimension wwnvp(8), dvt(4), trans(8,maxt), chuv(2), vpal(4,256)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      character*1 image(lenb)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c gs = temporary array for hex output
      character*1 gs(2*lzsm)
c lbl = character variable for numerical label
      character*8 lbl
      save zero,csize,alx
   91 format(1x,i3,' scalefont setfont')
   92 format ('%%Page: ? ',i4)
   93 format (i3,1x,i3,' translate')
   94 format (i3,1x,i3,' scale')
   95 format (i3,1x,i3,1x,i1,' [',i4,' 0 0 ',i3,1x,i3,' 0]')
   96 format ('{currentfile ',i3,' string readhexstring pop} image')
   97 format (80a1)
   98 format (f8.6,1x,f8.6,' scale')
   99 format(1x,f4.2,' setgray')
   81 format ('{currentfile ',i3,' string readhexstring pop} false 3 col
     1orimage')
   82 format(1x,3(f4.2,1x),'setrgbcolor')
   83 format ('%%BeginDocument')
   84 format ('%!PS-Adobe-2.1 EPSF-2.0')
   85 format ('%%BoundingBox: 0 0 612 792')
   86 format ('%%Extensions: CMYK')
   87 format ('%%EndDocument')
   88 format ('%!PS-Adobe-3.0 EPSF-3.0')
   89 format ('%%LanguageLevel: 2')
c csize = vertical size of characters
c alx = x coordinate for numerical labels
      data zero,csize,alx /0.,.02,.86/
c nbit = the number of colors, pixel depth
      nbit = 8/npix
c lx = size of unreduced image
      lx = (8*lz)/nbit
c normalize characters and location
      lafx = alx*real(isx)
      lafy = zero*real(isy)
      nfs = csize*real(isy)
c monochrome output
      if ((iwt.eq.1).or.(iwt.eq.2)) then
         mbit = 1
c grey-scale or color output if possible
      else
         mbit = min(8,nbit)
      endif
      mpix = nbit/mbit
c nbc = number of bytes per color
      nbc = 1
c color postscript colorimage command does not use a palette
      if ((mbit.eq.8).and.((iwt.eq.5).or.(iwt.eq.6))) nbc = 3
c l2 = (0,1) = (no,yes) use level2 commands
      l2 = 0
      if ((mbit.eq.8).and.((iwt.eq.7).or.(iwt.eq.8))) l2 = 1
c ntc = number of colors possible
      ntc = 2**nbit
c ls = shift operator
      ls = 2**(8/mpix)
c lzo = number of bytes in output scan line
      lzo = (lz - 1)/mpix + 1
c process lines in groups of lzs, default is 40
      lzs = lzsm/nbc
      if (lzs.gt.lzo) lzs = lzo
c lg = number of full or partial groups of lzs
      lg = (lzo - 1)/lzs + 1
c lgr = size of last group
      lgr = lzo - lzs*(lg - 1)
c location of base characters in hex conversion
      i0 = ichar('0')
      ia = ichar('A') - 10
c set scale and origin in landscape mode
      if (iwt.eq.(2*(iwt/2))) then
         lsm = 1
         iscx = 720
         iscy = 540
         igx = 576
         igy = 756 - iscx
         scl = 1.0
c set scale and origin in portrait mode
      else
         lsm = 0
         iscx = 528
         iscy = 396
         igx = 36
         igy = 756 - iscy
         scl = 11./15.
      endif
c restore and resave context
      write (21,*) 'grestore restore save'
c write out new page comments
      write (21,92) nf
c begin encapsulated document
      write (21,83)
      if (l2.eq.1) then
         write (21,88)
      else
         write (21,84)
      endif
      write (21,85)
c flag colorimage command as a level 1 extension
      if (nbc.eq.3) write (21,86)
      if (l2.eq.1) write (21,89)
c write out origin
      write (21,93) igx, igy
c rotate if landscape
      if (lsm.eq.1) write (21,*) '90 rotate'
c save graphics context
      write (21,*) 'gsave'
c write out scale
      write (21,94) lx, ly
c write out image operator command
      if (l2.eq.1) then
         call indimg(lx,ly)
      else
c write out arguments for image operator
         write (21,95) lx, ly, mbit, lx, ly, 0
c grey-scale image
         if (nbc.eq.1) then
            write (21,96) lzo
c color image
         else
            write (21,81) lzo
         endif
      endif
c outer loop over number of lines
      do 60 k = 1, ly
      joff = lz*(ly - k)
      lzb = lzs
c loop over number of lzs blocks in row
      do 50 j = 1, lg
      ioff = lzs*(j - 1)
c adjust size of last group
      if (j.eq.lg) lzb = lgr
c loop over lzs items in row
      do 40 i = 1, lzb
      iioff = nbc*(i - 1)
      do 30 ii = 1, nbc
c no packing required
      if (mpix.eq.1) then
         it = ichar(image(i+ioff+joff))
         if (nbit.eq.8) then
            if (nbc.eq.1) then
c map to grey-scale
               if ((iwt.eq.3).or.(iwt.eq.4)) then
                  il = 255.*vpal(4,it+1)
c map to indexed color image
               else
                  il = it
               endif
c map to color image
            else
               il = 255.*vpal(ii,it+1)
            endif
         else
            il = 255 - it
         endif
c pack pixels into bytes
      else
         moff = mpix*(i + ioff - 1) + joff
         is = 0
         do 20 m = 1, mpix
         it = ichar(image(m+moff))
c one pixel per byte
         if (npix.eq.1) then
c use lookup table
            if (lupt.ne.0) it = ipal(it+1)
            is = 2*is + it
c multiple pixels per byte
         else
            its = 0
            lfs = 1
c unpack pixels from byte
            do 10 n = 1, npix
            itt = it/ntc
            it = it - itt*ntc
            if (lupt.ne.0) it = ipal(it+1)
            its = its + lfs*it
            it = itt
            lfs = lfs + lfs
   10       continue
            is = lfs*is + its
         endif
   20    continue
c invert image
         il = 255 - is
      endif
c ih = 4 high bits
      ih = il/16
c il = 4 low bits
      il = il - 16*ih
c find ascii code for high bits
      if (ih.gt.9) then
         ih = ih + ia
      else
         ih = ih + i0
      endif
c find ascii code for low bits
      if (il.gt.9) then
         il = il + ia
      else
         il = il + i0
      endif
c convert to printable character
      gs(2*(ii+iioff)-1) = char(ih)
      gs(2*(ii+iioff)) = char(il)
   30 continue
   40 continue
c write characters in hex
      write (21,97) (gs(i),i=1,2*nbc*lzb)
   50 continue
   60 continue
c add EOF code for ASCIIHexDecode filter
      if (l2.eq.1) write (21,*) '>'
c add label
c restore and resave context
      write (21,*) 'grestore gsave'
c write out scale
      write (21,98) scl, scl
c set current font
      write (21,*) '/Courier findfont'
c first find how many digits in nf
      id = 0
      n = nf
   70 id = id + 1
      n = n/10
      if (n.gt.0) go to 70
c create label template
      lbl = '#        '
c create left justfied label
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
c scale and select font
      write (21,91) nfs
c set text color
      it = ifrg + 1
      cgrs = vpal(4,it)
c grayscale postscript format
      if ((iwt.eq.3).or.(iwt.eq.4)) then
         write (21,99) cgrs
c color postscript format
      elseif ((iwt.ge.5).and.(iwt.le.8)) then
         write (21,82) vpal(1,it), vpal(2,it), vpal(3,it)
      endif
c draw character string
      write (21,*) '('//lbl//')', lafx, lafy, ' moveto show'
c end endcapsulated document
      write (21,87)
c write image
      write (21,*) 'showpage'
      return
      end
      integer function kywait()
c special function to request keystroke
c returns: ascii code for keystroke
      kywait = 0
      return
      end
