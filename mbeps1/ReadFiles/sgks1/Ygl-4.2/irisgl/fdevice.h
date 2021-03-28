C     -*- fortran -*-
C     
C     Ygl: Run GL programs with standard X11 and/or OpenGL routines.
C     (C) Fred Hucht 1993-2006
C     EMail: fred<at>thp.Uni-Duisburg.de
C
C     $Id: Yfgl.h,v 4.3 2005-02-08 16:31:04+01 fred Exp fred $

C     FORTRAN includes for every routine
C     #ifndef	_YFGL_INCLUDED_
C     #define _YFGL_INCLUDED_
      
C     The colors
      integer*4  BLACK
      parameter( BLACK  = 0 )
      integer*4  WHITE
      parameter( WHITE  = 1 )
      integer*4  GREEN
      parameter( GREEN  = 2 )
      integer*4  YELLOW
      parameter( YELLOW = 3 )
      integer*4  BLUE
      parameter( BLUE   = 4 )
      integer*4  MAGENT
      parameter( MAGENT = 5 )
      integer*4  CYAN
      parameter( CYAN   = 6 )
      integer*4  RED
      parameter( RED    = 7 )
      
C     The Keys
      integer*4  BREAKK
      parameter( BREAKK = 1 )
      integer*4  SETUPK
      parameter( SETUPK = 2 )
      integer*4  LEFTCT
      parameter( LEFTCT = 3 )
      integer*4  CAPSLO
      parameter( CAPSLO = 4 )
      integer*4  RIGHTS
      parameter( RIGHTS = 5 )
      integer*4  LEFTSH
      parameter( LEFTSH = 6 )
      integer*4  ESCKEY
      parameter( ESCKEY = 7 )
      integer*4  ONEKEY
      parameter( ONEKEY = 8 )
      integer*4  TABKEY
      parameter( TABKEY = 9 )
      integer*4  QKEY
      parameter( QKEY   = 10 )
      integer*4  AKEY
      parameter( AKEY   = 11 )
      integer*4  SKEY
      parameter( SKEY   = 12 )
      integer*4  NOSCRL
      parameter( NOSCRL = 13 )
      integer*4  TWOKEY
      parameter( TWOKEY = 14 )
      integer*4  THREEK
      parameter( THREEK = 15 )
      integer*4  WKEY
      parameter( WKEY   = 16 )
      integer*4  EKEY
      parameter( EKEY   = 17 )
      integer*4  DKEY
      parameter( DKEY   = 18 )
      integer*4  FKEY
      parameter( FKEY   = 19 )
      integer*4  ZKEY
      parameter( ZKEY   = 20 )
      integer*4  XKEY
      parameter( XKEY   = 21 )
      integer*4  FOURKE
      parameter( FOURKE = 22 )
      integer*4  FIVEKE
      parameter( FIVEKE = 23 )
      integer*4  RKEY
      parameter( RKEY   = 24 )
      integer*4  TKEY
      parameter( TKEY   = 25 )
      integer*4  GKEY
      parameter( GKEY   = 26 )
      integer*4  HKEY
      parameter( HKEY   = 27 )
      integer*4  CKEY
      parameter( CKEY   = 28 )
      integer*4  VKEY
      parameter( VKEY   = 29 )
      integer*4  SIXKEY
      parameter( SIXKEY = 30 )
      integer*4  SEVENK
      parameter( SEVENK = 31 )
      integer*4  YKEY
      parameter( YKEY   = 32 )
      integer*4  UKEY
      parameter( UKEY   = 33 )
      integer*4  JKEY
      parameter( JKEY   = 34 )
      integer*4  KKEY
      parameter( KKEY   = 35 )
      integer*4  BKEY
      parameter( BKEY   = 36 )
      integer*4  NKEY
      parameter( NKEY   = 37 )
      integer*4  EIGHTK
      parameter( EIGHTK = 38 )
      integer*4  NINEKE
      parameter( NINEKE = 39 )
      integer*4  IKEY
      parameter( IKEY   = 40 )
      integer*4  OKEY
      parameter( OKEY   = 41 )
      integer*4  LKEY
      parameter( LKEY   = 42 )
      integer*4  SEMICO
      parameter( SEMICO = 43 )
      integer*4  MKEY
      parameter( MKEY   = 44 )
      integer*4  COMMAK
      parameter( COMMAK = 45 )
      integer*4  ZEROKE
      parameter( ZEROKE = 46 )
      integer*4  MINUSK
      parameter( MINUSK = 47 )
      integer*4  PKEY
      parameter( PKEY   = 48 )
      integer*4  LEFTBR
      parameter( LEFTBR = 49 )
      integer*4  QUOTEK
      parameter( QUOTEK = 50 )
      integer*4  RETKEY
      parameter( RETKEY = 51 )
      integer*4  PERIOD
      parameter( PERIOD = 52 )
      integer*4  VIRGUL
      parameter( VIRGUL = 53 )
      integer*4  EQUALK
      parameter( EQUALK = 54 )
      integer*4  ACCENT
      parameter( ACCENT = 55 )
      integer*4  RIGHTB
      parameter( RIGHTB = 56 )
      integer*4  BACKSL
      parameter( BACKSL = 57 )
      integer*4  PAD1
      parameter( PAD1   = 58 )
      integer*4  PAD0
      parameter( PAD0   = 59 )
      integer*4  LINEFE
      parameter( LINEFE = 60 )
      integer*4  BACKSP
      parameter( BACKSP = 61 )
      integer*4  DELKEY
      parameter( DELKEY = 62 )
      integer*4  PAD4
      parameter( PAD4   = 63 )
      integer*4  PAD2
      parameter( PAD2   = 64 )
      integer*4  PAD3
      parameter( PAD3   = 65 )
      integer*4  PADPER
      parameter( PADPER = 66 )
      integer*4  PAD7
      parameter( PAD7   = 67 )
      integer*4  PAD8
      parameter( PAD8   = 68 )
      integer*4  PAD5
      parameter( PAD5   = 69 )
      integer*4  PAD6
      parameter( PAD6   = 70 )
      integer*4  PADPF2
      parameter( PADPF2 = 71 )
      integer*4  PADPF1
      parameter( PADPF1 = 72 )
      integer*4  LEFTAR
      parameter( LEFTAR = 73 )
      integer*4  DOWNAR
      parameter( DOWNAR = 74 )
      integer*4  PAD9
      parameter( PAD9   = 75 )
      integer*4  PADMIN
      parameter( PADMIN = 76 )
      integer*4  PADCOM
      parameter( PADCOM = 77 )
      integer*4  PADPF4
      parameter( PADPF4 = 78 )
      integer*4  PADPF3
      parameter( PADPF3 = 79 )
      integer*4  RIGHTA
      parameter( RIGHTA = 80 )
      integer*4  UPARRO
      parameter( UPARRO = 81 )
      integer*4  PADENT
      parameter( PADENT = 82 )
      integer*4  SPACEK
      parameter( SPACEK = 83 )
      integer*4  LEFTAL
      parameter( LEFTAL = 143 )
      integer*4  RGHTAL
      parameter( RGHTAL = 144 )
      integer*4  RIGHTC
      parameter( RIGHTC = 145 )
      integer*4  F1KEY
      parameter( F1KEY  = 146 )
      integer*4  F2KEY
      parameter( F2KEY  = 147 )
      integer*4  F3KEY
      parameter( F3KEY  = 148 )
      integer*4  F4KEY
      parameter( F4KEY  = 149 )
      integer*4  F5KEY
      parameter( F5KEY  = 150 )
      integer*4  F6KEY
      parameter( F6KEY  = 151 )
      integer*4  F7KEY
      parameter( F7KEY  = 152 )
      integer*4  F8KEY
      parameter( F8KEY  = 153 )
      integer*4  F9KEY
      parameter( F9KEY  = 154 )
      integer*4  F10KEY
      parameter( F10KEY = 155 )
      integer*4  F11KEY
      parameter( F11KEY = 156 )
      integer*4  F12KEY
      parameter( F12KEY = 157 )
      integer*4  PRINTS
      parameter( PRINTS = 158 )
      integer*4  SCROLL
      parameter( SCROLL = 159 )
      integer*4  PAUSEK
      parameter( PAUSEK = 160 )
      integer*4  INSERT
      parameter( INSERT = 161 )
      integer*4  HOMEKE
      parameter( HOMEKE = 162 )
      integer*4  PAGEUP
      parameter( PAGEUP = 163 )
      integer*4  ENDKEY
      parameter( ENDKEY = 164 )
      integer*4  PAGEDO
      parameter( PAGEDO = 165 )
      integer*4  NUMLOC
      parameter( NUMLOC = 166 )
      integer*4  PADVIR
      parameter( PADVIR = 167 )
      integer*4  PADAST
      parameter( PADAST = 168 )
      integer*4  PADPLU
      parameter( PADPLU = 169 )
      
C     The Mouse stuff
      integer*4  MOUSE1
      parameter( MOUSE1 = 101 )
      integer*4  MOUSE2
      parameter( MOUSE2 = 102 )
      integer*4  MOUSE3
      parameter( MOUSE3 = 103 )
      integer*4  LEFTMO
      parameter( LEFTMO = 103 )
      integer*4  MIDDLE
      parameter( MIDDLE = 102 )
      integer*4  RIGHTM
      parameter( RIGHTM = 101 )
      integer*4  MENUBU
      parameter( MENUBU = 101 )
      
      integer*4  WHEELU
      parameter( WHEELU = 200 )
      integer*4  WHEELD
      parameter( WHEELD = 201 )
      
      integer*4  MOUSEX
      parameter( MOUSEX = 266 )
      integer*4  MOUSEY
      parameter( MOUSEY = 267 )
      
C     Other Devices
      integer*4  KEYBD
      parameter( KEYBD  = 513 )
      integer*4  REDRAW
      parameter( REDRAW = 528 )
      integer*4  INPUTC
      parameter( INPUTC = 534 )
C     Device WINCLO conflict with routine winclo()
      integer*4  WINCLOSE
      parameter( WINCLOSE = 537 )
      integer*4  WINFRE
      parameter( WINFRE = 539 )
      integer*4  WINTHA
      parameter( WINTHA = 540 )
      integer*4  WINQUI
      parameter( WINQUI = 542 )
      integer*4  DEPTHC
      parameter( DEPTHC = 543 )
      integer*4  SRCAUT
      parameter( SRCAUT = 0 )
      integer*4  SRCFRO
      parameter( SRCFRO = 1 )
      integer*4  SRCBAC
      parameter( SRCBAC = 2 )
      integer*4  DMRGB
      parameter( DMRGB  = 0 )
      integer*4  DMSING
      parameter( DMSING = 1 )
      integer*4  DMDOUB
      parameter( DMDOUB = 2 )
      integer*4  DMRGBD
      parameter( DMRGBD = 5 )
      integer*4  GDXPMA
      parameter( GDXPMA = 1 )
      integer*4  GDYPMA
      parameter( GDYPMA = 2 )
      integer*4  PUPNON
      parameter( PUPNON = 0 )
      integer*4  PUPGRE
      parameter( PUPGRE = 1 )
      integer*4  XMAXSC
      parameter( XMAXSC = 1279 )
      integer*4  YMAXSC
      parameter( YMAXSC = 1023 )

C     The logicops
      integer*4  LOZERO
      parameter( LOZERO = 0 )
      integer*4  LOAND
      parameter( LOAND = 1 )
      integer*4  LOANDR
      parameter( LOANDR = 2 )
      integer*4  LOSRC
      parameter( LOSRC = 3 )
      integer*4  LOANDI
      parameter( LOANDI = 4 )
      integer*4  LODST
      parameter( LODST = 5 )
      integer*4  LOXOR
      parameter( LOXOR = 6 )
      integer*4  LOOR
      parameter( LOOR = 7 )
      integer*4  LONOR
      parameter( LONOR = 8 )
      integer*4  LOXNOR
      parameter( LOXNOR = 9 )
      integer*4  LONDST
      parameter( LONDST = 10 )
      integer*4  LOORR
      parameter( LOORR = 11 )
      integer*4  LONSRC
      parameter( LONSRC = 12 )
      integer*4  LOORI
      parameter( LOORI = 13 )
      integer*4  LONAND
      parameter( LONAND = 14 )
      integer*4  LOONE
      parameter( LOONE = 15 )
      integer*4  LOMIN
      parameter( LOMIN = 16 )
      integer*4  LOMAX
      parameter( LOMAX = 17 )
      integer*4  LOAVG
      parameter( LOAVG = 18 )
      integer*4  LODMS
      parameter( LODMS = 19 )
      integer*4  LOSMD
      parameter( LOSMD = 20 )
      integer*4  LOSUM
      parameter( LOSUM = 21 )
      
C     The matrix modes
      integer*4   MSINGL
      parameter ( MSINGL = 0 )
      integer*4   MPROJE
      parameter ( MPROJE = 1 )
      integer*4   MVIEWI
      parameter ( MVIEWI = 2 )

C     The blend functions
      integer*4   BFZERO
      parameter ( BFZERO = 0 )
      integer*4   BFON
      parameter ( BFONE = 1 )
      integer*4   BFDC
      parameter ( BFDC = 2 )
      integer*4   BFSC
      parameter ( BFSC = 2 )
      integer*4   BFMDC
      parameter ( BFMDC = 3 )
      integer*4   BFMSC
      parameter ( BFMSC = 3 )
      integer*4   BFSA
      parameter ( BFSA = 4 )
      integer*4   BFMSA
      parameter ( BFMSA = 5 )
      integer*4   BFDA
      parameter ( BFDA = 6 )
      integer*4   BFMDA
      parameter ( BFMDA = 7 )
      
C     The z functions
      integer*4   ZFNEVE
      parameter ( ZFNEVE = 0 )
      integer*4   ZFLESS
      parameter ( ZFLESS = 1 )
      integer*4   ZFEQUA
      parameter ( ZFEQUA = 2 )
      integer*4   ZFLEQU
      parameter ( ZFLEQU = 3 )
      integer*4   ZFGREA
      parameter ( ZFGREA = 4 )
      integer*4   ZFNOTE
      parameter ( ZFNOTE = 5 )
      integer*4   ZFGEQU
      parameter ( ZFGEQU = 6 )
      integer*4   ZFALWA
      parameter ( ZFALWA = 7 )

C     The lighting stuff
      real*4      LMNULL
      parameter ( LMNULL = 0.0 )
      integer*4   MAXLIG
      parameter ( MAXLIG = 8 )
      integer*4   MAXRES
      parameter ( MAXRES = 4 )
      
      integer*4   DEFMAT
      parameter ( DEFMAT = 0 )
      integer*4   EMISSI
      parameter ( EMISSI = 1 )
      integer*4   AMBIEN
      parameter ( AMBIEN = 2 )
      integer*4   DIFFUS
      parameter ( DIFFUS = 3 )
      integer*4   SPECUL
      parameter ( SPECUL = 4 )
      integer*4   SHININ
      parameter ( SHININ = 5 )
      integer*4   COLORI
      parameter ( COLORI = 6 )
      integer*4   ALPHA
      parameter ( ALPHA = 7 )
      
      integer*4   DEFLIG
      parameter ( DEFLIG = 100 )
      integer*4   LCOLOR
      parameter ( LCOLOR = 101 )
      integer*4   POSITI
      parameter ( POSITI = 102 )
      integer*4   SPOTDI
      parameter ( SPOTDI = 103 )
      integer*4   SPOTLI
      parameter ( SPOTLI = 104 )
      
      integer*4   DEFLMO
      parameter ( DEFLMO = 200 )
      integer*4   LOCALV
      parameter ( LOCALV = 201 )
      integer*4   ATTENU
      parameter ( ATTENU = 202 )
      integer*4   ATTEN2
      parameter ( ATTEN2 = 203 )
      integer*4   TWOSID
      parameter ( TWOSID = 204 )
      
      integer*4   MATERI
      parameter ( MATERI = 1000 )
      integer*4   LIGHT0
      parameter ( LIGHT0 = 1100 )
      integer*4   LIGHT1
      parameter ( LIGHT1 = 1101 )
      integer*4   LIGHT2
      parameter ( LIGHT2 = 1102 )
      integer*4   LIGHT3
      parameter ( LIGHT3 = 1103 )
      integer*4   LIGHT4
      parameter ( LIGHT4 = 1104 )
      integer*4   LIGHT5
      parameter ( LIGHT5 = 1105 )
      integer*4   LIGHT6
      parameter ( LIGHT6 = 1106 )
      integer*4   LIGHT7
      parameter ( LIGHT7 = 1107 )
      integer*4   LMODEL
      parameter ( LMODEL = 1200 )
      
      integer*4   LMCCOL
      parameter ( LMCCOL = 0 )
      integer*4   LMCEMI
      parameter ( LMCEMI = 1 )
      integer*4   LMCAMB
      parameter ( LMCAMB = 2 )
      integer*4   LMCDIF
      parameter ( LMCDIF = 3 )
      integer*4   LMCSPE
      parameter ( LMCSPE = 4 )
      integer*4   LMCAD
      parameter ( LMCAD = 5 )
      integer*4   LMCNUL
      parameter ( LMCNUL = 6 )
      
C     The routines
      integer*4 isqueu
      external  isqueu
      integer*4 qtest
      external  qtest
      integer*4 qread
      external  qread
C     getXdpy and getXwid are >6 chars long
      integer*4 getxdpy
      external  getxdpy
      integer*4 getxwid
      external  getxwid
      integer*4 getxdid
      external  getxdid
      integer*4 getxgc
      external  getxgc
      integer*4 winget
      external  winget
      integer*4 getpla
      external  getpla
      integer*4 getval
      external  getval
      integer*4 getbut
      external  getbut
      integer*4 gversi
      external  gversi
      integer*4 windep
      external  windep
      integer*4 getlwi
      external  getlwi
      integer*4 getlst
      external  getlst
      integer*4 getlsr
      external  getlsr
      integer*4 getdis
      external  getdis
      integer*4 getgde
      external  getgde
      integer*4 getfon
      external  getfon
      integer*4 gethei
      external  gethei
      integer*4 getdes
      external  getdes
      integer*4 strwid
      external  strwid
      integer*4 getcol
      external  getcol
      integer*4 crectr
      external  crectr
      integer*4 rectre
      external  rectre
      integer*4 lrectr
      external  lrectr
      integer*4 readpi
      external  readpi
      integer*4 readrg
      external  readrg
      integer*4 dopup
      external  dopup
      integer*4 newpup
      external  newpup
      integer*4 winope
      external  winope
      integer*4 swinop
      external  swinop
      integer*4 winx
      external  winx
      integer*4 gl2ppm
      external  gl2ppm
      integer*4 getmmo
      external  getmmo
      integer*4 genobj
      external  genobj
      integer*4 isobj
      external  isobj
      integer*4 getope
      external  getope
      integer*4 endpic
      external  endpic
C
C     Not implemented (yet)
C
C     integer*4 blkqre
C     external  blkqre
C     integer*4 endsel
C     external  endsel
C     integer*4 gentag
C     external  gentag
C     integer*4 getbac
C     external  getbac
C     integer*4 getbuf
C     external  getbuf
C     integer*4 getcmm
C     external  getcmm
C     integer*4 getdcm
C     external  getdcm
C     integer*4 getdra
C     external  getdra
C     integer*4 getmap
C     external  getmap
C     integer*4 getpat
C     external  getpat
C     integer*4 getsha
C     external  getsha
C     integer*4 getsm
C     external  getsm
C     integer*4 getwri
C     external  getwri
C     integer*4 getzbu
C     external  getzbu
C     integer*4 istag
C     external  istag

C     #endif /* _YFGL_INCLUDED_ */
