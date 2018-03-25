/*
 *    Ygl: Run GL programs with standard X11 routines.
 *    (C) Fred Hucht 1993-2006
 *    EMail: fred@thp.Uni-Duisburg.de
 *
 *    $Id: Ygl.h,v 4.5 2005-02-08 17:01:45+01 fred Exp fred $
 */

#ifndef	_YGL_INCLUDED_
#define _YGL_INCLUDED_

#include <sys/types.h>

#ifdef _AUX_SOURCE
# include <X11/Yglprefix.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
  
#define BLACK			0
#define WHITE			1
#define GREEN			2
#define YELLOW			3
#define BLUE			4
#define MAGENTA			5
#define CYAN			6
#define RED			7
  
  /* for queue */
#define NULLDEV			0
#define BREAKKEY		1
#define SETUPKEY		2
#define LEFTCTRLKEY		3
#define CAPSLOCKKEY		4
#define RIGHTSHIFTKEY		5
#define LEFTSHIFTKEY		6
#define ESCKEY			7
#define ONEKEY			8
#define TABKEY			9
#define QKEY			10
#define AKEY			11
#define SKEY			12
#define NOSCRLKEY		13
#define TWOKEY			14
#define THREEKEY		15
#define WKEY			16
#define EKEY			17
#define DKEY			18
#define FKEY			19
#define ZKEY			20
#define XKEY			21
#define FOURKEY			22
#define FIVEKEY			23
#define RKEY			24
#define TKEY			25
#define GKEY			26
#define HKEY			27
#define CKEY			28
#define VKEY			29
#define SIXKEY			30
#define SEVENKEY		31
#define YKEY			32
#define UKEY			33
#define JKEY			34
#define KKEY			35
#define BKEY			36
#define NKEY			37
#define EIGHTKEY		38
#define NINEKEY			39
#define IKEY			40
#define OKEY			41
#define LKEY			42
#define SEMICOLONKEY		43
#define MKEY			44
#define COMMAKEY		45
#define ZEROKEY			46
#define MINUSKEY		47
#define PKEY			48
#define LEFTBRACKETKEY		49
#define QUOTEKEY		50
#define RETKEY			51
#define PERIODKEY		52
#define VIRGULEKEY		53
#define EQUALKEY		54
#define ACCENTGRAVEKEY		55
#define RIGHTBRACKETKEY		56
#define BACKSLASHKEY		57
#define PAD1			58
#define PAD0			59
#define LINEFEEDKEY		60
#define BACKSPACEKEY		61
#define DELKEY			62
#define PAD4			63
#define PAD2			64
#define PAD3			65
#define PADPERIOD		66
#define PAD7			67
#define PAD8			68
#define PAD5			69
#define PAD6			70
#define PADPF2			71
#define PADPF1			72
#define LEFTARROWKEY		73
#define DOWNARROWKEY		74
#define PAD9			75
#define PADMINUS		76
#define PADCOMMA		77
#define PADPF4			78
#define PADPF3			79
#define RIGHTARROWKEY		80
#define UPARROWKEY		81
#define PADENTER		82
#define SPACEKEY		83
#define LEFTALTKEY	 	143
#define	RIGHTALTKEY	 	144
#define	RIGHTCTRLKEY	 	145
#define	F1KEY	 		146
#define	F2KEY	 		147
#define	F3KEY	 		148
#define	F4KEY	 		149
#define	F5KEY	 		150
#define	F6KEY	 		151
#define	F7KEY	 		152
#define	F8KEY	 		153
#define	F9KEY	 		154
#define	F10KEY			155
#define	F11KEY			156
#define	F12KEY			157
#define	PRINTSCREENKEY		158
#define	SCROLLLOCKKEY		159
#define	PAUSEKEY		160
#define	INSERTKEY		161
#define	HOMEKEY			162
#define	PAGEUPKEY	 	163
#define	ENDKEY			164
#define	PAGEDOWNKEY		165
#define	NUMLOCKKEY		166
#define	PADVIRGULEKEY	 	167
#define PADASTERKEY	 	168
#define PADPLUSKEY	 	169
  
#define MOUSE1			101
#define MOUSE2			102
#define MOUSE3			103
#define LEFTMOUSE		103
#define MIDDLEMOUSE		102
#define RIGHTMOUSE		101
#define MENUBUTTON		101
  
  /* Wheel mouse (not in GL) */
#define WHEELUP			200
#define WHEELDOWN		201
  
#define MOUSEX			266
#define MOUSEY			267
  
#define ANYKEY			512
#define KEYBD			513
#define TIMER0			515
#define TIMER1			516
#define TIMER2			517
#define TIMER3			518
#define REDRAW			528
#define INPUTCHANGE		534
#define WINCLOSE		537
#define WINFREEZE		539
#define WINTHAW			540
#define WINQUIT			542
#define DEPTHCHANGE		543
  
#define MAXYGLDEVICE		544
  
  /* for readsource(): */
#define SRC_AUTO		0
#define SRC_FRONT		1
#define SRC_BACK		2
  
  /* for getdisplaymode(): */
#define DMRGB			0L
#define DMSINGLE		1L
#define DMDOUBLE		2L
#define DMRGBDOUBLE		5L
  
  /* for getgdesc(): */
#define GD_XPMAX 		1L
#define GD_YPMAX 		2L
  
  /* {XY}MAXSCREEN are dynamic */
#define XMAXSCREEN (getgdesc(GD_XPMAX) - 1)
#define YMAXSCREEN (getgdesc(GD_YPMAX) - 1)

  /* for setpup() */
#define PUP_NONE		0
#define PUP_GREY		1
  
  /* for logicop() */
#define LO_ZERO			0x0
#define LO_AND			0x1
#define LO_ANDR			0x2
#define LO_SRC			0x3
#define LO_ANDI			0x4
#define LO_DST			0x5
#define LO_XOR			0x6
#define LO_OR			0x7
#define LO_NOR			0x8
#define LO_XNOR			0x9
#define LO_NDST			0xa
#define LO_ORR			0xb
#define LO_NSRC			0xc
#define LO_ORI			0xd
#define LO_NAND			0xe
#define LO_ONE			0xf
#define LO_MIN			0x10
#define LO_MAX			0x11
#define LO_AVG			0x12
#define LO_DMS			0x13
#define LO_SMD			0x14
#define LO_SUM			0x15
  
  /* for mmode() */
#define MSINGLE			0
#define MPROJECTION		1
#define MVIEWING		2

  /* for blendfunction() */
#define BF_ZERO			0
#define BF_ONE			1
#define BF_SC			2
#define BF_MSC			3
#define BF_SA			4
#define BF_MSA			5
#define BF_DA			6
#define BF_MDA			7
#define BF_DC			8
#define BF_MDC			9
#define BF_MIN_SA_MDA		10
  
  /* for zfunction() */
#define ZF_NEVER		0
#define ZF_LESS			1
#define ZF_EQUAL		2
#define ZF_LEQUAL		3
#define ZF_GREATER		4
#define ZF_NOTEQUAL		5
#define ZF_GEQUAL		6
#define ZF_ALWAYS		7

  /* for lmdef(): MATERIAL properties */
#define DEFMATERIAL		0
#define EMISSION		1
#define AMBIENT			2
#define DIFFUSE			3
#define SPECULAR		4
#define SHININESS		5
#define COLORINDEXES		6
#define ALPHA			7
  /* for lmdef(): LIGHT properties */
#define DEFLIGHT		100
#define LCOLOR			101
#define POSITION		102
#define SPOTDIRECTION		103
#define SPOTLIGHT		104

  /* LIGHTINGMODEL properties */
#define DEFLMODEL		200
#define LOCALVIEWER		201
#define ATTENUATION		202
#define ATTENUATION2		203 /* used by SGI */
#define TWOSIDE			204

  /* TARGET constants */
#define MATERIAL		1000
#define BACKMATERIAL		1001
#define LIGHT0			1100
#define LIGHT1			1101
#define LIGHT2			1102
#define LIGHT3			1103
#define LIGHT4			1104
#define LIGHT5			1105
#define LIGHT6			1106
#define LIGHT7			1107
#define LMODEL			1200

  /* for lmcolor(): modes */
#define LMC_COLOR       0
#define LMC_EMISSION    1
#define LMC_AMBIENT     2
#define LMC_DIFFUSE     3
#define LMC_SPECULAR    4
#define LMC_AD          5
#define LMC_NULL        6

  /* for lmdef(): constants     */
#define LMNULL           0.0

  /* for shademodel() */
#define FLAT			0
#define GOURAUD			1

  /* Types */
#include <X11/Ygltypes.h>
  typedef char			Char8;
  typedef char			Void;
  
  typedef Uint8			Byte;
  typedef Uint8			RGBvalue;
  
  typedef Uint16		Colorindex;
  typedef Uint16		Device;
  typedef Uint16		Linestyle;
  
  typedef Int16			Angle;
  typedef Int16			Scoord;
  typedef Int16			Screencoord;
  
  typedef Int32			Icoord;
  
  typedef Float32		Coord;
  
  typedef Float32		Matrix[4][4];

  /********************* draw.c */
  extern void  clear		( void );
  
  /* Points */
  extern void  pnt2 		(  Coord,  Coord );
  extern void  pnt2i		( Icoord, Icoord );
  extern void  pnt2s		( Scoord, Scoord );
  
  /* Lines */
  extern void  move2 		(  Coord,  Coord );
  extern void  move2i		( Icoord, Icoord );
  extern void  move2s		( Scoord, Scoord );
  
  extern void  rmv2 		(  Coord,  Coord );
  extern void  rmv2i		( Icoord, Icoord );
  extern void  rmv2s		( Scoord, Scoord );
  
  extern void  draw2 		(  Coord,  Coord );
  extern void  draw2i		( Icoord, Icoord );
  extern void  draw2s		( Scoord, Scoord );
  
  extern void  rdr2 		(  Coord,  Coord );
  extern void  rdr2i		( Icoord, Icoord );
  extern void  rdr2s		( Scoord, Scoord );
  
  /* Arcs & Circles */
  extern void  arc 		(  Coord,  Coord,  Coord, Angle, Angle );
  extern void  arci		( Icoord, Icoord, Icoord, Angle, Angle );
  extern void  arcs		( Scoord, Scoord, Scoord, Angle, Angle );
  
  extern void  arcf 		(  Coord,  Coord,  Coord, Angle, Angle );
  extern void  arcfi		( Icoord, Icoord, Icoord, Angle, Angle );
  extern void  arcfs		( Scoord, Scoord, Scoord, Angle, Angle );
  
  extern void  circ 		(  Coord,  Coord,  Coord );
  extern void  circi		( Icoord, Icoord, Icoord );
  extern void  circs		( Scoord, Scoord, Scoord );
  
  extern void  circf 		(  Coord,  Coord,  Coord );
  extern void  circfi		( Icoord, Icoord, Icoord );
  extern void  circfs		( Scoord, Scoord, Scoord );
  
  /* Rects & Boxes */
  extern void  rect 		(  Coord,  Coord,  Coord,  Coord );
  extern void  recti		( Icoord, Icoord, Icoord, Icoord );
  extern void  rects		( Scoord, Scoord, Scoord, Scoord );
  
  extern void  rectf 		(  Coord,  Coord,  Coord,  Coord );
  extern void  rectfi		( Icoord, Icoord, Icoord, Icoord );
  extern void  rectfs		( Scoord, Scoord, Scoord, Scoord );
  
  extern void  sbox 		(  Coord,  Coord,  Coord,  Coord );
  extern void  sboxi		( Icoord, Icoord, Icoord, Icoord );
  extern void  sboxs		( Scoord, Scoord, Scoord, Scoord );
  
  extern void  sboxf 		(  Coord,  Coord,  Coord,  Coord );
  extern void  sboxfi		( Icoord, Icoord, Icoord, Icoord );
  extern void  sboxfs		( Scoord, Scoord, Scoord, Scoord );
  
  /* Filled Polygons */
  extern void  concave		( Int32 );
  
  extern void  pmv2 		(  Coord,  Coord );
  extern void  pmv2i		( Icoord, Icoord );
  extern void  pmv2s		( Scoord, Scoord );
  
  extern void  rpmv2 		(  Coord,  Coord );
  extern void  rpmv2i		( Icoord, Icoord );
  extern void  rpmv2s		( Scoord, Scoord );
  
  extern void  pdr2 		(  Coord,  Coord );
  extern void  pdr2i		( Icoord, Icoord );
  extern void  pdr2s		( Scoord, Scoord );
  
  extern void  rpdr2 		(  Coord,  Coord );
  extern void  rpdr2i		( Icoord, Icoord );
  extern void  rpdr2s		( Scoord, Scoord );
  
  extern void  pclos		( void );
  
  extern void  poly2 		( Int32,  Coord[][2] );
  extern void  poly2i		( Int32, Icoord[][2] );
  extern void  poly2s		( Int32, Scoord[][2] );
  
  extern void  polf2 		( Int32,  Coord[][2] );
  extern void  polf2i		( Int32, Icoord[][2] );
  extern void  polf2s		( Int32, Scoord[][2] );
  
  /* Vertex graphics */
  extern void  bgnpoint		( void );
  extern void  bgnline		( void );
  extern void  bgnclosedline	( void );
  extern void  bgnpolygon	( void );
  extern void  bgntmesh 	( void );
  
  extern void  endpoint		( void );
  extern void  endline		( void );
  extern void  endclosedline	( void );
  extern void  endpolygon	( void );
  extern void  endtmesh 	( void );
  
  extern void  v2s		( Int16[2] );
  extern void  v2i		( Int32[2] );
  extern void  v2f		( Float32[2] );
  extern void  v2d		( Float64[2] );
  
  /* Text */
  extern void  cmov2 		(  Coord,  Coord );
  extern void  cmov2i		( Icoord, Icoord );
  extern void  cmov2s		( Scoord, Scoord );
  
  extern void  getcpos		( Screencoord *, Screencoord * );
  
  /* Extensions: Routines not in gl by MiSt (michael@thp.Uni-Duisburg.de) */
#ifdef X11
  extern void  arcx 		(  Coord,  Coord,  Coord,  Coord, Angle, Angle );
  extern void  arcxi		( Icoord, Icoord, Icoord, Icoord, Angle, Angle );
  extern void  arcxs		( Scoord, Scoord, Scoord, Scoord, Angle, Angle );
  
  extern void  arcxf 		(  Coord,  Coord,  Coord,  Coord, Angle, Angle );
  extern void  arcxfi		( Icoord, Icoord, Icoord, Icoord, Angle, Angle );
  extern void  arcxfs		( Scoord, Scoord, Scoord, Scoord, Angle, Angle );
#endif
  
  /********************* queue.c */
  extern void  tie		( Device, Device, Device );
  extern void  noise 		( Device, Int16 );
  extern Int32 isqueued 	( Int16 );
  extern void  qdevice		( Device );
  extern void  unqdevice	( Device );
  extern void  qreset		( void );
  extern Int32 qtest		( void );
  extern Int32 qread		( Int16 * );
  extern void  qenter		( Int16, Int16 );

  extern void pick 		( Int16 *, Int32 );
  extern Int32 endpick 		( Int16[] );
  extern void picksize 		( Int16, Int16 );
  
  extern void initnames 	( void );
  extern void loadname 		( Int16 );
  extern void pushname 		( Int16 );
  extern void popname 		( void );
  
  /********************* misc.c */
  extern void  singlebuffer 	( void );
  extern void  doublebuffer 	( void );
  extern void  swapbuffers 	( void );
  extern void  frontbuffer	( Int32 );
  extern void  backbuffer	( Int32 );
  
  extern void  gflush		( void );
  extern void  gsync		( void );
  
#ifdef _XLIB_H_ 		/* Declare if <X11/Xlib.h> is included */
  extern Display *getXdpy 	( void );
  extern Window   getXwid 	( void );
#ifdef X11
  extern Window   getXdid 	( void );
  extern GC       getXgc 	( void );
#endif
#endif /* _XLIB_H_ */
  
  extern void  wintitle		( Char8 * );
  extern void  winset		( Int32 );
  extern Int32 winget		( void );
  extern Int32 getplanes	( void );
  extern Int32 getvaluator 	( Device );
  extern Int32 getbutton	( Device );
  extern Int32 gversion		( Char8[12] );
  
  extern void  ortho2		( Coord, Coord, Coord, Coord );
  extern void  viewport		( Screencoord, Screencoord, Screencoord, Screencoord );
  extern void  getviewport	( Screencoord *, Screencoord *, Screencoord *, Screencoord * );
  extern void  reshapeviewport 	( void );
  extern void  pushviewport 	( void );
  extern void  popviewport 	( void );
  
  extern void  winpop		( void );
  extern void  winpush		( void );
  extern Int32 windepth		( Int32 );
  
  extern void  linewidth	( Int16 );
  extern Int32 getlwidth	( void );
  extern void  deflinestyle	( Int32, Linestyle );
  extern void  setlinestyle	( Int32 );
  extern Int32 getlstyle	( void );
  extern void  lsrepeat		( Int32 );
  extern Int32 getlsrepeat	( void );
  extern Int32 getdisplaymode	( void );
  
  extern void  setbell		( Char8 );
  extern void  ringbell		( void );
  
  extern Int32 getgdesc		( Int32 );
  
  extern void  foreground	( void );
  
  extern void  logicop		( Int32 );
  
  extern void  getmatrix	( Matrix );
  
  /********************* font.c */
  extern void  loadXfont	( Int32 , Char8 * );
  extern void  font		( Int16 );
  extern Int32 getfont		( void );
  extern void  getfontencoding 	( Char8 * );
  extern Int32 getheight		( void );
  extern Int32 getdescender 	( void );
  extern Int32 strwidth		( Char8 * );
  extern void  charstr		( Char8 * );
  
  /********************* color.c */
  extern void  mapcolor		( Colorindex, Int16, Int16, Int16 );
  extern void  RGBcolor		( Int16, Int16, Int16 );
  extern void  cpack		( Uint32 );
  
  extern void  c3s		( Int16[3] );
  extern void  c3i		( Int32[3] );
  extern void  c3f		( Float32[3] );
  
  extern Int32 getcolor		( void );
  extern void  getmcolor	( Colorindex, Int16 *, Int16 *, Int16 * );
  extern void  getmcolors	( Colorindex, Colorindex, Int16 *, Int16 *, Int16 * );
  extern void  gRGBcolor	( Int16 *, Int16 *, Int16 * );
  
  extern void  color		( Colorindex );
  extern void  readsource	( Int32 );
  
  extern void  rectzoom 	( Float32, Float32 );
  
  extern Int32 crectread	( Screencoord, Screencoord, Screencoord, Screencoord, Uint8 * );
  extern Int32 rectread		( Screencoord, Screencoord, Screencoord, Screencoord, Int16 * );
  extern Int32 lrectread	( Screencoord, Screencoord, Screencoord, Screencoord, Int32 * );
  
  extern void  crectwrite	( Screencoord, Screencoord, Screencoord, Screencoord, Uint8 * );
  extern void  rectwrite	( Screencoord, Screencoord, Screencoord, Screencoord, Int16 * );
  extern void  lrectwrite	( Screencoord, Screencoord, Screencoord, Screencoord, Int32 * );
  
  extern void  rectcopy		( Screencoord, Screencoord, Screencoord, Screencoord, Screencoord, Screencoord );
  
  extern Int32 readpixels	( Int16, Colorindex[] );
  extern void  writepixels	( Int16, Colorindex[] );
  extern Int32 readRGB		( Int16, RGBvalue[], RGBvalue[], RGBvalue[] );
  extern void  writeRGB		( Int16, RGBvalue[], RGBvalue[], RGBvalue[] );

  /* for/from Bill Bishop */
  extern void blendfunction 	( Int32, Int32 );
  
  /********************* menu.c */
  extern void  addtopup		( Int32, Char8 *, ... );
  extern Int32 defpup		( Char8 *, ... );
  extern Int32 dopup		( Int32 );
  extern void  freepup		( Int32 );
  extern Int32 newpup		( void );
  extern void  setpup		( Int32, Int32, Int32 );
  
  /********************* ygl.c */
  /* Contraints */
  extern void  minsize		( Int32, Int32 );
  extern void  maxsize		( Int32, Int32 );
  extern void  prefsize		( Int32, Int32 );
  extern void  prefposition 	( Int32, Int32, Int32, Int32 );
  extern void  stepunit		( Int32, Int32 );
  extern void  keepaspect	( Int32, Int32 );
  extern void  noport		( void );
  extern void  noborder		( void );
  
  extern void  ginit		( void );
  extern void  winconstraints 	( void );
  extern Int32 winopen		( Char8 * );
  extern Int32 swinopen		( Int32 );
  
  extern void  winposition 	( Int32, Int32, Int32, Int32 );
  extern void  winmove		( Int32, Int32 );
  
  extern void  getsize		( Int32 *, Int32 * );
  extern void  getorigin	( Int32 *, Int32 * );
  
  extern void  RGBmode		( void );
  extern void  cmode		( void );
  
  extern void  gconfig		( void );
  extern void  winclose		( Int32 );
  extern void  gexit		( void );
  
#ifdef _XLIB_H_ 		/* Declare if <X11/Xlib.h> is included */
  extern Int32 winX 		( Display *, Window );
#endif /* _XLIB_H_ */
  
  /* gl2ppm.c */
  extern int   gl2ppm		( const char * );
  
  /* 3d.c */
  extern void  cmov 		( Coord, Coord, Coord );
  extern void  cmovi 		( Icoord, Icoord, Icoord );
  extern void  cmovs 		( Scoord, Scoord, Scoord );
  
  extern void  pnt 		( Coord, Coord, Coord );
  extern void  pnti 		( Icoord, Icoord, Icoord );
  extern void  pnts 		( Scoord, Scoord, Scoord );
  
  extern void  move 		( Coord, Coord, Coord );
  extern void  movei 		( Icoord, Icoord, Icoord );
  extern void  moves 		( Scoord, Scoord, Scoord );
  
  extern void  rmv 		( Coord, Coord, Coord );
  extern void  rmvi 		( Icoord, Icoord, Icoord );
  extern void  rmvs 		( Scoord, Scoord, Scoord );
  
  extern void  draw 		( Coord, Coord, Coord );
  extern void  drawi 		( Icoord, Icoord, Icoord );
  extern void  draws 		( Scoord, Scoord, Scoord );
  
  extern void  rdr 		( Coord, Coord, Coord );
  extern void  rdri 		( Icoord, Icoord, Icoord );
  extern void  rdrs 		( Scoord, Scoord, Scoord );

  extern void  pmv 		( Coord, Coord, Coord );
  extern void  pmvi 		( Icoord, Icoord, Icoord );
  extern void  pmvs 		( Scoord, Scoord, Scoord );
  
  extern void  rpmv 		( Coord, Coord, Coord );
  extern void  rpmvi 		( Icoord, Icoord, Icoord );
  extern void  rpmvs 		( Scoord, Scoord, Scoord );
  
  extern void  pdr 		( Coord, Coord, Coord );
  extern void  pdri 		( Icoord, Icoord, Icoord );
  extern void  pdrs 		( Scoord, Scoord, Scoord );
  
  extern void  rpdr 		( Coord, Coord, Coord );
  extern void  rpdri 		( Icoord, Icoord, Icoord );
  extern void  rpdrs 		( Scoord, Scoord, Scoord );
  
  
  
  
  
  extern void  polf 		( Int32, Coord[][3] );
  extern void  polfi 		( Int32, Icoord[][3] );
  extern void  polfs 		( Int32, Scoord[][3] );

  extern void  v3s 		( Int16[3] );
  extern void  v3i 		( Int32[3] );
  extern void  v3f 		( Float32[3] );
  extern void  v3d 		( Float64[3] );
  
  extern void  v4s 		( Int16[4] );
  extern void  v4i 		( Int32[4] );
  extern void  v4f 		( Float32[4] );
  extern void  v4d 		( Float64[4] );

  extern void  swaptmesh 	( void );

  extern void  ortho 		( Coord, Coord, Coord, Coord, Coord, Coord );
  extern void  lookat 		( Coord, Coord, Coord, Coord, Coord, Coord, Angle );
  extern void  window 		( Coord, Coord, Coord, Coord, Coord, Coord );
  extern void  perspective 	( Angle, Float32, Coord, Coord );
  extern void  polarview 	( Coord, Angle, Angle, Angle );
  extern void  rot 		( Float32, char );
  extern void  rotate 		( Angle, char );
  extern void  scale 		( Float32, Float32, Float32 );
  extern void  translate 	( Coord, Coord, Coord );

  extern void  loadmatrix 	( Matrix );
  extern void  multmatrix 	( Matrix );
  extern void  pushmatrix 	( void );
  extern void  popmatrix 	( void );

  extern void  shademodel 	( Int32 );
  
  extern void  c4f 		( Float32[4] );
  extern void  c4i 		( Int32[4] );
  extern void  c4s 		( Int16[4] );

  extern void  n3f 		( Float32[3] );
  extern void  normal 		( Coord[3] );

  extern void  backface		( Int32 );
  extern void  frontface 	( Int32 );

  extern Int32 getmmode 	( void );
  extern void  mmode 		( Int16 );
  
  extern void  zbuffer 		( Int32 );
  extern void  zclear 		( void );
  extern void  zdraw 		( Int32 );
  extern void  zfunction 	( Int32 );
  extern void  czclear 		( Int32, Int32 );
  extern void  depthcue 	( Int32 );
  extern void  lRGBrange 	( Int16, Int16, Int16, Int16, Int16, Int16, Int32, Int32 );
  extern void  lsetdepth 	( Int32, Int32 );
  extern void  lshaderange 	( Colorindex, Colorindex, Int32, Int32 );

  /* Display lists */
  extern Int32 genobj 		( void );
  extern Int32 isobj 		( Int32 );
  extern void  makeobj 		( Int32 );
  extern Int32 getopenobj 	( void );
  extern void  closeobj		( void );
  extern void  callobj 		( Int32 );
  extern void  delobj 		( Int32 );
  
  /* Lighting */
  extern void  lmbind 		( Int32, Int32 );
  extern void  lmcolor 		( Int32 );
  extern void  lmdef 		( Int16, Int16, Int16, Float32[] );

  extern void  RGBwritemask 	( Int16, Int16, Int16 );
  
  /* for Pete Riley */
  extern void  drawmode 	( Int32 );
  extern void  iconsize 	( Int32, Int32 );
  extern void  overlay 		( Int32 );
  extern void  pushattributes 	( void );
  extern void  popattributes 	( void );
  extern void  fullscrn 	( void );
  extern void  endfullscrn 	( void );
  extern void  scrmask 		( Screencoord , Screencoord , Screencoord , Screencoord  );
  
  /* not implemented (yet) */
#if 0
  extern void attachcursor ( Device, Device );
  extern void bbox2 ( Screencoord, Screencoord, Coord, Coord, Coord, Coord );
  extern void bbox2i ( Screencoord, Screencoord, Icoord, Icoord, Icoord, Icoord );
  extern void bbox2s ( Screencoord, Screencoord, Scoord, Scoord, Scoord, Scoord );
  extern void bgnsurface ( void );
  extern void bgntrim ( void );
  extern void blankscreen ( Int32 );
  extern void blanktime ( Int32 );
  extern void blink ( Int16, Colorindex, Int16, Int16, Int16 );
  extern Int32 blkqread ( Int16 *, Int16 );
  extern void chunksize ( Int32 );
  extern void clkoff ( void );
  extern void clkon ( void );
  extern void colorf ( Float32 );
  extern void compactify ( Int32 );
  extern void crv ( Coord[4][3] );
  extern void crvn ( Int32, Coord[][3] );
  extern void curorigin ( Int16, Int16, Int16 );
  extern void cursoff ( void );
  extern void curson ( void );
  extern void curstype ( Int32 );
  extern void curvebasis ( Int16 );
  extern void curveit ( Int16 );
  extern void curveprecision ( Int16 );
  extern void cyclemap ( Int16, Int16, Int16 );
  extern void defbasis ( Int32, Matrix );
  extern void defcursor ( Int16, Uint16 * );
  extern void defpattern ( Int16, Int16, Int16 * );
  /*extern void defrasterfont ( Int16, Int16, Int16, Fontchar[], Int16, Int16[] );*/
  extern void deltag ( Int32 );
  extern void editobj ( Int32 );
  extern Int32 endselect ( Int16[] );
  extern void endsurface ( void );
  extern void endtrim ( void );
  extern void finish ( void );
  extern void fudge ( Int32, Int32 );
  extern void gRGBmask ( Int16 *, Int16 *, Int16 * );
  extern void gammaramp ( Int16[256], Int16[256], Int16[256] );
  extern void gbegin ( void );
  extern Int32 gentag ( void );
  extern Int32 getbackface ( void );
  extern Int32 getbuffer ( void );
  extern Int32 getcmmode ( void );
  extern void getcursor ( Int16 *, Colorindex *, Colorindex *, Int32 * );
  extern Int32 getdcm ( void );
  extern void getdev ( Int32, Device[], short[] );
  extern Int32 getdrawmode ( void );
  extern void getgpos ( Coord *, Coord *, Coord *, Coord * );
  extern Int32 getmap ( void );
  extern void getnurbsproperty (Int32, Float32 *);
  extern Int32 getpattern ( void );
  extern void getscrmask ( Screencoord *, Screencoord *, Screencoord *, Screencoord * );
  extern Int32 getshade ( void );
  extern Int32 getsm ( void );
  extern Int32 getwritemask ( void );
  extern Int32 getzbuffer ( void );
  extern void greset ( void );
  extern void gselect ( Int16 *, Int32 );
  extern void icontitle ( char * );
  extern void imakebackground ( void );
  extern Int32 istag ( Int32 );
  extern void lampoff ( char );
  extern void lampon ( char );
  extern void lgetdepth ( Int32 *, Int32 * );
  extern void linesmooth ( Int32 );
  extern void maketag ( Int32 );
  extern void mapcolors ( Colorindex, Colorindex, Int16 *, Int16 *, Int16 *);
  extern void mapw ( Int32, Screencoord, Screencoord, Coord *, Coord *, Coord *, Coord *, Coord *, Coord * );
  extern void mapw2 ( Int32, Screencoord, Screencoord, Coord *, Coord * );
  extern void multimap ( void );
  /*extern void newtag ( Int32, Int32, Offset );*/
  extern void nurbscurve (Int32, Float64 *,Int32, Float64 *, Int32, Int32);
  extern void nurbssurface (Int32, Float64 *,Int32, Float64 *, Int32, Int32, Float64 *, Int32, Int32, Int32);
  extern void objdelete ( Int32, Int32 );
  extern void objinsert ( Int32 );
  extern void objreplace ( Int32 );
  extern void onemap ( void );
  extern void patch ( Matrix, Matrix, Matrix );
  extern void patchbasis ( Int32, Int32 );
  extern void patchcurves ( Int32, Int32 );
  extern void patchprecision ( Int32, Int32 );
  extern void pixmode (Int32,Int32);
  extern void pntsmooth ( Int32 );

  extern void poly 		( Int32, Coord[][3] );
  extern void polyi 		( Int32, Icoord[][3] );
  extern void polys 		( Int32, Scoord[][3] );
  
  extern void polygonlist ( Int32, Int32, void * );
  extern void polylinelist ( Int32, Int32, void * );
  extern void pwlcurve (Int32, Float64 *,Int32, Int32);
  extern void rcrv ( Coord[4][4] );
  extern void rcrvn ( Int32, Coord[4][4] );
  extern void rpatch ( Matrix, Matrix, Matrix, Matrix );
  extern void screenspace ( void );
  extern void set_dither ( Int32 );
  extern void setcursor ( Int16, Colorindex, Colorindex );
  extern void setdblights ( Int32 );
  extern void setmap ( Int16 );
  extern void setnurbsproperty (Int32, Float32);
  extern void setpattern ( Int16 );
  extern void setvaluator ( Device, Int16, Int16, Int16 );
  extern void splf ( Int32, Coord[][3], Colorindex[] );
  extern void splf2 ( Int32, Coord[][2], Colorindex[] );
  extern void splf2i ( Int32, Icoord[][2], Colorindex[] );
  extern void splf2s ( Int32, Scoord[][2], Colorindex[] );
  extern void splfi ( Int32, Icoord[][3], Colorindex[] );
  extern void splfs ( Int32, Scoord[][3], Colorindex[] );
  extern void subpixel ( Int32 );
  extern void swapinterval ( Int16 );
  extern void textport ( Screencoord, Screencoord, Screencoord, Screencoord );
  extern void tpoff ( void );
  extern void tpon ( void );
  extern void underlay ( Int32 );
  extern void wmpack ( Uint32 );
  extern void writemask ( Colorindex );
  extern void zsource ( Int32 );
  extern void zwritemask ( Int32 );
#endif
#ifdef __cplusplus
}
#endif

#endif /* _YGL_INCLUDED_ */
