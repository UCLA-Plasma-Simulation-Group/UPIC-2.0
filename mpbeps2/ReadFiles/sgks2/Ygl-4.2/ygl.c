/*
 *    Ygl: Run GL programs with standard X11 and/or OpenGL routines.
 *    (C) Fred Hucht 1993-2007
 *    EMail: fred<at>thp.Uni-Duisburg.de
 */

#define VERSION "4.2e" /* gversion() only allows 12 chars (3 of VERSION) */

static const char vcid[] = "$Id: ygl.c,v 4.17 2007-05-11 20:27:39+02 fred Exp $";

static const char * const WhatString[]= {
  "@(#)Ygl V" VERSION " by Fred Hucht (C) 1993-2007",
  "@(#)http://www.thp.Uni-Duisburg.de/Ygl/",
  "@(#)EMail: fred <at> thp.Uni-Duisburg.de"
};

#include "header.h"
#include <X11/keysym.h>
#include <X11/cursorfont.h>

#ifdef HASXMU
#include <X11/Xatom.h>
#include <X11/Xmu/StdCmap.h>
#endif

static XSizeHints YglSizeHints;       /* for winconstraints() */
static int        YglMapNext = True,  /* flag to tell winconstraints() to map window */
		  YglRGB     = False; /* flag to tell gconfig() to switch to RGB mode */
static Cursor     YglCursor;          /* our special cursor */

static void       noinit(const char *caller);
static int        notimplemented(const char *caller);
static int        warn(const char *name, const char *format, ...);
static int        XRscale(double x);
static int        YRscale(double x);
static int        Xscale (double x);
static int        Yscale (double x);
static void       await_windowevent(Window win, long mask, int type);
static const char*devicename(Device dev);

static void       init_keymap(void);
static Int32      x2gl_wid(Window win, int which);
static int        is_wid(Int32 wid);
static int        ginit_badWindowHandler(Display *dpy, XErrorEvent *error);
static void       get_environment(const char *caller);
static void       create_gc(const char *caller, Int32);
static void       free_gc(Int32);
static void       clear_children(Int32);
static void       set_cmapwins(const char *caller, const YglWindow *cw, int install);
static void       create_main_window(const char *caller, const char *Title, int x, int y, Uint width, Uint height);
#ifdef RGBWIN
static void       create_RGB_window(const char *caller);
static void       destroy_RGB_window(const char *caller);
#endif
#if 0
static int        errorHandler(Display *, XErrorEvent *);
static int        unused_colormaps(void);
#endif
static void       signal_handler(int);
static Int32      win_open(const char *caller, const char *Title, Int32 pid, Window xpid);
static Ulong      parent_height(void);
static int        setup_visuals(const char *caller);
static void       setup_colormaps(const char *caller);
static void       setup_shift_and_mask(void);
static YglWindow* find_new_window_id(const char *caller);
static int        winX_badWindowHandler(Display *, XErrorEvent *);
static Window     get_real_parent(Window win);
static Window     get_parent(Window win);
static int        translate_xy(int xi, int yi, int *xf, int *yf);

#ifdef DEBUG
static const char *XEventNames[] = {
  "Unused0",
  "Unused1",
  "KeyPress",
  "KeyRelease",
  "ButtonPress",
  "ButtonRelease",
  "MotionNotify",
  "EnterNotify",
  "LeaveNotify",
  "FocusIn",
  "FocusOut",
  "KeymapNotify",
  "Expose",
  "GraphicsExpose",
  "NoExpose",
  "VisibilityNotify",
  "CreateNotify",
  "DestroyNotify",
  "UnmapNotify",
  "MapNotify",
  "MapRequest",
  "ReparentNotify",
  "ConfigureNotify",
  "ConfigureRequest",
  "GravityNotify",
  "ResizeRequest",
  "CirculateNotify",
  "CirculateRequest",
  "PropertyNotify",
  "SelectionClear",
  "SelectionRequest",
  "SelectionNotify",
  "ColormapNotify",
  "ClientMessage",
  "MappingNotify"
};
#endif

static const struct Device_ {
  const char *name;
  Device dev;
} device[] = {
  {"NULLDEV",		NULLDEV},
  {"BREAKKEY",		BREAKKEY},
  {"SETUPKEY",		SETUPKEY},
  {"LEFTCTRLKEY",	LEFTCTRLKEY},
  {"CAPSLOCKKEY",	CAPSLOCKKEY},
  {"RIGHTSHIFTKEY",	RIGHTSHIFTKEY},
  {"LEFTSHIFTKEY",	LEFTSHIFTKEY},
  {"ESCKEY",		ESCKEY},
  {"ONEKEY",		ONEKEY},
  {"TABKEY",		TABKEY},
  {"QKEY",		QKEY},
  {"AKEY",		AKEY},
  {"SKEY",		SKEY},
  {"NOSCRLKEY",		NOSCRLKEY},
  {"TWOKEY",		TWOKEY},
  {"THREEKEY",		THREEKEY},
  {"WKEY",		WKEY},
  {"EKEY",		EKEY},
  {"DKEY",		DKEY},
  {"FKEY",		FKEY},
  {"ZKEY",		ZKEY},
  {"XKEY",		XKEY},
  {"FOURKEY",		FOURKEY},
  {"FIVEKEY",		FIVEKEY},
  {"RKEY",		RKEY},
  {"TKEY",		TKEY},
  {"GKEY",		GKEY},
  {"HKEY",		HKEY},
  {"CKEY",		CKEY},
  {"VKEY",		VKEY},
  {"SIXKEY",		SIXKEY},
  {"SEVENKEY",		SEVENKEY},
  {"YKEY",		YKEY},
  {"UKEY",		UKEY},
  {"JKEY",		JKEY},
  {"KKEY",		KKEY},
  {"BKEY",		BKEY},
  {"NKEY",		NKEY},
  {"EIGHTKEY",		EIGHTKEY},
  {"NINEKEY",		NINEKEY},
  {"IKEY",		IKEY},
  {"OKEY",		OKEY},
  {"LKEY",		LKEY},
  {"SEMICOLONKEY",	SEMICOLONKEY},
  {"MKEY",		MKEY},
  {"COMMAKEY",		COMMAKEY},
  {"ZEROKEY",		ZEROKEY},
  {"MINUSKEY",		MINUSKEY},
  {"PKEY",		PKEY},
  {"LEFTBRACKETKEY",	LEFTBRACKETKEY},
  {"QUOTEKEY",		QUOTEKEY},
  {"RETKEY",		RETKEY},
  {"PERIODKEY",		PERIODKEY},
  {"VIRGULEKEY",	VIRGULEKEY},
  {"EQUALKEY",		EQUALKEY},
  {"ACCENTGRAVEKEY",	ACCENTGRAVEKEY},
  {"RIGHTBRACKETKEY",	RIGHTBRACKETKEY},
  {"BACKSLASHKEY",	BACKSLASHKEY},
  {"PAD1",		PAD1},
  {"PAD0",		PAD0},
  {"LINEFEEDKEY",	LINEFEEDKEY},
  {"BACKSPACEKEY",	BACKSPACEKEY},
  {"DELKEY",		DELKEY},
  {"PAD4",		PAD4},
  {"PAD2",		PAD2},
  {"PAD3",		PAD3},
  {"PADPERIOD",		PADPERIOD},
  {"PAD7",		PAD7},
  {"PAD8",		PAD8},
  {"PAD5",		PAD5},
  {"PAD6",		PAD6},
  {"PADPF2",		PADPF2},
  {"PADPF1",		PADPF1},
  {"LEFTARROWKEY",	LEFTARROWKEY},
  {"DOWNARROWKEY",	DOWNARROWKEY},
  {"PAD9",		PAD9},
  {"PADMINUS",		PADMINUS},
  {"PADCOMMA",		PADCOMMA},
  {"PADPF4",		PADPF4},
  {"PADPF3",		PADPF3},
  {"RIGHTARROWKEY",	RIGHTARROWKEY},
  {"UPARROWKEY",	UPARROWKEY},
  {"PADENTER",		PADENTER},
  {"SPACEKEY",		SPACEKEY},
  {"LEFTALTKEY",	LEFTALTKEY},
  {"RIGHTALTKEY",	RIGHTALTKEY},
  {"RIGHTCTRLKEY",	RIGHTCTRLKEY},
  {"F1KEY",		F1KEY},
  {"F2KEY",		F2KEY},
  {"F3KEY",		F3KEY},
  {"F4KEY",		F4KEY},
  {"F5KEY",		F5KEY},
  {"F6KEY",		F6KEY},
  {"F7KEY",		F7KEY},
  {"F8KEY",		F8KEY},
  {"F9KEY",		F9KEY},
  {"F10KEY",		F10KEY},
  {"F11KEY",		F11KEY},
  {"F12KEY",		F12KEY},
  {"PRINTSCREENKEY",	PRINTSCREENKEY},
  {"SCROLLLOCKKEY",	SCROLLLOCKKEY},
  {"PAUSEKEY",		PAUSEKEY},
  {"INSERTKEY",		INSERTKEY},
  {"HOMEKEY",		HOMEKEY},
  {"PAGEUPKEY",		PAGEUPKEY},
  {"ENDKEY",		ENDKEY},
  {"PAGEDOWNKEY",	PAGEDOWNKEY},
  {"NUMLOCKKEY",	NUMLOCKKEY},
  {"PADVIRGULEKEY",	PADVIRGULEKEY},
  {"PADASTERKEY",	PADASTERKEY},
  {"PADPLUSKEY",	PADPLUSKEY},
  {"MOUSEX",		MOUSEX},
  {"MOUSEY",		MOUSEY},
  {"ANYKEY",		ANYKEY},
  {"KEYBD",		KEYBD},
  {"TIMER0",		TIMER0},
  {"TIMER1",		TIMER1},
  {"TIMER2",		TIMER2},
  {"TIMER3",		TIMER3},
  {"REDRAW",		REDRAW},
  {"INPUTCHANGE",	INPUTCHANGE},
  {"WINCLOSE",		WINCLOSE},
  {"WINFREEZE",		WINFREEZE},
  {"WINTHAW",		WINTHAW},
  {"WINQUIT",		WINQUIT},
  {"DEPTHCHANGE",	DEPTHCHANGE},
  {"Unknown",		MAXYGLDEVICE}
};

YglControl Ygl = {
  NULL,			/* display */
  NULL, NULL,		/* Windows, active */
  0, 1,			/* ActiveWindow, NextWindow */
  /*---------------------- Bitfields: */
#ifdef AUTOFLUSH
  False,		/* flush */
#endif
  False,		/* wm_dw_flag */
  False,		/* stat */
  False,		/* PCM */	/* use private colormap? */
  False,		/* BS */	/* YGL_BS (backing store ) set */
  False,		/* GC */	/* Use single GC in CMap windows (slower) */
  False,		/* cmapmode uses TrueColor or DirectColor visual */
#ifdef OGL
  False,	       	/* Use OpenGL routines */
  False,		/* in GL_SELECT mode */
#endif
  /*---------------------- End of Bitfields */
  NULL,			/* Fonts */
  0,			/* LastFont */
  0,			/* Eventmask */
  0,			/* Atom wm_dw */

  100,			/* FT */	/* autoflush time in msec */
  {-2,-2},		/* V */ 	/* User requested visuals */
  0,			/* PWID */	/* Parent window ID of winopen wins */
  4,			/* DSZ */	/* Size of RGBmode dither matrix */
  
  0, 0,			/* CCmap, RCmap */
  
  NULL,			/* Colors */
  NULL,			/* ColorsInv */
  
  0, 0, 0,		/* rs, gs, bs */
  0, 0, 0,		/* rb, gb, bb */
#ifdef DEBUG
  XEventNames,		/* For debugging */
#else
  NULL,
#endif
  noinit,
  notimplemented,
  warn,
  x2gl_wid,
  XRscale,
  YRscale,
  Xscale,
  Yscale,
  await_windowevent,
  devicename
};

/* Static functions */
static const char* devicename(Device dev) {
  int i;
  for(i = 0; device[i].dev != dev && device[i].dev != MAXYGLDEVICE; i++);
  return device[i].name;
}

static void init_keymap(void) {
  static int inited = False;
  
  if(!inited) {
    int code;
    unsigned char wieoft[MAXYGLDEVICE];
    inited = True;
    
    memset(wieoft, 0, sizeof(wieoft));
    
    for(code = 0; code < KEYMAP_LEN; code++) {
      KeySym ks = XKeycodeToKeysym(D, code, 0);
      switch(ks) {
	/* main keys */
      case XK_a:		Ygl.keymap[code] = AKEY;		break;
      case XK_b:		Ygl.keymap[code] = BKEY;		break;
      case XK_c:		Ygl.keymap[code] = CKEY;		break;
      case XK_d:		Ygl.keymap[code] = DKEY;		break;
      case XK_e:		Ygl.keymap[code] = EKEY;		break;
      case XK_f:		Ygl.keymap[code] = FKEY;		break;
      case XK_g:		Ygl.keymap[code] = GKEY;		break;
      case XK_h:		Ygl.keymap[code] = HKEY;		break;
      case XK_i:		Ygl.keymap[code] = IKEY;		break;
      case XK_j:		Ygl.keymap[code] = JKEY;		break;
      case XK_k:		Ygl.keymap[code] = KKEY;		break;
      case XK_l:		Ygl.keymap[code] = LKEY;		break;
      case XK_m:		Ygl.keymap[code] = MKEY;		break;
      case XK_n:		Ygl.keymap[code] = NKEY;		break;
      case XK_o:		Ygl.keymap[code] = OKEY;		break;
      case XK_p:		Ygl.keymap[code] = PKEY;		break;
      case XK_q:		Ygl.keymap[code] = QKEY;		break;
      case XK_r:		Ygl.keymap[code] = RKEY;		break;
      case XK_s:		Ygl.keymap[code] = SKEY;		break;
      case XK_t:		Ygl.keymap[code] = TKEY;		break;
      case XK_u:		Ygl.keymap[code] = UKEY;		break;
      case XK_v:		Ygl.keymap[code] = VKEY;		break;
      case XK_w:		Ygl.keymap[code] = WKEY;		break;
      case XK_x:		Ygl.keymap[code] = XKEY;		break;
      case XK_y:		Ygl.keymap[code] = YKEY;		break;
      case XK_z:		Ygl.keymap[code] = ZKEY;		break;
      case XK_0:		Ygl.keymap[code] = ZEROKEY;		break;
      case XK_1:		Ygl.keymap[code] = ONEKEY;		break;
      case XK_2:		Ygl.keymap[code] = TWOKEY;		break;
      case XK_3:		Ygl.keymap[code] = THREEKEY;		break;
      case XK_4:		Ygl.keymap[code] = FOURKEY;		break;
      case XK_5:		Ygl.keymap[code] = FIVEKEY;		break;
      case XK_6:		Ygl.keymap[code] = SIXKEY;		break;
      case XK_7:		Ygl.keymap[code] = SEVENKEY;		break;
      case XK_8:		Ygl.keymap[code] = EIGHTKEY;		break;
      case XK_9:		Ygl.keymap[code] = NINEKEY;		break;
      case XK_apostrophe:	Ygl.keymap[code] = QUOTEKEY;		break;
      case XK_backslash:	Ygl.keymap[code] = BACKSLASHKEY;	break;
      case XK_bracketleft: 	Ygl.keymap[code] = LEFTBRACKETKEY;	break;
      case XK_bracketright:	Ygl.keymap[code] = RIGHTBRACKETKEY;	break;
      case XK_comma:		Ygl.keymap[code] = COMMAKEY;		break;
      case XK_equal:		Ygl.keymap[code] = EQUALKEY;		break;
      case XK_grave:		Ygl.keymap[code] = ACCENTGRAVEKEY;	break;
      case XK_minus:		Ygl.keymap[code] = MINUSKEY;		break;
      case XK_period:		Ygl.keymap[code] = PERIODKEY;		break;
      case XK_semicolon:	Ygl.keymap[code] = SEMICOLONKEY;	break;
      case XK_slash:		Ygl.keymap[code] = VIRGULEKEY;		break;
      case XK_space:		Ygl.keymap[code] = SPACEKEY;		break;
      
      case XK_BackSpace:	Ygl.keymap[code] = BACKSPACEKEY;	break;
      case XK_Escape:		Ygl.keymap[code] = ESCKEY;		break;
      case XK_Return:		Ygl.keymap[code] = RETKEY;		break;
      case XK_Tab:		Ygl.keymap[code] = TABKEY;		break;
	/* modifiers */
      case XK_Shift_L:		Ygl.keymap[code] = LEFTSHIFTKEY;	break;
      case XK_Shift_R:		Ygl.keymap[code] = RIGHTSHIFTKEY;	break;
      case XK_Alt_L:
      case XK_Meta_L:		Ygl.keymap[code] = LEFTALTKEY;		break;
      case XK_Alt_R:
      case XK_Meta_R:
      case XK_ISO_Level3_Shift:
      case XK_Mode_switch:	Ygl.keymap[code] = RIGHTALTKEY;		break;
      case XK_Control_L:	Ygl.keymap[code] = LEFTCTRLKEY;		break;
      case XK_Control_R:
      case XK_Multi_key:
      case XK_Execute:		Ygl.keymap[code] = RIGHTCTRLKEY;	break;
      case XK_Num_Lock:		Ygl.keymap[code] = NUMLOCKKEY;		break;
      case XK_Caps_Lock:
      case XK_Shift_Lock:	Ygl.keymap[code] = CAPSLOCKKEY;		break;
	/* cursor control */
      case XK_Delete:		Ygl.keymap[code] = DELKEY;		break;
      case XK_Down:		Ygl.keymap[code] = DOWNARROWKEY;	break;
      case XK_End:		Ygl.keymap[code] = ENDKEY;		break;
      case XK_Home:		Ygl.keymap[code] = HOMEKEY;		break;
      case XK_Insert:		Ygl.keymap[code] = INSERTKEY;		break;
      case XK_Left:		Ygl.keymap[code] = LEFTARROWKEY;	break;
      case XK_Next:		Ygl.keymap[code] = PAGEDOWNKEY;		break;
      case XK_Prior:		Ygl.keymap[code] = PAGEUPKEY;		break;
      case XK_Right:		Ygl.keymap[code] = RIGHTARROWKEY;	break;
      case XK_Up:		Ygl.keymap[code] = UPARROWKEY;		break;
	/* keypad */
      case XK_KP_0:
      case XK_KP_Insert:	Ygl.keymap[code] = PAD0;		break;
      case XK_KP_1:
      case XK_KP_End:		Ygl.keymap[code] = PAD1;		break;
      case XK_KP_2:
      case XK_KP_Down:		Ygl.keymap[code] = PAD2;		break;
      case XK_KP_3:
      case XK_KP_Next:		Ygl.keymap[code] = PAD3;		break;
      case XK_KP_4:
      case XK_KP_Left:		Ygl.keymap[code] = PAD4;		break;
      case XK_KP_5:
      case XK_KP_Begin:		Ygl.keymap[code] = PAD5;		break;
      case XK_KP_6:
      case XK_KP_Right:		Ygl.keymap[code] = PAD6;		break;
      case XK_KP_7:
      case XK_KP_Home:		Ygl.keymap[code] = PAD7;		break;
      case XK_KP_8:
      case XK_KP_Up:		Ygl.keymap[code] = PAD8;		break;
      case XK_KP_9:
      case XK_KP_Prior:		Ygl.keymap[code] = PAD9;		break;
      case XK_KP_F1:		Ygl.keymap[code] = PADPF1;		break;
      case XK_KP_F2:		Ygl.keymap[code] = PADPF2;		break;
      case XK_KP_F3:		Ygl.keymap[code] = PADPF3;		break;
      case XK_KP_F4:		Ygl.keymap[code] = PADPF4;		break;
      case XK_KP_Add:		Ygl.keymap[code] = PADPLUSKEY;		break;
      case XK_KP_Decimal:
      case XK_KP_Delete:	Ygl.keymap[code] = PADPERIOD;		break;
      case XK_KP_Enter:		Ygl.keymap[code] = PADENTER;		break;
      case XK_KP_Divide: 	Ygl.keymap[code] = PADVIRGULEKEY;	break;
      case XK_KP_Multiply: 	Ygl.keymap[code] = PADASTERKEY;		break;
      case XK_KP_Separator:	Ygl.keymap[code] = PADCOMMA;		break;
      case XK_KP_Subtract: 	Ygl.keymap[code] = PADMINUS;		break;
	/* function keys */
      case XK_F1:		Ygl.keymap[code] = F1KEY;		break;
      case XK_F2:		Ygl.keymap[code] = F2KEY;		break;
      case XK_F3:		Ygl.keymap[code] = F3KEY;		break;
      case XK_F4:		Ygl.keymap[code] = F4KEY;		break;
      case XK_F5:		Ygl.keymap[code] = F5KEY;		break;
      case XK_F6:		Ygl.keymap[code] = F6KEY;		break;
      case XK_F7:		Ygl.keymap[code] = F7KEY;		break;
      case XK_F8:		Ygl.keymap[code] = F8KEY;		break;
      case XK_F9:		Ygl.keymap[code] = F9KEY;		break;
      case XK_F10:		Ygl.keymap[code] = F10KEY;		break;
      case XK_F11:		Ygl.keymap[code] = F11KEY;		break;
      case XK_F12:		Ygl.keymap[code] = F12KEY;		break;
      case XK_F13:
      case XK_Print:		Ygl.keymap[code] = PRINTSCREENKEY;	break;
      case XK_F14:
      case XK_Scroll_Lock:
      case XK_Cancel:		Ygl.keymap[code] = SCROLLLOCKKEY;	break;
      case XK_F15:
      case XK_Pause:		Ygl.keymap[code] = PAUSEKEY;		break;
      
      default:			Ygl.keymap[code] = 0;			break;
      }
      wieoft[Ygl.keymap[code]]++;
#if DEBUG > 2
      if(ks != 0 || Ygl.keymap[code] != 0) 
	fprintf(stderr, "code = %d, XK = 0x%x, dev = %d (%s)\n",
		code, ks, Ygl.keymap[code], Ygl.devicename(Ygl.keymap[code]));
#endif
    }
#ifdef DEBUG
    for(code = 1; code < sizeof(wieoft)/sizeof(wieoft[0]); code++) if(wieoft[code] > 1) {
      fprintf(stderr, "device %d (%s) is %d times defined\n",
	      code, Ygl.devicename(code), wieoft[code]);
    }
#endif
    /* Try to fix PAD keys, as Linux, Solaris and maybe others define them as
     * simple numeric keysyms (XK_0 instead of XK_KP_0)
     * Here we assume that the keypad keycodes are > than the normal codes */
    {
      Device FixUs[][2] = {
	{ZEROKEY,   PAD0},
	{ONEKEY,    PAD1},
	{TWOKEY,    PAD2},
	{THREEKEY,  PAD3},
	{FOURKEY,   PAD4},
	{FIVEKEY,   PAD5},
	{SIXKEY,    PAD6},
	{SEVENKEY,  PAD7},
	{EIGHTKEY,  PAD8},
	{NINEKEY,   PAD9},
	{MINUSKEY,  PADMINUS},
	{RETKEY,    PADENTER},
	{PERIODKEY, PADCOMMA}
      };
      int i;
      for(i = 0; i < sizeof(FixUs)/sizeof(FixUs[0]); i++)
	if(wieoft[FixUs[i][0]] == 2) {
	  for(code = KEYMAP_LEN; code > 0 && Ygl.keymap[code] != FixUs[i][0]; code--);
	  if(code > 0) {
#ifdef DEBUG
	    fprintf(stderr, "fixing keymap[%d] from %d (%s) -> %d (%s)\n",
		    code,
		    Ygl.keymap[code], Ygl.devicename(Ygl.keymap[code]),
		    FixUs[i][1], Ygl.devicename(FixUs[i][1]));
#endif
	    Ygl.keymap[code] = FixUs[i][1];
	  }
	}
    }
  }
}

static void noinit(const char *caller) {
  Yprintf(caller, "Ygl not initialized, call winopen or ginit first.\n");
  exit(1);
}

static int notimplemented(const char *caller) {
  char * mode;
  IFOGL(mode = "OpenGL", mode = "X11");
  return Yprintf(caller, "not implemented in %s mode.\n", mode);
}

static int warn(const char *name, const char *format, ...) {
  int r;
  va_list Argp;
  va_start(Argp, format);
  r  = fprintf(stderr, "Ygl: %s: ", name);
  r += vfprintf(stderr, format, Argp);
  va_end(Argp);
  return r;
}

static Int32 x2gl_wid(Window win, int which) { /* Computes Ygl wid from X wid */
  Int32 i = Ygl.NextWindow - 1;
  switch(which) {
  case X2GL_MAIN: while(i > 0 && win != Ygl.Windows[i].main) i--; break;
#ifdef RGBWIN
  case X2GL_TOP:  while(i > 0 && win != Ygl.Windows[i].top ) i--; break;
  case X2GL_WIN:  while(i > 0 && win != Ygl.Windows[i].win ) i--; break;
  case X2GL_DRAW: while(i > 0 && win != Ygl.Windows[i].draw) i--; break;
#endif
  default:        abort();
  }
  return i; /* 0 if not found */
}

/* Routines to scale Ygl coords to X11 coords. */
#define IFLOOR(x)	(((int)(x)) - ((x) < 0 && (x) != ((double)(int)(x)) ? 1 : 0))
static int XRscale(double x) { x =  x          * W->xf + 0.5; return IFLOOR(x)            ; }
static int YRscale(double y) { y =  y          * W->yf + 0.5; return IFLOOR(y)            ; }
static int Xscale (double x) { x = (x - W->xo) * W->xf + 0.5; return IFLOOR(x)            ; }
static int Yscale (double y) { y = (W->yo - y) * W->yf + 0.5; return IFLOOR(y) + W->ym - 1; }

static void await_windowevent(Window win, long mask, int type) {
  XEvent ev;
#ifdef DEBUG
  fprintf(stderr,
	  "await_windowevent: window 0x%x, waiting for %d(%s) with mask 0x%x\n",
	  win, type, Ygl.XEventNames[type], mask);
#endif
  do {
    XWindowEvent(D, win, mask, &ev);
#ifdef DEBUG
    fprintf(stderr, "await_windowevent: received %d(%s)\n", ev.type, Ygl.XEventNames[ev.type]);
#endif
  } while(ev.type != type);
}

static int is_wid(Int32 wid) {
  return wid > 0 && wid < Ygl.NextWindow && Ygl.Windows[wid].main != 0;
}

#if 0
static int errorHandler(Display *dpy, XErrorEvent *error) {
  char errortext[1024];
  XGetErrorText(dpy, error->error_code, errortext, 1024);
  fprintf(stderr, "Ygl: X Error: %s on 0x%x.\n", errortext, error->resourceid);
  return False;
}
#endif

#if 0
static int unused_colormaps(void) { /* returns number of unused colormaps */
  int num;
  Colormap *cmaps = XListInstalledColormaps(D, Ygl.PWID, &num);
#ifdef DEBUG
  fprintf(stderr,"unused_colormaps: %d of %d installed Colormaps\n",
	 num, MaxCmapsOfScreen(DefaultScreenOfDisplay(D)));
#endif
  XFree((char*) cmaps);
  return MaxCmapsOfScreen(DefaultScreenOfDisplay(D)) - num;
}
#endif /* 0 */

#ifdef AUTOFLUSH
static void signal_handler(int sig_no) {
#ifdef DEBUG
  fprintf(stderr, "In signal_handler, sig_no = %d\n", sig_no);
#endif
  Ygl.flush = True; /* XFlush() at next possibility. Don't try to XFlush() 
		     * here, as it will totally screw up XLibs 
		     * internal buffers */
}
#endif /* AUTOFLUSH */

static int ginit_badWindowHandler(Display *dpy, XErrorEvent *error) {
  Yprintf("ginit", "invalid window id 0x%x in YGL_PARENTWID.\n",
	  error->resourceid);
  exit (1);
}

#define NUMVARS 9

static void get_environment(const char *caller) {  /* Parse environment */
  char *env;
  int v, i;
  int useogl = False;
  const char *names[NUMVARS][2] = { /* valid names of environment vars */
    { "YGL_FLUSHTIME", 	     "YGL_FT"  },
    { "YGL_PRIVATECOLORMAP", "YGL_PCM" },
    { "YGL_BACKINGSTORE",    "YGL_BS"  },
    { "YGL_SINGLEGC",        "YGL_GC"  },
    { "YGL_CMAPVISUAL",      "YGL_CV"  },
    { "YGL_RGBVISUAL",       "YGL_RV"  },
    { "YGL_PARENTWID",       "YGL_PWID"},
    { "YGL_DITHERSIZE",      "YGL_DSZ" },
    { "YGL_USEOGL",          "YGL_OGL" }
  };
  
  Ygl.PWID = RootWindow(D, YglScreen);
  
  for(v = 0; v < NUMVARS; v++) { /* loop over vars */
    for(i = 0; i < 2; i++) { /* loop over names */
      if((env = getenv(names[v][i])) != NULL && env[0] != '\0') {
	i = 4711; /* Found, end of i loop */
	switch(v) {
	case 0: Ygl.FT   = 1000 * atof(env); break;
	case 1: Ygl.PCM  = (strcmp(env,"1") == 0); break; 	/* Use private colormap */
	case 2: Ygl.BS   = (strcmp(env,"1") == 0); break; 	/* Use backing store */
	case 3: Ygl.GC   = (strcmp(env,"1") == 0); break; 	/* Use single graphics context */
	case 4: Ygl.V[0] = (int) strtol(env, NULL, 0); break;	/* Colormap Visual */
	case 5: Ygl.V[1] = (int) strtol(env, NULL, 0); break;	/* RGBmode  Visual */
	case 6: Ygl.PWID = (Window) strtol(env, NULL, 0); break;/* Parent window ID */
	case 7: Ygl.DSZ  = (int) strtol(env, NULL, 0); break; 	/* Size of dither matrix */
	case 8: useogl   = (int) strtol(env, NULL, 0); break;	/* Use OpenGL routines */
	}
      }
    }
  }
  
#ifdef DEBUG
  fprintf(stderr, 
	  "get_environment: FT = %d, PCM = %d, BS = %d, GC = %d,\n"
	  "                 V = {0x%x,0x%x}, PWID = 0x%x, DSZ = %d, UseOGL = %d\n",
	  Ygl.FT, Ygl.PCM, Ygl.BS, Ygl.GC,
	  Ygl.V[0], Ygl.V[1], Ygl.PWID, Ygl.DSZ, useogl);
#endif
  { /* Check if PWID is valid window */
    XWindowAttributes wa;
    int (*old_handler)(Display*, XErrorEvent*);
    old_handler = XSetErrorHandler(ginit_badWindowHandler);
    XGetWindowAttributes(D, Ygl.PWID, &wa);
    XSync(D, False);
    XSetErrorHandler(old_handler);
  }
  if(Ygl.DSZ < 0 || Ygl.DSZ > 6) {
    Yprintf(caller, "Invalid value of YGL_DSZ. Must be >= 0 and <= 6.\n");
    exit(1);
  }
  Ygl.V[0] = (Ygl.V[0] == -1) ? DefaultVisual(D, YglScreen)->visualid : -1;
  Ygl.V[1] = (Ygl.V[1] == -1) ? DefaultVisual(D, YglScreen)->visualid : -1;
#ifdef OGL
  Ygl.UseOGL = useogl;
#else
  if(useogl == 1) {
    Yprintf(caller, "Ygl not compiled for OpenGL support, ignoring YGL_OGL.\n");
  }
#endif
#ifndef X11
  if(Ygl.UseOGL == 0) {
    Yprintf(caller, "Ygl not compiled for X11 support, ignoring YGL_OGL.\n");
    Ygl.UseOGL = 1;
  }
#endif
}

static void create_gc(const char *caller, Int32 id) {
  int i;
  Ulong      mask = GCGraphicsExposures | GCLineWidth | GCForeground;
  XGCValues  values;
  YglWindow *w = &Ygl.Windows[id];
  
  values.graphics_exposures = False; 	/* so rectcopy won't generate 
					 * {Graphics|No}Expose events */
  values.line_width = 0; 		/* default line width */
  values.foreground = 0; 		/* fg color */
  
  w->linewidth  = 0;
  w->linestyle  = 0;
  w->color      = 0;
  
  if(w->rgb || Ygl.GC) { /* RGB windows have only one GC */
    w->gc     = XCreateGC(D, w->draw, mask, &values);
    w->chargc = w->gc;
  } else {
    w->chargc = XCreateGC(D, w->draw, mask, &values);
    if(NULL == (w->gclist = (GC*) calloc(CMapSize, sizeof(GC)))) {
      Yprintf(caller, "can't allocate memory.\n");
      exit(1);
    }
    for(i = 0; i < CMapSize; i++) {
      values.foreground = i;
      w->gclist[i] = XCreateGC(D, w->draw, mask, &values);
    }
    w->gc = w->gclist[0];
  }
  if(w->rgb && Ygl.DSZ > 1) { /* Initialize dithering */
    Ulong mask = 0;
    XGCValues values;
    w->pm   = XCreatePixmap(D, w->draw, Ygl.DSZ, Ygl.DSZ, YglDepth());
    w->pmgc = XCreateGC(D, w->draw, mask, &values);
    w->pmi  = XCreateImage(D, YglVisual(), YglDepth(), ZPixmap, 
			   0, NULL, Ygl.DSZ, Ygl.DSZ, 8, 0);
    w->pmi->data = (char*) malloc(w->pmi->bytes_per_line * Ygl.DSZ);
    if(w->pmi->data == NULL) {
      Yprintf(caller, "can't allocate memory.\n");
      exit(1);
    }
    w->red = w->green = w->blue = -1;
  } else {
    w->pm = 0; /* Used as flag for dithering */
  }
#ifdef DEBUG
  XSync(D, False);
  fprintf(stderr, "%s: create_gc: id=%d rgb=%d main=0x%x draw=0x%x gc=0x%x\n",
	  caller, id, w->rgb, w->main, w->draw, w->gc);
#endif
}

static void free_gc(Int32 id) {
  int i;
  YglWindow *w = &Ygl.Windows[id];
  if(w->rgb || Ygl.GC) {
    XFreeGC(D, w->gc);
    w->gc = NULL;
  } else {
    for(i = 0; i < CMapSize; i++) {
      XFreeGC(D, w->gclist[i]);
    }
    XFreeGC(D, w->chargc);
    free(w->gclist);
    w->gclist = NULL;
    w->gc     = NULL;
    w->chargc = NULL;
  }
  if(w->rgb && w->pm) {
    XFreeGC(D, w->pmgc);
    XFreePixmap(D, w->pm);
    free(w->pmi->data);
    XDestroyImage(w->pmi);
    w->pm = 0;
  }
}

static void clear_children(Int32 id) {
  Window root, parent, *children;
  Uint nchildren, n;
  Status r;
#ifdef DEBUG
  fprintf(stderr, "clear_children: id = %d\n", id);
#endif
  r = XQueryTree(D, Ygl.Windows[id].main, &root, &parent, &children, &nchildren);
  
  for(n = 0; n < nchildren; n++) {
    Int32 cid = x2gl_wid(children[n], X2GL_MAIN);
    if(cid) {
      clear_children(cid);
      Ygl.Windows[cid].main = 0;
    }
  }
  
  if(r) XFree((char*) children);
}

static Ulong parent_height(void) {
  XWindowAttributes pa;
  
  if(W->xpid == RootWindow(D, YglScreen)) {
    pa.height = DisplayHeight(D, YglScreen);
  } else {
    XGetWindowAttributes(D, W->xpid, &pa);
  }
  return pa.height;
}

static int translate_xy(int xi, int yi, int *xf, int *yf) {
  int r;
  Window junkwin;
  XTranslateCoordinates(D, W->xpid, W->main,
			xi, yi, xf, yf, &junkwin);
  r = (*xf != 0 || *yf != 0) ? 1 : 0;
#ifdef DEBUG
  fprintf(stderr, 
	  "translate_xy: active=%d, xi=%d, yi=%d, xf=%d, yf=%d, r=%d\n",
	  Ygl.ActiveWindow, xi, yi, *xf, *yf, r);
#endif
  *xf += xi;
  *yf += yi;
#ifdef DEBUG
  fprintf(stderr, "              xf=%d, yf=%d\n", *xf, *yf);
#endif
  return r;
}

static void set_cmapwins(const char *caller, const YglWindow *cw, int install) {
  Window *wins = (Window*)malloc((Ygl.NextWindow+2) * sizeof(Window));
  Int32 id, num = 0;
#ifdef RGBWIN
  if(cw->top == cw->main) {
    /* cw is top-level window, install RGB map first */
    wins[num++] = cw->win;
    wins[num++] = cw->top;
  } else {
    /* We are subwindow, install RGB map last */
    wins[num++] = cw->top;
    wins[num++] = cw->win;
  }
#ifdef DEBUG
  fprintf(stderr, "%s: CMAP_WINS(0x%x): 0x%x 0x%x",
	  caller, cw->top, wins[0], wins[1]);
#endif
#else /* RGBWIN */
  wins[num++] = cw->top;
#ifdef DEBUG
  fprintf(stderr, "%s: CMAP_WINS(0x%x): 0x%x",
	  caller, cw->top, wins[0]);
#endif
#endif /* RGBWIN */
  
  /* If any window is subwindow of me and has different color mode,
   * add to window list. One window is enough, as we only use two 
   * colormaps (cmap mode, RGB mode) */
  for(id = 1; id < Ygl.NextWindow; id++) {
    YglWindow *w = &Ygl.Windows[id];
    if(w->top == cw->top			/* is subwindow */
       && wins[0] != IF_RGBWIN(w->win, w->main) /* not already in list */
       ) {
      wins[num++] = IF_RGBWIN(w->win, w->main);
#ifdef DEBUG
      fprintf(stderr, " 0x%x", wins[num-1]);
#endif
    }
  }
  if((!install) || (num == 1 && cw->top == wins[0])) {
    /* Don't install redundant window info */
    /* Changed last arg from True to False due to problems with XServer
     * in RedHat 6.0. Thx to
     * Frederick W. Koehler IV <fwkoehl@sandia.gov> &
     * Kristina L. Klinkner <kklinkne@euclid.math.usfca.edu> */
    Atom WM_C_W = XInternAtom(D, "WM_COLORMAP_WINDOWS", False);
#ifdef DEBUG
    fprintf(stderr, " clearing\n");
#endif
    XDeleteProperty(D, cw->top, WM_C_W);
  } else {
#ifdef DEBUG
    fprintf(stderr, " installing\n");
#endif
    if(XSetWMColormapWindows(D, cw->top, wins, num) == 0) {
      Yprintf(caller, "cannot change WM_COLORMAP_WINDOWS property.\n");
    }
  }
  free(wins);
}

#ifdef OGL
static void create_main_window_ogl(const char *caller, YglWindow *w) {
  GLXContext oldcx = w->cx; /* If saved from gconfig() */
  if(NULL == (w->cx = glXCreateContext(D, w->rgb ? &Ygl.RV : &Ygl.CV,
				       NULL, True))) {
    Yprintf(caller, "cannot create context.\n");
    exit(1);
  }
  if(oldcx) {
    if (Ygl.CV.visualid != Ygl.RV.visualid) {
      /* copy state from window with other mode */
      glXCopyContext(D, oldcx, w->cx,
		     GL_COLOR_BUFFER_BIT |
		     GL_DEPTH_BUFFER_BIT |
		     GL_ENABLE_BIT |
		     /* GL_LIGHTING_BIT lighting fails when copied | */
		     GL_LINE_BIT |
		     GL_POINT_BIT |
		     GL_POLYGON_BIT |
		     GL_POLYGON_STIPPLE_BIT |
		     GL_TRANSFORM_BIT |
		     GL_VIEWPORT_BIT);
      glXDestroyContext(D, oldcx);
    } else {
      glXDestroyContext(D, w->cx);
      w->cx = oldcx;
    }
  }
  glXMakeCurrent(D, w->main, w->cx); /* Before anything else */
  /* Initialize */
  glDrawBuffer(w->dbuf ? GL_BACK : GL_FRONT);
  if(!oldcx) {
#if 0
    w->circ  = gluNewQuadric(); /* circ() */
    gluQuadricDrawStyle(w->circ,  GLU_SILHOUETTE);
    w->circf = gluNewQuadric(); /* circf() */
    gluQuadricDrawStyle(w->circf, GLU_FILL);
#endif
    ortho2  (0.0, w->xm - 1.0, 0.0, w->ym - 1.0);	/* set X11 default scale */
    viewport(0  , w->xm - 1  , 0  , w->ym - 1  );	/* set X11 default viewport */
  }
  if(glXIsDirect(D, w->cx))
    Yprintf(caller, "direct rendering enabled.\n");

  /*
   * Some optimization settings
   * Moved to here as glColorMask() fails in setup_visuals under MacOSX
   */
  /* glEnable(GL_NORMALIZE); GL sets this */
  /* glEnable(GL_ALPHA_TEST); */
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_FALSE); /* from OpenGL README */
  /* Yprintf(caller, "alpha writes disabled.\n"); */
}
#endif

static void create_main_window(const char *caller, const char *Title,
			       int x, int y, Uint width, Uint height) {
  XSetWindowAttributes swa;
  Ulong swa_mask = 0;
  YglWindow *w = W;
  int we_are_top;
  XWMHints   xwmh;
  XClassHint xch;
  
#ifdef DEBUG
  fprintf(stderr, "create_main_window: '%s': (x,y) = (%d,%d), (w,h) = (%d,%d)\n",
	  Title, x, y, width, height);
#endif
  
  we_are_top = w->pid == 0 && Ygl.PWID == RootWindow(D, YglScreen);
  
  swa.background_pixel = WhitePixel(D,YglScreen);
  swa.border_pixel     = WhitePixel(D,YglScreen);
  swa.colormap         = YglColormap();
  swa.cursor           = YglCursor;
  swa_mask = CWBackPixel | CWBorderPixel | CWColormap | CWCursor;
  
  if(!we_are_top) {
    /* Used with subwindows and with YGL_PWID,
     * as window positions count from leftlower corner in Ygl. 
     * Don't do this with normal WM-controlled windows. */
    swa.win_gravity    = SouthWestGravity;
    swa_mask          |= CWWinGravity;
  }
  
  if(Ygl.BS) {
    swa.backing_store  = Always;
    swa_mask          |= CWBackingStore;
  }
  
  w->main = 
    XCreateWindow(D, w->xpid,
		  x, y,
		  width, height,
		  0,		/* BorderWidth */
		  YglDepth(),	/* Depth */
		  InputOutput,	/* Class */
		  YglVisual(),
		  swa_mask,
		  &swa);
  
  w->draw = IF_RGBWIN(w->win = w->main, w->main);
  
  if(we_are_top) {         /* We are top-level window created by winopen() */
    w->top = w->main;
  } else if(w->pid == 0) { /* Created by winopen() and YGL_PWID was set */
    w->top = Ygl.PWID;
  } else {                 /* Copy from parent */
    w->top = Ygl.Windows[w->pid].top;
  }
  set_cmapwins(caller, w, 1);
  create_gc(caller, Ygl.ActiveWindow);
  
  XStoreName(D, w->main, Title);
  
  xwmh.flags = StateHint | InputHint;
  xwmh.input = True;
  xwmh.initial_state = NormalState;
  XSetWMHints(D, w->main, &xwmh);
  
  xch.res_name  = (char*)Title;
  xch.res_class = "Ygl";
  XSetClassHint(D, w->main, &xch);

  IFOGL(create_main_window_ogl(caller, w),
	ortho2  (0.0, w->xm - 1.0, 0.0, w->ym - 1.0);	/* set X11 default scale */
	viewport(0  , w->xm - 1  , 0  , w->ym - 1  )	/* set X11 default viewport */
	);
  
  winconstraints();
  
  font(0); /* Set default font to "fixed" */
  
#if 1
  if(MaxCmapsOfScreen(DefaultScreenOfDisplay(D)) > 1) {
    /*  if(1 || unused_colormaps() > 0) { */
    XInstallColormap(D, YglColormap());
#ifdef DEBUG
    fprintf(stderr,"create_main_window: Installing Colormap\n");
#endif
  }
#endif
}

static YglWindow *find_new_window_id(const char *caller) {
  Ygl.ActiveWindow = 1;
  while(is_wid(Ygl.ActiveWindow)) Ygl.ActiveWindow++;
  if(Ygl.ActiveWindow == Ygl.NextWindow) {
    Ygl.NextWindow++;
    if(NULL == (Ygl.Windows = (YglWindow*) realloc(Ygl.Windows, Ygl.NextWindow * sizeof(YglWindow)))) {
      Yprintf(caller, "can't allocate memory.\n", caller);
      exit(1);
    }
  }
  return W = &Ygl.Windows[Ygl.ActiveWindow];
}

static Int32 win_open(const char *caller, const char *Title, Int32 pid, Window xpid) {
  YglWindow *w;
  
  w = find_new_window_id(caller); /* Get new W */
  
  /* Setup window structure */
  w->pid  = pid;
  w->xpid = xpid;
  
  if(YglSizeHints.flags & PMinSize) {
    w->xm = YglSizeHints.min_width;
    w->ym = YglSizeHints.min_height;
  } else {
    w->xm = 100;
    w->ym = 100;
  }
  
  /* flip-fix y relative to parent window */
  YglSizeHints.y = parent_height() - w->ym - YglSizeHints.y; 
  
  w->xp = w->yp = w->xc = w->yc = 0;
  YglRGB = w->rgb = False; /* New windows are always in cmode */
  w->mapped  = False; /* window is not mapped yet */
  w->dbuf    = False; /* New windows are always in singlebuffer mode */
  w->clipped = False; /* No clipping activated */
  w->vmode   = VertexNone;
  w->pm      = 0;     /* Used in draw.c:pnt2() */
  w->vw      = w->xm; /* else set_scales() crashes the first time */
  w->vh      = w->ym;
#ifdef OGL
  w->zp = w->zc = 0;
  w->cx      = NULL;
#endif

  create_main_window(caller, Title, YglSizeHints.x, YglSizeHints.y, w->xm, w->ym);
  
  return Ygl.ActiveWindow;
}

static Window get_parent(Window win) {
  Window root, parent, *children;
  Uint nchildren;
  if(XQueryTree(D, win, &root, &parent, &children, &nchildren))
    XFree((char*) children);
  return parent;
}  

static Window get_real_parent(Window win) {
  Window root, parent, *children;
  Uint nchildren;
  Status r;
  int (*old_handler)(Display*, XErrorEvent*);
  Atom WM_STATE = XInternAtom(D, "WM_STATE", True);
  Atom type = None;
  
  old_handler = XSetErrorHandler(winX_badWindowHandler);
  r = XQueryTree(D, win, &root, &parent, &children, &nchildren);
  XSync (D, False);
  XSetErrorHandler(old_handler);
  if(r) XFree((char*) children);
  
  while(1) {
    int format;
    unsigned long nitems, after;
    unsigned char *data;
    
    if(parent == root) break; /* win is maindow, no WM */
    /* Find a window with WM_STATE as per ICCCM */
    if(XGetWindowProperty(D, parent, WM_STATE, 0, 0,
			  False, AnyPropertyType,
			  &type, &format, &nitems, &after, &data))
      XFree(data);
    if(type) break;
    win = parent;
    if(XQueryTree(D, win, &root, &parent, &children, &nchildren))
      XFree((char*) children);
#ifdef DEBUG
    fprintf(stderr, "get_real_parent: root=0x%x, win=0x%x, parent=0x%x\n",
	    /**/                      root     , win     , parent);
#endif
  }
  return parent;
}

/* Initialization (constraints) */

void minsize(Int32 x, Int32 y) {
  YglSizeHints.min_width  = x;
  YglSizeHints.min_height = y;
  YglSizeHints.flags |= PMinSize;
}

void maxsize(Int32 x, Int32 y) {
  YglSizeHints.max_width  = x;
  YglSizeHints.max_height = y;
  YglSizeHints.flags |= PMaxSize;
}

void prefsize(Int32 x, Int32 y) {
  YglSizeHints.width  = x;
  YglSizeHints.height = y;
  YglSizeHints.flags |= PSize;
  YglSizeHints.flags |= USSize;
  minsize(x, y);
  maxsize(x, y);
}

void prefposition(Int32 x1, Int32 x2, Int32 y1, Int32 y2) {
  YglSizeHints.x = x1;
  YglSizeHints.y = y1;
  YglSizeHints.flags |= PPosition;
  YglSizeHints.flags |= USPosition;
  prefsize(x2 - x1 + 1, y2 - y1 + 1);
}

void stepunit(Int32 x, Int32 y) {
  YglSizeHints.width_inc  = x;
  YglSizeHints.height_inc = y;
  YglSizeHints.flags |= PResizeInc;
}

void keepaspect(Int32 x, Int32 y) {
  YglSizeHints.min_aspect.x = YglSizeHints.max_aspect.x = x;
  YglSizeHints.min_aspect.y = YglSizeHints.max_aspect.y = y;
  YglSizeHints.flags |= PAspect;
}

void noport(void) {
  YglMapNext = False;
}

void noborder(void) {
  /* Do nothing... */
}

static void setup_shift_and_mask(void) {
  /* Determine shift value */
  Ygl.rs = 0; while (!((1 << Ygl.rs) & Ygl.RMask)) Ygl.rs++;
  Ygl.gs = 0; while (!((1 << Ygl.gs) & Ygl.GMask)) Ygl.gs++;
  Ygl.bs = 0; while (!((1 << Ygl.bs) & Ygl.BMask)) Ygl.bs++;
  /* Determine mask widths */
  Ygl.rb = 0; while ( (1 << (Ygl.rs + Ygl.rb)) & Ygl.RMask) Ygl.rb++;
  Ygl.gb = 0; while ( (1 << (Ygl.gs + Ygl.gb)) & Ygl.GMask) Ygl.gb++;
  Ygl.bb = 0; while ( (1 << (Ygl.bs + Ygl.bb)) & Ygl.BMask) Ygl.bb++;
#ifdef DEBUG
  fprintf(stderr,"setup_visuals: Shift = (%d,%d,%d), width = (%d,%d,%d)\n",
	  Ygl.rs, Ygl.gs, Ygl.bs, Ygl.rb, Ygl.gb, Ygl.bb);
#endif
  /* We need it in this form in RGBcolor ... */
  Ygl.rb = 8-Ygl.rb; Ygl.gb = 8-Ygl.gb; Ygl.bb = 8-Ygl.bb;
}

#define BETTER_VISUAL(oldi,newi) ( ((oldi) == -1 || xv[oldi].depth < xv[newi].depth) ? (newi) : (oldi) )

#ifdef OGL
static int setup_visuals_ogl(const char *caller) { /* find Visuals */
  int junk[2];
  int rgblist[]  = {GLX_RGBA,
		    GLX_DOUBLEBUFFER, 
		    GLX_RED_SIZE, 1, 
		    GLX_GREEN_SIZE, 1, 
		    GLX_BLUE_SIZE, 1,
		    GLX_DEPTH_SIZE, 16,
		    None};
  int cmaplist[] = {GLX_DOUBLEBUFFER,
		    GLX_RED_SIZE, 1, 
		    GLX_GREEN_SIZE, 1, 
		    GLX_BLUE_SIZE, 1,
		    None};
  XVisualInfo *xv;
  
  if(!glXQueryExtension(D, junk, junk + 1)) {
    Yprintf(caller, "X server has no OpenGL GLX extension.\n");
    return False;
  }
  
  if((xv = glXChooseVisual(D, YglScreen, rgblist))) {
    memcpy(&Ygl.RV, xv, sizeof(XVisualInfo));
    XFree((char*)xv);
    setup_shift_and_mask();
  } else {
    Yprintf(caller, "no OpenGL RGB visual found.\n");
    return False;
  }
  
  if((xv = glXChooseVisual(D, YglScreen, cmaplist))) {
    memcpy(&Ygl.CV, xv, sizeof(XVisualInfo));
    XFree((char*)xv);
  } else {
    Yprintf(caller,
	    "no OpenGL colormap visual found, emulating using visual 0x%x.\n",
	    Ygl.RV.visualid);
    memcpy(&Ygl.CV, &Ygl.RV, sizeof(XVisualInfo));
    Ygl.EmulateCmap = True;
    CMapSize = 1 << EMULATE_CMAP_DEPTH;
    /*return False;*/
  }
  
  /*Ygl.RCmap = XCreateColormap(D, Ygl.PWID, Ygl.RV.visual, AllocNone);*/
  /*Ygl.CCmap = XCreateColormap(D, Ygl.PWID, Ygl.CV.visual, AllocAll);*/
  
#ifdef DEBUG
  fprintf(stderr,
	  "setup_visuals: CMap VisualID = 0x%x, Depth = %d, CMEntries = %d, Ygl.stat = %d\n"
	  "               RGB  VisualID = 0x%x, Depth = %d, CMEntries = %d\n",
	  Ygl.CV.visualid, Ygl.CV.depth, Ygl.CV.colormap_size, Ygl.stat,
	  Ygl.RV.visualid, Ygl.RV.depth, Ygl.RV.colormap_size);
#endif
  
  return True;
}
#endif /* OGL */

#ifdef X11
static int setup_visuals_x11(const char *caller) { /* find best X11 Visuals */
  XVisualInfo xvproto, *xv;
  /* Note: Some implementations have VisualID == 0 */
  int Pseudoi = -1, Grayi = -1, Truei = -1, Directi = -1, SColori = -1, SGrayi = -1, Useri0 = -1, Useri1 = -1;
  int xvn, i = 0;
  
  xvproto.screen = YglScreen;
  xv = XGetVisualInfo(D, VisualScreenMask, &xvproto, &xvn);
  
  /* Find best visual of each class */
  for(i = 0; i < xvn; i++) switch(xv[i].class) {
  case StaticGray:  SGrayi  = BETTER_VISUAL(SGrayi , i); break;
  case GrayScale:   Grayi   = BETTER_VISUAL(Grayi  , i); break;
  case StaticColor: SColori = BETTER_VISUAL(SColori, i); break;
  case PseudoColor: Pseudoi = BETTER_VISUAL(Pseudoi, i); break;
  case TrueColor:   Truei   = BETTER_VISUAL(Truei  , i); break;
  case DirectColor: Directi = BETTER_VISUAL(Directi, i); break;
  default:
    Yprintf(caller, "unknown visual class: %d.\n", xv[i].class); break;
  }
  
  /* Find best colormap mode visual */
  if(Ygl.V[0] != -1) {
    for(i = 0; i < xvn; i++) if(Ygl.V[0] == xv[i].visualid) Useri0 = i;
    if(Useri0 == -1)
      Yprintf(caller, "VisualID 0x%x not found, using default.\n", Ygl.V[0]);
    else {
      switch(xv[Useri0].class) {
      case DirectColor: /* Consider as static */
      case TrueColor:	Ygl.GC = True; Ygl.EmulateCmap = True;
      case StaticGray: 
      case StaticColor: Ygl.stat = True; break;
      case GrayScale: 
      case PseudoColor: break;
      default:
	Yprintf(caller,	"VisualID 0x%x invalid for colormap mode, using default.\n", Ygl.V[0]);
	Useri0 = -1;
	break;
      }
    }
  }
  
  /* Find best RGBmode visual */
  if(Ygl.V[1] != -1) {
    for(i = 0; i < xvn; i++) if(Ygl.V[1] == xv[i].visualid) Useri1 = i;
    if(Useri1 == -1)
      Yprintf(caller, "VisualID 0x%x not found, using default.\n", Ygl.V[1]);
    else {
      switch(xv[Useri1].class) {
      case TrueColor:
      case DirectColor: break;
      default:
	Yprintf(caller, "VisualID 0x%x invalid for RGBmode, using default.\n", Ygl.V[1]);
	Useri1 = -1;
	break;
      }
    }
  }
  
  if      (Useri0  >= 0) i = Useri0;
  else if (Pseudoi >= 0) i = Pseudoi;
  else if (SColori >= 0){i = SColori; Ygl.stat = True; }
  else if (Grayi   >= 0) i = Grayi;
  else if (SGrayi  >= 0){i = SGrayi;  Ygl.stat = True; }
  else if (Truei   >= 0){i = Truei;   Ygl.GC = Ygl.stat = Ygl.EmulateCmap = True; }
  else if (Directi >= 0){i = Directi; Ygl.GC = Ygl.stat = Ygl.EmulateCmap = True; }
  else {
    Yprintf(caller, "no appropiate visual found.\n");
    return False;
  }
  
  if(i >= 0) {
    memcpy(&Ygl.CV, xv+i, sizeof(XVisualInfo));
    if(Ygl.EmulateCmap) {
      CMapSize  = 1 << EMULATE_CMAP_DEPTH;
    }
    if(Ygl.GC == False && CMapSize > MAX_MULTIGC_CMAPSIZE) {
      Ygl.GC = True;
#ifdef DEBUG
      fprintf(stderr, "setup_visuals: CMapSize > %d, using single GC.\n",
	      MAX_MULTIGC_CMAPSIZE);
#endif
    }
  }
  
  if      (Useri1  >= 0) i = Useri1;
  else if (Truei   >= 0) i = Truei;
  else if (Directi >= 0) i = Directi;
  else                   i = -1;
  
  if(i >= 0) {
    memcpy(&Ygl.RV, xv+i, sizeof(XVisualInfo));
    setup_shift_and_mask();
  }
  
  XFree((char*) xv);
  
#ifdef DEBUG
  fprintf(stderr,
	  "setup_visuals: CMap VisualID = 0x%x, Depth = %d, CMEntries = %d, Ygl.stat = %d\n"
	  "               RGB  VisualID = 0x%x, Depth = %d, CMEntries = %d\n",
	  Ygl.CV.visualid, Ygl.CV.depth, Ygl.CV.colormap_size, Ygl.stat,
	  Ygl.RV.visualid, Ygl.RV.depth, Ygl.RV.colormap_size);
#endif
  return True;
}
#endif /* X11 */

static int setup_visuals(const char *caller) { /* find Visuals */
#ifdef OGL
  if(Ygl.UseOGL) Ygl.UseOGL = setup_visuals_ogl(caller);
  if(Ygl.UseOGL) return True; /* OK, use OpenGL */
#endif
#ifdef X11
  return setup_visuals_x11(caller);
#else
  return False;
#endif
}

Int32 gversion(Char8 *v) {
  const char * MyName = "gversion";
  const char * version_x11 = "Ygl:X11-" VERSION;
  const char * version_ogl = "Ygl:OGL-" VERSION;
  const char * version;
  int r, useogl;
  
  if(D == NULL) {
    YglControl tmp;
    memcpy(&tmp, &Ygl, sizeof(Ygl)); /* Save a copy of fresh struct */
    D = XOpenDisplay(NULL);
    if(D == NULL) {
      r = False;
    } else {
      get_environment(MyName);
      r = setup_visuals(MyName);
      XCloseDisplay(D);
    }
#ifdef OGL
    useogl = Ygl.UseOGL;
#endif
    memcpy(&Ygl, &tmp, sizeof(Ygl)); /* Restore */
#ifdef OGL
    Ygl.UseOGL = useogl;
#endif
  } else {
    r = True;
  }
  
  IFOGL(version = version_ogl, version = version_x11);
  
  strncpy(v, version, 12);
  v[11] = '\0';
  
  return r ? 0 : -1;
}

#ifdef HASXMU
static void find_colormap_xmu(void) {
  XStandardColormap *stdcmap;
  int num, r, i, found = False;
  Atom map;
  
  switch(Ygl.CV.class) {
  case StaticGray:  map = XA_RGB_GRAY_MAP; break;
  case StaticColor: map = XA_RGB_BEST_MAP; break;
  case GrayScale:   map = XA_RGB_GRAY_MAP; break;
  case PseudoColor: map = XA_RGB_DEFAULT_MAP; break;
  }
  
  r = XmuLookupStandardColormap(D,
				YglScreen,
				Ygl.CV.visualid,
				Ygl.CV.depth,
				map,
				False, 
				True);
  
#ifdef DEBUG
  fprintf(stderr, "ginit: XmuLookupStandardColormap = %d\n", r);
#endif
  
  r = XGetRGBColormaps(D, Ygl.PWID, &stdcmap, &num, map);
  
  for(i = 0; i < num; i++) if(stdcmap[i].visualid == Ygl.CV.visualid) {
    Ygl.CCmap = stdcmap[i].colormap;
    found = True;
  }
  
  if(!found) {
    Yprintf(caller, "no appropiate standard colormap installed.\n");
    exit(1);
  }
}
#endif /* HASXMU */

static void setup_colormaps(const char *caller) {
  struct colors_ {
    Ushort red, green, blue;
    const char *name;
  } colors[] = {
    {    0,     0,     0, "black"  },
    {65535, 65535, 65535, "white"  },
    {    0, 65535,     0, "green"  },
    {65535, 65535,     0, "yellow" },
    {    0,     0, 65535, "blue"   },
    {65535,     0, 65535, "magenta"},
    {    0, 65535, 65535, "cyan"   },
    {65535,     0,     0, "red"    }
  };
  int i;
  
  if(!Ygl.PCM) { /* use global colormap */
#ifdef HASXMU
    find_colormap_xmu(void);
#else /* HASXMU */
    if(Ygl.CV.visual == DefaultVisual(D, YglScreen)) {
      Ygl.CCmap = DefaultColormap(D, YglScreen);
    } else {
      Ygl.CCmap = XCreateColormap(D, Ygl.PWID, Ygl.CV.visual, AllocNone);
    }
#endif /* HASXMU */
    
    if(NULL == (Ygl.Colors = (Ulong*) calloc(CMapSize, sizeof(Ulong)))) {
      Yprintf(caller, "can't allocate memory.\n");
      exit(1);
    }
    
    if(!Ygl.EmulateCmap) {
      if(NULL == (Ygl.ColorsInv = (Ulong*) calloc(CMapSize, sizeof(Ulong)))) {
	Yprintf(caller, "can't allocate memory.\n");
	exit(1);
      }
    }
    
    for(i = 0; i < MIN(8,CMapSize); i++) { /* define only allocated colors,
					    * others may change */
      XColor xc;
      xc.red   = colors[i].red;
      xc.green = colors[i].green;
      xc.blue  = colors[i].blue;
      xc.flags = DoRGB;
      
      if(XAllocColor(D, Ygl.CCmap, &xc)) {
	Ygl.Colors[i] = xc.pixel;
	if(!Ygl.EmulateCmap) Ygl.ColorsInv[xc.pixel] = i;
#ifdef DEBUG
	fprintf(stderr, "alloc color '%s', i = %d, pixel = 0x%x (%d,%d,%d)\n",
		colors[i].name, i, xc.pixel,  xc.red,  xc.green,  xc.blue);
#endif
      } else {
	Yprintf(caller, "can't allocate color %s = {%d,%d,%d}\n",
		colors[i].name, colors[i].red, colors[i].green, colors[i].blue);
      }
    }
  } else { /* PCM */
    XColor *xcs;
    if(Ygl.stat) {
      Yprintf(caller, "can't use private colormap on static visuals.\n");
      exit(1);
    }
    Ygl.CCmap = XCreateColormap(D, Ygl.PWID, Ygl.CV.visual, AllocAll);
    
    /* copy default colormap to reduce technicolor effects */
    
    if(NULL == (xcs = (XColor*) calloc(CMapSize, sizeof(XColor)))) {
      Yprintf(caller, "can't allocate memory.\n");
      exit(1);
    }
    
    for(i = 0; i < CMapSize; i++) xcs[i].pixel = i;
    XQueryColors(D, DefaultColormap(D, YglScreen), xcs, CMapSize);
    
    for(i = 0; i < MIN(8,CMapSize); i++) {
      xcs[CMapSize - 1 - i].red   = colors[i].red;
      xcs[CMapSize - 1 - i].green = colors[i].green;
      xcs[CMapSize - 1 - i].blue  = colors[i].blue;
    }
    XStoreColors(D, Ygl.CCmap, xcs, CMapSize);
    free(xcs);
  } /* PCM */

  if(Ygl.RV.visual) {
    Ygl.RCmap = XCreateColormap(D, Ygl.PWID, Ygl.RV.visual, AllocNone);
  }
}

#if DEBUG > 10
static int syncfn(Display *dpy) {
  fprintf(stderr, "Sync\n");
  return XSync(dpy, False);
}
#endif

void ginit(void) {
  const char * MyName = "ginit";
  if ((D = XOpenDisplay(NULL)) == NULL) {
    Yprintf(MyName, "can\'t open display \"%s\".\n", XDisplayName(NULL));
    exit(1);
  }

  /* allocate one window structure. The first window(0) is never used...*/
  if(NULL == (Ygl.Windows = (YglWindow*) calloc(1, sizeof(YglWindow)))) {
    Yprintf(MyName, "can't allocate memory.\n");
    exit(1);
  }
  
  Ygl.ActiveWindow = 0;
  Ygl.NextWindow   = 1;
  /* GL seems to set these... (see header.h) */
  Ygl.EventMask    = EnterLeaveMask | RedrawMask;

  get_environment(MyName);
  
#ifndef DEBUG
  /* XSetErrorHandler(errorHandler); */
#endif
  
  if(setup_visuals(MyName) == False) exit(1);
  
  setup_colormaps(MyName);
  
  if(Ygl.FT == -1000) {
    Yprintf(MyName, "switching to synchronous mode.\n");
    XSynchronize(D, True);
  }
  else if(Ygl.FT > 0) {
#ifdef AUTOFLUSH
    struct itimerval itimer, otimer;
    struct sigaction action, oaction;
    itimer.it_interval.tv_sec  =  Ygl.FT / 1000;
    itimer.it_value.tv_sec     =  Ygl.FT / 1000;
    itimer.it_interval.tv_usec = (Ygl.FT % 1000) * 1000;
    itimer.it_value.tv_usec    = (Ygl.FT % 1000) * 1000;
    Ygl.flush = False;
    action.sa_handler = signal_handler;
    sigfillset(&action.sa_mask);
    action.sa_flags   = 0;
    if(sigaction(SIGVTALRM, &action, &oaction)) {
      perror(MyName);
      exit(1);
    }
    if(oaction.sa_handler) Yprintf(MyName, "warning: SIGVTALRM already used.\n");
    setitimer(ITIMER_VIRTUAL, &itimer, &otimer);
#endif /* AUTOFLUSH */
  }

#if DEBUG > 10
  fprintf(stderr, "DEBUG: setting XSetAfterFunction(D, syncfn)\n");
  XSetAfterFunction(D, syncfn);
#elif DEBUG > 1
  fprintf(stderr, "DEBUG: setting XSynchronize(D, True)\n");
  XSynchronize(D, True);
#endif
  
  loadXfont(0, "fixed");
  
  Ygl.wm_dw = XInternAtom(D, "WM_DELETE_WINDOW", False);

  YglCursor = XCreateFontCursor(D, XC_left_ptr);
  {
    XColor fg, bg;
    bg.red   = 0     ; bg.green = 0     ; bg.blue  = 0xffff;
    fg.red   = 0xffff; fg.green = 0xffff; fg.blue  = 0     ;
    
    XRecolorCursor(D, YglCursor, &fg, &bg);
  }
  
  if(Ygl.RV.depth >= 24) {
#ifdef DEBUG
    fprintf(stderr, "Won't dither, depth >= 24\n");
#endif
    Ygl.DSZ = 0; /* Don't need to dither if we have 24 bits */
  }
  
  init_keymap();
}

void winconstraints(void) {
  const char * MyName = "winconstraints";
  YglWindow  *w = W;
  
  I(MyName, "");
  XSetWMNormalHints(D, w->main, &YglSizeHints);
  
  if(YglMapNext && !w->mapped) { /* Map window? */
#ifdef DEBUG
    fprintf(stderr, "winconstraints: Mapping window %d\n", Ygl.ActiveWindow);
#endif
    XSelectInput(D, w->main, StructureNotifyMask | ExposureMask);
    /* So we get the MapNotify and the first Exposure(REDRAW) event */
    w->mapped = True;
    XMapWindow(D, w->main);
#ifdef RGBWIN
    if(w->rgb) XMapWindow(D, w->win);
#endif
    Ygl.await_windowevent(w->main, StructureNotifyMask, MapNotify);
    
    if(YglSizeHints.flags & PPosition && get_parent(w->main) != w->xpid) {
      /* Fix window position gobbled by window manager */
      int xf, yf;
#ifdef DEBUG
      fprintf(stderr, "winconstraints: Moving to (%d,%d)\n", YglSizeHints.x, YglSizeHints.y);
#endif
      XMoveWindow(D, w->main, YglSizeHints.x, YglSizeHints.y);
      Ygl.await_windowevent(w->main, StructureNotifyMask, ConfigureNotify);
      if(translate_xy(YglSizeHints.x, YglSizeHints.y, &xf, &yf)) {
#ifdef DEBUG
	fprintf(stderr, "winconstraints: Moving to (%d,%d)\n", xf, yf);
#endif
	XMoveWindow(D, w->main, xf, yf);
	Ygl.await_windowevent(w->main, StructureNotifyMask, ConfigureNotify);
      }
    }
    
    XSetWMProtocols(D, w->main, &Ygl.wm_dw, Ygl.wm_dw_flag ? 1 : 0);
    XSelectInput(D, w->main, Ygl.EventMask);
#ifdef RGBWIN
    if(w->rgb) /* Expose events don't propagate */
      XSelectInput(D, w->win, Ygl.EventMask & ExposureMask);
#endif
    qenter(REDRAW, Ygl.ActiveWindow);
  } else if(w->mapped && !YglMapNext) { /* unmap window */
#ifdef DEBUG
    fprintf(stderr, "winconstraints: Unmapping window %d\n", Ygl.ActiveWindow);
#endif
    XUnmapWindow(D, w->main);
    w->mapped = False;
  }
  
  /* Reinit constraints */
  memset(&YglSizeHints, 0, sizeof(XSizeHints)); /* clear YglSizeHints */
  YglMapNext = True;
  /* YglSizeHints.flags |= PSize;
     YglSizeHints.width  = 100;
     YglSizeHints.height = 100; */
}

Int32 winopen(char *Title) {
  const char * MyName = "winopen";
  if(D == NULL) ginit(); /* first window */
#if 0  
  if(0 == (YglSizeHints.flags & PMinSize)) minsize(40,30);
  if(0 == (YglSizeHints.flags & PMaxSize)) maxsize(YglScreenWidth, YglScreenHeight);
#endif
  return win_open(MyName, Title, 0, Ygl.PWID);
}

Int32 swinopen(Int32 pid) {
  const char * MyName = "swinopen";
  I(MyName, "%d", pid);
  if(is_wid(pid)) {
    return win_open(MyName, NULL, pid, Ygl.Windows[pid].main);
  } else {
    Yprintf(MyName, "parent does not exist: %d.\n", pid);
    return 0;
  }
}

void winposition(Int32 x1, Int32 x2, Int32 y1, Int32 y2) {
  const char * MyName = "winposition";
  YglWindow *w = W;
  int x, y;
  I(MyName, "%d,%d,%d,%d", x1, x2, y1, y2 );
  x = x1;
  y = parent_height() - 1 - y2;
  
  if(get_parent(w->main) != w->xpid) {
    /* Fix window position gobbled by window manager */
    int xf, yf;
    XSelectInput(D, w->main, StructureNotifyMask);
    XMoveResizeWindow(D, w->main, x,  y,  x2 - x1 + 1, y2 - y1 + 1);
#ifdef DEBUG
    XSync(D, False);
#else
    Ygl.await_windowevent(w->main, StructureNotifyMask, ConfigureNotify);
#endif
    XSelectInput(D, w->main, Ygl.EventMask);
    if(translate_xy(x, y, &xf, &yf))
      XMoveWindow(D, w->main, xf, yf);
  } else {
    XMoveResizeWindow(D, w->main, x,  y,  x2 - x1 + 1, y2 - y1 + 1);
  }
}

void winmove(Int32 x, Int32 y) {
  const char * MyName = "winmove";
  Int32 xs, ys;
  I(MyName, "%d,%d", x, y);
  getsize(&xs, &ys);
  XMoveWindow(D, W->main, x, parent_height() - y - ys);
}

void getsize(Int32 *x, Int32 *y) { 
  const char * MyName = "getsize";
  XWindowAttributes wa;
  I(MyName, "*,*");
  XGetWindowAttributes(D, W->main, &wa);
  *x = wa.width;
  *y = wa.height;
}

void getorigin(Int32 *x, Int32 *y) {
  /* Return the position of 
   * the lower left corner of the drawing area relative to 
   * the lower left corner of the parent window
   * WRONG! relative to Ygl.PWID (normally root) */
  const char * MyName = "getorigin";
  XWindowAttributes wa, pa;
  Window junkwin;
  int rx, ry;
  I(MyName, "*,*");
  XGetWindowAttributes(D, W->main, &wa);
  XTranslateCoordinates(D, W->main, Ygl.PWID,
			-wa.border_width,
			-wa.border_width,
			&rx, &ry, &junkwin);
  
  if(Ygl.PWID == RootWindow(D, YglScreen)) {
    pa.height = DisplayHeight(D, YglScreen);
  } else {
    XGetWindowAttributes(D, Ygl.PWID, &pa);
  }
  
  *x = rx;
  *y = pa.height - (wa.height + ry);
#ifdef DEBUG
  fprintf(stderr, "getorigin: %d %d, returning %d %d\n",
	  rx, ry, *x, *y);
#endif
}

void RGBmode(void) {
  /* switch the window active on next gconfig() to RGB mode */
  YglRGB = True;
}

void cmode(void) {
  /* switch the window active on next gconfig() to colormap mode */
  YglRGB = False;
}

#ifdef RGBWIN
static void create_RGB_window(const char *caller) {
  XSetWindowAttributes swa;
  Ulong swa_mask = 0;
  YglWindow *w = W;
  
#ifdef DEBUG
  fprintf(stderr, "create_RGB_window: (w,h) = (%d,%d)\n", w->xm, w->ym);
#endif
  
  /*swa.background_pixel = WhitePixel(D,YglScreen);*/
  swa.colormap         = Ygl.RCmap;
  swa.cursor           = YglCursor;
  swa_mask = /*CWBackPixel|*/CWColormap | CWCursor;
  
  if(Ygl.BS) {
    swa.backing_store  = Always;
    swa_mask          |= CWBackingStore;
  }
  
  w->draw = w->win = 
    XCreateWindow(D, w->main,
		  0, 0,
		  w->xm, w->ym,
		  0,			/* BorderWidth */
		  Ygl.RV.depth,		/* Depth */
		  InputOutput,		/* Class */
		  Ygl.RV.visual,
		  swa_mask,
		  &swa);

  XStoreName(D, w->win, "RGB");
  XSelectInput(D, w->win, StructureNotifyMask);
  XMapWindow(D, w->win);
  Ygl.await_windowevent(w->win, StructureNotifyMask, MapNotify);
  XSelectInput(D, w->win, Ygl.EventMask & ExposureMask);
  XLowerWindow(D, w->win); /* So it does not obscure any subwindows */
}

static void destroy_RGB_window(const char *caller) {
  YglWindow *w = W;
  XEvent ev;
  
  XDestroyWindow(D, w->win);
  XSync(D, False);
  while(XCheckWindowEvent(D, w->win, Ygl.EventMask, &ev)) {
#ifdef DEBUG
    fprintf(stderr, "destroy_RGB_window: dropped %d(%s) for window 0x%x\n",
	    ev.type, Ygl.XEventNames[ev.type], w->win);
#endif
  }
  
  w->draw = w->win = w->main;
}
#endif /* RGBWIN */

#ifdef OGL
static int destroycontext = True;
#endif

void gconfig(void) {
  const char * MyName = "gconfig";
  YglWindow *w = W;
  short dbuf;
#ifndef RGBWIN
  XWindowAttributes wa;
  Window junkwin;
  int x, y;
  char *Title;
#endif
  
  I(MyName, "");
  if(YglRGB == w->rgb) return; /* no mode change */
  if(Ygl.RV.visual == NULL) {
    Yprintf(MyName, "RGBmode requires a TrueColor or DirectColor visual.\n");
    exit(1);
  }
  
  dbuf = w->dbuf;
#ifdef RGBWIN
  singlebuffer();
  free_gc(Ygl.ActiveWindow);
#else /* RGBWIN */
  /* Save state of active window */
  XGetWindowAttributes(D, w->main, &wa);
  XTranslateCoordinates(D, w->main, w->xpid,
			-wa.border_width/*-wa.x*/,
			-wa.border_width/*-wa.y*/,
			&x, &y, &junkwin);
  XFetchName(D, w->main, &Title);
  XGetNormalHints(D, w->main, &YglSizeHints);  /* ??? */
#ifdef DEBUG
  fprintf(stderr,"gconfig: Title = '%s', x = %d, y = %d, xt = %d, yt = %d\n", 
	  Title, wa.x, wa.y, x, y);
#endif
  /* And close it */
#ifdef OGL
  destroycontext = False;
#endif
  winclose(Ygl.ActiveWindow);
#endif /* RGBWIN */
  
  w->rgb = YglRGB; /* Switch to new mode */
  
#ifdef RGBWIN
  if(YglRGB) { /* switch to RGB */
    create_RGB_window(MyName);
  } else { /* switch to cmap */
    destroy_RGB_window(MyName);
  }
  create_gc(MyName, Ygl.ActiveWindow);
  font(0); /* Set default font to "fixed" */
#else /* RGBWIN */
  /* And open a new one with other visual at same place */
  YglSizeHints.x = x;
  YglSizeHints.y = y;
  YglSizeHints.flags |= PPosition | USPosition;
  YglMapNext = w->mapped; /* map if window was mapped */
  w->mapped  = False;     /* window is not mapped yet */
  create_main_window(MyName, Title, x, y, wa.width, wa.height);
  XFree(Title);
#endif /* RGBWIN */
  /*set_cmapwins(MyName, w, 1);*/
  
  if(dbuf) doublebuffer();
}

void winclose(Int32 id) {
  const char * MyName = "winclose";
  XEvent ev;
  
  I(MyName, "%d", id);
  
  if(is_wid(id)) {
    Int32 active = winget();
    YglWindow *w = &Ygl.Windows[id];
    set_cmapwins(MyName, w, 0); /* Remove from list */
    free_gc(id);
    clear_children(id); /* marks all children as closed */

    winset(id);
    singlebuffer();
    if (Ygl.Windows[active].main != 0)
      winset(active);
    else
      winset(1);

#ifdef RGBWIN
    if(w->rgb) destroy_RGB_window(MyName);
#endif
#ifdef OGL
    if(Ygl.UseOGL) {
      glXMakeCurrent(D, None, NULL);
#if 0
      gluDeleteQuadric(w->circ);
      gluDeleteQuadric(w->circf);
#endif
      if(destroycontext) glXDestroyContext(D, w->cx);
      destroycontext = True; /* Reinit */
    }
#endif
    XDestroyWindow(D, w->main);
    XSync(D, False);
    while(XCheckWindowEvent(D, w->main, Ygl.EventMask, &ev)) {
#ifdef DEBUG
      fprintf(stderr, "winclose: dropped %d(%s) for window 0x%x\n",
	      ev.type, Ygl.XEventNames[ev.type], w->main);
#endif
    }
    w->main = 0;
  } else {
    Yprintf(MyName, "invalid window id: %d\n", id);
  }
}

void gexit(void) {
  const char * MyName = "gexit";
  int i;
  I(MyName, "");
  for(i = 1; i < Ygl.NextWindow; i++) if(is_wid(i)) winclose(i);

#ifdef HASXMU
  if(Ygl.PCM) {
    XFreeColormap(D, Ygl.CCmap);
  }
#else
  if(Ygl.PCM || Ygl.CCmap != DefaultColormap(D, YglScreen)) {
    XFreeColormap(D, Ygl.CCmap);
  }
#endif
  
  if(!Ygl.PCM) {
    free(Ygl.Colors);
    if(!Ygl.EmulateCmap)
      free(Ygl.ColorsInv);
  }

  if(Ygl.RCmap != 0) XFreeColormap(D, Ygl.RCmap);

  free(Ygl.Windows);

  XFreeCursor(D, YglCursor);
  
  XCloseDisplay(D);
  D = NULL;
}

static int winX_badWindowHandler(Display *dpy, XErrorEvent *error) {
  Yprintf("winX", "no window with id 0x%x found.\n", error->resourceid);
  exit (1);
}

Int32 winX(Display *dpy, Window win) {
  const char * MyName = "winX";
  XWindowAttributes wa;
  XSetWindowAttributes swa;
  Ulong swa_mask = 0;
  YglWindow *w;
  
  if(D == NULL) ginit(); /* first window. Here we assume that dpy is taken from envvar DISPLAY... */
  
  if(win == RootWindow(D, YglScreen)) {
    Yprintf(MyName, "Cannot convert root window.\n");
    exit (1);
  }
  
  w = find_new_window_id(MyName);
  
  w->xpid = get_real_parent(win); /* Also checks if win is valid */
  
  XGetWindowAttributes(D, win, &wa);
  
  /* Setup window structure */
  w->top  = win;
  w->main = win;
#ifdef RGBWIN
  w->win  = win;
#endif
  w->draw = win;
  w->xm   = wa.width;
  w->ym   = wa.height;
  w->vw   = w->xm; /* else set_scales() crashes the firat time */
  w->vh   = w->ym;
  w->xp   = w->yp = w->xc = w->yc = 0;
  
  if(wa.visual->visualid == Ygl.CV.visualid) {
    w->rgb = False;
  } else if(wa.visual->visualid == Ygl.RV.visualid) {
    w->rgb = True;
  } else {
    Yprintf(MyName,
	    "unsupported visual type, must be 0x%x (cmap) or 0x%x (RGB).\n",
	    Ygl.CV.visualid, Ygl.RV.visualid);
    exit(1);
  }
  
  w->mapped  = wa.map_state != IsUnmapped;
  w->pid     = -1;
  w->dbuf    = False; /* New windows are always in singlebuffer mode */
  w->clipped = False; /* No clipping activated */
  w->vmode   = VertexNone;
  
  /* Set up window attributes */
  swa.background_pixel = WhitePixel(D,YglScreen);
  swa.border_pixel     = WhitePixel(D,YglScreen);
  swa.colormap         = YglColormap();
  swa.cursor           = YglCursor;
  swa_mask = CWBackPixel | CWBorderPixel | CWColormap | CWCursor;
  
  if(Ygl.BS) {
    swa.backing_store = Always;
    swa_mask |= CWBackingStore;
  }
  
  XChangeWindowAttributes(D, w->main, swa_mask, &swa);
  
  create_gc(MyName, Ygl.ActiveWindow);
  font(0); /* Set default font to "fixed" */
  
  ortho2  (0.0, w->xm - 1.0, 0.0, w->ym - 1.0);	/* set default scale */
  viewport(0  , w->xm - 1  , 0  , w->ym - 1  );	/* set default viewport */
  
#ifdef DEBUG
  fprintf(stderr, "winX: id = %d, done.\n", Ygl.ActiveWindow);
#endif
  
  qenter(REDRAW, Ygl.ActiveWindow);
  
  return Ygl.ActiveWindow;
}
