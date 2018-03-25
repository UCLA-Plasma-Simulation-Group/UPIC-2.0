/*
 *    Ygl: Run GL programs with standard X11 and/or OpenGL routines.
 *    (C) Fred Hucht 1993-2007
 *    EMail: fred<at>thp.Uni-Duisburg.de
 */

static const char vcid[] = "$Id: misc.c,v 4.14 2007-05-08 13:25:45+02 fred Exp $";

#ifdef OGL
# define GLX_GLXEXT_PROTOTYPES /* we want prototypes for glX*VideoSyncSGI() */
#endif

#include "header.h"

#ifdef _AIX
/* # include <sys/select.h> fails under AIX 3.2 */
#endif

#ifdef X11
# ifdef DOUBLEBUF
#  include <X11/extensions/Xdbe.h>
# endif

# ifdef MULTIBUF
#  ifndef NO_MULTIBUF_H /* we have both multibuf.h and libXext.a */
#   include <X11/extensions/multibuf.h>
#  else /* we don't have the include file, but the library (HP/UX 8.0) */
#   define MultibufferUpdateActionUndefined	0
#   define MultibufferUpdateHintFrequent		0
#  endif
# endif
#endif /* X11 */

#define YglClipMask 	1L
#define YglDashMask 	2L
#define YglGCMask	4L

typedef struct YglGCValues_ {
  Ulong		gcmask;		/* The GC mask */
  XGCValues 	values;		/* The GC values */
  XRectangle	rect;		/* The clipping rectangle */
  char		*dashes;	/* The dash list */
  int		ndashes;	/* Number of elements in dashes */
} YglGCValues;

static void change_gc(Ulong, YglGCValues *);
#ifdef X11
static void set_scales(YglWindow *w);
#endif
static int  linestyle_id2idx(Int32 id);
#ifdef MULTIBUF
static int  mbuf_errorhandler(Display *dpy, XErrorEvent *error);
#endif

static void change_gc(Ulong mask, YglGCValues *yv) {
  YglWindow *w = W;
  if(w->rgb || Ygl.GC) { /* RGB windows have only one GC */
    if(mask & YglClipMask) XSetClipRectangles(D, w->gc, 0, 0, &yv->rect, 1, Unsorted);
    if(mask & YglDashMask) XSetDashes        (D, w->gc, 0, yv->dashes, yv->ndashes);
    if(mask & YglGCMask)   XChangeGC         (D, w->gc, yv->gcmask, &yv->values);
  } else {
    int i;
    for(i = 0; i < CMapSize; i++) {
      if(mask & YglClipMask) XSetClipRectangles(D, w->gclist[i], 0, 0, &yv->rect, 1, Unsorted);
      if(mask & YglDashMask) XSetDashes        (D, w->gclist[i], 0, yv->dashes, yv->ndashes);
      if(mask & YglGCMask)   XChangeGC         (D, w->gclist[i], yv->gcmask, &yv->values);
    }
  }
  if(mask & YglClipMask) XSetClipRectangles(D, w->chargc, 0, 0, &yv->rect, 1, Unsorted);
  if(mask & YglGCMask)   XChangeGC         (D, w->chargc, yv->gcmask, &yv->values); /* Added 960130 */
}

enum buffers { NoBuf, DBuf, MBuf, OBuf };
static int HasBuf = NoBuf;

#ifdef MULTIBUF
static int mbuf_errorhandler(Display *dpy, XErrorEvent *error) {
  char errortext[1024];
  XGetErrorText(dpy, error->error_code, errortext, 1024);
  Yprintf("doublebuffer", "X Error: %s in XmbufCreateBuffers, strange...\n",
	  errortext);
  return False;
}
#endif

void doublebuffer(void) {
  const char * MyName = "doublebuffer";
  YglWindow *w = W;
  int junk[2];
  I(MyName, "");
#ifdef OGL
  if(Ygl.UseOGL) {
    HasBuf = OBuf;
    /*if(W == NULL) { lmodtest
      prefsize(1280, 1024);
      winopen("YglFullscreen");
      w = W;
      }*/
    glDrawBuffer(GL_BACK);
    w->dbuf = True;
    return;
  }
#endif /* OGL */
#ifdef DOUBLEBUF
  if(HasBuf == NoBuf)
    HasBuf = XdbeQueryExtension(D, junk, junk+1) ? DBuf : NoBuf;
  if(HasBuf == NoBuf) {
    Yprintf(MyName, "extension 'DOUBLE-BUFFER' not found in X-Server.\n");
  } else if(HasBuf == DBuf && !w->dbuf) { /* If not inited */
    w->dbufs[0] = w->main;
    w->dbufs[1] = XdbeAllocateBackBufferName(D, w->main, XdbeUndefined);
    w->dbuf     = True;
    w->draw     = w->dbufs[1];
  }
#endif /* DOUBLEBUF */
#ifdef MULTIBUF
  if(HasBuf == NoBuf)
    HasBuf = XmbufQueryExtension(D, junk, junk+1) ? MBuf : NoBuf;
  if(HasBuf == NoBuf) {
    Yprintf(MyName, "extension 'Multi-Buffering' not found in X-Server.\n");
  } else if(HasBuf == MBuf && !w->dbuf) { /* If not inited */
    int (*old_handler)(Display*, XErrorEvent*);
    int r;
    
    old_handler = XSetErrorHandler(mbuf_errorhandler);
    r = XmbufCreateBuffers(D, IF_RGBWIN(w->win, w->main), 2,
			   MultibufferUpdateActionUndefined,
			   MultibufferUpdateHintFrequent, w->dbufs);
    XSync(D, False);
    XSetErrorHandler(old_handler);
    
    if(r != 2) {
      Yprintf(MyName, "unable to create 2 buffers.\n");
      HasBuf = NoBuf;
    } else {
      w->dbuf    = True;	/* Window is in dbuf mode */
      w->dispbuf = 0;		/* Displayed buffer */
      w->draw    = w->dbufs[1];	/* draw to backbuffer */
    }
  }
#endif /* MULTIBUF */
  if(HasBuf == NoBuf) {
    Yprintf(MyName, "Warning: cannot do doublebuffering.\n");
  }
#if !(defined(DOUBLEBUF) || defined(MULTIBUF))
  Yprintf(MyName, "Warning: Ygl is not compiled for doublebuffering.\n");
#endif
}

void singlebuffer(void) {
  const char * MyName = "singlebuffer";
  YglWindow *w = W;
  I(MyName, "");
  if(w->dbuf) switch(HasBuf) {
#ifdef OGL
  case OBuf:
    glDrawBuffer(GL_FRONT);
    break;
#endif
#ifdef DOUBLEBUF
  case DBuf:
    XdbeDeallocateBackBufferName(D, w->dbufs[1]);
    W->draw = IF_RGBWIN(W->win, W->main);
    break;
#endif
#ifdef MULTIBUF
  case MBuf:
    XmbufDestroyBuffers(D, IF_RGBWIN(w->win, w->main));
    W->draw = IF_RGBWIN(W->win, W->main);
    break;
#endif
  case NoBuf:
    /* Do nothing */
    break;
  }
  W->dbuf = False; /* Window is not in dbuf mode */
}

void swapbuffers(void) {
  const char * MyName = "swapbuffers";
  YglWindow *w = W;
  I(MyName, "");

  if(HasBuf != NoBuf && !w->dbuf) {
    Yprintf(MyName, "window %d is not in doublebuffer mode.\n", winget());
    return;
  }
  switch(HasBuf) {
#ifdef OGL
  case OBuf:
    glXSwapBuffers(D, w->main);
    break;
#endif
#ifdef DOUBLEBUF
  case DBuf:
    {
      XdbeSwapInfo xdbesi;
      xdbesi.swap_window = w->main;
      xdbesi.swap_action = XdbeUndefined;
      XdbeSwapBuffers(D, &xdbesi, 1);
      XFlush(D); /* Need it? */
      break;
    }
#endif
#ifdef MULTIBUF
  case MBuf:
    w->draw = w->dbufs[w->dispbuf];
    w->dispbuf = 1 - w->dispbuf;
    XmbufDisplayBuffers(D, 1, w->dbufs + w->dispbuf, 0, 0);
    XFlush(D); /* XmbufDisplayBuffers() seems not to flush */
    break;
#endif
  case NoBuf:
    /* Do nothing */
    break;
  }
}

void frontbuffer(Int32 bool) {
  const char * MyName = "frontbuffer";
  YglWindow *w = W;
  I(MyName, "%d", bool);
  if(HasBuf != NoBuf && !w->dbuf) {
    /* Ignore if not in doublebuffer mode */
    /* Yprintf(MyName, "window %d is not in doublebuffer mode.\n", winget()); */
    return;
  }
  switch(HasBuf) {
#ifdef OGL
  case OBuf:
    glDrawBuffer(bool ? GL_FRONT : GL_BACK);
    break;
#endif
#ifdef DOUBLEBUF
  case DBuf:
    w->draw = bool ? w->dbufs[0] : w->dbufs[1];
    break;
#endif
#ifdef MULTIBUF
  case MBuf:
    w->draw = bool ? w->dbufs[w->dispbuf] : w->dbufs[1 - w->dispbuf];
    break;
#endif
  case NoBuf:
    /* Do nothing */
    break;
  }
}

void backbuffer(Int32 bool) {
  const char * MyName = "backbuffer";
  I(MyName, "%d", bool);
  frontbuffer(!bool);
}

void gflush(void) {
  const char * MyName = "gflush";
  I(MyName, "");
  IFOGL(glFlush(), XFlush(D));
}

#ifdef OGL
static void gsync_ogl(const char *caller) {
  unsigned int retraceCount;
  static int hasExtension = -1;
  if (hasExtension == -1) {
    const char * extensionsString, * pos;
    const char * extension = "GLX_SGI_video_sync";
    extensionsString = glXQueryExtensionsString(D, YglScreen);
    pos = strstr(extensionsString, extension);
    hasExtension = pos != NULL
      && (pos == extensionsString || pos[-1] == ' ')
      && (pos[strlen(extension)] == ' ' || pos[strlen(extension)] == '\0');
    if (!hasExtension) {
      Yprintf(caller, "extension %s not present\n", extension);
    }
  }
  if (hasExtension) {
#ifdef GSvs
    glXGetVideoSyncSGI(&retraceCount);
    glXWaitVideoSyncSGI(2, (retraceCount+1)%2, &retraceCount);
#else
    Yprintf(caller, "glX*VideoSyncSGI() not enabled in Ygl\n");
#endif
  }
}
#endif

void gsync(void) {
  const char * MyName = "gsync";
  I(MyName, "");
  IFOGL(gsync_ogl(MyName),/**/);
}

Display *getXdpy(void) {
  const char * MyName = "getXdpy";
  I(MyName, "");
  return D;
}

Window getXwid(void) {
  /* Return the main window. Usable for move/resize, map/unmap,
   * event stuff */
  const char * MyName = "getXwid";
  I(MyName, "");
  return W->main;
}

#ifdef X11
Window getXdid(void) {
  /* Return the drawable. Usable for drawing. */
  const char * MyName = "getXdid";
  I(MyName, "");
  return W->draw;
}

GC getXgc(void) {
  const char * MyName = "getXgc";
  I(MyName, "");
  return W->gc;
}
#endif

void wintitle(Char8 *Title) {
  const char * MyName = "wintitle";
  I(MyName, "'%s'", Title);
  XStoreName(D, W->main, Title);
}

void winset(Int32 wid) {
  const char * MyName = "winset";
  I(MyName, "%d", wid);
  if(wid > 0 && wid <= Ygl.NextWindow && Ygl.Windows[wid].main != 0) {
    W = &Ygl.Windows[Ygl.ActiveWindow = wid];
    IFOGL(glXMakeCurrent(D, W->main, W->cx), /**/);
  } else {
    Yprintf(MyName, "invalid window id: %d\n", wid);
  }
}

Int32 winget(void) {
  const char * MyName = "winget";
  I(MyName, "");
  return Ygl.ActiveWindow;
}

Int32 getplanes(void) {
  const char * MyName = "getplanes";
  I(MyName, "");
  return W->rgb ? Ygl.RV.depth : Ygl.EmulateCmap ? EMULATE_CMAP_DEPTH : Ygl.CV.depth;
}

Int32 getvaluator(Device dev) {
  Window junkwin;
  int rx, ry, cx, cy;
  Uint mask;
  Int32 r = -1;
  const char * MyName = "getvaluator";
  I(MyName, "%d", dev);
  XQueryPointer(D, W->main, &junkwin, &junkwin, &rx, &ry, &cx, &cy, &mask);
#ifdef DEBUG
  fprintf(stderr, "getvaluator: root = (%d,%d), child = (%d,%d) mask=0x%x\n", rx, ry, cx, cy, mask);
#endif
  switch(dev) {
  case MOUSEX: r = rx; break;
  case MOUSEY: r = YglScreenHeight - ry - 1; break;
  default:     Yprintf(MyName, "unknown device: %d.\n", dev); break;
  }
  return r;
}

Int32 getbutton(Device dev) {
  Window junkwin;
  int junk;
  Uint mask, bmask;
  Int32 r = -1;
  int i, j, code;
  char keys[32];
  const char * MyName = "getbutton";
  I(MyName, "%d", dev);
  switch(dev) {
  case 0:
    break;
  case   LEFTMOUSE: bmask = Button1Mask; goto query;
  case MIDDLEMOUSE: bmask = Button2Mask; goto query;
  case  RIGHTMOUSE: bmask = Button3Mask; goto query;
  query:
    XQueryPointer(D, W->main, &junkwin, &junkwin, 
		  &junk, &junk, &junk, &junk, &mask);
#ifdef DEBUG
    fprintf(stderr, "getbutton: mask = 0x%x\n", mask);
#endif
    r = 0 != (mask & bmask);
    break;
  default:
    XQueryKeymap(D, keys);
    for(i = code = 0; i < 32; i++) for(j = 1; j < 256; j <<= 1, code++) {
      if(dev == (Ygl.keymap[code] & (KEYMAP_BIT-1))) { /* dev found */
	r = 0 != (keys[i] & j);
#ifdef DEBUG
	fprintf(stderr, "getbutton: key %d %s, device = %d\n",
		code,
		r ? "pressed" : "released",
		Ygl.keymap[code] & (KEYMAP_BIT-1));
#endif
      }
    }
    break;
  }
  if(r == -1) {
    Yprintf(MyName, "unknown device: %d.\n", dev);
  }
  return r;
}

static void set_scales(YglWindow *w) {
  w->xf = (w->vw - 1) / (w->or - w->ol);
  w->yf = (w->vh - 1) / (w->ot - w->ob);
  w->xo = w->ol - w->vl / w->xf;
  w->yo = w->ob - w->vb / w->yf;
}

void ortho2(Coord left, Coord right, Coord bottom, Coord top) {
  const char * MyName = "ortho2";
  YglWindow *w = W;
  I(MyName, "%g,%g,%g,%g", left, right, bottom, top);
  
  if(left == right || bottom == top) {
    Yprintf(MyName, "x-range or y-range is empty.\n");
    return;
  }
  
  w->ol = left;
  w->or = right;
  w->ob = bottom;
  w->ot = top;
  
  set_scales(w);
  
  IFOGL(ortho(left, right, bottom, top, -1.0, 1.0),
	;/**/);
}

#ifdef X11
static void viewport_x11(Screencoord left,   Screencoord right,
			 Screencoord bottom, Screencoord top, int clip) {
  YglWindow *w = W;
  
  if ( clip ) {
    YglGCValues yv;
    /* Clip window */
    w->clipped     = True;
    yv.rect.x      =             MIN(left, right);
    yv.rect.y      = w->ym - 1 - MAX(bottom, top);
    yv.rect.width  = w->vw;
    yv.rect.height = w->vh;
#ifdef DEBUG
    fprintf(stderr, "viewport: Clipping window with rect=(%d,%d,%d,%d).\n",
	    yv.rect.x, yv.rect.y, yv.rect.width, yv.rect.height);
#endif
    change_gc(YglClipMask, &yv);
  } else if(w->clipped) {
    /* Unclip window */
    YglGCValues yv;
    
    yv.values.clip_mask = None;
    yv.gcmask = GCClipMask;
    change_gc(YglGCMask, &yv);
    w->clipped = False;
#ifdef DEBUG
    fprintf(stderr, "viewport: Unclipping window.\n");
#endif
  }
}
#endif

void viewport(Screencoord left,   Screencoord right,
	      Screencoord bottom, Screencoord top) {
  const char * MyName = "viewport";
  YglWindow *w = W;
  int clip;
  
  I(MyName, "%d,%d,%d,%d", left, right, bottom, top);
  
  w->vl = left;
  w->vr = right;
  w->vb = bottom;
  w->vt = top;
  w->vw = ABS(right - left) + 1;
  w->vh = ABS(top - bottom) + 1;
  
  set_scales(w);

  clip = (w->vl > 0 || w->vr < w->xm - 1 || w->vb > 0 || w->vt < w->ym - 1);
  
  IFOGL(
	glViewport(w->vl, w->vb, w->vw, w->vh);
	if (clip) {
	  w->clipped = True;
	  glScissor(w->vl, w->vb, w->vw, w->vh);
	  glEnable(GL_SCISSOR_TEST);
	} else if (w->clipped) {
	  w->clipped = False;
	  glDisable(GL_SCISSOR_TEST);
	},
	viewport_x11(w->vl, w->vr, w->vb, w->vt, clip)
	);
}

void getviewport(Screencoord*left, Screencoord*right, Screencoord*bottom, Screencoord*top) {
  YglWindow *w = W;
  const char * MyName = "getviewport";
  I(MyName, "...");
  *left   = w->vl;
  *right  = w->vr;
  *bottom = w->vb;
  *top    = w->vt;
}

void reshapeviewport(void) {
  YglWindow *w = W;
  Int32 x, y;
  const char * MyName = "reshapeviewport";
  I(MyName, "");
  
  getsize(&x, &y);
  w->xm = x;
  w->ym = y;
  
  viewport(0, x - 1, 0, y - 1);
#ifdef RGBWIN
  if(w->rgb) {
    /* Ignore ExposeEvent */
    XSelectInput(D, w->win,  NoEventMask);
    XResizeWindow(D, w->win, x, y);
    XSelectInput(D, w->win,  Ygl.EventMask &  ExposureMask);
  }
#endif
}

void pushviewport(void) {
  const char * MyName = "pushviewport";
  I(MyName, "");
  IFOGL(glPushAttrib(GL_VIEWPORT_BIT),
	NI(MyName));
}

void popviewport(void) {
  const char * MyName = "popviewport";
  I(MyName, "");
  IFOGL(glPopAttrib(),
	NI(MyName));
}

void winpop(void) {
  const char * MyName = "winpop";
  I(MyName, "");
  XRaiseWindow(D, W->main);
}

void winpush(void) {
  const char * MyName = "winpush";
  I(MyName, "");
  XLowerWindow(D, W->main);
}

Int32 windepth(Int32 wid) {
  Window root, parent, *children, mainwin;
  Uint nchildren;
  Int32 n = 0;
  const char * MyName = "windepth";
  I(MyName, "%d", wid);
  
  if(wid <= 0 || wid > Ygl.NextWindow || (mainwin = Ygl.Windows[wid].main) == 0) {
    Yprintf(MyName, "invalid window id: %d\n", wid);
    return 0;
  }
  
  if(XQueryTree(D, Ygl.PWID, &root, &parent, &children, &nchildren)) {
    for(n = 0; n < nchildren && children[n] != mainwin; n++);
    n++;
    XFree((char*) children);
  }
  return n;  
}

#ifdef X11
static void linewidth_x11(Int16 w) {
  YglGCValues  yv;
  if(w == W->linewidth) return; /* already set */
  yv.gcmask = GCLineWidth; W->linewidth = yv.values.line_width = w;
  change_gc(YglGCMask, &yv);
}
#endif

void linewidth(Int16 w) {
  const char * MyName = "linewidth";
  I(MyName, "%d", w);
  IFOGL(glLineWidth(w), linewidth_x11(w));
}

Int32 getlwidth(void) {
  const char * MyName = "getlwidth";
  int r;
  I(MyName, "");
  IFOGL({GLint w;glGetIntegerv(GL_LINE_WIDTH, &w);r = w;},r = W->linewidth);
  return r;
}

typedef struct YglLineStyle_ {
  Int32 id;
  int   ndashes;
  int   repeat;
  char 	dashes[17];
} YglLineStyle;

static YglLineStyle *linestyles = NULL;
static int laststyle = -1;

static int linestyle_id2idx(Int32 id) {
  int i = laststyle;
  while(i >= 0 && linestyles[i].id != id) i--;
  return i;
}

static void init_linestyles(void) {
  YglLineStyle *ls;
  laststyle  = 0;
  linestyles = (YglLineStyle*)malloc(2 * sizeof(YglLineStyle));
  /* setup element 0 */
  ls = &linestyles[0];
  ls->id        = 0;
  ls->ndashes   = 1;
  ls->repeat    = 1;
  ls->dashes[0] = 1;
}

void deflinestyle(Int32 id, Linestyle style) {
  int i, bit, imax = 16;
  YglLineStyle *ls;
  const char * MyName = "deflinestyle";
  I(MyName, "%d,%d", id, style);
  
  if(id == 0) {
    Yprintf(MyName, "cannot modify linestyle 0.\n");
    return;
  }
  
  if(linestyles == NULL) { /* initialize */
    init_linestyles();
    i = ++laststyle;
  } else {
    i = linestyle_id2idx(id);
    if(i < 0) { /* not found */
      i = ++laststyle;
      linestyles = (YglLineStyle*)realloc(linestyles,
					  (laststyle + 1) * sizeof(YglLineStyle));
    }
  }
  
  if(linestyles == NULL) {
    Yprintf(MyName, "can't allocate memory for linestyle %d'.\n", id);
    exit(-1);
  }
  
  ls = &linestyles[i];
  ls->id        = id;
  ls->ndashes   = 0;
  ls->repeat    = 1;
  ls->dashes[0] = 0;
  bit = 1;
  
  if((style & 0xf) == ((style >> 8) & 0xf)) { /* Periodic pattern 8 bit */
    imax /= 2;
    if((style & 0x7) == ((style >> 4) & 0x7)) { /* 4 bit */
      imax /= 2;
      if((style & 0x3) == ((style >> 2) & 0x3)) { /* 2 bit */
	imax /= 2;
      }
    }
  }
  
  for(i = 0; i < imax; i++) {
    if(bit == ((style >> i) & 1)) {
      ls->dashes[ls->ndashes]++;	/* Same value */
    } else {
      ls->dashes[++ls->ndashes] = 1;
      bit ^= 1;
    }
  }
  ls->ndashes++;
  
  if(ls->dashes[0] == 0) {
    /* Pattern starts with 0s, rotate, it must start with 1s */
    int tmp = ls->dashes[1];
    ls->ndashes--;
    for(i = 0; i < ls->ndashes - 1; i++) 
      ls->dashes[i] = ls->dashes[i+2];
    ls->dashes[i] = tmp;
  }
  
#ifdef DEBUG
  fprintf(stderr, "%s: id = %d, linestyle = 0x%x, dashes = {",
	  MyName, id, style);
  for(i = 0; i < ls->ndashes; i++) fprintf(stderr, " %d", ls->dashes[i]);
  fprintf(stderr, "}\n");
#endif
}

void setlinestyle(Int32 id) {
  int idx;
  YglWindow *w = W;
  YglLineStyle *ls;
  const char * MyName = "setlinestyle";
  I(MyName, "%d", id);

  if(linestyles == NULL) { /* initialize */
    init_linestyles();
  }
  
  idx = linestyle_id2idx(id);
#ifdef DEBUG
  fprintf(stderr, "%s: id=%d, idx=%d\n", MyName, id, idx);  
#endif
  if(idx < 0) {
    Yprintf(MyName, "invalid linestyle %d.\n", id);
    return;
  }
  
  if(idx == w->linestyle) return; /* already set */
  w->linestyle = idx;
  
  if(idx == 0) {
    YglGCValues yv;
    yv.gcmask = GCLineStyle; yv.values.line_style = LineSolid;
    change_gc(YglGCMask, &yv);
  } else {
    YglGCValues yv;
    ls = &linestyles[idx];
    yv.gcmask  = GCLineStyle; yv.values.line_style = LineOnOffDash;
    yv.dashes  = ls->dashes;
    yv.ndashes = ls->ndashes;
    change_gc(YglDashMask | YglGCMask, &yv);
  }
}

Int32 getlstyle(void) {
  const char * MyName = "getlstyle";
  I(MyName, "");
  return linestyles[W->linestyle].id;
}

void lsrepeat(Int32 f) {
  int i;
  YglLineStyle *ls;
  YglGCValues yv;
  const char * MyName = "lsrepeat";
  I(MyName, "%d", f);
  
  if(f < 1 || f > 255) {
    Yprintf(MyName, "argument %d must be >= 1 and <= 255.\n", f);
    return;
  }
  
  if(W->linestyle != 0) {
    ls  = &linestyles[W->linestyle];
    for(i = 0; i < ls->ndashes; i++) 
      ls->dashes[i] = f * ls->dashes[i] / ls->repeat;
    ls->repeat = f;
    
    yv.dashes  = ls->dashes;
    yv.ndashes = ls->ndashes;
    change_gc(YglDashMask, &yv);
  }
  
#ifdef DEBUG
  fprintf(stderr, "%s: id = %d, factor = %d, dashes = {",
	  MyName, ls->id, f);
  for(i = 0; i < ls->ndashes; i++) fprintf(stderr, " %d", ls->dashes[i]);
  fprintf(stderr, "}\n");
#endif
}  

Int32 getlsrepeat(void) {
  const char * MyName = "getlsrepeat";
  I(MyName, "");
  return linestyles[W->linestyle].repeat;
}

Int32 getdisplaymode(void) {
  const char * MyName = "getdisplaymode";
  I(MyName, "");
  if(W->rgb)
    if(W->dbuf) return DMRGBDOUBLE;
    else        return DMRGB;
  else
    if(W->dbuf) return DMDOUBLE;
    else        return DMSINGLE;
}

void setbell(Char8 t) {
  XKeyboardControl xkbc;
  const char * MyName = "setbell";
  I(MyName, "%d", t);
  switch(t) {
  case 0: xkbc.bell_duration = 0;   break;
  case 1: xkbc.bell_duration = 100; break;
  case 2: xkbc.bell_duration = 400; break;
  default:
    Yprintf(MyName, "invalid value: %d.\n", t);
    return; break;
  }
  XChangeKeyboardControl(D, KBBellDuration, &xkbc);
  F;
}

void ringbell(void) {
  const char * MyName = "ringbell";
  I(MyName, "");
  XBell(D, 0);
  F;
}

Int32 getgdesc(Int32 what) {
  const char * MyName = "getgdesc";
  Int32 r = -1;
  Display *dpy;
  
  if(D == NULL) {
    if ((dpy = XOpenDisplay(NULL)) == NULL) {
      Yprintf(MyName, "can\'t open display \"%s\".\n", XDisplayName(NULL));
      exit(1);
    }
  } else {
    dpy = D;
  }
  
  switch(what) {
  case GD_XPMAX: r = DisplayWidth (dpy, DefaultScreen(dpy)); break;
  case GD_YPMAX: r = DisplayHeight(dpy, DefaultScreen(dpy)); break;
  default:
    Yprintf(MyName, "unsupported or unknown argument: %d.\n", what);
    break;
  }

  if(D == NULL) {
    XCloseDisplay(dpy);
  }

  return r;
}

void foreground(void) {
  /* Do nothing */
  Yprintf("foreground", "ignored.\n");
}

#ifdef X11
static void logicop_x11(Int32 op, const char *caller) {
  int xop;
  YglGCValues yv;
  switch(op) {
  case LO_ZERO: xop = GXclear; 		break;
  case LO_AND:  xop = GXand;   		break;
  case LO_ANDR: xop = GXandReverse; 	break;
  case LO_SRC:  xop = GXcopy;		break;
  case LO_ANDI: xop = GXandInverted;	break;
  case LO_DST:  xop = GXnoop;		break;
  case LO_XOR:  xop = GXxor;		break;
  case LO_OR:   xop = GXor;		break;
  case LO_NOR:  xop = GXnor;		break;
  case LO_XNOR: xop = GXequiv;		break;
  case LO_NDST: xop = GXinvert;		break;
  case LO_ORR:  xop = GXorReverse;	break;
  case LO_NSRC: xop = GXcopyInverted;	break;
  case LO_ORI:  xop = GXorInverted;	break;
  case LO_NAND: xop = GXnand;		break;
  case LO_ONE:  xop = GXset;		break;
    
  case LO_MIN:
  case LO_MAX:
  case LO_AVG:
  case LO_DMS:
  case LO_SMD:
  case LO_SUM:
    Yprintf(caller, "unsupported argument: %d.\n", op);
    return;
  default:
    Yprintf(caller, "unknown argument: %d.\n", op);
    return;
  }
  
  yv.gcmask = GCFunction; yv.values.function = xop;  
  change_gc(YglGCMask, &yv);
}
#endif

#ifdef OGL
static void logicop_ogl(Int32 op, const char *caller) {
  GLenum oop;
  switch(op) {
  case LO_ZERO: oop = GL_CLEAR; 	break;
  case LO_AND:  oop = GL_AND;   	break;
  case LO_ANDR: oop = GL_AND_REVERSE; 	break;
  case LO_SRC:  oop = GL_COPY;		break;
  case LO_ANDI: oop = GL_AND_INVERTED;	break;
  case LO_DST:  oop = GL_NOOP;		break;
  case LO_XOR:  oop = GL_XOR;		break;
  case LO_OR:   oop = GL_OR;		break;
  case LO_NOR:  oop = GL_NOR;		break;
  case LO_XNOR: oop = GL_EQUIV;		break;
  case LO_NDST: oop = GL_INVERT;	break;
  case LO_ORR:  oop = GL_OR_REVERSE;	break;
  case LO_NSRC: oop = GL_COPY_INVERTED;	break;
  case LO_ORI:  oop = GL_OR_INVERTED;	break;
  case LO_NAND: oop = GL_NAND;		break;
  case LO_ONE:  oop = GL_SET;		break;
    
  case LO_MIN:
  case LO_MAX:
  case LO_AVG:
  case LO_DMS:
  case LO_SMD:
  case LO_SUM:
    Yprintf(caller, "unsupported argument: %d.\n", op);
    return;
  default:
    Yprintf(caller, "unknown argument: %d.\n", op);
    return;
  }
  glLogicOp(oop);
  if(oop != GL_COPY) {
    glEnable(GL_LOGIC_OP);
  } else {
    glDisable(GL_LOGIC_OP);
  }
}
#endif

void logicop(Int32 op) {
  const char * MyName = "logicop";
  I(MyName, "%d", op);
  IFOGL(logicop_ogl(op, MyName),
	logicop_x11(op, MyName));
}

#ifdef OGL
static void getmatrix_ogl(const char *caller, Matrix m) {
  GLint mode, what;
  GLfloat glm[16];
  short i, j;
  
  glGetIntegerv(GL_MATRIX_MODE, &mode);
  switch(mode) {
  case GL_MODELVIEW:  what = GL_MODELVIEW_MATRIX; break;
  case GL_PROJECTION: what = GL_PROJECTION_MATRIX; break;
  default:
    Yprintf(caller, "strange OpenGL matrix mode %d\n", mode);
    return;
  }
  glGetFloatv(what, glm);
  /* Under OpenGL matrices are in column order, transpose... */
  for(i = 0; i < 4; i++) for(j = 0; j < 4; j++) m[j][i] = glm[4*i + j];
}
#endif

#ifdef X11
static void getmatrix_x11(Matrix m) {
  /* We only have orthographic views, so... */
  YglWindow *w = W;
  memset(m, 0, sizeof(Matrix));
  m[0][0] =  2.0 / (w->or - w->ol);
  m[1][1] =  2.0 / (w->ot - w->ob);
  m[2][2] = -1.0;
  m[3][0] = -(w->or + w->ol)/(w->or - w->ol);
  m[3][1] = -(w->ot + w->ob)/(w->ot - w->ob);
  m[3][3] =  1.0;
}
#endif

void getmatrix(Matrix m) {
  const char * MyName = "getmatrix";
  I(MyName, "Matrix");
  IFOGL(getmatrix_ogl(MyName, m), getmatrix_x11(m));
}

#ifdef COVERSLEEP

#if 0

/* Uncomment next line if you have problems with usleep() in old (prior 1998) versions of Alpha Linux */
/* #define ALPHA_LINUX_USLEEP_BUG */

#ifdef __useconds_t_defined
# define USLEEP_RET_TYPE int
# define USLEEP_RET_VAL  0
# define USLEEP_ARG_TYPE useconds_t
#else /* __useconds_t_defined */

# if defined(__linux)

#  if defined(__i386) || defined(__alpha)
#   if defined(ALPHA_LINUX_USLEEP_BUG)
/* Under old DEC Alpha Linux, it was "unsigned int usleep(unsigned int)"...
 * Thanks to Kai Froese <froese@traf8.uni-duisburg.de> */
#    define USLEEP_RET_TYPE unsigned int
#    define USLEEP_RET_VAL  0
#    define USLEEP_ARG_TYPE unsigned int
#   elif defined(__GLIBC__)
/* In new GNU libc, usleep is "void usleep(unsigned int)"
 * Thanks to Juergen Holm <holm@Theorie.Physik.UNI-Goettingen.DE> */
#    define USLEEP_RET_TYPE void
#    define USLEEP_RET_VAL
#    define USLEEP_ARG_TYPE unsigned int
#   else
/* In old libc, usleep was "void usleep(unsigned long)" */
#    define USLEEP_RET_TYPE void
#    define USLEEP_RET_VAL
#    define USLEEP_ARG_TYPE unsigned long
#   endif
#  else /* defined(__i386) || defined(__alpha) */
/* Unknown linux derivate */
#   error "Unknown Linux derivate, please report to author"
#  endif

# elif defined(__FreeBSD__)

/* Under FreeBSD, usleep is "void usleep(unsigned int)"
 * Thanx to Pedro Giffuni <pgiffuni@fps.biblos.unal.edu.co>
 * and      Stephen Kennedy <steve@maths.tcd.ie> */
#  define USLEEP_RET_TYPE void
#  define USLEEP_RET_VAL
#  define USLEEP_ARG_TYPE unsigned int

# else /* defined(__linux) || defined(__FreeBSD__) */

/* On all (?) other platforms, usleep is "int usleep(unsigned int)" */
#  define USLEEP_RET_TYPE int
#  define USLEEP_RET_VAL  0
#  define USLEEP_ARG_TYPE unsigned int

# endif /* defined (__linux) */

#endif /* __useconds_t_defined */

#endif /* 0 */

#include "usleep.h"
#ifndef USLEEP_RET_TYPE
# define USLEEP_RET_TYPE int
# define USLEEP_RET_VAL  0
#endif
#ifndef USLEEP_ARG_TYPE
# define USLEEP_ARG_TYPE unsigned int
#endif

USLEEP_RET_TYPE usleep(USLEEP_ARG_TYPE Useconds) {
  struct timeval tmout;
  if(D != NULL) IFOGL(glFlush(), XFlush(D));
  tmout.tv_usec = Useconds % 1000000;
  tmout.tv_sec  = Useconds / 1000000;
  (void) select(0, NULL, NULL, NULL, &tmout);
  return USLEEP_RET_VAL;
}

Uint sleep(Uint seconds) {
  struct timeval tmout;
  if(D != NULL) IFOGL(glFlush(), XFlush(D));
  tmout.tv_usec = 0;
  tmout.tv_sec  = seconds;
  (void) select(0, NULL, NULL, NULL, &tmout);
  return 0;
}
#endif
