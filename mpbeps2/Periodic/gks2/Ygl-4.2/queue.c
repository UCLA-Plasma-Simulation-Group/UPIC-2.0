/*
 *    Ygl: Run GL programs with standard X11 and/or OpenGL routines.
 *    (C) Fred Hucht 1993-2006
 *    EMail: fred<at>thp.Uni-Duisburg.de
 */

static const char vcid[] = "$Id: queue.c,v 4.8 2007-05-08 13:27:54+02 fred Exp $";

#include "header.h"
#include <setjmp.h>

static Ulong YglButtonMask = 0;

static int   change_map(Device dev, int on_off);
static void  set_YglEventMask(void);
static Bool  more_redraw(Display*, XEvent*, char *);
static int   check_more_redraw(Int32 wid);
static Int32 q_next(int parent, Int16 *val);

/* XChangeActivePointerGrab only understands these: */
#define XCAPG_MASK (ButtonPressMask|ButtonReleaseMask|\
		    EnterWindowMask|LeaveWindowMask|\
		    Button1MotionMask|Button2MotionMask|\
		    Button3MotionMask|Button4MotionMask|\
		    Button5MotionMask|PointerMotionHintMask|\
		    PointerMotionMask|ButtonMotionMask|KeymapStateMask)

#define KeyMask    (KeyPressMask | KeyReleaseMask)
#define ButtonMask (ButtonPressMask | ButtonReleaseMask)

static char IsQueued[MAXYGLDEVICE];

Int32 isqueued(Int16 dev) {
  return dev >= 0 && dev < MAXYGLDEVICE && IsQueued[dev];
}

static int change_map(Device dev, int on_off) {
  int code, found = False;
  for (code = 0; code < KEYMAP_LEN; code++) /* several codes may return same keysym (eg MacX) */
    if (dev == (Ygl.keymap[code] & (KEYMAP_BIT-1))) { /* found */
      found = True;
      if (on_off) Ygl.keymap[code] |=  KEYMAP_BIT;
      else        Ygl.keymap[code] &= ~KEYMAP_BIT;
#ifdef DEBUG
      fprintf(stderr, "change_map: dev = %d (%s), code = %d, new Ygl.keymap[code] = %d\n",
	      dev, Ygl.devicename(dev), code, Ygl.keymap[code]);
#endif
    }
  return found;
}

static void set_YglEventMask(void) {
  int i;
  
  /* To change event mask even when a button is pressed: */
  XChangeActivePointerGrab(D, Ygl.EventMask & XCAPG_MASK, None, CurrentTime);
  
#ifdef DEBUG
  fprintf(stderr, "set_YglEventMask: EventMask = 0x%x, IsQueued[ANYKEY] = %d\n",
	  Ygl.EventMask, IsQueued[ANYKEY]);
#endif
  
  for (i = 1; i < Ygl.NextWindow; i++) {
    YglWindow *w = &Ygl.Windows[i];
    if (w->main != 0) {
      XSetWMProtocols(D, w->main, &Ygl.wm_dw, Ygl.wm_dw_flag ? 1 : 0);
      XSelectInput(D, w->main, Ygl.EventMask);
#ifdef RGBWIN
      if (w->rgb) /* Expose events don't propagate */
	XSelectInput(D, w->win,  Ygl.EventMask &  ExposureMask);
#endif
    }
  }
}

static int YglTie[3][2] = {{0,0},{0,0},{0,0}}; /* ties to left, middle, right button */

void tie(Device button, Device val1, Device val2) {
  const char * MyName = "tie";
  int tieb;
  I(MyName, "%d,%d,%d", button, val1, val2);
  
  switch (button) {
  case LEFTMOUSE:   tieb = 0; break;
  case MIDDLEMOUSE: tieb = 1; break;
  case RIGHTMOUSE:  tieb = 2; break;
  default:
    Yprintf(MyName, "first argument is not a button: %d.\n", button); return;
  }

  switch (val1) {
  case MOUSEX: break;
  case MOUSEY: break;
  default:
    Yprintf(MyName, "second argument is not a valuator: %d.\n", val1); return;
  }
  
  switch (val2) {
  case MOUSEX: break;
  case MOUSEY: break;
  default:
    Yprintf(MyName, "third argument is not a valuator: %d.\n", val2); return;
  }
  
  YglTie[tieb][0] = val1;
  YglTie[tieb][1] = val2;
}

static int got_sigalrm     = 0;
static int in_nextevent    = False;
static jmp_buf alrm_context;

static void sigalrm_handler(int sig_no) {
#ifdef DEBUG
  struct timeval tv;
  struct timezone tz;
  static double tm, otm;
  gettimeofday(&tv, NULL);
  tm = tv.tv_sec + 1e-6 * tv.tv_usec;
  fprintf(stderr, "In sigalrm_handler after %.6f sec\n", tm - otm);
  otm = tm;
#endif
  got_sigalrm++;
  /* Bump out of XNextEvent() */
  if (in_nextevent) {
    in_nextevent = False;
    longjmp(alrm_context, 1);
  }
}

void noise(Device dev, Int16 delta) {
  const char * MyName = "noise";
  I(MyName, "%d,%d", dev, delta);
  if (isqueued(dev)) {
    switch (dev) {
    case TIMER0:
      {
	Ulong usec = 20000L * delta; /* default delta: 1 -> 50Hz */
	struct itimerval itimer, otimer;
	itimer.it_interval.tv_sec  = usec / 1000000;
	itimer.it_interval.tv_usec = usec % 1000000;
	itimer.it_value.tv_sec     = usec / 1000000;
	itimer.it_value.tv_usec    = usec % 1000000;
#ifdef DEBUG
	fprintf(stderr, "ALRM time: %g sec\n",
		itimer.it_value.tv_sec + 1e-6 * itimer.it_value.tv_usec);
#endif
	setitimer(ITIMER_REAL, &itimer, &otimer);
      }
      break;
    default:
      Yprintf(MyName, "device %d (%s) not implemented.\n", dev, Ygl.devicename(dev));
      break;
    }
  } else {
    Yprintf(MyName, "device %d (%s) not queued.\n", dev, Ygl.devicename(dev));
  }
}

void qdevice(Device dev) {
  const char * MyName = "qdevice";
  
  I(MyName, "%d", dev);
  switch (dev) {
  case WINCLOSE:    
    Yprintf(MyName, "device WINCLOSE not longer supported. Using WINQUIT instead...\n");
  case WINQUIT:      /* This is the event GL sends for WM_DELETE_WINDOW... */
                     Ygl.wm_dw_flag = True;              break;
  case REDRAW:       Ygl.EventMask |= RedrawMask;        break;
  case INPUTCHANGE:  Ygl.EventMask |= EnterLeaveMask;    break;
  case MOUSEX:
  case MOUSEY:       Ygl.EventMask |= PointerMotionMask; break;
  case LEFTMOUSE:    YglButtonMask |= Button1Mask;       break;
  case MIDDLEMOUSE:  YglButtonMask |= Button2Mask;       break;
  case RIGHTMOUSE:   YglButtonMask |= Button3Mask;       break;
  case WHEELUP:      YglButtonMask |= Button4Mask;       break;
  case WHEELDOWN:    YglButtonMask |= Button5Mask;       break;
  case KEYBD:        Ygl.EventMask |= KeyMask;           break;
  case TIMER0:
    {
      struct sigaction action, oaction;
      action.sa_handler = sigalrm_handler;
      sigfillset(&action.sa_mask);
      action.sa_flags   = 0;
      if (sigaction(SIGALRM, &action, &oaction)) {
	perror(MyName);
	exit(1);
      }
      if (oaction.sa_handler)
	Yprintf(MyName, "TIMER0: warning: SIGALRM already used.\n");
      IsQueued[TIMER0] = True;
      noise(TIMER0, 1);
    }
    break;
  default:
    if (change_map(dev, True)) {
      IsQueued[ANYKEY] = True;
      Ygl.EventMask |= KeyMask;
    } else {
      Yprintf(MyName, "unsupported device: %d (%s).\n", dev, Ygl.devicename(dev));
      return;
    }
    break;
  }
  
  if (YglButtonMask) Ygl.EventMask |= ButtonMask;
  
  set_YglEventMask();
  IsQueued[dev] = True;
}

void unqdevice(Device dev) {
  int i;
  const char * MyName = "unqdevice";
  
  I(MyName, "%d", dev);
  switch (dev) {
  case WINCLOSE:    
    Yprintf(MyName, "device WINCLOSE not longer supported. Using WINQUIT instead...\n");
  case WINQUIT:      /* This is the event GL sends for WM_DELETE_WINDOW... */
                     Ygl.wm_dw_flag = False;              break;
  case REDRAW:       Ygl.EventMask &= ~RedrawMask;        break;
  case INPUTCHANGE:  Ygl.EventMask &= ~EnterLeaveMask;    break;
  case MOUSEX:
  case MOUSEY:       Ygl.EventMask &= ~PointerMotionMask; break;
  case LEFTMOUSE:    YglButtonMask &= ~Button1Mask;       break;
  case MIDDLEMOUSE:  YglButtonMask &= ~Button2Mask;       break;
  case RIGHTMOUSE:   YglButtonMask &= ~Button3Mask;       break;
  case WHEELUP:      YglButtonMask &= ~Button4Mask;       break;
  case WHEELDOWN:    YglButtonMask &= ~Button5Mask;       break;
  case KEYBD: if (!IsQueued[ANYKEY]) Ygl.EventMask &= ~KeyMask; break;
  case TIMER0:
    {
      struct sigaction action;
      action.sa_handler = SIG_DFL;
      sigfillset(&action.sa_mask);
      action.sa_flags   = 0;
      noise(TIMER0, 0);
      if (sigaction(SIGALRM, &action, NULL)) {
	perror(MyName);
	exit(1);
      }
    }
    break;
  default:
    if (!change_map(dev, False)) {
      Yprintf(MyName, "unsupported device: %d (%s).\n", dev, Ygl.devicename(dev));
      return;
    }
    IsQueued[ANYKEY] = False;
    for (i = 0; i < KEYMAP_LEN; i++)
      IsQueued[ANYKEY] = IsQueued[ANYKEY] || (Ygl.keymap[i] & KEYMAP_BIT);
    if (!IsQueued[ANYKEY] && !IsQueued[KEYBD])
      Ygl.EventMask &= ~KeyMask;
    break;
  }
  
  if (!YglButtonMask) Ygl.EventMask &= ~ButtonMask;
  
  set_YglEventMask();
  IsQueued[dev] = False;
}

void qreset(void) {
  const char * MyName = "qreset";
  
  I(MyName, "");
  XSync(D, True);
  /* Don't remove saved mouse/key events, else we may never get 
   * MOUSEY/ANYKEY events. */
  /*
    mousex = -1; mousey = -1;
    mousexchanged = False;
    mouseychanged = False;
    keyevent = 0;
    */
}

struct r_info {
  Window main;
#ifdef RGBWIN
  Window win;
#endif
  int    num;
  int    ret;
};

static Bool more_redraw(Display *dpy, XEvent *e, char *rc) {
  struct r_info *r = (struct r_info*) rc;
  int ret;
  if ((e->type == Expose || e->type == ConfigureNotify) &&
      (e->xany.window == r->main || IF_RGBWIN(e->xany.window == r->win, False))) {
    r->ret = True; /* found Expose or ConfigureNotify event for same window */
  }
  ret = r->ret || --r->num == 0; /* return True if found or queue empty */
#ifdef DEBUG
#ifdef RGBWIN
  fprintf(stderr, 
	  "more_redraw: type=%s, r={main=0x%x,win=0x%x,num=%d,ret=%d}, ret=%d\n",
	  Ygl.XEventNames[e->type], r->main, r->win, r->num, r->ret, ret);
#else
  fprintf(stderr, 
	  "more_redraw: type=%s, r={main=0x%x,num=%d,ret=%d}, ret=%d\n",
	  Ygl.XEventNames[e->type], r->main, r->num, r->ret, ret);
#endif
#endif
  return ret;
}

static int check_more_redraw(Int32 wid) {
  XEvent p;
  struct r_info info;
  info.main = Ygl.Windows[wid].main;
#ifdef RGBWIN
  info.win  = Ygl.Windows[wid].win;
#endif
  info.num  = XPending(D);
  info.ret  = False;
  if (info.num) {
#ifdef DEBUG
    fprintf(stderr, "check_more_redraw: Calling XPeekIfEvent()\n");
#endif
    XPeekIfEvent(D, &p, more_redraw, (char*)&info);
  }
  return info.ret;
}

#define QREAD 1
#define QTEST 2

static Int32 q_next(int parent, Int16 *val) {
  XEvent e;
  Int32 dev = -1, dev2;
  int tieb = -1;
  char key;
  static Int16 keyevent = 0;
  static Int16 mousexchanged = False;
  static Int16 mouseychanged = False;
  static Int16 mousex = -1;
  static Int16 mousey = -1;
  char *caller = parent == QTEST ? "qtest" : "qread";
  
  I(caller, "");
  
  /* Let TIMER0, i.e. SIGALRM, bump out of XNextEvent() */
  if (parent == QREAD && IsQueued[TIMER0]) if (setjmp(alrm_context)) {
    /* if we are here, we came from longjmp() */
#ifdef DEBUG
    fprintf(stderr, "SIGALRM received in XNextEvent.\n");
#endif
  }
  
  while (dev == -1) {
    if (got_sigalrm) {
      dev  = TIMER0;
      *val = got_sigalrm;
      if (parent == QREAD) got_sigalrm = 0;
    } else if (mousexchanged) {	/* Hack for mouse events */
      dev  = MOUSEX;
      *val = mousex;
      if (parent == QREAD) mousexchanged = False;
    } else if (mouseychanged) {	/* Hack for mouse events */
      dev  = MOUSEY;
      *val = YglScreenHeight - mousey - 1;
      if (parent == QREAD) mouseychanged = False;
    } else if (keyevent > 0) {	/* Hack for reporting both KEYBD and ANYKEY device */
      dev  = keyevent;
      *val = 1; /* Must be 1 (KeyPress) */
      if (parent == QREAD) keyevent = 0;
    } else if (parent == QTEST && XPending(D) == 0) {
      dev = 0;
    } else {
      
      if (parent == QREAD) {
	IFOGL(glFlush(), XFlush(D));
	in_nextevent = True;
	XNextEvent(D, &e);
	in_nextevent = False;
      } else {
	XPeekEvent(D, &e);
      }
#ifdef DEBUG
      fprintf(stderr, 
	      "%s: e.type = %s for window 0x%x, XPending returns %d\n",
	      caller, Ygl.XEventNames[e.type], e.xany.window, XPending(D));
#endif
      switch (e.type) {
	Int32 wid;
      case ConfigureNotify:
      case Expose:
#ifdef RGBWIN
	wid = Ygl.x2gl_wid(e.xany.window, X2GL_WIN);
	if (wid == 0) wid = Ygl.x2gl_wid(e.xany.window, X2GL_MAIN);
#else
	wid = Ygl.x2gl_wid(e.xany.window, X2GL_MAIN);
#endif
	
	if (wid == 0) Yprintf(caller, "REDRAW received for unknown window\n");
	
	if (parent == QTEST) {
	  dev = REDRAW;
	} else if (e.type == Expose && e.xexpose.count != 0) { /* ignore */
#ifdef DEBUG
	  fprintf(stderr,
		  "%s: e.xexpose.count = %d\n",
		  caller, e.xexpose.count);
#endif
	} else if (check_more_redraw(wid)) {
	  /* more Expose events for the same window in queue? */
#ifdef DEBUG
	  fprintf(stderr,
		  "More Expose or ConfigureNotify events for same window in queue, ignoring this one...\n",
		  caller);
#endif
	} else {
	  dev  = REDRAW;
	  *val = wid;
	}
	break;
	
      case EnterNotify:
	dev  = INPUTCHANGE;
	*val = Ygl.x2gl_wid(e.xany.window, X2GL_MAIN);
	break;
	
      case LeaveNotify:
	dev  = INPUTCHANGE;
	*val = 0;
	break;
	
      case MotionNotify:
#ifdef DEBUG
	fprintf(stderr, "%s: e.xmotion.x_root = %d, e.xmotion.y_root = %d\n",
		/**/ caller, e.xmotion.x_root,      e.xmotion.y_root);
#endif
	if (parent == QTEST) {
	  if (mousex != e.xmotion.x_root)
	    dev = MOUSEX;	/* if both values have changed, report x first... */
	  else
	    dev = MOUSEY;
	} else {
	  if (mousey != e.xmotion.y_root) {
	    mousey = e.xmotion.y_root;
	    mouseychanged = True;
	  }
	  if (mousex != e.xmotion.x_root) {
	    dev  = MOUSEX;
	    *val = mousex = e.xmotion.x_root;
	  } 
	}
	break;
	
      case KeyPress:
      case KeyRelease: /* Note: KeyReleases are not reported as KEYBD events */
#ifdef DEBUG
	fprintf(stderr, "%s: e.xkey.keycode = %d, Ygl.keymap[%d] = %d\n",
		caller, e.xkey.keycode, e.xkey.keycode, Ygl.keymap[e.xkey.keycode]);
#endif
	if (IsQueued[KEYBD] &&
	    e.type == KeyPress &&
	    1 == XLookupString(&e.xkey, &key, 1, NULL, NULL)) {
	  /* return code */
	  dev  = KEYBD;
	  *val = (Int16)key;
	}
	if (IsQueued[ANYKEY] &&
	    (dev2 = Ygl.keymap[e.xkey.keycode]) & KEYMAP_BIT ) {
	  /* Individual keys requested? */
	  dev2 &= (KEYMAP_BIT-1);
	  if (dev == KEYBD) {	/* first report KEYBD event */
	    if (parent == QREAD) keyevent = dev2;
	  } else {
	    dev  = dev2;
	    *val = e.type == KeyPress ? 1 : 0;
	  }
	}
      	break;
	
      case ButtonPress:
      case ButtonRelease:
	switch (e.xbutton.button) {
	case Button1:
	  if (YglButtonMask & Button1Mask) {
	    dev = LEFTMOUSE;
	    if (parent == QREAD) tieb = 0;
	  }
	  break;
	case Button2:
	  if (YglButtonMask & Button2Mask) {
	    dev = MIDDLEMOUSE;
	    if (parent == QREAD) tieb = 1;
	  }
	  break;
	case Button3:
	  if (YglButtonMask & Button3Mask) {
	    dev = RIGHTMOUSE;
	    if (parent == QREAD) tieb = 2;
	  }
	  break;
	case Button4:
	  if (YglButtonMask & Button4Mask) {
	    dev = WHEELUP;
	  }
	  break;
	case Button5:
	  if (YglButtonMask & Button5Mask) {
	    dev = WHEELDOWN;
	  }
	  break;
	default:
	  Yprintf(caller, "unknown button: %d.\n", e.xbutton.button);
	  break;
	}
	
	if (parent == QREAD && dev != -1) {
	  *val = e.type == ButtonPress ? 1 : 0;
	  
	  if (tieb != -1) {
	    switch (YglTie[tieb][0]) {
	    case MOUSEX: mousexchanged = True; mousex = e.xbutton.x_root; break;
	    case MOUSEY: mouseychanged = True; mousey = e.xbutton.y_root; break;
	    }
	    switch (YglTie[tieb][1]) {
	    case MOUSEX: mousexchanged = True; mousex = e.xbutton.x_root; break;
	    case MOUSEY: mouseychanged = True; mousey = e.xbutton.y_root; break;
	    }
	  }
	}
	break;
	
      case ClientMessage: 
#ifdef DEBUG
	fprintf(stderr, "%s: e.xclient.message_type = %d\n",
		caller, e.xclient.message_type);
#endif
	dev  = WINQUIT;
	*val = Ygl.x2gl_wid(e.xany.window, X2GL_MAIN);
	break;
	
      case CirculateNotify:
      case CreateNotify:
      case DestroyNotify:
      case GravityNotify:
      case ReparentNotify:
      case MapNotify:
      case UnmapNotify:
	/* Ignore these Events... */
	break;
      case MappingNotify:
	/* We should reload Ygl.keymap[] here, but who cares... */
	break;
      default:
	Yprintf(caller,
		"unknown event (type=%d) for window %d(0x%x) received.\n",
		e.type,	Ygl.x2gl_wid(e.xany.window, X2GL_MAIN), e.xany.window);
	break;
      } /* end switch */
      
      if (dev == -1 && parent == QTEST) {
#ifdef DEBUG
	fprintf(stderr, "qtest: eating event: e.type = %d(%s)\n",
		e.type, Ygl.XEventNames[e.type]);
#endif
	XNextEvent(D, &e); /* read and ignore unreported events so they don't block the queue */
      }
      
      if (dev != -1 && parent == QREAD) Ygl.lastreadevent = e;
    }
  }
#ifdef DEBUG
  if (dev != 0) fprintf(stderr, "%s: returning dev = %d (%s), *val = %d\n",
			caller, dev, Ygl.devicename(dev), *val);
#endif
  return dev;
}

Int32 qtest(void) { 
  Int16 val;
  return q_next(QTEST, &val);
}

Int32 qread(Int16 *val) {
  return q_next(QREAD, val);
}

void qenter(Int16 dev, Int16 val) {
  XEvent e;
  const char * MyName = "qenter";
  
  I(MyName, "%d,%d", dev, val);

  switch (dev) {
  case REDRAW:
    e.type = Expose;
    /* if RGBWIN, send to .win, not .main */
    e.xexpose.window = (val > 0 && val < Ygl.NextWindow) ? Ygl.Windows[val].IF_RGBWIN(win,main) : 0;
    e.xexpose.count = 0;
    break;
#if 0
  case INPUTCHANGE:  Yprintf(MyName, "unsupported device: INPUTCHANGE.\n") ; return; break;
    /* e.type = EnterNotify; e.xcrossing.window = val; break; */
  case KEYBD:        Yprintf(MyName, "unsupported device: KEYBD.\n")       ; return; break;
  case MENUBUTTON:   Yprintf(MyName, "unsupported device: MENUBUTTON.\n")  ; return; break;
#endif
  default:
    Yprintf(MyName, "unsupported device: %d (%s).\n", dev, Ygl.devicename(dev));
    return;
    break;
  }
  XPutBackEvent(D, &e);
}

#ifdef OGL

/* picking/names stuff */
static Int32 pick_bufferlen, mySelectBufferLen;
static GLuint *mySelectBuffer;
static Int16 picksize_deltax = 10;
static Int16 picksize_deltay = 10;

void picksize(Int16 deltax, Int16 deltay) {
  const char * MyName = "picksize";
  I(MyName, "%d,%d", deltax, deltay);
  if(!Ygl.UseOGL) NI(MyName);
  picksize_deltax = deltax;
  picksize_deltay = deltay;
}

void pick(Int16 buffer[], Int32 bufferlen) {
  const char * MyName = "pick";
  I(MyName, "0x%x,%d", buffer, bufferlen);
  if(!Ygl.UseOGL) NI(MyName);
  
  if (Ygl.SelectMode) {
    Yprintf(MyName, "already in picking mode\n");
    return;
  }
  
  pick_bufferlen    =     bufferlen; /* save for endpick() */
  mySelectBufferLen = 4 * bufferlen; /* should be sufficient */
  
  mySelectBuffer = (GLuint*) calloc(mySelectBufferLen, sizeof(GLuint));
  if (mySelectBuffer == NULL) {
    Yprintf(MyName, "can't allocate memory.\n");
    exit(1);
  }
  
  glSelectBuffer(mySelectBufferLen, mySelectBuffer);

#ifndef VIEWPICK
  glRenderMode(GL_SELECT);
#endif
  Ygl.SelectMode = True;
  initnames();
  
  {
    Window junkwin;
    int rx, ry, cx, cy;
    Uint mask;
    GLint viewport[4];
    GLdouble x, y, delx, dely;
    
    glGetIntegerv(GL_VIEWPORT, viewport);
    
    XQueryPointer(D, W->main, &junkwin, &junkwin, &rx, &ry, &cx, &cy, &mask);
#ifdef DEBUG
    fprintf(stderr, "pick: root = (%d,%d), child = (%d,%d) mask=0x%x\n",
	    rx, ry, cx, cy, mask);
#endif
    x    = (double)cx;
    y    = (double)viewport[3] - cy;
    delx = (double)picksize_deltax;
    dely = (double)picksize_deltay;
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
#ifdef DEBUG
    fprintf(stderr, "pick: gluPickMatrix(%g, %g, %g, %g, (%d,%d,%d,%d))\n",
	    x, y, delx, dely,
	    viewport[0],viewport[1],viewport[2],viewport[3]);
#endif
    gluPickMatrix(x, y, delx, dely, viewport);
    glMatrixMode(GL_MODELVIEW);
  }
}

Int32 endpick(Int16 buffer[]) {
  const char * MyName = "endpick";
  Int32 hits, i, bi, mi;
  I(MyName, "0x%x", buffer);
  if(!Ygl.UseOGL) NIR(MyName, 0);

  if (!Ygl.SelectMode) {
    Yprintf(MyName, "not in picking mode\n");
    return 0;
  }
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix(); /* pushed in pick() */
  glMatrixMode(GL_MODELVIEW);
  
  Ygl.SelectMode = False;
  hits = glRenderMode(GL_RENDER);
#ifdef VIEWPICK
  swapbuffers();
#endif
  
#ifdef DEBUG
  Yprintf(MyName, "mySelectBuffer[0..%d] : ", MIN(mySelectBufferLen, 20) - 1);
  for(i = 0; i < MIN(mySelectBufferLen, 20); i++) {
    fprintf(stderr, "%d ", mySelectBuffer[i]);
  }
  fprintf(stderr, "\n");
#endif
  
#define BUF_CHK() (bi < pick_bufferlen && mi < mySelectBufferLen)
  /* copy select buffer from OpenGL to GL format */
  for(i = bi = mi = 0; i < hits; i++) {
    GLuint depth = 0; /* stack depth */
    if (BUF_CHK()) depth = buffer[bi++] = mySelectBuffer[mi++];
    else hits = -i; /* return -valid hits if buffer too small */
    mi += 2; /* skip depth info */
    for (; depth > 0; depth--) {
      if (BUF_CHK()) buffer[bi++] = mySelectBuffer[mi++];
      else hits = -i;
    }
  }
#undef BUF_CHK
  if (mi >= mySelectBufferLen) { /* should not happen... */
    Yprintf(MyName, "mySelectBuffer too small, this should not happen\n");
  }
  
  free(mySelectBuffer);

#ifdef DEBUG
  Yprintf(MyName, "hits = %d, buffer[0..%d] : ", hits, bi - 1);
  for(i = 0; i < bi; i++) {
    fprintf(stderr, "%d ", buffer[i]);
  }
  fprintf(stderr, "\n");
#endif
  
  return hits;
}

void initnames(void) {
  const char * MyName = "initnames";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  glInitNames();
  glPushName(-1); /* Why do we have to push someting onto the stack? */
}

void loadname(Int16 name) {
  const char * MyName = "loadname";
  I(MyName, "%d", name);
  if(!Ygl.UseOGL) NI(MyName);
  glLoadName(name);
}

void pushname(Int16 name){
  const char * MyName = "pushname";
  I(MyName, "%d", name);
  if(!Ygl.UseOGL) NI(MyName);
  glPushName(name);
}

void popname(void){
  const char * MyName = "popname";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  glPopName();
}

#endif  /* OGL */
