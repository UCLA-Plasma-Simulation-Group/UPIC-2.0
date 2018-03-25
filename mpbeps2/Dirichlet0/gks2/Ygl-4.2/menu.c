/*
 *    Ygl: Run GL programs with standard X11 routines.
 *    (C) Fred Hucht 1993-2007
 *    EMail: fred@thp.Uni-Duisburg.de
 */

static const char vcid[] = "$Id: menu.c,v 3.6 2007-05-08 13:25:35+02 fred Exp $";

#include "header.h"
/* #include <X11/bitmaps/gray1> */
#include "X11/gray1.h"
#include <X11/cursorfont.h>

#define MENUPADW 	14
#define MENUPADH 	4
#define MENUSHADE 	8
#define SUBMENUOFFSET	10
#define ON		1
#define OFF		0

typedef Int32(*YglMenuCallback)(Int32, ...);

typedef struct YglMenu_ {
  Window		win;
  int			nitem;
  int 			x, y;
  Uint			w, h;
  struct YglMenuItem_ 	*longest;
  struct YglMenuItem_ 	**item;
} YglMenu;

typedef struct YglMenuItem_ {
  char 			*name;
  YglMenuCallback	callback;
  struct YglMenu_ 	*submenu;
  Int32			retval;
  Int32			mode;
  Uint 			callF:1;
  Uint			line:1;
} YglMenuItem;

static YglMenu 	**YglMenus = NULL;
static int 	YglLastMenu;
static GC	menugc = NULL;
static Cursor 	menucursor;
static int 	fontascent;
static int	menuheight;
static int 	white, black;

static YglMenu 	*mid_to_m(const char *caller, Int32 mid);
static void 	add_to_pup(const char *caller, Int32 mid, const char *string, va_list argp);
static void 	create_menu(YglMenu *m);
static void 	map_menu(YglMenu *m, YglMenu *parent, int item);
static int 	do_pup(YglMenu *m, int issub);
static Int32 	displayitem(YglMenu *m, int i, int on, int dosub);

static void create_menu(YglMenu *m) {
  if(menugc == NULL) {
    XColor fg, bg;
    menucursor = XCreateFontCursor(D, XC_right_ptr);
    bg.red = 0     ; bg.green = 0     ; bg.blue  = 0xffff;
    fg.red = 0xffff; fg.green = 0xffff; fg.blue  = 0     ;
    XRecolorCursor(D, menucursor, &fg, &bg);
    black = YGL_COLORS(BLACK);
    white = YGL_COLORS(WHITE);
  }
  {
    Ulong mask = 0;
    XSetWindowAttributes swa;
    
    mask |= CWBackPixmap;       swa.background_pixmap = None;
    mask |= CWSaveUnder;        swa.save_under        = True;
    mask |= CWEventMask;        swa.event_mask        = ExposureMask /*StructureNotifyMask*/;
    mask |= CWOverrideRedirect; swa.override_redirect = True;
    mask |= CWColormap;         swa.colormap          = Ygl.CCmap;
    mask |= CWCursor;           swa.cursor            = menucursor;
  
    m->win = XCreateWindow(D, RootWindow(D, YglScreen),
			   0, 0,
			   m->w + MENUSHADE,
			   m->h + MENUSHADE,
			   0,
			   Ygl.CV.depth,
			   InputOutput,
			   Ygl.CV.visual,
			   mask,
			   &swa);
  }
  if(menugc == NULL) {
    Ulong mask = 0;
    XGCValues gcv;
    Pixmap stipple;
    
    stipple = XCreatePixmapFromBitmapData(D, m->win,
					  gray1_bits, 
					  gray1_width, 
					  gray1_height,
					  black,
					  white,
					  1);
    
    mask |= GCGraphicsExposures; gcv.graphics_exposures = False;
    mask |= GCLineWidth;         gcv.line_width = 0;
    mask |= GCForeground;        gcv.background = black;
    mask |= GCBackground;        gcv.foreground = white;
    mask |= GCStipple;           gcv.stipple = stipple;
    
    menugc = XCreateGC(D, m->win, mask, &gcv);
  }
}

static YglMenu *mid_to_m(const char *caller, Int32 mid) {
  if(mid < 1 || mid > YglLastMenu || YglMenus[mid] == NULL) {
    Yprintf(caller, "invalid menu id: %d\n", mid);
    exit(1);
  }
  return YglMenus[mid];
}

static void add_to_pup(const char *caller, Int32 mid, const char *s, va_list argp) {
  YglMenu *m;
  int i;
  
  m = mid_to_m(caller, mid);
  
  for(i = 0; s[i];) {
    YglMenuItem *next;
    char buf[1024], *bufp = buf;
    int tflag = False;
    
    next = m->item[m->nitem];
    next->callF    = True;
    next->callback = NULL;
    next->submenu  = NULL;
    next->retval   = m->nitem;
    next->mode     = PUP_NONE;
    next->line     = False;
#ifdef DEBUG
    fprintf(stderr, "add_to_pup: caller = '%s', mid = %d, i = %d, string = '%s'\n",
	    caller, mid, i, &s[i]);
#endif
    for(;s[i] && s[i] != '|'; i++) {
      switch(s[i]) {
      case '%':
	i++;
	switch(s[i]) {
	case '%':
	  *bufp++ = s[i];
	  break;
	case 't':
	  tflag = True;
	  next = m->item[0];
	  break;
	case 'F':
	  /* m->item[0]->callback = (Int32(*)(Int32, ...)) va_arg(argp, void*);*/
	  m->item[0]->callback = va_arg(argp, YglMenuCallback);
#ifdef DEBUG
	  fprintf(stderr, "added %%F callback = 0x%x\n",
		  m->item[0]->callback);
#endif
	  break;
	case 'n':
	case 'm':
	  next->callF = False;
	  /* No break */
	case 'f':
	case 'M': /* %M is not available in SGI's gl. See examples/popup.c. */
	  if(tflag) {
	    Yprintf(caller, 
		    "invalid combination of '%%t' and '%%%c' in '%s', position %d\n",
		    s[i], s, i);
	    exit(1);
	  }
	  if(s[i] == 'm' || s[i] == 'M') {
	    Int32 mid = va_arg(argp, Int32);
	    next->submenu = mid_to_m(caller, mid);
#ifdef DEBUG
	    fprintf(stderr, "added %%m submenu = %d\n", mid);
#endif
	  } else { /* 'f' or 'n' */
	    next->callback = va_arg(argp, YglMenuCallback);
	    /* (Int32(*)(Int32, ...)) va_arg(argp, void*); */
#ifdef DEBUG
	    fprintf(stderr, "added %%fn callback = 0x%x\n", next->callback);
#endif
	  }
	  break;
	case 'x': {
	  char *endp;
	  i++;
	  next->retval = (Int32) strtol(&s[i], &endp, 0);
	  i = endp - s - 1;
	}
	  break;
	case 'l':
	  next->line = True;
	  break;
	default:
	  Yprintf(caller,
		  "unknown token '%%%c' in '%s', position %d\n", 
		  s[i], s, i);
	  /* exit(1); */
	  break;
	}
	break;
      default:
	*bufp++ = s[i];
	break;
      }
    } /* for(;...) */
    *bufp = '\0';
    
    if(NULL == (next->name = (char*)malloc(strlen(buf) + 1))) {
      Yprintf(caller, "can't allocate memory.\n");
      exit(-1);
    }
    
    strcpy(next->name, buf);
    
    if (strlen(next->name) > strlen(m->longest->name)) m->longest = next;
    
    if(!tflag) {
      m->nitem++;
      if(NULL == (m->item = (YglMenuItem**)realloc(m->item, (m->nitem + 1) * sizeof(YglMenuItem*))) ||
	 NULL == (m->item[m->nitem] = (YglMenuItem*)calloc(1, sizeof(YglMenuItem)))) {
	Yprintf(caller, "can't allocate memory.\n");
	exit(-1);
      }
    }
    
    if(s[i] == '|') i++;
  } /* for */
}

static Int32 displayitem(YglMenu *m, int i, int on, int dosub) {
  int fg = 0, bg = 0, r = -1;
  YglMenuItem *item;
  
  if(i == -1) return -1;
  
  item = m->item[i];
  
  switch(item->mode) {
  case PUP_GREY:
    XSetFillStyle(D, menugc, FillStippled);
    fg = bg = black;
    break;
  case PUP_NONE:
    if(on) {
      fg = white; bg = black;
    } else {
      bg = white; fg = black;
    }
    break;
  }
  
  XSetForeground(D, menugc, bg);
  XFillRectangle(D, m->win, menugc,
		 1, 
		 menuheight * i + 5,
		 m->w - 1, 
		 menuheight);
  
  XSetForeground(D, menugc, fg);
  XSetFillStyle(D, menugc, FillSolid);
  XDrawString(D, m->win, menugc,
	      MENUPADW, 
	      menuheight * i + (menuheight + fontascent)/2 + 3,
	      item->name, strlen(item->name));
  if(item->line) 
    XDrawLine(D, m->win, menugc, 
	      1, 
	      menuheight * (i + 1) + 4,
	      m->w - 1,
	      menuheight * (i + 1) + 4);
  
  if(item->submenu) {
    XPoint arrow[] = {{0, 0}, {-6, -3}, {0, 6}, {6, -3}};
    arrow[0].x = m->w - SUBMENUOFFSET + 6;
    arrow[0].y = menuheight * i + menuheight/2 + 4;
    XFillPolygon(D, m->win, menugc, arrow, 4, Convex, CoordModePrevious);
    
    if(dosub && item->mode != PUP_GREY) {
      if(on) {
	map_menu(item->submenu, m, i);
	r = do_pup(item->submenu, True);
      } else {
	XUnmapWindow(D, item->submenu->win);
      }
    }
  }
  return r;
}


static void map_menu(YglMenu *m, YglMenu *parent, int item) {
  int i;
  
  if(parent == NULL) { /* Mainmenu */
    m->x = Ygl.lastreadevent.xbutton.x_root;
    m->y = Ygl.lastreadevent.xbutton.y_root - menuheight;
  } else { /* Submenu */
    m->x = parent->x + parent->w - SUBMENUOFFSET;
    m->y = parent->y + item * menuheight + 4;
  }
  
  m->w = strwidth(m->longest->name) + 2 * MENUPADW;
  m->h = m->nitem * menuheight + 5;
  
#ifdef DEBUG
  fprintf(stderr, "map_menu: width = %d height = %d longest = %s\n",
	  m->w, m->h, m->longest->name);
#endif
    
  XMoveResizeWindow(D, 
		    m->win,
		    MIN(m->x, YglScreenWidth  - (m->w + MENUSHADE)) - 2,
		    MIN(m->y, YglScreenHeight - (m->h + MENUSHADE)) - 2,
		    m->w + MENUSHADE,
		    m->h + MENUSHADE);
  
  XMapRaised(D, m->win);
  XSync(D, False);
  
  /* SaveUnder doesn't work properly when we only wait for MapNotify. Strange */
  /*Ygl.await_windowevent(m->win, StructureNotifyMask, MapNotify);*/
  Ygl.await_windowevent(m->win, ExposureMask, Expose);
  
  XSetForeground(D, menugc, white);
  XFillRectangle(D, m->win, menugc, 0, 0, m->w, m->h);
  
  XSetForeground(D, menugc, black);
  XDrawRectangle(D, m->win, menugc, 0, 0, m->w, m->h);
  
  XSetFillStyle(D, menugc, FillStippled);
  XFillRectangle(D, m->win, menugc, 
		 m->w, MENUSHADE, MENUSHADE, m->h);
  XFillRectangle(D, m->win, menugc, 
		 MENUSHADE, m->h, m->w - MENUSHADE, MENUSHADE);
  XSetFillStyle(D, menugc, FillSolid);
  
  XDrawString(D, m->win, menugc,
	      (m->w - strwidth(m->item[0]->name)) / 2, 
	      (menuheight + fontascent)/2 + 2,
	      m->item[0]->name, strlen(m->item[0]->name));
  
  XDrawRectangle(D, m->win, menugc, 
		 2, 2,
		 m->w - 4,
		 menuheight);
  
  XDrawLine(D, m->win, menugc, 
	    0,    menuheight + 4,
	    m->w, menuheight + 4);
  
  for (i = 1; i < m->nitem; i++) {
    displayitem(m, i, OFF, 0);
  }
}

static int do_pup(YglMenu *m, int issub) {
  int menuval = -1, oldmenuval = -1, r, subr = -1, olddosub = -1;
  Uint mask;
  YglMenuItem *item = NULL;
  
  do {
    Window root, win;
    int rx, ry, wx, wy, dosub;
    subr = -1;
#if DEBUG > 1
    fprintf(stderr, "XQ: win = 0x%x\n", m->win);
#endif
    XQueryPointer(D, m->win, &root, &win, &rx, &ry, &wx, &wy, &mask);
    
    dosub = wx > m->w - SUBMENUOFFSET;
    
    if(issub && wx < 0) { /* Close submenus on left leave */
      menuval = -1;
      break;
    }
    
    if(wx > 0 &&
       wx < m->w &&
       wy > menuheight + 4 &&
       wy < menuheight * m->nitem + 4) {
      menuval = (wy - 4) / menuheight;
      item    = m->item[menuval];
      
      if(oldmenuval != menuval || olddosub != dosub) { /* Something changed */
#ifdef DEBUG
	fprintf(stderr, "menuval = %d old = %d\n", 
		menuval, oldmenuval);
#endif
	/*  */ displayitem(m, oldmenuval, OFF, dosub);
	subr = displayitem(m,    menuval,  ON, dosub);
	oldmenuval = menuval;
	olddosub   = dosub;
      }
    } else if(oldmenuval != -1) {
      displayitem(m, oldmenuval, OFF, dosub);
      oldmenuval = menuval = -1;
    }
  } while(subr == -1 && mask & Button3Mask);
  
#ifdef DEBUG
  if(subr >= 0)
    fprintf(stderr, "while done, subr = %d, menuval = %d, retval = %d, "
	    "item[0]->callback = 0x%x, callback = 0x%x\n", 
	    subr, menuval, item->retval,
	    m->item[0]->callback, item->callback);
#endif
  
  XUnmapWindow(D, m->win);
  
  if(subr != -1) {
    if(item->callF && m->item[0]->callback) {
#ifdef DEBUG
      Int32 prev = subr;
#endif
      subr = m->item[0]->callback(menuval, subr);
#ifdef DEBUG
      fprintf(stderr, "do_pupx: %%F(%d, %d) = %d\n", 
	      menuval, prev, subr);
#endif
    }
    return subr;
  }
  
  if(menuval != -1 && item->submenu == NULL && item->mode != PUP_GREY) {
    /* Valid item selected */
    r = item->retval;
    if(item->callback) {
#ifdef DEBUG
      Int32 prev = r;
#endif
      r = item->callback(r);
#ifdef DEBUG
      fprintf(stderr, "do_pup: %%f(%d) = %d\n", prev, r);
#endif
    }
    if(item->callF && m->item[0]->callback) {
#ifdef DEBUG
      Int32 prev = r;
#endif
      r = m->item[0]->callback(r);
#ifdef DEBUG
      fprintf(stderr, "do_pup: %%F(%d) = %d\n", prev, r);
#endif
    }
  } else {
    r = -1;
  }      
  return r;
}


Int32 newpup(void) {
  YglMenu *m;
  int mid;
  const char * MyName = "newpup";
  I(MyName, "");
  
  if(YglMenus == NULL) { /* initialize */
    mid = YglLastMenu = 1;
    YglMenus = (YglMenu**)malloc(2 * sizeof(YglMenu*)); /* First el. never used */
  } else {
    for(mid = 1; mid <= YglLastMenu && YglMenus[mid] != NULL; mid++);
    if(mid > YglLastMenu) {
      YglLastMenu++;
      YglMenus = (YglMenu**)realloc(YglMenus, (YglLastMenu + 1) * sizeof(YglMenu*));
    }
  }
  if(YglMenus == NULL ||
     NULL == (m = YglMenus[mid] = (YglMenu*)calloc(1, sizeof(YglMenu)))) {
    Yprintf(MyName, "can't allocate memory.\n");
    exit(-1);
  }
  if(NULL == (m->item = (YglMenuItem**)malloc(2 * sizeof(YglMenuItem*))) ||
     NULL == (m->item[0] = (YglMenuItem*)calloc(1, sizeof(YglMenuItem))) ||
     NULL == (m->item[1] = (YglMenuItem*)calloc(1, sizeof(YglMenuItem))) ||
     NULL == (m->item[0]->name = malloc(16))) {
    Yprintf(MyName, "can't allocate memory.\n");
    exit(-1);
  }
  
  sprintf(m->item[0]->name, "Menu %d", mid);
  m->item[0]->callback = NULL;
  m->longest	= m->item[0];
  m->nitem	= 1;
  m->w 		= 10;
  m->h 		= 10;
  
  create_menu(m);
#ifdef DEBUG
  fprintf(stderr, "newpup: ready, mid = %d, m = 0x%x\n", mid, m);
#endif
  return mid;
}

Int32 defpup(Char8 *string, ...) {
  Int32 mid;
  va_list Argp;
  const char * MyName = "defpup";
  I(MyName, "'%s',...", string);
  mid = newpup();
  va_start(Argp, string);
  add_to_pup(MyName, mid, string, Argp);
  va_end(Argp);
  return mid;
}

void addtopup(Int32 mid, Char8 *string, ...) {
  va_list Argp;
  const char * MyName = "addtopup";
  I(MyName, "%d,'%s',...", mid, string);
  va_start(Argp, string);
  add_to_pup(MyName, mid, string, Argp);
  va_end(Argp);
}

void freepup(Int32 mid) {
  YglMenu *m;
  int i;
  const char * MyName = "freepup";
  I(MyName, "%d", mid);
  m = mid_to_m(MyName, mid);
  XDestroyWindow(D, m->win);
  for(i = 0; i <= m->nitem; i++) {
    if(m->item[i]->name) free(m->item[i]->name);
    free(m->item[i]);
  }
  free(m->item);
  free(m);
  YglMenus[mid] = NULL;
}

void setpup(Int32 mid, Int32 entry, Int32 mode) {
  YglMenu *m;
  const char * MyName = "setpup";
  I(MyName, "%d,%d,%d", mid, entry, mode);
  m = mid_to_m(MyName, mid);
  if(entry < 1 || entry > m->nitem - 1) {
    Yprintf(MyName, "entry %d in menu %d does not exist.\n",
	    entry, mid);
    return;
  }
  switch(mode) {
  case PUP_NONE:
  case PUP_GREY:
    m->item[entry]->mode = mode;
    break;
  default:
    Yprintf(MyName, "invalid mode: %d.\n", mode);
    return;
  }
}

Int32 dopup(Int32 mid) {
  YglMenu *m;
  Int32 r;
  const char * MyName = "dopup";
  Colormap *cmaps;
  int i, numcmaps, install = 1;
  /*Window CMapWindows[3];*/
  
  I(MyName, "%d", mid);
  m = mid_to_m(MyName, mid);
  if(Ygl.lastreadevent.type != ButtonPress) {
    Yprintf(MyName, "last event was not ButtonPress.\n");
    return -1;
  }
  XSetFont(D, menugc, Ygl.Fonts[W->font].fs->fid);
  fontascent = Ygl.Fonts[W->font].fs->ascent;
  menuheight = MENUPADH + Ygl.Fonts[W->font].fs->descent + fontascent;
  
  cmaps = XListInstalledColormaps(D, m->win, &numcmaps);
  
  /* Do we have to install menu cmap? */
  for(i = 0; i < numcmaps && (install &= cmaps[i] != Ygl.CCmap); i++);
  
#if 0
  /* Any better idea how to get the cmap installed when popping 
     up the menu in RGB mode??? */
  
  CMapWindows[0] = m->win;
  CMapWindows[1] = W->win;
  CMapWindows[2] = W->top;
  
  if(XSetWMColormapWindows(D, W->top, CMapWindows, 3) == 0) {
    Yprintf(MyName, "cannot change WM_COLORMAP_WINDOWS property.\n");
  }
  
  /*XUngrabPointer(D, 0);*/
  
#endif
  if(install) {
#ifdef DEBUG
    fprintf(stderr, "dopup: installing colormap\n");
#endif
    XInstallColormap(D, Ygl.CCmap); /* As we may be in RGBmode */
  }
  
  /* Switch cursor on */
  XChangeActivePointerGrab(D, 0, menucursor, CurrentTime);
  
  map_menu(m, NULL, 0);
  r = do_pup(m, False);
  
  /* Restore installed colormaps */
  if(install) for(i = 0; i < numcmaps; i++) XInstallColormap(D, cmaps[i]);
  XFree((char*) cmaps);
  
  return r;
}
