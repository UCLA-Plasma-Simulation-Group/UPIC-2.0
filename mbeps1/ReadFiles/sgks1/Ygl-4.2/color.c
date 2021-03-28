/*
 *    Ygl: Run GL programs with standard X11 and/or OpenGL routines.
 *    (C) Fred Hucht 1993-2007
 *    EMail: fred<at>thp.Uni-Duisburg.de
 */

static const char vcid[] = "$Id: color.c,v 4.9 2007-05-08 11:01:05+02 fred Exp $";

#include "header.h"

static int   rect_read_errh(Display *, XErrorEvent *);
static Int32 rect_read(const char*, Screencoord, Screencoord, Screencoord, Screencoord, int, void *);
static void rect_write(const char*, Screencoord, Screencoord, Screencoord, Screencoord, int, void *);
static Ulong emuColorsInv(Ulong x);
#ifdef X11
static void create_dither(int r, int g, int b);
static void set_color(const char *caller, int r, int g, int b);
#endif

#define YGL_COLORSINV(x) (Ygl.PCM ? CMapSize-1-(x) : Ygl.EmulateCmap ? emuColorsInv(x) : Ygl.ColorsInv[x])

static Ulong emuColorsInv(Ulong x) {
  Ulong i;
  for(i = 0; i < CMapSize && Ygl.Colors[i] != x; i++);
  return i < CMapSize ? i : 0;
}

void mapcolor(Colorindex ind, Int16 r, Int16 g, Int16 b) {
  XColor xc;
  const char * MyName = "mapcolor";
  I(MyName, "%d,%d,%d,%d", ind, r, g, b);
  if(W->rgb) { Yprintf(MyName, "not in CMap mode.\n"); return; }
  if(ind >= CMapSize) {
    Yprintf(MyName, "can't map cell %d (must be between 0 and %d).\n",
	    ind, CMapSize-1);
    return;
  }
  xc.red   = r << 8;
  xc.green = g << 8;
  xc.blue  = b << 8;
  if(!Ygl.PCM) { /* use default colormap */
    if (Ygl.Colors[ind]) XFreeColors(D, Ygl.CCmap, Ygl.Colors + ind, 1, 0);
    
    if(XAllocColor(D, Ygl.CCmap, &xc)) {
      Ygl.Colors[ind] = xc.pixel;
      if(!Ygl.EmulateCmap)
	Ygl.ColorsInv[xc.pixel] = ind;
    } else {
      Yprintf(MyName,
	      "can't allocate color (%d, %d, %d) for cell %d, colormap full.\n",
	      r, g, b, ind);
    }
  }
  else { /* use private colormap */
    xc.pixel = CMapSize - 1 - ind;
    xc.flags = DoRGB;
    XStoreColor(D, Ygl.CCmap, &xc);
  }
#ifdef DEBUG
  fprintf(stderr, "mapcolor: index = %d, xc = %d (%d,%d,%d)\n",
	  ind, xc.pixel, xc.red, xc.green, xc.blue);
#endif
  if(!Ygl.GC) { /* if GC list in CMap windows */
    int i;
    YglWindow *w = &Ygl.Windows[1]; /* first window */
    for(i = 1; i < Ygl.NextWindow; i++, w++) {
      if(w->main != 0 && !w->rgb) {
	XSetForeground(D, w->gclist[ind], ind);
      }
    }
  }
  F;
}

#ifdef X11
static const char *DI[] = { 	/* index vector for dither matrices */
  "\000", 			/* DI[sz][c] >> 3 : x-coordinate */
  "\000", 			/* DI[sz][c] &  7 : y-coordinate */
  				/* 2x2  ***+-0--1-*/
  "\000\011"				/*0| 0  3 */
  "\001\010", 				/*1| 2  1 */
    				/* 3x3  ***+-0--1--2-*/
  "\000\011\012"			/*0| 0  5  3 */
  "\020\001\010"			/*1| 4  1  7 */
  "\022\021\002", 			/*2| 8  2  6 */
    				/* 4x4  ***+-0--1--2--3-*/
  "\000\022\002\020"			/*0| 0  8  3 11 */
  "\011\033\013\031"			/*1|12  4 15  7 */
  "\010\032\012\030"			/*2| 2 10  1  9 */
  "\021\003\023\001", 			/*3|14  6 13  5 */
    				/* 5x5  ***+-0--1--2--3--4-*/
  "\000\012\024\031\043"		/*0| 0  5 10 15 20 */
  "\010\022\034\041\003"		/*1|13 18 23  3  8 */
  "\020\032\044\001\013"		/*2|21  1  6 11 16 */
  "\030\042\004\011\023"		/*3| 9 14 19 24  4 */
  "\040\002\014\021\033", 		/*4|17 22  2  7 12 */
				/* 6x6  ***+-0--1--2--3--4--5-*/
  "\000\033\002\035\004\031" 		/*0| 0 18 16 34  9 27 */
  "\011\044\013\040\015\042" 		/*1|33  6 24  5 23 15 */
  "\022\055\024\051\020\053" 		/*2| 2 20 12 30 11 29 */
  "\010\043\012\045\014\041" 		/*3|35  8 26  1 19 17 */
  "\021\054\023\050\025\052" 		/*4| 4 22 14 32  7 25 */
  "\032\005\034\001\030\003" 		/*5|31 10 28  3 21 13 */
};

static void create_dither(int r, int g, int b) {
  YglWindow *w = W;
  int i, c;
  int sz  = Ygl.DSZ, szq = sz * sz;
  int ri, gi, bi, rx, gx, bx, ry, gy, by, red, green, blue;
  static int init = 0;
  static int rn, gn, bn;
  
  if(init == 0) {		/* Initialize static vars */
    init = 1;
    rn = ((256 >> Ygl.rb) - 1) * szq + 1; /* Number of dither patterns red   */
    gn = ((256 >> Ygl.gb) - 1) * szq + 1; /* Number of dither patterns green */
    bn = ((256 >> Ygl.bb) - 1) * szq + 1; /* Number of dither patterns blue  */
#ifdef DEBUG
    fprintf(stderr, "create_dither: sz = %d, {rn, gn, bn} = {%d,%d,%d}\n",
	    sz, rn, gn, bn);
#endif
  }
  
  /* Index of dither pat., base color, # elements with color rx + 1 */
  ri = (r * rn) >> 8; rx = ri / szq; ry = ri % szq;  
  gi = (g * gn) >> 8; gx = gi / szq; gy = gi % szq;
  bi = (b * bn) >> 8; bx = bi / szq; by = bi % szq;
  
  red   = (ri << 8) / rn;
  green = (gi << 8) / gn;
  blue  = (bi << 8) / bn;
  
  if(red != w->red || green != w->green || blue != w->blue) { /*Already set?*/
    /* No, we must change pattern/color */
    Ulong mask;
    XGCValues values;
    
    w->red = red; w->green = green; w->blue = blue;
    if(ry + gy + by == 0) {
      /* We hit pure color, don't dither but set foreground in w->gc */
      mask = GCFillStyle | GCForeground;
      values.fill_style = FillSolid;
      values.foreground =
	((rx << Ygl.rs) & Ygl.RMask) +
	((gx << Ygl.gs) & Ygl.GMask) +
	((bx << Ygl.bs) & Ygl.BMask);
    } else {
      /* Dither. */
      mask = GCTile | GCFillStyle;
      values.fill_style = FillTiled;
      values.tile       = w->pm;
      
      for(c = szq; c > 0; c--) {
	if(c == ry) rx++;	/* Go to next color */
	if(c == gy) gx++;
	if(c == by) bx++;
	
	i = DI[sz][c-1];	/* Next position in dither matrix */
	
	XPutPixel(w->pmi, i >> 3, i & 7, /* Simplified RGB_TO_RGBVisual */
		  ((rx << Ygl.rs) & Ygl.RMask) +
		  ((gx << Ygl.gs) & Ygl.GMask) +
		  ((bx << Ygl.bs) & Ygl.BMask));
      }
      
      XPutImage(D, w->pm, w->pmgc, w->pmi, 0, 0, 0, 0, sz, sz);
    }
    XChangeGC(D, w->gc, mask, &values);
  }
}
#endif /* X11 */
/* converts 24 bit RGB values ( 0xBBGGRR ) to value appropiate for RGBVisual */
#define RGB24_TO_RGBVisual(rgb24) ( \
	((((rgb24) >> (Ygl.rb +  0)) << Ygl.rs) & Ygl.RMask) \
      + ((((rgb24) >> (Ygl.gb +  8)) << Ygl.gs) & Ygl.GMask) \
      + ((((rgb24) >> (Ygl.bb + 16)) << Ygl.bs) & Ygl.BMask) )

#define RGBVisual_TO_RGB24(pixel) ( \
	((((pixel) & Ygl.RMask) >> Ygl.rs) << (Ygl.rb +  0)) \
      + ((((pixel) & Ygl.GMask) >> Ygl.gs) << (Ygl.gb +  8)) \
      + ((((pixel) & Ygl.BMask) >> Ygl.bs) << (Ygl.bb + 16)) )

#define RGB_TO_RGBVisual(r,g,b) ( \
	((((r) >> Ygl.rb) << Ygl.rs) & Ygl.RMask) \
      + ((((g) >> Ygl.gb) << Ygl.gs) & Ygl.GMask) \
      + ((((b) >> Ygl.bb) << Ygl.bs) & Ygl.BMask) )
#ifdef X11
static void set_color(const char *caller, int r, int g, int b) {
  if(!W->rgb) { Yprintf(caller, "not in RGB mode.\n"); return; }
  if(W->pm) {
    create_dither(r, g, b);
  } else {
    XSetForeground(D, W->gc, RGB_TO_RGBVisual(r, g, b));
  }
}
#endif /* X11 */

void RGBcolor(Int16 r, Int16 g, Int16 b) {
  const char * MyName = "RGBcolor";
  I(MyName, "%d,%d,%d", r, g, b);
  IFOGL(glColor3ub(r, g, b),
	set_color(MyName, r, g, b)
	);
}

void cpack(Uint32 rgb) {
  const char * MyName = "cpack";
  I(MyName, "%0x%x", rgb);
  IFOGL(glColor4ubv((GLubyte *)&rgb),
	set_color(MyName,
		  (rgb >>  0) & 0xFF,
		  (rgb >>  8) & 0xFF,
		  (rgb >> 16) & 0xFF)
	);
}

void c3s(Int16 c[3]) { IFOGL(glColor3usv((GLushort *)c),set_color("c3s", c[0], c[1], c[2])); }
void c4s(Int16 c[4]) { IFOGL(glColor4usv((GLushort *)c),set_color("c4s", c[0], c[1], c[2])); }

void c3i(Int32 c[3]) { IFOGL(glColor3uiv((GLuint *)c),set_color("c3i", c[0], c[1], c[2])); }
void c4i(Int32 c[4]) { IFOGL(glColor4uiv((GLuint *)c),set_color("c4i", c[0], c[1], c[2])); }

void c3f(Float32 c[3]) { IFOGL(glColor3fv(c), float f = 255.9999999; set_color("c3f", (int)(f * c[0]), (int)(f * c[1]), (int)(f * c[2]))); }
void c4f(Float32 c[4]) { IFOGL(glColor4fv(c), float f = 255.9999999; set_color("c4f", (int)(f * c[0]), (int)(f * c[1]), (int)(f * c[2]))); }

#if 0 /* This version is veeeryyy slooowwww !!! */
void RGBcolor (Int16 r, Int16 g, Int16 b) {
  XColor xc;
  const char * MyName = "RGBcolor";
  I(MyName);
  if(!W->rgb) { Yprintf(MyName, "not in RGB mode.\n"); return; }
  xc.red   = r << 8;
  xc.green = g << 8;
  xc.blue  = b << 8;
  if(XAllocColor(D, Ygl.RCmap, &xc))
    XSetForeground(D, W->gc, xc.pixel);
  else
    Yprintf(MyName, "can't allocate color (%d, %d, %d), colormap full.\n",
	    r, g, b);
}
#endif

Int32 getcolor(void) {
  const char * MyName = "getcolor";
  I(MyName, "");
  if(W->rgb) { Yprintf(MyName, "not in CMap mode.\n"); return 0; }
  return W->color;
}

void getmcolor (Colorindex ind, Int16 *r, Int16 *g, Int16 *b) {
  XColor xc;
  const char * MyName = "getmcolor";
  I(MyName, "%d,*,*,*", ind);
  if(W->rgb) { Yprintf(MyName, "not in CMap mode.\n"); return; }
  xc.pixel = YGL_COLORS(ind);
  XQueryColor(D, Ygl.CCmap, &xc);
  *r = xc.red   >> 8;
  *g = xc.green >> 8;
  *b = xc.blue  >> 8;
}

void getmcolors (Colorindex ind1, Colorindex ind2,
		 Int16 *r, Int16 *g, Int16 *b) {
  XColor *xc;
  Colorindex i;
  int n = ind2 - ind1 + 1;
  const char * MyName = "getmcolors";
  I(MyName, "%d,%d,*,*,*", ind1, ind2);
  if(W->rgb) {
    Yprintf(MyName, "not in CMap mode.\n");
    return;
  }

  if(n < 0) {
    Yprintf(MyName, "2nd argument < 1st argument.\n");
    return;
  }
  
  if(NULL == (xc = (XColor*) malloc(n * sizeof(XColor)))) {
    Yprintf(MyName, "can't allocate memory.\n");
    exit(1);
  }

  for(i = 0; i < n; i++) xc[i].pixel = YGL_COLORS(i);
  XQueryColors(D, Ygl.CCmap, xc, n);
  for(i = 0; i < n; i++) {
    r[i] = xc[i].red   >> 8;
    g[i] = xc[i].green >> 8;
    b[i] = xc[i].blue  >> 8;
  }
  free(xc);
}

void gRGBcolor(Int16 *r, Int16 *g, Int16 *b) {
  XGCValues ret;
  const char * MyName = "gRGBcolor";
  I(MyName, "*,*,*");
  if(!W->rgb) {
    Yprintf(MyName, "not in RGB mode.\n");
    return;
  }
  if(W->pm) {
    *r = W->red;
    *g = W->green;
    *b = W->blue;
  } else {
    if(XGetGCValues(D, W->gc, GCForeground, &ret)) {
      *r = ((ret.foreground & Ygl.RMask) >> Ygl.rs) << Ygl.rb;
      *g = ((ret.foreground & Ygl.GMask) >> Ygl.gs) << Ygl.gb;
      *b = ((ret.foreground & Ygl.BMask) >> Ygl.bs) << Ygl.bb;
    } else {
      Yprintf(MyName, "cannot get GCValues.\n");
    }
  }
#ifdef DEBUG
  fprintf(stderr, "gRGBcolor: foreground = %d, (r,g,b) = (%d,%d,%d).\n",
	  ret.foreground,*r,*g,*b);
#endif
}

void color(Colorindex ind) {
  const char * MyName = "color";
  I(MyName, "%d", ind);
  if(W->rgb) {
    Yprintf(MyName, "not in CMap mode.\n");
  } else if(ind >= CMapSize) {
    Yprintf(MyName, "invalid color %d (must be between 0 and %d).\n",
	    ind, CMapSize-1);
  }
  else {
    Ulong col = YGL_COLORS(ind);
    W->color = ind;
    IFOGL(
	  if(Ygl.EmulateCmap) {
	    glColor3ub((col >> 16) & 0xff, (col >> 8) & 0xff, col & 0xff);
	  } else {
	    /*fprintf(stderr, "glIndexi(0x%x)\n", col);*/
	    glIndexi(col);
	  }
	  ,
	  if(Ygl.GC) XSetForeground(D, W->gc, col); /* single GC */
	  else       W->gc = W->gclist[col]
	  );
  }
}

#ifndef YGL_PREFIX
void Xcolor(unsigned int ind) {
  const char * MyName = "Xcolor";
  I(MyName, "%d", ind);
  if(W->rgb) {
    Yprintf(MyName, "not in CMap mode.\n");
  } /*else if(ind >= CMapSize)
      Yprintf(MyName, "invalid color %d (must be between 0 and %d).\n",
      ind, CMapSize-1);*/
  else {
    W->color = YGL_COLORSINV(ind);
    IFOGL(/*fprintf(stderr, "glIndexi(0x%x)\n", ind);*/
	  if(Ygl.EmulateCmap) {
	    glColor3ub((ind >> 16) & 0xff, (ind >> 8) & 0xff, ind & 0xff);
	  } else {
	    /*fprintf(stderr, "glIndexi(0x%x)\n", ind);*/
	    glIndexi(ind);
	  }
	  ,
	  if(Ygl.GC) XSetForeground(D, W->gc, ind); /* single GC */
	  else       W->gc = W->gclist[ind]
	  );
  }
}

void getmXcolor (Colorindex ind, Int16 *r, Int16 *g, Int16 *b) {
  XColor xc;
  const char * MyName = "getmXcolor";
  I(MyName, "%d,*,*,*", ind);
  if(W->rgb)
    Yprintf(MyName, "not in CMap mode.\n");
  else {
    xc.pixel = ind;
    XQueryColor(D, Ygl.CCmap, &xc);
    *r = xc.red   >> 8;
    *g = xc.green >> 8;
    *b = xc.blue  >> 8;
  }
}
#endif

/* block pixel transfer routines */
static int rect_read_front = True;

void readsource(Int32 source) {
  const char * MyName = "readsource";
  I(MyName, "%d", source);
  switch(source) {
  case SRC_AUTO: 
  case SRC_BACK: 
    rect_read_front = False;
    break;
  case SRC_FRONT: 
    rect_read_front = True;
    break;
  default:
    Yprintf(MyName, "unknown mode %d.\n", source);
    break;
  }
}

static int rect_read_err;
static int rect_read_errh(Display *dpy, XErrorEvent *error) {
  rect_read_err = error->error_code;
  return 0;
}

static Float32 Zoom[2] = {1.0, 1.0};
static int DoZoom = False;

void rectzoom(Float32 xfactor, Float32 yfactor) {
  const char * MyName = "readsource";
  I(MyName, "%g,%g", xfactor, yfactor);
  if(xfactor < 0.0 || yfactor < 0.0) {
    Yprintf(MyName, "xfactor of yfactor negative.\n");
    return;
  }
  Zoom[0] = xfactor;
  Zoom[1] = yfactor;
  DoZoom = xfactor != 1.0 || yfactor != 1.0;
}

static Int32 rect_read(const char *caller,
		       Screencoord x1, Screencoord y1,
		       Screencoord x2, Screencoord y2,
		       int size, void *data) {
  XImage *XI;
  long x, y;
  Int32 *data32 = (Int32*) data;
  Int16 *data16 = (Int16*) data;
  Uint8 *data8  = (Uint8*) data;
  Screencoord width  = x2 - x1 + 1;
  Screencoord height = y2 - y1 + 1;
  int (*old_handler)(Display*, XErrorEvent*);
  
  rect_read_err = 0; /* reset */
  old_handler = XSetErrorHandler(rect_read_errh);
  
  XI = XGetImage(D, 
		 rect_read_front ? IF_RGBWIN(W->win, W->main) : W->draw,
		 x1, W->ym - 1 - y2,
		 width, height,
		 (1 << YglDepth())-1,
		 ZPixmap);
  
  XSync (D, False);
  XSetErrorHandler(old_handler); /* reset */
  
#ifdef DEBUG
  fprintf(stderr, "rect_read: (w, h) = (%d, %d), XI->(w, h) = (%d, %d), err = %d\n",
	  width, height, XI->width, XI->height, rect_read_err);
#endif
  if(rect_read_err == BadMatch) {
    Yprintf(caller, "window not completely visible.\n");
    return 0;
  }
  
  for(y = XI->height-1; y >= 0; y--) {
#if DEBUG > 1
    fprintf(stderr, "y = %d\n", y);
#endif
    for(x = 0; x < XI->width; x++) {
      Ulong pixel;
      long val;
      
      pixel = XGetPixel(XI, x, y);
      
      if(W->rgb) {
	val = RGBVisual_TO_RGB24(pixel);
      } else {
	val = YGL_COLORSINV(pixel);
      }
      
      switch(size) {
      case 1: *(data8++)  = val; break;
      case 2: *(data16++) = val; break;
      case 4: *(data32++) = val; break;
      }
#if DEBUG > 2
      fprintf(stderr, "%x,%x ", pixel, val);
#endif
    }
#if DEBUG > 2
    fputc(13, stderr);
#endif
  }
  x = XI->width * XI->height;
#ifdef DEBUG
  fprintf(stderr, "loop end.");
#endif
  XDestroyImage(XI);
  F;
  return x;
}

static void rect_write(const char *caller,
		       Screencoord x1, Screencoord y1,
		       Screencoord x2, Screencoord y2,
		       int size, void *data) {
  XImage *XI;
  long x, y, ys;
  Int32 *data32 = (Int32*) data;
  Int16 *data16 = (Int16*) data;
  Uint8 *data8  = (Uint8*) data;
  Screencoord width  = x2 - x1 + 1;
  Screencoord height = y2 - y1 + 1;
  Ulong xwidth  = Zoom[0] * width;
  Ulong xheight = Zoom[1] * height;
  float xf = 1.0 / Zoom[0];
  float yf = 1.0 / Zoom[1];
  
  XI = XCreateImage(D,
		    YglVisual(),
		    YglDepth(),
		    ZPixmap,	/* Format */
		    0,		/* Offset */
		    NULL,	/* Data */
		    xwidth,	/* Width */
		    xheight,	/* Height */
		    8, 		/* BitmapPad */
		    0);		/* BytesPerLine */

  if(NULL == (XI->data = (char*) malloc(XI->bytes_per_line * xheight))) {
    Yprintf(caller, "can't allocate memory.\n");
    exit(1);
  }
  
  for(ys = 0, y = xheight - 1; y >= 0; ys++, y--) for(x = 0; x < xwidth; x++) {
    Ulong idx, pixel;
    long val = -1;
    switch(size) {
#define ZOOMDATA(p, x_, y_) (p)[idx = (long)(xf * (x_)) + width * (long)(yf * (y_))]
    case 1: val = DoZoom ? ZOOMDATA(data8,  x, ys) : *(data8++);  break;
    case 2: val = DoZoom ? ZOOMDATA(data16, x, ys) : *(data16++); break;
    case 4: val = DoZoom ? ZOOMDATA(data32, x, ys) : *(data32++); break;
#undef ZOOMDATA
    }
    if(W->rgb) {
      pixel = RGB24_TO_RGBVisual(val);
    } else { /* CMap */
      if(val < 0 || val >= CMapSize) {
	Yprintf(caller, "invalid color %d at position (%d,%d) in CMap mode "
		"(must be between 0 and %d).\n", val, x, y, CMapSize-1);
	val = 0;
      }
      pixel = YGL_COLORS(val);
    }
#if USE_XPUTPIXEL
    XPutPixel(XI, x, y, pixel);
#else
    {
      void * addr;
      addr = XI->data + y * XI->bytes_per_line + XI->bits_per_pixel/8 * x;
      switch(XI->bits_per_pixel) {
      case  8: *(Int8 *)addr = pixel; break;
      case 16: *(Int16*)addr = pixel; break;
      case 32: *(Int32*)addr = pixel; break;
      }
    }
#endif
  }
  XPutImage(D, W->draw, W->gc, XI, 0, 0,
	    x1, W->ym - (y1 + xheight),
	    xwidth, xheight);
  
  /* free(XI->data) done by XDestroyImage() */
  XDestroyImage(XI);
  F;
}

Int32 lrectread(Screencoord x1, Screencoord y1,	Screencoord x2, Screencoord y2,	Int32 *data) { const char * MyName = "lrectread"; I(MyName, "%d,%d,%d,%d,*",x1,y1,x2,y2); return rect_read(MyName, x1,y1,x2,y2, 4, data);}
Int32  rectread(Screencoord x1, Screencoord y1,	Screencoord x2, Screencoord y2,	Int16 *data) { const char * MyName = "rectread" ; I(MyName, "%d,%d,%d,%d,*",x1,y1,x2,y2); return rect_read(MyName, x1,y1,x2,y2, 2, data);}
Int32 crectread(Screencoord x1, Screencoord y1,	Screencoord x2, Screencoord y2,	Uint8 *data) { const char * MyName = "crectread"; I(MyName, "%d,%d,%d,%d,*",x1,y1,x2,y2); return rect_read(MyName, x1,y1,x2,y2, 1, data);}

void lrectwrite(Screencoord x1, Screencoord y1, Screencoord x2, Screencoord y2, Int32 *data) { const char * MyName = "lrectwrite"; I(MyName, "%d,%d,%d,%d,*",x1,y1,x2,y2); rect_write(MyName, x1,y1,x2,y2, 4, data);}
void  rectwrite(Screencoord x1, Screencoord y1,	Screencoord x2, Screencoord y2,	Int16 *data) { const char * MyName = "rectwrite" ; I(MyName, "%d,%d,%d,%d,*",x1,y1,x2,y2); rect_write(MyName, x1,y1,x2,y2, 2, data);}
void crectwrite(Screencoord x1, Screencoord y1,	Screencoord x2, Screencoord y2,	Uint8 *data) { const char * MyName = "crectwrite"; I(MyName, "%d,%d,%d,%d,*",x1,y1,x2,y2); rect_write(MyName, x1,y1,x2,y2, 1, data);}

void rectcopy(Screencoord x1, Screencoord y1,
	      Screencoord x2, Screencoord y2,
	      Screencoord xn, Screencoord yn) {
  const char * MyName = "rectcopy";
  I(MyName, "%d,%d,%d,%d,%d,%d", x1, y1, x2, y2, xn, yn);
  if(DoZoom) {
    int n = (x2 - x1 + 1) * (y2 - y1 + 1);
    Int32 *buf;
    if(NULL == (buf = (Int32*) calloc(n, sizeof (Int32)))) {
      Yprintf(MyName, "can't allocate memory.\n");
      exit(1);
    }
    if(rect_read(MyName, x1, y1, x2, y2, 4, buf))
      rect_write(MyName, xn, yn, xn+x2-x1, yn+y2-y1, 4, buf);
    free((char *)buf);
  } else {
    XCopyArea(D, W->draw, W->draw, W->gc,
	      x1, W->ym - 1 - y2,
	      x2 - x1 + 1,
	      y2 - y1 + 1,
	      xn, W->ym - 1 - yn - (y2-y1));
  }
  F;
}

Int32 readpixels(Int16 n, Colorindex data[]) {
  Screencoord x, y;
  Int32 r;
  const char * MyName = "readpixels";
  I(MyName, "%d,*", n);
  if(W->rgb) {
    Yprintf(MyName, "not in CMap mode.\n");
    exit(1);
  }
  x = X(W->xc);
  y = Y(W->yc);
  r = rect_read(MyName, x, y, x + n, y + 1, 2, data);
  W->xc += n / W->xf;
  return r;
}

void writepixels(Int16 n, Colorindex data[]) {
  Screencoord x, y;
  const char * MyName = "writepixels";
  I(MyName, "%d,*", n);
  if(W->rgb) {
    Yprintf(MyName, "not in CMap mode.\n");
    exit(1);
  }
  x = X(W->xc);
  y = Y(W->yc);
  rect_write(MyName, x, y, x + n, y + 1, 2, data);
  W->xc += n / W->xf;
}

Int32 readRGB(Int16 n, RGBvalue r[], RGBvalue g[], RGBvalue b[]) {
  Screencoord x, y;
  Int32 ret, *data, i;
  const char * MyName = "readRGB";
  I(MyName, "%d,*,*,*", n);
  if(!W->rgb) {
    Yprintf(MyName, "not in RGB mode.\n");
    exit(1);
  }
  x = X(W->xc);
  y = Y(W->yc);
  W->xc += n / W->xf;
  if(NULL == (data = (Int32*)malloc(n * sizeof(Int32)))) {
    Yprintf(MyName, "out of memory.\n");
    exit(1);
  }
  ret = rect_read(MyName, x, y, x + n, y + 1, 4, data);
  for(i = 0; i < n; i++) {
    r[i] = (data[i] >>  0) & 0xFF;
    g[i] = (data[i] >>  8) & 0xFF;
    b[i] = (data[i] >> 16) & 0xFF;
  }
  free(data);
  return ret;
}

void writeRGB(Int16 n, RGBvalue r[], RGBvalue g[], RGBvalue b[]) {
  Screencoord x, y;
  Int32 *data, i;
  const char * MyName = "writeRGB";
  I(MyName, "%d,*,*,*", n);
  if(!W->rgb) {
    Yprintf(MyName, "not in RGB mode.\n");
    exit(1);
  }
  x = X(W->xc);
  y = Y(W->yc);
  W->xc += n / W->xf;
  if(NULL == (data = (Int32*)malloc(n * sizeof(Int32)))) {
    Yprintf(MyName, "out of memory.\n");
    exit(1);
  }
  for(i = 0; i < n; i++) {
    data[i] = (b[i] << 16) + (g[i] << 8) + r[i];
  }
  rect_write(MyName, x, y, x + n, y + 1, 4, data);  
  free(data);
}

#ifdef OGL
static Int32 blendfunction_translate(const char *caller, Int32 factor) {
  Int32 r;
  switch (factor) {
  case BF_ZERO:       r = GL_ZERO;                break;
  case BF_ONE:        r = GL_ONE;                 break;
  case BF_SC:         r = GL_SRC_COLOR;           break;
  case BF_MSC:        r = GL_ONE_MINUS_SRC_COLOR; break;
  case BF_SA:         r = GL_SRC_ALPHA;           break;
  case BF_MSA:        r = GL_ONE_MINUS_SRC_ALPHA; break;
  case BF_DA:         r = GL_DST_ALPHA;           break;
  case BF_MDA:        r = GL_ONE_MINUS_DST_ALPHA; break;
  case BF_DC:         r = GL_DST_COLOR;           break;
  case BF_MDC:        r = GL_ONE_MINUS_DST_COLOR; break;
  case BF_MIN_SA_MDA: r = GL_SRC_ALPHA_SATURATE;  break;
  default:
    Yprintf(caller, "invalid argument %d\n", factor);
    r = 0; break;
  }
  return r;
}
#endif

void blendfunction(Int32 sfactor, Int32 dfactor) {
  const char * MyName = "blendfunction";
  I(MyName, "%d,%d", sfactor, dfactor);
  IFOGL(glBlendFunc(blendfunction_translate(MyName, sfactor),
		    blendfunction_translate(MyName, dfactor));
	if (sfactor == BF_ONE && dfactor == BF_ZERO) glDisable(GL_BLEND);
	else                                         glEnable (GL_BLEND),
	NI(MyName));
}

/* $Log: color.c,v $
 * Revision 4.9  2007-05-08 11:01:05+02  fred
 * Version 4.2d
 *
 * Revision 4.8  2006-01-10 10:15:20+01  fred
 * Added blendfunction()
 *
 * Revision 4.7  2005-02-15 11:30:39+01  fred
 * *** empty log message ***
 *
 * Revision 4.6  1998-10-27 17:30:33+01  fred
 * bug in rect_write(), picture written one pixel too high
 *
 * Revision 4.5  1997-09-17 14:16:25+02  fred
 * *** empty log message ***
 *
 * Revision 4.4  1997-07-07 11:09:39+02  fred
 * *** empty log message ***
 *
 * Revision 4.3  1996-10-22 19:00:06+02  fred
 * rectzoom added
 * */
