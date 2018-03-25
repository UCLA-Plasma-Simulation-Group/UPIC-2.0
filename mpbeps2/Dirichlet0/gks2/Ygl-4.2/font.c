/*
 *    Ygl: Run GL programs with standard X11 and/or OpenGL routines.
 *    (C) Fred Hucht 1993-2007
 *    EMail: fred<at>thp.Uni-Duisburg.de
 */

static const char vcid[] = "$Id: font.c,v 4.3 2007-05-08 11:01:38+02 fred Exp $";

#include "header.h"

#ifdef OGL
static GLuint font_base = 0;
#endif

void loadXfont(Int32 id, Char8 *name) {
  int i;
  XFontStruct *fs;
  const char * MyName = "loadXfont";
  I(MyName, "%d,'%s'", id, name);
  if(NULL == (fs = XLoadQueryFont(D, name))) {
    Yprintf(MyName, "can't find font '%s'.\n", name);
    return;
  }

  if(Ygl.Fonts == NULL) { /* initialize */
    i = Ygl.LastFont = 0;
    Ygl.Fonts = (YglFont*)malloc(sizeof(YglFont));
  } else {
    for(i = Ygl.LastFont; i >= 0 && Ygl.Fonts[i].id != id; i--);
    if(i < 0) { /* not found */
      i = ++Ygl.LastFont;
      Ygl.Fonts = (YglFont*)realloc(Ygl.Fonts,
				    (Ygl.LastFont + 1) * sizeof(YglFont));
    }
  }
  
  if(Ygl.Fonts == NULL) {
    Yprintf(MyName, "can't allocate memory for font '%s'.\n", name);
    exit(-1);
  }
  
  Ygl.Fonts[i].fs = fs;
  Ygl.Fonts[i].id = id;
#ifdef DEBUG
  fprintf(stderr, 
	  "loadXfont: name = '%s', fs = 0x%x, id = %d.\n", 
	  name, fs, id);
#endif
}

#ifdef OGL
static void oglSetFont(const char *caller, int i) {
  int first;
  int last;
  if (font_base == 0) {
    /* initialize */
    font_base = glGenLists(256);
    if (!glIsList(font_base)) {
      Yprintf(caller, "out of display lists (%d).\n", font_base);
      exit (1);
    }
  }
  first = Ygl.Fonts[i].fs->min_char_or_byte2;
  last  = Ygl.Fonts[i].fs->max_char_or_byte2;
  glXUseXFont(Ygl.Fonts[i].fs->fid, first, last-first+1, font_base+first);
}
#endif

void font(Int16 id) {
  int i = Ygl.LastFont;
  const char * MyName = "font";
  I(MyName, "%d", id);
  while(i > 0 && Ygl.Fonts[i].id != id) i--;
  W->font = i;
#ifdef DEBUG
  fprintf(stderr, "font: id = %d, W->font = %d, fid = 0x%x.\n",
	  id, i, Ygl.Fonts[i].fs->fid);
#endif
  
  IFOGL(oglSetFont(MyName, i),
	XSetFont(D, W->chargc, Ygl.Fonts[i].fs->fid));
}

Int32 getfont(void) {
  const char * MyName = "getfont";
  I(MyName, "");
  return Ygl.Fonts[W->font].id;
}

void getfontencoding(char *r) {
  XFontStruct *fs;
  XFontProp *fp;
  int i;
  Atom fontatom;
  char *name, *np = NULL, *rp = r;
  
  const char * MyName = "getfontencoding";
  I(MyName, "'%s'", r);
  fs = Ygl.Fonts[W->font].fs;
  fontatom = XInternAtom(D, "FONT", False);
  
  for (i = 0, fp = fs->properties; i < fs->n_properties; i++, fp++) {
    if (fp->name == fontatom) {
      np = name = XGetAtomName(D, fp->card32);
      i = 0;
      while(i < 13 && *np != 0) if(*np++ == '-') i++;
      do {
	if(*np != '-') *rp++ = *np;
      }
      while(*np++ != 0);
      XFree(name);
    }
  }
  if(np == NULL) {
    Yprintf(MyName, "can't determine fontencoding.\n");
    *r = '\0';
  }
}

Int32 getheight(void) {
  const char * MyName = "getheight";
  I(MyName, "");
  return(Ygl.Fonts[W->font].fs->ascent +
	 Ygl.Fonts[W->font].fs->descent);
}

Int32 getdescender(void) {
  const char * MyName = "getdescender";
  I(MyName, "");
  return Ygl.Fonts[W->font].fs->descent;
}

Int32 strwidth(Char8 *string) {
  const char * MyName = "strwidth";
  I(MyName, "'%s'", string);
  return XTextWidth(Ygl.Fonts[W->font].fs, string, strlen(string));
}

void charstr(Char8 *Text) {
  const char * MyName = "charstr";
  I(MyName, "'%s'", Text);
  IFOGL(/* OpenGL */
	glRasterPos3f(W->xc, W->yc, W->zc);
	glPushAttrib(GL_LIST_BIT);
	glListBase(font_base);
	glCallLists(strlen(Text), GL_UNSIGNED_BYTE, (GLubyte *)Text);
	glPopAttrib(),
	/* X11 */
	if(!(W->rgb || Ygl.GC)) { /* set text color to active color */
	  XSetForeground(D, W->chargc, YGL_COLORS(W->color));
	}
	XDrawString(D, W->draw, W->chargc, X(W->xc), Y(W->yc),
		    Text, strlen(Text))
	);
  W->xc += strwidth(Text) / W->xf;
  F;
}

void cmov2 (Coord  x, Coord  y ) { I("cmov2" , "%g,%g", x, y); W->xc = x; W->yc = y; IFOGL(W->zc = 0,/**/); }
void cmov2i(Icoord x, Icoord y ) { I("cmov2i", "%d,%d", x, y); W->xc = x; W->yc = y; IFOGL(W->zc = 0,/**/); }
void cmov2s(Scoord x, Scoord y ) { I("cmov2s", "%d,%d", x, y); W->xc = x; W->yc = y; IFOGL(W->zc = 0,/**/); }
#ifdef OGL
void cmov  (Coord  x, Coord  y, Coord  z) { I("cmov" , "%g,%g,%g", x, y, z); W->xc = x; W->yc = y; W->zc = z; }
void cmovi (Icoord x, Icoord y, Icoord z) { I("cmovi", "%d,%d,%d", x, y, z); W->xc = x; W->yc = y; W->zc = z; }
void cmovs (Scoord x, Scoord y, Scoord z) { I("cmovs", "%d,%d,%d", x, y, z); W->xc = x; W->yc = y; W->zc = z; }
#endif

void getcpos(Screencoord *x, Screencoord *y) {
  /* Note: Due to a bug in GL getcpos() returns coordinates
   * relative to the screen, not to the window. We emulate this. :-((
   */
  const char *MyName = "getcpos";
  Int32 xo, yo;
  I(MyName, "*,*");
  getorigin(&xo, &yo);
  *x = xo + XR(W->xc);
  *y = yo + YR(W->yc);
#ifdef DEBUG
  fprintf(stderr, "getcpos: *x = %d, *y = %d\n", *x, *y);
#endif
}
