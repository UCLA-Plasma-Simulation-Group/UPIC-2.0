/*
 *    Ygl: Run GL programs with standard X11 and/or OpenGL routines.
 *    (C) Fred Hucht 1993-2007
 *    EMail: fred<at>thp.Uni-Duisburg.de
 */

static const char vcid[] = "$Id: draw.c,v 4.8 2007-05-08 11:01:24+02 fred Exp $";

#include "header.h"
#ifdef DEBUG /* see arc_ogl() */
# include <math.h>
#endif

#ifdef X11
#define DWGC D,W->draw,W->gc
#endif
#define Coo Coord
#define Ico Icoord
#define Sco Scoord

/* Drawing */
void clear(void) {
  I("clear", "");
  IFOGL(
	if (W->rgb || Ygl.EmulateCmap) {
	  GLfloat color[4]; /* Ask OpenGL for the current color */
	  glGetFloatv(GL_CURRENT_COLOR, color);
	  glClearColor(color[0], color[1], color[2], color[3]);
	} else {
	  GLfloat color;
	  glGetFloatv(GL_CURRENT_INDEX, &color);
	  glClearIndex(color);
	}
	glClear(GL_COLOR_BUFFER_BIT);
	,
	XFillRectangle(DWGC, 0, 0, W->xm, W->ym)
	);
  F;
}

#define  PNT2_FN(t, x, y) IFOGL(\
glBegin(GL_POINTS); glVertex2##t(W->xp = (x), W->yp = (y)); glEnd(),\
if(W->pm) XFillRectangle(DWGC, XS(x), YS(y), 1, 1);else XDrawPoint(DWGC, XS(x), YS(y)))

#define  PNT_FN(t, x, y, z) IFOGL(\
glBegin(GL_POINTS); glVertex3##t(W->xp = (x), W->yp = (y), W->zp = (z)); glEnd(),\
if(W->pm) XFillRectangle(DWGC, XS(x), YS(y), 1, 1);else XDrawPoint(DWGC, XS(x), YS(y)))
     
#define DRAW2_FN(t, x, y) IFOGL(\
glBegin(GL_LINES); glVertex2##t(W->xp, W->yp); glVertex2##t(x, y); glEnd(),\
XDrawLine(DWGC, X(W->xp), Y(W->yp), X(x), Y(y)))

#define DRAW_FN(t, x, y, z) IFOGL(\
glBegin(GL_LINES); glVertex3##t(W->xp, W->yp, W->zp); glVertex3##t(x, y, z); glEnd(),\
XDrawLine(DWGC, X(W->xp), Y(W->yp), X(x), Y(y)))

#define  RECT_FN(t, x1, y1, x2, y2) IFOGL(\
glBegin(GL_LINE_LOOP);glVertex2##t(x1, y1);glVertex2##t(x2, y1);glVertex2##t(x2, y2);glVertex2##t(x1, y2);glEnd(),\
XDrawRectangle(DWGC, X(MIN(x1,x2)), Y(MAX(y1,y2)), XR(ABS((x2)-(x1))), YR(ABS((y2)-(y1)))))

#define RECTF_FN(t, x1, y1, x2, y2) IFOGL(\
glBegin(GL_POLYGON  );glVertex2##t(x1, y1);glVertex2##t(x2, y1);glVertex2##t(x2, y2);glVertex2##t(x1, y2);glEnd(),\
XFillRectangle(DWGC, X(MIN(x1,x2)), Y(MAX(y1,y2)), XR(ABS((x2)-(x1))), YR(ABS((y2)-(y1)))))
     
#define ARC_FN(x, y, r, s, e) IFOGL(\
arc_ogl(x, y, r, s, e, 0),\
r=MAX(0,r); XDrawArc(DWGC, X((x)-(r)), Y((y)+(r)), XR(2*(r)), YR(2*(r)), 64*(s)/10, 64*((e)-(s))/10))
#define ARCF_FN(x, y, r, s, e) IFOGL(\
arc_ogl(x, y, r, s, e, 1),\
r=MAX(0,r); XFillArc(DWGC, X((x)-(r)), Y((y)+(r)), XR(2*(r)), YR(2*(r)), 64*(s)/10, 64*((e)-(s))/10))

/* Points */
void pnt2  (Coo x, Coo y       ) { I("pnt2" , "args"); PNT2_FN(f, x, y); F; }
void pnt2i (Ico x, Ico y       ) { I("pnt2i", "args"); PNT2_FN(i, x, y); F; }
void pnt2s (Sco x, Sco y       ) { I("pnt2s", "args"); PNT2_FN(s, x, y); F; }
void pnt   (Coo x, Coo y, Coo z) { I("pnt"  , "args"); PNT_FN(f, x, y, z); F; }
void pnti  (Ico x, Ico y, Ico z) { I("pnti" , "args"); PNT_FN(i, x, y, z); F; }
void pnts  (Sco x, Sco y, Sco z) { I("pnts" , "args"); PNT_FN(s, x, y, z); F; }

/* Lines */
void move2 (Coo x, Coo y       ) { SCP( =x, =y, =0);}
void move2i(Ico x, Ico y       ) { SCP( =x, =y, =0);}
void move2s(Sco x, Sco y       ) { SCP( =x, =y, =0);}
void move  (Coo x, Coo y, Coo z) { SCP( =x, =y, =z);}
void movei (Ico x, Ico y, Ico z) { SCP( =x, =y, =z);}
void moves (Sco x, Sco y, Sco z) { SCP( =x, =y, =z);}

void rmv2  (Coo x, Coo y       ) { SCP(+=x,+=y,+=0);}
void rmv2i (Ico x, Ico y       ) { SCP(+=x,+=y,+=0);}
void rmv2s (Sco x, Sco y       ) { SCP(+=x,+=y,+=0);}
void rmv   (Coo x, Coo y, Coo z) { SCP(+=x,+=y,+=z);}
void rmvi  (Ico x, Ico y, Ico z) { SCP(+=x,+=y,+=z);}
void rmvs  (Sco x, Sco y, Sco z) { SCP(+=x,+=y,+=z);}

void draw2 (Coo x, Coo y       ) { I("draw2" , "args"); DRAW2_FN(f, x, y); SCP( =x, =y, =0); F;}
void draw2i(Ico x, Ico y       ) { I("draw2i", "args"); DRAW2_FN(i, x, y); SCP( =x, =y, =0); F;}
void draw2s(Sco x, Sco y       ) { I("draw2s", "args"); DRAW2_FN(s, x, y); SCP( =x, =y, =0); F;}
void draw  (Coo x, Coo y, Coo z) { I("draw"  , "args");  DRAW_FN(f, x, y, z); SCP( =x, =y, =z); F;}
void drawi (Ico x, Ico y, Ico z) { I("drawi" , "args");  DRAW_FN(i, x, y, z); SCP( =x, =y, =z); F;}
void draws (Sco x, Sco y, Sco z) { I("draws" , "args");  DRAW_FN(s, x, y, z); SCP( =x, =y, =z); F;}

void rdr2  (Coo x, Coo y       ) { I("rdr2"  , "args"); DRAW2_FN(f, W->xp+x, W->yp+y); SCP(+=x,+=y,+=0); F;}
void rdr2i (Ico x, Ico y       ) { I("rdr2i" , "args"); DRAW2_FN(i, W->xp+x, W->yp+y); SCP(+=x,+=y,+=0); F;}
void rdr2s (Sco x, Sco y       ) { I("rdr2s" , "args"); DRAW2_FN(s, W->xp+x, W->yp+y); SCP(+=x,+=y,+=0); F;}
void rdr   (Coo x, Coo y, Coo z) { I("rdr"   , "args");  DRAW_FN(f, W->xp+x, W->yp+y, W->zp+z); SCP(+=x,+=y,+=z); F;}
void rdri  (Ico x, Ico y, Ico z) { I("rdri"  , "args");  DRAW_FN(i, W->xp+x, W->yp+y, W->zp+z); SCP(+=x,+=y,+=z); F;}
void rdrs  (Sco x, Sco y, Sco z) { I("rdrs"  , "args");  DRAW_FN(s, W->xp+x, W->yp+y, W->zp+z); SCP(+=x,+=y,+=z); F;}

/* Arcs & Circles */
#ifdef OGL
static void arc_ogl(double x, double y, double r,
		    int s, int e, int filled) {
#if 0
  if(r > 0) {
    glPushMatrix();
    glTranslatef(x, y, 0.0);
    /* Mesa requires 2 loops -----\ */
    gluPartialDisk(filled ? W->circf : W->circ,
		   0.0, r, 12, 2, 0.1 * s, 0.1 * (e - s));
    glPopMatrix();
  }
#else
  int i;
#define NARC 32
#define SIN(i) ((1 - 2 * (((i) / (NARC/2)) & 1)) * Sin[(i) % (NARC/2)])
#define COS(i) SIN(i + NARC/4)
#ifdef DEBUG_NO
# define DEBUGPRINT(i) fprintf(stderr, "%g,%g %g %f,%f %f,%f\n",\
			       x, y, r, glx, gly,\
			       x + r * cos(2 * M_PI * i / NARC),\
			       y + r * sin(2 * M_PI * i / NARC));
#else
# define DEBUGPRINT(i)
#endif
  const double Sin[NARC] = {
    0.0000000000000000, 0.1950903220161283, 0.3826834323650898, 0.5555702330196022,
    0.7071067811865476, 0.8314696123025452, 0.9238795325112867, 0.9807852804032304,
    1.0000000000000000, 0.9807852804032304, 0.9238795325112867, 0.8314696123025455,
    0.7071067811865476, 0.5555702330196022, 0.3826834323650899, 0.1950903220161286
  };
  
  if(r > 0) {
    GLdouble glx, gly;
    glBegin(filled ? GL_POLYGON : GL_LINE_STRIP);
    if(s == 0 && e == 3600) { /* circ */
      for(i = 0; i <= NARC; i++) {
	glx = x + r * COS(i);
	gly = y + r * SIN(i);
	glVertex2d(glx, gly);
      }
    } else { /* arc */
      double sd = (double)s/3600*NARC;
      double ed = (double)e/3600*NARC;
      double d;
      i = sd;
      d = sd - i;
      /* interpolate first point */
      glx = x + r * ((1 - d) * COS(i) + d * COS(i+1));
      gly = y + r * ((1 - d) * SIN(i) + d * SIN(i+1));
      DEBUGPRINT(sd);
      glVertex2d(glx, gly);
      for(i++; i <= (int)ed; i++) {
	glx = x + r * COS(i);
	gly = y + r * SIN(i);
	DEBUGPRINT(i);
	glVertex2d(glx, gly);
      }
      i = ed;
      d = ed - i;
      /* interpolate last point */
      glx = x + r * ((1 - d) * COS(i) + d * COS(i+1));
      gly = y + r * ((1 - d) * SIN(i) + d * SIN(i+1));
      DEBUGPRINT(ed);
      glVertex2d(glx, gly);
    }
    glEnd();
  }
#endif
}
#endif
void arc   (Coo x, Coo y, Coo r, Angle s, Angle e) { I("arc"   , "args");  ARC_FN(x, y, r, s, e); F;}
void arci  (Ico x, Ico y, Ico r, Angle s, Angle e) { I("arci"  , "args");  ARC_FN(x, y, r, s, e); F;}
void arcs  (Sco x, Sco y, Sco r, Angle s, Angle e) { I("arcs"  , "args");  ARC_FN(x, y, r, s, e); F;}
void arcf  (Coo x, Coo y, Coo r, Angle s, Angle e) { I("arcf"  , "args"); ARCF_FN(x, y, r, s, e); F;}
void arcfi (Ico x, Ico y, Ico r, Angle s, Angle e) { I("arcfi" , "args"); ARCF_FN(x, y, r, s, e); F;}
void arcfs (Sco x, Sco y, Sco r, Angle s, Angle e) { I("arcfs" , "args"); ARCF_FN(x, y, r, s, e); F;}

void circ  (Coo x, Coo y, Coo r) { I("circ"  , "args");  ARC_FN(x, y, r, 0, 3600); F;}
void circi (Ico x, Ico y, Ico r) { I("circi" , "args");  ARC_FN(x, y, r, 0, 3600); F;}
void circs (Sco x, Sco y, Sco r) { I("circs" , "args");  ARC_FN(x, y, r, 0, 3600); F;}
void circf (Coo x, Coo y, Coo r) { I("circf" , "args"); ARCF_FN(x, y, r, 0, 3600); F;}
void circfi(Ico x, Ico y, Ico r) { I("circfi", "args"); ARCF_FN(x, y, r, 0, 3600); F;}
void circfs(Sco x, Sco y, Sco r) { I("circfs", "args"); ARCF_FN(x, y, r, 0, 3600); F;}

/* Rects & Boxes */
void rect  (Coo x1, Coo y1, Coo x2, Coo y2) { I("rect"  , "args");  RECT_FN(f, x1, y1, x2, y2); F;}
void sbox  (Coo x1, Coo y1, Coo x2, Coo y2) { I("sbox"  , "args");  RECT_FN(f, x1, y1, x2, y2); F;}
void recti (Ico x1, Ico y1, Ico x2, Ico y2) { I("recti" , "args");  RECT_FN(i, x1, y1, x2, y2); F;}
void sboxi (Ico x1, Ico y1, Ico x2, Ico y2) { I("sboxi" , "args");  RECT_FN(i, x1, y1, x2, y2); F;}
void rects (Sco x1, Sco y1, Sco x2, Sco y2) { I("rects" , "args");  RECT_FN(s, x1, y1, x2, y2); F;}
void sboxs (Sco x1, Sco y1, Sco x2, Sco y2) { I("sboxs" , "args");  RECT_FN(s, x1, y1, x2, y2); F;}
void rectf (Coo x1, Coo y1, Coo x2, Coo y2) { I("rectf" , "args"); RECTF_FN(f, x1, y1, x2, y2); F;}
void sboxf (Coo x1, Coo y1, Coo x2, Coo y2) { I("sboxf" , "args"); RECTF_FN(f, x1, y1, x2, y2); F;}
void rectfi(Ico x1, Ico y1, Ico x2, Ico y2) { I("rectfi", "args"); RECTF_FN(i, x1, y1, x2, y2); F;}
void sboxfi(Ico x1, Ico y1, Ico x2, Ico y2) { I("sboxfi", "args"); RECTF_FN(i, x1, y1, x2, y2); F;}
void rectfs(Sco x1, Sco y1, Sco x2, Sco y2) { I("rectfs", "args"); RECTF_FN(s, x1, y1, x2, y2); F;}
void sboxfs(Sco x1, Sco y1, Sco x2, Sco y2) { I("sboxfs", "args"); RECTF_FN(s, x1, y1, x2, y2); F;}

static Int32 PMode = Convex;

#ifdef X11
static Ulong q_len = 0, q_max = 0;
static XPoint *queue = NULL;

static void q_add(const char *caller, int x, int y) {
  if(q_max <= q_len) {
    q_max += 0x1000;
    if(queue == NULL)	/* Make new queue */
      queue = (XPoint*)malloc(q_max * sizeof(XPoint));
    else		/* Resize queue */
      queue = (XPoint*)realloc(queue, q_max * sizeof(XPoint));
    if(queue == NULL) {
      Yprintf(caller, "out of memory.\n");
      exit(1);
    }
  }
  /* Add point */
  queue[q_len].x = x;
  queue[q_len].y = y;
  q_len++;
}
#endif

void concave(Int32 bool) { PMode = bool ? Complex : Convex;}

void  pmv2 (Coo x, Coo y) { const char *MyName = "pmv2" ; I(MyName, "args"); SCP(=x, =y, =0); IFOGL(NI(MyName),q_len = 0; q_add(MyName, X(x), Y(y)));}
void  pmv2i(Ico x, Ico y) { const char *MyName = "pmv2i"; I(MyName, "args"); SCP(=x, =y, =0); IFOGL(NI(MyName),q_len = 0; q_add(MyName, X(x), Y(y)));}
void  pmv2s(Sco x, Sco y) { const char *MyName = "pmv2s"; I(MyName, "args"); SCP(=x, =y, =0); IFOGL(NI(MyName),q_len = 0; q_add(MyName, X(x), Y(y)));}

void rpmv2 (Coo x, Coo y) { const char *MyName = "rpmv2" ; I(MyName, "args"); SCP(+=x, +=y, =0); IFOGL(NI(MyName),q_len = 0; q_add(MyName, X(x), Y(y)));}
void rpmv2i(Ico x, Ico y) { const char *MyName = "rpmv2i"; I(MyName, "args"); SCP(+=x, +=y, =0); IFOGL(NI(MyName),q_len = 0; q_add(MyName, X(x), Y(y)));}
void rpmv2s(Sco x, Sco y) { const char *MyName = "rpmv2s"; I(MyName, "args"); SCP(+=x, +=y, =0); IFOGL(NI(MyName),q_len = 0; q_add(MyName, X(x), Y(y)));}

void  pdr2 (Coo x, Coo y) { const char *MyName = "pdr2" ; I(MyName, "args"); IFOGL(NI(MyName),q_add(MyName, X(x), Y(y)));}
void  pdr2i(Ico x, Ico y) { const char *MyName = "pdr2i"; I(MyName, "args"); IFOGL(NI(MyName),q_add(MyName, X(x), Y(y)));}
void  pdr2s(Sco x, Sco y) { const char *MyName = "pdr2s"; I(MyName, "args"); IFOGL(NI(MyName),q_add(MyName, X(x), Y(y)));}

void rpdr2 (Coo x, Coo y) { const char *MyName = "rpdr2" ; I(MyName, "args"); IFOGL(NI(MyName),q_add(MyName, queue[q_len-1].x + X(x), queue[q_len-1].y - Y(y)));}
void rpdr2i(Ico x, Ico y) { const char *MyName = "rpdr2i"; I(MyName, "args"); IFOGL(NI(MyName),q_add(MyName, queue[q_len-1].x + X(x), queue[q_len-1].y - Y(y)));}
void rpdr2s(Sco x, Sco y) { const char *MyName = "rpdr2s"; I(MyName, "args"); IFOGL(NI(MyName),q_add(MyName, queue[q_len-1].x + X(x), queue[q_len-1].y - Y(y)));}

void pclos(void) { const char *MyName = "pclos";I(MyName, ""); IFOGL(NI(MyName),if(q_len) XFillPolygon(DWGC, queue, q_len, PMode, CoordModeOrigin)); F;}

void  pmv (Coo x, Coo y, Coo z) { const char *MyName = "pmv" ; NI(MyName);}
void  pmvi(Ico x, Ico y, Ico z) { const char *MyName = "pmvi"; NI(MyName);}
void  pmvs(Sco x, Sco y, Sco z) { const char *MyName = "pmvs"; NI(MyName);}

void rpmv (Coo x, Coo y, Coo z) { const char *MyName = "rpmv" ; NI(MyName);}
void rpmvi(Ico x, Ico y, Ico z) { const char *MyName = "rpmvi"; NI(MyName);}
void rpmvs(Sco x, Sco y, Sco z) { const char *MyName = "rpmvs"; NI(MyName);}

void  pdr (Coo x, Coo y, Coo z) { const char *MyName = "pdr" ; NI(MyName);}
void  pdri(Ico x, Ico y, Ico z) { const char *MyName = "pdri"; NI(MyName);}
void  pdrs(Sco x, Sco y, Sco z) { const char *MyName = "pdrs"; NI(MyName);}

void rpdr (Coo x, Coo y, Coo z) { const char *MyName = "rpdr" ; NI(MyName);}
void rpdri(Ico x, Ico y, Ico z) { const char *MyName = "rpdri"; NI(MyName);}
void rpdrs(Sco x, Sco y, Sco z) { const char *MyName = "rpdrs"; NI(MyName);}

#define MKP(n, p) { Int32 i; q_len = 0; for(i = 0; i < n; i++) { q_add(MyName, X(p[i][0]), Y(p[i][1]));}}

void poly2 (Int32 n, Coo p[][2]) { const char *MyName = "poly2" ; I(MyName, "args"); IFOGL(NI(MyName),MKP(n, p); q_add(MyName, XS(p[0][0]), YS(p[0][1])); XDrawLines(DWGC, queue, q_len, CoordModeOrigin)); F;}
void poly2i(Int32 n, Ico p[][2]) { const char *MyName = "poly2i"; I(MyName, "args"); IFOGL(NI(MyName),MKP(n, p); q_add(MyName, XS(p[0][0]), YS(p[0][1])); XDrawLines(DWGC, queue, q_len, CoordModeOrigin)); F;}
void poly2s(Int32 n, Sco p[][2]) { const char *MyName = "poly2s"; I(MyName, "args"); IFOGL(NI(MyName),MKP(n, p); q_add(MyName, XS(p[0][0]), YS(p[0][1])); XDrawLines(DWGC, queue, q_len, CoordModeOrigin)); F;}

void polf2 (Int32 n, Coo p[][2]) { const char *MyName = "polf2" ; I(MyName, "args"); IFOGL(NI(MyName),MKP(n, p); SCP(=p[0][0],=p[0][1], =0); XFillPolygon(DWGC, queue, q_len, PMode, CoordModeOrigin)); F;}
void polf2i(Int32 n, Ico p[][2]) { const char *MyName = "polf2i"; I(MyName, "args"); IFOGL(NI(MyName),MKP(n, p); SCP(=p[0][0],=p[0][1], =0); XFillPolygon(DWGC, queue, q_len, PMode, CoordModeOrigin)); F;}
void polf2s(Int32 n, Sco p[][2]) { const char *MyName = "polf2s"; I(MyName, "args"); IFOGL(NI(MyName),MKP(n, p); SCP(=p[0][0],=p[0][1], =0); XFillPolygon(DWGC, queue, q_len, PMode, CoordModeOrigin)); F;}

#ifdef OGL
static void polf_ogl (Int32 n, Coo p[][3]) {int i;glBegin(GL_POLYGON);for(i = 0; i < n; i++) glVertex3fv(p[i]);glEnd();}
static void polfi_ogl(Int32 n, Ico p[][3]) {int i;glBegin(GL_POLYGON);for(i = 0; i < n; i++) glVertex3iv(p[i]);glEnd();}
static void polfs_ogl(Int32 n, Sco p[][3]) {int i;glBegin(GL_POLYGON);for(i = 0; i < n; i++) glVertex3sv(p[i]);glEnd();}
#endif

void polf  (Int32 n, Coo p[][3]) { const char *MyName = "polf" ; I(MyName, "args"); IFOGL(polf_ogl (n, p),MKP(n, p); SCP(=p[0][0],=p[0][1], =0); XFillPolygon(DWGC, queue, q_len, PMode, CoordModeOrigin)); F;}
void polfi (Int32 n, Ico p[][3]) { const char *MyName = "polfi"; I(MyName, "args"); IFOGL(polfi_ogl(n, p),MKP(n, p); SCP(=p[0][0],=p[0][1], =0); XFillPolygon(DWGC, queue, q_len, PMode, CoordModeOrigin)); F;}
void polfs (Int32 n, Sco p[][3]) { const char *MyName = "polfs"; I(MyName, "args"); IFOGL(polfs_ogl(n, p),MKP(n, p); SCP(=p[0][0],=p[0][1], =0); XFillPolygon(DWGC, queue, q_len, PMode, CoordModeOrigin)); F;}

/* Vertex graphics */
void bgnpoint     (void) { const char *MyName = "bgnpoint"     ; I(MyName, ""); W->vmode = VertexPoint;IFOGL(glBegin(GL_POINTS    ), q_len = 0);}
void bgnline      (void) { const char *MyName = "bgnline"      ; I(MyName, ""); W->vmode = VertexLine; IFOGL(glBegin(GL_LINE_STRIP), q_len = 0);}
void bgnclosedline(void) { const char *MyName = "bgnclosedline"; I(MyName, ""); W->vmode = VertexCLine;IFOGL(glBegin(GL_LINE_LOOP ), q_len = 0);}
void bgnpolygon   (void) { const char *MyName = "bgnpolygon"   ; I(MyName, ""); W->vmode = VertexPoly; IFOGL(glBegin(GL_POLYGON   ), q_len = 0);}
void bgntmesh     (void) { const char *MyName = "bgntmesh"     ; I(MyName, ""); W->vmode = VertexTMesh;IFOGL(glBegin(GL_TRIANGLE_STRIP), NI(MyName));}

void endpoint     (void) { const char *MyName = "endpoint"     ; I(MyName, ""); if(W->vmode != VertexPoint) {Yprintf(MyName, "missing bgnpoint().\n"     ); exit(1);} IFOGL(glEnd(), if(q_len) {                                        XDrawPoints (DWGC, queue, q_len,        CoordModeOrigin); });F;}
void endline      (void) { const char *MyName = "endline"      ; I(MyName, ""); if(W->vmode != VertexLine ) {Yprintf(MyName, "missing bgnline().\n"      ); exit(1);} IFOGL(glEnd(), if(q_len) {                                        XDrawLines  (DWGC, queue, q_len,        CoordModeOrigin); });F;}
void endclosedline(void) { const char *MyName = "endclosedline"; I(MyName, ""); if(W->vmode != VertexCLine) {Yprintf(MyName, "missing bgnclosedline().\n"); exit(1);} IFOGL(glEnd(), if(q_len) { q_add(MyName, queue[0].x, queue[0].y); XDrawLines  (DWGC, queue, q_len,        CoordModeOrigin); });F;}
void endpolygon   (void) { const char *MyName = "endpolygon"   ; I(MyName, ""); if(W->vmode != VertexPoly ) {Yprintf(MyName, "missing bgnpolygon().\n"   ); exit(1);} IFOGL(glEnd(), if(q_len) {                                        XFillPolygon(DWGC, queue, q_len, PMode, CoordModeOrigin); });F;}
void endtmesh     (void) { const char *MyName = "endtmesh"     ; I(MyName, ""); if(W->vmode != VertexTMesh) {Yprintf(MyName, "missing bgntmesh().\n"     ); exit(1);} IFOGL(glEnd(), NI(MyName)); F;}

#if 0
#define V_FN(t, gt, v) if(W->vmode == VertexNone) { Yprintf(MyName, "not in vertex mode\n"); exit(1);}\
IFOGL(glVertex##t##v((gt *)v), q_add(MyName, X(v[0]), Y(v[1])))
#endif

float vo[2][3];
int vi;

#ifdef OGL
Float32 no[2][3];

void n3f(Float32 v[3]) {
  const char * MyName = "n3f";
  I(MyName, "{%g,%g,%g}", v[0], v[1], v[2]);
  if(W->vmode == VertexTMesh) 
    memcpy(no[vi], v, sizeof(no[0]));
  IFOGL(glNormal3fv(v),NI(MyName));
}

void normal(Coord v[3]) {
  const char * MyName = "normal";
  I(MyName, "{%g,%g,%g}", v[0], v[1], v[2]);
  IFOGL(glNormal3fv(v),NI(MyName));
}
#endif

#if 1
#define V_FN(t, gt, v)						\
switch (W->vmode) {						\
case VertexNone: 						\
  Yprintf(MyName, "not in vertex mode.\n");			\
  exit(1);							\
case VertexTMesh:						\
  vo[vi][0] = v[0];						\
  vo[vi][1] = v[1];						\
  vo[vi][2] = v[2];						\
  vi = 1 - vi;							\
}								\
IFOGL(glVertex##t##v((gt *)v), q_add(MyName, X(v[0]), Y(v[1])))
#endif

void swaptmesh(void) {
  const char *MyName = "swaptmesh";
  if(W->vmode != VertexTMesh) {
    Yprintf(MyName, "missing bgntmesh.\n");
    exit(1);
  }
  /* To swap the tmesh we simply resend the 2nd-last vertex */
  /*glNormal3fv(no[vi]);*/
  /*glVertex3fv(vo[vi]);*/
  /* v3f(vo[vi]); */
}

void v2s(Int16   v[2]) { const char *MyName = "v2s"; V_FN(2s, GLshort,  v);}
void v2i(Int32   v[2]) { const char *MyName = "v2i"; V_FN(2i, GLint,    v);}
void v2f(Float32 v[2]) { const char *MyName = "v2f"; V_FN(2f, GLfloat,  v);}
void v2d(Float64 v[2]) { const char *MyName = "v2d"; V_FN(2d, GLdouble, v);}
void v3s(Int16   v[3]) { const char *MyName = "v3s"; V_FN(3s, GLshort,  v);}
void v3i(Int32   v[3]) { const char *MyName = "v3i"; V_FN(3i, GLint,    v);}
void v3f(Float32 v[3]) { const char *MyName = "v3f"; V_FN(3f, GLfloat,  v);}
void v3d(Float64 v[3]) { const char *MyName = "v3d"; V_FN(3d, GLdouble, v);}
void v4s(Int16   v[4]) { const char *MyName = "v4s"; V_FN(4s, GLshort,  v);}
void v4i(Int32   v[4]) { const char *MyName = "v4i"; V_FN(4i, GLint,    v);}
void v4f(Float32 v[4]) { const char *MyName = "v4f"; V_FN(4f, GLfloat,  v);}
void v4d(Float64 v[4]) { const char *MyName = "v4d"; V_FN(4d, GLdouble, v);}

/* Extensions: Routines from X not in gl. Contributed by MiSt (michael<at>thp.Uni-Duisburg.de) */
#ifdef X11
void arcx  (Coo x,Coo y,Coo rx,Coo ry,Angle s,Angle e) { I("arcx" , "args"); XDrawArc(DWGC, X(x-rx), Y(y+ry), XR(2*MAX(0,rx)), YR(2*MAX(0,ry)), (64*s)/10, (64*(e-s))/10); F;}
void arcxi (Ico x,Ico y,Ico rx,Ico ry,Angle s,Angle e) { I("arcxi", "args"); XDrawArc(DWGC, X(x-rx), Y(y+ry), XR(2*MAX(0,rx)), YR(2*MAX(0,ry)), (64*s)/10, (64*(e-s))/10); F;}
void arcxs (Sco x,Sco y,Sco rx,Sco ry,Angle s,Angle e) { I("arcxs", "args"); XDrawArc(DWGC, X(x-rx), Y(y+ry), XR(2*MAX(0,rx)), YR(2*MAX(0,ry)), (64*s)/10, (64*(e-s))/10); F;}

void arcxf (Coo x,Coo y,Coo rx,Coo ry,Angle s,Angle e) { I("arcxf" , "args"); XFillArc(DWGC, X(x-rx), Y(y+ry), XR(2*MAX(0,rx)), YR(2*MAX(0,ry)), (64*s)/10, (64*(e-s))/10); F;}
void arcxfi(Ico x,Ico y,Ico rx,Ico ry,Angle s,Angle e) { I("arcxfi", "args"); XFillArc(DWGC, X(x-rx), Y(y+ry), XR(2*MAX(0,rx)), YR(2*MAX(0,ry)), (64*s)/10, (64*(e-s))/10); F;}
void arcxfs(Sco x,Sco y,Sco rx,Sco ry,Angle s,Angle e) { I("arcxfs", "args"); XFillArc(DWGC, X(x-rx), Y(y+ry), XR(2*MAX(0,rx)), YR(2*MAX(0,ry)), (64*s)/10, (64*(e-s))/10); F;}
#endif
