/*
 * xmap.c: display the active colormap in a [Y]gl window.
 * 
 * (C) Fred Hucht 1994, EMail: fred@hal6000.Uni-Duisburg.DE
 * 
 * Feel free to copy and redistribute in terms of the
 * GNU public license.
 */

#include <stdio.h>
/*#include <gl/gl.h>
  #include <gl/device.h>*/
#include <X11/Ygl.h>

#define SZ 160
#define SI (SZ/m)

#ifdef XCOLOR
#define COLOR		Xcolor
#define GETMCOLOR 	getmXcolor
#else
#define COLOR		color
#define GETMCOLOR 	getmcolor
#endif


void set_header(char *text) {
  color(BLACK); 
  rectfi( 0, SZ, SZ, SZ+20);
  color(WHITE); 
  cmov2i(10, SZ + 5); 
  charstr(text);
}

int main(int argc, char *argv[]) {
  Device dev;
  short val;
  int i, j, m, p, left = 0;
  Int16 x,y;
  Int32 win, index;
  Int32 xo, yo, xs, ys;
  char version[20], header[80], buf[80], *name = argv[0], *title, 
       *titles[]={"X11 ColorMap (<ESC> quit)",
		  "YGL ColorMap (<ESC> quit)",
		   "GL ColorMap (<ESC> quit)"};

  if(-1 == gversion(version)) { 
    fprintf(stderr, "%s: won't run on this display.\n", version); 
    exit(1);
  }
  
  minsize   (SZ, SZ + 20);
  /* stepunit  (SZ, SZ + 20); */
  keepaspect(SZ, SZ + 20);
  
  if(name[0] == 'g') title = titles[2];        /* glmap  */
  else if(name[1] == 'g') title = titles[1];   /* yglmap */
  else title = titles[0];                      /* xmap   */
  
  win = winopen(title);winconstraints();

  qdevice(WINQUIT);
  qdevice(KEYBD);
  qdevice(REDRAW);
  qdevice(INPUTCHANGE);
  qenter(REDRAW, win);
  qdevice(MOUSEX);
  qdevice(MOUSEY);
  qdevice(LEFTMOUSE);
  
  p = getplanes();
  m = 1 << (p/2);
  
  sprintf(header, "%s, %d bpp", version, p);
  
  while (1) {
    dev = qread(&val);
#ifdef DEBUG
    printf("%d %d\n", dev, val);
#endif
    switch(dev) {
    case INPUTCHANGE:
      if(val==0) set_header(header);
      break;
    case MOUSEX: x = val; goto mouse;
    case MOUSEY: y = val; goto mouse;
    mouse: {
      Int32 xw, yw, q;
      Int16 r,g,b;
      
      q = qtest();
      if(q == MOUSEX || q == MOUSEY) break;
      
      xw = x - xo;
      yw = y - yo;
      index = m * (yw*m/xs) + (xw*m/xs);
      if(index >= 0 && index < (1<<p)) {
	GETMCOLOR(index, &r, &g, &b);
#ifdef DEBUG
	printf("%d %d %d %d %d %d %d=(%d,%d,%d)\n", xs, ys, 
	       xw*m/xs, yw*m/xs, x, y, index, r, g, b);
#endif
	sprintf(buf, "%d=(%d,%d,%d)", index, r, g, b);
	set_header(buf);
	if(left) {
	  COLOR(index);
	  rectf(0, 0, SZ-1, SZ-1);
	}
      }
    }
      break;
    case KEYBD:
      if(val == '\033' || val == 'q') {
	gexit();
	exit(0);
      }
      break;
    case WINQUIT:
      gexit();
      exit(0);
      break;
    case LEFTMOUSE:
      if(val == 1) {
	left = 1;
	if(index >= 0 && index < (1<<p)) {
	  COLOR(index);
	  rectf(0, 0, SZ-1, SZ-1);
	}
	break;
      }
      /* else redraw */
      left = 0;
    case REDRAW:
      reshapeviewport();
      getorigin(&xo, &yo);
      getsize  (&xs, &ys);
      color(BLACK); clear();
      set_header(header);
      for(i = 0; i < m; i++) for(j = 0; j < m; j++) {
	COLOR(i * m + j);
	rectfi(j     * SI,
	       i     * SI,
	       (j+1) * SI - 1,
	       (i+1) * SI - 1);
      }
      break;
    }
  }
}
