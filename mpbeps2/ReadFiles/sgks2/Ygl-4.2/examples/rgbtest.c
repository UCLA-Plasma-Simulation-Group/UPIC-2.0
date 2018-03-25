/* rgbtest.c by Fred Hucht (C) 1993-2002.
 * Example program with two windows,
 * one in RGBmode and one in Colormap mode */

static const char vcid[] = "$Id: rgbtest.c,v 3.3 1996/07/18 16:35:57 fred Exp $";

#include <X11/Ygl.h>
#include <stdio.h>
#include <unistd.h>

#define STEP 2
#define SIZE (256 * STEP)

Int32 RGBwin, cwin;

void drawit(int g) {
  short r,b;
  reshapeviewport();
  for(r = 0; r < 256; r++) for(b = 0; b < 256; b++) {
    RGBcolor(r, g, b);
    rectfi(r*STEP, b*STEP, (r+1)*STEP, (b+1)*STEP);
  }
  sleep(1);
  RGBcolor(255,255,0);
  font(4711);
  cmov2i(10, 10);
  charstr("Hello in Yellow");
}

void drawitgrey(void) {
  short i;
  reshapeviewport();
  for(i = 0; i < 256; i++) {
    RGBcolor(i, i, i);
    rectfi(i*STEP, 0, (i+1)*STEP, 255);
  }
  RGBcolor(255,255,0);
  font(4711);
  cmov2i(10, 10);
  charstr("Hey in Grey");
}

void drawit2(void) {
  winset(cwin);
  reshapeviewport();
  color(BLACK); clear();
  color(RED); circf(0.0, 0.0, 0.8);
  sleep(1);
  winset(RGBwin);
}

int main() {
  Device dev;
  short val;
  int green = 128, openwins=0;
  Int32 xo, yo;
  Char8 *Title = "RGBmode Window";
  
  /* putenv("YGL_BS=1");*/
  /* setenv("YGL_BS","1",1); *//* use this ones under BSD unixes */
  
  minsize(512, 512);
  cwin = winopen("Colormap Window");
  openwins++;
  
  winconstraints(); /* remove all constraints */
  ortho2(-1.0, 1.0, -1.0, 1.0);
  color(BLUE); clear();
  
  /* noport(); */
  minsize(SIZE, SIZE);
#ifdef SUBWIN
  RGBwin = swinopen(cwin);
#else
  RGBwin = winopen("RGB mode");
#endif
  openwins++;
  RGBmode();
  gconfig();
  
  /*ortho2(0.0, SIZE - 1, 0.0, SIZE - 1);*/
  
  qdevice(REDRAW);
  qdevice(KEYBD);
  qdevice(WINQUIT);
  qdevice(LEFTMOUSE);
  tie(LEFTMOUSE, MOUSEX, MOUSEY); /* report mouse position on buttonpress */
  
  winset(RGBwin);
  getsize(&xo,&yo);
  
  loadXfont(4711, "-*-times-medium-r-*-*-*-140-*-*-*-*-iso8859-1");
  
  while (1) {
    dev = qread(&val);
    switch(dev) {
    case LEFTMOUSE:
      if(val==1) { /* If button pressed */
	Int32 x, y, xs, ys;
	Int16 r, g, b, mx, my;
	char buf[80];
	getorigin(&x, &y);
	getsize(&xs, &ys);
	if(MOUSEX != qread(&mx) || MOUSEY != qread(&my))
	  printf("tie doesn't work, strange...\n");
	x = (mx - x) * xo / xs;
	y = (my - y) * yo / ys;
	/* Next line is to get the real color (not the desired)
	 * on <24 bit Visuals */
	RGBcolor(x, green, y); gRGBcolor(&r, &g, &b);
	sprintf(buf,"Color at (%d,%d) is (%d,%d,%d).",x,y,r,g,b);
	wintitle(buf);
      } else { /* Restore original title */
	wintitle(Title);
      }
      break;
    case KEYBD:
      switch(val) {
      case '\033' :
	gexit();
	return(0);
	break;
      case '0':	green  =   0; drawit(green); break;
      case '9':	green  = 255; drawit(green); break;
      case '+':	green +=  16; drawit(green); break;
      case '-':	green -=  16; drawit(green); break;
      case 'g': drawitgrey(); break;
      case 'p':
	printf("Piping RGB window through 'ppmtogif > RGB.gif'\n");
	gl2ppm("| ppmtogif > RGB.gif");
	ringbell();
	break;
      }
      break;
    case WINQUIT:
      winclose(val);
      if(0 == --openwins) {
	gexit();
	return(0);
      }
      break;
    case REDRAW:
      if(val == cwin)
	drawit2();
      else if(val == RGBwin)
	drawit(green);
      break;
    }
  }
}
