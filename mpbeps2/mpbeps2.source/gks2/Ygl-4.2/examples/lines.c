/* lines.c by Fred Hucht (C) 1993-2005.
 * Example program to compare performance */

static const char vcid[] = "$Id: lines.c,v 3.3 2005-02-15 11:29:09+01 fred Exp $";

#include <X11/Ygl.h>
#include <stdio.h>

#define SIZE 400	/* Size of main window */
#define STEP 4		/* Width of lines */

#ifndef TRUE
# define TRUE  (0==0)
# define FALSE (!TRUE)
#endif

int ncolors;
Int32 win, swin;

void drawit(Int32 w) {
  int i, j, col;
  winset(w);
  if(w == win) {
    char *text = "This is font Times 14. Press ESC to exit...";
    color(WHITE); clear();
    col = 0;
    for(i=0; i<SIZE; i+=STEP) {
      color(col);
      /* The next line is a loop for performance comparsion. */
      for(j=0; j<SIZE; j+=STEP) rectfi(i, j, i+STEP, j+STEP);
      
      /* Of course it's much faster to use this... */
      /* rectfi(i, 0, i+STEP, SIZE); */
      col = (col+1) % ncolors;
    }
    color(WHITE); rectfi(5, 5, 10+strwidth(text), 30);
    cmov2i(10, 10); color(BLUE); charstr(text);
  } else if(w == swin) {
    color(WHITE); clear();
    color(RED); circf(50, 50, 10);
  } else {
    fprintf(stderr, "Unknown window: %d\n", w);
  }
}

int main() {
  char buf1[20], buf2[80];
  int running = TRUE;

  gversion(buf1);
  sprintf(buf2, "Lines drawn with %s",buf1);

  minsize(SIZE, SIZE);
  win  = winopen(buf2);
  
  prefposition(100, 199, 100, 199);
  swin = swinopen(win);
  winposition(100, 199, 100, 199);
  
  qdevice(KEYBD);
  qdevice(UPARROWKEY);
  unqdevice(DOWNARROWKEY);
  qdevice(WINQUIT);
  qdevice(REDRAW);
  unqdevice(INPUTCHANGE);

  winset(win);
  loadXfont(4711, "-*-times-medium-r-*-*-*-140-*-*-*-*-iso8859-1");
  getfontencoding(buf2);
#ifdef DEBUG
  fprintf(stderr, "font = %d, encoding = '%s'\n", getfont(), buf2);
#endif
  font(4711);
  getfontencoding(buf2);
#ifdef DEBUG
  fprintf(stderr, "font = %d, encoding = '%s'\n", getfont(), buf2);
#endif
  
  ncolors = 1 << getplanes(); /* Number of colors */
  ncolors = (ncolors<8) ? ncolors : 8;    /* Number of predefined colors */

  while (running) {
    Device dev;
    short val;
    dev = qread(&val);
#ifdef DEBUG
    fprintf(stderr, "dev = %d, *val = %d\n", dev, val);
#endif
    switch(dev) {
    case REDRAW:
      winset(val);
      reshapeviewport();
      drawit(val);
      break;
    case KEYBD:
      switch(val) {
      case 'u':
	unqdevice(UPARROWKEY);
	printf("UPARROWKEY disabled\n");
	break;
      case 'U':
	qdevice(UPARROWKEY);
	printf("UPARROWKEY enabled\n");
	break;
      case 'p':
	printf("Piping main window through 'ppmtogif > Window1.gif'\n");
	winset( win); gl2ppm("| ppmtogif > Window1.gif");
	printf("Piping subwindow through 'ppmtogif > Window2.gif'\n");
	winset(swin); gl2ppm("| ppmtogif > Window2.gif");
	ringbell();
	break;
      case 'q':
      case '\033':
	running = FALSE;
	break;
      }
      break;
    case UPARROWKEY:
      printf("UpArrow %s\n", (val==1) ? "pressed" : "released");
      break;
    case WINQUIT:
      running = FALSE;
      break;
    }
  }
  gexit();
  return(0);
}
