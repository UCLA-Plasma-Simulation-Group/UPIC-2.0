/* coltest.c by Fred Hucht (C) 1993-2002.
 * Example for color animation with private colormap */

static const char vcid[] = "$Id: coltest.c,v 3.3 1996/07/18 16:35:57 fred Exp $";

#include <X11/Ygl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265359
#endif

void map_hue(Colorindex i, double phi) {
  mapcolor(i,
	   (int)(127.5 * (1.0+cos(phi           ))),
	   (int)(127.5 * (1.0+cos(phi-2./3.*M_PI))),
	   (int)(127.5 * (1.0+cos(phi-4./3.*M_PI)))
	   );
}

int main() {
  int i=0;
  double phi;
  Int16 val;
  Int32 win;
  
  minsize(400,400);
  
  puts("Using private colormap");

  putenv("YGL_PCM=1");
  
  /* setenv("YGL_PCM","1",1); *//* use this ones under BSD unixes */
  
  win = winopen("Coltest, <ESC> to quit");

  qdevice(ESCKEY);
  qdevice(REDRAW);
  qdevice(WINQUIT);
  unqdevice(INPUTCHANGE);

  while(1) {
    phi = 2.0 * M_PI * (double)(++i) / 64.0;
    map_hue(12, phi             );
    map_hue(13, phi + 0.5 * M_PI);
    map_hue(14, phi + 1.0 * M_PI);
    map_hue(15, phi + 1.5 * M_PI);
    
    if(qtest()) switch(qread(&val)) { /* process events */
    case REDRAW:
      reshapeviewport();
      color(WHITE); clear();
      color(12); circfi(100,100,80);
      color(13); circfi(100,300,80);
      color(14); circfi(300,300,80);
      color(15); circfi(300,100,80);
      break;
    case ESCKEY:
    case WINQUIT:
      gexit();
      return(0);
      break;
    }
    usleep(100000); /* Too short sleeps may stress your XServer... */
  }
}
