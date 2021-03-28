/* smile.c by Fred Hucht (C) 1993-2002.
 * Short and simple example program */

static const char vcid[] = "$Id: smile.c,v 3.2 1996/07/18 16:35:57 fred Exp $";

#include <X11/Ygl.h>
#include <stdio.h>
#include <unistd.h>

int main() {
  Int32 win;
  
  prefsize(100, 100);
  win = winopen("Smile!");

  /*printf("%d\n", win);*/

  /* The background */
  color(BLACK); clear();

  /* The face */
  color(RED); circfi(50, 50, 40);
  
  /* The eyes */
  color(WHITE); circfi(30, 60, 10); circfi(70, 60, 10);
  color(BLACK); circfi(30, 60,  5); circfi(70, 60,  5);
  
  /* The smile */
  arci(50, 50, 25, 2000, 3400);
  
  /* The sleep */
  sleep(2);

  /* The twinkle */
  color(RED);circfi(30, 65, 10);
  sleep(1);
  color(WHITE); circfi(30, 60, 10);
  color(BLACK); circfi(30, 60,  5);

  /* The end */
  sleep(2);
  gexit();
  return(0);
}
