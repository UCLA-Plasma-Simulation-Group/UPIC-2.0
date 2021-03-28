/*
 * yup.c: display a graph of CPU uptime in a [Y]gl window.
 * This program may not run on non-AIX.
 * 
 * (C) Fred Hucht 1994, EMail: fred@hal6000.Uni-Duisburg.DE
 * 
 * Feel free to copy and redistribute in terms of the
 * GNU public license.
 */

#include <gl/gl.h>
#include <stdio.h>

#define FACT 2
#define X 600
#define Y 200

main(int argc, char *argv[]) {
  FILE *p;
  char cmd[80];
  
  sprintf(cmd, "rup %s", argc == 1 ? "localhost" : argv[1]);
  
  putenv("YGL_BS=1");
  prefsize(X+1, Y+1);
  winopen(cmd);
  color(BLACK);
  clear();
  
  while(1) {
    double a1,a2,a3;
    short i;
    char x[80];
    p = (FILE*) popen(cmd, "r");
    fscanf(p, "%51c %lg%*c%lg%*c%lg", x, &a1, &a2, &a3);
    pclose(p);
    /*    printf("%s-%g %g %g\n", x, a1,a2,a3); */
    rectcopy(1, 0, X, Y, 0, 0);
    color(BLACK); move2(X, 0);draw2(X, Y);
    color(YELLOW);
    for(i=0; i<=10; i++) pnt2(X, 10*i*FACT);
    color(BLUE);  pnt2(X, a3*FACT);
    color(GREEN); pnt2(X, a2*FACT);
    color(RED);   pnt2(X, a1*FACT);
    sleep(5);
  }
}
