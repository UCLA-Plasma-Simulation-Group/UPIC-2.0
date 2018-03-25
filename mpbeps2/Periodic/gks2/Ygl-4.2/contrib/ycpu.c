/*
 * ycpu.c: display a graph of CPU usage obtained from vmstat
 * in a [Y]gl window. This program may not run on non-AIX.
 * 
 * (C) Fred Hucht 1994,1995, EMail: fred@hal6000.Uni-Duisburg.DE
 * 
 * Feel free to copy and redistribute in terms of the
 * GNU public license.
 */

#include <gl/gl.h>
#include <gl/device.h>
#include <stdio.h>
#include <stdlib.h>

int goon() {
  int r = 1;
  if(qtest()) {
    Int32 dev;
    Int16 val;
    dev = qread(&val);
    switch(dev) {
    case WINQUIT:
      r = 0; break;
    case KEYBD:
      switch(val) {
      case '\033':
      case 'q':
	r = 0; break;
      }
      break;
    }
  }
  return(r);
}

int main(int argc, char *argv[]) {
  FILE *p;
  char hostname[256], cmd[128];
  int X = 200, Y = 200, sec = 1;
  
  putenv("YGL_BS=1");
  
  if(argc > 1) X   = atoi(argv[1]);
  if(argc > 2) sec = atoi(argv[2]);

  prefsize(X+1, Y+1);
  gethostname(hostname, sizeof(hostname));
  winopen(hostname);
  
  qdevice(KEYBD);
  qdevice(WINQUIT);
  
  color(BLACK);
  clear();
  
  sprintf(cmd, "vmstat %d", sec);
  
  if(NULL == (p = (FILE*) popen(cmd, "r"))) {
    perror("pipe error");
    exit(-1);
  }
  
  while(goon()) {
    short us, sy, id, wa;
    short i;
    char line[160];

    if(NULL == fgets(line, 160, p)) {
      perror("fgets error");
      exit(-2);
    }
    if(4 == sscanf(line,
		   "%*d %*d %*d %*d %*d %*d %*d %*d %*d %*d %*d %*d %*d %hd %hd %hd %hd",
		   &us, &sy, &id, &wa)) {
      /*printf("%d %d %d %d\n", us, sy, id, wa);*/
      rectcopy(1, 0, X, Y, 0, 0);
      move2(X, 0);
      color(GREEN);  rdr2(0, Y/100.0 * us);
      color(YELLOW); rdr2(0, Y/100.0 * sy);
      color(BLUE);   rdr2(0, Y/100.0 * id);
      color(RED);    rdr2(0, Y/100.0 * wa);
      for(i=0; i<=10; i++) pnt2(X, Y/100.0 * 10 * i);
      sleep(0);
    }
  }
  gexit();
  pclose(p);
  exit(0);
}
