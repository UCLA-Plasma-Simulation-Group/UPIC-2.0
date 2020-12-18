/* lmbind.c by Fred Hucht (C) 1996-2002.
 * Example program for 3D stuff */

static const char vcid[] = "$Id: lmbind.c,v 3.3 1998-10-28 19:23:01+01 fred Exp fred $";

#include <X11/Ygl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
# define M_PI 3.1415926536
#endif

#define NSPHERE 10

float lmodel[] = {
  TWOSIDE,   0.00,
  LMNULL
};

float material1[] = {
  EMISSION,  0.10, 0.10, 0.10,
  AMBIENT,   0.00, 0.00, 0.00,
  DIFFUSE,   0.80, 0.80, 0.80,
  SPECULAR,  1.00, 1.00, 1.00,
  SHININESS, 10.0,
  ALPHA,     1.0,
  LMNULL
};

float material2[] = {
  EMISSION,  1.00, 1.00, 1.00,
  AMBIENT,   0.00, 0.00, 0.00,
  DIFFUSE,   0.00, 0.00, 0.00,
  SPECULAR,  0.00, 0.00, 0.00,
  LMNULL
};

float light1[] = {
  POSITION,      0.0, 0.0, 1.0, 0.0, /* z +inf */
  /* AMBIENT,       0.0, 0.0, 0.0,
     SPOTDIRECTION, 0.0, 0.0,-1.0,*/
  LCOLOR,        1.0, 0.0, 0.0,
  LMNULL
};

float light2[] = {
  POSITION,      0.0, 1.0, 0.0, 0.0, /* y +inf */
  LCOLOR,        0.0, 1.0, 0.0,
  LMNULL
};

float light3[] = {
  POSITION,      1.0, 0.0, 0.0, 0.0, /* x +inf */
  LCOLOR,        0.0, 0.0, 1.0,
  LMNULL
};

Matrix idmat = {
  {1.0,0.0,0.0,0.0},  /* identity matrix */
  {0.0,1.0,0.0,0.0},
  {0.0,0.0,1.0,0.0},
  {0.0,0.0,0.0,1.0}
};

int main(int argc, char *argv[]) {
  int i;
  Int32 dev;
  Int16 val;
  int mx = 0, my = 0, mz = 0, smi = 0, bf = 0, ff = 0, zb = 1;
  int l1 = 1, l2 = 2, l3 = 3;
  int sm[] = {FLAT, GOURAUD};
  Int32 Kugel;
  double d = 0.0;
  int nsphere = NSPHERE;
  
  if (argc > 1) {
    nsphere = atoi(argv[1]);
  }
  
  minsize(400,400);
  
  putenv("YGL_OGL=1");
  /* setenv("YGL_OGL","1",1); *//* use this ones under BSD unixes */
  
  winopen("lmbind");
  RGBmode();
  doublebuffer();
  gconfig();
  
  qdevice(MOUSEX);
  qdevice(MOUSEY);
  qdevice(KEYBD);
  qdevice(UPARROWKEY);
  qdevice(DOWNARROWKEY);  
  
  Kugel = genobj();
  makeobj(Kugel);

  bgntmesh();
  for (i = nsphere / 6; i < nsphere; i++) {
    double theta  = (M_PI / nsphere) * i;
    double theta2 = (M_PI / nsphere) * (i+1);
    int j;
    
    for (j = 0; j <= 2 * nsphere; j++) {
      double phi = (M_PI / nsphere) * j;
      Float32 v[3], n[3];
      n[0] = sin(theta) * cos(phi);
      n[1] = sin(theta) * sin(phi);
      n[2] = cos(theta);
      v[0] = n[0];
      v[1] = n[1];
      v[2] = n[2];
      n3f(n);
      v3f(v);
      
      n[0] = sin(theta2) * cos(phi);
      n[1] = sin(theta2) * sin(phi);
      n[2] = cos(theta2);
      v[0] = n[0];
      v[1] = n[1];
      v[2] = n[2];
      n3f(n);
      v3f(v);
    }
  }
  endtmesh();
  closeobj();
  
  lsetdepth(0, 0x7FFFFF);
  zbuffer(zb);
#if 1
  mmode(MVIEWING);
  perspective(600, 1.0, 1.0, 30.0);
  
  loadmatrix(idmat); /* redundant here */
  lookat(3.0, 3.0, 3.0,
	 0.0, 0.0, 0.0,
	 0);
#else
  ortho(-2.0,2.0,-2.0,2.0,-2.0,2.0);
  /*  lookat(0,0,0, 3,3,2, 0);*/
#endif
  
  lmdef(DEFLMODEL,   1, sizeof(lmodel)/sizeof(lmodel[0]), lmodel);
  
  lmdef(DEFMATERIAL, 1, sizeof(material1)/sizeof(material1[0]), material1);
  lmdef(DEFMATERIAL, 2, sizeof(material2)/sizeof(material2[0]), material2);
  
  /* lmdef(DEFLIGHT,    1,  0, NULL); */
  lmdef(DEFLIGHT, 1, sizeof(light1)/sizeof(light1[0]), light1);
  lmdef(DEFLIGHT, 2, sizeof(light2)/sizeof(light2[0]), light2);
  lmdef(DEFLIGHT, 3, sizeof(light3)/sizeof(light3[0]), light3);
  
  lmbind(LMODEL, 1);
  lmbind(MATERIAL, 1);
  lmbind(BACKMATERIAL, 2);
  lmbind(LIGHT3, l1);
  lmbind(LIGHT1, l2);
  lmbind(LIGHT2, l3);
  
  while ((dev = qread(&val))) {
    d = 0;
    switch (dev) {
    case MOUSEX: mx = val; break;
    case MOUSEY: my = val; break;
    case KEYBD:
      switch (val) {
      case '\033':
      case 'q': exit(0);
      case '1': l1 = 1 - l1; lmbind(LIGHT3, l1); break;
      case '2': l2 = 2 - l2; lmbind(LIGHT1, l2); break;
      case '3': l3 = 3 - l3; lmbind(LIGHT2, l3); break;
      case 'x': mx += 5; break;
      case 'y': my += 5; break;
      case 'z': mz += 5; break;
      case 'X': mx -= 5; break;
      case 'Y': my -= 5; break;
      case 'Z': mz -= 5; break;
      case 'r': mx = my = mz = 0; break;
      case 'g': shademodel(sm[smi]); smi = 1 - smi; break;
      case 'b': bf = 1 - bf;  backface(bf); break;
      case 'f': ff = 1 - ff; frontface(ff); break;
      case 'u': zb = 1 - zb; zbuffer(zb); break;
      case 'p': gl2ppm("| ppmtogif > lmbind.gif"); break;
      }
      break;
    case UPARROWKEY:   if (val) d =  0.5; break;
    case DOWNARROWKEY: if (val) d = -0.5; break;
    }
    
    /* No qreset(), drops important events */
    
    /* Move eye */
    mmode(MPROJECTION);
    translate(0.0, 0.0, d);
    mmode(MVIEWING);
    
    pushmatrix();
    
    rotate(30 * mx, 'x');
    rotate(30 * my, 'y');
    rotate(30 * mz, 'z');
    
    RGBcolor(0,0,0);
    clear();
    zclear();
    
    callobj(Kugel);
#if 1
    translate(0.0, 0.0, 2.0);
    callobj(Kugel);
    
    translate(0.0, 2.0,-2.0);
    callobj(Kugel);
    
    translate(2.0,-2.0, 0.0);
    callobj(Kugel);
#endif
    popmatrix();
    swapbuffers();
  }
  return 0;
}
