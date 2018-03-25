/*
 *    Ygl: Run GL programs with standard X11 and/or OpenGL routines.
 *    (C) Fred Hucht 1996-2007
 *    EMail: fred<at>thp.Uni-Duisburg.de
 */

static const char vcid[] = "$Id: 3d.c,v 4.9 2007-05-08 13:28:43+02 fred Exp $";

#include "header.h"

#ifdef OGL
static int msingle = 1; /* IrisGL is in single matrix mode per default */
#endif

/*
 * Projection transformations (apply to projection matrix stack)
 */
void ortho(Coord left, Coord right, Coord bottom,
	   Coord top,  Coord near,  Coord far) {
  const char * MyName = "ortho";
  I(MyName, "%g,%g,%g,%g,%g,%g", left, right, bottom, top, near, far);
  IFOGL(
	GLint old;
	glGetIntegerv(GL_MATRIX_MODE, &old);	/* save old mmode */
	glMatrixMode(GL_PROJECTION);		/*switch to PROJECTION mode*/
	if (!Ygl.SelectMode) glLoadIdentity();	/* reset matrix */
	glOrtho(left, right, bottom, top, near, far);
	if (msingle || Ygl.SelectMode) {
	  glMatrixMode(GL_MODELVIEW);
	  glLoadIdentity();
	}  
	glMatrixMode(old),			/* restore mmode */
	/* In X11 2d mode we allow ortho() calls and ignore the z-component */
	ortho2(left, right, bottom, top)
	);	
}

#ifdef OGL
/* the remaining stuff is only available in OpenGL mode */

void perspective(Angle ang, Float32 aspect, Coord near, Coord far) {
  const char * MyName = "perspective";
  GLint old;
  I(MyName, "%d,%g,%g,%g", ang, aspect, near, far);
  if(!Ygl.UseOGL) NI(MyName);
  
  glGetIntegerv(GL_MATRIX_MODE, &old);
  glMatrixMode(GL_PROJECTION);
  if (!Ygl.SelectMode) glLoadIdentity();
  gluPerspective(0.1 * ang, aspect, near, far);
  if (msingle || Ygl.SelectMode) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  }  
  glMatrixMode(old);
}

void window(Coord left, Coord right, Coord bottom,
	    Coord top,  Coord near,  Coord far) {
  const char * MyName = "window";
  GLint old;
  I(MyName, "%g,%g,%g,%g,%g,%g", left, right, bottom, top, near, far);
  if(!Ygl.UseOGL) NI(MyName);
  
  glGetIntegerv(GL_MATRIX_MODE, &old);
  glMatrixMode(GL_PROJECTION);
  if (!Ygl.SelectMode) glLoadIdentity();
  glFrustum(left, right, bottom, top, near, far);
  if (msingle || Ygl.SelectMode) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  }  
  glMatrixMode(old);
}

/*
 * Coordinate transformations 
 */
void lookat(Coord eyex, Coord eyey, Coord eyez,
	    Coord cenx, Coord ceny, Coord cenz,
	    Angle twist) {
  const char * MyName = "lookat";
  GLdouble rx = cenx - eyex;
  GLdouble ry = ceny - eyey;
  GLdouble rz = cenz - eyez;
  GLdouble upx, upy, upz;
  I(MyName, "%g,%g,%g,%g,%g,%g,%d", eyex, eyey, eyez, cenx, ceny, cenz, twist);
  if(!Ygl.UseOGL) NI(MyName);
  if(rx != 0.0 || rz != 0.0) {
    /* OK, view line is not y-axis, so up direction is y-axis */
    upx = 0.0; upy = 1.0; upz = 0.0;
  } else {
    /* Special case: looking along y-axis, up = -z */
    upx = 0.0; upy = 0.0; upz = -1.0;
  }
  gluLookAt(eyex, eyey, eyez, cenx, ceny, cenz, upx, upy, upz);
  glRotatef(0.1 * twist, rx, ry, rz);
}

void polarview(Coord distance, Angle azim, Angle inc, Angle twist) {
  const char * MyName = "polarview";
  I(MyName, "%g,%d,%d,%d", distance, azim, inc, twist);
  if(!Ygl.UseOGL) NI(MyName);
  glTranslatef(0.0, 0.0, -distance);
  glRotatef(-0.1 * twist, 0.0, 0.0, 1.0);
  glRotatef(-0.1 * inc,   1.0, 0.0, 0.0);
  glRotatef(-0.1 * azim,  0.0, 0.0, 1.0);
}

void rotate(Angle angle, Char8 axis) {
  const char * MyName = "rotate";
  I(MyName, "%d,'%c'", angle, axis);
  if(!Ygl.UseOGL) NI(MyName);
  rot(0.1 * angle, axis);
}

void rot(Float32 angle, Char8 axis) {
  const char * MyName = "rot";
  GLfloat x = 0.0F, y = 0.0F, z = 0.0F;
  I(MyName, "%g,'%c'", angle, axis);
  if(!Ygl.UseOGL) NI(MyName);
  switch(axis) {
  case 'x': case 'X': x = 1.0F; break;
  case 'y': case 'Y': y = 1.0F; break;
  case 'z': case 'Z': z = 1.0F; break;
  default:
    Yprintf(MyName, "invalid axis '%c'.\n", axis);
    return;
  }
  glRotatef(angle, x, y, z);
}

void translate(Coord x, Coord y, Coord z) {
  const char * MyName = "translate";
  I(MyName, "%g,%g,%g", x, y, z);
  if(!Ygl.UseOGL) NI(MyName);
  glTranslatef(x, y, z);
}

void scale(Float32 x, Float32 y, Float32 z) {
  const char * MyName = "scale";
  I(MyName, "%g,%g,%g", x, y, z);
  if(!Ygl.UseOGL) NI(MyName);
  glScalef(x, y, z);
#if 0
  /* done in setup_visuals */
  if(x != 1.0 || y != 1.0 || z != 1.0) {
    glEnable(GL_NORMALIZE); /* Could be better... */
  }
#endif
}

/*
 * Matrix ops
 */
void pushmatrix(void) {
  const char * MyName = "pushmatrix";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  glPushMatrix();
}

void popmatrix(void) {
  const char * MyName = "popmatrix";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  glPopMatrix();
}

void loadmatrix(Matrix m) {
  const char * MyName = "loadmatrix";
  GLfloat glm[16];
  short i, j;
  I(MyName, "Matrix");
  if(!Ygl.UseOGL) NI(MyName);
  /* Under OpenGL matrices are in column order, transpose... */
  for(i = 0; i < 4; i++) for(j = 0; j < 4; j++) glm[4*j + i] = m[i][j];
  glLoadMatrixf(glm);
}

/* getmatrix: see misc.c */

void multmatrix(Matrix m) {
  const char * MyName = "multmatrix";
  GLfloat glm[16];
  short i, j;
  I(MyName, "Matrix");
  if(!Ygl.UseOGL) NI(MyName);
  /* Under OpenGL matrices are in column order, transpose... */
  for(i = 0; i < 4; i++) for(j = 0; j < 4; j++) glm[4*j + i] = m[i][j];
  glMultMatrixf(glm);
}

void mmode(Int16 mode) {
  const char * MyName = "mmode";
  I(MyName, "%d", mode);
  if(!Ygl.UseOGL) NI(MyName);
  switch(mode) {
  case MSINGLE:
    msingle = 1;
    glMatrixMode(GL_MODELVIEW);
    break;
  case MPROJECTION:
    msingle = 0;
    glMatrixMode(GL_PROJECTION);
    break;
  case MVIEWING:
    msingle = 0;
    glMatrixMode(GL_MODELVIEW);
    break;
  default:
    Yprintf(MyName, "invalid mode %d.\n", mode);
    break;
  }
}

Int32 getmmode(void) {
  const char * MyName = "getmmode";
  GLint mode;
  I(MyName, "");
  if(!Ygl.UseOGL) NIR(MyName, -1);
  
  glGetIntegerv(GL_MATRIX_MODE, &mode);
  
  switch(mode) {
  case GL_MODELVIEW:  return MVIEWING;
  case GL_PROJECTION: return MPROJECTION;
  default:
    Yprintf(MyName, "strange OpenGL matrix mode %d\n", mode);
    return -1;
  }
}

/*
 * z-buffer stuff
 */
void depthcue(Int32 mode) {
  const char * MyName = "depthcue";
  I(MyName, "%d", mode);
  if(!Ygl.UseOGL) NI(MyName);
  if(!Ygl.PCM && !W->rgb) {
    Yprintf(MyName, "private colormap required.\n");
    return;
  }
  switch (mode) {
  case True:
    glEnable (GL_FOG);
    break;
  case False:
    glDisable(GL_FOG);
    break;
  default:
    Yprintf(MyName, "invalid mode %d.\n", mode);
    return;
  }
}

#define GL2OGL_Z(z) (((double)(z) + 0x800000)/0xFFFFFF)

void lshaderange(Colorindex low, Colorindex high, Int32 near, Int32 far) {
  const char * MyName = "lshaderange";
  I(MyName, "%d,%d,%d,%d", low, high, near, far);
  if(!Ygl.UseOGL) NI(MyName);
  glFogi(GL_FOG_MODE, GL_LINEAR);
  glFogi(GL_FOG_START, near);
  glFogi(GL_FOG_END,   far);
  glFogi(GL_FOG_INDEX, YGL_COLORS(low));
}
  
void lRGBrange(Int16 rmin, Int16 gmin, Int16 bmin,
	       Int16 rmax, Int16 gmax, Int16 bmax, 
	       Int32 near, Int32 far) {
  const char * MyName = "lRGBrange";
  GLfloat color[3];
  I(MyName, "%d,%d,%d,%d,%d,%d,%d,%d", rmin, gmin, bmin, rmax, gmax, bmax, near, far);
  if(!Ygl.UseOGL) NI(MyName);
  color[0] = rmax;
  color[1] = gmax;
  color[2] = bmax;
  glFogi(GL_FOG_MODE, GL_LINEAR);
  glFogi(GL_FOG_START, near);
  glFogi(GL_FOG_END,   far);
  glFogfv(GL_FOG_COLOR, color);
}

void zbuffer(Int32 mode) {
  const char * MyName = "zbuffer";
  I(MyName, "%d", mode);
  if(!Ygl.UseOGL) NI(MyName);
  switch (mode) {
  case True:
    glEnable (GL_DEPTH_TEST);
    break;
  case False:
    glDisable(GL_DEPTH_TEST);
    break;
  default:
    Yprintf(MyName, "invalid mode %d.\n", mode);
    return;
  }
}

void lsetdepth(Int32 near, Int32 far) {
  const char * MyName = "lsetdepth";
  I(MyName, "%d,%d", near, far);
  if(!Ygl.UseOGL) NI(MyName);
  glDepthRange(GL2OGL_Z(near),
	       GL2OGL_Z(far));
}

void zclear(void) {
  const char * MyName = "zclear";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  /*glClearDepth(1.0);*/
  glClear(GL_DEPTH_BUFFER_BIT);
}

void zdraw(Int32 bool) {
  const char * MyName = "zdraw";
  I(MyName, "%d", bool);
  if(!Ygl.UseOGL) NI(MyName);
  glDepthMask(bool);
}

void czclear(Int32 cval, Int32 zval) {
  const char * MyName = "czclear";
  I(MyName, "%d,%d", cval, zval);
  if(!Ygl.UseOGL) NI(MyName);
  glClearDepth(GL2OGL_Z(zval));
  /*if(Ygl.RGB) glClearColor(glColor4ubv((GLubyte *)&cval);
    else */
  glClearIndex(cval);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void zfunction(Int32 func) {
  const char * MyName = "zfunction";
  GLenum glfunc;
  I(MyName, "%d", func);
  if(!Ygl.UseOGL) NI(MyName);
  switch(func) {
  case ZF_NEVER:    glfunc = GL_NEVER;    break;
  case ZF_LESS:     glfunc = GL_LESS;     break;
  case ZF_LEQUAL:   glfunc = GL_LEQUAL;   break;
  case ZF_GREATER:  glfunc = GL_GREATER;  break;
  case ZF_NOTEQUAL: glfunc = GL_NOTEQUAL; break;
  case ZF_GEQUAL:   glfunc = GL_GEQUAL;   break;
  case ZF_ALWAYS:   glfunc = GL_ALWAYS;   break;
  default:
    Yprintf(MyName, "invalid function %d.\n", func);
    return;
  }
  glDepthFunc(glfunc);
}

/* Display lists */
Int32 genobj(void) {
  const char * MyName = "genobj";
  I(MyName, "");
  if(!Ygl.UseOGL) NIR(MyName, -1);
  return glGenLists(1);
}

Int32 isobj(Int32 object) {
  const char * MyName = "isobj";
  I(MyName, "%d", object);
  if(!Ygl.UseOGL) NIR(MyName, -1);
  return glIsList(object);
}

void makeobj(Int32 object) {
  const char * MyName = "makeobj";
  I(MyName, "%d", object);
  if(!Ygl.UseOGL) NI(MyName);
  glNewList(object, GL_COMPILE);
}

Int32 getopenobj(void) {
  const char * MyName = "getopenobj";
  GLint object;
  I(MyName, "");
  if(!Ygl.UseOGL) NIR(MyName, -1);
  glGetIntegerv(GL_LIST_INDEX, &object);
  return object;
}

void closeobj(void) {
  const char * MyName = "closeobj";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  glEndList();
}

void callobj(Int32 object) {
  const char * MyName = "callobj";
  I(MyName, "%d", object);
  if(!Ygl.UseOGL) NI(MyName);
  glCallList(object);
}

void delobj(Int32 object) {
  const char * MyName = "delobj";
  I(MyName, "%d", object);
  if(!Ygl.UseOGL) NI(MyName);
  glDeleteLists(object, 1);
}

/* Lighting */

/*
 * Material
 */
static int MaterialNum = 0;
static struct Material {
  Uint    index;		/* lmdef index */
  GLuint  list;			/* 1st of 2 OpenGL List IDs for MATERIAL */
  Float32 emission[4];		/* and BACKMATERIAL */
  Float32 ambient[4];
  Float32 diffuse[4];
  Float32 specular[4];
  Float32 shininess[1];
  Float32 alpha[1];
  Float32 colorindexes[4];
} *Materials = NULL, Materialdefault = {
  0, 0,
  {0.0F, 0.0F, 0.0F, 1.0F},	/* emission */
  {0.2F, 0.2F, 0.2F, 1.0F},	/* ambient */
  {0.8F, 0.8F, 0.8F, 1.0F},	/* diffuse */
  {0.0F, 0.0F, 0.0F, 1.0F},	/* specular */
  {0.0F},			/* shininess */
  {1.0F},			/* alpha */
  {0.0F, 127.5F, 255.0F}};	/* colorindexes */

static void defmaterial(const char *caller, 
			Int16 index, Int16 numpoints, Float32 prop[]) {
  int i;
  int newmaterial = False;
  struct Material *material;
  
  if(Materials == NULL) { /* initialize */
    i = MaterialNum++;
    Materials = (struct Material*)malloc(sizeof(struct Material));
    newmaterial = True;
  } else {
    for(i = MaterialNum - 1; i >= 0 && Materials[i].index != index; i--);
    if(i < 0) { /* not found, allocate new */
      i = MaterialNum++;
      Materials = (struct Material*)
	realloc(Materials, MaterialNum * sizeof(struct Material));
      newmaterial = True;
    }
  }
  
  if(Materials == NULL) {
    Yprintf(caller, "can't allocate memory.\n");
    exit(-1);
  }
  
  material = &Materials[i];
  
  if(newmaterial) { /* initialize with defaults */
    memcpy(material, &Materialdefault, sizeof(struct Material));
    material->index = index;
    material->list  = glGenLists(2);
  } else {
    glDeleteLists(material->list, 2);
  }
  
  for(i = 0; i < numpoints && prop[i];) switch((int)(prop[i++])) {
    int j;
  case EMISSION:
    for(j = 0; j < 3; j++) material->emission[j] = prop[i++];
    break;
  case AMBIENT:
    for(j = 0; j < 3; j++) material->ambient[j] = prop[i++];
    break;
  case DIFFUSE:
    for(j = 0; j < 3; j++) material->diffuse[j] = prop[i++];
    break;
  case SPECULAR:
    for(j = 0; j < 3; j++) material->specular[j] = prop[i++];
    break;
  case SHININESS:
    material->shininess[0] = prop[i++];
    break;
  case ALPHA:
    material->alpha[0] = prop[i++];
    break;
  case COLORINDEXES:
    for(j = 0; j < 3; j++) material->colorindexes[j] = prop[i++];
    break;
  default:
    Yprintf(caller, "invalid MATERIAL property %g.\n", prop[i-1]);
    return;
  }
  /* Set alpha components */
  material->emission[3] = material->alpha[0];
  material->ambient[3]  = material->alpha[0];
  material->diffuse[3]  = material->alpha[0];
  material->specular[3] = material->alpha[0];
  
  /* Put the complete definition into a display list */
  glNewList(material->list, GL_COMPILE);
  glMaterialfv(GL_FRONT, GL_EMISSION,  material->emission);
  glMaterialfv(GL_FRONT, GL_AMBIENT,   material->ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   material->diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR,  material->specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, material->shininess);
  glMaterialfv(GL_FRONT, GL_COLOR_INDEXES, material->colorindexes);
  glEndList();

  /* and again for BACKMATERIAL */
  glNewList(material->list + 1, GL_COMPILE);
  glMaterialfv(GL_BACK, GL_EMISSION,  material->emission);
  glMaterialfv(GL_BACK, GL_AMBIENT,   material->ambient);
  glMaterialfv(GL_BACK, GL_DIFFUSE,   material->diffuse);
  glMaterialfv(GL_BACK, GL_SPECULAR,  material->specular);
  glMaterialfv(GL_BACK, GL_SHININESS, material->shininess);
  glMaterialfv(GL_BACK, GL_COLOR_INDEXES, material->colorindexes);
  glEndList();
}

/*
 * LightModel
 */
static int LightModelNum = 0;
static struct LightModel {
  Uint    index;		/* lmdef index */
  GLuint  list;			/* OpenGL List ID */
  Float32 ambient[4];
  Float32 localviewer[1];
  Float32 twoside[1];
} *LightModels = NULL, LightModeldefault = {
  0, 0,
  {0.2F, 0.2F, 0.2F, 1.0F},	/* ambient */
  {0.0F},			/* localviewer */
  {0.0F}			/* twoside */
};

static struct Attenuation_ {
  Float32 constant;
  Float32 linear;
  Float32 quadratic;
} Attenuation = {
  1.0F,
  0.0F,
  0.0F
};
  
static void deflightmodel(const char *caller, 
			  Int16 index, Int16 numpoints, Float32 prop[]) {
  int i;
  int newlightmodel = False;
  struct LightModel *lightmodel;
  
  if(LightModels == NULL) { /* initialize */
    i = LightModelNum++;
    LightModels = (struct LightModel*)malloc(sizeof(struct LightModel));
    newlightmodel = True;
  } else {
    for(i = LightModelNum - 1; i >= 0 && LightModels[i].index != index; i--);
    if(i < 0) { /* not found, allocate new */
      i = LightModelNum++;
      LightModels = (struct LightModel*)
	realloc(LightModels, LightModelNum * sizeof(struct LightModel));
      newlightmodel = True;
    }
  }
  
  if(LightModels == NULL) {
    Yprintf(caller, "can't allocate memory.\n");
    exit(-1);
  }
  
  lightmodel = &LightModels[i];
  
  if(newlightmodel) { /* initialize with defaults */
    memcpy(lightmodel, &LightModeldefault, sizeof(struct LightModel));
    lightmodel->index = index;
    lightmodel->list  = glGenLists(1);
  } else {
    glDeleteLists(lightmodel->list, 1);
  }
  
  for(i = 0; i < numpoints && prop[i];) switch((int)(prop[i++])) {
    int j;
  case AMBIENT:
    for(j = 0; j < 3; j++) lightmodel->ambient[j] = prop[i++];
    break;
  case LOCALVIEWER:
    lightmodel->localviewer[0] = prop[i++];
    break;
  case ATTENUATION:
    /* used in lmbind */
    Attenuation.constant = prop[i++];
    Attenuation.linear   = prop[i++];
    /*Yprintf(caller, "LMODEL ATTENUATION experimental.\n");*/
    break;
  case ATTENUATION2:
    /* used in lmbind */
    Attenuation.quadratic = prop[i++];
    /*Yprintf(caller, "LMODEL ATTENUATION2 experimental.\n");*/
    break;
  case TWOSIDE:
    lightmodel->twoside[0] = prop[i++];
    break;
  default:
    Yprintf(caller, "invalid LMODEL property %g.\n", prop[i-1]);
    return;
  }
  
  /* Put the complete definition into a display list */
  glNewList(lightmodel->list, GL_COMPILE);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,
		 lightmodel->ambient);
  glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER,
		 lightmodel->localviewer);
  glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE,
		 lightmodel->twoside);
  glEndList();
}

/*
 * Light
 */
static int LightNum = 0;
static struct Light {
  Uint    index;
  Float32 ambient[4];
  Float32 lcolor[4];
  Float32 position[4];
  Float32 spotdir[3];
  Float32 spotexp[1];
  Float32 spotcut[1];
} *Lights = NULL, Lightdefault = {
  0,
  {0.0F, 0.0F, 0.0F, 1.0F},	/* ambient */
  {1.0F, 1.0F, 1.0F, 1.0F},	/* lcolor */
  {0.0F, 0.0F, 1.0F, 0.0F},	/* position */
  {0.0F, 0.0F,-1.0F},		/* spotdir */
  {0.0F},			/* spotexp */
  {180.0F}};			/* spotcut */

/* Which definition is already bound to light i? */
static Int32 boundLightindex[8];

static void deflight(const char *caller, 
		     Int16 index, Int16 numpoints, Float32 prop[]) {
  int i;
  int newlight = False;
  struct Light *light;
  
  /* Reset entry in boundLightindex so lmbind reload definition */
  for(i = 0; i < 8; i++) if(boundLightindex[i] == index)
    boundLightindex[i] = 0;
  
  if(Lights == NULL) { /* initialize */
    i = LightNum++;
    Lights = (struct Light*)malloc(sizeof(struct Light));
    newlight = True;
    
  } else {
    for(i = LightNum - 1; i >= 0 && Lights[i].index != index; i--);
    if(i < 0) { /* not found, allocate new */
      i = LightNum++;
      Lights = (struct Light*)
	realloc(Lights, LightNum * sizeof(struct Light));
      newlight = True;
    }
  }
  
  if(Lights == NULL) {
    Yprintf(caller, "can't allocate memory.\n");
    exit(-1);
  }
  
  light = &Lights[i];
  
  if(newlight) { /* initialize with defaults */
    memcpy(light, &Lightdefault, sizeof(struct Light));
    light->index = index;
  }
  
  for(i = 0; i < numpoints && prop[i];) switch((int)(prop[i++])) {
    int j;
  case AMBIENT:
    for(j = 0; j < 3; j++) light->ambient[j] = prop[i++];
    break;
  case LCOLOR:
    for(j = 0; j < 3; j++) light->lcolor[j] = prop[i++];
    break;
  case POSITION:
    for(j = 0; j < 4; j++) light->position[j] = prop[i++];
    break;
  case SPOTDIRECTION:
    for(j = 0; j < 3; j++) light->spotdir[j] = prop[i++];
    break;
  case SPOTLIGHT:
    light->spotexp[0] = prop[i++];
    light->spotcut[0] = prop[i++];
    break;
  default:
    Yprintf(caller, "invalid LIGHT property %g.\n", prop[i-1]);
    return;
  }
}

void lmdef(Int16 deftype, Int16 index, Int16 numpoints, Float32 prop[]) {
  const char * MyName = "lmdef";
  I(MyName, "%d,%d,%d,*", deftype, index, numpoints);
  if(!Ygl.UseOGL) NI(MyName);
  
#ifdef DEBUG
  {
    int i;
    Yprintf(MyName, "deftype=%d, index=%d, numpoints=%d: ",
	    deftype, index, numpoints);
    for(i = 0; i < numpoints; i++) fprintf(stderr, "%g ", prop[i]);
    fprintf(stderr, "\n");
  }
#endif
  
  if(index == 0) {
    Yprintf(MyName, "can't redefine index 0.\n");
    return;
  }
  
  switch(deftype) {
  case DEFMATERIAL:
    defmaterial(MyName, index, numpoints, prop);
    break;
  case DEFLIGHT:
    deflight(MyName, index, numpoints, prop);
    break;
  case DEFLMODEL:
    deflightmodel(MyName, index, numpoints, prop);
    break;
  default:
    Yprintf(MyName, "invalid type %d.\n", deftype);
    return;
  }
}

void lmbind(Int32 target, Int32 index) {
  const char * MyName = "lmbind";
  GLenum gllight;
  int i, boundi;
  I(MyName, "%d,%d", target, index);
  if(!Ygl.UseOGL) NI(MyName);
  
  if(target == MATERIAL || target == LMODEL) {
    if(index == 0) {
      glDisable(GL_LIGHTING);
      glDisable(GL_NORMALIZE);
      return;
    } else {
      glEnable(GL_LIGHTING);
      glEnable(GL_NORMALIZE); /* GL sets this */
    }
  }
  
  switch(target) {
  case MATERIAL:
    for(i = MaterialNum - 1; i >= 0 && Materials[i].index != index; i--);
    if(i < 0) {
      Yprintf(MyName, "MATERIAL %d not defined.\n", index);
      return;
    }
    glCallList(Materials[i].list);
    return;
  case BACKMATERIAL:
    for(i = MaterialNum - 1; i >= 0 && Materials[i].index != index; i--);
    if(i < 0) {
      Yprintf(MyName, "MATERIAL %d not defined.\n", index);
      return;
    }
    glCallList(Materials[i].list + 1);
    return;
  case LMODEL:
    for(i = LightModelNum - 1; i >= 0 && LightModels[i].index != index; i--);
    if(i < 0) {
      Yprintf(MyName, "MATERIAL %d not defined.\n", index);
      return;
    }
    glCallList(LightModels[i].list);
    return;
  case LIGHT0: boundi = 0; gllight = GL_LIGHT0; break;
  case LIGHT1: boundi = 1; gllight = GL_LIGHT1; break;
  case LIGHT2: boundi = 2; gllight = GL_LIGHT2; break;
  case LIGHT3: boundi = 3; gllight = GL_LIGHT3; break;
  case LIGHT4: boundi = 4; gllight = GL_LIGHT4; break;
  case LIGHT5: boundi = 5; gllight = GL_LIGHT5; break;
  case LIGHT6: boundi = 6; gllight = GL_LIGHT6; break;
  case LIGHT7: boundi = 7; gllight = GL_LIGHT7; break;
  default:
    Yprintf(MyName, "invalid target %d.\n", target);
    return;
  }
  
  if(index == 0) {
    glDisable(gllight);
  } else {
    glEnable(gllight);
    if(index != boundLightindex[boundi]) {
      /* load definition for light if changed or not already loaded */
      boundLightindex[boundi] = index;
      for(i = LightNum - 1; i >= 0 && Lights[i].index != index; i--);
      if(i < 0) {
	Yprintf(MyName, "LIGHT %d not defined.\n", index);
	return;
      }
      glLightfv(gllight, GL_AMBIENT,        Lights[i].ambient);
      glLightfv(gllight, GL_DIFFUSE,        Lights[i].lcolor);
      glLightfv(gllight, GL_SPECULAR,       Lights[i].lcolor);
      glLightfv(gllight, GL_POSITION,       Lights[i].position);
      glLightfv(gllight, GL_SPOT_DIRECTION, Lights[i].spotdir);
      glLightfv(gllight, GL_SPOT_EXPONENT,  Lights[i].spotexp);
      glLightfv(gllight, GL_SPOT_CUTOFF,    Lights[i].spotcut);
      glLightf(gllight, GL_CONSTANT_ATTENUATION,  Attenuation.constant);
      glLightf(gllight, GL_LINEAR_ATTENUATION,    Attenuation.linear);
      glLightf(gllight, GL_QUADRATIC_ATTENUATION, Attenuation.quadratic);
    }
  }
}

void lmcolor(Int32 mode) {
  const char * MyName = "lmcolor";
  GLenum glmode;
  I(MyName, "%d", mode);
  if(!Ygl.UseOGL) NI(MyName);

  switch(mode) {
  case LMC_COLOR:
    glDisable(GL_COLOR_MATERIAL);
    return;
  case LMC_NULL:
    Yprintf(MyName, "unimplemented mode %d.\n", mode);
    return;
  case LMC_EMISSION: glmode = GL_EMISSION; break;
  case LMC_AMBIENT:  glmode = GL_AMBIENT;  break;
  case LMC_DIFFUSE:  glmode = GL_DIFFUSE;  break;
  case LMC_SPECULAR: glmode = GL_SPECULAR; break;
  case LMC_AD:       glmode = GL_AMBIENT_AND_DIFFUSE; break;
  default:
    Yprintf(MyName, "invalid mode %d.\n", mode);
    return;
  }
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, glmode);
}

void shademodel(Int32 mode) {
  const char * MyName = "shademodel";
  GLenum glmode;
  I(MyName, "%d", mode);
  if(!Ygl.UseOGL) NI(MyName);
  
  switch(mode) {
  case FLAT:    glmode = GL_FLAT;   break;
  case GOURAUD: glmode = GL_SMOOTH; break;
  default:
    Yprintf(MyName, "invalid mode %d.\n", mode);
    return;
  }
  glShadeModel(glmode);
}

static int Face[] = {False, False};

static void front_back(const char *caller, GLenum glcull, Int32 mode) {
  I(caller, "%d", mode);
  if(!Ygl.UseOGL) NI(caller);
  switch(mode) {
  case True:
    Face[0] = True;
    glCullFace(glcull);
    break;
  case False:
    Face[0] = False;
    break;
  default:
    Yprintf(caller, "invalid mode %d.\n", mode);
    return;
  }
  if(Face[0] || Face[1]) {
    glEnable(GL_CULL_FACE);
  } else {
    glDisable(GL_CULL_FACE);
  }
}

void frontface(Int32 mode) {
  const char * MyName = "frontface";
  front_back(MyName, GL_FRONT, mode);
}

void backface(Int32 mode) {
  const char * MyName = "backface";
  front_back(MyName, GL_BACK, mode);
}

void RGBwritemask(Int16 red, Int16 green, Int16 blue) {
  const char * MyName = "RGBwritemask";
  I(MyName, "%d,%d,%d", red, green, blue);
  if(!Ygl.UseOGL) NI(MyName);
  glColorMask(red != 0, green != 0, blue != 0, GL_FALSE);
}

/* for Pete Riley */
void drawmode(Int32 mode) {
  const char * MyName = "drawmode";
  I(MyName, "%d", mode);
  if(!Ygl.UseOGL) NI(MyName);
  NI(MyName);
}

void iconsize(Int32 sx, Int32 sy) {
  const char * MyName = "iconsize";
  I(MyName, "%d,%d", sx, sy);
  if(!Ygl.UseOGL) NI(MyName);
  NI(MyName);
}

void overlay(Int32 arg) {
  const char * MyName = "overlay";
  I(MyName, "%d", arg);
  if(!Ygl.UseOGL) NI(MyName);
  NI(MyName);
}

void pushattributes(void) {
  const char * MyName = "pushattributes";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  glPushAttrib(GL_ALL_ATTRIB_BITS);
}

void popattributes(void) {
  const char * MyName = "popattributes";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  glPopAttrib();
}

void fullscrn(void) {
  const char * MyName = "fullscrn";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  NI(MyName);
}

void endfullscrn(void) {
  const char * MyName = "endfullscrn";
  I(MyName, "");
  if(!Ygl.UseOGL) NI(MyName);
  NI(MyName);
}

void scrmask(Screencoord l, Screencoord r, Screencoord b, Screencoord t) {
  const char * MyName = "scrmask";
  I(MyName, "%d,%d,%d,%d", l, r, b, t);
  if(!Ygl.UseOGL) NI(MyName);
  glScissor(l, b, r - l, t - b);
  glEnable(GL_SCISSOR_TEST);
}

#endif /* OGL */
