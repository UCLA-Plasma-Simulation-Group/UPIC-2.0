/*
 *    Ygl: Run GL programs with standard X11 and/or OpenGL routines.
 *    (C) Fred Hucht 1993-2006
 *    EMail: fred<at>thp.Uni-Duisburg.de
 */

static const char vcid[] = "$Id: fortran.c,v 4.2 2007-05-08 13:23:57+02 fred Exp $";

#include "header.h" /* for COVERSLEEP & string.h */

Display *getXdpy(void);
Window getXwid(void);
Window getXdid(void);
GC getXgc(void);
Int32 winX(Display *dpy, Window win);

#ifndef YGL_PRE
# define YGL_PRE(x) x
#endif

/* Fortran types */
#define REAL  	Float32
#define REAL8 	Float64
#define INT2  	Int16
#define INT4  	Int32
#define LOGICAL	Int32
#define CHAR	char

/* Append '\0' to FORTRAN string.
 * Thx to Frank Scott <fscott@frasco.demon.co.uk> for reporting bug and fix */
#define L_TBUF 255
static char Tbuf[L_TBUF+1];
static char *apnd_0(CHAR * str, INT4 * len) {
  int l = *len > L_TBUF ? L_TBUF : *len;
  strncpy(Tbuf, str, l);
  Tbuf[l] = '\0';
  return Tbuf;
}

void setval_(INT4*x1, INT4*x2, INT4*x3, INT4*x4) {}
void curson_(void) {}
void cursof_(void) {}
void chunks_(void) {}
void greset_(void) {}

/********************* draw.c */
void YGL_PRE(clear_) (void) { clear();}

/* Points */
void YGL_PRE(pnt2_)  (REAL*x1, REAL*x2) { pnt2 (*x1, *x2);}
void YGL_PRE(pnt2i_) (INT4*x1, INT4*x2) { pnt2i(*x1, *x2);}
void YGL_PRE(pnt2s_) (INT2*x1, INT2*x2) { pnt2s(*x1, *x2);}

/* Lines */
void YGL_PRE(move2_) (REAL*x1, REAL*x2) { move2 (*x1, *x2);}
void YGL_PRE(move2i_)(INT4*x1, INT4*x2) { move2i(*x1, *x2);}
void YGL_PRE(move2s_)(INT2*x1, INT2*x2) { move2s(*x1, *x2);}

void YGL_PRE(rmv2_)  (REAL*x1, REAL*x2) { rmv2 (*x1, *x2);}
void YGL_PRE(rmv2i_) (INT4*x1, INT4*x2) { rmv2i(*x1, *x2);}
void YGL_PRE(rmv2s_) (INT2*x1, INT2*x2) { rmv2s(*x1, *x2);}

void YGL_PRE(draw2_) (REAL*x1, REAL*x2) { draw2 (*x1, *x2);}
void YGL_PRE(draw2i_)(INT4*x1, INT4*x2) { draw2i(*x1, *x2);}
void YGL_PRE(draw2s_)(INT2*x1, INT2*x2) { draw2s(*x1, *x2);}

void YGL_PRE(rdr2_)  (REAL*x1, REAL*x2) { rdr2 (*x1, *x2);}
void YGL_PRE(rdr2i_) (INT4*x1, INT4*x2) { rdr2i(*x1, *x2);}
void YGL_PRE(rdr2s_) (INT2*x1, INT2*x2) { rdr2s(*x1, *x2);}

/* Arcs & Circles */
void YGL_PRE(arc_)   (REAL*x1, REAL*x2, REAL*x3, INT4*x4, INT4*x5) { arc  (*x1, *x2, *x3, *x4, *x5);}
void YGL_PRE(arci_)  (INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT4*x5) { arci (*x1, *x2, *x3, *x4, *x5);}
void YGL_PRE(arcs_)  (INT2*x1, INT2*x2, INT2*x3, INT4*x4, INT4*x5) { arcs (*x1, *x2, *x3, *x4, *x5);}

void YGL_PRE(arcf_)  (REAL*x1, REAL*x2, REAL*x3, INT4*x4, INT4*x5) { arcf (*x1, *x2, *x3, *x4, *x5);}
void YGL_PRE(arcfi_) (INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT4*x5) { arcfi(*x1, *x2, *x3, *x4, *x5);}
void YGL_PRE(arcfs_) (INT2*x1, INT2*x2, INT2*x3, INT4*x4, INT4*x5) { arcfs(*x1, *x2, *x3, *x4, *x5);}

void YGL_PRE(circ_)  (REAL*x1, REAL*x2, REAL*x3) { circ  (*x1, *x2, *x3);}
void YGL_PRE(circi_) (INT4*x1, INT4*x2, INT4*x3) { circi (*x1, *x2, *x3);}
void YGL_PRE(circs_) (INT2*x1, INT2*x2, INT2*x3) { circs (*x1, *x2, *x3);}

void YGL_PRE(circf_) (REAL*x1, REAL*x2, REAL*x3) { circf (*x1, *x2, *x3);}
void YGL_PRE(circfi_)(INT4*x1, INT4*x2, INT4*x3) { circfi(*x1, *x2, *x3);}
void YGL_PRE(circfs_)(INT2*x1, INT2*x2, INT2*x3) { circfs(*x1, *x2, *x3);}

/* Rects & Boxes */
void YGL_PRE(rect_)  (REAL*x1, REAL*x2, REAL*x3, REAL*x4) { rect  (*x1, *x2, *x3, *x4);}
void YGL_PRE(recti_) (INT4*x1, INT4*x2, INT4*x3, INT4*x4) { recti (*x1, *x2, *x3, *x4);}
void YGL_PRE(rects_) (INT2*x1, INT2*x2, INT2*x3, INT2*x4) { rects (*x1, *x2, *x3, *x4);}

void YGL_PRE(rectf_) (REAL*x1, REAL*x2, REAL*x3, REAL*x4) { rectf (*x1, *x2, *x3, *x4);}
void YGL_PRE(rectfi_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4) { rectfi(*x1, *x2, *x3, *x4);}
void YGL_PRE(rectfs_)(INT2*x1, INT2*x2, INT2*x3, INT2*x4) { rectfs(*x1, *x2, *x3, *x4);}

void YGL_PRE(sbox_)  (REAL*x1, REAL*x2, REAL*x3, REAL*x4) { sbox  (*x1, *x2, *x3, *x4);}
void YGL_PRE(sboxi_) (INT4*x1, INT4*x2, INT4*x3, INT4*x4) { sboxi (*x1, *x2, *x3, *x4);}
void YGL_PRE(sboxs_) (INT2*x1, INT2*x2, INT2*x3, INT2*x4) { sboxs (*x1, *x2, *x3, *x4);}

void YGL_PRE(sboxf_) (REAL*x1, REAL*x2, REAL*x3, REAL*x4) { sboxf (*x1, *x2, *x3, *x4);}
void YGL_PRE(sboxfi_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4) { sboxfi(*x1, *x2, *x3, *x4);}
void YGL_PRE(sboxfs_)(INT2*x1, INT2*x2, INT2*x3, INT2*x4) { sboxfs(*x1, *x2, *x3, *x4);}

/* Filled Polygons */
void YGL_PRE(concav_)(LOGICAL*x1) { concave(*x1);}

void YGL_PRE(pmv2_)  (REAL*x1, REAL*x2) { pmv2 (*x1, *x2);}
void YGL_PRE(pmv2i_) (INT4*x1, INT4*x2) { pmv2i(*x1, *x2);}
void YGL_PRE(pmv2s_) (INT2*x1, INT2*x2) { pmv2s(*x1, *x2);}

void YGL_PRE(rpmv2_) (REAL*x1, REAL*x2) { rpmv2 (*x1, *x2);}
void YGL_PRE(rpmv2i_)(INT4*x1, INT4*x2) { rpmv2i(*x1, *x2);}
void YGL_PRE(rpmv2s_)(INT2*x1, INT2*x2) { rpmv2s(*x1, *x2);}

void YGL_PRE(pdr2_)  (REAL*x1, REAL*x2) { pdr2 (*x1, *x2);}
void YGL_PRE(pdr2i_) (INT4*x1, INT4*x2) { pdr2i(*x1, *x2);}
void YGL_PRE(pdr2s_) (INT2*x1, INT2*x2) { pdr2s(*x1, *x2);}

void YGL_PRE(rpdr2_) (REAL*x1, REAL*x2) { rpdr2 (*x1, *x2);}
void YGL_PRE(rpdr2i_)(INT4*x1, INT4*x2) { rpdr2i(*x1, *x2);}
void YGL_PRE(rpdr2s_)(INT2*x1, INT2*x2) { rpdr2s(*x1, *x2);}

void YGL_PRE(pclos_) (void) { pclos();}

void YGL_PRE(poly2_) (INT4*x1, REAL x2[][2]) { polf2 (*x1, x2); }
void YGL_PRE(poly2i_)(INT4*x1, INT4 x2[][2]) { polf2i(*x1, x2); }
void YGL_PRE(poly2s_)(INT4*x1, INT2 x2[][2]) { polf2s(*x1, x2); }

void YGL_PRE(polf2_) (INT4*x1, REAL x2[][2]) { polf2 (*x1, x2); }
void YGL_PRE(polf2i_)(INT4*x1, INT4 x2[][2]) { polf2i(*x1, x2); }
void YGL_PRE(polf2s_)(INT4*x1, INT2 x2[][2]) { polf2s(*x1, x2); }

/* Vertex graphics */
void YGL_PRE(bgnpoi_)(void) { bgnpoint     ();}
void YGL_PRE(bgnlin_)(void) { bgnline      ();}
void YGL_PRE(bgnclo_)(void) { bgnclosedline();}
void YGL_PRE(bgnpol_)(void) { bgnpolygon   ();}
void YGL_PRE(bgntme_)(void) { bgntmesh     ();}

void YGL_PRE(endpoi_)(void) { endpoint     ();}
void YGL_PRE(endlin_)(void) { endline      ();}
void YGL_PRE(endclo_)(void) { endclosedline();}
void YGL_PRE(endpol_)(void) { endpolygon   ();}
void YGL_PRE(endtme_)(void) { endtmesh     ();}

void YGL_PRE(v2s_)   (INT2  x1[2]) { v2s(x1);}
void YGL_PRE(v2i_)   (INT4  x1[2]) { v2i(x1);}
void YGL_PRE(v2f_)   (REAL  x1[2]) { v2f(x1);}
void YGL_PRE(v2d_)   (REAL8 x1[2]) { v2d(x1);}

/* Text */
void YGL_PRE(cmov2_) (REAL*x1, REAL*x2) { cmov2 (*x1, *x2);}
void YGL_PRE(cmov2i_)(INT4*x1, INT4*x2) { cmov2i(*x1, *x2);}
void YGL_PRE(cmov2s_)(INT2*x1, INT2*x2) { cmov2s(*x1, *x2);}

void YGL_PRE(getcpo_)(INT2*x1, INT2*x2) { getcpos(x1, x2);}

/* Extensions: Routines not in gl by MiSt (michael@hal6000.thp.Uni-Duisburg.DE) */
#ifdef X11
void YGL_PRE(arcx_)  (REAL*x1, REAL*x2, REAL*x3, REAL*x4, INT4*x5, INT4*x6) { arcx  (*x1,*x2,*x3,*x4,*x5,*x6); }
void YGL_PRE(arcxi_) (INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT4*x5, INT4*x6) { arcxi (*x1,*x2,*x3,*x4,*x5,*x6); }
void YGL_PRE(arcxs_) (INT2*x1, INT2*x2, INT2*x3, INT2*x4, INT4*x5, INT4*x6) { arcxs (*x1,*x2,*x3,*x4,*x5,*x6); }

void YGL_PRE(arcxf_) (REAL*x1, REAL*x2, REAL*x3, REAL*x4, INT4*x5, INT4*x6) { arcxf (*x1,*x2,*x3,*x4,*x5,*x6); }
void YGL_PRE(arcxfi_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT4*x5, INT4*x6) { arcxfi(*x1,*x2,*x3,*x4,*x5,*x6); }
void YGL_PRE(arcxfs_)(INT2*x1, INT2*x2, INT2*x3, INT2*x4, INT4*x5, INT4*x6) { arcxfs(*x1,*x2,*x3,*x4,*x5,*x6); }
#endif

/********************* queue.c */
void YGL_PRE(tie_)   (INT4*x1, INT4*x2, INT4*x3) { tie(*x1, *x2, *x3); }
void YGL_PRE(noise_) (INT4*x1, INT4*x2) { noise(*x1, *x2); }
INT4 YGL_PRE(isqueu_)(INT4*x1) { return isqueued(*x1); }
void YGL_PRE(qdevic_)(INT4*x1) { qdevice(*x1); }
void YGL_PRE(unqdev_)(INT4*x1) { unqdevice(*x1); }
void YGL_PRE(qreset_)(void)    { qreset(); }
INT4 YGL_PRE(qtest_) (void)    { return qtest(); }
INT4 YGL_PRE(qread_) (INT2*x1) { return qread(x1); }
void YGL_PRE(qenter_)(INT4*x1, INT4*x2) { qenter(*x1, *x2); }
#ifdef OGL
void YGL_PRE(pick_)  (INT2*x1, INT4*x2) { pick(x1, *x2); }
INT4 YGL_PRE(endpic_)(INT2*x1) { return (INT4)endpick(x1); }
void YGL_PRE(picksi_)(INT4*x1, INT4*x2) { picksize(*x1, *x2); }

void YGL_PRE(initna_)(void)    { initnames(); }
void YGL_PRE(loadna_)(INT4*x1) { loadname(*x1); }
void YGL_PRE(pushna_)(INT4*x1) { pushname(*x1); }
void YGL_PRE(popnam_)(void)    { popname(); }
#endif /* OGL */
/********************* misc.c */
void YGL_PRE(single_)(void)       { singlebuffer(); }
void YGL_PRE(doublebuffer_)(void) { doublebuffer(); } /* > 6 */
void YGL_PRE(double_)(void)       { doublebuffer(); }
void YGL_PRE(swapbu_)(void)       { swapbuffers(); }
void YGL_PRE(frontb_)(LOGICAL*x1) { frontbuffer(*x1); }
void YGL_PRE(backbu_)(LOGICAL*x1) { backbuffer (*x1); }

void YGL_PRE(gflush_)(void) { gflush(); }
void YGL_PRE(gsync_) (void) { gsync(); }


INT4 YGL_PRE(getxdpy_)(void) { return (INT4)getXdpy(); } /* > 6 */
INT4 YGL_PRE(getxwid_)(void) { return (INT4)getXwid(); } /* > 6 */
#ifdef X11
INT4 YGL_PRE(getxdid_)(void) { return (INT4)getXdid(); } /* > 6 */
INT4 YGL_PRE(getxgc_) (void) { return (INT4)getXgc (); }
#endif


void YGL_PRE(wintit_)(CHAR*x1, INT4*len) { wintitle(apnd_0(x1,len));}
void YGL_PRE(winset_)(INT4*x1) { winset(*x1);}
INT4 YGL_PRE(winget_)(void)    { return winget();}
INT4 YGL_PRE(getpla_)(void)    { return getplanes();}
INT4 YGL_PRE(getval_)(INT4*x1) { return getvaluator(*x1);}
INT4 YGL_PRE(getbut_)(INT4*x1) { return getbutton(*x1);}
INT4 YGL_PRE(gversi_)(CHAR*x1, INT4*len) { INT4 r = gversion(x1);if(len) *len = strlen(x1);return r;}

void YGL_PRE(ortho2_)(REAL*x1, REAL*x2, REAL*x3, REAL*x4) { ortho2(*x1, *x2, *x3, *x4);}
void YGL_PRE(viewpo_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4) { viewport(*x1, *x2, *x3, *x4);}
void YGL_PRE(getvie_)(INT2*x1, INT2*x2, INT2*x3, INT2*x4) { getviewport(x1, x2, x3, x4);}
void YGL_PRE(reshap_)(void)    { reshapeviewport();}
void YGL_PRE(pushvi_)(void)    { pushviewport();}
void YGL_PRE(popvie_)(void)    { popviewport();}

void YGL_PRE(winpop_)(void)    { winpop();}
void YGL_PRE(winpus_)(void)    { winpush();}
INT4 YGL_PRE(windep_)(INT4*x1) { return windepth(*x1);}

void YGL_PRE(linewi_)(INT4*x1) { linewidth(*x1);}
INT4 YGL_PRE(getlwi_)(void)    { return getlwidth();}
void YGL_PRE(deflin_)(INT4*x1, INT4*x2) { deflinestyle(*x1, *x2);}
void YGL_PRE(setlin_)(INT4*x1) { setlinestyle(*x1);}
INT4 YGL_PRE(getlst_)(void)    { return getlstyle();}
void YGL_PRE(lsrepe_)(INT4*x1) { lsrepeat(*x1);}
INT4 YGL_PRE(getlsr_)(void)    { return getlsrepeat();}
INT4 YGL_PRE(getdis_)(void)    { return getdisplaymode();}

void YGL_PRE(setbel_)(CHAR*x1) { setbell(*x1);}
void YGL_PRE(ringbe_)(void)    { ringbell();}

INT4 YGL_PRE(getgde_)(INT4*x1) { return getgdesc(*x1);}

void YGL_PRE(foregr_)(void)    { }

void YGL_PRE(logico_)(INT4*x1) { logicop(*x1);}

void YGL_PRE(getmat_)(Matrix x1) { getmatrix(x1);}

/********************* font.c */
void YGL_PRE(loadxf_)(INT4*x1, CHAR*x2, INT4*len) { loadXfont(*x1, apnd_0(x2,len)); }
void YGL_PRE(font_)  (INT4*x1) { font(*x1);}
INT4 YGL_PRE(getfon_)(void) { return getfont();}
void YGL_PRE(getfontencoding_)(CHAR*x1, INT4*len) { getfontencoding(x1);if(len) *len = strlen(x1);} /* > 6 */
INT4 YGL_PRE(gethei_)(void) { return getheight();}
INT4 YGL_PRE(getdes_)(void) { return getdescender();}
INT4 YGL_PRE(strwid_)(CHAR*x1, INT4*len) { return strwidth(apnd_0(x1,len));}
void YGL_PRE(charst_)(CHAR*x1, INT4*len) { charstr(apnd_0(x1,len));}

/********************* color.c */
void YGL_PRE(mapcol_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4) { mapcolor(*x1, *x2, *x3, *x4);}
void YGL_PRE(rgbcol_)(INT4*x1, INT4*x2, INT4*x3) { RGBcolor(*x1, *x2, *x3);}
void YGL_PRE(cpack_) (INT4*x1) { cpack(*x1);}

void YGL_PRE(c3s_)   (INT2 x1[3]) { c3s(x1);}
void YGL_PRE(c3i_)   (INT4 x1[3]) { c3i(x1);}
void YGL_PRE(c3f_)   (REAL x1[3]) { c3f(x1);}

INT4 YGL_PRE(getcol_)(void) { return getcolor();}
void YGL_PRE(getmco_)(INT4*x1, INT2*x2, INT2*x3, INT2*x4) { getmcolor (*x1, x2, x3, x4); }
void YGL_PRE(getmcolors_)(INT4*x1, INT4*x2, INT2*x3, INT2*x4, INT2*x5) { getmcolors(*x1, *x2, x3, x4, x5); }
void YGL_PRE(grgbco_)(INT2*x1, INT2*x2, INT2*x3) { gRGBcolor(x1, x2, x3);}

void YGL_PRE(color_) (INT4*x1) { color(*x1);}
void YGL_PRE(readso_)(INT4*x1) { readsource(*x1);}

void YGL_PRE(rectzo_)(REAL*x1, REAL*x2) { rectzoom(*x1, *x2);}

INT4 YGL_PRE(crectr_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4,Uint8*x5) { return crectread(*x1,*x2,*x3,*x4,x5); }
INT4 YGL_PRE(rectre_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT2*x5) { return  rectread(*x1,*x2,*x3,*x4,x5); }
INT4 YGL_PRE(lrectr_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT4*x5) { return lrectread(*x1,*x2,*x3,*x4,x5); }

void YGL_PRE(crectw_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4,Uint8*x5) { crectwrite(*x1,*x2,*x3,*x4,x5); }
void YGL_PRE(rectwr_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT2*x5) {  rectwrite(*x1,*x2,*x3,*x4,x5); }
void YGL_PRE(lrectw_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT4*x5) { lrectwrite(*x1,*x2,*x3,*x4,x5); }

void YGL_PRE(rectco_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT4*x5, INT4*x6) { rectcopy(*x1,*x2,*x3,*x4,*x5,*x6); }

INT4 YGL_PRE(readpi_)(INT4*x1, INT2 x2[]) { return readpixels(*x1, (Colorindex*)x2); }
void YGL_PRE(writep_)(INT4*x1, INT2 x2[]) {       writepixels(*x1, (Colorindex*)x2); }
INT4 YGL_PRE(readrg_)(INT4*x1, CHAR*x2, CHAR*x3, CHAR*x4) { return readRGB(*x1,(RGBvalue*)x2,(RGBvalue*)x3,(RGBvalue*)x4); }
void YGL_PRE(writer_)(INT4*x1, CHAR*x2, CHAR*x3, CHAR*x4) {       writeRGB(*x1,(RGBvalue*)x2,(RGBvalue*)x3,(RGBvalue*)x4); }


void YGL_PRE(blendf_)(INT4*x1, INT4*x2) { blendfunction(*x1,*x2); }

/********************* menu.c */
void YGL_PRE(addtop_)(INT4*x1, CHAR*x2, INT4*len, INT4*x3) { addtopup(*x1, apnd_0(x2,len), *x3);}
/* defpup() not available in FORTRAN */
INT4 YGL_PRE(dopup_) (INT4*x1) { return dopup(*x1);}
void YGL_PRE(freepu_)(INT4*x1) { freepup(*x1);}
INT4 YGL_PRE(newpup_)(void)    { return newpup();}
void YGL_PRE(setpup_)(INT4*x1, INT4*x2, INT4*x3) { setpup(*x1, *x2, *x3);}

/********************* ygl.c */
/* Contraints */
void YGL_PRE(minsiz_)(INT4*x1, INT4*x2) { minsize(*x1, *x2);}
void YGL_PRE(maxsiz_)(INT4*x1, INT4*x2) { maxsize(*x1, *x2);}
void YGL_PRE(prefsi_)(INT4*x1, INT4*x2) { prefsize(*x1, *x2);}
void YGL_PRE(prefpo_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4) { prefposition(*x1, *x2, *x3, *x4);}
void YGL_PRE(stepun_)(INT4*x1, INT4*x2) { stepunit(*x1, *x2);}
void YGL_PRE(keepas_)(INT4*x1, INT4*x2) { keepaspect(*x1, *x2);}
void YGL_PRE(noport_)(void) { noport();}
void YGL_PRE(nobord_)(void) { noborder();}

void YGL_PRE(ginit_) (void) { ginit();}
void YGL_PRE(wincon_)(void) { winconstraints();}
INT4 YGL_PRE(winope_)(CHAR*x1, INT4*len) { return winopen(apnd_0(x1,len));}
INT4 YGL_PRE(swinop_)(INT4*x1) { return swinopen(*x1);}

void YGL_PRE(winpos_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4) { winposition(*x1, *x2, *x3, *x4);}
void YGL_PRE(winmov_)(INT4*x1, INT4*x2) { winmove(*x1, *x2);}

void YGL_PRE(getsiz_)(INT4*x1, INT4*x2) { getsize  (x1,  x2);}
void YGL_PRE(getori_)(INT4*x1, INT4*x2) { getorigin(x1,  x2);}

void YGL_PRE(rgbmod_)(void) { RGBmode();}
void YGL_PRE(cmode_) (void) { cmode();}

void YGL_PRE(gconfi_)(void) { gconfig();}
void YGL_PRE(winclo_)(INT4*x1) { winclose(*x1);}
void YGL_PRE(gexit_) (void) { gexit();}


INT4 YGL_PRE(winx_)  (INT4*dpy, INT4*win) { return winX((Display*)dpy, (Window)*win);}


/* gl2ppm.c */
INT4 YGL_PRE(gl2ppm_)(CHAR*x1, INT4*len) { return gl2ppm(apnd_0(x1,len)); }

/* 3d.c */
#ifdef OGL
void YGL_PRE(cmov_)  (REAL*x1, REAL*x2, REAL*x3) { cmov (*x1, *x2, *x3);}
void YGL_PRE(cmovi_) (INT4*x1, INT4*x2, INT4*x3) { cmovi(*x1, *x2, *x3);}
void YGL_PRE(cmovs_) (INT2*x1, INT2*x2, INT2*x3) { cmovs(*x1, *x2, *x3);}

void YGL_PRE(pnt_)   (REAL*x1, REAL*x2, REAL*x3) { pnt (*x1, *x2, *x3);}
void YGL_PRE(pnti_)  (INT4*x1, INT4*x2, INT4*x3) { pnti(*x1, *x2, *x3);}
void YGL_PRE(pnts_)  (INT2*x1, INT2*x2, INT2*x3) { pnts(*x1, *x2, *x3);}

void YGL_PRE(move_)  (REAL*x1, REAL*x2, REAL*x3) { move (*x1, *x2, *x3);}
void YGL_PRE(movei_) (INT4*x1, INT4*x2, INT4*x3) { movei(*x1, *x2, *x3);}
void YGL_PRE(moves_) (INT2*x1, INT2*x2, INT2*x3) { moves(*x1, *x2, *x3);}

void YGL_PRE(rmv_)   (REAL*x1, REAL*x2, REAL*x3) { rmv (*x1, *x2, *x3);}
void YGL_PRE(rmvi_)  (INT4*x1, INT4*x2, INT4*x3) { rmvi(*x1, *x2, *x3);}
void YGL_PRE(rmvs_)  (INT2*x1, INT2*x2, INT2*x3) { rmvs(*x1, *x2, *x3);}

void YGL_PRE(draw_)  (REAL*x1, REAL*x2, REAL*x3) { draw (*x1, *x2, *x3);}
void YGL_PRE(drawi_) (INT4*x1, INT4*x2, INT4*x3) { drawi(*x1, *x2, *x3);}
void YGL_PRE(draws_) (INT2*x1, INT2*x2, INT2*x3) { draws(*x1, *x2, *x3);}

void YGL_PRE(rdr_)   (REAL*x1, REAL*x2, REAL*x3) { rdr (*x1, *x2, *x3);}
void YGL_PRE(rdri_)  (INT4*x1, INT4*x2, INT4*x3) { rdri(*x1, *x2, *x3);}
void YGL_PRE(rdrs_)  (INT2*x1, INT2*x2, INT2*x3) { rdrs(*x1, *x2, *x3);}

void YGL_PRE(pmv_)   (REAL*x1, REAL*x2, REAL*x3) { pmv (*x1, *x2, *x3);}
void YGL_PRE(pmvi_)  (INT4*x1, INT4*x2, INT4*x3) { pmvi(*x1, *x2, *x3);}
void YGL_PRE(pmvs_)  (INT2*x1, INT2*x2, INT2*x3) { pmvs(*x1, *x2, *x3);}

void YGL_PRE(rpmv_)  (REAL*x1, REAL*x2, REAL*x3) { rpmv (*x1, *x2, *x3);}
void YGL_PRE(rpmvi_) (INT4*x1, INT4*x2, INT4*x3) { rpmvi(*x1, *x2, *x3);}
void YGL_PRE(rpmvs_) (INT2*x1, INT2*x2, INT2*x3) { rpmvs(*x1, *x2, *x3);}

void YGL_PRE(pdr_)   (REAL*x1, REAL*x2, REAL*x3) { pdr (*x1, *x2, *x3);}
void YGL_PRE(pdri_)  (INT4*x1, INT4*x2, INT4*x3) { pdri(*x1, *x2, *x3);}
void YGL_PRE(pdrs_)  (INT2*x1, INT2*x2, INT2*x3) { pdrs(*x1, *x2, *x3);}

void YGL_PRE(rpdr_)  (REAL*x1, REAL*x2, REAL*x3) { rpdr (*x1, *x2, *x3);}
void YGL_PRE(rpdri_) (INT4*x1, INT4*x2, INT4*x3) { rpdri(*x1, *x2, *x3);}
void YGL_PRE(rpdrs_) (INT2*x1, INT2*x2, INT2*x3) { rpdrs(*x1, *x2, *x3);}

void YGL_PRE(poly_)  (INT4*x1, REAL x2[][3]) { polf (*x1, x2); }
void YGL_PRE(polyi_) (INT4*x1, INT4 x2[][3]) { polfi(*x1, x2); }
void YGL_PRE(polys_) (INT4*x1, INT2 x2[][3]) { polfs(*x1, x2); }

void YGL_PRE(polf_)  (INT4*x1, REAL x2[][3]) { polf (*x1, x2); }
void YGL_PRE(polfi_) (INT4*x1, INT4 x2[][3]) { polfi(*x1, x2); }
void YGL_PRE(polfs_) (INT4*x1, INT2 x2[][3]) { polfs(*x1, x2); }

void YGL_PRE(v3s_)   (INT2  x1[3]) { v3s(x1);}
void YGL_PRE(v3i_)   (INT4  x1[3]) { v3i(x1);}
void YGL_PRE(v3f_)   (REAL  x1[3]) { v3f(x1);}
void YGL_PRE(v3d_)   (REAL8 x1[3]) { v3d(x1);}

void YGL_PRE(v4s_)   (INT2  x1[4]) { v4s(x1);}
void YGL_PRE(v4i_)   (INT4  x1[4]) { v4i(x1);}
void YGL_PRE(v4f_)   (REAL  x1[4]) { v4f(x1);}
void YGL_PRE(v4d_)   (REAL8 x1[4]) { v4d(x1);}

void YGL_PRE(swaptm_)(void) { swaptmesh();}

void YGL_PRE(ortho_) (REAL*x1, REAL*x2, REAL*x3, REAL*x4, REAL*x5, REAL*x6) { ortho(*x1, *x2, *x3, *x4, *x5, *x6);}
void YGL_PRE(lookat_)(REAL*x1, REAL*x2, REAL*x3, REAL*x4, REAL*x5, REAL*x6, INT4*x7) { lookat(*x1, *x2, *x3, *x4, *x5, *x6, *x7);}
void YGL_PRE(window_)(REAL*x1, REAL*x2, REAL*x3, REAL*x4, REAL*x5, REAL*x6) { window(*x1, *x2, *x3, *x4, *x5, *x6);}
void YGL_PRE(perspe_)(INT4*x1, REAL*x2, REAL*x3, REAL*x4) { perspective(*x1, *x2, *x3, *x4);}
void YGL_PRE(polarv_)(REAL*x1, INT4*x2, INT4*x3, INT4*x4) { polarview(*x1, *x2, *x3, *x4);}
void YGL_PRE(rot_)   (REAL*x1, CHAR*x2) { rot(*x1, *x2);}
void YGL_PRE(rotate_)(INT4*x1, CHAR*x2) { rotate(*x1, *x2);}
void YGL_PRE(scale_) (REAL*x1, REAL*x2, REAL*x3) { scale(*x1, *x2, *x3);}
void YGL_PRE(transl_)(REAL*x1, REAL*x2, REAL*x3) { translate(*x1, *x2, *x3);}

void YGL_PRE(loadma_)(Matrix x1) { loadmatrix(x1);}
void YGL_PRE(multma_)(Matrix x1) { multmatrix(x1);}
void YGL_PRE(pushma_)(void) { pushmatrix();}
void YGL_PRE(popmat_)(void) { popmatrix();}

void YGL_PRE(shadem_)(INT4*x1) { shademodel(*x1);}
  
void YGL_PRE(c4s_)   (INT2 x1[4]) { c4s(x1);}
void YGL_PRE(c4i_)   (INT4 x1[4]) { c4i(x1);}
void YGL_PRE(c4f_)   (REAL x1[4]) { c4f(x1);}

void YGL_PRE(n3f_)   (REAL x1[3]) { n3f(x1);}
void YGL_PRE(normal_)(REAL x1[3]) { normal(x1);}

void YGL_PRE(backfa_)(LOGICAL*x1) { backface(*x1);}
void YGL_PRE(frontf_)(LOGICAL*x1) { frontface(*x1);}

INT4 YGL_PRE(getmmo_)(void) { return getmmode();}
void YGL_PRE(mmode_) (INT4*x1) { mmode(*x1);}

void YGL_PRE(zbuffe_)(LOGICAL*x1) { zbuffer(*x1);}
void YGL_PRE(zclear_)(void)    { zclear();}
void YGL_PRE(zdraw_) (INT4*x1) { zdraw(*x1);}
void YGL_PRE(zfunct_)(INT4*x1) { zfunction(*x1);}
void YGL_PRE(czclea_)(INT4*x1, INT4*x2) { czclear(*x1, *x2);}
void YGL_PRE(depthc_)(INT4*x1) { depthcue(*x1);}
void YGL_PRE(lrgbra_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4, INT4*x5, INT4*x6, INT4*x7, INT4*x8) { lRGBrange(*x1, *x2, *x3, *x4, *x5, *x6, *x7, *x8);}
void YGL_PRE(lsetde_)(INT4*x1, INT4*x2) { lsetdepth(*x1, *x2);}
void YGL_PRE(lshade_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4) { lshaderange(*x1, *x2, *x3, *x4);}

/* Display lists */
INT4 YGL_PRE(genobj_)(void) { return genobj();}
INT4 YGL_PRE(isobj_) (INT4*x1) { return isobj(*x1);}
void YGL_PRE(makeob_)(INT4*x1) { makeobj(*x1);}
INT4 YGL_PRE(getope_)(void) { return getopenobj();}
void YGL_PRE(closeo_)(void) { closeobj();}
void YGL_PRE(callob_)(INT4*x1) { callobj(*x1);}
void YGL_PRE(delobj_)(INT4*x1) { delobj(*x1);}

/* Lighting */
void YGL_PRE(lmbind_)(INT4*x1, INT4*x2) { lmbind(*x1, *x2);}
void YGL_PRE(lmcolo_)(INT4*x1) { lmcolor(*x1);}
void YGL_PRE(lmdef_) (INT4*x1, INT4*x2, INT4*x3, REAL x4[]) { lmdef(*x1, *x2, *x3, x4);}

void YGL_PRE(rgbwri_)(INT4*x1, INT4*x2, INT4*x3) { RGBwritemask(*x1, *x2, *x3);}

/* for Pete Riley */
void YGL_PRE(drawmo_)(INT4*x1) { drawmode(*x1);}
void YGL_PRE(iconsi_)(INT4*x1, INT4*x2) { iconsize(*x1, *x2);}
void YGL_PRE(overla_)(INT4*x1) { overlay(*x1);}
void YGL_PRE(pushat_)(void)    { pushattributes();}
void YGL_PRE(popatt_)(void)    { popattributes();}
void YGL_PRE(fullsc_)(void)    { fullscrn();}
void YGL_PRE(endful_)(void)    { endfullscrn();}
void YGL_PRE(scrmas_)(INT4*x1, INT4*x2, INT4*x3, INT4*x4) { scrmask(*x1,*x2,*x3,*x4);}

#endif /* OGL */

#ifdef COVERSLEEP
void YGL_PRE(usleep_)(INT4*x1) {usleep(*x1);}
void YGL_PRE(sleep_) (INT4*x1) { sleep(*x1);}
#endif /* COVERSLEEP */
