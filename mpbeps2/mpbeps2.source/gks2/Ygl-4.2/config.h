/*
 *    Ygl: Run GL programs with standard X11 routines.
 *    (C) Fred Hucht 1993-96
 *    EMail: fred@thp.Uni-Duisburg.DE
 *
 *    $Id: config.h,v 3.5 1997-07-07 11:09:39+02 fred Exp $
 */

#define AUTOFLUSH	/* define to generate code for flushing the
			   server via timer */ 

#define COVERSLEEP	/* define to cover the system commands sleep
			   and usleep with own versions that flushes
			   Xlibs output buffer */ 

/*#define RGBWIN*/

#ifdef RGBWIN
# define IF_RGBWIN(x, y) x
#else
# define IF_RGBWIN(x, y) y
#endif
