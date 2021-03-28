#
#    Ygl: Run 3d-GL programs with standard X11 routines.
#    (C) Fred Hucht 1993-97
#    EMail: fred@thp.Uni-Duisburg.DE
#
#    $Id: Makefile.std,v 3.3 1997-03-13 16:19:21+01 fred Exp fred $

# Uncomment next line to use X11 bindings

#X11LIB	= -lX11
X11LIB = -L/usr/X11/lib -LX11

# Uncomment next line to use OpenGL bindings
#OGLIB	= -lGL -lGLU

# Uncomment next line to use DoubleBuffer or MultiBuffer extension with X11
MBUFLIB	= -lXext

CDEBUGFLAGS = -O

COPTS =
# COPTS = -Aa -D_HPUX_SOURCE -L/usr/lib/X11R4 -I/usr/include/X11R4 # For HP-UX 8.0x
# COPTS = -Ae -L/usr/lib/X11R5 -I/usr/include/X11R5 # For HP-UX 9.0x

# *************************** End of configuration *************************

.SUFFIXES: .gl

CFLAGS= $(CDEBUGFLAGS) $(COPTS) -I..

REQUIREDLIBS 	= -L.. -lYgl $(MBUFLIB) $(X11LIB) $(OGLIB)

GLTARGETS	= coltest.gl lines.gl rgbtest.gl smile.gl popup.gl lmbind.gl
YGLTARGETS	= coltest    lines    rgbtest    smile    popup    lmbind

all: 	$(YGLTARGETS)
gl: 	$(GLTARGETS)

.c:
	$(CC) $(CFLAGS) -o $@ $< $(REQUIREDLIBS) -lm

.c.gl:
	$(CC) $(CFLAGS) -o $@ $< ../gl2ppm.o -lgl -lm

clean:
	/bin/rm -f *.o *~ $(YGLTARGETS) $(GLTARGETS) *.gif smile_xlf smile_f77 smile_f2c smile_f77.c

smile_f2c: smile_f77.f
	f2c smile_f77.f
	$(CC) $(CFLAGS) -o $@ smile_f77.c -lf2c $(REQUIREDLIBS)

smile_f77: smile_f77.f
	$(FC) $(FFLAGS) -o $@ smile_f77.f $(REQUIREDLIBS)

smile_xlf: smile_xlf.f
	xlf $(FFLAGS) -o $@ smile_xlf.f $(REQUIREDLIBS)
