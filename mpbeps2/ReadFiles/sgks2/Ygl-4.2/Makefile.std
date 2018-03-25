#
#    Ygl: Run 2d-GL programs with standard X11 routines.
#    (C) Fred Hucht 1993-2002
#    EMail: fred@thp.Uni-Duisburg.DE
#
#    $Id: Makefile.std,v 3.4 1996-07-18 18:38:53+02 fred Exp fred $

# Uncomment next line to use X11 bindings (faster, but only 2d)
X11		= -DX11

# Uncomment next line to use OpenGL bindings (slower than X11, but also 3d)
OGL		= -DOGL

# Uncomment next two lines to include FORTRAN bindings
FBO		= fortran.o
FBH		= X11/Yfgl.h

# Uncomment next line to use DoubleBuffer extension with X11
DOUBLEBUF 	= -DDOUBLEBUF

# Uncomment next line to use MultiBuffer extension with X11
MULTIBUF 	= -DMULTIBUF

# Uncomment next two lines to prepend "ygl_" to all function names 
#YGL_PREFIX	= -DYGL_PREFIX
#PH		= X11/Yglprefix.h
#TARGET		= libYglp.a

CDEBUGFLAGS	= -O

OBJS		= ygl.o draw.o misc.o font.o queue.o color.o menu.o gl2ppm.o $(FBO)
TARGET          = libYgl.a

COPTS = $(X11) $(OGL) $(DOUBLEBUF) $(MULTIBUF)
# COPTS = -Aa -D_HPUX_SOURCE -DMULTIBUF -DNO_MULTIBUF_H -L/usr/lib/X11R4 -I/usr/include/X11R4 # For HP-UX 8.0x
# COPTS = -Ae $(DOUBLEBUF) $(MULTIBUF) -L/usr/lib/X11R5 -I/usr/include/X11R5 # For HP-UX 9.0x

# End of configuration

CFLAGS		= -I. $(CDEBUGFLAGS) $(COPTS) $(YGL_PREFIX)

all: 	$(TARGET)

.c.o:	
	$(CC) -c $(CFLAGS) $<

X11/Yglprefix.h:	makeYglprefix X11/Ygl.h
	./makeYglprefix > $@

X11/Ygltypes.h:		makeYgltypes
	./makeYgltypes > $@

makeYgltypes:		makeYgltypes.c
	$(CC) $(CFLAGS) -o makeYgltypes makeYgltypes.c

usleep.i:	header.h
	$(RM) usleep_tst.c
	/bin/ln -s header.h usleep_tst.c
	$(CC) $(CFLAGS) -E usleep_tst.c > usleep.i
	$(RM) usleep_tst.c

usleep.h:	usleep.i
	sed -nf makeusleep.sed usleep.i > usleep.h

misc.o:		usleep.h

$(OBJS):	header.h config.h X11/Ygl.h X11/Ygltypes.h $(PH)

$(TARGET):	$(OBJS)
	/bin/rm -f $@
	ar rv $@ $(OBJS)
	ranlib $@

install:	$(TARGET)
	/bin/cp $(TARGET) /usr/lib/$(TARGET) 
	ranlib /usr/lib/$(TARGET)
	/bin/cp X11/Ygl.h X11/Ygltypes.h $(PH) $(FBH) /usr/include/X11

clean:
	/bin/rm -f *.o *~ $(TARGET) X11/Yglprefix.h X11/Ygltypes.h makeYgltypes usleep.h usleep.i

etags:
	etags *.[ch]
