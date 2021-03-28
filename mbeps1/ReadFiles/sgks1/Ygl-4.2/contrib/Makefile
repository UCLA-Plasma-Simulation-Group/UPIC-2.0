
.SUFFIXES: .gl

TARGETS	= ycpu yup xmap yglmap

all: $(TARGETS)

.c.gl:
	$(CC) $(CFLAGS) -o $@ $< -lgl

.c:
	$(CC) $(CFLAGS) -o $@ $< -lYgl -lXext -lX11

glmap:	xmap.c
	$(CC)          $(CFLAGS) -o $@ $? -lgl

yglmap:	xmap.c
	$(CC)          $(CFLAGS) -o $@ $? -lYgl -lXext -lX11

xmap:	xmap.c
	$(CC) -DXCOLOR $(CFLAGS) -o $@ $? -lYgl -lXext -lX11

clean:
	/bin/rm -f *.o *~ $(TARGETS)
