#
#    Ygl: Run 2d-GL programs with standard X11 routines.
#    (C) Fred Hucht 1993-2002
#    EMail: fred@thp.Uni-Duisburg.DE
#
#    $Id: makeusleep.sed,v 1.2 2002-04-02 15:10:10+02 fred Exp fred $
#
/^[^#].*[ 	_]usleep *(/{
# Remove YGL_PREFIX if used
s/ygl_usleep/usleep/
i\
#if 0
p
i\
#endif
h
s/^\(.*\)[ 	]usleep.*$/#define USLEEP_RET_TYPE \1/
s/extern//
p
/void/{
i\
#define USLEEP_RET_VAL
}
/void/!{
i\
#define USLEEP_RET_VAL 0
}
g
s/^.*usleep *( *\(unsigned\) \([^ )]*\).*$/#define USLEEP_ARG_TYPE \1 \2/p
s/^.*usleep *( *\(signed\) \([^ )]*\).*$/#define USLEEP_ARG_TYPE \1 \2/p
s/^.*usleep *( *\([^ )]*\).*$/#define USLEEP_ARG_TYPE \1/p
}
