/* popup.c by Fred Hucht (C) 1993-2002.
 * Example program for popup menus */

static const char vcid[] = "$Id: popup.c,v 3.2 1996/07/18 16:35:57 fred Exp $";

#include <X11/Ygl.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

/* Callback functions */
Int32 fn1(Int32 x) { printf("---> fn1(%d)\n", x); return(100 * x);}
Int32 fn2(Int32 x) { printf("---> fn2(%d)\n", x); return(10  + x);}
Int32 fe0(Int32 m) { return(m-1); }
Int32 fe1(Int32 m, Int32 s) { return(10 *(m-1) + s); }
Int32 fe2(Int32 m, Int32 s) { return(100*(m-1) + s); }
Int32 End(Int32 x) { return(4711); }

int main() {
  Int32 menu, sub, subsub, dev, r, sub1, sub10, sub100;
  Int16 val;
  minsize(300, 300);
  winopen("Menu");
  ortho2(-1.0, 1.0, -1.0, 1.0);
  
  subsub  = defpup("Sub%%sub!%t%F|#s4%x34|#s5%x35",
		   fn1);	/* %F */
  
  sub     = defpup("Sub!%t%F|#s4%x24%m|#s5%x25", 
		   fn1,		/* %F */
		   subsub);	/* %m */
  
  sub1    = defpup("0%t%F|0|1|2|3|4|5|6|7|8|9",
		   fe0);	/* %F: fe0 is applied to default return code */
  
  sub10   = defpup("00%t%F|00%M|10%M|20%M|30%M|40%M|50%M|60%M|70%M|80%M|90%M",
		   fe1,		/* %F */
		   /* All submenus are defined with %M,
		    * hence fe1(retval, subretval) is returned.
		    * Example: 30 (item #4) and 6 (item #7) selected 
		    *          => fe1(4, fe0(7)) == 36 */
		   sub1,sub1,sub1,sub1,sub1,sub1,sub1,sub1,sub1,sub1);
  
  sub100  = defpup("000%t%F|000%M|100%M|200%M|300%M|400%M|500%M|600%M|700%M|800%M|900%M",
		   fe2,		/* %F */
		   /* All submenus are defined with %M,
		    * hence fe2(retval, subretval) is returned.
		    * Example: 100 (item #2), 30 (item #4) and 6 (item #7) selected
		    *          => fe2(2, fe1(4, fe0(7))) == 136 */
		   sub10,sub10,sub10,sub10,sub10,sub10,sub10,sub10,sub10,sub10);
  
  menu    = defpup("The Title%t%F|Numbers%m|#2%l%f|#3%n|SubSub%m|#5|Sub%m|Exit%n",
		   fn1,		/* %F: overall menu callback */
		   /*		#1 */
		   sub100,	/* %m: item #1 has submenu 'sub100' */
		   /*		#2 */
		   fn2,		/* %f: return fn1(fn2(2)) == 1200 */
		   /*		#3 */
		   fn2,		/* %n: don't call fn1(), return fn2(3) == 13*/
		   /*		#4 */
		   subsub,	/* %m: item #4 has submenu 'subsub'  */
		   /*		#5 */
		   /*		#6 */
		   sub,		/* %m: item #6 has submenu 'sub' */
		   /*		#7 */
		   End);	/* %n: don't call fn1(), return End(7) == 4711 */
  
  setpup(menu, 5, PUP_GREY);	/* Greyout item #5 in menu menu */
  setpup(menu, 4, PUP_GREY);	/* Greyout item #4 in menu menu */
  
  qdevice(MENUBUTTON);
  
  loadXfont(2, "-*-times-medium-r-*-*-*-140-*-*-*-*-iso8859-1");
  font(2);
  
  while((dev = qread(&val))) {
    switch(dev) {
    case REDRAW:
      color(BLACK);
      clear();
      color(WHITE);
      circf(0.0, 0.0, 0.8);
      sleep(0);
      break;
    case MENUBUTTON: /* Right button */
      if(val == 1) {
	r = dopup(menu);
	printf("----> dopup returned %d\n", r);
	if(r == 4711) {
	  gexit();
	  exit(0);
	}
      }
      break;
    }
  }
  return 0;
}
