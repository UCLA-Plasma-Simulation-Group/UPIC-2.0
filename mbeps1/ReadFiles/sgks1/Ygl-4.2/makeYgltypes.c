/*
 *    Ygl: Run GL programs with standard X11 routines.
 *    (C) Fred Hucht 1993-97
 *    EMail: fred@thp.Uni-Duisburg.DE
 *
 *    $Id: makeYgltypes.c,v 4.4 1997-09-15 21:03:06+02 fred Exp $
 */

#include <stdio.h>

int main() {
  int i;
  int done[] = {0,0};
  char *prefix[] = {" Int", "Float"};
  struct {
    int len;
    int i_f;
    char *name;
  } sz[] = {
    {sizeof(char),	0, "char  "},
    {sizeof(short),	0, "short "},
    {sizeof(int),	0, "int   "},
    {sizeof(long),	0, "long  "},
    {sizeof(float),	1, "float "},
    {sizeof(double),	1, "double"},
    {-1,		-1,""}
  };
  
  printf("/*\n"
	 " * Automagically created by makeYgltypes. Do not change!\n"
	 " */\n");
  
  for(i = 0; sz[i].len > 0; i++) {
    int i_f = sz[i].i_f;
    if(!(done[i_f] & sz[i].len)) {
      done[i_f] |= sz[i].len;
      printf("typedef %s%s %s%d;\n",
	     i_f ? "         " : "signed   ",
	     sz[i].name, prefix[i_f], 8*sz[i].len);
      if(i_f == 0) 
	printf("typedef unsigned %s Uint%d;\n",
	       sz[i].name, 8*sz[i].len);
    }
  }
  return 0;
}
