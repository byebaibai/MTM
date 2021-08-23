#include <stdio.h>
/* kern and Ritch page 26*/
int getline(FILE *input,char s[],int lim)
  {
  int c,i;

  for(i=0; i<lim-1 && (c=fgetc(input))!=EOF && c!= '\n'; ++i)
     s[i] = c;
     if( c=='\n') {
       s[i] = c;
       ++i;
       }
       s[i] = '\0';
       return(i);
  }
