#ifndef PCWIN_H
#define PCWIN_H

#define TRUE True
#define FALSE False

#ifndef PI   /* Avoid Linux Warnings! */
#define PI   3.14159265358979323846
#endif

#define Rad2Deg      (180.0/PI)
#define Deg2Rad      (PI/180.0)
#define AbsFun(a)    (((a)<0)? -(a) : (a))
#define MinFun(a,b)  (((a)<(b))? (a) : (b) )
#define MaxFun(a,b)  (((a)>(b))? (a) : (b) )


#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

struct FileInfoStruct {
  int  ftype ;
  char   path[255] ;
  char   fname[255] ;
} ;

typedef struct FileInfoStruct Boxstruct;
Boxstruct Openbox ,Savebox;

 char * fexists_error_message ( void ) ;
void * malloc_filename ( char * , char * ) ;
int fexists_path ( char * , char * ) ;
FILE * fopen_path ( char * , char * , char * ) ;
#endif
