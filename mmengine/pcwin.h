# Copyright (C) 2004-2016 Kevin E. Gilbert, dba Serena Software
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Kevin E. Gilbert
# Serena Software
# Box 3076
# Bloomington, IN 47402-3076
#
# gilbert@serenasoft.com

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
