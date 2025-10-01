// Copyright (C) 2004-2016 Kevin E. Gilbert, dba Serena Software
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Kevin E. Gilbert
// Serena Software
// Box 3076
// Bloomington, IN 47402-3076
//
// gilbert@serenasoft.com

void densty(int, float **x, float **y, float);
void printarray(float a[], char *s);
void printmatrix(float **, char *s);
double sign(double a,double b );
void linear(float **, float *, int);
void square(float *, float **, int);
void scfchk(int, int, int *,float *);
void diag(int,float **);
void frag_diag(float **, float **, float *);
void xdiag(int,float **, float **, float *);
void fofgen(float **,float **,float **, float **, int *, float **);
void otable(void);
void message_alert(char *,char *);
void vescf(void);
void printcont(float **,int,int,char *);
void randnum(long int *,float *);
