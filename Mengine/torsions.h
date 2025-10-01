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

#ifndef TORSIONS_H
#define TORSIONS_H

struct t_torsions {
        int ntor, **i14;
        float *v1,*v2,*v3,*v4,*v5,*v6;
        int   *ph1,*ph2,*ph3,*ph4,*ph5,*ph6;
        } torsions;

#endif

int isbond(int, int);
int isangle(int,int);
int istorsion(int, int);
void get_torsions(void);
int is_metal(char *);
int is_linear(int);
double dihdrl(int,int,int,int);
void etorsion(void);
void etorsion1(void);
void etorsion2(int);
int isbond(int,int);
int is_delocalbond(int,int);
void numeral(int,char *,int);
void angle(int,int,int,float *);
void four(char *,char *,char *,char *, char *);
int is_ring42(int,int);
int is_ring54(int, int, int, int);
void message_alert(char *, char *);
void  transomeg(float *, float *, float *, int *, int *, char *, char *, char *, char *,
         char *, char *, char *, char *, char *, char *, char *);
