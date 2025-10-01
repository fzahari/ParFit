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

int find_bond(int,int);
int find_fixed_bond(int,int);
void bnd_corr(int,int, int,int,double *);
double dihdrl(int,int,int,int);
double deltaks(double,double);
void ebond(void);
void ebond1(void);
void ebond2(int);

FILE * fopen_path ( char * , char * , char * ) ;
void find_pibond(int *,char *, float *,float *,float *,float *,float *);
void find_metalbond(int *, int, int,char *,float *, float *);
void angle(int ,int ,int ,float *);
int is_ring31(int);
int is_ring41(int);
int is_ring51(int);
int is_ring61(int);
int is_ring42(int, int);
int is_ring52(int, int);
int is_ring62(int,int);
void pcmfout(int);
double arcos(float );
double angs(float,float,float,float,float ,float ,float ,float);
double powi(double, int);
void numeral(int, char *, int);
int is_delocalbond(int ia, int ib);
void message_alert(char *, char *);
double sign(double ,double  ); 
void transbond(int *,char *,float *,float *);
int is_bond(int,int);
double deltaks(double,double);
float get_bond(int, int);
void set_missing_constants(int);        
int get_missing_constants(void);
int kbond(void);
void kstrtor(void);
void kstrbnd(void);
void kvdw(void);
