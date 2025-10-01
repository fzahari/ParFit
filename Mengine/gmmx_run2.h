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

FILE * fopen_path ( char * , char * , char * ) ;
void rdpcm(int,int);
void pcmfout(int);
int get_filename(void);
void hster(int, int *);
int lepimer(int);
int check_id(void);  
int check_stop(void);
void minimize(void);
void mmp22mod(int,int);
void pcmfin(int,int);
void initialize(void);
void generate_bonds(void);
void check_numfile(int);
int bad15(int);
void set_active(void);
void pimarkselection(void);
void hbondreset(void);
void pireset(void);
void mmxsub(int);
double distance(int,int);
void angle(int,int,int,float *);
double couple(int,int,int,int);
double dihdrl(int,int,int,int);
void dipole_dipole(void);
void charge_dipole(void);
void eheat(void);
void sort_file(char *, char *);
int  bstereo(int);
void message_alert(char *, char *);
void rdfile(int);
void bkboltz(void);
int make_atom(int,float,float,float,char *);
void make_bond(int,int,int);
void metalcheck(char *, int *, float *,float *,int *);
void metal(void);
void gettoken(void);
void write_gmmx(void);
void gmmx2(void);
void end_calculation(void);
double energy(void);
int setup_calculation(void);
float get_heat_of_formation(void);
void jacobi(int,int,double **,double *,double **,double *,double *);
