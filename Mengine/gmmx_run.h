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

int FetchRecord(FILE *, char *);
void sort_file(char *, char *);
void gettoken(void);
int is_bond(int,int);
double compute_rms(int, int **, double **, double **);
void fndumt( int, int **, double **, double ** );
int fast_compare(double **, double **);
void slow_compare(int, int **, double **, double **);
void message_alert(char *, char *);
void gradient(void);
void eigen(int norder,double (*)[3],double *, double (*)[3]);
int check_structure(void);
int check_id(void);  
double distance(int,int);
void rotb(int,int,int,int, float);
void make_bond(int , int , int );
void deletebond(int, int);
int setup_calculation(void);
void end_calculation(void);
void hster(int, int *);
int lepimer(int);
double dihdrl(int,int,int,int);
void mvatom(void);
void mvbond(void);
void minimize(void);
void initial_picalc(void);
int findcd(int);
int findtor(int);
int bad15(int);
int close_contact(int);
void gmmx2(void);
int check_stop(void);
void hbondreset(void);
void setup_output(char *);
void pimarkselection(void);
void pireset(void);
void pcmfout(int);
void center(int);
void pcmfin(int,int);
void initialize(void);
void generate_bonds(void);
int  bstereo(int);
void type(void);
void mmxsub(int);
void eheat(void);
void charge_dipole(void);
void check_numfile(int);
void dipole_dipole(void);
void set_active(void);
double energy(void);
void write_restart(void);
int  read_restart(void);
void rdfile(int);
void montc(void);
void include_min(void);
void jacobi(int,int,double **,double *,double **,double *,double *);
void quatfit(int,int **,double **,double **);
FILE * fopen_path ( char * , char * , char * ) ;
void read_gmmxinp(void);
void run_gmmx(void);
int rmsck(int, int *);
void write_gmmx(void);
void read_sdf(int,int);
float ran1(int *);
void find_rings(void);
void hcoord(void);
void initialize_gmmx(void);
int check_bond(int,int,int *);
