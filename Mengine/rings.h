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

#ifndef RINGS_H
#define RINGS_H

struct t_rings {
                int nring3, **r13; //[MAXATOM][3];
                int nring4, **r14; //[MAXATOM][4];
                int nring5, **r15; //[MAXATOM][5];
                int nring6, **r16; //[MAXATOM/2][6];
                } rings;

#endif

int isbond(int, int);
int isangle(int,int);
int is_ring31(int);
int is_ring41(int);
int is_ring42(int,int);
int is_ring43(int, int, int);
int is_ring44(int *);
int is_ring51(int);
int is_ring52(int, int);
int is_ring53(int, int, int);
int is_ring54(int, int, int, int);
int is_ring55(int *);
int is_ring61(int);
int is_ring62(int, int);
int is_ring63(int,int,int);
int is_ring66(int *);
void ksort(int, int *);
void message_alert(char *, char *);
int is_cyclo6(int, int *);
int find_rsize(int, int);
void get_rsize(int,int,int,int *);
int aromatic_5(int,int *);
int aromatic_6(int *);
void get_rings(void);
void allocate_rings(void);
void free_rings(void);
int icompare(int, int *, int *);
int has_double_bond(int);
int check_connectivity(int ia);
int check_connectivity_yatom(int ia);
int get_bondorder(int,int);
