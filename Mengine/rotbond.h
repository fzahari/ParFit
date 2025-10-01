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

#ifndef ROTBOND_H
#define ROTBOND_H

struct t_rotbond {
    int nbonds, bond[200][2], bond_res[200];
    int incl_amide, incl_alkene;
    } rotbond;

#endif

void find_rings(void);
void find_bonds(void);
void hdel(int);
void hadd(void);
int check_bond(int,int,int *);
void message_alert(char *, char *);
void matident(float (*)[4]);
void matrotsc(float (*)[4], int*, float*, float*);
void matrotat(float (*)[4], int*, float*);
void mattrans(float (*)[4], float*, float*, float*);
void matriply(float (*)[4], float (*)[4], float (*)[4]);
void matxform(float (*)[4], int);
void rotb(int, int, int, int, float);
void chain(int, int, int);
double dihdrl(int,int,int,int);
