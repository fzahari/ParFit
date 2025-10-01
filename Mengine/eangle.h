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

void eangle2a(int);
void eangle2a_fixed(int);
void eangle2b(int);
void eangle3(int, double (*)[3]);
void eangle(void);
void eangle1(void);
void eangle2(int);
//double acos(double);
int is_ring43(int, int, int);
int is_ring53(int, int, int);
int is_ring63(int, int, int);
int is_delocalbond(int, int);
int is_TS(int, int, int);
int isbond(int,int);
int find_bond(int,int);
void transtheta(int *, float *, float *, char *, char *, char *, char *);
void ts_four(int *, float *, float *,char *, char *,char *,char *);
long ipack(int ,int ,int ,int);
long kang(int,int,int);
void numeral(int,char *,int);
