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

// fix distantce data
#ifndef FIX_H
#define FIX_H

#define MAXDIS 300

struct t_fixdis {
        int FixBond;
        int nfxstr, ifxstr[MAXDIS][2];
        float fxstrd[MAXDIS], fxstrc[MAXDIS];
        int fixatm;
        }       fixdis;

// fix angle data
struct t_fixangle {
        int nfxangle, ifxangle[10][3];
        float anat[10], acon[10];
        } fixangle;
        
// fix torsion data
struct t_fxtor {
        int nfxtor, iatom[10][4];
        int start_ang[10], curr_ang[10], final_ang[10];
        int step[10];
        } fxtor;
#endif


int find_bond(int,int);
int find_fixed_bond(int,int);
void egeom(void);
void egeom1(void);
void egeom2(int);
