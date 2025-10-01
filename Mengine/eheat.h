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

#ifndef EHEAT_P
#define EHEAT_P

struct t_ehpara {
        int   neheat, nevhf;
        int nn[220];
        float ee[220], estr[220], et[220];   // bond data
        char cc[220][21], mm[155][21];
        int   cent[155];
        char  kheat[155][13];
        float ss[155], ssless[155];          // structural data
        } ehpara;

struct t_epiheat {
        int  npihf;
        long int nnp[125];
        float eep[125], aa[125], bb[125], ccc[125];
        char ccp[125][21];} epiheat;

struct t_piheat {
        float tepi, tesig, ecpi, etpi;
        }       piheat;

struct t_elowht {
        float epi;
        int mult;
        int uhf;
        }       elowht;

#endif

double dihdrl(int,int,int,int);
void numeral(int, char *, int);
void eheatpi(float *,float *,float *,float *,int *);
float get_heat_of_formation(void);
void eheat(void);
