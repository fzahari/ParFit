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

#ifndef DERIV_H
#define DERIV_H

/* EXTERN struct t_deriv {
        double d1[MAXATOM][3],
               deb[MAXATOM][3],  dea[MAXATOM][3],  destb[MAXATOM][3],deopb[MAXATOM][3],
               detor[MAXATOM][3],de14[MAXATOM][3], devdw[MAXATOM][3],deqq[MAXATOM][3],
               deaa[MAXATOM][3], destor[MAXATOM][3], dehb[MAXATOM][3], deimprop[MAXATOM][3];
       } deriv;   */

struct t_deriv {
        double **d1,
               **deb,  **dea,  **destb, **deopb, **detor,**de14, **devdw, **deqq,
               **deaa, **destor, **dehb, **deimprop, **deub, **desolv, *drb;
        } deriv;  
        
#endif
