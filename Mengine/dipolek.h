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

#ifndef DIPOLEK_H
#define DIPOLEK_H

struct  t_dipole_k {
        int ndipole,ndipole3,ndipole4,ndipole5;
        char   kb[MAXBONDCONST][7], kb3[MAXBOND3CONST][7], kb4[MAXBOND4CONST][7], kb5[MAXBOND5CONST][7];
        float bmom[MAXBONDCONST],bmom3[MAXBOND3CONST],bmom4[MAXBOND4CONST],bmom5[MAXBOND5CONST];
         } dipole_k;

#endif
