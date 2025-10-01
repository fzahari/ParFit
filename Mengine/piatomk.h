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

#ifndef PIATOM_K_H
#define PIATOM_K_H

struct t_piatomk {
        int   npiatom;
        int   kat[MAXPIATOM], qp[MAXPIATOM];
        float q[MAXPIATOM], ion[MAXPIATOM], emz[MAXPIATOM], zor[MAXPIATOM], zz[MAXPIATOM];
        float w1[MAXPIATOM], w2[MAXPIATOM], w3[MAXPIATOM];
        }  piatomk;

#endif
