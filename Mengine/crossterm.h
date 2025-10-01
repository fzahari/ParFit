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

#ifndef CROSSTERM_H
#define CROSSTERM_H

struct t_crossterm_k {
        int nangang, nstrbnd, nstrtor;
        int ang_ang[MAXAA], stbnindex[MAXSTBN] ;
        char str_tor[MAXSTRTOR][7], stbn[MAXSTBN][10];
        float aacon[MAXAA][3], stbncon[MAXSTBN][3], str_torcon[MAXSTRTOR];
        } crossterm_k;

#endif
