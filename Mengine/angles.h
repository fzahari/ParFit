# Copyright (C) 2004-2016 Kevin E. Gilbert, dba Serena Software
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Kevin E. Gilbert
# Serena Software
# Box 3076
# Bloomington, IN 47402-3076
#
# gilbert@serenasoft.com

#ifndef ANGLES_H
#define ANGLES_H

#define HARMONIC  1
#define FOURIER   2

#define MAXMETALANG 2000

struct t_angles {
        int nang, i13[MAXANG][4], index[MAXANG];
        int nopb, angin[MAXANG],angtype[MAXANG];
        float acon[MAXANG], anat[MAXANG], copb[MAXANG];
        float fcon[MAXANG],c0[MAXANG],c1[MAXANG],c2[MAXANG],c3[MAXANG],c4[MAXANG],c5[MAXANG],c6[MAXANG];
        } angles;

struct t_high_coord {
        int ncoord, i13[MAXMETALANG][3];
        } high_coord;
#endif

void message_alert(char *, char *);
void get_angles(void);
int isangle(int, int);
