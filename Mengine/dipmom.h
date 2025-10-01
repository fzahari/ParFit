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

#ifndef DIPMOM_H
#define DIPMOM_H

struct t_dipolemom {
        double total, xdipole, ydipole, zdipole;
       }  dipolemom;
struct t_dmomv {
        float xn,yn,zn,xp,yp,zp;
        }       dmomv;
struct t_quadmom {
        double xx,xy,xz,yx,yy,yz,zx,zy,zz;
       } quadmom;
 #endif

void charge_dipole(void);
void dipole_dipole(void);
void vibcharge_dipole(double **);
void vibdipole_dipole(double **);
