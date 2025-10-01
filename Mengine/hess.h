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

#ifndef HESS_H
#define HESS_H

double hesscut;
    
struct t_hess {
    float **hessx, **hessy, **hessz;
    } hess; 
   
#endif 

void piseq(int, int);
void hessian(double *,int *, int *, int *,double *);
void ebond2(int);
void eangle2(int);
void estrbnd2(int);
void eangang2(int);
void eopbend2(int);
void etorsion2(int);
void estrtor2(int);
void ecoord_bond2(int);
void echarge2(int);
void edipole2(int);
void ehal2(int);
void ebufcharge2(int);
void eopbend_wilson2(int);      
void esolv2(int);
void ehigh_coord2(int);
void message_alert(char *,char *);


