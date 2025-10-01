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

#ifndef ANGK_H
#define ANGK_H

struct t_angk1 {
        int use_ang3, use_ang4, use_ang5;
        int nang, nang3, nang4, nang5;
        int ndel, ndel3, ndel4;
        char  ktype[MAXANGCONST][10],ktype3[MAXANG3CONST][10],ktype4[MAXANG4CONST][10],ktype5[MAXANG5CONST][10];
        char  kdel[MAXANGDEL][10],kdel3[MAXANG3DEL][10],kdel4[MAXANG4DEL][10];
        float con[MAXANGCONST], ang[MAXANGCONST][3];
        float con3[MAXANG3CONST], ang3[MAXANG3CONST][3];        
        float con4[MAXANG4CONST], ang4[MAXANG4CONST][3];        
        float con5[MAXANG5CONST], ang5[MAXANG5CONST][3];
        float condel[MAXANGDEL],  angdel[MAXANGDEL][3];      
        float condel3[MAXANG3DEL], angdel3[MAXANG3DEL][3];      
        float condel4[MAXANG4DEL], angdel4[MAXANG4DEL][3];      
         } angk1;

struct t_angkp1 {
        int  npiang;
        char kpa[MAXPIANGCONST][10];
        float pacon[MAXPIANGCONST], panat[MAXPIANGCONST][3], piobp[MAXPIANGCONST];
        } angkp1;

struct t_angf {
        int use_angf, nfang;
        char kftype[MAXANGCONST][10];
        float fcon[MAXANGCONST],fc0[MAXANGCONST],fc1[MAXANGCONST],fc2[MAXANGCONST];
        float fc3[MAXANGCONST],fc4[MAXANGCONST],fc5[MAXANGCONST],fc6[MAXANGCONST];
        } angf;

#endif
