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

#ifndef TORSIONK_H
#define TORSIONK_H

struct t_allene {
    int nallene, ntor[10];
     } allene;
        
        
struct t_torkn1 {
        int  use_tor4, use_tor5;
        int  ntor, ntor4, ntor5, ntordel, torindex[MAXTORDEL];
        char  kv[MAXTORCONST][13], kv4[MAXTOR4CONST][13], kv5[MAXTOR5CONST][13];
        char  kvdel[MAXTORDEL][13];
        float tv1[MAXTORCONST], tv2[MAXTORCONST], tv3[MAXTORCONST];
        float tv4[MAXTORCONST], tv5[MAXTORCONST], tv6[MAXTORCONST];
        int   phase1[MAXTORCONST], phase2[MAXTORCONST], phase3[MAXTORCONST];
        int   phase4[MAXTORCONST], phase5[MAXTORCONST], phase6[MAXTORCONST];
        float tv41[MAXTOR4CONST], tv42[MAXTOR4CONST], tv43[MAXTOR4CONST];
        int   phase41[MAXTOR4CONST], phase42[MAXTOR4CONST], phase43[MAXTOR4CONST];
        float tv51[MAXTOR5CONST], tv52[MAXTOR5CONST], tv53[MAXTOR5CONST];
        int   phase51[MAXTOR5CONST], phase52[MAXTOR5CONST], phase53[MAXTOR5CONST];
        float tvdel1[MAXTORDEL], tvdel2[MAXTORDEL], tvdel3[MAXTORDEL];
        int   phasedel1[MAXTORDEL], phasedel2[MAXTORDEL], phasedel3[MAXTORDEL];
        } torkn1;

struct t_torknp {
        int npitor;
        char kv[MAXPITORCONST][13];
        int ph1[MAXPITORCONST], ph2[MAXPITORCONST], ph3[MAXPITORCONST];
        float   tv1[MAXPITORCONST], tv2[MAXPITORCONST], tv3[MAXPITORCONST];
        } torknp;

#endif
