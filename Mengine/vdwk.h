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

#ifndef VDWK_H
#define VDWK_H

struct t_vdw1 {
        int  nvdw;
        float rad[MAXVDWCONST], eps[MAXVDWCONST];
        int lpd[MAXVDWCONST], ihtyp[MAXVDWCONST], ihdon[MAXVDWCONST];
        float alpha[MAXVDWCONST],n[MAXVDWCONST],a[MAXVDWCONST],g[MAXVDWCONST];
        char da[MAXVDWCONST][2];
        } vdw1;

struct t_vdwpr_k {
        int nvdwpr;
        int  ia1[MAXBONDCONST],ia2[MAXBONDCONST];
        char kv[MAXBONDCONST][7];
        float radius[MAXBONDCONST], eps[MAXBONDCONST];
        } vdwpr_k;

#endif
