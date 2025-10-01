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

#ifndef CHARGEK_H
#define CHARGEK_H

struct t_charge_k {
         int ncharge, nbndchrg, nbndchrgdel;
         int type[MAXATOMTYPE], btype[MAXBONDCONST], btypedel[MAXBONDCONST];
         float charge[MAXATOMTYPE], bcharge[MAXBONDCONST], formchrg[MAXATOMTYPE], bchargedel[MAXBONDCONST];
         float typechrg[MAXATOMTYPE];
         } charge_k;

#endif
