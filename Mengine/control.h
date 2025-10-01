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

#ifndef CONTROL_H
#define CONTROL_H

struct t_control {
        int mode, nconfout,ffield,iprint,jprint,addconst;
        int optmet,electro,getcharge,iuv,pipl,gbsa;
        int inftype,outftype;
        char infile[128],outfile[128],addparfile[128],errorfile[128];
        }control;

#endif
