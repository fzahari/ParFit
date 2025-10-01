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

#ifndef PIDATA_H
#define PIDATA_H

struct t_penint {
        float **pen, *ed;
        }       penint;
struct t_ed123 {
        int nscft;
        float *ed1, *ed2, *ed3, *ed4, *ed5;
        }       ed123;
struct t_densit {
        float **den, **dena, **denb;
        }       densit;
struct t_local {
        float *da, *adnmo, *adnmt, *d12, *d13, *fn;
        float *f22, *f27, *f29, *f1, *fa1;
        }  local;
struct t_uvnpi {
        int ifirst, kfirst;
        }       uvnpi;
struct t_pigroup {
        int nfrag;
        int fraglist[MAXORB],grlist[MAXORB],nao[MAXORB/4],nelecs[MAXORB/4];
        } pigroup;       

#endif
