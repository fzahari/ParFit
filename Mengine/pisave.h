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

#ifndef PISAVE_H
#define PISAVE_H

struct t_pibondsave {
        int *ibsave, *kbsave, *nbsave;
        float *bksave, *blsave, *pksave, *plsave, *plsave2;
        }       pibondsave;

struct t_pitorsave {
        int *itsave, *ktsave, *ncsave;
        float    *tcsave;
        }       pitorsave;

struct t_pisave {
        float **fjsave, **fksave, **gjsave,**gksave;
        }       pisave;

struct t_mtape {
        float **vmtape, **vbtape, **pensave;
        }       mtape;

#endif
