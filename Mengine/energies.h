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

#ifndef ENERGIES_H
#define ENERGIES_H

struct  t_energies {
        double total, estr, ebend, estrbnd, e14, evdw, etor, eu, eopb,
               eangang, estrtor, ehbond, efix, eimprop, eimptors, eurey, esolv; } energies;

struct t_virial {
        double virx, viry, virz;
         } virial;
 #endif
