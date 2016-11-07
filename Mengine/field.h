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

#ifndef FIELD_H
#define FIELD_H

struct t_field {
  char    name[30];
  int     type;
  float   bondunit, bond_cubic, bond_quartic;
  float   angleunit, angle_cubic, angle_quartic, angle_pentic, angle_sextic;
  float   str_bndunit, ang_angunit, str_torunit, torsionunit;
  char    vdwtype[20], radiusrule[20], radiustype[20], radiussize[20], epsrule[20];
  float   vdw_14scale, a_expterm, b_expterm, c_expterm;
  float   chg_14scale, dielectric;
  float pos_12scale,pos_13scale;
        } field;
 #endif
void set_field(void);
void potoff(void);
