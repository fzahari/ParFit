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

#ifndef POT_H
#define POT_H

struct t_pot {
        int use_bond, use_angle, use_strbnd, use_urey, use_angang;
        int use_opbend, use_improp, use_imptor, use_tors, use_strtor;
        int use_tortor, use_vdw, use_lj, use_buck, use_hal, use_gauss;
        int use_charge, use_bufcharge, use_chrgdpl, use_dipole, use_polar, use_solv;
        int use_geom, use_extra, use_picalc, use_hbond, use_coordb, use_opbend_wilson;
        int use_bounds, use_image, use_replica, use_deform, use_ewald, use_highcoord;
        }  pot;

#endif

