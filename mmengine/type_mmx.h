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

void set_atomtype(int,int,int,int,int,int);
void set_atomtypes(int);
int is_ring31(int);
int is_ring41(int);
int is_ring51(int);
int is_cyclo5(int, int *);
int isbond(int,int);
int aromatic_5(int,int *);
int find_rsize(int,int);
void get_rsize(int,int,int, int *);
int icompare(int, int *, int *);
void adjust_mmfftypes();
void deletebond(int,int);
void pireset(void);

