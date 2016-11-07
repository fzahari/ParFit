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

void message_alert(char *,char *);
void find_rings1(int *, int **);
void del_atom(int, int **);
void log_ring1(int, int *);
int icompare(int, int *, int *);
int in_ring(int,int);
void add_all_rings(int,int *);
void trim(int ,int *,int *);
void check_N2(int,int *,int *);
void getUnion(int *, int *,int *);
int getRing(int,int,int *, int **,int *);
void delete_bond(int,int, int **);
void breakBond(int,int **);
void repack(int,int **);
void addq(int,int *,int);
int deleteq(int *,int);
void add_pair(int ,int ,int,int);
void find_rings(void);
int get_intersection(int *, int *);
