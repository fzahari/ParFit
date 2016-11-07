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

/* utility routines */
extern float *vector( int, int);
extern  int *ivector( int, int);
extern  long int *ilvector( int, int);
extern double *dvector( int, int);
extern float **matrix( int,  int,  int,  int);
extern  int **imatrix( int,  int,  int,  int);
extern double **dmatrix( int,  int,  int,  int);
extern void free_vector(float *,  int,  int);
extern void free_ivector( int *,  int,  int);
extern void free_ilvector( long int *,  int,  int);
extern void free_dvector(double *,  int,  int);
extern void free_matrix(float **, int, int, int, int);
extern void free_imatrix( int **, int, int, int, int);
extern void free_dmatrix(double **, int, int, int, int);
