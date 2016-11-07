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

double energy(void);
void tncg(int,int,int *,double *,double *, double, double (*)(),  void (*)());
void mqn(int , int, int *,double *,double *, double, double (*)() );
double minimize1(double *, double *);
void search(int,double *,double *,double *,double *,double,double *,int *,double (*)(),int *);
double newton1(double *, double *);
void newton2(int, double *,double *,int *, int *, int *, double *);
void hessian(double *, int *, int *, int *,double *);
void gradient(void);
void piseq(int, int);
void piden(void);
void inesc(char *);
void update(void);
void hdel(int);
void hadd(void);
void minimize(void);
double minimiz1(double *, double *);
