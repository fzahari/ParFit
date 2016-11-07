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

int  make_atom(int , float , float , float,char * );
void make_bond(int , int , int );
void generate_mm2lists(int *, int *, int **, int **);
void gettoken(void);
void mm3type(void);
void message_alert(char *, char *);
void set_atomdata(int,int,int,int,int,int);
int mm3_mmxtype(int);
void setbondorder(int,int);
void set_atomtypes(int);
void hdel(int);
void hadd(void);
FILE * fopen_path ( char * , char * , char * ) ;
int ReadMMInt(char *,int,int);
double ReadMMFloat(char *,int,int);
void FetchRecord(FILE *, char *);
void write_mm3(void);
void mm32mod(int);

