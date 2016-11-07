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

FILE * fopen_path ( char * , char * , char * ) ;
void InitialTransform(void);
void zero_data(void);
void read_datafiles(char *);
void initialize_pcmodel(void);
void fixdisreset(void);
void coordreset(void);
void ddrivereset(void);
void hbondreset(void);
void pireset(void);
void generate_bonds(void);
void set_field(void);
void fixangle_reset(void);
void message_alert(char *, char *);
int strmessage_alert(char *);
void gettoken(void);
void free_dotmem(void);
void init_pov(void);
void reset_fixtype(void);
void set_window_title(void);
void remove_file(char *,char *);
void reset_calc_parameters(void);
void reset_atom_data(void);
void initialize(void);
void resetsubstrmem(void);
