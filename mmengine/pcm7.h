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

void set_missing_constants(int);
int get_missing_constants(void);
float get_heat_of_formation(void);
FILE * fopen_path ( char * , char * , char * ) ;
void message_alert(char *,char *);
void get_added_const(void);
void assign_gast_charge();
void lattice(void);
void bounds(void);
void hadd(void);
void pcm7_min(void);
int iscoord_bond(int, int);
void get_coordbonds(void);
int isbond(int, int);
int isangle(int,int);
int ishbond(int, int);
int istorsion(int, int);
void get_hbond(void);
void mark_hbond(void);
void reset_atom_data(void);
int initialise(void);
void get_bonds(void);
void get_angles(void);
void get_torsions(void);
void get_rings(void);
void set_field(void);
double energy(void);
double print_energy(void);
void pcmfout(int);
void orbital(void);
void free_orbit(void);
void set_active(void);
void estrtor(void);
void echarge(void);
void echrg_dipole(void);
void echrg_dipole1(void);
void ebuck(void);
void ebuckmm3(void);
void estrbnd(void);
void edipole(void);
void etorsion(void);
void eangang(void);
void eangle(void);
void ebond(void);
void eopbend(void);
void ehbond(void);
void ecoord_bond(void);
void ebuck_charge(void);
void ehigh_coord(void);

void estrtor1(void);
void echarge1(void);
void ebuck1(void);
void ebuckmm31(void);
void estrbnd1(void);
void edipole1(void);
void etorsion1(void);
void eangang1(void);
void eangle1(void);
void ebond1(void);
void eopbend1(void);
void ehbond1(void);
void elj1(void);
void ecoord_bond1(void);
void ebuck_charge1(void);
void ehigh_coord1(void);

void kangle(void);
void korbit(void);
void ktorsion(void);
void kdipole(void);
void kstrbnd(void);
void kcharge(void);
void ksolv(void);
void kangang(void);
void kopbend(void);
void kstrtor(void);
void khbond(void);
void piseq(int, int);
void kcoord_bonds(void);
void kvdw(void);
int  kbond(void);

void get_memory(void);
void free_memory(void);
void gradient(void);
void minimize(void);
void dynamics(void);
void read_datafiles(char *);
void attach(void);
void eheat(void);
void pirite(void);
void charge_dipole(void);
void dipole_dipole(void);
void piden(void);
void hdel(int);
void set_atomtypes(int);
void type(void);
void generate_bonds(void);
void zero_data(void);
void ehal(void);
void ehal1(void);
void ebufcharge(void);
void ebufcharge1(void);
void eopbend_wilson(void);
void eopbend_wilson1(void);
void hbondreset(void);
void vibrate(void);
void pireset(void);
void adjust_mmfftypes(void);
int setup_calculation(void);
void end_calculation(void);
void egeom1(void);
void esolv1(void);
void egeom(void);
void esolv(void);
void xlogp(float *);
void mmxsub(int);
