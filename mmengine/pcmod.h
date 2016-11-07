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

#ifndef PCMOD_H
#define PCMOD_H

#ifndef True
#define True 1
#define False 0
#endif 

#include "pcmsiz.h"


int     LPTYPE;

#define radian 57.29577951308


#define         MMX             1
#define         MM2             2
#define         MM3             3
#define         MM4             4
#define         AMBER           5
#define         CHARMM          6
#define         MMFF94          7
#define         OPLS            8
#define         OPLSAA          9
#define         UNKNOWN         10

// File Information
#define FTYPE_MMX                   102
#define FTYPE_MM2                   103
#define FTYPE_MM3                   104
#define FTYPE_ALC                   106
#define FTYPE_SYBYL                 107
#define FTYPE_MMOD                  108
#define FTYPE_PCM                   109
#define FTYPE_XRA                   110
#define FTYPE_C3D                   111
#define FTYPE_MOL                   112
#define FTYPE_PDB                   115
#define FTYPE_CSD                   118
#define FTYPE_SDF                   123
#define FTYPE_TINKER                124
#define FTYPE_MOL2                  125

#define FTYPE_MOP                   105
#define FTYPE_ARC                   114
#define FTYPE_GAU                   113   // read gaussian output
#define FTYPE_GAUSOUT               119   // write gaussian job file
#define FTYPE_GAU_IRC               122
#define FTYPE_GAUSFCHK              128
#define FTYPE_PSGVBIN               116
#define FTYPE_PSGVBOUT              117
#define FTYPE_GAMES                 120
#define FTYPE_GAMESOUT              121
#define FTYPE_EHT                   126
#define FTYPE_HONDO                 127
#define FTYPE_HONDOPUN              131
#define FTYPE_TURBOMOLE             129
#define FTYPE_ADF                   130

#define FTYPE_SMILES                131
#define FTYPE_CHEMDRAW              132
#define FTYPE_XML                   133
#define FTYPE_CML                   134

#define SUB_MOVE                        0
#define SUB_HIDE                        1
#define SUB_MINIMIZE                    2

//  flags definitions
#define PI_MASK                 0
#define HBOND_MASK              1
#define AROMATIC_MASK           2
//  metal flags
#define METCOORD_MASK           3
#define SATMET_MASK             4
#define GT18e_MASK              5
#define LOWSPIN_MASK            6
#define SQPLAN_MASK             7

// type rules
#define NO_RETYPE               8
//  invisible
#define VIS_MASK                9
//  minimize 
#define MIN_MASK                10
// cpk surface
#define CPK_SURF                11
//  dotsurf
#define DOT_SURF                12
// Nterm, CNterm, Oterm, COterm, DUMMY
#define NTERM                   13
#define CNTERM                  14
#define OTERM                   15
#define COTERM                  16
#define DUMMY                   17
#define P5                      18
#define P3                      19
// Ring Size
#define RING3                   20
#define RING4                   21
#define RING5                   22
#define RING6                   23

/*  PCMODEL specific definitions */
#ifndef ATOMTYPE
typedef struct {
                        double x,y,z;
                        int type;
                        int tclass;
                        int mmx_type;
                        int mm3_type;
                        int amber_type;
                        int mmff_type;
                        int charm_type;
                        int opls_type;
                        int atomnum;
                        int serno;
                        int molecule;
                        int residue;
                        int biotype;
                        double atomwt;
                        float energy;
                        int use;
                        int color;
                        int chrg_color;
                        int iat[MAXIAT];
                        int bo[MAXIAT];
                        char name[3];
                        double charge;
                        float formal_charge;
                        float sigma_charge;
                        float radius;
                        float vdw_radius;
                        long int flags;
                        long int substr[MAXSSCLASS];
                        } ATOMTYPE;
#endif

struct {
                int numbonds;
                int ia1[MAXATOM], ia2[MAXATOM], bondorder[MAXATOM];
                } bonds;

struct {
                int   istereo;
                int   noh;
                int   immx;
                int   minimized;
                int   modified;
                float dielc;
                float avleg;
                int   rescale;
                int   nohyd; } flags;


ATOMTYPE         atom[MAXATOM];
int              natom;
char             Struct_Title[100];
FILE             *pcmoutfile, *errfile;
char             pcwindir[80];
int              Mopac_Charges;
int              use_external_chrg;
int              use_gast_chrg;
int              default_intype, default_outtype;

#endif
