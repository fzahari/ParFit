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

#include "pcwin.h"
#include "pcmod.h"
#include "pot.h"
#include "field.h"
#include "cutoffs.h"
#include "energies.h"
#include "utility.h"
#include "fix.h"
#include "optimize.h"
#include "files.h"
#include "minim_control.h"
#include "minim_values.h"
#include "pistuf.h"
#include "units.h"
#include "pcmfile.h"
#include "dipmom.h"
#include "user.h"
#include "initiali.h"
#include "tsbond.h"

#define NL gettoken(); if (pcmfile.head == 1000) goto L_30;
#define ALPHABETIC      2
#define NUMERIC 1

// =======================================
void message_alert(char *astring, char *title)
{
//    MessageBox(NULL, astring, title, MB_ICONEXCLAMATION);
      printf("%s\n",astring);
}
// ===================================
void initialize_pcmodel()
{
    char string[256];
        strcpy(pcwindir,"");
        pcmoutfile = fopen_path(pcwindir,"pcmod.out","w");
        user.dielec = FALSE;
        units.dielec = 1.0;
        dipolemom.total = 0.0;
        use_external_chrg = FALSE;
        use_gast_chrg = FALSE;
        zero_data();
        strcpy(string,"mmxconst.prm");
        read_datafiles(string);  
        set_field();
        LPTYPE = 20;
        if (field.type == MMX ||field.type == MM2 ||field.type == MM3  )
          LPTYPE = 20;
        else
          LPTYPE = 0;
        Openbox.ftype = FTYPE_PCM;
        Savebox.ftype = FTYPE_PCM;
        default_intype = MMX;
        default_outtype = MMX;
        strcpy(Openbox.path,pcwindir);
        strcpy(Savebox.path,pcwindir);
        minim_values.ndc = 4;
        reset_calc_parameters();
        natom = 0;
        bonds.numbonds = -1;
/*  minimizer control */
        minim_control.method = 3;
        minim_control.field = MMX;
        minim_control.added_const = FALSE;
// printout
        minim_values.iprint = FALSE;
        pistuf.jprint = FALSE;
/*  cutoffs    */
        cutoffs.vdwcut = 10.0;
        cutoffs.pmecut = 9.0;
        cutoffs.chrgcut = 100.0;
        cutoffs.dipcut = 8.0;
/*  solvation stuff  */
        pot.use_bounds = FALSE;
        pot.use_image = FALSE;
        pot.use_deform = FALSE;
        pot.use_solv = FALSE;
}

void initialize(void)
{
    natom = 0;
    bonds.numbonds = -1;
    reset_calc_parameters();
    pot.use_bounds = FALSE;
    pot.use_image = FALSE;
    pot.use_deform = FALSE;
    dipolemom.total = 0.0;
// optimization flags
    optimize_data.param_avail = FALSE;
    optimize_data.converge = FALSE;
    optimize_data.initial_energy = 0.0;
    optimize_data.final_energy = 0.0;
    optimize_data.initial_heat = 0.0;
    optimize_data.final_heat = 0.0;
//
    reset_atom_data();
}
/* ============================================== */
void reset_calc_parameters(void)
{
   flags.noh = False;
   flags.nohyd = False;
   flags.immx = True;

/* default to hydrogen bonding on  */
//   minim_values.ndc = 4;
   minim_values.nconst= 0;
/* pi stuff */
   pistuf.nplane = 1;
   pistuf.iuv = 1;
   pistuf.ihuck = 0;
}
/* =============================================== */
void pireset()
{
        long int i,mask;
        mask = 1L << PI_MASK;
        for (i = 1; i <= natom; i++)
                atom[i].flags &= ~mask;
}
/* =============================================== */
void reset_fixtype()
{
        long int i,mask;
        mask = 1L << NO_RETYPE;
        for (i = 1; i <= natom; i++)
                atom[i].flags &= ~mask;
}
/* ------------------------- */
void hbondreset()
{
}
/* ------------------------- */
void resetsubstrmem()
{
}
/* ------------------------- */
void coordreset()
{
        long int i,j, mask4,mask5,mask6,mask7,mask8;

        mask4 = 1L << METCOORD_MASK;
        mask5 = 1L << SATMET_MASK;
        mask6 = 1L << GT18e_MASK;
        mask7 = 1L << LOWSPIN_MASK;
        mask8 = 1L << SQPLAN_MASK;

        for (i=1; i <= natom; i++)
        {
                atom[i].flags &= ~mask4;
                atom[i].flags &= ~mask5;
                atom[i].flags &= ~mask6;
                atom[i].flags &= ~mask7;
                atom[i].flags &= ~mask8;

                for (j=0; j < MAXIAT; j++)
                {
                        if (atom[i].bo[j] == 9)
                        {
                                atom[i].bo[j] = 0;
                                atom[i].iat[j] = 0;
                        }
                }
        }
        generate_bonds();
}               
/* ------------------------- */
void fixdisreset()
{
}
/* ============================================== */
void fixangle_reset()
{
}
/* ============================================== */       
void reset_atom_data(void)
{
        int i, j;

    for (i = 0; i < MAXATOM; i++) 
    {
        atom[i].use = 0;
        atom[i].x = 0.0F; atom[i].y = 0.0F; atom[i].z = 0.0F;
        atom[i].color = 0;
        atom[i].chrg_color = 7;
        atom[i].charge = 0.0F;
        atom[i].energy = 0.0F;
        atom[i].atomnum = 0;
        atom[i].serno = 0;
        atom[i].molecule = 0;
        atom[i].residue = 0;
        atom[i].formal_charge = 0;
        atom[i].atomwt = 0.0;
        atom[i].type = 0;
        atom[i].mmx_type = 0;
        atom[i].mm3_type = 0;
        atom[i].mmff_type = 0;
        atom[i].tclass = 0;
        *atom[i].name = '\0';
        atom[i].flags = 0;
        for (j = 0; j < MAXIAT; j++) 
        {
           atom[i].bo[j] = 0;
           atom[i].iat[j] = 0;
        }
     }
}
// =======================================
#ifdef PCM_WIN
int strcasecmp(char *str1,char *str2)
{
     char c1,c2;
     while(1)
     {
         c1 = tolower(*str1++);
         c2 = tolower(*str2++);
         if (c1 < c2) return -1;
         if (c1 > c2) return 1;
         if (c1 == 0) return 0;
     }
}
// ========================================
int strncasecmp(char *s1,char *s2,int n)
{
    int i;
    char c1, c2;
    for (i=0; i<n; i++)
    {
        c1 = tolower(*s1++);
        c2 = tolower(*s2++);
        if (c1 < c2) return -1;
        if (c1 > c2) return 1;
        if (!c1) return 0;
    }
    return 0;
} 
#endif
// void itoa(int n,char *s, int radix)
//{
//    sprintf(s,"%d",n);
//}
// ========================================
void remove_file(char *path,char *name)
{
    char tempname[255];
    int ix;

    strcpy(tempname,"");
    if ( (ix = strlen(path)) != 0)
    {
        strcpy(tempname,path);
        strcat(tempname,"\\");
        strcat(tempname,name);
    } else
    {
        strcpy(tempname,name);
    }
    remove(tempname);
}
