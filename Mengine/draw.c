// Copyright (C) 2004-2016 Kevin E. Gilbert, dba Serena Software
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Kevin E. Gilbert
// Serena Software
// Box 3076
// Bloomington, IN 47402-3076
//
// gilbert@serenasoft.com

#include "pcwin.h"
#include "pcmod.h"
#include "field.h"
#include "atom_k.h"
#include "pot.h"
#include "energies.h"
#include "fix.h"
#include "minim_control.h"
#include "minim_values.h"
#include "units.h"
#include "draw.h"
#include "elements.h"


// ===============================================
void set_atomdata(int ia, int mmxtype, int mm3type, int mmfftype,int ambertype,int oplstype)
{
    int type = 0;
    
        if (mmxtype > 0)
           atom[ia].mmx_type = mmxtype;
        if (mm3type > 0)
           atom[ia].mm3_type = mm3type;
        if (mmfftype > 0)
           atom[ia].mmff_type = mmfftype;
           
        if (field.type == MMX || field.type == MM2)
        {
           type = mmxtype;
        } else if (field.type == MM3)
        {
           type = mm3type;
        } else if (field.type == MMFF94)
        {
            type = mmfftype;
        }

        if ( type < 300) 
        {
           atom[ia].tclass = atom_k.tclass[mmxtype];
           atom[ia].atomnum = atom_k.number[mmxtype];
           atom[ia].atomwt = atom_k.weight[mmxtype];
           strcpy(atom[ia].name, atom_k.symbol[mmxtype]);
        }        
}
// ===============================================
void set_atomtype(int ia, int mmxtype, int mm3type, int mmfftype,int ambertype,int oplstype)
{
        if (mmxtype > 0)
           atom[ia].mmx_type = mmxtype;
        if (mm3type > 0)
           atom[ia].mm3_type = mm3type;
        if (mmfftype > 0)
           atom[ia].mmff_type = mmfftype;
           
}
// ======================================================
int make_atom(int type, float x, float y, float z,char *name)
{
     int i,iz;
     natom++;
     if (natom >= MAXATOM - 10)
     {
         natom--;
         message_alert("Max atoms exceeded in makeatom","PCMODEL Error");
         return -1;
     }
        iz = strlen(name);
        if (iz > 0) // got an atom name use it 
        {
            for (i=1; i <= MAXATOMTYPE; i++)
            {
                if (strcmp(name,atom_k.symbol[i]) == 0) // names match
                {
                    strcpy(atom[natom].name,name);
                    atom[natom].tclass = atom_k.tclass[i];
                    atom[natom].atomwt = atom_k.weight[i];
                    atom[natom].atomnum = atom_k.number[i];
                    atom[natom].type = atom_k.type[i];
                    if (type == 0) type = atom[natom].type;
                    if (field.type == MMX || field.type == MM2)
                          atom[natom].mmx_type = type;
                    else if (field.type == MM3)
                          atom[natom].mm3_type = type;
                    else if (field.type == MMFF94)
                          atom[natom].mmff_type = type;
                    if (type != 0)
                    {
                        atom[natom].type = type;
                        atom[natom].tclass = atom_k.tclass[type];
                    }
                    atom[natom].x = x; atom[natom].y = y; atom[natom].z = z;
                    atom[natom].vdw_radius = Elements[atom[natom].atomnum-1].covradius;
                    return natom;
                } 
            }
           for( i = 0; i < 103; i++ )
           {
              if( strcasecmp(Elements[i].symbol,name) == 0 )
              {
                  strcpy(atom[natom].name,Elements[i].symbol);
                  atom[natom].atomwt = Elements[i].weight;
                  atom[natom].atomnum = Elements[i].atomnum;
                  atom[natom].type = Elements[i].type;    // set generic type
                  atom[natom].tclass = Elements[i].type;    // set generic type
                  atom[natom].mmx_type = Elements[i].type;    // set generic type
                  if (type > 0) 
                  {
                       atom[natom].type = type; // overwrite type if type is passed in
                       if (field.type == MMX || field.type == MM2)
                          atom[natom].mmx_type = type;
                       else if (field.type == MM3)
                          atom[natom].mm3_type = type;
                       else if (field.type == MMFF94)
                          atom[natom].mmff_type = type;
                  }
                  if (type != 0) atom[natom].type = type;
                  atom[natom].x = x; atom[natom].y = y; atom[natom].z = z;
                  atom[natom].vdw_radius = Elements[i].covradius;
                  return natom;
              }
            }
//  got a name but it is not recognized set default values
            strcpy(atom[natom].name,name);
            atom[natom].atomwt = 1.0;
            atom[natom].atomnum = 1;
            if (type == 20)
            {
                atom[natom].atomwt = 0.0;
                atom[natom].atomnum = 0;
            }
            atom[natom].vdw_radius = 2.0;
            atom[natom].x = x; atom[natom].y = y; atom[natom].z = z;
            if (type > 0) 
            {
                atom[natom].type = type;
                if (field.type == MMX || field.type == MM2)
                    atom[natom].mmx_type = type;
                else if (field.type == MM3)
                    atom[natom].mm3_type = type;
                else if (field.type == MMFF94)
                    atom[natom].mmff_type = type;
             } else
            {
               atom[natom].mmx_type = 60;  // set type to Dummy atom
               atom[natom].mm3_type = 299;
            }
            return natom;
        } else if (type > 0)
        {
            atom[natom].type = type;
            if (field.type == MMX || field.type == MM2)
                atom[natom].mmx_type = type;
            else if (field.type == MM3)
                atom[natom].mm3_type = type;
            else if (field.type == MMFF94)
                atom[natom].mmff_type = type;
            atom[natom].tclass = atom_k.tclass[type];
            if (atom[natom].tclass == 0) atom[natom].tclass = atom[natom].type;
            atom[natom].atomnum = atom_k.number[type];
            atom[natom].atomwt = atom_k.weight[type];
            strcpy(atom[natom].name, atom_k.symbol[type]);
        
            atom[natom].x = x; atom[natom].y = y; atom[natom].z = z;
            return natom;
        } else
        {
           message_alert("Input atom does not have valid atomic symbol or atom type","Error");
           return -1;
        }
}
// ====================================================
void make_bond(int ia1, int ia2, int bo)
{
    int loop,loop1;

   if (ia1 < 1 || ia1 > MAXATOM) return;
   if (ia2 < 1 || ia2 > MAXATOM) return;
   /*   find if bond exists        */
   for (loop = 0; loop < MAXIAT; loop++)
   {
     if (atom[ia1].iat[loop] == ia2)
     {
        if (bo == 1)
           atom[ia1].bo[loop] += bo;
        else
           atom[ia1].bo[loop] = bo;
        
        for (loop1 = 0; loop1 < MAXIAT; loop1++)
        {
         if (atom[ia2].iat[loop1] == ia1)
         {
             if ( bo == 1)
             {
                atom[ia2].bo[loop1] += bo; /* bond exists return */
                return;
             } else
             {
                atom[ia2].bo[loop1] = bo; /* bond exists return */
                return;
             }
         }
        }
      }
   }
/*  bond does not exist create it               */     
   for (loop = 0; loop < MAXIAT; loop++) 
   {
      if (atom[ia1].iat[loop] == 0) 
      {
          atom[ia1].iat[loop] = ia2;
          atom[ia1].bo[loop] = bo;
          break;
      }
   }
   for (loop = 0; loop < MAXIAT; loop++) 
   {
      if (atom[ia2].iat[loop] == 0) 
      {
          atom[ia2].iat[loop] = ia1;
          atom[ia2].bo[loop] = bo;
          break;
      }
   }
}
// ===============================================
void deleteatom(int i)
{
   int j,k, ibo;
   atom[i].type = 0;
   atom[i].x = 0.0F;
   atom[i].y = 0.0F;
   atom[i].z = 0.0F;
   atom[i].flags = 0;
   for (j=0; j<MAXIAT; j++)
   {
      if (atom[i].iat[j] != 0 )
      {
         if (atom[atom[i].iat[j]].atomnum == 1 || atom[atom[i].iat[j]].atomnum == 0 )
         {
             atom[atom[i].iat[j]].type = 0;
             atom[atom[i].iat[j]].x = 0.0F;
             atom[atom[i].iat[j]].y = 0.0F;
             atom[atom[i].iat[j]].z = 0.0F;
             atom[atom[i].iat[j]].flags = 0;
         }
         ibo = atom[i].bo[j];
         for (k=1; k <= ibo; k++)
            deletebond(i,atom[i].iat[j]);
         atom[i].iat[j] = 0;
     }
   }
   natom--;
}
// ==================================================
void deletebond(int i, int j) 
{
   int i1,j1;
   for (i1=0; i1 < MAXIAT; i1++) 
   {
      if (atom[i].iat[i1] == j) 
      {
        if (atom[i].bo[i1] == 9)  /* if deleting a coordinated bond - delete it completely */
        {
            atom[i].bo[i1] = 0;
            atom[i].iat[i1] = 0;
        }
        atom[i].bo[i1] -= 1;
        if (atom[i].bo[i1] <= 0 )
            atom[i].iat[i1] = 0;
        break;
      }
   }
   for (j1=0; j1 < MAXIAT; j1++) 
   {
      if (atom[j].iat[j1] == i) 
      {
        if (atom[j].bo[j1] == 9)  /* if deleting a coordinated bond - delete it completely */
        {
            atom[j].bo[j1] = 0;
            atom[j].iat[j1] = 0;
        }
        atom[j].bo[j1] -= 1;
        if (atom[j].bo[j1] <= 0)
            atom[j].iat[j1] = 0;
        break;
      }
   }
}
/* --------------------------- */
void generate_bonds(void)
{
   int i,j;

   bonds.numbonds = 0;
   for (i= 0; i <= natom; i++)
   {
           bonds.ia1[i] = 0;
           bonds.ia2[i] = 0;
           bonds.bondorder[i] = 0;
   }

   for (i = 1; i <= natom; i++)
   {
      for (j = 0; j < MAXIAT; j++) 
      {
         if (atom[i].iat[j] > i && atom[i].iat[j] != 0) 
         {
            bonds.ia1[bonds.numbonds] = i;
            bonds.ia2[bonds.numbonds] = atom[i].iat[j];
            bonds.bondorder[bonds.numbonds] = atom[i].bo[j];
            bonds.numbonds++;
         }
     }
   }
}
