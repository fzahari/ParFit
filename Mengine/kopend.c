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

#include "angles.h"
#include "pot.h"
#include "field.h"
#include "ooplane.h"
#include "opbend.h"

void kopbend()
{
   int i, j, ia,ib,ic,id, jj, ierr;
   int ita,itb, itc,itd, itd_class;
   int it, it1;
   long mask5, mask6, mask8;
   char pa[4],pb[4],pc[4],pd[4],pdc[4],pt[13],pt1[13],pt2[13];
 
   mask5 = 1L << SATMET_MASK;
   mask6 = 1L << GT18e_MASK;
   mask8 = 1L << SQPLAN_MASK;

   if (ooplane_k.nopbend <= 0)
   {
       pot.use_opbend = FALSE;
       return;
   }
   angles.nopb = 0;  
   for (i=0; i < angles.nang; i++)
   {
      ia = angles.i13[i][0];
      ib = angles.i13[i][1];
      ic = angles.i13[i][2];
      ita = atom[ia].type;
      itb = atom[ib].type;
      itc = atom[ic].type;
      ierr = FALSE;
      if (angles.i13[i][3] != 0 )
      {
             // angle ia-ib-id-ic
             id = angles.i13[i][3];
             itd = atom[id].type;
             itd_class = atom[id].tclass;
             numeral(ita,pa,3);
             numeral(itb,pb,3);
             numeral(itc,pc,3);
             numeral(itd,pd,3);
             numeral(itd_class,pdc,3);
             if (field.type == MMX || field.type == MM3)
             {
                 strcpy(pt,"  0"); strcat(pt,pb); strcat(pt,pd); strcat(pt,"  0");
             }else if (field.type == MMFF94)
             {
               if (ita <= itc && ita <= itd)
               {
                 if (itc <= itd)
                 {
                   strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 } else
                 {
                   strcpy(pt,pa); strcat(pt,pb); strcat(pt,pd); strcat(pt,pc);
                 }
               } else if ( itc <= ita && itc <= itd)
               {
                 if (ita <= itd)
                 {
                    strcpy(pt,pc); strcat(pt,pb); strcat(pt,pa); strcat(pt,pd);
                 } else
                 {
                    strcpy(pt,pc); strcat(pt,pb); strcat(pt,pd); strcat(pt,pa);
                 }
               } else if (itd <= ita && itd <= itc)
               {
                 if ( ita <= itc)
                 {
                   strcpy(pt,pd); strcat(pt,pb); strcat(pt,pa); strcat(pt,pc);
                 } else
                 {
                   strcpy(pt,pd); strcat(pt,pb); strcat(pt,pc); strcat(pt,pa);
                 }
               }
             }
             strcpy(pt1,"  0"); strcat(pt1,pb); strcat(pt1,pdc); strcat(pt1,"  0");
             strcpy(pt2,"  0"); strcat(pt2,pb); strcat(pt2,"  0"); strcat(pt2,"  0");
             it = itb*100 + itd;
             it1 = itb*100;
             ierr = FALSE;
             for (j=0; j < ooplane_k.nopbend; j++)
             {
                if (strcmp(pt,ooplane_k.iopb[j]) == 0 || strcmp(pt1,ooplane_k.iopb[j]) == 0 || strcmp(pt2,ooplane_k.iopb[j]) == 0)
                {
                     angles.copb[i] = ooplane_k.copb[j];
                     angles.angin[i] = TRUE;
                     angles.nopb++;
                     ierr = TRUE;
                     break;
                }
             }
      }
   }
//  look for sq planar metals
   for (i=1; i <= natom; i++)
   {
       if (atom[i].type >= 300)
       {
           jj = 0;
           for (j=0; j < MAXIAT; j++)
           {
               if (atom[i].iat[j] != 0)
                jj++;
           }
           if (jj == 4)
           {
               strcpy(pt1,"  0"); strcat(pt1,"300"); strcat(pt1,"  0"); strcat(pt1,"  0");
               if ( (atom[i].flags & mask5) && (atom[i].flags & mask6) && (atom[i].flags & mask8) )
               {
                  for (j=0; j < ooplane_k.nopbend; j++)
                  {
                     if (strcmp(pt,ooplane_k.iopb[j]) == 0 || strcmp(pt1,ooplane_k.iopb[j]) == 0)
                     {
                        angles.i13[angles.nang][0] = atom[i].iat[0];
                        angles.i13[angles.nang][1] = i;
                        angles.i13[angles.nang][2] = atom[i].iat[1];                        
                        angles.i13[angles.nang][3] = atom[i].iat[2];                        
                        angles.copb[angles.nang] = ooplane_k.copb[j];
                        angles.angin[angles.nang] = TRUE;
                        angles.nang++;
                        angles.nopb++;
                        angles.i13[angles.nang][0] = atom[i].iat[0];
                        angles.i13[angles.nang][1] = i;
                        angles.i13[angles.nang][2] = atom[i].iat[2];                        
                        angles.i13[angles.nang][3] = atom[i].iat[1];                        
                        angles.copb[angles.nang] = ooplane_k.copb[j];
                        angles.angin[angles.nang] = TRUE;
                        angles.nang++;
                        angles.nopb++;
                        angles.i13[angles.nang][0] = atom[i].iat[0];
                        angles.i13[angles.nang][1] = i;
                        angles.i13[angles.nang][2] = atom[i].iat[3];                        
                        angles.i13[angles.nang][3] = atom[i].iat[1];                        
                        angles.copb[angles.nang] = ooplane_k.copb[j];
                        angles.angin[angles.nang] = TRUE;
                        angles.nang++;
                        angles.nopb++;

                        angles.i13[angles.nang][0] = atom[i].iat[1];
                        angles.i13[angles.nang][1] = i;
                        angles.i13[angles.nang][2] = atom[i].iat[2];                        
                        angles.i13[angles.nang][3] = atom[i].iat[0];                        
                        angles.copb[angles.nang] = ooplane_k.copb[j];
                        angles.angin[angles.nang] = TRUE;
                        angles.nang++;
                        angles.nopb++;
                        angles.i13[angles.nang][0] = atom[i].iat[1];
                        angles.i13[angles.nang][1] = i;
                        angles.i13[angles.nang][2] = atom[i].iat[3];                        
                        angles.i13[angles.nang][3] = atom[i].iat[0];                        
                        angles.copb[angles.nang] = ooplane_k.copb[j];
                        angles.angin[angles.nang] = TRUE;
                        angles.nang++;
                        angles.nopb++;
                        
                        angles.i13[angles.nang][0] = atom[i].iat[2];
                        angles.i13[angles.nang][1] = i;
                        angles.i13[angles.nang][2] = atom[i].iat[3];                        
                        angles.i13[angles.nang][3] = atom[i].iat[0];                        
                        angles.copb[angles.nang] = ooplane_k.copb[j];
                        angles.angin[angles.nang] = TRUE;
                        angles.nang++;
                        angles.nopb++;
                        break;
                     }
                  }
               }
           }
       }
   }
//
   if (angles.nopb == 0)
     pot.use_opbend = FALSE; 
}
