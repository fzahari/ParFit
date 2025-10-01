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
#include "bonds_ff.h"
#include "field.h"
#include "crossterm.h"
#include "ebond.h"

struct t_strbnd {
        int nstrbnd, isb[MAXANG][3];
        float ksb1[MAXANG], ksb2[MAXANG];
        } strbnd;
struct stbndat {
        int irow, krow, lrow;
        float con1, con2;
      } stbn[] =
      {
        0,    1,    0,      0.15,      0.15,
        0,    1,    1,      0.10,      0.30,
        0,    1,    2,      0.05,      0.35,
        0,    1,    3,      0.05,      0.35,
        0,    1,    4,      0.05,      0.35,
        0,    2,    0,      0.00,      0.00,
        0,    2,    1,      0.00,      0.15,
        0,    2,    2,      0.00,      0.15,
        0,    2,    3,      0.00,      0.15,
        0,    2,    4,      0.00,      0.15,
        1,    1,    1,      0.30,      0.30,
        1,    1,    2,      0.30,      0.50,
        1,    1,    3,      0.30,      0.50,
        1,    1,    4,      0.30,      0.50,
        2,    1,    2,      0.50,      0.50,
        2,    1,    3,      0.50,      0.50,
        2,    1,    4,      0.50,      0.50,
        3,    1,    3,      0.50,      0.50,
        3,    1,    4,      0.50,      0.50,
        4,    1,    4,      0.50,      0.50,
        1,    2,    1,      0.30,      0.30,
        1,    2,    2,      0.25,      0.25,
        1,    2,    3,      0.25,      0.25,
        1,    2,    4,      0.25,      0.25,
        2,    2,    2,      0.25,      0.25,
        2,    2,    3,      0.25,      0.25,
        2,    2,    4,      0.25,      0.25,
        3,    2,    3,      0.25,      0.25,
        3,    2,    4,      0.25,      0.25,
        4,    2,    4,      0.25,      0.25  };

void kstrbnd()
{
   int i, j, iz, nb1, nb2, irow, krow, lrow;
   int ia, ib, ic, itb, ita, itc;
   int k, nstbn, index;
   char pa[4],pb[4],pc[4],pt[10],pt1[10];

   irow = krow = lrow = 0;
   strbnd.nstrbnd = 0;
   nb1 = -1;
   nb2 = -1;

   for (i=0; i < angles.nang; i++)
   {
       ia = angles.i13[i][0];
       ib = angles.i13[i][1];
       ic = angles.i13[i][2];
       ita = atom[ia].type;
       itb = atom[ib].type;
       itc = atom[ic].type;
       if (field.type == MMX && (atom[ia].type >= 300 || atom[ib].type >= 300 || atom[ic].type >= 300))
          goto L_10;
       if (field.type == MMX && (atom[ia].type >= 20 || atom[ic].type >= 20))
          goto L_10;
       if (field.type == MM3 && (atom[ia].type >= 300 || atom[ib].type >= 300 || atom[ic].type >= 300) )
          goto L_10;  
       numeral(ita,pa,3);
       numeral(itb,pb,3);
       numeral(itc,pc,3);
       if (ita <= itc)
       {
           strcpy(pt,pa);     strcat(pt,pb);  strcat(pt,pc);
       }else
       {
           strcpy(pt,pc);     strcat(pt,pb);  strcat(pt,pa);
       }
       itb = atom[ib].tclass;
       numeral(itb,pb,3);
       strcpy(pt1,"  0"); strcat(pt1,pb); strcat(pt1,"  0");
       
	   nstbn = 0;
       for (j=0; j < crossterm_k.nstrbnd; j++)
       {
           if (strcmp(crossterm_k.stbn[j],pt) == 0 || strcmp(crossterm_k.stbn[j],pt1) == 0 )
           {
               nstbn = j;
               break;
           }
       }
       for (j=0; j < bonds_ff.nbnd; j++)
       {
           iz = find_bond(ia, ib);
           if (iz != -1)
           {
              nb1 = iz;
              break;
           }
       }
      for (j=0; j < bonds_ff.nbnd; j++)
       {
           iz = find_bond(ib, ic);
           if (iz != -1)
           {
              nb2 = iz;
              break;
           }
       }
	   index = 0;
       if (angles.index[i] == 0)
          index = 0;
       else if (angles.index[i] == 1)
       {
          if (bonds_ff.index[nb1] == 1)
             index = 1;
          else
             index = 2;
       } else if (angles.index[i] == 2)
          index = 3;
       else if (angles.index[i] == 3)
          index = 5;
       else if (angles.index[i] == 4)
          index = 4;
       else if (angles.index[i] == 5)
       {
          if (bonds_ff.index[nb1] == 1)
             index = 6;
          else
             index = 7;
       } else if (angles.index[i] == 6)
          index = 8;
       else if (angles.index[i] == 7)
       {
          if (bonds_ff.index[nb1] == 1)
             index = 9;
          else
             index = 10;
       } else if (angles.index[i] == 8)
          index = 11;
       
       if (crossterm_k.stbnindex[nstbn] != index)
          nstbn = -1;
          
       if (field.type == MMFF94)
       {
		   irow = 0;
		   krow = 0;
           if (angles.anat[i] < 175.0)
           {
              if (nstbn == -1)  // constant not found  assign constants base on periodic table
              {
                 if (atom[ia].atomnum <= 1)
                  irow = 0;
                 else if (atom[ia].atomnum <= 10)
                  irow = 1;
                 else if (atom[ia].atomnum <= 18)
                  irow = 2;
                 else if (atom[ia].atomnum <= 36)
                  irow = 3;
                 else if (atom[ia].atomnum <= 54)
                  irow = 4;
                  
                 if (atom[ib].atomnum <= 1)
                  krow = 0;
                 else if (atom[ib].atomnum <= 10)
                  krow = 1;
                 else if (atom[ib].atomnum <= 18)
                  krow = 2;
                 else if (atom[ib].atomnum <= 36)
                  krow = 3;
                 else if (atom[ib].atomnum <= 54)
                  krow = 4;
                  
                 if (atom[ic].atomnum <= 1)
                  lrow = 0;
                 else if (atom[ic].atomnum <= 10)
                  lrow = 1;
                 else if (atom[ic].atomnum <= 18)
                  lrow = 2;
                 else if (atom[ic].atomnum <= 36)
                  lrow = 3;
                 else if (atom[ic].atomnum <= 54)
                  lrow = 4;

                for (j=0; j < 30; j++)
                {
                  if ( irow == stbn[j].irow && krow == stbn[j].krow && lrow == stbn[j].lrow )
                  {
                      strbnd.isb[strbnd.nstrbnd][0] = i;
                      strbnd.isb[strbnd.nstrbnd][1] = nb1;
                      strbnd.isb[strbnd.nstrbnd][2] = nb2;
                      strbnd.ksb1[strbnd.nstrbnd] = stbn[j].con1;
                      strbnd.ksb2[strbnd.nstrbnd] = stbn[j].con2;
                      strbnd.nstrbnd++;
                      break;
                  } else if (lrow == stbn[j].irow && krow == stbn[j].krow && irow == stbn[j].lrow)
                  {
                      strbnd.isb[strbnd.nstrbnd][0] = i;
                      strbnd.isb[strbnd.nstrbnd][1] = nb1;
                      strbnd.isb[strbnd.nstrbnd][2] = nb2;
                      strbnd.ksb1[strbnd.nstrbnd] = stbn[j].con2;
                      strbnd.ksb2[strbnd.nstrbnd] = stbn[j].con1;
                      strbnd.nstrbnd++;
                      break;
                  }
                }   
             }else
             {
                strbnd.isb[strbnd.nstrbnd][0] = i;
                strbnd.isb[strbnd.nstrbnd][1] = nb1;
                strbnd.isb[strbnd.nstrbnd][2] = nb2;
                if (angles.anat[i] > 175.0)
                {
                   strbnd.ksb1[strbnd.nstrbnd] = 0.0;
                   strbnd.ksb2[strbnd.nstrbnd] = 0.0;
                } else
                { 
                 if (ita < itc)
                 {
                    strbnd.ksb1[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][0];
                    strbnd.ksb2[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][1];
                 } else
                 {
                    strbnd.ksb1[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][1];
                    strbnd.ksb2[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][0];
                 }
                }
                strbnd.nstrbnd++;
             }
           }
       } else
       {         
             k = 0;
             if ( atom[ia].atomnum <= 1)
             {
                 k++;
                 if (fabs(crossterm_k.stbncon[nstbn][2]) < 0.1)
                   nb1 = -1;
             }
             if (atom[ic].atomnum <= 1)
             {
                 k++;
                 if (fabs(crossterm_k.stbncon[nstbn][2]) < 0.1)
                    nb2 = -1;
             }
             if (crossterm_k.stbncon[nstbn][k] != 0.0)
             {
                 strbnd.isb[strbnd.nstrbnd][0] = i;
                 strbnd.isb[strbnd.nstrbnd][1] = nb1;
                 strbnd.isb[strbnd.nstrbnd][2] = nb2;
                 strbnd.ksb1[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][k];           
                 strbnd.ksb2[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][k];           
                 strbnd.nstrbnd++;
             }
       }
L_10:
    continue;
   }
}
           
       
