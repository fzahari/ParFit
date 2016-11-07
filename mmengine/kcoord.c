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
#include "field.h"
#include "coordb.h"
#include "vdwk.h"

void kcoord_bonds()
{
     int i, ia, ib,ierr, ita,itb;
     int jj, nusi;
     long int mask5,mask6,mask7,mask8;
     float rad,eps;

     mask5 = 1L << SATMET_MASK;
     mask6 = 1L << GT18e_MASK;
     mask7 = 1L << LOWSPIN_MASK;
     mask8 = 1L << SQPLAN_MASK;

     for (i=0; i < cbonds.nbnd; i++)
     {
         ia = cbonds.icbonds[i][0];
         ib = cbonds.icbonds[i][1];
         ita = atom[ia].type;
         itb = atom[ib].type;
         nusi = 0;
         rad = vdw1.rad[ita];
         eps = vdw1.eps[ita];
         ierr = FALSE;
         if (strcmp(field.radiusrule,"ARITHMETIC") == 0)
         {
            cbonds.vrad[i] = rad + vdw1.rad[atom[ib].type];
         } else if (strcmp(field.radiusrule,"GEOMETRIC") == 0)
         {
            cbonds.vrad[i] = 2.0*( sqrt(rad*vdw1.rad[atom[ib].type]));
         } else
            cbonds.vrad[i] = rad + vdw1.rad[atom[ib].type]; 
      
         if (strcmp(field.epsrule,"ARITHMETIC") == 0)
         {
             cbonds.veps[i] = 0.5*(eps + vdw1.eps[atom[ib].type]);
         } else if (strcmp(field.epsrule,"GEOMETRIC") == 0)
         {
             cbonds.veps[i] = sqrt( eps*vdw1.eps[atom[ib].type] );
         }else
             cbonds.veps[i] = sqrt( eps*vdw1.eps[atom[ib].type] );

         if ( (atom[ia].flags & mask6) && !( atom[ia].flags & mask5 ))
               nusi = 2;
         else if( ( atom[ia].flags &mask5 ) && ( atom[ia].flags & mask6 ) )
         {
             if( ( atom[ia].flags & mask7 ) && !( atom[ia].flags & mask8 ) )
               nusi = 1;
             else if(( atom[ia].flags & mask8 ) && !( atom[ia].flags & mask7 ) )
               nusi = 3;
          }

     // metal to pi atom
         if (atom[ib].type == 2 || atom[ib].type == 4 || atom[ib].type == 29 ||
             atom[ib].type == 30 || atom[ib].type == 40 || atom[ib].type == 48 )
         {
             cbonds.ia1[i] = 0;
             cbonds.ia2[i] = 0;
             for (jj = 0 ; jj < MAXIAT; jj++)
             {
                 if (atom[ib].iat[jj] != 0 && atom[ib].iat[jj] != ia && cbonds.ia1[i] == 0)
                     cbonds.ia1[i] = atom[ib].iat[jj];
                 else if (atom[ib].iat[jj] != 0 && atom[ib].iat[jj] != ia && cbonds.ia2[i] == 0)
                     cbonds.ia2[i] = atom[ib].iat[jj];
             }
             cbonds.vrad[i] = rad + 0.66; 
             if (nusi > 0)
                 cbonds.vrad[i] += 0.15;
             cbonds.lpd[i] = vdw1.lpd[atom[ib].type];
             if (nusi > 0)
                cbonds.lpd[i] *= 0.85;
         } else
         {     // metal to nitrogen and oxygen lone pairs
             cbonds.vrad[i] = rad + vdw1.rad[atom[ib].type]; // use radius of lone pair 
             cbonds.vrad[i] *= 0.6; // cut radius in half
             
             cbonds.lpd[i] = vdw1.lpd[atom[ib].type];
             if (atom[ib].type == 7)
                cbonds.lpd[i] *= 3.0;
             else if (atom[ib].type == 10)
                cbonds.lpd[i] *= 4.0;

             if( ( atom[ia].flags & mask6 ) && !( atom[ia].flags & mask5 ) )
             {
                 cbonds.vrad[i] += .15;
                 cbonds.lpd[i] *= .85;
             }

             if( nusi == 1 || (nusi == 3 && (atom[atom[ib].iat[0]].type == 6 || atom[atom[ib].iat[0]].type == 15)) )
                 cbonds.vrad[i] -= .1;
         }
     }

}
