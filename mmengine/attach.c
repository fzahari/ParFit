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
#include "attached.h"
#include "minim_control.h"

// ===================    
void attach()
{
    int i,j,k,m,p;
    int jj,kk,mm, jji;

    for (i=1; i <= natom; i++)
    {
        attached.n13[i] = 0;
        for (j=0; j < MAXIAT; j++)
        {
            jj = atom[i].iat[j];
            jji = 0;
            for (mm = 0; mm < MAXIAT; mm++)
            {
                if (atom[jj].iat[mm] != 0 && atom[jj].bo[mm] != 9)
                   jji++;
            }
            
            if ( atom[jj].type < 300 && jji < 5 )
            {
               if (jj != 0 && atom[i].bo[j] != 9)
               {
                  for(k=0; k < MAXIAT; k++)
                  {
                      kk = atom[jj].iat[k];
                      if (kk != i && kk != 0 && atom[jj].bo[k] != 9)
                      {
                         for (m=0; m < MAXIAT; m++)
                         {
                           if (kk == atom[i].iat[m])
                            break;
                         }
                         attached.i13[attached.n13[i]][i] = kk;
                         attached.n13[i]++;
                      }
                  }
               }
            }
        }
    }
  // find 14 relations
    for (i=1; i <= natom; i++)
    {
        attached.n14[i] = 0;
        for (j=0; j < MAXIAT; j++)
        {
            jj = atom[i].iat[j];
            if (jj != 0 && atom[i].bo[j] != 9 && atom[jj].type < 300)
            {
               for(k=0; k < MAXIAT; k++)
               {
                   kk = atom[jj].iat[k];
                   if (kk != 0 && atom[jj].bo[k] != 9 && atom[kk].type < 300)
                   {
                      for (m = 0; m < MAXIAT; m++)
                      {
                        mm = atom[kk].iat[m];
                        if (mm != i && mm != 0 && atom[kk].bo[m] != 9)
                        {
                          for (p=0; p < MAXIAT; p++)
                          {
                            if (mm == atom[i].iat[p])
                              goto L_20;
                          }
                          for (p=0; p < attached.n13[i]; p++)
                          {
                             if (mm == attached.i13[p][i])
                               goto L_20;
                          }
                          attached.i14[attached.n14[i]][i] = mm;
                          attached.n14[i]++;
                        }
L_20:
                       continue;
                     }
                   }
               }
            }
        }
    }
}


