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
#include "angang.h"
#include "crossterm.h"
#include "eangang.h"                

void kangang()
{
    int i, j, k, m, l, jj, ierr;
    int it, ia, ic, itc;
    int nang;
    float fa, faa;
    angang.nangang = 0;
    ierr = FALSE;
    
    for (i=1; i <= natom; i++)
    {
        jj = 0;
        for (j=0; j < MAXIAT; j++)
           if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
             jj++;
        nang = (jj * (jj-1))/2;
        if (nang > 0)
        {
            for (j=0; j < angles.nang -1; j++)
            {
                if (angles.i13[j][1] == i && atom[angles.i13[j][0]].type < 200 && atom[angles.i13[j][2]].type < 200)
                {
                    ia = angles.i13[j][0];
                    ic = angles.i13[j][2];
                    it = atom[angles.i13[j][1]].type;
                    itc = atom[angles.i13[j][1]].tclass;
                    m = 0;
                    if (atom[ia].atomnum <= 1) m++;
                    if (atom[ic].atomnum <= 1) m++;
                    fa = 0.0;
                    ierr = FALSE;
                    for (k = 0; k < crossterm_k.nangang; k++)
                    {
                        if (crossterm_k.ang_ang[k] == it || crossterm_k.ang_ang[k] == itc)
                        {
                            fa = crossterm_k.aacon[k][m];
                            ierr = TRUE;
                            break;
                        }
                    }
                    // now find second angle
                    for (k = j+1; k < angles.nang; k++)
                    {
                        if(angles.i13[k][1] == i)
                        {
                            ia = angles.i13[k][0];
                            ic = angles.i13[k][2];
                            it = atom[angles.i13[k][1]].type;
                            m = 0;
                            if (atom[ia].atomnum <= 1) m++;
                            if (atom[ic].atomnum <= 1) m++;
                            faa = 0.0;
                            ierr = FALSE;
                            if (atom[ia].type < 200 && atom[ic].type < 200)
                            {
                               for (l = 0; l < crossterm_k.nangang; l++)
                               {
                                  if (crossterm_k.ang_ang[l] == it || crossterm_k.ang_ang[l] == itc)
                                  {
                                     faa = crossterm_k.aacon[l][m];
                                     ierr = TRUE;
                                     break;
                                  }
                                }
                             }
                             faa *= fa;
                             if (faa != 0.0)
                             {
                                 angang.iaa[angang.nangang][0] = j;
                                 angang.iaa[angang.nangang][1] = k;
                                 angang.fa[angang.nangang] = faa;
                                 angang.nangang++;
                             }
                        }
                    }
                }
            }
        }
    }
}
