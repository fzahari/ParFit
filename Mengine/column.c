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

#include "energies.h"
#include "derivs.h"
#include "column.h"


void column(int nvar, int *h_init, int *h_stop,int *h_index,
             int *c_init,int *c_stop,int *c_index, int *c_value)
{
    int i,j,k,m;

    for (i=0; i < nvar; i++)
    {
        c_init[i] = 0;
        c_stop[i] = 0;
    }
// count number of elements in each column
    for (i=0; i < nvar; i++)
    {
        for (j=h_init[i]; j < h_stop[i]; j++)
        {
            k = h_index[j];
            c_stop[k]++;
        }
    }
// set each beginning marker past the last element for its column
    c_init[0] = c_stop[0] + 1;
    for (i=1; i <nvar; i++)
    {
        c_init[i] = c_init[i-1] + c_stop[i];
    }
// set column index scanning rows in reverse order
    for (i = nvar-1; i >= 0; i--)
    {
        for (j=h_init[i]; j < h_stop[i]; j++)
        {
            k = h_index[j];
            m = c_init[k]-1;
            c_init[k] = m;
            c_index[m] = i;
            c_value[m] = j;
        }
    }
    for(i=0; i < nvar; i++)
    {
       c_stop[i] = c_init[i] + c_stop[i]-1;
    }
}
            
        
