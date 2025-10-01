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

#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"

#include "angles.h"
#include "field.h"
#include "pot.h"

// number of atoms with 5 or greater real bonds
void message_alert(char *, char *);
void get_angles(void);
int isangle(int, int);
int isbond(int,int);

void get_angles()
{
    int i,j,k,icount;
    int ia, ib, ic, id, ie;
    int jj, icoord;
    int itemp[5];

    ia = ib = ic = id = ie = 0;

    angles.nang = 0;
    high_coord.ncoord = 0;
    for (i=0; i < MAXANG; i++)
    {
        angles.i13[i][0] = 0;
        angles.i13[i][1] = 0;
        angles.i13[i][2] = 0;
        angles.i13[i][3] = 0;
    }
    for (i=0; i < MAXMETALANG; i++)
    {
        high_coord.i13[i][0] = 0;   
        high_coord.i13[i][1] = 0;   
        high_coord.i13[i][2] = 0;   
    }
//
    for (i = 1; i <= natom; i++)
    {
        jj = 0;
        icoord = 0;
        for (j=0; j < MAXIAT; j++)
        {
            if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
              jj++;
            if (atom[i].bo[j] == 9)
              icoord++;
        }

        if (jj == 2)
        {
            icount = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                {
                    icount++;
                    if (icount == 1)
                       ia = atom[i].iat[j];
                    else if (icount == 2)
                    {
                        ib = atom[i].iat[j];
                        break;
                    }
                }
            }
            if (atom[i].type < 300)
            {
               angles.i13[angles.nang][0] = ia;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = ib;
               angles.i13[angles.nang][3] = 0;
               angles.nang++;
               if (angles.nang > MAXANG)
               {
                    message_alert("Too many angles","Error in Get Angles");
                    exit(0);
                }
            } else
            {
                high_coord.i13[high_coord.ncoord][0] = ia;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = ib;
                high_coord.ncoord++;
               if (high_coord.ncoord > MAXMETALANG)
               {
                    message_alert("Too many metal angles","Error in Get Angles");
                    exit(0);
                }
            }         
        } else if (jj == 3)
        {
            if (icoord == 0)
            {
             ia = atom[i].iat[0];
             ib = atom[i].iat[1];
             ic = atom[i].iat[2];
            } else
            {
              icoord = 0;
              for(j=0; j < MAXIAT; j++)
              {
                 if (atom[i].iat[j] !=0 && atom[i].bo[j] != 9)
                 {
                     itemp[icoord] = atom[i].iat[j];
                     icoord++;
                 }
              }
              ia = itemp[0];
              ib = itemp[1];
              ic = itemp[2];
            }
            if (atom[i].type < 300)
            {
               angles.i13[angles.nang][0] = ia;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = ib;
               angles.i13[angles.nang][3] = ic;
               angles.nang++;            
               angles.i13[angles.nang][0] = ia;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = ic;
               angles.i13[angles.nang][3] = ib;
               angles.nang++;            
               angles.i13[angles.nang][0] = ib;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = ic;
               angles.i13[angles.nang][3] = ia;
               angles.nang++;
               if (angles.nang > MAXANG)
               {
                    message_alert("Too many angles","Error in Get Angles");
                    exit(0);
                }
            } else
            {
                high_coord.i13[high_coord.ncoord][0] = ia;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = ib;
                high_coord.ncoord++;
                high_coord.i13[high_coord.ncoord][0] = ia;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = ic;
                high_coord.ncoord++;
                high_coord.i13[high_coord.ncoord][0] = ib;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = ic;
                high_coord.ncoord++;
               if (high_coord.ncoord > MAXMETALANG)
               {
                    message_alert("Too many metal angles","Error in Get Angles");
                    exit(0);
                }
            }          
        } else if (jj == 4) 
        {
            if (icoord == 0)
            {
              ia = atom[i].iat[0];
              ib = atom[i].iat[1];
              ic = atom[i].iat[2];
              id = atom[i].iat[3];
            } else
            {
                icoord = 0;
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                    {
                        itemp[icoord] = atom[i].iat[j];
                        icoord++;
                    }
                }
                ia = itemp[0];
                ib = itemp[1];
                ic = itemp[2];
                id = itemp[3];
            }
            if (atom[i].type < 300)
            {
               angles.i13[angles.nang][0] = ia;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = ib;
               angles.i13[angles.nang][3] = 0;
               angles.nang++;            
               angles.i13[angles.nang][0] = ia;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = ic;
               angles.i13[angles.nang][3] = 0;
               angles.nang++;            
               angles.i13[angles.nang][0] = ia;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = id;
               angles.i13[angles.nang][3] = 0;
               angles.nang++;            
               angles.i13[angles.nang][0] = ib;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = ic;
               angles.i13[angles.nang][3] = 0;
               angles.nang++;            
               angles.i13[angles.nang][0] = ib;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = id;
               angles.i13[angles.nang][3] = 0;
               angles.nang++;            
               angles.i13[angles.nang][0] = ic;
               angles.i13[angles.nang][1] = i;
               angles.i13[angles.nang][2] = id;
               angles.i13[angles.nang][3] = 0;
               angles.nang++;
               if (angles.nang > MAXANG)
               {
                    message_alert("Too many angles","Error in Get Angles");
                    exit(0);
                }
            } else
            {
                high_coord.i13[high_coord.ncoord][0] = ia;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = ib;
                high_coord.ncoord++;
                high_coord.i13[high_coord.ncoord][0] = ia;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = ic;
                high_coord.ncoord++;
                high_coord.i13[high_coord.ncoord][0] = ia;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = id;
                high_coord.ncoord++;
                high_coord.i13[high_coord.ncoord][0] = ib;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = ic;
                high_coord.ncoord++;
                high_coord.i13[high_coord.ncoord][0] = ib;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = id;
                high_coord.ncoord++;
                high_coord.i13[high_coord.ncoord][0] = ic;          
                high_coord.i13[high_coord.ncoord][1] = i;          
                high_coord.i13[high_coord.ncoord][2] = id;
                high_coord.ncoord++;
               if (high_coord.ncoord > MAXMETALANG)
               {
                    message_alert("Too many metal angles","Error in Get Angles");
                    exit(0);
                }
            }                     
        } else if (jj == 5)  // must be POS metal ??
        {
            if (icoord == 0)
            {
              ia = atom[i].iat[0];
              ib = atom[i].iat[1];
              ic = atom[i].iat[2];
              id = atom[i].iat[3];
              ie = atom[i].iat[4];
            } else
            {
                icoord = 0;
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                    {
                        itemp[icoord] = atom[i].iat[j];
                        icoord++;
                    }
                }
                ia = itemp[0];
                ib = itemp[1];
                ic = itemp[2];
                id = itemp[3];
                ie = itemp[4];
            }
              // ia angles
              high_coord.i13[high_coord.ncoord][0] = ia;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = ib;
              high_coord.ncoord++;
              high_coord.i13[high_coord.ncoord][0] = ia;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = ic;
              high_coord.ncoord++;
              high_coord.i13[high_coord.ncoord][0] = ia;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = id;
              high_coord.ncoord++;
              high_coord.i13[high_coord.ncoord][0] = ia;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = ie;
              high_coord.ncoord++;
              // ib angles
              high_coord.i13[high_coord.ncoord][0] = ib;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = ic;
              high_coord.ncoord++;
              high_coord.i13[high_coord.ncoord][0] = ib;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = id;
              high_coord.ncoord++;
              high_coord.i13[high_coord.ncoord][0] = ib;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = ie;
              high_coord.ncoord++;
              // ic angles
              high_coord.i13[high_coord.ncoord][0] = ic;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = id;
              high_coord.ncoord++;
              high_coord.i13[high_coord.ncoord][0] = ic;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = ie;
              high_coord.ncoord++;
              // id angles
              high_coord.i13[high_coord.ncoord][0] = id;
              high_coord.i13[high_coord.ncoord][1] = i;
              high_coord.i13[high_coord.ncoord][2] = ie;
              high_coord.ncoord++;
               if (high_coord.ncoord > MAXMETALANG)
               {
                    message_alert("Too many metal angles","Error in Get Angles");
                    exit(0);
                }
        } else if (jj > 5)
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                {
                    ia = atom[i].iat[j];
                    for (k=j+1; k < MAXIAT; k++)
                    {
                        if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
                        {
                            high_coord.i13[high_coord.ncoord][0] = ia;
                            high_coord.i13[high_coord.ncoord][1] = i;
                            high_coord.i13[high_coord.ncoord][2] = atom[i].iat[k];
                            high_coord.ncoord++;
                            if (high_coord.ncoord > MAXMETALANG)
                            {
                                message_alert("Too many metal angles","Error in Get Angles");
                                exit(0);
                            }
                        }
                    }
                }
            }
        }

        if (angles.nang > MAXANG)
        {
            message_alert("Too many angles. PCMODEL will now exit","Angle Setup");
            fprintf(pcmoutfile,"Too many angles. PCMODEL will now exit\n");
            exit(0);
        }
    }
    if (high_coord.ncoord > 0)
       pot.use_highcoord = TRUE;
       
    for (i=0; i < angles.nang; i++)
       angles.angin[i] = FALSE;
}

int isangle(int i, int j)
{
    int k,l,katm;
    for (k=0; k < MAXIAT; k++)
    {
        if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
        {
            katm = atom[i].iat[k];
	    if (atom[katm].type < 300)  // test for POS metals
	    {
				for (l=0; l < MAXIAT; l++)
				{
					if (atom[katm].iat[l] == j)
						return TRUE;
				}
	    }
        }
    }
    return FALSE;
}

