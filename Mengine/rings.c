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
#include "rings.h"
#include "utility.h"

// ===========================
void allocate_rings()
{
  rings.nring3 = 0;
  rings.nring4 = 0;
  rings.nring5 = 0;
  rings.nring6 = 0;
  rings.r13 = imatrix(0,(natom+10)/3, 0,3);
  rings.r14 = imatrix(0,(natom+10)/4, 0,4);
  rings.r15 = imatrix(0,(natom+10)/5, 0,5);
  rings.r16 = imatrix(0,(natom+10)/2, 0,6);
}
// =========================
void free_rings()
{
  free_imatrix(rings.r13 ,0,(natom+10)/3, 0,3);
  free_imatrix(rings.r14 ,0,(natom+10)/4, 0,4);
  free_imatrix(rings.r15 ,0,(natom+10)/5, 0,5);
  free_imatrix(rings.r16 ,0,(natom+10)/2, 0,6);
}
/* --------------------------------------------------------- */
// iatom is atom number to search on
// isize = ring size
// num = number of ring if more than one ring containing atom iatom
// array = array of ring atoms
//
void get_rsize(int iatom, int isize, int num, int *array)
{
    int i, j, k, icount;
    icount = -1;
    
    if (isize == 3)
    {
       for (i=0; i < rings.nring3; i++)
       {
           for (j=0; j < 3; j++)
           {
               if (iatom == rings.r13[i][j])
                  icount++;
               if (icount == num)
               {
                  for (k=0; k < isize; k++)
                      array[k] = rings.r13[i][k];
                  return;
               }
           }
        }
    } else if (isize == 4)
    {
       for (i=0; i < rings.nring4; i++)
       {
           for (j=0; j < 4; j++)
           {
               if (iatom == rings.r14[i][j])
                  icount++;
               if (icount == num)
               {
                  for (k=0; k < isize; k++)
                      array[k] = rings.r14[i][k];
                  return;
               }
           }
        }
    } else if (isize == 5)
    {
       for (i=0; i < rings.nring5; i++)
       {
           for (j=0; j < 5; j++)
           {
               if (iatom == rings.r15[i][j])
                  icount++;
               if (icount == num)
               {
                  for (k=0; k < isize; k++)
                      array[k] = rings.r15[i][k];
                  return;
               }
           }
        }
    } else if (isize == 6)
    {
       for (i=0; i < rings.nring6; i++)
       {
           for (j=0; j < 6; j++)
           {
               if (iatom == rings.r16[i][j])
                  icount++;
               if (icount == num)
               {
                  for (k=0; k < isize; k++)
                      array[k] = rings.r16[i][k];
                  return;
               }
           }
        }
    } 
}
/* --------------------------------------------------------- */
int find_rsize(int isize,int iatom)
{
    int i, j, icount;
    
    icount = 0;
    if (isize == 3)
    {
        for (i=0; i < rings.nring3; i++)
        {
            for(j=0; j < 3; j++)
            {
                if (rings.r13[i][j] == iatom)
                    icount++;
            }
        }
        return(icount);
    }else if (isize == 4)
    {
        for (i=0; i < rings.nring4; i++)
        {
            for(j=0; j < 4; j++)
            {
                if (rings.r14[i][j] == iatom)
                    icount++;
            }
        }
        return(icount);
    }else if (isize == 5)
    {
        for (i=0; i < rings.nring5; i++)
        {
            for(j=0; j < 5; j++)
            {
                if (rings.r15[i][j] == iatom)
                    icount++;
            }
        }
        return(icount);
    }else if (isize == 6)
    {
        for (i=0; i < rings.nring6; i++)
        {
            for(j=0; j < 6; j++)
            {
                if (rings.r16[i][j] == iatom)
                    icount++;
            }
        }
        return(icount);
    } 
    return(icount);
}
// ===============================
int is_ring31(int ia)
{
    int i,j;

    for (i=0; i < rings.nring3; i++)
    {
        for (j=0; j < 3; j++)
        {
            if (ia == rings.r13[i][j])
            {
               return(TRUE);
            }
        }
    }
    return FALSE;
}
/* ============================ */
int is_ring41(int ia)
{
    int i,j, found;

    found = FALSE;
    for (i=0; i < rings.nring4; i++)
    {
        for (j=0; j < 4; j++)
        {
            if (ia == rings.r14[i][j])
            {
                found = TRUE;
                return (found);
            }
        }
    }
    return (found);
}
/* -------------------------------------------------------- */            
int is_ring42(int ia, int ib)
{
    int i,j, k, itmp, jtmp;
    if (ia < ib )
    {
        itmp = ia;
        jtmp = ib;
    }else
    {
        itmp = ib;
        jtmp = ia;
    }
    for (i=0; i < rings.nring4; i++)
    {
        for (j = 0; j < 3; j++)
        {
            if ( itmp == rings.r14[i][j])
            {
                for (k=j+1; k < 4; k++)
                {
                    if (jtmp == rings.r14[i][k])
                       return(TRUE);
                }
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */            
int is_ring43(int ia, int ib, int ic)
{
    int i, found;
    int array[5];

    found = FALSE;
    array[0] = ia;
    array[1] = ib;
    array[2] = ic;
    ksort(3,array);
    for (i=0; i < rings.nring4; i++)
    {
        if ( array[0] == rings.r14[i][0] && rings.r14[i][1] == array[1] && rings.r14[i][2] == array[2])
        {
            found = TRUE;
            return(found);
        }
        if ( rings.r14[i][0] == array[0] && rings.r14[i][1] == array[1] && rings.r14[i][3] == array[2])
        {
            found = TRUE;
            return(found);
        }
        if ( rings.r14[i][0] == array[0] && rings.r14[i][2] == array[1] && rings.r14[i][3] == array[2])
        {
            found = TRUE;
            return(found);
        }
        if ( rings.r14[i][1] == array[0] && rings.r14[i][2] == array[1] && rings.r14[i][3] == array[2])
        {
            found = TRUE;
            return(found);
        }
    }
    return(found);
}
/* -------------------------------------------------------- */    
int is_ring44(int *array)
{
    int i;
    for (i=0; i < rings.nring4; i++)
    {
        if (array[0] == rings.r14[i][0] && array[1] == rings.r14[i][1] &&
            array[2] == rings.r14[i][2] && array[3] == rings.r14[i][3])
            return(TRUE);
    }
    return FALSE;
}
/* -------------------------------------------------------- */    
int is_ring51(int ia)
{
    int i,j;

    for(i=0; i < rings.nring5; i++)
    {
        for (j=0; j < 5; j++)
        {
            if (rings.r15[i][j] == ia)
            {
                return(TRUE);
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */
int is_ring52(int ia, int ib)
{
    int i, j,k, itmp, jtmp;
    if (ia < ib)
    {
        itmp = ia;
        jtmp = ib;
    }else
    {
        itmp = ib;
        jtmp = ia;
    }
    for (i=0; i < rings.nring5; i++)
    {
       for (j=0; j < 4; j++)
       {
           if (itmp == rings.r15[i][j])
           {
               for (k=j+1; k < 5; k++)
               {
                   if (jtmp == rings.r15[i][k])
                       return(TRUE);
               }
           }
       }
    }
    return(FALSE);
}     
/* -------------------------------------------------------- */    
int is_ring53(int ia, int ib, int ic)
{
    int i, j,k,l;
    int array[5];

    array[0] = ia;
    array[1] = ib;
    array[2] = ic;
    ksort(3,array);
    for (i=0; i < rings.nring5; i++)
    {
        for (j=0; j < 3; j++)
        {
            if (array[0] == rings.r15[i][j])
            {
                for (k=j+1; k < 4; k++)
                {
                    if (array[1] == rings.r15[i][k])
                    {
                        for (l=j+1; l < 5; l++)
                        {
                            if (array[2] == rings.r15[i][l])
                                return(TRUE);
                        }
                    }
                }
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */
int is_ring54(int ia, int ib, int ic, int id)
{
    int i;
    int array[5];
    array[0] = ia;
    array[1] = ib;
    array[2] = ic;
    array[3] = id;
    ksort(4,array);
    for (i=0; i < rings.nring5; i++)
    {
        if (array[0] == rings.r15[i][0])
        {
            if (array[1] == rings.r15[i][1])
            {
                if (array[2] == rings.r15[i][2])
                {
                    if (array[3] == rings.r15[i][3])
                       return(TRUE);
                    else if (array[3] == rings.r15[i][4])
                       return(TRUE);
                    else
                       return(FALSE);
                } else if (array[2] == rings.r15[i][3] && array[3] == rings.r15[i][4])
                     return(TRUE);
            } else if (array[1] == rings.r15[i][2] &&
                       array[2] == rings.r15[i][3] && array[3] == rings.r15[i][4])
                       return(TRUE);
            
        } else if (array[0] == rings.r15[i][1] && array[1] == rings.r15[i][2] &&
            array[2] == rings.r15[i][3] && array[3] == rings.r15[i][4])
            return(TRUE); 
    }
    return (FALSE);
}
/* -------------------------------------------------------- */    
int is_ring55(int *array)
{
    int i;
    for(i=0; i < rings.nring5; i++)
    {
        if ( array[0] == rings.r15[i][0] && array[1] == rings.r15[i][1] && array[2] == rings.r15[i][2] &&
             array[3] == rings.r15[i][3] && array[4] == rings.r15[i][4])
             return(TRUE);
    }
    return (FALSE);
}
/* -------------------------------------------------------- */    
int is_ring61(int ia)
{
    int i, j;
    for(i=0; i < rings.nring6; i++)
    {
        for(j=0; j < 6; j++)
        {
            if (ia == rings.r16[i][j])
               return (TRUE);
        }
    }
    return (FALSE);
}
/* -------------------------------------------------------- */
int is_ring62(int ia, int ib)
{
    int i, j,k, itmp, jtmp;
/*
    if (ia < ib)
    {
        itmp = ia;
        jtmp = ib;
    }else
    {
        itmp = ib;
        jtmp = ia;
    }
*/
    itmp = ia;
    jtmp = ib;    
        
    for (i=0; i < rings.nring6; i++)
    {
       for (j=0; j < 6; j++)
       {
           if (itmp == rings.r16[i][j])
           {
               for (k=0; k < 6; k++)
               {
                   if(k!=j){
                       if (jtmp == rings.r16[i][k])
                           return(TRUE);
                   }
               }
           }
       }
    }
    return(FALSE);
}     
/* -------------------------------------------------------- */    
int is_ring63(int ia, int ib, int ic)
{
    int i, j, k, l;
    int array[5];

    array[0] = ia;
    array[1] = ib;
    array[2] = ic;
/*
    ksort(3,array);
*/
    for (i=0; i < rings.nring6; i++)
    {
        for (j=0; j < 6; j++)
        {
            if (array[0] == rings.r16[i][j])
            {
                for (k=0; k < 6; k++)
                {
                    if(k!=j){
                        if (array[1] == rings.r16[i][k])
                        {
                            for (l=0; l < 6; l++)
                            {
                                if(l!=j && l!=k){
                                    if (array[2] == rings.r16[i][l])
                                        return(TRUE);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */    
int is_ring66(int *array)
{
    int i;
    for(i=0; i < rings.nring6; i++)
    {
        if ( array[0] == rings.r16[i][0] && array[1] == rings.r16[i][1] && array[2] == rings.r16[i][2] &&
             array[3] == rings.r16[i][3] && array[4] == rings.r16[i][4] && array[5] == rings.r16[i][5])
             return(TRUE);
    }
    return (FALSE);
}
/* -------------------------------------------------------- */    
void get_rings()
{
	int i,j, add_ring, k, l,m, n, jj,kk,isix;
	int jatm,katm,latm,matm;
	int array[6], array1[6],array2[6];
	long int aromatic_mask;
	
	aromatic_mask = (1L << AROMATIC_MASK);

   rings.nring3 = 0;
   rings.nring4 = 0;
   rings.nring5 = 0;
   rings.nring6 = 0;
   // allocate memory
   allocate_rings();
// get three membered rings
   for (i=1; i <= natom; i++)
   {
     for (j=0; j < MAXIAT; j++)
       {
	 if (atom[i].iat[j] != 0)
	   {
	     jatm = atom[i].iat[j];
	     for (k=0; k < MAXIAT; k++)
	       {
		 if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i)
		   {
		     if (isbond(i,atom[jatm].iat[k]))
		       {
			 array[0] = i;
			 array[1] = jatm;
			 array[2] = atom[jatm].iat[k];
			 add_ring = TRUE;
			 for (jj=0; jj < 3; jj++)
			   array1[jj] = array[jj];

			 for (jj = 0; jj < rings.nring3; jj++)
			   {
			     for (kk = 0; kk < 3; kk++)
			       array2[kk] = rings.r13[jj][kk];
			     if (icompare(3,array1,array2))
			       add_ring = FALSE;
			   }
			 if (add_ring == TRUE)
			   { 
				   ksort(3, array);
			     rings.r13[rings.nring3][0] = array[0];
			     rings.r13[rings.nring3][1] = array[1];
			     rings.r13[rings.nring3][2] = array[2];
			     rings.nring3++;
			   }
		       }
		   }
	       }
	   }
       }
   }
// get four membered rings
   for (i=1; i <= natom; i++)
   {
     for (j=0; j < MAXIAT; j++)
       {
	 if (atom[i].iat[j] != 0)
	   {
	     jatm = atom[i].iat[j];
	     for (k=0; k < MAXIAT; k++)
	       {
		 if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i)
		   {
		     katm = atom[jatm].iat[k];
		     for (l = 0; l < MAXIAT; l++)
		       {
			 if (atom[katm].iat[l] != 0 && atom[katm].iat[l] != i && atom[katm].iat[l] != jatm)
			   {
			     if (isbond(i,atom[katm].iat[l]))
			       {
				 array[0] = i;
				 array[1] = jatm;
				 array[2] = katm;
				 array[3] = atom[katm].iat[l];
				 add_ring = TRUE;
				 for (jj=0; jj < 4; jj++)
				   array1[jj] = array[jj];

				 for (jj = 0; jj < rings.nring4; jj++)
				   {
				     for (kk = 0; kk < 4; kk++)
				       array2[kk] = rings.r14[jj][kk];
				     if (icompare(4,array1,array2))
				       add_ring = FALSE;
				   }
				 if (add_ring == TRUE)
				   {  
					   ksort(4, array);
				     rings.r14[rings.nring4][0] = array[0];
				     rings.r14[rings.nring4][1] = array[1];
				     rings.r14[rings.nring4][2] = array[2];
				     rings.r14[rings.nring4][3] = array[3];
				     rings.nring4++;
				   }
			       }
			   }
		       }
		   }
	       }
	   }
       }
   }
// get five membered rings
   for (i=1; i <= natom; i++)
   {
     for (j=0; j < MAXIAT; j++)
       {
	 if (atom[i].iat[j] != 0)
	   {
	     jatm = atom[i].iat[j];
	     for (k=0; k < MAXIAT; k++)
	       {
		 if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i)
		   {
		     katm = atom[jatm].iat[k];
		     for (l = 0; l < MAXIAT; l++)
		       {
			 if (atom[katm].iat[l] != 0 && atom[katm].iat[l] != i && atom[katm].iat[l] != jatm)
			   {
			     latm = atom[katm].iat[l];
			     for (m=0; m < MAXIAT; m++)
			       {
				 if (atom[latm].iat[m] != 0 && atom[latm].iat[m] != i && atom[latm].iat[m] != jatm && atom[latm].iat[m] != katm)
				   {
				     if (isbond(i,atom[latm].iat[m]))
				       {
					 array[0] = i;
					 array[1] = jatm;
					 array[2] = katm;
					 array[3] = latm;
					 array[4] = atom[latm].iat[m];
					 add_ring = TRUE;
					 for (jj=0; jj < 5; jj++)
					   array1[jj] = array[jj];

					 for (jj = 0; jj < rings.nring5; jj++)
					   {
					     for (kk = 0; kk < 5; kk++)
					       array2[kk] = rings.r15[jj][kk];
					     if (icompare(5,array1,array2) == TRUE)
					       add_ring = FALSE;
					   }
					 if (add_ring == TRUE)
					   {   
						   ksort(5, array);
					     rings.r15[rings.nring5][0] = array[0];
					     rings.r15[rings.nring5][1] = array[1];
					     rings.r15[rings.nring5][2] = array[2];
					     rings.r15[rings.nring5][3] = array[3];
					     rings.r15[rings.nring5][4] = array[4];
					     rings.nring5++;
					   }
				       }
				   }
			       }
			   }
		       }
		   }
	       }
	   }
       }
   }
// get six membered rings
   for (i=1; i <= natom; i++)
   {
     for (j=0; j < MAXIAT; j++)
       {
	 if (atom[i].iat[j] != 0)
	   {
	     jatm = atom[i].iat[j];
	     for (k=0; k < MAXIAT; k++)
	       {
		 if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i)
		   {
		     katm = atom[jatm].iat[k];
		     for (l = 0; l < MAXIAT; l++)
		       {
			 if (atom[katm].iat[l] != 0 && atom[katm].iat[l] != i && atom[katm].iat[l] != jatm)
			   {
			     latm = atom[katm].iat[l];
			     for (m=0; m < MAXIAT; m++)
			       {
				 if (atom[latm].iat[m] != 0 && atom[latm].iat[m] != i && atom[latm].iat[m] != jatm && atom[latm].iat[m] != katm)
				   {
				     matm = atom[latm].iat[m];
				     for (n=0; n < MAXIAT; n++)
				       {
					 if (atom[matm].iat[n] != 0 && atom[matm].iat[n] != i && atom[matm].iat[n] != jatm && atom[matm].iat[n] != katm && atom[matm].iat[n] != latm)
					   {
					     if (isbond(i,atom[matm].iat[n]))
					       {
						 array[0] = i;
						 array[1] = jatm;
						 array[2] = katm;
						 array[3] = latm;
						 array[4] = matm;
						 array[5] = atom[matm].iat[n];
						 add_ring = TRUE;
						 for (jj=0; jj < 6; jj++)
						   array1[jj] = array[jj];

						 for (jj = 0; jj < rings.nring6; jj++)
						   {
						     for (kk = 0; kk < 6; kk++)
						       array2[kk] = rings.r16[jj][kk];
						     if (icompare(6,array1,array2))
						       add_ring = FALSE;
						   }
						 if (add_ring == TRUE)
						   {  
							   ksort(6, array);
						     rings.r16[rings.nring6][0] = array[0];
						     rings.r16[rings.nring6][1] = array[1];
						     rings.r16[rings.nring6][2] = array[2];
						     rings.r16[rings.nring6][3] = array[3];
						     rings.r16[rings.nring6][4] = array[4];
						     rings.r16[rings.nring6][5] = array[5];
						     rings.nring6++;
						   }
					       }
					   }
				       }
				   }
			       }
			   }
		       }
		   }
	       }
	   }
       }
   }
	// find aromatic rings
	for (i=0; i < rings.nring6; i++)
	{
		for(j=0; j < 6; j++)
			array[j] = rings.r16[i][j];
		isix = aromatic_6(array);
		if (isix == TRUE)
		{
			for(j=0; j < 6; j++)
				atom[rings.r16[i][j]].flags |= aromatic_mask;
		}
	}
    for (i=0; i < rings.nring5; i++)
    {
		for(j=0; j < 5; j++)
			array[j] = rings.r15[i][j];
		isix = aromatic_5(i,array);
		if (isix == TRUE)
		{
			for(j=0; j < 5; j++)
				atom[rings.r15[i][j]].flags |= aromatic_mask;
		}
    }
}
/* =================================================== */
void ksort(int num, int array[])
{
    int i,temp;
    int found;

L_1:
    found = FALSE;
    for (i=0; i < num-1; i++)
    {
       if (array[i+1] < array[i])
       {
          temp = array[i];
          array[i] = array[i+1];
          array[i+1] = temp;
          found = TRUE;
       }
    }
    if (found == TRUE)
       goto L_1;
}
/* =======================================  */
int icompare(int iatom, int *iarray1, int *iarray2)
{
    int i;
    int scratch1[50],scratch2[50];

    for (i=0; i < iatom; i++)
    {
        scratch1[i] = iarray1[i];
        scratch2[i] = iarray2[i];
    }
    ksort(iatom, scratch1);
    ksort(iatom, scratch2);
    for (i=0; i < iatom; i++)
    {
        if (scratch1[i] != scratch2[i])
            return FALSE;
    }
    return TRUE;
}
/* --------------------------------------------------------- */
int aromatic_5(int ij, int *array)
{
  int i,j,npi,ihetero,found;
  int jatm,katm;
  int ia,ib,ic,idbl,jdbl;
  int db[5];
    
    npi = 0;
    ihetero = 0;
    for (i=0; i < 5; i++)
      db[i] = 0;
    // printf("array: %d %d %d %d %d\n",array[0],array[1],array[2],array[3],array[4]);
    for (i=0; i < 5; i++)
      {
	if (has_double_bond(array[i]) == TRUE)
	  {
	    if (check_connectivity(array[i]) == FALSE)
	      {
		return FALSE;
	      }
	  } else
	  {
	    if (check_connectivity_yatom(array[i]) == FALSE)
	      {
		return FALSE;
	      }
	  }
      }
    // printf("all aromatic\n");
    //
    for (i=0; i < 5; i++)
      {
	jatm = array[i];
	if (atom[jatm].atomnum > 6)
	  {
	    found = FALSE;
	    for (j=0; j < MAXIAT;  j++)
	      {
		if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] == 2)
		  found = TRUE;
	      }
	    if (found == FALSE) ihetero = jatm;
	  }
      }
    for (i=0; i < 4; i++)
      {
	jatm = array[i];
	for (j=i+1; j < 5; j++)
	  {
	    katm = array[j];
	    if (isbond(jatm,katm) && get_bondorder(jatm,katm) == 2)
	      {
		npi++;
		if (i==0 && j == 1)
		  db[0] = 1;
		else if (i == 0 && j == 4)
		  db[4] = 1;
		else 
		  db[i] = 1;
	      }
	  }
      }
    //printf("aro 5: %d %d - %d %d %d %d %d\n",npi,ihetero,db[0],db[1],db[2],db[3],db[4]);
    ia = ib = ic = 0;
    if (npi == 2 && ihetero != 0) return TRUE;  // two double bonds and one hetero atom - aromatic
    if (npi == 0) return FALSE; // no double bonds in ring not aromatic
    if (npi == 1) // 
      {
        if (db[0] == 1)
        {
            ia = array[2];
            ib = array[3];
            ic = array[4];
        } else if (db[1] == 1)
        {
            ia = array[0];
            ib = array[3];
            ic = array[4];
        } else if (db[2] == 1)
        {
            ia = array[0];
            ib = array[1];
            ic = array[4];
        } else if (db[3] == 1)
        {
            ia = array[0];
            ib = array[1];
            ic = array[2];
        } else if (db[4] == 1)
        {
            ia = array[1];
            ib = array[2];
            ic = array[3];
        }
	if (ia == ihetero)
	  {
	    idbl = FALSE;
	    jdbl = FALSE;
	    for (i=0; i < MAXIAT; i++)
	      {
		if (atom[ib].iat[i] != 0 && atom[ib].bo[i] == 2 && is_ring61(atom[ib].iat[i]))
		  idbl = TRUE;
	      }
	    for (i=0; i < MAXIAT; i++)
	      {
		if (atom[ic].iat[i] != 0 && atom[ic].bo[i] == 2 && is_ring61(atom[ic].iat[i]))
		  jdbl = TRUE;
	      }
	    if (idbl && jdbl) return TRUE;
	  } else if (ib == ihetero)
	  {
	    idbl = FALSE;
	    jdbl = FALSE;
	    for (i=0; i < MAXIAT; i++)
	      {
		if (atom[ia].iat[i] != 0 && atom[ia].bo[i] == 2 && is_ring61(atom[ia].iat[i]))
		  idbl = TRUE;
	      }
	    for (i=0; i < MAXIAT; i++)
	      {
		if (atom[ic].iat[i] != 0 && atom[ic].bo[i] == 2 && is_ring61(atom[ic].iat[i]))
		  jdbl = TRUE;
	      }
	    if (idbl && jdbl) return TRUE;
	  } else if (ic == ihetero)
	  {
	    idbl = FALSE;
	    jdbl = FALSE;
	    for (i=0; i < MAXIAT; i++)
	      {
		if (atom[ib].iat[i] != 0 && atom[ib].bo[i] == 2 && is_ring61(atom[ib].iat[i]))
		  idbl = TRUE;
	      }
	    for (i=0; i < MAXIAT; i++)
	      {
		if (atom[ia].iat[i] != 0 && atom[ia].bo[i] == 2 && is_ring61(atom[ia].iat[i]))
		  jdbl = TRUE;
	      }
	    if (idbl && jdbl) return TRUE;
	  }
      }
    return FALSE;                  
}
/* --------------------------------------------------------- */
int aromatic_6(int *array)
{
    int i,j,k, ia,ib, idbl,kdbl;
    int num, inarray[10];
    long int aromatic_mask, mask6;
	
    aromatic_mask = (1 << AROMATIC_MASK);
    mask6 = (1L << RING6);
    num = 0;
    for (i=0; i < 10; i++)
		inarray[i] = FALSE;
	
    for (i=0; i < 5; i++)
    {
        ia = array[i];
        for (j=i+1; j < 6; j++)
        {
            for (k=0; k < MAXIAT; k++)
            {
                if (atom[ia].iat[k] == array[j])
                {
                    if (atom[ia].bo[k] == 2)
                    {
						num++;
						inarray[i] = TRUE;
						inarray[j] = TRUE;
                    }
                }
            }
        }
    }
    if (num == 3)
		return TRUE;
	
    if (num >= 1)
    {
        for(i=0; i < 5; i++)
        {
            if (inarray[i] == FALSE)
            {
                ia = array[i];
                for(j=i+1;j < 6; j++)
                {
                    if (inarray[j] == FALSE)
                    {
                        ib = array[j];
                        if ( (atom[ia].flags & aromatic_mask) && (atom[ib].flags & aromatic_mask))
                        {
							num++;
							if (num == 3)
								return TRUE;
                        } else if (isbond(ia,ib))
                        {
							idbl = FALSE;
							kdbl = FALSE;
							for(k=0; k < MAXIAT; k++)
							{
								if (atom[ia].iat[k] !=0 && atom[ia].bo[k] == 2)
								{
									if (atom[atom[ia].iat[k]].atomnum == 6 || atom[atom[ia].iat[k]].atomnum == 7)
									{
										if (atom[atom[ia].iat[k]].flags & mask6)
											idbl = TRUE;
									}
								}
								if (atom[ib].iat[k] !=0 && atom[ib].bo[k] == 2)
								{
									if (atom[atom[ib].iat[k]].atomnum == 6 || atom[atom[ib].iat[k]].atomnum == 7)
									{
										if (atom[atom[ib].iat[k]].flags & mask6)
											kdbl = TRUE;
									}
								}
							}
							if (idbl == TRUE && kdbl == TRUE)
							{
								num++;
								if (num == 3)
									return TRUE;
							}
                        }
                        inarray[i] = TRUE;
                        inarray[j] = TRUE;
                    }
                }
            }
        }
    }
    return FALSE;
}
// ==============================
int check_connectivity(int ia)
{
	int nconnect,i,atomnum,input_charge;

  atomnum = atom[ia].atomnum;
  input_charge = atom[ia].formal_charge;
  nconnect = 0;
  for (i=0; i < MAXIAT; i++)
    {
      if (atom[ia].iat[i] != 0)
	nconnect++;
    }
  // C (x3)
  if (atomnum == 6 && nconnect == 3)
    return TRUE;
  // N (x2)
  else if (atomnum == 7 && nconnect == 2)
    return TRUE;
  // P (x2)
  else if (atomnum == 15 && nconnect == 2)
    return TRUE;
  // N+ (x2)
  else if (atomnum == 7 && nconnect == 3 && input_charge == 1)
    return TRUE;
  // P+ (x2)
  else if (atomnum == 15 && nconnect == 3 && input_charge == 1)
    return TRUE;
  // O+ (x2)
  else if (atomnum == 8 && nconnect == 2 && input_charge == 1)
    return TRUE;
  // S+ (x2)
  else if (atomnum == 16 && nconnect == 2 && input_charge == 1)
    return TRUE;

  return FALSE;
}
// ==============================
int check_connectivity_yatom(int ia)
{
  int nconnect,i,atomnum,input_charge;

  atomnum = atom[ia].atomnum;
  input_charge = atom[ia].formal_charge;
  nconnect = 0;
  for (i=0; i < MAXIAT; i++)
    {
      if (atom[ia].iat[i] != 0)
	nconnect++;
    }
  // C- (x3) cyclopentadiene anion
    if (atomnum == 6 && input_charge == -1 && nconnect == 3)
    return TRUE;
  // N- (x2)
  else if (atomnum == 7 && input_charge == -1 && nconnect == 2)
    return TRUE;
  // N  (x3)
  else if (atomnum == 7 && nconnect == 3)
    return TRUE;
  // O (x2)
  else if (atomnum == 8 && nconnect == 2)
    return TRUE;
  // S (x2)
  else if (atomnum == 16 && nconnect == 2)
    return TRUE;
  // P (x2)
  else if (atomnum == 15 && nconnect == 3)
    return TRUE;
  return FALSE;
}
//  ====================================
int has_double_bond(int ia)
 {
   int i;
   for (i=0; i < MAXIAT; i++)
     {
       if (atom[ia].iat[i] != 0 && atom[ia].bo[i] == 2)
	 return TRUE;
     }
   return FALSE;
 }


