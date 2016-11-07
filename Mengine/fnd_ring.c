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
#include "utility.h"
#include "gmmx.h"
#include "find_ring.h"

#define MAX_RINGSIZE 30
  
void find_rings()
{
    int i,j,k,l,kk, jji;
    int found_terminal;
    int ia,ib,icount,found,tmp_size, tmp_ring[50];
    int *iexist, **itemiat;


    iexist = ivector(1, natom+1);
    itemiat = imatrix(1,natom+1, 0,MAXIAT);
    gmmxring.nrings = 0;
    for (i=0; i < 200; i++)
    {
        gmmxring.nring_atoms[i] = 0;
        gmmxring.fused[i] = -1;
        for (j=0; j < 30; j++)
           gmmxring.ring_atoms[i][j] = 0;
    }

    for (i=1; i <= natom; i++)
    {
        iexist[i] = atom[i].type;
        for(j=0; j < MAXIAT; j++)
            itemiat[i][j] = atom[i].iat[j];
    }

L_1:
    found_terminal = FALSE;
    for (i=1; i <= natom; i++)
    {
        if (iexist[i] != 0)
        {
            jji = 0;
            for(j=0; j < MAXIAT; j++)
            {
                if (itemiat[i][j] != 0)
                   jji++;
            }
            if (jji <= 1)
            {
                del_atom(i,itemiat);
                iexist[i] = 0;
                found_terminal = TRUE;
            }
        }
    }
    if (found_terminal == TRUE)
       goto L_1;
    find_rings1(iexist, itemiat);
    free_ivector(iexist, 1, natom+1);
    free_imatrix(itemiat, 1, natom+1, 0, MAXIAT);
// sort on ring size
L_2:
    found = FALSE;
    for (i=0; i < gmmxring.nrings-1; i++)
    {
        for (j=i+1; j < gmmxring.nrings; j++)
        {
            if (gmmxring.nring_atoms[j] < gmmxring.nring_atoms[i])
            {
                tmp_size = gmmxring.nring_atoms[i];
                for (k=0; k < tmp_size; k++)
                   tmp_ring[k] = gmmxring.ring_atoms[i][k];
                   
                gmmxring.nring_atoms[i] = gmmxring.nring_atoms[j];
                for (k=0; k < gmmxring.nring_atoms[j]; k++)
                   gmmxring.ring_atoms[i][k] = gmmxring.ring_atoms[j][k];
                   
                gmmxring.nring_atoms[j] = tmp_size;
                for (k=0; k < tmp_size; k++)
                   gmmxring.ring_atoms[j][k] = tmp_ring[k];
                found = TRUE;
            }
        }
    }
    if (found == TRUE)  goto L_2;
// find fused rings
//  now find branch atom pairs and what rings they occur in
     gmmxring.npairs = 0;
     for (i=0; i < 200; i++)
     {
         gmmxring.pair[i][0] = -1;
         gmmxring.pair[i][1] = -1;
         for (j=0; j < 10; j++)
             gmmxring.in_rings[i][j] = -1;
     }              
     for (i=0; i < gmmxring.nbranch; i++)
     {
         icount = 0;
         ia = gmmxring.branch_atoms[i];
         for (j=0; j < gmmxring.nrings; j++)
         {
             if (in_ring(ia,j))            
             {
                 gmmxring.in_rings[i][icount] = j;
                 icount++;
             }
         }
         // if icount == 1 we should remove - not a branch atom
     }
//  found which rings each branch atom belongs to
     for (i=0; i < gmmxring.nbranch; i++)
     {
         for (k=0; gmmxring.in_rings[i][k] != -1; k++)
         {
             ia = gmmxring.in_rings[i][k];
             for (kk=k+1; gmmxring.in_rings[i][kk] != -1; kk++)
             {
                ib = gmmxring.in_rings[i][kk];
                for (j=i+1; j < gmmxring.nbranch; j++)
                {
                   for (l=0; gmmxring.in_rings[j][l] != -1; l++)
                   {
                     if (ia == gmmxring.in_rings[j][l] && ib == gmmxring.in_rings[j][l+1])
                         add_pair(i,j,ia,ib);
                   }
                }
             }
         }
     }
// now group fused systems
     gmmxring.nfused = 0;
     for (i=0; i < gmmxring.npairs; i++)
     {
         found = FALSE;
         ia = gmmxring.pair[i][2];
         ib = gmmxring.pair[i][3];
         for (j=i+1; j < gmmxring.npairs; j++)
         {
             if (ia == gmmxring.pair[j][2] || ia == gmmxring.pair[j][3] || ib == gmmxring.pair[j][2] || ib == gmmxring.pair[j][3])
             {
                 if (gmmxring.fused[ia] == -1)
                 {
                     found = TRUE;
                     gmmxring.fused[ia] = gmmxring.nfused;
                     gmmxring.fused[ib] = gmmxring.nfused;
                     gmmxring.fused[gmmxring.pair[j][2]] = gmmxring.nfused;
                     gmmxring.fused[gmmxring.pair[j][3]] = gmmxring.nfused;
                     gmmxring.nfused++;
                 } else
                 {
                     found = TRUE;
                     gmmxring.fused[ib] = gmmxring.fused[ia];
                     gmmxring.fused[gmmxring.pair[j][2]] = gmmxring.fused[ia];
                     gmmxring.fused[gmmxring.pair[j][3]] = gmmxring.fused[ia];
                 }
             }
         }
         if (found == FALSE && gmmxring.fused[ia] == -1)
         {
             gmmxring.fused[ia] = gmmxring.nfused;
             gmmxring.fused[ib] = gmmxring.nfused;
             gmmxring.nfused++;
         }
     }
}
/* ==================================== */
void del_atom(int i, int **itemiat)
{
    int j, itiat, k;

    for (j=0; j < MAXIAT; j++)
    {
        if (itemiat[i][j] != 0)
        {
            itiat = itemiat[i][j];
            for (k=0; k < MAXIAT; k++)
            {
                if (itemiat[itiat][k] == i)
                  itemiat[itiat][k] = 0;
            }
        }
    }
    for (j=0; j < MAXIAT; j++)
       itemiat[i][j] = 0;
}
/* ==================================== */
void add_all_rings(int j,int *tmp)
{
    int i,ia,k;
    for (i=0; i < 10; i++)
    {
        if (gmmxring.in_rings[j][i] != -1)
        {
            ia = gmmxring.in_rings[j][i];
            for (k=0; k < 10; k++)
            {
                if (tmp[k] == ia)
                   break;
                else if (tmp[k] == -1)
                {
                    tmp[k] = ia;
                    break;
                }
            }
        }
    }
}
/* =================================================== */
int in_ring(int ia,int k)
{
    int i;
    for (i=0; i < gmmxring.nring_atoms[k]; i++)
    {
        if (ia == gmmxring.ring_atoms[k][i])
          return TRUE;
    }
    return FALSE;
}
/* =================================================== */
void add_pair(int ia,int ib,int k,int l)
{
    int tmpa, tmpb;

    if (ia < ib)
    {
        tmpa = ia;
        tmpb = ib;
    } else
    {
        tmpa = ib;
        tmpb = ia;
    }

// get here pair is not found create new pair
    gmmxring.pair[gmmxring.npairs][0] = tmpa;
    gmmxring.pair[gmmxring.npairs][1] = tmpb;
    gmmxring.pair[gmmxring.npairs][2] = k;
    gmmxring.pair[gmmxring.npairs][3] = l;
    gmmxring.npairs++;
    return;
}
// ============================================
void find_rings1(int *iexist,int **itemiat)
{
    int i,j,k,ia,ib;
    int rsize,ncount,numbond,nring,ringsize;
    int *newatom,*atnum,*jj;
//    int atbond[50][6];
    int **atbond;
    int *ring;
    int n2count,ntrimset;
    int *nodeN2,*trimSet,*nodesDone;
    int smallestdegree,smallest;
    int nodesToBreak;
    int newring;

    n2count = 0;
    gmmxring.nrings = 0;

// count size of ring system
    ncount = 0;
    for (i=1; i <= natom; i++)
       if (iexist[i] != 0)
          ncount++;

// allocate space
    rsize = ncount;
    newatom = ivector(0,natom+1);
    atnum = ivector(0,rsize);
    jj = ivector(0,rsize);
    ring = ivector(0,rsize);
    nodeN2 = ivector(0,rsize);
    trimSet = ivector(0,rsize);
    nodesDone = ivector(0,rsize);
    atbond = imatrix(0,rsize+2,0,6);
 
    ncount = 0;
    for (i=1; i <= natom; i++)
    {
        if (iexist[i] != 0)
        {
            atnum[ncount] = i;
            newatom[i] = ncount;
            ncount++;
        }
    }
//  count attachments and get bridgehead atoms
    gmmxring.nbranch = 0;
    numbond = 0;
    for (i=0; i < rsize; i++)
    {
       nodeN2[i] = -1;
       for (j=0; j < 6;j++)
         atbond[i][j] = -1;
    }
    for (i=0; i < rsize; i++)
    {
        ia = atnum[i];
        ncount = 0;
        n2count = 0;
        for (j=0; j < MAXIAT; j++)
        {
            if (itemiat[ia][j] != 0)
            {
                atbond[i][ncount] = newatom[itemiat[ia][j]];
                ncount++;
                if (ia < itemiat[ia][j])
                   numbond++;
            }
        }
        jj[i] = ncount;
        if (ncount >= 3)
        {
            gmmxring.branch_atoms[gmmxring.nbranch] = ia;
            gmmxring.nbranch++;
        } else if (ncount == 2)
        {
            nodeN2[n2count] = i;
            n2count++;
        }
            
    }
// number of rings
    nring = numbond - rsize + 1;
// loop over N2
    ntrimset = 0;
    for (i=0; i < rsize; i++)
       trimSet[i] = -1;
    do {

     // check current connectivities
        smallestdegree = 7;
        smallest = -1;
        for (i=0; i < rsize; i++)
        {
            ncount = 0;
            for (j=0; j < 6; j++)
            {
                if (atbond[i][j] != -1)
                   ncount++;
            }
            jj[i] = ncount;
            
            if (jj[i] == 1)
            {
               trim(i,&ntrimset,trimSet);
               if (ntrimset > rsize)
                  message_alert("Error in TrimSet","Error");
               for (j=0; j < 6; j++)
               {
                   if (atbond[i][j] != -1)
                   {
                       ia = atbond[i][j];
                       for (k=0; j < 6; k++)
                       {
                           if (atbond[ia][k] == i)
                           {
                               atbond[i][j] = -1;
                               atbond[ia][k] = -1;
                               break;
                           }
                       }
                   }
               }
            } else if (jj[i] == 2)
            {
                check_N2(i,&n2count,nodeN2);
                if (n2count > rsize)
                   message_alert("Error in check_N2","Error");
            }
           
            if (jj[i] < smallestdegree && jj[i] > 0)
            {
                smallest = i;
                smallestdegree = jj[i];
            }
        }
        if (smallest == -1)
           break;
        for(i=0; i < rsize; i++)
          ring[i] = -1;
        if (smallestdegree == 2)
        {
            nodesToBreak = 0;
            for (i=0; i < n2count; i++)
            {
                newring = getRing(nodeN2[i],rsize,jj,atbond,ring);
                if (newring)
                {
                    for (j=0; ring[j] != -1 && j < rsize; j++)
                           ring[j] = atnum[ring[j]];
                    ringsize = j;
                    if (ringsize > rsize)
                      message_alert("Error in getRing","Error");
                    log_ring1(ringsize,ring);
                }
                nodesDone[nodesToBreak] = nodeN2[i];
                nodesToBreak++;
            }
            for (i=0; i < nodesToBreak; i++)
            {
                breakBond(nodesDone[i],atbond);
                nodeN2[i] = -1;
            }
            n2count = 0;
        } else if (smallestdegree == 3)
        {
                newring = getRing(smallest,rsize,jj,atbond,ring);
                if (newring)
                {
                    ia = ring[0];
                    ib = ring[1];
                    
                    for (j=0; ring[j] != -1 && j < rsize; j++)
                           ring[j] = atnum[ring[j]];
                    ringsize = j;
                    log_ring1(ringsize,ring);
                    delete_bond(ia,ib,atbond);
                }
        }
        repack(rsize,atbond);       
    } while (ntrimset < rsize);
// free space
    free_ivector(newatom ,0,natom+1);
    free_ivector(atnum ,0,rsize);
    free_ivector(jj ,0,rsize);
    free_ivector(ring ,0,rsize);
    free_ivector(nodeN2 ,0,rsize);
    free_ivector(trimSet ,0,rsize);
    free_ivector(nodesDone ,0,rsize);
    free_imatrix(atbond ,0,rsize+2,0,6);
}
// =============  queque stuff =======================
#define MAX_QUE_SIZE 100
int  queue[MAX_QUE_SIZE];

void addq(int front,int *rear,int item)
{
    *rear = (*rear+1)%MAX_QUE_SIZE;
    if (front == *rear)
    {
        message_alert("Error in Queue","Error");
        return;
    }
    queue[*rear] = item;
}
int deleteq(int *front,int rear)
{
    if (*front == rear)  return(-1);
    *front = (*front+1)%MAX_QUE_SIZE;
    return queue[*front];
}
// ===================================================
int getRing(int ia,int rsize,int *jj, int **atbond,int *ring)
{
    int i,j,ib,ic,intersect,itmp;
    int front,rear;
    int **path;

    path = imatrix(0,rsize,0,rsize);
    for (i=0; i < rsize; i++)
    {
        for (j=0; j < rsize; j++)
          path[i][j] = -1;
    }
    front = rear = -1;
// get all neighbors of root node, put them on queque
// and update path
    for (i=0; i < jj[ia]; i++)
    {
        itmp = atbond[ia][i];
        if (jj[atbond[ia][i]] > 0)
        {
            addq(front,&rear,atbond[ia][i]);
            path[ia][0] = ia;
            path[atbond[ia][i]][0] = ia;
            path[atbond[ia][i]][1] = atbond[ia][i];
        }
    }
    while (front >= -1)
    {
        ib = deleteq(&front,rear); // get new node from queue
        if (ib == -1)  // empty queue
        {
            free_imatrix(path,0,rsize,0,rsize);
            return FALSE;
        }
        for (i=0; i < jj[ib]; i++)
        {
            ic = atbond[ib][i];
            if (path[ic][0] != -1) // already visited
            {
                for (j=0; path[ib][j] != -1; j++)
                        ;
                if (path[ib][j-2] != ic) // found collision
                {
                    intersect = get_intersection(path[ib],path[ic]);
                    if (intersect == 1)
                    {
                        for(j=0; j < rsize; j++)
                          ring[j] = -1;
                       getUnion(path[ib],path[ic],ring);
                       free_imatrix(path,0,rsize,0,rsize);
                       return TRUE;
                    }
                }
            } else  // put on queue
            {
                addq(front,&rear,ic);
                for (j=0;path[ib][j] != -1; j++)
                   path[ic][j] = path[ib][j];
                path[ic][j] = ic;
            }
        }
            
    }
    free_imatrix(path,0,rsize,0,rsize);
    return FALSE;   
}            
// ======================================
int get_intersection(int *path1, int *path2)
{
    int i,j,ia,icount;
    icount = 0;
    for (i=0; path1[i] != -1; i++)
    {
        ia = path1[i];
        for (j=0; path2[j] != -1; j++)
        {
            if (path2[j] == ia)
              icount++;
        }
    }
    return icount;
}   
// ======================================
void getUnion(int *path1, int *path2,int *ring)
{
    int i,j,k,icount;
    icount = 0;

    for (i=0; path1[i] != -1; i++)
    {
        ring[icount] = path1[i];
        icount++;
    }
    for (j=0; path2[j] != -1; j++)
        ;
    for (k=j-1; k > 0; k--)
    {
        ring[icount] = path2[k];
        icount++;
    }
}
// ======================================
void log_ring1(int iatom, int *path)
{
    int i, isame, j;
    int iarray1[50], iarray2[50];

    if (iatom > MAX_RINGSIZE)
    {
        message_alert("Ring too large in find ring","Error");
        return;
    }
    for (i=0; i < gmmxring.nrings; i++)
    {
        if (gmmxring.nring_atoms[i] == iatom)
        {
            isame = FALSE;
            for (j=0; j < 50; j++)
            {
                iarray1[j] = 0;
                iarray2[j] = 0;
            }
            for(j=0; j < iatom; j++)
               iarray1[j] = path[j];
            for(j=0; j < gmmxring.nring_atoms[i]; j++)
               iarray2[j] = gmmxring.ring_atoms[i][j];

            isame = icompare(iatom,iarray1, iarray2);
            if (isame == TRUE)
               return;
        }
    }
    gmmxring.nring_atoms[gmmxring.nrings] = iatom;
    for (i=0; i < iatom; i++)
        gmmxring.ring_atoms[gmmxring.nrings][i] = path[i];
    gmmxring.nrings++;
}
// =================================
void delete_bond(int ia,int ib, int **atbond)
{
    int i;

    for (i=0; i < 6; i++)
    {
      if (atbond[ia][i] == ib)
                atbond[ia][i] = -1;
      if (atbond[ib][i] == ia)
                atbond[ib][i] = -1;
    }
}
// =================================
void check_N2(int ia,int *n2count,int *nodes)
{
    int i;

    for (i=0; i < *n2count; i++)
    {
        if (nodes[i] == ia) // already there
           return;
    }
    nodes[*n2count] = ia;
    *n2count += 1;
}
// =================================
void trim(int ia,int *j,int *trimset)
{
    int i;

    for (i=0; i < *j; i++)
    {
        if (trimset[i] == ia)
           return;
    }
    trimset[*j] = ia;
    *j += 1;
}
// =================================
void breakBond(int ia,int **atbond)
{
    int i,j,ib;

    for (i=0; i < 6; i++)
    {
        if (atbond[ia][i] != -1)
        {
            ib = atbond[ia][i];
            for (j=0; j < 6; j++)
            {
                if (atbond[ib][j] == ia)
                {
                    atbond[ia][i] = -1;
                    atbond[ib][j] = -1;
                    return;
                }
            }
        }
    }
}
// =====================================
void repack(int rsize,int **atbond)
{
    int i,j,k;

    for (i=0; i < rsize; i++)
    {
        for (j=0; j < 5; j++)
        {
            if (atbond[i][j] == -1)
            {
                for (k=j+1; k < 6; k++)
                    atbond[i][k-1] = atbond[i][k];
                atbond[i][5] = -1;
            }
        }
    }
}
                    
