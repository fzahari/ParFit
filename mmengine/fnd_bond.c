#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"
#include "gmmx.h"
#include "rotbond.h"
    
void find_bonds()
{
     int i,j,jji, ia1, ia2, k, l, found1, found2, ibres, nRc;
     int jj1, bondorder, bondres;
     
     rotbond.nbonds = 0;
     for (i=0; i < 200; i++)
     {
        rotbond.bond[i][0] = 0;
        rotbond.bond[i][1] = 0;
     }
     
     find_rings();

     bondorder = 1;
     if (rotbond.incl_alkene == TRUE)
        bondorder = 2;
     bondres = 180;
     if (rotbond.incl_amide == TRUE)
        bondres = 1000;
        
//      hdel(0);
     for (i=1; i <= natom; i++)
     {
         jji = 0;
         for (j=0; j < MAXIAT; j++)
         {
             if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].atomnum != 1 && atom[atom[i].iat[j]].mmx_type != 20)
               jji++;
         }
         if (jji > 1)
         {
             for (j=0; j < MAXIAT; j++)
             {
                 if (atom[i].iat[j] != 0 && atom[i].bo[j] <= bondorder)
                 {
                   ia1 = atom[i].iat[j];
                   jj1 = 0;
                   for (k=0; k < MAXIAT; k++)
                   {
                       if (atom[ia1].iat[k] != 0 && atom[ia1].bo[k] != 9 && atom[atom[ia1].iat[k]].atomnum != 1 && atom[atom[ia1].iat[k]].mmx_type != 20)
                          jj1++;
                   }
                   if (jj1 > 1)
                   {
                     if (i < atom[i].iat[j])
                     {
                         ia1 = i;
                         ia2 = atom[i].iat[j];
                     }else
                     {
                         ia1 = atom[i].iat[j];
                         ia2 = i;
                     }
                     if (atom[ia1].atomnum == 1)
                        goto L_10;
                     if (atom[ia2].atomnum == 1)
                        goto L_10;
                     // check to see if we have this bond
                     for (k=0; k < rotbond.nbonds; k++)
                     {
                         if (ia1 == rotbond.bond[k][0] && ia2 == rotbond.bond[k][1])
                            goto L_10;
                     }
                     // check to see if part of ring
                     for (k=0; k < gmmxring.nrings; k++)
                     {
                         found1 = FALSE;
                         found2 = FALSE;
                         for (l=0; l < gmmxring.nring_atoms[k]; l++)
                         {
                             if (ia1 == gmmxring.ring_atoms[k][l])
                                found1 = TRUE;
                             if (ia2 == gmmxring.ring_atoms[k][l])
                                found2 = TRUE;
                         }
                         if (found1 == TRUE && found2 == TRUE)
                            goto L_10;
                     }
                     // got here with new bond
                     nRc = check_bond(ia1,ia2,&ibres);
                     if (nRc == TRUE && ibres != bondres)
                     {
                        rotbond.bond[rotbond.nbonds][0] = ia1;
                        rotbond.bond[rotbond.nbonds][1] = ia2;
                        rotbond.bond_res[rotbond.nbonds] = ibres;
                        rotbond.nbonds++;
                        if (rotbond.nbonds > 200)
                        {
                            message_alert("Too many rotatable bonds","Find Bond Setup");
                            return;
                        }
                     }
                   }                    
                 }
L_10:
                 continue;
             }         
         }
     }
//     hadd();    
}
