#include "pcwin.h"
#include "pcmod.h"
#include "bonds_ff.h"
#include "field.h"
#include "utility.h"
#include "atom_k.h"
#include "pot.h"
#include "elements.h"
#include "bondk.h"
#include "chargek.h"
#include "dipolek.h"
#include "gastype.h"
#include "echarge.h"

void kcharge()
{
   int i, j, ia, ib, it, iit, kit, jji;
   long int mask;
   char pa[4],pb[4],pt[7];
   float bmomj;
   
   mask = 1L << 0;
   
   for (i=1; i <= natom; i++)
   {
      atom[i].sigma_charge = 0.0;
      if (atom[i].type < 300)
        atom[i].charge = 0.0;
   }

   if (field.type == MMFF94)
   {
       compute_fcharge();
       for(i=0; i < bonds_ff.nbnd; i++)
       {
           ia = bonds_ff.i12[i][0];
           ib = bonds_ff.i12[i][1];
           iit = atom[ia].type;
           kit = atom[ib].type;
           if (iit < kit)
             it = iit*100+kit;
           else
             it = kit*100+iit;
           if (bonds_ff.index[i] == 1)
           {
              for(j=0; j < charge_k.nbndchrgdel; j++)
              {
                 if (it == charge_k.btypedel[j])
                 {
                   if (iit < kit)
                   {
                     atom[ia].charge -= charge_k.bchargedel[j];
                     atom[ib].charge += charge_k.bchargedel[j];
                   }else
                   {
                     atom[ib].charge -= charge_k.bchargedel[j];
                     atom[ia].charge += charge_k.bchargedel[j];
                   }
                   break;
                 }
              }
           } else
           {
              for(j=0; j < charge_k.nbndchrg; j++)
              {
                 if (it == charge_k.btype[j])
                 {
                   if (iit < kit)
                   {
                     atom[ia].charge -= charge_k.bcharge[j];
                     atom[ib].charge += charge_k.bcharge[j];
                   }else
                   {
                     atom[ib].charge -= charge_k.bcharge[j];
                     atom[ia].charge += charge_k.bcharge[j];
                   }
                   break;
                 }
             }
           }
       }
       // add in formal charges
       for (i=1; i <= natom; i++)
       {
           if (atom[i].atomnum == 15)
               atom[i].formal_charge = 0.0;
           if (atom[i].atomnum == 16 && atom[i].type != 72)
               atom[i].formal_charge = 0.0;             
       }              
       for (i=1; i <= natom; i++)
       {
           atom[i].charge += (1.0 - atom_k.ligands[atom[i].type]*charge_k.formchrg[atom[i].type])*atom[i].formal_charge ;
           for(j=0; j < MAXIAT; j++)
           {
               if (atom[i].iat[j] != 0)
               {
                  if (atom[atom[i].iat[j]].formal_charge < 0.0)
                           atom[i].charge += charge_k.formchrg[atom[atom[i].iat[j]].type]*atom[atom[i].iat[j]].formal_charge;
               }
           }
       }
   } else if ((field.type == MM3 && pot.use_dipole == FALSE) )
   {
      // compute charges from bond moments
      for (i=0; i < bonds_ff.nbnd; i++)
      {
          ia = bonds_ff.i12[i][0];
          ib = bonds_ff.i12[i][1];
          iit = atom[ia].type;
          kit = atom[ib].type;
         if( field.type == MMX)
         {
             if (iit >= 300)
                iit = 300;
             if( kit >= 300)
                kit = 300;
             if (iit == 40 && ( kit != 2 && kit != 3 && kit != 4 && kit != 40) )
                iit = 2;
             if (kit == 40 && ( iit != 2 && iit != 3 && iit != 4 && iit != 40) )
                kit = 2;
         }
          numeral(iit,pa,3);
          numeral(kit,pb,3);
          
          if (iit < kit)
          {
              strcpy(pt,pa);
              strcat(pt,pb);
          } else
          {
              strcpy(pt,pb);
              strcat(pt,pa);
          }
          if( (atom[ia].flags & mask) && (atom[ib].flags & mask) )
          {  // search through pi constants
             bmomj = 0.0;
             for (j=0; j < bondpk.npibond; j++)
             {
                 if (strcmp(pt, bondpk.kb[j]) == 0)
                 {
                    bmomj = bondpk.bmom[j];
                    break;
                 }
              }
          } else
          {
              bmomj = 0.0;
              for (j=0; j < dipole_k.ndipole; j++)
              {
                  if (strcmp(pt,dipole_k.kb[j]) == 0)
                  {
                      bmomj = dipole_k.bmom[j];
                      break;
                  }
              }
          }
          if (bmomj != 0.0)
          {
              bmomj = 0.42*bmomj/bonds_ff.bl[i];
              if (iit <= kit)
              {
                  atom[ia].charge += 0.5*bmomj;
                  atom[ib].charge -= 0.5*bmomj;
              } else
              {
                  atom[ia].charge -= 0.5*bmomj;
                  atom[ib].charge += 0.5*bmomj;
              }
          }
          if (kit == 66)
          {
              if (iit == 3)
                atom[ib].charge -= 0.5;
              if (iit == 18)
                atom[ib].charge -= 0.33;
          } else if (iit == 66)
          {
              if (kit == 3)
                atom[ia].charge -= 0.5;
              if (kit == 18)
                atom[ia].charge -= 0.33;
          }
      }
      // check for isolated charged atoms
      for (i=1; i <= natom; i++)
      {
          if ( !(atom[i].flags & mask))
          {
              if (atom[i].type == 16 || atom[i].type == 30)
                  atom[i].charge += 1.0;
              if (atom[i].type == 41)  // N+
              {
                  for(j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0)
                      {
                          if (atom[atom[i].iat[j]].type == 42 || atom[atom[i].iat[j]].type == 4
                          || atom[atom[i].iat[j]].type == 66)
                             goto L_10;
                      }
                  }
                  atom[i].charge += 1.0;
              }
              if (atom[i].type == 46) // O+
              {
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].bo[j] == 3) // co in metal carbonyl
                        goto L_10;
                  }                 
                  atom[i].charge += 1.0;
              }
              if (atom[i].type == 27 || atom[i].type == 42 || atom[i].type == 48)
                 atom[i].charge -= 1.0;
              if (atom[i].type >= 11 && atom[i].type <= 14)
              {
                  jji = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0)
                        jji++;
                  }
                  if (jji == 0)
                    atom[i].charge -= 1.0;
              }
              if (atom[i].type == 25 || atom[i].type == 47)
              {
                  jji = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0)
                        jji++;
                  }
                  if (jji == 6)
                    atom[i].charge -= 1.0;
              }
              if (atom[i].type == 58) // Al
                 atom[i].charge += 1.0;
                     
          }
L_10:
          continue;
      }
      // search list of charged atoms read in
      for (i=1; i <= natom; i++)
      {
                for(j=0; j < charge_k.ncharge; j++)
        {
                        if (atom[i].type == charge_k.type[j])
                        {
              atom[i].charge = charge_k.charge[j];
              break;
                        }                   
        }
      }     
   } else
   {
      // search list of charges
       for (i=1; i <= natom; i++)
       {
           for(j=0; j < charge_k.ncharge; j++)
           {
               if (atom[i].type == charge_k.type[j])
               {
                   atom[i].charge = charge_k.charge[j];
                   break;
               }                   
           }
       }
   }
   for (i=1;i <= natom; i++)
      atom[i].sigma_charge = atom[i].charge; 
}
/* -------------------------------------------------------- */
void compute_fcharge()
{
    int i, j, jjbo;
    int jatm, k, jji, attached[4];
    int *used;
    float charge;
    
    used = ivector(1,natom+1);
    for (i=1; i <= natom; i++)
       used[i] = FALSE;
       
    for (i=1; i <= natom; i++)
    {
        atom[i].formal_charge = charge_k.typechrg[atom[i].type];
    }
//  find variable charges: types 32, 72, 76, 81
    for (i=1; i <= natom; i++)
    {
        if (used[i] == FALSE)
        {
            if (atom[i].formal_charge != 0.0)
            {
               if (atom[i].type == 32)
               {
                   jatm = atom[i].iat[0];
                   jji = 0;
                   jjbo = 0;
                   for (k=0; k < MAXIAT; k++)
                   {
                       if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9)
                       {
                           if (atom[atom[jatm].iat[k]].type == 32)
                           {
                               if (atom[jatm].bo[k] == 1)
                               {
                                   jjbo++;
                                   attached[jji] = atom[jatm].iat[k];
                                   jji++;
                               } else
                               {
                                   attached[jji] = atom[jatm].iat[k];
                                   jji++;
                               }
                           }
                       }
                   }
                   charge = -(double)jjbo/jji;
                   if (jji == 1)  // only 1 type 32
                     atom[i].formal_charge = 0.0;
                   else if (jji == 2)
                   {
                       if (atom[jatm].atomnum == 15)  // p-o  treat as single negative charge;
                          charge = -1.00/jji;
                       if (atom[jatm].type == 73)
                       {
                           charge = -0.500;
                           atom[attached[0]].formal_charge = charge;
                           atom[attached[1]].formal_charge = charge;
                           used[attached[0]] = TRUE;                           
                           used[attached[1]] = TRUE;
                       } else if (atom[jatm].atomnum != 7)
                       {
                           atom[attached[0]].formal_charge = charge;
                           atom[attached[1]].formal_charge = charge;
                           used[attached[0]] = TRUE;                           
                           used[attached[1]] = TRUE;
                       }else
                       {
                           atom[attached[0]].formal_charge = 0.0;
                           atom[attached[1]].formal_charge = 0.0;
                           used[attached[0]] = TRUE;                           
                           used[attached[1]] = TRUE;
                       }                         
                   } else if (jji == 3)
                   {
                        atom[attached[0]].formal_charge = charge;
                        atom[attached[1]].formal_charge = charge;
                        atom[attached[2]].formal_charge = charge;
                        used[attached[0]] = TRUE;                           
                        used[attached[1]] = TRUE;                           
                        used[attached[2]] = TRUE;                           
                   } else if (jji == 4)
                   {
                        atom[attached[0]].formal_charge = charge;
                        atom[attached[1]].formal_charge = charge;
                        atom[attached[2]].formal_charge = charge;
                        atom[attached[3]].formal_charge = charge;
                        used[attached[0]] = TRUE;                           
                        used[attached[1]] = TRUE;                           
                        used[attached[2]] = TRUE;                           
                        used[attached[3]] = TRUE;                           
                   }
               }else if (atom[i].type == 72)
               {
                   jatm = atom[i].iat[0];
                   jji = 0;
                   jjbo = 0;
                   for (k=0; k < MAXIAT; k++)
                   {
                       if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9)
                       {
                           if (atom[atom[jatm].iat[k]].type == 72)
                           {
                               if (atom[jatm].bo[k] == 1)
                               {
                                   jjbo++;
                                   attached[jji] = atom[jatm].iat[k];
                                   jji++;
                               } else
                               {
                                   attached[jji] = atom[jatm].iat[k];
                                   jji++;
                               }
                           }
                       }
                   }
                   if (jji == 1)
                   {
                       if (atom[jatm].atomnum == 15)
                          atom[i].formal_charge = 0.0;
                       else
                          atom[i].formal_charge = -1.00;
                   }else
                   {
                       atom[attached[0]].formal_charge = -0.5;
                       atom[attached[1]].formal_charge = -0.5;
                       used[attached[0]] = TRUE;
                       used[attached[1]] = TRUE;
                   }
               } else if (atom[i].type == 81)
               {
                   jji = 0;
                   for (k=0; k < MAXIAT ; k++)
                   {
                       if (atom[i].iat[k] != 0 && atom[atom[i].iat[k]].type == 80)
                       {
                           jatm = atom[i].iat[k];
                           jji = 0;
                           for (j=0; j < MAXIAT; j++)
                           {
                               if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] != 9)
                               {
                                   if (atom[atom[jatm].iat[j]].atomnum == 7)
                                   {
                                       attached[jji] = atom[jatm].iat[j];
                                       jji++;
                                   }
                               }
                           }
                       }
                   }
                   if (jji == 0)
                   {
                     charge = 1.00;
                     atom[i].formal_charge = charge;
                     used[i] = TRUE;
                   } else
                   {
                      charge = 1.00/jji;
                      for (j=0; j < jji; j++)
                      {
                          atom[attached[j]].formal_charge = charge;
                          used[attached[j]] = TRUE;
                      }
                   }
               }
            }
        }
    }
    free_ivector(used,1, natom+1);    
}
/* ====================================== */
void assign_gast_charge()
{
    int i,j,t,cycle;
    double d1,d2,z1,q1;
    double *xx, *atomchrg;
    Gastype *atm_par;


    atomchrg = dvector(0,natom+1);
    xx = dvector(0,natom+1);
    atm_par = (Gastype *)malloc( (natom+1)* sizeof(Gastype));

    for (i=1; i <= natom; i++)
    {
       if (atom[i].mmx_type < 300)
       {
         atm_par[i].a = Gastchrg[atom[i].mmx_type-1].a;
         atm_par[i].b = Gastchrg[atom[i].mmx_type-1].b;
         atm_par[i].c = Gastchrg[atom[i].mmx_type-1].c;
         atm_par[i].d = atm_par[i].a + atm_par[i].b + atm_par[i].c;
         d2 = atm_par[i].d;
         if (fabs(atm_par[i].d) < 0.01)
           atm_par[i].d = 1.0;
       }
    }
    
    for (i=1; i <= natom; i++)
    {
        atomchrg[i] = 0.0;
        if (atom[i].mmx_type == 30 ) // C+
           atomchrg[i] += 1.0;
        else if (atom[i].mmx_type == 41) // N+
           atomchrg[i] += 1.0;
        else if (atom[i].mmx_type == 42) // O-
           atomchrg[i] -= 1.0;
        else if (atom[i].mmx_type == 48) // C-
           atomchrg[i] -= 1.0;
        else if (atom[i].mmx_type == 66) // O-
           atomchrg[i] -= 0.5;
        else if (atom[i].mmx_type >= 300) // Metals 
           atomchrg[i] = atom[i].charge;
        xx[i] = atm_par[i].a;
    }
    z1 = 1.0;

    cycle = 0;
    do {
        z1 *= 0.5;
        d1 = 0.0;
        for (i=1; i <= natom; i++)
        {
            if (atm_par[i].d != 1.0 && atom[i].mmx_type < 300)
            {
                q1 = atomchrg[i];
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 )
                    {
                        t = atom[i].iat[j];
                        if (atm_par[t].d != 1.0 && atom[t].mmx_type < 300)
                        {
                           d2 = atm_par[t].d;
                           if (xx[t] > xx[i])
                             d2 = atm_par[i].d;

                           if (atom[atom[i].iat[0]].atomnum == 1)
                              d2 = 20.02;
                           if (atom[atom[t].iat[0]].atomnum == 1)
                              d2 = 20.02;

                           atomchrg[i] += (xx[t] - xx[i])/d2*z1;
                        }
                    }
                }
                q1 = fabs(atomchrg[i] - q1);
                if (q1 > d1)
                   d1 = q1;
            }
        }
        if (d1 >= 0.001)
        {
            for (i=1; i <= natom; i++)
               if (atom[i].mmx_type < 300)
                 xx[i] = atm_par[i].a + atm_par[i].b*atomchrg[i] + atm_par[i].c*atomchrg[i]*atomchrg[i];
        }
        cycle++;
    }while ( (d1 > 0.001) && cycle <= 5);

    for (i=1; i <= natom; i++)
       if (atom[i].mmx_type < 300)
           atom[i].charge = atomchrg[i];
      
    free(atm_par);
    free_dvector(atomchrg ,0,natom+1);
    free_dvector(xx ,0,natom+1);
    
}
// sort in reverse order on first array, keep other arrays in sync
//
void ksortr(int num, int array[], int array1[],int array2[])
{
    int i,temp;
    int found;

L_1:
    found = FALSE;
    for (i=0; i < num-1; i++)
    {
       if (array[i+1] > array[i])
       {
          temp = array[i];
          array[i] = array[i+1];
          array[i+1] = temp;
          
          temp = array1[i];
          array1[i] = array1[i+1];
          array1[i+1] = temp;

          temp = array2[i];
          array2[i] = array2[i+1];
          array2[i+1] = temp;

          found = TRUE;
       }
    }
    if (found == TRUE)
       goto L_1;
}
/* ==================================================== */
void ksort_ring(int ia1,int size,int *array)
{
    int i,j,ipoint;
    int tmp1[6], tmp2[6];
    int best;

    ipoint = 0;
    for (i=0; i < size; i++)
    {
        if (ia1 == array[i])
        {
            ipoint = i;
            break;
        }
    }

    for (i=0; i < size; i++)
    {
        j = ipoint + i;
        if (j >= size) j -= size;
        
        tmp1[i] = array[j];
    }
//    
    tmp2[0] = tmp1[0];
    for (i=1; i < size; i++)
        tmp2[i] = tmp1[size-i];
//
    best = 1;
    for (i=0; i < size; i++)
    {
        if (atom[tmp2[i]].atomnum > atom[tmp1[i]].atomnum)
        {
            best = 2;
            break;
        } else if (atom[tmp2[i]].atomnum < atom[tmp1[i]].atomnum)
        {
            best = 1;
            break;
        }
    }
    if (best == 1)
    {
        for (i=0; i < size; i++)
           array[i] = tmp1[i];
    } else
    {
        for (i=0; i < size; i++)
           array[i] = tmp2[i];
    }
}
