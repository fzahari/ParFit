#include "pcwin.h"
#include "pcmod.h"
#include "bonds_ff.h"
#include "rings.h"
#include "pot.h"
#include "field.h"
#include "bondk.h"
#include "dipolek.h"
#include "dipole.h"

void kdipole()
{
    int i, j, ia, ib, ita, itb, it;
    int ierr, use_dipole;
    long int mask;
    char pa[4],pb[4],pt[7];

    mask = 1L << 0;
    use_dipole = FALSE;
    
    for (i=0; i < bonds_ff.nbnd; i++)
    {
        bonds_ff.bmom[i] = 0.0;
        ia = bonds_ff.i12[i][0];
        ib = bonds_ff.i12[i][1];
        bonds_ff.idpl[i][0] = -1;
        bonds_ff.idpl[i][1] = -1;
        ita = atom[ia].type;
        itb = atom[ib].type;
        numeral(ita,pa,3);
        numeral(itb,pb,3);
        if(ita < itb )
        {
            it = ita*100+itb;
            strcpy(pt,pa);
            strcat(pt,pb);
        } else
        {
            it = itb*100+ita;
            strcpy(pt,pb);
            strcat(pt,pa);
        }

        ierr = FALSE;
//  Do we need to check for transition state atoms
        if (bondpk.use_pibond == TRUE)
        {
            if ( (atom[ia].flags & mask) && (atom[ib].flags & mask) )
            {
                for (j=0; j < bondpk.npibond; j++)
                {
                    if (strcmp(pt,bondpk.kb[j]) == 0)
                    {
                        if (ita <= itb )
                        {
                            bonds_ff.idpl[i][0] = ia;
                            bonds_ff.idpl[i][1] = ib;
                        } else
                        {
                            bonds_ff.idpl[i][0] = ib;
                            bonds_ff.idpl[i][1] = ia;
                        }
                        bonds_ff.bmom[i] = bondpk.bmom[j];
                        ierr = TRUE;
                        use_dipole = TRUE;
                        break;
                    }
                }
            }
        }
        if (ierr == TRUE)
          goto L_10;
//   check for metals
//   check 5 membered rings
        if (rings.nring5 > 0 && bondk1.use_ring5 == TRUE)
        {
         for (j=0; j < dipole_k.ndipole5; j++)
         {
          if (strcmp(dipole_k.kb5[j],pt) == 0)
          {
             if (ita <= itb )
             {
                 bonds_ff.idpl[i][0] = ia;
                 bonds_ff.idpl[i][1] = ib;
              } else
              {
                 bonds_ff.idpl[i][0] = ib;
                 bonds_ff.idpl[i][1] = ia;
              }
              bonds_ff.bmom[i] = dipole_k.bmom5[j];
              ierr = TRUE;
              use_dipole = TRUE;
              break;
          }
         }
       }
      if (ierr == TRUE)
          goto L_10;
//   check 4 membered rings
        if (rings.nring4 > 0 && bondk1.use_ring4 == TRUE)
        {
         for (j=0; j < dipole_k.ndipole4; j++)
         {
          if (strcmp(dipole_k.kb4[j],pt) == 0)
          {
             if (ita <= itb )
             {
                 bonds_ff.idpl[i][0] = ia;
                 bonds_ff.idpl[i][1] = ib;
              } else
              {
                 bonds_ff.idpl[i][0] = ib;
                 bonds_ff.idpl[i][1] = ia;
              }
              bonds_ff.bmom[i] = dipole_k.bmom4[j];
              ierr = TRUE;
              use_dipole = TRUE;
              break;
          }
         }
       }
      if (ierr == TRUE)
          goto L_10;
//   check 3 membered rings
        if (rings.nring3 > 0 && bondk1.use_ring3 == TRUE)
        {
         for (j=0; j < dipole_k.ndipole3; j++)
         {
          if (strcmp(dipole_k.kb3[j],pt) == 0)
          {
             if (ita <= itb )
             {
                 bonds_ff.idpl[i][0] = ia;
                 bonds_ff.idpl[i][1] = ib;
              } else
              {
                 bonds_ff.idpl[i][0] = ib;
                 bonds_ff.idpl[i][1] = ia;
              }
              bonds_ff.bmom[i] = dipole_k.bmom3[j];
              ierr = TRUE;
              use_dipole = TRUE;
              break;
          }
         }
       }
      if (ierr == TRUE)
          goto L_10;

 // regular constants
         for (j=0; j < dipole_k.ndipole; j++)
         {
          if (strcmp(dipole_k.kb[j],pt) ==0)
          {
             if (ita <= itb )
             {
                 bonds_ff.idpl[i][0] = ia;
                 bonds_ff.idpl[i][1] = ib;
              } else
              {
                 bonds_ff.idpl[i][0] = ib;
                 bonds_ff.idpl[i][1] = ia;
              }
              bonds_ff.bmom[i] = dipole_k.bmom[j];
              ierr = TRUE;
              use_dipole = TRUE;
              break;
          }
        }
L_10:
        continue;
    }
//    if (use_dipole == FALSE)
//      pot.use_dipole = FALSE;

}

