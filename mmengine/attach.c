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


