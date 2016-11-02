#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "cutoffs.h"
#include "pot.h"
#include "units.h"
#include "minim_values.h"
#include "skip.h"
#include "echrgdpl.h"

void echrg_dipole()
{
    int i,j,k,k1,k2;
    double cc, e,rk2,rkr3,dotk;
    double debye, f,fi,fik,cutoff;
    double xk,yk,zk;
    double xr,yr,zr,r2;
    
    cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
    cc = 4.80298*4.80298*14.39418;
    debye = 4.8033324;
    f = -cc/(debye*units.dielec);
    for (i=1; i <= natom; i++)
       skip[i] = 0;
      
    for (i=1; i <= natom; i++)
    {
        if (atom[i].charge != 0.0)
        {
            skip[i] = i;
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                   skip[atom[i].iat[j]] = i;
            }
            fi = f*atom[i].charge;
            for (k=0; k < bonds_ff.nbnd; k++)
            {
               if (bonds_ff.bmom[k] != 0.0)
               {
                  k1 = bonds_ff.idpl[k][0];
                  k2 = bonds_ff.idpl[k][1];
                  if ( (atom[i].use || atom[k1].use || atom[k2].use) && skip[i] != i)
                  {
                        xk = atom[k2].x - atom[k1].x;
                        yk = atom[k2].y - atom[k1].y;
                        zk = atom[k2].z - atom[k1].z;
                        xr = atom[i].x - atom[k1].x - xk*0.5;
                        yr = atom[i].y - atom[k1].y - yk*0.5;
                        zr = atom[i].z - atom[k1].z - zk*0.5;
                        r2 = xr*xr + yr*yr + zr*zr;
                        if (r2 <= cutoff)
                        {
                            fik = fi*bonds_ff.bmom[k];
                            rk2 = xk*xk + yk*yk + zk*zk;
                            rkr3 = sqrt(rk2*r2)*r2;
                            dotk = xk*xr + yk*yr + zk*zr;
                            e = fik*dotk/rkr3;
                            energies.eu += e;
                            atom[i].energy += e;
                            atom[k1].energy += e;
                            atom[k2].energy += e;
                            if (minim_values.iprint)
                           fprintf(pcmoutfile,"Q-DD: %2s(%-3d) -- %2s(%-3d)- %2s(%-3d) :  %-8.3f  %-8.3f   %-8.3f  = %8.4f\n",
                            atom[i].name, i,atom[k1].name,k1,atom[k2].name,k2,atom[i].charge,bonds_ff.bmom[k], rkr3, e);
                        }
                  }
               }
            }
        }
    }
}
// ===============================================================
void echrg_dipole1()
{
    int i,j,k,k1,k2;
    double cc, e,rk2,rkr3,dotk;
    double debye, f,fi,fik,cutoff;
    double xk,yk,zk;
    double xr,yr,zr,r2;
    double term,term2,term3;
    double termx,termy,termz;
    double termxk,termyk,termzk;
    double dedxi1,dedyi1,dedzi1;
    double dedxk1,dedyk1,dedzk1;
    double dedxk2,dedyk2,dedzk2;
    
    cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
    cc = 4.80298*4.80298*14.39418;
    debye = 4.8033324;
    f = -cc/(debye*units.dielec);
    for (i=1; i <= natom; i++)
       skip[i] = 0;
      
    for (i=1; i <= natom; i++)
    {
        if (atom[i].charge != 0.0)
        {
            skip[i] = i;
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                   skip[atom[i].iat[j]] = i;
            }
            fi = f*atom[i].charge;
            for (k=0; k < bonds_ff.nbnd; k++)
            {
               if (bonds_ff.bmom[k] != 0.0)
               {
                  k1 = bonds_ff.idpl[k][0];
                  k2 = bonds_ff.idpl[k][1];
                  if ( (atom[i].use || atom[k1].use || atom[k2].use) && skip[i] != i)
                  {
                        xk = atom[k2].x - atom[k1].x;
                        yk = atom[k2].y - atom[k1].y;
                        zk = atom[k2].z - atom[k1].z;
                        xr = atom[i].x - atom[k1].x - xk*0.5;
                        yr = atom[i].y - atom[k1].y - yk*0.5;
                        zr = atom[i].z - atom[k1].z - zk*0.5;
                        r2 = xr*xr + yr*yr + zr*zr;
                        if (r2 <= cutoff)
                        {
                            fik = fi*bonds_ff.bmom[k];
                            rk2 = xk*xk + yk*yk + zk*zk;
                            rkr3 = sqrt(rk2*r2)*r2;
                            dotk = xk*xr + yk*yr + zk*zr;
                            e = fik*dotk/rkr3;
                            term = fik/rkr3;
                            term2 = -3.0 * dotk / r2;
                            term3 = dotk / rk2;
                            termx = term * (xk+xr*term2);
                            termy = term * (yk+yr*term2);
                            termz = term * (zk+zr*term2);
                            termxk = term * (xr-xk*term3);
                            termyk = term * (yr-yk*term3);
                            termzk = term * (zr-zk*term3);
                            dedxi1 = termx;
                            dedyi1 = termy;
                            dedzi1 = termz;
                            dedxk1 = -0.5*termx - termxk;
                            dedyk1 = -0.5*termy - termyk;
                            dedzk1 = -0.5*termz - termzk;
                            dedxk2 = -0.5*termx + termxk;
                            dedyk2 = -0.5*termy + termyk;
                            dedzk2 = -0.5*termz + termzk;
                            
                            energies.eu += e;
                            deriv.deqq[i][0] += dedxi1;
                            deriv.deqq[i][1] += dedyi1;
                            deriv.deqq[i][2] += dedzi1;
                            deriv.deqq[k1][0] += dedxk1;
                            deriv.deqq[k1][1] += dedyk1;
                            deriv.deqq[k1][2] += dedzk1;
                            deriv.deqq[k2][0] += dedxk2;
                            deriv.deqq[k2][1] += dedyk2;
                            deriv.deqq[k2][2] += dedzk2;
                        }
                  }
               }
            }
        }
    }
}
                                                     
