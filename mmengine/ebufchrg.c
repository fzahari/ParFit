#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "angles.h"
#include "derivs.h"
#include "hess.h"
#include "utility.h"
#include "field.h"
#include "attached.h"
#include "cutoffs.h"
#include "pot.h"
#include "units.h"
#include "minim_values.h"
#include "skip.h"
#include "ebufcharge.h"
        
void ebufcharge()
{
   int i, j, k,temp;
   double xi, yi, zi, xk,yk,zk;
   double xr, yr, zr;
   double e, rik, cc, rik2;
   double rdielc;
   double rdn;
   double cutoff;

   cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
   
   energies.eu = 0.0;
   cc = 332.0716;
   rdielc = 1.0/units.dielec;

   rdn = 1.0;
   if (field.type == MMX || field.type == MM2 )
     rdn = .915;
   else if (field.type == MM3)
     rdn = .923;

   if (minim_values.iprint)
   {
       fprintf(pcmoutfile,"\nCharge-Charge Interactions\n");
       fprintf(pcmoutfile,"      At1      At2      q1        q2         Rik          Eqq\n");
   }
   for (i=1; i <= natom; i++)
      skip[i] = 0;

   for (i=1; i < natom; i++)
   {
      if (atom[i].charge != 0.0 )
      {
         for (k=0; k < MAXIAT; k++)
         {
            if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
              skip[atom[i].iat[k]] = i;
         }
         for(k=0; k < attached.n13[i]; k++)
            skip[attached.i13[k][i]] = i;
         for(k=0; k < attached.n14[i]; k++)
         {
            temp = attached.i14[k][i];
            skip[attached.i14[k][i]] = -i;
         }
         
         for(j=i+1; j <= natom; j++)
         {
           if (atom[i].use || atom[j].use)
           {
             if ( skip[j] != i)
             {
               if (atom[j].charge != 0.0)
               { 
                  if (atom[i].type == 5 || atom[i].type == 36)
                  {
                   // compute reduced distance
                    xi = rdn*(atom[i].x-atom[atom[i].iat[0]].x) + atom[atom[i].iat[0]].x;
                    yi = rdn*(atom[i].y-atom[atom[i].iat[0]].y) + atom[atom[i].iat[0]].y;
                    zi = rdn*(atom[i].z-atom[atom[i].iat[0]].z) + atom[atom[i].iat[0]].z;
                  } else
                  {
                    xi = atom[i].x;
                    yi = atom[i].y;
                    zi = atom[i].z;
                  }
                  if (atom[j].type == 5 || atom[j].type == 36)
                  {
                   // compute reduced distance
                     xk = rdn*(atom[j].x-atom[atom[j].iat[0]].x) + atom[atom[j].iat[0]].x;
                     yk = rdn*(atom[j].y-atom[atom[j].iat[0]].y) + atom[atom[j].iat[0]].y;
                     zk = rdn*(atom[j].z-atom[atom[j].iat[0]].z) + atom[atom[j].iat[0]].z;
                  } else
                  {
                     xk = atom[j].x;
                     yk = atom[j].y;
                     zk = atom[j].z;
                  }
                  xr = xi - xk;
                  yr = yi - yk;
                  zr = zi - zk;
                  rik2 =  xr*xr + yr*yr + zr*zr;
                  if (rik2 < cutoff)
                  {
                     rik = sqrt( xr*xr + yr*yr + zr*zr);
                     e = cc*rdielc*atom[i].charge*atom[j].charge/(rik+0.05);
                     if (skip[j] == -i)
                       e *= units.chgscale;  
                     energies.eu += e;
                     atom[i].energy += e;
                     atom[j].energy += e;
                     if (minim_values.iprint)
                      fprintf(pcmoutfile,"QQ: %2s(%-3d)- %2s(%-3d)   %-8.3f  %-8.3f   %-8.3f = %8.4f\n",
                       atom[i].name, i, atom[j].name, j,atom[i].charge, atom[j].charge,rik,e);
                  }
                }
             }
           }
         }
      }
   }

}

void ebufcharge1()
{
   int i,j, k;
   double xi, yi, zi;
   double xr, yr, zr;
   double xk,yk,zk,rdn;
   double e, rik, cc, rik2, st;
   double rdielc;
   double cutoff;

   cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
   
   energies.eu = 0.0;
   cc = 332.0716;
   rdielc = 1.0/units.dielec;
   rdn = 1.0;
   if (field.type == MMX || field.type == MM2 )
     rdn = .915;
   else if (field.type == MM3)
     rdn = .923;

   for (i=0; i <= natom; i++)
   {
       deriv.deqq[i][0] = 0.0;
       deriv.deqq[i][1] = 0.0;
       deriv.deqq[i][2] = 0.0;
   }
      
   for (i=1; i <= natom; i++)
      skip[i] = 0;

   for (i=1; i < natom; i++)
   {
       if (atom[i].charge != 0.0 )
       {
         for (k=0; k < MAXIAT; k++)
         {
            if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
              skip[atom[i].iat[k]] = i;
         }
         for(k=0; k < attached.n13[i]; k++)
            skip[attached.i13[k][i]] = i;
         for(k=0; k < attached.n14[i]; k++)
            skip[attached.i14[k][i]] = -i;

           for(j=i+1; j <= natom; j++)
           {
             if (atom[i].use || atom[j].use)
             {
               if (atom[j].charge != 0.0 &&  skip[j] != i  )
               {
                if (atom[i].type == 5 || atom[i].type == 36)
                {
                   // compute reduced distance
                   xi = rdn*(atom[i].x-atom[atom[i].iat[0]].x) + atom[atom[i].iat[0]].x;
                   yi = rdn*(atom[i].y-atom[atom[i].iat[0]].y) + atom[atom[i].iat[0]].y;
                   zi = rdn*(atom[i].z-atom[atom[i].iat[0]].z) + atom[atom[i].iat[0]].z;
               } else
               {
                  xi = atom[i].x;
                  yi = atom[i].y;
                  zi = atom[i].z;
               }
                if (atom[j].type == 5 || atom[j].type == 36)
                {
                   // compute reduced distance
                   xk = rdn*(atom[j].x-atom[atom[j].iat[0]].x) + atom[atom[j].iat[0]].x;
                   yk = rdn*(atom[j].y-atom[atom[j].iat[0]].y) + atom[atom[j].iat[0]].y;
                   zk = rdn*(atom[j].z-atom[atom[j].iat[0]].z) + atom[atom[j].iat[0]].z;
               } else
               {
                  xk = atom[j].x;
                  yk = atom[j].y;
                  zk = atom[j].z;
               }
                xr = xi - xk;
                yr = yi - yk;
                zr = zi - zk;
                rik2 = xr*xr +yr*yr + zr*zr;
                if (rik2 < cutoff)
                {
                   rik = sqrt( rik2);
                   st = -cc*rdielc*atom[i].charge*atom[j].charge/(rik*(rik+0.05)*(rik+0.05));
                   if (skip[j] == -i)
                     st *= units.chgscale;
                   e = cc*rdielc*atom[i].charge*atom[j].charge/(rik+0.05);
                   if (skip[j] == -i)
                     e *= units.chgscale;
                  
                   deriv.deqq[i][0] += st*xr;
                   deriv.deqq[i][1] += st*yr;
                   deriv.deqq[i][2] += st*zr;
         
                   deriv.deqq[j][0] -= st*xr;
                   deriv.deqq[j][1] -= st*yr;
                   deriv.deqq[j][2] -= st*zr;
         
                   virial.virx += xr*st*xr;
                   virial.viry += yr*st*yr;
                   virial.virz += zr*st*zr;
                   
                   energies.eu += e;
                }
               }
             }
           }
       }
   }
}

void ebufcharge2(int i)
{
    int j,k, temp;
    double fik,de,d2e,d2edx,d2edy,d2edz;
    double xi,yi,zi,xr,yr,zr,term[3][3];
    double xk,yk,zk,rdn;
    double r,r2;
    double cc, rdielc;
    double cutoff;
    
    if (fabs(atom[i].charge) < 0.01)
       return;

    cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
    
    cc = 332.0716;
    rdielc = 1.0/units.dielec;
    rdn = 1.0;
    if (field.type == MMX || field.type == MM2 )
     rdn = .915;
    else if (field.type == MM3)
     rdn = .923;
    
    for (j=1; j <= natom; j++)
      skip[j] = 0;
      
    skip[i] = i;
    for (k=0; k < MAXIAT; k++)
    {
        if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
          skip[atom[i].iat[k]] = i;
    }
    for(k=0; k < attached.n13[i]; k++)
       skip[attached.i13[k][i]] = i;
    for(k=0; k < attached.n14[i]; k++)
    {
       temp = attached.i14[k][i];
       skip[attached.i14[k][i]] = -i;
    }
    
    if (atom[i].type == 5 || atom[i].type == 36)
    {
          // compute reduced distance
          xi = rdn*(atom[i].x-atom[atom[i].iat[0]].x) + atom[atom[i].iat[0]].x;
          yi = rdn*(atom[i].y-atom[atom[i].iat[0]].y) + atom[atom[i].iat[0]].y;
          zi = rdn*(atom[i].z-atom[atom[i].iat[0]].z) + atom[atom[i].iat[0]].z;
    } else
    {
         xi = atom[i].x;
         yi = atom[i].y;
         zi = atom[i].z;
     }
    
    for (k=1; k <= natom; k++)
    {
        if (atom[k].charge != 0.0  &&  (skip[k] != i) )
        {
                if (atom[j].type == 5 || atom[j].type == 36)
                {
                   // compute reduced distance
                   xk = rdn*(atom[k].x-atom[atom[k].iat[0]].x) + atom[atom[k].iat[0]].x;
                   yk = rdn*(atom[k].y-atom[atom[k].iat[0]].y) + atom[atom[k].iat[0]].y;
                   zk = rdn*(atom[k].z-atom[atom[k].iat[0]].z) + atom[atom[k].iat[0]].z;
               } else
               {
                  xk = atom[k].x;
                  yk = atom[k].y;
                  zk = atom[k].z;
               }
               xr = xi - xk;
               yr = yi - yk;
               zr = zi - zk;
              r2 = xr*xr + yr*yr + zr*zr;
              if (r2 < cutoff)
             {
              r = sqrt(r2);
              fik = cc*rdielc*atom[i].charge*atom[k].charge;
              if (skip[k] == -i)
                 fik *= units.chgscale;
              de = -fik/((r+0.05)*(r+0.05));
              d2e = -2.0*de/r;
              de = de / r;
              d2e = (d2e-de) / r2;
              d2edx = d2e * xr;
              d2edy = d2e * yr;
              d2edz = d2e * zr;
              term[0][0]= d2edx*xr + de;
              term[1][0] = d2edx*yr;
              term[2][0] = d2edx*zr;
              term[0][1] = term[1][0];
              term[1][1] = d2edy*yr + de;
              term[2][1] = d2edy*zr;
              term[0][2] = term[2][0];
              term[1][2] = term[2][1];
              term[2][2] = d2edz*zr + de;

              for (j=0; j < 3; j++)
              {
                  hess.hessx[i][j] += term[j][0];
                  hess.hessy[i][j] += term[j][1];
                  hess.hessz[i][j] += term[j][2];
                  hess.hessx[k][j] -= term[j][0];
                  hess.hessy[k][j] -= term[j][1];
                  hess.hessz[k][j] -= term[j][2];
              }
             }
        }
    }
}
