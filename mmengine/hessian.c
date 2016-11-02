#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "derivs.h"
#include "hess.h"
#include "utility.h"
#include "pot.h"

       
void hessian(double *h, int *hinit, int *hstop,int *hindex, double *hdiag)
{
   int i,j,k,nhess;
   double hmax;
   char hatext[60];
   
   int *keep, maxhess;  //  keep[MAXATOM];

   keep = ivector(0,natom);

   if (natom < 300)
      maxhess = (3*natom*(3*natom-1))/2;
   else if (natom < 800)
      maxhess = (3*natom*(3*natom-1))/4;
   else
      maxhess = (3*natom*(3*natom-1))/20;
   

   nhess = 0;
   for (i=1; i <= natom; i++)
   {

       hinit[nhess] = 0;
       hstop[nhess] = 0;
       hdiag[nhess] = 0.0;
       nhess++;
       hinit[nhess] = 0;
       hstop[nhess] = 0;
       hdiag[nhess] = 0.0;
       nhess++;
       hinit[nhess] = 0;
       hstop[nhess] = 0;
       hdiag[nhess] = 0.0;
       nhess++;       
   }
/*   periodic boundary conditions set here  */

 //   if (pot.use_picalc) piseq(0,0);

    nhess= 0;
    for (i=1; i <= natom; i++)
    {
       for(k=1; k <= natom; k++)
       {
           for(j=0; j < 3; j++)
           {
              hess.hessx[k][j] = 0.0;
              hess.hessy[k][j] = 0.0;
              hess.hessz[k][j] = 0.0;
           }
       }
          
      if (atom[i].use)
      {
        if (pot.use_bond)ebond2(i);
        if (pot.use_angle)eangle2(i); 
        if (pot.use_strbnd)estrbnd2(i);
        if (pot.use_angang)eangang2(i);
        if (pot.use_opbend)eopbend2(i);
        if (pot.use_opbend_wilson) eopbend_wilson2(i);
        if (pot.use_tors)etorsion2(i); 
        if (pot.use_strtor) estrtor2(i);
        if (pot.use_coordb)ecoord_bond2(i);
        
        if (pot.use_charge)echarge2(i);
        if (pot.use_hal)ehal2(i);  
        if (pot.use_bufcharge)ebufcharge2(i);
        if (pot.use_dipole)edipole2(i);
        if (pot.use_solv)esolv2(i);
        if (pot.use_highcoord) ehigh_coord2(i);

        hdiag[(i-1)*3]   += hess.hessx[i][0];
        hdiag[(i-1)*3+1] += hess.hessy[i][1];
        hdiag[(i-1)*3+2] += hess.hessz[i][2];

        for (k=i+1; k <= natom; k++)
        {
          keep[k] = FALSE;
          if (atom[k].use)
          {
              hmax = fabs(hess.hessx[k][0]);
              for (j=0; j < 3;j++)
              {
                  if (fabs(hess.hessx[k][j]) > hmax)
                     hmax = fabs(hess.hessx[k][j]);
              }
              for (j=0; j < 3;j++)
              {
                  if (fabs(hess.hessy[k][j]) > hmax)
                     hmax = fabs(hess.hessy[k][j]);
              }
              for (j=0; j < 3;j++)
              {
                  if (fabs(hess.hessz[k][j]) > hmax)
                     hmax = fabs(hess.hessz[k][j]);
              }
              if (hmax >= hesscut) keep[k] = TRUE;
          }
        }
      
        hinit[(i-1)*3] = nhess;  //  x component of new atom
      
        hindex[nhess] = 3*(i-1) + 1;      
        h[nhess] = hess.hessx[i][1];
        nhess++;
        hindex[nhess] = 3*(i-1) + 2;       
        h[nhess] = hess.hessx[i][2];
        nhess++;
      
        for (k=i+1; k <= natom; k++)
        {
          if (keep[k])
          {
              for (j=0; j < 3; j++)
              {
                  hindex[nhess] = 3*(k-1) +j;
                  h[nhess] = hess.hessx[k][j];
                  nhess++;
              }
          }
        }
        hstop[(i-1)*3] = nhess;
      
        hinit[(i-1)*3+1] = nhess;  // y component of atom i
        hindex[nhess] = 3*(i-1) + 2;
        h[nhess] = hess.hessy[i][2];
        nhess++;
        for (k=i+1; k <= natom; k++)
        {
          if (keep[k])
          {
              for (j=0; j < 3; j++)
              {
                  hindex[nhess] = 3*(k-1) +j;
                  h[nhess] = hess.hessy[k][j];
                  nhess++;
              }
          }
        }
        hstop[(i-1)*3+1] = nhess;
        hinit[(i-1)*3+2] = nhess;
        for (k=i+1; k <= natom; k++)
        {
          if (keep[k])
          {
              for (j=0; j < 3; j++)
              {
                  hindex[nhess] = 3*(k-1) +j;
                  h[nhess] = hess.hessz[k][j];
                  nhess++;
              }
          }
        }
        hstop[(i-1)*3+2] = nhess;

        if (nhess > maxhess)
        {
          sprintf(hatext,"Too many hessian elements: %d elements  %d max\n",nhess,maxhess);
          message_alert(hatext,"Error");
          printf("Too many hessian elements !\n");
          exit(0);
        }
      }
    }
    free_ivector(keep, 0,natom);
}



