#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "derivs.h"
#include "utility.h"
#include "pot.h"
#include "gmmx.h"
#include "minim_values.h"
#include "minim_control.h"
#include "minvar.h"
#include "minimize.h"

// mode
#define NONE      0
#define NEWTON    1
#define TNCG      2
#define DTNCG     3
#define AUTO      4
// method
#define DIAG      5
#define SSOR      6
#define ICCG      7
#define BLOCK     8


#define  Success     0
#define  ReSearch    1
#define  WideAngle   2
#define  BadIntpln   3
#define  IntplnErr   4
#define  blank       5
#define  Failure     6
#define  SmallFct    7
#define  SmallGrad   8


         
double scale2;
static int status;

void minimize()
{
   int i,nvar, iter, icount;
   double minimum,grdmin;
   double minimiz1(), newton1();
   double *xx;  //  xx[maxvar];
   int method, maxvar;
   
   maxvar = 3*natom;
   xx = dvector(0,maxvar);
      
   grdmin = 1.0;
   scale2 = 12.0; // bfgs  12.0

   if (minim_control.method == 1 || minim_control.method == 3 || minim_control.method == 4)
   {
     nvar = 0;
     icount = 0;
     for (i=1; i <= natom; i++)
     {
       if (atom[i].use)
       {
         xx[nvar] = atom[i].x*scale2;
         nvar++;
         xx[nvar] = atom[i].y*scale2;
         nvar++;
         xx[nvar] = atom[i].z*scale2;
         nvar++;
       }
     }

     method = 1;
     grdmin = 0.5;
     if (minim_control.method == 1)
        grdmin = 0.1;
     else
        grdmin = 0.5;
     if (pot.use_picalc)
     {
        while (TRUE)
        {
         iter = 0;
         piseq(0,0);
         if (pot.use_charge) piden();
         mqn(nvar,method, &iter,xx, &minimum,grdmin, minimiz1 );
         icount++;
         if (iter < 10 && icount > 1)
           break;
                 if (icount > 10)
                 {
                        status = Failure;
                        break;
                 }
        }
     } else
        mqn(nvar,method, &iter,xx, &minimum,grdmin, minimiz1 );
   }
   
   if (status != Failure)
   {
      if (minim_control.method == 2 || minim_control.method == 3 || minim_control.method == 4)
      {
        scale2 = 1.0;   // tcng
        grdmin = 0.0001;
        nvar = 0;
        icount = 0;
        for (i=1; i <= natom; i++)
        {
          if (atom[i].use)
          {
            xx[nvar] = atom[i].x*scale2;
            nvar++;
            xx[nvar] = atom[i].y*scale2;
            nvar++;
            xx[nvar] = atom[i].z*scale2;
            nvar++;
          }
        }
        if (pot.use_picalc)
        {
           while(TRUE)
           {
             iter=0;
             piseq(0,0);
             if (pot.use_charge) piden();
             tncg(nvar,method,&iter, xx, &minimum, grdmin,newton1, newton2);
             icount++;
             if(iter < 10 && icount > 1)
               break;
                        if (icount > 10)
                        {
                                status = Failure;
                                break;
                        }
           }
        }else   
           tncg(nvar,method,&iter, xx, &minimum, grdmin,newton1, newton2); 
   
        nvar = 0;
        for (i=1; i <= natom; i++)
        {
          if (atom[i].use)
          {
            atom[i].x =  xx[nvar]/scale2;
            nvar++;
            atom[i].y =  xx[nvar]/scale2;
            nvar++;
            atom[i].z =  xx[nvar]/scale2;
            nvar++;
          }
        }
      }
   }
   free_dvector(xx, 0, maxvar);
}
// =============================
void mqn(int nvar,int method,int *iter, double *x,  double *minimum, double grdmin, double (*fgvalue) ())
{
      int i,ncalls,nerror;
      int niter,period,nstart;
      double fast,slow,epsln,d1temp,d2temp;
      double f,f_old,f_new,f_move;
      double rms,beta,x_move,g_norm,g_rms;
      double gg,gg_old;
      double sg,dg,sd,dd,angle;
      double *g, *p;                         //   g[maxvar];
      double *x_old, *g_old;             //   x_old[maxvar],g_old[maxvar];
      double  *s, *d;                 //   p[maxvar],s[maxvar],d[maxvar];
      double fctmin;
      int restart, terminate;
      int maxiter, nextiter, maxvar;

      maxvar = 3*natom;
      x_old = dvector(0,maxvar);
      g_old = dvector(0,maxvar);
      g = dvector(0,maxvar);
      p = dvector(0,maxvar);
      s = dvector(0,maxvar);
      d = dvector(0,maxvar);
      
      ncalls = 0;
      rms = sqrt((float)nvar)/ sqrt(3.0);
      restart = TRUE;
      terminate = FALSE;
      status = blank;
      nerror = 0;

      fctmin = -10000.0;
      maxiter = 1000;
      nextiter = 1;
      fast = 0.5;
      slow = 0.0;
      epsln = 1.0e-16;
      if (nvar > 200)
         period = nvar;
      else
         period = 200;
      minvar.cappa = .1;
      minvar.stpmin = 1.0e-20;
      minvar.stpmax = 5.0;
      minvar.angmax = 100.0;
      minvar.intmax = 5;

      niter = nextiter -1;
      ncalls = ncalls + 1;
      f = fgvalue(x, g); // get function and first deriv at original point
      g_norm = 0.0;
      for (i=0; i < nvar; i++)
      {
          x_old[i] = x[i];
          g_old[i] = g[i];
          g_norm += g[i]*g[i];
      }
      g_norm = sqrt(g_norm);
      f_move = 0.5*minvar.stpmax*g_norm;
      g_rms = g_norm*scale2/rms;

     if (niter > maxiter)
        terminate = TRUE;
     if (f < fctmin)
        terminate = TRUE;
     if (g_rms < grdmin)
         terminate = TRUE;

     while ( ! terminate)
     {
         niter++;
         status = Success;   
         if (niter > maxiter)
            terminate = TRUE;


         if (restart || method == 0)
         {
             for (i=0; i < nvar; i++)
                 p[i] = -g[i];
             nstart = niter;
             restart = FALSE;
         } else if (method == 1) // BFGS method
         {
             sg = 0.0;
             dg = 0.0;
             dd = 0.0;
             sd = 0.0;
             for (i=0; i < nvar; i++)
             {
               sg += s[i]*g[i];
               dg += d[i]*g[i];
               dd += d[i]*d[i];
               sd += s[i]*d[i];
             }
             for (i=0; i < nvar; i++)
             {
                 d1temp = (d[i]*sg + s[i]*dg)/sd;
                 d2temp = (1.0+dd/sd)*(s[i]*sg/sd);
                 p[i] = -g[i] + d1temp - d2temp;
             }
         } else if (method == 2)  // Fletcher Reeves
         {
             gg = 0.0;
             gg_old = 0.0;
             for (i=0; i < nvar; i++)
             {
                 gg += g[i]*g[i];
                 gg_old += g_old[i]*g_old[i];
             }
             beta = gg/gg_old;
             for (i=0; i < nvar; i++)
                p[i] = -g[i] + beta*p[i];
         } else if (method == 3)  // Polak Ribere
         {
             dg = 0.0;
             gg_old = 0.0;
             for (i=0; i < nvar; i++)
             {
                 dg += d[i]*g[i];
                 gg_old += g_old[i]*g_old[i];
             }
             beta = dg/gg_old;
             for (i=0; i < nvar; i++)
                p[i] = -g[i] + beta*p[i];
        }

         // do a line search
         f_old = f;
         search(nvar,&f,g,x,p,f_move,&angle,&ncalls,fgvalue,&status);
         f_new = f;
         
         f_move = f_old - f_new;
         x_move = 0.0;
         g_norm = 0.0;
         for (i=0; i < nvar; i++)
         {
            s[i] = x[i] - x_old[i];
            d[i] = g[i] - g_old[i];
            x_move += s[i]*s[i];
            g_norm += g[i]*g[i];
            x_old[i] = x[i];
            g_old[i] = g[i];
         }
         x_move = sqrt(x_move) / (scale2 * rms);
         g_norm = sqrt(g_norm);
         g_rms = g_norm * scale2/rms;
         
// function increase
         if (f_move <= 0.0)
         {
            status = Failure;
            nerror = nerror + 1;
            if (nerror == 3)
               terminate = TRUE;
            else
               restart = TRUE;

            for(i=0; i < nvar; i++)
            {
               x[i] = x_old[i];
               g[i] = g_old[i];
            }
         }
         if (x_move < epsln)
         {
             nerror++;
             if (nerror > 3)
               terminate = TRUE;
             else
               restart = TRUE;
         }
// normal termination
         if (f < fctmin)
         {
            status = Success;
            terminate = TRUE;
         }
         if (g_rms < grdmin)
         {
            status = Success;
            terminate = TRUE;
         }

         if (terminate || (niter%10 == 0) )
         {
            
            if (fabs(f_move) < 0.00001)
            { 
                terminate = TRUE;
                status = Failure;
            }
            if (fabs(x_move) < 0.00001) 
            {   
                terminate = TRUE;
                status = Failure;
            }            

            if ( !gmmx.run )
             fprintf(pcmoutfile," %d   %f    %f   %f   %f  \n",niter,
             f,g_rms,f_move,x_move);
             
         }
     }
     *minimum = f;
     *iter = niter;
     free_dvector(x_old ,0,maxvar);
     free_dvector(g_old,0,maxvar);
     free_dvector( g ,0,maxvar);
     free_dvector( p ,0,maxvar);
     free_dvector( s ,0,maxvar);
     free_dvector( d ,0,maxvar);             
}

double minimiz1(double *xx, double *g)
{
    int i,nvar;
    double e_min;

    nvar = 0;
    for (i = 1; i <= natom; i++)
    {
        if (atom[i].use)
        {
            atom[i].x = xx[nvar]/scale2;
            nvar++;
            atom[i].y = xx[nvar]/scale2;
            nvar++;
            atom[i].z = xx[nvar]/scale2;
            nvar++;
        }
    }

    gradient();
    e_min = energies.total;

    nvar = 0;
    for (i=1; i <= natom; i++)
    {
        if (atom[i].use)
        {
            xx[nvar] = atom[i].x*scale2;
            g[nvar] = deriv.d1[i][0]/scale2;
            nvar++;
            xx[nvar] = atom[i].y*scale2;
            g[nvar] = deriv.d1[i][1]/scale2;
            nvar++;
            xx[nvar] = atom[i].z*scale2;
            g[nvar] = deriv.d1[i][2]/scale2;
            nvar++;
        }
    }
    return(e_min);
}

double newton1(double *xx, double *g)
{
    int i, nvar;
    double e;

    nvar = 0;
    for (i=1; i <= natom; i++)
    {
        if (atom[i].use)
        {
            atom[i].x = xx[nvar];
            nvar++;
            atom[i].y = xx[nvar];
            nvar++;
            atom[i].z = xx[nvar];
            nvar++;
        }
    }

    gradient();
    e = energies.total;

    nvar = 0;
    for (i=1; i <= natom; i++)
    {
        if (atom[i].use)
        {
            xx[nvar] = atom[i].x;
            g[nvar] = deriv.d1[i][0];
            nvar++;
            xx[nvar] = atom[i].y;
            g[nvar] = deriv.d1[i][1];
            nvar++;
            xx[nvar] = atom[i].z;
            g[nvar] = deriv.d1[i][2];
            nvar++;
        }
    }
    return(e);
}

void newton2(int mode, double *xx,double *h,int *hinit,
               int *hstop, int *hindex, double *hdiag)
{
    int i,j,k,nvar, maxvar, nuse;
    int *hvar, *huse;   //  hvar[maxvar],huse[maxvar];

    if (mode == NONE)
       return;

    maxvar = 3*natom;
    hvar = ivector(0,maxvar);
    huse = ivector(0,maxvar);
    
    nvar = 0;
    nuse = TRUE;
    for (i=1; i <= natom; i++)
    {
        if (atom[i].use)
        {
           atom[i].x = xx[nvar];
           nvar++;
           atom[i].y = xx[nvar];
           nvar++;
           atom[i].z = xx[nvar];
           nvar++;
        } else
           nuse = FALSE;
    }

    hessian(h,hinit,hstop,hindex,hdiag);

    nvar = 0;
    if (nuse == FALSE)
    {
       for (i=1; i <= natom; i++)
       {
           k = 3*(i-1);
           if (atom[i].use)
           {
              for (j=0; j < 3; j++)
              {
                hvar[nvar] = j+k;
                huse[j+k] = nvar;
                nvar++;
              }
           } else
           {
              for (j=0; j < 3; j++)
                 huse[j+k] = 0;
           }
       }
       for (i=0; i < nvar; i++)
       {
          k = hvar[i];
          hinit[i] = hinit[k];
          hstop[i] = hstop[k];
          hdiag[i] = hdiag[k];
          for (j=hinit[i]; j < hstop[i]; j++)
            hindex[j] = huse[hindex[j]];
       }
    }
//
    nvar = 0;
    for (i=1; i <= natom; i++)
    {
        if (atom[i].use)
        {
            xx[nvar] = atom[i].x;
            nvar++;
            xx[nvar] = atom[i].y;
            nvar++;
            xx[nvar] = atom[i].z;
            nvar++;
        }
    }         
    free_ivector(hvar ,0,maxvar);
    free_ivector(huse ,0,maxvar);

}
            
            
            
