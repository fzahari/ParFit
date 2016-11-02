#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"
#include "bonds_ff.h"
#include "torsions.h"
#include "pimatx.h"
#include "pistuf.h"
#include "pisave.h"
#include "minim_values.h"
#include "minim_control.h"
#include "pibonds.h"
#include "pidata.h"
#include "k.h"
        
void piread(int iset)
{
    int i,k,ks;
    
    ks = 0;

    huckel();
    
    for (i=0; i < pistuf.norb; i++)
    {
        for (k=0; k < pistuf.norb; k++)
            pisave.fjsave[i][k] = pimatx.f[i][k];
    }
    if (iset != 1)
    {
        if (pistuf.nplane == 1 || ks == 0)
        {
            ks = 1;
            for (i=0; i < pistuf.norb; i++)
            {
                for(k=0; k < pistuf.norb; k++)
                    pisave.fksave[i][k] = pimatx.f[i][k];
            }
        }
    }  
}

void huckel()
{
    int i,ii,ik,il,kl,no;
    int ia, ib;
    long int mask;

    no = pistuf.norb-1;
    mask = (1L << 0);
    kl = il = 0;

    for (i=0; i < bonds_ff.nbnd; i++)
    {
        ia = bonds_ff.i12[i][0];
        ib = bonds_ff.i12[i][1];
        if (atom[ia].flags & mask && atom[ib].flags & mask)
        {
            for (ii = 0 ; ii < pistuf.norb; ii++)
            {
                if (pistuf.listpi[ii] == ia)
                   il = ii;
                if (pistuf.listpi[ii] == ib)
                   kl = ii;
            } 
            ii = atom[ia].mmx_type;
            ik = atom[ib].mmx_type;
            pimatx.f[il][kl] = -1.0;
            
            if (ii == 6 || ii == 9 || ii == 42 || ii == 46)
            {
                pimatx.f[il][kl] = -0.8;
            } else if (ii == 15)
            {
                pimatx.f[il][kl] = -0.4;               
            } else if ( ii == 7)
            {
                pimatx.f[il][kl] = -1.0;
            } else if ( ii == 10 || ii == 41)
            {
                if (ik == 6 || ik == 8 || ik == 9)
                    pimatx.f[il][kl] = -0.8;
                else
                    pimatx.f[il][kl] = -1.0;
            } else
            {
                if (ik == 6 || ik == 9 || ik == 46)
                    pimatx.f[il][kl] = -0.8;
                else if (ik == 15)
                    pimatx.f[il][kl] = -0.4;
                else
                    pimatx.f[il][kl] = -1.0;
            }
            pimatx.f[kl][il] = pimatx.f[il][kl];
        }
    }
    
    for (i=0; i < pistuf.norb; i++)
    {
        ia = pistuf.listpi[i];
        ii = atom[ia].mmx_type;
        pimatx.f[i][i] = 0.0;
        if (ii == 6 || ii == 41 || ii == 46)
           pimatx.f[i][i] = 2.0;
        else if (ii == 7 || ii == 42 || ii == 66)
           pimatx.f[i][i] = 1.0;
        else if (ii == 8)
           pimatx.f[i][i] = 1.5;
        else if (ii == 10)
           pimatx.f[i][i] = 0.5;
        else if (ii == 15)
           pimatx.f[i][i] = 0.5;
    }
}
// ================================
void dircos()
{
    int i,ii,ik,k,l,list[8],no;
    float a1, a2, b1, b2, c1, c2, rik, rnorm,
           s, x2, x3, xr[8], y2, y3, yr[8], z2, z3, zr[8];

    no = pistuf.norb -1;
    pibonds.npisave = 0;
    for (i=0; i < no; i++)
    {
        list[0] = pistuf.listpi[i];
        for (k = i+1; k < pistuf.norb; k++)
        {
            list[1] = pistuf.listpi[k];
            for (l=0; l < 3; l++)
            {
                list[l+2] = pistuf.piplane[i][l];
                list[l+5] = pistuf.piplane[k][l];
            }
            rotat(list, xr, yr, zr);
            ii = atom[pistuf.listpi[i]].mmx_type;
            ik = atom[pistuf.listpi[k]].mmx_type;
            if ( ((ii==4 || ii==10) && (ik==4 || ik==10)) ||
                 ((ii==4 || ii==41) && (ik==4 || ik==10)) )
            {
                a1 = a2 = 0.0;
                b1 = b2 = 1.0;
                c1 = c2 = 0.0;
            } else
            {
               x2 = xr[3] - xr[2];
               y2 = yr[3] - yr[2];
               z2 = zr[3] - zr[2];
               x3 = xr[4] - xr[2];
               y3 = yr[4] - yr[2];
               z3 = zr[4] - zr[2];
/*  DIRECTION NUMBERS AND COSINES OF THE FIRST ATOM.  */
               a1 = y2*z3 - y3*z2;
               b1 = x3*z2 - x2*z3;
               c1 = x2*y3 - x3*y2;
	       if ( fabs(a1) < 0.001 && fabs(b1) < 0.001 && fabs(c1) < 0.001)
		 {
                   a1 = 0.0;
                   b1 = 1.0;
                   c1 = 0.0;
		 } else
		 {
		   rnorm = sqrt( a1*a1 + b1*b1 + c1*c1 );
		   if (rnorm > 0.1)
		     {
		       a1 /= rnorm;
		       b1 /= rnorm;
		       c1 /= rnorm;
		     }else
		     {
		       a1 = 0.0;
		       b1 = 1.0;
		       c1 = 0.0;
		     }
		 }
/*  THE SECOND ATOM, BASED ON ATOM IPL(I,K), IS TRANSLATED TO THE ORIGIN. */
               x2 = xr[6] - xr[5];
               y2 = yr[6] - yr[5];
               z2 = zr[6] - zr[5];
               x3 = xr[7] - xr[5];
               y3 = yr[7] - yr[5];
               z3 = zr[7] - zr[5];
/*  DIRECTION NUMBERS AND COSINES OF THE SECOND ATOM.  */
               a2 = y2*z3 - y3*z2;
               b2 = x3*z2 - x2*z3;
               c2 = x2*y3 - x3*y2;
	       if ( (fabs(a2) < 0.001) && (fabs(b2) < 0.001) && (fabs(c2) < 0.001) )
		 {
                   a2 = 0.0;
                   b2 = 1.0;
                   c2 = 0.0;
		 } else
		 {
		   rnorm = sqrt( a2*a2 + b2*b2 + c2*c2 );
		   /*  THE A2 COSINE CHANGES SIGN TO CORRESPOND TO THE INTERNUCLEAR AXIS.  */
		   if (rnorm > 0.1)
		     {
		       a2 = -a2/rnorm;
		       b2 /= rnorm;
		       c2 /= rnorm;
		     } else
		     {
		       a2 = 0.0;
		       b2 = 1.0;
		       c2 = 0.0;
		     }
		   s = a1*a2 + b1*b2 + c1*c2;
                  if( s < 0 )
                  {
                     a2 = -a2;
                     b2 = -b2;
                     c2 = -c2;
                  }
               }
            }
            rik = fabs( xr[1] )/0.5292;
            pibonds.riks[pibonds.npisave] = rik;
            pibonds.a1s[pibonds.npisave] = a1;
            pibonds.b1s[pibonds.npisave] = b1;
            pibonds.c1s[pibonds.npisave] = c1;
            pibonds.a2s[pibonds.npisave] = a2;
            pibonds.b2s[pibonds.npisave] = b2;
            pibonds.c2s[pibonds.npisave] = c2;
            pibonds.npisave += 1;
        }
    }
}
// ==============================================
void picon()
{
    int i,ia,ib,ic,id,ii,k,kk,nc;
    int iti,itk,ikno;
    long int mask;
    float bkon, blen, factor, pbo, pbw, slpk, slpl,slpl2,
           sumcos, sumdif, tc, ang, dihed;
    double temp;

    sumdif = 0.0;
    sumcos = 0.0;
    factor = 1.0;
    mask = 1L << 0;

    for (i=0; i < pistuf.norb; i++)
    {
        for(k=0; k < pistuf.norb; k++)
          pimatx.hc[i][k] = pisave.gjsave[i][k];
    }
    for (i=0; i < pistuf.norb; i++)
    {
        for(k=0; k < pistuf.norb; k++)
          pimatx.v[i][k] = pisave.gksave[i][k];
    }


    if (pistuf.nplane != 0)
    {
        for (i=0; i < pistuf.ntpi; i++)
        {
            ii = pitorsave.ncsave[i];
            ia = torsions.i14[ii][0];
            ib = torsions.i14[ii][1];
            ic = torsions.i14[ii][2];
            id = torsions.i14[ii][3];
            dihed = dihdrl(ia,ib,ic,id);
            ang = .01745329F*fabs( dihed );
            ang = cos( ang );

            ii = pitorsave.itsave[i];
            kk = pitorsave.ktsave[i];
            pbw = pimatx.hc[kk][ii];
            pbo = pimatx.v[kk][ii];
            sumdif += pbo - pbw;
            sumcos += pbo*(1.0 - ang*ang);
        }
        sumdif = fabs(sumdif);
        if ( (sumdif > .001) && (sumcos > 0.001) )
          factor = sumdif/sumcos;
        if (factor > 1.0)
          factor = 1.0;
        if (minim_control.field == MM3)
          factor = 1.0;
    }
    for (i=0; i < pistuf.ntpi; i++)
    {
        ii = pitorsave.itsave[i];
        kk = pitorsave.ktsave[i];
        nc = pitorsave.ncsave[i];
        tc = pitorsave.tcsave[i];
        torsions.v2[nc] = factor*pimatx.v[kk][ii]*tc;
    }

    for (i=0; i < pistuf.nbpi; i++)
    {
        ii = pibondsave.ibsave[i];
        kk = pibondsave.kbsave[i];
        nc = pibondsave.nbsave[i];
        bkon = pibondsave.bksave[i];
        blen = pibondsave.blsave[i];
        slpk = pibondsave.pksave[i];
        slpl = pibondsave.plsave[i];
        slpl2 = pibondsave.plsave2[i];
        pbo = pimatx.hc[ii][kk];
        bonds_ff.bk[nc] = bkon - slpk + slpk*pbo;
        if (pbo < 0.0) pbo = 0.01;
        temp = pow(pbo,(1.0-slpl2));
        bonds_ff.bl[nc] = blen + slpl - slpl*temp;

        if (minim_control.field == MMX)
        {
           iti = atom[pistuf.listpi[ii]].mmx_type;
           itk = atom[pistuf.listpi[kk]].mmx_type;
           if (iti <= itk)
           {
             ikno = 100*iti + itk;
             ia = pistuf.listpi[kk];
           } else
           {
             ikno = 100*itk + iti;
             ia = pistuf.listpi[ii];
           }
           if (ikno == 215)
               bonds_ff.bl[nc] = blen + slpl - slpl*pbo;         
           continue;
        }
    }
}

void pirite()
{
    int i,j,k,l;
    
    for( i = 0; i < pistuf.norb; i++ )
    {
        for( j = 0; j < pistuf.norb; j++ )
        {
            pimatx.g[i][j] = pisave.gjsave[i][j];
            pimatx.f[i][j] = pisave.fjsave[i][j];
        }
    }

    if( pistuf.jprint )
        fprintf( pcmoutfile, "\nEIGENVALUES-Most to least stable \n");
    if( pistuf.nplane == 1 )
    {
      for( i = 0; i < pistuf.norb; i++ )
      {
          for( k = 0; k < pistuf.norb; k++ )
             pimatx.f[i][k] = pisave.fksave[i][k];
      }
    }else if( pistuf.nplane > 1 )
    {
      for( i = 0; i < pistuf.norb; i++ )
      {
         for( k = 0; k < pistuf.norb; k++ )
             pimatx.f[i][k] = pisave.fksave[i][k];
      }
    }else
    {
      for( i = 0; i < pistuf.norb; i++ )
      {
         for( k = 0; k < pistuf.norb; k++ )
             pimatx.f[i][k] = pisave.fksave[i][k];
      }
    }
    
    if (pistuf.jprint) printmatrix(pimatx.f," F Matrix ");

    l = 0;
    if( pistuf.jprint )
    {
       for( i = 0; i < pistuf.norb; i++ )
       {
          fprintf( pcmoutfile, " %8.4f", pimatx.en[i] );
          l++;
          if (l >=10 )
          {
              l = 0;
              fprintf(pcmoutfile,"\n");
          }
        }
        fprintf( pcmoutfile, "\n" );
     }

     for( i = 0; i < pistuf.norb; i++ )
     {
        for( k = 0; k < pistuf.norb; k++ )
             pimatx.v[i][k] = mtape.vmtape[i][k];
     }
     if (pistuf.jprint) printmatrix(pimatx.v," Molecular Orbitals ");
    
 /* use pen as a variable (real*4) to read densities from ktape */
     for( i = 0; i < pistuf.norb; i++ )
     {
         for( k = 0; k < pistuf.norb; k++ )
              penint.pen[i][k] = mtape.pensave[i][k];
         pimatx.g[i][i] = penint.pen[i][i];
     }
     if (pistuf.jprint)  printmatrix(densit.den," Electron Densit- Bond Order Matrix");
}
// ===============================
void rotat(int list[], float xr[], float yr[], float zr[])
{
        int i, l, i_, l_;
        float cose, denom, sine, xsav, xt, yt, zt;

        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         *
         *  THIS SUBROUTINE ROTATES THE I,K VECTOR SO THAT I IS AT ORIGIN AND K 
         *  ALONG THE X-AXIS.  THE ATOMS DEFINING THE RESPECTIVE PLANES ARE ALSO 
         *  ROTATED AND THE BOND LENGTHS NORMALIZED. 
         *
         * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         * */

        l = list[0];
        xt = atom[l].x;
        yt = atom[l].y;
        zt = atom[l].z;
        for( i = 0; i < 8; i++ )
        {
            l = list[i];
            xr[i] = atom[l].x - xt;
            yr[i] = atom[l].y - yt;
            zr[i] = atom[l].z - zt;
        }
        denom = sqrt( xr[1]*xr[1] + yr[1]*yr[1] );
        if( denom > 0.001 )
        {
           sine = yr[1]/denom;
           cose = xr[1]/denom;
           for( i = 0; i < 8; i++ )
           {
               xsav = xr[i];
               xr[i] = xr[i]*cose + yr[i]*sine;
               yr[i] = yr[i]*cose - xsav*sine;
           }
        }
        denom = sqrt( xr[1]*xr[1] + zr[1]*zr[1] );
        if( denom > 0.001 )
        {
           sine = zr[1]/denom;
           cose = xr[1]/denom;
           for( i = 0; i < 8; i++ )
           {
               xsav = xr[i];
               xr[i] = xr[i]*cose + zr[i]*sine;
               zr[i] = zr[i]*cose - xsav*sine;
            }
        }
        /*  THE NEXT LOOP 'NORMALIZES' THE COORDINATES (IE. MAKES ALL ATOMS ATTAC
         *  TO THE CENTRAL TWO ATOMS HAVE EQUAL BOND LEGENTHS) SO THAT THE ORBITA
         *  MAKES EQUAL ANGLES WITH THE BONDS RATHER THAN SIMPLY BEING PERPENDICU
         *  TO THE COMMON PLANE OF THE ATOMS.  */
        for( i = 3; i <= 8; i++ )
        {
                l = i/3;
                i_ = i-1;
                l_ = l-1;
                denom = (xr[i_] - xr[l_])*(xr[i_] - xr[l_]) + (yr[i_] - yr[l_])*(yr[i_] - yr[l_])
                 + (zr[i_] - zr[l_])*(zr[i_] - zr[l_]);
                /*  YOU DON'T WANT TO NORMALIZE NON-BONDED DISTANCES  */
                if( !((denom < 0.001) || (denom > 2.0)) )
                {
                        denom = sqrt( denom );
                        xr[i_] = (xr[i_] - (l - 1)*xr[l_])/denom + (l - 1)*xr[l_];
                        /*  NOTE:  THE NORMALIZATION OF THE X-COORDINATE IS COMPLICATED BY THE FA
                         *  THAT FOR ATOMS ATTACHED TO ATOM=2 YOU MUST FIRST TRANSLATE BACK TO TH
                         *  ORIGIN AND THEN NORMALIZE AND THEN RETRANSLATE BACK ALONG THE X-AXIS. */
                        yr[i_] /= denom;
                        zr[i_] /= denom;
                }
        }
        return;
}
       
        
      
                               
    
