#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"
#include "pimatx.h"
#include "minim_values.h"
#include "pistuf.h"
#include "pidata.h"
#include "pisave.h"
#include "eheat.h"
#include "u.h"

int nocca, noccb;

static float **flocal,**vlocal,*elocal;
static int npiel = 0;

void  vescf()
{
        static int i, iter, ik, irand, j, k, line, nctr, no, nocc, npi, ncc;
        static int icore, iexcit, iuhf, npass;
        float **fb,**vb, *eig, *eigb;
        static float beta, betaij, btx, dbii, dbiigii, dbij, dbijhij, dbjj, denom, 
         dif0, diff, dii, diigii, dij, dijhij, djj, ecore, ecoul,  
         exch, fact, gamii, gii, hcii, pbiipjj, pii, 
         piipjj, pjj, qi, sumii, sumij, sumij2, sx, 
         vn, vo, vold, wi, wx, xzt, sctest;
        static float chg, pl;
        static int uhf = False;
        static double shift = 25.0;
 
         fb = matrix(0, pistuf.norb, 0, pistuf.norb);
         vb = matrix(0, pistuf.norb, 0, pistuf.norb);
         eig = vector(0, pistuf.norb);
         eigb = vector(0, pistuf.norb);

         /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         *
         *  CALCULATES THE MOLECULAR ORBITALS ACCORDING TO VESCF PROCEDURE.
         *  Extensively modified by John McKelvey to do open shell pi systems
         *  and includes improved convergence routines.
         *
         * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         * */

         iter = 80;

        /*  BTX, WX AND SX ARE THE STANDARD BETA, STANDARD IONIZATION POTENTIAL  AND STANDARD OVERLAP. */
        if( uvnpi.ifirst == TRUE )
        {
                chg = 0.;
                icore = 0;
                iexcit = 0;
                npi = pimatx.npiel*2;
                irand = 0;
                noccb = 0;
                sctest = 1/10000.;
// uhf calc
                if( pistuf.iuv == 2 )
                {
                        nocca = npi/2;
                        uhf = True;
                        noccb = (npiel - (elowht.mult - 1))/2;
                        nocca = npiel - noccb;
                        iuhf = 1;
                } 
// rhf calc
               else
                { 
                        nocca = npi/2;
                        elowht.mult = 1;
                        iuhf = 0;
                        uhf = False;
                } 
                elowht.uhf = uhf;
                npass = 0;
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( j = 0; j < pistuf.norb; j++ )
                        {
                                pimatx.g[i][j] = 0.;
                                densit.dena[i][j] = 0.;
                                densit.denb[i][j] = 0.;
                                densit.den[i][j] = 0.;
                        }
                        penint.ed[i] = pimatx.q[i];
                        pimatx.zt[i] = pimatx.zor[i];
                        if( uhf )
                        {
                                densit.denb[i][i] = pimatx.q[i]*(float)( noccb )/
                                 (float)( nocca + noccb );
                                densit.dena[i][i] = pimatx.q[i]*(float)( nocca )/
                                 (float)( nocca + noccb );
                        }
                        else
                        {
                                densit.dena[i][i] = 0.5*pimatx.q[i];
                        }
                        densit.den[i][i] = pimatx.q[i];
                        pimatx.g[i][i] = pimatx.q[i];
                }
        }
        pl = 1.0;
                
        if( icore )
                pimatx.q[icore] = chg;
        if( !pimatx.npbmx )
        {
            if( pistuf.jprint )
            {
                if( minim_values.iprint)
                    {
                        fprintf( pcmoutfile,"MOLECULAR ORBITALS BY VESCF CALCULATION (PLANAR)\n" );
                        fprintf( pcmoutfile,"          UNIT = ELECTRON VOLT (EV)\n" );
                    }
            }
        }
        else
        {
            if( pistuf.jprint )
            {
                if( minim_values.iprint)
                {
                        fprintf( pcmoutfile,"MOLECULAR ORBITALS BY VESCF CALCULATION (ACTUAL)\n" );
                        fprintf( pcmoutfile,"          UNIT = ELECTRON VOLT (EV)\n" );
                }
           }
        }
        btx = -2.5;
        line = ((pistuf.norb*pistuf.norb) + pistuf.norb)/2;
        wx = -8.52893;
        sx = 0.32818;
        beta = btx/(wx*sx);

        no = pistuf.norb - 1;
        vold = 0.;
        nctr = 0;
        npass += 1;
        for( nctr = 0; nctr < iter; nctr++ )
        {
//                printmatrix(pimatx.f,"f at start");
                if( pl > 0.1 || nctr == 0 )
                {
                       frag_diag(pimatx.f, pimatx.v, eig );
//                       xdiag( pimatx.f, pimatx.v, eig );
                       if( nctr == 1 )
                        {
                                for( i = 0; i < pistuf.norb; i++ )
                                {
                                        for( j = 0; j < pistuf.norb; j++ )
                                        {
                                                vb[i][j] = pimatx.v[i][j];
                                        }
                                }
                        }
//                        if( uhf && nctr > 1 )
//                                xdiag( fb, vb, eigb );
                }
                else if( pl <= 0.1 && nctr > 1 )
                {
                        frag_diag(pimatx.f, pimatx.v, eig );
//                        xdiag( pimatx.f, pimatx.v, eig );
//                        if( uhf && nctr > 1 )
//                                xdiag( fb, vb, eigb );
                }
//                printmatrix(pimatx.f,"f aft frag_diag");
//                printmatrix(pimatx.v,"v aft frag_diag");
                /*  THE NEXT LOOP CALCULATES ELECTRON DENSITIES AND USES THEM TO ARRIVE
                 *  AT A NEW ZETA FOR EACH ORBITAL.
                 *  DENSITY MATRIX */
                fact = 1.0;
                densty( nocca, pimatx.v, densit.dena, fact );
//                printmatrix(pimatx.v,"v aft densty");
//                printmatrix(densit.dena,"dena aft densty");
                if( uhf )
                {
                        densty( noccb, vb, densit.denb, fact );
                        for( i = 0; i < pistuf.norb; i++ )
                        {
                                for( j = 0; j <= i; j++ )
                                {
                                   densit.den[j][i] = densit.dena[j][i] + densit.denb[j][i];
                                   densit.den[i][j] = densit.den[j][i];
                                }
                        }
                }
                else
                {
                        for( i = 0; i < pistuf.norb; i++ )
                        {
                                for( j = 0; j <= i; j++ )
                                {
                                        densit.den[i][j] = 2.0*densit.dena[i][j];
                                        densit.den[j][i] = densit.den[i][j];
                                }
                        }
                }
                /*..CALL ZERNER SCFCHK HERE... */
                ncc = npi;
                if (uhf)
                   ncc = nocca;
                scfchk( nctr,ncc, &iuhf, &pl );
//                printmatrix(densit.den,"den aft schchk");
                
//                printmatrix(densit.den,"den aft1");

                if( uhf )
                {
                        for( i = 0; i < pistuf.norb; i++ )
                        {
                                for( j = 0; j <= i; j++ )
                                {
                                        densit.den[i][j] = densit.dena[i][j] + densit.denb[i][j];
                                        pimatx.g[j][i] = densit.den[i][j];
                                        densit.den[j][i] = densit.den[i][j];
                                }
                        }
                }
                else
                {
                        for( i = 0; i < pistuf.norb; i++ )
                        {
                                for( j = 0; j <= i; j++ )
                                {
                                        pimatx.g[j][i] = densit.den[i][j];
                                        densit.dena[i][j] = 0.5*densit.den[i][j];
                                        densit.dena[j][i] = densit.dena[i][j];
                                }
                        }
                }

                for( i = 0; i <  pistuf.norb; i++ )
                {
                        ed123.ed5[i] = ed123.ed4[i];
                        ed123.ed4[i] = ed123.ed3[i];
                        ed123.ed3[i] = ed123.ed2[i];
                        ed123.ed2[i] = ed123.ed1[i];
                        ed123.ed1[i] = densit.den[i][i];
                        if (ed123.nscft >= 4)
                        {
                            denom = 15.0;
                            pimatx.edd[i] = (5*ed123.ed1[i] + 4*ed123.ed2[i] + 
                               3*ed123.ed3[i] + 2*ed123.ed4[i] + ed123.ed5[i])/denom;
                            densit.den[i][i]  = pimatx.edd[i];
                            densit.dena[i][i] = 0.5*densit.den[i][i];
                            densit.denb[i][i] = 0.5*densit.den[i][i];                            
                            pimatx.g[i][i] = densit.den[i][i];
                        } else if (ed123.nscft == 3)
                        {
                            denom = 10.0;
                            pimatx.edd[i] = (4.0*ed123.ed1[i] + 3.0*ed123.ed2[i] + 2.0*ed123.ed3[i] + ed123.ed4[i])/denom;
                            densit.den[i][i]  = pimatx.edd[i];
                            densit.dena[i][i] = 0.5*densit.den[i][i];
                            densit.denb[i][i] = 0.5*densit.den[i][i];                            
                            pimatx.g[i][i] = densit.den[i][i];
                        } else if( ed123.nscft == 2 ) 
                        {
                            denom = 6.0;
                            pimatx.edd[i] = (3.0*ed123.ed1[i] + 2.0*ed123.ed2[i] + ed123.ed3[i])/denom;
                            densit.den[i][i]  = pimatx.edd[i];
                            densit.dena[i][i] = 0.5*densit.den[i][i];
                            densit.denb[i][i] = 0.5*densit.den[i][i];                            
                            pimatx.g[i][i] = densit.den[i][i];
                        } else if( ed123.nscft == 1 )
                        {
                            denom = 3.0;
                            pimatx.edd[i] = (2*ed123.ed1[i] + ed123.ed2[i])/denom;
                            densit.den[i][i]  = pimatx.edd[i];
                            densit.dena[i][i] = 0.5*densit.den[i][i];
                            densit.denb[i][i] = 0.5*densit.den[i][i];                            
                            pimatx.g[i][i] = densit.den[i][i];
                        } else
                        {
                            densit.dena[i][i] = 0.5*densit.den[i][i];
                            densit.denb[i][i] = 0.5*densit.den[i][i];                            
                            pimatx.edd[i] = densit.den[i][i];
                            pimatx.g[i][i] = densit.den[i][i];
                        }
//                      printmatrix(densit.den,"den aft2");
//                      printmatrix(pimatx.g,"g aft2");
                        /*..THE FOLLOWING STATEMENT IS True FOR ONLY SECOND ROW ATOMS (LI-F): */
                      // change for sulfur by jjg
                        if (atom[pistuf.listpi[i]].type == 15 || atom[pistuf.listpi[i]].type == 38)
                           pimatx.zt[i] = pimatx.zor[i];
                        else
                          pimatx.zt[i] = (2*pimatx.zor[i] - (pimatx.edd[i] - pimatx.q[i])*.35)/2.;
                }
                /*  FROM HERE TO STATEMENT 230 SYMMETRY RELATED ORBITALS ARE FORCED TO
                 *  HAVE THE SAME CHARGE (WHICH IS JUST THE AVERAGE OF THE ORBITALS). */
                xzt = 0;
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( k = 0; k < pistuf.norb; k++ )
                        {
                                pimatx.hc[i][k] = 0.;
                        }
                }

                /*  SUBROUTINE OTABLE CALCULATES THE OVERLAP INTEGRALS (RETURNED IN F),
                 *  THE THEORETICAL REPULSION INTEGRALS (RETURNED IN G), THE OVERLAP
                 *  PORTION OF THE CORE INTEGRALS (RETURNED IN HC) AND THE PENETRATION
                 *  INTEGRALS (RETURNED IN PEN).
                 *         CALL OTABLE(natom,ZT,F,G,HC,QP,PEN,ITAPE) */
//                printmatrix(pimatx.g,"g bef otable");
//                printmatrix(pimatx.hc,"hc bef otable");
                otable();
//                printmatrix(pimatx.f,"f aft otable");
//                printmatrix(pimatx.g,"g aft otable");
//                printmatrix(pimatx.hc,"hc aft otable");

                for( i = 0; i < pistuf.norb; i++ )
                {
                        pimatx.gamma[i][i] = 2.*pimatx.zt[i]*pimatx.emz[i]/pimatx.zz[i];
                }
                for( i = 0; i < no; i++ )
                {
                        gii = pimatx.g[i][i];
                        gamii = pimatx.gamma[i][i];
                        ik = i + 1;
                        for( k = ik; k < pistuf.norb; k++ )
                        {
                                /*  THE OVERLAP IS NOT SCALED DOWN IN CALCULATING THE REPULSION INTEGRAL.
                                 *  THE REPULSION INTEGRALS ARE SCALED DOWN BY THE PRODUCT OF OVERLAP AND
                                 *  THE AVERAGE OF THE DIFFERENCES BETWEEN THE THEORETICAL AND EMPIRICAL
                                 *  ONE-CENTER REPULSION INTEGRALS. */
                                pimatx.gamma[i][k] = pimatx.g[i][k] - .5*pimatx.f[i][k]*
                                 (gii + pimatx.g[k][k] - gamii - pimatx.gamma[k][k]);
                                pimatx.gamma[k][i] = pimatx.gamma[i][k];
                        }
                }
                 for( i = 0; i < pistuf.norb; i++ )
                {
                        /*  IN CALCULATING ALPHA THE ONE CENTER REPULSION INTEGRAL IS ADDED FOR
                         *  ORBITALS WITH TWO ELECTRONS. THE FIRST IONIZATION POTENTIAL FOR THE
                         *  ORBITAL MUST BE USED. */
                        pimatx.w[i] = -((2.*pimatx.zt[i])*(2.*pimatx.zt[i])*pimatx.w1[i] + 
                                2.*pimatx.zt[i]*pimatx.w2[i] + pimatx.w3[i]);
                        wi = pimatx.w[i];
                        hcii = pimatx.w[i] - (pimatx.q[i] - 1.)*pimatx.gamma[i][i];
                        for( k = 0; k < pistuf.norb; k++ )
                        {
                                if( i != k )
                                        hcii += -pimatx.q[k]*pimatx.gamma[i][k] - penint.pen[i][k];
                        }
                        pimatx.hc[i][i] = hcii;
                }
                for( i = 0; i < no; i++ )
                {
                        wi = pimatx.w[i];
                        qi = pimatx.q[i];
                        gii = pimatx.gamma[i][i];
                        ik = i + 1;
                        for( k = ik; k < pistuf.norb; k++ )
                        {
                                pimatx.hc[i][k] *= .5*beta*(wi + pimatx.w[k] - ((qi - 1.)*
                                 gii + (pimatx.q[k] - 1.)*pimatx.gamma[k][k]));
                                pimatx.hc[k][i] = pimatx.hc[i][k];
                        }
                }
               /*..CALCULATE FOCK MATRIX ELEMENTS */
                fofgen( pimatx.f, pimatx.hc, densit.den, densit.dena, &pistuf.norb, 
                 pimatx.gamma );
                if( uhf )
                        fofgen( fb, pimatx.hc, densit.den, densit.denb, &pistuf.norb, 
                         pimatx.gamma );
// shift levels
                if (nctr >= 3)
                {
                    shift = shift/(pow((1.0+(float)nctr),1.25));
                    for (i=0; i < pistuf.norb; i++)
                    {
                        for (j=0; j < i; j++)
                        {
                            pimatx.f[i][j] += shift*densit.dena[i][j];
                            pimatx.f[j][i] = pimatx.f[i][j];
                        }
                        pimatx.f[i][i] -= shift;
                    }
                }
//                printmatrix(pimatx.f,"f after shift");
                /*...
                 *...CALCULATE TOTAL ELECTRONIC ENERGY AND ADD NUCLEAR REPULSION ENERGY:
                 *... */
                ecore = 0.;
                exch = 0.;
                vo = 0.e0;
                vn = 0.e0;
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( j = 0; j < pistuf.norb; j++ )
                        {
                                ecore += pimatx.hc[i][j]*densit.den[i][j];
                                exch += -pimatx.gamma[i][j]*(densit.den[i][j]*densit.den[i][j]);
                                vo += 0.5*densit.dena[i][j]*(pimatx.f[i][j] + pimatx.hc[i][j]);
                         }
                        if( uhf )
                        {
                                for( j = 0; j < pistuf.norb; j++ )
                                {
                                        vo += 0.5*densit.denb[i][j]*(fb[i][j] + pimatx.hc[i][j]);
                                }
                        }
                        for( j = 0; j <= i; j++ )
                        {
                                vn += pimatx.q[i]*pimatx.q[j]*pimatx.gamma[i][j];
                        }
                        /*..SUBTRACT OFF THE 1-CENTER NUCLEAR REPULSION TERM ADDED FOR I=J IN THE
                         *..ABOVE LOOP!! */
                        vn += -pimatx.q[i]*pimatx.q[i]*pimatx.gamma[i][i];
                }
                if( !uhf )
                        vo *= 2.e0;
                exch *= 0.25;
                ecoul = vo - exch - ecore;
                vo += vn;
                diff = 0.;
                if( nctr > 0 )
                        diff = fabs( vold - vo );
                dif0 = vo - vold;
                vold = vo;
                if( pistuf.jprint )
                        fprintf( pcmoutfile,"Pi iteration %d V0 = %g EV    DIFF =%g  PL = %g\n", nctr, vo, dif0,pl );

                if( nctr != 0 )
                {
                        sctest = 0.00001;
                        piheat.tepi = diff*23.061;
                       if( diff <= sctest || (pl <= sctest) )
//                        if( diff <= sctest )
                                goto L_9877;
                }
        }

        /*  IF THIS POINT IS REACHED THE SYSTEM HAS NOT REACHED SELF-CONSISTENCY. */
        if( pistuf.jprint )
              fprintf( pcmoutfile,"THE VESCF CALC. HAS NOT REACHED SELF-CONSISTENCY\n" );
L_9877:
        if( iexcit == 1 && !iuhf )
        {
                if( pistuf.jprint )
                    fprintf( pcmoutfile," EXCITED STATE CORRECTION TO DENSITY\n" );
                nocc = nocca;
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( j = 0; j <= i; j++ )
                        {
                                densit.den[j][i] += pimatx.v[nocc][i]*pimatx.v[nocc][j] - 
                                 pimatx.v[nocc-1][i]*pimatx.v[nocc-1][j];
                                densit.den[i][j] = densit.den[j][i];
                                pimatx.g[i][j] = densit.den[i][j];
                        }
                }
        }
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( k = 0; k < pistuf.norb; k++ )
                        {
                                mtape.vmtape[i][k] = pimatx.v[i][k];
                        }
                }
                if( uhf )
                {
                        for( i = 0; i < pistuf.norb; i++ )
                        {
                                for( k = 0; k < pistuf.norb; k++ )
                                {
                                        mtape.vbtape[i][k] = vb[i][k];
                                }
                        }
                }

        uvnpi.ifirst = FALSE;
        fofgen( pimatx.f, pimatx.hc, densit.den, densit.dena, &pistuf.norb, 
         pimatx.gamma );
        if( uhf )
                fofgen( fb, pimatx.hc, densit.den, densit.denb, &pistuf.norb, 
                 pimatx.gamma );
        for( i = 0; i < pistuf.norb; i++ )
        {
                for( j = 0; j < pistuf.norb; j++ )
                {
                        if( uhf )
                                pimatx.f[i][j] += fb[i][j];
                        pimatx.g[i][j] = densit.den[i][j];
                }
        }
        elowht.epi = 0.;
        sumii = 0.;
        sumij = 0.;
        sumij2 = 0.;
        for( i = 0; i < pistuf.norb; i++ )
        {
                if( uhf )
                {
                        dii = densit.dena[i][i];
                        dbii = densit.denb[i][i];
                        dbiigii = dii*dbii*pimatx.gamma[i][i];
                        elowht.epi -= dbiigii;
                }
                else
                {
                        dii = densit.den[i][i];
                        diigii = dii*dii*pimatx.gamma[i][i]*.25;
                        sumii += diigii;
                        elowht.epi -= diigii;
                }
                if( i == pistuf.norb )
                        break;
                for( j = i + 1; j < pistuf.norb; j++ )
                {
                        betaij = pimatx.hc[i][j];
                        if( uhf )
                        {
                                djj = densit.dena[j][j];
                                dij = densit.dena[i][j];
                                dbjj = densit.denb[j][j];
                                dbij = densit.denb[i][j];
                                dbijhij = 2.*(dbij + dij)*betaij;
                                pii = dii + dbii;
                                pjj = djj + dbjj;
                                pbiipjj = ((pii - pimatx.q[i])*(pjj - pimatx.q[j]) - (dij*dij + 
                                 dbij*dbij))*pimatx.gamma[i][j];
                                elowht.epi += -dbijhij - pbiipjj;
                        }
                        else
                        {
                                djj = densit.den[j][j];
                                dij = densit.den[i][j];
                                dijhij = 2.*dij*betaij;
                                sumij += dijhij;
                                piipjj = ((dii - pimatx.q[i])*(djj - pimatx.q[j]) - .5*dij*dij)*
                                 pimatx.gamma[i][j];
                                sumij2 += piipjj;
                                elowht.epi += -dijhij - piipjj;
                        }
                }
        }
         free_matrix(fb ,0, pistuf.norb, 0, pistuf.norb);
         free_matrix(vb ,0, pistuf.norb, 0, pistuf.norb);
         free_vector(eig ,0, pistuf.norb);
         free_vector(eigb ,0, pistuf.norb);

        return;
} 

void fofgen(float **f,float **hc,float **pt,float **px,int *norb,float **gamma)
{
        int ip, jp;

        for( ip = 0; ip < *norb; ip++ )
        {
                /*..CORE HAMILTONIAN + DIAGONAL EXCHANGE */
                f[ip][ip] = hc[ip][ip] + (pt[ip][ip] - px[ip][ip])* gamma[ip][ip];
                /*..GET ELECTROSTATIC INTERACTIONS FROM NEIGHBORS */
                for( jp = 0; jp <= (ip - 1); jp++ )
                {
                        /*..RESONANCE AND OFF DIAGONAL EXCHANGE */
                        f[ip][jp] = hc[ip][jp] - px[ip][jp]*gamma[ip][jp];
                        f[jp][ip] = f[ip][jp];
                        f[ip][ip] += pt[jp][jp]*gamma[ip][jp];
                }
                if( ip != *norb )
                {
                        for( jp = ip + 1; jp < *norb; jp++ )
                        {
                                f[ip][ip] += pt[jp][jp]*gamma[ip][jp];
                        }
                }
        }
        return;
} 

void densty(int nocc,float **v,float **p,float fact)
{
        int i, j, k;
        float s;
        for( i = 0; i < pistuf.norb; i++ )
        {
                for( j = 0; j <= i; j++ )
                {
                        s = 0.;
                        for( k = 0; k < nocc; k++ )
                        {
                                s += v[k][i]*v[k][j]*fact;
                        }
                        p[i][j] = s;
                        p[j][i] = s;                        
                }
        }
        return;
}
// =========================================
void frag_diag(float **f,float **v,float *eig)
{
    int i,nocc_all,nocc2,nvirt2,nocc,nvirt1;
    int nocci,naoi,nvirti,i1,i2,nao1,nao2,nocc1;
    int j,k,jj,kk,m,l;
    float fact;
    
// local fock matrix
    flocal = matrix(0, pistuf.norb, 0, pistuf.norb);
    vlocal = matrix(0, pistuf.norb, 0, pistuf.norb);
    elocal = vector(0, pistuf.norb);
//
    for (i=0; i < pistuf.norb; i++)
    {
        for (j=0; j < pistuf.norb; j++)
        {
            vlocal[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
// 
    nocc_all = 0;
    for (i=0; i < pigroup.nfrag; i++)
       nocc_all += pigroup.nelecs[i]/2;
    fact = 2.0;
    i2 = 0;
    nocc2 = 0;
    nvirt2 = nocc_all;
    nocc = 0;
    nvirt1 = 0;
    nao2 = 0;
//
    for (i=0; i < pigroup.nfrag; i++)
    {
        nocci = pigroup.nelecs[i]/2;
        naoi = pigroup.nao[i];
        nvirti = naoi - nocci;
        i1 = i2;
        i2 = i2 + naoi;
        nao1 = nao2;
        nao2 = nao2 + naoi;
        nocc1 = nocc2;
        nocc2 = nocc2 + nocci;
        nvirt1 = nvirt2;
        nvirt2 = nvirt2 + nvirti;

        jj = 0;
        for (j = i1; j < i2; j++)
        {
            kk = 0;
            for (k=i1; k < i2; k++)
            {
               flocal[jj][kk] = f[j][k];
               kk++;
            }
            jj++;
        }
        xdiag(naoi,flocal,vlocal,elocal);
// copy back eigenvalues
        k = 0;
        for (j=nocc1; j < nocc2; j++)
        {
            eig[j] = elocal[k];
            k++;
        }
        for (j=nvirt1; j < nvirt2; j++)
        {
            eig[j] = elocal[k];
            k++;
        }
// copy back eigenvectors
        k = 0;
        for (j = nocc1; j < nocc2; j++)
        {
            m = 0;
            for (l=nao1; l < nao2; l++)
            {
                v[j][l] = vlocal[k][m];
                m++;
            }
            k++;
        }
        for (j = nvirt1; j < nvirt2; j++)
        {
            m = 0;
            for (l=nao1; l < nao2; l++)
            {
                v[j][l] = vlocal[k][m];
                m++;
            }
            k++;
        }
    }
    free_matrix(flocal ,0, pistuf.norb, 0, pistuf.norb);
    free_matrix(vlocal ,0, pistuf.norb, 0, pistuf.norb);
    free_vector(elocal ,0, pistuf.norb);
}              
// =========================================
void xdiag(int norb,float **f,float **v,float *eig)
{
        int i, k, k1, n1;
        float hcii,ed[MAXORB];

        diag(norb,f);
        for( i = 0; i < norb; i++ )
        {
                f[i][i] = -f[i][i];
                ed[i] = 1.00001;
        }
        for( i = 0; i < (norb - 1); i++ )
        {
                hcii = f[i][i];
                k1 = i + 1;
                for( k = k1; k < norb; k++ )
                {
                        if( hcii >= f[k][k] )
                        {
                                ed[i] += 1.00001;
                        }
                        else
                        {
                                ed[k] += 1.00001;
                        }
                }
        }
        n1 = norb + 1;
        for( i = 0; i < norb; i++ )
        {
                k1 = n1 - ed[i] + 0.1 - 1;
                pimatx.en[k1] = -f[i][i];
                eig[k1] = -f[i][i];
                for( k = 0; k < norb; k++ )
                {
                        v[k1][k] = pimatx.hc[k][i];
                }
        }
        return;
} 
// ==================================================
void scfchk(int nctr, int nael,int *annisw,float *pl)
{
        int i, ij, itt, j, k, l, na, nrr;
        static float fac, sa, slp, sum, tden, test;

        na = pistuf.norb;
        nrr = (pistuf.norb*pistuf.norb + pistuf.norb)/2;
        test = 1.e-8;
        if( uvnpi.kfirst == TRUE )
        {
                for( i = 0; i < pistuf.norb; i++ )
                {
                        local.d12[i] = 0.;
                        local.d13[i] = 0.;
                }
                uvnpi.kfirst = FALSE;
        }

        if( *annisw == 0 )
        {
                linear( densit.den, local.f1, na );
                for( i = 0; i < nrr; i++ )
                {
                        local.fa1[i] = 0.;
                }
        }
        else
        {
                linear( densit.dena, local.f1, na );
                linear( densit.denb, local.fa1,na );
        }


        na = pistuf.norb;
        for( i = 0; i < na; i++ )
        {
                local.fn[i] = 0.0;
                local.d13[i] = local.d12[i];
                local.d12[i] = 0.0;
        }
        l = 0;
        for( i = 1; i <= pistuf.norb; i++ )
        {
                k = i-1;
                l += i;
                local.d12[k] += local.f1[l-1] + local.fa1[l-1];
        }


        if( nctr >= 3 )                // old nscft -1 > = 2
        {
                if( *annisw != 0)
                {
                        for( i = 0; i < nrr; i++ )
                        {
                                local.f22[i] = local.fa1[i];
                        }
                }
                for( i = 0; i < na; i++ )
                {
                        tden = local.d13[i] - local.adnmt[i];
                        if( fabs( tden ) >= test )
                        {
                                slp = (local.d12[i] - local.adnmo[i])/tden;
                                if( slp <= 0.0 && slp >= -19.0 )
                                        local.fn[i] = slp/(slp - 1.0);
                        }
                        local.adnmt[i] = local.d13[i];
                        local.adnmo[i] = local.d12[i];
                }
                sum = 0.0;
                for( i = 0; i < na; i++ )
                {
                        sum += local.fn[i];
                }
                fac = sum/ pistuf.norb;
                if( fac > 1.0 )
                        fac = 1.0;
                if( fac < 0.0 )
                        fac = 0.0;
        }
        else
        {
                itt = nctr-1 + 1;
                if( itt == 2 )
                {
                        if( *annisw != 0)
                        {
                                for( i = 0; i < nrr; i++ )
                                {
                                        local.f22[i] = local.fa1[i];
                                }
                        }
                        for( i = 0; i < na; i++ )
                        {
                                local.adnmo[i] = local.d12[i];
                        }
                        fac = 0.50;
                }
                else
                {
                        test = 1.e-8;
                        for( i = 0; i < na; i++ )
                        {
                                local.adnmt[i] = local.d12[i];
                        }
                        for( i = 0; i < nrr; i++ )
                        {
                                local.f29[i] = local.f1[i];
                        }
                        if( *annisw != 0)
                        {
                                for( i = 0; i < nrr; i++ )
                                {
                                        local.f22[i] = local.fa1[i];
                                        local.f27[i] = local.fa1[i];
                                }
                        }
                        goto L_14002;
                }
        }
        for( i = 0; i < nrr; i++ )
        {
                local.fa1[i] = local.f29[i];
        }
        sa = 0.;
        ij = 0;
        *pl = 0.0;
        for( i = 0; i < pistuf.norb; i++ )
        {
                for( j = 0; j <= i; j++ )
                {
                        if (fabs( local.f1[ij]-local.fa1[ij]) > *pl )
                                *pl = fabs( local.f1[ij] - local.fa1[ij] );
                        local.f1[ij] = fac*local.fa1[ij] + (1.0 - fac)*local.f1[ij];
                        ij += 1;
                }
        }
/* =========================================  */
        j = 0;
        for( i = 1; i <= pistuf.norb; i++ )
        {
                j += i;
                local.da[i-1] = local.f1[j-1];
        }
        for( i = 0; i < nrr; i++ )
        {
                local.f29[i] = local.f1[i];
        }
        if( *annisw != 0 )
        {
                /*..retrieve current beta den matrix */
                for( i = 0; i < nrr; i++ )
                {
                        local.fa1[i] = local.f22[i];
                }
                /*..retrieve previous beta den matrix */
                for( i = 0; i < nrr; i++ )
                {
                        local.f1[i] = local.f27[i];
                }
                ij = 0;
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( j = 0; j <= i; j++ )
                        {
                                if (fabs( local.f1[ij]-local.fa1[ij]) > *pl)
                                        *pl = fabs( local.f1[ij] - local.fa1[ij] );
                                local.fa1[ij] = fac*local.f1[ij] + (1.0 - fac)*local.fa1[ij];
                                ij += 1;
                        }
                }
/* ============================================= */
                j = 0;
                for( i = 1; i <= pistuf.norb; i++ )
                {
                        j += i;
                        local.da[i-1] += local.fa1[j-1];
                }
                for( i = 0; i < nrr; i++ )
                {
                        local.f27[i] = local.fa1[i];
                }
                /*..restore alpha den to f */
                for( i = 0; i < nrr; i++ )
                {
                        local.f1[i] = local.f29[i];
                }
        }
L_14002:
        if( *annisw == 0 )
        {
                square( local.f1, densit.den, na );
        }
        else
        {
                square( local.f1, densit.dena, na );
                square( local.fa1, densit.denb, na );
        }
        return;
}
/* -------------------------------- */
void square(float *a,float **b,int nb)
{
        int i, ix, j, jx, k;

        k = (nb*(nb + 1))/2 -1;
        for( j = 0; j < nb; j++ )
        {
                jx = nb - j - 1;
                for( i = 0; i <= jx; i++ )
                {
                        ix = jx - i;
                        b[ix][jx] = a[k];
                        k -= 1;
                }
        }
        for( j = 0; j < nb; j++ )
        {
                for( i = 0; i <= j; i++ )
                {
                        b[j][i] = b[i][j];
                }
        }
        return;
} 


void  linear(float **a,float *b,int n)
{
        int i, j, k;

        k = 0;
        for( j = 0; j < n; j++ )
        {
                for( i = 0; i <= j; i++ )
                {
                        b[k] = a[i][j];
                        k += 1;
                }
        }
        return;
} 



void randnum(long int *i,float *x)
{
        float c, y;
        double a, b, d, e, z;

        d = 48828125.;
        e = 2097152.;
        y = (float)( *i );
        z = (double)( y );
        a = z*d;
        b = fmod( a, e );
        c = (float)( b );
        *i = c;
        *x = c/e;
        return;
} /* end of function */
/*      ---------------------- */
void printarray(float *a, char *s)
{
    int i;
    fprintf(pcmoutfile," %s \n",s);
    for (i=0; i < pistuf.norb; i++)
        fprintf(pcmoutfile," %8.4f",a[i]);
    fprintf(pcmoutfile,"\n");
    return;
}
// ================================        
void printmatrix(float **a, char *s)
{
    int i,j;
    fprintf(pcmoutfile," %s \n",s);
    for (i=0; i < pistuf.norb; i++)
    {
        for (j = 0; j < pistuf.norb; j++)
        {
            fprintf(pcmoutfile," %8.4f",a[i][j]);
        }
        fprintf(pcmoutfile,"\n");
    }
    fprintf(pcmoutfile,"\n");
    return;
}
// =============================================
void printcont(float **a,int nx, int ny, char *s)
{
    int i,j;
    fprintf(pcmoutfile," %s \n",s);
    for (i=0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf(pcmoutfile," %5.2f",a[i][j]);
        }
        fprintf(pcmoutfile,"\n");
    }
    fprintf(pcmoutfile,"\n");
    return;
}
