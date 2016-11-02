#include "pcwin.h"
#include "pcmod.h"
#include "derivs.h"
#include "hess.h"
#include "rings.h"
#include "utility.h"
#include "nonbond.h"
#include "torsions.h"
#include "pot.h"
#include "attached.h"
#include "solv.h"
#include "get_mem.h"
#include "skip.h"


#define STILL 1
#define HCT   2

/* ================================================================ */
void get_memory()
{
    int i, j;
    int ntor;
    int ntypes, found, itype[MAXATOMTYPE];
    
    skip = ivector(0,natom+1);

    deriv.d1     = dmatrix(0,natom+1, 0,3);
    deriv.deb    = dmatrix(0,natom+1, 0,3);
    deriv.dea    = dmatrix(0,natom+1, 0,3);
    deriv.destb  = dmatrix(0,natom+1, 0,3);
    deriv.deopb  = dmatrix(0,natom+1, 0,3);
    deriv.detor  = dmatrix(0,natom+1, 0,3);
    deriv.de14  = dmatrix(0,natom+1, 0,3);
    deriv.devdw  = dmatrix(0,natom+1, 0,3);
    deriv.deqq   = dmatrix(0,natom+1, 0,3);
    deriv.deaa   = dmatrix(0,natom+1, 0,3);
    deriv.destor = dmatrix(0,natom+1, 0,3);   
    deriv.dehb = dmatrix(0,natom+1, 0,3);   
    deriv.deimprop = dmatrix(0,natom+1, 0,3);   
    deriv.deub = dmatrix(0,natom+1, 0,3);   
    deriv.desolv = dmatrix(0,natom+1, 0,3);
    deriv.drb = dvector(0,natom+1);  
    
    ntypes = 0;
    nonbond.maxnbtype = 0;
    
    for (i=1; i <= natom; i++)
    {
        found = FALSE;
        if (atom[i].type > nonbond.maxnbtype)
           nonbond.maxnbtype = atom[i].type;
        for (j=0; j < ntypes; j++)
        {
            if (atom[i].type == itype[j])
            {
                found = TRUE;
                break;
            }
        }
        if (found == FALSE)
        {
            itype[ntypes] = atom[i].type;
            ntypes++;
        }
    }

    j=0;
    for (i=0; i <= ntypes; i++)
       j+= i;
       
    nonbond.nonbond = j;    
    nonbond.iNBtype = imatrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.ipif =    imatrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.vrad = matrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.veps = matrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.vrad14 = matrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.veps14 = matrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);

    max_torsions();
    ntor = torsions.ntor;
    if (ntor > 0)
    {
       torsions.i14 = imatrix(0,ntor,0,4);
       torsions.v1 = vector(0,ntor);
       torsions.v2 = vector(0,ntor);
       torsions.v3 = vector(0,ntor);
       torsions.v4 = vector(0,ntor);
       torsions.v5 = vector(0,ntor);
       torsions.v6 = vector(0,ntor);
       torsions.ph1 = ivector(0,ntor);
       torsions.ph2 = ivector(0,ntor);
       torsions.ph3 = ivector(0,ntor);
       torsions.ph4 = ivector(0,ntor);
       torsions.ph5 = ivector(0,ntor);
       torsions.ph6 = ivector(0,ntor);
       for (i=0; i < ntor; i++)
       {
           torsions.v1[i] = 0.0;
           torsions.v2[i] = 0.0;
           torsions.v3[i] = 0.0;
           torsions.v4[i] = 0.0;
           torsions.v5[i] = 0.0;
           torsions.v6[i] = 0.0;
           torsions.ph1[i] = 0;
           torsions.ph2[i] = 0;
           torsions.ph3[i] = 0;
           torsions.ph4[i] = 0;
           torsions.ph5[i] = 0;
           torsions.ph6[i] = 0;
       }
    }    
    hess.hessx = matrix(0,natom+1, 0,3);
    hess.hessy = matrix(0,natom+1, 0,3);
    hess.hessz = matrix(0,natom+1, 0,3);

/*
    rings.r13 = imatrix(0,natom+1, 0,3);
    rings.r14 = imatrix(0,natom+1, 0,4);
    rings.r15 = imatrix(0,natom+1, 0,5);
    rings.r16 = imatrix(0,natom+1, 0,6);
*/

    attached.n13 = ivector(0,natom+1);
    attached.n14 = ivector(0,natom+1);
    attached.i13 = imatrix(0, 20, 0,natom+1);
    attached.i14 = imatrix(0,144, 0,natom+1);

    for(i=0; i <= natom; i++)
    {
        for (j=0; j < 3; j++)
        {
            deriv.d1[i][j] = 0.0;
            deriv.deb[i][j] = 0.0;
            deriv.dea[i][j] = 0.0;
            deriv.destb[i][j] = 0.0;
            deriv.deopb[i][j] = 0.0;
            deriv.detor[i][j] = 0.0;
            deriv.de14[i][j] = 0.0;
            deriv.devdw[i][j] = 0.0;
            deriv.deqq[i][j] = 0.0;
            deriv.deaa[i][j] = 0.0;
            deriv.destor[i][j] = 0.0;
            deriv.dehb[i][j] = 0.0;
            deriv.deimprop[i][j] = 0.0;
            deriv.deub[i][j] = 0.0;

            hess.hessx[i][j] = 0.0;
            hess.hessy[i][j] = 0.0;
            hess.hessz[i][j] = 0.0;
        }
    }
    if (pot.use_solv)
    {
        solvent.asolv = dvector(0,natom+1);
        solvent.rsolv = dvector(0,natom+1);
        solvent.rborn = dvector(0,natom+1);
        if (solvent.type == STILL)
        {
             solvent.vsolv = dvector(0,natom+1);
             solvent.gpol = dvector(0,natom+1);
        } else if (solvent.type == HCT)
        {
             solvent.shct = dvector(0,natom+1);
        }
    }  
}

void free_memory()
{
    int ntor;
    free_ivector(skip, 0, natom+1);      
    free_dmatrix(deriv.d1,      0,natom+1, 0,3);
    free_dmatrix(deriv.deb,     0,natom+1, 0,3);
    free_dmatrix(deriv.dea,     0,natom+1, 0,3);
    free_dmatrix(deriv.destb,   0,natom+1, 0,3);
    free_dmatrix(deriv.deopb,   0,natom+1, 0,3);
    free_dmatrix(deriv.detor,   0,natom+1, 0,3);
    free_dmatrix(deriv.de14,    0,natom+1, 0,3);
    free_dmatrix(deriv.devdw,   0,natom+1, 0,3);
    free_dmatrix(deriv.deqq,    0,natom+1, 0,3);
    free_dmatrix(deriv.deaa,    0,natom+1, 0,3);
    free_dmatrix(deriv.destor,  0,natom+1, 0,3);
    free_dmatrix(deriv.dehb,    0,natom+1, 0,3);
    free_dmatrix(deriv.deimprop,0,natom+1, 0,3);
    free_dmatrix(deriv.deub,    0,natom+1, 0,3);
    free_dmatrix(deriv.desolv ,0,natom+1, 0,3);   
    free_dvector(deriv.drb ,0,natom+1);  

    free_imatrix(nonbond.iNBtype ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_imatrix(nonbond.ipif ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_matrix(nonbond.vrad ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_matrix(nonbond.veps ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_matrix(nonbond.vrad14 ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_matrix(nonbond.veps14 ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);

//    ntor = 9*natom;
    ntor = torsions.ntor;
    if (ntor > 0)
    {
       free_imatrix(torsions.i14 ,0,ntor,0,4);
       free_vector(torsions.v1 ,0,ntor);
       free_vector(torsions.v2 ,0,ntor);
       free_vector(torsions.v3 ,0,ntor);
       free_vector(torsions.v4 ,0,ntor);
       free_vector(torsions.v5 ,0,ntor);
       free_vector(torsions.v6 ,0,ntor);
       free_ivector(torsions.ph1 ,0,ntor);
       free_ivector(torsions.ph2 ,0,ntor);
       free_ivector(torsions.ph3 ,0,ntor);
       free_ivector(torsions.ph4 ,0,ntor);
       free_ivector(torsions.ph5 ,0,ntor);
       free_ivector(torsions.ph6 ,0,ntor);
    }

    free_matrix(hess.hessx, 0,natom+1, 0,3);
    free_matrix(hess.hessy, 0,natom+1, 0,3);
    free_matrix(hess.hessz, 0,natom+1, 0,3);
    free_rings();
/*
    free_imatrix(rings.r13,0,natom+1, 0,3);
    free_imatrix(rings.r14,0,natom+1, 0,4);
    free_imatrix(rings.r15,0,natom+1, 0,5);
    free_imatrix(rings.r16,0,natom+1, 0,6);
*/

    free_ivector(attached.n13 ,0,natom+1);
    free_ivector(attached.n14 ,0,natom+1);
    free_imatrix(attached.i13 ,0, 20, 0,natom+1);
    free_imatrix(attached.i14 ,0,144, 0,natom+1);

    if (pot.use_solv)
    {
        free_dvector(solvent.asolv ,0,natom+1);
        free_dvector(solvent.rsolv ,0,natom+1);
        free_dvector(solvent.rborn ,0,natom+1);
        if (solvent.type == STILL)
        {
             free_dvector(solvent.vsolv ,0,natom+1);
             free_dvector(solvent.gpol ,0,natom+1);
        } else if (solvent.type == HCT)
        {
             free_dvector(solvent.shct ,0,natom+1);
        }
    }
}

