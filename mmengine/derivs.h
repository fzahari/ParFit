#ifndef DERIV_H
#define DERIV_H

/* EXTERN struct t_deriv {
        double d1[MAXATOM][3],
               deb[MAXATOM][3],  dea[MAXATOM][3],  destb[MAXATOM][3],deopb[MAXATOM][3],
               detor[MAXATOM][3],de14[MAXATOM][3], devdw[MAXATOM][3],deqq[MAXATOM][3],
               deaa[MAXATOM][3], destor[MAXATOM][3], dehb[MAXATOM][3], deimprop[MAXATOM][3];
       } deriv;   */

struct t_deriv {
        double **d1,
               **deb,  **dea,  **destb, **deopb, **detor,**de14, **devdw, **deqq,
               **deaa, **destor, **dehb, **deimprop, **deub, **desolv, *drb;
        } deriv;  
        
#endif
