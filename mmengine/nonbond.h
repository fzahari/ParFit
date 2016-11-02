#ifndef NONBOND_H
#define NONBOND_H

struct t_nonbond {
        int nonbond, npair, maxnbtype, **iNBtype, **ipif;
        float **vrad, **veps, **vrad14,**veps14;
        } nonbond;
          
 /* EXTERN struct t_nonbond {
        int nonbond, npair, ntype, *iNBtype;
        float *vrad, *veps, *ipif;
        } nonbond;  */


/* EXTERN struct t_nonbond {
        int nonbond, iNB[MAXNB][3];
        float vrad[MAXNB], veps[MAXNB];
        } nonbond;  */
#endif
