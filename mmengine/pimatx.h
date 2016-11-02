#ifndef PIMATX_H
#define PIMATX_H

struct t_pimatx {
       int npbmx, npiel;
       int *qp;
       float *eedd, *edd;
       float **hc, **v, **f,**g, **gamma, *en, *w, 
         *zor, *zt, *zz, *emz, *w1, *w2, *w3, *q;
        }       pimatx;

#endif

