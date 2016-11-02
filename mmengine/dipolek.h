#ifndef DIPOLEK_H
#define DIPOLEK_H

struct  t_dipole_k {
        int ndipole,ndipole3,ndipole4,ndipole5;
        char   kb[MAXBONDCONST][7], kb3[MAXBOND3CONST][7], kb4[MAXBOND4CONST][7], kb5[MAXBOND5CONST][7];
        float bmom[MAXBONDCONST],bmom3[MAXBOND3CONST],bmom4[MAXBOND4CONST],bmom5[MAXBOND5CONST];
         } dipole_k;

#endif
