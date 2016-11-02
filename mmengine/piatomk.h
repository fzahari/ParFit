#ifndef PIATOM_K_H
#define PIATOM_K_H

struct t_piatomk {
        int   npiatom;
        int   kat[MAXPIATOM], qp[MAXPIATOM];
        float q[MAXPIATOM], ion[MAXPIATOM], emz[MAXPIATOM], zor[MAXPIATOM], zz[MAXPIATOM];
        float w1[MAXPIATOM], w2[MAXPIATOM], w3[MAXPIATOM];
        }  piatomk;

#endif
