#ifndef BONDS_FF_H
#define BONDS_FF_H

struct t_bonds_ff {
        int nbnd, i12[MAXBND][2], idpl[MAXBND][2], index[MAXBND];
        double bk[MAXBND],bl[MAXBND], bmom[MAXBND];
        } bonds_ff;

#endif

