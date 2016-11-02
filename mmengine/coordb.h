#ifndef CBONDS_H
#define CBONDS_H

struct t_cbonds {
        int nbnd, icbonds[MAXCBND][3], lpd[MAXCBND], ia1[MAXCBND],ia2[MAXCBND];
        double vrad[MAXCBND],veps[MAXCBND];
        } cbonds;

#endif

void coscalc(int,int,int,int,double *);
void ecoord_bond(void);
void ecoord_bond1(void);
void ecoord_bond2(int);
