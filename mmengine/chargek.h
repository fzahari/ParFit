#ifndef CHARGEK_H
#define CHARGEK_H

struct t_charge_k {
         int ncharge, nbndchrg, nbndchrgdel;
         int type[MAXATOMTYPE], btype[MAXBONDCONST], btypedel[MAXBONDCONST];
         float charge[MAXATOMTYPE], bcharge[MAXBONDCONST], formchrg[MAXATOMTYPE], bchargedel[MAXBONDCONST];
         float typechrg[MAXATOMTYPE];
         } charge_k;

#endif
