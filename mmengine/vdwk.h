#ifndef VDWK_H
#define VDWK_H

struct t_vdw1 {
        int  nvdw;
        float rad[MAXVDWCONST], eps[MAXVDWCONST];
        int lpd[MAXVDWCONST], ihtyp[MAXVDWCONST], ihdon[MAXVDWCONST];
        float alpha[MAXVDWCONST],n[MAXVDWCONST],a[MAXVDWCONST],g[MAXVDWCONST];
        char da[MAXVDWCONST][2];
        } vdw1;

struct t_vdwpr_k {
        int nvdwpr;
        int  ia1[MAXBONDCONST],ia2[MAXBONDCONST];
        char kv[MAXBONDCONST][7];
        float radius[MAXBONDCONST], eps[MAXBONDCONST];
        } vdwpr_k;

#endif
