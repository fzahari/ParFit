#ifndef BONDK_H
#define BONDK_H

struct  t_bondk1 {
        int use_ring3, use_ring4, use_ring5;
        int nbnd, nbnd3, nbnd4, nbnd5, ndeloc;
        char kb[MAXBONDCONST][7], kb3[MAXBOND3CONST][7],kb4[MAXBOND4CONST][7],kb5[MAXBOND5CONST][7];
        char kbdel[MAXBONDDELOC][7];
        float s[MAXBONDCONST], t[MAXBONDCONST];
        float s3[MAXBOND3CONST], t3[MAXBOND3CONST];
        float s4[MAXBOND4CONST], t4[MAXBOND4CONST];
        float s5[MAXBOND5CONST], t5[MAXBOND5CONST];
        float sdel[MAXBONDDELOC], tdel[MAXBONDDELOC];
        }  bondk1;
                
        
struct t_bondpk {
        int     npibond;
        int     use_pibond;
        char    kb[MAXPIBONDCONST][7];
        float   bk[MAXPIBONDCONST], bl[MAXPIBONDCONST], bmom[MAXPIBONDCONST],sslop[MAXPIBONDCONST], tslop[MAXPIBONDCONST], tslop2[MAXPIBONDCONST];
        } bondpk;

        
struct t_electroneg {
        int nelecti, nelectj;
        int itype[200],ibond[200],iattach[200];
        int jtype[200],jbond[200],jattach[200];              
        float icorr[200],jcorr[200];
        }  electroneg;

#endif
