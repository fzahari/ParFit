#ifndef ATOM_K_H
#define ATOM_K_H

struct t_atom_k {
        int natomtype;
        int type[MAXATOMTYPE], valency[MAXATOMTYPE];
        int tclass[MAXATOMTYPE], tclass1[MAXATOMTYPE], tclass2[MAXATOMTYPE];
        char symbol[MAXATOMTYPE][5], description[MAXATOMTYPE][21];
        int  number[MAXATOMTYPE];
        float weight[MAXATOMTYPE];
        int   ligands[MAXATOMTYPE];
        }  atom_k;

#endif
