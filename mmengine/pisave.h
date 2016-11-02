#ifndef PISAVE_H
#define PISAVE_H

struct t_pibondsave {
        int *ibsave, *kbsave, *nbsave;
        float *bksave, *blsave, *pksave, *plsave, *plsave2;
        }       pibondsave;

struct t_pitorsave {
        int *itsave, *ktsave, *ncsave;
        float    *tcsave;
        }       pitorsave;

struct t_pisave {
        float **fjsave, **fksave, **gjsave,**gksave;
        }       pisave;

struct t_mtape {
        float **vmtape, **vbtape, **pensave;
        }       mtape;

#endif
