#ifndef PIDATA_H
#define PIDATA_H

struct t_penint {
        float **pen, *ed;
        }       penint;
struct t_ed123 {
        int nscft;
        float *ed1, *ed2, *ed3, *ed4, *ed5;
        }       ed123;
struct t_densit {
        float **den, **dena, **denb;
        }       densit;
struct t_local {
        float *da, *adnmo, *adnmt, *d12, *d13, *fn;
        float *f22, *f27, *f29, *f1, *fa1;
        }  local;
struct t_uvnpi {
        int ifirst, kfirst;
        }       uvnpi;
struct t_pigroup {
        int nfrag;
        int fraglist[MAXORB],grlist[MAXORB],nao[MAXORB/4],nelecs[MAXORB/4];
        } pigroup;       

#endif
