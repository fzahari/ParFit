#ifndef CONTROL_H
#define CONTROL_H

struct t_control {
        int mode, nconfout,ffield,iprint,jprint,addconst;
        int optmet,electro,getcharge,iuv,pipl,gbsa;
        int inftype,outftype;
        char infile[128],outfile[128],addparfile[128],errorfile[128];
        }control;

#endif
