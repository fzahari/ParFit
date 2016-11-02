#ifndef EHEAT_P
#define EHEAT_P

struct t_ehpara {
        int   neheat, nevhf;
        int nn[220];
        float ee[220], estr[220], et[220];   // bond data
        char cc[220][21], mm[155][21];
        int   cent[155];
        char  kheat[155][13];
        float ss[155], ssless[155];          // structural data
        } ehpara;

struct t_epiheat {
        int  npihf;
        long int nnp[125];
        float eep[125], aa[125], bb[125], ccc[125];
        char ccp[125][21];} epiheat;

struct t_piheat {
        float tepi, tesig, ecpi, etpi;
        }       piheat;

struct t_elowht {
        float epi;
        int mult;
        int uhf;
        }       elowht;

#endif

double dihdrl(int,int,int,int);
void numeral(int, char *, int);
void eheatpi(float *,float *,float *,float *,int *);
float get_heat_of_formation(void);
void eheat(void);
