void set_atomtype(int,int,int,int,int,int);
void set_atomtypes(int);
int is_ring31(int);
int is_ring41(int);
int is_ring51(int);
int is_cyclo5(int, int *);
int isbond(int,int);
int aromatic_5(int,int *);
int find_rsize(int,int);
void get_rsize(int,int,int, int *);
int icompare(int, int *, int *);
void adjust_mmfftypes();
void deletebond(int,int);
void pireset(void);

