int  make_atom(int , float , float , float,char * );
void make_bond(int , int , int );
void generate_mm2lists(int *, int *, int **, int **);
void gettoken(void);
void mm3type(void);
void message_alert(char *, char *);
void set_atomdata(int,int,int,int,int,int);
int mm3_mmxtype(int);
void setbondorder(int,int);
void set_atomtypes(int);
void hdel(int);
void hadd(void);
FILE * fopen_path ( char * , char * , char * ) ;
int ReadMMInt(char *,int,int);
double ReadMMFloat(char *,int,int);
void FetchRecord(FILE *, char *);
void write_mm3(void);
void mm32mod(int);

