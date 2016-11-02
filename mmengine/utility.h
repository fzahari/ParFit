/* utility routines */
extern float *vector( int, int);
extern  int *ivector( int, int);
extern  long int *ilvector( int, int);
extern double *dvector( int, int);
extern float **matrix( int,  int,  int,  int);
extern  int **imatrix( int,  int,  int,  int);
extern double **dmatrix( int,  int,  int,  int);
extern void free_vector(float *,  int,  int);
extern void free_ivector( int *,  int,  int);
extern void free_ilvector( long int *,  int,  int);
extern void free_dvector(double *,  int,  int);
extern void free_matrix(float **, int, int, int, int);
extern void free_imatrix( int **, int, int, int, int);
extern void free_dmatrix(double **, int, int, int, int);
