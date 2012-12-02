double *combin ( double alpha, double beta, int n );
double *combin_inverse ( double alpha, double beta, int n );
double condition_hager ( int n, double a[] );
double condition_linpack ( int n, double a[] );
double condition_sample1 ( int n, double a[], int m );
double *conex1 ( double alpha );
double *conex1_inverse ( double alpha );
double *conex2 ( double alpha );
double *conex2_inverse ( double alpha );
double *conex3 ( int n );
double *conex3_inverse ( int n );
double *conex4 ( );
double *conex4_inverse ( );
double *kahan ( double alpha, int m, int n );
double *kahan_inverse ( double alpha, int n );
int r8ge_fa ( int n, double a[], int pivot[] );
double *r8ge_inverse ( int n, double a[], int pivot[] );
void r8ge_sl ( int n, double a_lu[], int pivot[], double b[], int job );

