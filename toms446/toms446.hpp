double *binom ( double x[], double xx[], int npl, int m, int nt );
double *cheby ( int nf, int npl, double *functn ( double x ) );
double *dfrnt ( double xx[], int npl );
double echeb ( double x, double coef[], int npl );
double edcheb ( double x, double coef[], int npl );
double *invert ( double x[], double xx[], int npl, int net );
void mltply ( double xx[], double x2[], int npl, double x3[] );
double *mltply_new ( double xx[], double x2[], int npl );
double *ntgrt ( double xx[], int npl );
double r8_abs ( double x );
void r8vec_add ( int n, double a1[], double a2[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double *r8vec_copy_new ( int n, double a1[] );
void r8vec_scale ( double s, int n, double a[] );
double *r8vec_zero_new ( int n );
void timestamp ( );
void xalfa2 ( double x[], double xx[], int npl, int m, int maxet, double epsln, 
  int &net );
void xalfa3 ( double x[], double xx[], int npl, int m, int maxet, double epsln, 
  int &net );
