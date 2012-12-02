void chkder ( int m, int n, double x[], double fvec[], double fjac[], 
  int ldfjac, double xp[], double fvecp[], int mode, double err[] );
void dogleg ( int n, double r[], int lr, double diag[], double qtb[],
  double delta, double x[], double wa1[], double wa2[] );
double enorm ( int n, double x[] );
void fdjac1 ( void fcn ( int n, double x[], double f[], int *iflag ), 
  int n, double x[], double fvec[], double fjac[], int ldfjac, int *iflag,
  int ml, int mu, double epsfcn, double wa1[], double wa2[] );
void fdjac2 ( void fcn ( int m, int n, double x[], double fvec[], int *iflag ), 
  int m, int n, double x[], double fvec[], double fjac[], int ldfjac, 
  int *iflag, double epsfcn, double wa[] );
int hybrd ( void fcn ( int n, double x[], double fvec[], int *iflag ), 
  int n, double x[], double fvec[], double xtol, int maxfev, int ml, 
  int mu, double epsfcn, double diag[], int mode, double factor, int nprint, 
  int nfev, double fjac[], int ldfjac, double r[], int lr, double qtf[], 
  double wa1[], double wa2[], double wa3[], double wa4[] );
int hybrd1 ( void fcn ( int n, double x[], double fvec[], int *iflag ), int n, 
  double x[], double fvec[], double tol, double wa[], int lwa );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void qform ( int m, int n, double q[], int ldq, double wa[] );
void qrfac ( int m, int n, double a[], int lda, bool pivot, int ipvt[],
  int lipvt, double rdiag[], double acnorm[], double wa[] );
void r1mpyq ( int m, int n, double a[], int lda, double v[], double w[] );
bool r1updt ( int m, int n, double s[], int ls, double u[], double v[], double w[] );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_tiny ( );
double r8_uniform_01 ( int *seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );
