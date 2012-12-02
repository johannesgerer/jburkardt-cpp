double *cawiq ( int nt, double t[], int mlt[], int nwts, int ndx[], int key, 
  int nst, double aj[], double bj[], int *jdf, double zemu );
void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
  double wts[] );
double cegqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double f ( double x, int i ) );
double cegqfs ( int nt, int kind, double alpha, double beta, 
  double f ( double x, int i ) );
double ceiqf ( int nt, double t[], int mlt[], int kind, double alpha, 
  double beta, double a, double b, double f ( double x, int i ) );
double ceiqfs ( int nt, double t[], int mlt[], int kind, double alpha, 
  double beta, double f ( double x, int i ) );
void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double t[], double wts[] );
void cgqfs ( int nt, int kind, double alpha, double beta, int lo, double t[], 
  double wts[] );
void chkqf ( double t[], double wts[], int mlt[], int nt, int nwts, int ndx[], 
  int key, int mop, int mex, int kind, double alpha, double beta, int lo, 
  double a, double b );
void chkqfs ( double t[], double wts[], int mlt[], int nt, int nwts, int ndx[], 
  int key, double w[], int mop, int mex, int kind, double alpha, double beta, 
  int lo );
double *ciqf ( int nt, double t[], int mlt[], int nwts, int ndx[], int key, 
  int kind, double alpha, double beta, double a, double b, int lo );
double *ciqfs ( int nt, double t[], int mlt[], int nwts, int ndx[], int key, 
  int kind, double alpha, double beta, int lo );
double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
  double bj[] );
double *cliqf ( int nt, double t[], int kind, double alpha, double beta, 
  double a, double b, int lo );
double *cliqfs ( int nt, double t[], int kind, double alpha, double beta, 
  int lo );
double *cwiqd ( int m, int nm, int l, double v, double xk[], int nstar, 
  double phi[], double a[], double r[] );
double eiqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  int key, double f ( double x, int i ) );
double eiqfs ( int nt, double t[], double wts[], double f ( double x, int i ) );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_sign ( int i );
void imtqlx ( int n, double d[], double e[], double z[] );
void parchk ( int kind, int m, double alpha, double beta );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_sign ( double x );
void r8vec_print ( int n, double a[], string title );
double *scmm ( int m, int kind, double alpha, double beta, double a, double b );
void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  double swts[], double st[], int kind, double alpha, double beta, double a, 
  double b );
double *sct ( int nt, double t[], int kind, double a, double b );
void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
  double wts[] );
void timestamp ( );
double *wm ( int m, int kind, double alpha, double beta );
double *wtfn ( double t[], int nt, int kind, double alpha, double beta );
