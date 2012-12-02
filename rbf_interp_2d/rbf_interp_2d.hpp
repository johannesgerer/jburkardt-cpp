void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
double dnrm2 ( int n, double x[], int incx );
void drot ( int n, double x[], int incx, double y[], int incy, double c,
  double s );
void drotg ( double *sa, double *sb, double *c, double *s );
void dscal ( int n, double sa, double x[], int incx );
int dsvdc ( double a[], int lda, int m, int n, double s[], double e[], 
  double u[], int ldu, double v[], int ldv, double work[], int job );
void dswap ( int n, double x[], int incx, double y[], int incy );
void phi1 ( int n, double r[], double r0, double v[] );
void phi2 ( int n, double r[], double r0, double v[] );
void phi3 ( int n, double r[], double r0, double v[] );
void phi4 ( int n, double r[], double r0, double v[] );
double *r8mat_solve_svd ( int m, int n, double a[], double b[] );
double *rbf_interp ( int m, int nd, double xd[], double r0, 
  void phi ( int n, double r[], double r0, double v[] ), double w[], 
  int ni, double xi[] );
double *rbf_weight ( int m, int nd, double xd[], double r0, 
  void phi ( int n, double r[], double r0, double v[] ), 
  double fd[] );
