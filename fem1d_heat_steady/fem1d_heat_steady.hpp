double *fem1d_heat_steady ( int n, double a, double b, double ua, double ub,
  double k ( double x ), double f ( double x ), double x[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4vec_zero_new ( int n );
double r8_abs ( double x );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_solve2 ( int n, double a[], double b[], int *ierror );
double *r8mat_zero_new ( int m, int n );
double *r8vec_even_new ( int n, double alo, double ahi );
double *r8vec_zero_new ( int n );
void timestamp ( );
