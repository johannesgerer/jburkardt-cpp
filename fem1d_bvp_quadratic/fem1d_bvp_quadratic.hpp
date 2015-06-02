double *fem1d_bvp_quadratic ( int n, double a ( double x ), double c ( double x ), 
  double f ( double x ), double x[] );
double h1s_error_quadratic ( int n, double x[], double u[], 
  double exact_ux ( double x ) );
int *i4vec_zero_new ( int n );
double l1_error ( int n, double x[], double u[], double exact ( double x ) );
double l2_error_quadratic ( int n, double x[], double u[], 
  double exact ( double x ) );
double *r8mat_solve2 ( int n, double a[], double b[], int *ierror );
double *r8mat_zero_new ( int m, int n );
double *r8vec_even_new ( int n, double alo, double ahi );
double *r8vec_zero_new ( int n );
void timestamp ( );
