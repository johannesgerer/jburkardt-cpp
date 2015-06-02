double *fem2d_bvp_quadratic ( int nx, int ny, double a ( double x, double y ), 
  double c ( double x, double y ), double f ( double x, double y ), 
  double x[], double y[] );
double fem2d_h1s_error_quadratic ( int nx, int ny, double x[], double y[], 
  double u[], double exact_ux ( double x, double y ), 
  double exact_uy ( double x, double y ) );
double fem2d_l1_error ( int nx, int ny, double x[], double y[], double u[], 
  double exact ( double x, double y ) );
double fem2d_l2_error_quadratic ( int nx, int ny, double x[], double y[], 
  double u[], double exact ( double x, double y ) );
int *i4vec_zero_new ( int n );
double *r8mat_solve2 ( int n, double a[], double b[], int &ierror );
double *r8mat_zero_new ( int m, int n );
double *r8vec_even_new ( int n, double alo, double ahi );
double *r8vec_zero_new ( int n );
void timestamp ( );

