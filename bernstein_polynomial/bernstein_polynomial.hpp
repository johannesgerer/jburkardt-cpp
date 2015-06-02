double *bernstein_matrix ( int n );
double *bernstein_matrix_inverse ( int n );
double *bernstein_poly_01 ( int n, double x );
void bernstein_poly_01_values ( int *n_data, int *n, int *k, double *x, 
  double *b );
double *bernstein_poly_ab ( int n, double a, double b, double x );
double *bernstein_poly_ab_approx ( int n, double a, double b, double ydata[], 
  int nval, double xval[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8_choose ( int n, int k );
double r8_max ( double x, double y );
double r8_mop ( int i );
double r8_uniform_01 ( int *seed );
double r8mat_is_identity ( int n, double a[] );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mv_new ( int m, int n, double a[], double x[] );
double r8mat_norm_fro ( int m, int n, double a[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double *r8vec_linspace_new ( int n, double a_first, double a_last );
double r8vec_sum ( int n, double a[] );
void timestamp ( );
