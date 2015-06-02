void jacobi_eigenvalue ( int n, double a[], int it_max, double v[], 
  double d[], int &it_num, int &rot_num );
void moment_method ( int n, double moment[], double x[], double w[] );
double *moments_laguerre ( int m );
double *moments_legendre ( int m, double a, double b );
double *moments_normal_01 ( int m );
double *moments_normal ( int m, double mu, double sigma );
double *moments_truncated_normal_ab ( int m, double mu, double sigma, 
  double a, double b );
double *moments_truncated_normal_a ( int m, double mu, double sigma, 
  double a );
double *moments_truncated_normal_b ( int m, double mu, double sigma, 
  double b );
double normal_01_cdf ( double x );
double normal_01_pdf ( double x );
double r8_choose ( int n, int k );
double r8_factorial ( int n );
double r8_factorial2 ( int n );
double r8_mop ( int i );
double *r8mat_cholesky_factor_upper ( int n, double a[], int &flag );
double *r8mat_copy_new ( int m, int n, double a1[] );
void r8mat_diag_get_vector ( int n, double a[], double v[] );
void r8mat_identity  ( int n, double a[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8vec_print_dupe ( int n, double a[], string title );
void r8vec2_print ( int n, double a1[], double a2[], string title );
void timestamp_dupe ( );
double truncated_normal_ab_moment ( int order, double mu, double s, double a,
  double b );
double truncated_normal_a_moment ( int order, double mu, double s, double a );
double truncated_normal_b_moment ( int order, double mu, double s, double b );
