int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8_abs ( double x );
double r8_normal_01 ( int *seed );
double r8_sign ( double x );
double r8_uniform_01 ( int *seed );
void r8bin_print ( int bin_num, int bin[], double bin_limit[], string title );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
double *r8mat_house_axh_new ( int n, double a[], double v[] );
void r8mat_identity ( int n, double a[] );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
void r8mat_orth_uniform ( int n, int *seed, double a[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8symm_test ( int n, double lambda_mean, double lambda_dev, int *seed, 
  double a[], double q[], double lambda[] );
void r8vec_bin ( int n, double x[], int bin_num, double bin_min, double bin_max,
  int bin[], double bin_limit[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double *r8vec_house_column ( int n, double a[], int k );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
double r8vec_norm_l2 ( int n, double a[] );
void r8vec_normal ( int n, double b, double c, int *seed, double x[] );
void r8vec_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int *seed );
double *r8vec_zero_new ( int n );
void r8vec2_print ( int n, double a1[], double a2[], string title );
void timestamp ( );
