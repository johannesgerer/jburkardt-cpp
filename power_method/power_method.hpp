double cpu_time ( void );
double *fibonacci2 ( int n );
void power_method ( int n, double a[], double y[], int it_max, double tol,
  double *lambda, int *it_num );
void power_method2 ( int n, double a[], double x_init[], int it_max, 
  double tol, complex <double> *lambda, complex <double> v[], int *it_num );
double r8_abs ( double x );
double r8_epsilon ( void );
void r8mat_mv ( int m, int n, double a[], double x[], double ax[] );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_divide ( int n, double a[], double s );
double r8vec_dot ( int n, double a1[], double a2[] );
double r8vec_norm_l2 ( int n, double a[] );
double *r8vec_uniform_01 ( int n, int *seed );
void timestamp ( void );
double *tris ( int m, int n, double x, double y, double z );
complex <double> *tris_eigenvalues ( int n, double x, double y, double z );
