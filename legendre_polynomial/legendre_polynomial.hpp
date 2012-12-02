int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void imtqlx ( int n, double d[], double e[], double z[] );
double *p_exponential_product ( int p, double b );
double p_integral ( int n );
double *p_polynomial ( int m, int n, double x[] );
double *p_polynomial_coefficients ( int n );
double *p_polynomial_prime ( int m, int n, double x[] );
void p_polynomial_values ( int &n_data, int &n, double &x, double &fx );
double *p_polynomial_zeros ( int nt );
double *p_power_product ( int p, int e );
void p_quadrature_rule ( int nt, double t[], double wts[] );
double *pm_polynomial ( int mm, int n, int m, double x[] );
void pm_polynomial_values ( int &n_data, int &n, int &m, double &x, double &fx );
double *pmn_polynomial ( int mm, int n, int m, double x[] );
void pmn_polynomial_values ( int &n_data, int &n, int &m, double &x, 
  double &fx );
double *pmns_polynomial ( int mm, int n, int m, double x[] );
void pmns_polynomial_values ( int &n_data, int &n, int &m, double &x, 
  double &fx );
double *pn_polynomial ( int m, int n, double x[] );
double *pn_pair_product ( int p );
double r8_abs ( double x );
double r8_add ( double x, double y );
double r8_epsilon ( );
double r8_factorial ( int n );
double r8_sign ( double x );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void r8vec_print ( int n, double a[], string title );
void r8vec2_print ( int n, double a1[], double a2[], string title );
void timestamp ( );
