int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void imtqlx ( int n, double d[], double e[], double z[] );
double *l_exponential_product ( int p, double b );
double l_integral ( int n );
double *l_polynomial ( int m, int n, double x[] );
double *l_polynomial_coefficients ( int n );
void l_polynomial_values ( int &n_data, int &n, double &x, double &fx );
double *l_polynomial_zeros ( int n );
double *l_power_product ( int p, int e );
void l_quadrature_rule ( int n, double x[], double w[] );
double *lf_function ( int mm, int n, double alpha, double x[] );
void lf_function_values ( int &n_data, int &n, double &a, double &x,
  double &fx );
double *lf_function_zeros ( int n, double alpha );
double lf_integral ( int n, double alpha );
void lf_quadrature_rule ( int n, double alpha, double x[], double w[] );
double lm_integral ( int n, int m );
double *lm_polynomial ( int mm, int n, int m, double x[] );
double *lm_polynomial_coefficients ( int n, int m );
void lm_polynomial_values ( int &n_data, int &n, int &m, double &x,
  double &fx );
double *lm_polynomial_zeros ( int n, int m );
void lm_quadrature_rule ( int n, int m, double x[], double w[] );
double r8_epsilon ( );
double r8_factorial ( int n );
double r8_gamma ( double x );
double r8_sign ( double x );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void r8vec_print ( int n, double a[], string title );
void r8vec2_print ( int n, double a1[], double a2[], string title );
void timestamp ( );
