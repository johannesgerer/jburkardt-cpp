int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
string i4_to_string ( int i4 );
void imtqlx ( int n, double d[], double e[], double z[] );
double j_double_product_integral ( int i, int j, double a, double b );
double j_integral ( int n );
double *j_polynomial ( int m, int n, double alpha, double beta, double x[] );
void j_polynomial_values ( int &n_data, int &n, double &a, double &b, double &x,
  double &fx );
double *j_polynomial_zeros ( int n, double alpha, double beta );
void j_quadrature_rule ( int n, double alpha, double beta, double x[], 
  double w[] );
double r8_choose ( int n, int k );
double r8_epsilon ( );
double r8_factorial ( int n );
double r8_max ( double x, double y );
double r8_sign ( double x );
string r8_to_string ( double r8, string format );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double *r8vec_linspace_new ( int n, double a_first, double a_last );
void r8vec_print ( int n, double a[], string title );
void r8vec2_print ( int n, double a1[], double a2[], string title );
void timestamp ( );
