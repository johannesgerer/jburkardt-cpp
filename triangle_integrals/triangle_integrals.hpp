void i4_to_pascal ( int k, int &i, int &j );
int i4_to_pascal_degree ( int k );
int pascal_to_i4 ( int i, int j );
double *poly_power ( int d1, double p1[], int n );
double *poly_power_linear ( int d1, double p1[], int n );
void poly_print ( int d, double p[], string title );
double *poly_product ( int d1, double p1[], int d2, double p2[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void rs_to_xy_map ( double t[], double &a, double &b, double &c, double &d, 
  double &e, double &f );
void timestamp ( );
double triangle_area ( double t[] );
double triangle_monomial_integral ( int i, int j, double t[] );
double triangle_poly_integral ( int d, double p[], double t[] );
double triangle_xy_integral ( double x1, double y1, double x2, double y2, 
  double x3, double y3 );
double triangle01_monomial_integral ( int i, int j );
double triangle01_poly_integral ( int d, double p[] );
int trinomial ( int i, int j, int k );
void xy_to_rs_map ( double t[], double &a, double &b, double &c, double &d, 
  double &e, double &f );
