void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
double cube_monomial ( double a[], double b[], int expon[3] );
void cube_monomial_test ( int degree_max );
void cube_quad_test ( int degree_max );;
void cube_rule ( double a[], double b[], int order[], double w[], double xyz[] );
double cube_volume ( double a[], double b[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
bool i4vec_odd_any ( int n, int a[] );
int i4vec_product ( int n, int a[] );
void line_unit_o01 ( double w[], double x[] );
void line_unit_o02 ( double w[], double x[] );
void line_unit_o03 ( double w[], double x[] );
void line_unit_o04 ( double w[], double x[] );
void line_unit_o05 ( double w[], double x[] );
void line_unit_quad_test ( int degree_max );;
double *monomial_value ( int dim_num, int point_num, int expon[], double x[] );
double r8_abs ( double x );
double r8_choose ( int n, int k );
double r8_mop ( int i );
double r8mat_det_4d ( double a[4*4] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_direct_product ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double x[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void subcomp_next ( int n, int k, int a[], bool *more, int *h, int *t );
void timestamp ( );

