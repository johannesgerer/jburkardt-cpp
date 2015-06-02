void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
void line_unit_o01 ( double w[], double x[] );
void line_unit_o02 ( double w[], double x[] );
void line_unit_o03 ( double w[], double x[] );
void line_unit_o04 ( double w[], double x[] );
void line_unit_o05 ( double w[], double x[] );
double *monomial_value ( int dim_num, int point_num, int expon[], double x[] );
double square_monomial ( double a[], double b[], int expon[2] );
void square_monomial_test ( int degree_max );
void square_quad_test ( int degree_max );;
void square_rule ( double a[], double b[], int order[], double w[], double xy[] );
double square_volume ( double a[], double b[] );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_direct_product ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double x[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void subcomp_next ( int n, int k, int a[], bool *more, int *h, int *t );
void timestamp ( );

