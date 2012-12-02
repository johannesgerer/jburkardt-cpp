double *cc_compute_points ( int n );
int i4_power ( int i, int j );
int i4vec_product ( int n, int a[] );
double *lagrange_base_1d ( int nd, double xd[], int ni, double xi[] );
double *lagrange_interp_nd_grid ( int m, int n_1d[], double a[], double b[], int nd );
double *lagrange_interp_nd_grid2 ( int m, int ind[], double a[], double b[], int nd );
int lagrange_interp_nd_size ( int m, int ind[] );
int lagrange_interp_nd_size2 ( int m, int ind[] );
double *lagrange_interp_nd_value ( int m, int n_1d[], double a[], double b[], int nd, 
  double zd[], int ni, double xi[] );
double *lagrange_interp_nd_value2 ( int m, int ind[], double a[], double b[], int nd, 
  double zd[], int ni, double xi[] );
int order_from_level_135 ( int l );
double r8_abs ( double x );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
void r8vec_direct_product ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double x[] );
void r8vec_direct_product2 ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double w[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_norm_affine ( int n, double v0[], double v1[] );
void timestamp ( );
