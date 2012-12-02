double arc_sine ( double s );
double ball_f1_nd ( double func ( int n, double x[] ), int n, double center[], 
  double r );
double ball_f3_nd ( double func ( int n, double x[] ), int n, double center[], 
  double r );
double ball_monomial_nd ( int n, int p[], double r );
double ball_unit_07_3d ( double func ( double x, double y, double z ) );
double ball_unit_14_3d ( double func ( double x, double y, double z ) );
double ball_unit_15_3d ( double func ( double x, double y, double z ) );
double ball_unit_f1_nd ( double func ( int n, double x[] ), int n );
double ball_unit_f3_nd ( double func ( int n, double x[] ), int n );
double ball_unit_volume_3d ( );
double ball_unit_volume_nd ( int n );
double ball_volume_3d ( double r );
double ball_volume_nd ( int n, double r );
double c1_geg_monomial_integral ( double alpha, int expon );;
double c1_jac_monomial_integral ( double alpha, double beta, int expon );
double c1_leg_monomial_integral ( int expon );
double circle_annulus ( double func ( double x, double y ), double center[2], 
  double radius1, double radius2, int nr );
double circle_annulus_area_2d ( double radius1, double radius2 );
double circle_annulus_sector ( double func ( double x, double y ), 
  double center[2], double radius1, double radius2, double theta1, 
  double theta2, int nr );
double circle_annulus_sector_area_2d ( double radius1, double radius2, 
  double theta1, double theta2 );
double circle_area_2d ( double r );
double circle_cap_area_2d ( double r, double h );
double circle_cum ( double func ( double x, double y ), double center[2], 
  double radius, int order );
void circle_rt_set ( int rule, int nr, int nt, int nc, double ra[], 
  double rw[], double ta[], double tw[], double *cw );
void circle_rt_size ( int rule, int *nr, int *nt, int *nc );
double circle_rt_sum ( double func ( double x, double y ), double center[2], 
  double radius, int nr, double ra[], double rw[], int nt, double ta[], 
  double tw[], double zw );
double circle_sector ( double func ( double x, double y ), double center[2], 
  double radius, double theta1, double theta2, int nr );
double circle_sector_area_2d ( double r, double theta1, double theta2 );
double circle_triangle_area_2d ( double r, double theta1, double theta2 );
void circle_xy_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] );
int circle_xy_size ( int rule );
double circle_xy_sum ( double func ( double x, double y ), double center[2], 
  double r, int order, double xtab[], double ytab[], double weight[] );
void cn_geg_00_1 ( int n, double alpha, int o, double x[], double w[] );
int cn_geg_00_1_size ( int n, double alpha );
void cn_geg_01_1 ( int n, double alpha, int o, double x[], double w[] );
int cn_geg_01_1_size ( int n, double alpha );
void cn_geg_02_xiu ( int n, double alpha, int o, double x[], double w[] );
int cn_geg_02_xiu_size ( int n, double alpha );
void cn_geg_03_xiu ( int n, double alpha, int o, double x[], double w[] );
int cn_geg_03_xiu_size ( int n, double alpha );
double cn_geg_monomial_integral ( int n, double alpha, int expon[] );
void cn_jac_00_1 ( int n, double alpha, double beta, int o, double x[], 
  double w[] );
int cn_jac_00_1_size ( int n, double alpha, double beta );
void cn_jac_01_1 ( int n, double alpha, double beta, int o, double x[], 
  double w[] );
int cn_jac_01_1_size ( int n, double alpha, double beta );
double cn_jac_monomial_integral ( int n, double alpha, double beta, 
  int expon[] );
void cn_jac_02_xiu ( int n, double alpha, double beta, int o, double x[], 
  double w[] );
int cn_jac_02_xiu_size ( int n, double alpha, double beta );
void cn_leg_01_1 ( int n, int o, double x[], double w[] );
int cn_leg_01_1_size ( int n );
void cn_leg_02_xiu ( int n, int o, double x[], double w[] );
int cn_leg_02_xiu_size ( int n );
void cn_leg_03_1 ( int n, int o, double x[], double w[] );
int cn_leg_03_1_size ( int n );
void cn_leg_03_xiu ( int n, int o, double x[], double w[] );
int cn_leg_03_xiu_size ( int n );
void cn_leg_05_1 ( int n, int option, int o, double x[], double w[] );
int cn_leg_05_1_size ( int n );
void cn_leg_05_2 ( int n, int o, double x[], double w[] );
int cn_leg_05_2_size ( int n );
double cn_leg_monomial_integral ( int n, int expon[] );
double cone_unit_3d ( double func ( double x, double y, double z ) );
double cone_volume_3d ( double r, double h );
double cube_shell_nd ( double func ( int n, double x[] ), int n, double r1, 
  double r2 );
double cube_shell_volume_nd ( int n, double r1, double r2 );
double cube_unit_3d ( double func ( double x, double y, double z ) );
void cube_unit_nd ( double func ( int n, double x[] ), double qa[], 
  double qb[], int n, int k );
double cube_unit_volume_nd ( int n );
double ellipse_area_2d ( double r1, double r2 );
double ellipse_circumference_2d ( double r1, double r2 );
double ellipse_eccentricity_2d ( double r1, double r2 );
double ellipsoid_volume_3d ( double r1, double r2, double r3 );
int en_r2_01_1_size ( int n );
void en_r2_01_1 ( int n, int o, double x[], double w[] );
int en_r2_02_xiu_size ( int n );
void en_r2_02_xiu ( int n, int o, double x[], double w[] );
int en_r2_03_1_size ( int n );
void en_r2_03_1 ( int n, int o, double x[], double w[] );
int en_r2_03_2_size ( int n );
void en_r2_03_2 ( int n, int o, double x[], double w[] );
void en_r2_03_xiu ( int n, int o, double x[], double w[] );
int en_r2_03_xiu_size ( int n );
int en_r2_05_1_size ( int n );
void en_r2_05_1 ( int n, int option, int o, double x[], double w[] );
int en_r2_05_2_size ( int n );
void en_r2_05_2 ( int n, int o, double x[], double w[] );
int en_r2_05_3_size ( int n );
void en_r2_05_3 ( int n, int o, double x[], double w[] );
int en_r2_05_4_size ( int n );
void en_r2_05_4 ( int n, int o, double x[], double w[] );
int en_r2_05_5_size ( int n );
void en_r2_05_5 ( int n, int o, double x[], double w[] );
int en_r2_05_6_size ( int n );
void en_r2_05_6 ( int n, int o, double x[], double w[] );
int en_r2_07_1_size ( int n );
void en_r2_07_1 ( int n, int option, int o, double x[], double w[] );
int en_r2_07_2_size ( int n );
void en_r2_07_2 ( int n, int o, double x[], double w[] );
int en_r2_07_3_size ( int n );
void en_r2_07_3 ( int n, int option, int o, double x[], double w[] );
int en_r2_09_1_size ( int n );
void en_r2_09_1 ( int n, int option, int o, double x[], double w[] );
int en_r2_11_1_size ( int n );
void en_r2_11_1 ( int n, int option, int o, double x[], double w[] );
double en_r2_monomial_integral ( int n, int alpha[] );
double ep1_glg_monomial_integral ( int expon, double alpha );
double ep1_lag_monomial_integral ( int expon );
void epn_glg_00_1 ( int n, double alpha, int o, double x[], double w[] );
int epn_glg_00_1_size ( int n, double alpha );
void epn_glg_01_1 ( int n, double alpha, int o, double x[], double w[] );
int epn_glg_01_1_size ( int n, double alpha );
void epn_glg_02_xiu ( int n, double alpha, int o, double x[], double w[] );
int epn_glg_02_xiu_size ( int n, double alpha );
double epn_glg_monomial_integral ( int n, int expon[], double alpha );
void epn_lag_00_1 ( int n, int o, double x[], double w[] );
int epn_lag_00_1_size ( int n );
void epn_lag_01_1 ( int n, int o, double x[], double w[] );
int epn_lag_01_1_size ( int n );
void epn_lag_02_xiu ( int n, int o, double x[], double w[] );
int epn_lag_02_xiu_size ( int n );
double epn_lag_monomial_integral ( int n, int expon[] );
void gw_02_xiu ( int n, int o, double gamma0, double delta0, double c1, 
  double volume_1d, double x[], double w[] );
int gw_02_xiu_size ( int n );
double hexagon_area_2d ( double r );
double hexagon_sum ( double func ( double x, double y ), double center[2], 
  double r, int order, double xtab[], double ytab[], double weight[] );
double hexagon_unit_area_2d ( );
void hexagon_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] );
int hexagon_unit_size ( int rule );
int i4_factorial ( int n );
int i4_factorial2 ( int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
int i4vec_sum ( int n, int a[] );
void i4vec_zero ( int n, int a[] );
void ksub_next2 ( int n, int k, int a[], int *in, int *iout );
void legendre_set ( int order, double xtab[], double weight[] );
void legendre_set_x1 ( int order, double xtab[], double weight[] );
void legendre_set_x2 ( int order, double xtab[], double weight[] );
double lens_half_2d ( double func ( double x, double y ), double center[2], 
  double r, double theta1, double theta2, int order );
double lens_half_area_2d ( double r, double theta1, double theta2 );
double lens_half_h_area_2d ( double r, double h );
double lens_half_w_area_2d ( double r, double w );
double *monomial_value ( int dim_num, int point_num, double x[], int expon[] );
double octahedron_unit_nd ( double func ( int n, double x[] ), int n );
double octahedron_unit_volume_nd ( int n );
double parallelipiped_volume_3d ( double x[4], double y[4], double z[4] );
double parallelipiped_volume_nd ( int n, double v[] );
double polygon_1_2d ( int n, double x[], double y[] );
double polygon_x_2d ( int n, double x[], double y[] );
double polygon_xx_2d ( int n, double x[], double y[] );
double polygon_xy_2d ( int n, double x[], double y[] );
double polygon_y_2d ( int n, double x[], double y[] );
double polygon_yy_2d ( int n, double x[], double y[] );
double pyramid_unit_o01_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o05_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o06_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o08_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o08b_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o09_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o13_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o18_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o27_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_o48_3d ( double func ( double x, double y, double z ) );
double pyramid_unit_monomial_3d ( int alpha, int beta, int gamma );
double pyramid_unit_volume_3d ( );
double pyramid_volume_3d ( double r, double h );
double qmdpt ( double func ( int n, double x[] ), int n, int nsub );
double qmult_1d ( double func ( double x ), double a, double b );
double qmult_2d ( double func ( double x, double y ), double a, double b, 
  double fup ( double x ), double flo ( double x ) );
double qmult_3d ( double func ( double x, double y, double z ), double a, 
  double b, double fup1 ( double x ), double flo1 ( double x ), 
  double fup2 ( double x, double y ), double flo2 ( double x, double y ) );
double r8_abs ( double x );
double r8_choose ( int n, int k );
double r8_epsilon ( );
double r8_factorial ( int n );
double r8_gamma ( double x );
double r8_gamma_log ( double x );
double r8_huge ( );
double r8_hyper_2f1 ( double a, double b, double c, double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_mop ( int i );
double r8_psi ( double xx );
void r8_swap ( double *x, double *y );
void r8_swap3 ( double *x, double *y, double *z );
double r8_uniform_01 ( int *seed );
double r8ge_det ( int n, double a_lu[], int pivot[] );
int r8ge_fa ( int n, double a[], int pivot[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_even_select ( int n, double xlo, double xhi, int ival );
bool r8vec_mirror_next ( int n, double a[] );
void r8vec_print ( int n, double a[], char *title );
void r8vec_zero ( int n, double a[] );
double rectangle_3d ( double func ( double x, double y, double z ),
  double a[3], double b[3] );
double rectangle_sub_2d ( double func ( double x, double y ), double xval[2], 
  double yval[2], int nsub[2], int order, double xtab[], double ytab[], 
  double weight[] );
void rule_adjust ( double a, double b, double c, double d, int order, 
  double x[], double w[] );
int s_len_trim ( char *s );
double simplex_nd ( double func ( int n, double x[] ), int n, double v[] );
double simplex_unit_01_nd ( double func ( int n, double x[] ), int n );
double simplex_unit_03_nd ( double func ( int n, double x[] ), int n );
double simplex_unit_05_nd ( double func ( int n, double x[] ), int n );
double simplex_unit_05_2_nd ( double func ( int n, double x[] ), int n );
double simplex_unit_volume_nd ( int n );
double simplex_volume_nd ( int n, double v[] );
double sin_power_int ( double a, double b, int n );
double sphere_05_nd ( double func ( int n, double x[] ), int n, double center[], 
  double r );
double sphere_07_1_nd ( double func ( int n, double x[] ), int n, 
  double center[], double r );
double sphere_area_3d ( double r );
double sphere_area_nd ( int n, double r );
double sphere_cap_area_2d ( double r, double h );
double sphere_cap_area_3d ( double r, double h );
double sphere_cap_area_nd ( int dim_num, double r, double h );
double sphere_cap_volume_2d ( double r, double h );
double sphere_cap_volume_3d ( double r, double h );
double sphere_cap_volume_nd ( int dim_num, double r, double h );
double sphere_k ( int n );
double sphere_monomial_int_nd ( int n, double r, int e[] );
double sphere_shell_03_nd ( double func ( int n, double x[] ), int n, 
  double center[], double r1, double r2 );
double sphere_shell_volume_nd ( int n, double r1, double r2 );
double sphere_unit_03_nd ( double func ( int n, double x[] ), int n );
double sphere_unit_04_nd ( double func ( int n, double x[] ), int n );
double sphere_unit_05_nd ( double func ( int n, double x[] ), int n );
double sphere_unit_07_3d ( double func ( double x, double y, double z ) );
double sphere_unit_07_1_nd ( double func ( int n, double x[] ), int n );
double sphere_unit_07_2_nd ( double func ( int n, double x[] ), int n );
double sphere_unit_11_3d ( double func ( double x, double y, double z ) );
double sphere_unit_11_nd ( double func ( int n, double x[] ), int n );
double sphere_unit_14_3d ( double func ( double x, double y, double z ) );
double sphere_unit_15_3d ( double func ( double x, double y, double z ) );
double sphere_unit_area_3d ( );
double sphere_unit_area_nd ( int n );
void sphere_unit_area_values ( int *n_data, int *n, double *area );
double sphere_unit_monomial_nd ( int n, int p[] );
double sphere_unit_volume_nd ( int dim_num );
void sphere_unit_volume_values ( int *n_data, int *n, double *volume );
double sphere_volume_2d ( double r );
double sphere_volume_3d ( double r );
double sphere_volume_nd ( int dim_num, double r );
double square_sum ( double func ( double x, double y ), double center[2], 
  double r, int order, double xtab[], double ytab[], double weight[] );
void square_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] );
int square_unit_size ( int rule );
double square_unit_sum ( double func ( double x, double y ), int order, 
  double xtab[], double ytab[], double weight[] );
void subset_gray_next ( int n, int a[], bool *more, int *ncard, int *iadd );
double tetra_07 ( double func ( double x, double y, double z ), double x[], 
  double y[], double z[] );
double tetra_sum ( double func ( double x, double y, double z ), double x[4], 
  double y[4], double z[4], int order, double xtab[], double ytab[], 
  double ztab[], double weight[] );
double tetra_tproduct ( double func ( double x, double y, double z ), 
  int order, double x[4], double y[4], double z[4] );
int tetra_unit_size ( int rule );
void tetra_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double ztab[], double weight[] );
double tetra_unit_sum ( double func ( double x, double y, double z ), 
  int order, double xtab[], double ytab[], double ztab[], double weight[] );
double tetra_unit_volume ( );
double tetra_volume ( double x[4], double y[4], double z[4] );
void timestamp ( );
double torus_1 ( double func ( double x, double y, double z ), double r1, 
  double r2, int n );
double torus_14s ( double func ( double x, double y, double z ), double r1, 
  double r2 );
double torus_5s2 ( double func ( double x, double y, double z ), double r1, 
  double r2 );
double torus_6s2 ( double func ( double x, double y, double z ), double r1, 
  double r2 );
double torus_area_3d ( double r1, double r2 );
double torus_square_14c ( double func ( double x, double y, double z ), 
  double r1, double r2 );
double torus_square_5c2 ( double func ( double x, double y, double z ), 
  double r1, double r2 );
double torus_square_area_3d ( double r1, double r2 );
double torus_square_volume_3d ( double r1, double r2 );
double torus_volume_3d ( double r1, double r2 );
void triangle_rule_adjust ( double xval[3], double yval[3], int order, 
  double xtab[], double ytab[], double weight[], double xtab2[], double ytab2[], 
  double weight2[] );
double triangle_sub ( double func ( double x, double y ), double xval[3], 
  double yval[], int nsub, int order, double xtab[], double ytab[], 
  double weight[] );
double triangle_sum ( double func ( double x, double y ), 
  double xval[3], double yval[3], int order, double xtab[], double ytab[], 
  double weight[] );
double triangle_sum_adjusted ( double func ( double x, double y ), 
  int order, double xtab[], double ytab[], double weight[] );
void triangle_unit_product_set ( int rule, int order, double xtab[], 
  double ytab[], double weight[] );
int triangle_unit_product_size ( int rule );
void triangle_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] );
int triangle_unit_size ( int rule );
double triangle_unit_sum ( double func ( double x, double y ), int order, 
  double xtab[], double ytab[], double weight[] );
double triangle_unit_volume ( );
double triangle_volume ( double x[3], double y[3] );
double *tvec_even ( int nt );
double *tvec_even2 ( int nt );
double *tvec_even3 ( int nt );
double *tvec_even_bracket ( int nt, double theta1, double theta2 );
double *tvec_even_bracket2 ( int nt, double theta1, double theta2 );
double *tvec_even_bracket3 ( int nt, double theta1, double theta2 );
void vec_lex_next ( int dim_num, int base, int a[], bool *more );
