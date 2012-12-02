int *abscissa_level_closed_nd ( int level_max, int dim_num, int test_num, 
  int test_val[] );
int *abscissa_level_open_nd ( int level_max, int dim_num, int test_num, 
  int test_val[] );
double cc_abscissa ( int order, int i );
double *cc_weights ( int n );
void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
double f1_abscissa ( int order, int i );
double *f1_weights ( int order );
double f2_abscissa ( int order, int i );
double *f2_weights ( int order );
void gh_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] );
double *gh_weights ( int order );
void gl_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] );
double *gl_weights ( int order );
double gp_abscissa ( int order, int i );
double *gp_weights ( int order );
int i4_log_2 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
int i4vec_product ( int n, int a[] );
int *index_level_own ( int level, int level_max, int dim_num, int point_num, 
  int grid_index[], int grid_base[] );
int index_to_level_closed ( int dim_num, int t[], int order, int level_max );
int index_to_level_open ( int dim_num, int t[], int order, int level_max );
void level_to_order_closed ( int dim_num, int level[], int order[] );
void level_to_order_open ( int dim_num, int level[], int order[] );
void levels_index ( int dim_num, int level_max, int rule, int point_num, 
  int grid_index[], int grid_base[] );
void levels_index_cfn ( int dim_num, int level_max, int point_num, 
  int grid_index[], int grid_base[] );
void levels_index_ofn ( int dim_num, int level_max, int point_num, 
  int grid_index[], int grid_base[] );
void levels_index_onn ( int dim_num, int level_max, int point_num, 
  int grid_index [], int grid_base[] );
void levels_index_own ( int dim_num, int level_max, int point_num, 
  int grid_index [], int grid_base[] );
int levels_index_size ( int dim_num, int level_max, int rule );
int levels_index_size_onn ( int dim_num, int level_max );
int levels_index_size_own ( int dim_num, int level_max );
void lg_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] );
double *lg_weights ( int order );
double monomial_integral_hermite ( int dim_num, int expon[] );
double monomial_integral_laguerre ( int dim_num, int expon[] );
double monomial_integral_legendre ( int dim_num, int expon[] );
double monomial_quadrature ( int dim_num, int expon[], int point_num, 
  double weight[], double x[], int rule );
double *monomial_value ( int dim_num, int point_num, double x[], int expon[] );
int *multigrid_index_cfn ( int dim_num, int order_1d[], int order_nd );
int *multigrid_index_ofn ( int dim_num, int order_1d[], int order_nd );
int *multigrid_index_onn ( int dim_num, int order_1d[], int order_nd );
int *multigrid_index_own ( int dim_num, int order_1d[], int order_nd );
void multigrid_scale_closed ( int dim_num, int order_nd, int level_max, 
  int level_1d[], int grid_index[] );
void multigrid_scale_open ( int dim_num, int order_nd, int level_max, 
  int level_1d[], int grid_index[] );
double *product_weights ( int dim_num, int order_1d[], int order_nd, int rule );
double r8_abs ( double x );
double r8_choose ( int n, int k );
double r8_factorial ( int n );
double r8_factorial2 ( int n );
double r8_huge ( );
double r8_mop ( int i );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
void sparse_grid ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] );
int sparse_grid_cc_size ( int dim_num, int level_max );
void sparse_grid_cfn ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] );
void sparse_grid_ofn ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] );
int sparse_grid_ofn_size ( int dim_num, int level_max );
void sparse_grid_onn ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] );
void sparse_grid_own ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] );
void sparse_grid_weights_cfn ( int dim_num, int level_max, int rule, 
  int point_num, int grid_index[], double grid_weight[] );
void sparse_grid_weights_ofn ( int dim_num, int level_max, int rule, 
  int point_num, int grid_index[], double grid_weight[] );
void timestamp ( );
void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more );

