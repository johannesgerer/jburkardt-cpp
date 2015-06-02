int *abscissa_level_closed_nd ( int level_max, int dim_num, int test_num, 
  int test_val[] );
double cc_abscissa ( int order, int i );
double *cc_weights ( int n );
void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
int i4_choose ( int n, int k );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_mop ( int i );
int i4_power ( int i, int j );
string i4_to_string ( int i4, string format );
int i4vec_product ( int n, int a[] );
int index_to_level_closed ( int dim_num, int t[], int order, int level_max );
void level_to_order_closed ( int dim_num, int level[], int order[] );
double monomial_int01 ( int dim_num, int expon[] );
double monomial_quadrature ( int dim_num, int expon[], int point_num, 
  double w[], double x[] );
double *monomial_value ( int dim_num, int point_num, double x[], int expon[] );
int *multigrid_index0 ( int dim_num, int order_1d[], int order_nd );
void multigrid_scale_closed ( int dim_num, int order_nd, int level_max, 
  int level_1d[], int grid_index[] );
double *product_weights_cc ( int dim_num, int order_1d[], int order_nd );
double r8_epsilon ( );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
void sparse_grid_cc ( int dim_num, int level_max, int point_num, 
  double grid_weight[], double grid_point[] );
int *sparse_grid_cc_index ( int dim_num, int level_max, int point_num );
int sparse_grid_cfn_size ( int dim_num, int level_max );
int sparse_grid_cc_size_old ( int dim_num, int level_max );
void sparse_grid_cc_weights ( int dim_num, int level_max, int point_num, 
  int grid_index[], double grid_weight[] );
int sparse_grid_ccs_size ( int dim_num, int level_max );
void timestamp ( );
void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more );

