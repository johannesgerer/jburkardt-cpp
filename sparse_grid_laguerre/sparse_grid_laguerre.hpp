int choose ( int n, int k );
void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
int i4_log_2 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
string i4_to_string ( int i4, string format );
int i4vec_product ( int n, int a[] );
void laguerre_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] );
double laguerre_integral_nd ( int dim_num, int expon[] );
void laguerre_weights ( int order, double weight[] );
void level_to_order_open ( int dim_num, int level[], int order[] );
double monomial_quadrature ( int dim_num, int expon[], int point_num, 
  double weight[], double x[] );
double *monomial_value ( int dim_num, int point_num, double x[], int expon[] );
int *multigrid_index_one ( int dim_num, int order_1d[], int order_nd );
double *product_weight_laguerre ( int dim_num, int order_1d[], int order_nd );
double r8_epsilon ( );
double r8_factorial ( int n );
double r8_huge ( );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
int s_len_trim ( string s );
void sparse_grid_laguerre ( int dim_num, int level_max, int point_num, 
  double grid_weight[], double grid_point[] );
void sparse_grid_laguerre_index ( int dim_num, int level_max, int point_num, 
  int grid_index [], int grid_base[] );
int sparse_grid_laguerre_size ( int dim_num, int level_max );
void timestamp ( );
void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more );
