void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
void gl_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] );
void gl_weights ( int order, double weight[] );
int i4_choose ( int n, int k );
int i4_log_2 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
string i4_to_string ( int i4, string format );
int i4vec_product ( int n, int a[] );
int *index_level_gl ( int level, int level_max, int dim_num, int point_num, 
  int grid_index[], int grid_base[] );
void level_to_order_open ( int dim_num, int level[], int order[] );
double monomial_int01 ( int dim_num, int expon[] );
double monomial_quadrature ( int dim_num, int expon[], int point_num, 
  double weight[], double x[] );
double *monomial_value ( int dim_num, int point_num, double x[], int expon[] );
int *multigrid_index_z ( int dim_num, int order_1d[], int order_nd );
double *product_weight_gl ( int dim_num, int order_1d[], int order_nd );
double r8_epsilon ( );
double r8_huge ( );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
int s_len_trim ( char *s );
void sparse_grid_gl ( int dim_num, int level_max, int point_num, 
  double grid_weight[], double grid_point[] );
void sparse_grid_gl_index ( int dim_num, int level_max, int point_num, 
  int grid_index [], int grid_base[] );
int sparse_grid_gl_size ( int dim_num, int level_max );
void timestamp ( );
void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more );
