int *abscissa_level_open_nd ( int level_max, int dim_num, int test_num, 
  int test_val[] );
void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
double f2_abscissa ( int order, int i );
void gl_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] );
double gp_abscissa ( int order, int i );
int i4_choose ( int n, int k );
int i4_log_2 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
string i4_to_string ( int i4, string format );
int i4vec_product ( int n, int a[] );
int index_to_level_open ( int dim_num, int t[], int order, int level_max );
void level_to_order_open ( int dim_num, int level[], int order[] );
int *levels_open_index ( int dim_num, int level_max, int point_num );
int *multigrid_index1 ( int dim_num, int order_1d[], int order_nd );
void multigrid_scale_open ( int dim_num, int order_nd, int level_max, 
  int level_1d[], int grid_index[] );
double nco_abscissa ( int order, int i );
double r8_epsilon ( );
double r8_huge ( );
void r8mat_write ( string output_filename, int m, int n, double table[] );
int sparse_grid_f2s_size ( int dim_num, int level_max );
int sparse_grid_gps_size ( int dim_num, int level_max );
int sparse_grid_ofn_size ( int dim_num, int level_max );
int sparse_grid_onn_size ( int dim_num, int level_max );
int sparse_grid_own_size ( int dim_num, int level_max );
void timestamp ( );
void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more );
