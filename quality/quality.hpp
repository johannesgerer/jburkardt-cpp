double alpha_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[] );
double arc_cosine ( double c );
double area_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[] );
void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m );
double beta_measure ( int dim_num, int n, double z[] );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
double chi_measure ( int dim_num, int n, double z[], int ns, 
  double *sample_routine ( int dim_num, int n, int *seed ), 
  int seed_init );
double d_measure ( int dim_num, int n, double z[], int ns, 
  double *sample_routine ( int dim_num, int n, int *seed ), 
  int seed_init );
double dge_det ( int n, double a[] );
int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2, 
  double x3, double y3 );
double *dtable_data_read ( char *input_filename, int m, int n );
void dtable_header_read ( char *input_filename, int *m, int *n );
int dtris2 ( int point_num, double point_xy[], int *tri_num, 
  int tri_vert[], int tri_nabe[] );
double e_measure ( int dim_num, int n, double z[], int ns, 
  double *sample_routine ( int dim_num, int n, int *seed ), 
  int seed_init );
int file_column_count ( char *input_filename );
int file_row_count ( char *input_filename );
void find_closest ( int dim_num, int n, int sample_num, double s[], double r[],
  int nearest[] );
double gamma_measure ( int dim_num, int n, double z[] );
double h_measure ( int dim_num, int n, double z[], int ns, 
  double *sample_routine ( int dim_num, int n, int *seed ), 
  int seed_init );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_sign ( int i );
int i4_wrap ( int ival, int ilo, int ihi );
int *i4vec_indicator ( int n );
void i4vec_print ( int n, int a[], char *title );
double lambda_measure ( int dim_num, int n, double z[] );
int lrline ( double xu, double yu, double xv1, double yv1, double xv2, 
  double yv2, double dv );
double mu_measure ( int dim_num, int n, double z[], int ns, 
  double *sample_routine ( int dim_num, int n, int *seed ), 
  int seed_init );
double nu_measure ( int dim_num, int n, double z[], int ns, 
  double *sample_routine ( int dim_num, int n, int *seed ), 
  int seed_init );
void perm_inv ( int n, int p[] );
double *pointset_spacing ( int dim_num, int n, double z[] );
double q_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[] );
double r0_measure ( int dim_num, int n, double z[] );
double r8_epsilon ( void );
double r8_huge ( void );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int *seed );
void r82vec_permute ( int n, double a[], int p[] );
int *r82vec_sort_heap_index_a ( int n, double a[] );
bool r8mat_in_01 ( int m, int n, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
double *r8mat_uniform_01 ( int m, int n, int *seed );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
double *r8vec_normal_01 ( int n, int *seed );
double *r8vec_uniform_01 ( int n, int *seed );
double *radius_maximus ( int dim_num, int n, double z[], bool walls );
int s_len_trim ( char *s );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
double *sample_hypercube_uniform ( int dim_num, int n, int *seed );
double *sample_sphere_uniform ( int m, int n, int *seed );
double sphere_measure ( int dim_num, int n, double z[] );
double sphere_volume_nd ( int dim_num, double r );
int swapec ( int i, int *top, int *btri, int *bedg, int point_num, 
  double point_xy[], int tri_num, int tri_vert[], int tri_nabe[], 
  int stack[] );
double tau_measure ( int dim_num, int n, double z[], int ns, 
  double *sample_routine ( int dim_num, int n, int *seed ), 
  int seed_init );
void timestamp ( void );
char *timestring ( void );
void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num, 
  int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg );

