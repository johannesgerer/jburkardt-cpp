int *adj_set_q4_mesh ( int node_num, int element_num,
  int element_node[], int element_neighbor[], int adj_num, int adj_col[] );
int adj_size_q4_mesh ( int node_num, int element_num, int element_node[], 
  int element_neighbor[], int adj_col[] );
void area_q4_mesh ( int node_num, int element_num, double node_xy[], 
  int element_node[], double element_area[], double *mesh_area );
double area_quad ( double quad_xy[2*4] );
void bandwidth ( int element_order, int element_num, int element_node[], 
  int *ml, int *mu, int *m );
int boundary_edge_count_q4_mesh ( int element_num, int element_node[] );
int boundary_edge_count_euler_q4_mesh ( int node_num, int element_num, 
  int hole_num );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void example1_q4_mesh ( int node_num, int element_num, double node_xy[], 
  int element_node[], int element_neighbor[] );
void example1_q4_mesh_size ( int *node_num, int *element_num, int *hole_num );
void example2_q4_mesh ( int node_num, int element_num, double node_xy[], 
  int element_node[], int element_neighbor[] );
void example2_q4_mesh_size ( int *node_num, int *element_num, int *hole_num );
int file_column_count ( string input_filename );
int file_row_count ( string input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
int i4col_sorted_unique_count ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
void i4mat_copy ( int m, int n, int a1[], int a2[] );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
int i4row_compare ( int m, int n, int a[], int i, int j );
void i4row_sort_a ( int m, int n, int a[] );
void i4row_swap ( int m, int n, int a[], int irow1, int irow2 );
void i4vec_heap_d ( int n, int a[] );
void i4vec_print ( int n, int a[], string title );
void i4vec_sort_heap_a ( int n, int a[] );
void mesh_base_zero ( int node_num, int element_order, 
  int element_num, int element_node[] );
int *neighbor_elements_q4_mesh ( int element_num, int element_node[] );
int *node_order_q4_mesh ( int element_num, int element_node[], int node_num );
void plot_q4_mesh ( int node_num, int element_num, double node_xy[], 
  int element_node[], int node_show, int element_show, string output_filename );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_huge ( );
int r8_nint ( double x );
double r8_uniform_01 ( int *seed );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_bracket ( int n, double x[], double xval, int *left, 
  int *right );
void r8vec_print ( int n, double a[], string title );
double r8vec_sum ( int n, double a[] );
void reference_to_physical_q4 ( double q4[2*4], int n, double rs[], 
  double xy[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void sample_q4_mesh ( int node_num, double node_xy[], int element_num, 
  int element_node[], int sample_num, int *seed, double sample_xy[], 
  int sample_element[] );
void sample_quad ( double quad_xy[2*4], int n, int *seed, double xy[] );
double *sample_quad_new ( double quad_xy[2*4], int n, int *seed );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void timestamp ( );
double triangle_area ( double t[2*3] );
void triangle_sample ( double t[2*3], int n, int *seed, double p[] );
