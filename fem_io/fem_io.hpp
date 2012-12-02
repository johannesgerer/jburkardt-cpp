char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
void fem_data_read ( string node_coord_file_name, string element_file_name, 
  string node_data_file_name, int dim_num, int node_num, int element_num,
  int element_order, int node_data_num, double **node_coord, int **element_node,
  double **node_data );
void fem_header_print ( int dim_num, int node_num, int element_num, 
  int element_order, int node_data_num );
void fem_header_read ( string node_coord_file_name, string element_file_name, 
  string node_data_file_name, int *dim_num, int *node_num, int *element_num,
  int *element_order, int *node_data_num );
void fem_write ( string node_coord_file_name, string element_file_name, 
  string node_data_file_name, int dim_num, int node_num, int element_num, 
  int element_order, int node_data_num, double node_coord[], 
  int element_node[], double node_data[] );
int file_column_count ( string input_filename );
int file_row_count ( string input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
double r8_epsilon ( void );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( void );
