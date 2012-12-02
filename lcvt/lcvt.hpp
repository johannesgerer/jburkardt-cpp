char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
double cluster_energy ( int dim_num, int n, double cell_generator[],
  int sample_num_cvt, int sample_function_cvt, int *seed );
void cvt ( int m, int n, int sample_function_init, 
  int sample_function_cvt, int sample_num_cvt, int maxit, int *seed, 
  double generator[] );
void cvt_iteration ( int m, int n, double generator[], int sample_num_cvt,
 int sample_function_cvt, int *seed, double *change_l2 );
void cvt_write ( int dim_num, int n, int batch, int seed_init, int seed, 
  char *init_string, int it_max, int it_fixed, int it_num, 
  double it_diff, double energy, char *sample_string, int sample_num, 
  double r[], char *file_out_name, bool comment );
int file_column_count ( char *input_file_name );
int file_row_count ( char *input_file_name );
int find_closest ( int m, int n, double x[], double generator[] );
int get_seed ( void );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4_to_halton ( int seed, int base[], int ndim, double r[] );
void lcvt_write ( int dim_num, int n, int seed_start, int sample_function_init,
  char* file_in_name, int sample_function_cvt, int sample_num_cvt, int cvt_it,
  double cvt_energy, int latin_it, double latin_energy, double cell_generator[],
  char *file_out_name );
void param_print ( int dim_num, int n, int cvt_it, int latin_it, int seed, 
  int seed_start, int sample_function_cvt, int sample_function_init, 
  int sample_num_cvt );
int prime ( int n );
double r8_epsilon ( void );
double r8_uniform_01 ( int *seed );
void r8mat_latinize ( int m, int n, double table[] );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
double *r8table_data_read ( char *input_filename, int m, int n );
int *r8vec_sort_heap_index_a ( int n, double a[] );
void region_sampler ( int m, int n, int n_total, double x[], 
  int sample_function, bool reset, int *seed );
bool s_eqi ( char *s1, char *s2 );
int s_len_trim ( char* s );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
int test_region ( double x[], int dim_num );
void timestamp ( void );
char *timestring ( void );
void tuple_next_fast ( int m, int n, int rank, int x[] );
