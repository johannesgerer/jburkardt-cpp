char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
double chi_measure ( int ndim, int n, double z[], int ns, int seed_init );
void cvt_iterate ( int ndim, int n, int batch, int sample, bool reset, 
  int sample_num, int *seed, double r[], double *energy );
void cvt_sample ( int ndim, int n, int n_now, int sample, bool reset, 
  int *seed, double r[] );
void cvt_write ( int ndim, int n, int batch, int seed_init, int seed, 
  char *init_string, int it_max, int it_num, char *sample_string, 
  int sample_num, double r[], char *file_out_name );
double d_measure ( int ndim, int n, double z[], int ns, int seed_init );
double dge_det ( int n, double a[] );
char digit_to_ch ( int i );
double *dtable_data_read ( char *input_filename, int m, int n );
void dtable_data_write ( int m, int n, double table[], ofstream &output );
void dtable_header_read ( char *input_filename, int *m, int *n );
void dtable_header_write ( int m, int n, char *output_filename,
   ofstream &output );
int file_column_count ( char *input_filename );
void file_name_ext_get ( char *file_name, int *i, int *j );
char *file_name_ext_swap ( char *file_name, char *ext );
int file_row_count ( char *input_filename );
void find_closest ( int ndim, int n, int sample_num, double s[], double r[],
  int nearest[] );
void fsu_cvt ( int ndim, int n, int batch, int init, int sample, int sample_num, 
  int it_max, int *seed, double r[], int *it_num );
void fsu_halton ( int ndim, int n, int step, int seed[], int leap[],
  int base[], double r[] );
void fsu_hammersley  ( int ndim, int n, int step, int seed[], int leap[],
  int base[], double r[] );
void fsu_latinize ( int m, int n, double table[] );
int get_seed ( void );
double h_measure ( int ndim, int n, double z[], int ns, int seed_init );
bool halham_dim_num_check ( int dim_num );
bool halham_leap_check ( int ndim, int leap[] );
bool halham_n_check ( int n );
bool halham_ndim_check ( int ndim );
bool halham_seed_check ( int ndim, int seed[] );
bool halham_step_check ( int step );
void halham_write ( int ndim, int n, int step, int seed[], int leap[], 
  int base[], double r[], char *file_out_name );
bool halton_base_check ( int ndim, int base[] );
bool hammersley_base_check ( int ndim, int base[] );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
char *i4_to_s ( int i );
void i4vec_transpose_print ( int n, int a[], char *title );
double *pointset_spacing ( int ndim, int n, double z[] );
int prime ( int n );
double r8_epsilon ( void );
double r8_huge ( void );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
bool r8mat_in_01 ( int m, int n, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
int *r8vec_sort_heap_index_a ( int n, double a[] );
void r8vec_uniform_01 ( int n, int *seed, double r[] );
unsigned long random_initialize ( int seed );
bool s_eqi ( char *s1, char *s2 );
int s_index_last_c ( char *s, char c );
int s_len_trim ( char *s );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
double tau_measure ( int ndim, int n, double z[], int ns, int seed_init );
void timestamp ( void );
char *timestring ( void );
void tuple_next_fast ( int m, int n, int rank, int x[] );

