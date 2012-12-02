char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
void cvt ( int dim_num, int n, int batch, int init, int sample, int sample_num, 
  int it_max, int it_fixed, int *seed, double r[], int *it_num, 
  double *it_diff, double *energy );
double cvt_energy ( int dim_num, int n, int batch, int sample, bool initialize,
  int sample_num, int *seed, double r[] );
void cvt_iterate ( int dim_num, int n, int batch, int sample, bool initialize, 
  int sample_num, int *seed, double r[], double *it_diff, double *energy );
void cvt_sample ( int dim_num, int n, int n_now, int sample, bool initialize, 
  int *seed, double r[] );
void data_read ( char *file_in_name, int dim_num, int n, double r[] );
char digit_to_ch ( int i );
void find_closest ( int dim_num, int n, int sample_num, double s[], double r[],
  int nearest[] );
int get_seed ( void );
bool halham_leap_check ( int dim_num, int leap[] );
bool halham_n_check ( int n );
bool halham_dim_num_check ( int dim_num );
bool halham_seed_check ( int dim_num, int seed[] );
bool halham_step_check ( int step );
bool halton_base_check ( int dim_num, int base[] );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4_to_halton_sequence ( int dim_num, int n, int step, int seed[], int leap[],
  int base[], double r[] );
char *i4_to_s ( int i );
int prime ( int n );
double r8_epsilon ( void );
double r8_huge ( void );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
unsigned long random_initialize ( int seed );
void s_blank_delete ( char *s );
void s_cap ( char *s );
bool s_eqi ( char *s1, char *s2 );
int s_len_trim ( char* s );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double rvec[] );
void timestamp ( void );
void tuple_next_fast ( int m, int n, int rank, int x[] );
void user ( int dim_num, int n, int *seed, double r[] );
