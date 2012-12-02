void covariance ( int m, int n, int x[], double *average, double *std,
  double *covc );
char digit_to_ch ( int i );
int get_seed ( );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
char *i4_to_s ( int i );
int i4_uniform ( int b, int c, int *seed );
void ihs ( int m, int n, int d, int *seed, int x[] );
void ihs_write ( int m, int n, int d, int seed_init, int seed,
  int r[], char *file_out_name );
float r4_abs ( float x );
int r4_nint ( float x );
double r8_epsilon ( );
double r8_huge ( );
double r8_uniform_01 ( int *seed );
double r8vec_average ( int n, double a[] );
double r8vec_std ( int n, double a[] );
void timestamp ( );
