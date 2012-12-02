double *asset_path ( double s0, double mu, double sigma, double t1, int n, 
  int *seed );
double binomial ( double s0, double e, double r, double sigma, double t1, 
  int m );
double bsf ( double s0, double t0, double e, double r, double sigma, double t1 );
double *forward ( double e, double r, double sigma, double t1, int nx, 
  int nt, double smax );
double *mc ( double s0, double e, double r, double sigma, double t1, int m, 
  int *seed );
double r8_max ( double x, double y );
double *r8vec_normal_01_new ( int n, int *seed );
void r8vec_print_part ( int n, double a[], int max_print, string title );
double *r8vec_uniform_01_new ( int n, int *seed );
void r8vec_write ( string output_filename, int n, double x[] );
void timestamp ( );
