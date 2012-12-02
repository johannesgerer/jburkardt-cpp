int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int *seed );
void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] );
int s_len_trim ( string s );
void timestamp ( );
double triangle_area ( double t[2*3] );
double *triangle_integrand_01 ( int p_num, double p[], int f_num );
double *triangle_integrand_02 ( int p_num, double p[], int f_num );
double *triangle_integrand_03 ( int p_num, double p[], int f_num );
double *triangle_integrand_04 ( int p_num, double p[], int f_num );
double *triangle_integrand_05 ( int p_num, double p[], int f_num );
double *triangle_monte_carlo ( double t[], int p_num, int f_num, 
  double *triangle_unit_sample ( int p_num, int *seed ), 
  double *triangle_integrand ( int p_num, double p[], int f_num ), int *seed );
double *triangle_unit_sample_01 ( int p_num, int *seed );
double *triangle_unit_sample_02 ( int p_num, int *seed );
double *triangle_unit_sample_03 ( int p_num, int *seed );
double *triangle_unit_sample_04 ( int p_num, int *seed );
