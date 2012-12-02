int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8mat_det_4d ( double a[4*4] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int *seed );
double *reference_to_physical_tet4 ( double t[], int n, double ref[] );
int s_len_trim ( string s );
double *tetrahedron_integrand_01 ( int p_num, double p[], int f_num );
double *tetrahedron_integrand_02 ( int p_num, double p[], int f_num );
double *tetrahedron_integrand_03 ( int p_num, double p[], int f_num );
double *tetrahedron_integrand_04 ( int p_num, double p[], int f_num );
double *tetrahedron_integrand_05 ( int p_num, double p[], int f_num );
double *tetrahedron_monte_carlo ( double t[], int p_num, int f_num, 
  double *tetrahedron_unit_sample ( int p_num, int *seed ), 
  double *tetrahedron_integrand ( int p_num, double p[], int f_num ), int *seed );
double *tetrahedron_unit_sample_01 ( int p_num, int *seed );
double *tetrahedron_unit_sample_02 ( int p_num, int *seed );
double *tetrahedron_unit_sample_03 ( int p_num, int *seed );
double *tetrahedron_unit_sample_04 ( int p_num, int *seed );
double tetrahedron_volume ( double tetra[3*4] );
void timestamp ( );
