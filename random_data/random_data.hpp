double arc_cosine ( double c );
double *bad_in_simplex01 ( int dim_num, int point_num, int *seed );
double *brownian ( int dim_num, int n, int *seed );
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
double *dge_mxv ( int m, int n, double a[], double x[] );
void direction_uniform_nd ( int dim_num, int *seed, double w[] );
int dpofa ( double a[], int lda, int n );
void dposl ( double a[], int lda, int n, double b[] );
double *dut_mxv ( int m, int n, double a[], double x[] );
unsigned long get_seed ( );
double *grid_in_cube01 ( int dim_num, int n, int center, int *seed );
int grid_side ( int dim_num, int n );
bool halham_dim_num_check ( int dim_num );
bool halham_leap_check ( int dim_num, int leap[] );
bool halham_n_check ( int n );
bool halham_seed_check ( int dim_num, int seed[] );
bool halham_step_check ( int step );
bool halton_base_check ( int dim_num, int base[] );
double *halton_in_circle01_accept ( int m, int n, int *seed );
double *halton_in_circle01_map ( int m, int n, int *seed );
double *halton_in_cube01 ( int m, int n, int *seed );
bool hammersley_base_check ( int dim_num, int base[] );
double *hammersley_in_cube01 ( int m, int n, int *seed );
int i4_factorial ( int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
void i4_to_halton ( int dim_num, int step, int seed[], int leap[], int base[],
  double r[] );
void i4_to_halton_sequence ( int dim_num, int n, int step, int seed[], int leap[],
  int base[], double r[] );
void i4_to_hammersley ( int dim_num, int step, int seed[], int leap[], int base[],
  double r[] );
void i4_to_hammersley_sequence ( int dim_num, int n, int step, int seed[], int leap[],
  int base[], double r[] );
int i4_uniform ( int b, int c, int *seed );
void i4vec_transpose_print ( int n, int a[], string title );
void ksub_random2 ( int n, int k, int *seed, int a[] );
double *normal ( int m, int n, double r[], double mu[], int *seed );
double *normal_circular ( int m, int n, int *seed );
double *normal_multivariate ( int m, int n, double r[], double mu[],
  int *seed );
double *normal_simple ( int m, int n, int *seed );
double *polygon_centroid_2d ( int n, double v[] );
int prime ( int n );
float r4_abs ( float x );
int r4_nint ( float x );
double r8_epsilon ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_nint ( double x );
double r8_normal_01 ( int *seed );
double r8_pi ( );
double r8_uniform_01 ( int *seed );
double *r8mat_normal_01_new ( int m, int n, int *seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_norm ( int n, double a[] );
void r8vec_normal_01 ( int n, int *seed, double x[] );
double *r8vec_normal_01_new ( int n, int *seed );
void r8vec_print ( int n, double a[], string title );
double r8vec_sum ( int n, double a[] );
void r8vec_uniform_01 ( int n, int *seed, double r[] );
double *r8vec_uniform_01_new ( int n, int *seed );
double *r8vec_zero_new ( int n );
unsigned long random_initialize ( unsigned long seed );
int s_len_trim ( string s );
void scale_from_simplex01 ( int dim_num, int n, double t[], double x[] );
void scale_to_ball01 ( int dim_num, int n, double x[] );
void scale_to_block01 ( int dim_num, int n, double x[] );
void scale_to_cube01 ( int dim_num, int n, double x[] );
double stri_angles_to_area ( double r, double a, double b, double c );
void stri_sides_to_angles ( double r, double as, double bs, double cs,
  double *a, double *b, double *c );
void stri_vertices_to_sides ( double r, double v1[3], double v2[3],
  double v3[3], double *as, double *bs, double *cs );
void timestamp ( );
double triangle_area_2d ( double t[2*3] );
void tuple_next_fast ( int m, int n, int rank, int x[] );
double *uniform_in_annulus ( double pc[], double r1, double r2, int n, int *seed );
double *uniform_in_annulus_accept ( double pc[], double r1, double r2, int n,
  int *seed );
double *uniform_in_annulus_sector ( double pc[], double r1, double r2, double theta1,
  double theta2, int n, int *seed );
double *uniform_in_circle01_map ( int n, int *seed );
double *uniform_in_cube01 ( int m, int n, int *seed );
double *uniform_in_ellipsoid_map ( int dim_num, int n, double a[], double r,
  int *seed );
double *uniform_in_parallelogram_map ( double v1[2], double v2[2],
  double v3[2], int n, int *seed );
double *uniform_in_polygon_map ( int nv, double v[], int n, int *seed );
double *uniform_in_sector_map ( double r1, double r2, double t1,
  double t2, int n, int *seed );
double *uniform_in_simplex01_map ( int dim_num, int n, int *seed );
double *uniform_in_sphere01_map ( int dim_num, int n, int *seed );
double *uniform_in_tetrahedron ( double v[], int n, int *seed );
double *uniform_in_triangle_map1 ( double v1[2], double v2[2], double v3[2],
  int n, int *seed );
double *uniform_in_triangle_map2 ( double v1[2], double v2[2], double v3[2],
  int n, int *seed );
double *uniform_in_triangle01_map ( int n, int *seed );
double *uniform_on_ellipsoid_map ( int dim_num, int n, double a[],
  double r, int *seed );
double *uniform_on_hemisphere01_phong ( int n, int m, int *seed );
double *uniform_on_simplex01_map ( int dim_num, int n, int *seed );
double *uniform_on_sphere01_map ( int dim_num, int n, int *seed );
double *uniform_on_sphere01_patch_tp ( int n, double phi1, double phi2,
  double theta1, double theta2, int *seed );
double *uniform_on_sphere01_patch_xyz ( int n, double phi1, double phi2,
  double theta1, double theta2, int *seed );
double *uniform_on_sphere01_triangle_xyz ( int n, double v1[], double v2[],
  double v3[], int *seed );
double *uniform_walk ( int dim_num, int n, int *seed );

