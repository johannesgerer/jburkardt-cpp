int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void plane_normal_basis_3d ( double pp[3], double pn[3], double pq[3],
  double pr[3] );
double r8_abs ( double x );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
double *r8vec_any_normal ( int dim_num, double v1[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double *r8vec_cross_product_3d ( double v1[3], double v2[3] );
double r8vec_norm ( int n, double a[] );
double r8vec_norm_affine ( int n, double v0[], double v1[] );
void r8vec_normal_01 ( int n, int *seed, double x[] );
void r8vec_transpose_print ( int n, double x[], string title );
void r8vec_uniform_01 ( int n, int *seed, double r[] );
double *r8vec_uniform_01_new ( int n, int *seed );
void r8vec_zero ( int n, double a[] );
double *sphere_stereograph ( int m, int n, double p[] );
double *sphere_stereograph_inverse ( int m, int n, double q[] );
double *sphere_stereograph2 ( int m, int n, double p[], double focus[],
  double center[] );
double *sphere_stereograph2_inverse ( int m, int n, double q[], double focus[],
  double center[] );
void timestamp ( );
double *uniform_on_sphere01_map ( int dim_num, int n, int *seed );
