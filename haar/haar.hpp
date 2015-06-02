void haar_1d ( int n, double x[] );
void haar_1d_inverse ( int n, double x[] );
void haar_2d ( int m, int n, double u[] );
void haar_2d_inverse ( int m, int n, double u[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double *r8mat_copy_new ( int m, int n, double a1[] );
double r8mat_dif_fro ( int m, int n, double a[], double b[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
double *r8vec_copy_new ( int n, double a[] );
double r8vec_diff_norm ( int n, double a[], double b[] );
double *r8vec_linspace_new ( int n, double a_lo, double a_hi );
double *r8vec_ones_new ( int n );
void r8vec_transpose_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
