double cpu_time ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform ( int a, int b, int *seed );
void i4vec_print ( int n, int a[], string title );
int point_radial_unique_count ( int m, int n, double a[], int *seed );
int point_radial_tol_unique_count ( int m, int n, double a[], double tol,
  int *seed );
int point_radial_tol_unique_index ( int m, int n, double a[], double tol,
  int *seed, int undx[], int xdnu[] );
int point_tol_unique_count ( int m, int n, double a[], double tol );
int point_tol_unique_index ( int m, int n, double a[], double tol, int xdnu[] );
int point_unique_count ( int m, int n, double a[] );
void point_unique_index ( int m, int n, double a[], int unique_num, int undx[],
  int xdnu[] );
int r4_nint ( float x );
double r8_abs ( double x );
double r8_max ( double x, double y );
double *r8col_duplicates ( int m, int n, int n_unique, int *seed );
int *r8col_sort_heap_index_a ( int m, int n, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
int r8vec_compare ( int n, double a[], double b[] );
double r8vec_norm_l2 ( int n, double a[] );
void r8vec_print ( int n, double a[], string title );
int *r8vec_sort_heap_index_a ( int n, double a[] );
double r8vec_sum ( int n, double a[] );
void r8vec_uniform_01 ( int n, int *seed, double r[] );
double *r8vec_uniform_01_new ( int n, int *seed );
void timestamp ( );
