int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8_abs ( double x );
double r8_factorial ( int n );
double r8mat_det ( int n, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
double *r8mat_zero_new ( int m, int n );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_norm ( int n, double a[] );
double r8vec_sum ( int n, double a[] );
double *simplex_coordinates1 ( int n );
double *simplex_coordinates2 ( int n );
double simplex_volume ( int n, double x[] );
void timestamp ( );
