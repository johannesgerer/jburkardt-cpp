double *monomial_value ( int m, int n, int e[], double x[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
double wedge01_integral ( int e[] );
double *wedge01_sample ( int n, int &seed );
double wedge01_volume ( );
