double *ellipsoid_sample ( int m, int n, double a[], double v[], double r, 
  int &seed );
double ellipsoid_volume ( int m, double a[], double v[], double r );
double hypersphere_unit_volume ( int m );
double *monomial_value ( int m, int n, int e[], double x[] );
void r8_print ( double r, string title );
double r8_uniform_01 ( int &seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
double *r8po_fa ( int n, double a[] );
double *r8po_sl ( int n, double a_lu[], double b[] );
double r8vec_norm ( int n, double a[] );
double *r8vec_normal_01_new ( int n, int &seed );
void r8vec_print ( int n, double a[], string title );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
double *uniform_in_sphere01_map ( int dim_num, int n, int &seed );

