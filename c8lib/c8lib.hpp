double c8_abs ( complex <double> x );
double c8_argument ( complex <double> x );
complex <double> c8_cube_root ( complex <double> x );
complex <double> c8_i ( );
bool c8_le_l1 ( complex <double> x, complex <double> y );
bool c8_le_l2 ( complex <double> x, complex <double> y );
bool c8_le_li ( complex <double> x, complex <double> y );
double c8_magnitude ( complex <double> x );
double c8_norm_l1 ( complex <double> x );
double c8_norm_l2 ( complex <double> x );
double c8_norm_li ( complex <double> x );
complex <double> c8_normal_01 ( int *seed );
complex <double> c8_one ( );
void c8_print ( complex <double> a, string title );
complex <double> c8_sqrt ( complex <double> x );
void c8_swap ( complex <double> *x, complex <double> *y );
complex <double> c8_uniform_01 ( int *seed );
complex <double> c8_zero ( );
void c8mat_copy ( int m, int n, complex <double> a1[], complex <double> a2[] );
complex <double> *c8mat_copy_new ( int m, int n, complex <double> a1[] );
complex <double> *c8mat_identity ( int n );
complex <double> *c8mat_indicator_new ( int m, int n );
complex <double> *c8mat_mm_new ( int n1, int n2, int n3, complex <double> a[], 
  complex <double> b[] );
void c8mat_nint ( int m, int n, complex <double> a[] );
void c8mat_print ( int m, int n, complex <double> a[], string title );
void c8mat_print_some ( int m, int n, complex <double> a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void c8mat_uniform_01 ( int m, int n, int *seed, complex <double> c[] );
complex <double> *c8mat_uniform_01_new ( int m, int n, int *seed );
void c8vec_copy ( int n, complex <double> a1[], complex <double> a2[] );
complex <double> *c8vec_copy_new ( int n, complex <double> a1[] );
complex <double> *c8vec_indicator_new ( int n );
double c8vec_norm_l2 ( int n, complex <double> a[] );
void c8vec_print ( int n, complex <double> a[], string title );
void c8vec_print_part ( int n, complex <double> a[], int max_print, 
  string title );
void c8vec_print_some ( int n, complex <double> a[], int i_lo, int i_hi, 
  string title );
complex <double> *c8vec_uniform_01_new ( int n, int *seed );
complex <double> *c8vec_unity ( int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8_abs ( double x );
complex <double> r8_csqrt ( double x );
double r8_max ( double x, double y );
int r8_nint ( double x );
double r8_sign ( double x );
double r8_uniform_01 ( int *seed );
void r8poly2_root ( double a, double b, double c, complex <double> *r1,
  complex <double> *r2 );
void r8poly3_root ( double a, double b, double c, double d, complex <double> *r1,
  complex <double> *r2, complex <double> *r3 );
void r8poly4_root ( double a, double b, double c, double d, double e,
  complex <double> *r1, complex <double> *r2, complex <double> *r3,
  complex <double> *r4 );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void timestamp ( );
