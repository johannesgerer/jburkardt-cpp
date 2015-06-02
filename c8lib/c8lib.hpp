double            c8_abs ( complex <double> x );
complex <double>  c8_acos ( complex <double> c1 );
complex <double>  c8_acosh ( complex <double> c1 );
complex <double>  c8_add ( complex <double> c1, complex <double> c2 );
double            c8_arg ( complex <double> x );
complex <double>  c8_asin ( complex <double> c1 );
complex <double>  c8_asinh ( complex <double> c1 );
complex <double>  c8_atan ( complex <double> c1 );
complex <double>  c8_atanh ( complex <double> c1 );
complex <double>  c8_conj ( complex <double> c1 );
void              c8_copy ( complex <double> c1, complex <double> c2 );
complex <double>  c8_cos ( complex <double> c1 );
complex <double>  c8_cosh ( complex <double> c1 );
complex <double>  c8_cube_root ( complex <double> x );
complex <double>  c8_div ( complex <double> c1, complex <double> c2 );
complex <double>  c8_div_r8 ( complex <double> c1, double r );
complex <double>  c8_exp ( complex <double> c1 );
complex <double>  c8_i ( );
double            c8_imag ( complex <double> c );
complex <double>  c8_inv ( complex <double> c1 );
bool              c8_le_l1 ( complex <double> x, complex <double> y );
bool              c8_le_l2 ( complex <double> x, complex <double> y );
bool              c8_le_li ( complex <double> x, complex <double> y );
complex <double>  c8_log ( complex <double> c1 );
double            c8_mag ( complex <double> x );
complex <double>  c8_mul ( complex <double> c1, complex <double> c2 );
complex <double>  c8_neg ( complex <double> c1 );
complex <double>  c8_nint ( complex <double> c1 );
double            c8_norm_l1 ( complex <double> x );
double            c8_norm_l2 ( complex <double> x );
double            c8_norm_li ( complex <double> x );
complex <double>  c8_normal_01 ( int &seed );
complex <double>  c8_one ( );
void              c8_print ( complex <double> a, string title );
double            c8_real ( complex <double> c );
complex <double>  c8_sin ( complex <double> c1 );
complex <double>  c8_sinh ( complex <double> c1 );
complex <double>  c8_sqrt ( complex <double> x );
complex <double>  c8_sub ( complex <double> c1, complex <double> c2 );
void              c8_swap ( complex <double> &x, complex <double> &y );
complex <double>  c8_tan ( complex <double> c1 );
complex <double>  c8_tanh ( complex <double> c1 );
void              c8_to_cartesian ( complex <double> c, double &x, double &y );
void              c8_to_polar ( complex <double> c, double &r, double &theta );
complex <double>  c8_uniform_01 ( int &seed );
complex <double>  c8_zero ( );

void              c8mat_add ( int m, int n, complex <double> alpha, complex <double> a[],
                  complex <double> beta, complex <double> b[], complex <double> c[] );
void              c8mat_add_r8 ( int m, int n, double alpha, complex <double> a[],
                    double beta, complex <double> b[], complex <double> c[] );
void              c8mat_copy ( int m, int n, complex <double> a1[], 
                    complex <double> a2[] );
complex <double> *c8mat_copy_new ( int m, int n, complex <double> a1[] );
void                c8mat_fss ( int n, complex <double> a[], int nb, complex <double> x[] );
complex <double> *c8mat_fss_new ( int n, complex <double> a[], int nb, 
                    complex <double> b[] );
complex <double> *c8mat_identity_new ( int n );
complex <double> *c8mat_indicator_new ( int m, int n );
void              c8mat_minvm ( int n1, int n2, complex <double> a[], 
                    complex <double> b[], complex <double> e[] );
complex <double> *c8mat_minvm_new ( int n1, int n2, complex <double> a[], 
                    complex <double> b[] );
void              c8mat_mm ( int n1, int n2, int n3, complex <double> a[], 
                    complex <double> b[], complex <double> c[] );
complex <double> *c8mat_mm_new ( int n1, int n2, int n3, complex <double> a[], 
                    complex <double> b[] );
void              c8mat_nint ( int m, int n, complex <double> a[] );
double            c8mat_norm_fro ( int m, int n, complex <double> a[] );
double            c8mat_norm_l1 ( int m, int n, complex <double> a[] );
double            c8mat_norm_li ( int m, int n, complex <double> a[] );
void              c8mat_print ( int m, int n, complex <double> a[], string title );
void              c8mat_print_some ( int m, int n, complex <double> a[], int ilo, 
                    int jlo, int ihi, int jhi, string title );
void              c8mat_scale ( int m, int n, complex <double> alpha, complex <double> a[] );
void              c8mat_scale_r8 ( int m, int n, double alpha, complex <double> a[] );
void              c8mat_uniform_01 ( int m, int n, int &seed, 
                    complex <double> c[] );
complex <double> *c8mat_uniform_01_new ( int m, int n, int &seed );
complex <double> *c8mat_zero_new ( int m, int n );

void              c8vec_copy ( int n, complex <double> a1[], 
                    complex <double> a2[] );
complex <double> *c8vec_copy_new ( int n, complex <double> a1[] );
complex <double> *c8vec_indicator_new ( int n );
void              c8vec_nint ( int n, complex <double> a[] );
double            c8vec_norm_l1 ( int n, complex <double> a[] );
double            c8vec_norm_l2 ( int n, complex <double> a[] );
double            c8vec_norm_li ( int n, complex <double> a[] );
void              c8vec_print ( int n, complex <double> a[], string title );
void              c8vec_print_part ( int n, complex <double> a[], int max_print, 
                    string title );
void              c8vec_print_some ( int n, complex <double> a[], int i_lo, 
                    int i_hi, string title );
void              c8vec_sort_a_l1 ( int n, complex <double> x[] );
void              c8vec_sort_a_l2 ( int n, complex <double> x[] );
void              c8vec_sort_a_li ( int n, complex <double> x[] );
complex <double> *c8vec_spiral_new ( int n, int m, complex <double> c1, 
                    complex <double> c2 );
complex <double> *c8vec_uniform_01_new ( int n, int &seed );
complex <double> *c8vec_unity_new ( int n );

complex <double>  cartesian_to_c8 ( double x, double y );

complex <double>  polar_to_c8 ( double r, double theta );

double            r8_abs ( double x );
double            r8_atan ( double y, double x );
complex <double>  r8_csqrt ( double x );
double            r8_floor ( double x );
double            r8_max ( double x, double y );
double            r8_sign ( double x );
double            r8_uniform_01 ( int &seed );
void              r8poly2_root ( double a, double b, double c, 
                    complex <double> &r1, complex <double> &r2 );
void              r8poly3_root ( double a, double b, double c, double d, 
                    complex <double> &r1, complex <double> &r2, 
                    complex <double> &r3 );
void              r8poly4_root ( double a, double b, double c, double d, double e,
                    complex <double> &r1, complex <double> &r2, complex <double> &r3,
                    complex <double> &r4 );
void              sort_heap_external ( int n, int &indx, int &i, int &j, int isgn );
void              timestamp ( );
