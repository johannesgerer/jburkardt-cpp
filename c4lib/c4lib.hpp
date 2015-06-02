float            c4_abs ( complex <float> x );
complex <float>  c4_acos ( complex <float> c1 );
complex <float>  c4_acosh ( complex <float> c1 );
complex <float>  c4_add ( complex <float> c1, complex <float> c2 );
float            c4_arg ( complex <float> x );
complex <float>  c4_asin ( complex <float> c1 );
complex <float>  c4_asinh ( complex <float> c1 );
complex <float>  c4_atan ( complex <float> c1 );
complex <float>  c4_atanh ( complex <float> c1 );
complex <float>  c4_conj ( complex <float> c1 );
void             c4_copy ( complex <float> c1, complex <float> c2 );
complex <float>  c4_cos ( complex <float> c1 );
complex <float>  c4_cosh ( complex <float> c1 );
complex <float>  c4_cube_root ( complex <float> x );
complex <float>  c4_div ( complex <float> c1, complex <float> c2 );
complex <float>  c4_div_r4 ( complex <float> c1, float r );
complex <float>  c4_exp ( complex <float> c1 );
complex <float>  c4_i ( );
float            c4_imag ( complex <float> c );
complex <float>  c4_inv ( complex <float> c1 );
bool             c4_le_l1 ( complex <float> x, complex <float> y );
bool             c4_le_l2 ( complex <float> x, complex <float> y );
bool             c4_le_li ( complex <float> x, complex <float> y );
complex <float>  c4_log ( complex <float> c1 );
float            c4_mag ( complex <float> x );
complex <float>  c4_mul ( complex <float> c1, complex <float> c2 );
complex <float>  c4_neg ( complex <float> c1 );
complex <float>  c4_nint ( complex <float> c1 );
float            c4_norm_l1 ( complex <float> x );
float            c4_norm_l2 ( complex <float> x );
float            c4_norm_li ( complex <float> x );
complex <float>  c4_normal_01 ( int *seed );
complex <float>  c4_one ( );
void             c4_print ( complex <float> a, string title );
float            c4_real ( complex <float> c );
complex <float>  c4_sin ( complex <float> c1 );
complex <float>  c4_sinh ( complex <float> c1 );
complex <float>  c4_sqrt ( complex <float> x );
complex <float>  c4_sub ( complex <float> c1, complex <float> c2 );
void             c4_swap ( complex <float> *x, complex <float> *y );
complex <float>  c4_tan ( complex <float> c1 );
complex <float>  c4_tanh ( complex <float> c1 );
void             c4_to_cartesian ( complex <float> c, float *x, float *y );
void             c4_to_polar ( complex <float> c, float *r, float *theta );
complex <float>  c4_uniform_01 ( int *seed );
complex <float>  c4_zero ( );
void             c4mat_add ( int m, int n, complex <float> alpha, complex <float> a[],
                 complex <float> beta, complex <float> b[], complex <float> c[] );
void             c4mat_add_r4 ( int m, int n, float alpha, complex <float> a[],
                 float beta, complex <float> b[], complex <float> c[] );
void             c4mat_copy ( int m, int n, complex <float> a1[], complex <float> a2[] );
complex <float> *c4mat_copy_new ( int m, int n, complex <float> a1[] );
void             c4mat_fss ( int n, complex <float> a[], int nb, complex <float> x[] );
complex <float> *c4mat_fss_new ( int n, complex <float> a[], int nb, 
                 complex <float> b[] );
complex <float> *c4mat_identity_new ( int n );
complex <float> *c4mat_indicator_new ( int m, int n );
void             c4mat_minvm ( int n1, int n2, complex <float> a[], 
                 complex <float> b[], complex <float> c[] );
complex <float> *c4mat_minvm_new ( int n1, int n2, complex <float> a[], 
                 complex <float> b[] );
void             c4mat_mm ( int n1, int n2, int n3, complex <float> a[], 
                 complex <float> b[], complex <float> c[] );
complex <float> *c4mat_mm_new ( int n1, int n2, int n3, complex <float> a[], 
                 complex <float> b[] );
void             c4mat_nint ( int m, int n, complex <float> a[] );
float            c4mat_norm_fro ( int m, int n, complex <float> a[] );
float            c4mat_norm_l1 ( int m, int n, complex <float> a[] );
float            c4mat_norm_li ( int m, int n, complex <float> a[] );
void             c4mat_print ( int m, int n, complex <float> a[], string title );
void             c4mat_print_some ( int m, int n, complex <float> a[], int ilo, int jlo, 
                 int ihi, int jhi, string title );
void             c4mat_scale ( int m, int n, complex <float> alpha, complex <float> a[] );
void             c4mat_scale_r4 ( int m, int n, float alpha, complex <float> a[] );
void             c4mat_uniform_01 ( int m, int n, int *seed, complex <float> c[] );
complex <float> *c4mat_uniform_01_new ( int m, int n, int *seed );
complex <float> *c4mat_zero_new ( int m, int n );
void             c4vec_copy ( int n, complex <float> a1[], complex <float> a2[] );
complex <float> *c4vec_copy_new ( int n, complex <float> a1[] );
complex <float> *c4vec_indicator_new ( int n );
void             c4vec_nint ( int n, complex <float> a[] );
float            c4vec_norm_l2 ( int n, complex <float> a[] );
void             c4vec_print ( int n, complex <float> a[], string title );
void             c4vec_print_part ( int n, complex <float> a[], int max_print, 
                 string title );
void             c4vec_print_some ( int n, complex <float> a[], int i_lo, int i_hi, 
                 string title );
void             c4vec_sort_a_l2 ( int n, complex <float> x[] );
complex <float> *c4vec_spiral ( int n, int m, complex <float> c1, 
                 complex <float> c2 );
void             c4vec_uniform_01 ( int n, int *seed, complex <float> c[] );
complex <float> *c4vec_uniform_01_new ( int n, int *seed );
complex <float>  cartesian_to_c4 ( float x, float y );
int              i4_max ( int i1, int i2 );
int              i4_min ( int i1, int i2 );
complex <float>  polar_to_c4 ( float r, float theta );
float            r4_abs ( float x );
complex <float>  r4_csqrt ( float x );
float            r4_floor ( float x );
float            r4_max ( float x, float y );
int              r4_nint ( float x );
float            r4_sign ( float x );
float            r4_uniform_01 ( int *seed );
void             r4poly2_root ( float a, float b, float c, complex <float> *r1,
                 complex <float> *r2 );
void             r4poly3_root ( float a, float b, float c, float d,
                 complex <float> *r1, complex <float> *r2, complex <float> *r3 );
void             r4poly4_root ( float a, float b, float c, float d, float e,
                 complex <float> *r1, complex <float> *r2, complex <float> *r3,
                 complex <float> *r4 );
void             sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void             timestamp ( );
