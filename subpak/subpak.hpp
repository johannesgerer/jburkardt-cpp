double angle_shift ( double alpha, double beta );
double angle_shift_deg ( double alpha, double beta );
double *angle_to_rgb ( double angle );
void axis_limits ( double xmin, double xmax, int ndivs, double *pxmin, 
  double *pxmax, double *pxdiv, int *nticks );
int bar_check ( int digit[12] );
char *bar_code ( int digit[] );
char *bar_digit_code_left ( int digit );
char *bar_digit_code_right ( int digit );
void bin_to_d_even ( int nbin, int bin, double a, double b, double *cmin, 
  double *cmax );
double bmi_english ( double w_lb, double h_ft, double h_in );
double bmi_metric ( double w_kg, double h_m );
double c8_argument ( complex <double> x );
double c8_magnitude ( complex <double> x );
complex <double> c8_sqrt ( complex <double> x );
char ch_cap ( char c );
bool ch_eqi ( char ch1, char ch2 );
bool ch_is_digit ( char c );
int ch_to_digit ( char ch );
int charstar_len_trim ( char *s );
double degrees_to_radians ( double degrees );
double e_constant ( );
double euler_constant ( );
void fac_div ( int prime_num, int npower1[], int npower2[], int npower3[] );
void fac_gcd ( int prime_num, int npower1[], int npower2[], int npower3[] );
void fac_lcm ( int prime_num, int npower1[], int npower2[], int npower3[] );
void fac_mul ( int prime_num, int npower1[], int npower2[], int npower3[] );
void fac_print ( int prime_num, int npower[] );
int fac_to_i4 ( int prime_num, int npower[] );
void fac_to_rat ( int prime_num, int npower[], int *top, int *bot );
double feet_to_meters ( double ft );
double gauss_sum ( int dim_num, int n, double amplitude[], double center[], 
  double width[], double x[] );
unsigned long get_seed ( );
double *grid1 ( int dim_num, int nstep, double x1[], double x2[] );
double *grid1n ( int j, int dim_num, int nstep, double x1[], double x2[] );
double *grid2 ( int j1, int j2, int dim_num, int nstep, double x1[], 
  double x2[] );
double *grid2n ( int j, int j1, int j2, int dim_num, double x1[], double x2[] );
double *grid3 ( int dim_num, int nstep1, int nstep2, double x1[], double x2[],
  double x3[] );
double *grid3n ( int j, int k, int dim_num, int nstep1, int nstep2, 
  double x1[], double x2[], double x3[] );
double *grid4 ( int j1, int j2, int k1, int k2, int dim_num, int nstep1, 
  int nstep2, double x1[], double x2[], double x3[] );
double *grid4n ( int j, int j1, int j2, int k, int k1, int k2, int dim_num, 
  int nstep1, int nstep2, double x1[], double x2[], double x3[] );
short int i2_reverse_bytes ( short int x );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
int i4_sign ( int i );
void i4_swap ( int *i, int *j );
int *i4_to_digits_decimal ( int i, int n );
int *i4_to_fac ( int i, int prime_num );
char i4_to_isbn ( int i );
int i4_uniform ( int a, int b, int *seed );
double i4int_to_r8int ( int imin, int imax, int i, double rmin, double rmax );
void i4vec_copy ( int n, int a1[], int a2[] );
int *i4vec_indicator_new ( int n );
int i4vec_min ( int n, int a[] );
void i4vec_permute ( int n, int p[], int base, int a[] );
void i4vec_print ( int n, int a[], string title );
int i4vec_sum ( int n, int a[] );
void i4vec_zero ( int n, int a[] );
int *i4vec_zero_new ( int n );
void ij_next ( int *i, int *j, int n );
void ij_next_gt ( int *i, int *j, int n );
void index_box2_next_2d ( int n1, int n2, int ic, int jc, int *i, int *j, 
  int *more );
void index_box2_next_3d ( int n1, int n2, int n3, int ic, int jc, int kc, 
  int *i, int *j, int *k, bool *more );
int index1_col ( int i_min, int i, int i_max, int index_min );
int index1_row ( int i_min, int i, int i_max, int index_min );
int index2_col ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int index_min );
int index2_row ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int index_min );
int index3_col ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int k_min, int k, int k_max, int index_min );
int index3_row ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int k_min, int k, int k_max, int index_min );
int index4_col ( int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max, 
  int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max, 
  int index_min );
int index4_row ( int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max, 
  int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max, 
  int index_min );
int indexn_col ( int n, int i_min[], int i[], int i_max[], int index_min );
int indexn_row ( int n, int i_min[], int i[], int i_max[], int index_min );
int isbn_check ( string isbn );
void isbn_fill ( string isbn );
int isbn_to_i4 ( char c );
int iset2_compare ( int x1, int y1, int x2, int y2 );
int lcm_12n ( int n );
void lmat_print ( int m, int n, bool a[], string title );
void lmat_print_some ( int m, int n, bool a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void lmat_transpose_print ( int m, int n, bool a[], string title );
void lmat_transpose_print_some ( int m, int n, bool a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
int luhn_check ( int digit_num, int digit[] ) ;
void lvec_print ( int n, bool a[], string title );
int pause_input ( );
void pause_seconds ( int seconds );
bool perm_check ( int n, int p[], int base );
void perm_cycle ( int n, int p[], int *isgn, int *ncycle, int iopt );
int *perm_free ( int npart, int ipart[], int nfree );
void perm_inverse ( int n, int p[] );
void perm_print ( int n, int p[], string title );
int *perm_uniform_new ( int n, int base, int *seed );
double pounds_to_kilograms ( double lb );
int prime ( int n );
int prime_ge ( int n );
int r4_nint ( float x );
double r8_abs ( double x );
double r8_huge ( );
double r8_log_10 ( double x );
double r8_max ( double x, double y );
double r8_modp ( double x, double y );
double r8_uniform ( double b, double c, int *seed );
double r8_uniform_01 ( int *seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
int r8poly_degree ( int na, double a[] );
void r8poly_print ( int n, double a[], string title );
double *r8vec_indicator_new ( int n );
double r8vec_max ( int n, double r8vec[] );
double r8vec_mean ( int n, double x[] );
double r8vec_min ( int n, double r8vec[] );
void r8vec_print ( int n, double a[], string title );
double r8vec_variance ( int n, double x[] );
double *r8vec_zero_new ( int n );
double radians_to_degrees ( double angle );
unsigned long rand_initialize ( unsigned long seed );
unsigned long random_initialize ( unsigned long seed );
void rat_factor ( int m, int n, int factor_max, int *factor_num, int factor[], 
  int power[], int *mleft, int *nleft );
double rickey ( int ab, int bb, int er, double f, int h, int hb, int hp, 
  int r, int so, int tb );
int *roots_to_i4poly ( int n, int x[] );
double *roots_to_r8poly ( int n, double x[] );
bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
int s_to_i4 ( char *s, int *last, bool *error );
bool s_to_i4vec ( char *s, int n, int ivec[] );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double rvec[] );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void timestamp ( );
void tuple_next2 ( int n, int xmin[], int xmax[], int x[], int *rank );
double *tvec_even ( int nt );
double *tvec_even2 ( int nt );
double *tvec_even3 ( int nt );
double *tvec_even_bracket ( int nt, double theta1, double theta2 );
double *tvec_even_bracket2 ( int nt, double theta1, double theta2 );
double *tvec_even_bracket3 ( int nt, double theta1, double theta2 );
int upc_check_digit ( int p, int l, int r );
double versine_pulse ( double t, double ta, double tb, double v1, double amp );
