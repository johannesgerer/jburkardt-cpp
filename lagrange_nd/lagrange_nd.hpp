int *comp_unrank_grlex ( int kc, int rank );
int i4_choose ( int n, int k );
int i4_min ( int i1, int i2 );
int i4_max ( int i1, int i2 );
void i4mat_print ( int m, int n, int a[], string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void i4vec_concatenate ( int n1, int a[], int n2, int b[], int c[] );
void i4vec_permute ( int n, int p[], int a[] );
void i4vec_print ( int n, int a[], string title );
int *i4vec_sort_heap_index_a ( int n, int a[] );
int i4vec_sum ( int n, int a[] );
double *interpolant_value ( int d, int r, int pn, int po[], double pc[], 
  int pe[], double pd[], int ni, double xi[] );
void lagrange_complete ( int d, int n, int r, int nd, double xd[], int po[], 
  double pc[], int pe[] );
void lagrange_complete2 ( int d, int n, int r, int nd, double xd[], int po[], 
  double pc[], int pe[] );
void lagrange_partial ( int d, int n, int r, int nd, double xd[], int po[], 
  double pc[], int pe[] );
void lagrange_partial2 ( int d, int n, int r, int nd, double xd[], int po[], 
  double pc[], int pe[] );
void lagrange_partial3 ( int d, int n, int nd, double xd[], int option, 
  int po[], double **pc, int **pe, int &n2 );
bool lagrange_partial4 ( int d, int n, int r, int nd, double xd[], 
  int option, double tol, int po[], double pc[], int pe[] );
void lp_coefficients ( int n, int &o, double c[], int f[] );
void lpp_to_polynomial ( int m, int l[], int o_max, int &o, double c[], int e[] );
int mono_between_enum ( int d, int n1, int n2 );
void mono_between_next_grlex ( int d, int n1, int n2, int x[] );
void mono_next_grlex ( int d, int x[] );
int mono_rank_grlex ( int m, int x[] );
int mono_total_enum ( int d, int n );
void mono_total_next_grlex ( int d, int n, int x[] );
int *mono_unrank_grlex ( int d, int rank );
int mono_upto_enum ( int d, int n );
double *mono_value ( int d, int nx, int f[], double x[] );
void perm_check0 ( int n, int p[] );
void polynomial_axpy ( double s, int o1, double c1[], int e1[], int o2, 
  double c2[], int e2[], int &o, double c[], int e[] );
void polynomial_compress ( int o1, double c1[], int e1[], int &o2, double c2[], 
  int e2[] );
void polynomial_print ( int d, int o, double c[], int e[], string title );
void polynomial_sort ( int o, double c[], int e[] );
double *polynomial_value ( int d, int o, double c[], int e[], int nx, 
  double x[] );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void r8col_separation ( int m, int n, double a[], double &d_min, double &d_max );
double r8mat_is_identity ( int n, double a[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8vec_concatenate ( int n1, double a[], int n2, double b[], double c[] );
void r8vec_permute ( int n, int p[], double a[] );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );

