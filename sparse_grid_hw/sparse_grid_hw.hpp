int cce_order ( int l );
int ccl_order ( int l );
int ccs_order ( int l );
void cc ( int n, double x[], double w[] );
double cpu_time ( );
double fn_integral ( int d );
double *fn_value ( int d, int n, double x[] );
double fu_integral ( int d );
double *fu_value ( int d, int n, double x[] );
int *get_seq ( int d, int norm, int seq_num );
void gqn ( int n, double x[], double w[] );
int gqn_order ( int l );
int gqn2_order ( int l );
void gqu ( int n, double x[], double w[] );
int gqu_order ( int l );
int i4_choose ( int n, int k );
int i4_factorial2 ( int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_mop ( int i );
int i4_power ( int i, int j );
void i4mat_print ( int m, int n, int a[], string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title );
int *i4vec_cum0_new ( int n, int a[] );
void i4vec_print ( int n, int a[], string title );
int i4vec_product ( int n, int a[] );
int i4vec_sum ( int n, int a[] );
void i4vec_transpose_print ( int n, int a[], string title );
void kpn ( int n, double x[], double w[] );
int kpn_order ( int l );
void kpu ( int n, double x[], double w[] );
int kpu_order ( int l );
int num_seq ( int n, int k );
void nwspgr ( void rule ( int n, double x[], double w[] ), 
  int rule_order ( int l ), int dim, int k, int r_size, int &s_size, 
  double nodes[], double weights[] );
int nwspgr_size ( int rule_order ( int l ), int dim, int k );
void quad_rule_print ( int m, int n, double x[], double w[], string title );
double r8_abs ( double x );
double r8_uniform_01 ( int &seed );
int *r8cvv_offset ( int m, int nr[] );
void r8cvv_print ( int mn, double a[], int m, int roff[], string title );
double *r8cvv_rget_new ( int mn, double a[], int m, int roff[], int i );
void r8cvv_rset ( int mn, double a[], int m, int roff[], int i, double ai[] );
double *r8mat_normal_01_new ( int m, int n, int &seed );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_direct_product ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double x[] );
void r8vec_direct_product2 ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double w[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double *r8vec_normal_01_new ( int n, int &seed );
void r8vec_print ( int n, double a[], string title );
double r8vec_sum ( int n, double a[] );
void r8vec_transpose_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
void rule_adjust ( double a, double b, double c, double d, int n, double x[], 
  double w[] );
void rule_sort ( int m, int n, double x[], double w[] );
void sort_heap_external ( int n, int &indx, int &i, int &j, int isgn );
int symmetric_sparse_size ( int nr, int dim, double nodes[], double x0 );
void tensor_product ( int d, int order1d[], int n1d, double x1d[], 
  double w1d[], int n, double xnd[], double wnd[] );
void tensor_product_cell ( int nc, double xc[], double wc[], int dim, int nr[], 
  int roff[], int np, double xp[], double wp[] );
void timestamp ( );

