int i4_min ( int i1, int i2 );
double *orth_random ( int n, int &seed );
double *pds_random ( int n, int &seed );
double r8_normal_01 ( int &seed );
double r8_sign ( double x );
double r8_uniform_01 ( int &seed );
void r83_cg ( int n, double a[], double b[], double x[] );
double *r83_dif2 ( int m, int n );
double *r83_mv ( int m, int n, double a[], double x[] );
double *r83_resid ( int m, int n, double a[], double x[], double b[] );
void r83s_cg ( int n, double a[], double b[], double x[] );
double *r83s_dif2 ( int m, int n );
double *r83s_mv ( int m, int n, double a[], double x[] );
double *r83s_resid ( int m, int n, double a[], double x[], double b[] );
void r83t_cg ( int n, double a[], double b[], double x[] );
double *r83t_dif2 ( int m, int n );
double *r83t_mv ( int m, int n, double a[], double x[] );
double *r83t_resid ( int m, int n, double a[], double x[], double b[] );
void r8ge_cg ( int n, double a[], double b[], double x[] );
double *r8ge_dif2 ( int m, int n );
double *r8ge_mv ( int m, int n, double a[], double x[] );
double *r8ge_resid ( int m, int n, double a[], double x[], double b[] );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
void r8mat_house_axh ( int n, double a[], double v[] );
double *r8mat_identity_new ( int n );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_zero_new ( int m, int n );
void r8pbu_cg ( int n, int mu, double a[], double b[], double x[] );
double *r8pbu_dif2 ( int m, int n, int mu );
double *r8pbu_mv ( int m, int n, int mu, double a[], double x[] );
double *r8pbu_resid ( int m, int n, int mu, double a[], double x[], double b[] );
void r8sd_cg ( int n, int ndiag, int offset[], double a[], double b[], 
  double x[] );
double *r8sd_dif2 ( int m, int n, int ndiag, int offset[] );
double *r8sd_mv ( int m, int n, int ndiag, int offset[], double a[], double x[] );
double *r8sd_resid ( int m, int n, int ndiag, int offset[], double a[], 
  double x[], double b[] );
void r8sp_cg ( int n, int nz_num, int row[], int col[], double a[], 
  double b[], double x[] );
double *r8sp_dif2 ( int m, int n, int nz_num, int row[], int col[] );
double *r8sp_mv ( int m, int n, int nz_num, int row[], int col[], 
  double a[], double x[] );
double *r8sp_resid ( int m, int n, int nz_num, int row[], int col[], double a[], 
  double x[], double b[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_diff_norm ( int n, double a[], double b[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double *r8vec_house_column ( int n, double a[], int k );
double r8vec_norm ( int n, double a[] );
void r8vec_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
double *r8vec_zero_new ( int n );
void timestamp ( );

