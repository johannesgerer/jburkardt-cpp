double *bartlett_sample ( int m, int df, double sigma[] );
double *bartlett_unit_sample ( int m, int df );
void jacobi_eigenvalue ( int n, double a[], int it_max, double v[], 
  double d[], int &it_num, int &rot_num );
void r8mat_add ( int m, int n, double a[], double b[] );
double *r8mat_cholesky_factor_upper ( int n, double a[], int &flag );
double *r8mat_copy_new ( int m, int n, double a1[] );
void r8mat_diag_get_vector ( int n, double a[], double v[] );
double *r8mat_diagonal_new ( int n, double diag[] );
void r8mat_divide ( int m, int n, double s, double a[] );
void r8mat_identity  ( int n, double a[] );
double *r8mat_identity_new ( int n );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] );
double r8mat_norm_fro_affine ( int m, int n, double a1[], double a2[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_zero_new ( int m, int n );
double *r8ut_inverse ( int n, double a[] );
void r8vec_print ( int n, double a[], string title );
double *wishart_sample ( int m, int df, double sigma[] );
double *wishart_sample_inverse ( int m, int df, double sigma[] );
double *wishart_unit_sample ( int m, int df );
double *wishart_unit_sample_inverse ( int m, int df );;
