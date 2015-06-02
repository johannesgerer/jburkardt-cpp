//
//  Utilities for Complex double precision vectors.
//
void c8vec_print ( int n, double a[], string title );
void c8vec_sort_a2 ( int n, double x[] );
double *c8vec_unity ( int n );
//
//  Utility to add a multiple of one real vector to another.
//
void daxpy ( int n, double sa, double x[], int incx, double y[], int incy );
//
//  File utilities.
//
int file_delete ( string file_name );
bool file_exist ( string file_name );
//
//  Utility to get a seed for the random number generator.
//
int get_seed ( );
//
//  Utility for one of the tests.
//
double *hilbert_inverse ( int n );
//
//  Utilities for integers.
//
int i4_huge ( );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
int i4_uniform ( int ilo, int ihi, int *seed );
//
//  Utilities for integer vectors.
//
void i4vec_print ( int n, int a[], string title );
int i4vec_search_binary_a ( int n, int a[], int b );
//
//  R4 Utilities.
//
float r4_abs ( float x );
int r4_nint ( float x );
float r4_uniform ( float b, float c, int *seed );
float r4_uniform_01 ( int *seed );
//
//  R8 Utilities.
//
double r8_abs ( double x );
bool r8_is_int ( double r );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_sign ( double x );
double r8_sign2 ( double x, double y );
void r8_swap ( double *x, double *y );
double r8_uniform ( double rlo, double rhi, int *seed );
double r8_uniform_01 ( int *seed );
//
//  Real double precision Tridiagonal.
//
double *r83_cr_fa ( int n, double a[] );
double *r83_cr_sl ( int n, double a_cr[], double b[] );
double *r83_cr_sls ( int n, double a_cr[], int nb, double b[] );
void r83_gs_sl ( int n, double a[], double b[], double x[], int it_max, int job );
double *r83_indicator ( int n );
void r83_jac_sl ( int n, double a[], double b[], double x[], int it_max, int job );
double *r83_mxv ( int n, double a[], double x[] );
double r83_np_det ( int n, double a[] );
int r83_np_fa ( int n, double a[] );
double *r83_np_fs ( int n, double a[], double b[] );
double *r83_np_ml ( int n, double a[], double x[], int job );
double *r83_np_sl ( int n, double a[], double b[], int job );
void r83_print ( int n, double a[], string title );
void r83_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi, 
  string title );
double *r83_random ( int n, int *seed );
double *r83_to_r8ge ( int n, double a[] );
double *r83_vxm ( int n, double a[], double x[] );
double *r83_zero ( int n );
//
//  Real double precision Tridiagonal No Pivoting.
//
double *r83np_fs ( int n, double a[], double b[] );
//
//  Real double precision Tridiagonal Periodic.
//
double r83p_det ( int n, double a[], double work4 );
int r83p_fa ( int n, double a[], double work2[], double work3[], double *work4 );
double *r83p_indicator ( int n );
double *r83p_ml ( int n, double a[], double x[], int job );
double *r83p_mxv ( int n, double a[], double x[] );
void r83p_print ( int n, double a[], string title );
void r83p_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r83p_random ( int n, int *seed );
double *r83p_sl ( int n, double a[], double b[], int job, double work2[], 
  double work3[], double work4 );
double *r83p_to_r8ge ( int n, double a[] );
double *r83p_vxm ( int n, double a[], double x[] );
double *r83p_zero ( int n );
//
//  Real double precision Tridiagonal Scalar.
//
double *r83s_mxv ( int n, double a[], double x[] );
void r83s_print ( int n, double a[], string title );
void r83s_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r83s_random ( int n, int &seed );
//
//  Real double precision Pentagonal.
//
double *r85_indicator ( int n );
double *r85_mxv ( int n, double a[], double x[] );
double *r85_np_fs ( int n, double a[], double b[] );
void r85_print ( int n, double a[], string title );
void r85_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r85_random ( int n, int *seed );
double *r85_to_r8ge ( int n, double a[] );
double *r85_vxm ( int n, double a[], double x[] );
double *r85_zero ( int n );
//
//  Real double precision Border Banded.
//
void r8bb_add ( int n1, int n2, int ml, int mu, double a[], int i, int j, 
  double value );
int r8bb_fa ( int n1, int n2, int ml, int mu, double a[], int pivot[] );
double r8bb_get ( int n1, int n2, int ml, int mu, double a[], int i, int j );
double *r8bb_indicator ( int n1, int n2, int ml, int mu );
double *r8bb_mxv ( int n1, int n2, int ml, int mu, double a[], double x[] );
void r8bb_print ( int n1, int n2, int ml, int mu, double a[], string title );
void r8bb_print_some ( int n1, int n2, int ml, int mu, double a[], int ilo, 
  int jlo, int ihi, int jhi, string title );
double *r8bb_random ( int n1, int n2, int ml, int mu, int *seed );
void r8bb_set ( int n1, int n2, int ml, int mu, double a[], int i, int j, 
  double value );
double *r8bb_sl ( int n1, int n2, int ml, int mu, double a[], int pivot[], 
  double b[] );
double *r8bb_to_r8ge ( int n1, int n2, int ml, int mu, double a[] );
double *r8bb_vxm ( int n1, int n2, int ml, int mu, double a[], double x[] );
double *r8bb_zero ( int n1, int n2, int ml, int mu );
//
//  Real double precision Banded Lower Triangular.
//
double r8blt_det ( int n, int ml, double a[] );
double *r8blt_indicator ( int n, int ml );
double *r8blt_mxv ( int n, int ml, double a[], double x[] );
void r8blt_print ( int n, int ml, double a[], string title );
void r8blt_print_some ( int n, int ml, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8blt_random ( int n, int ml, int *seed );
double *r8blt_sl ( int n, int ml, double a[], double b[], int job );
double *r8blt_to_r8ge ( int n, int ml, double a[] );
double *r8blt_vxm ( int n, int ml, double a[], double x[] );
double *r8blt_zero ( int n, int ml );
//
//  Real double precision Block Toeplitz.
//
double *r8bto_indicator ( int m, int l );
double *r8bto_mxv ( int m, int l, double a[], double x[] );
void r8bto_print ( int m, int l, double a[], string title );
void r8bto_print_some ( int m, int l, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8bto_random ( int m, int l, int *seed );
double *r8bto_sl ( int m, int l, double a[], double b[] );
double *r8bto_to_r8ge ( int m, int l, double a[] );
double *r8bto_vxm ( int m, int l, double a[], double x[] );
double *r8bto_zero ( int m, int l );
//
//  Real double precision Banded Upper Triangular.
//
double r8but_det ( int n, int mu, double a[] );
double *r8but_indicator ( int n, int mu );
double *r8but_mxv ( int n, int mu, double a[], double x[] );
void r8but_print ( int n, int mu, double a[], string title );
void r8but_print_some ( int n, int mu, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8but_random ( int n, int mu, int *seed );
double *r8but_sl ( int n, int mu, double a[], double b[], int job );
double *r8but_to_r8ge ( int n, int mu, double a[] );
double *r8but_vxm ( int n, int mu, double a[], double x[] );
//
//  Real double precision Compact Banded.
//
double r8cb_det ( int n, int ml, int mu, double a[] );
double *r8cb_indicator ( int m, int n, int ml, int mu );
double *r8cb_ml ( int n, int ml, int mu, double a[], double x[], int job );
double *r8cb_mxv ( int n, int ml, int mu, double a[], double x[] );
int r8cb_np_fa ( int n, int ml, int mu, double a[] );
double *r8cb_np_sl ( int n, int ml, int mu, double a[], double b[], int job );
void r8cb_print ( int m, int n, int ml, int mu, double a[], string title );
void r8cb_print_some ( int m, int n, int ml, int mu, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
double *r8cb_random ( int n, int ml, int mu, int *seed );
double *r8cb_to_r8vec ( int m, int n, int ml, int mu, double *a );
double *r8cb_to_r8ge ( int n, int ml, int mu, double a[] );
double *r8cb_vxm ( int n, int ml, int mu, double a[], double x[] );
double *r8cb_zero ( int n, int ml, int mu );
//
//  Real double precision Compact Border Banded.
//
void r8cbb_add ( int n1, int n2, int ml, int mu, double a[], int i, int j, 
  double value );
bool r8cbb_error ( int n1, int n2, int ml, int mu );
int r8cbb_fa ( int n1, int n2, int ml, int mu, double a[] );
double r8cbb_get ( int n1, int n2, int ml, int mu, double a[], int i, int j );
double *r8cbb_indicator ( int n1, int n2, int ml, int mu );
double *r8cbb_mxv ( int n1, int n2, int ml, int mu, double a[], double x[] );
void r8cbb_print ( int n1, int n2, int ml, int mu, double a[], string title );
void r8cbb_print_some ( int n1, int n2, int ml, int mu, double a[], int ilo, 
  int jlo, int ihi, int jhi, string title );
double *r8cbb_random ( int n1, int n2, int ml, int mu, int *seed );
void r8cbb_set ( int n1, int n2, int ml, int mu, double a[], int i, int j, 
  double value );
double *r8cbb_sl ( int n1, int n2, int ml, int mu, double a[], double b[] );
double *r8cbb_to_r8ge ( int n1, int n2, int ml, int mu, double a[] );
double *r8cbb_vxm ( int n1, int n2, int ml, int mu, double a[], double x[] );
double *r8cbb_zero ( int n1, int n2, int ml, int mu );
//
//  Real double precision (Sparse) Compressed Column.
//
double r8cc_get ( int m, int n, int nnzero, int colptr[], int rowind[], 
  double a[], int i, int j );
int r8cc_ijk ( int m, int n, int nnzero, int colptr[], int rowind[], 
  int i, int j );
void r8cc_inc ( int m, int n, int nnzero, int colptr[], int rowind[], 
  double a[], int i, int j, double aij );
double *r8cc_indicator ( int m, int n, int nnzero, int colptr[], int rowind[] );
void r8cc_kij ( int m, int n, int nnzero, int colptr[], int rowind[], int k, 
  int *i, int *j );
double *r8cc_mxv ( int m, int n, int nnzero, int colptr[], int rowind[], 
  double a[], double x[] );
void r8cc_print ( int m, int n, int nnzero, int colptr[], int rowind[], 
  double a[], string title );
void r8cc_print_some ( int m, int n, int nnzero, int colptr[], int rowind[], 
  double a[], int ilo, int jlo, int ihi, int jhi, string title );
double *r8cc_random ( int m, int n, int nnzero, int colptr[], int rowind[], 
  int *seed );
void r8cc_read ( string col_file, string row_file, string a_file, int m, 
  int n, int nz_num, int col[], int row[], double a[] );
void r8cc_read_size ( string col_file, string row_file, int *m, int *n, 
  int *nz_num, int *base );
void r8cc_set ( int m, int n, int nnzero, int colptr[], int rowind[], 
  double a[], int i, int j, double aij );
double *r8cc_to_r8ge ( int m, int n, int nnzero, int colptr[], int rowind[], 
  double a[] );
double *r8cc_vxm ( int m, int n, int nnzero, int colptr[], int rowind[], 
  double a[], double x[] );
void r8cc_write ( string col_file, string row_file, string a_file, int m, int n,
  int nz_num, int col[], int row[], double a[] );
double *r8cc_zero ( int m, int n, int nnzero, int colptr[], int rowind[] );
//
//  Real double precision Circulant.
//
double *r8ci_eval ( int n, double a[] );
double *r8ci_indicator ( int n );
double *r8ci_mxv ( int n, double a[], double x[] );
void r8ci_print ( int n, double a[], string title );
void r8ci_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8ci_random ( int n, int *seed );
double *r8ci_sl ( int n, double a[], double b[], int job );
double *r8ci_to_r8ge ( int n, double a[] );
double *r8ci_vxm ( int n, double a[], double x[] );
double *r8ci_zero ( int n );
//
//  Real double precision General Banded.
//
double r8gb_det ( int n, int ml, int mu, double a[], int pivot[] );
int r8gb_fa ( int n, int ml, int mu, double a[], int pivot[] );
double *r8gb_indicator ( int m, int n, int ml, int mu );
double *r8gb_ml ( int n, int ml, int mu, double a[], int pivot[], double x[], 
  int job );
double *r8gb_mu ( int n, int ml, int mu, double a[], int pivot[], double x[], 
  int job );
double *r8gb_mxv ( int m, int n, int ml, int mu, double a[], double x[] );
int r8gb_nz_num ( int m, int n, int ml, int mu, double a[] );
void r8gb_print ( int m, int n, int ml, int mu, double a[], string title );
void r8gb_print_some ( int m, int n, int ml, int mu, double a[], int ilo, 
  int jlo, int ihi, int jhi, string title );
double *r8gb_random ( int m, int n, int ml, int mu, int *seed );
double *r8gb_sl ( int n, int ml, int mu, double a[], int pivot[], 
  double b[], int job );
void r8gb_to_r8s3 ( int m, int n, int ml, int mu, double a[], int nz_num, 
  int *isym, int row[], int col[], double b[] );
void r8gb_to_r8sp ( int m, int n, int ml, int mu, double a[], int nz_num, 
  int row[], int col[], double b[] );
double *r8gb_to_r8vec ( int m, int n, int ml, int mu, double *a );
double *r8gb_to_r8ge ( int m, int n, int ml, int mu, double a[] );
void r8gb_to_r8s3 ( int m, int n, int ml, int mu, double a[], int nz_num, 
  int *isym, int row[], int col[], double b[] );
int r8gb_trf ( int m, int n, int ml, int mu, double a[], int pivot[] );
double *r8gb_trs ( int n, int ml, int mu, int nrhs, char trans, double a[], 
  int pivot[], double b[] );
double *r8gb_vxm ( int m, int n, int ml, int mu, double a[], double x[] );
double *r8gb_zero ( int m, int n, int ml, int mu );
//
//  Real double precision General Diagonal.
//
bool r8gd_error ( int n, int ndiag );
double *r8gd_indicator ( int n, int ndiag, int offset[] );
double *r8gd_mxv ( int n, int ndiag, int offset[], double a[], double x[] );
void r8gd_print ( int n, int ndiag, int offset[], double a[], string title );
void r8gd_print_some ( int n, int ndiag, int offset[], double a[], int ilo, 
  int jlo, int ihi, int jhi, string title );
double *r8gd_random ( int n, int ndiag, int offset[], int *seed );
double *r8gd_to_r8ge ( int n, int ndiag, int offset[], double a[] );
double *r8gd_vxm ( int n, int ndiag, int offset[], double a[], double x[] );
double *r8gd_zero ( int n, int ndiag );
//
//  Real double precision General.
//
double r8ge_co ( int n, double a[], int pivot[] );
double r8ge_det ( int n, double a[], int pivot[] );
double *r8ge_dilu ( int m, int n, double a[] );
int r8ge_fa ( int n, double a[], int pivot[] );
void r8ge_fs ( int n, double a[], double x[] );
double *r8ge_fs_new ( int n, double a[], double b[] );
void r8ge_fss ( int n, double a[], int nb, double b[] );
double *r8ge_fss_new ( int n, double a[], int nb, double b[] );
double *r8ge_identity ( int n );
void r8ge_ilu ( int m, int n, double a[], double l[], double u[] );
double *r8ge_indicator ( int m, int n );
double *r8ge_inverse ( int n, double a[], int pivot[] );
double *r8ge_ml ( int n, double a[], int pivot[], double x[], int job );
double *r8ge_mu ( int m, int n, double a[], char trans, int pivot[], double x[] );
double *r8ge_mxm ( int n, double a[], double b[] );
double *r8ge_mxv ( int m, int n, double a[], double x[] );
double r8ge_np_det ( int n, double a[] );
int r8ge_np_fa ( int n, double a[] );
double *r8ge_np_inverse ( int n, double a[] );
double *r8ge_np_ml ( int n, double a[], double x[], int job );
double *r8ge_np_sl ( int n, double a[], double b[], int job );
double *r8ge_np_trm ( int m, int n, double a[], double x[], int job );
int r8ge_np_trf ( int m, int n, double a[] );
double *r8ge_np_trs ( int n, int nrhs, char trans, double a[], double b[] );
void r8ge_plu ( int m, int n, double a[], double p[], double l[], double u[] );
double *r8ge_poly ( int n, double a[] );
void r8ge_print ( int m, int n, double a[], string title );
void r8ge_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8ge_random ( int m, int n, int *seed );
double *r8ge_res ( int m, int n, double a[], double x[], double b[] );
void r8ge_sl ( int n, double a[], int pivot[], double b[], int job );
double *r8ge_sl_it ( int n, double a[], double alu[], int pivot[], double b[], 
  int job, double x[] );
double *r8ge_sl_new ( int n, double a[], int pivot[], double b[], int job );
double *r8ge_to_r8gb ( int m, int n, int ml, int mu, double a[] );
void r8ge_to_r8ri ( int n, double a[], int nz, int ija[], double sa[] );
int r8ge_to_r8ri_size ( int n, double a[] );
double *r8ge_to_r8vec ( int m, int n, double *a );
int r8ge_trf ( int m, int n, double a[], int pivot[] );
double *r8ge_trs ( int n, int nrhs, char trans, double a[], int pivot[], double b[] );
double *r8ge_vxm ( int m, int n, double a[], double x[] );
double *r8ge_zero ( int m, int n );
//
//  Real double precision Lower Triangular, Full Storage.
//
double r8lt_det ( int n, double a[] );
double *r8lt_indicator ( int m, int n );
double *r8lt_inverse ( int n, double a[] );
double *r8lt_mxm ( int n, double a[], double b[] );
double *r8lt_mxv ( int m, int n, double a[], double x[] );
void r8lt_print ( int m, int n, double a[], string title );
void r8lt_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8lt_random ( int m, int n, int *seed );
double *r8lt_sl ( int n, double a[], double b[], int job );
double *r8lt_vxm ( int m, int n, double a[], double x[] );
double *r8lt_zero ( int m, int n );
//
//  R8MAT utilities.
//
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
//
//  Real double precision Nonsymmetric Coordinate Format.
//
double *r8ncf_indicator ( int m, int n, int nz_num, int rowcol[] );
void r8ncf_print ( int m, int n, int nz_num, int rowcol[], double a[], 
  string title );
void r8ncf_print_some ( int m, int n, int nz_num, int rowcol[], 
  double a[], int ilo, int jlo, int ihi, int jhi, string title );
//
//  Real double precision Positive Definite Symmetric Band, Lower.
//
double r8pbl_det ( int n, int mu, double a[] );
double *r8pbl_indicator ( int n, int mu );
void r8pbl_print ( int n, int mu, double a[], string title );
void r8pbl_print_some ( int n, int mu, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8pbl_random ( int n, int mu, int *seed );
double *r8pbl_to_r8ge ( int n, int mu, double a[] );
double *r8pbl_zero ( int n, int mu );
//
//  Real double precision Positive Definite Symmetric Band, Upper.
//
double *r8pbu_cg ( int n, int mu, double a[], double b[], double x[] );
double r8pbu_det ( int n, int mu, double a[] );
double *r8pbu_fa ( int n, int mu, double a[] );
double *r8pbu_indicator ( int n, int mu );
double *r8pbu_ml ( int n, int mu, double a[], double x[] );
double *r8pbu_mxv ( int n, int mu, double a[], double x[] );
void r8pbu_print ( int n, int mu, double a[], string title );
void r8pbu_print_some ( int n, int mu, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8pbu_random ( int n, int mu, int *seed );
double *r8pbu_sl ( int n, int mu, double a[], double b[] );
double *r8pbu_sor ( int n, int mu, double a[], double b[], double eps, int itchk, 
  int itmax, double omega, double x_init[] );
double *r8pbu_to_r8ge ( int n, int mu, double a[] );
double *r8pbu_zero ( int n, int mu );
//
//  Real double precision Positive Definite Symmetric.
//
double r8po_det ( int n, double a[] );
double *r8po_fa ( int n, double a[] );
double *r8po_indicator ( int n );
double *r8po_inverse ( int n, double a[] );
double *r8po_ml ( int n, double a[], double x[] );
double *r8po_mxm ( int n, double a[], double b[] );
double *r8po_mxv ( int n, double a[], double x[] );
void r8po_print ( int n, double a[], string title );
void r8po_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8po_random ( int n, int *seed );
double *r8po_sl ( int n, double a[], double b[] );
double *r8po_to_r8ge ( int n, double a[] );
double *r8po_zero ( int n );
//
//  Real double precision Positive Definite Symmetric Packed.
//
double r8pp_det ( int n, double a[] );
double *r8pp_fa ( int n, double a[] );
double *r8pp_indicator ( int n );
double *r8pp_mxv ( int n, double a[], double x[] );
void r8pp_print ( int n, double a[], string title );
void r8pp_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8pp_random ( int n, int *seed );
double *r8pp_sl ( int n, double a[], double b[] );
double *r8pp_to_r8ge ( int n, double a[] );
double *r8pp_zero ( int n );
//
//  Real double precision Row Indexed.
//
double *r8ri_to_r8ge ( int nz, int ija[], double sa[], int n );
//
//  Utilities for double precision matrices, considered as arrays of rows.
//
void r8row_swap ( int m, int n, double a[], int irow1, int irow2 );
//
//  Real double precision SLAP Triad.
//
double *r8s3_indicator ( int n, int nz_num, int isym, int row[], int col[] );
void r8s3_print ( int m, int n, int nz_num, int isym, int row[], int col[], 
  double a[], string title );
void r8s3_print_some ( int m, int n, int nz_num, int isym, int row[], int col[], 
  double a[], int ilo, int jlo, int ihi, int jhi, string title );
void r8s3_read ( string input_file, int n, int nz_num, int row[], int col[], 
  double a[] );
void r8s3_read_size ( string input_file, int *n, int *nz_num );
void r8s3_write ( int n, int nz_num, int isym, int row[], int col[], 
  double a[], string output_file );
//
//  Real double precision Symmetric Diagonal
//
double *r8sd_cg ( int n, int ndiag, int offset[], double a[], double b[], 
  double x[] );
double *r8sd_indicator ( int n, int ndiag, int offset[] );
double *r8sd_mxv ( int n, int ndiag, int offset[], double a[], double x[] );
void r8sd_print ( int n, int ndiag, int offset[], double a[], string title );
void r8sd_print_some ( int n, int ndiag, int offset[], double a[], int ilo, 
  int jlo, int ihi, int jhi, string title );
double *r8sd_random ( int n, int ndiag, int offset[], int *seed );
double *r8sd_to_r8ge ( int n, int ndiag, int offset[], double a[] );
double *r8sd_zero ( int n, int ndiag );
//
//  Real double precision Sherman Morrison.
//
double *r8sm_ml ( int n, double a[], double u[], double v[], int pivot[], 
  double x[], int job );
double *r8sm_mxv ( int m, int n, double a[], double u[], double v[], double x[] );
void r8sm_print ( int m, int n, double a[], double u[], double v[], string title );
void r8sm_print_some ( int m, int n, double a[], double u[], double v[], int ilo, 
  int jlo, int ihi, int jhi, string title );
void r8sm_random ( int m, int n, double a[], double u[], double v[], int *seed );
double *r8sm_sl ( int n, double a[], double u[], double v[], double b[], 
  int pivot[], int job );
double *r8sm_to_r8ge ( int m, int n, double a[], double u[], double v[] );
double *r8sm_vxm ( int m, int n, double a[], double u[], double v[], double x[] );
void r8sm_zero ( int m, int n, double a[], double u[], double v[] );
//
//  Real double precision Sparse storage.
//
bool r8sp_check ( int m, int n, int nz_num, int row[], int col[] );
int r8sp_ij_to_k ( int nz_num, int row[], int col[], int i, int j );
double *r8sp_indicator ( int m, int n, int nz_num, int row[], int col[] );
double *r8sp_mxv ( int m, int n, int nz_num, int row[], int col[], 
  double a[], double x[] );
void r8sp_print ( int m, int n, int nz_num, int row[], int col[], 
  double a[], string title );
void r8sp_print_some ( int m, int n, int nz_num, int row[], int col[], 
  double a[], int ilo, int jlo, int ihi, int jhi, string title );
double *r8sp_random ( int m, int n, int nz_num, int row[], int col[], 
  int *seed );
void r8sp_read ( string input_file, int m, int n, int nz_num, int row[], 
  int col[], double a[] );
void r8sp_read_size ( string input_file, int *m, int *n, int *nz_num );
double *r8sp_to_r8ge ( int m, int n, int nz_num, int row[], int col[], 
  double a[] );
void r8sp_to_r8ncf ( int m, int n, int nz_num, int row[], int col[], 
  double a[], int rowcol[] );
double *r8sp_vxm ( int m, int n, int nz_num, int row[], int col[], 
  double a[], double x[] );
void r8sp_write ( int m, int n, int nz_num, int row[], int col[], double a[], 
  string output_file );
double *r8sp_zero ( int m, int n, int nz_num, int row[], int col[] );
//
//  Real double precision Sparse Row storage.
//
void r8sr_indicator ( int n, int nz, int row[], int col[], double diag[], 
  double off[] );
double *r8sr_mxv ( int n, int nz, int row[], int col[], double diag[], 
  double off[], double x[] );
void r8sr_print ( int n, int nz, int row[], int col[], double diag[], 
  double off[], string title );
void r8sr_print_some ( int n, int nz, int row[], int col[], double diag[], 
  double off[], int ilo, int jlo, int ihi, int jhi, string title );
void r8sr_random ( int n, int nz, int row[], int col[], double diag[], 
  double off[], int *seed );
double *r8sr_to_r8ge ( int n, int nz, int row[], int col[], double diag[], 
  double off[] );
double *r8sr_vxm ( int n, int nz, int row[], int col[], double diag[], 
  double off[], double x[] );
void r8sr_zero ( int n, int nz, int row[], int col[], double diag[], 
  double off[] );
//
//  Real double precision Symmetric Skyline.
//
bool r8ss_error ( int diag[], int n, int na );
double *r8ss_indicator ( int n, int *na, int diag[]);
double *r8ss_mxv ( int n, int na, int diag[], double a[], double x[] );
void r8ss_print ( int n, int na, int diag[], double a[], string title );
void r8ss_print_some ( int n, int na, int diag[], double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8ss_random ( int n, int *na, int diag[], double a[], int *seed );
double *r8ss_to_r8ge ( int n, int na, int diag[], double a[] );
double *r8ss_zero ( int n, int na, int diag[] );
//
//  Real double precision Symmetric Toeplitz.
//
double *r8sto_indicator ( int n );
double *r8sto_inverse ( int n, double a[] );
double *r8sto_mxv ( int n, double a[], double x[] );
void r8sto_print ( int n, double a[], string title );
void r8sto_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8sto_random ( int n, int *seed );
double *r8sto_sl ( int n, double a[], double b[] );
double *r8sto_to_r8ge ( int n, double a[] );
double *r8sto_yw_sl ( int n, double b[] );
double *r8sto_zero ( int n );
//
//  Real double precision Toeplitz.
//
double *r8to_indicator ( int n );
double *r8to_mxv ( int n, double a[], double x[] );
void r8to_print ( int n, double a[], string title );
void r8to_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8to_random ( int n, int *seed );
double *r8to_sl ( int n, double a[], double b[], int job );
double *r8to_to_r8ge ( int n, double a[] );
double *r8to_vxm ( int n, double a[], double x[] );
double *r8to_zero ( int n );
//
//  Real double precision Upper Triangular.
//
double r8ut_det ( int n, double a[] );
double *r8ut_indicator ( int m, int n );
double *r8ut_inverse ( int n, double a[] );
double *r8ut_mxm ( int n, double a[], double b[] );
double *r8ut_mxv ( int m, int n, double a[], double x[] );
void r8ut_print ( int m, int n, double a[], string title );
void r8ut_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8ut_random ( int m, int n, int *seed );
double *r8ut_sl ( int n, double a[], double b[], int job );
double *r8ut_vxm ( int m, int n, double a[], double x[] );
double *r8ut_zero ( int m, int n );
//
//  Real double precision Upper Triangular Packed.
//
void r8utp_print ( int n, double a[], string title );
void r8utp_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
//
//  R8VEC utilities.
//
double r8vec_dot_product ( int n, double x[], double y[] );
double *r8vec_indicator_new ( int n );
void r8vec_print ( int n, double a[], string title );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title );
double *r8vec_random ( int n, double alo, double ahi, int *seed );
double *r8vec_read ( string input_file, int n );
int r8vec_read_size ( string input_file );
double *r8vec_to_r8cb ( int m, int n, int ml, int mu, double *x );
double *r8vec_to_r8gb ( int m, int n, int ml, int mu, double *x );
double *r8vec_to_r8ge ( int m, int n, double *x );
double *r8vec_uniform_01_new ( int n, int &seed );
void r8vec_write ( string output_filename, int n, double x[] );
//
//  R8VEC2 Utilities.
//
void r8vec2_print_some ( int n, double x1[], double x2[], int max_print, 
  string title );
//
//  Real double precision Vandermonde.
//
double r8vm_det ( int n, double a[] );
double *r8vm_mxv ( int m, int n, double a[], double x[] );
void r8vm_print ( int m, int n, double a[], string title );
void r8vm_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8vm_random ( int m, int n, int &seed );
void r8vm_sl ( int n, double a[], double b[], int job, double x[], int *info );
double *r8vm_sl_new ( int n, double a[], double b[], int job, int *info );
double *r8vm_to_r8ge ( int m, int n, double a[] );
double *r8vm_vxm ( int m, int n, double a[], double x[] );
double *r8vm_zero ( int m, int n );
//
//  Utilities for strings.
//
int s_len_trim ( string s );
//
//  Utility routine for external heap sort.
//
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
//
//  Utility to print the time.
//
void timestamp ( );
