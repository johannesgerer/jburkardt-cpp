float c4_argument ( complex <float> x );
float c4_magnitude ( complex <float> x );
complex <float> c4_sqrt ( complex <float> x );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
int i4_sign ( int i );
int i4_uniform ( int a, int b, int &seed );
int i4_wrap ( int ival, int ilo, int ihi );
float i4int_to_r4int ( int imin, int imax, int i, float rmin, float rmax );
void i4vec_copy ( int n, int a1[], int a2[] );
int *i4vec_indicator0_new ( int n );
int *i4vec_indicator1_new ( int n );
void i4vec_permute ( int n, int p[], int base, int a[] );
void i4vec_print ( int n, int a[], string title );
void i4vec_zero ( int n, int a[] );
int *i4vec_zero_new ( int n );
bool perm_check ( int n, int p[], int base );
int *perm_uniform_new ( int n, int base, int &seed );
float r4_abs ( float x );
float r4_add ( float y, float x );
float r4_atan ( float y, float x );
float r4_big ( );
float r4_cas ( float x );
float r4_ceiling ( float x );
float r4_choose ( int n, int k );
float r4_chop ( int place, float x );
complex <float> r4_csqrt ( float x );
float r4_cube_root ( float x );
float r4_diff ( float x, float y, int n );
int r4_digit ( float x, int idigit );
float r4_epsilon ( );
float r4_epsilon_compute ( );
float r4_exp ( float x );
float r4_factorial ( int n );
float r4_factorial2 ( int n );
float r4_floor ( float x );
float r4_fraction ( int i, int j );
float r4_fractional ( float x );
float r4_huge ( );
bool r4_in_01 ( float a );
bool r4_is_int ( float r );
float r4_log_10 ( float x );
float r4_log_2 ( float x );
float r4_log_b ( float x, float b );
void r4_mant ( float x, int *s, float *r, int *l );
float r4_max ( float x, float y );
float r4_min ( float x, float y );
float r4_mod ( float x, float y );
float r4_modp ( float x, float y );
float r4_mop ( int i );
int r4_nint ( float x );
float r4_normal_01 ( int &seed );
float r4_normal_ab ( float a, float b, int &seed );
float r4_pi ( );
float r4_power ( float r, int p );
float r4_power_fast ( float r, int p, int *mults );
float r4_pythag ( float a, float b );
float r4_reverse_bytes ( float x );
int r4_round_i4 ( float x );
float r4_round2 ( int nplace, float x );
float r4_roundb ( int base, int nplace, float x );
float r4_roundx ( int nplace, float x );
float r4_sign ( float x );
float r4_sign3 ( float x );
bool r4_sign_opposite ( float r1, float r2 );
bool r4_sign_opposite_strict ( float r1, float r2 );
void r4_swap ( float &x, float &y );
void r4_swap3 ( float &x, float &y, float &z );
float r4_tiny ( );
float r4_to_r4_discrete ( float r, float rmin, float rmax, int nr );
void r4_to_dhms ( float r, int *d, int *h, int *m, int *s );
int r4_to_i4 ( float x, float xmin, float xmax, int ixmin, int ixmax );
float r4_uniform_01 ( int &seed );
float r4_uniform_ab ( float a, float b, int &seed );
void r4_unswap3 ( float &x, float &y, float &z );
float r4_walsh_1d ( float x, int digit );
float *r42_cheby ( int n, float alo, float ahi );
void r42_print ( float a[2], string title );
void r42_uniform ( float b, float c, int &seed, float r[] );
void r42poly2_print ( float a, float b, float c, float d, float e, float f );
int r42poly2_type ( float a, float b, float c, float d, float e, float f );
void r42poly2_type_print ( int type );
float *r42vec_max ( int n, float a[] );
float *r42vec_min ( int n, float a[] );
int r42vec_order_type ( int n, float a[] );
void r42vec_part_quick_a ( int n, float a[], int *l, int *r );
void r42vec_permute ( int n, int p[], int base, float a[] );
void r42vec_print ( int n, float a[], string title );
int *r42vec_sort_heap_index_a ( int n, int base, float a[] );
void r42vec_sort_quick_a ( int n, float a[] );
float *r43vec_max ( int n, float a[] );
float *r43vec_min ( int n, float a[] );
void r43vec_part_quick_a ( int n, float a[], int *l, int *r );
void r43vec_sort_quick_a ( int n, float a[] );
float *r4block_expand_linear ( int l, int m, int n, float x[], int lfat, int mfat,
  int nfat );
void r4block_print ( int l, int m, int n, float a[], string title );
float *r4block_zero_new ( int l, int m, int n );
int r4col_compare ( int m, int n, float a[], int i, int j );
float *r4col_duplicates ( int m, int n, int n_unique, int &seed );
int r4col_find ( int m, int n, float a[], float x[] );
int *r4col_first_index ( int m, int n, float a[], float tol );
int r4col_insert ( int n_max, int m, int n, float a[], float x[] );
float *r4col_max ( int m, int n, float a[] );
int *r4col_max_index ( int m, int n, float a[] );
void r4col_max_one ( int m, int n, float a[] );
float *r4col_mean ( int m, int n, float a[] );
float *r4col_min ( int m, int n, float a[] );
int *r4col_min_index ( int m, int n, float a[] );
void r4col_part_quick_a ( int m, int n, float a[], int *l, int *r );
void r4col_permute ( int m, int n, int p[], int base, float a[] );
void r4col_sort_heap_a ( int m, int n, float a[] );
int *r4col_sort_heap_index_a ( int m, int n, int base, float a[] );
void r4col_sort_quick_a ( int m, int n, float a[] );
void r4col_sorted_tol_undex ( int m, int n, float a[], int unique_num,
  float tol, int undx[], int xdnu[] );
int r4col_sorted_tol_unique ( int m, int n, float a[], float tol );
int r4col_sorted_tol_unique_count ( int m, int n, float a[], float tol );
void r4col_sorted_undex ( int m, int n, float a[],
  int unique_num, int undx[], int xdnu[] );
int r4col_sorted_unique ( int m, int n, float a[] );
int r4col_sorted_unique_count ( int m, int n, float a[] );
void r4col_sortr_a ( int m, int n, float a[], int key );
float *r4col_sum ( int m, int n, float a[] );
void r4col_swap ( int m, int n, float a[], int j1, int j2 );
float *r4col_to_r4vec ( int m, int n, float a[] );
void r4col_tol_undex ( int x_dim, int x_num, float x_val[], int x_unique_num,
  float tol, int undx[], int xdnu[] );
int r4col_tol_unique_count ( int m, int n, float a[], float tol );
int *r4col_tol_unique_index ( int m, int n, float a[], float tol );
void r4col_undex ( int x_dim, int x_num, float x_val[], int x_unique_num,
  int undx[], int xdnu[] );
int r4col_unique_count ( int m, int n, float a[] );
int *r4col_unique_index ( int m, int n, float a[] );
float *r4col_variance ( int m, int n, float a[] );
float r4int_to_r4int ( float rmin, float rmax, float r, float r2min,
  float r2max );
int r4int_to_i4int ( float rmin, float rmax, float r, int imin, int imax );
float *r4mat_border_add ( int m, int n, float table[] );
float *r4mat_border_cut ( int m, int n, float table[] );;
float *r4mat_cholesky_factor ( int n, float a[] );
float *r4mat_cholesky_solve ( int n, float a[], float b[] );
float *r4mat_choresky_factor ( int n, float a[] );
void r4mat_copy ( int m, int n, float a1[], float a2[] );
float *r4mat_copy_new ( int m, int n, float a1[] );
void r4mat_delete ( float **a, int m, int n );
float r4mat_det ( int n, float a[] );
float r4mat_det_2d ( float a[] );
float r4mat_det_3d ( float a[] );
float r4mat_det_4d ( float a[] );
float r4mat_det_5d ( float a[] );
void r4mat_diag_add_scalar ( int n, float a[], float s );
void r4mat_diag_add_vector ( int n, float a[], float v[] );
float *r4mat_diag_get_vector ( int n, float a[] );
void r4mat_diag_set_scalar ( int n, float a[], float s );
void r4mat_diag_set_vector ( int n, float a[], float v[] );
float *r4mat_expand_linear ( int m, int n, float x[], int mfat, int nfat );
float *r4mat_expand_linear2 ( int m, int n, float a[], int m2, int n2 );
void r4mat_flip_cols ( int m, int n, float a[] );
void r4mat_flip_rows ( int m, int n, float a[] );
float *r4mat_givens_post ( int n, float a[], int row, int col );
float *r4mat_givens_pre ( int n, float a[], int row, int col );
float *r4mat_hess ( float (*fx )( int n, float x[] ), int n, float x[] );
void r4mat_house_axh ( int n, float a[], float v[] );
float *r4mat_house_axh_new ( int n, float a[], float v[] );
float *r4mat_house_form ( int n, float v[] );
float *r4mat_house_hxa ( int n, float a[], float v[] );
float *r4mat_house_post ( int n, float a[], int row, int col );
float *r4mat_house_pre ( int n, float a[], int row, int col );
float *r4mat_identity ( int n );
bool r4mat_in_01 ( int m, int n, float a[] );
float *r4mat_indicator_new ( int m, int n );
float *r4mat_inverse_2d ( float a[] );
float *r4mat_inverse_3d ( float a[] );
float *r4mat_inverse_4d ( float a[] );
float *r4mat_jac ( int m, int n, float eps,
  float *(*fx) ( int m, int n, float x[] ), float x[] );
float *r4mat_l_inverse ( int n, float a[] );
void r4mat_l_print ( int m, int n, float a[], string title );
float *r4mat_l_solve ( int n, float a[], float b[] );
float *r4mat_l1_inverse ( int n, float a[] );
float *r4mat_lt_solve ( int n, float a[], float b[] );
void r4mat_lu ( int m, int n, float a[], float l[], float p[], float u[] );
float r4mat_max ( int m, int n, float a[] );
void r4mat_max_index ( int m, int n, float a[], int *i_max, int *j_max );
float r4mat_maxcol_minrow ( int m, int n, float a[] );
float r4mat_maxrow_mincol ( int m, int n, float a[] );
float r4mat_min ( int m, int n, float a[] );
void r4mat_min_index ( int m, int n, float a[], int *i_min, int *j_min );
float r4mat_mincol_maxrow ( int m, int n, float a[] );
float r4mat_minrow_maxcol ( int m, int n, float a[] );
void r4mat_mm ( int n1, int n2, int n3, float a[], float b[], float c[] );
float *r4mat_mm_new ( int n1, int n2, int n3, float a[], float b[] );
float *r4mat_mtv ( int m, int n, float a[], float x[] );
void r4mat_mtxv ( int m, int n, float a[], float x[], float atx[] );
float *r4mat_mv ( int m, int n, float a[], float x[] );
void r4mat_mxm ( int n1, int n2, int n3, float a[], float b[], float c[] );
float *r4mat_mxm_new ( int n1, int n2, int n3, float a[], float b[] );
void r4mat_mxv ( int m, int n, float a[], float x[], float ax[] );
float **r4mat_new ( int m, int n );
void r4mat_nint ( int m, int n, float a[] );
float r4mat_norm_eis ( int m, int n, float a[] );
float r4mat_norm_fro ( int m, int n, float a[] );
float *r4mat_nullspace ( int m, int n, float a[], int nullspace_size );
int r4mat_nullspace_size ( int m, int n, float a[] );
float *r4mat_orth_uniform_new ( int n, int &seed );
void r4mat_plot ( int m, int n, float a[], string title );
char r4mat_plot_symbol ( float r );
float *r4mat_poly_char ( int n, float a[] );
float *r4mat_power ( int n, float a[], int npow );
void r4mat_power_method ( int n, float a[], float *r, float v[] );
void r4mat_print ( int m, int n, float a[], string title );
void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r4mat_ref ( int m, int n, float a[] );
void r4mat_rref ( int m, int n, float a[] );
int r4mat_solve ( int n, int nrhs, float a[] );
float *r4mat_solve_2d ( float a[], float b[], float *det );
float *r4mat_solve_3d ( float a[], float b[], float *det );
float *r4mat_solve2 ( int n, float a[], float b[], int *ierror );
float *r4mat_symm_eigen ( int n, float x[], float q[] );
void r4mat_symm_jacobi ( int n, float a[] );
int r4mat_to_r4plu ( int n, float a[], int pivot[], float lu[] );
float r4mat_trace ( int n, float a[] );
float *r4mat_transpose ( int m, int n, float a[] );
void r4mat_transpose_in_place ( int n, float a[] );
void r4mat_transpose_print ( int m, int n, float a[], string title );
void r4mat_transpose_print_some ( int m, int n, float a[], int ilo, int jlo,
  int ihi, int jhi, string title );
float *r4mat_u_inverse ( int n, float a[] );
float *r4mat_u1_inverse ( int n, float a[] );
float *r4mat_uniform_new ( int m, int n, float b, float c, int &seed );
void r4mat_uniform_01 ( int m, int n, int &seed, float r[] );
float *r4mat_uniform_01_new ( int m, int n, int &seed );
float *r4mat_vand2 ( int n, float x[] );
void r4mat_zero ( int m, int n, float a[] );
float *r4mat_zero_new ( int m, int n );
float r4plu_det ( int n, int pivot[], float lu[] );
void r4plu_inverse ( int n, int pivot[], float lu[], float a_inverse[] );
void r4plu_mul ( int n, int pivot[], float lu[], float x[], float b[] );
void r4plu_sol ( int n, int pivot[], float lu[], float b[], float x[] );
void r4plu_to_r4mat ( int n, int pivot[], float lu[], float a[] );
int r4poly_degree ( int na, float a[] );
float *r4poly_deriv ( int n, float c[], int p );
float r4poly_lagrange_0 ( int npol, float xpol[], float xval );
float r4poly_lagrange_1 ( int npol, float xpol[], float xval );
float r4poly_lagrange_2 ( int npol, float xpol[], float xval );
float *r4poly_lagrange_coef ( int npol, int ipol, float xpol[] );
void r4poly_lagrange_factor ( int npol, float xpol[], float xval,
  float *wval, float *dwdx );
int r4poly_lagrange_val ( int npol, int ipol, float xpol[], float xval,
  float *pval, float *dpdx );
int r4poly_order ( int na, float a[] );
void r4poly_print ( int n, float a[], string title );
void r4poly_shift ( float scale, float shift, int n, float poly_cof[] );
float r4poly_val_horner ( int n, float c[], float x );
float r4poly_value ( int n, float a[], float x );
int r4poly2_ex ( float x1, float y1, float x2, float y2, float x3, float y3,
  float *x, float *y );
int r4poly2_ex2 ( float x1, float y1, float x2, float y2, float x3, float y3,
  float *x, float *y, float *a, float *b, float *c );
void r4poly2_root ( float a, float b, float c, complex <float> *r1,
  complex <float> *r2 );
void r4poly2_rroot ( float a, float b, float c, float *r1, float *r2 );
void r4poly2_val ( float x1, float y1, float x2, float y2,
  float x3, float y3, float x, float *y, float *yp, float *ypp );
void r4poly2_val2 ( int ndata, float tdata[],
  float ydata[], int left, float tval, float *yval );
void r4poly3_root ( float a, float b, float c, float d, complex <float> *r1,
  complex <float> *r2, complex <float> *r3 );
void r4poly4_root ( float a, float b, float c, float d, float e,
  complex <float> *r1, complex <float> *r2, complex <float> *r3,
  complex <float> *r4 );
void r4pp_delete ( float **a, int m, int n );
float **r4pp_new ( int m, int n );
int r4r4_compare ( float x1, float y1, float x2, float y2 );
void r4r4_print ( float a1, float a2, string title );
int r4r4r4_compare ( float x1, float y1, float z1, float x2, float y2,
  float z2 );
void r4r4r4vec_index_insert_unique ( int maxn, int *n, float x[], float y[],
  float z[], int indx[], float xval, float yval, float zval, int *ival,
  int *ierror );
void r4r4r4vec_index_search ( int n, float x[], float y[], float z[],
  int indx[], float xval, float yval, float zval, int *less, int *equal,
  int *more );
void r4r4vec_index_insert_unique ( int maxn, int *n, float x[], float y[],
  int indx[], float xval, float yval, int *ival, int *ierror );
void r4r4vec_index_search ( int n, float x[], float y[], int indx[],
  float xval, float yval, int *less, int *equal, int *more );
float *r4row_max ( int m, int n, float a[] );
float *r4row_mean ( int m, int n, float a[] );
float *r4row_min ( int m, int n, float a[] );
float *r4row_sum ( int m, int n, float a[] );
void r4row_swap ( int m, int n, float a[], int irow1, int irow2 );
float *r4row_to_r4vec ( int m, int n, float a[] );
float *r4row_variance ( int m, int n, float a[] );
void r4slmat_print ( int m, int n, float a[], string title );
void r4vec_01_to_ab ( int n, float a[], double amax, double amin );
void r4vec_ab_to_01 ( int n, float a[] );
float *r4vec_ab_to_cd ( int n, float a[], float bmin, float bmax );
float r4vec_amax ( int n, float a[] );
int r4vec_amax_index ( int n, float a[] );
float r4vec_amin ( int n, float a[] );
int r4vec_amin_index ( int n, float a[] );
float *r4vec_any_normal ( int dim_num, float v1[] );
void r4vec_bracket ( int n, float x[], float xval, int *left,
  int *right );
void r4vec_bracket2 ( int n, float x[], float xval, int start, int *left,
  int *right );
void r4vec_bracket3 ( int n, float t[], float tval, int *left );
void r4vec_bracket4 ( int nt, float t[], int ns, float s[], int left[] );
float r4vec_circular_variance ( int n, float x[] );
int r4vec_compare ( int n, float a[], float b[] );
float *r4vec_convolve_circ ( int n, float x[], float y[] );
void r4vec_copy ( int n, float a1[], float a2[] );
float *r4vec_copy_new ( int n, float a1[] );
float r4vec_correlation ( int n, float x[], float y[] );
float r4vec_covar ( int n, float x[], float y[] );
float r4vec_cross_product_2d ( float v1[2], float v2[2] );
float r4vec_cross_product_affine_2d ( float v0[2], float v1[2], float v2[2] );
float *r4vec_cross_product_3d ( float v1[3], float v2[3] );
float *r4vec_cross_product_affine_3d ( float v0[3], float v1[3], float v2[3] );
float *r4vec_dif ( int n, float h );
float r4vec_diff_norm ( int n, float a[], float b[] );
float r4vec_diff_norm_l1 ( int n, float a[], float b[] );
float r4vec_diff_norm_l2 ( int n, float a[], float b[] );
float r4vec_diff_norm_li ( int n, float a[], float b[] );
float r4vec_diff_norm_squared ( int n, float a[], float b[] );
void r4vec_direct_product ( int factor_index, int factor_order,
  float factor_value[], int factor_num, int point_num, float x[] );
void r4vec_direct_product2 ( int factor_index, int factor_order,
  float factor_value[], int factor_num, int point_num, float w[] );
float r4vec_distance ( int dim_num, float v1[], float v2[] );
bool r4vec_distinct ( int n, float x[] );
void r4vec_divide ( int n, float a[], float s );
float r4vec_dot_product ( int n, float v1[], float v2[] );
float r4vec_dot_product_affine ( int n, float v0[], float v1[], float v2[] );
bool r4vec_eq ( int n, float a1[], float a2[] );
void r4vec_even ( int n, float alo, float ahi, float a[] );
float *r4vec_even_new ( int n, float alo, float ahi );
float r4vec_even_select ( int n, float xlo, float xhi, int ival );
void r4vec_even2 ( int maxval, int nfill[], int nold, float xold[],
  int *nval, float xval[] );
void r4vec_even3 ( int nold, int nval, float xold[], float xval[] );
float *r4vec_expand_linear ( int n, float x[], int fat );
int *r4vec_first_index ( int n, float a[], float tol );
float r4vec_frac ( int n, float a[], int k );
float *r4vec_fraction ( int n, float x[] );
bool r4vec_gt ( int n, float a1[], float a2[] );
void r4vec_heap_a ( int n, float a[] );
void r4vec_heap_d ( int n, float a[] );
int *r4vec_histogram ( int n, float a[], float a_lo, float a_hi, int histo_num );
float *r4vec_house_column ( int n, float a[], int k );
float r4vec_i4vec_dot_product ( int n, float r4vec[], int i4vec[] );
bool r4vec_in_01 ( int n, float a[] );
void r4vec_index_delete_all ( int n, float x[], int indx[], float xval,
  int *n2, float x2[], int indx2[] );
void r4vec_index_delete_dupes ( int n, float x[], int indx[],
  int *n2, float x2[], int indx2[] );
void r4vec_index_delete_one ( int n, float x[], int indx[], float xval,
  int *n2, float x2[], int indx2[] );
void r4vec_index_insert ( int *n, float x[], int indx[], float xval );
void r4vec_index_insert_unique ( int *n, float x[], int indx[], float xval );
void r4vec_index_order ( int n, float x[], int indx[] );
void r4vec_index_search ( int n, float x[], int indx[], float xval, int *less,
  int *equal, int *more );
void r4vec_index_sort_unique ( int n, float x[], int *n2, float x2[],
  int indx2[] );
void r4vec_index_sorted_range ( int n, float r[], int indx[], float r_lo,
  float r_hi, int *i_lo, int *i_hi );
void r4vec_indexed_heap_d ( int n, float a[], int indx[] );
int r4vec_indexed_heap_d_extract ( int *n, float a[], int indx[] );
void r4vec_indexed_heap_d_insert ( int *n, float a[], int indx[],
  int indx_insert );
int r4vec_indexed_heap_d_max ( int n, float a[], int indx[] );
void r4vec_indicator0 ( int n, float a[] );
float *r4vec_indicator0_new ( int n );
void r4vec_indicator1 ( int n, float a[] );
float *r4vec_indicator1_new ( int n );
void r4vec_insert ( int n, float a[], int pos, float value );
bool r4vec_is_int ( int n, float a[] );
bool r4vec_is_nonnegative ( int n, float x[] );
bool r4vec_is_zero ( int n, float x[] );
float *r4vec_linspace_new ( int n, float a_lo, float a_hi );
bool r4vec_lt ( int n, float a1[], float a2[] );
void r4vec_mask_print ( int n, float a[], int mask_num, int mask[],
  string title );
float r4vec_max ( int n, float x[] );
int r4vec_max_index ( int n, float a[] );
float r4vec_mean ( int n, float x[] );
float r4vec_median ( int n, float a[] );
float r4vec_min ( int n, float r4vec[] );
int r4vec_min_index ( int n, float a[] );
float r4vec_min_pos ( int n, float a[] );
bool r4vec_mirror_next ( int n, float a[] );
bool r4vec_negative_strict ( int n, float a[] );
float *r4vec_nint ( int n, float a[] );
float r4vec_norm ( int n, float v[] );
float r4vec_norm_affine ( int n, float v0[], float v1[] );
float r4vec_norm_l1 ( int n, float a[] );
float r4vec_norm_l2 ( int n, float a[] );
float r4vec_norm_li ( int n, float a[] );
float r4vec_norm_lp ( int n, float a[], float p );
float *r4vec_normal_01 ( int n, int &seed );
void r4vec_normalize ( int n, float a[] );
void r4vec_normalize_l1 ( int n, float a[] );
float r4vec_normsq ( int n, float a[] );
float r4vec_normsq_affine ( int n, float v0[], float v1[] );
float *r4vec_ones_new ( int n );
int r4vec_order_type ( int n, float x[] );
void r4vec_part_quick_a ( int n, float a[], int *l, int *r );
void r4vec_permute ( int n, int p[], int base, float a[] );
void r4vec_permute_cyclic ( int n, int k, float a[] );
void r4vec_permute_uniform ( int n, float a[], int &seed );
void r4vec_polarize ( int n, float a[], float p[], float a_normal[],
  float a_parallel[] );
bool r4vec_positive_strict ( int n, float a[] );
void r4vec_print ( int n, float a[], string title );
void r4vec_print_part ( int n, float a[], int max_print, string title );
void r4vec_print_some ( int n, float a[], int i_lo, int i_hi, string title );
float r4vec_product ( int n, float a[] );
void r4vec_range ( int n, float x[], float xmin, float xmax, float y[],
  float *ymin, float *ymax );
void r4vec_range_2 ( int n, float a[], float *amin, float *amax );
void r4vec_reverse ( int n, float a[] );
void r4vec_rotate ( int n, float a[], int m );
float r4vec_scalar_triple_product ( float v1[3], float v2[3], float v3[3] );
int r4vec_search_binary_a ( int n, float a[], float aval );
void r4vec_shift ( int shift, int n, float x[] );
void r4vec_shift_circular ( int shift, int n, float x[] );
void r4vec_sort_bubble_a ( int n, float a[] );
void r4vec_sort_bubble_d ( int n, float a[] );
void r4vec_sort_heap_a ( int n, float a[] );
void r4vec_sort_heap_d ( int n, float a[] );
void r4vec_sort_heap_index_a ( int n, float a[], int indx[] );
int *r4vec_sort_heap_index_a_new ( int n, float a[] );
void r4vec_sort_heap_index_d ( int n, float a[], int indx[] );
int *r4vec_sort_heap_index_d_new ( int n, float a[] );
int *r4vec_sort_heap_mask_a ( int n, float a[], int mask_num, int mask[] );
void r4vec_sort_insert_a ( int n, float a[] );
int *r4vec_sort_insert_index_a ( int n, float a[] );
void r4vec_sort_quick_a ( int n, float a[] );
void r4vec_sort_shell_a ( int n, float a[] );
float *r4vec_sorted_merge_a ( int na, float a[], int nb, float b[], int *nc );
int r4vec_sorted_nearest ( int n, float a[], float value );
void r4vec_sorted_range ( int n, float r[], float r_lo, float r_hi,
  int *i_lo, int *i_hi );
void r4vec_sorted_split ( int n, float a[], float split, int *i_lt, int *i_gt );
void r4vec_sorted_undex ( int x_num, float x_val[], int x_unique_num,
  float tol, int undx[], int xdnu[] );
float *r4vec_sorted_unique ( int n, float a[], float tol, int *unique_num );
int r4vec_sorted_unique_count ( int n, float a[], float tol );
void r4vec_sorted_unique_hist ( int n, float a[], float tol, int maxuniq,
  int *unique_num, float auniq[], int acount[] );
int r4vec_split ( int n, float a[], float split );
float r4vec_std ( int n, float a[] );
void r4vec_stutter ( int n, float a[], int m, float am[] );
float *r4vec_stutter_new ( int n, float a[], int m );
float r4vec_sum ( int n, float a[] );
void r4vec_swap ( int n, float a1[], float a2[] );
void r4vec_transpose_print ( int n, float a[], string title );
void r4vec_undex ( int x_num, float x_val[], int x_unique_num, float tol,
  int undx[], int xdnu[] );
void r4vec_uniform ( int n, float b, float c, int &seed, float r[] );
float *r4vec_uniform_new ( int n, float b, float c, int &seed );
void r4vec_uniform_01 ( int n, int &seed, float r[] );
float *r4vec_uniform_01_new ( int n, int &seed );
int r4vec_unique_count ( int n, float a[], float tol );
int *r4vec_unique_index ( int n, float a[], float tol );
void r4vec_unit_sum ( int n, float a[] );
float r4vec_variance ( int n, float x[] );
float *r4vec_vector_triple_product ( float v1[3], float v2[3], float v3[3] );
void r4vec_write ( int n, float r[], string output_file );
void r4vec_zero ( int n, float a1[] );
float *r4vec_zero_new ( int n );
int r4vec2_compare ( int n, float a1[], float a2[], int i, int j );
void r4vec2_print ( int n, float a1[], float a2[], string title );
void r4vec2_print_some ( int n, float x1[], float x2[], int max_print,
  string title );
void r4vec2_sort_a ( int n, float a1[], float a2[] );
void r4vec2_sort_d ( int n, float a1[], float a2[] );
int *r4vec2_sort_heap_index_a ( int n, int base, float x[], float y[] );
void r4vec2_sorted_unique ( int n, float a1[], float a2[], int *unique_num );
void r4vec2_sorted_unique_index ( int n, float a1[], float a2[],
  int *unique_num, int indx[] );
int r4vec2_sum_max_index ( int n, float a[], float b[] );
void r4vec3_print ( int n, float a1[], float a2[], float a3[], string title );
float *roots_to_r4poly ( int n, float x[] );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void timestamp ( );
