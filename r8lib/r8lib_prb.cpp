# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "r8lib.hpp"

int main ( );

void i4int_to_r8int_test ( );

void perm0_check_test ( );
void perm0_uniform_test ( );

void perm1_check_test ( );
void perm1_uniform_test ( );

void r8_abs_test ( );
void r8_acos_test ( );
void r8_acosh_test ( );
void r8_asinh_test ( );
void r8_atan_test ( );
void r8_atanh_test ( );
void r8_big_test ( );
void r8_cas_test ( );
void r8_ceiling_test ( );
void r8_choose_test ( );
void r8_cosd_test ( );
void r8_cotd_test ( );
void r8_cscd_test ( );
void r8_cube_root_test ( );
void r8_diff_test ( );
void r8_digit_test ( );
void r8_e_test ( );
void r8_epsilon_test ( );
void r8_epsilon_compute_test ( );
void r8_factorial_test ( );
void r8_factorial2_test ( );
void r8_fall_test ( );
void r8_fractional_test ( );
void r8_gamma_test ( );
void r8_gamma_log_test ( );
void r8_huge_test ( );
void r8_log_2_test ( );
void r8_log_b_test ( );
void r8_mant_test ( );
void r8_max_test ( );
void r8_min_test ( );
void r8_mod_test ( );
void r8_modp_test ( );
void r8_mop_test ( );
void r8_nint_test ( );
void r8_normal_01_test ( );
void r8_pi_test ( );
void r8_power_test ( );
void r8_power_fast_test ( );
void r8_rise_test ( );
void r8_round2_test ( );
void r8_roundb_test ( );
void r8_roundx_test ( );
void r8_secd_test ( );
void r8_sign_test ( );
void r8_sign3_test ( );
void r8_sind_test ( );
void r8_swap_test ( );
void r8_swap3_test ( );
void r8_tand_test ( );
void r8_to_i4_test ( );
void r8_to_r8_discrete_test ( );
void r8_uniform_01_test ( );
void r8_uniform_ab_test ( );
void r8_walsh_1d_test ( );
void r8_wrap_test ( );

void r82col_print_part_test ( );

void r82poly2_type_test ( );

void r82row_order_type_test( );
void r82row_part_quick_a_test ( );
void r82row_print_part_test ( );
void r82row_sort_heap_index_a_test ( );
void r82row_sort_quick_a_test ( );

void r83col_print_part_test ( );

void r83row_print_part_test ( );

void r8block_expand_linear_test ( );
void r8block_new_test ( );
void r8block_print_test ( );

void r8cmat_to_r8mat_new_test ( );

void r8col_find_test ( );
void r8col_insert_test ( );
void r8col_sort_heap_a_test ( );
void r8col_sort_heap_index_a_test ( );
void r8col_sort_quick_a_test ( );
void r8col_sorted_tol_unique_test ( );
void r8col_sorted_unique_count_test ( );
void r8col_sorted_tol_undex_test ( );
void r8col_max_test ( );
void r8col_mean_test ( );
void r8col_min_test ( );
void r8col_permute_test ( );
void r8col_sortr_a_test ( );
void r8col_sum_test ( );
void r8col_swap_test ( );
void r8col_to_r8vec_test ( );
void r8col_tol_undex_test ( );
void r8col_undex_test ( );
void r8col_unique_count_test ( );
void r8col_variance_test ( );

void r8int_to_i4int_test ( );

void r8mat_cholesky_inverse_test ( );
void r8mat_cholesky_solve_test ( );
void r8mat_cholesky_solve_upper_test ( );
void r8mat_det_2d_test ( );
void r8mat_det_3d_test ( );
void r8mat_det_4d_test ( );
void r8mat_det_5d_test ( );
void r8mat_expand_linear_test ( );
void r8mat_expand_linear2_test ( );
void r8mat_fs_new_test ( );
void r8mat_fss_new_test ( );
void r8mat_givens_post_test ( );
void r8mat_givens_pre_test ( );
void r8mat_hess_test ( );
double r8mat_hess_f ( int n, double x[] );
double *r8mat_hess_exact ( int n, double x[] );
void r8mat_house_axh_test ( );
void r8mat_house_form_test ( );
void r8mat_house_post_test ( );
void r8mat_house_pre_test ( );
void r8mat_indicator_new_test ( );
void r8mat_inverse_2d_test ( );
void r8mat_inverse_3d_test ( );
void r8mat_inverse_4d_test ( );
void r8mat_jac_test ( );
double *r8mat_jac_f ( int m, int n, double x[] );
double *r8mat_jac_exact ( int m, int n, double x[] );
void r8mat_kronecker_test ( );
void r8mat_l_inverse_test ( );
void r8mat_l_print_test ( );
void r8mat_l1_inverse_test ( );
void r8mat_lu_test ( );
void r8mat_max_test ( );
void r8mat_max_index_test ( );
void r8mat_maxcol_minrow_test ( );
void r8mat_maxrow_mincol_test ( );
void r8mat_min_test ( );
void r8mat_min_index_test ( );
void r8mat_mincol_maxrow_test ( );
void r8mat_minrow_maxcol_test ( );
void r8mat_mm_test ( );
void r8mat_mm_new_test ( );
void r8mat_mv_new_test ( );
void r8mat_mv_test ( );
void r8mat_mtv_new_test ( );
void r8mat_mtv_test ( );
void r8mat_nint_test ( );
void r8mat_nonzeros_test ( );
void r8mat_norm_fro_test ( );
void r8mat_norm_l1_test ( );
void r8mat_nullspace_test ( );
void r8mat_nullspace_size_test ( );
void r8mat_orth_uniform_new_test ( );
void r8mat_plot_test ( );
void r8mat_power_method_test ( );
void r8mat_print_test ( );
void r8mat_print_some_test ( );
void r8mat_ref_test ( );
void r8mat_rref_test ( );
void r8mat_solve_test ( );
void r8mat_solve_2d_test ( );
void r8mat_solve_3d_test ( );
void r8mat_solve2_test ( );
void r8mat_sub_new_test ( );
void r8mat_symm_jacobi_test ( );
void r8mat_to_r8cmat_new_test ( );
void r8mat_to_r8plu_test ( );
void r8mat_to_r8rmat_test ( );
void r8mat_trace_test ( );
void r8mat_transpose_new_test ( );
void r8mat_transpose_print_test ( );
void r8mat_u_inverse_test ( );
void r8mat_u1_inverse_test ( );
void r8mat_uniform_ab_new_test ( );

void r8plu_det_test ( );
void r8plu_inverse_test ( );
void r8plu_mul_test ( );
void r8plu_sol_test ( );
void r8plu_to_r8mat_test ( );

void r8poly_degree_test ( );
void r8poly_deriv_test ( );
void r8poly_lagrange_coef_test ( );
void r8poly_lagrange_0_test ( );
void r8poly_lagrange_1_test ( );
void r8poly_lagrange_2_test ( );
void r8poly_lagrange_factor_test ( );
void r8poly_lagrange_val_test ( );
void r8poly_print_test ( );
void r8poly_value_horner_test ( );
void r8poly_values_horner_test ( );

void r8poly2_ex_test ( );
void r8poly2_ex2_test ( );
void r8poly2_val_test ( );
void r8poly2_val_f ( double x, double *y, double *yp, double *ypp );
void r8poly2_val2_test ( );

void r8rmat_new_test ( );
void r8rmat_to_r8mat_test ( );

void r8row_max_test ( );
void r8row_mean_test ( );
void r8row_min_test ( );
void r8row_sum_test ( );
void r8row_swap_test ( );
void r8row_to_r8vec_test ( );
void r8row_variance_test ( );

void r8r8vec_index_insert_unique_test ( );

void r8r8r8vec_index_insert_unique_test ( );

void r8slmat_print_test ( );

void r8vec_amax_test ( );
void r8vec_amin_test ( );
void r8vec_bracket_test ( );
void r8vec_bracket2_test ( );
void r8vec_bracket3_test ( );
void r8vec_bracket5_test ( );
void r8vec_chebyspace_new_test ( );
void r8vec_concatenate_new_test ( );
void r8vec_convolution_test ( );
void r8vec_convolution_circ_test ( );
void r8vec_dif_test ( );
double r8vec_dif_f ( double x );
void r8vec_direct_product_test ( );
void r8vec_direct_product2_test ( );
void r8vec_even_test ( );
void r8vec_even2_test ( );
void r8vec_expand_linear_test ( );
void r8vec_frac_test ( );
void r8vec_histogram_test ( );
void r8vec_house_column_test ( );
void r8vec_index_delete_all_test ( );
void r8vec_index_delete_dupes_test ( );
void r8vec_index_delete_one_test ( );
void r8vec_index_insert_test ( );
void r8vec_index_insert_unique_test ( );
void r8vec_index_order_test ( );
void r8vec_index_search_test ( );
void r8vec_index_sorted_range_test ( );
void r8vec_indexed_heap_d_test ( );
void r8vec_indexed_heap_d_extract_test ( );
void r8vec_indexed_heap_d_insert_test ( );
void r8vec_indexed_heap_d_max_test ( );
void r8vec_indicator0_new_test ( );
void r8vec_legendre_test ( );
void r8vec_linspace_new_test ( );
void r8vec_max_test ( );
void r8vec_max_index_test ( );
void r8vec_mean_test ( );
void r8vec_median_test ( );
void r8vec_midspace_new_test ( );
void r8vec_min_test ( );
void r8vec_min_index_test ( );
void r8vec_nint_test ( );
void r8vec_norm_l0_test ( );
void r8vec_norm_l1_test ( );
void r8vec_norm_l2_test ( );
void r8vec_norm_li_test ( );
void r8vec_normal_01_test ( );
void r8vec_normalize_l1_test ( );
void r8vec_order_type_test ( );
void r8vec_permute_test ( );
void r8vec_permute_uniform_test ( );
void r8vec_polarize_test ( );
void r8vec_print_test ( );
void r8vec_rotate_test ( );
void r8vec_reverse_test ( );
void r8vec_search_binary_a_test ( );
void r8vec_sort_bubble_a_test ( );
void r8vec_sort_heap_a_test ( );
void r8vec_sort_heap_d_test ( );
void r8vec_sort_heap_index_a_new_test ( );
void r8vec_sort_heap_index_d_new_test ( );
void r8vec_sort_heap_mask_a_test ( );
void r8vec_sort_insert_a_test ( );
void r8vec_sort_insert_index_a_test ( );
void r8vec_sort_quick_a_test ( );
void r8vec_sorted_merge_a_test ( );
void r8vec_sorted_nearest_test ( );
void r8vec_sorted_range_test ( );
void r8vec_sorted_split_test ( );
void r8vec_sorted_undex_test ( );
void r8vec_sorted_unique_test ( );
void r8vec_sorted_unique_count_test ( );
void r8vec_sorted_unique_hist_test ( );
void r8vec_split_test ( );
void r8vec_transpose_print_test ( );
void r8vec_undex_test ( );
void r8vec_uniform_01_new_test ( );
void r8vec_uniform_ab_new_test ( );
void r8vec_variance_test ( );

void r8vec2_sort_a_test ( );
void r8vec2_sort_d_test ( );
void r8vec2_sort_heap_index_a_test ( );
void r8vec2_sorted_unique_test ( );
void r8vec2_sorted_unique_index_test ( );
void r8vec2_sum_max_index_test ( );

void roots_to_r8poly_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for R8LIB_PRB.
//
//  Discussion:
//
//    R8LIB_PRB tests the R8LIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "R8LIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the R8LIB library.\n";

  i4int_to_r8int_test ( );

  perm0_check_test ( );
  perm0_uniform_test ( );

  perm1_check_test ( );
  perm1_uniform_test ( );

  r8_abs_test ( );
  r8_acos_test ( );
  r8_acosh_test ( );
  r8_asinh_test ( );
  r8_atan_test ( );
  r8_atanh_test ( );
  r8_big_test ( );
  r8_cas_test ( );
  r8_ceiling_test ( );
  r8_choose_test ( );
  r8_cosd_test ( );
  r8_cotd_test ( );
  r8_cscd_test ( );
  r8_cube_root_test ( );
  r8_diff_test ( );
  r8_digit_test ( );
  r8_e_test ( );
  r8_epsilon_test ( );
  r8_epsilon_compute_test ( );
  r8_factorial_test ( );
  r8_factorial2_test ( );
  r8_fall_test ( );
  r8_fractional_test ( );
  r8_gamma_test ( );
  r8_gamma_log_test ( );
  r8_huge_test ( );
  r8_log_2_test ( );
  r8_log_b_test ( );
  r8_mant_test ( );
  r8_max_test ( );
  r8_min_test ( );
  r8_mod_test ( );
  r8_modp_test ( );
  r8_mop_test ( );
  r8_nint_test ( );
  r8_normal_01_test ( );
  r8_pi_test ( );
  r8_power_test ( );
  r8_power_fast_test ( );
  r8_rise_test ( );
  r8_round2_test ( );
  r8_roundb_test ( );
  r8_roundx_test ( );
  r8_secd_test ( );
  r8_sign_test ( );
  r8_sign3_test ( );
  r8_sind_test ( );
  r8_swap_test ( );
  r8_swap3_test ( );
  r8_tand_test ( );
  r8_to_i4_test ( );
  r8_to_r8_discrete_test ( );
  r8_uniform_01_test ( );
  r8_uniform_ab_test ( );
  r8_walsh_1d_test ( );
  r8_wrap_test ( );

  r82col_print_part_test ( );

  r82poly2_type_test ( );

  r82row_order_type_test ( );
  r82row_part_quick_a_test ( );
  r82row_print_part_test ( );
  r82row_sort_heap_index_a_test ( );
  r82row_sort_quick_a_test ( );

  r83col_print_part_test ( );

  r83row_print_part_test ( );

  r8block_expand_linear_test ( );
  r8block_new_test ( );
  r8block_print_test ( );

  r8cmat_to_r8mat_new_test ( );

  r8col_find_test ( );
  r8col_insert_test ( );
  r8col_sort_heap_a_test ( );
  r8col_sort_heap_index_a_test ( );
  r8col_sort_quick_a_test ( );
  r8col_sorted_tol_unique_test ( );
  r8col_sorted_unique_count_test ( );
  r8col_sorted_tol_undex_test ( );
  r8col_max_test ( );
  r8col_mean_test ( );
  r8col_min_test ( );
  r8col_permute_test ( );
  r8col_sortr_a_test ( );
  r8col_sum_test ( );
  r8col_swap_test ( );
  r8col_to_r8vec_test ( );
  r8col_tol_undex_test ( );
  r8col_undex_test ( );
  r8col_unique_count_test ( );
  r8col_variance_test ( );

  r8int_to_i4int_test ( );

  r8mat_cholesky_inverse_test ( );
  r8mat_cholesky_solve_test ( );
  r8mat_cholesky_solve_upper_test ( );
  r8mat_det_2d_test ( );
  r8mat_det_3d_test ( );
  r8mat_det_4d_test ( );
  r8mat_det_5d_test ( );
  r8mat_expand_linear_test ( );
  r8mat_expand_linear2_test ( );
  r8mat_fs_new_test ( );
  r8mat_fss_new_test ( );
  r8mat_givens_post_test ( );
  r8mat_givens_pre_test ( );
  r8mat_hess_test ( );
  r8mat_house_axh_test ( );
  r8mat_house_form_test ( );
  r8mat_house_post_test ( );
  r8mat_house_pre_test ( );
  r8mat_indicator_new_test ( );
  r8mat_inverse_2d_test ( );
  r8mat_inverse_3d_test ( );
  r8mat_inverse_4d_test ( );
  r8mat_jac_test ( );
  r8mat_kronecker_test ( );
  r8mat_l_inverse_test ( );
  r8mat_l_print_test ( );
  r8mat_l1_inverse_test ( );
  r8mat_lu_test ( );
  r8mat_max_test ( );
  r8mat_max_index_test ( );
  r8mat_maxcol_minrow_test ( );
  r8mat_maxrow_mincol_test ( );
  r8mat_min_test ( );
  r8mat_min_index_test ( );
  r8mat_mincol_maxrow_test ( );
  r8mat_minrow_maxcol_test ( );
  r8mat_mm_test ( );
  r8mat_mm_new_test ( );
  r8mat_mv_new_test ( );
  r8mat_mv_test ( );
  r8mat_mtv_new_test ( );
  r8mat_mtv_test ( );
  r8mat_nint_test ( );
  r8mat_nonzeros_test ( );
  r8mat_norm_fro_test ( );
  r8mat_norm_l1_test ( );
  r8mat_nullspace_test ( );
  r8mat_nullspace_size_test ( );
  r8mat_orth_uniform_new_test ( );
  r8mat_plot_test ( );
  r8mat_power_method_test ( );
  r8mat_print_test ( );
  r8mat_print_some_test ( );
  r8mat_ref_test ( );
  r8mat_rref_test ( );
  r8mat_solve_test ( );
  r8mat_solve_2d_test ( );
  r8mat_solve_3d_test ( );
  r8mat_solve2_test ( );
  r8mat_sub_new_test ( );
  r8mat_symm_jacobi_test ( );
  r8mat_to_r8cmat_new_test ( );
  r8mat_to_r8plu_test ( );
  r8mat_to_r8rmat_test ( );
  r8mat_trace_test ( );
  r8mat_transpose_new_test ( );
  r8mat_transpose_print_test ( );
  r8mat_u_inverse_test ( );
  r8mat_u1_inverse_test ( );
  r8mat_uniform_ab_new_test ( );

  r8plu_det_test ( );
  r8plu_inverse_test ( );
  r8plu_mul_test ( );
  r8plu_sol_test ( );
  r8plu_to_r8mat_test ( );

  r8poly_degree_test ( );
  r8poly_deriv_test ( );
  r8poly_lagrange_coef_test ( );
  r8poly_lagrange_0_test ( );
  r8poly_lagrange_1_test ( );
  r8poly_lagrange_2_test ( );
  r8poly_lagrange_factor_test ( );
  r8poly_lagrange_val_test ( );
  r8poly_print_test ( );
  r8poly_value_horner_test ( );
  r8poly_values_horner_test ( );

  r8poly2_ex_test ( );
  r8poly2_ex2_test ( );
  r8poly2_val_test ( );
  r8poly2_val2_test ( );

  r8r8vec_index_insert_unique_test ( );

  r8r8r8vec_index_insert_unique_test ( );

  r8rmat_new_test ( );
  r8rmat_to_r8mat_test ( );

  r8row_max_test ( );
  r8row_mean_test ( );
  r8row_min_test ( );
  r8row_sum_test ( );
  r8row_swap_test ( );
  r8row_to_r8vec_test ( );
  r8row_variance_test ( );

  r8slmat_print_test ( );

  r8vec_amax_test ( );
  r8vec_amin_test ( );
  r8vec_bracket_test ( );
  r8vec_bracket2_test ( );
  r8vec_bracket3_test ( );
  r8vec_bracket5_test ( );
  r8vec_chebyspace_new_test ( );
  r8vec_concatenate_new_test ( );
  r8vec_convolution_test ( );
  r8vec_convolution_circ_test ( );
  r8vec_dif_test ( );
  r8vec_direct_product_test ( );
  r8vec_direct_product2_test ( );
  r8vec_even_test ( );
  r8vec_even2_test ( );
  r8vec_expand_linear_test ( );
  r8vec_frac_test ( );
  r8vec_histogram_test ( );
  r8vec_house_column_test ( );
  r8vec_index_delete_all_test ( );
  r8vec_index_delete_dupes_test ( );
  r8vec_index_delete_one_test ( );
  r8vec_index_insert_test ( );
  r8vec_index_insert_unique_test ( );
  r8vec_index_order_test ( );
  r8vec_index_search_test ( );
  r8vec_index_sorted_range_test ( );
  r8vec_indexed_heap_d_test ( );
  r8vec_indexed_heap_d_extract_test ( );
  r8vec_indexed_heap_d_insert_test ( );
  r8vec_indexed_heap_d_max_test ( );
  r8vec_indicator0_new_test ( );
  r8vec_legendre_test ( );
  r8vec_linspace_new_test ( );
  r8vec_max_test ( );
  r8vec_max_index_test ( );
  r8vec_mean_test ( );
  r8vec_median_test ( );
  r8vec_midspace_new_test ( );
  r8vec_min_test ( );
  r8vec_min_index_test ( );
  r8vec_nint_test ( );
  r8vec_norm_l0_test ( );
  r8vec_norm_l1_test ( );
  r8vec_norm_l2_test ( );
  r8vec_norm_li_test ( );
  r8vec_normal_01_test ( );
  r8vec_normalize_l1_test ( );
  r8vec_order_type_test ( );
  r8vec_permute_test ( );
  r8vec_permute_uniform_test ( );
  r8vec_polarize_test ( );
  r8vec_print_test ( );
  r8vec_rotate_test ( );
  r8vec_reverse_test ( );
  r8vec_search_binary_a_test ( );
  r8vec_sort_bubble_a_test ( );
  r8vec_sort_heap_a_test ( );
  r8vec_sort_heap_d_test ( );
  r8vec_sort_heap_index_a_new_test ( );
  r8vec_sort_heap_index_d_new_test ( );
  r8vec_sort_heap_mask_a_test ( );
  r8vec_sort_insert_a_test ( );
  r8vec_sort_insert_index_a_test ( );
  r8vec_sort_quick_a_test ( );
  r8vec_sorted_merge_a_test ( );
  r8vec_sorted_nearest_test ( );
  r8vec_sorted_range_test ( );
  r8vec_sorted_split_test ( );
  r8vec_sorted_undex_test ( );
  r8vec_sorted_unique_test ( );
  r8vec_sorted_unique_count_test ( );
  r8vec_sorted_unique_hist_test ( );
  r8vec_split_test ( );
  r8vec_transpose_print_test ( );
  r8vec_undex_test ( );
  r8vec_uniform_01_new_test ( );
  r8vec_uniform_ab_new_test ( );
  r8vec_variance_test ( );

  r8vec2_sort_a_test ( );
  r8vec2_sort_d_test ( );
  r8vec2_sort_heap_index_a_test ( );
  r8vec2_sorted_unique_test ( );
  r8vec2_sorted_unique_index_test ( );
  r8vec2_sum_max_index_test ( );

  roots_to_r8poly_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "R8LIB_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void i4int_to_r8int_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4INT_TO_R8INT_TEST tests I4INT_TO_R8INT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi = 11;
  int ilo = 1;
  int ir;
  double r;
  double r2;
  double rhi = 200.0;
  double rhi2;
  double rlo = 100.0;
  double rlo2;
  int seed;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "I4INT_TO_R8INT_TEST\n";
  cout << "  For data in an interval,\n";
  cout << "  I4INT_TO_R8INT converts an integer to a real;\n";
  cout << "\n";
  cout << "  Integer interval: [" << ilo << ", " << ihi << "]\n";
  cout << "  Real interval:    [" << rlo << ", " << rhi << "]\n";
  cout << "\n";
  cout << "         R         I(R)        R(I(R))\n";
  cout << "\n";

  seed = 123456789;

  rlo2 = rlo - 15.0;
  rhi2 = rhi + 15.0;

  for ( test = 1; test <= test_num; test++ )
  {
    r = r8_uniform_ab ( rlo2, rhi2, seed );
    ir = r8int_to_i4int ( rlo, rhi, r, ilo, ihi );
    r2 = i4int_to_r8int ( ilo, ihi, ir, rlo, rhi );
    cout << "  " << setw(12) << r
         << "  " << setw(6)  << ir
         << "  " << setw(12) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void perm0_check_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PERM0_CHECK_TEST tests PERM0_CHECK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  bool check;
  int n = 5;
  int p1[5] = { 5, 2, 3, 4, 1 };
  int p2[5] = { 4, 1, 3, 0, 2 };
  int p3[5] = { 0, 2, 1, 3, 2 };

  cout << "\n";
  cout << "PERM0_CHECK_TEST\n";
  cout << "  PERM0_CHECK checks a permutation of 0, ..., N-1.\n";
  cout << "\n";

  i4vec_transpose_print ( n, p1, "  Permutation 1:" );
  check = perm0_check( n, p1 );

  i4vec_transpose_print ( n, p2, "  Permutation 2:" );
  check = perm0_check( n, p2 );

  i4vec_transpose_print ( n, p3, "  Permutation 3:" );
  check = perm0_check( n, p3 );

  return;
}
//****************************************************************************80

void perm0_uniform_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PERM0_UNIFORM_TEST tests PERM0_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 10;
  int *p;
  int seed;
  int test;

  cout << "\n";
  cout << "PERM0_UNIFORM_TEST\n";
  cout << "  PERM0_UNIFORM randomly selects a permutation of 0,...,N-1.\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    p = perm0_uniform_new ( n, seed );
    cout << "  ";
    for ( i = 0; i < n; i++ )
    {
      cout << setw(4) << p[i];
    }
    cout << "\n";
    delete [] p;
  }
  return;
}
//****************************************************************************80

void perm1_check_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PERM1_CHECK_TEST tests PERM1_CHECK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  bool check;
  int n = 5;
  int p1[5] = { 5, 2, 3, 4, 1 };
  int p2[5] = { 4, 1, 3, 0, 2 };
  int p3[5] = { 0, 2, 1, 3, 2 };

  cout << "\n";
  cout << "PERM1_CHECK_TEST\n";
  cout << "  PERM1_CHECK checks a permutation of 1, ..., N.\n";
  cout << "\n";

  i4vec_transpose_print ( n, p1, "  Permutation 1:" );
  check = perm1_check( n, p1 );

  i4vec_transpose_print ( n, p2, "  Permutation 2:" );
  check = perm1_check( n, p2 );

  i4vec_transpose_print ( n, p3, "  Permutation 3:" );
  check = perm1_check( n, p3 );

  return;
}
//****************************************************************************80

void perm1_uniform_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PERM1_UNIFORM_TEST tests PERM1_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 10;
  int *p;
  int seed;
  int test;

  cout << "\n";
  cout << "PERM1_UNIFORM_TEST\n";
  cout << "  PERM1_UNIFORM randomly selects a permutation of 1,...,N.\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    p = perm1_uniform_new ( n, seed );
    cout << "  ";
    for ( i = 0; i < n; i++ )
    {
      cout << setw(4) << p[i];
    }
    cout << "\n";
    delete [] p;
  }
  return;
}
//****************************************************************************80

void r8_abs_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS_TEST tests R8_ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double r8;
  double r8_absolute;
  double r8_hi = 5.0;
  double r8_lo = -3.0;
  int seed;
  int test;
  int test_num = 10;

  seed = 123456789;

  cout << "\n";
  cout << "R8_ABS_TEST\n";
  cout << "  R8_ABS returns the absolute value of an R8.\n";
  cout << "\n";
  cout << "      X         R8_ABS(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    r8 = r8_uniform_ab ( r8_lo, r8_hi, seed );
    r8_absolute = r8_abs ( r8 );
    cout << "  " << setw(10) << r8
         << "  " << setw(10) << r8_absolute << "\n";
  }

  return;
}
//****************************************************************************80

void r8_acos_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ACOS_TEST tests R8_ACOS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  int test;

  cout << "\n";
  cout << "R8_ACOS_TEST\n";
  cout << "  R8_ACOS computes the arc-cosine of an angle.\n"; 
  cout << "\n";
  cout << "       C            R8_ACOS(C)        ACOS(C)\n";
  cout << "\n";

  for ( test = -1; test <= 13; test++ )
  {
    c = ( double ) ( test - 6 ) / ( double ) ( 6 );

    cout << setw(14) << c << "  "
         << setw(14) << r8_acos ( c );

    if ( -1.0 <= c && c <= 1.0 )
    {
      cout << "  " << setw(14) << acos ( c );
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void r8_acosh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ACOSH_TEST tests R8_ACOSH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int test;
  double x;
  double x2;

  cout << "\n";
  cout << "R8_ACOSH_TEST\n";
  cout << "  R8_ACOSH computes the arc-hyperbolic-cosine of an angle.\n";
  cout << "\n";
  cout << "       X            A=R8_ACOSH(X)    COSH(A)\n";
  cout << "\n";

  for ( test = 0; test <= 8; test++ )
  {
    x = 1.0 + ( double ) ( test ) / 2.0;
    a = r8_acosh ( x );
    x2 = cosh ( a );
    cout << setw(14) << x << "  "
         << setw(14) << a << "  "
         << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_asinh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ASINH_TEST tests R8_ASINH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double x;
  double x2;

  cout << "\n";
  cout << "R8_ASINH_TEST\n";
  cout << "  R8_ASINH computes the inverse hyperbolic sine\n";
  cout << "  of a given value.\n";
  cout << "\n";
  cout << "         X   R8_ASINH(X)     SINH(R8_ASINH(X))\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    a = r8_asinh ( x );
    x2 = sinh ( a );

    cout                   << "  "
         << setw(10) << x  << "  "
         << setw(10) << a  << "  "
         << setw(10) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_atan_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ATAN_TEST tests R8_ATAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 8

  int test;
  double x;
  double xtest[TEST_NUM] = {
     1.0,  1.0,  0.0, -1.0,
    -1.0, -1.0,  0.0,  1.0 };
  double y;
  double ytest[TEST_NUM] = {
     0.0,  1.0,  1.0,  1.0,
     0.0, -1.0, -1.0, -1.0 };

  cout << "\n";
  cout << "R8_ATAN_TEST\n";
  cout << "  R8_ATAN computes the arc-tangent given Y and X;\n";
  cout << "  ATAN2 is the system version of this routine.\n";
  cout << "\n";
  cout << "       X             Y          ATAN2(Y,X)    R8_ATAN(Y,X)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = xtest[test];
    y = ytest[test];
    cout << "  " << setw(14) << x
         << "  " << setw(14) << y
         << "  " << setw(14) << atan2 ( y, x )
         << "  " << setw(14) << r8_atan ( y, x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8_atanh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ATANH_TEST tests R8_ATANH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double x;
  double x2;

  cout << "\n";
  cout << "R8_ATANH_TEST\n";
  cout << "  R8_ATANH computes the inverse hyperbolic tangent\n";
  cout << "  of a given value.\n";
  cout << "\n";
  cout << "         X     R8_ATANH(X)     TANH(R8_ATANH(X))\n";
  cout << "\n";

  for ( i = -2; i <= 9; i++ )
  {
    x = ( ( double ) i ) / 10.0;
    a = r8_atanh ( x );
    x2 = tanh ( a );

    cout                   << "  "
         << setw(10) << x  << "  "
         << setw(10) << a  << "  "
         << setw(10) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_big_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BIG_TEST tests R8_BIG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "R8_BIG_TEST\n";
  cout << "  R8_BIG returns a 'big' R8 value;\n";
  cout << "\n";
  cout << "  R8_BIG =   " << r8_big ( ) << "\n";

  return;
}
//****************************************************************************80

void r8_cas_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CAS_TEST tests R8_CAS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 12

  int test;
  double x;

  cout << "\n";
  cout << "R8_CAS_TEST\n";
  cout << "  R8_CAS evaluates the casine of a number.\n";
  cout << "\n";
  cout << "          X           R8_CAS ( X )\n";
  cout << "\n";

  for ( test = 0; test <= TEST_NUM; test++ )
  {
    x = r8_pi ( ) * ( double ) ( test ) / ( double ) ( TEST_NUM );
    cout << "  " << setw(14) << x
         << "  " << setw(14) << r8_cas ( x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8_ceiling_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CEILING_TEST tests R8_CEILING.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double rval;
  double rval_rounded;

  cout << "\n";
  cout << "R8_CEILING_TEST\n";
  cout << "  R8_CEILING rounds a value up.\n";
  cout << "\n";
  cout << "           X       R8_CEILING(X)\n";
  cout << "\n";

  for ( i = -6; i <= 6; i++ )
  {
    rval = ( double ) ( i ) / 5.0;
    rval_rounded = r8_ceiling ( rval );
    cout << "  " << setw(14) << rval
         << "  " << setw(14) << rval_rounded << "\n";
  }

  return;
}
//****************************************************************************80

void r8_choose_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE_TEST tests R8_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double cnk;
  int k;
  int n;

  cout << "\n";
  cout << "R8_CHOOSE_TEST\n";
  cout << "  R8_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "         N         K       CNK\n";
 
  for ( n = 0; n <= 5; n++ )
  {
    cout << "\n";
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );
      cout << setw(10) << n << "  "
           << setw(8) << k << "  "
           << setw(14) << cnk << "\n";
    }
  }
 
  return;
}
//****************************************************************************80

void r8_cosd_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_COSD_TEST tests R8_COSD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double angle;
  int i;

  cout << "\n";
  cout << "R8_COSD_TEST\n";
  cout << "  R8_COSD computes the cosine of an angle\n";
  cout << "  given in degrees.\n";
  cout << "\n";
  cout << "  ANGLE    R8_COSD(ANGLE)\n";
  cout << "\n";
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    cout << "  " << setw(8) << angle
         << "  " << setw(14) <<  r8_cosd ( angle ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void r8_cotd_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_COTD_TEST tests R8_COTD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double angle;
  int i;

  cout << "\n";
  cout << "R8_COTD_TEST\n";
  cout << "  R8_COTD computes the cotangent of an angle\n";
  cout << "  given in degrees.\n";
  cout << "\n";
  cout << "  ANGLE    R8_COTD(ANGLE)\n";
  cout << "\n";
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    if ( i % 180 == 0 )
    {
      cout << "  " << setw(8) << angle
           << "    Undefined\n";
    }
    else
    {
      cout << "  " << setw(8) << angle
           << "  " << setw(14) <<  r8_cotd ( angle ) << "\n";
    }
  }
 
  return;
}
//****************************************************************************80

void r8_cscd_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CSCD_TEST tests R8_CSCD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double angle;
  int i;

  cout << "\n";
  cout << "R8_CSCD_TEST\n";
  cout << "  R8_CSCD computes the cosecant of an angle\n";
  cout << "  given in degrees.\n";
  cout << "\n";
  cout << "  ANGLE    R8_CSCD(ANGLE)\n";
  cout << "\n";
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    if ( i % 180 == 0 )
    {
      cout << "  " << setw(8) << angle
           << "    Undefined\n";
    }
    else
    {
      cout << "  " << setw(8) << angle
           << "  " << setw(14) <<  r8_cscd ( angle ) << "\n";
    }
  }
 
  return;
}
//****************************************************************************80

void r8_cube_root_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CUBE_ROOT_TEST tests R8_CUBE_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int seed;
  double x1;
  double y;
  double x2;

  cout << "\n";
  cout << "R8_CUBE_ROOT_TEST\n";
  cout << "  R8_CUBE_ROOT computes the cube root of an R8.\n";
  cout << "\n";
  cout << "       X               Y               Y^3\n";
  cout << "\n";

  a = -10.0;
  b = +10.0;
  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x1 = r8_uniform_ab ( a, b, seed );
    y = r8_cube_root ( x1 );
    x2 = pow ( y, 3 );
    cout << setw(14) << x1 << "  "
         << setw(14) << y << "  "
         << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_diff_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_DIFF_TEST tests R8_DIFF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 15

  int ndig = 3;
  int test;
  double x = 1.0;
  double y;
  double y_test[TEST_NUM] = {
    0.0625, 0.125, 0.25, 0.50,  0.874,
    0.876,  0.90,  0.95, 0.99,  1.0,
    1.01,   1.05,  1.10, 3.0,  10.0 };

  cout << "\n";
  cout << "R8_DIFF_TEST\n";
  cout << "  R8_DIFF computes a difference X-Y to a given\n";
  cout << "  number of binary places.\n";
  cout << "\n";
  cout << "  For this test, we use " << ndig << " binary places.\n";
  cout << "\n";
  cout << "       X       Y       X-Y     R8_DIFF(X,Y)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    y = y_test[test];
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y
         << "  " << setw(10) << x-y
         << "  " << setw(10) << r8_diff ( x, y, ndig ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8_digit_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_DIGIT_TEST tests R8_DIGIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define MAXDIG 20

  int idigit;
  double x;

  x = r8_pi ( );

  cout << "\n";
  cout << "R8_DIGIT_TEST\n";
  cout << "  R8_DIGIT extracts decimal digits.\n";
  cout << "\n";
  cout << "  Here, we get digits of " << x << "\n";
  cout << "\n";

  cout << "  ";
  for ( idigit = -2; idigit <= MAXDIG; idigit++ )
  {
    cout << setw(3) << idigit;
  }
  cout << "\n";

  cout << "  ";
  for ( idigit = -2; idigit <= MAXDIG; idigit++ )
  {
    cout << setw(3) << r8_digit ( x, idigit );
  }
  cout << "\n";

  return;
# undef MAXDIG
}
//****************************************************************************80

void r8_e_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_E_TEST tests R8_E.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double value1;
  double value2;

  cout << "\n";
  cout << "R8_E_TEST\n";
  cout << "  R8_E returns the value of E.\n";
  cout << "  Compare E to (1+1/n)^n\n";
  value1 = r8_e ( );
  cout << "  R8_E =      " << value1 << "\n";
  cout << "\n";
  cout << "         N     Estimate      Error\n";
  cout << "\n";

  n = 1;
  for ( i = 0; i <= 20; i++ )
  {
    value2 = pow ( ( double ) ( n + 1 ) / ( double ) ( n ), n );
    cout << "  " << setw(8) << n
         << "  " << setw(14) << value2
         << "  " << setw(14) << fabs ( value1 - value2 ) << "\n";
    n = n * 2;
  }

  return;
}
//****************************************************************************80

void r8_epsilon_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON_TEST tests R8_EPSILON.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double r;
  double s;
  streamsize ss;
  double t;

  cout << "\n";
  cout << "R8_EPSILON_TEST\n";
  cout << "  R8_EPSILON produces the R8 roundoff unit.\n";
  cout << "\n";
//
//  Save the current precision.
//
  ss = cout.precision ( );

  r = r8_epsilon ( );
  cout << "  R = R8_EPSILON()         = " 
       << setprecision(16) << setw(24) << r << "\n";

  s = 1.0 + r;
  t = s - 1.0;
  cout << "  ( 1 + R ) - 1           =     " 
       << setprecision(16) << setw(24) << t << "\n";

  s = 1.0 + ( r / 2.0 );
  t = s - 1.0;
  cout << "  ( 1 + (R/2) ) - 1       = " << setw(10) << t << "\n";
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_epsilon_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON_COMPUTE_TEST tests R8_EPSILON_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double r;
  double s;
  streamsize ss;
  double t;

  cout << "\n";
  cout << "R8_EPSILON_COMPUTE_TEST\n";
  cout << "  R8_EPSILON_COMPUTE computes the R8 roundoff unit.\n";
  cout << "\n";
//
//  Save the current precision.
//
  ss = cout.precision ( );

  r = r8_epsilon_compute ( );
  cout << "  R = R8_EPSILON_COMPUTE() =  " 
       << setprecision(16) << setw(24) << r << "\n";

  s = 1.0 + r;
  t = s - 1.0;
  cout << "  ( 1 + R ) - 1           =     " 
       << setprecision(16) << setw(24) << t << "\n";

  s = 1.0 + ( r / 2.0 );
  t = s - 1.0;
  cout << "  ( 1 + (R/2) ) - 1       = " << setw(10) << t << "\n";
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_factorial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_TEST tests R8_FACTORIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double f1;
  double f2;
  int n;
  int n_data;
  streamsize ss;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  cout << "\n";
  cout << "R8_FACTORIAL_TEST\n";
  cout << "  R8_FACTORIAL evaluates the factorial function.\n";
  cout << "\n";
  cout << "    N                Exact                  Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( n_data, n, f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_factorial ( n );

    cout << "  "
         << setw(4) << n << "  "
         << setprecision(16) << setw(24) << f1 << "  "
         << setprecision(16) << setw(24) << f2 << "\n";
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_factorial2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL2_TEST tests R8_FACTORIAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double f1;
  double f2;
  int n;
  int n_data;
  streamsize ss;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  cout << "\n";
  cout << "R8_FACTORIAL2_TEST\n";
  cout << "  R8_FACTORIAL2 evaluates the double factorial function.\n";
  cout << "\n";
  cout << "    N                Exact                  Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial2_values ( n_data, n, f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_factorial2 ( n );

    cout << "  "
         << setw(4) << n << "  "
         << setprecision(16) << setw(24) << f1 << "  "
         << setprecision(16) << setw(24) << f2 << "\n";
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_fall_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FALL_TEST tests R8_FALL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double f1;
  double f2;
  int n;
  int n_data;
  streamsize ss;
  double x;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  cout << "\n";
  cout << "R8_FALL_TEST\n";
  cout << "  R8_FALL evaluates the falling factorial Fall(X,N).\n";
  cout << "\n";
  cout << "    X          N                Exact                  Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_fall_values ( n_data, x, n, f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_fall ( x, n );

    cout << "  "
         << setw(8) << x << "  "
         << setw(4) << n << "  "
         << setprecision(16) << setw(24) << f1 << "  "
         << setprecision(16) << setw(24) << f2 << "\n";
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_fractional_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FRACTIONAL_TEST tests R8_FRACTIONAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  double fractional;
  double r8;
  double r8_hi = 5.0;
  double r8_lo = -3.0;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "R8_FRACTIONAL_TEST\n";
  cout << "  R8_FRACTIONAL returns the fractional part of an R8.\n";
  cout << "\n";
  cout << "          X           R8_FRACTIONAL ( X )\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    r8 = r8_uniform_ab ( r8_lo, r8_hi, seed );
    fractional = r8_fractional ( r8 );
    cout << "  " << setw(14) << r8
         << "  " << setw(14) << fractional << "\n";
  }

  return;
}
//****************************************************************************80

void r8_gamma_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_TEST tests R8_GAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int n_data;
  streamsize ss;
  double x;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  cout << "\n";
  cout << "R8_GAMMA_TEST:\n";
  cout << "   R8_GAMMA evaluates the Gamma function.\n";
  cout << "\n";
  cout << "      X            GAMMA(X)     R8_GAMMA(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = r8_gamma ( x );

    cout << "  " << setw(12)                     << x  
         << "  " << setw(24) << setprecision(16) << fx1 
         << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_gamma_log_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG_TEST tests R8_GAMMA_LOG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int n_data;
  streamsize ss;
  double x;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  cout << "\n";
  cout << "R8_GAMMA_LOG_TEST:\n";
  cout << "   R8_GAMMA_LOG evaluates the logarithm of the Gamma function.\n";
  cout << "\n";
  cout << "      X            GAMMA_LOG(X)     R8_GAMMA_LOG(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = r8_gamma_log ( x );
    cout << "  " << setw(12)                     << x
         << "  " << setw(24) << setprecision(16) << fx1
         << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }

//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_huge_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE_TEST tests R8_HUGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "R8_HUGE_TEST\n";
  cout << "  R8_HUGE returns a 'huge' R8 value;\n";
  cout << "\n";
  cout << "  R8_HUGE =   " << r8_huge ( ) << "\n";

  return;
}
//****************************************************************************80

void r8_log_2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LOG_2_TEST tests R8_LOG_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 18

  int test;
  double x;
  double x_test[TEST_NUM] = {
    0.0,  1.0,  2.0,   3.0,  9.0,
   10.0, 11.0, 99.0, 101.0, -1.0,
   -2.0, -3.0, -9.0,   0.5,  0.33,
    0.25, 0.20, 0.01 };

  cout << "\n";
  cout << "R8_LOG_2_TEST\n";
  cout << "  R8_LOG_2: computes the logarithm base 2.\n";
  cout << "\n";
  cout << "        X       R8_LOG_2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    cout << "  " << setw(12) << x
         << "  " << setw(12) << r8_log_2 ( x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8_log_b_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LOG_B_TEST tests R8_LOG_B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 10

  double b;
  double b_test[TEST_NUM] = {
    2.0, 3.0, 4.0, 5.0, 6.0,
    7.0, 8.0, 16.0, 32.0, 256.0 };
  int test;
  double x;

  x = 16.0;

  cout << "\n";
  cout << "R8_LOG_B_TEST\n";
  cout << "  R8_LOG_B computes the logarithm base B.\n";
  cout << "\n";
  cout << "       X       B      R8_LOG_B\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    b = b_test[test];

    cout << "  " << setw(12) << x
         << "  " << setw(12) << b
         << "  " << setw(12) << r8_log_b ( x, b ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8_mant_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MANT_TEST tests R8_MANT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int l;
  double r;
  int s;
  double x;

  x = -314.159;

  cout << "\n";
  cout << "R8_MANT_TEST\n";
  cout << "  R8_MANT decomposes a value.\n";
  cout << "\n";
  cout << "  Number to be decomposed: X = " << x << "\n";

  r8_mant ( x, s, r, l );

  cout << "\n";
  cout << "  X = " << s << " * " << r << " * 2 ^ " << l << "\n";

  return;
}
//****************************************************************************80

void r8_max_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX_TEST tests R8_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  int i;
  double r8_hi;
  double r8_lo;
  int seed;

  cout << "\n";
  cout << "R8_MAX_TEST\n";
  cout << "  R8_MAX returns the maximum of two R8's.\n";
  cout << "\n";
  cout << "       A       B      C=R8_MAX(A,B)\n";
  cout << "\n";

  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;

  for ( i = 0; i < 10; i++ )
  {
    a = r8_uniform_ab ( r8_lo, r8_hi, seed );
    b = r8_uniform_ab ( r8_lo, r8_hi, seed );
    c = r8_max ( a, b );
    cout << "  " << setw(8) << a
         << "  " << setw(8) << b
         << "  " << setw(8) << c << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void r8_min_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN_TEST tests R8_MIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  int i;
  double r8_hi;
  double r8_lo;
  int seed;

  cout << "\n";
  cout << "R8_MIN_TEST\n";
  cout << "  R8_MIN returns the minimum of two R8's.\n";
  cout << "\n";
  cout << "       A       B      C=R8_MIN(A,B)\n";
  cout << "\n";

  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;

  for ( i = 0; i < 10; i++ )
  {
    a = r8_uniform_ab ( r8_lo, r8_hi, seed );
    b = r8_uniform_ab ( r8_lo, r8_hi, seed );
    c = r8_min ( a, b );
    cout << "  " << setw(8) << a
         << "  " << setw(8) << b
         << "  " << setw(8) << c << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void r8_mod_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOD_TEST tests R8_MOD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  double x;
  double x_hi = 10.0;
  double x_lo = -10.0;
  double y;
  double z1;
  double z2;

  cout << "\n";
  cout << "R8_MOD_TEST\n";
  cout << "  R8_MOD returns the remainder after division.\n";
  cout << "  R8_MOD ( X, Y ) has the same sign as X.\n";
  cout << "\n";
  cout << "      X         Y    FMOD(X,Y)    R8_MOD(X,Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r8_uniform_ab ( x_lo, x_hi, seed );
    y = r8_uniform_ab ( x_lo, x_hi, seed );

    z1 =   fmod ( x, y );
    z2 = r8_mod ( x, y );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << y
         << "  " << setw(12) << z1
         << "  " << setw(12) << z2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_modp_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MODP_TEST tests R8_MODP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  double x;
  double x_hi = 10.0;
  double x_lo = -10.0;
  double y;
  double z1;
  double z2;

  cout << "\n";
  cout << "R8_MODP_TEST\n";
  cout << "  R8_MODP returns the remainder after division.\n";
  cout << "  R8_MODP ( X, Y ) is positive.\n";
  cout << "\n";
  cout << "      X       Y     FMOD(X,Y)  R8_MODP(X,Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r8_uniform_ab ( x_lo, x_hi, seed );
    y = r8_uniform_ab ( x_lo, x_hi, seed );

    z1 =   fmod  ( x, y );
    z2 = r8_modp ( x, y );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << y
         << "  " << setw(12) << z1
         << "  " << setw(12) << z2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_mop_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP_TEST tests R8_MOP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i4;
  int i4_max;
  int i4_min;
  double r8;
  int seed = 123456789;
  int test;

  cout << "\n";
  cout << "R8_MOP_TEST\n";
  cout << "  R8_MOP evaluates (-1.0)^I4 as an R8.\n";
  cout << "\n";
  cout << "    I4  R8_MOP(I4)\n";
  cout << "\n";

  i4_min = -100;
  i4_max = +100;

  for ( test = 1; test <= 10; test++ )
  {
    i4 = i4_uniform_ab ( i4_min, i4_max, seed );
    r8 = r8_mop ( i4 );
    cout << "  "
         << setw(4) << i4 << "  "
         << setw(4) <<r8 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_nint_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT_TEST tests R8_NINT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double c;
  int seed = 123456789;
  int test;
  int test_num = 10;
  double x;

  cout << "\n";
  cout << "R8_NINT_TEST\n";
  cout << "  R8_NINT produces the nearest integer.\n";
  cout << "\n";
  cout << "      X      R8_NINT(X)\n";
  cout << "\n";

  b = -10.0;
  c = +10.0;

  for ( test = 1; test <= test_num; test++ )
  {
    x = r8_uniform_ab ( b, c, seed );
    cout << setw(10) << x << "  "
         << setw(6)  << r8_nint ( x ) << "\n";
  }

  return;
}
//****************************************************************************80

void r8_normal_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01_TEST tests R8_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 20

  int seed = 123456789;
  int test;
  double x;

  cout << "\n";
  cout << "R8_NORMAL_01_TEST\n";
  cout << "  R8_NORMAL_01 generates normally distributed random values.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = r8_normal_01 ( seed );
    cout                  << "  "
         << setw(10) << x << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8_pi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PI_TEST tests R8_PI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double four;
  double one;
  double v1;
  double v2;

  four = ( double ) ( 4 );
  one = ( double ) ( 1 );

  cout << "\n";
  cout << "R8_PI_TEST\n";
  cout << "  R8_PI returns the value of PI.\n";
  cout << "\n";
  v1 = r8_pi ( );
  cout << "  R8_PI =     " << v1 << "\n";
  v2 = four * atan ( one );
  cout << "  4*atan(1) = " << v2 << "\n";

  return;
}
//****************************************************************************80

void r8_power_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_POWER_TEST tests R8_POWER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int p;
  double r;
  double value;

  cout << "\n";
  cout << "R8_POWER_TEST\n";
  cout << "  R8_POWER computes R^P\n";
  cout << "\n";
  cout << "      R          P       R^P\n";
  cout << "\n";

  for ( p = -5; p <= 5; p++ )
  {
    r = 2.0;
    value = r8_power ( r, p );
    cout << "  " << setw(12) << r
         << "  " << setw(6)  << p
         << "  " << setw(12) << value << "\n";
  }

  return;
}
//****************************************************************************80

void r8_power_fast_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_POWER_FAST_TEST tests R8_POWER_FAST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int mults;
  int p;
  double r;
  double rp;

  cout << "\n";
  cout << "R8_POWER_FAST_TEST\n";
  cout << "  R8_POWER_FAST computes R^P, economizing on\n";
  cout << "  multiplications.\n";
  cout << "\n";
  cout << "      R          P       R^P         Mults\n";
  cout << "\n";

  for ( i = -10; i <= 40; i++ )
  {
    r = 2.0;
    p = i;
    rp = r8_power_fast ( r, p, mults );
    cout << "  " << setw(12) << r
         << "  " << setw(6)  << p
         << "  " << setw(12) << rp
         << "  " << setw(6)  << mults << "\n";
  }

  return;
}
//****************************************************************************80

void r8_rise_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_RISE_TEST tests R8_RISE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double f1;
  double f2;
  int n;
  int n_data;
  streamsize ss;
  double x;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  cout << "\n";
  cout << "R8_RISE_TEST\n";
  cout << "  R8_RISE evaluates the rising factorial Fall(X,N).\n";
  cout << "\n";
  cout << "    X          N                Exact                  Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_rise_values ( n_data, x, n, f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_rise ( x, n );

    cout << "  "
         << setw(8) << x << "  "
         << setw(4) << n << "  "
         << setprecision(16) << setw(24) << f1 << "  "
         << setprecision(16) << setw(24) << f2 << "\n";
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_round2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ROUND2_TEST tests R8_ROUND2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nplace;
  streamsize ss;
  double x;
  double xround;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  x = r8_pi ( );

  cout << "\n";
  cout << "R8_ROUND2_TEST\n";
  cout << "  R8_ROUND2 rounds a number to a\n";
  cout << "  specified number of base 2 digits.\n";
  cout << "\n";
  cout << "  Test effect on PI:\n";
  cout << "  X = " << x << "\n";
  cout << "\n";
  cout << "  NPLACE  XROUND\n";
  cout << "\n";

  for ( i = 0; i <= 20; i++ )
  {
    nplace = i;
    xround = r8_round2 ( nplace, x );
    cout << "  " << setw(8) << i
         << "  " << setprecision(16) << setw(24) << xround << "\n";
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_roundb_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ROUNDB_TEST tests R8_ROUNDB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int base;
  int i;
  int nplace;
  streamsize ss;
  double x;
  double xround;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  base = 3;
  x = r8_pi ( );

  cout << "\n";
  cout << "R8_ROUNDB_TEST\n";
  cout << "  R8_ROUNDB rounds a number to a \n";
  cout << "  specified number of base IBASE digits.\n";
  cout << "\n";
  cout << "  Here, we will use IBASE = " << base << "\n";
  cout << "\n";
  cout << "  Test effect on PI:\n";
  cout << "  X = " << setprecision(16) << setw(24) << x << "\n";
  cout << "\n";
  cout << "  NPLACE  XROUND\n";
  cout << "\n";

  for ( i = 0; i <= 20; i++ )
  {
    nplace = i;
    xround = r8_roundb ( base, nplace, x );
    cout << "  " << setw(8) << i
         << "  " << setprecision(16) << setw(24) << xround << "\n";
  }

  cout << "\n";
  cout << "  Try with a negative base:\n";
  x = 121.0;
  base = -3;
  nplace = 3;
  cout << "\n";
  cout << "  Input quantity is X = " << x << "\n";
  cout << "  to be rounded in base " << base << "\n";

  for ( nplace = 1; nplace <= 5; nplace++ )
  {
    xround = r8_roundb ( base, nplace, x );

    cout << "\n";
    cout << "  Output value to " << nplace << " places is " << xround << "\n";
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_roundx_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ROUNDX_TEST tests R8_ROUNDX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nplace;
  int seed;
  streamsize ss;
  double x;
  double xround;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  seed = 123456789;
  x = r8_pi ( );

  cout << "\n";
  cout << "R8_ROUNDX_TEST\n";
  cout << "  R8_ROUNDX rounds a number to a \n";
  cout << "  specified number of decimal digits.\n";
  cout << "\n";
  cout << "  Test effect on PI:\n";
  cout << "  X = " << setprecision(16) << x << "\n";
  cout << "\n";
  cout << "  NPLACE  XROUND\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    nplace = i;
    xround = r8_roundx ( nplace, x );
    cout << "  " << setw(6) << i
         << "  " << setw(20) << xround << "\n";
  }

  cout << "\n";
  cout << "  Test effect on random values:\n";
  cout << "\n";
  cout << "  NPLACE  X     XROUND\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    x = r8_uniform_01 ( seed );

    cout << "\n";

    for ( nplace = 0; nplace <= 10; nplace = nplace + 2 )
    {
      xround = r8_roundx ( nplace, x );

      cout << "  " << setw(6)  << nplace
           << "  " << setw(16) << x
           << "  " << setw(20) << xround << "\n";
    }
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_secd_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SECD_TEST tests R8_SECD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double angle;
  int i;

  cout << "\n";
  cout << "R8_SECD_TEST\n";
  cout << "  R8_SECD computes the secant of an angle\n";
  cout << "  given in degrees.\n";
  cout << "\n";
  cout << "  ANGLE    R8_SECD(ANGLE)\n";
  cout << "\n";
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    if ( ( i + 90 ) % 180 == 0 )
    {
      cout << "  " << setw(8) << angle
           << "    Undefined\n";
    }
    else
    {
      cout << "  " << setw(8) << angle
           << "  " << setw(14) <<  r8_secd ( angle ) << "\n";
    }
  }
 
  return;
}
//****************************************************************************80

void r8_sign_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN_TEST tests R8_SIGN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  double r8;
  double r8_test[5] = { -1.25, -0.25, 0.0, +0.5, +9.0 };
  double s;
  int test;
  const int test_num = 5;

  cout << "\n";
  cout << "R8_SIGN_TEST\n";
  cout << "  R8_SIGN returns the sign of an R8.\n";
  cout << "\n";
  cout << "      R8      R8_SIGN(R8)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    r8 = r8_test[test];
    s = r8_sign ( r8 );
    cout << setw(10) << r8 << "  "
         << setw(10) << s << "\n";
  }

  return;
}
//****************************************************************************80

void r8_sign3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN3_TEST tests R8_SIGN3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  double r8;
  double r8_test[5] = { -1.25, -0.25, 0.0, +0.5, +9.0 };
  double s;
  int test;
  const int test_num = 5;

  cout << "\n";
  cout << "R8_SIGN3_TEST\n";
  cout << "  R8_SIGN3 returns the three-way sign of an R8.\n";
  cout << "\n";
  cout << "      R8      R8_SIGN3(R8)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    r8 = r8_test[test];
    s = r8_sign3 ( r8 );
    cout << setw(10) << r8 << "  "
         << setw(10) << s << "\n";
  }

  return;
}
//****************************************************************************80

void r8_sind_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIND_TEST tests R8_SIND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double angle;
  int i;

  cout << "\n";
  cout << "R8_SIND_TEST\n";
  cout << "  R8_SIND computes the sine of an angle\n";
  cout << "  given in degrees.\n";
  cout << "\n";
  cout << "  ANGLE    R8_SIND(ANGLE)\n";
  cout << "\n";
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    cout << "  " << setw(8) << angle
         << "  " << setw(14) <<  r8_sind ( angle ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void r8_swap_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP_TEST tests R8_SWAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double x;
  double y;

  cout << "\n";
  cout << "R8_SWAP_TEST\n";
  cout << "  R8_SWAP swaps two reals.\n";

  x = 1.0;
  y = 3.14159;

  cout << "\n";
  cout << "  Before swapping: \n";
  cout << "\n";
  cout << "    X = " << x << "\n";
  cout << "    Y = " << y << "\n";

  r8_swap ( x, y );

  cout << "\n";
  cout << "  After swapping: \n";
  cout << "\n";
  cout << "    X = " << x << "\n";
  cout << "    Y = " << y << "\n";

  return;
}
//****************************************************************************80

void r8_swap3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP3_TEST tests R8_SWAP3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double x;
  double y;
  double z;

  cout << "\n";
  cout << "R8_SWAP3_TEST\n";
  cout << "  R8_SWAP3 swaps three reals.\n";

  x = 1.0;
  y = 3.14159;
  z = 1952.0;

  cout << "\n";
  cout << "              X       Y       Z\n";
  cout << "\n";
  cout << "  Start: " << x << "  " << y << "  " << z << "\n";

  for ( i = 1; i <= 3; i++ )
  {
    r8_swap3 ( x, y, z );
    cout << "  Swap " << i << "  " << x << "  " << y << "  " << z << "\n";
  }

  return;
}
//****************************************************************************80

void r8_tand_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TAND_TEST tests R8_TAND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double angle;
  int i;

  cout << "\n";
  cout << "R8_TAND_TEST\n";
  cout << "  R8_TAND computes the tangent of an angle\n";
  cout << "  given in degrees.\n";
  cout << "\n";
  cout << "  ANGLE    R8_TAND(ANGLE)\n";
  cout << "\n";
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    if ( ( i + 90 ) % 180 == 0 )
    {
      cout << "  " << setw(8) << angle
           << "    Undefined\n";
    }
    else
    {
      cout << "  " << setw(8) << angle
           << "  " << setw(14) <<  r8_tand ( angle ) << "\n";
    }
  }
 
  return;
}
//****************************************************************************80

void r8_to_r8_discrete_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_R8_DISCRETE_TEST tests R8_TO_R8_DISCRETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int ndx = 19;
  double r;
  double rd;
  double rhi = 10.0;
  double rhi2;
  double rlo = 1.0;
  double rlo2;
  int seed;
  int test;
  int test_num = 15;

  cout << "\n";
  cout << "R8_TO_R8_DISCRETE_TEST\n";
  cout << "  R8_TO_R8_DISCRETE maps numbers to a discrete set\n";
  cout << "  of equally spaced numbers in an interval.\n";
  cout << "\n";
  cout << "  Number of discrete values = " << ndx << "\n";
  cout << "  Real interval: " << rlo << "  " << rhi << "\n";
  cout << "\n";
  cout << "      R       RD\n";
  cout << "\n";

  seed = 123456789;

  rlo2 = rlo - 2.0;
  rhi2 = rhi + 2.0;

  for ( test = 0; test < test_num; test++ )
  {
    r = r8_uniform_ab ( rlo2, rhi2, seed );
    rd = r8_to_r8_discrete ( r, rlo, rhi, ndx );
    cout << "  " << setw(14) << r
         << "  " << setw(14) << rd << "\n";
  }

  return;
}
//****************************************************************************80

void r8_to_i4_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_I4_TEST tests R8_TO_I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  int ix;
  int ixmax;
  int ixmin;
  double x;
  double xmax;
  double xmin;

  cout << "\n";
  cout << "R8_TO_I4_TEST\n";
  cout << "  R8_TO_I4 finds an integer IX in [IXMIN,IXMAX]\n";
  cout << "  corresponding to X in [XMIN,XMAX].\n";

  xmin = 2.5;
  x = 3.5;
  xmax = 5.5;

  ixmin = 10;
  ixmax = 40;

  ix = r8_to_i4 ( xmin, xmax, x, ixmin, ixmax );

  cout << "\n";
  cout << "   XMIN " <<  xmin << "   X = " <<  x << "  XMAX = "
    <<  xmax << "\n";
  cout << "  IXMIN " << ixmin << "  IX = " << ix << " IXMAX = "
    << ixmax << "\n";

  return;
}
//****************************************************************************80

void r8_uniform_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 1000

  int i;
  double max;
  double mean;
  double min;
  int n;
  int seed = 123456789;
  double x[N];
  double variance;

  cout << "\n";
  cout << "R8_UNIFORM_01_TEST\n";
  cout << "  R8_UNIFORM_01 samples a uniform random distribution in [0,1].\n";
  cout << "  distributed random numbers.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = r8_uniform_01 ( seed );
  }

  cout << "\n";
  cout << "  First few values:\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  min = r8vec_min ( N, x );
  max = r8vec_max ( N, x );
  mean = r8vec_mean ( N, x );
  variance = r8vec_variance ( N, x );

  cout << "\n";
  cout << "  Number of samples was " << N << "\n";
  cout << "  Minimum value was " << min << "\n";
  cout << "  Maximum value was " << max << "\n";
  cout << "  Average value was " << mean << "\n";
  cout << "  Variance was      " << variance << "\n";

  return;
# undef N
}
//****************************************************************************80

void r8_uniform_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB_TEST tests R8_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  int i;
  int seed;

  b = 10.0;
  c = 25.0;
  seed = 17;

  cout << "\n";
  cout << "R8_UNIFORM_AB_TEST\n";
  cout << "  R8_UNIFORM_AB produces a random real in a given range.\n";
  cout << "\n";
  cout << "  Using range " << b << " <= A <= " << c << ".\n";
  cout << "\n";

  cout << "\n";
  cout << "      I       A\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    a = r8_uniform_ab ( b, c, seed );
    cout << setw ( 6 )  << i << " "
         << setw ( 10 ) << a << "\n";
  }

  return;
}
//****************************************************************************80

void r8_walsh_1d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_WALSH_1D_TEST tests R8_WALSH_1D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double w0;
  double wm1;
  double wm2;
  double wm3;
  double wp1;
  double wp2;
  double x;

  cout << "\n";
  cout << "R8_WALSH_1D_TEST\n";
  cout << "  R8_WALSH_1D evaluates 1D Walsh functions:\n";
  cout << "\n";
  cout << "  X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)\n";
  cout << "\n";

  for ( i = 0; i <= 32; i++ )
  {
    x = ( double ) ( i ) / 4.0;

    wp2 = r8_walsh_1d ( x,  2 );
    wp1 = r8_walsh_1d ( x,  1 );
    w0  = r8_walsh_1d ( x,  0 );
    wm1 = r8_walsh_1d ( x, -1 );
    wm2 = r8_walsh_1d ( x, -2 );
    wm3 = r8_walsh_1d ( x, -3 );

    cout << "  " << setw(10) << x
         << "  " << setw(2)  << wp2
         << "  " << setw(2)  << wp1
         << "  " << setw(2)  << w0
         << "  " << setw(2)  << wm1
         << "  " << setw(2)  << wm2
         << "  " << setw(2)  << wm3 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_wrap_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_WRAP_TEST tests R8_WRAP;
//
//  Discussion:
//
//    Apparently if you turn on high precision somewhere in COUT, you're
//    stuck with it forever...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a = - 2.0;
  double b = 12.0;
  double r;
  double r2;
  double rhi = 6.5;
  double rlo = 3.0;
  int seed;
  int test;
  int test_num = 20;

  cout << "\n";
  cout << "R8_WRAP_TEST\n";
  cout << "  R8_WRAP \"wraps\" an R8 to lie within an interval:\n";
  cout << "\n";
  cout << "  Wrapping interval is " << rlo << ", " << rhi << "\n";
  cout << "\n";
  cout << "      R      R8_WRAP ( R )\n";
  cout << "\n";
  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    r = r8_uniform_ab ( a, b, seed );
    r2 = r8_wrap ( r, rlo, rhi );
    cout << "  " << setw(10) << r
         << "  " << setw(10) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void r82col_print_part_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R82COL_PRINT_PART_TEST tests R82COL_PRINT_PART.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int max_print;
  int n = 10;
  double v[10*2] = {
    11.0,  21.0, 31.0, 41.0, 51.0, 61.0, 71.0, 81.0, 91.0, 101.0, 
    12.0,  22.0, 32.0, 42.0, 52.0, 62.0, 72.0, 82.0, 92.0, 102.0 };

  cout << "\n";
  cout << "R82COL_PRINT_PART_TEST\n";
  cout << "  R82COL_PRINT_PART prints part of an R82COL.\n";

  max_print = 2;
  r82col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 2" );

  max_print = 5;
  r82col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 5" );

  max_print = 25;
  r82col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 25" );

  return;
}
//****************************************************************************80

void r82poly2_type_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R82POLY2_TYPE_TEST tests R82POLY2_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 12

  double a;
  double a_test[TEST_NUM] = {
    9.0, 4.0, 9.0,  1.0, 0.0,
    1.0, 0.0, 0.0,  0.0, 0.0,
    0.0, 0.0 };
  double b;
  double b_test[TEST_NUM] = {
    -4.0, 1.0,  16.0,  1.0,  0.0,
     2.0, 1.0,   1.0,  1.0,  0.0,
     0.0, 0.0 };
  double c;
  double c_test[TEST_NUM] = {
     0.0, -4.0,   0.0,   0.0, 1.0,
     0.0,  0.0,   0.0,  0.0,  0.0,
     0.0,  0.0 };
  double d;
  double r8_test[TEST_NUM] = {
    -36.0,  3.0,  36.0,  -6.0, 3.0,
    -2.0,   0.0,   0.0,  0.0,  2.0,
     0.0, 0.0 };
  double e;
  double e_test[TEST_NUM] = {
    -24.0, -4.0, -32.0, -10.0, -1.0,
     16.0, -6.0, -6.0, -2.0, -1.0,
     0.0, 0.0 };
  double f;
  double f_test[TEST_NUM] = {
    -36.0,  1.0, -92.0, 115.0, -3.0,
     33.0, +8.0, 10.0,  +1.0,  1.0,
      0.0, 1.0 };
  int test;
  int type;

  cout << "\n";
  cout << "R82POLY2_TYPE_TEST\n";
  cout << "  R82POLY2_TYPE determines the type of a second order\n";
  cout << "  equation in two variables.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];
    e = e_test[test];
    f = f_test[test];

    cout << "\n";

    r82poly2_print ( a, b, c, d, e, f );

    type = r82poly2_type ( a, b, c, d, e, f );

    cout << "  Type = " << type << "\n";

    r82poly2_type_print ( type );
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r82row_order_type_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R82ROW_ORDER_TYPE_TEST tests R82ROW_ORDER_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define TEST_NUM 10

  int i;
  int j;
  int order;
  int seed = 123456789;
  int test;
  double *x;

  cout << "\n";
  cout << "R82ROW_ORDER_TYPE_TEST\n";
  cout << "  R82ROW_ORDER_TYPE classifies an R8VEC as\n";
  cout << "  -1: no order\n";
  cout << "   0: all equal;\n";
  cout << "   1: ascending;\n";
  cout << "   2: strictly ascending;\n";
  cout << "   3: descending;\n";
  cout << "   4: strictly descending.\n";
  cout << "\n";

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    x = r8mat_uniform_01_new ( 2, N, seed );

    for ( j = 0; j < N; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        x[i+j*2] = ( double ) ( r8_nint ( 3.0 * x[i+j*2] ) );
      }
    }
    order = r82row_order_type ( N, x );

    cout << "  Order type = " << order << "\n";

    r82row_print ( N, x, " " );

    delete [] x;
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void r82row_part_quick_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R82ROW_PART_QUICK_A_TEST tests R82ROW_PART_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  double *a;
  double b = 0.0E+00;
  double c = 2.0E+00;
  int i;
  int l;
  int r;
  int seed = 123456789;

  cout << "\n";
  cout << "R82ROW_PART_QUICK_A_TEST\n";
  cout << "  R82ROW_PART_QUICK_A reorders an R82ROW\n";
  cout << "  as part of a quick sort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8mat_uniform_ab_new ( 2, N, b, c, seed );

  r82row_print ( N, a, "  Before rearrangment:" );

  r82row_part_quick_a ( N, a, l, r );

  cout << "\n";
  cout << "  Rearranged array\n";
  cout << "  Left index =  " << l << "\n";
  cout << "  Key index =   " << l+1 << "\n";
  cout << "  Right index = " << r << "\n";

  r82row_print ( l,     a,         "  Left half:" );
  r82row_print ( 1,     a+2*l,     "  Key:" );
  r82row_print ( N-l-1, a+2*(l+1), "  Right half:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r82row_print_part_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R82ROW_PRINT_PART_TEST tests R82ROW_PRINT_PART.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int max_print;
  int n = 10;
  double v[2*10] = {
     11.0,  21.0, 
     12.0,  22.0, 
     13.0,  23.0, 
     14.0,  24.0, 
     15.0,  25.0, 
     16.0,  26.0, 
     17.0,  27.0, 
     18.0,  28.0, 
     19.0,  29.0, 
     20.0,  30.0 };

  cout << "\n";
  cout << "R82ROW_PRINT_PART_TEST\n";
  cout << "  R82ROW_PRINT_PART prints part of an R82ROW.\n";

  max_print = 2;
  r82row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 2" );

  max_print = 5;
  r82row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 5" );

  max_print = 25;
  r82row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 25" );

  return;
}
//****************************************************************************80

void r82row_sort_heap_index_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R82ROW_SORT_HEAP_INDEX_A_TEST tests R82ROW_SORT_HEAP_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int *indx;
  int seed = 123456789;

  cout << "\n";
  cout << "R82ROW_SORT_HEAP_INDEX_A_TEST\n";
  cout << "  R82ROW_SORT_HEAP_INDEX_A index sorts an R82ROW\n";
  cout << "  using heapsort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8mat_uniform_ab_new ( 2, N, b, c, seed );
//
//  Give a few elements the same first component.
//
  a[0+2*2] = a[0+4*2];
  a[0+3*2] = a[0+11*2];
//
//  Give a few elements the same second component.
//
  a[1+5*2] = a[1+0*2];
  a[1+1*2] = a[1+8*2];
//
//  Make two entries equal.
//
  a[0+6*2] = a[0+10*2];
  a[1+6*2] = a[1+10*2];

  r82row_print ( N, a, "  Before rearrangement:" );

  indx = r82row_sort_heap_index_a ( N, a );

  cout << "\n";
  cout << "         I     Index   A(Index)\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << indx[i]
         << "  " << setw(12) << a[0+indx[i]*2]
         << "  " << setw(12) << a[1+indx[i]*2] << "\n";
  }

  r82row_permute ( N, indx, a );

  r82row_print ( N, a, "  After rearrangement by R82ROW_PERMUTE:" );

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void r82row_sort_quick_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R82ROW_SORT_QUICK_A_TEST tests R82ROW_SORT_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "R82ROW_SORT_QUICK_A_TEST\n";
  cout << "  R82ROW_SORT_QUICK_A sorts an R82ROW\n";
  cout << "  as part of a quick sort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8mat_uniform_ab_new ( 2, N, b, c, seed );
//
//  For better testing, give a few elements the same first component.
//
  a[2*(3-1)+0] = a[2*(5-1)+0];
  a[2*(4-1)+0] = a[2*(12-1)+0];
//
//  Make two entries equal.
//
  a[2*(7-1)+0] = a[2*(11-1)+0];
  a[2*(7-1)+1] = a[2*(11-1)+1];

  r82row_print ( N, a, "  Before sorting:" );

  r82row_sort_quick_a ( N, a );

  r82row_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r83col_print_part_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R83COL_PRINT_PART_TEST tests R83COL_PRINT_PART.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int max_print;
  int n = 10;
  double v[10*3] = {
    11.0,  21.0, 31.0, 41.0, 51.0, 61.0, 71.0, 81.0, 91.0, 101.0, 
    12.0,  22.0, 32.0, 42.0, 52.0, 62.0, 72.0, 82.0, 92.0, 102.0,
    13.0,  23.0, 33.0, 43.0, 53.0, 63.0, 73.0, 83.0, 93.0, 103.0 };

  cout << "\n";
  cout << "R83COL_PRINT_PART_TEST\n";
  cout << "  R83COL_PRINT_PART prints part of an R83COL.\n";

  max_print = 2;
  r83col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 2" );

  max_print = 5;
  r83col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 5" );

  max_print = 25;
  r83col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 25" );

  return;
}
//****************************************************************************80

void r83row_print_part_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R83ROW_PRINT_PART_TEST tests R83ROW_PRINT_PART.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int max_print;
  int n = 10;
  double v[3*10] = {
     11.0,  21.0,  31.0,
     12.0,  22.0,  32.0,
     13.0,  23.0,  33.0,
     14.0,  24.0,  34.0,
     15.0,  25.0,  35.0,
     16.0,  26.0,  36.0,
     17.0,  27.0,  37.0,
     18.0,  28.0,  38.0,
     19.0,  29.0,  39.0,
     20.0,  30.0,  40.0, };

  cout << "\n";
  cout << "R83ROW_PRINT_PART_TEST\n";
  cout << "  R83ROW_PRINT_PART prints part of an R83ROW.\n";

  max_print = 2;
  r83row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 2" );

  max_print = 5;
  r83row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 5" );

  max_print = 25;
  r83row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 25" );

  return;
}
//****************************************************************************80

void r8block_expand_linear_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLOCK_EXPAND_LINEAR_TEST tests R8BLOCK_EXPAND_LINEAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define L 4
# define M 3
# define N 2

  int l2;
  int lfat = 1;
  int m2;
  int mfat = 2;
  int n2;
  int nfat = 1;
  double x[L*M*N] = {
        1.0,  2.0,  3.0,   4.0,  1.0,
        4.0,  9.0, 16.0,   1.0,  8.0,
       27.0, 64.0,  2.0,   4.0,  6.0,
        8.0,  2.0,  8.0,  18.0, 32.0,
        2.0, 16.0, 54.0, 128.0 };
  double *xfat;

  l2 = ( L - 1 ) * ( lfat + 1 ) + 1;
  m2 = ( M - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( N - 1 ) * ( nfat + 1 ) + 1;

  cout << "\n";
  cout << "R8BLOCK_EXPAND_LINEAR_TEST\n";
  cout << "  R8BLOCK_EXPAND_LINEAR linearly interpolates new data\n";
  cout << "  between old values in a 3D block.\n";

  r8block_print ( L, M, N, x, "  Original block:" );

  cout << "\n";
  cout << "  LFAT = " << lfat << "\n";
  cout << "  MFAT = " << mfat << "\n";
  cout << "  NFAT = " << nfat << "\n";

  xfat = r8block_expand_linear ( L, M, N, x, lfat, mfat, nfat );

  r8block_print ( l2, m2, n2, xfat, "  Fattened block:" );

  delete [] xfat;

  return;
# undef L
# undef M
# undef N
}
//****************************************************************************80

void r8block_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0363 tests R8BLOCK_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double ***a;
  double ***b;
  int i;
  int j;
  int k;
  int l;
  int m;
  int n;

  cout << "\n";
  cout << "R8BLOCK_NEW_TEST:\n";
  cout << "  R8BLOCK_NEW dynamically creates a 3D array.\n";
  cout << "  Array entries can be addressed using the\n";
  cout << "  notation \"a[i][j][k]\".\n";
//
//  These dimensions could be entered by the user; they could depend on
//  some other calculation; or they could be changed repeatedly during this
//  computation, as long as old memory is deleted by R8BLOCK_DELETE and new memory
//  requested by R8BLOCK_NEW.
//
  l = 2;
  m = 3;
  n = 2;
//
//  Allocate memory.
//
  cout << "\n";
  cout << "  Allocating memory for array A of size " << l << " by " << m << " by " << n << ".\n";

  a = r8block_new ( l, m, n );

  cout << "\n";
  cout << "  Assigning values to A.\n";
//
//  Store values in A.
//
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        a[i][j][k] = ( double ) ( 100 * i + 10 * j + k );
      }
    }
  }
//
//  Print A.
//
  cout << "\n";
  cout << "  Dynamically allocated matrix A:\n";
  cout << "\n";
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        cout << "  " << setw(8) << a[i][j][k];
      }
      cout << "\n";
    }
    cout << "\n";
  }
//
//  Free memory.
//
  r8block_delete ( l, m, n, a );

  return;
}
//****************************************************************************80

void r8block_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLOCK_PRINT_TEST tests R8BLOCK_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2012
//
//  Author:
//
//    John Burkardt
//
{
# define L 4
# define M 3
# define N 2

  double x[L*M*N] = {
        1.0,  2.0,  3.0,   4.0,  1.0,
        4.0,  9.0, 16.0,   1.0,  8.0,
       27.0, 64.0,  2.0,   4.0,  6.0,
        8.0,  2.0,  8.0,  18.0, 32.0,
        2.0, 16.0, 54.0, 128.0 };

  cout << "\n";
  cout << "R8BLOCK_PRINT_TEST\n";
  cout << "  R8BLOCK_PRINT prints an R8BLOCK.\n";

  r8block_print ( L, M, N, x, "  The 3D array:" );

  return;
# undef L
# undef M
# undef N
}
//****************************************************************************80

void r8cmat_to_r8mat_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8CMAT_TO_R8MAT_NEW_TEST tests R8CMAT_TO_R8MAT_NEW;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double **b;
  double *c;
  int m = 5;
  int n = 4;

  cout << "\n";
  cout << "R8CMAT_TO_R8MAT_NEW_TEST\n";
  cout << "  R8CMAT_TO_R8MAT_NEW converts an R8CMAT to an R8MAT.\n";
  cout << "\n";
  cout << "  Data is of order (" << m << "," << n << ".\n";
//
//  Set the R8MAT.
//
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  The R8MAT A:" );
//
//  Convert.
//
  b = r8mat_to_r8cmat_new ( m, n, a );
  r8cmat_print ( m, n, b, "  The R8CMAT B:" );
//
//  Recover the matrix.
//
  c = r8cmat_to_r8mat_new ( m, n, b );
  r8mat_print ( m, n, c, "  The R8MAT C:" );
//
//  Free memory.
//
  delete [] a;
  r8cmat_delete ( m, n, b );
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void r8col_find_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_FIND_TEST tests R8COL_FIND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  int col;
  double dtab[M*N];
  double r8vec[M];
  int i;
  int j;
  int k;

  k = 1;

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      dtab[i+j*M] = ( double ) k;
      if ( j == 2 )
      {
        r8vec[i] = ( double ) k;
      }
      k = k + 1;
    }
  }

  col = r8col_find ( M, N, dtab, r8vec );

  cout << "\n";
  cout << "R8COL_FIND_TEST\n";
  cout << "  R8COL_FIND finds a column in a table matching\n";
  cout << "  a given set of data.\n";
  cout << "\n";
  cout << "  R8COL_FIND returns COL = " << col << "\n";

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_insert_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_INSERT_TEST tests R8COL_INSERT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N_MAX 10

  double a[M*N_MAX] = {
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0 };
  int col;
  double r8vec1[M] = { 3.0, 7.0, 11.0 };
  double r8vec2[M] = { 3.0, 4.0, 18.0 };
  int n;

  cout << "\n";
  cout << "R8COL_INSERT_TEST\n";
  cout << "  R8COL_INSERT inserts new columns into a sorted R8COL.\n";

  n = 4;

  r8mat_print ( M, n, a, "  The unsorted matrix:" );

  r8col_sort_heap_a ( M, n, a );

  r8mat_print ( M, n, a, "  The sorted matrix:" );

  r8vec_print ( M, r8vec1, "  New column:" );

  col = r8col_insert ( N_MAX, M, n, a, r8vec1 );

  if ( col < 0 )
  {
    cout << "\n";
    cout << "  The data was already in column " << abs ( col ) << "\n";
  }
  else
  {
    r8mat_print ( M, n, a, "  The updated matrix:" );
  }

  r8vec_print ( M, r8vec2, "  New column:" );

  col = r8col_insert ( N_MAX, M, n, a, r8vec2 );

  if ( col < 0 )
  {
    cout << "\n";
    cout << "  The data was already in column " << abs ( col ) << "\n";
  }
  else
  {
    r8mat_print ( M, n, a, "  The updated matrix:" );
  }

  return;
# undef M
# undef N_MAX
}
//****************************************************************************80

void r8col_sort_heap_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SORT_HEAP_A_TEST tests R8COL_SORT_HEAP_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N_MAX 10

  double a[M*N_MAX] = {
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0 };
  int col;
  double r8vec1[M] = { 3.0, 7.0, 11.0 };
  double r8vec2[M] = { 3.0, 4.0, 18.0 };
  int n;

  cout << "\n";
  cout << "R8COL_SORT_HEAP_A_TEST\n";
  cout << "  R8COL_SORT_HEAP_A ascending heap sorts a table of columns.\n";

  n = 4;

  r8mat_print ( M, n, a, "  The unsorted matrix:" );

  r8col_sort_heap_a ( M, n, a );

  r8mat_print ( M, n, a, "  The sorted matrix:" );

  return;
# undef M
# undef N_MAX
}
//****************************************************************************80

void r8col_sort_heap_index_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SORT_HEAP_INDEX_A_TEST tests R8COL_SORT_HEAP_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 15

  double a[M*N] = {
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    3.0,  4.0, 18.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int i;
  int *indx;
  int j;
  int j2;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8COL_SORT_HEAP_INDEX_A_TEST\n";
  cout << "  R8COL_SORT_HEAP_INDEX_A computes an index vector which\n";
  cout << "  ascending sorts an R8COL.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  indx = r8col_sort_heap_index_a ( m, n, a );

  cout << "\n";
  cout << "  The implicitly sorted R8COL (transposed)\n";
  cout << "\n";

  for ( j = 0; j < n; j++ )
  {
    j2 = indx[j];
    cout << "  " << setw(4) << j2 << ":";
    for ( i = 0; i < m; i++ )
    {
      cout << "  " << setw(10) << a[i+j2*m];
    }
    cout << "\n";
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_sort_quick_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SORT_QUICK_A_TEST tests R8COL_SORT_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 10

  double *a;
  double b = 0.0;
  double c = 10.0;
  int seed;

  cout << "\n";
  cout << "R8COL_SORT_QUICK_A_TEST\n";
  cout << "  R8COL_SORT_QUICK_A sorts a table of columns.\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  The unsorted matrix:" );

  r8col_sort_quick_a ( M, N, a );

  r8mat_print ( M, N, a, "  The sorted matrix:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_sorted_tol_unique_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SORTED_TOL_UNIQUE_TEST tests R8COL_SORTED_TOL_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  cout << "\n";
  cout << "R8COL_SORTED_TOL_UNIQUE_TEST\n";
  cout << "  R8COL_SORTED_TOL_UNIQUE finds tolerably unique columns \n";
  cout << "  in a sorted R8COL.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  cout << "\n";
  cout << "  Using tolerance = " << tol << "\n";

  unique_num = r8col_sorted_tol_unique ( m, n, a, tol );

  cout << "\n";
  cout << "  Number of tolerably unique columns is " << unique_num << "\n";

  r8mat_transpose_print ( m, unique_num, a,
    "  The sorted tolerably unique R8COL (transposed):" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_sorted_unique_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SORTED_UNIQUE_COUNT_TEST tests R8COL_SORTED_UNIQUE_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  cout << "\n";
  cout << "R8COL_SORTED_UNIQUE_COUNT_TEST\n";
  cout << "  R8COL_SORTED_UNIQUE_COUNT counts tolerably unique columns \n";
  cout << "  in a sorted R8COL.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  cout << "\n";
  cout << "  Using tolerance = " << tol << "\n";

  unique_num = r8col_sorted_tol_unique_count ( m, n, a, tol );

  cout << "\n";
  cout << "  Number of tolerably unique columns is " << unique_num << "\n";

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_sorted_tol_undex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SORTED_TOL_UNDEX_TEST tests R8COL_SORTED_TOL_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  double tol;
  int *undx;
  int unique_num;
  int *xdnu;

  cout << "\n";
  cout << "R8COL_SORTED_TOL_UNDEX_TEST\n";
  cout << "  R8COL_SORTED_TOL_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the tolerably unique columns of a sorted R8COL,\n";
  cout << "  and a map from the original R8COL to the (implicit)\n";
  cout << "  R8COL of sorted tolerably unique elements.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  cout << "\n";
  cout << "  Using tolerance = " << tol << "\n";

  n_unique = r8col_sorted_tol_unique_count ( m, n, a, tol );

  cout << "\n";
  cout << "  Number of tolerably unique columns is " << n_unique << "\n";

  au = new double[m*n_unique];
  undx = new int[n_unique];
  xdnu = new int[n];

  r8col_sorted_tol_undex ( m, n, a, n_unique, tol, undx, xdnu );

  cout << "\n";
  cout << "  XDNU points to the representative for each item.\n";
  cout << "  UNDX selects the representatives.\n";
  cout << "\n";
  cout << "     I  XDNU  UNDX\n";
  cout << "\n";
  for ( i = 0; i < n_unique; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << undx[i] << "\n";
  }
  for ( i = n_unique; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i] << "\n";
  }
  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au,
    "  The tolerably unique R8COL (transposed):" );

  delete [] au;
  delete [] undx;
  delete [] xdnu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_max_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_MAX_TEST tests R8COL_MAX;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *amax;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "R8COL_MAX_TEST\n";
  cout << "  R8COL_MAX computes maximums of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  amax = r8col_max ( M, N, a );

  r8vec_print ( N, amax, "  Column maximums:" );

  delete [] amax;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_mean_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_MEAN_TEST tests R8COL_MEAN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *colsum;
  int i;
  int j;
  int k;
  double *mean;

  cout << "\n";
  cout << "R8COL_MEAN_TEST\n";
  cout << "  R8COL_MEAN computes means of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  mean = r8col_mean ( M, N, a );

  r8vec_print ( N, mean, "  The column means:" );

  delete [] mean;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_min_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_MIN_TEST tests R8COL_MIN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *amin;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "R8COL_MIN_TEST\n";
  cout << "  R8COL_MIN computes minimums of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  amin = r8col_min ( M, N, a );

  r8vec_print ( N, amin, "  Column minimums:" );

  delete [] amin;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_permute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_PERMUTE_TEST tests R8COL_PERMUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 5

  double a[M*N] = {
    11.0, 21.0, 31.0,
    12.0, 22.0, 32.0,
    13.0, 23.0, 33.0,
    14.0, 24.0, 34.0,
    15.0, 25.0, 35.0 };
  int perm[N] = { 1, 3, 4, 0, 2 };


  cout << "\n";
  cout << "R8COL_PERMUTE_TEST\n";
  cout << "  R8COL_PERMUTE permutes an R8COL in place.\n";

  r8mat_print ( M, N, a, "  A (unpermuted):" );

  i4vec_print ( N, perm, "  The (column) permutation vector:" );

  r8col_permute ( M, N, perm, a );

  r8mat_print ( M, N, a, "  A (permuted):" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_sortr_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SORTR_A_TEST tests R8COL_SORTR_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define N 3

  double *a;
  double b = 0.0;
  double c = 10.0;
  int key;
  int seed;

  cout << "\n";
  cout << "R8COL_SORTR_A_TEST\n";
  cout << "  R8COL_SORTR_A is given an array, and reorders\n";
  cout << "  it so that a particular column is sorted.\n";

  key = 2;
  cout << "\n";
  cout << "  Here, the special column is " << key << "\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  Unsorted array:" );

  r8col_sortr_a ( M, N, a, key );

  r8mat_print ( M, N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_sum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SUM_TEST tests R8COL_SUM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *colsum;
  int i;
  int j;
  int k;
  double *mean;

  cout << "\n";
  cout << "R8COL_SUM_TEST\n";
  cout << "  R8COL_SUM computes sums of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  colsum = r8col_sum ( M, N, a );

  r8vec_print ( N, colsum, "  The column sums" );

  delete [] colsum;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_swap_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SWAP_TEST tests R8COL_SWAP;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int icol1;
  int icol2;
  int m = 3;
  int n = 4;

  cout << "\n";
  cout << "R8COL_SWAP_TEST\n";
  cout << "  R8COL_SWAP swaps two columns of an R8COL;\n";

  a = r8mat_indicator_new ( m, n );

  r8mat_print ( m, n, a, "  The array:" );

  icol1 = 1;
  icol2 = 3;

  cout << "\n";
  cout << "  Swap columns " << icol1 << " and " << icol2 << ":\n";

  r8col_swap ( m, n, a, icol1, icol2 );

  r8mat_print ( m, n, a, "  The updated matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void r8col_to_r8vec_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_TO_R8VEC_TEST tests R8COL_TO_R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  double *x;

  cout << "\n";
  cout << "R8COL_TO_R8VEC_TEST\n";
  cout << "  R8COL_TO_R8VEC converts an array of columns to a vector.\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = ( double ) ( 10 * i + j );
    }
  }

  r8mat_print ( M, N, a, "  The array of columns:" );

  x = r8col_to_r8vec ( M, N, a );

  r8vec_print ( M*N, x, "  The resulting vector of columns:" );

  delete [] x;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_tol_undex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_TOL_UNDEX_TEST tests R8COL_TOL_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  double tol;
  int *undx;
  int unique_num;
  int *xdnu;

  cout << "\n";
  cout << "R8COL_TOL_UNDEX_TEST\n";
  cout << "  R8COL_TOL_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the tolerably unique columns of an R8COL,\n";
  cout << "  and a map from the original R8COL to the (implicit)\n";
  cout << "  R8COL of sorted tolerably unique elements.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  tol = 0.25;

  cout << "\n";
  cout << "  Using tolerance = " << tol << "\n";

  n_unique = r8col_tol_unique_count ( m, n, a, tol );

  cout << "\n";
  cout << "  Number of tolerably unique columns is " << n_unique << "\n";

  au = new double[m*n_unique];
  undx = new int[n_unique];
  xdnu = new int[n];

  r8col_tol_undex ( m, n, a, n_unique, tol, undx, xdnu );

  cout << "\n";
  cout << "  XDNU points to the representative for each item.\n";
  cout << "  UNDX selects the representatives.\n";
  cout << "\n";
  cout << "     I  XDNU  UNDX\n";
  cout << "\n";
  for ( i = 0; i < n_unique; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << undx[i] << "\n";
  }
  for ( i = n_unique; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i] << "\n";
  }

  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au,
    "  The tolerably unique R8COL (transposed):" );

  delete [] au;
  delete [] undx;
  delete [] xdnu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_undex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_UNDEX_TEST tests R8COL_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  int *undx;
  int unique_num;
  int *xdnu;

  cout << "\n";
  cout << "R8COL_UNDEX_TEST\n";
  cout << "  R8COL_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the unique columns of an (unsorted) R8COL,\n";
  cout << "  and a map from the original R8COL to the (implicit)\n";
  cout << "  R8COL of sorted unique elements.\n";

  r8mat_transpose_print ( m, n, a, "  The R8COL (transposed):" );

  n_unique = r8col_unique_count ( m, n, a );

  cout << "\n";
  cout << "  Number of unique columns is " << n_unique << "\n";

  au = new double[m*n_unique];
  undx = new int[n_unique];
  xdnu = new int[n];

  r8col_undex ( m, n, a, n_unique, undx, xdnu );

  cout << "\n";
  cout << "  XDNU points to the representative for each item.\n";
  cout << "  UNDX selects the representatives.\n";
  cout << "\n";
  cout << "     I  XDNU  UNDX\n";
  cout << "\n";
  for ( i = 0; i < n_unique; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << undx[i] << "\n";
  }
  for ( i = n_unique; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i] << "\n";
  }

  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au, "  The Unique R8COL (transposed):" );

  delete [] au;
  delete [] undx;
  delete [] xdnu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_unique_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_UNIQUE_COUNT_TEST tests R8COL_UNIQUE_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  cout << "\n";
  cout << "R8COL_UNIQUE_COUNT_TEST\n";
  cout << "  R8COL_UNIQUE_COUNT counts unique columns.\n";

  r8mat_transpose_print ( m, n, a, "  The R8COL (transposed):" );

  unique_num = r8col_unique_count ( m, n, a );

  cout << "\n";
  cout << "  Number of unique columns is " << unique_num << "\n";

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8col_variance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_VARIANCE_TEST tests R8COL_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *variance;

  cout << "\n";
  cout << "R8COL_VARIANCE_TEST\n";
  cout << "  R8COL_VARIANCE computes variances of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) ( k );
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  variance = r8col_variance ( M, N, a );

  cout << "\n";
  cout << "  Column  variance:\n";
  cout << "\n";

  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(6)  << j
         << "  " << setw(10) << variance[j] << "\n";
  }

  delete [] variance;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8int_to_i4int_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8INT_TO_I4INT_TEST tests R8INT_TO_I4INT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi = 11;
  int ilo = 1;
  int ir;
  double r;
  double r2;
  double rhi = 200.0;
  double rhi2;
  double rlo = 100.0;
  double rlo2;
  int seed;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "R8INT_TO_I4INT_TEST\n";
  cout << "  For data in an interval,\n";
  cout << "  R8INT_TO_I4INT converts a real to an integer.\n";
  cout << "\n";
  cout << "  Integer interval: [" << ilo << ", " << ihi << "]\n";
  cout << "  Real interval:    [" << rlo << ", " << rhi << "]\n";
  cout << "\n";
  cout << "  R   I(R)  R(I(R))\n";
  cout << "\n";

  seed = 123456789;

  rlo2 = rlo - 15.0;
  rhi2 = rhi + 15.0;

  for ( test = 1; test <= test_num; test++ )
  {
    r = r8_uniform_ab ( rlo2, rhi2, seed );
    ir = r8int_to_i4int ( rlo, rhi, r, ilo, ihi );
    r2 = i4int_to_r8int ( ilo, ihi, ir, rlo, rhi );
    cout << "  " << setw(12) << r
         << "  " << setw(6)  << ir
         << "  " << setw(12) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8r8vec_index_insert_unique_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8R8_INDEX_INSERT_UNIQUE_TEST tests R8R8VEC_INDEX_INSERT_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 30

  int i;
  int ierror;
  int indx[N_MAX];
  int ival;
  int n;
  int seed;
  double x[N_MAX];
  double x_max = 4.0;
  double x_min = 1.0;
  double xval;
  double y[N_MAX];
  double y_max = 3.0;
  double y_min = 1.0;
  double yval;

  n = 0;

  cout << "\n";
  cout << "R8R8_INDEX_INSERT_UNIQUE_TEST\n";
  cout << "  R8R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an\n";
  cout << "  index sorted array.\n";
  cout << "\n";
  cout << "  Generate " << N_MAX << " random values:\n";
  cout << "\n";
  cout << "    XVAL    YVAL   Index\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( x_min, x_max, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    yval = r8_uniform_ab ( y_min, y_max, seed );
    yval = ( double ) ( r8_nint ( yval ) );

    r8r8vec_index_insert_unique ( N_MAX, n, x, y, indx, xval, yval,
      ival, ierror );

    cout << "  " << setw(6)  << ival
         << "  " << setw(12) << xval
         << "  " << setw(12) << yval << "\n";
  }

  cout << "\n";
  cout << "  Vector of unique X Y values:\n";
  cout << "\n";
  cout << "  I  X(I)   Y(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6)  << i+1
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << y[i] << "\n";
  }

  cout << "\n";
  cout << "  X, Y sorted by index\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(INDX(I))  Y(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i+1
         << "  " << setw(6) << indx[i]
         << "  " << setw(12) << x[indx[i]-1]
         << "  " << setw(12) << y[indx[i]-1] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8r8r8vec_index_insert_unique_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8R8R8VEC_INDEX_INSERT_UNIQUE_TEST tests R8R8R8VEC_INDEX_INSERT_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 30

  int i;
  int ierror;
  int indx[N_MAX];
  int ival;
  int n;
  int seed;
  double x[N_MAX];
  double xval;
  double y[N_MAX];
  double yval;
  double z[N_MAX];
  double zval;

  n = 0;

  cout << "\n";
  cout << "R8R8R8VEC_INDEX_INSERT_UNIQUE_TEST\n";
  cout << "  R8R8R8VEC_INDEX_INSERT_UNIQUE inserts unique values into\n";
  cout << "  an index sorted array.\n";
  cout << "\n";
  cout << "  Number of random values to generate = " << N_MAX << "\n";
  cout << "\n";
  cout << "    XVAL    YVAL  ZVAL  Index\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( 1.0, 4.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    yval = r8_uniform_ab ( 1.0, 3.0, seed );
    yval = ( double ) ( r8_nint ( yval ) );
    zval = r8_uniform_ab ( 1.0, 4.0, seed );
    zval = ( double ) ( r8_nint ( zval ) );

    r8r8r8vec_index_insert_unique ( N_MAX, n, x, y, z, indx,
      xval, yval, zval, ival, ierror );

    cout << "  " << setw(6) << xval
         << "  " << setw(6) << yval
         << "  " << setw(6) << zval
         << "  " << setw(6) << ival << "\n";
  }

  cout << "\n";
  cout << "  Vector of unique X Y Z values:\n";
  cout << "\n";
  cout << "  I  X(I)   Y(I)    Z(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << y[i]
         << "  " << setw(6) << z[i] << "\n";
  }

  cout << "\n";
  cout << "  X Y Z sorted by index:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[indx[i]-1]
         << "  " << setw(6) << y[indx[i]-1]
         << "  " << setw(6) << z[indx[i]-1] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8mat_cholesky_inverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_CHOLESKY_INVERSE_TEST tests R8MAT_CHOLESKY_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a2;
  double *a3;
  int i;
  int j;
  int n = 5;
  int test;

  cout << "\n";
  cout << "R8MAT_CHOLESKY_INVERSE_TEST\n";
  cout << "  For a positive definite symmetric matrix,\n";
  cout << "  R8MAT_CHOLESKY_INVERSE computes the inverse.\n";

  a = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a[i+j*n] = -1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  r8mat_print ( n, n, a, "  Matrix to be inverted:" );

  a2 = r8mat_copy_new ( n, n, a );

  r8mat_cholesky_inverse ( n, a2 );

  r8mat_print ( n, n, a2, "  Inverse matrix:" );

  a3 = r8mat_mm_new ( n, n, n, a2, a );
  
  r8mat_print ( n, n, a3, "  Product inv(A) * A:" );

  delete [] a;
  delete [] a2;
  delete [] a3;

  return;
}
//****************************************************************************80

void r8mat_cholesky_solve_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_CHOLESKY_SOLVE_TEST tests R8MAT_CHOLESKY_SOLVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a2;
  double *a3;
  double *b;
  double *d;
  int flag;
  int i;
  int j;
  double *l;
  int n = 5;
  double *r;
  int test;
  double *x;

  cout << "\n";
  cout << "R8MAT_CHOLESKY_SOLVE_TEST\n";
  cout << "  For a positive definite symmetric matrix,\n";
  cout << "  R8MAT_CHOLESKY_SOLVE solves a linear system\n";
  cout << "  using the lower Cholesky factorization.\n";

  a = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a[i+j*n] = -1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  r8mat_print ( n, n, a, "  Matrix to be factored:" );
//
//  Compute the Cholesky factor.
//
  l = r8mat_cholesky_factor ( n, a, flag );
  if ( flag != 0 )
  {
    cerr << "\n";
    cerr << "  R8MAT_CHOLESKY_FACTOR failed.\n";
    return;
  }
  r8mat_print ( n, n, l, "  Cholesky factor L:" );
  d = r8mat_mmt_new ( n, n, n, l, l );
  r8mat_print ( n, n, d, "  Product L * L':" );
//
//  Solve a system.
//
  b = new double[n];

  for ( i = 0; i < n - 1; i++ )
  {
    b[i] = 0.0;
  }
  b[n-1] = ( double ) ( n + 1 );

  r8vec_print ( n, b, "  Right hand side:" );

  x = r8mat_cholesky_solve ( n, l, b );

  r8vec_print ( n, x, "  Computed solution:" );

  delete [] a;
  delete [] b;
  delete [] d;
  delete [] l;
  delete [] x;

  return;
}
//****************************************************************************80

void r8mat_cholesky_solve_upper_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_CHOLESKY_SOLVE_UPPER_TEST tests the Cholesky routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a2;
  double *a3;
  double *b;
  double *d;
  int flag;
  int i;
  int j;
  double *l;
  int n = 5;
  double *r;
  int test;
  double *x;

  cout << "\n";
  cout << "R8MAT_CHOLESKY_SOLVE_UPPER_TEST\n";
  cout << "  For a positive definite symmetric matrix,\n";
  cout << "  R8MAT_CHOLESKY_SOLVE_UPPER solves a linear system\n";
  cout << "  using the upper Cholesky factorization.\n";

  a = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a[i+j*n] = -1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  r8mat_print ( n, n, a, "  Matrix to be factored:" );
//
//  Compute the Cholesky factor.
//
  r = r8mat_cholesky_factor_upper ( n, a, flag );
  if ( flag != 0 )
  {
    cerr << "\n";
    cerr << "  R8MAT_CHOLESKY_FACTOR_UPPER failed.\n";
    return;
  }
  r8mat_print ( n, n, r, "  Cholesky factor R:" );
  d = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, d, "  Product R' * R:" );
//
//  Solve a system.
//
  b = new double[n];

  for ( i = 0; i < n - 1; i++ )
  {
    b[i] = 0.0;
  }
  b[n-1] = ( double ) ( n + 1 );

  r8vec_print ( n, b, "  Right hand side:" );

  x = r8mat_cholesky_solve_upper ( n, r, b );

  r8vec_print ( n, x, "  Computed solution:" );

  delete [] a;
  delete [] b;
  delete [] d;
  delete [] r;
  delete [] x;

  return;
}
//****************************************************************************80

void r8mat_det_2d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_2D_TEST tests R8MAT_DET_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0 };

  cout << "\n";
  cout << "R8MAT_DET_2D_TEST\n";
  cout << "  R8MAT_DET_2D: determinant of a 2 by 2 matrix;\n";

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_2d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  cout << "\n";
  cout << "  R8MAT_DET_2D computes determinant:" << det << "\n";
//
//  Special formula for the determinant of a Vandermonde matrix:
//
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  cout << "  Exact determinant is " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8mat_det_3d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_3D_TEST tests R8MAT_DET_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0 };

  cout << "\n";
  cout << "R8MAT_DET_3D_TEST\n";
  cout << "  R8MAT_DET_3D: determinant of a 3 by 3 matrix;\n";

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_3d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  cout << "\n";
  cout << "  R8MAT_DET_3D computes determinant:" << det << "\n";
//
//  Special formula for the determinant of a Vandermonde matrix:
//
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  cout << "  Exact determinant is " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8mat_det_4d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_4D_TEST tests R8MAT_DET_4D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0, 2.0 };

  cout << "\n";
  cout << "R8MAT_DET_4D_TEST\n";
  cout << "  R8MAT_DET_4D determinant of a 4 by 4 matrix;\n";

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_4d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  cout << "\n";
  cout << "  R8MAT_DET_4D computes determinant:" << det << "\n";
//
//  Special formula for the determinant of a Vandermonde matrix:
//
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  cout << "  Exact determinant is " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8mat_det_5d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_5D_TEST tests R8MAT_DET_5D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0, 2.0, 3.0 };

  cout << "\n";
  cout << "R8MAT_DET_5D_TEST\n";
  cout << "  R8MAT_DET_5D determinant of a 5 by 5 matrix;\n";

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_5d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  cout << "\n";
  cout << "  R8MAT_DET_5D computes determinant:" << det << "\n";
//
//  Special formula for the determinant of a Vandermonde matrix:
//
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  cout << "  Exact determinant is " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8mat_expand_linear_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXPAND_LINEAR_TEST tests R8MAT_EXPAND_LINEAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 3

  int m2;
  int mfat = 2;
  int n2;
  int nfat = 1;
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double x[M*N] = {
    1.0, 2.0, 3.0, 4.0, 1.0,
    4.0, 9.0, 16.0, 1.0, 8.0,
    27.0, 64.0 };
  double *xfat;

  m2 = ( M - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( N - 1 ) * ( nfat + 1 ) + 1;

  cout << "\n";
  cout << "R8MAT_EXPAND_LINEAR_TEST\n";
  cout << "  R8MAT_EXPAND_LINEAR linearly interpolates new data\n";
  cout << "  between old values in a matrix.\n";

  r8mat_print ( M, N, x, "  Original matrix:" );

  cout << "\n";
  cout << "  MFAT = " << mfat << "\n";
  cout << "  NFAT = " << nfat << "\n";

  xfat = r8mat_expand_linear ( M, N, x, mfat, nfat );

  r8mat_print ( m2, n2, xfat, "  Fattened matrix:" );

  delete [] xfat;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_expand_linear2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXPAND_LINEAR2_TEST tests R8MAT_EXPAND_LINEAR2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 2

  double a[M*N];
  double *a2;
  int i;
  int j;
  int m2 = 10;
  int n2 = 5;

  cout << "\n";
  cout << "R8MAT_EXPAND_LINEAR2_TEST\n";
  cout << "  R8MAT_EXPAND_LINEAR2 fills in a large array by\n";
  cout << "  interpolating data from a small array.\n";
  cout << "\n";
  cout << "  Original matrix has dimensions:\n";
  cout << "\n";
  cout << "  M = " << M << ", N = " << N << "\n";
  cout << "\n";
  cout << "  Expanded matrix has dimensions:\n";
  cout << "\n";
  cout << "  M2 = " << m2 << ", N2 = " << n2 << "\n";

  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      a[i-1+(j-1)*M] = 10.0 * ( double ) ( i ) + ( double ) ( j );
    }
  }

  r8mat_print ( M, N, a, "  The little matrix A:" );

  a2 = r8mat_expand_linear2 ( M, N, a, m2, n2 );

  r8mat_print ( m2, n2, a2, "  Expanded array A2:" );

  delete [] a2;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_fs_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_FS_NEW_TEST tests R8MAT_FS_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int k;
  int n = N;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "R8MAT_FS_NEW_TEST\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8MAT_FS_NEW factors and solves a linear system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << n << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( n, n, seed );
//
//  Set the desired solutions.
//
  b = new double[n];

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*n] * x[j];
    }
  }
//
//  Factor and solve the system.
//
  delete [] x;

  x = r8mat_fs_new ( n, a, b );
  
  r8vec_print ( n, x, "  Solution:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void r8mat_fss_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_FSS_NEW_TEST tests R8MAT_FSS_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define NB 3

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int k;
  int n = N;
  int nb = NB;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "R8MAT_FSS_NEW_TEST\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8MAT_FSS_NEW factors and solves multiple linear systems.\n";
  cout << "\n";
  cout << "  Matrix order N = " << n << "\n";
  cout << "  Number of systems NB = " << nb << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( n, n, seed );
//
//  Set the desired solutions.
//
  b = new double[n * nb];

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0;
  }
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  k = 1;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( i % 3 ) + 1;
  }
  k = 2;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
//
//  Factor and solve the system.
//
  delete [] x;

  x = r8mat_fss_new ( n, a, nb, b );
  
  r8mat_print ( n, nb, x, "  Solutions:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
# undef NB
}
//****************************************************************************80

void r8mat_givens_post_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_GIVENS_POST_TEST tests R8MAT_GIVENS_POST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N*N];
  double *ag;
  int col;
  double *g;
  int i;
  int j;
  int row;

  cout << "\n";
  cout << "R8MAT_GIVENS_POST_TEST\n";
  cout << "  R8MAT_GIVENS_POST computes a Givens postmultiplier rotation matrix.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) i4_power ( i + 1, j );
    }
  }

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 3;
  col = 2;

  cout << "\n";
  cout << "  I = " << row << "  J = " << col << "\n";

  g = r8mat_givens_post ( N, a, row, col );

  r8mat_print ( N, N, g, "  G" );

  ag = r8mat_mm_new ( N, N, N, a, g );

  r8mat_print ( N, N, ag, "  A*G" );

  delete [] ag;
  delete [] g;

  return;
# undef N
}
//****************************************************************************80

void r8mat_givens_pre_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_GIVENS_PRE_TEST tests R8MAT_GIVENS_PRE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N*N];
  int col;
  double *g;
  double *ga;
  int i;
  int j;
  int row;

  cout << "\n";
  cout << "R8MAT_GIVENS_PRE_TEST\n";
  cout << "  R8MAT_GIVENS_PRE computes a Givens premultiplier rotation matrix.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) i4_power ( i + 1, j );
    }
  }

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 3;
  col = 2;

  cout << "\n";
  cout << "  I = " << row << "  J = " << col << "\n";

  g = r8mat_givens_pre ( N, a, row, col );

  r8mat_print ( N, N, g, "  G" );

  ga = r8mat_mm_new ( N, N, N, g, a );

  r8mat_print ( N, N, ga, "  G*A" );

  delete [] g;
  delete [] ga;

  return;
# undef N
}
//****************************************************************************80

void r8mat_hess_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HESS_TEST tests R8MAT_HESS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double *h;
  double x[N] = { 1.0, 2.0, 3.0 };

  cout << "\n";
  cout << "R8MAT_HESS_TEST\n";
  cout << "  R8MAT_HESS estimates the Hessian matrix\n";
  cout << "  of a scalar function.\n";

  h = r8mat_hess ( r8mat_hess_f, N, x );

  r8mat_print ( N, N, h, "  Estimated jacobian:" );

  delete [] h;

  h = r8mat_hess_exact ( N, x );

  r8mat_print ( N, N, h, "  Exact jacobian:" );

  delete [] h;

  return;
# undef N
}
//****************************************************************************80

double r8mat_hess_f ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HESS_F is a sample nonlinear function for treatment by R8MAT_HESS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of parameters.
//
//    Input, double X[N], the parameter values.
//
//    Output, double R8MAT_HESS_F, the function value.
//
{
  double f;

  f = x[0] * x[0] + x[0] * x[1] + x[1] * cos ( 10.0 * x[2] );

  return f;
}
//****************************************************************************80

double *r8mat_hess_exact ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HESS_EXACT is the exact Hessian of R8MAT_HESS_F.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of parameters.
//
//    Input, double X[N], the parameter values.
//
//    Output, double R8MAT_HESS_EXACT[N*N], the Hessian values.
//
{
  double *h;

  h = new double[n*n];

  h[0+0*3] = 2.0;
  h[0+1*3] = 1.0;
  h[0+2*3] = 0.0;

  h[1+0*3] = 1.0;
  h[1+1*3] = 0.0;
  h[1+2*3] = -10.0 * sin ( 10.0 * x[2] );

  h[2+0*3] = 0.0;
  h[2+1*3] = -10.0 * sin ( 10.0 * x[2] );
  h[2+2*3] = -100.0 * x[1] * cos ( 10.0 * x[2] );

  return h;
}
//****************************************************************************80

void r8mat_house_axh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HOUSE_AXH_TEST tests R8MAT_HOUSE_AXH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_col;
  double *ah;
  double *h;
  double *ha;
  int k;
  int n = 5;
  double r8_hi;
  double r8_lo;
  int seed;
  double *v;

  cout << "\n";
  cout << "R8MAT_HOUSE_AXH_TEST\n";
  cout << "  R8MAT_HOUSE_AXH multiplies a matrix A times a\n";
  cout << "  compact Householder matrix.\n";

  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;

  a = r8mat_uniform_ab_new ( n, n, r8_lo, r8_hi, seed );

  r8mat_print ( n, n, a, "  Matrix A:" );
//
//  Request V, the compact form of the Householder matrix H
//  such that H*A packs column 3 of A.
//
  k = 3;
  a_col = ( a + ( k - 1 ) * n );
  v = r8vec_house_column ( n, a_col, k );

  r8vec_print ( n, v, "  Compact vector V so H*A packs column 3:" );

  h = r8mat_house_form ( n, v );

  r8mat_print ( n, n, h, "  Householder matrix H:" );
//
//  Compute A*H.
//
  ah = r8mat_house_axh_new ( n, a, v );

  r8mat_print ( n, n, ah, "  Indirect product A*H:" );

  delete [] ah;
//
//  Compare with a direct calculation.
//
  ah = r8mat_mm_new ( n, n, n, a, h );

  r8mat_print ( n, n, ah, "  Direct product A*H:" );
//
//  Verify that H*A packs column 3:
//
  ha = r8mat_mm_new ( n, n, n, h, a );

  r8mat_print ( n, n, ha, "  H*A should have packed column 3:" );
//
//  Free memory.
//
  delete [] a;
  delete [] ah;
  delete [] h;
  delete [] ha;
  delete [] v;

  return;
# undef N
}
//****************************************************************************80

void r8mat_house_form_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HOUSE_FORM_TEST tests R8MAT_HOUSE_FORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *h;
  double v[N] = { 0.0, 0.0, 1.0, 2.0, 3.0 };

  cout << "\n";
  cout << "R8MAT_HOUSE_FORM_TEST\n";
  cout << "  R8MAT_HOUSE_FORM forms a Householder\n";
  cout << "  matrix from its compact form.\n";

  r8vec_print ( N, v, "  Compact vector form V:" ) ;

  h = r8mat_house_form ( N, v );

  r8mat_print ( N, N, h, "  Householder matrix H:" );

  delete [] h;

  return;
# undef N
}
//****************************************************************************80

void r8mat_house_post_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HOUSE_POST_TEST tests R8MAT_HOUSE_POST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *ah;
  double b = 0.0;
  double c = 5.0;
  int col;
  double *h;
  int n = 5;
  int row;
  int seed;

  cout << "\n";
  cout << "R8MAT_HOUSE_POST_TEST\n";
  cout << "  R8MAT_HOUSE_POST computes a Householder postmultiplier;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( n, n, b, c, seed );

  r8mat_print ( n, n, a, "  Matrix A:" );

  row = 2;
  col = 3;

  cout << "\n";
  cout << "  I = " << row << "  J = " << col << "\n";

  h = r8mat_house_post ( n, a, row, col );

  r8mat_print ( n, n, h, "  Householder matrix H:" );

  ah = r8mat_mm_new ( n, n, n, a, h );

  r8mat_print ( n, n, ah, "  Product A*H:" );

  delete [] a;
  delete [] ah;
  delete [] h;

  return;
}
//****************************************************************************80

void r8mat_house_pre_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HOUSE_PRE_TEST tests R8MAT_HOUSE_PRE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double b = 0.0;
  double c = 5.0;
  int col;
  double *h;
  double *ha;
  int row;
  int seed;

  cout << "\n";
  cout << "R8MAT_HOUSE_PRE_TEST\n";
  cout << "  R8MAT_HOUSE_PRE computes a Householder premultiplier;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( N, N, b, c, seed );

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 2;
  col = 3;

  cout << "\n";
  cout << "  I = " << row << "  J = " << col << "\n";

  h = r8mat_house_pre ( N, a, row, col );

  r8mat_print ( N, N, h, "  Householder matrix H:" );

  ha = r8mat_mm_new ( N, N, N, h, a );

  r8mat_print ( N, N, ha, "  Product H*A:" );

  delete [] a;
  delete [] h;
  delete [] ha;

  return;
# undef N
}
//****************************************************************************80

void r8mat_indicator_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_INDICATOR_NEW_TEST tests R8MAT_INDICATOR_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int m = 5;
  int n = 4;

  cout << "\n";
  cout << "R8MAT_INDICATOR_NEW_TEST\n";
  cout << "  R8MAT_INDICATOR_NEW returns an indicator matrix.\n";

  a = r8mat_indicator_new ( m, n );

  r8mat_print ( m, n, a, "  The indicator matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_inverse_2d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_INVERSE_2D_TEST tests R8MAT_INVERSE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "R8MAT_INVERSE_2D_TEST\n";
  cout << "  R8MAT_INVERSE_2D inverts a 2 by 2 matrix.\n";

  a[0+0*N] = 1.0;
  a[0+1*N] = 2.0;

  a[1+0*N] = 3.0;
  a[1+1*N] = 4.0;

  r8mat_print ( 2, 2, a, "  Matrix A:" );

  b = r8mat_inverse_2d ( a );

  r8mat_print ( 2, 2, b, "  Inverse matrix B:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 2, 2, c, "  C = A * B:" );

  delete [] b;

  return;

# undef N
}
//****************************************************************************80

void r8mat_inverse_3d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_INVERSE_3D_TEST tests R8MAT_INVERSE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "R8MAT_INVERSE_3D_TEST\n";
  cout << "  R8MAT_INVERSE_3D inverts a 3 by 3 matrix.\n";

  a[0+0*N] = 3.0;
  a[0+1*N] = 2.0;
  a[0+2*N] = 1.0;

  a[1+0*N] = 2.0;
  a[1+1*N] = 2.0;
  a[1+2*N] = 1.0;

  a[2+0*N] = 0.0;
  a[2+1*N] = 1.0;
  a[2+2*N] = 1.0;

  r8mat_print ( 3, 3, a, "  Matrix A:" );

  b = r8mat_inverse_3d ( a );

  r8mat_print ( 3, 3, b, "  Inverse matrix B:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 3, 3, c, "  C = A * B:" );

  delete [] b;

  return;

# undef N
}
//****************************************************************************80

void r8mat_inverse_4d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_INVERSE_4D_TEST tests R8MAT_INVERSE_4D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "R8MAT_INVERSE_4D_TEST\n";
  cout << "  R8MAT_INVERSE_4D inverts a 4 x 4 matrix.\n";


  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i <= j )
      {
        a[i+j*N] = ( double ) ( N - j );
      }
      else if ( j == i - 1 )
      {
        a[i+j*N] = ( double ) ( N - j - 1 );
      }
      else {
        a[i+j*N] = 0.0;
      }
    }
  }

  r8mat_print ( 4, 4, a, "  Matrix A:" );

  b = r8mat_inverse_4d ( a );

  r8mat_print ( 4, 4, b, "  Inverse matrix B:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 4, 4, c, "  C = A * B:" );

  delete [] b;

  return;

# undef N
}
//****************************************************************************80

void r8mat_jac_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_JAC_TEST tests R8MAT_JAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double eps = 0.00001;
  double *fprime;
  int m = 3;
  double x[N] = { 1.0, 2.0, 3.0, 4.0 };

  cout << "\n";
  cout << "R8MAT_JAC_TEST\n";
  cout << "  R8MAT_JAC estimates the M by N jacobian matrix\n";
  cout << "  of a nonlinear function.\n";

  fprime = r8mat_jac ( m, N, eps, r8mat_jac_f, x );

  r8mat_print ( m, N, fprime, "  Estimated jacobian:" );

  delete [] fprime;

  fprime = r8mat_jac_exact ( m, N, x );

  r8mat_print (  m, N, fprime, "  Exact jacobian:" );

  delete [] fprime;

  return;
# undef N
}
//****************************************************************************80

double *r8mat_jac_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_JAC_F is a sample nonlinear function for treatment by R8MAT_JAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of functions.
//
//    Input, int N, the number of parameters.
//
//    Input, double X[N], the parameter values.
//
//    Output, double R8MAT_JAC_F[M], the function values.
//
{
  double *f;

  f = new double[m];

  f[0] = sin ( x[0] * x[1] );
  f[1] = sqrt ( 1.0 + x[0] * x[0] ) + x[2];
  f[2] = x[0] + 2.0 * x[1] + 3.0 * x[2] + 4.0 * x[3];

  return f;
}
//****************************************************************************80

double *r8mat_jac_exact ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_JAC_EXACT is the exact jacobian of R8MAT_JAC_F.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of functions.
//
//    Input, int N, the number of parameters.
//
//    Input, double X[N], the parameter values.
//
//    Output, double R8MAT_JAC_EXACT[M*N], the jacobian values.
//
{
  double *fprime;

  fprime = new double[m*n];

  fprime[0+0*3] = cos ( x[0] * x[1] ) * x[1];
  fprime[0+1*3] = cos ( x[0] * x[1] ) * x[0];
  fprime[0+2*3] = 0.0;
  fprime[0+3*3] = 0.0;

  fprime[1+0*3] = x[0] / sqrt ( 1.0 + x[0] * x[0] );
  fprime[1+1*3] = 0.0;
  fprime[1+2*3] = 1.0;
  fprime[1+3*3] = 0.0;

  fprime[2+0*3] = 1.0;
  fprime[2+1*3] = 2.0;
  fprime[2+2*3] = 3.0;
  fprime[2+3*3] = 4.0;

  return fprime;
}
//****************************************************************************80

void r8mat_kronecker_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_KRONECKER_TEST tests R8MAT_KRONECKER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m1 = 2;
  int m2 = 3;
  int m = m1 * m2;
  int n1 = 3;
  int n2 = 2;
  int n = n1 * n2;

  double a[2*3] = {
    1.0, 4.0, 
    2.0, 5.0, 
    3.0, 6.0 };
  double b[3*2] = {
    7.0,  9.0, 11.0, 
    8.0, 10.0, 12.0 };
  double *c;

  cout << "\n";
  cout << "R8MAT_KRONECKER_TEST\n";
  cout << "  R8MAT_KRONECKER computes the Kronecker product\n";
  cout << "  of two R8MAT's.\n";

  r8mat_print ( m1, n1, a, "  Factor matrix A:" );
  r8mat_print ( m2, n2, b, "  Factor matrix B:" );

  c = r8mat_kronecker ( m1, n1, a, m2, n2, b );

  r8mat_print ( m, n, c, "  Kronecker product C = kron(A,B)" );

  delete [] c;

  return;
}
//****************************************************************************80

void r8mat_l_inverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_L_INVERSE_TEST tests R8MAT_L_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N*N] = {
    1.0, 2.0, 4.0,  7.0,
    0.0, 3.0, 5.0,  8.0,
    0.0, 0.0, 6.0,  9.0,
    0.0, 0.0, 0.0, 10.0 };
  double *b;
  double *c;

  cout << "\n";
  cout << "R8MAT_L_INVERSE_TEST\n";
  cout << "  R8MAT_L_INVERSE inverts a lower triangular R8MAT.\n";

  r8mat_print ( N, N, a, "  Matrix A to be inverted:" );

  b = r8mat_l_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void r8mat_l_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_L_PRINT_TEST tests R8MAT_L_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double a1[28] = {
    11.0, 21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          22.0, 32.0, 42.0, 52.0, 62.0, 72.0,
                33.0, 43.0, 53.0, 63.0, 73.0,
                      44.0, 54.0, 64.0, 74.0,
                            55.0, 65.0, 75.0,
                                  66.0, 76.0,
                                        77.0 };
  double a2[18] = {
    11.0, 21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          22.0, 32.0, 42.0, 52.0, 62.0, 72.0,
                33.0, 43.0, 53.0, 63.0, 73.0 };
  double a3[10] = {
    11.0, 21.0, 31.0, 41.0,
          22.0, 32.0, 42.0,
                33.0, 43.0,
                      44.0 };
  int m1 = 7;
  int m2 = 7;
  int m3 = 4;
  int n1 = 7;
  int n2 = 3;
  int n3 = 7;

  cout << "\n";
  cout << "R8MAT_L_PRINT_TEST\n";
  cout << "  R8MAT_L_PRINT prints a lower triangular R8MAT\n";
  cout << "  stored compactly.  Only the (possibly) nonzero\n";
  cout << "  elements are printed.\n";

  r8mat_l_print ( m1, n1, a1, "  A 7 by 7 matrix." );

  r8mat_l_print ( m2, n2, a2, "  A 7 by 3 matrix." );

  r8mat_l_print ( m3, n3, a3, "  A 4 by 7 matrix." );

  return;
}
//****************************************************************************80

void r8mat_l1_inverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_L1_INVERSE_TEST tests R8MAT_L1_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 6
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N*N] = {
     1.0, 2.0, 0.0, 5.0, 0.0, 75.0,
     0.0, 1.0, 0.0, 0.0, 0.0,  0.0,
     0.0, 0.0, 1.0, 3.0, 0.0,  0.0,
     0.0, 0.0, 0.0, 1.0, 0.0,  6.0,
     0.0, 0.0, 0.0, 0.0, 1.0,  4.0,
     0.0, 0.0, 0.0, 0.0, 0.0,  1.0 };
  double *b;
  double *c;

  cout << "\n";
  cout << "R8MAT_L1_INVERSE_TEST\n";
  cout << "  R8MAT_L1_INVERSE inverts a unit lower triangular R8MAT.\n";

  r8mat_print ( N, N, a, "  Matrix A to be inverted:" );

  b = r8mat_l1_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void r8mat_lu_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_LU_TEST tests R8MAT_LU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5

  double *a;
  double l[M*M];
  double *lu;
  double p[M*M];
  double *plu;
  double u[M*N];
  double x[N] = { 1.0, 10.0, 4.0, 2.0, 3.0 };

  cout << "\n";
  cout << "R8MAT_LU_TEST\n";
  cout << "  R8MAT_LU computes the LU factors of an R8MAT.\n";

  a = r8mat_vand2 ( N, x );

  r8mat_print ( M, N, a, "  Matrix to be factored:" );

  r8mat_lu ( M, N, a, l, p, u );

  r8mat_print ( M, M, p, "  P factor:" );

  r8mat_print ( M, M, l, "  L factor:" );

  r8mat_print ( M, N, u, "  U factor:" );

  lu = r8mat_mm_new ( M, M, N, l, u );

  plu = r8mat_mm_new ( M, M, N, p, lu );

  r8mat_print ( M, N, plu, "  P*L*U:" );

  delete [] a;
  delete [] lu;
  delete [] plu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_max_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MAX_TEST tests R8MAT_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp;

  cout << "\n";
  cout << "R8MAT_MAX_TEST\n";
  cout << "  R8MAT_MAX computes the maximum value of an R8MAT.\n";
 
  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, seed );

  r8mat_print ( m, n, a, "  Random array:" );

  temp = r8mat_max ( m, n, a );

  cout << "\n";
  cout << "  Maximum value = " << temp << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_max_index_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MAX_INDEX_TEST tests R8MAT_MAX_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 3

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "R8MAT_MAX_INDEX_TEST\n";
  cout << "  R8MAT_MAX_INDEX locates the maximum entry of an R8MAT;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  Random array:" );

  r8mat_max_index ( M, N, a, i, j );

  cout << "\n";
  cout << "  Maximum I,J indices            " << i << "  " << j << "\n";

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_maxcol_minrow_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MAXCOL_MINROW_TEST tests R8MAT_MAXCOL_MINROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;

  cout << "\n";
  cout << "R8MAT_MAXCOL_MINROW_TEST\n";
  cout << "  R8MAT_MAXCOL_MINROW computes the maximum over\n";
  cout << "  columns of the mininum over rows;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, seed );

  r8mat_print ( m, n, a, "  Random array:" );

  cout << "  MAXCOL_MINROW = " << r8mat_maxcol_minrow ( m, n, a ) << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_maxrow_mincol_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MAXROW_MINCOL_TEST tests R8MAT_MAXROW_MINCOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  cout << "\n";
  cout << "R8MAT_MAXROW_MINCOL_TEST\n";
  cout << "  R8MAT_MAXROW_MINCOL computes the maximum over\n";
  cout << "  rows of the mininum over columns;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, seed );

  r8mat_print ( m, n, a, "  Random array:" );

  cout << "  MAXROW_MINCOL = " << r8mat_maxrow_mincol ( m, n, a ) << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_min_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MIN_TEST tests R8MAT_MIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp;

  cout << "\n";
  cout << "R8MAT_MIN_TESTn";
  cout << "  For a real matrix,\n";
  cout << "  R8MAT_MIN computes the minimum value;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, seed );

  r8mat_print ( m, n, a, "  Random array:" );

  temp = r8mat_min ( m, n, a );

  cout << "\n";
  cout << "  Minimum value = " << temp << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_min_index_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MIN_INDEX_TEST tests R8MAT_MIN_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 3

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "R8MAT_MIN_INDEX_TEST\n";
  cout << "  R8MAT_MIN_INDEX locates the minimum entry of an R8MAT;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  Random array:" );

  cout << "\n";
  r8mat_min_index ( M, N, a, i, j );
  cout << "  Minimum I,J indices            " << i << "  " << j << "\n";

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_mincol_maxrow_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MINCOL_MAXROW_TEST tests R8MAT_MINCOL_MAXROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  cout << "\n";
  cout << "R8MAT_MINCOL_MAXROW_TEST\n";
  cout << "  R8MAT_MINCOL_MAXROW computes the minimum over\n";
  cout << "  columns of the maxinum over rows;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, seed );

  r8mat_print ( m, n, a, "  Random array:" );

  cout << "  MINCOL_MAXROW = " << r8mat_mincol_maxrow ( m, n, a ) << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_minrow_maxcol_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MINROW_MAXCOL_TEST tests R8MAT_MINROW_MAXCOL_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  cout << "\n";
  cout << "R8MAT_MINROW_MAXCOL_TEST\n";
  cout << "  R8MAT_MINROW_MAXCOL computes the minimum over\n";
  cout << "  rows of the maxinum over columns;\n";
  cout << "\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, seed );

  r8mat_print ( m, n, a, "  Random array:" );

  cout << "  MINROW_MAXCOL = " << r8mat_minrow_maxcol ( m, n, a ) << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_mm_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM_NEW_TEST tests R8MAT_MM_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
# define N3 4
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2*N3] = {
     1.0,  2.0,  3.0,
     4.0,  5.0,  6.0,
     0.0,  0.0,  1.0,
    -1.0,  2.0, -1.0 };
  double *c;

  cout << "\n";
  cout << "R8MAT_MM_NEW_TEST\n";
  cout << "  R8MAT_MM_NEW multiplies two (rectangular) matrices\n";
  cout << "  and returns the result as the function value.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8mat_print ( N2, N3, b, "  Matrix B:" );

  c = r8mat_mm_new ( N1, N2, N3, a, b );

  r8mat_print ( N1, N3, c, "  Product C = A * B:" );

  delete [] c;

  return;
# undef N1
# undef N2
# undef N3
}
//****************************************************************************80

void r8mat_mm_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM_TEST tests R8MAT_MM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
# define N3 4
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2*N3] = {
     1.0,  2.0,  3.0,
     4.0,  5.0,  6.0,
     0.0,  0.0,  1.0,
    -1.0,  2.0, -1.0 };
  double c[N1*N3];

  cout << "\n";
  cout << "R8MAT_MM_TEST\n";
  cout << "  R8MAT_MM multiplies two (rectangular) matrices\n";
  cout << "  and returns the result as an argument.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8mat_print ( N2, N3, b, "  Matrix B:" );

  r8mat_mm ( N1, N2, N3, a, b, c );

  r8mat_print ( N1, N3, c, "  Product C = A * B:" );

  return;
# undef N1
# undef N2
# undef N3
}
//****************************************************************************80

void r8mat_mv_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_NEW_TEST tests R8MAT_MV_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2] = {
     1.0,  2.0,  3.0 };
  double *c;

  cout << "\n";
  cout << "R8MAT_MV_NEW_TEST\n";
  cout << "  R8MAT_MV_NEW multiplies a (rectangular) matrix times a vector,\n";
  cout << "  and returns the result as the function value.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N2, b, "  Vector B:" );

  c = r8mat_mv_new ( N1, N2, a, b );

  r8vec_print ( N1, c, "  Product C = A * B:" );

  delete [] c;

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void r8mat_mv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_TEST tests R8MAT_MV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2] = {
     1.0,  2.0,  3.0 };
  double c[N1];

  cout << "\n";
  cout << "R8MAT_MV_TEST\n";
  cout << "  R8MAT_MV multiplies a (rectangular) matrix times a vector,\n";
  cout << "  and returns the result as an argument.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N2, b, "  Vector B:" );

  r8mat_mv ( N1, N2, a, b, c );

  r8vec_print ( N1, c, "  Product C = A * B:" );

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void r8mat_mtv_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MTV_NEW_TEST tests R8MAT_MTV_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N1] = {
     1.0,  2.0 };
  double *c;

  cout << "\n";
  cout << "R8MAT_MTV_NEW_TEST\n";
  cout << "  R8MAT_MTV_NEW multiplies a transposed matrix times a vector,\n";
  cout << "  and returns the result as the function value.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N1, b, "  Vector B:" );

  c = r8mat_mtv_new ( N1, N2, a, b );

  r8vec_print ( N2, c, "  Product C = A' * B:" );

  delete [] c;

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void r8mat_mtv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MTV_TEST tests R8MAT_MTV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N1] = {
     1.0,  2.0 };
  double c[N2];

  cout << "\n";
  cout << "R8MAT_MTV_TEST\n";
  cout << "  R8MAT_MTV multiplies a transposed matrix times a vector,\n";
  cout << "  and returns the result as an argument.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N1, b, "  Vector B:" );

  r8mat_mtv ( N1, N2, a, b, c );

  r8vec_print ( N2, c, "  Product C = A' * B:" );

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void r8mat_nint_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NINT_TEST tests R8MAT_NINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int m;
  int n;
  int seed;
  double x1;
  double x2;

  cout << "\n";
  cout << "R8MAT_NINT_TEST\n";
  cout << "  R8MAT_NINT rounds an R8MAT.\n";

  m = 5;
  n = 4;
  x1 = -5.0;
  x2 = +5.0;
  seed = 123456789;
  a = r8mat_uniform_ab_new ( m, n, x1, x2, seed );
  r8mat_print ( m, n, a, "  Matrix A:" );
  r8mat_nint ( m, n, a );
  r8mat_print ( m, n, a, "  Rounded matrix A:" );

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_nonzeros_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NONZEROS_TEST tests R8MAT_NONZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int j;
  int m = 5;
  int n = 4;
  int c1;
  int c2;

  cout << "\n";
  cout << "R8MAT_NONZEROS_TEST\n";
  cout << "  R8MAT_NONZEROS counts nonzeros in an R8MAT.\n";

  a = new double[m*n];

  c1 = 0;
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( ( i % 2 ) == 0 && ( j % 2 ) == 0 )
      {
        a[i+j*m] = 1;
        c1 = c1 + 1;
      }
      else
      {
        a[i+j*m] = 0;
      }
    }
  }

  r8mat_print ( m, n, a, "  Matrix A:" );

  c2 = r8mat_nonzeros ( m, n, a );

  cout << "\n";
  cout << "  Expected nonzeros = " << c1 << "\n";
  cout << "  Computed nonzeros = " << c2 << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_norm_fro_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORM_FRO_TEST tests R8MAT_NORM_FRO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int j;
  int k;
  int m = 5;
  int n = 4;
  double t1;
  double t2;

  cout << "\n";
  cout << "R8MAT_NORM_FRO_TEST\n";
  cout << "  R8MAT_NORM_FRO computes the Frobenius norm of a matrix.\n";

  a = new double[m*n];

  k = 0;
  t1 = 0.0;
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      k = k + 1;
      a[i+j*m] = ( double ) ( k );
      t1 = t1 + k * k;
    }
  }
  t1 = sqrt ( t1 );

  r8mat_print ( m, n, a, "  Matrix A:" );

  t2 = r8mat_norm_fro ( m, n, a );

  cout << "\n";
  cout << "  Expected Frobenius norm = " << t1 << "\n";
  cout << "  Computed Frobenius norm = " << t2 << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_norm_l1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORM_L1_TEST tests R8MAT_NORM_L1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int m;
  int n;
  int seed;
  double t;
  double x1;
  double x2;

  cout << "\n";
  cout << "R8MAT_NORM_L1_TEST\n";
  cout << "  R8MAT_NORM_L1 computes the L1 norm of a matrix.\n";

  m = 5;
  n = 4;
  x1 = -5.0;
  x2 = +5.0;
  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, x1, x2, seed );
  r8mat_nint ( m, n, a );

  r8mat_print ( m, n, a, "  Matrix A:" );

  t = r8mat_norm_l1 ( m, n, a );

  cout << "\n";
  cout << "  L1 norm = " << t << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8mat_nullspace_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NULLSPACE_TEST tests R8MAT_NULLSPACE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  double *ax;
  int m = M;
  int n = N;
  double *nullspace;
  int nullspace_size;

  cout << "\n";
  cout << "R8MAT_NULLSPACE_TEST\n";
  cout << "  R8MAT_NULLSPACE computes the nullspace of a matrix.\n";

  r8mat_print ( m, n, a, "  Input A:" );

  nullspace_size = r8mat_nullspace_size ( m, n, a );

  cout << "\n";
  cout << "  Nullspace size is " << nullspace_size << "\n";

  nullspace = r8mat_nullspace ( m, n, a, nullspace_size );

  r8mat_print ( n, nullspace_size, nullspace, "  Nullspace vectors:" );

  ax = r8mat_mm_new ( m, n, nullspace_size, a, nullspace );

  r8mat_print ( m, nullspace_size, ax, "  Product A * Nullspace vectors:" );

  delete [] ax;
  delete [] nullspace;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_nullspace_size_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NULLSPACE_SIZE_TEST tests R8MAT_NULLSPACE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  int m = M;
  int n = N;
  int nullspace_size;

  cout << "\n";
  cout << "R8MAT_NULLSPACE_SIZE_TEST\n";
  cout << "  R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.\n";

  r8mat_print ( m, n, a, "  Input A:" );

  nullspace_size = r8mat_nullspace_size ( m, n, a );

  cout << "\n";
  cout << "  Nullspace size is " << nullspace_size << "\n";

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_orth_uniform_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ORTH_UNIFORM_NEW_TEST tests R8MAT_ORTH_UNIFORM_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *at;
  double *ata;
  int n = 5;
  int seed;

  cout << "\n";
  cout << "R8MAT_ORTH_UNIFORM_NEW_TEST\n";
  cout << "  R8MAT_ORTH_UNIFORM_NEW computes a random orthogonal matrix.\n";

  seed = 123456789;

  a = r8mat_orth_uniform_new ( n, seed );

  r8mat_print ( n, n, a, "  Random orthogonal matrix A" );

  at = r8mat_transpose_new ( n, n, a );

  ata = r8mat_mm_new ( n, n, n, at, a );

  r8mat_print ( n, n, ata, "  AT*A" );

  delete [] a;
  delete [] at;
  delete [] ata;

  return;
}
//****************************************************************************80

void r8mat_plot_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PLOT_TEST tests R8MAT_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define N 100

  double a[M*N];
  int i;
  int im1;
  int ip1;

  r8mat_zero ( M, N, a );

  for ( i = 0; i < M; i++ )
  {
    a[i+i*M] = -2.0;

    if ( i+1 < N )
    {
      ip1 = i+1;
    }
    else
    {
      ip1 = 0;
    }

    a[i+ip1*M] = 1.0;

    if ( 0 <= i-1 )
    {
      im1 = i-1;
    }
    else
    {
      im1 = N-1;
    }
    a[i+im1*M] = 1.0;
  }

  cout << "\n";
  cout << "R8MAT_PLOT_TEST\n";
  cout << "  R8MAT_PLOT prints a symbolic picture of a matrix.\n";
  cout << "  Typically,\n";
  cout << "\n";
  cout << "    - for negative, \n";
  cout << "    0 for zero, and\n";
  cout << "    + for positive entries\n";
  cout << "\n";
  cout << "  or\n";
  cout << "\n";
  cout << "    X for nonzero and\n";
  cout << "    0 for zero.\n";
  cout << "\n";

  r8mat_plot ( M, N, a, "  A plot of the matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_power_method_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_POWER_METHOD_TEST tests R8MAT_POWER_METHOD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N*N];
  double *av;
  int i;
  int j;
  double r;
  double v[N];

  cout << "\n";
  cout << "R8MAT_POWER_METHOD_TEST\n";
  cout << "  R8MAT_POWER_METHOD applies the power method\n";
  cout << "  to a matrix.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( j == i - 1 || j == i + 1 )
      {
        a[i+j*N] = -1.0;
      }
      else if ( j == i )
      {
        a[i+j*N] = 2.0;
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }
  r8vec_zero ( N, v );

  r8mat_power_method ( N, a, &r, v );

  cout << "\n";
  cout << "  Estimated eigenvalue = " << r << "\n";

  r8vec_print ( N, v, "  Estimated eigenvector V:" );

  av = r8mat_mv_new ( N, N, a, v );

  r8vec_print ( N, av, "  Value of A*V:" );

  delete [] av;

  return;
# undef N
}
//****************************************************************************80

void r8mat_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_TEST tests R8MAT_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8MAT_PRINT_TEST\n";
  cout << "  R8MAT_PRINT prints an R8MAT.\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print ( m, n, a, "  The R8MAT:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_print_some_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME_TEST tests R8MAT_PRINT_SOME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8MAT_PRINT_SOME_TEST\n";
  cout << "  R8MAT_PRINT_SOME prints some of an R8MAT.\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print_some ( m, n, a, 2, 1, 4, 2, "  The R8MAT, rows 2:4, cols 1:2:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_ref_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_REF_TEST tests R8MAT_REF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8MAT_REF_TEST\n";
  cout << "  R8MAT_REF computes the row echelon form of a matrix.\n";

  r8mat_print ( m, n, a, "  Input A:" );

  r8mat_ref ( m, n, a );

  r8mat_print ( m, n, a, "  REF form:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_rref_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_RREF_TEST tests R8MAT_RREF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8MAT_RREF_TEST\n";
  cout << "  R8MAT_RREF computes the reduced row echelon form of a matrix.\n";

  r8mat_print ( m, n, a, "  Input A:" );

  r8mat_rref ( m, n, a );

  r8mat_print ( m, n, a, "  REF form:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_solve_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE_TEST tests R8MAT_SOLVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define RHS_NUM 2

  double a[N*(N+RHS_NUM)] = {
     1.0,  4.0,  7.0,
     2.0,  5.0,  8.0,
     3.0,  6.0,  0.0,
    14.0, 32.0, 23.0,
     7.0, 16.0,  7.0 };
  int i;
  int info;
  int j;

  cout << "\n";
  cout << "R8MAT_SOLVE_TEST\n";
  cout << "  R8MAT_SOLVE solves linear systems.\n";
//
//  Print out the matrix to be inverted.
//
  r8mat_print ( N, N+RHS_NUM, a, "  The linear system:" );
//
//  Solve the systems.
//
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  The input matrix was singular.\n";
    cout << "  The solutions could not be computed.\n";
    return;
  }

  cout << "\n";
  cout << "  The computed solutions:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    for ( j = N; j < N+RHS_NUM; j++ )
    {
      cout << setw(10) << a[i+j*N] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef RHS_NUM
}
//****************************************************************************80

void r8mat_solve_2d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE_2D_TEST tests R8MAT_SOLVE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double det;
  int i;
  int n = 2;
  int seed;
  int test;
  int test_num = 5;
  double *x;
  double *x2;

  cout << "\n";
  cout << "R8MAT_SOLVE_2D_TEST\n";
  cout << "  R8MAT_SOLVE_2D solves 2D linear systems.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8mat_uniform_01_new ( n, n, seed );
    x = r8vec_uniform_01_new ( n, seed );
    b = r8mat_mv_new ( n, n, a, x );

    x2 = r8mat_solve_2d ( a, b, &det );

    cout << "\n";
    cout << "  Solution / Computed:\n";
    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(14) << x[i]
           << "  " << setw(14) << x2[i] << "\n";
    }

    delete [] a;
    delete [] b;
    delete [] x;
    delete [] x2;
  }

  return;
}
//****************************************************************************80

void r8mat_solve_3d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE_3D_TEST tests R8MAT_SOLVE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double det;
  int i;
  int n = 3;
  int seed;
  int test;
  int test_num = 5;
  double *x;
  double *x2;

  cout << "\n";
  cout << "R8MAT_SOLVE_3D_TEST\n";
  cout << "  R8MAT_SOLVE_3D solves 3D linear systems.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8mat_uniform_01_new ( n, n, seed );
    x = r8vec_uniform_01_new ( n, seed );
    b = r8mat_mv_new ( n, n, a, x );

    x2 = r8mat_solve_3d ( a, b, &det );

    cout << "\n";
    cout << "  Solution / Computed:\n";
    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(14) << x[i]
           << "  " << setw(14) << x2[i] << "\n";
    }

    delete [] a;
    delete [] b;
    delete [] x;
    delete [] x2;
  }

  return;
}
//****************************************************************************80

void r8mat_solve2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE2_TEST tests R8MAT_SOLVE2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 February 2014
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double *a;
  double a1[2*2] = {
    1.0, 3.0,
    2.0, 4.0 };
  double a2[3*3] = {
    2.0, 1.0, 1.0,
    1.0, 1.0, 0.0,
    1.0, 0.0, 1.0 };
  double a3[4*4] = {
    1.0, 2.0, 1.0, 3.0,
    0.0, 1.0, 2.0, 1.0,
    0.0, 0.0, 3.0, 2.0,
    1.0, 3.0, 0.0, 1.0 };
  double a4[3*3] = {
    2.0, 1.0, 3.0,
    4.0, 2.0, 6.0,
    1.0, 4.0, 5.0 };
  double *b;
  double b1[2] = { 5.0, 11.0 };
  double b2[3] = { 4.0, 2.0, 2.0 };
  double b3[4] = { 5.0, 11.0, 16.0, 15.0 };
  double b4[3] = { 13.0, 17.0, 20.0 };
  int ierror;
  int n;
  int n_test[TEST_NUM] = { 2, 3, 4, 3 };
  int test;
  double *x;

  cout << "\n";
  cout << "R8MAT_SOLVE2_TEST\n";
  cout << "  R8MAT_SOLVE2 is a linear solver.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = n_test[test];

    if ( test == 0 )
    {
      a = a1;
      b = b1;
    }
    else if ( test == 1 )
    {
      a = a2;
      b = b2;
    }
    else if ( test == 2 )
    {
      a = a3;
      b = b3;
    }
    else if ( test == 3 )
    {
      a = a4;
      b = b4;
    }

    r8vec_print ( n, b, "  Right hand side:" );

    x = r8mat_solve2 ( n, a, b, ierror );

    cout << "\n";
    if ( ierror == 0 )
    {
      cout << "  The system is nonsingular.\n";
    }
    else if ( ierror == 1 )
    {
      cout << "  The system is singular, but consistent.\n";
    }
    else if ( ierror == 2 )
    {
      cout << "  The system is singular and inconsistent.\n";
    }

    r8vec_print ( n, x, "  Computed solution:" );

    delete [] x;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8mat_sub_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SUB_NEW_TEST tests R8MAT_SUB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double *c;
  int m = 4;
  int n = 4;

  cout << "\n";
  cout << "R8MAT_SUB_NEW_TEST\n";
  cout << "  R8MAT_SUB_NEW computes C = A - B for R8MAT's\n";

  a = r8mat_indicator_new ( m, n );

  b = r8mat_transpose_new ( m, n, a );

  c = r8mat_sub_new ( m, n, a, b );

  r8mat_print ( m, n, a, "  A:" );
  r8mat_print ( m, n, b, "  B:" );
  r8mat_print ( m, n, c, "  C = A-B:" );

  delete [] a;
  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void r8mat_symm_jacobi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SYMM_JACOBI_TEST tests R8MAT_SYMM_JACOBI;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int n = 5;
  double *q;
  int seed;
  double *x;

  cout << "\n";
  cout << "R8MAT_SYMM_JACOBI_TEST\n";
  cout << "  For a symmetric R8MAT:\n";
  cout << "  R8MAT_SYMM_JACOBI diagonalizes;\n";
//
//  Choose eigenvalues.
//
  x = r8vec_indicator1_new ( n );
//
//  Choose eigenvectors.
//
  seed = 123456789;

  q = r8mat_orth_uniform_new ( n, seed );
//
//  Now get A = Q*X*Q.
//
  a = r8mat_symm_eigen ( n, x, q );

  r8mat_print ( n, n, a, "  Matrix to diagonalize:" );

  r8mat_symm_jacobi ( n, a );

  cout << "\n";
  cout << "  Computed Eigenvalues:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(12) << a[i+i*n] << "\n";
  }

  delete [] a;
  delete [] q;
  delete [] x;

  return;
}
//****************************************************************************80

void r8mat_to_r8cmat_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TO_R8CMAT_NEW_TEST tests R8MAT_TO_R8CMAT_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double **b;
  double *c;
  int m = 5;
  int n = 4;

  cout << "\n";
  cout << "R8MAT_TO_R8CMAT_NEW_TEST\n";
  cout << "  R8MAT_TO_R8CMAT_NEW converts an R8MAT to an R8CMAT.\n";
  cout << "\n";
  cout << "  Data is of order (" << m << "," << n << ".\n";
//
//  Set the R8MAT.
//
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  The R8MAT A:" );
//
//  Convert.
//
  b = r8mat_to_r8cmat_new ( m, n, a );
  r8cmat_print ( m, n, b, "  The R8CMAT B:" );
//
//  Recover the matrix.
//
  c = r8cmat_to_r8mat_new ( m, n, b );
  r8mat_print ( m, n, c, "  The R8MAT C:" );
//
//  Free memory.
//
  delete [] a;
  r8cmat_delete ( m, n, b );
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void r8mat_to_r8plu_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TO_R8PLU_TEST tests R8MAT_TO_R8PLU;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double a2[N*N];
  double b = 0.0;
  double c = 1.0;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "R8MAT_TO_R8PLU_TEST\n";
  cout << "  R8MAT_TO_R8PLU determines the compressed PLU factors\n";
  cout << "  of a real general matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_ab_new ( N, N, b, c, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
  }
//
//  Display the gory details.
//
  i4vec_print ( N, pivot, "  The pivot vector P:" );

  r8mat_print ( N, N, lu, "  The compressed LU factors:" );
//
//  Recover the matrix from the PLU factors.
//
  r8plu_to_r8mat ( N, pivot, lu, a2 );

  r8mat_print ( N, N, a2, "  The recovered matrix A2:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8mat_to_r8rmat_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TO_R8RMAT_TEST tests R8MAT_TO_R8RMAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double **b;
  double *c;
  int m = 5;
  int n = 4;

  cout << "\n";
  cout << "R8MAT_TO_R8RMAT_TEST\n";
  cout << "  R8MAT_TO_R8RMAT converts an R8MAT to an R8RMAT.\n";
  cout << "\n";
  cout << "  Data is of order (" << m << "," << n << ").\n";
//
//  Set the R8MAT.
//
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  The R8MAT A:" );
//
//  Convert.
//
  b = r8mat_to_r8rmat ( m, n, a );
  r8rmat_print ( m, n, b, "  The R8RMAT B:" );
//
//  Recover the matrix.
//
  c = r8rmat_to_r8mat ( m, n, b );
  r8mat_print ( m, n, c, "  The R8MAT C:" );
//
//  Free memory.
//
  free ( a );
  r8rmat_delete ( m, n, b );
  free ( c );

  return;
}
//****************************************************************************80

void r8mat_trace_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRACE_TEST tests R8MAT_TRACE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N];
  int i;
  int j;
  double trace;

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i <= j )
      {
        a[i+j*N] = ( double ) ( N - j );
      }
      else if ( j == i - 1 )
      {
        a[i+j*N] = ( double ) ( N - j - 1 );
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }

  cout << "\n";
  cout << "R8MAT_TRACE_TEST\n";
  cout << "  R8MAT_TRACE computes the trace of an R8MAT\n";

  r8mat_print ( N, N, a, "  Matrix:" );

  trace = r8mat_trace ( N, a );

  cout << "\n";
  cout << "  Trace is " << trace << "\n";

  return;
# undef N
}
//****************************************************************************80

void r8mat_transpose_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_NEW_TEST tests R8MAT_TRANSPOSE_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *at;
  int m;
  int n;

  cout << "\n";
  cout << "R8MAT_TRANSPOSE_NEW_TEST\n";
  cout << "  R8MAT_TRANSPOSE_NEW transposes an R8MAT.\n";

  m = 5;
  n = 4;
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  Matrix A:" );

  at = r8mat_transpose_new ( m, n, a );
  r8mat_print ( n, m, at, "  Transposed matrix At:" );

  delete [] a;
  delete [] at;

  return;
}
//****************************************************************************80

void r8mat_transpose_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_TEST tests R8MAT_TRANSPOSE_PRINT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 7
# define N 12

  double a[M*N];
  int i;
  int j;

  cout << "\n";
  cout << "R8MAT_TRANSPOSE_PRINT_TEST\n";
  cout << "  R8MAT_TRANSPOSE_PRINT prints an R8MAT,\n";
  cout << "  transposed.\n";
  cout << "\n";
  cout << "  Matrix row order M =    " << M << "\n";
  cout << "  Matrix column order N = " << N << "\n";
//
//  Set the matrix.
//
  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      a[i-1+(j-1)*M] = ( double ) ( i * 100 + j );
    }
  }

  r8mat_transpose_print ( M, N, a, "  The transposed matrix A:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_u_inverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_U_INVERSE_TEST tests R8MAT_U_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N*N] = {
    1.0, 0.0, 0.0,  0.0,
    2.0, 3.0, 0.0,  0.0,
    4.0, 5.0, 6.0,  0.0,
    7.0, 8.0, 9.0, 10.0 };
  double *b;
  double *c;
  int i;

  cout << "\n";
  cout << "R8MAT_U_INVERSE_TEST\n";
  cout << "  R8MAT_U_INVERSE inverts an upper triangular R8MAT.\n";

  r8mat_print ( N, N, a, "  Input matrix A" );

  b = r8mat_u_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void r8mat_u1_inverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_U1_INVERSE_TEST tests R8MAT_U1_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 6
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N*N] = {
    1.0, 0.0, 0.0, 0.0, 0.0,  0.0,
    2.0, 1.0, 0.0, 0.0, 0.0,  0.0,
    0.0, 0.0, 1.0, 0.0, 0.0,  0.0,
    5.0, 0.0, 3.0, 1.0, 0.0,  0.0,
    0.0, 0.0, 0.0, 0.0, 1.0,  0.0,
   75.0, 0.0, 0.0, 6.0, 4.0,  1.0 };
  double *b;
  double *c;

  cout << "\n";
  cout << "R8MAT_U1_INVERSE_TEST\n";
  cout << "  R8MAT_U1_INVERSE inverts a unit upper triangular R8MAT.\n";

  r8mat_print ( N, N, a, "  Input matrix A" );

  b = r8mat_u1_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void r8mat_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_AB_NEW_TEST tests R8MAT_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 4

  double *a;
  double b = 2.0E+00;
  double c = 10.0E+00;
  int seed = 123456789;

  cout << "\n";
  cout << "R8MAT_UNIFORM_AB_NEW_TEST\n";
  cout << "  R8MAT_UNIFORM_AB_NEW returns a random R8MAT in [A,B].\n";
  cout << "\n";

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  The random R8MAT:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8plu_det_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8PLU_DET_TEST tests R8PLU_DET;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double b = 0.0;
  double c = 1.0;
  double det;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "R8PLU_DET_TEST\n";
  cout << "  R8PLU_DET determines the determinant of an R8MAT from its\n";
  cout << "  compressed PLU factors.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_ab_new ( N, N, b, c, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
  }
//
//  Compute the determinant.
//
  det = r8plu_det ( N, pivot, lu );

  cout << "\n";
  cout << "  The determinant = " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8plu_inverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8PLU_INVERSE_TEST tests R8PLU_INVERSE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double b[N*N];
  double *c;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "R8PLU_INVERSE_TEST\n";
  cout << "  R8PLU_INVERSE determines the inverse of an R8MAT from its\n";
  cout << "  compressed PLU factors.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( N, N, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
  }
//
//  Compute the inverse.
//
  r8plu_inverse ( N, pivot, lu, b );

  r8mat_print ( N, N, b, "  The inverse B:" );
//
//  Compute the product C = A * B.
//
  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] a;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void r8plu_mul_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8PLU_MUL_TEST tests R8PLU_MUL;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;
  double x[N];

  cout << "\n";
  cout << "R8PLU_MUL_TEST\n";
  cout << "  R8PLU_MUL computes the product A*x\n";
  cout << "  using the compressed PLU factors of A.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( N, N, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Set the right hand side B1.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) (i+1);
  }

  b = r8mat_mv_new ( N, N, a, x );

  r8vec_print ( N, b, "  The right hand side B (computed from A):" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );
//
//  Solve the system.
//
  r8plu_mul ( N, pivot, lu, x, b );

  r8vec_print ( N, b, "  The right hand side B (computed from PLU):" );

  delete [] a;
  delete [] b;

  return;
# undef N
}
//****************************************************************************80

void r8plu_sol_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8PLU_SOL_TEST tests R8PLU_SOL;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;
  double x[N];

  cout << "\n";
  cout << "R8PLU_SOL_TEST\n";
  cout << "  R8PLU_SOL solves the linear system A*x=b\n";
  cout << "  using the compressed PLU factors of A.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( N, N, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Set the desired solution.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) (i+1);
  }
//
//  Set the right hand side.
//
  b = r8mat_mv_new ( N, N, a, x );

  r8vec_print ( N, b, "  The right hand side B (computed from A):" );
//
//  Destroy the desired solution (no cheating now!)
//
  r8vec_zero ( N, x );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  r8plu_sol ( N, pivot, lu, b, x );

  r8vec_print ( N, x, "  The computed solution X:" );

  delete [] a;
  delete [] b;

  return;
# undef N
}
//****************************************************************************80

void r8plu_to_r8mat_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8PLU_TO_R8MAT_TEST tests R8PLU_TO_R8MAT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double a2[N*N];
  double b = 0.0;
  double c = 1.0;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "R8PLU_TO_R8MAT_TEST\n";
  cout << "  R8PLU_TO_R8MAT determines the original matrix from\n";
  cout << "  the compressed PLU factors.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_ab_new ( N, N, b, c, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
  }
//
//  Display the gory details.
//
  i4vec_print ( N, pivot, "  The pivot vector P:" );

  r8mat_print ( N, N, lu, "  The compressed LU factors:" );
//
//  Recover the matrix from the PLU factors.
//
  r8plu_to_r8mat ( N, pivot, lu, a2 );

  r8mat_print ( N, N, a2, "  The recovered matrix A2:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8poly_degree_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DEGREE_TEST tests R8POLY_DEGREE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c1[4] = { 1.0, 2.0, 3.0, 4.0 }; 
  double c2[4] = { 1.0, 2.0, 3.0, 0.0 };
  double c3[4] = { 1.0, 2.0, 0.0, 4.0 };
  double c4[4] = { 1.0, 0.0, 0.0, 0.0 };
  double c5[4] = { 0.0, 0.0, 0.0, 0.0 };
  int d;
  int m;
 
  cout << "\n";
  cout << "R8POLY_DEGREE_TEST\n";
  cout << "  R8POLY_DEGREE determines the degree of an R8POLY.\n";

  m = 3;

  r8poly_print ( m, c1, "  The R8POLY:" );
  d = r8poly_degree ( m, c1 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  r8poly_print ( m, c2, "  The R8POLY:" );
  d = r8poly_degree ( m, c2 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  r8poly_print ( m, c3, "  The R8POLY:" );
  d = r8poly_degree ( m, c3 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  r8poly_print ( m, c4, "  The R8POLY:" );
  d = r8poly_degree ( m, c4 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  r8poly_print ( m, c5, "  The R8POLY:" );
  d = r8poly_degree ( m, c5 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  return;
}
//****************************************************************************80

void r8poly_deriv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DERIV_TEST tests R8POLY_DERIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double *c;
  double *cp;
  int d;
  double *x;

  cout << "\n";
  cout << "R8POLY_DERIV_TEST\n";
  cout << "  R8POLY_DERIV computes the coefficients of\n";
  cout << "  the derivative of a polynomial.\n";

  x = r8vec_indicator1_new ( N );

  c = roots_to_r8poly ( N, x );

  r8poly_print ( N, c, "  The initial polynomial" );

  for ( d = 0; d <= N; d++ )
  {
    cp = r8poly_deriv ( N, c, d );
    cout << "\n";
    cout << "  The derivative of order " << d << "\n";
    cout << "\n";
    r8poly_print ( N-d, cp, " " );
    delete [] cp;
  }

  delete [] c;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void r8poly_lagrange_coef_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_LAGRANGE_COEF_TEST tests R8POLY_LAGRANGE_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 3

  int i;
  int ipol;
  double *pcof;
  double *xpol;

  cout << "\n";
  cout << "R8POLY_LAGRANGE_COEF_TEST\n";
  cout << "  R8POLY_LAGRANGE_COEF returns the coefficients for\n";
  cout << "  a Lagrange basis polynomial.\n";

  xpol = r8vec_indicator1_new ( NPOL );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );

  for ( ipol = 1; ipol <= NPOL; ipol++ )
  {
    pcof = r8poly_lagrange_coef ( NPOL, ipol, xpol );

    cout << "\n";
    cout << "  Lagrange basis polynomial " << setw(4) << ipol << ":\n";
    cout << "\n";

    for ( i = 0; i < NPOL; i++ )
    {
      cout << setw(10) << pcof[i] << "  "
           << setw(4)  << i << "\n";
    }
    delete [] pcof;

  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void r8poly_lagrange_0_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_LAGRANGE_0_TEST tests R8POLY_LAGRANGE_0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  int ival;
  int nx;
  double wval;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  cout << "\n";
  cout << "R8POLY_LAGRANGE_0_TEST\n";
  cout << "  R8POLY_LAGRANGE_0 evaluates the Lagrange\n";
  cout << "  factor W(X) at a point.\n";
  cout << "\n";
  cout << "  The number of data points is " << NPOL << "\n";
//
//  Set the abscissas of the polynomials.
//
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );

  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
//
//  Evaluate W(X), W'(X).
//
  cout << "\n";
  cout << "      X          W(X)\n";
  cout << "\n";

  nx = 4 * NPOL - 1;

  for ( ival = 1; ival <= nx; ival++ )
  {
    xval = r8vec_even_select ( nx, xlo, xhi, ival );

    wval = r8poly_lagrange_0 ( NPOL, xpol, xval );

    cout << setw(12) << xval   << "  "
         << setw(12) << wval   << "\n";
  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void r8poly_lagrange_1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_LAGRANGE_1_TEST tests R8POLY_LAGRANGE_1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  double dwdx;
  int ival;
  int nx;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  cout << "\n";
  cout << "R8POLY_LAGRANGE_1_TEST\n";
  cout << "  R8POLY_LAGRANGE_1 evaluates the Lagrange\n";
  cout << "  factor W'(X) at a point.\n";
  cout << "\n";
  cout << "  The number of data points is " << NPOL << "\n";
//
//  Set the abscissas of the polynomials.
//
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );

  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
//
//  Evaluate W(X), W'(X).
//
  cout << "\n";
  cout << "      X          W'(X)\n";
  cout << "\n";

  nx = 4 * NPOL - 1;

  for ( ival = 1; ival <= nx; ival++ )
  {
    xval = r8vec_even_select ( nx, xlo, xhi, ival );

    dwdx = r8poly_lagrange_1 ( NPOL, xpol, xval );

    cout << setw(12) << xval   << "  "
         << setw(12) << dwdx   << "\n";
  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void r8poly_lagrange_2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_LAGRANGE_2_TEST tests R8POLY_LAGRANGE_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  double dw2dx2;
  int ival;
  int nx;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  cout << "\n";
  cout << "R8POLY_LAGRANGE_2_TEST\n";
  cout << "  R8POLY_LAGRANGE_2 evaluates the Lagrange\n";
  cout << "  factor W''(X) at a point.\n";
  cout << "\n";
  cout << "  The number of data points is " << NPOL << "\n";
//
//  Set the abscissas of the polynomials.
//
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );

  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
//
//  Evaluate W''.
//
  cout << "\n";
  cout << "      X          W''(X)\n";
  cout << "\n";

  nx = 4 * NPOL - 1;

  for ( ival = 1; ival <= nx; ival++ )
  {
    xval = r8vec_even_select ( nx, xlo, xhi, ival );

    dw2dx2 = r8poly_lagrange_2 ( NPOL, xpol, xval );

    cout << setw(12) << xval   << "  "
         << setw(12) << dw2dx2 << "\n";
  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void r8poly_lagrange_factor_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_LAGRANGE_FACTOR_TEST tests R8POLY_LAGRANGE_FACTOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  double dwdx;
  int i;
  double wval;
  double xhi;
  double xlo;
  double xpol[NPOL];
  double xval;

  cout << "\n";
  cout << "R8POLY_LAGRANGE_FACTOR_TEST\n";
  cout << "  R8POLY_LAGRANGE_FACTOR evaluates the Lagrange\n";
  cout << "  factor W(X) at a point.\n";
  cout << "\n";
  cout << "  For this test, we use " << NPOL << " functions.\n";
//
//  Set the abscissas of the polynomials.
//
  xlo = 0.0;
  xhi = ( double ) ( NPOL - 1 );
  for ( i = 0; i < NPOL; i++ )
  {
    xpol[i] = ( ( double ) ( NPOL - i ) * xlo + ( double ) i * xhi )
      / ( double ) ( NPOL );
  }

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
//
//  Evaluate W(X) and W'(X).
//
  cout << "\n";
  cout << "      X          W(X)          W'(X)\n";
  cout << "\n";

  for ( i = 0; i < 2 * NPOL - 2; i++ )
  {
    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, i );

    r8poly_lagrange_factor ( NPOL, xpol, xval, &wval, &dwdx );

    cout << setw ( 10 ) << xval << " "
         << setw ( 10 ) << wval << " "
         << setw ( 10 ) << dwdx << "\n";
  }

  return;
# undef NPOL
}
//****************************************************************************80

void r8poly_lagrange_val_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_LAGRANGE_VAL_TEST tests R8POLY_LAGRANGE_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  int i;
  int ipol;
  int ival;
  double p1;
  double p2;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  cout << "\n";
  cout << "R8POLY_LAGRANGE_VAL_TEST\n";
  cout << "  R8POLY_LAGRANGE_VAL evaluates a Lagrange\n";
  cout << "  interpolating polynomial at a point.\n";
  cout << "\n";
  cout << "  For this test, we use " << NPOL << " functions.\n";
//
//  Set the abscissas of the polynomials.
//
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );
  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
//
//  Evaluate the polynomials.
//
  cout << "\n";
  cout << "  Here are the values of the functions at\n";
  cout << "  several points:\n";
  cout << "\n";
  cout << "      X          L1          L2          L3      L4          L5\n";
  cout << "\n";

  for ( ival = 0; ival < 2 * NPOL - 1; ival++ )
  {

    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, ival );
    cout << setw(10) << xval << "  ";

    for ( ipol = 0; ipol < NPOL; ipol++ )
    {
      r8poly_lagrange_val ( NPOL, ipol, xpol, xval, &p1, &p2 );
      cout << setw(10) << p1 << "  ";
    }
    cout << "\n";
  }
  cout << "\n";
  cout << "  And the derivatives:\n";
  cout << "\n";
  cout << "      X          L'1         L'2         L'3     L'4         L'5\n";
  cout << "\n";

  for ( ival = 0; ival < 2 * NPOL - 1; ival++ )
  {
    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, ival );
    cout << setw(10) << xval << " ";

    for ( ipol = 0; ipol < NPOL; ipol++ )
    {
      r8poly_lagrange_val ( NPOL, ipol, xpol, xval, &p1, &p2 );
      cout << setw(10) << p2 << " ";
    }
    cout << "\n";
  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void r8poly_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PRINT_TEST tests R8POLY_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = { 2.0, -3.4, 56.0, 0.0, 0.78, 9.0 };
  int m = 5;

  cout << "\n";
  cout << "R8POLY_PRINT_TEST\n";
  cout << "  R8POLY_PRINT prints an R8POLY.\n";

  r8poly_print ( m, c, "  The R8POLY:" );

  return;
}
//****************************************************************************80

void r8poly_value_horner_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_VALUE_HORNER_TEST tests R8POLY_VALUE_HORNER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c[5] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  int i;
  int m = 4;
  int n = 16;
  double p;
  double *x;
  double x_hi;
  double x_lo;

  cout << "\n";
  cout << "R8POLY_VALUE_HORNER_TEST\n";
  cout << "  R8POLY_VALUE_HORNER evaluates a polynomial at\n";
  cout << "  one point, using Horner's method.\n";

  r8poly_print ( m, c, "  The polynomial coefficients:" );

  x_lo = 0.0;
  x_hi = 5.0;
  x = r8vec_linspace_new ( n, x_lo, x_hi );

  cout << "\n";
  cout << "   I    X    P(X)\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    p = r8poly_value_horner ( m, c, x[i] );
    cout << "  " << setw(2) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(14) << p << "\n";
  }

  delete [] x;

  return;
}
//****************************************************************************80

void r8poly_values_horner_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_VALUES_HORNER_TEST tests R8POLY_VALUES_HORNER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  double c[5] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  int i;
  int m = 4;
  int n = 16;
  double *p;
  double *x;
  double x_hi;
  double x_lo;

  cout << "\n";
  cout << "R8POLY_VALUES_HORNER_TEST\n";
  cout << "  R8POLY_VALUES_HORNER evaluates a polynomial at a\n";
  cout << "  point, using Horner's method.\n";

  r8poly_print ( m, c, "  The polynomial:" );

  x_lo = 0.0;
  x_hi = 5.0;
  x = r8vec_linspace_new ( n, x_lo, x_hi );

  p = r8poly_values_horner ( m, c, n, x );

  r8vec2_print ( n, x, p, "  X, P(X):" );

  delete [] p;
  delete [] x;

  return;
}
//****************************************************************************80

void r8poly2_ex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY2_EX_TEST tests R8POLY2_EX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  int ierror;
  double x1;
  double x2;
  double x3;
  double xmin;
  double y1;
  double y2;
  double y3;
  double ymin;

  cout << "\n";
  cout << "R8POLY2_EX_TEST\n";
  cout << "  R8POLY2_EX finds the extreme value\n";
  cout << "  of a parabola determined by three points.\n";

  a =  2.0;
  b = -4.0;
  c = 10.0;

  x1 = 1.0;
  y1 = a * x1 * x1 + b * x1 + c;
  x2 = 2.0;
  y2 = a * x2 * x2 + b * x2 + c;
  x3 = 3.0;
  y3 = a * x3 * x3 + b * x3 + c;

  cout << "\n";
  cout << "  Parabolic coefficients A = "
    << a << ", B = " << b << ", c = " << c << "\n";
  cout << "\n";
  cout << "  X, Y data:\n";
  cout << "\n";
  cout << "  " << x1 << "  " << y1;
  cout << "  " << x2 << "  " << y2;
  cout << "  " << x3 << "  " << y3;

  a = 0.0;
  b = 0.0;
  c = 0.0;

  ierror = r8poly2_ex ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  if ( ierror == 0 )
  {
    cout << "\n";
    cout << "  R8POLY2_EX returns XMIN = "
      << xmin << ", YMIN = " << ymin << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  R8POLY2_EX returns error code " << ierror << ".\n";
  }

  return;
}
//****************************************************************************80

void r8poly2_ex2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY2_EX2_TEST tests R8POLY2_EX and R8POLY2_EX2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  int ierror;
  double x1;
  double x2;
  double x3;
  double xmin;
  double y1;
  double y2;
  double y3;
  double ymin;

  cout << "\n";
  cout << "R8POLY2_EX2_TEST\n";
  cout << "  R8POLY2_EX2 finds the extreme value\n";
  cout << "  of a parabola determined by three points.\n";

  a =  2.0;
  b = -4.0;
  c = 10.0;

  x1 = 1.0;
  y1 = a * x1 * x1 + b * x1 + c;
  x2 = 2.0;
  y2 = a * x2 * x2 + b * x2 + c;
  x3 = 3.0;
  y3 = a * x3 * x3 + b * x3 + c;

  cout << "\n";
  cout << "  Parabolic coefficients A = "
    << a << ", B = " << b << ", c = " << c << "\n";
  cout << "\n";
  cout << "  X, Y data:\n";
  cout << "\n";
  cout << "  " << x1 << "  " << y1;
  cout << "  " << x2 << "  " << y2;
  cout << "  " << x3 << "  " << y3;

  a = 0.0;
  b = 0.0;
  c = 0.0;

  ierror = r8poly2_ex2 ( x1, y1, x2, y2, x3, y3, &xmin, &ymin, &a, &b, &c );

  if ( ierror == 0 )
  {
    cout << "\n";
    cout << "  R8POLY2_EX2 returns XMIN = "
      << xmin << ", YMIN = " << ymin << "\n";
    cout << "  and A = " << a << ", B = " << b << ", c = " << c << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  R8POLY2_EX2 returns error code " << ierror << ".\n";
  }

  return;
}
//****************************************************************************80

void r8poly2_val_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY2_VAL_TEST tests R8POLY2_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double x;
  double x1;
  double x2;
  double x3;
  double y;
  double y1;
  double y2;
  double y3;
  double yp;
  double ypp;

  cout << "\n";
  cout << "R8POLY2_VAL_TEST\n";
  cout << "  R8POLY2_VAL evaluates a parabola given\n";
  cout << "  3 data points.\n";
  cout << "\n";
  cout << "  Our parabola will be 2*x^2 + 3 * x + 1.\n";
  cout << "\n";
  cout << "  Case 1: 3 distinct data points:\n";
  cout << "\n";

  x1 = -1.0;
  x2 = 1.0;
  x3 = 3.0;

  r8poly2_val_f ( x1, &y1, &yp, &ypp );
  r8poly2_val_f ( x2, &y2, &yp, &ypp );
  r8poly2_val_f ( x3, &y3, &yp, &ypp );

  cout << "  " << x1 << " " << y1 << "\n";
  cout << "  " << x2 << " " << y2 << "\n";
  cout << "  " << x3 << " " << y3 << "\n";

  cout << "\n";
  cout << "  Sampled data:\n";
  cout << "\n";
  cout << "  X, Y, Y', Y''\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    cout << "  " << x << "  " << y << "  " << yp << "  " << ypp << "\n";
  }

  cout << "\n";
  cout << "  Case 2: X1=X2, X3 distinct:\n";
  cout << "\n";

  x1 = -1.0;
  x2 = -1.0;
  x3 = 3.0;

  r8poly2_val_f ( x1, &y1, &y2, &ypp );
  r8poly2_val_f ( x3, &y3, &yp, &ypp );

  cout << "  " << x1 << "  " << y1 << "\n";
  cout << "  " << x2 << "  " << y2 << "\n";
  cout << "  " << x3 << "  " << y3 << "\n";

  cout << "\n";
  cout << "  Sampled data:\n";
  cout << "\n";
  cout << "   X, Y, Y', Y''\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    cout << "  " << x << "  " << y << "  " << yp << "  " << ypp << "\n";
  }

  cout << "\n";
  cout << "  Case 3: X1=X2=X3:\n";
  cout << "\n";

  x1 = -1.0;
  x2 = -1.0;
  x3 = -1.0;

  r8poly2_val_f ( x1, &y1, &y2, &y3 );

  cout << "  " << x1 << "  " << y1 << "\n";
  cout << "  " << x2 << "  " << y2 << "\n";
  cout << "  " << x3 << "  " << y3 << "\n";

  cout << "\n";
  cout << "  Sampled data:\n";
  cout << "\n";
  cout << "  X, Y, Y', Y''\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    cout << "  " << x << "  " << y << "  " << yp << "  " << ypp << "\n";
  }

  return;
}
//****************************************************************************80

void r8poly2_val_f ( double x, double *y, double *yp, double *ypp )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY2_VAL_F evaluates a parabola for us.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  *y = 2.0 * x * x + 3.0 * x + 1.0;
  *yp = 4.0 * x + 3.0;
  *ypp = 4.0;

  return;
}
//****************************************************************************80

void r8poly2_val2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY2_VAL2_TEST tests R8POLY2_VAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 5

  int i;
  int left;
  double xdata[NDATA];
  double xval;
  double ydata[NDATA];
  double yval;
  double zdata[NDATA];
  double zval;

  cout << "\n";
  cout << "R8POLY2_VAL2_TEST\n";
  cout << "  R8POLY2_VAL2 evaluates parabolas through\n";
  cout << "  3 points in a table\n";
  cout << "\n";
  cout << "  Our data tables will actually be parabolas:\n";
  cout << "    A: 2*x^2 + 3 * x + 1.\n";
  cout << "    B: 4*x^2 - 2 * x + 5.\n";
  cout << "\n";

  for ( i = 0; i < NDATA; i++ )
  {
    xval = 2.0 * ( double ) ( i + 1 );
    xdata[i] = xval;
    ydata[i] = 2.0 * xval * xval + 3.0 * xval + 1.0;
    zdata[i] = 4.0 * xval * xval - 2.0 * xval + 5.0;
    cout << setw(6)  << i << " "
         << setw(10) << xdata[i] << "  "
         << setw(10) << ydata[i] << "  "
         << setw(10) << zdata[i] << "\n";
  }

  cout << "\n";
  cout << "  Interpolated data:\n";
  cout << "\n";
  cout << "  LEFT, X, Y1, Y2\n";
  cout << "\n";

  for ( i = 0; i <= 4; i++ )
  {
    xval = ( double ) ( 2 * i + 1 );
    left = i;
    if ( NDATA - 3 < left )
    {
      left = NDATA - 3;
    }
    if ( left < 0 )
    {
      left = 0;
    }
    r8poly2_val2 ( NDATA, xdata, ydata, left, xval, &yval );
    r8poly2_val2 ( NDATA, xdata, zdata, left, xval, &zval );

    cout << setw(6)  << left << "  "
         << setw(10) << xval << "  "
         << setw(10) << yval << "  "
         << setw(10) << zval << "\n";
  }

  return;
# undef NDATA
}
//****************************************************************************80

void r8rmat_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_NEW_TEST tests R8RMAT_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double **a;
  double **b;
  int i;
  int j;
  int k;
  int m;
  int n;

  cout << "\n";
  cout << "R8RMAT_NEW_TEST:\n";
  cout << "  R8RMAT_NEW dynamically creates a 2D array.\n";
  cout << "  Array entries can be addressed using the\n";
  cout << "  notation \"a[i][j]\".\n";
//
//  These dimensions could be entered by the user; they could depend on
//  some other calculation; or they could be changed repeatedly during this
//  computation, as long as old memory is deleted by R8MAT_DELETE and new memory
//  requested by R8RMAT_NEW.
//
  m = 4;
  n = 5;
//
//  Allocate memory.
//
  cout << "\n";
  cout << "  Allocating memory for array A of size " << m << " by " << n << ".\n";

  a = r8rmat_new ( m, n );

  cout << "\n";
  cout << "  Assigning values to A.\n";
//
//  Store values in A.
//
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = ( double ) ( 10 * i + j );
    }
  }
//
//  Print A.
//
  cout << "\n";
  cout << "  Dynamically allocated matrix A:\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << a[i][j];
    }
    cout << "\n";
  }
//
//  Create a new matrix B to store A' * A.
//
  b = r8rmat_new ( n, n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        b[i][j] = b[i][j] + a[k][i] * a[k][j];
      }
    }
  }
//
//  Print the matrix.
//
  cout << "\n";
  cout << "  Dynamically allocated matrix B = A' * A:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << b[i][j];
    }
    cout << "\n";
  }
//
//  Free memory.
//
  r8rmat_delete ( m, n, a );
  r8rmat_delete ( n, n, b );

  return;
}
//****************************************************************************80

void r8rmat_to_r8mat_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_TO_R8MAT_TEST tests R8RMAT_TO_R8MAT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double **b;
  double *c;
  int m = 5;
  int n = 4;

  cout << "\n";
  cout << "R8RMAT_TO_R8MAT_TEST\n";
  cout << "  R8RMAT_TO_R8MAT converts an R8RMAT to an R8MAT.\n";
  cout << "\n";
  cout << "  Data is of order (" << m << "," << n << ").\n";
//
//  Set the R8MAT.
//
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  The R8MAT A:" );
//
//  Convert.
//
  b = r8mat_to_r8rmat ( m, n, a );
  r8rmat_print ( m, n, b, "  The R8RMAT B:" );
//
//  Recover the matrix.
//
  c = r8rmat_to_r8mat ( m, n, b );
  r8mat_print ( m, n, c, "  The R8MAT C:" );
//
//  Free memory.
//
  free ( a );
  r8rmat_delete ( m, n, b );
  free ( c );

  return;
}
//****************************************************************************80

void r8row_max_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_MAX_TEST tests R8ROW_MAX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *amax;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "R8ROW_MAX_TEST\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_MAX computes maximums;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  amax = r8row_max ( M, N, a );

  r8vec_print ( M, amax, "  The row maximums:" );

  delete [] amax;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8row_mean_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_MEAN_TEST tests R8ROW_MEAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *mean;
  double *rowsum;

  cout << "\n";
  cout << "R8ROW_MEAN_TEST\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_MEAN computes means;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  mean = r8row_mean ( M, N, a );

  r8vec_print ( M, mean, "  Row means:" );

  delete [] mean;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8row_min_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_MIN_TEST tests R8ROW_MIN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *amin;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "R8ROW_MIN_TEST\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_MIN computes minimums;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  amin = r8row_min ( M, N, a );

  r8vec_print ( M, amin, "  The row minimums:" );

  delete [] amin;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8row_sum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_SUM_TEST tests R8ROW_SUM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *rowsum;

  cout << "\n";
  cout << "R8ROW_SUM_TEST\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_SUM computes sums;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  rowsum = r8row_sum ( M, N, a );

  r8vec_print ( M, rowsum, "  The row sums:" );

  delete [] rowsum;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8row_swap_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_SWAP_TEST tests R8ROW_SWAP;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int row1;
  int row2;
  int j;
  int k;

  cout << "\n";
  cout << "R8ROW_SWAP_TEST\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_SWAP swaps two rows;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  row1 = 1;
  row2 = 3;

  cout << "\n";
  cout << "  Swap rows " << row1 << " and " << row2 << "\n";

  r8row_swap ( M, N, a, row1, row2 );

  r8mat_print ( M, N, a, "  The modified matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8row_to_r8vec_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_TO_R8VEC_TEST tests R8ROW_TO_R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *x;

  cout << "\n";
  cout << "R8ROW_TO_R8VEC_TEST\n";
  cout << "  R8ROW_TO_R8VEC converts an array of rows into a vector.\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) ( 10 * i + j );
    }
  }

  r8mat_print ( M, N, a, "  The array of rows:" );

  x = r8row_to_r8vec ( M, N, a );

  r8vec_print ( M*N, x, "  The resulting vector of rows:" );

  delete [] x;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8row_variance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_VARIANCE_TEST tests R8ROW_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *variance;

  cout << "\n";
  cout << "R8ROW_VARIANCE_TEST\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_VARIANCE computes variances;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  variance = r8row_variance ( M, N, a );

  cout << "\n";
  cout << "  Row variances:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(10) << variance[i] << "\n";
  }

  delete [] variance;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8slmat_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8SLMAT_PRINT_TEST tests R8SLMAT_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double *a;
  double a1[21] = {
    21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          32.0, 42.0, 52.0, 62.0, 72.0,
                43.0, 53.0, 63.0, 73.0,
                      54.0, 64.0, 74.0,
                            65.0, 75.0,
                                  76.0 };
  double a2[15] = {
    21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          32.0, 42.0, 52.0, 62.0, 72.0,
                43.0, 53.0, 63.0, 73.0 };
  double a3[6] = {
    21.0, 31.0, 41.0,
          32.0, 42.0,
                43.0 };
  int m;
  int m_test[TEST_NUM] = { 7, 7, 4 };
  int n;
  int n_test[TEST_NUM] = { 7, 3, 7 };
  int size;
  int size_test[TEST_NUM] = { 21, 15, 6 };
  int test;

  cout << "\n";
  cout << "R8SLMAT_PRINT_TEST\n";
  cout << "  R8SLMAT_PRINT prints a strictly lower triangular matrix\n";
  cout << "  stored compactly.  Only the (possibly) nonzero \n";
  cout << "  elements are printed.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    m = m_test[test];
    n = n_test[test];
    size = size_test[test];
    a = new double[size];

    if ( test == 0 )
    {
      r8vec_copy ( size, a1, a );
    }
    else if ( test == 1 )
    {
      r8vec_copy ( size, a2, a );
    }
    else if ( test == 2 )
    {
      r8vec_copy ( size, a3, a );
    }

    r8slmat_print ( m, n, a, "  R8SLMAT:" );

    delete [] a;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8vec_amax_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_AMAX_TEST tests R8VEC_AMAX;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_AMAX_TEST\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_AMAX:      maximum magnitude entry;\n";

  b = - ( double ) N;
  c =  ( double ) N;

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  aval = r8vec_amax ( N, a );
  cout << "  Maximum absolute:         " << aval << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_amin_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_AMIN_TEST tests R8VEC_AMIN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_AMIN_TEST\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_AMIN:      minimum magnitude entry.\n";

  b = - ( double ) N;
  c =  ( double ) N;

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  aval = r8vec_amin ( N, a );
  cout << "  Minimum absolute:         " << aval << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_bracket_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET_TEST tests R8VEC_BRACKET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define TEST_NUM 6

  int i;
  int left;
  int right;
  int test;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "R8VEC_BRACKET_TEST\n";
  cout << "  R8VEC_BRACKET finds a pair of entries in a\n";
  cout << "  sorted real array which bracket a value.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!)" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    xval = xtest[test];

    cout << "\n";
    cout << "  Search for XVAL = " << xval << "\n";

    r8vec_bracket ( N, x, xval, left, right );

    cout << "  X[" << left  << "-1] = " << x[left-1]  << "\n";
    cout << "  X[" << right << "-1] = " << x[right-1] << "\n";
  }

  return;

# undef N
# undef TEST_NUM
}
//****************************************************************************80

void r8vec_bracket2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET2_TEST tests R8VEC_BRACKET2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define TEST_NUM 6

  int i;
  int left;
  int right;
  int start;
  int test;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "R8VEC_BRACKET2_TEST\n";
  cout << "  R8VEC_BRACKET2 finds a pair of entries in a\n";
  cout << "  sorted R8VEC which bracket a value.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!)" );

  for ( test = 0; test < TEST_NUM; test++ )
  {

    xval = xtest[test];

    cout << "\n";
    cout << "  Search for XVAL = " << xval << "\n";

    if ( 0 < left )
    {
      start = left;
    }
    else
    {
      start = ( N + 1 ) / 2;
    }

    cout << "  Start = " << start << "\n";

    r8vec_bracket2 ( N, x, xval, start, left, right );

    cout << "  Left =  " << left  << "\n";
    cout << "  Right = " << right << "\n";

    if ( 1 <= left )
    {
      cout << "  X[" << left  << "-1] = " << x[left-1]  << "\n";
    }

    if ( 1 <= right )
    {
      cout << "  X[" << right << "-1] = " << x[right-1] << "\n";
    }
  }
  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void r8vec_bracket3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET3_TEST tests R8VEC_BRACKET3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define TEST_NUM 6

  int i;
  int itest;
  int left;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "R8VEC_BRACKET3_TEST\n";
  cout << "  R8VEC_BRACKET3 finds a pair of entries in a\n";
  cout << "  sorted real array which bracket a value.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!):" );

  left = ( N - 1 ) / 2;

  for ( itest = 0; itest < TEST_NUM; itest++ )
  {
    xval = xtest[itest];

    cout << "\n";
    cout << "  Search for XVAL = " << xval << "\n";

    cout << "  Starting guess for interval is = " << left << "\n";

    r8vec_bracket3 ( N, x, xval, left );

    cout << "  Nearest interval:\n";
    cout << "   X[" << left   << "]= " << x[left  ] << "\n";
    cout << "   X[" << left+1 << "]= " << x[left+1] << "\n";
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void r8vec_bracket5_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET5_TEST tests R8VEC_BRACKET5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int left;
  int n = 10;
  int right;
  int test;
  int test_num = 6;
  double *x;
  double xtest[6] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "R8VEC_BRACKET5_TEST\n";
  cout << "  R8VEC_BRACKET5 finds a pair of entries in a\n";
  cout << "  sorted R8VEC which bracket a value.\n";

  x = r8vec_indicator1_new ( n );
  x[5] = x[4];

  r8vec_print ( n, x, "  Sorted array:" );

  cout << "\n";
  cout << "        LEFT                   RIGHT\n";
  cout << "      X(LEFT)       XVAL     X(RIGHT)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    xval = xtest[test];

    left = r8vec_bracket5 ( n, x, xval );

    if ( left == -1 )
    {
      cout << "  " << setw(10) <<left << "\n";
      cout << "              " << setw(10) << xval << "  (Not bracketed!)\n";
    }
    else
    {
      right = left + 1;
      cout << "  " << setw(10) << left
           << "  " << "          "
           << "  " << setw(10) << right << "\n";
      cout << "  " << setw(10) << x[left]
           << "  " << setw(10) << xval
           << "  " << setw(10) << x[right] << "\n";
    }
  }

  delete [] x;

  return;
}
//****************************************************************************80

void r8vec_chebyspace_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CHEBYSPACE_NEW_TEST tests R8VEC_CHEBYSPACE_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double *r;
  double r1;
  double r2;

  cout << "\n";
  cout << "R8VEC_CHEBYSPACE_NEW_TEST\n";
  cout << "  R8VEC_CHEBYSPACE_NEW computes N Chebyshev points in [R1,R2].\n";

  r1 = -1.0;
  r2 = +1.0;
  n = 5;

  r = r8vec_chebyspace_new ( n, r1, r2 );

  cout << "\n";
  cout << "  N = " << n << ", R1 = " << r1 << ", R2 = " << r2 << "\n";

  r8vec_print ( n, r, "  Chebyshev points:" );

  delete [] r;

  r1 =   0.0;
  r2 = +10.0;
  n = 7;

  r = r8vec_chebyspace_new ( n, r1, r2 );

  cout << "\n";
  cout << "  N = " << n << ", R1 = " << r1 << ", R2 = " << r2 << "\n";

  r8vec_print ( n, r, "  Chebyshev points:" );

  delete [] r;

  return;
}
//****************************************************************************80

void r8vec_concatenate_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CONCATENATE_NEW_TEST tests R8VEC_CONCATENATE_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n1 = 5;
  int n2 = 3;
  int n3 = n1 + n2;

  double a1[5] = { 91.1, 31.2, 71.3, 51.4, 31.5 };
  double a2[3] = { 42.6, 22.7, 12.8 };
  double *a3;

  cout << "\n";
  cout << "R8VEC_CONCATENATE_NEW_TEST\n";
  cout << "  R8VEC_CONCATENATE_NEW concatenates two R8VECs\n";

  r8vec_print ( n1, a1, "  Array 1:" );
  r8vec_print ( n2, a2, "  Array 2:" );
  a3 = r8vec_concatenate_new ( n1, a1, n2, a2 );
  r8vec_print ( n3, a3, "  Array 3 = Array 1 + Array 2:" );

  delete [] a3;

  return;
}
//****************************************************************************80

void r8vec_convolution_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CONVOLUTION_TEST tests R8VEC_CONVOLUTION
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2012
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 3

  int m = M;
  int n = N;

  double x[M] = { 1.0, 2.0, 3.0, 4.0 };
  double y[N] = { -1.0, 5.0, 3.0 };
  double *z;
  double z_correct[M+N-1] = { -1.0, 3.0, 10.0, 17.0, 29.0, 12.0 };

  cout << "\n";
  cout << "R8VEC_CONVOLUTION_TEST\n";
  cout << "  R8VEC_CONVOLUTION computes the convolution\n";
  cout << "  of two vectors.\n";

  r8vec_print ( m, x, "  The factor X:" );
  r8vec_print ( n, y, "  The factor Y:" );

  z = r8vec_convolution ( m, x, n, y );

  r8vec_print ( m + n - 1, z, "  The convolution z = x star y:" );

  r8vec_print ( m + n - 1, z_correct, "  Correct answer:" );

  delete [] z;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8vec_convolution_circ_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CONVOLUTION_CIRC_TEST tests R8VEC_CONVOLUTION_CIRC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double x[N] = { 1.0, 2.0, 3.0, 4.0 };
  double y[N] = { 1.0, 2.0, 4.0, 8.0 };
  double *z;
  double z_correct[N] = { 37.0, 44.0, 43.0, 26.0 };

  cout << "\n";
  cout << "R8VEC_CONVOLUTION_CIRC_TEST\n";
  cout << "  R8VEC_CONVOLUTION_CIRC computes the circular convolution\n";
  cout << "  of two vectors.\n";

  r8vec_print ( N, x, "  The factor X:" );
  r8vec_print ( N, y, "  The factor Y:" );

  z = r8vec_convolution_circ ( N, x, y );

  r8vec_print ( N, z, "  The circular convolution z = xCCy:" );

  r8vec_print ( N, z_correct, "  Correct answer:" );

  delete [] z;

  return;
# undef N
}
//****************************************************************************80

void r8vec_dif_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIF_TEST tests R8VEC_DIF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *cof;
  double fdif;
  double h = 0.01;
  int i;
  int n = 4;
  double x = 1.0;
  double xi;

  cout << "\n";
  cout << "R8VEC_DIF_TEST\n";
  cout << "  R8VEC_DIF estimates derivatives.\n";
  cout << "\n";
  cout << "  Estimate the derivative of order N = " << n << "\n";
  cout << "  Using H = " << h << "\n";
  cout << "  at argument X = " << x << "\n";
//
//  Get the coefficients.
//
  cof = r8vec_dif ( n, h );

  r8vec_print ( n+1, cof, "  The difference coefficients:" );

  fdif = 0.0;
  for ( i = 0; i <= n; i++ )
  {
    xi = x + ( double ) ( 2 * i - n ) * h;
    fdif = fdif + cof[i] * r8vec_dif_f ( xi );
  }

  cout << "\n";
  cout << "  Estimate is FDIF = " << fdif << "\n";

  delete [] cof;

  return;
}
//****************************************************************************80

double r8vec_dif_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIF_F evaluates the function used in R8VEC_DIF_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double value;

  value = exp ( x );

  return value;
}
//****************************************************************************80

void r8vec_direct_product_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIRECT_PRODUCT_TEST tests R8VEC_DIRECT_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int factor_num = 3;
  int point_num = 24;

  int factor_index;
  int factor_order;
  double *factor_value;
  int i;
  int j;
  double x[factor_num*point_num];

  cout << "\n";
  cout << "R8VEC_DIRECT_PRODUCT_TEST\n";
  cout << "  R8VEC_DIRECT_PRODUCT forms the entries of a\n";
  cout << "  direct product of a given number of R8VEC factors.\n";

  for ( j = 0; j < point_num; j++ )
  {
    for ( i = 0; i < factor_num; i++ )
    {
      x[i+j*factor_num] = 0.0;
    }
  }

  for ( factor_index = 0; factor_index < factor_num; factor_index++ )
  {
    if ( factor_index == 0 )
    {
      factor_order = 4;
      factor_value = new double[factor_order];
      factor_value[0] = 1.0;
      factor_value[1] = 2.0;
      factor_value[2] = 3.0;
      factor_value[3] = 4.0;
    }
    else if ( factor_index == 1 )
    {
      factor_order = 3;
      factor_value = new double[factor_order];
      factor_value[0] = 50.0;
      factor_value[1] = 60.0;
      factor_value[2] = 70.0;
    }
    else if ( factor_index == 2 )
    {
      factor_order = 2;
      factor_value = new double[factor_order];
      factor_value[0] = 800.0;
      factor_value[1] = 900.0;
    }

    r8vec_direct_product ( factor_index, factor_order, factor_value,
      factor_num, point_num, x );

    delete [] factor_value;
  }

  cout << "\n";
  cout << "     J         X(1)      X(2)      X(3)\n";
  cout << "\n";

  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j << "  ";
    for ( i = 0; i < factor_num; i++ )
    {
      cout << "  " << setw(8) << x[i+j*factor_num];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_direct_product2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIRECT_PRODUCT2_TEST tests R8VEC_DIRECT_PRODUCT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int factor_num = 3;
  int point_num = 24;

  int factor_index;
  int factor_order;
  double *factor_value;
  int i;
  int j;
  double w[point_num];

  cout << "\n";
  cout << "R8VEC_DIRECT_PRODUCT2_TEST\n";
  cout << "  R8VEC_DIRECT_PRODUCT2 forms the entries of a\n";
  cout << "  direct product of a given number of R8VEC factors.\n";

  for ( j = 0; j  < point_num; j++ )
  {
    w[j] = 1.0;
  }

  for ( factor_index = 0; factor_index < factor_num; factor_index++ )
  {
    if ( factor_index == 0 )
    {
      factor_order = 4;
      factor_value = new double[factor_order];
      factor_value[0] = 2.0;
      factor_value[1] = 3.0;
      factor_value[2] = 5.0;
      factor_value[3] = 7.0;
    }
    else if ( factor_index == 1 )
    {
      factor_order = 3;
      factor_value = new double[factor_order];
      factor_value[0] = 11.0;
      factor_value[1] = 13.0;
      factor_value[2] = 17.0;
    }
    else if ( factor_index == 2 )
    {
      factor_order = 2;
      factor_value = new double[factor_order];
      factor_value[0] = 19.0;
      factor_value[1] = 21.0;
    }

    r8vec_direct_product2 ( factor_index, factor_order, factor_value,
      factor_num, point_num, w );

    delete [] factor_value;
  }

  cout << "\n";
  cout << "     J         W(J)\n";
  cout << "\n";

  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j << "  "
         << "  " << setw(8) << w[j] << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_even_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN_TEST tests R8VEC_EVEN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *x;
  double xhi = 99.0;
  double xlo = 0.0;

  cout << "\n";
  cout << "R8VEC_EVEN_TEST\n";
  cout << "  R8VEC_EVEN computes N evenly spaced values\n";
  cout << "  between XLO and XHI.\n";
  cout << "\n";
  cout << "  XLO = " << xlo << "\n";
  cout << "  XHI = " << xhi << "\n";
  cout << "  while N = " << N << "\n";

  x = r8vec_even_new ( N, xlo, xhi );

  r8vec_print ( N, x, "  Resulting array:" );

  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void r8vec_even2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN2_TEST tests R8VEC_EVEN2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NOLD 5
# define MAXVAL 20

  int i;
  int istar;
  int jstar;
  int nfill[NOLD-1] = { 4, 3, 5, 0 };
  int nval;
  double xold[NOLD] = { 0.0, 1.0, 5.0, 2.0, 0.0 };
  double xval[MAXVAL];

  cout << "\n";
  cout << "R8VEC_EVEN2_TEST:\n";
  cout << "  R8VEC_EVEN2 interpolates a specified number of\n";
  cout << "  points pairs of values in a vector.\n";
  cout << "\n";
  cout << "  Input data:\n";
  cout << "\n";
  for ( i = 1; i <= NOLD; i++ )
  {
    cout << "  " << setw(12) << xold[i-1] << "\n";
    if ( i < NOLD )
    {
      cout << "  (" << nfill[i-1] << ")\n";
    }
  }

  r8vec_even2 ( MAXVAL, nfill, NOLD, xold, nval, xval );

  cout << "\n";
  cout << "  Resulting vector:\n";
  cout << "\n";

  istar = 1;
  jstar = 1;
  for ( i = 1; i <= nval; i++ )
  {
    if ( i == istar )
    {
      cout << "  " << '*' << "  " << xval[i-1] << "\n";

      if ( jstar < NOLD )
      {
        istar = istar + nfill[jstar-1] + 1;
        jstar = jstar + 1;
      }
    }
    else
    {
      cout << "     " << setw(12) << xval[i-1] << "\n";
    }
  }

  return;
# undef MAXVAL
# undef NOLD
}
//****************************************************************************80

void r8vec_expand_linear_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EXPAND_LINEAR_TEST tests R8VEC_EXPAND_LINEAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int fat = 3;
  int nfat;
  double x[N] = { 16.0, 4.0, 0.0, 4.0, 16.0, 36.0 };
  double *xfat;

  cout << "\n";
  cout << "R8VEC_EXPAND_LINEAR_TEST\n";
  cout << "  R8VEC_EXPAND_LINEAR linearly interpolates new data\n";
  cout << "  between old values.\n";
  cout << "\n";

  r8vec_print ( N, x, "  Original vector:" );

  cout << "\n";
  cout << "  Expansion factor is " << fat << "\n";

  xfat = r8vec_expand_linear ( N, x, fat );

  nfat = ( N - 1 ) * ( fat + 1 ) + 1;

  r8vec_print ( nfat, xfat, "  Fattened vector:" );

  delete [] xfat;

  return;
# undef N
}
//****************************************************************************80

void r8vec_frac_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_FRAC_TEST tests R8VEC_FRAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double afrac;
  int k;
  int seed;

  cout << "\n";
  cout << "R8VEC_FRAC_TEST\n";
  cout << "  R8VEC_FRAC: K-th smallest R8VEC entry;\n";

  seed = 123456789;

  a = r8vec_uniform_01_new ( N, seed );

  r8vec_print ( N, a, "  Array to search:" );

  cout << "\n";
  cout << "  Fractile  Value\n";
  cout << "\n";

  for ( k = 1; k < N; k = k + N/2 )
  {
    afrac = r8vec_frac ( N, a, k );
    cout << "  " << setw(6) << k
         << "  " << setw(14) << afrac << "\n";
  }

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_histogram_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_HISTOGRAM_TEST tests R8VEC_HISTOGRAM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define HISTO_NUM 20
# define N 1000

  double *a;
  double a_hi;
  double a_lo;
  double bin_hi;
  double bin_lo;
  int *histo_gram;
  int i;
  int seed = 123456789;
  int test;
  int test_num = 2;

  cout << "\n";
  cout << "R8VEC_HISTOGRAM_TEST\n";
  cout << "  R8VEC_HISTOGRAM histograms a real vector.\n";

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      cout << "\n";
      cout << "  Uniform data:\n";

      a_lo =  0.0;
      a_hi = +1.0;
      a = r8vec_uniform_01_new ( N, seed );
    }
    else if ( test == 2 )
    {
      cout << "\n";
      cout << "  Normal data:\n";
      a_lo = -3.0;
      a_hi = +3.0;
      a = r8vec_normal_01_new ( N, seed );
    }

    histo_gram = r8vec_histogram ( N, a, a_lo, a_hi, HISTO_NUM );

    cout << "\n";
    cout << "  Histogram of data:\n";
    cout << "\n";

    for ( i = 0; i <= HISTO_NUM+1; i++ )
    {
      if ( i == 0 )
      {
        cout << "  " << "          "
             << "  " << setw(10) << a_lo
             << "  " << setw(6)  << histo_gram[i] << "\n";
      }
      else if ( i <= HISTO_NUM )
      {
        bin_lo = ( ( double ) ( HISTO_NUM - i + 1 ) * a_lo
                 + ( double ) (             i - 1 ) * a_hi )
                 / ( double ) ( HISTO_NUM         );

        bin_hi = ( ( double ) ( HISTO_NUM - i     ) * a_lo
                 + ( double ) (             i     ) * a_hi )
                 / ( double ) ( HISTO_NUM         );

        cout << "  " << setw(10) << bin_lo
             << "  " << setw(10) << bin_hi
             << "  " << setw(6)  << histo_gram[i] << "\n";
      }
      else if ( i == HISTO_NUM+1 )
      {
        cout << "  " << setw(10) << a_hi
             << "  " << "          "
             << "  " << setw(6)  << histo_gram[i] << "\n";
      }
    }
    delete [] a;
    delete [] histo_gram;
  }

  return;
# undef HISTO_NUM
# undef N
}
//****************************************************************************80

void r8vec_house_column_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_HOUSE_COLUMN_TEST tests R8VEC_HOUSE_COLUMN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 5.0;
  double *h;
  double *ha;
  int i;
  int j;
  int k;
  int n = 4;
  int seed;
  double *v;

  cout << "\n";
  cout << "R8VEC_HOUSE_COLUMN_TEST\n";
  cout << "  R8VEC_HOUSE_COLUMN returns the compact form of\n";
  cout << "  a Householder matrix that packs a column\n";
  cout << "  of a matrix.\n";
//
//  Get a random matrix.
//
  seed = 123456789;

  a = r8mat_uniform_ab_new ( n, n, b, c, seed );

  r8mat_print ( n, n, a, "  Matrix A:" );

  for ( k = 1; k <= n-1; k++ )
  {
    cout << "\n";
    cout << "  Working on column K = " << k << "\n";

    v = r8vec_house_column ( n, a+(k-1)*n, k );

    h = r8mat_house_form ( n, v );

    r8mat_print ( n, n, h, "  Householder matrix H:" );

    ha = r8mat_mm_new ( n, n, n, h, a );

    r8mat_print ( n, n, ha, "  Product H*A:" );
//
//  If we set A := HA, then we can successively convert A to upper
//  triangular form.
//
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        a[i+j*n] = ha[i+j*n];
      }
    }

    delete [] h;
    delete [] ha;
    delete [] v;
  }

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_index_delete_all_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEX_DELETE_ALL_TEST tests R8VEC_INDEX_DELETE_ALL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "R8VEC_INDEX_DELETE_ALL_TEST\n";
  cout << "  R8VEC_INDEX_DELETE_ALL deletes all copies of a\n";
  cout << "  particular value.\n";

  xval = 8.0;
  r8vec_index_insert ( n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "  " << xval << "\n";
    r8vec_index_insert ( n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( n, x, indx, xval );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call R8VEC_INDEX_DELETE_ALL to delete all values of 7:\n";

  xval = 7.0;
  r8vec_index_delete_all ( n, x, indx, xval, n, x, indx );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8vec_index_delete_dupes_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEX_DELETE_DUPES_TEST tests R8VEC_INDEX_DELETE_DUPES
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "R8VEC_INDEX_DELETE_DUPES_TEST\n";
  cout << "  R8VEC_INDEX_DELETE_DUPES deletes duplicates.\n";

  xval = 8.0;
  r8vec_index_insert ( n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "  " << xval << "\n";
    r8vec_index_insert ( n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( n, x, indx, xval );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call R8VEC_INDEX_DELETE_DUPES to delete duplicates:\n";

  r8vec_index_delete_dupes ( n, x, indx, n, x, indx );

  cout << "\n";
  cout << "  Indexed list of unique entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8vec_index_delete_one_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEX_DELETE_ONE_TEST tests R8VEC_DELETE_ONE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "R8VEC_INDEX_DELETE_ONE_TEST\n";
  cout << "  R8VEC_INDEX_DELETE_ONE deletes one copies of a\n";
  cout << "  particular value.\n";

  xval = 8.0;
  r8vec_index_insert ( n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "  " << xval << "\n";
    r8vec_index_insert ( n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( n, x, indx, xval );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call R8VEC_INDEX_DELETE_ONE to delete one value of 8:\n";

  xval = 8.0;
  r8vec_index_delete_one ( n, x, indx, xval, n, x, indx );

  cout << "\n";
  cout << "  Indexed list of unique entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8vec_index_insert_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEX_INSERT_TEST tests R8VEC_INDEX_INSERT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "R8VEC_INDEX_INSERT_TEST\n";
  cout << "  R8VEC_INDEX_INSERT inserts values into an\n";
  cout << "  index sorted array.\n";

  xval = 8.0;
  r8vec_index_insert ( n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "  " << xval << "\n";
    r8vec_index_insert ( n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( n, x, indx, xval );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8vec_index_insert_unique_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEX_INSERT_UNIQUE_TEST tests R8VEC_INDEX_INSERT_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double b;
  double c;
  int equal;
  int i;
  int indx[N_MAX];
  int less;
  int more;
  int n;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "R8VEC_INDEX_INSERT_UNIQUE_TEST\n";
  cout << "  R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an\n";
  cout << "  index sorted array.\n";

  b = 0.0;
  c = 20.0;
  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( b, c, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "    " << setw(6) << xval << "\n";
    r8vec_index_insert_unique ( n, x, indx, xval );
  }

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Results of search for given XVAL:\n";
  cout << "\n";
  cout << "  XVAL  Less Equal More\n";
  cout << "\n";

  for ( i = 0; i <= N_MAX; i++ )
  {
    xval = ( double ) ( i );
    r8vec_index_search ( n, x, indx, xval, less, equal, more );
    cout << "  " << setw(6) << xval
         << "  " << setw(3) << less
         << "  " << setw(3) << equal
         << "  " << setw(3) << more << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8vec_index_order_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEX_ORDER_TEST tests R8VEC_INDEX_ORDER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  int i;
  int indx[N_MAX];
  int n;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "R8VEC_INDEX_ORDER_TEST\n";
  cout << "  R8VEC_INDEX_ORDER sorts an index sorted array.\n";
  cout << "\n";
  cout << "  Generate some random values:\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "  " << setw(6) << xval << "\n";
    r8vec_index_insert_unique ( n, x, indx, xval );
  }

  cout << "\n";
  cout << "  Indexed list of unique entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Now call R8VEC_INDEX_ORDER to carry out the sorting:\n";

  r8vec_index_order ( n, x, indx );

  r8vec_print ( n, x, "  X:" );

  return;
# undef N_MAX
}
//****************************************************************************80

void r8vec_index_search_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEX_SEARCH_TEST tests R8VEC_INDEX_SEARCH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double b;
  double c;
  int equal;
  int i;
  int indx[N_MAX];
  int less;
  int more;
  int n;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "R8VEC_INDEX_SEARCH_TEST\n";
  cout << "  R8VEC_INDEX_SEARCH searches for an entry \n";
  cout << "  with a given value.\n";
  cout << "\n";
  cout << "  Generate some random values:\n";
  cout << "\n";

  b = 0.0;
  c = 20.0;
  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( b, c, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "    " << setw(6) << xval << "\n";
    r8vec_index_insert_unique ( n, x, indx, xval );
  }

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Results of search for given XVAL:\n";
  cout << "\n";
  cout << "  XVAL  Less Equal More\n";
  cout << "\n";

  for ( i = 0; i <= N_MAX; i++ )
  {
    xval = ( double ) ( i );
    r8vec_index_search ( n, x, indx, xval, less, equal, more );
    cout << "  " << setw(6) << xval
         << "  " << setw(3) << less
         << "  " << setw(3) << equal
         << "  " << setw(3) << more << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8vec_index_sorted_range_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEX_SORTED_RANGE_TEST tests R8VEC_INDEX_SORTED_RANGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_hi;
  int i_lo;
  int *indx;
  int n = 20;
  double r[20];
  double r_lo;
  double r_hi;
  int seed;
  double t;
  int test;

  cout << "\n";
  cout << "R8VEC_INDEX_SORTED_RANGE_TEST\n";
  cout << "  R8VEC_INDEX_SORTED_RANGE seeks the range I_LO:I_HI\n";
  cout << "  of entries of sorted indexed R so that\n";
  cout << "  R_LO <= R(INDX(I)) <= R_HI for I_LO <= I <= I_HI.\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r8vec_uniform_01 ( n, seed, r );

    r8vec_print ( n, r, "  Array" );

    indx = r8vec_sort_heap_index_a_new ( n, r );

    cout << "\n";
    cout << "     I  INDX    R(INDX(I))\n";
    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << i
           << "  " << setw(4) << indx[i]
           << "  " << setw(14) << r[indx[i]] << "\n";
    }

    r_lo = r8_uniform_01 ( seed );
    r_hi = r8_uniform_01 ( seed );

    if ( r_hi < r_lo )
    {
      t = r_lo;
      r_lo = r_hi;
      r_hi = t;
    }

    r8vec_index_sorted_range ( n, r, indx, r_lo, r_hi, i_lo, i_hi );

    cout << "\n";
    if ( i_hi < i_lo )
    {
      cout << "  " << "R_LO" << "      "
           << "  " << setw(14) << r_lo << "\n";
      cout << "  " << "R_HI" << "      "
           << "  " << setw(14) << r_hi << "\n";
      cout << "  Empty range in R.\n";
    }
    else
    {
      cout << "  " << "R_LO" << "      "
           << "  " << setw(14) << r_lo << "\n";
      for ( i = i_lo; i <= i_hi; i++ )
      {
        cout << "  " << setw(4) << i
             << "  " << setw(4) << indx[i]
             << "  " << setw(14) << r[indx[i]] << "\n";
      }
      cout << "  " << "R_HI" << "      "
           << "  " << setw(14) << r_hi << "\n";
    }
    delete [] indx;
  }

  return;
}
//****************************************************************************80

void r8vec_indexed_heap_d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEXED_HEAP_D_TEST tests R8VEC_INDEXED_HEAP_D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double a[20] = {
    101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0,
    111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0 };
  int i;
  int indx[10] = {
    0, 10, 16, 4, 6, 12, 14, 2, 18, 8 };
  int m = 20;
  int n = 10;

  cout << "\n";
  cout << "R8VEC_INDEXED_HEAP_D_TEST\n";
  cout << "  R8VEC_INDEXED_HEAP_D creates a descending heap\n";
  cout << "  from an indexed vector.\n";
//
//  Print before.
//
  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  cout << "\n";
  cout << "  A(INDX):\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Heap the data.
//
  r8vec_indexed_heap_d ( n, a, indx );
//
//  Print afterwards.  Only INDX should change.
//
  r8vec_print ( m, a, "  The data vector (should NOT change):" );
  i4vec_print ( n, indx, "  The index vector (may change):" );
  cout << "\n";
  cout << "  A(INDX) is now a descending heap:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_indexed_heap_d_extract_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEXED_HEAP_D_EXTRACT_TEST tests R8VEC_INDEXED_HEAP_D_EXTRACT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  cout << "\n";
  cout << "R8VEC_INDEXED_HEAP_D_EXTRACT_TEST\n";
  cout << "  For an indexed R8VEC,\n";
  cout << "  R8VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;\n";
//
//  Set the data array.  To keep things easy, we will use the indicator vector.
//
  a = r8vec_indicator1_new ( m );
//
//  The index array will initially be a random subset of the numbers 1 to M,
//  in random order.
//
  n = 5;
  indx[0]  =  8;
  indx[1]  =  1;
  indx[2]  =  7;
  indx[3]  = 13;
  indx[4]  =  4;
  indx[5]  =  6;
  indx[6]  = 14;
  indx[7]  =  0;
  indx[8]  = 18;
  indx[9]  = 19;
  indx[10] =  2;

  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  cout << "\n";
  cout << "  A(INDX):\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Create a descending heap from the indexed array.
//
  r8vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  cout << "\n";
  cout << "  A(INDX) after heaping:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Insert five entries, and monitor the maximum.
//
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    cout << "\n";
    cout << "  Inserting value " << a[indx_insert] << "\n";

    r8vec_indexed_heap_d_insert ( n, a, indx, indx_insert );

    indx_max = r8vec_indexed_heap_d_max ( n, a, indx );

    cout << "  Current maximum is " << a[indx_max] << "\n";
  }
  r8vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  cout << "\n";
  cout << "  A(INDX) after insertions:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Extract the first 5 largest elements.
//
  cout << "\n";
  cout << "  Now extract the maximum several times.\n";
  cout << "\n";

  for ( i = 0; i < 5; i++ )
  {
    indx_extract = r8vec_indexed_heap_d_extract ( n, a, indx );
    cout << "  Extracting maximum element A[" << indx_extract
         << "] = " << a[indx_extract] << "\n";
  }
  r8vec_print ( m, a, "  The data vector after extractions:" );
  i4vec_print ( n, indx, "  The index vector after extractions:" );
  cout << "\n";
  cout << "  A(INDX) after extractions:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_indexed_heap_d_insert_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEXED_HEAP_D_INSERT_TEST tests R8VEC_INDEXED_HEAP_D_INSERT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  cout << "\n";
  cout << "R8VEC_INDEXED_HEAP_D_INSERT_TEST\n";
  cout << "  For an indexed R8VEC,\n";
  cout << "  R8VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.\n";
//
//  Set the data array.  To keep things easy, we will use the indicator vector.
//
  a = r8vec_indicator1_new ( m );
//
//  The index array will initially be a random subset of the numbers 1 to M,
//  in random order.
//
  n = 5;
  indx[0]  =  8;
  indx[1]  =  1;
  indx[2]  =  7;
  indx[3]  = 13;
  indx[4]  =  4;
  indx[5]  =  6;
  indx[6]  = 14;
  indx[7]  =  0;
  indx[8]  = 18;
  indx[9]  = 19;
  indx[10] =  2;

  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  cout << "\n";
  cout << "  A(INDX):\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Create a descending heap from the indexed array.
//
  r8vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  cout << "\n";
  cout << "  A(INDX) after heaping:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Insert five entries, and monitor the maximum.
//
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    cout << "\n";
    cout << "  Inserting value " << a[indx_insert] << "\n";

    r8vec_indexed_heap_d_insert ( n, a, indx, indx_insert );

    indx_max = r8vec_indexed_heap_d_max ( n, a, indx );

    cout << "  Current maximum is " << a[indx_max] << "\n";
  }
  r8vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  cout << "\n";
  cout << "  A(INDX) after insertions:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_indexed_heap_d_max_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDEXED_HEAP_D_MAX_TEST tests R8VEC_INDEXED_HEAP_D_MAX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  cout << "\n";
  cout << "R8VEC_INDEXED_HEAP_D_MAX_TEST\n";
  cout << "  For an indexed R8VEC,\n";
  cout << "  R8VEC_INDEXED_HEAP_D_MAX reports the maximum value.\n";
//
//  Set the data array.  To keep things easy, we will use the indicator vector.
//
  a = r8vec_indicator1_new ( m );
//
//  The index array will initially be a random subset of the numbers 1 to M,
//  in random order.
//
  n = 5;
  indx[0]  =  8;
  indx[1]  =  1;
  indx[2]  =  7;
  indx[3]  = 13;
  indx[4]  =  4;
  indx[5]  =  6;
  indx[6]  = 14;
  indx[7]  =  0;
  indx[8]  = 18;
  indx[9]  = 19;
  indx[10] =  2;

  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  cout << "\n";
  cout << "  A(INDX):\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Create a descending heap from the indexed array.
//
  r8vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  cout << "\n";
  cout << "  A(INDX) after heaping:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Insert five entries, and monitor the maximum.
//
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    cout << "\n";
    cout << "  Inserting value " << a[indx_insert] << "\n";

    r8vec_indexed_heap_d_insert ( n, a, indx, indx_insert );

    indx_max = r8vec_indexed_heap_d_max ( n, a, indx );

    cout << "  Current maximum is " << a[indx_max] << "\n";
  }
  r8vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  cout << "\n";
  cout << "  A(INDX) after insertions:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_indicator0_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDICATOR0_NEW_TEST tests R8VEC_INDICATOR0_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double *v;

  cout << "\n";
  cout << "TEST012555\n";
  cout << "  R8VEC_INDICATOR0_NEW returns an indicator vector.\n";

  n = 10;
  v = r8vec_indicator0_new ( n );
  r8vec_print ( n, v, "  Indicator0 vector:" );
  delete [] v;

  return;
}
//****************************************************************************80

void r8vec_legendre_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LEGENDRE_TEST tests R8VEC_LEGENDRE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double *r;
  double r1;
  double r2;

  cout << "\n";
  cout << "R8VEC_LEGENDRE_TEST\n";
  cout << "  R8VEC_LEGENDRE computes N Legendre points in [R1,R2].\n";

  r1 = -1.0;
  r2 = +1.0;
  n = 5;

  r = r8vec_legendre_new ( n, r1, r2 );

  cout << "\n";
  cout << "  N = " << n
       << "  R1 = " << r1
       << "  R2 = " << r2 << "\n";

  r8vec_print ( n, r, "  Legendre points:" );

  delete [] r;

  r1 =   0.0;
  r2 = +10.0;
  n = 7;

  r = r8vec_legendre_new ( n, r1, r2 );

  cout << "\n";
  cout << "  N = " << n
       << "  R1 = " << r1
       << "  R2 = " << r2 << "\n";

  r8vec_print ( n, r, "  Legendre points:" );

  delete [] r;

  return;
}
//****************************************************************************80

void r8vec_linspace_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW_TEST tests R8VEC_LINSPACE_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int n = 5;
  double *x;

  cout << "\n";
  cout << "R8VEC_LINSPACE_NEW_TEST\n";
  cout << "  For a R8VEC:\n";
  cout << "  R8VEC_LINSPACE_NEW: evenly spaced points between A and B;\n";

  a = 10.0;
  b = 20.0;

  x = r8vec_linspace_new ( n, a, b );
  r8vec_print ( n, x, "  r8vec_linspace ( 5, 10, 20 )" );
  delete [] x;

  return;
}
//****************************************************************************80

void r8vec_max_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX_TEST tests R8VEC_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int n;
  double rmax;
  int seed;

  cout << "\n";
  cout << "R8VEC_MAX_TEST\n";
  cout << "  R8VEC_MAX produces the maximum entry in a real array.\n";

  n = 10;
  seed = 123456789;

  a = r8vec_uniform_01_new ( n, seed );

  r8vec_print ( n, a, "  The array:" );

  rmax = r8vec_max ( n, a );

  cout << "\n";
  cout << "  R8VEC_MAX reports the maximum value is " << rmax << ".\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_max_index_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX_INDEX_TEST tests R8VEC_MAX_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int ival;
  int seed;

  cout << "\n";
  cout << "R8VEC_MAX_INDEX_TEST\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_MAX_INDEX: index of maximum entry;\n";

  b = - ( double ) ( N );
  c =   ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  ival = r8vec_max_index ( N, a );
  cout << "  Maximum index:           " << ival << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_mean_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEAN_TEST tests R8VEC_MEAN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double mean;
  int n = 10;
  double r8_hi;
  double r8_lo;
  int seed;

  cout << "\n";
  cout << "R8VEC_MEAN_TEST\n";
  cout << "  R8VEC_MEAN computes the mean of an R8VEC;\n";

  r8_lo = - 5.0;
  r8_hi = + 5.0;
  seed = 123456789;
  a = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );

  r8vec_print ( n, a, "  Input vector:" );

  mean = r8vec_mean ( n, a );

  cout << "\n";
  cout << "  Mean:    " << mean << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_median_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEDIAN_TEST tests R8VEC_MEDIAN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double b;
  double c;
  double median;
  int seed;

  cout << "\n";
  cout << "R8VEC_MEDIAN_TEST\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_MEDIAN:    median value;\n";

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  median = r8vec_median ( N, a );
  cout << "  Median:  " << median << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_midspace_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIDSPACE_NEW_TEST tests R8VEC_MIDSPACE_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int n = 5;
  double *x;

  cout << "\n";
  cout << "R8VEC_MIDSPACE_NEW_TEST\n";
  cout << "  For a R8VEC:\n";
  cout << "  R8VEC_MIDSPACE_NEW: evenly spaced midpoints between A and B\n";

  a = 10.0;
  b = 20.0;

  x = r8vec_midspace_new ( n, a, b );
  r8vec_print ( n, x, "  r8vec_midspace ( 5, 10, 20 )" );
  delete [] x;

  return;
}
//****************************************************************************80

void r8vec_min_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN_TEST tests R8VEC_MIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int n;
  double rmin;
  int seed;

  cout << "\n";
  cout << "R8VEC_MIN_TEST\n";
  cout << "  R8VEC_MIN produces the minimum entry.\n";

  n = 10;
  seed = 123456789;

  a = r8vec_uniform_01_new ( n, seed );

  r8vec_print ( n, a, "  The array:" );

  rmin = r8vec_min ( n, a );

  cout << "\n";
  cout << "  R8VEC_MIN reports the minimum value is " << rmin << ".\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_min_index_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN_INDEX_TEST tests R8VEC_MIN_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int ival;
  int seed;

  cout << "\n";
  cout << "R8VEC_MIN_INDEX_TEST\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_MIN_INDEX: index of minimum entry;\n";

  b = - ( double ) ( N );
  c =   ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  ival = r8vec_min_index ( N, a );
  cout << "  Minimum index:           " << ival << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_nint_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NINT_TEST tests R8VEC_NINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int n;
  int seed;
  double x1;
  double x2;

  cout << "\n";
  cout << "R8VEC_NINT_TEST\n";
  cout << "  R8VEC_NINT rounds an R8VEC.\n";

  n = 5;
  x1 = -5.0;
  x2 = +5.0;
  seed = 123456789;
  a = r8vec_uniform_ab_new ( n, x1, x2, seed );
  r8vec_print ( n, a, "  Vector A:" );
  r8vec_nint ( n, a );
  r8vec_print ( n, a, "  Rounded vector A:" );

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_norm_l0_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_L0_TEST tests R8VEC_NORM_L0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double a_hi;
  double a_lo;
  int n;
  int seed;

  cout << "\n";
  cout << "R8VEC_NORM_L0_TEST\n";
  cout << "  R8VEC_NORM_L0 computes the L0 'norm' of an R8VEC.\n";

  n = 10;
  a_lo = - 2.0;
  a_hi = + 2.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, a_lo, a_hi, seed );

  r8vec_nint ( n, a );

  r8vec_print ( n, a, "  Input vector:" );

  cout << "\n";
  cout << "  L0 norm:           " << r8vec_norm_l0 ( n, a ) << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_norm_l1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_L1_TEST tests R8VEC_NORM_L1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_NORM_L1_TEST\n";
  cout << "  R8VEC_NORM_L1 computes the L1 norm of an R8VEC.\n";

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";
  cout << "  L1 norm:           " << r8vec_norm_l1 ( N, a ) << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_norm_l2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_L2_TEST tests R8VEC_NORM_L2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_NORM_L2_TEST\n";
  cout << "  R8VEC_NORM_L2 computes the L2 norm of an R8VEC.\n";

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";
  cout << "  L2 norm:           " << r8vec_norm_l2 ( N, a ) << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_norm_li_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_LI_TEST tests R8VEC_NORM_LI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_NORM_LI_TEST\n";
  cout << "  R8VEC_NORM_LI computes the Loo norm of an R8VEC.\n";

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";
  cout << "  L-Infinity norm:   " << r8vec_norm_li ( N, a ) << "\n";;

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_normal_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_TEST tests R8VEC_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 1000

  int i;
  int n;
  int seed = 123456789;
  double *x;
  double x_max;
  double x_mean;
  double x_min;
  double x_var;

  cout << "\n";
  cout << "R8VEC_NORMAL_01_TEST\n";
  cout << "  R8VEC_NORMAL_01 computes a vector of normally\n";
  cout << "  distributed random numbers.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
//
//  Test 1:
//  Simply call 5 times for 1 value, and print.
//
  cout << "\n";
  cout << "  Test #1: Call 5 times, 1 value each time.\n";
  cout << "\n";

  n = 1;
  for ( i = 0; i < 5; i++ )
  {
    x = r8vec_normal_01_new ( n, seed );
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[0] << "\n";
    delete [] x;
  }
//
//  Test 2:
//  Restore the random number seed, and repeat.
//
  cout << "\n";
  cout << "  Test #2: Restore the random number seed.\n";
  cout << "  Call 5 times, 1 value each time.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  seed = 123456789;

  n = 1;
  for ( i = 0; i < 5; i++ )
  {
    x = r8vec_normal_01_new ( n, seed );
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[0] << "\n";
    delete [] x;
  }
//
//  Test 3:
//  Restore the random number seed, compute all 5 values at once.
//
  cout << "\n";
  cout << "  Test #3: Restore the random number seed.\n";
  cout << "  Call 1 time for 5 values.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  seed = 123456789;

  n = 5;
  x = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  delete [] x;
//
//  Test 4:
//  Restore the random number seed, compute all 5 values at once.
//
  cout << "\n";
  cout << "  Test #4: Restore the random number seed.\n";
  cout << "  Call for 2, 1, and 2 values.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  seed = 123456789;

  n = 2;
  x = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  delete [] x;

  n = 1;
  x = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  delete [] x;

  n = 2;
  x = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  delete [] x;
//
//  Test 5:
//  Determine the minimum, maximum, mean and variance.
//
  n = N_MAX;
  x = r8vec_normal_01_new ( n, seed );
  x_min = r8vec_min ( n, x );
  x_max = r8vec_max ( n, x );
  x_mean = r8vec_mean ( n, x );
  x_var = r8vec_variance ( n, x );
  delete [] x;

  cout << "\n";
  cout << "  Test #5:\n";
  cout << "  Number of samples was " << n << "\n";
  cout << "  Minimum value was " << x_min << "\n";
  cout << "  Maximum value was " << x_max << "\n";
  cout << "  Average value was " << x_mean << "\n";
  cout << "  Variance was      " << x_var << "\n";
  cout << "  Expected average  = 0.0\n";
  cout << "  Expected variance = 1.0\n";

  return;
# undef N_MAX
}
//****************************************************************************80

void r8vec_normalize_l1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMALIZE_L1_TEST tests R8VEC_NORMALIZE_L1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_NORMALIZE_L1_TEST\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_NORMALIZE_L1:  make unit sum;\n";

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  r8vec_normalize_l1 ( N, a );

  r8vec_print ( N, a, "  After calling R8VEC_NORMALIZE_L1:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_order_type_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ORDER_TYPE_TEST tests R8VEC_ORDER_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define TEST_NUM 6

  int itest;
  int j;
  int order;
  double x[N];

  cout << "\n";
  cout << "R8VEC_ORDER_TYPE_TEST\n";
  cout << "  R8VEC_ORDER_TYPE classifies a real vector as\n";
  cout << "  -1: no order\n";
  cout << "   0: all equal;\n";
  cout << "   1: ascending;\n";
  cout << "   2: strictly ascending;\n";
  cout << "   3: descending;\n";
  cout << "   4: strictly descending.\n";
  cout << "\n";

  for ( itest = 1; itest <= TEST_NUM; itest++ )
  {
    if ( itest == 1 )
    {
      x[0] = 1.0;
      x[1] = 3.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 2 )
    {
      x[0] = 2.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 2.0;
    }
    else if ( itest == 3 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 4 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 3.0;
      x[3] = 4.0;
    }
    else if ( itest == 5 )
    {
      x[0] = 4.0;
      x[1] = 4.0;
      x[2] = 3.0;
      x[3] = 1.0;
    }
    else if ( itest == 6 )
    {
      x[0] = 9.0;
      x[1] = 7.0;
      x[2] = 3.0;
      x[3] = 0.0;
    }

    order = r8vec_order_type ( N, x );

    cout << "\n";
    cout << "The following vector has order type " << order << ".\n";
    cout << "\n";

    r8vec_print ( N, x, "" );
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void r8vec_permute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PERMUTE_TEST tests R8VEC_PERMUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int p[5] = { 1, 3, 4, 0, 2 };
  double x[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

  cout << "\n";
  cout << "R8VEC_PERMUTE_TEST\n";
  cout << "  R8VEC_PERMUTE permutes an R8VEC in place.\n";

  r8vec_print ( n, x, "  Original Array X[]:" );

  i4vec_print ( n, p, "  Permutation Vector P[]:" );

  r8vec_permute ( n, p, x );

  r8vec_print ( n, x, "  Permuted array X[P[]]:" );

  return;
}
//****************************************************************************80

void r8vec_permute_uniform_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PERMUTE_UNIFORM_TEST tests R8VEC_PERMUTE_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int n = 10;
  int seed;

  cout << "\n";
  cout << "R8VEC_PERMUTE_UNIFORM_TEST\n";
  cout << "  R8VEC_PERMUTE_UNIFORM randomly reorders an R8VEC.\n";

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( double ) ( 101 + i );
  }
  seed = 123456789;

  r8vec_print ( n, a, "  A, before rearrangement:" );

  r8vec_permute_uniform ( n, a, seed );

  r8vec_print ( n, a, "  A, after random permutation:" );

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_polarize_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_POLARIZE_TEST tests R8VEC_POLARIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N] = { 1.0, 2.0,  3.0 };
  double a2[N];
  double a_normal[N];
  double a_parallel[N];
  double ap_norm;
  int i;
  double p[N] = { 3.0, 1.0, -2.0 };
  double p_norm;
  double pan;
  double pap;

  cout << "\n";
  cout << "R8VEC_POLARIZE_TEST\n";
  cout << "  R8VEC_POLARIZE decomposes a vector into\n";
  cout << "  components parallel and normal to a direction.\n";

  r8vec_print ( N, a, "  Original vector:" );

  r8vec_print ( N, p, "  Direction vector:" );

  r8vec_polarize ( N, a, p, a_normal, a_parallel );

  r8vec_print ( N, a_normal, "  Normal component:" );

  r8vec_print ( N, a_parallel, "  Parallel component:" );

  pan = r8vec_dot_product ( N, p, a_normal );
  p_norm = r8vec_norm ( N, p );
  ap_norm = r8vec_norm ( N, a_parallel );

  pap = r8vec_dot_product ( N, p, a_parallel ) / ( p_norm * ap_norm );

  cout << "\n";
  cout << "  Dot product of P and A_normal (should be 0) " << pan << "\n";
  cout << "  Cosine of angle between P and A_parallel (should be 1 or -1) "
       << pap << "\n";

  for ( i = 0; i < N; i++ )
  {
    a2[i] = a_normal[i] + a_parallel[i];
  }

  r8vec_print ( N, a2, "  Sum of components (should equal A):" );

  return;
# undef N
}
//****************************************************************************80

void r8vec_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_TEST tests R8VEC_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[4] = { 123.456, 0.000005, -1.0E+06, 3.14159265 };
  int n = 4;

  cout << "\n";
  cout << "TEST1335\n";
  cout << "  R8VEC_PRINT prints an R8VEC.\n";

  r8vec_print ( n, a, "  The R8VEC:" );

  return;
}
//****************************************************************************80

void r8vec_rotate_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ROTATE_TEST tests R8VEC_ROTATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  int m = 2;

  cout << "\n";
  cout << "R8VEC_ROTATE_TEST\n";
  cout << "  R8VEC_ROTATE rotates an R8VEC in place.\n";
  cout << "\n";
  cout << "  Rotate entries " << m << " places to the right.\n";

  r8vec_print ( N, a, "  Original array:" );

  r8vec_rotate ( N, a, m );

  r8vec_print ( N, a, "  Rotated array:" );

  return;
# undef N
}
//****************************************************************************80

void r8vec_reverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_REVERSE_TEST tests R8VEC_REVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;

  cout << "\n";
  cout << "R8VEC_REVERSE_TEST\n";
  cout << "  R8VEC_REVERSE reverses an R8VEC.\n";

  a = r8vec_indicator1_new ( N );

  r8vec_print ( N, a, "  Original array:" );

  r8vec_reverse ( N, a );

  r8vec_print ( N, a, "  Reversed array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_search_binary_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SEARCH_BINARY_A_TEST tests R8VEC_SEARCH_BINARY_A;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  int index;
  int nc;
  double search_val;
  int seed = 123456789;

  cout << "\n";
  cout << "R8VEC_SEARCH_BINARY_A_TEST\n";
  cout << "  For ascending order:\n";
  cout << "  R8VEC_SEARCH_BINARY_A searches a sorted array;\n";
  cout << "\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8vec_uniform_01_new ( N, seed );

  search_val = a[0];

  r8vec_sort_heap_a ( N, a );

  r8vec_print ( N, a, "  Sorted vector A:" );
//
//  Now search the sorted array for a given value.
//
  cout << "\n";
  cout << "  Search the array for the value " << search_val << "\n";

  index = r8vec_search_binary_a ( N, a, search_val );

  cout << "\n";
  cout << "  SEARCH RESULT:\n";
  cout << "\n";

  if ( 0 < index )
  {
    cout << "    The value occurs in index " << index << "\n";
  }
  else
  {
    cout << "    The value does not occur in the array.\n";
  }

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sort_bubble_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_BUBBLE_A_TEST tests R8VEC_SORT_BUBBLE_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  int i;
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_BUBBLE_A_TEST\n";
  cout << "  R8VEC_SORT_BUBBLE_A sorts a real array.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = r8vec_uniform_01_new ( N, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_bubble_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sort_heap_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_A_TEST tests R8VEC_SORT_HEAP_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_HEAP_A_TEST\n";
  cout << "  R8VEC_SORT_HEAP_A ascending sorts an R8VEC.\n";
  cout << "\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print_some ( N, a, 1, 10, "  Original array:" );

  r8vec_sort_heap_a ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Ascending sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sort_heap_d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_D_TEST tests R8VEC_SORT_HEAP_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_HEAP_D_TEST\n";
  cout << "  R8VEC_SORT_HEAP_D descending sorts an R8VEC.\n";
  cout << "\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print_some ( N, a, 1, 10, "  Original array:" );

  r8vec_sort_heap_d ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Descending sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sort_heap_index_a_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_INDEX_A_NEW_TEST tests R8VEC_SORT_HEAP_INDEX_A_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int i;
  int *indx;
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_HEAP_INDEX_A_NEW_TEST\n";
  cout << "  R8VEC_SORT_HEAP_INDEX_A_NEW creates an ascending\n";
  cout << "  sort index for an R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  indx = r8vec_sort_heap_index_a_new ( N, a );

  cout << "\n";
  cout << "  After indexed ascending sort:\n";
  cout << "\n";
  cout << "         I INDX(I)       A(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << indx[i]
         << "  " << setw(12) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  Now use the index array to carry out the\n";
  cout << "  permutation implicitly.\n";
  cout << "\n";
  cout << "   INDX(I)  A(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << indx[i]
         << "  " << setw(12) << a[indx[i]] << "\n";
  }
  cout << "\n";
  cout << "  Call R8VEC_PERMUTE to carry out the permutation explicitly.\n";
  cout << "\n";

  r8vec_permute ( N, indx, a );

  r8vec_print ( N, a, "  I, A(I)" );

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sort_heap_index_d_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_INDEX_D_NEW_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int i;
  int *indx;
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_HEAP_INDEX_D_NEW_TEST\n";
  cout << "  R8VEC_SORT_HEAP_INDEX_D_NEW creates a descending\n";
  cout << "  sort index for an R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  indx = r8vec_sort_heap_index_d_new ( N, a );

  cout << "\n";
  cout << "  After indexed descending sort:\n";
  cout << "\n";
  cout << "         I INDX(I)       A(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << indx[i]
         << "  " << setw(12) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  Now use the index array to carry out the\n";
  cout << "  permutation implicitly.\n";
  cout << "\n";
  cout << "   INDX(I)  A(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << indx[i]
         << "  " << setw(12) << a[indx[i]] << "\n";
  }
  cout << "\n";
  cout << "  Call R8VEC_PERMUTE to carry out the permutation explicitly.\n";
  cout << "\n";

  r8vec_permute ( N, indx, a );

  r8vec_print ( N, a, "  I, A(I)" );

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sort_heap_mask_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_MASK_A_TEST tests R8VEC_SORT_HEAP_MASK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define MASK_NUM 10
# define N 20

  double *a;
  double b;
  double c;
  int i;
  int *indx;
  int mask[MASK_NUM] = { 1, 3, 6, 7, 8, 11, 12, 15, 17, 18 };
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_HEAP_MASK_A_TEST\n";
  cout << "  R8VEC_SORT_HEAP_MASK_A creates an ascending\n";
  cout << "  sort index for a masked R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  i4vec_print ( MASK_NUM, mask, "  The mask array:" );

  r8vec_mask_print ( N, a, MASK_NUM, mask, "  The masked unsorted array:" );

  indx = r8vec_sort_heap_mask_a ( N, a, MASK_NUM, mask );

  cout << "\n";
  cout << "  After masked indexed ascending sort:\n";
  cout << "\n";
  cout << "  I, INDX(I), MASK(INDX(I)), A(MASK(INDX(I)))\n";
  cout << "\n";
  for ( i = 0; i < MASK_NUM; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << indx[i]
         << "  " << setw(6) << mask[indx[i]]
         << "  " << setw(14) << a[mask[indx[i]]] << "\n";
  }

  cout << "\n";
  cout << "  Call I4VEC_PERMUTE to carry out the index permutation\n";
  cout << "  explicitly on the MASK vector.\n";
  cout << "\n";

  i4vec_permute ( MASK_NUM, indx, mask );
//
//  Essentially, INDX becomes the identity vector now.
//
  delete [] indx;

  indx = i4vec_indicator1_new ( MASK_NUM );

  i4vec_print ( MASK_NUM, mask, "  The reordered mask array:" );

  r8vec_mask_print ( N, a, MASK_NUM, mask,
    "  The reordered masked sorted array:" );

  delete [] a;
  delete [] indx;

  return;
# undef MASK_NUM
# undef N
}
//****************************************************************************80

void r8vec_sort_insert_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_INSERT_A_TEST tests R8VEC_SORT_INSERT_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_INSERT_A_TEST\n";
  cout << "  R8VEC_SORT_INSERT_A ascending sorts an R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print_some ( N, a, 1, 10, "  Unsorted array:" );

  r8vec_sort_insert_a ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sort_insert_index_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_INSERT_INDEX_A_TEST tests R8VEC_SORT_INSERT_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int i;
  int *indx;
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_INSERT_INDEX_A_TEST\n";
  cout << "  R8VEC_SORT_INSERT_INDEX_A creates an ascending\n";
  cout << "  sort index for an R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print_some ( N, a, 1, 10, "  Unsorted array:" );

  indx = r8vec_sort_insert_index_a ( N, a );

  cout << "\n";
  cout << "  After indexed ascending sort:\n";
  cout << "\n";
  cout << "       I  INDX(I), A(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << indx[i]
         << "  " << setw(12) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  Now use the index array to carry out the\n";
  cout << "  permutation implicitly.\n";
  cout << "\n";
  cout << "       I  INDX(I), A(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << indx[i]
         << "  " << setw(12) << a[indx[i]] << "\n";
  }

  cout << "\n";
  cout << "  Call R8VEC_PERMUTE to carry out the permutation explicitly.\n";
  cout << "\n";

  r8vec_permute ( N, indx, a );

  r8vec_print_some ( N, a, 1, 10, "  Permuted data" );

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sort_quick_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_QUICK_A_TEST tests R8VEC_SORT_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "R8VEC_SORT_QUICK_A_TEST\n";
  cout << "  R8VEC_SORT_QUICK_A sorts an R8VEC\n";
  cout << "  using quick sort.\n";

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_quick_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sorted_merge_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTEC_MERGE_A_TEST tests R8VEC_SORTED_MERGE_A;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double *c;
  int na = 10;
  int nb = 10;
  int nc;
  int seed = 123456789;

  cout << "\n";
  cout << "R8VEC_SORTED_MERGE_A_TEST\n";
  cout << "  For ascending order:\n";
  cout << "  R8VEC_SORTED_MERGE_A merges two sorted R8VEC's;\n";
  cout << "\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8vec_uniform_01_new ( na, seed );
  b = r8vec_uniform_01_new ( nb, seed );

  r8vec_sort_heap_a ( na, a );

  r8vec_sort_heap_a ( nb, b );

  r8vec_print ( na, a, "  Sorted vector A:" );

  r8vec_print ( nb, b, "  Sorted vector B:" );

  c = r8vec_sorted_merge_a ( na, a, nb, b, nc );

  r8vec_print ( nc, c, "  Merged vector C:" );

  delete [] a;
  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void r8vec_sorted_nearest_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTED_NEAREST_TEST tests R8VEC_SORTED_NEAREST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double b;
  double c;
  int i;
  int j;
  int seed = 123456789;
  double *x;
  double xval;

  cout << "\n";
  cout << "R8VEC_SORTED_NEAREST_TEST\n";
  cout << "  R8VEC_SORTED_NEAREST finds the nearest entry\n";
  cout << "  in a sorted real array.\n";

  b = 0.0;
  c = 10.0;

  x = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_sort_heap_a ( N, x );

  r8vec_print ( N, x, "  Sorted array:" );

  cout << "\n";
  cout << "     Test        Nearest\n";
  cout << "     Value    Index   Value\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    xval = r8_uniform_01 ( seed );

    j = r8vec_sorted_nearest ( N, x, xval );

    cout << "  "
         << setw(8) << xval   << "    "
         << setw(6) << j      << "  "
         << setw(8) << x[j-1] << "\n";
  }

  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sorted_range_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTED_RANGE_TEST tests R8VEC_SORTED_RANGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_hi;
  int i_lo;
  int n = 10;
  double r[10];
  double r_lo;
  double r_hi;
  int seed;
  double t;
  int test;

  cout << "\n";
  cout << "R8VEC_SORTED_RANGE_TEST\n";
  cout << "  R8VEC_SORTED_RANGE seeks the range of indices\n";
  cout << "  in a sorted vector R so that\n";
  cout << "  R_LO <= R(I_LO:I_HI) <= R_HI.\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r8vec_uniform_01 ( n, seed, r );

    r8vec_sort_heap_a ( n, r );

    r8vec_print ( n, r, "  Sorted array R:" );

    r_lo = r8_uniform_01 ( seed );
    r_hi = r8_uniform_01 ( seed );

    if ( r_hi < r_lo )
    {
      t = r_lo;
      r_lo = r_hi;
      r_hi = t;
    }

    r8vec_sorted_range ( n, r, r_lo, r_hi, i_lo, i_hi );

    cout << "\n";
    if ( i_hi < i_lo )
    {
      cout << "  R_LO  " << setw(14) << r_lo << "\n";
      cout << "  R_HI  " << setw(14) << r_hi << "\n";
      cout << "  Empty range in R.\n";
    }
    else
    {
      cout << "  R_LO  " << setw(14) << r_lo << "\n";
      for ( i = i_lo; i <= i_hi; i++)
      {
        cout << "  " << setw(4) << i << "  " << setw(14) << r[i] << "\n";
      }
      cout << "  R_HI  " << setw(14) << r_hi << "\n";
    }
  }
  return;
}
//****************************************************************************80

void r8vec_sorted_split_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTED_SPLIT_TEST tests R8VEC_SORTED_SPLIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b;
  double c;
  int i;
  int i_gt;
  int i_lt;
  int isplit;
  int n = 25;
  int seed;
  double split;

  cout << "\n";
  cout << "R8VEC_SORTED_SPLIT_TEST\n";
  cout << "  R8VEC_SORTED_SPLIT splits a sorted vector into\n";
  cout << "  entries less than and greater than a\n";
  cout << "  splitting value.\n";

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, b, c, seed );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.5 * ( double ) ( r8_nint ( a[i] ) );
  }

  r8vec_sort_heap_a ( n, a );

  split = 0.5 * ( a[0] + a[n-1] );

  r8vec_print ( n, a, "  The sorted array:" );

  cout << "\n";
  cout << "  Splitting value is " << split << "\n";
  cout << "\n";

  r8vec_sorted_split ( n, a, split, i_lt, i_gt );

  cout << "  Lower index I_LT = " << i_lt << "\n";
  cout << "  Upper index I_GT = " << i_gt << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_sorted_undex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTED_UNDEX_TEST tests R8VEC_SORTED_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2008
//
//  Author:
//
//    John Burkardt
//
{
# define X_NUM 9

  int i;
  double tol;
  int *undx;
  int x_num = X_NUM;
  int x_unique_num;
  double x_val[X_NUM] = { 11.0, 11.0, 11.0, 22.0, 22.0, 33.0, 33.0, 55.0, 55.0 };
  int *xdnu;
  double *xu_val;

  cout << "\n";
  cout << "R8VEC_SORTED_UNDEX_TEST\n";
  cout << "  R8VEC_SORTED_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the unique elements of a sorted R8VEC,\n";
  cout << "  and a map from the original vector to the (implicit)\n";
  cout << "  vector of sorted unique elements.\n";

  r8vec_print ( x_num, x_val, "  The vector X:" );

  tol = r8_epsilon ( );
  x_unique_num = r8vec_sorted_unique_count ( x_num, x_val, tol );

  undx = new int[x_unique_num];
  xu_val = new double[x_unique_num];

  xdnu = new int[x_num];

  cout << "\n";
  cout << "  Number of unique entries in X is " << x_unique_num << "\n";

  r8vec_sorted_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu );

  cout << "\n";
  cout << "  UNDX can be used to list the unique elements of X\n";
  cout << "  in sorted order.\n";
  cout << "\n";
  cout << "     I  UNDX   X(UNDX)\n";
  cout << "\n";

  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(8) << x_val[undx[i]] << "\n";
  }

  for ( i = 0; i < x_unique_num; i++ )
  {
    xu_val[i] = x_val[undx[i]];
  }
  cout << "\n";
  cout << "  UNDX can be used to created XU, a copy of X\n";
  cout << "  containing only the unique elements, in sorted order.\n";
  cout << "\n";
  cout << "     I  UNDX     XU(I)\n";
  cout << "\n";
  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(8) << xu_val[i] << "\n";
  }

  cout << "\n";
  cout << "  XDNU can be used to match each element of X with one of the\n";
  cout << "  unique elements\n";
  cout << "\n";
  cout << "     I  XDNU  X(I)   XU(XDNU(I))\n";
  cout << "\n";

  for ( i = 0; i < x_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << x_val[i]
         << "  " << setw(12) << xu_val[xdnu[i]] << "\n";
  }

  delete [] undx;
  delete [] xdnu;
  delete [] xu_val;

  return;
# undef X_NUM
}
//****************************************************************************80

void r8vec_sorted_unique_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTED_UNIQUE_TEST tests R8VEC_SORTED_UNIQUE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_unique;
  double b;
  double c;
  int i;
  int n = 20;
  int seed;
  double tol = 0.25;
  int unique_num;

  cout << "\n";
  cout << "R8VEC_SORTED_UNIQUE_TEST\n";
  cout << "  R8VEC_SORTED_UNIQUE finds unique entries in a sorted R8VEC;\n";

  b = 0.0;
  c = ( double ) ( n );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, b, c, seed );

  r8vec_nint ( n, a );

  r8vec_print_some ( n, a, 1, 10, "  Unsorted array:" );

  r8vec_sort_heap_a ( n, a );

  a_unique = r8vec_sorted_unique ( n, a, tol, unique_num );

  r8vec_print ( unique_num, a_unique, "  Unique entries" );

  delete [] a;
  delete [] a_unique;

  return;
}
//****************************************************************************80

void r8vec_sorted_unique_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTED_UNIQUE_COUNT_TEST tests R8VEC_SORTED_UNIQUE_COUNT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 30

  double *a;
  double b;
  double c;
  int i;
  int unique_num;
  int seed;
  double tol = 0.25;

  cout << "\n";
  cout << "R8VEC_SORTED_UNIQUE_COUNT_TEST\n";
  cout << "  R8VEC_SORTED_UNIQUE_COUNT counts unique entries in a sorted R8VEC;\n";

  b = 0.0;
  c = ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  for ( i = 0; i < N; i++ )
  {
    a[i] = ( double ) r8_nint ( a[i] );
  }

  unique_num = r8vec_sorted_unique_count ( N, a, tol );

  cout << "\n";
  cout << "  Using a tolerance of " << tol << "\n";
  cout << "  R8VEC_SORTED_UNIQUE_COUNT counts " << unique_num
       << " unique entries in A.\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void r8vec_sorted_unique_hist_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTED_UNIQUE_HIST_TEST tests R8VEC_SORTED_UNIQUE_HIST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define MAXUNIQ 30
# define N 30

  double *a;
  int acount[MAXUNIQ];
  double auniq[MAXUNIQ];
  double b;
  double c;
  int i;
  int unique_num;
  int seed;
  double tol = 0.25;

  cout << "\n";
  cout << "R8VEC_SORTED_UNIQUE_HIST_TEST\n";
  cout << "  R8VEC_SORTED_UNIQUE_HIST stores the unique entries\n";
  cout << "  and their multiplicities.\n";

  b = 0.0;
  c = ( double ) N;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  for ( i = 0; i < N; i++ )
  {
    a[i] = ( ( double ) ( ( int ) a[i] ) ) + 0.5;
  }

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_bubble_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  r8vec_sorted_unique_hist ( N, a, tol, MAXUNIQ, unique_num, auniq, acount );

  cout << "\n";
  cout << "  R8VEC_SORTED_UNIQUE_HIST counts " << unique_num << " unique entries.\n";
  cout << "\n";
  cout << "  Value  Multiplicity\n";
  cout << "\n";
  for ( i = 0; i < unique_num; i++ )
  {
    cout << setw(6)  << i         << "  "
         << setw(12) << auniq[i]  << "  "
         << setw(6)  << acount[i] << "\n";
  }

  delete [] a;

  return;
# undef MAXUNIQ
# undef N
}
//****************************************************************************80

void r8vec_split_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SPLIT_TEST tests R8VEC_SPLIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b;
  double c;
  int i;
  int i_gt;
  int i_lt;
  int isplit;
  int n = 25;
  int seed;
  double split;

  cout << "\n";
  cout << "R8VEC_SPLIT_TEST\n";
  cout << "  R8VEC_SPLIT splits a vector into\n";
  cout << "  entries less than and greater than a\n";
  cout << "  splitting value.\n";

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, b, c, seed );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.5 * ( double ) ( r8_nint ( a[i] ) );
  }

  split = 0.5 * ( a[0] + a[n-1] );

  cout << "\n";
  cout << "  Splitting value is " << split << "\n";
  cout << "\n";

  isplit = r8vec_split ( n, a, split );

  r8vec_print ( n, a, "  The split array:" );

  cout << "\n";
  cout << "  Array entries <= SPLIT up to index " << isplit << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec_transpose_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TRANSPOSE_PRINT_TEST tests R8VEC_TRANSPOSE_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  int n = 12;
  int seed;
  double *x;

  seed = 123456789;

  cout << "\n";
  cout << "R8VEC_TRANSPOSE_PRINT_TEST\n";
  cout << "  R8VEC_TRANSPOSE_PRINT prints an R8VEC \"transposed\",\n";
  cout << "  that is, placing multiple entries on a line.\n";

  x = r8vec_uniform_01_new ( n, seed );

  r8vec_transpose_print ( n, x, "  The vector X:" );

  delete [] x;

  return;
}
//****************************************************************************80

void r8vec_undex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNDEX_TEST tests R8VEC_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define X_NUM 9

  int i;
  double tol;
  int *undx;
  int x_num = X_NUM;
  int x_unique_num;
  double x_val[X_NUM] = { 33.0, 55.0, 11.0, 11.0, 55.0, 33.0, 22.0, 22.0, 11.0 };
  int *xdnu;
  double *xu_val;

  cout << "\n";
  cout << "R8VEC_UNDEX_TEST\n";
  cout << "  R8VEC_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the unique elements of an (unsorted) R8VEC,\n";
  cout << "  and a map from the original vector to the (implicit)\n";
  cout << "  vector of sorted unique elements.\n";

  r8vec_print ( x_num, x_val, "  The vector X:" );

  tol = r8_epsilon ( );
  x_unique_num = r8vec_unique_count ( x_num, x_val, tol );

  undx = new int[x_unique_num];
  xu_val = new double[x_unique_num];

  xdnu = new int[x_num];

  cout << "\n";
  cout << "  Number of unique entries in X is " << x_unique_num << "\n";

  r8vec_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu );

  cout << "\n";
  cout << "  UNDX can be used to list the unique elements of X\n";
  cout << "  in sorted order.\n";
  cout << "\n";
  cout << "     I  UNDX   X(UNDX)\n";
  cout << "\n";

  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(8) << x_val[undx[i]] << "\n";
  }

  for ( i = 0; i < x_unique_num; i++ )
  {
    xu_val[i] = x_val[undx[i]];
  }
  cout << "\n";
  cout << "  UNDX can be used to created XU, a copy of X\n";
  cout << "  containing only the unique elements, in sorted order.\n";
  cout << "\n";
  cout << "     I  UNDX     XU(I)\n";
  cout << "\n";
  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(8) << xu_val[i] << "\n";
  }

  cout << "\n";
  cout << "  XDNU can be used to match each element of X with one of the\n";
  cout << "  unique elements\n";
  cout << "\n";
  cout << "     I  XDNU  X(I)   XU(XDNU(I))\n";
  cout << "\n";

  for ( i = 0; i < x_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << x_val[i]
         << "  " << setw(12) << xu_val[xdnu[i]] << "\n";
  }

  delete [] undx;
  delete [] xdnu;
  delete [] xu_val;

  return;
# undef X_NUM
}
//****************************************************************************80

void r8vec_uniform_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW_TEST tests R8VEC_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int j;
  double *r;
  int seed;

  cout << "\n";
  cout << "R8VEC_UNIFORM_01_NEW_TEST\n";
  cout << "  R8VEC_UNIFORM returns a random R8VEC\n";
  cout << "  with entries in [ 0.0, 1.0 ]\n";
  cout << "\n";

  seed = 123456789;

  for ( j = 1; j <= 3; j++ )
  {
    cout << "\n";
    cout << "  Input SEED = " << seed << "\n";
    cout << "\n";

    r = r8vec_uniform_01_new ( N, seed );

    r8vec_print ( N, r, "  Random R8VEC:" );

    delete [] r;
  }

  return;
# undef N
}
//****************************************************************************80

void r8vec_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_AB_NEW_TEST tests R8VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a = 10.0;
  double b = 20.0;
  int j;
  double *r;
  int seed;

  cout << "\n";
  cout << "R8VEC_UNIFORM_AB_NEW_TEST\n";
  cout << "  R8VEC_UNIFORM returns a random R8VEC\n";
  cout << "  with entries in a given range [ A, B ]\n";
  cout << "\n";
  cout << "  For this problem:\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "\n";

  seed = 123456789;

  for ( j = 1; j <= 3; j++ )
  {
    cout << "\n";
    cout << "  Input SEED = " << seed << "\n";
    cout << "\n";

    r = r8vec_uniform_ab_new ( N, a, b, seed );

    r8vec_print ( N, r, "  Random R8VEC:" );

    delete [] r;
  }

  return;
# undef N
}
//****************************************************************************80

void r8vec_variance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_VARIANCE_TEST tests R8VEC_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int n = 10;
  double r8_hi;
  double r8_lo;
  int seed;
  double variance;

  cout << "\n";
  cout << "R8VEC_VARIANCE_TEST\n";
  cout << "  R8VEC_VARIANCE computes the variance of an R8VEC;\n";

  r8_lo = - 5.0;
  r8_hi = + 5.0;
  seed = 123456789;
  a = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );

  r8vec_print ( n, a, "  Input vector:" );

  variance = r8vec_variance ( n, a );

  cout << "\n";
  cout << "  Variance:    " << variance << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void r8vec2_sort_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_SORT_A_TEST tests R8VEC2_SORT_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a1;
  double *a2;
  double b;
  double c;
  int n = 10;
  int seed;

  cout << "\n";
  cout << "R8VEC2_SORT_A_TEST\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SORT_A ascending sorts;\n";

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( n, a1, a2, "  The pair of arrays:" );

  r8vec2_sort_a ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after ascending sort:" );

  delete [] a1;
  delete [] a2;

  return;
}
//****************************************************************************80

void r8vec2_sort_d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_SORT_D_TEST tests R8VEC2_SORT_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a1;
  double *a2;
  double b;
  double c;
  int n = 10;
  int seed;

  cout << "\n";
  cout << "R8VEC2_SORT_D_TEST\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SORT_D descending sorts;\n";

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( n, a1, a2, "  The pair of arrays:" );

  r8vec2_sort_d ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after descending sort:" );

  delete [] a1;
  delete [] a2;

  return;
}
//****************************************************************************80

void r8vec2_sort_heap_index_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_SORT_HEAP_INDEX_A_TEST tests R8VEC2_SORT_HEAP_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int i;
  int *indx;
  int seed = 123456789;
  double x[N];
  double y[N];

  cout << "\n";
  cout << "R8VEC2_SORT_HEAP_INDEX_A_TEST\n";
  cout << "  R8VEC2_SORT_HEAP_INDEX_A creates a sort index\n";
  cout << "  for an (X,Y) array.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) i4_uniform_ab ( 0, N, seed ) / ( double ) N;
    y[i] = ( double ) i4_uniform_ab ( 0, N, seed ) / ( double ) N;
  }

  cout << "\n";
  cout << "  The unsorted array:\n";
  cout << "\n";
  cout << "         I  X(I), Y(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << y[i] << "\n";
  }

  indx = r8vec2_sort_heap_index_a ( N, x, y );

  cout << "\n";
  cout << "  After sorting:\n";
  cout << "\n";
  cout << "         I  INDX(I), X(I), Y(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << indx[i]
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << y[i] << "\n";
  }

  cout << "\n";
  cout << "  Now use the index array to carry out the\n";
  cout << "  permutation implicitly.\n";
  cout << "\n";
  cout << "         I  INDX(I), X(INDX(I)), Y(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << indx[i]
         << "  " << setw(12) << x[indx[i]]
         << "  " << setw(12) << y[indx[i]] << "\n";
  }

  cout << "\n";
  cout << "  R8VEC_PERMUTE carries out the permutation.\n";

  r8vec_permute ( N, indx, x );
  r8vec_permute ( N, indx, y );

  cout << "\n";
  cout << "         I X(I), Y(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << y[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void r8vec2_sorted_unique_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_SORTED_UNIQUE_TEST tests R8VEC2_SORTED_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a1;
  double *a2;
  double b;
  double c;
  int n = 10;
  int seed;
  int unique_num;

  cout << "\n";
  cout << "R8VEC2_SORTED_UNIQUE_TEST\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SORTED_UNIQUE counts unique entries.\n";

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( n, a1, a2, "  The pair of arrays:" );

  r8vec2_sort_a ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after ascending sort:" );

  r8vec2_sorted_unique ( n, a1, a2, unique_num );

  r8vec2_print ( unique_num, a1, a2, "  UNIQed array:" );

  delete [] a1;
  delete [] a2;

  return;
}
//****************************************************************************80

void r8vec2_sorted_unique_index_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_SORTED_UNIQUE_INDEX_TEST tests R8VEC2_SORTED_UNIQUE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a1;
  double *a2;
  double b;
  double c;
  int indx[N];
  int seed;
  int unique_num;

  cout << "\n";
  cout << "R8VEC2_SORTED_UNIQUE_INDEX_TEST\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SORTED_UNIQUE_INDEX indexes unique entries.\n";

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( N, b, c, seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( N, b, c, seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( N, a1, a2, "  The pair of arrays:" );

  r8vec2_sorted_unique_index ( N, a1, a2, unique_num, indx );

  cout << "\n";
  cout << "  The number of unique elements is " << unique_num << "\n";

  i4vec_print ( unique_num, indx, "  Index of Unique Elements:" );

  r8vec_index_order ( unique_num, a1, indx );
  r8vec_index_order ( unique_num, a2, indx );

  r8vec2_print ( unique_num, a1, a2, "  After Indexed Nonunique Deletion." );

  delete [] a1;
  delete [] a2;

  return;
# undef N
}
//****************************************************************************80

void r8vec2_sum_max_index_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_SUM_MAX_INDEX_TEST tests R8VEC2_SUM_MAX_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a1;
  double *a2;
  double b;
  double c;
  int ival;
  int n = 10;
  int seed;

  cout << "\n";
  cout << "R8VEC2_SUM_MAX_INDEX_TEST\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SUM_MAX_INDEX: index of the sum vector\n";
  cout << "  with maximum value.\n";

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, seed );

  b = 0.0;
  c = 5.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, seed );

  r8vec2_print ( n, a1, a2, "  The pair of vectors:" );

  ival = r8vec2_sum_max_index ( n, a1, a2 );

  cout << "\n";
  cout << "  Index of maximum in A+B: " << ival << "\n";

  delete [] a1;
  delete [] a2;

  return;
}
//****************************************************************************80

void roots_to_r8poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ROOTS_TO_R8POLY_TEST tests ROOTS_TO_R8POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int n = 5;
  double x[5] = { 1.0, -4.0, 3.0, 0.0, 3.0 };

  cout << "\n";
  cout << "ROOTS_TO_R8POLY_TEST:\n";
  cout << "  ROOTS_TO_R8POLY is given N real roots,\n";
  cout << "  and constructs the coefficient vector\n";
  cout << "  of the corresponding polynomial.\n";

  r8vec_print ( n, x, "  N real roots:" );

  c = roots_to_r8poly ( n, x );

  r8poly_print ( n, c, "  Corresponding polynomial:" );

  delete [] c;

  return;
}
