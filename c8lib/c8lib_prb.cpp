# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "c8lib.hpp"

int main ( );

void c8_abs_test ( );
void c8_acos_test ( );
void c8_acosh_test ( );
void c8_add_test ( );
void c8_arg_test ( );
void c8_asin_test ( );
void c8_asinh_test ( );
void c8_atan_test ( );
void c8_atanh_test ( );
void c8_conj_test ( );
void c8_cos_test ( );
void c8_cosh_test ( );
void c8_cube_root_test ( );
void c8_div_test ( );
void c8_div_r8_test ( );
void c8_exp_test ( );
void c8_imag_test ( );
void c8_i_test ( );
void c8_inv_test ( );
void c8_le_l1_test ( );
void c8_le_l2_test ( );
void c8_le_li_test ( );
void c8_log_test ( );
void c8_mag_test ( );
void c8_mul_test ( );
void c8_nint_test ( );
void c8_norm_l1_test ( );
void c8_norm_l2_test ( );
void c8_norm_li_test ( );
void c8_normal_01_test ( );
void c8_one_test ( );
void c8_print_test ( );
void c8_real_test ( );
void c8_sin_test ( );
void c8_sinh_test ( );
void c8_sqrt_test ( );
void c8_sub_test ( );
void c8_tan_test ( );
void c8_tanh_test ( );
void c8_to_cartesian_test ( );
void c8_to_polar_test ( );
void c8_uniform_01_test ( );
void c8_zero_test ( );

void c8mat_identity_new_test ( );
void c8mat_indicator_new_test ( );
void c8mat_norm_fro_test ( );
void c8mat_norm_l1_test ( );
void c8mat_norm_li_test ( );
void c8mat_uniform_01_new_test ( );

void c8vec_indicator_new_test ( );
void c8vec_nint_test ( );
void c8vec_norm_l1_test ( );
void c8vec_norm_l2_test ( );
void c8vec_norm_li_test ( );
void c8vec_print_test ( );
void c8vec_sort_a_l1_test ( );
void c8vec_sort_a_l2_test ( );
void c8vec_sort_a_li_test ( );
void c8vec_spiral_new_test ( );
void c8vec_uniform_01_new_test ( );
void c8vec_unity_new_test ( );

void cartesian_to_c8_test ( );

void polar_to_c8_test ( );

void r8_atan_test ( );
void r8_uniform_01_test ( );

void r8poly2_root_test ( );
void r8poly3_root_test ( );
void r8poly4_root_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for C8LIB_PRB.
//
//  Discussion:
//
//    C8LIB_PRB tests the C8LIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "C8LIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the C8LIB library.\n";

  c8_abs_test ( );
  c8_acos_test ( );
  c8_acosh_test ( );
  c8_add_test ( );
  c8_arg_test ( );
  c8_asin_test ( );
  c8_asinh_test ( );
  c8_atan_test ( );
  c8_atanh_test ( );
  c8_conj_test ( );
  c8_cos_test ( );
  c8_cosh_test ( );
  c8_cube_root_test ( );
  c8_div_test ( );
  c8_div_r8_test ( );
  c8_exp_test ( );
  c8_i_test ( );
  c8_imag_test ( );
  c8_inv_test ( );
  c8_le_l1_test ( );
  c8_le_l2_test ( );
  c8_le_li_test ( );
  c8_log_test ( );
  c8_mag_test ( );
  c8_mul_test ( );
  c8_nint_test ( );
  c8_norm_l1_test ( );
  c8_norm_l2_test ( );
  c8_norm_li_test ( );
  c8_normal_01_test ( );
  c8_one_test ( );
  c8_print_test ( );
  c8_real_test ( );
  c8_sin_test ( );
  c8_sinh_test ( );
  c8_sqrt_test ( );
  c8_sub_test ( );
  c8_tan_test ( );
  c8_tanh_test ( );
  c8_to_cartesian_test ( );
  c8_to_polar_test ( );
  c8_uniform_01_test ( );
  c8_zero_test ( );

  c8mat_identity_new_test ( );
  c8mat_indicator_new_test ( );
  c8mat_norm_fro_test ( );
  c8mat_norm_l1_test ( );
  c8mat_norm_li_test ( );
  c8mat_uniform_01_new_test ( );

  c8vec_indicator_new_test ( );
  c8vec_nint_test ( );
  c8vec_norm_l1_test ( );
  c8vec_norm_l2_test ( );
  c8vec_norm_li_test ( );
  c8vec_print_test ( );
  c8vec_sort_a_l1_test ( );
  c8vec_sort_a_l2_test ( );
  c8vec_sort_a_li_test ( );
  c8vec_spiral_new_test ( );
  c8vec_uniform_01_new_test ( );
  c8vec_unity_new_test ( );

  cartesian_to_c8_test ( );

  polar_to_c8_test ( );

  r8_atan_test ( );
  r8_uniform_01_test ( );

  r8poly2_root_test ( );
  r8poly3_root_test ( );
  r8poly4_root_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "C8LIB_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void c8_abs_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ABS_TEST tests C8_ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  int i;
  double r2;
  double r3;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_ABS_TEST\n";
  cout << "  C8_ABS computes the absolute value of a C8.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          R2=C8_ABS(C1)    R3=CABS(C1)\n";
  cout << "     ---------------------     ---------------  ---------------\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    r2 = c8_abs ( c1 );
    r3 = abs ( c1 );
    cout << "  " << setw(24) << c1
         << "  " << setw(12) << r2
         << "  " << setw(12) << r3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_acos_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ACOS_TEST tests C8_ACOS.
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
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_ACOS_TEST\n";
  cout << "  C8_ACOS computes the inverse cosine;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_ACOS(C1)           C3 = C8_COS(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_acos ( c1 );
    c3 = c8_cos ( c2 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_acosh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ACOSH_TEST tests C8_ACOSH.
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
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_ACOSH_TEST\n";
  cout << "  C8_ACOSH computes the inverse hyperbolic cosine;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_ACOSH(C1)          C3 = C8_COSH(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_acosh ( c1 );
    c3 = c8_cosh ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_add_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ADD_TEST tests C8_ADD.
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
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_ADD_TEST\n";
  cout << "  C8_ADD computes C3 = C1 + C2.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2=C8_UNIFORM_01          C3 = C8_ADD(C1,C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_uniform_01 ( seed );
    c3 = c8_add ( c1, c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_arg_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ARG_TEST tests C8_ARG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double argument;
  complex <double> c1;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "C8_ARG_TEST\n";
  cout << "  C8_ARG computes the argument of a C8.\n";
  cout << "\n";
  cout << 
    "            C1=random            ARG=C8_ARG(C1)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c8_uniform_01 ( seed );
    argument = c8_arg ( c1 );

    cout << "  " << setw(24) << c1 << ""
         << "  " << setw(12) << argument << "\n";
  }

  return;
}
//****************************************************************************80

void c8_asin_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ASIN_TEST tests C8_ASIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_ASIN_TEST\n";
  cout << "  C8_ASIN computes the inverse sine;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_ASIN(C1)           C3 = C8_SIN(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_asin ( c1 );
    c3 = c8_sin ( c2 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_asinh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ASINH_TEST tests C8_ASINH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_ASINH_TEST\n";
  cout << "  C8_ASINH computes the inverse hyperbolic sine;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_ASINH(C1)          C3 = C8_SINH(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_asinh ( c1 );
    c3 = c8_sinh ( c2 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_atan_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ATAN_TEST tests C8_ATAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_ATAN_TEST\n";
  cout << "  C8_ATAN computes the inverse tangent;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_ATAN(C1)          C3 = C8_TAN(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_atan ( c1 );
    c3 = c8_tan ( c2 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_atanh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ATANH_TEST tests C8_ATANH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_ATANH_TEST\n";
  cout << "  C8_ATANH computes the inverse hyperbolic tangent;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_ATANH(C1)         C3 = C8_TANH(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_atanh ( c1 );
    c3 = c8_tanh ( c2 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_conj_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_CONJ_TEST tests C8_CONJ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_CONJ_TEST\n";
  cout << "  C8_CONJ computes C2 = conj ( C1 ).\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2=C8_CONJ(C1)             C3=C8_CONJ(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_conj ( c1 );
    c3 = c8_conj ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2 
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_cos_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_COS_TEST tests C8_COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_COS_TEST\n";
  cout << "  C8_COS computes the cosine of a C8;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_COS(C1)           C3 = C8_COS(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_cos ( c1 );
    c3 = c8_acos ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_cosh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_COSH_TEST tests C8_COSH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_COSH_TEST\n";
  cout << "  C8_COSH computes the hyperbolic cosine of a C8;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_COSH(C1)           C3 = C8_COSH(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_cosh ( c1 );
    c3 = c8_acosh ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2 
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_cube_root_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_CUBE_ROOT_TEST tests C8_CUBE_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "C8_CUBE_ROOT_TEST\n";
  cout << "  C8_CUBE_ROOT computes the principal cube root of a C8.\n";
  cout << "\n";
  cout << 
    "            C1=random            C2=C8_CUBE_ROOT(C1)         C3=C2*C2*C2\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_cube_root ( c1 );
    c3 = c2 * c2 * c2;

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_div_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_DIV_TEST tests C8_DIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  complex <double> c4;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_DIV_TEST\n";
  cout << "  C8_DIV computes C3 = C1 / C2.\n";
  cout << "\n";
  cout << "        C1=C8_UNIFORM_01          C2=C8_UNIFORM_01  ";
  cout << "       C3=C8_DIV(C1,C2)           C4 = C1 / C2\n";
  cout << "     ---------------------     ---------------------";
  cout << "     ---------------------    ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_uniform_01 ( seed );
    c3 = c8_div ( c1, c2 );
    c4 = c1 / c2;

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3
         << "  " << setw(24) << c4 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_div_r8_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_DIV_R8_TEST tests C8_DIV_R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c3;
  complex <double> c4;
  int i;
  double r2;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_DIV_R8_TEST\n";
  cout << "  C8_DIV_R8 computes C3 = C1 / R2.\n";
  cout << "\n";
  cout << "        C1=C8_UNIFORM_01          R2=R8_UNIFORM_01  ";
  cout << "       C3=C8_DIV_R8(C1,R2)           C4 = C1 / R2\n";
  cout << "     ---------------------     ---------------------";
  cout << "     ---------------------    ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    r2 = r8_uniform_01 ( seed );
    c3 = c8_div_r8 ( c1, r2 );
    c4 = c1 / r2;

    cout << "  " << setw(24) << c1
         << "  " << setw(12) << r2 << "            "
         << "  " << setw(24) << c3
         << "  " << setw(24) << c4 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_exp_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_EXP_TEST tests C8_EXP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_EXP_TEST\n";
  cout << "  C8_EXP computes C2 = e ^ C1.\n";
  cout << "\n";
  cout << "        C1=C8_UNIFORM_01          C2=C8_EXP(C1)           C3=C8_LOG(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_exp ( c1 );
    c3 = c8_log ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2 
         << "  " << setw(24) << c3 << "\n";
  }
  return;
}
//****************************************************************************80

void c8_i_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_I_TEST tests C8_I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;

  cout << "\n";
  cout << "C8_I_TEST\n";
  cout << "  C8_I returns the imaginary unit as a C8.\n";
 
  c1 = c8_i ( );
  c8_print ( c1, "  C1 = C8_I ( ) = " );

  c2 = c1 * c1;
  c8_print ( c2, "  C2 = C1 * C1 = " );

  return;
}
//****************************************************************************80

void c8_imag_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_IMAG_TEST tests C8_IMAG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  int i;
  double r2;
  double r3;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_IMAG_TEST\n";
  cout << "  C8_IMAG computes the imaginary part of a C8.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          R2=C8_IMAG(C1)    R3=imag(C1)\n";
  cout << "     ---------------------     ---------------  ---------------\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    r2 = c8_imag ( c1 );
    r3 = imag ( c1 );
    cout << "  " << setw(24) << c1
         << "  " << setw(12) << r2
         << "  " << setw(12) << r3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_inv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_INV_TEST tests C8_INV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_INV_TEST\n";
  cout << "  C8_INV computes C2 = 1 / C1.\n";
  cout << "  Check by C3 = 1 / C2.\n";
  cout << "\n";
  cout << "        C1=C8_UNIFORM_01          C2=C8_INV(C1)             C3=C8_INV(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_inv ( c1 );
    c3 = c8_inv ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_le_l1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_LE_L1_TEST tests C8_LE_L1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  int l3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_LE_L1_TEST\n";
  cout << "  C8_LE_L1 evalues (C1 <= C2) using the L1 norm.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01()      C2 = C8_UNIFORM_01()      L3 = C8_LE_L1(C1,C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_uniform_01 ( seed );
    l3 = c8_le_l1 ( c1, c2 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "          " << l3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_le_l2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_LE_L2_TEST tests C8_LE_L2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  int l3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_LE_L2_TEST\n";
  cout << "  C8_LE_L2 evalues (C1 <= C2) using the L2 norm.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01()      C2 = C8_UNIFORM_01()      L3 = C8_LE_L2(C1,C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_uniform_01 ( seed );
    l3 = c8_le_l2 ( c1, c2 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "          " << l3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_le_li_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_LE_LI_TEST tests C8_LE_LI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  int l3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_LE_LI_TEST\n";
  cout << "  C8_LE_LI evalues (C1 <= C2) using the Loo norm.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01()      C2 = C8_UNIFORM_01()      L3 = C8_LE_LI(C1,C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_uniform_01 ( seed );
    l3 = c8_le_li ( c1, c2 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "          " << l3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_log_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_LOG_TEST tests C8_LOG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_LOG_TEST\n";
  cout << "  C8_LOG computes log ( Z ).\n";
  cout << "\n";
  cout << "        C1=C8_UNIFORM_01          C2=C8_LOG(C1)             C3=C8_EXP(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_log ( c1 );
    c3 = c8_exp ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_mag_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_MAG_TEST tests C8_MAG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  double magnitude;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "C8_MAG_TEST\n";
  cout << "  C8_MAG computes the magnitude of a C8.\n";
  cout << "\n";
  cout << 
    "            C1=random            MAG=C8_MAG(C1)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c8_uniform_01 ( seed );
    magnitude = c8_mag ( c1 );

    cout << "  " << setw(24) << c1
         << "  " << setw(12) << magnitude << "\n";
  }

  return;
}
//****************************************************************************80

void c8_mul_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_MUL_TEST tests C8_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  complex <double> c4;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_MUL_TEST\n";
  cout << "  C8_MUL computes C3 = C1 * C2.\n";
  cout << "\n";
  cout << "        C1=C8_UNIFORM_01          C2=C8_UNIFORM_01  ";
  cout << "        C3=C8_MUL(C1,C2)          C4=C1*C2\n";
  cout << "     ---------------------     ---------------------";
  cout << "     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_uniform_01 ( seed );
    c3 = c8_mul ( c1, c2 );
    c4 = c1 * c2;

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3
         << "  " << setw(24) << c4 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_nint_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NINT_TEST tests C8_NINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_NINT_TEST\n";
  cout << "  C8_NINT computes the nearest integer to a C8.\n";
  cout << "\n";
  cout << "       C1=10*C8_UNIFORM_01      C2=C8_NINT(C1)\n";
  cout << "     ---------------------     ---------------\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    c1 = 10.0 * c8_uniform_01 ( seed );
    c2 = c8_nint ( c1 );
    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_norm_l1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORM_L1_TEST tests C8_NORM_L1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  int i;
  double r2;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_NORM_L1_TEST\n";
  cout << "  C8_NORM_L1 computes the L1 norm of a C8.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01        R2=C8_NORM_L1(C1)\n";
  cout << "     ---------------------     -----------------\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    r2 = c8_norm_l1 ( c1 );
    cout << "  " << setw(24) << c1
         << "  " << setw(12) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_norm_l2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORM_L2_TEST tests C8_NORM_L2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  int i;
  double r2;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_NORM_L2_TEST\n";
  cout << "  C8_NORM_L2 computes the L2 norm of a C8.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01        R2=C8_NORM_L2(C1)\n";
  cout << "     ---------------------     -----------------\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    r2 = c8_norm_l2 ( c1 );
    cout << "  " << setw(24) << c1
         << "  " << setw(12) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_norm_li_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORM_LI_TEST tests C8_NORM_LI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  int i;
  double r2;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_NORM_LI_TEST\n";
  cout << "  C8_NORM_LI computes the Loo norm of a C8.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01        R2=C8_NORM_LI(C1)\n";
  cout << "     ---------------------     -----------------\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    r2 = c8_norm_li ( c1 );
    cout << "  " << setw(24) << c1
         << "  " << setw(12) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_normal_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORMAL_01_TEST tests C8_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  int seed = 123456789;
  int test;
  int test_num = 20;

  cout << "\n";
  cout << "C8_NORMAL_01_TEST\n";
  cout << "  C8_NORMAL_01 generates unit pseudonormal\n";
  cout << "  complex values.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    c1 = c8_normal_01 ( seed );
    cout << "  " << setw(24) << c1 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_one_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ONE_TEST tests C8_ONE.
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
  complex <double> c1;
  complex <double> c2;

  cout << "\n";
  cout << "C8_ONE_TEST\n";
  cout << "  C8_ONE returns one as a C8.\n";
 
  c1 = c8_one ( );
  c8_print ( c1, "  C1 = C8_ONE ( ) = " );

  c2 = c1 + c1;
  c8_print ( c2, "  C2 = C1 + C1 = " );

  return;
}
//****************************************************************************80

void c8_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_PRINT_TEST tests C8_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  complex <double> c4;
  complex <double> c5;
  complex <double> c6;
  complex <double> c7;
  complex <double> c8;

  cout << "\n";
  cout << "C8_PRINT_TEST\n";
  cout << "  C8_PRINT prints a C8.\n";
  cout << "\n";

  c1 = complex <double> (  0.0,               0.0 );
  c2 = complex <double> (  1.0,               0.0 );
  c3 = complex <double> (  3.141592653589793, 0.0 );
  c4 = complex <double> (  0.0,               1.0 );
  c5 = complex <double> (  1.0,               2.0 );
  c6 = complex <double> ( -12.34,            56.78 );
  c7 = complex <double> (  0.001,             0.000002 );
  c8 = complex <double> (  3.0E+08,          -4.5E+09 );

  c8_print ( c1, "  Zero:" );
  c8_print ( c2, "  One:" );
  c8_print ( c3, "  Pi:" );
  c8_print ( c4, "  i:" );
  c8_print ( c5, "  1+2i:" );
  c8_print ( c6, " -12.34 + 56.78i:" );
  c8_print ( c7, "  1E-3 + 2E-6i:" );
  c8_print ( c8, "  3E8 - 4.5E9i:" );

  return;
}
//****************************************************************************80

void c8_real_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_REAL_TEST tests C8_REAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  int i;
  double r2;
  double r3;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_REAL_TEST\n";
  cout << "  C8_REAL computes the real part of a C8.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          R2=C8_REAL(C1)    R3=REAL(C1)_\n";
  cout << "     ---------------------     ---------------  ---------------\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    r2 = c8_real ( c1 );
    r3 = real ( c1 );
    cout << "  " << setw(24) << c1
         << "  " << setw(12) << r2 
         << "  " << setw(12) << r3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_sin_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_SIN_TEST tests C8_SIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_SIN_TEST\n";
  cout << "  C8_SIN computes the sine of a C8;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_SIN(C1)            C3 = C8_ASIN(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_sin ( c1 );
    c3 = c8_asin ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2 
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_sinh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_SINH_TEST tests C8_SINH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_SINH_TEST\n";
  cout << "  C8_SINH computes the hyperbolic sine of a C8;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_SINH(C1)           C3 = C8_ASINH(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_sinh ( c1 );
    c3 = c8_asinh ( c1 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_sqrt_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_SQRT_TEST tests C8_SQRT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "C8_SQRT_TEST\n";
  cout << "  C8_SQRT computes the principal square root of a C8.\n";
  cout << "\n";
  cout << 
    "            C1=random            C2=C8_SQRT(C1)              C3=C2*C2\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_sqrt ( c1 );
    c3 = c2 * c2;

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_sub_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_SUB_TEST tests C8_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_SUB_TEST\n";
  cout << "  C8_SUB computes C3 = C1 - C2.\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2=C8_UNIFORM_01          C3 = C8_SUB(C1,C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_uniform_01 ( seed );
    c3 = c8_sub ( c1, c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_tan_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_TAN_TEST tests C8_TAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_TAN_TEST\n";
  cout << "  C8_TAN computes the tangent of a C8;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01          C2 = C8_TAN(C1)           C3 = C8_ATAN(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_tan ( c1 );
    c3 = c8_atan ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_tanh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_TANH_TEST tests C8_TANH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_TANH_TEST\n";
  cout << "  C8_TANH computes the hyperbolic tangent of a C8;\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01         C2 = C8_TANH(C1)            C3 = C8_ATANH(C2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c2 = c8_tanh ( c1 );
    c3 = c8_atanh ( c2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(24) << c2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_to_cartesian_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_TO_CARTESIAN_TEST tests C8_TO_CARTESIAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c3;
  int i;
  int seed;
  double x2;
  double y2;

  seed = 123456789;

  cout << "\n";
  cout << "C8_TO_CARTESIAN_TEST\n";
  cout << "  C8_TO_CARTESIAN computes C8 -> ( X, Y ).\n";
  cout << "\n";
  cout << "        C1=C8_UNIFORM_01    (X2,Y2)=C8_TO_CARTESIAN(C1)     C3=CARTESIAN_TO_C8(X2,Y2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c8_to_cartesian ( c1, x2, y2 );
    c3 = cartesian_to_c8 ( x2, y2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(12) << x2
         << "  " << setw(12) << y2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_to_polar_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_TO_POLAR_TEST tests C8_TO_POLAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c3;
  int i;
  double r2;
  int seed;
  double t2;

  seed = 123456789;

  cout << "\n";
  cout << "C8_TO_POLAR_TEST\n";
  cout << "  C8_TO_POLAR computes C8 -> ( R, T ).\n";
  cout << "\n";
  cout << "        C1=C8_UNIFORM_01     (X2,Y2)=C8_TO_POLAR(C1)     C3=POLAR_TO_C8(X2,Y2)\n";
  cout << "     ---------------------     ---------------------     ---------------------\n";

  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( seed );
    c8_to_polar ( c1, r2, t2 );
    c3 = polar_to_c8 ( r2, t2 );

    cout << "  " << setw(24) << c1
         << "  " << setw(12) << r2
         << "  " << setw(12) << t2
         << "  " << setw(24) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void c8_uniform_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_UNIFORM_01_TEST tests C8_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c;
  int i;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "C8_UNIFORM_01_TEST\n";
  cout << "  C8_UNIFORM_01 returns a uniformly random \"unit\" C8\n";
  cout << "\n";
  cout << "       C1=C8_UNIFORM_01(SEED)\n";
  cout << "     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    c = c8_uniform_01 ( seed );

    cout << "  " << setw(24) << c << "\n";
  }

  return;
}
//****************************************************************************80

void c8_zero_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ZERO_TEST tests C8_ZERO.
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
  complex <double> c1;

  cout << "\n";
  cout << "C8_ZERO_TEST\n";
  cout << "  C8_ZERO returns zero as a C8.\n";
 
  c1 = c8_zero ( );
  c8_print ( c1, "  C1 = C8_ZERO ( ) = " );

  return;
}
//****************************************************************************80

void c8mat_identity_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_IDENTITY_NEW_TEST tests C8MAT_IDENTITY_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int m = 4;
  int n = 4;

  cout << "\n";
  cout << "C8MAT_IDENTITY_NEW_TEST\n";
  cout << "  C8MAT_IDENTITY_NEW returns the complex identity matrix.\n";

  a = c8mat_identity_new ( n );

  c8mat_print ( m, n, a, "  The identity matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void c8mat_indicator_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_INDICATOR_NEW_TEST tests C8MAT_INDICATOR_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int m = 5;
  int n = 3;

  cout << "\n";
  cout << "C8MAT_INDICATOR_NEW_TEST\n";
  cout << "  C8MAT_INDICATOR_NEW returns the complex indicator matrix.\n";

  a = c8mat_indicator_new ( m, n );

  c8mat_print ( m, n, a, "  The indicator matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void c8mat_norm_fro_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_NORM_FRO_TEST tests C8MAT_NORM_FRO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int m = 5;
  int n = 4;
  double value;

  cout << "\n";
  cout << "C8MAT_NORM_FRO_TEST\n";
  cout << "  C8MAT_NORM_FRO returns the Frobenius norm of a matrix.\n";

  a = c8mat_indicator_new ( m, n );

  c8mat_print ( m, n, a, "  The indicator matrix:" );
 
  value = c8mat_norm_fro ( m, n, a );

  cout << "\n";
  cout << "  Frobenius norm = " << value << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void c8mat_norm_l1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_NORM_L1_TEST tests C8MAT_NORM_L1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int m = 5;
  int n = 4;
  double value;

  cout << "\n";
  cout << "C8MAT_NORM_L1_TEST\n";
  cout << "  C8MAT_NORM_L1 returns the L1 norm of a matrix.\n";

  a = c8mat_indicator_new ( m, n );

  c8mat_print ( m, n, a, "  The indicator matrix:" );
 
  value = c8mat_norm_l1 ( m, n, a );

  cout << "\n";
  cout << "  L1 norm = " << value << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void c8mat_norm_li_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_NORM_LI_TEST tests C8MAT_NORM_LI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int m = 5;
  int n = 4;
  double value;

  cout << "\n";
  cout << "C8MAT_NORM_LI_TEST\n";
  cout << "  C8MAT_NORM_LI returns the Loo norm of a matrix.\n";

  a = c8mat_indicator_new ( m, n );

  c8mat_print ( m, n, a, "  The indicator matrix:" );
 
  value = c8mat_norm_li ( m, n, a );

  cout << "\n";
  cout << "  Loo norm = " << value << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void c8mat_uniform_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_UNIFORM_01_NEW_TEST tests C8MAT_UNIFORM_01_NEW.
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
  complex <double> *a;
  int m = 5;
  int n = 4;
  int seed = 123456789;

  cout << "\n";
  cout << "C8MAT_UNIFORM_01_NEW_TEST\n";
  cout << "  C8MAT_UNIFORM_01_NEW computes a random complex matrix.\n";

  a = c8mat_uniform_01_new ( m, n, seed );

  c8mat_print ( m, n, a, "  The matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_indicator_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_INDICATOR_NEW_TEST tests C8VEC_INDICATOR_NEW;
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
  complex <double> *a;
  int n = 10;

  cout << "\n";
  cout << "C8VEC_INDICATOR_NEW_TEST\n";
  cout << "  C8VEC_INDICATOR_NEW sets A = (1-1i,2-2i,...,N-Ni)\n";

  a = c8vec_indicator_new ( n );
 
  c8vec_print ( n, a, "  The indicator vector:" );

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_nint_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_NINT_TEST tests C8VEC_NINT_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int i;
  int n = 5;
  complex <double> s;
  int seed;

  cout << "\n";
  cout << "C8VEC_NINT_TEST\n";
  cout << "  C8VEC_NINT rounds a complex vector.\n";

  seed = 123456789;

  a = c8vec_uniform_01_new ( n, seed );

  s = complex <double> ( 5.0, 3.0 );

  for ( i = 0; i < n; i++ )
  {
    a[i] = s * a[i];
  }

  c8vec_print ( n, a, "  The initial vector:" );

  c8vec_nint ( n, a );

  c8vec_print ( n, a, "  The random C8VEC:" );

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_norm_l1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_NORM_L1_TEST tests C8VEC_NORM_L1;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int n = 5;
  double value;

  cout << "\n";
  cout << "C8VEC_NORM_L1_TEST\n";
  cout << "  C8VEC_NORM_L1 computes the L1 norm of a C8VEC\n";

  a = c8vec_indicator_new ( n );
 
  c8vec_print ( n, a, "  The indicator vector:" );

  value = c8vec_norm_l1 ( n, a );
  cout << "\n";
  cout << "  L1 norm = " << value << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_norm_l2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_NORM_L2_TEST tests C8VEC_NORM_L2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int n = 5;
  double value;

  cout << "\n";
  cout << "C8VEC_NORM_L2_TEST\n";
  cout << "  C8VEC_NORM_L2 computes the L2 norm of a C8VEC\n";

  a = c8vec_indicator_new ( n );
 
  c8vec_print ( n, a, "  The indicator vector:" );

  value = c8vec_norm_l2 ( n, a );
  cout << "\n";
  cout << "  L2 norm = " << value << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_norm_li_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_NORM_LI_TEST tests C8VEC_NORM_LI;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int n = 5;
  double value;

  cout << "\n";
  cout << "C8VEC_NORM_LI_TEST\n";
  cout << "  C8VEC_NORM_LI computes the Loo norm of a C8VEC\n";

  a = c8vec_indicator_new ( n );
 
  c8vec_print ( n, a, "  The indicator vector:" );

  value = c8vec_norm_li ( n, a );
  cout << "\n";
  cout << "  Loo norm = " << value << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_PRINT_TEST tests C8VEC_PRINT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int n = 5;

  cout << "\n";
  cout << "C8VEC_PRINT_TEST\n";
  cout << "  C8VEC_PRINT prints a C8VEC\n";

  a = c8vec_indicator_new ( n );
 
  c8vec_print ( n, a, "  The indicator vector:" );

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_sort_a_l1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_SORT_A_L1_TEST tests C8VEC_SORT_A_L1;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int i;
  int n = 10;
  int seed;
  double value;

  cout << "\n";
  cout << "C8VEC_SORT_A_L1_TEST\n";
  cout << "  C8VEC_SORT_A_L1 sorts a C8VEC by L1 norm\n";

  seed = 123456789;
  a = c8vec_uniform_01_new ( n, seed );
 
  c8vec_print ( n, a, "  The unsorted vector:" );

  c8vec_sort_a_l1 ( n, a );

  cout << "\n";
  cout << "   I             A(I)                   L1(A(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(28) << a[i]
         << "  " << setw(14) <<  c8_norm_l1 ( a[i] ) << "\n";
  }

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_sort_a_l2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_SORT_A_L2_TEST tests C8VEC_SORT_A_L2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int i;
  int n = 10;
  int seed;
  double value;

  cout << "\n";
  cout << "C8VEC_SORT_A_L2_TEST\n";
  cout << "  C8VEC_SORT_A_L2 sorts a C8VEC by L2 norm\n";

  seed = 123456789;
  a = c8vec_uniform_01_new ( n, seed );
 
  c8vec_print ( n, a, "  The unsorted vector:" );

  c8vec_sort_a_l2 ( n, a );

  cout << "\n";
  cout << "   I             A(I)                   L2(A(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(28) << a[i]
         << "  " << setw(14) <<  c8_norm_l2 ( a[i] ) << "\n";
  }

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_sort_a_li_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_SORT_A_LI_TEST tests C8VEC_SORT_A_LI;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int i;
  int n = 10;
  int seed;
  double value;

  cout << "\n";
  cout << "C8VEC_SORT_A_LI_TEST\n";
  cout << "  C8VEC_SORT_A_LI sorts a C8VEC by Loo norm\n";

  seed = 123456789;
  a = c8vec_uniform_01_new ( n, seed );
 
  c8vec_print ( n, a, "  The unsorted vector:" );

  c8vec_sort_a_li ( n, a );

  cout << "\n";
  cout << "   I             A(I)                   Loo(A(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(28) << a[i]
         << "  " << setw(14) <<  c8_norm_li ( a[i] ) << "\n";
  }

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_spiral_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_SPIRAL_NEW_TEST tests C8VEC_SPIRAL_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *c;
  complex <double> c1;
  complex <double> c2;
  int m;
  int n = 13;

  cout << "\n";
  cout << "C8VEC_SPIRAL_NEW_TEST\n";
  cout << "  C8VEC_SPIRAL_NEW returns N points on a spiral\n";
  cout << "  which includes M complete turns.\n";

  m = 1;
  c1 = complex <double> ( 5.0, 0.0 );
  c2 = complex <double> ( 3.0, 0.0 );

  c = c8vec_spiral_new ( n, m, c1, c2 );

  c8vec_print ( n, c, "  The spiral points:" );

  delete [] c;

  return;
}
//****************************************************************************80

void c8vec_uniform_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_UNIFORM_01_NEW_TEST tests C8VEC_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int n = 5;
  int seed = 123456789;

  cout << "\n";
  cout << "C8VEC_UNIFORM_01_NEW_TEST\n";
  cout << "  C8VEC_UNIFORM_01_NEW computes a random complex vector.\n";

  a = c8vec_uniform_01_new ( n, seed );

  c8vec_print ( n, a, "  The random C8VEC:" );

  delete [] a;

  return;
}
//****************************************************************************80

void c8vec_unity_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_UNITY_NEW_TEST tests C8VEC_UNITY_NEW;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int n = 12;

  cout << "\n";
  cout << "C8VEC_UNITY_NEW_TEST\n";
  cout << "  C8VEC_UNITY_NEW sets A to the N roots of unity.\n";

  a = c8vec_unity_new ( n );
 
  c8vec_print ( n, a, "  The N roots of unity:" );

  delete [] a;

  return;
}
//****************************************************************************80

void cartesian_to_c8_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CARTESIAN_TO_C8_TEST tests CARTESIAN_TO_C8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c2;
  int i;
  int seed;
  double x1;
  double y1;
  double x3;
  double y3;

  seed = 123456789;

  cout << "\n";
  cout << "CARTESIAN_TO_C8_TEST\n";
  cout << "  CARTESIAN_TO_C8 computes ( X, Y ) -> C8.\n";
  cout << "\n";
  cout << "        X1,Y1=R8_UNIFORM_01    C2=CARTESIAN_TO_C8(X1,Y1)   (X3,Y3)=C8_TO_CARTESIAN(C2)     \n";
  cout << "     ---------------------     ---------------------     ---------------------\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    x1 = r8_uniform_01 ( seed );
    y1 = r8_uniform_01 ( seed );
    c2 = cartesian_to_c8 ( x1, y1 );
    c8_to_cartesian ( c2, x3, y3 );

    cout << "  " << setw(12) << x1
         << "  " << setw(12) << y1
         << "  " << setw(24) << c2 
         << "  " << setw(12) << x3
         << "  " << setw(12) << y3 << "\n";
  }

  return;
}
//****************************************************************************80

void polar_to_c8_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLAR_TO_C8_TEST tests POLAR_TO_C8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c2;
  int i;
  double r1;
  double r3;
  const double r8_pi = 3.141592653589793;
  int seed;
  double t1;
  double t3;

  seed = 123456789;

  cout << "\n";
  cout << "POLAR_TO_C8_TEST\n";
  cout << "  POLAR_TO_C8 computes ( R, T ) -> C8.\n";
  cout << "\n";
  cout << "        X1,Y1=R8_UNIFORM_01    C2=POLAR_TO_C8(X1,Y1)   (X3,Y3)=C8_TO_POLAR(C2)     \n";
  cout << "     ---------------------     ---------------------     ---------------------\n";

  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    r1 = r8_uniform_01 ( seed );
    t1 = r8_uniform_01 ( seed ) * 2.0 * r8_pi;
    c2 = polar_to_c8 ( r1, t1 );
    c8_to_polar ( c2, r3, t3 );

    cout << "  " << setw(12) << r1
         << "  " << setw(12) << t1
         << "  " << setw(24) << c2
         << "  " << setw(12) << r3
         << "  " << setw(12) << t3 << "\n";
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
  min = x[0];
  max = x[0];
  mean = 0.0;
  for ( i = 0; i < N; i++ )
  {
    if ( x[i] < min )
    {
      min = x[i];
    }
    if ( max < x[i] )
    {
      max = x[i];
    }
    mean = mean + x[i];
  }
  mean = mean / ( double ) ( N );
  variance = 0.0;
  for ( i = 0; i < N; i++ )
  {
    variance = variance + pow ( x[i] - mean, 2 );
  }
  variance =  variance / ( double ) ( N );


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

void r8poly2_root_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY2_ROOT_TEST tests R8POLY2_ROOT.
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
# define TEST_NUM 3

  double a;
  double a_test[TEST_NUM] = { 2.0, 1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = { -2.0, -20.0, -2.0 };
  double c;
  double c_test[TEST_NUM] = { -24.0, 100.0, 10.0 };
  complex <double> r1;
  complex <double> r2;
  int test;

  cout << "\n";
  cout << "R8POLY2_ROOT_TEST\n";
  cout << "  R8POLY2_ROOT finds quadratic equation roots.\n";
  cout << "\n";
  cout << "         A         B         C     R1         R2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];

    r8poly2_root ( a, b, c, r1, r2 );

    cout << "  " << setw(8) << a
         << "  " << setw(8) << b
         << "  " << setw(8) << c
         << "  " << setw(24) << r1
         << "  " << setw(24) << r2 << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8poly3_root_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY3_ROOT_TEST tests R8POLY3_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double a;
  double a_test[TEST_NUM] = { 1.0, 9.0, 1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = { -6.0, -36.0, -5.0, -8.0 };
  double c;
  double c_test[TEST_NUM] = { 11.0, 54.0, 8.0, 25.0 };
  double d;
  double r8_test[TEST_NUM] = { -6.0, -27.0, -4.0, -26.0 };
  complex <double> r1;
  complex <double> r2;
  complex <double> r3;
  int test;
//
//  1: Three distinct real roots, 1, 2, 3.
//  2: One repeated real root, 1.5, 1.5, 1.5.
//  3: Two real roots, one repeated, 1, 2, 2.
//  4: One real root, a complex conjugate pair, 2, 3+2I, 3-2I.
//
  cout << "\n";
  cout << "R8POLY3_ROOT_TEST\n";
  cout << "  R8POLY3_ROOT finds roots of cubic equations.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];

    cout << "\n";
    cout << "  Polynomial coefficients A, B, C, D:\n";
    cout << "\n";

    cout << "  A = " << setw(8) << a << "\n";
    cout << "  B = " << setw(8) << b << "\n";
    cout << "  C = " << setw(8) << c << "\n";
    cout << "  D = " << setw(8) << d << "\n";

    r8poly3_root ( a, b, c, d, r1, r2, r3 );

    cout << "\n";
    cout << "  Roots:\n";
    cout << "\n";
    cout << "  " << r1 << "\n";
    cout << "  " << r2 << "\n";
    cout << "  " << r3 << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void r8poly4_root_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY4_ROOT_TEST tests R8POLY4_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 7

  double a;
  double a_test[TEST_NUM] = {
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = {
    -10.0, -5.0, -22.0, -16.0, -20.0,
    2.0, 0.0 };
  double c;
  double c_test[TEST_NUM] = {
    35.0, 1.0, 141.0, 72.0, 150.0,
    1.0, 13.0 };
  double d;
  double r8_test[TEST_NUM] = {
    -50.0, 21.0, -220.0, -128.0, -500.0,
    8.0, 0.0 };
  double e;
  double e_test[TEST_NUM] = {
    24.0, -18.0, +100.0, 80.0, 625.0,
    -12.0, 36.0 };
  complex <double> r1;
  complex <double> r2;
  complex <double> r3;
  complex <double> r4;
  int test;
//
//  1: Four distinct real roots, 1, 2, 3, 4.
//  2: Three distinct real roots, 1, -2, 3, 3
//  3: Two distinct real roots, 1, 1, 10, 10.
//  4: Two distinct real roots, 2, 2, 2, 10
//  5: One real root, 5, 5, 5, 5
//  6: Two distinct real roots, one complex conjugate pair.
//  7: Two distinct complex conjugate pairs.
//
  cout << "\n";
  cout << "R8POLY4_ROOT_TEST\n";
  cout << "  R8POLY4_ROOT finds roots of quartic equations.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];
    e = e_test[test];

    cout << "\n";
    cout << "  A = " << setw(8) << a << "\n";
    cout << "  B = " << setw(8) << b << "\n";
    cout << "  C = " << setw(8) << c << "\n";
    cout << "  D = " << setw(8) << d << "\n";
    cout << "  E = " << setw(8) << e << "\n";

    r8poly4_root ( a, b, c, d, e, r1, r2, r3, r4 );

    cout << "\n";
    cout << "  Roots:\n";
    cout << "\n";
    cout << "  " << r1 << "\n";
    cout << "  " << r2 << "\n";
    cout << "  " << r3 << "\n";
    cout << "  " << r4 << "\n";
  }

  return;
# undef TEST_NUM
}
