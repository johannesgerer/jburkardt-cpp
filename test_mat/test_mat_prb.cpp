# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <cstring>

using namespace std;

# include "test_mat.hpp"

int main ( );
void test_cond ( );
void test_determinant ( );
void test_eigen ( );
void test_inverse ( );
void test_null ( );
void test_plu ( );
void test_solution ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_MAT_PRB.
//
//  Discussion:
//
//    TEST_MAT_PRB calls the TEST_MAT test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TEST_MAT_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_MAT library.\n";

  test_cond ( );
  test_determinant ( );
  test_eigen ( );
  test_inverse ( );
  test_null ( );
  test_plu ( );
  test_solution ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_MAT_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test_cond ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_COND tests the condition number computations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  double cond;
  int n;
  int seed;

  cout << "\n";
  cout << "TEST_COND\n";
  cout << "  Compute the condition number of an example of each\n";
  cout << "  test matrix\n";
  cout << "\n";
  cout << "  Matrix title             N      COND\n";
  cout << "\n";
//
//  AEGERTER matrix.
//
  n = 5;
  cond = aegerter_condition ( n );
  cout << "  " << setw(20) << "AEGERTER            "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  BAB matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  cond = bab_condition ( n, alpha, beta );
  cout << "  " << setw(20) << "BAB                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  BODEWIG matrix.
//
  n = 4;
  cond = bodewig_condition ( );
  cout << "  " << setw(20) << "BODEWIG             "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  COMBIN matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  cond = combin_condition ( alpha, beta, n );
  cout << "  " << setw(20) << "COMBIN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  CONEX3 matrix.
//
  n = 5;
  cond = conex3_condition ( n );
  cout << "  " << setw(20) << "CONEX3              "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  RUTIS5 matrix.
//
  n = 4;
  cond = rutis5_condition ( );
  cout << "  " << setw(20) << "RUTIS5              "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  SUMMATION matrix.
//
  n = 5;
  cond = summation_condition ( n );
  cout << "  " << setw(20) << "SUMMATION           "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  TRI_UPPER matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  cond = tri_upper_condition ( alpha, n );
  cout << "  " << setw(20) << "TRI_UPPER           "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  WILK03 matrix.
//
  n = 3;
  cond = wilk03_condition ( );
  cout << "  " << setw(20) << "WILK03              "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";
//
//  WILSON matrix.
//
  n = 4;
  cond = wilson_condition ( );
  cout << "  " << setw(20) << "WILSON              "
       << "  " << setw(4) << n
       << "  " << setw(14) << cond << "\n";

  return;
}
//****************************************************************************80

void test_determinant ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_DETERMINANT tests the determinant computations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double b;
  double beta;
  int col_num;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double da;
  double determ1;
  double determ2;
  double di;
  double gamma;
  int i;
  int ii;
  int jj;
  int k;
  double *l;
  int m;
  int n;
  double norm_frobenius;
  double *p;
  double perturb;
  int *pivot;
  double prob;
  int rank;
  int row_num;
  int seed;
  int seed_save;
  string title;
  double *u;
  double *v1;
  double *v2;
  double *v3;
  double *w;
  double *x;
  int x_n;
  double *y;
  int y_n;
  double y_sum;
  double *z;

  cout << "\n";
  cout << "TEST_DETERMINANT\n";
  cout << "  Compute the determinants of an example of each\n";
  cout << "  test matrix; compare with the determinant routine,\n";
  cout << "  if available.  Print the matrix Frobenius norm\n";
  cout << "  for an estimate of magnitude.\n";
  cout << "\n";
  cout << "  Matrix title             N          "
       << "Determ          Determ          ||A||\n";
  cout << "\n";
//
//  AEGERTER matrix.
//
  title = "AEGERTER            ";
  n = 5;
  a = aegerter ( n );
  determ1 = aegerter_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << title
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ANTICIRCULANT matrix.
//
  title = "ANTICIRCULANT       ";
  n = 3;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( r8_nint ( 50.0 * x[i] - 25.0 ) ) / 5.0;
  }
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << title
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  ANTICIRCULANT matrix.
//
  title = "ANTICIRCULANT       ";
  n = 4;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( r8_nint ( 50.0 * x[i] - 25.0 ) ) / 5.0;
  }
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << title
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  ANTICIRCULANT matrix.
//
  title = "ANTICIRCULANT       ";
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( r8_nint ( 50.0 * x[i] - 25.0 ) ) / 5.0;
  }
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << title
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  ANTIHADAMARD matrix.
//
  title = "ANTIHADAMARD        ";
  n = 5;
  a = antihadamard ( n );
  determ1 = antihadamard_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << title
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ANTISYMM_RANDOM matrix.
//
  title = "ANTISYMM_RANDOM     ";
  n = 5;
  seed = 123456789;
  a = antisymm_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << title
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ANTISYMM_RANDOM matrix.
//
  title = "ANTISYMM_RANDOM     ";
  n = 6;
  seed = 123456789;
  a = antisymm_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << title
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BAB matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  a = bab ( n, alpha, beta );
  determ1 = bab_determinant ( n, alpha, beta );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "BAB                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BIMARKOV_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  a = bimarkov_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "BIMARKOV_RANDOM     "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BIS matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  a = bis ( alpha, beta, n, n );
  determ1 = bis_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "BIS                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BODEWIG matrix.
//
  n = 4;
  a = bodewig ( );
  determ1 = bodewig_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "BODEWIG             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BOOTHROYD matrix.
//
  n = 5;
  a = boothroyd ( n );
  determ1 = boothroyd_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "BOOTHROYD           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BORDERBAND matrix.
//
  n = 5;
  a = borderband ( n );
  determ1 = borderband_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "BORDERBAND          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CARRY matrix.
//
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 2, 20, &seed );
  a = carry ( k, n );
  determ1 = carry_determinant ( k, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CARRY               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CAUCHY matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  y = r8vec_uniform_01_new ( n, &seed );
  a = cauchy ( n, x, y );
  determ1 = cauchy_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CAUCHY              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  CHEBY_DIFF1 matrix.
//
  n = 5;
  a = cheby_diff1 ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CHEBY_DIFF1         "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHEBY_DIFF1 matrix.
//
  n = 6;
  a = cheby_diff1 ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CHEBY_DIFF1         "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHEBY_T matrix.
//
  n = 5;
  a = cheby_t ( n );
  determ1 = cheby_t_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CHEBY_T             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHEBY_U matrix.
//
  n = 5;
  a = cheby_u ( n );
  determ1 = cheby_u_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CHEBY_U             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHEBY_VAN1 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  a = cheby_van1 ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CHEBY_VAN1          "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  CHEBY_VAN2 matrix.
//
  for ( n = 2; n <= 10; n++ )
  {
    a = cheby_van2 ( n );
    determ1 = cheby_van2_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << "CHEBY_VAN2          "
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
//
//  CHEBY_VAN3 matrix.
//
  n = 5;
  a = cheby_van3 ( n );
  determ1 = cheby_van3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CHEBY_VAN3          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHOW matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = chow ( alpha, beta, n, n );
  determ1 = chow_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CHOW                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CIRCULANT matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = circulant ( n, n, x );
  determ1 = circulant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CIRCULANT           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  CIRCULANT2 matrix.
//
  n = 3;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CIRCULANT2          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CIRCULANT2 matrix.
//
  n = 4;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CIRCULANT2          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CIRCULANT2 matrix.
//
  n = 5;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CIRCULANT2          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CLEMENT1 matrix.
//
  n = 5;
  a = clement1 ( n );
  determ1 = clement1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CLEMENT1            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CLEMENT1 matrix.
//
  n = 6;
  a = clement1 ( n );
  determ1 = clement1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CLEMENT1            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CLEMENT2 matrix.
//
  n = 5;
  a = clement2 ( n );
  determ1 = clement2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CLEMENT2            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CLEMENT2 matrix.
//
  n = 6;
  a = clement2 ( n );
  determ1 = clement2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CLEMENT2            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CLEMENT3 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  y = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    y[i] = r8_nint ( 10.0 * y[i] - 5.0 );
  }
  a = clement3 ( n, x, y );
  determ1 = clement3_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CLEMENT3            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  CLEMENT3 matrix.
//
  n = 6;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  y = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    y[i] = r8_nint ( 10.0 * y[i] - 5.0 );
  }
  a = clement3 ( n, x, y );
  determ1 = clement3_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CLEMENT3            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  COMBIN matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = combin ( alpha, beta, n );
  determ1 = combin_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "COMBIN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  COMPANION matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  a = companion ( n, x );
  determ1 = companion_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "COMPANION           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  COMPLEX_I matrix.
//
  n = 2;
  a = complex_i ( );
  determ1 = complex_i_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "COMPLEX_I           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONEX1 matrix.
//
  n = 4;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = conex1 ( alpha );
  determ1 = conex1_determinant ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CONEX1              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONEX2 matrix.
//
  n = 3;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = conex2 ( alpha );
  determ1 = conex2_determinant ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CONEX2              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONEX3 matrix.
//
  n = 5;
  a = conex3 ( n );
  determ1 = conex3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CONEX5              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONFERENCE matrix.
//
  n = 6;
  a = conference ( n );
  determ1 = conference_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CONFERENCE          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CREATION matrix.
//
  n = 5;
  a = creation ( n, n );
  determ1 = creation_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CREATION            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB2 matrix.
//
  n = 4;
  a = daub2 ( n );
  determ1 = daub2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DAUB2               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB4 matrix.
//
  n = 8;
  a = daub4 ( n );
  determ1 = daub4_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DAUB4               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB6 matrix.
//
  n = 12;
  a = daub6 ( n );
  determ1 = daub6_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DAUB6               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB8 matrix.
//
  n = 16;
  a = daub8 ( n );
  determ1 = daub8_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DAUB8               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB10 matrix.
//
  n = 20;
  a = daub10 ( n );
  determ1 = daub10_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DAUB10              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB12 matrix.
//
  n = 24;
  a = daub12 ( n );
  determ1 = daub12_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DAUB12              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIAGONAL matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  a = diagonal ( n, n, x );
  determ1 = diagonal_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DIAGONAL            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  DIF1 matrix.
//
  n = 5;
  a = dif1 ( n, n );
  determ1 = dif1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DIF1                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIF1CYCLIC matrix.
//
  n = 5;
  a = dif1cyclic ( n );
  determ1 = dif1cyclic_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DIF1CYCLIC          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIF2 matrix.
//
  n = 5;
  a = dif2 ( n, n );
  determ1 = dif2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DIF2                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIF2CYCLIC matrix.
//
  n = 5;
  a = dif2cyclic ( n );
  determ1 = dif2cyclic_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DIF2CYCLIC          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DORR matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = dorr ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DORR                "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DOWNSHIFT matrix.
//
  n = 5;
  a = downshift ( n );
  determ1 = downshift_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DOWNSHIFT           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  EBERLEIN matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = eberlein ( alpha, n );
  determ1 = eberlein_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "EBERLEIN            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  EULERIAN matrix.
//
  n = 5;
  a = eulerian ( n, n );
  determ1 = eulerian_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "EULERIAN            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  EXCHANGE matrix.
//
  n = 5;
  a = exchange ( n, n );
  determ1 = exchange_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "EXCHANGE            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FIBONACCI1 matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = fibonacci1 ( n, alpha, beta );
  determ1 = fibonacci1_determinant ( n, alpha, beta );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "FIBONACCI1          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FIBONACCI2 matrix.
//
  n = 5;
  a = fibonacci2 ( n );
  determ1 = fibonacci2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "FIBONACCI2          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FIBONACCI3 matrix.
//
  n = 5;
  a = fibonacci3 ( n );
  determ1 = fibonacci3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "FIBONACCI3          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FIEDLER matrix.
//
  n = 7;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = fiedler ( n, n, x );
  determ1 = fiedler_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "FIEDLER             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  FORSYTHE matrix.
//
  n = 7;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = forsythe ( alpha, beta, n );
  determ1 = forsythe_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "FORSYTHE            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FOURIER_COSINE matrix.
//
  n = 5;
  a = fourier_cosine ( n );
  determ1 = fourier_cosine_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "FOURIER_COSINE      "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FOURIER_SINE matrix.
//
  n = 5;
  a = fourier_sine ( n );
  determ1 = fourier_sine_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "FOURIER_SINE        "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FRANK matrix.
//
  n = 5;
  a = frank ( n );
  determ1 = frank_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "FRANK               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GEAR matrix.
//
  for ( n = 4; n <= 8; n++ )
  {
    seed = 123456789;
    ii = i4_uniform ( -n, n, &seed );
    jj = i4_uniform ( -n, n, &seed );
    a = gear ( ii, jj, n );
    determ1 = gear_determinant ( ii, jj, n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << "GEAR                "
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
//
//  GFPP matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = gfpp ( n, alpha );
  determ1 = gfpp_determinant ( n, alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "GFPP                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GIVENS matrix.
//
  n = 5;
  a = givens ( n, n );
  determ1 = givens_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "GIVENS              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GK316 matrix.
//
  n = 5;
  a = gk316 ( n );
  determ1 = gk316_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "GK316               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GK323 matrix.
//
  n = 5;
  a = gk323 ( n, n );
  determ1 = gk323_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "GK323               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GK324 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = gk324 ( n, n, x );
  determ1 = gk324_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "GK324               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  GRCAR matrix.
//
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 1, n - 1, &seed );
  a = grcar ( n, n, k );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "GRCAR               "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HADAMARD matrix.
//
  n = 5;
  a = hadamard ( n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HADAMARD            "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HANKEL matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( 2 * n - 1, &seed );
  for ( i = 0; i < 2 * n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = hankel ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HANKEL              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  HANOWA matrix.
//
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = hanowa ( alpha, n );
  determ1 = hanowa_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HANOWA              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HARMAN matrix.
//
  n = 8;
  a = harman ( );
  determ1 = harman_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HARMAN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HARTLEY matrix.
//
  for ( n = 5; n <= 8; n++ )
  {
    a = hartley ( n );
    determ1 = hartley_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << "HARTLEY             "
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
//
//  HELMERT matrix.
//
  n = 5;
  a = helmert ( n );
  determ1 = helmert_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HELMERT             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HELMERT2 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = helmert2 ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HELMERT2            "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  HERMITE matrix.
//
  n = 5;
  a = hermite ( n );
  determ1 = hermite_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HERMITE             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HERNDON matrix.
//
  n = 5;
  a = herndon ( n );
  determ1 = herndon_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HERNDON             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HILBERT matrix.
//
  n = 5;
  a = hilbert ( n, n );
  determ1 = hilbert_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HILBERT             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HOUSEHOLDER matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = householder ( n, x );
  determ1 = householder_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "HOUSEHOLDER         "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  IDEM_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  rank = i4_uniform ( 0, n, &seed );
  a = idem_random ( n, rank, &seed );
  determ1 = idem_random_determinant ( n, rank );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "IDEM_RANDOM         "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  IDENTITY matrix.
//
  n = 5;
  a = identity ( n, n );
  determ1 = identity_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "IDENTITY            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  IJFACT1 matrix.
//
  n = 5;
  a = ijfact1 ( n );
  determ1 = ijfact1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "IJFACT1             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  IJFACT2 matrix.
//
  n = 5;
  a = ijfact2 ( n );
  determ1 = ijfact2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "IJFACT2             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ILL3 matrix.
//
  n = 3;
  a = ill3 ( );
  determ1 = ill3_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ILL3                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  INTEGRATION matrix.
//
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = integration ( alpha, n );
  determ1 = integration_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "INTEGRATION         "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  INVOL matrix.
//
  n = 5;
  a = invol ( n );
  determ1 = invol_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "INVOL               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  INVOL_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  rank = i4_uniform ( 0, n, &seed );
  a = invol_random ( n, rank, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "INVOL_RANDOM        "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  JACOBI matrix.
//
  n = 5;
  a = jacobi ( n, n );
  determ1 = jacobi_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "JACOBI              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  JACOBI matrix.
//
  n = 6;
  a = jacobi ( n, n );
  determ1 = jacobi_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "JACOBI              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  JORDAN matrix.
//
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = jordan ( alpha, n, n );
  determ1 = jordan_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "JORDAN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  KAHAN matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = kahan ( alpha, n, n );
  determ1 = kahan_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "KAHAN               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  KERSHAW matrix.
//
  n = 4;
  a = kershaw ( );
  determ1 = kershaw_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "KERSHAW             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  KERSHAWTRI matrix.
//
  n = 5;
  x_n = ( n + 1 ) / 2;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = kershawtri ( n, x );
  determ1 = kershawtri_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "KERSHAWTRI          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  KMS matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = kms ( alpha, n, n );
  determ1 = kms_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "KMS                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LAGUERRE matrix.
//
  n = 5;
  a = laguerre ( n );
  determ1 = laguerre_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LAGUERRE            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LEHMER matrix.
//
  n = 5;
  a = lehmer ( n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LEHMER              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LESLIE matrix.
//
  n = 4;
  b = 0.025;
  di = 0.010;
  da = 0.100;
  a = leslie ( b, di, da );
  determ1 = leslie_determinant ( b, di, da );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LESLIE              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LESP matrix.
//
  n = 5;
  a = lesp ( n, n );
  determ1 = lesp_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LESP                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LIETZKE matrix.
//
  n = 5;
  a = lietzke ( n );
  determ1 = lietzke_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LIETZKE             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LIGHTS_OUT matrix.
//
  if ( false )
  {
    row_num = 5;
    col_num = 5;
    n = row_num * col_num;
//  a = lights_out ( row_num, col_num, n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << "LIGHTS_OUT          "
         << "  " << setw(4) << n
         << "  " << "              "
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
  else
  {
    cout << "  " << setw(20) << "LIGHTS_OUT          -----Not ready----\n";
  }
//
//  LINE_ADJ matrix.
//
  n = 5;
  a = line_adj ( n );
  determ1 = line_adj_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LINE_ADJ            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LINE_LOOP_ADJ matrix.
//
  n = 5;
  a = line_loop_adj ( n );
  determ1 = line_loop_adj_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LINE_LOOP_ADJ       "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LOEWNER matrix.
//
  n = 5;
  seed = 123456789;
  w = r8vec_uniform_01_new ( n, &seed );
  x = r8vec_uniform_01_new ( n, &seed );
  y = r8vec_uniform_01_new ( n, &seed );
  z = r8vec_uniform_01_new ( n, &seed );
  a = loewner ( w, x, y, z, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LOEWNER             "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] w;
  delete [] x;
  delete [] y;
  delete [] z;
//
//  LOTKIN matrix.
//
  n = 5;
  a = lotkin ( n, n );
  determ1 = lotkin_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "LOTKIN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MARKOV_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  a = markov_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "MARKOV_RANDOM       "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MAXIJ matrix.
//
  n = 5;
  a = maxij ( n, n );
  determ1 = maxij_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "MAXIJ               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MILNES matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = milnes ( n, n, x );
  determ1 = milnes_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "MILNES              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  MINIJ matrix.
//
  n = 5;
  a = minij ( n, n );
  determ1 = minij_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "MINIJ               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MOLER1 matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = moler1 ( alpha, n, n );
  determ1 = moler1_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "MOLER1              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MOLER2 matrix.
//
  n = 5;
  a = moler2 ( );
  determ1 = moler2_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "MOLER2              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MOLER3 matrix.
//
  n = 5;
  a = moler3 ( n, n );
  determ1 = moler3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "MOLER3              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  NEUMANN matrix.
//
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  a = neumann ( row_num, col_num );
  determ1 = neumann_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "NEUMANN             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ONE matrix.
//
  n = 5;
  a = one ( n, n );
  determ1 = one_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ONE                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ORTEGA matrix.
//
  n = 5;
  seed = 123456789;
  v1 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v1[i] = r8_nint ( 50.0 * v1[i] - 25.0 ) / 5.0;
  }
  v2 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v2[i] = r8_nint ( 50.0 * v2[i] - 25.0 ) / 5.0;
  }
  v3 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v3[i] = r8_nint ( 50.0 * v3[i] - 25.0 ) / 5.0;
  }
  a = ortega ( n, v1, v2, v3 );
  determ1 = ortega_determinant ( n, v1, v2, v3 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ORTEGA              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] v1;
  delete [] v2;
  delete [] v3;
//
//  ORTH_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  a = orth_random ( n, &seed );
  determ1 = orth_random_determinant ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ORTH_RANDOM         "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ORTH_SYMM matrix.
//
  n = 5;
  a = orth_symm ( n );
  determ1 = orth_symm_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ORTH_SYMM           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  OTO matrix.
//
  n = 5;
  a = oto ( n, n );
  determ1 = oto_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "OTO                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PARTER matrix.
//
  n = 5;
  a = parter ( n, n );
  determ1 = parter_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PARTER              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PASCAL1 matrix.
//
  n = 5;
  a = pascal1 ( n );
  determ1 = pascal1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PASCAL1             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PASCAL2 matrix.
//
  n = 5;
  a = pascal2 ( n );
  determ1 = pascal2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PASCAL2             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PASCAL3 matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = pascal3 ( n, alpha );
  determ1 = pascal3_determinant ( n, alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PASCAL3             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PDS_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  seed_save = seed;
  a = pds_random ( n, &seed );
  seed = seed_save;
  determ1 = pds_random_determinant ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PDS_RANDOM          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PEI matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = pei ( alpha, n );
  determ1 = pei_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PEI                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PERMUTATION_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  seed_save = seed;
  a = permutation_random ( n, &seed );
  seed = seed_save;
  determ1 = permutation_random_determinant ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PERMUTATION_RANDOM  "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PLU matrix.
//
  n = 5;
  l = new double[n*n];
  p = new double[n*n];
  pivot = new int[n];
  u = new double[n*n];
  for ( i = 0; i < n; i++ )
  {
    pivot[i] = i + 1;
  }
  a = plu ( n, pivot, p, l, u );
  determ1 = plu_determinant ( n, p, l, u );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PLU                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] pivot;
  delete [] u;
//
//  POISSON matrix.
//
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  a = poisson ( row_num, col_num, n );
  determ1 = poisson_determinant ( row_num, col_num, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "POISSON             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PROLATE matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = prolate ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PROLATE             "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RECTANGLE_ADJ matrix.
//
  if ( false )
  {
    row_num = 5;
    col_num = 5;
    n = row_num * col_num;
//  a = rectangle_adj ( row_num, col_num, n );
    determ1 = rectangle_adj_determinant ( row_num, col_num );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << "RECTANGLE_ADJ       "
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
  else
  {
    cout << "  " << setw(20) << "RECTANGLE_ADJ       -----Not ready-----\n";
  }
//
//  REDHEFFER matrix.
//
  n = 5;
  a = redheffer ( n );
  determ1 = redheffer_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "REDHEFFER           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  REF_RANDOM matrix.
//
  n = 5;
  prob = 0.65;
  seed_save = 123456789;
  seed = seed_save;
  a = ref_random ( n, n, prob, &seed );
  seed = seed_save;
  determ1 = ref_random_determinant ( n, prob, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "REF_RANDOM          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  REF_RANDOM matrix.
//
  n = 5;
  prob = 0.85;
  seed_save = 123456789;
  seed = seed_save;
  a = ref_random ( n, n, prob, &seed );
  seed = seed_save;
  determ1 = ref_random_determinant ( n, prob, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "REF_RANDOM          "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RIEMANN matrix.
//
  n = 5;
  a = riemann ( n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RIEMANN             "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RING_ADJ matrix.
//
  for ( n = 1; n <= 8; n++ )
  {
    a = ring_adj ( n );
    determ1 = ring_adj_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << "RING_ADJ            "
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
//
//  RIS matrix.
//
  n = 5;
  a = ris ( n );
  determ1 = ris_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RIS                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RODMAN matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = rodman ( alpha, n, n );
  determ1 = rodman_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RODMAN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ROSSER1 matrix.
//
  n = 8;
  a = rosser1 ( );
  determ1 = rosser1_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ROSSER1             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ROUTH matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] ) / 5.0;
  }
  a = routh ( n, x );
  determ1 = routh_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ROUTH               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  RUTIS1 matrix.
//
  n = 4;
  a = rutis1 ( );
  determ1 = rutis1_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RUTIS1              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RUTIS2 matrix.
//
  n = 4;
  a = rutis2 ( );
  determ1 = rutis2_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RUTIS2              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RUTIS3 matrix.
//
  n = 4;
  a = rutis3 ( );
  determ1 = rutis3_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RUTIS3              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RUTIS4 matrix.
//
  n = 4;
  a = rutis4 ( n );
  determ1 = rutis4_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RUTIS4              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RUTIS5 matrix.
//
  n = 4;
  a = rutis5 ( );
  determ1 = rutis5_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RUTIS5              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SCHUR_BLOCK matrix.
//
  n = 5;
  x_n = ( n + 1 ) / 2;
  y_n = n / 2;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( y_n, &seed );
  for ( i = 0; i < y_n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  a = schur_block ( n, x, y );
  determ1 = schur_block_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SCHUR_BLOCK         "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  SKEW_CIRCULANT matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = skew_circulant ( n, n, x );
  determ1 = skew_circulant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SKEW_CIRCULANT      "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  SPLINE matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = spline ( n, x );
  determ1 = spline_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SPLINE              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  STIRLING matrix.
//
  n = 5;
  a = stirling ( n, n );
  determ1 = stirling_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "STIRLING            "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  STRIPE matrix.
//
  n = 5;
  a = stripe ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "STRIPE              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SUMMATION matrix.
//
  n = 5;
  a = summation ( n, n );
  determ1 = summation_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SUMMATION           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SWEET1 matrix.
//
  n = 6;
  perturb = 0.0;
  a = sweet1 ( perturb );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SWEET1              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SWEET2 matrix.
//
  n = 6;
  perturb = 0.0;
  a = sweet2 ( perturb );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SWEET2              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SWEET3 matrix.
//
  n = 6;
  perturb = 0.0;
  a = sweet3 ( perturb );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SWEET3              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SWEET4 matrix.
//
  n = 13;
  perturb = 0.0;
  a = sweet4 ( perturb );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SWEET4              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SYLVESTER matrix.
//
  n = 5;
  x_n = 3 + 1;
  y_n = 2 + 1;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( y_n, &seed );
  for ( i = 0; i < y_n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  a = sylvester ( n, x_n - 1, x, y_n - 1, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SYLVESTER           "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  SYMM_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = symm_random ( n, x, &seed );
  determ1 = symm_random_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "SYMM_RANDOM         "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  TOEPLITZ matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( 2 * n - 1, &seed );
  for ( i = 0; i < 2 * n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = toeplitz ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TOEPLITZ            "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  TOEPLITZ_5DIAG matrix.
//
  n = 5;
  seed = 123456789;
  d1 = r8_uniform_01 ( &seed );
  d1 = r8_nint ( 50.0 * d1 - 25.0 ) / 5.0;
  d2 = r8_uniform_01 ( &seed );
  d2 = r8_nint ( 50.0 * d2 - 25.0 ) / 5.0;
  d3 = r8_uniform_01 ( &seed );
  d3 = r8_nint ( 50.0 * d3 - 25.0 ) / 5.0;
  d4 = r8_uniform_01 ( &seed );
  d4 = r8_nint ( 50.0 * d4 - 25.0 ) / 5.0;
  d5 = r8_uniform_01 ( &seed );
  d5 = r8_nint ( 50.0 * d5 - 25.0 ) / 5.0;
  a = toeplitz_5diag ( n, d1, d2, d3, d4, d5 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TOEPLITZ_5DIAG      "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TOEPLITZ_5S matrix.
//
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta - 25.0 ) / 5.0;
  gamma = r8_uniform_01 ( &seed );
  gamma = r8_nint ( 50.0 * gamma - 25.0 ) / 5.0;
  a = toeplitz_5s ( row_num, col_num, alpha, beta, gamma, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TOEPLITZ_5S         "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TOEPLITZ_PDS matrix.
//
  m = 3;
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( m, &seed );
  y = r8vec_uniform_01_new ( m, &seed );
  y_sum = r8vec_sum ( m, y );
  for ( i = 0; i < m; i++ )
  {
    y[i] = y[i] / y_sum;
  }
  a = toeplitz_pds ( m, n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TOEPLITZ_PDS        "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  TOURNAMENT_RANDOM matrix.
//
  n = 5;
  seed_save = 123456789;
  seed = seed_save;
  a = tournament_random ( n, &seed );
  seed = seed_save;
  determ1 = tournament_random_determinant ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TOURNAMENT_RANDOM   "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRANSITION_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  a = transition_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TRANSITION_RANDOM   "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRENCH matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = trench ( alpha, n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TRENCH              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRI_UPPER matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = tri_upper ( alpha, n );
  determ1 = tri_upper_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TRI_UPPER           "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRIS matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta - 25.0 ) / 5.0;
  gamma = r8_uniform_01 ( &seed );
  gamma = r8_nint ( 50.0 * gamma - 25.0 ) / 5.0;
  a = tris ( n, n, alpha, beta, gamma );
  determ1 = tris_determinant ( n, alpha, beta, gamma );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TRIS                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRIV matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  z = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    z[i] = r8_nint ( 50.0 * z[i] - 25.0 ) / 5.0;
  }
  a = triv ( n, x, y, z );
  determ1 = triv_determinant ( n, x, y, z );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TRIV                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
  delete [] z;
//
//  TRIW matrix.
//
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 0, n - 1, &seed );
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = triw ( alpha, k, n );
  determ1 = triw_determinant ( alpha, k, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "TRIW                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  UPSHIFT matrix.
//
  n = 5;
  a = upshift ( n );
  determ1 = upshift_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "UPSHIFT             "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  VAND1 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = vand1 ( n, x );
  determ1 = vand1_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "VAND1               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  VAND2 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = vand2 ( n, x );
  determ1 = vand2_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "VAND2               "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  WATHEN matrix.
//
  if ( false )
  {
  row_num = 5;
  col_num = 5;
  n = wathen_order ( row_num, col_num );
  a = wathen ( row_num, col_num, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WATHEN             "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  }
  else
  {
  cout << "  " << setw(20) << "WATHEN             -----Not ready-----\n";
  }
//
//  WILK03 matrix.
//
  n = 3;
  a = wilk03 ( );
  determ1 = wilk03_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILK03              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK04 matrix.
//
  n = 4;
  a = wilk04 ( );
  determ1 = wilk04_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILK04              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK05 matrix.
//
  n = 5;
  a = wilk05 ( );
  determ1 = wilk05_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILK05              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK12 matrix.
//
  n = 12;
  a = wilk12 ( );
  determ1 = wilk12_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILK12              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK20 matrix.
//
  n = 20;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50 * alpha - 25.0 ) / 5.0;
  a = wilk20 ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILK20              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK21 matrix.
//
  n = 21;
  a = wilk21 ( n );
  determ1 = wilk21_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILK21              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILSON matrix.
//
  n = 4;
  a = wilson ( );
  determ1 = wilson_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILSON              "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ZERO matrix.
//
  n = 5;
  a = zero ( n, n );
  determ1 = zero_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ZERO                "
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ZIELKE matrix.
//
  n = 5;
  seed = 123456789;
  d1 = r8_uniform_01 ( &seed );
  d1 = r8_nint ( 50.0 * d1 - 25.0 ) / 5.0;
  d2 = r8_uniform_01 ( &seed );
  d2 = r8_nint ( 50.0 * d2 - 25.0 ) / 5.0;
  d3 = r8_uniform_01 ( &seed );
  d3 = r8_nint ( 50.0 * d3 - 25.0 ) / 5.0;
  a = zielke ( n, d1, d2, d3 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ZIELKE              "
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;

  return;
}
//****************************************************************************80

void test_eigen ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_EIGEN tests the eigenvalue computations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double beta;
  double error_frobenius;
  double gamma;
  int i;
  int i1;
  int k;
  double *lambda;
  int n;
  double norm_frobenius;
  int rank;
  int seed;
  int seed_save;
  string title;
  double *v1;
  double *v2;
  double *v3;
  double *x;

  cout << "\n";
  cout << "TEST_EIGEN\n";
  cout << "  Compute the Frobenius norm of the eigenvalue error:\n";
  cout << "    A * X - X * LAMBDA\n";
  cout << "  given a set of K eigenvectors X and eigenvalues LAMBDA.\n";
  cout << "\n";
  cout << "  Matrix title             N     K          ||A||       ||(A-Lambda*I)*X||\n";
  cout << "\n";
//
//  BODEWIG matrix.
//
  n = 4;
  k = 4;
  a = bodewig ( );
  lambda = bodewig_eigenvalues ( );
  x = bodewig_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "BODEWIG             "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  CARRY matrix.
//
  n = 5;
  k = 5;
  seed = 123456789;
  i1 = i4_uniform ( 2, 20, &seed );
  a = carry ( i1, n );
  lambda = carry_eigenvalues ( i1, n );
  x = carry_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CARRY               "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  CHOW matrix.
//
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = chow ( alpha, beta, n, n );
  lambda = chow_eigenvalues ( alpha, beta, n );
  x = chow_right ( alpha, beta, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "CHOW                "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  COMBIN matrix.
//
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = combin ( alpha, beta, n );
  lambda = combin_eigenvalues ( alpha, beta, n );
  x = combin_right ( alpha, beta, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "COMBIN              "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  DIF2 matrix.
//
  n = 5;
  k = 5;
  a = dif2 ( n, n );
  lambda = dif2_eigenvalues ( n );
  x = dif2_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "DIF2                "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  EXCHANGE matrix.
//
  n = 5;
  k = 5;
  a = exchange ( n, n );
  lambda = exchange_eigenvalues ( n );
  x = exchange_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "EXCHANGE            "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  IDEM_RANDOM matrix.
//
  n = 5;
  k = 5;
  rank = 3;
  seed_save = 987654321;
  seed = seed_save;
  a = idem_random ( n, rank, &seed );
  lambda = idem_random_eigenvalues ( n, rank );
  seed = seed_save;
  x = idem_random_right ( n, rank, &seed );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "IDEM_RANDOM         "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  IDENTITY matrix.
//
  n = 5;
  k = 5;
  a = identity ( n, n );
  lambda = identity_eigenvalues ( n );
  x = identity_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "IDENTITY            "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ILL3 matrix.
//
  n = 3;
  k = 3;
  a = ill3 ( );
  lambda = ill3_eigenvalues ( );
  x = ill3_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ILL3                "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  KMS matrix.
//
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  a = kms ( alpha, n, n );
  lambda = kms_eigenvalues ( alpha, n );
  x = kms_right ( alpha, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "KMS                 "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ONE matrix.
//
  n = 5;
  k = 5;
  a = one ( n, n );
  lambda = one_eigenvalues ( n );
  x = one_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ONE                 "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ORTEGA matrix.
//
  n = 5;
  k = 5;
  seed = 123456789;
  v1 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v1[i] = r8_nint ( 50.0 * v1[i] - 25.0 ) / 5.0;
  }
  v2 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v2[i] = r8_nint ( 50.0 * v2[i] - 25.0 ) / 5.0;
  }
  v3 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v3[i] = r8_nint ( 50.0 * v3[i] - 25.0 ) / 5.0;
  }
  a = ortega ( n, v1, v2, v3 );
  lambda = ortega_eigenvalues ( n, v1, v2, v3 );
  x = ortega_right ( n, v1, v2, v3 );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ORTEGA              "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] v1;
  delete [] v2;
  delete [] v3;
  delete [] x;
//
//  OTO matrix.
//
  n = 5;
  k = 5;
  a = oto ( n, n );
  lambda = oto_eigenvalues ( n );
  x = oto_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "OTO                 "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  PDS_RANDOM matrix.
//
  n = 5;
  k = 5;
  seed_save = 123456789;
  seed = seed_save;
  a = pds_random ( n, &seed );
  seed = seed_save;
  lambda = pds_random_eigenvalues ( n, &seed );
  seed = seed_save;
  x = pds_random_right ( n, &seed );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PDS_RANDOM          "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  PEI matrix.
//
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = pei ( alpha, n );
  lambda = pei_eigenvalues ( alpha, n );
  x = pei_right ( alpha, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "PEI                 "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  RODMAN matrix.
//
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = rodman ( alpha, n, n );
  lambda = rodman_eigenvalues ( alpha, n );
  x = rodman_right ( alpha, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RODMAN              "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ROSSER1 matrix.
//
  n = 8;
  k = 8;
  a = rosser1 ( );
  lambda = rosser1_eigenvalues ( );
  x = rosser1_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ROSSER1             "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  RUTIS1 matrix.
//
  n = 4;
  k = 4;
  a = rutis1 ( );
  lambda = rutis1_eigenvalues ( );
  x = rutis1_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RUTIS1              "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  RUTIS2 matrix.
//
  n = 4;
  k = 4;
  a = rutis2 ( );
  lambda = rutis2_eigenvalues ( );
  x = rutis2_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RUTIS2              "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  RUTIS5 matrix.
//
  n = 4;
  k = 4;
  a = rutis5 ( );
  lambda = rutis5_eigenvalues ( );
  x = rutis5_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "RUTIS5              "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  WILK12 matrix.
//
  n = 12;
  k = 12;
  a = wilk12 ( );
  lambda = wilk12_eigenvalues ( );
  x = wilk12_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILK12              "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  WILSON matrix.
//
  n = 4;
  k = 4;
  a = wilson ( );
  lambda = wilson_eigenvalues ( );
  x = wilson_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "WILSON              "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ZERO matrix.
//
  n = 5;
  k = 5;
  a = zero ( n, n );
  lambda = zero_eigenvalues ( n );
  x = zero_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << "ZERO                "
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;

  return;
}
//****************************************************************************80

void test_inverse ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_INVERSE tests the inverse computations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double *b;;
  double beta;
  double *c;
  double error_ab;
  double error_ac;;
  double gamma;
  int i;
  int ii;
  int jj;
  int k;
  double *l;
  int n;
  double norma_frobenius;
  double normc_frobenius;
  double *p;
  int *pivot;
  int seed;
  int seed_save;
  string title;
  double *u;
  double *v1;
  double *v2;
  double *v3;
  double *w;
  double *x;
  int x_n;
  double *y;
  int y_n;
  double *z;

  cout << "\n";
  cout << "TEST_INVERSE\n";
  cout << "  A = a test matrix of order N;\n";
  cout << "  B = inverse as computed by a routine.\n";
  cout << "  C = inverse as computed by R8MAT_INVERSE.\n";
  cout << "\n";
  cout << "  ||A||    = Frobenius norm of A.\n";
  cout << "  ||C||    = Frobenius norm of C.\n";
  cout << "  ||I-AC|| = Frobenius norm of I-A*C.\n";
  cout << "  ||I-AB|| = Frobenius norm of I-A*B.\n";
  cout << "\n";
  cout << "  Matrix title             N        "
       << "   ||A||          ||C||      ||I-AC||        ||I-AB||\n";
  cout << "\n";
//
//  AEGERTER matrix.
//
  n = 5;
  a = aegerter ( n );
  b = aegerter_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "AEGERTER            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BAB matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  a = bab ( n, alpha, beta );
  b = bab_inverse ( n, alpha, beta );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "BAB                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BERNSTEIN matrix.
//
  n = 5;
  a = bernstein ( n );
  b = bernstein_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "BERNSTEIN           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BIS matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  a = bis ( alpha, beta, n, n );
  b = bis_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "BIS                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BODEWIG matrix.
//
  n = 4;
  a = bodewig ( );
  b = bodewig_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "BODEWIG             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BOOTHROYD matrix.
//
  n = 5;
  a = boothroyd ( n );
  b = boothroyd_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "BOOTHROYD           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BORDERBAND matrix.
//
  n = 5;
  a = borderband ( n );
  b = borderband_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "BORDERBAND          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CARRY matrix.
//
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 2, 20, &seed );
  a = carry ( k, n );
  b = carry_inverse ( k, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CARRY               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CAUCHY matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  y = r8vec_uniform_01_new ( n, &seed );
  a = cauchy ( n, x, y );
  b = cauchy_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CAUCHY              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] y;
//
//  CHEBY_T matrix.
//
  n = 5;
  a = cheby_t ( n );
  b = cheby_t_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CHEBY_T             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CHEBY_U matrix.
//
  n = 5;
  a = cheby_u ( n );
  b = cheby_u_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CHEBY_U             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CHEBY_VAN2 matrix.
//
  n = 5;
  a = cheby_van2 ( n );
  b = cheby_van2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CHEBY_VAN2          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CHEBY_VAN3 matrix.
//
  n = 5;
  a = cheby_van3 ( n );
  b = cheby_van3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CHEBY_VAN3          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CHOW matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = chow ( alpha, beta, n, n );
  b = chow_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CHOW                "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CIRCULANT matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = circulant ( n, n, x );
  b = circulant_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CIRCULANT           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  CIRCULANT2 matrix.
//
  n = 5;
  a = circulant2 ( n );
  b = circulant2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CIRCULANT2          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CLEMENT1 matrix.
//
  n = 6;
  a = clement1 ( n );
  b = clement1_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CLEMENT1            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CLEMENT2 matrix.
//
  n = 6;
  a = clement2 ( n );
  b = clement2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CLEMENT2            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CLEMENT3.
//
  n = 6;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  a = clement3 ( n, x, y );
  b = clement3_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CLEMENT3            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] y;
//
//  COMBIN matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = combin ( alpha, beta, n );
  b = combin_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "COMBIN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  COMPANION.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  a = companion ( n, x );
  b = companion_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "COMPANION           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  COMPLEX_I
//
  n = 2;
  a = complex_i ( );
  b = complex_i_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "COMPLEX_I           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CONEX1 matrix.
//
  n = 4;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = conex1 ( alpha );
  b = conex1_inverse ( alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CONEX1              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CONEX2 matrix.
//
  n = 3;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = conex2 ( alpha );
  b = conex2_inverse ( alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CONEX2              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CONEX3 matrix.
//
  n = 5;
  a = conex3 ( n );
  b = conex3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CONEX3              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CONFERENCE matrix.
//
  n = 6;
  a = conference ( n );
  b = conference_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "CONFERENCE          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB2 matrix.
//
  n = 4;
  a = daub2 ( n );
  b = daub2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DAUB2               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB4 matrix.
//
  n = 8;
  a = daub4 ( n );
  b = daub4_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DAUB4               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB6 matrix.
//
  n = 12;
  a = daub6 ( n );
  b = daub6_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DAUB6               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB8 matrix.
//
  n = 16;
  a = daub8 ( n );
  b = daub8_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DAUB8               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB10 matrix.
//
  n = 20;
  a = daub10 ( n );
  b = daub10_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DAUB10              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB12 matrix.
//
  n = 24;
  a = daub12 ( n );
  b = daub12_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DAUB12              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DIAGONAL.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = diagonal ( n, n, x );
  b = diagonal_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DIAGONAL            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  DIF2 matrix.
//
  n = 5;
  a = dif2 ( n, n );
  b = dif2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DIF2                "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DOWNSHIFT matrix.
//
  n = 5;
  a = downshift ( n );
  b = downshift_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "DOWNSHIFT           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DRMAC
//

//
//  EULERIAN matrix.
//
  n = 5;
  a = eulerian ( n, n );
  b = eulerian_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "EULERIAN            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  EXCHANGE matrix.
//
  n = 5;
  a = exchange ( n, n );
  b = exchange_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "EXCHANGE            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FIBONACCI2 matrix.
//
  n = 5;
  a = fibonacci2 ( n );
  b = fibonacci2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "FIBONACCI2          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FIBONACCI3 matrix.
//
  n = 5;
  a = fibonacci3 ( n );
  b = fibonacci3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "FIBONACCI3          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FIEDLER.
//  The FIEDLER_INVERSE routine assumes the X vector is sorted.
//
  n = 7;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  r8vec_sort_bubble_a ( n, x );
  a = fiedler ( n, n, x );
  b = fiedler_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "FIEDLER             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  FORSYTHE matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = forsythe ( alpha, beta, n );
  b = forsythe_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "FORSYTHE            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FOURIER_COSINE matrix.
//
  n = 5;
  a = fourier_cosine ( n );
  b = fourier_cosine_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "FOURIER_COSINE      "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FOURIER_SINE matrix.
//
  n = 5;
  a = fourier_sine ( n );
  b = fourier_sine_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "FOURIER_SINE        "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FRANK matrix.
//
  n = 5;
  a = frank ( n );
  b = frank_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "FRANK               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GFPP matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  a = gfpp ( n, alpha );
  b = gfpp_inverse ( n, alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "GFPP                "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GIVENS matrix.
//
  n = 5;
  a = givens ( n, n );
  b = givens_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "GIVENS              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GK316 matrix.
//
  n = 5;
  a = gk316 ( n );
  b = gk316_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "GK316               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GK323 matrix.
//
  n = 5;
  a = gk323 ( n, n );
  b = gk323_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "GK323               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GK324 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = gk324 ( n, n, x );
  b = gk324_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "GK324               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  HANOWA matrix.
//
  n = 8;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  a = hanowa ( alpha, n );
  b = hanowa_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HANOWA              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HARMAN matrix.
//
  n = 8;
  a = harman (  );
  b = harman_inverse (  );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HARMAN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HARTLEY matrix.
//
  n = 5;
  a = hartley ( n );
  b = hartley_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HARTLEY             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HELMERT matrix.
//
  n = 5;
  a = helmert ( n );
  b = helmert_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HELMERT             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HELMERT2 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = helmert2 ( n, x );
  b = helmert2_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HELMERT2            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  HERMITE matrix.
//
  n = 5;
  a = hermite ( n );
  b = hermite_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HERMITE             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HERNDON matrix.
//
  n = 5;
  a = herndon ( n );
  b = herndon_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HERNDON             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HILBERT matrix.
//
  n = 5;
  a = hilbert ( n, n );
  b = hilbert_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HILBERT             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HOUSEHOLDER matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = householder ( n, x );
  b = householder_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "HOUSEHOLDER         "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  IDENTITY matrix.
//
  n = 5;
  a = identity ( n, n );
  b = identity_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "IDENTITY            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  ILL3 matrix.
//
  n = 3;
  a = ill3 ( );
  b = ill3_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "ILL3                "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  INTEGRATION matrix.
//
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = integration ( alpha, n );
  b = integration_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "INTEGRATION         "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  INVOL matrix.
//
  n = 5;
  a = invol ( n );
  b = invol_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "INVOL               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  JORDAN matrix.
//
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = jordan ( alpha, n, n );
  b = jordan_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "JORDAN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  KAHAN matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = kahan ( alpha, n, n );
  b = kahan_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "KAHAN               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  KERSHAW matrix.
//
  n = 4;
  a = kershaw ( );
  b = kershaw_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "KERSHAW             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  KERSHAWTRI matrix.
//
  n = 5;
  x_n = ( n + 1 ) / 2;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = kershawtri ( n, x );
  b = kershawtri_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "KERSHAWTRI          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  KMS matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = kms ( alpha, n, n );
  b = kms_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "KMS                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LAGUERRE matrix.
//
  n = 5;
  a = laguerre ( n );
  b = laguerre_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "LAGUERRE            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LEGENDRE matrix.
//
  n = 5;
  a = legendre ( n );
  b = legendre_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "LEGENDRE            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LEHMER matrix.
//
  n = 5;
  a = lehmer ( n, n );
  b = lehmer_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "LEHMER              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LIETZKE matrix.
//
  n = 5;
  a = lietzke ( n );
  b = lietzke_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "LIETZKE             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LOTKIN matrix.
//
  n = 5;
  a = lotkin ( n, n );
  b = lotkin_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "LOTKIN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  MAXIJ matrix.
//
  n = 5;
  a = maxij ( n, n );
  b = maxij_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "MAXIJ               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  MILNES matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = milnes ( n, n, x );
  b = milnes_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "MILNES              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  MINIJ matrix.
//
  n = 5;
  a = minij ( n, n );
  b = minij_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "MINIJ               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  MOLER1 matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = moler1 ( alpha, n, n );
  b = moler1_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "MOLER1              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  MOLER3 matrix.
//
  n = 5;
  a = moler3 ( n, n );
  b = moler3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "MOLER3              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  ORTEGA matrix.
//
  n = 5;
  seed = 123456789;
  v1 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v1[i] = r8_nint ( 50.0 * v1[i] - 25.0  ) / 5.0;
  }
  v2 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v2[i] = r8_nint ( 50.0 * v2[i] - 25.0  ) / 5.0;
  }
  v3 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v3[i] = r8_nint ( 50.0 * v3[i] - 25.0  ) / 5.0;
  }
  a = ortega ( n, v1, v2, v3 );
  b = ortega_inverse ( n, v1, v2, v3 );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "ORTEGA              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] v1;
  delete [] v2;
  delete [] v3;
//
//  ORTH_SYMM matrix.
//
  n = 5;
  a = orth_symm ( n );
  b = orth_symm_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "ORTH_SYMM           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  OTO matrix.
//
  n = 5;
  a = oto ( n, n );
  b = oto_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "OTO                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PARTER matrix.
//
  n = 5;
  a = parter ( n, n );
  b = parter_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "PARTER              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PASCAL1 matrix.
//
  n = 5;
  a = pascal1 ( n );
  b = pascal1_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "PASCAL1             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PASCAL2 matrix.
//
  n = 5;
  a = pascal2 ( n );
  b = pascal2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "PASCAL2             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PASCAL3 matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = pascal3 ( n, alpha );
  b = pascal3_inverse ( n, alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "PASCAL3             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PDS_RANDOM matrix.
//
  n = 5;
  seed_save = 123456789;
  seed = seed_save;
  a = pds_random ( n, &seed );
  seed = seed_save;
  b = pds_random_inverse ( n, &seed );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "PDS_RANDOM          "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PEI matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = pei ( alpha, n );
  b = pei_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "PEI                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PERMUTATION_RANDOM matrix.
//
  n = 5;
  seed = 123456789;
  seed_save = seed;
  a = permutation_random ( n, &seed );
  seed = seed_save;
  b = permutation_random_inverse ( n, &seed );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "PERMUTATION_RANDOM  "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PLU matrix.
//
  n = 5;
  pivot = new int[n];
  p = new double[n*n];
  l = new double[n*n];
  u = new double[n*n];
  for ( i = 0; i < n; i++ )
  {
    pivot[i] = i + 1;
  }
  a = plu ( n, pivot, p, l, u );
  b = plu_inverse ( n, p, l, u );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "PLU                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] l;
  delete [] p;
  delete [] pivot;
  delete [] u;
//
//  RIS matrix.
//
  n = 5;
  a = ris ( n );
  b = ris_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "RIS                 "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RODMAN matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = rodman ( alpha, n, n );
  b = rodman_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "RODMAN              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS1 matrix.
//
  n = 4;
  a = rutis1 ( );
  b = rutis1_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "RUTIS1              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS2 matrix.
//
  n = 4;
  a = rutis2 ( );
  b = rutis2_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "RUTIS2              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS3 matrix.
//
  n = 4;
  a = rutis3 ( );
  b = rutis3_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "RUTIS3              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS4 matrix.
//
  n = 5;
  a = rutis4 ( n );
  b = rutis4_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "RUTIS4              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS5 matrix.
//
  n = 4;
  a = rutis5 ( );
  b = rutis5_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "RUTIS5              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SCHUR_BLOCK matrix.
//
  n = 5;
  x_n = ( n + 1 ) / 2;
  y_n = n / 2;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( y_n, &seed );
  for ( i = 0; i < y_n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  a = schur_block ( n, x, y );
  b = schur_block_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "SCHUR_BLOCK         "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] y;
//
//  SPLINE matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = spline ( n, x );
  b = spline_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "SPLINE              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  STIRLING matrix.
//
  n = 5;
  a = stirling ( n, n );
  b = stirling_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "STIRLING            "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SUMMATION matrix.
//
  n = 5;
  a = summation ( n, n );
  b = summation_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "SUMMATION           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  TRI_UPPER matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = tri_upper ( alpha, n );
  b = tri_upper_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "TRI_UPPER           "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  TRIS matrix.
//
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta - 25.0 ) / 5.0;
  gamma = r8_uniform_01 ( &seed );
  gamma = r8_nint ( 50.0 * gamma - 25.0 ) / 5.0;
  a = tris ( n, n, alpha, beta, gamma );
  b = tris_inverse ( n, alpha, beta, gamma );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "TRIS                "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  TRIV matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  z = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    z[i] = r8_nint ( 50.0 * z[i] - 25.0 ) / 5.0;
  }
  a = triv ( n, x, y, z );
  b = triv_inverse ( n, x, y, z );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "TRIV                "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] y;
  delete [] z;
//
//  TRIW matrix.
//
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 0, n - 1, &seed );
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = triw ( alpha, k, n );
  b = triw_inverse ( alpha, k, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "TRIW                "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  UPSHIFT matrix.
//
  n = 5;
  a = upshift ( n );
  b = upshift_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "UPSHIFT             "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  VAND1 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = vand1 ( n, x );
  b = vand1_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "VAND1               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  VAND2 matrix.
//
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = vand2 ( n, x );
  b = vand2_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "VAND2               "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILK03 matrix.
//
  n = 3;
  a = wilk03 ( );
  b = wilk03_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "WILK03              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILK04 matrix.
//
  n = 4;
  a = wilk04 ( );
  b = wilk04_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "WILK04              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILK05 matrix.
//
  n = 5;
  a = wilk05 ( );
  b = wilk05_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "WILK05              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILK21 matrix.
//
  n = 21;
  a = wilk21 ( n );
  b = wilk21_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "WILK21              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILSON matrix.
//
  n = 4;
  a = wilson ( );
  b = wilson_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << "WILSON              "
       << "  " << setw(4) << n
       << "  " << setw(14) << norma_frobenius
       << "  " << setw(14) << normc_frobenius
       << "  " << setw(14) << error_ac
       << "  " << setw(14) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void test_null ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_NULL tests the null vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *at;
  double alpha;
  int col_num;
  double error_l2;
  double f1;
  double f2;
  int m;
  int mt;
  int n;
  int nt;
  double norm_a_frobenius;
  double norm_x_l2;
  int row_num;
  int seed;
  string title;
  double *x;

  cout << "\n";
  cout << "TEST_NULL\n";
  cout << "  A = a test matrix of order M by N\n";
  cout << "  x = an N vector, candidate for a null vector.\n";
  cout << "\n";
  cout << "  ||A|| = Frobenius norm of A.\n";
  cout << "  ||x|| = L2 norm of x.\n";
  cout << "  ||A*x||/||x|| = L2 norm of A*x over L2 norm of x.\n";
  cout << "\n";
  cout << "  Matrix title	           M     N      "
       << "||A||            ||x||        ||A*x||/||x||\n";
  cout << "\n";
//
//  ARCHIMEDES matrix.
//
  m = 7;
  n = 8;
  a = archimedes ( );
  x = archimedes_null ( );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "ARCHIMEDES          "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  CHEBY_DIFF1 matrix.
//
  m = 5;
  n = 5;
  a = cheby_diff1 ( n );
  x = cheby_diff1_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "CHEBY_DIFF1         "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  CREATION matrix.
//
  m = 5;
  n = 5;
  a = creation ( m, n );
  x = creation_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "CREATION            "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF1 matrix.
//  Only has null vectors for N odd.
//
  m = 5;
  n = 5;
  a = dif1 ( m, n );
  x = dif1_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "DIF1                "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF1CYCLIC matrix.
//
  m = 5;
  n = 5;
  a = dif1cyclic ( n );
  x = dif1cyclic_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "DIF1CYCLIC          "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF2CYCLIC matrix.
//
  m = 5;
  n = 5;
  a = dif2cyclic ( n );
  x = dif2cyclic_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "DIF2CYCLIC          "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  EBERLEIN matrix.
//  We have a LEFT null vector.
//
  m = 5;
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = eberlein ( alpha, n );
  mt = n;
  nt = m;
  at = r8mat_transpose_new ( m, n, a );
  x = eberlein_null_left ( n );
  error_l2 = r8mat_is_null_vector ( mt, nt, at, x );
  norm_a_frobenius = r8mat_norm_fro ( mt, nt, at );
  norm_x_l2 = r8vec_norm_l2 ( nt, x );
  cout << "  " << setw(20) << "EBERLEIN (left)     "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] at;
  delete [] x;
//
//  FIBONACCI1 matrix.
//
  m = 5;
  n = 5;
  seed = 123456789;
  f1 = r8_uniform_01 ( &seed );
  f1 = r8_nint ( 50.0 * f1 - 25.0 ) / 5.0;
  f2 = r8_uniform_01 ( &seed );
  f2 = r8_nint ( 50.0 * f2 - 25.0 ) / 5.0;
  a = fibonacci1 ( n, f1, f2 );
  x = fibonacci1_null ( n, f1, f2 );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "FIBONACCI1          "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  LAUCHLI matrix.
//  We have a LEFT null vector of a RECTANGULAR matrix.
//
  m = 6;
  n = m - 1;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = lauchli ( alpha, m, n );
  mt = n;
  nt = m;
  at = r8mat_transpose_new ( m, n, a );
  x = lauchli_null_left ( alpha, m, n );
  error_l2 = r8mat_is_null_vector ( mt, nt, at, x );
  norm_a_frobenius = r8mat_norm_fro ( mt, nt, at );
  norm_x_l2 = r8vec_norm_l2 ( nt, x );
  cout << "  " << setw(20) << "LAUCHLI (left)      "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] at;
  delete [] x;
//
//  LINE_ADJ matrix.
//
  m = 7;
  n = 7;
  a = line_adj ( n );
  x = line_adj_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "LINE_ADJ            "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  MOLER2 matrix.
//
  m = 5;
  n = 5;
  a = moler2 ( );
  x = moler2_null ( );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "MOLER2              "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  NEUMANN matrix.
//
  row_num = 5;
  col_num = 5;
  m = row_num * col_num;
  n = row_num * col_num;
  a = neumann ( row_num, col_num );
  x = neumann_null ( row_num, col_num );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "NEUMANN             "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ONE matrix.
//
  m = 5;
  n = 5;
  a = one ( n, n );
  x = one_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "ONE                 "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  RING_ADJ matrix.
//  N must be a multiple of 4 for there to be a null vector.
//
  m = 12;
  n = 12;
  a = ring_adj ( n );
  x = ring_adj_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "RING_ADJ            "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ROSSER1 matrix.
//
  m = 8;
  n = 8;
  a = rosser1 ( );
  x = rosser1_null ( );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "ROSSER1             "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ZERO matrix.
//
  m = 5;
  n = 5;
  a = zero ( m, n );
  x = zero_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << "ZERO                "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;

  return;
}
//****************************************************************************80

void test_plu ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_PLU tests the PLU factors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double error_frobenius;
  double *l;
  int m;
  int n;
  double norm_a_frobenius;
  double *p;
  int seed;
  string title;
  double *u;

  cout << "\n";
  cout << "TEST_PLU\n";
  cout << "  A = a test matrix of order M by N\n";
  cout << "  P, L, U are the PLU factors.\n";
  cout << "\n";
  cout << "  ||A|| = Frobenius norm of A.\n";
  cout << "  ||A-PLU|| = Frobenius norm of A-P*L*U.\n";
  cout << "\n";
  cout << "  Matrix title	           M     N      ";
  cout << "||A||            ||A-PLU||\n";
  cout << "\n";
//
//  BODEWIG matrix.
//
  m = 4;
  n = 4;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = bodewig ( );
  bodewig_plu ( p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  BODEWIG             "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  BORDERBAND matrix.
//
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = borderband ( n );
  borderband_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  BORDERBAND          "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  DIF2 matrix.
//
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = dif2 ( m, n );
  dif2_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  DIF2                "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  GFPP matrix.
//
  m = 5;
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = gfpp ( n, alpha );
  gfpp_plu ( n, alpha, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  GFPP                "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  GIVENS matrix.
//
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = givens ( m, n );
  givens_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  GIVENS              "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  KMS matrix.
//
  m = 5;
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = kms ( alpha, m, n );
  kms_plu ( alpha, n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  KMS                 "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  MAXIJ matrix.
//
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = maxij ( m, n );
  maxij_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  MAXIJ               "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  MINIJ matrix.
//
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = minij ( m, n );
  minij_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  MINIJ               "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  MOLER1 matrix.
//
  m = 5;
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = moler1 ( alpha, m, n );
  moler1_plu ( alpha, n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  MOLER1              "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  MOLER3 matrix.
//
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = moler3 ( m, n );
  moler3_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  MOLER3              "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  OTO matrix.
//
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = oto ( m, n );
  oto_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  OTO                 "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  PASCAL2 matrix.
//
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = pascal2 ( n );
  pascal2_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  PASCAL2             "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  WILSON matrix.
//
  m = 4;
  n = 4;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = wilson ( );
  wilson_plu ( p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  WILSON              "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;

  return;
}
//****************************************************************************80

void test_solution ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_SOLUTION tests the linear solution computations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double *b;
  double beta;
  double error_frobenius;
  double gamma;
  int i1;
  int k;
  int m;
  int n;
  int ncol;
  double norm_frobenius;
  int nrow;
  int seed;
  int seed_save;
  string title;
  double *x;

  cout << "\n";
  cout << "TEST_SOLUTION\n";
  cout << "  Compute the Frobenius norm of the solution error:\n";
  cout << "    A * X - B\n";
  cout << "  given MxN matrix A, NxK solution X, MxK right hand side B.\n";
  cout << "\n";
  cout << "  Matrix title             M     N     K      ||A||         ||A*X-B||\n";
  cout << "\n";
//
//  BODEWIG matrix.
//
  m = 4;
  n = 4;
  k = 1;
  a = bodewig ( );
  b = bodewig_rhs ( );
  x = bodewig_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << "BODEWIG             "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  DIF2 matrix.
//
  m = 10;
  n = 10;
  k = 2;
  a = dif2 ( m, n );
  b = dif2_rhs ( m, k );
  x = dif2_solution ( n, k );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << "DIF2                "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  FRANK matrix.
//
  m = 10;
  n = 10;
  k = 2;
  a = frank ( n );
  b = frank_rhs ( m, k );
  x = frank_solution ( n, k );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << "FRANK               "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  POISSON matrix.
//
  nrow = 4;
  ncol = 5;
  m = nrow * ncol;
  n = nrow * ncol;
  k = 1;
  a = poisson ( nrow, ncol, n );
  b = poisson_rhs ( nrow, ncol, n );
  x = poisson_solution ( nrow, ncol, n );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << "POISSON             "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  WILK03 matrix.
//
  m = 3;
  n = 3;
  k = 1;
  a = wilk03 ( );
  b = wilk03_rhs ( );
  x = wilk03_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << "WILK03              "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  WILK04 matrix.
//
  m = 4;
  n = 4;
  k = 1;
  a = wilk04 ( );
  b = wilk04_rhs ( );
  x = wilk04_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << "WILK04              "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  WILSON matrix.
//
  m = 4;
  n = 4;
  k = 1;
  a = wilson ( );
  b = wilson_rhs ( );
  x = wilson_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << "WILSON              "
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;

  return;
}
