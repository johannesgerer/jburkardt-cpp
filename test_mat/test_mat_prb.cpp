# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <cstring>

using namespace std;

# include "test_mat.hpp"

int main ( );
void bvec_next_grlex_test ( );
void legendre_zeros_test ( );
void mertens_test ( );
void moebius_test ( );
void r8mat_is_eigen_left_test ( );
void r8mat_is_llt_test ( );
void r8mat_is_null_left_test ( );
void r8mat_is_null_right_test ( );
void r8mat_is_solution_test ( );
void r8mat_norm_fro_test ( );
void test_condition ( );
void test_determinant ( );
void test_eigen_left ( );
void test_eigen_right ( );
void test_inverse ( );
void test_llt ( );
void test_null_left ( );
void test_null_right ( );
void test_plu ( );
void test_solution ( );
void test_type ( );

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
//    TEST_MAT_PRB tests the TEST_MAT library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2015
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
//
//  Utilities.
//
  bvec_next_grlex_test ( );
  legendre_zeros_test ( );
  mertens_test ( );
  moebius_test ( );
  r8mat_is_eigen_left_test ( );
  r8mat_is_llt_test ( );
  r8mat_is_null_left_test ( );
  r8mat_is_null_right_test ( );
  r8mat_is_solution_test ( );
  r8mat_norm_fro_test ( );
//
//  Interesting stuff.
//
  test_condition ( );
  test_determinant ( );
  test_eigen_left ( );
  test_eigen_right ( );
  test_inverse ( );
  test_llt ( );
  test_null_left ( );
  test_null_right ( );
  test_plu ( );
  test_solution ( );
  test_type ( );
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

void bvec_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_NEXT_GRLEX_TEST tests BVEC_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2015
//
//  Author:
//
//    John Burkardt
//
{ 
  int *b;
  int i;
  int j;
  int n = 4;

  cout << "\n";
  cout << "BVEC_NEXT_GRLEX_TEST\n";
  cout << "  BVEC_NEXT_GRLEX computes binary vectors in GRLEX order.\n";
  cout << "\n";

  b = new int[n];

  for ( j = 0; j < n; j++ )
  {
    b[j] = 0;
  }

  for ( i = 0; i <= 16; i++ )
  {
    cout << "  " << setw(2) << i << ":  ";
    for ( j = 0; j < n; j++ )
    {
      cout << b[j];
    }
    cout << "\n";
    bvec_next_grlex ( n, b );
  }

  delete [] b;

  return;
}
//****************************************************************************80

void legendre_zeros_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ZEROS_TEST tests LEGENDRE_ZEROS.
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
  double *l;
  int n;

  cout << "\n";
  cout << "LEGENDRE_ZEROS_TEST:\n";
  cout << "  LEGENDRE_ZEROS computes the zeros of the N-th Legendre\n";
  cout << "  polynomial.\n";

  for ( n = 1; n <= 7; n++ )
  {
    l = legendre_zeros ( n );
    r8vec_print ( n, l, "  Legendre zeros" );
    delete [] l;
  }

  return;
}
//****************************************************************************80

void mertens_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MERTENS_TEST tests MERTENS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "MERTENS_TEST\n";
  cout << "  MERTENS computes the Mertens function.\n";
  cout << "\n";
  cout << "      N   Exact   MERTENS(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
     mertens_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                              << "  "
         << setw(8)  << n             << "  "
         << setw(10) << c             << "  "
         << setw(10) << mertens ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void moebius_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MOEBIUS_TEST tests MOEBIUS.
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
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "MOEBIUS_TEST\n";
  cout << "  MOEBIUS computes the Moebius function.\n";
  cout << "\n";
  cout << "      N   Exact   MOEBIUS(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
     moebius_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                              << "  "
         << setw(8)  << n             << "  "
         << setw(10) << c             << "  "
         << setw(10) << moebius ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void r8mat_is_eigen_left_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_EIGEN_LEFT_TEST tests R8MAT_IS_EIGEN_LEFT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 March 2015
//
//  Author:
//
//    John Burkardt
//
{
//
//  This is the CARRY ( 4.0, 4 ) matrix.
//
  double a[4*4] = {
   0.13671875,   0.05859375,   0.01953125,   0.00390625, 
   0.60546875,   0.52734375,   0.39453125,   0.25390625, 
   0.25390625,   0.39453125,   0.52734375,   0.60546875, 
   0.00390625,   0.01953125,   0.05859375,   0.13671875 };
  int k = 4;
  double lam[4] = {
     1.000000000000000, 
     0.250000000000000, 
     0.062500000000000, 
     0.015625000000000 };
  int n = 4;
  double value;
  double x[4*4] = {
       1.0, 11.0, 11.0,  1.0, 
       1.0,  3.0, -3.0, -1.0, 
       1.0, -1.0, -1.0,  1.0, 
       1.0, -3.0,  3.0, -1.0 };

  cout << "\n";
  cout << "R8MAT_IS_EIGEN_LEFT_TEST:\n";
  cout << "  R8MAT_IS_EIGEN_LEFT tests the error in the left eigensystem\n";
  cout << "    A' * X - X * LAMBDA = 0\n";

  r8mat_print ( n, n, a, "  Matrix A:" );
  r8mat_print ( n, k, x, "  Eigenmatrix X:" );
  r8vec_print ( n, lam, "  Eigenvalues LAM:" );

  value = r8mat_is_eigen_left ( n, k, a, x, lam );

  cout << "\n";
  cout << "  Frobenius norm of A'*X-X*LAMBDA is " << value << "\n";

  cout << "\n";
  cout << "R8MAT_IS_EIGEN_LEFT_TEST\n";
  cout << "  Normal end of execution.\n";

  return;
}
//****************************************************************************80

void r8mat_is_llt_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_LLT_TEST tests R8MAT_IS_LLT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 March 2015
//
//  Author:
//
//    John Burkardt
//
{
//
//  Matrix is listed by columns.
//
  double a[4*4] = { 
    2.0, 1.0, 0.0, 0.0, 
    1.0, 2.0, 1.0, 0.0, 
    0.0, 1.0, 2.0, 1.0, 
    0.0, 0.0, 1.0, 2.0 };
  double enorm;
  double l[4*4] = { 
    1.414213562373095, 0.707106781186547, 
    0.0,               0.0,               
    0.0,               1.224744871391589, 
    0.816496580927726, 0.0,               
    0.0,               0.0,               
    1.154700538379251, 0.866025403784439, 
    0.0,               0.0,               
    0.0,               1.118033988749895  };
  int m = 4;
  int n = 4;

  cout << "\n";
  cout << "R8MAT_IS_LLT_TEST:\n";
  cout << "  R8MAT_IS_LLT tests the error in a lower triangular\n";
  cout << "  Cholesky factorization A = L * L' by looking at A-L*L'\n";

  r8mat_print ( m, m, a, "  Matrix A:" );
  r8mat_print ( m, n, l, "  Cholesky factor L:" );

  enorm = r8mat_is_llt ( m, n, a, l );

  cout << "\n";
  cout << "  Frobenius norm of A-L*L' is " << enorm << "\n";

  cout << "\n";
  cout << "R8MAT_IS_LLT_TEST\n";
  cout << "  Normal end of execution.\n";

  return;
}
//****************************************************************************80

void r8mat_is_null_left_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_NULL_LEFT_TEST tests R8MAT_IS_NULL_LEFT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
//
//  Matrix is listed by columns.
//
  double a[3*3] = { 
    1.0, 4.0, 7.0, 
    2.0, 5.0, 8.0,
    3.0, 6.0, 9.0 };
  double enorm;
  int m = 3;
  int n = 3;
  double x[3] = { 1.0, -2.0, 1.0 };

  cout << "\n";
  cout << "R8MAT_IS_NULL_LEFT_TEST:\n";
  cout << "  R8MAT_IS_NULL_LEFT tests whether the M vector X\n";
  cout << "  is a left null vector of A, that is, x'*A=0.\n";

  r8mat_print ( m, n, a, "  Matrix A:" );
  r8vec_print ( m, x, "  Vector X:" );

  enorm = r8mat_is_null_left ( m, n, a, x );

  cout << "\n";
  cout << "  Frobenius norm of X'*A is " << enorm << "\n";

  return;
}
//****************************************************************************80

void r8mat_is_null_right_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_NULL_RIGHT_TEST tests R8MAT_IS_NULL_RIGHT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
//
//  Matrix is listed by columns.
//
  double a[3*3] = { 
    1.0, 4.0, 7.0, 
    2.0, 5.0, 8.0,
    3.0, 6.0, 9.0 };
  double enorm;
  int m = 3;
  int n = 3;
  double x[3] = { 1.0, -2.0, 1.0 };

  cout << "\n";
  cout << "R8MAT_IS_NULL_RIGHT_TEST:\n";
  cout << "  R8MAT_IS_NULL_RIGHT tests whether the N vector X\n";
  cout << "  is a right null vector of A, that is, A*x=0.\n";

  r8mat_print ( m, n, a, "  Matrix A:" );
  r8vec_print ( n, x, "  Vector X:" );

  enorm = r8mat_is_null_right ( m, n, a, x );

  cout << "\n";
  cout << "  Frobenius norm of A*x is " << enorm << "\n";

  return;
}
//****************************************************************************80

void r8mat_is_solution_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_SOLUTION_TEST tests R8MAT_IS_SOLUTION.
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
  double *a;
  double *b;
  double enorm;
  int i4_hi;
  int i4_lo;
  int k;
  int m;
  int n;
  double r8_hi;
  double r8_lo;
  int seed;
  double *x;

  cout << "\n";
  cout << "R8MAT_IS_SOLUTION_TEST:\n";
  cout << "  R8MAT_IS_SOLUTION tests whether X is the solution of\n";
  cout << "  A*X=B by computing the Frobenius norm of the residual.\n";
//
//  Get random shapes.
//
  i4_lo = 1;
  i4_hi = 10;
  seed = 123456789;
  m = i4_uniform_ab ( i4_lo, i4_hi, seed );
  n = i4_uniform_ab ( i4_lo, i4_hi, seed );
  k = i4_uniform_ab ( i4_lo, i4_hi, seed );
//
//  Get a random A.
//
  r8_lo = -5.0;
  r8_hi = +5.0;
  a = r8mat_uniform_ab_new ( m, n, r8_lo, r8_hi, seed );
//
//  Get a random X.
//
  r8_lo = -5.0;
  r8_hi = +5.0;
  x = r8mat_uniform_ab_new ( n, k, r8_lo, r8_hi, seed );
//
//  Compute B = A * X.
//
  b = r8mat_mm_new ( m, n, k, a, x );
//
//  Compute || A*X-B||
//
  enorm = r8mat_is_solution ( m, n, k, a, x, b );
  
  cout << "\n";
  cout << "  A is " << m << " by " << n << "\n";
  cout << "  X is " << n << " by " << k << "\n";
  cout << "  B is " << m << " by " << k << "\n";
  cout << "  Frobenius error in A*X-B is " << enorm << "\n";
//
//  Free memory.
//
  delete [] a;
  delete [] b;
  delete [] x;

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

void test_condition ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_CONDITION tests the condition number computations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double a_norm;
  double alpha;
  double *b;
  double b_norm;
  double beta;
  double cond1;
  double cond2;
  int i;
  int n;
  double r8_hi;
  double r8_lo;
  int seed;
  string title;
  double *x;

  cout << "\n";
  cout << "TEST_CONDITION\n";
  cout << "  Compute the L1 condition number of an example of each\n";
  cout << "  test matrix\n";
  cout << "\n";
  cout << "  Title                    N            COND            COND\n";
  cout << "\n";
//
//  AEGERTER
//
  title = "AEGERTER";
  n = 5;
  cond1 = aegerter_condition ( n );

  a = aegerter ( n );
  b = aegerter_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  BAB
//
  title = "BAB";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = bab_condition ( n, alpha, beta );

  a = bab ( n, alpha, beta );
  b = bab_inverse ( n, alpha, beta );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  BAUER
//
  title = "BAUER";
  n = 6;
  cond1 = bauer_condition ( );

  a = bauer ( );
  b = bauer_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  BIS
//
  title = "BIS";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = bis_condition ( alpha, beta, n );

  a = bis ( alpha, beta, n, n );
  b = bis_inverse ( alpha, beta, n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  BIW
//
  title = "BIW";
  n = 5;
  cond1 = biw_condition ( n );

  a = biw ( n );
  b = biw_inverse (  n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  BODEWIG
//
  title = "BODEWIG";
  n = 4;
  cond1 = bodewig_condition ( );

  a = bodewig ( );
  b = bodewig_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  BOOTHROYD
//
  title = "BOOTHROYD";
  n = 5;
  cond1 = boothroyd_condition ( n );

  a = boothroyd ( n );
  b = boothroyd_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  COMBIN
//
  title = "COMBIN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = combin_condition ( alpha, beta, n );

  a = combin ( alpha, beta, n );
  b = combin_inverse ( alpha, beta, n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  COMPANION
//
  title = "COMPANION";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  cond1 = companion_condition ( n, x );

  a = companion ( n, x );
  b = companion_inverse ( n, x );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  CONEX1
//
  title = "CONEX1";
  n = 4;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = conex1_condition ( alpha );

  a = conex1 ( alpha );
  b = conex1_inverse ( alpha );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  CONEX2
//
  title = "CONEX2";
  n = 3;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = conex2_condition ( alpha );

  a = conex2 ( alpha );
  b = conex2_inverse ( alpha );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  CONEX3
//
  title = "CONEX3";
  n = 5;
  cond1 = conex3_condition ( n );

  a = conex3 ( n );
  b = conex3_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  CONEX4
//
  title = "CONEX4";
  n = 4;
  cond1 = conex4_condition ( );

  a = conex4 ( );
  b = conex4_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  DAUB2
//
  title = "DAUB2";
  n = 4;
  cond1 = daub2_condition ( n );

  a = daub2 ( n );
  b = daub2_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  DAUB4
//
  title = "DAUB4";
  n = 8;
  cond1 = daub4_condition ( n );

  a = daub4 ( n );
  b = daub4_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  DAUB6
//
  title = "DAUB6";
  n = 12;
  cond1 = daub6_condition ( n );

  a = daub6 ( n );
  b = daub6_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  DAUB8
//
  title = "DAU8";
  n = 16;
  cond1 = daub8_condition ( n );

  a = daub8 ( n );
  b = daub8_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  DAUB10
//
  title = "DAUB10";
  n = 20;
  cond1 = daub10_condition ( n );

  a = daub10 ( n );
  b = daub10_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  DAUB12
//
  title = "DAUB12";
  n = 24;
  cond1 = daub12_condition ( n );

  a = daub12 ( n );
  b = daub12_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  DIAGONAL
//
  title = "DIAGONAL";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  cond1 = diagonal_condition ( n, x );

  a = diagonal ( n, n, x );
  b = diagonal_inverse ( n, x );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  DIF2
//
  title = "DIF2";
  n = 5;
  cond1 = dif2_condition ( n );

  a = dif2 ( n, n );
  b = dif2_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  DOWNSHIFT
//
  title = "DOWNSHIFT";
  n = 5;
  cond1 = downshift_condition ( n );

  a = downshift ( n );
  b = downshift_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  EXCHANGE
//
  title = "EXCHANGE";
  n = 5;
  cond1 = exchange_condition ( n );

  a = exchange ( n, n );
  b = exchange_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  FIBONACCI2
//
  title = "FIBONACCI2";
  n = 5;
  cond1 = fibonacci2_condition ( n );

  a = fibonacci2 ( n );
  b = fibonacci2_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  GFPP
//
  title = "GFPP";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = gfpp_condition ( n, alpha );

  a = gfpp ( n, alpha );
  b = gfpp_inverse ( n, alpha );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  GIVENS
//
  title = "GIVENS";
  n = 5;
  cond1 = givens_condition ( n );

  a = givens ( n, n );
  b = givens_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  HANKEL_N
//
  title = "HANKEL_N";
  n = 5;
  cond1 = hankel_n_condition ( n );

  a = hankel_n ( n );
  b = hankel_n_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  HARMAN
//
  title = "HARMAN";
  n = 8;
  cond1 = harman_condition ( );

  a = harman ( );
  b = harman_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  HARTLEY
//
  title = "HARTLEY";
  n = 5;
  cond1 = hartley_condition ( n );

  a = hartley ( n );
  b = hartley_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  IDENTITY
//
  title = "IDENTITY";
  n = 5;
  cond1 = identity_condition ( n );

  a = identity ( n, n );
  b = identity_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  ILL3
//
  title = "ILL3";
  n = 3;
  cond1 = ill3_condition ( );

  a = ill3 ( );
  b = ill3_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  JORDAN
//
  title = "JORDAN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = jordan_condition ( n, alpha );

  a = jordan ( n, n, alpha );
  b = jordan_inverse ( n, alpha );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  KERSHAW
//
  title = "KERSHAW";
  n = 4;
  cond1 = kershaw_condition ( );

  a = kershaw ( );
  b = kershaw_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  LIETZKE
//
  title = "LIETZKE";
  n = 5;
  cond1 = lietzke_condition ( n );

  a = lietzke ( n );
  b = lietzke_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  MAXIJ
//
  title = "MAXIJ";
  n = 5;
  cond1 = maxij_condition ( n );

  a = maxij ( n, n );
  b = maxij_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  MINIJ
//
  title = "MINIJ";
  n = 5;
  cond1 = minij_condition ( n );

  a = minij ( n, n );
  b = minij_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  ORTH_SYMM
//
  title = "ORTH_SYMM";
  n = 5;
  cond1 = orth_symm_condition ( n );

  a = orth_symm ( n );
  b = orth_symm_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  OTO
//
  title = "OTO";
  n = 5;
  cond1 = oto_condition ( n );

  a = oto ( n, n );
  b = oto_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  PASCAL1
//
  title = "PASCAL1";
  n = 5;
  cond1 = pascal1_condition ( n );

  a = pascal1 ( n );
  b = pascal1_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  PASCAL3
//
  title = "PASCAL3";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = pascal3_condition ( n, alpha );

  a = pascal3 ( n, alpha );
  b = pascal3_inverse ( n, alpha );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  PEI
//
  title = "PEI";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = pei_condition ( alpha, n );

  a = pei ( alpha, n );
  b = pei_inverse ( alpha, n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  RODMAN
//
  title = "RODMAN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = rodman_condition ( n, alpha );

  a = rodman ( n, n, alpha );
  b = rodman_inverse ( n, alpha );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  RUTIS1
//
  title = "RUTIS1";
  n = 4;
  cond1 = rutis1_condition ( );

  a = rutis1 ( );
  b = rutis1_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  RUTIS2
//
  title = "RUTIS2";
  n = 4;
  cond1 = rutis2_condition ( );

  a = rutis2 ( );
  b = rutis2_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  RUTIS3
//
  title = "RUTIS3";
  n = 4;
  cond1 = rutis3_condition ( );

  a = rutis3 ( );
  b = rutis3_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  RUTIS5
//
  title = "RUTIS5";
  n = 4;
  cond1 = rutis5_condition ( );

  a = rutis5 ( );
  b = rutis5_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  SUMMATION
//
  title = "SUMMATION";
  n = 5;
  cond1 = summation_condition ( n );

  a = summation ( n, n );
  b = summation_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  SWEET1
//
  title = "SWEET1";
  n = 6;
  cond1 = sweet1_condition ( );

  a = sweet1 ( );
  b = sweet1_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  SWEET2
//
  title = "SWEET2";
  n = 6;
  cond1 = sweet2_condition ( );

  a = sweet2 ( );
  b = sweet2_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  SWEET3
//
  title = "SWEET3";
  n = 6;
  cond1 = sweet3_condition ( );

  a = sweet3 ( );
  b = sweet3_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  SWEET4
//
  title = "SWEET4";
  n = 13;
  cond1 = sweet4_condition ( );

  a = sweet4 ( );
  b = sweet4_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  TRI_UPPER
//
  title = "TRI_UPPER";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  cond1 = tri_upper_condition ( alpha, n );

  a = tri_upper ( alpha, n );
  b = tri_upper_inverse ( alpha, n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  UPSHIFT
//
  title = "UPSHIFT";
  n = 5;
  cond1 = upshift_condition ( n );

  a = upshift ( n );
  b = upshift_inverse ( n );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  WILK03
//
  title = "WILK03";
  n = 3;
  cond1 = wilk03_condition ( );

  a = wilk03 ( );
  b = wilk03_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  WILK04
//
  title = "WILK04";
  n = 4;
  cond1 = wilk04_condition ( );

  a = wilk04 ( );
  b = wilk04_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  WILK05
//
  title = "WILK05";
  n = 5;
  cond1 = wilk05_condition ( );

  a = wilk05 ( );
  b = wilk05_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;
//
//  WILSON
//
  title = "WILSON";
  n = 4;
  cond1 = wilson_condition ( );

  a = wilson ( );
  b = wilson_inverse ( );
  a_norm = r8mat_norm_l1 ( n, n, a );
  b_norm = r8mat_norm_l1 ( n, n, b );
  cond2 = a_norm * b_norm;

  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << cond1
       << "  " << setw(14) << cond2 << "\n";
  delete [] a;
  delete [] b;

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
//    13 April 2015
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
  double *d;
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
  int i4_hi;
  int i4_lo;
  int ii;
  int jj;
  int k;
  int key;
  double *l;
  int m;
  int n;
  double norm_frobenius;
  double *p;
  int *pivot;
  double prob;
  double r8_hi;
  double r8_lo;
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
  double x_hi;
  double x_lo;
  int x_n;
  double x1;
  double x2;
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
  cout << "  Title                    N          "
       << "Determ          Determ          ||A||\n";
  cout << "\n";
//
//  A123
//
  title = "A123";
  n = 3;
  a = a123 ( );
  determ1 = a123_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  AEGERTER
//
  title = "AEGERTER";
  n = 5;
  a = aegerter ( n );
  determ1 = aegerter_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ANTICIRCULANT
//
  title = "ANTICIRCULANT       ";
  n = 3;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  ANTICIRCULANT
//
  title = "ANTICIRCULANT       ";
  n = 4;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  ANTICIRCULANT
//
  title = "ANTICIRCULANT       ";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  ANTIHADAMARD
//
  title = "ANTIHADAMARD        ";
  n = 5;
  a = antihadamard ( n );
  determ1 = antihadamard_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ANTISYMM_RANDOM
//
  title = "ANTISYMM_RANDOM     ";
  n = 5;
  key = 123456789;
  a = antisymm_random ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ANTISYMM_RANDOM
//
  title = "ANTISYMM_RANDOM     ";
  n = 6;
  key = 123456789;
  a = antisymm_random ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BAB
//
  title = "BAB";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = bab ( n, alpha, beta );
  determ1 = bab_determinant ( n, alpha, beta );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BAUER
//
  title = "BAUER";
  n = 6;
  a = bauer ( );
  determ1 = bauer_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BERNSTEIN
//
  title = "BERNSTEIN";
  n = 5;
  a = bernstein ( n );
  determ1 = bernstein_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BIMARKOV_RANDOM
//
  title = "BIMARKOV_RANDOM";
  n = 5;
  key = 123456789;
  a = bimarkov_random ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BIS
//
  title = "BIS";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = bis ( alpha, beta, n, n );
  determ1 = bis_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BIW
//
  title = "BIW";
  n = 5;
  a = biw ( n );
  determ1 = biw_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BODEWIG
//
  title = "BODEWIG";
  n = 4;
  a = bodewig ( );
  determ1 = bodewig_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BOOTHROYD
//
  title = "BOOTHROYD";
  n = 5;
  a = boothroyd ( n );
  determ1 = boothroyd_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  BORDERBAND
//
  title = "BORDERBAND";
  n = 5;
  a = borderband ( n );
  determ1 = borderband_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CARRY
//
  title = "CARRY";
  n = 5;
  i4_lo = 2;
  i4_hi = 20;
  seed = 123456789;
  k = i4_uniform_ab ( i4_lo, i4_hi, seed );
  a = carry ( n, k );
  determ1 = carry_determinant ( n, k );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CAUCHY
//
  title = "CAUCHY";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = cauchy ( n, x, y );
  determ1 = cauchy_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  CHEBY_DIFF1
//
  n = 5;
  a = cheby_diff1 ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHEBY_DIFF1
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
//  CHEBY_T
//
  title = "CHEBY_T";
  n = 5;
  a = cheby_t ( n );
  determ1 = cheby_t_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHEBY_U
//
  title = "CHEBY_U";
  n = 5;
  a = cheby_u ( n );
  determ1 = cheby_u_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHEBY_VAN1
//
  n = 5;
  r8_lo = -1.0;
  r8_hi = +1.0;
  x = r8vec_linspace_new ( n, r8_lo, r8_hi );
  a = cheby_van1 ( n, r8_lo, r8_hi, n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  CHEBY_VAN2
//
  for ( n = 2; n <= 10; n++ )
  {
    title = "CHEBY_VAN2";
    a = cheby_van2 ( n );
    determ1 = cheby_van2_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << left << title << right
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
//
//  CHEBY_VAN3
//
  title = "CHEBY_VAN3";
  n = 5;
  a = cheby_van3 ( n );
  determ1 = cheby_van3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CHOW
//
  title = "CHOW";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = chow ( alpha, beta, n, n );
  determ1 = chow_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CIRCULANT
//
  title = "CIRCULANT";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = circulant ( n, n, x );
  determ1 = circulant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  CIRCULANT2
//
  title = "CIRCULANT2";
  n = 3;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CIRCULANT2
//
  title = "CIRCULANT2";
  n = 4;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CIRCULANT2
//
  title = "CIRCULANT2";
  n = 5;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CLEMENT1
//
  title = "CLEMENT1";
  n = 5;
  a = clement1 ( n );
  determ1 = clement1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CLEMENT1
//
  title = "CLEMENT1";
  n = 6;
  a = clement1 ( n );
  determ1 = clement1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CLEMENT2
//
  title = "CLEMENT2";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  a = clement2 ( n, x, y );
  determ1 = clement2_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  CLEMENT2
//
  title = "CLEMENT2";
  n = 6;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  a = clement2 ( n, x, y );
  determ1 = clement2_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  COMBIN
//
  title = "COMBIN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = combin ( alpha, beta, n );
  determ1 = combin_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  COMPANION
//
  title = "COMPANION";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = companion ( n, x );
  determ1 = companion_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  COMPLEX_I
//
  title = "COMPLEX_I";
  n = 2;
  a = complex_i ( );
  determ1 = complex_i_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONEX1
//
  title = "CONEX1";
  n = 4;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = conex1 ( alpha );
  determ1 = conex1_determinant ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONEX2
//
  title = "CONEX2";
  n = 3;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = conex2 ( alpha );
  determ1 = conex2_determinant ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONEX3
//
  title = "CONEX3";
  n = 5;
  a = conex3 ( n );
  determ1 = conex3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONEX4
//
  title = "CONEX4";
  n = 4;
  a = conex4 ( );
  determ1 = conex4_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CONFERENCE
//  N-1 must be an odd prime or a power of an odd prime.
//
  title = "CONFERENCE";
  n = 6;
  a = conference ( n );
  determ1 = conference_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  CREATION
//
  title = "CREATION";
  n = 5;
  a = creation ( n, n );
  determ1 = creation_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB2
//
  title = "DAUB2";
  n = 4;
  a = daub2 ( n );
  determ1 = daub2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB4
//
  title = "DAUB4";
  n = 8;
  a = daub4 ( n );
  determ1 = daub4_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB6
//
  title = "DAUB6";
  n = 12;
  a = daub6 ( n );
  determ1 = daub6_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB8
//
  title = "DAUB8";
  n = 16;
  a = daub8 ( n );
  determ1 = daub8_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB10
//
  title = "DAUB10";
  n = 20;
  a = daub10 ( n );
  determ1 = daub10_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DAUB12
//
  title = "DAUB12";
  n = 24;
  a = daub12 ( n );
  determ1 = daub12_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIAGONAL
//
  title = "DIAGONAL";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = diagonal ( n, n, x );
  determ1 = diagonal_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  DIF1
//
  title = "DIF1";
  n = 5;
  a = dif1 ( n, n );
  determ1 = dif1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIF1
//
  title = "DIF1";
  n = 6;
  a = dif1 ( n, n );
  determ1 = dif1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIF1CYCLIC
//
  title = "DIF1CYCLIC";
  n = 5;
  a = dif1cyclic ( n );
  determ1 = dif1cyclic_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIF2
//
  title = "DIF2";
  n = 5;
  a = dif2 ( n, n );
  determ1 = dif2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DIF2CYCLIC
//
  title = "DIF2CYCLIC";
  n = 5;
  a = dif2cyclic ( n );
  determ1 = dif2cyclic_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DORR
//
  title = "DORR";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = dorr ( alpha, n );
  determ1 = dorr_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  DOWNSHIFT
//
  title = "DOWNSHIFT";
  n = 5;
  a = downshift ( n );
  determ1 = downshift_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  EBERLEIN
//
  title = "EBERLEIN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = eberlein ( alpha, n );
  determ1 = eberlein_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  EULERIAN
//
  title = "EULERIAN";
  n = 5;
  a = eulerian ( n, n );
  determ1 = eulerian_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  EXCHANGE
//
  title = "EXCHANGE";
  n = 5;
  a = exchange ( n, n );
  determ1 = exchange_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FIBONACCI1
//
  title = "FIBONACCI1";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = fibonacci1 ( n, alpha, beta );
  determ1 = fibonacci1_determinant ( n, alpha, beta );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FIBONACCI2
//
  title = "FIBONACCI2";
  n = 5;
  a = fibonacci2 ( n );
  determ1 = fibonacci2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FIBONACCI3
//
  title = "FIBONACCI3";
  n = 5;
  a = fibonacci3 ( n );
  determ1 = fibonacci3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FIEDLER
//
  title = "FIEDLER";
  n = 7;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = fiedler ( n, n, x );
  determ1 = fiedler_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  FORSYTHE
//
  title = "FORSYTHE";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = forsythe ( alpha, beta, n );
  determ1 = forsythe_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FORSYTHE
//
  title = "FORSYTHE";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = forsythe ( alpha, beta, n );
  determ1 = forsythe_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FOURIER_COSINE
//
  title = "FOURIER_COSINE";
  n = 5;
  a = fourier_cosine ( n );
  determ1 = fourier_cosine_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FOURIER_SINE
//
  title = "FOURIER_SINE";
  n = 5;
  a = fourier_sine ( n );
  determ1 = fourier_sine_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  FRANK
//
  title = "FRANK";
  n = 5;
  a = frank ( n );
  determ1 = frank_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GEAR
//
  for ( n = 4; n <= 8; n++ )
  {
    title = "GRCAR";
    i4_lo = -n;
    i4_hi = +n;
    seed = 123456789;
    ii = i4_uniform_ab ( i4_lo, i4_hi, seed );
    jj = i4_uniform_ab ( i4_lo, i4_hi, seed );
    a = gear ( ii, jj, n );
    determ1 = gear_determinant ( ii, jj, n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << left << title << right
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
//
//  GFPP
//
  title = "GFPP";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = gfpp ( n, alpha );
  determ1 = gfpp_determinant ( n, alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GIVENS
//
  title = "GIVENS";
  n = 5;
  a = givens ( n, n );
  determ1 = givens_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GK316
//
  title = "GK316";
  n = 5;
  a = gk316 ( n );
  determ1 = gk316_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GK323
//
  title = "GK323";
  n = 5;
  a = gk323 ( n, n );
  determ1 = gk323_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  GK324
//
  title = "GK324";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  a = gk324 ( n, n, x );
  determ1 = gk324_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  GRCAR
//
  title = "GRCAR";
  n = 5;
  i4_lo = 1;
  i4_hi = n - 1;
  seed = 123456789;
  k = i4_uniform_ab ( i4_lo, i4_hi, seed );
  a = grcar ( n, n, k );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HADAMARD
//
  title = "HADAMARD";
  n = 5;
  a = hadamard ( n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HANKEL
//
  title = "HANKEL";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( 2 * n - 1, r8_lo, r8_hi, seed );
  a = hankel ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  HANKEL_N
//
  title = "HANKEL_N";
  n = 5;
  a = hankel_n ( n );
  determ1 = hankel_n_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HANOWA
//
  title = "HANOWA";
  n = 6;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = hanowa ( alpha, n );
  determ1 = hanowa_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HARMAN
//
  title = "HARMAN";
  n = 8;
  a = harman ( );
  determ1 = harman_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HARTLEY
//
  title = "HARTLEY";
  for ( n = 5; n <= 8; n++ )
  {
    a = hartley ( n );
    determ1 = hartley_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << left << title << right
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
//
//  HELMERT
//
  title = "HELMERT";
  n = 5;
  a = helmert ( n );
  determ1 = helmert_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HELMERT2
//
  title = "HELMERT2";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = helmert2 ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  HERMITE
//
  title = "HERMITE";
  n = 5;
  a = hermite ( n );
  determ1 = hermite_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HERNDON
//
  title = "HERNDON";
  n = 5;
  a = herndon ( n );
  determ1 = herndon_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HILBERT
//
  title = "HILBERT";
  n = 5;
  a = hilbert ( n, n );
  determ1 = hilbert_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  HOUSEHOLDER
//
  title = "HOUSEHOLDER";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = householder ( n, x );
  determ1 = householder_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  IDEM_RANDOM
//
  title = "IDEM_RANDOM";
  n = 5;
  i4_lo = 0;
  i4_hi = n;
  seed = 123456789;
  rank = i4_uniform_ab ( i4_lo, i4_hi, seed );
  key = 123456789;
  a = idem_random ( n, rank, key );
  determ1 = idem_random_determinant ( n, rank, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  IDENTITY
//
  title = "IDENTITY";
  n = 5;
  a = identity ( n, n );
  determ1 = identity_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  IJFACT1
//
  title = "IJFACT1";
  n = 5;
  a = ijfact1 ( n );
  determ1 = ijfact1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  IJFACT2
//
  title = "IJFACT2";
  n = 5;
  a = ijfact2 ( n );
  determ1 = ijfact2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ILL3
//
  title = "ILL3";
  n = 3;
  a = ill3 ( );
  determ1 = ill3_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  INTEGRATION
//
  title = "INTEGRATION";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = integration ( alpha, n );
  determ1 = integration_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  INVOL
//
  title = "INVOL";
  n = 5;
  a = invol ( n );
  determ1 = invol_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  INVOL_RANDOM
//
  title = "INVOL_RANDOM";
  n = 5;
  i4_lo = 0;
  i4_hi = n;
  seed = 123456789;
  rank = i4_uniform_ab ( i4_lo, i4_hi, seed );
  key = 123456789;
  a = invol_random ( n, rank, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  JACOBI
//
  title = "JACOBI";
  n = 5;
  a = jacobi ( n, n );
  determ1 = jacobi_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  JACOBI
//
  title = "JACOBI";
  n = 6;
  a = jacobi ( n, n );
  determ1 = jacobi_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  JORDAN
//
  title = "JORDAN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = jordan ( n, n, alpha );
  determ1 = jordan_determinant ( n, alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  KAHAN
//
  title = "KAHAN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = kahan ( alpha, n, n );
  determ1 = kahan_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  KERSHAW
//
  title = "KERSHAW";
  n = 4;
  a = kershaw ( );
  determ1 = kershaw_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  KERSHAWTRI
//
  title = "KERSHAW_TRI";
  n = 5;
  x_n = ( n + 1 ) / 2;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( x_n, r8_lo, r8_hi, seed );
  a = kershawtri ( n, x );
  determ1 = kershawtri_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  KMS
//
  title = "KMS";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = kms ( alpha, n, n );
  determ1 = kms_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LAGUERRE
//
  title = "LAGUERRE";
  n = 5;
  a = laguerre ( n );
  determ1 = laguerre_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LEGENDRE
//
  title = "LEGENDRE";
  n = 5;
  a = legendre ( n );
  determ1 = legendre_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LEHMER
//
  title = "LEHMER";
  n = 5;
  a = lehmer ( n, n );
  determ1 = lehmer_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LESLIE
//
  title = "LESLIE";
  n = 4;
  b = 0.025;
  di = 0.010;
  da = 0.100;
  a = leslie ( b, di, da );
  determ1 = leslie_determinant ( b, di, da );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LESP
//
  title = "LESP";
  n = 5;
  a = lesp ( n, n );
  determ1 = lesp_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LIETZKE
//
  title = "LIETZKE";
  n = 5;
  a = lietzke ( n );
  determ1 = lietzke_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LIGHTS_OUT
//
  title = "LIGHTS_OUT";
  if ( false )
  {
    row_num = 5;
    col_num = 5;
    n = row_num * col_num;
//  a = lights_out ( row_num, col_num, n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << left << title << right
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
//  LINE_ADJ
//
  title = "LINE_ADJ";
  n = 5;
  a = line_adj ( n );
  determ1 = line_adj_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LINE_ADJ
//
  title = "LINE_ADJ";
  n = 6;
  a = line_adj ( n );
  determ1 = line_adj_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LINE_LOOP_ADJ
//
  title = "LINE_LOOP_ADJ";
  n = 5;
  a = line_loop_adj ( n );
  determ1 = line_loop_adj_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  LOEWNER
//
  title = "LOEWNER";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  w = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  z = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = loewner ( w, x, y, z, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
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
//  LOTKIN
//
  title = "LOTKIN";
  n = 5;
  a = lotkin ( n, n );
  determ1 = lotkin_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MARKOV_RANDOM
//
  title = "MARKOV_RANDOM";
  n = 5;
  key = 123456789;
  a = markov_random ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MAXIJ
//
  title = "MAXIJ";
  n = 5;
  a = maxij ( n, n );
  determ1 = maxij_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MILNES
//
  title = "MILNES";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = milnes ( n, n, x );
  determ1 = milnes_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  MINIJ
//
  title = "MINIJ";
  n = 5;
  a = minij ( n, n );
  determ1 = minij_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MOLER1
//
  title = "MOLER1";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = moler1 ( alpha, n, n );
  determ1 = moler1_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MOLER2
//
  title = "MOLER2";
  n = 5;
  a = moler2 ( );
  determ1 = moler2_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MOLER3
//
  title = "MOLER3";
  n = 5;
  a = moler3 ( n, n );
  determ1 = moler3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  MOLER4
//
  title = "MOLER4";
  n = 4;
  a = moler4 ( );
  determ1 = moler4_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  NEUMANN
//
  title = "NEUMANN";
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  a = neumann ( row_num, col_num );
  determ1 = neumann_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ONE
//
  title = "ONE";
  n = 5;
  a = one ( n, n );
  determ1 = one_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ORTEGA
//
  title = "ORTEGA";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  v1 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  v2 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  v3 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = ortega ( n, v1, v2, v3 );
  determ1 = ortega_determinant ( n, v1, v2, v3 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] v1;
  delete [] v2;
  delete [] v3;
//
//  ORTH_RANDOM
//
  title = "ORTH_RANDOM";
  n = 5;
  key = 123456789;
  a = orth_random ( n, key );
  determ1 = orth_random_determinant ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ORTH_SYMM
//
  title = "ORTH_SYMM";
  n = 5;
  a = orth_symm ( n );
  determ1 = orth_symm_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  OTO
//
  title = "OTO";
  n = 5;
  a = oto ( n, n );
  determ1 = oto_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PARTER
//
  title = "PARTER";
  n = 5;
  a = parter ( n, n );
  determ1 = parter_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PASCAL1
//
  title = "PASCAL1";
  n = 5;
  a = pascal1 ( n );
  determ1 = pascal1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PASCAL2
//
  title = "PASCAL2";
  n = 5;
  a = pascal2 ( n );
  determ1 = pascal2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PASCAL3
//
  title = "PASCAL3";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = pascal3 ( n, alpha );
  determ1 = pascal3_determinant ( n, alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PDS_RANDOM
//
  title = "PDS_RANDOM";
  n = 5;
  key = 123456789;
  a = pds_random ( n, key );
  determ1 = pds_random_determinant ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PEI
//
  title = "PEI";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = pei ( alpha, n );
  determ1 = pei_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PERMUTATION_RANDOM
//
  title = "PERMUTATION_RANDOM";
  n = 5;
  key = 123456789;
  a = permutation_random ( n, key );
  determ1 = permutation_random_determinant ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PLU
//
  title = "PLU";
  n = 5;
  pivot = new int[n];
  seed = 123456789;
  for ( i = 0; i < n; i++ )
  {
    i4_lo = i;
    i4_hi = n - 1;
    pivot[i] = i4_uniform_ab ( i4_lo, i4_hi, seed );
  }
  a = plu ( n, pivot );
  determ1 = plu_determinant ( n, pivot );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] pivot;
//
//  POISSON
//
  title = "POISSON";
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  a = poisson ( row_num, col_num );
  determ1 = poisson_determinant ( row_num, col_num );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  PROLATE
//
  title = "PROLATE";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = prolate ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RECTANGLE_ADJ
//
  if ( false )
  {
    title = "RECTANGLE_ADJ";
    row_num = 5;
    col_num = 5;
    n = row_num * col_num;
//  a = rectangle_adj ( row_num, col_num, n );
    determ1 = rectangle_adj_determinant ( row_num, col_num );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << left << title << right
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
//  REDHEFFER
//
  title = "REDHEFFER";
  n = 5;
  a = redheffer ( n );
  determ1 = redheffer_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  REF_RANDOM
//
  title = "REF_RANDOM";
  n = 5;
  prob = 0.65;
  key = 123456789;
  a = ref_random ( n, n, prob, key );
  determ1 = ref_random_determinant ( n, prob, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  REF_RANDOM
//
  title = "REF_RANDOM";
  n = 5;
  prob = 0.85;
  key = 123456789;
  a = ref_random ( n, n, prob, key );
  determ1 = ref_random_determinant ( n, prob, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RIEMANN
//
  title = "RIEMANN";
  n = 5;
  a = riemann ( n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RING_ADJ
//
  title = "RING_ADJ";
  for ( n = 1; n <= 8; n++ )
  {
    a = ring_adj ( n );
    determ1 = ring_adj_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    cout << "  " << setw(20) << left << title << right
         << "  " << setw(4) << n
         << "  " << setw(14) << determ1
         << "  " << setw(14) << determ2
         << "  " << setw(14) << norm_frobenius << "\n";
    delete [] a;
  }
//
//  RIS
//
  title = "RIS";
  n = 5;
  a = ris ( n );
  determ1 = ris_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RODMAN
//
  title = "RODMAN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = rodman ( n, n, alpha );
  determ1 = rodman_determinant ( n, alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ROSSER1
//
  title = "ROSSER1";
  n = 8;
  a = rosser1 ( );
  determ1 = rosser1_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ROUTH
//
  title = "ROUTH";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = routh ( n, x );
  determ1 = routh_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  RUTIS1
//
  title = "RUTIS1";
  n = 4;
  a = rutis1 ( );
  determ1 = rutis1_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RUTIS2
//
  title = "RUTIS2";
  n = 4;
  a = rutis2 ( );
  determ1 = rutis2_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RUTIS3
//
  title = "RUTIS3";
  n = 4;
  a = rutis3 ( );
  determ1 = rutis3_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RUTIS4
//
  title = "RUTIS4";
  n = 4;
  a = rutis4 ( n );
  determ1 = rutis4_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  RUTIS5
//
  title = "RUTIS5";
  n = 4;
  a = rutis5 ( );
  determ1 = rutis5_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SCHUR_BLOCK
//
  title = "SCHUR_BLOCK";
  n = 5;
  x_n = ( n + 1 ) / 2;
  y_n = n / 2;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( x_n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( y_n, r8_lo, r8_hi, seed );
  a = schur_block ( n, x, y );
  determ1 = schur_block_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  SKEW_CIRCULANT
//
  title = "SKEW_CIRCULANT";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = skew_circulant ( n, n, x );
  determ1 = skew_circulant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  SPLINE
//
  title = "SPLINE";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = spline ( n, x );
  determ1 = spline_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  STIRLING
//
  title = "STIRLING";
  n = 5;
  a = stirling ( n, n );
  determ1 = stirling_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  STRIPE
//
  title = "STRIPE";
  n = 5;
  a = stripe ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SUMMATION
//
  title = "SUMMATION";
  n = 5;
  a = summation ( n, n );
  determ1 = summation_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SWEET1
//
  title = "SWEET1";
  n = 6;
  a = sweet1 ( );
  determ1 = sweet1_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SWEET2
//
  title = "SWEET2";
  n = 6;
  a = sweet2 ( );
  determ1 = sweet2_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SWEET3
//
  title = "SWEET3";
  n = 6;
  a = sweet3 ( );
  determ1 = sweet3_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SWEET4
//
  title = "SWEET4";
  n = 13;
  a = sweet4 ( );
  determ1 = sweet4_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SYLVESTER
//
  title = "SYLVESTER";
  n = 5;
  x_n = 3 + 1;
  y_n = 2 + 1;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( x_n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( y_n, r8_lo, r8_hi, seed );
  a = sylvester ( n, x_n - 1, x, y_n - 1, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  SYLVESTER_KAC
//
  title = "SYLVESTER_KAC";
  n = 5;
  a = sylvester_kac ( n );
  determ1 = sylvester_kac_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SYLVESTER_KAC
//
  title = "SYLVESTER_KAC";
  n = 6;
  a = sylvester_kac ( n );
  determ1 = sylvester_kac_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  SYMM_RANDOM
//
  title = "SYMM_RANDOM";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  d = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  key = 123456789;
  a = symm_random ( n, d, key );
  determ1 = symm_random_determinant ( n, d, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] d;
//
//  TOEPLITZ
//
  title = "TOEPLITZ";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( 2 * n - 1, r8_lo, r8_hi, seed );
  a = toeplitz ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  TOEPLITZ_5DIAG
//
  title = "TOEPLITZ_5DIAG";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  d1 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  d2 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  d3 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  d4 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  d5 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = toeplitz_5diag ( n, d1, d2, d3, d4, d5 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TOEPLITZ_5S
//
  title = "TOEPLITZ_5S";
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  r8_lo = -5.0;
  r8_hi = +5.0; 
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  gamma = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = toeplitz_5s ( row_num, col_num, alpha, beta, gamma, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TOEPLITZ_PDS
//
  title = "TOEPLITZ_PDS";
  m = 3;
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( m, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( m, r8_lo, r8_hi, seed );
  y_sum = r8vec_sum ( m, y );
  for ( i = 0; i < m; i++ )
  {
    y[i] = y[i] / y_sum;
  }
  a = toeplitz_pds ( m, n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
//
//  TOURNAMENT_RANDOM
//
  title = "TOURNAMENT_RANDOM";
  n = 5;
  key = 123456789;
  a = tournament_random ( n, key );
  determ1 = tournament_random_determinant ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRANSITION_RANDOM
//
  title = "TRANSITION_RANDOM";
  n = 5;
  key = 123456789;
  a = transition_random ( n, key );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRENCH
//
  title = "TRENCH";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = trench ( alpha, n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRI_UPPER
//
  title = "TRI_UPPER";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = tri_upper ( alpha, n );
  determ1 = tri_upper_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRIS
//
  title = "TRIS";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  gamma = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = tris ( n, n, alpha, beta, gamma );
  determ1 = tris_determinant ( n, alpha, beta, gamma );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  TRIV
//
  title = "TRIV";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  z = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  a = triv ( n, x, y, z );
  determ1 = triv_determinant ( n, x, y, z );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
  delete [] y;
  delete [] z;
//
//  TRIW
//
  title = "TRIW";
  n = 5;
  i4_lo = 0;
  i4_hi = n - 1;
  seed = 123456789;
  k = i4_uniform_ab ( i4_lo, i4_hi, seed );
  r8_lo = -5.0;
  r8_hi = +5.0;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = triw ( alpha, k, n );
  determ1 = triw_determinant ( alpha, k, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  UPSHIFT
//
  title = "UPSHIFT";
  n = 5;
  a = upshift ( n );
  determ1 = upshift_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  VAND1
//
  title = "VAND1";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = vand1 ( n, x );
  determ1 = vand1_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  VAND2
//
  title = "VAND2";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = vand2 ( n, x );
  determ1 = vand2_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
  delete [] x;
//
//  WATHEN
//
  title = "WATHEN";
  if ( false )
  {
  row_num = 5;
  col_num = 5;
  n = wathen_order ( row_num, col_num );
  a = wathen ( row_num, col_num, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
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
//  WILK03
//
  title = "WILK03";
  n = 3;
  a = wilk03 ( );
  determ1 = wilk03_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK04
//
  title = "WILK04";
  n = 4;
  a = wilk04 ( );
  determ1 = wilk04_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK05
//
  title = "WILK05";
  n = 5;
  a = wilk05 ( );
  determ1 = wilk05_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK12
//
  title = "WILK12";
  n = 12;
  a = wilk12 ( );
  determ1 = wilk12_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK20
//
  title = "WILK20";
  n = 20;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = wilk20 ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILK21
//
  title = "WILK21";
  n = 21;
  a = wilk21 ( n );
  determ1 = wilk21_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  WILSON
//
  title = "WILSON";
  n = 4;
  a = wilson ( );
  determ1 = wilson_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ZERO
//
  title = "ZERO";
  n = 5;
  a = zero ( n, n );
  determ1 = zero_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(14) << determ1
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;
//
//  ZIELKE
//
  title = "ZIELKE";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  d1 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  d2 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  d3 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = zielke ( n, d1, d2, d3 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << "              "
       << "  " << setw(14) << determ2
       << "  " << setw(14) << norm_frobenius << "\n";
  delete [] a;

  return;
}
//****************************************************************************80

void test_eigen_left ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_EIGEN_LEFT tests left eigensystems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double beta;
  double *d;
  double error_frobenius;
  double gamma;
  int i;
  int i1;
  int i4_hi;
  int i4_lo;
  int k;
  int key;
  double *lambda;
  int n;
  double norm_frobenius;
  double r8_hi;
  double r8_lo;
  int rank;
  int seed;
  int seed_save;
  string title;
  double *v1;
  double *v2;
  double *v3;
  double *x;

  cout << "\n";
  cout << "TEST_EIGEN_LEFT\n";
  cout << "  Compute the Frobenius norm of the eigenvalue error:\n";
  cout << "    X * A - LAMBDA * X\n";
  cout << "  given K left eigenvectors X and eigenvalues LAMBDA.\n";
  cout << "\n";
  cout << "  Title                    N     K          ||A||       ||X*A-Lambda*X||\n";
  cout << "\n";
//
//  A123
//
  title = "A123";
  n = 3;
  k = 3;
  a = a123 ( );
  lambda = a123_eigenvalues ( );
  x = a123_eigen_left ( );
  error_frobenius = r8mat_is_eigen_left ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  CARRY
//
  title = "CARRY";
  n = 5;
  k = 5;
  i4_lo = 2;
  i4_hi = 20;
  seed = 123456789;
  i1 = i4_uniform_ab ( i4_lo, i4_hi, seed );
  a = carry ( n, i1 );
  lambda = carry_eigenvalues ( n, i1 );
  x = carry_eigen_left ( n, i1 );
  error_frobenius = r8mat_is_eigen_left ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  CHOW
//
  title = "CHOW";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = chow ( alpha, beta, n, n );
  lambda = chow_eigenvalues ( alpha, beta, n );
  x = chow_eigen_left ( alpha, beta, n );
  error_frobenius = r8mat_is_eigen_left ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  DIAGONAL
//
  title = "DIAGONAL";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  d = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = diagonal ( n, n, d );
  lambda = diagonal_eigenvalues ( n, d );
  x = diagonal_eigen_left ( n, d );
  error_frobenius = r8mat_is_eigen_left ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] d;
  delete [] lambda;
  delete [] x;
//
//  ROSSER1
//
  title = "ROSSER1";
  n = 8;
  k = 8;
  a = rosser1 ( );
  lambda = rosser1_eigenvalues ( );
  x = rosser1_eigen_left ( );
  error_frobenius = r8mat_is_eigen_left ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  SYMM_RANDOM
//
  title = "SYMM_RANDOM";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  d = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  key = 123456789;
  a = symm_random ( n, d, key );
  lambda = symm_random_eigenvalues ( n, d, key );
  x = symm_random_eigen_left ( n, d, key );
  error_frobenius = r8mat_is_eigen_left ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] d;
  delete [] lambda;
  delete [] x;

  return;
}
//****************************************************************************80

void test_eigen_right ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_EIGEN_RIGHT tests right eigensystems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double beta;
  double *d;
  double error_frobenius;
  double gamma;
  int i;
  int i1;
  int i4_hi;
  int i4_lo;
  int k;
  int key;
  double *lambda;
  int n;
  double norm_frobenius;
  double r8_hi;
  double r8_lo;
  int rank;
  int seed;
  int seed_save;
  string title;
  double *v1;
  double *v2;
  double *v3;
  double *x;

  cout << "\n";
  cout << "TEST_EIGEN_RIGHT\n";
  cout << "  Compute the Frobenius norm of the eigenvalue error:\n";
  cout << "    A * X - X * LAMBDA\n";
  cout << "  given K right eigenvectors X and eigenvalues LAMBDA.\n";
  cout << "\n";
  cout << "  Title                    N     K          ||A||       ||A*X-X*Lambda||\n";
  cout << "\n";
//
//  A123
//
  title = "A123";
  n = 3;
  k = 3;
  a = a123 ( );
  lambda = a123_eigenvalues ( );
  x = a123_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  BAB
//
  title = "BAB";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = bab ( n, alpha, beta );
  lambda = bab_eigenvalues ( n, alpha, beta );
  x = bab_eigen_right ( n, alpha, beta );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  BODEWIG
//
  title = "BODEWIG";
  n = 4;
  k = 4;
  a = bodewig ( );
  lambda = bodewig_eigenvalues ( );
  x = bodewig_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  CARRY
//
  title = "CARRY";
  n = 5;
  k = 5;
  i4_lo = 2;
  i4_hi = 20;
  seed = 123456789;
  i1 = i4_uniform_ab ( i4_lo, i4_hi, seed );
  a = carry ( n, i1 );
  lambda = carry_eigenvalues ( n, i1 );
  x = carry_eigen_right ( n, i1 );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  CHOW
//
  title = "CHOW";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = chow ( alpha, beta, n, n );
  lambda = chow_eigenvalues ( alpha, beta, n );
  x = chow_eigen_right ( alpha, beta, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  COMBIN
//
  title = "COMBIN";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = combin ( alpha, beta, n );
  lambda = combin_eigenvalues ( alpha, beta, n );
  x = combin_eigen_right ( alpha, beta, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  DIF2
//
  title = "DIF2";
  n = 5;
  k = 5;
  a = dif2 ( n, n );
  lambda = dif2_eigenvalues ( n );
  x = dif2_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  EXCHANGE
//
  title = "EXCHANGE";
  n = 5;
  k = 5;
  a = exchange ( n, n );
  lambda = exchange_eigenvalues ( n );
  x = exchange_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  IDEM_RANDOM
//
  title = "IDEM_RANDOM";
  n = 5;
  k = 5;
  rank = 3;
  key = 987654321;
  a = idem_random ( n, rank, key );
  lambda = idem_random_eigenvalues ( n, rank, key );
  x = idem_random_eigen_right ( n, rank, key );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  IDENTITY
//
  title = "IDENTITY";
  n = 5;
  k = 5;
  a = identity ( n, n );
  lambda = identity_eigenvalues ( n );
  x = identity_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ILL3
//
  title = "ILL3";
  n = 3;
  k = 3;
  a = ill3 ( );
  lambda = ill3_eigenvalues ( );
  x = ill3_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  KERSHAW
//
  title = "KERSHAW";
  n = 4;
  k = 4;
  a = kershaw ( );
  lambda = kershaw_eigenvalues ( );
  x = kershaw_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  KMS
//  Eigenvalue information requires 0 <= ALPHA <= 1.
//
  title = "KMS";
  n = 5;
  k = 5;
  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = kms ( alpha, n, n );
  lambda = kms_eigenvalues ( alpha, n );
  x = kms_eigen_right ( alpha, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  LINE_ADJ
//
  title = "LINE_ADJ";
  n = 5;
  k = 5;
  a = line_adj ( n );
  lambda = line_adj_eigenvalues ( n );
  x = line_adj_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  LINE_LOOP_ADJ
//
  title = "LINE_LOOP_ADJ";
  n = 5;
  k = 5;
  a = line_loop_adj ( n );
  lambda = line_loop_adj_eigenvalues ( n );
  x = line_loop_adj_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ONE
//
  title = "ONE";
  n = 5;
  k = 5;
  a = one ( n, n );
  lambda = one_eigenvalues ( n );
  x = one_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ORTEGA
//
  title = "ORTEGA";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  v1 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  v2 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  v3 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = ortega ( n, v1, v2, v3 );
  lambda = ortega_eigenvalues ( n, v1, v2, v3 );
  x = ortega_eigen_right ( n, v1, v2, v3 );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
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
//  OTO
//
  title = "OTO";
  n = 5;
  k = 5;
  a = oto ( n, n );
  lambda = oto_eigenvalues ( n );
  x = oto_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  PDS_RANDOM
//
  title = "PDS_RANDOM";
  n = 5;
  k = 5;
  key = 123456789;
  a = pds_random ( n, key );
  lambda = pds_random_eigenvalues ( n, key );
  x = pds_random_eigen_right ( n, key );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  PEI
//
  title = "PEI";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = pei ( alpha, n );
  lambda = pei_eigenvalues ( alpha, n );
  x = pei_eigen_right ( alpha, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  RODMAN
//
  title = "RODMAN";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = rodman ( n, n, alpha );
  lambda = rodman_eigenvalues ( n, alpha );
  x = rodman_eigen_right ( n, alpha );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ROSSER1
//
  title = "ROSSER1";
  n = 8;
  k = 8;
  a = rosser1 ( );
  lambda = rosser1_eigenvalues ( );
  x = rosser1_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  RUTIS1
//
  title = "RUTIS1";
  n = 4;
  k = 4;
  a = rutis1 ( );
  lambda = rutis1_eigenvalues ( );
  x = rutis1_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  RUTIS2
//
  title = "RUTIS2";
  n = 4;
  k = 4;
  a = rutis2 ( );
  lambda = rutis2_eigenvalues ( );
  x = rutis2_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  RUTIS5
//
  title = "RUTIS5";
  n = 4;
  k = 4;
  a = rutis5 ( );
  lambda = rutis5_eigenvalues ( );
  x = rutis5_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  SYLVESTER_KAC
//
  title = "SYLVESTER_KAC";
  n = 5;
  k = 5;
  a = sylvester_kac ( n );
  lambda = sylvester_kac_eigenvalues ( n );
  x = sylvester_kac_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  SYMM_RANDOM
//
  title = "SYMM_RANDOM";
  n = 5;
  k = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  d = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  key = 123456789;
  a = symm_random ( n, d, key );
  lambda = symm_random_eigenvalues ( n, d, key );
  x = symm_random_eigen_right ( n, d, key );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] d;
  delete [] lambda;
  delete [] x;
//
//  WILK12
//
  title = "WILK12";
  n = 12;
  k = 12;
  a = wilk12 ( );
  lambda = wilk12_eigenvalues ( );
  x = wilk12_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  WILSON
//
  title = "WILSON";
  n = 4;
  k = 4;
  a = wilson ( );
  lambda = wilson_eigenvalues ( );
  x = wilson_eigen_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] lambda;
  delete [] x;
//
//  ZERO
//
  title = "ZERO";
  n = 5;
  k = 5;
  a = zero ( n, n );
  lambda = zero_eigenvalues ( n );
  x = zero_eigen_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
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
//    17 April 2015
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
  double *d;
  double error_ab;
  double error_ac;;
  double gamma;
  int i;
  int i4_hi;
  int i4_lo;
  int ii;
  int jj;
  int k;
  int key;
  double *l;
  int n;
  double norma_frobenius;
  double normc_frobenius;
  double *p;
  int *pivot;
  double r8_hi;
  double r8_lo;
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
  cout << "  Title                     N    "
       << "   ||A||      ||C||  ||I-AC||    ||I-AB||\n";
  cout << "\n";
//
//  AEGERTER
//
  title = "AEGERTER";
  n = 5;
  a = aegerter ( n );
  b = aegerter_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BAB
//
  title = "BAB";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = bab ( n, alpha, beta );
  b = bab_inverse ( n, alpha, beta );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BAUER
//
  title = "BAUER";
  n = 6;
  a = bauer ( );
  b = bauer_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BERNSTEIN
//
  title = "BERNSTEIN";
  n = 5;
  a = bernstein ( n );
  b = bernstein_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BIS
//
  title = "BIS";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = bis ( alpha, beta, n, n );
  b = bis_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";;
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BIW
//
  title = "BIW";
  n = 5;
  a = biw ( n );
  b = biw_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BODEWIG
//
  title = "BODEWIG";
  n = 4;
  a = bodewig ( );
  b = bodewig_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BOOTHROYD
//
  title = "BOOTHROYD";
  n = 5;
  a = boothroyd ( n );
  b = boothroyd_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  BORDERBAND
//
  title = "BORDERBAND";
  n = 5;
  a = borderband ( n );
  b = borderband_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CARRY
//
  title = "CARRY";
  n = 5;
  i4_lo = 2;
  i4_hi = 20;
  seed = 123456789;
  k = i4_uniform_ab ( i4_lo, i4_hi, seed );
  a = carry ( n, k );
  b = carry_inverse ( n, k );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CAUCHY
//
  title = "CAUCHY";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = cauchy ( n, x, y );
  b = cauchy_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] y;
//
//  CHEBY_T
//
  title = "CHEBY_T";
  n = 5;
  a = cheby_t ( n );
  b = cheby_t_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CHEBY_U
//
  title = "CHEBY_U";
  n = 5;
  a = cheby_u ( n );
  b = cheby_u_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CHEBY_VAN2
//
  title = "CHEBY_VAN2";
  n = 5;
  a = cheby_van2 ( n );
  b = cheby_van2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CHEBY_VAN3
//
  title = "CHEBY_VAN3";
  n = 5;
  a = cheby_van3 ( n );
  b = cheby_van3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CHOW
//
  title = "CHOW";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = chow ( alpha, beta, n, n );
  b = chow_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CIRCULANT
//
  title = "CIRCULANT";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = circulant ( n, n, x );
  b = circulant_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  CIRCULANT2
//
  title = "CIRCULANT2";
  n = 5;
  a = circulant2 ( n );
  b = circulant2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CLEMENT1
//  N must be even.
//
  title = "CLEMENT1";
  n = 6;
  a = clement1 ( n );
  b = clement1_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CLEMENT2
//  N must be even.
//
  title = "CLEMENT2";
  n = 6;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  a = clement2 ( n, x, y );
  b = clement2_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] y;
//
//  COMBIN
//
  title = "COMBIN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = combin ( alpha, beta, n );
  b = combin_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  COMPANION
//
  title = "COMPANION";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = companion ( n, x );
  b = companion_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  COMPLEX_I
//
  title = "COMPLEX_I";
  n = 2;
  a = complex_i ( );
  b = complex_i_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CONEX1
//
  title = "CONEX1";
  n = 4;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = conex1 ( alpha );
  b = conex1_inverse ( alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CONEX2
//
  title = "CONEX2";
  n = 3;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = conex2 ( alpha );
  b = conex2_inverse ( alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CONEX3
//
  title = "CONEX3";
  n = 5;
  a = conex3 ( n );
  b = conex3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  CONFERENCE
//  N-1 must be an odd prime or a power of an odd prime.
//
  title = "CONFERENCE";
  n = 6;
  a = conference ( n );
  b = conference_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB2
//
  title = "DAUB2";
  n = 4;
  a = daub2 ( n );
  b = daub2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB4
//
  title = "DAUB4";
  n = 8;
  a = daub4 ( n );
  b = daub4_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB6
//
  title = "DAUB6";
  n = 12;
  a = daub6 ( n );
  b = daub6_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB8
//
  title = "DAUB8";
  n = 16;
  a = daub8 ( n );
  b = daub8_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB10
//
  title = "DAUB10";
  n = 20;
  a = daub10 ( n );
  b = daub10_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DAUB12
//
  title = "DAUB12";
  n = 24;
  a = daub12 ( n );
  b = daub12_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DIAGONAL
//
  title = "DIAGONAL";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = diagonal ( n, n, x );
  b = diagonal_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  DIF1
//  N must be even.
//
  title = "DIF1";
  n = 6;
  a = dif1 ( n, n );
  b = dif1_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DIF2
//
  title = "DIF2";
  n = 5;
  a = dif2 ( n, n );
  b = dif2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DORR
//
  title = "DORR";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = dorr ( alpha, n );
  b = dorr_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  DOWNSHIFT
//
  title = "DOWNSHIFT";
  n = 5;
  a = downshift ( n );
  b = downshift_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  EULERIAN
//
  title = "EULERIAN";
  n = 5;
  a = eulerian ( n, n );
  b = eulerian_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  EXCHANGE
//
  title = "EXCHANGE";
  n = 5;
  a = exchange ( n, n );
  b = exchange_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FIBONACCI2
//
  title = "FIBONACCI2";
  n = 5;
  a = fibonacci2 ( n );
  b = fibonacci2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FIBONACCI3
//
  title = "FIBONACCI3";
  n = 5;
  a = fibonacci3 ( n );
  b = fibonacci3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FIEDLER.
//  The FIEDLER_INVERSE routine assumes the X vector is sorted.
//
  title = "FIEDLER";
  n = 7;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  r8vec_sort_bubble_a ( n, x );
  a = fiedler ( n, n, x );
  b = fiedler_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  FORSYTHE
//
  title = "FORSYTHE";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = forsythe ( alpha, beta, n );
  b = forsythe_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FOURIER_COSINE
//
  title = "FOURIER_COSINE";
  n = 5;
  a = fourier_cosine ( n );
  b = fourier_cosine_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FOURIER_SINE
//
  title = "FOURIER_SINE";
  n = 5;
  a = fourier_sine ( n );
  b = fourier_sine_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  FRANK
//
  title = "FRANK";
  n = 5;
  a = frank ( n );
  b = frank_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GFPP
//
  title = "GFPP";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = gfpp ( n, alpha );
  b = gfpp_inverse ( n, alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GIVENS
//
  title = "GIVENS";
  n = 5;
  a = givens ( n, n );
  b = givens_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GK316
//
  title = "GK316";
  n = 5;
  a = gk316 ( n );
  b = gk316_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GK323
//
  title = "GK323";
  n = 5;
  a = gk323 ( n, n );
  b = gk323_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  GK324
//
  title = "GK324";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  a = gk324 ( n, n, x );
  b = gk324_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  HANKEL_N
//
  title = "HANKEL_N";
  n = 5;
  a = hankel_n ( n );
  b = hankel_n_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HANOWA
//  N must be even.
//
  title = "HANOWA";
  n = 6;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = hanowa ( alpha, n );
  b = hanowa_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HARMAN
//
  title = "HARMAN";
  n = 8;
  a = harman (  );
  b = harman_inverse (  );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HARTLEY
//
  title = "HARTLEY";
  n = 5;
  a = hartley ( n );
  b = hartley_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HELMERT
//
  title = "HELMERT";
  n = 5;
  a = helmert ( n );
  b = helmert_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HELMERT2
//
  title = "HELMERT2";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = helmert2 ( n, x );
  b = helmert2_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  HERMITE
//
  title = "HERMITE";
  n = 5;
  a = hermite ( n );
  b = hermite_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HERNDON
//
  title = "HERNDON";
  n = 5;
  a = herndon ( n );
  b = herndon_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HILBERT
//
  title = "HILBERT";
  n = 5;
  a = hilbert ( n, n );
  b = hilbert_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  HOUSEHOLDER
//
  title = "HOUSEHOLDER";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = householder ( n, x );
  b = householder_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  IDENTITY
//
  title = "IDENTITY";
  n = 5;
  a = identity ( n, n );
  b = identity_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  ILL3
//
  title = "ILL3";
  n = 3;
  a = ill3 ( );
  b = ill3_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  INTEGRATION
//
  title = "INTEGRATION";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = integration ( alpha, n );
  b = integration_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  INVOL
//
  title = "INVOL";
  n = 5;
  a = invol ( n );
  b = invol_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  JACOBI
//  N must be even.
//
  title = "JACOBI";
  n = 6;
  a = jacobi ( n, n );
  b = jacobi_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  JORDAN
//
  title = "JORDAN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = jordan ( n, n, alpha );
  b = jordan_inverse ( n, alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  KAHAN
//
  title = "KAHAN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = kahan ( alpha, n, n );
  b = kahan_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  KERSHAW
//
  title = "KERSHAW";
  n = 4;
  a = kershaw ( );
  b = kershaw_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  KERSHAWTRI
//
  title = "KERSHAWTRI";
  n = 5;
  x_n = ( n + 1 ) / 2;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( x_n, r8_lo, r8_hi, seed );
  a = kershawtri ( n, x );
  b = kershawtri_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  KMS
//
  title = "KMS";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = kms ( alpha, n, n );
  b = kms_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LAGUERRE
//
  title = "LAGUERRE";
  n = 5;
  a = laguerre ( n );
  b = laguerre_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LEGENDRE
//
  title = "LEGENDRE";
  n = 5;
  a = legendre ( n );
  b = legendre_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LEHMER
//
  title = "LEHMER";
  n = 5;
  a = lehmer ( n, n );
  b = lehmer_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LESP
//
  title = "LESP";
  n = 5;
  a = lesp ( n, n );
  b = lesp_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LIETZKE
//
  title = "LIETZKE";
  n = 5;
  a = lietzke ( n );
  b = lietzke_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LINE_ADJ
//  N must be even.
//
  title = "LINE_ADJ";
  n = 6;
  a = line_adj ( n );
  b = line_adj_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  LOTKIN
//
  title = "LOTKIN";
  n = 5;
  a = lotkin ( n, n );
  b = lotkin_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  MAXIJ
//
  title = "MAXIJ";
  n = 5;
  a = maxij ( n, n );
  b = maxij_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  MILNES
//
  title = "MILNES";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = milnes ( n, n, x );
  b = milnes_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  MINIJ
//
  title = "MINIJ";
  n = 5;
  a = minij ( n, n );
  b = minij_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  MOLER1
//
  title = "MOLER1";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = moler1 ( alpha, n, n );
  b = moler1_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  MOLER3
//
  title = "MOLER3";
  n = 5;
  a = moler3 ( n, n );
  b = moler3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  ORTEGA
//
  title = "ORTEGA";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  v1 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  v2 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  v3 = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = ortega ( n, v1, v2, v3 );
  b = ortega_inverse ( n, v1, v2, v3 );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] v1;
  delete [] v2;
  delete [] v3;
//
//  ORTH_SYMM
//
  title = "ORTH_SYMM";
  n = 5;
  a = orth_symm ( n );
  b = orth_symm_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  OTO
//
  title = "OTO";
  n = 5;
  a = oto ( n, n );
  b = oto_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PARTER
//
  title = "PARTER";
  n = 5;
  a = parter ( n, n );
  b = parter_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PASCAL1
//
  title = "PASCAL1";
  n = 5;
  a = pascal1 ( n );
  b = pascal1_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PASCAL2
//
  title = "PASCAL2";
  n = 5;
  a = pascal2 ( n );
  b = pascal2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PASCAL3
//
  title = "PASCAL3";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = pascal3 ( n, alpha );
  b = pascal3_inverse ( n, alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PDS_RANDOM
//
  title = "PDS_RANDOM";
  n = 5;
  key = 123456789;
  a = pds_random ( n, key );
  b = pds_random_inverse ( n, key );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PEI
//
  title = "PEI";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = pei ( alpha, n );
  b = pei_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PERMUTATION_RANDOM
//
  title = "PERMUTATION_RANDOM";
  n = 5;
  key = 123456789;
  a = permutation_random ( n, key );
  b = permutation_random_inverse ( n, key );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  PLU
//
  title = "PLU";
  n = 5;
  pivot = new int[n];
  seed = 123456789;
  for ( i = 0; i < n; i++ )
  {
    i4_lo = i;
    i4_hi = n - 1;
    pivot[i] = i4_uniform_ab ( i4_lo, i4_hi, seed );
  }
  a = plu ( n, pivot );
  b = plu_inverse ( n, pivot );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] pivot;
//
//  RIS
//
  title = "RIS";
  n = 5;
  a = ris ( n );
  b = ris_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RODMAN
//
  title = "RODMAN";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = rodman ( n, n, alpha );
  b = rodman_inverse ( n, alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS1
//
  title = "RUTIS1";
  n = 4;
  a = rutis1 ( );
  b = rutis1_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS2
//
  title = "RUTIS2";
  n = 4;
  a = rutis2 ( );
  b = rutis2_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS3
//
  title = "RUTIS3";
  n = 4;
  a = rutis3 ( );
  b = rutis3_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS4
//
  title = "RUTIS4";
  n = 5;
  a = rutis4 ( n );
  b = rutis4_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  RUTIS5
//
  title = "RUTIS5";
  n = 4;
  a = rutis5 ( );
  b = rutis5_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SCHUR_BLOCK
//
  title = "SCHUR_BLOCK";
  n = 5;
  x_n = ( n + 1 ) / 2;
  y_n = n / 2;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( x_n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( y_n, r8_lo, r8_hi, seed );
  a = schur_block ( n, x, y );
  b = schur_block_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] y;
//
//  SPLINE
//
  title = "SPLINE";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  a = spline ( n, x );
  b = spline_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  STIRLING
//
  title = "STIRLING";
  n = 5;
  a = stirling ( n, n );
  b = stirling_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SUMMATION
//
  title = "SUMMATION";
  n = 5;
  a = summation ( n, n );
  b = summation_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SWEET1
//
  title = "SWEET1";
  n = 6;
  a = sweet1 ( );
  b = sweet1_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SWEET2
//
  title = "SWEET2";
  n = 6;
  a = sweet2 ( );
  b = sweet2_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SWEET3
//
  title = "SWEET3";
  n = 6;
  a = sweet3 ( );
  b = sweet3_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SWEET4
//
  title = "SWEET4";
  n = 13;
  a = sweet4 ( );
  b = sweet4_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SYLVESTER_KAC
//  N must be even.
//
  title = "SYLVESTER_KAC";
  n = 6;
  a = sylvester_kac ( n );
  b = sylvester_kac_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  SYMM_RANDOM
//
  title = "SYMM_RANDOM";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  d = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  key = 123456789;
  a = symm_random ( n, d, key );
  b = symm_random_inverse ( n, d, key );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] d;
//
//  TRI_UPPER
//
  title = "TRI_UPPER";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = tri_upper ( alpha, n );
  b = tri_upper_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  TRIS
//
  title = "TRIS";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed );
  gamma = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = tris ( n, n, alpha, beta, gamma );
  b = tris_inverse ( n, alpha, beta, gamma );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  TRIV
//
  title = "TRIV";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  z = r8vec_uniform_ab_new ( n - 1, r8_lo, r8_hi, seed );
  a = triv ( n, x, y, z );
  b = triv_inverse ( n, x, y, z );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] y;
  delete [] z;
//
//  TRIW
//
  title = "TRIW";
  n = 5;
  i4_lo = 0;
  i4_hi = n - 1;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  k = i4_uniform_ab ( i4_lo, i4_hi, seed );
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = triw ( alpha, k, n );
  b = triw_inverse ( alpha, k, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  UPSHIFT
//
  title = "UPSHIFT";
  n = 5;
  a = upshift ( n );
  b = upshift_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  VAND1
//
  title = "VAND1";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = vand1 ( n, x );
  b = vand1_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
//
//  VAND2
//
  title = "VAND2";
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  a = vand2 ( n, x );
  b = vand2_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILK03
//
  title = "WILK03";
  n = 3;
  a = wilk03 ( );
  b = wilk03_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILK04
//
  title = "WILK04";
  n = 4;
  a = wilk04 ( );
  b = wilk04_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILK05
//
  title = "WILK05";
  n = 5;
  a = wilk05 ( );
  b = wilk05_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILK21
//
  title = "WILK21";
  n = 21;
  a = wilk21 ( n );
  b = wilk21_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;
//
//  WILSON
//
  title = "WILSON";
  n = 4;
  a = wilson ( );
  b = wilson_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << n
       << "  " << setw(10) << norma_frobenius
       << "  " << setw(10) << normc_frobenius
       << "  " << setw(10) << error_ac
       << "  " << setw(10) << error_ab << "\n";
  delete [] a;
  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void test_llt ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LLT tests LLT factors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2015
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
  double r8_hi;
  double r8_lo;
  int seed;
  string title;

  cout << "\n";
  cout << "TEST_LLT\n";
  cout << "  A = a test matrix of order M by M\n";
  cout << "  L is an M by N lower triangular Cholesky factor.\n";
  cout << "\n";
  cout << "  ||A|| = Frobenius norm of A.\n";
  cout << "  ||A-LLT|| = Frobenius norm of A-L*L'.\n";
  cout << "\n";
  cout << "  Title                    M     N      ||A||            ||A-LLT||\n";
  cout << "\n";
//
//  DIF2
//
  title = "DIF2";
  m = 5;
  n = 5;
  a = dif2 ( m, n );
  l = dif2_llt ( n );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  GIVENS
//
  title = "GIVENS";
  m = 5;
  n = 5;
  a = givens ( m, n );
  l = givens_llt ( n );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  KERSHAW
//
  title = "KERSHAW";
  m = 4;
  n = 4;
  a = kershaw ( );
  l = kershaw_llt ( );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  LEHMER
//
  title = "LEHMER";
  m = 5;
  n = 5;
  a = lehmer ( n, n );
  l = lehmer_llt ( n );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  MINIJ
//
  title = "MINIJ";
  m = 5;
  n = 5;
  a = minij ( n, n );
  l = minij_llt ( n );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  MOLER1
//
  title = "MOLER1";
  m = 5;
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = moler1 ( alpha, m, n );
  l = moler1_llt ( alpha, n );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  MOLER3
//
  title = "MOLER3";
  m = 5;
  n = 5;
  a = moler3 ( m, n );
  l = moler3_llt ( n );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  OTO
//
  title = "OTO";
  m = 5;
  n = 5;
  a = oto ( m, n );
  l = oto_llt ( n );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  PASCAL2
//
  title = "PASCAL2";
  m = 5;
  n = 5;
  a = pascal2 ( n );
  l = pascal2_llt ( n );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;
//
//  WILSON
//
  title = "WILSON";
  m = 4;
  n = 4;
  a = wilson ( );
  l = wilson_llt ( );
  error_frobenius = r8mat_is_llt ( m, n, a, l );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n",
  delete [] a;
  delete [] l;

  return;
}
//****************************************************************************80

void test_null_left ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_NULL_LEFT tests left null vectors.
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
  double alpha;
  double error_l2;
  double f1;
  double f2;
  int m;
  int n;
  double norm_a_frobenius;
  double norm_x_l2;
  double r8_hi;
  double r8_lo;
  int seed;
  string title;
  double *x;

  cout << "\n";
  cout << "TEST_NULL_LEFT\n";
  cout << "  A = a test matrix of order M by N\n";
  cout << "  x = an M vector, candidate for a left null vector.\n";
  cout << "\n";
  cout << "  ||A|| = Frobenius norm of A.\n";
  cout << "  ||x|| = L2 norm of x.\n";
  cout << "  ||A'*x||/||x|| = L2 norm of A'*x over L2 norm of x.\n";
  cout << "\n";
  cout << "  Title                    M     N      "
       << "||A||            ||x||        ||A'*x||/||x||\n";
  cout << "\n";
//
//  A123
//
  title = "A123";
  m = 3;
  n = 3;
  a = a123 ( );
  x = a123_null_left ( );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  CHEBY_DIFF1
//
  title = "CHEBY_DIFF1";
  m = 5;
  n = 5;
  a = cheby_diff1 ( n );
  x = cheby_diff1_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  CREATION
//
  title = "CREATION";
  m = 5;
  n = 5;
  a = creation ( m, n );
  x = creation_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF1
//  Only has null vectors for M odd.
//
  title = "DIF1";
  m = 5;
  n = 5;
  a = dif1 ( m, n );
  x = dif1_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF1CYCLIC
//
  title = "DIF1CYCLIC";
  m = 5;
  n = 5;
  a = dif1cyclic ( n );
  x = dif1cyclic_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF2CYCLIC
//
  title = "DIF2CYCLIC";
  m = 5;
  n = 5;
  a = dif2cyclic ( n );
  x = dif2cyclic_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  EBERLEIN
//
  title = "EBERLEIN";
  m = 5;
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = eberlein ( alpha, n );
  x = eberlein_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  FIBONACCI1
//
  title = "FIBONACCI1";
  m = 5;
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  f1 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  f2 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = fibonacci1 ( n, f1, f2 );
  x = fibonacci1_null_left ( m, n, f1, f2 );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  LAUCHLI
//
  title = "LAUCHLI";
  m = 6;
  n = m - 1;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = lauchli ( alpha, m, n );
  x = lauchli_null_left ( alpha, m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  LINE_ADJ
//
  title = "LINE_ADJ";
  m = 7;
  n = 7;
  a = line_adj ( n );
  x = line_adj_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  MOLER2
//
  title = "MOLER2";
  m = 5;
  n = 5;
  a = moler2 ( );
  x = moler2_null_left ( );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ONE
//
  title = "ONE";
  m = 5;
  n = 5;
  a = one ( n, n );
  x = one_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  RING_ADJ
//  M must be a multiple of 4 for there to be a null vector.
//
  title = "RING_ADJ";
  m = 12;
  n = 12;
  a = ring_adj ( n );
  x = ring_adj_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ROSSER1
//
  title = "ROSSER1";
  m = 8;
  n = 8;
  a = rosser1 ( );
  x = rosser1_null_left ( );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ZERO
//
  title = "ZERO";
  m = 5;
  n = 5;
  a = zero ( m, n );
  x = zero_null_left ( m, n );
  error_l2 = r8mat_is_null_left ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( m, x );
  cout << "  " << setw(20) << left << title << right
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

void test_null_right ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_NULL_RIGHT tests right null vectors.
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
  double alpha;
  int col_num;
  double error_l2;
  double f1;
  double f2;
  int m;
  int n;
  double norm_a_frobenius;
  double norm_x_l2;
  double r8_hi;
  double r8_lo;
  int row_num;
  int seed;
  string title;
  double *x;

  cout << "\n";
  cout << "TEST_NULL_RIGHT\n";
  cout << "  A = a test matrix of order M by N\n";
  cout << "  x = an N vector, candidate for a right null vector.\n";
  cout << "\n";
  cout << "  ||A|| = Frobenius norm of A.\n";
  cout << "  ||x|| = L2 norm of x.\n";
  cout << "  ||A*x||/||x|| = L2 norm of A*x over L2 norm of x.\n";
  cout << "\n";
  cout << "  Title                    M     N      "
       << "||A||            ||x||        ||A*x||/||x||\n";
  cout << "\n";
//
//  A123
//
  title = "A123";
  m = 3;
  n = 3;
  a = a123 ( );
  x = a123_null_right ( );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ARCHIMEDES
//
  title = "ARCHIMEDES";
  m = 7;
  n = 8;
  a = archimedes ( );
  x = archimedes_null_right ( );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  CHEBY_DIFF1
//
  title = "CHEBY_DIFF1";
  m = 5;
  n = 5;
  a = cheby_diff1 ( n );
  x = cheby_diff1_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  CREATION
//
  title = "CREATION";
  m = 5;
  n = 5;
  a = creation ( m, n );
  x = creation_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF1
//  Only has null vectors for N odd.
//
  title = "DIF1";
  m = 5;
  n = 5;
  a = dif1 ( m, n );
  x = dif1_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF1CYCLIC
//
  title = "DIF1CYCLIC";
  m = 5;
  n = 5;
  a = dif1cyclic ( n );
  x = dif1cyclic_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  DIF2CYCLIC
//
  title = "DIF2CYCLIC";
  m = 5;
  n = 5;
  a = dif2cyclic ( n );
  x = dif2cyclic_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  FIBONACCI1
//
  title = "FIBONACCI1";
  m = 5;
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  f1 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  f2 = r8_uniform_ab ( r8_lo, r8_hi, seed );
  a = fibonacci1 ( n, f1, f2 );
  x = fibonacci1_null_right ( m, n, f1, f2 );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  HAMMING
//
  title = "HAMMING";
  m = 5;
  n = i4_power ( 2, m ) - 1;
  a = hamming ( m, n );
  x = hamming_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  LINE_ADJ
//
  title = "LINE_ADJ";
  m = 7;
  n = 7;
  a = line_adj ( n );
  x = line_adj_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  MOLER2
//
  title = "MOLER2";
  m = 5;
  n = 5;
  a = moler2 ( );
  x = moler2_null_right ( );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  NEUMANN
//
  title = "NEUMANN";
  row_num = 5;
  col_num = 5;
  m = row_num * col_num;
  n = row_num * col_num;
  a = neumann ( row_num, col_num );
  x = neumann_null_right ( row_num, col_num );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ONE
//
  title = "ONE";
  m = 5;
  n = 5;
  a = one ( n, n );
  x = one_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  RING_ADJ
//  N must be a multiple of 4 for there to be a null vector.
//
  title = "RING_ADJ";
  m = 12;
  n = 12;
  a = ring_adj ( n );
  x = ring_adj_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ROSSER1
//
  title = "ROSSER1";
  m = 8;
  n = 8;
  a = rosser1 ( );
  x = rosser1_null_right ( );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << norm_x_l2
       << "  " << setw(14) << error_l2 << "\n";
  delete [] a;
  delete [] x;
//
//  ZERO
//
  title = "ZERO";
  m = 5;
  n = 5;
  a = zero ( m, n );
  x = zero_null_right ( m, n );
  error_l2 = r8mat_is_null_right ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  cout << "  " << setw(20) << left << title << right
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
//    07 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double error_frobenius;
  int i;
  int i4_hi;
  int i4_lo;
  double *l;
  int m;
  int n;
  double norm_a_frobenius;
  double *p;
  int *pivot;
  double r8_hi;
  double r8_lo;
  int seed;
  string title;
  double *u;
  double *x;

  cout << "\n";
  cout << "TEST_PLU\n";
  cout << "  A = a test matrix of order M by N\n";
  cout << "  P, L, U are the PLU factors.\n";
  cout << "\n";
  cout << "  ||A|| = Frobenius norm of A.\n";
  cout << "  ||A-PLU|| = Frobenius norm of A-P*L*U.\n";
  cout << "\n";
  cout << "  Title                  M     N      ";
  cout << "||A||            ||A-PLU||\n";
  cout << "\n";
//
//  A123
//
  title = "A123";
  m = 3;
  n = 3;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = a123 ( );
  a123_plu ( p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  BODEWIG
//
  title = "BODEWIG";
  m = 4;
  n = 4;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = bodewig ( );
  bodewig_plu ( p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  BORDERBAND
//
  title = "BORDERBAND";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = borderband ( n );
  borderband_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  DIF2
//
  title = "DIF2";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = dif2 ( m, n );
  dif2_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  GFPP
//
  title = "GFPP";
  m = 5;
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = gfpp ( n, alpha );
  gfpp_plu ( n, alpha, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  GIVENS
//
  title = "GIVENS";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = givens ( m, n );
  givens_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  KMS
//
  title = "KMS";
  m = 5;
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = kms ( alpha, m, n );
  kms_plu ( alpha, n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  LEHMER
//
  title = "LEHMER";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = lehmer ( m, n );
  lehmer_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  MAXIJ
//
  title = "MAXIJ";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = maxij ( m, n );
  maxij_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  MINIJ
//
  title = "MINIJ";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = minij ( m, n );
  minij_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  MOLER1
//
  title = "MOLER1";
  m = 5;
  n = 5;
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed );
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = moler1 ( alpha, m, n );
  moler1_plu ( alpha, n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  MOLER3
//
  title = "MOLER3";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = moler3 ( m, n );
  moler3_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  OTO
//
  title = "OTO";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = oto ( m, n );
  oto_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  PASCAL2
//
  title = "PASCAL2";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = pascal2 ( n );
  pascal2_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  PLU
//
  title = "PLU";
  m = 5;
  n = 5;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  pivot = new int[n];
  seed = 123456789;
  for ( i = 0; i < n; i++ )
  {
    i4_lo = i;
    i4_hi = n - 1;
    pivot[i] = i4_uniform_ab ( i4_lo, i4_hi, seed );
  }
  a = plu ( n, pivot );
  plu_plu ( n, pivot, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
//
//  VAND2
//
  title = "VAND2";
  m = 4;
  n = 4;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( m, r8_lo, r8_hi, seed );
  a = vand2 ( m, x );
  vand2_plu ( m, x, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_a_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] l;
  delete [] p;
  delete [] u;
  delete [] x;
//
//  WILSON
//
  title = "WILSON";
  m = 4;
  n = 4;
  p = new double[m*m];
  l = new double[m*m];
  u = new double[m*n];
  a = wilson ( );
  wilson_plu ( p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  cout << "  " << setw(20) << left << title << right
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
  cout << "  Title                    M     N     K      ||A||         ||A*X-B||\n";
  cout << "\n";
//
//  A123
//
  title = "A123";
  m = 3;
  n = 3;
  k = 1;
  a = a123 ( );
  b = a123_rhs ( );
  x = a123_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  BODEWIG
//
  title = "BODEWIG";
  m = 4;
  n = 4;
  k = 1;
  a = bodewig ( );
  b = bodewig_rhs ( );
  x = bodewig_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  DIF2
//
  title = "DIF2";
  m = 10;
  n = 10;
  k = 2;
  a = dif2 ( m, n );
  b = dif2_rhs ( m, k );
  x = dif2_solution ( n, k );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  FRANK
//
  title = "FRANK";
  m = 10;
  n = 10;
  k = 2;
  a = frank ( n );
  b = frank_rhs ( m, k );
  x = frank_solution ( n, k );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  POISSON
//
  title = "POISSON";
  nrow = 4;
  ncol = 5;
  m = nrow * ncol;
  n = nrow * ncol;
  k = 1;
  a = poisson ( nrow, ncol );
  b = poisson_rhs ( nrow, ncol );
  x = poisson_solution ( nrow, ncol );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  WILK03
//
  title = "WILK03";
  m = 3;
  n = 3;
  k = 1;
  a = wilk03 ( );
  b = wilk03_rhs ( );
  x = wilk03_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  WILK04
//
  title = "WILK04";
  m = 4;
  n = 4;
  k = 1;
  a = wilk04 ( );
  b = wilk04_rhs ( );
  x = wilk04_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(4) << k
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;
  delete [] b;
  delete [] x;
//
//  WILSON
//
  title = "WILSON";
  m = 4;
  n = 4;
  k = 1;
  a = wilson ( );
  b = wilson_rhs ( );
  x = wilson_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
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
//****************************************************************************80

void test_type ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_TYPE tests functions which test the type of a matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double error_frobenius;
  int key;
  int m;
  int n;
  double norm_frobenius;
  int seed;
  string title;

  cout << "\n";
  cout << "TEST_TYPE\n";
  cout << "  Demonstrate functions which test the type of a matrix.\n";
//
//  TRANSITION.
//
  cout << "\n";
  cout << "  Title                    M     N     ||A||";
  cout << "            ||Transition Error||\n";
  cout << "\n";

  title = "BODEWIG";
  m = 4;
  n = 4;
  a = bodewig ( );
  error_frobenius = r8mat_is_transition ( m, n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;

  title = "SNAKES";
  m = 101;
  n = 101;
  a = snakes ( );
  error_frobenius = r8mat_is_transition ( m, n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;

  title = "TRANSITION_RANDOM";
  m = 5;
  n = 5;
  key = 123456789;
  a = transition_random ( n, key );
  error_frobenius = r8mat_is_transition ( m, n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  cout << "  " << setw(20) << left << title << right
       << "  " << setw(4) << m
       << "  " << setw(4) << n
       << "  " << setw(14) << norm_frobenius
       << "  " << setw(14) << error_frobenius << "\n";
  delete [] a;

  return;
}
