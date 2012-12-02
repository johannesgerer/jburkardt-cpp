# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "condition.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CONDITION_PRB.
//
//  Discussion:
//
//    CONDITION_PRB tests the CONDITION library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "CONDITION_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CONDITION library.\n";
  cout << "  This test also requires the R8LIB library.\n";
 
  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CONDITION_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CONDITION_LINPACK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_inverse;
  double a_inverse_norm_l1;
  double *a_lu;
  double a_norm_l1;
  double alpha;
  double beta;
  double cond;
  double cond_l1;
  int i;
  int info;
  int n;
  string name;
  int *pivot;
  int seed;
  double *z;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  CONDITION_LINPACK estimates the L1 condition number.\n";
  cout << "\n";
  cout << "  Matrix               Order   Condition         Linpack\n";
  cout << "\n";
//
//  Combinatorial matrix.
//
  name = "Combinatorial";
  n = 4;
  alpha = 2.0;
  beta = 3.0;
  a = combin ( alpha, beta, n );
  a_inverse = combin_inverse ( alpha, beta, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  CONEX1
//
  name = "CONEX1";
  n = 4;
  alpha = 100.0;
  a = conex1 ( alpha );
  a_inverse = conex1_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  CONEX2
//
  name = "CONEX2";
  n = 3;
  alpha = 100.0;
  a = conex2 ( alpha );
  a_inverse = conex2_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  CONEX3
//
  name = "CONEX3";
  n = 5;
  a = conex3 ( n );
  a_inverse = conex3_inverse ( n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  CONEX4
//
  name = "CONEX4";
  n = 4;
  a = conex4 ( );
  a_inverse = conex4_inverse ( );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  KAHAN
//
  name = "KAHAN";
  n = 4;
  alpha = 0.25;
  a = kahan ( alpha, n, n );
  a_inverse = kahan_inverse ( alpha, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  Random
//
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    name = "RANDOM";
    n = 4;
    a = r8mat_uniform_01_new ( n, n, seed );
    a_lu = r8mat_copy_new ( n, n, a );
    pivot = new int[n];
    info = r8ge_fa ( n, a_lu, pivot );
    a_inverse = r8ge_inverse ( n, a_lu, pivot );
    a_norm_l1 = r8mat_norm_l1 ( n, n, a );
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
    cond_l1 = a_norm_l1 * a_inverse_norm_l1;
    cond = condition_linpack ( n, a );
    cout << "  " << setw(20) << name
         << "  " << setw(4) << n
         << "  " << setw(14) << cond_l1
         << "  " << setw(14) << cond << "\n";
    delete [] a;
    delete [] a_inverse;
    delete [] a_lu;
    delete [] pivot;
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CONDITION_SAMPLE1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_inverse;
  double a_inverse_norm_l1;
  double *a_lu;
  double a_norm_l1;
  double alpha;
  double beta;
  double cond;
  double cond_l1;
  int i;
  int info;
  int j;
  int m;
  int m_test[3] = { 10, 1000, 100000 };
  int n;
  string name;
  int *pivot;
  int seed;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  CONDITION_SAMPLE1 estimates the L1 condition number using sampling.\n";
  cout << "\n";
  cout << "  Matrix                 Samples Order   Condition        Estimate\n";
//
//  Combinatorial matrix.
//
  name = "Combinatorial";
  n = 4;
  alpha = 2.0;
  beta = 3.0;
  a = combin ( alpha, beta, n );
  a_inverse = combin_inverse ( alpha, beta, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    cout << "  " << setw(20) << name
         << "  " << setw(8) << m
         << "  " << setw(4) << n
         << "  " << setw(14) << cond_l1
         << "  " << setw(14) << cond << "\n";
  }
  delete [] a;
  delete [] a_inverse;
//
//  CONEX1
//
  name = "CONEX1";
  n = 4;
  alpha = 100.0;
  a = conex1 ( alpha );
  a_inverse = conex1_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    cout << "  " << setw(20) << name
         << "  " << setw(8) << m
         << "  " << setw(4) << n
         << "  " << setw(14) << cond_l1
         << "  " << setw(14) << cond << "\n";
  }
  delete [] a;
  delete [] a_inverse;
//
//  CONEX2
//
  name = "CONEX2";
  n = 3;
  alpha = 100.0;
  a = conex2 ( alpha );
  a_inverse = conex2_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    cout << "  " << setw(20) << name
         << "  " << setw(8) << m
         << "  " << setw(4) << n
         << "  " << setw(14) << cond_l1
         << "  " << setw(14) << cond << "\n";
  }
  delete [] a;
  delete [] a_inverse;
//
//  CONEX3
//
  name = "CONEX3";
  n = 5;
  a = conex3 ( n );
  a_inverse = conex3_inverse ( n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    cout << "  " << setw(20) << name
         << "  " << setw(8) << m
         << "  " << setw(4) << n
         << "  " << setw(14) << cond_l1
         << "  " << setw(14) << cond << "\n";
  }
  delete [] a;
  delete [] a_inverse;
//
//  CONEX4
//
  name = "CONEX4";
  n = 4;
  a = conex4 ( );
  a_inverse = conex4_inverse ( );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    cout << "  " << setw(20) << name
         << "  " << setw(8) << m
         << "  " << setw(4) << n
         << "  " << setw(14) << cond_l1
         << "  " << setw(14) << cond << "\n";
  }
  delete [] a;
  delete [] a_inverse;
//
//  KAHAN
//
  name = "KAHAN";
  n = 4;
  alpha = 0.25;
  a = kahan ( alpha, n, n );
  a_inverse = kahan_inverse ( alpha, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    cout << "  " << setw(20) << name
         << "  " << setw(8) << m
         << "  " << setw(4) << n
         << "  " << setw(14) << cond_l1
         << "  " << setw(14) << cond << "\n";
  }
  delete [] a;
  delete [] a_inverse;
//
//  Random
//
  seed = 123456789;

  for ( j = 1; j <= 5; j++ )
  {
    name = "RANDOM";
    n = 4;
    a = r8mat_uniform_01_new ( n, n, seed );
    a_lu = r8mat_copy_new ( n, n, a );
    pivot = new int[n];
    info = r8ge_fa ( n, a_lu, pivot );
    a_inverse = r8ge_inverse ( n, a_lu, pivot );
    a_norm_l1 = r8mat_norm_l1 ( n, n, a );
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
    cond_l1 = a_norm_l1 * a_inverse_norm_l1;
    cout << "\n";
    for ( i = 0; i < 3; i++ )
    {
      m = m_test[i];
      cond = condition_sample1 ( n, a, m );
      cout << "  " << setw(20) << name
           << "  " << setw(8) << m
           << "  " << setw(4) << n
           << "  " << setw(14) << cond_l1
           << "  " << setw(14) << cond << "\n";
    }
    delete [] a;
    delete [] a_inverse;
    delete [] a_lu;
    delete [] pivot;
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests CONDITION_HAGER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_inverse;
  double a_inverse_norm_l1;
  double *a_lu;
  double a_norm_l1;
  double alpha;
  double beta;
  double cond;
  double cond_l1;
  int i;
  int info;
  int n;
  string name;
  int *pivot;
  int seed;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  CONDITION_HAGER estimates the L1 condition number.\n";
  cout << "\n";
  cout << "  Matrix               Order   Condition         Hager\n";
  cout << "\n";
//
//  Combinatorial matrix.
//
  name = "Combinatorial";
  n = 4;
  alpha = 2.0;
  beta = 3.0;
  a = combin ( alpha, beta, n );
  a_inverse = combin_inverse ( alpha, beta, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  CONEX1
//
  name = "CONEX1";
  n = 4;
  alpha = 100.0;
  a = conex1 ( alpha );
  a_inverse = conex1_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  CONEX2
//
  name = "CONEX2";
  n = 3;
  alpha = 100.0;
  a = conex2 ( alpha );
  a_inverse = conex2_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  CONEX3
//
  name = "CONEX3";
  n = 5;
  a = conex3 ( n );
  a_inverse = conex3_inverse ( n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  CONEX4
//
  name = "CONEX4";
  n = 4;
  a = conex4 ( );
  a_inverse = conex4_inverse ( );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  KAHAN
//
  name = "KAHAN";
  n = 4;
  alpha = 0.25;
  a = kahan ( alpha, n, n );
  a_inverse = kahan_inverse ( alpha, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  cout << "  " << setw(20) << name
       << "  " << setw(4) << n
       << "  " << setw(14) << cond_l1
       << "  " << setw(14) << cond << "\n";
  delete [] a;
  delete [] a_inverse;
//
//  Random
//
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    name = "RANDOM";
    n = 4;
    a = r8mat_uniform_01_new ( n, n, seed );
    a_lu = r8mat_copy_new ( n, n, a );
    pivot = new int[n];
    info = r8ge_fa ( n, a_lu, pivot );
    a_inverse = r8ge_inverse ( n, a_lu, pivot );
    a_norm_l1 = r8mat_norm_l1 ( n, n, a );
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
    cond_l1 = a_norm_l1 * a_inverse_norm_l1;
    cond = condition_hager ( n, a );
    cout << "  " << setw(20) << name
         << "  " << setw(4) << n
         << "  " << setw(14) << cond_l1
         << "  " << setw(14) << cond << "\n";
    delete [] a;
    delete [] a_inverse;
    delete [] a_lu;
    delete [] pivot;
  }
  return;
}
