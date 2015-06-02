# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <cstring>

using namespace std;

# include "test_matrix_exponential.hpp"
# include "c8lib.hpp"
# include "r8lib.hpp"

int main ( );
void test_matrix_exponential_test01 ( );
void test_matrix_exponential_test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_MATRIX_EXPONENTIAL_PRB.
//
//  Discussion:
//
//    TEST_MATRIX_EXPONENTIAL_PRB tests the TEST_MATRIX_EXPONENTIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TEST_MATRIX_EXPONENTIAL_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_MATRIX_EXPONENTIAL library.\n";
  cout << "  The C8LIB and R8LIB libraries are needed.\n";

  test_matrix_exponential_test01 ( );
  test_matrix_exponential_test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_MATRIX_EXPONENTIAL_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test_matrix_exponential_test01 ( )

//*****************************************************************************80
//
//  Purpose:
//
//    TEST_MATRIX_EXPONENTIAL_TEST01 retrieves real test data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *expa;
  int n;
  int test;
  int test_num;

  cout << "\n";
  cout << "TEST_MATRIX_EXPONENTIAL_TEST01:\n";
  cout << "  Retrieve the data for each matrix exponential test.\n";

  test_num = r8mat_exp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    cout << "\n";
    cout << "  Test #" << test << "\n";

    n = r8mat_exp_n ( test );

    r8mat_exp_story ( test );

    cout << "\n";
    cout << "  Matrix order N = " << n << "\n";

    a = r8mat_exp_a ( test, n );
    r8mat_print ( n, n, a, "  Matrix A:" );

    expa = r8mat_exp_expa ( test, n );
    r8mat_print ( n, n, expa, "  Exact Exponential exp(A):" );

    delete [] a;
    delete [] expa;
  }
  return;
}
//****************************************************************************80

void test_matrix_exponential_test02 ( )

//*****************************************************************************80
//
//  Purpose:
//
//    TEST_MATRIX_EXPONENTIAL_TEST02 retrieves complex test data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  complex <double> *expa;
  int n;
  int test;
  int test_num;

  cout << "\n";
  cout << "TEST_MATRIX_EXPONENTIAL_TEST02:\n";
  cout << "  Retrieve the data for each matrix exponential test.\n";

  test_num = c8mat_exp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    cout << "\n";
    cout << "  Test #" << test << "\n";

    n = c8mat_exp_n ( test );

    c8mat_exp_story ( test );

    cout << "\n";
    cout << "  Matrix order N = " << n << "\n";

    a = c8mat_exp_a ( test, n );
    c8mat_print ( n, n, a, "  Matrix A:" );

    expa = c8mat_exp_expa ( test, n );
    c8mat_print ( n, n, expa, "  Exact Exponential exp(A):" );

    delete [] a;
    delete [] expa;
  }
  return;
}
