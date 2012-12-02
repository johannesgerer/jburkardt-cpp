# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "test_matrix_exponential.hpp"
# include "r8lib.hpp"

int main ( );
void test_matrix_exponential_test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_MATRIX_EXPONENTIAL_TEST tests some matrix exponential algorithms.
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
  timestamp ( );
  cout << "\n";
  cout << "TEST_MATRIX_EXPONENTIAL_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_MATRIX_EXPONENTIAL library.\n";
  cout << "  The R8LIB library is needed.\n";

  test_matrix_exponential_test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_MATRIX_EXPONENTIAL_TEST:\n";
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
//    TEST_MATRIX_EXPONENTIAL_TEST01 retrieves the test data.
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

  test_num = mexp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    cout << "\n";
    cout << "  Test #" << test << "\n";

    n = mexp_n ( test );

    mexp_story ( test );

    cout << "\n";
    cout << "  Matrix order N = " << n << "\n";

    a = mexp_a ( test, n );
    r8mat_print ( n, n, a, "  Matrix A:" );

    expa = mexp_expa ( test, n );
    r8mat_print ( n, n, expa, "  Exact Exponential exp(A):" );

    delete [] a;
    delete [] expa;
  }
  return;
}
