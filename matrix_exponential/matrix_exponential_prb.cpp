# include <cstdlib>
# include <iostream>
# include <cmath>

using namespace std;

# include "matrix_exponential.hpp"
# include "test_matrix_exponential.hpp"
# include "r8lib.hpp"

int main ( );
void matrix_exponential_test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MATRIX_EXPONENTIAL_TEST tests some matrix exponential algorithms.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "MATRIX_EXPONENTIAL_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the MATRIX_EXPONENTIAL library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test needs the TEST_MATRIX_EXPONENTIAL library.\n";

  matrix_exponential_test01 ( );
/*
  Terminate.
*/
  cout << "\n";
  cout << "MATRIX_EXPONENTIAL_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void matrix_exponential_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    MATRIX_EXPONENTIAL_TEST01 compares matrix exponential algorithms.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_exp;
  int n;
  int test;
  int test_num;

  cout << "\n";
  cout << "MATRIX_EXPONENTIAL_TEST01:\n";
  cout << "  EXPM is MATLAB's matrix exponential function\n";
  cout << "  EXPM11 is an equivalent to EXPM\n";
  cout << "  EXPM2 uses a Taylor series approach\n";
  cout << "  EXPM3 relies on an eigenvalue calculation.\n";

  test_num = mexp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    cout << "\n";
    cout << "  Test #" << test << "\n";

    mexp_story ( test );

    n = mexp_n ( test );

    cout << "  Matrix order N = " << n << "\n";

    a = mexp_a ( test, n );

    r8mat_print ( n, n, a, "  Matrix:" );

    a_exp = expm11 ( n, a );
    r8mat_print ( n, n, a_exp, "  EXPM1(A):" );
    delete [] a_exp;

    a_exp = expm2 ( n, a );
    r8mat_print ( n, n, a_exp, "  EXPM2(A):" );
    delete [] a_exp;
//
//    a_exp = expm3 ( n, a );
//    r8mat_print ( n, n, a_exp, "  EXPM3(A):" );
//    delete [] a_exp;
//
    a_exp = mexp_expa ( test, n );
    r8mat_print ( n, n, a_exp, "  Exact Exponential:" );
    delete [] a_exp;

    delete [] a;
  }

  return;
}
