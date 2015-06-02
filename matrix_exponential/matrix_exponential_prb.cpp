# include <cstdlib>
# include <iostream>
# include <cmath>
# include <complex>

using namespace std;

# include "matrix_exponential.hpp"
# include "test_matrix_exponential.hpp"
# include "c8lib.hpp"
# include "r8lib.hpp"

int main ( );
void matrix_exponential_test01 ( );
void matrix_exponential_test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MATRIX_EXPONENTIAL_PRB.
//
//  Discussion:
//
//    MATRIX_EXPONENTIAL_PRB tests the MATRIX_EXPONENTIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2013
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
  cout << "  The C8LIB and R8LIB libraries are needed.\n";
  cout << "  This test needs the TEST_MATRIX_EXPONENTIAL library.\n";

  matrix_exponential_test01 ( );
  matrix_exponential_test02 ( );
//
//  Terminate.
//
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
//    MATRIX_EXPONENTIAL_TEST01 compares real matrix exponential algorithms.
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
  cout << "  R8MAT_EXPM1 is an equivalent to EXPM\n";
  cout << "  R8MAT_EXPM2 uses a Taylor series approach\n";
  cout << "  R8MAT_EXPM3 relies on an eigenvalue calculation.\n";

  test_num = r8mat_exp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    cout << "\n";
    cout << "  Test #" << test << "\n";

    r8mat_exp_story ( test );

    n = r8mat_exp_n ( test );

    cout << "  Matrix order N = " << n << "\n";

    a = r8mat_exp_a ( test, n );

    r8mat_print ( n, n, a, "  Matrix:" );

    a_exp = r8mat_expm1 ( n, a );
    r8mat_print ( n, n, a_exp, "  R8MAT_EXPM1(A):" );
    delete [] a_exp;

    a_exp = r8mat_expm2 ( n, a );
    r8mat_print ( n, n, a_exp, "  R8MAT_EXPM2(A):" );
    delete [] a_exp;
//
//    a_exp = r8mat_expm3 ( n, a );
//    r8mat_print ( n, n, a_exp, "  R8MAT_EXPM3(A):" );
//    delete [] a_exp;
//
    a_exp = r8mat_exp_expa ( test, n );
    r8mat_print ( n, n, a_exp, "  Exact Exponential:" );
    delete [] a_exp;

    delete [] a;
  }

  return;
}
//****************************************************************************80

void matrix_exponential_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    MATRIX_EXPONENTIAL_TEST02 compares complex matrix exponential algorithms.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2013
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  complex <double> *a_exp;
  int n;
  int test;
  int test_num;

  cout << "\n";
  cout << "MATRIX_EXPONENTIAL_TEST02:\n";
  cout << "  EXPM is MATLAB's matrix exponential function\n";
  cout << "  C8MAT_EXPM1 is an equivalent to EXPM\n";
  cout << "  C8MAT_EXPM2 uses a Taylor series approach\n";
  cout << "  C8MAT_EXPM3 relies on an eigenvalue calculation.\n";

  test_num = c8mat_exp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    cout << "\n";
    cout << "  Test #" << test << "\n";

    c8mat_exp_story ( test );

    n = c8mat_exp_n ( test );

    cout << "  Matrix order N = " << n << "\n";

    a = c8mat_exp_a ( test, n );

    c8mat_print ( n, n, a, "  Matrix:" );

    a_exp = c8mat_expm1 ( n, a );
    c8mat_print ( n, n, a_exp, "  C8MAT_EXPM1(A):" );
    delete [] a_exp;

//  a_exp = c8mat_expm2 ( n, a );
//  c8mat_print ( n, n, a_exp, "  C8MAT_EXPM2(A):" );
//  delete [] a_exp;
//
//    a_exp = c8mat_expm3 ( n, a );
//    c8mat_print ( n, n, a_exp, "  C8MAT_EXPM3(A):" );
//    delete [] a_exp;
//
    a_exp = c8mat_exp_expa ( test, n );
    c8mat_print ( n, n, a_exp, "  Exact Exponential:" );
    delete [] a_exp;

    delete [] a;
  }

  return;
}
