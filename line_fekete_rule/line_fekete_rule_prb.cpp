# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "line_fekete_rule.hpp"
# include "qr_solve.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int m );
void test02 ( int m );
void test03 ( int m );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LINE_FEKETE_RULE_PRB.
//
//  Discussion:
//
//    LINE_FEKETE_RULE_PRB tests the LINE_FEKETE_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2014
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  int m;
  int m_test[TEST_NUM] = { 5, 11, 21 };
  int test;
  int test_num = TEST_NUM;

  timestamp ( );
  cout << "\n";
  cout << "LINE_FEKETE_RULE_PRB\n";
  cout << "  C version\n";
  cout << "  Test the LINE_FEKETE_RULE library.\n";

  for ( test = 0; test < test_num; test++ )
  {
    m = m_test[test];
    test01 ( m );
  }

  for ( test = 0; test < test_num; test++ )
  {
    m = m_test[test];
    test02 ( m );
  }

  for ( test = 0; test < test_num; test++ )
  {
    m = m_test[test];
    test03 ( m );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "LINE_FEKETE_RULE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
# undef TEST_NUM
}
//****************************************************************************80

void test01 ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 seeks Fekete points in [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alvise Sommariva, Marco Vianello,
//    Computing approximate Fekete points by QR factorizations of Vandermonde 
//    matrices,
//    Computers and Mathematics with Applications,
//    Volume 57, 2009, pages 1324-1336.
//
//  Parameters:
//
//    Input, int M, the dimension of the polynomial space.
//
{
# define N 5001

  double a;
  double b;
  int n = N;
  int nf;
  double *wf;
  double wf_sum;
  double *x;
  double *xf;

  a = -1.0;
  b = +1.0;
  x = r8vec_linspace_new ( n, a, b );

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Seek Fekete points in [" << a << "," << b << "]\n";
  cout << "  using " << n << " equally spaced sample points\n";
  cout << "  for polynomials of degree M = " << m << "\n";
  cout << "  using the monomial basis and uniform weight.\n";

  wf = new double[m];
  xf = new double[m];
  line_fekete_monomial ( m, a, b, n, x, nf, xf, wf );

  cout << "\n";
  cout << "  NF = " << nf << "\n";
  r8vec_print ( nf, xf, "  Estimated Fekete points XF:" );

  wf_sum = r8vec_sum ( nf, wf );
  cout << "\n";
  cout << "  Sum(WF) = " << wf_sum << "\n";

  delete [] wf;
  delete [] x;
  delete [] xf;

  return;
#  undef N
}
//****************************************************************************80

void test02 ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 seeks Fekete points in [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    L Bos, N Levenberg,
//    On the calculation of approximate Fekete points: the univariate case,
//    Electronic Transactions on Numerical Analysis,
//    Volume 30, pages 377-397, 2008.
//
//  Parameters:
//
//    Input, int M, the dimension of the polynomial space.
//
{
# define N 5001

  double a;
  double b;
  int n = N;
  int nf;
  double *wf;
  double wf_sum;
  double *x;
  double *xf;

  a = -1.0;
  b = +1.0;
  x = r8vec_linspace_new ( n, a, b );

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Seek Fekete points in [" << a << "," << b << "]\n";
  cout << "  using " << n << " equally spaced sample points\n";
  cout << "  for polynomials of degree M = " << m << "\n";
  cout << "  with the Chebyshev basis.\n";

  wf = new double[m];
  xf = new double[m];
  line_fekete_chebyshev ( m, a, b, n, x, nf, xf, wf );

  cout << "\n";
  cout << "  NF = " << nf << "\n";
  r8vec_print ( nf, xf, "  Estimated Fekete points XF:" );
  wf_sum = r8vec_sum ( nf, wf );
  cout << "\n";
  cout << "  Sum(WF) = " << wf_sum << "\n";

  delete [] wf;
  delete [] x;
  delete [] xf;

  return;
#  undef N
}
//****************************************************************************80

void test03 ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 seeks Fekete points in [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the polynomial space.
//
{
# define N 5001

  double a;
  double b;
  int n = N;
  int nf;
  double *wf;
  double wf_sum;
  double *x;
  double *xf;

  a = -1.0;
  b = +1.0;
  x = r8vec_linspace_new ( n, a, b );

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Seek Fekete points in [" << a << "," << b << "]\n";
  cout << "  using " << n << " equally spaced sample points\n";
  cout << "  for polynomials of degree M = " << m << "\n";
  cout << "  with the Legendre basis and uniform weight.\n";

  wf = new double[m];
  xf = new double[m];
  line_fekete_legendre ( m, a, b, n, x, nf, xf, wf );

  cout << "\n";
  cout << "  NF = " << nf << "\n";
  r8vec_print ( nf, xf, "  Estimated Fekete points XF:" );
  wf_sum = r8vec_sum ( nf, wf );
  cout << "\n";
  cout << "  Sum(WF) = " << wf_sum << "\n";

  delete [] wf;
  delete [] x;
  delete [] xf;

  return;
#  undef N
}
