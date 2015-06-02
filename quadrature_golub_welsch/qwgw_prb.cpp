# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "qwgw.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QWGW_PRB.
//
//  Discussion:
//
//    QWGW_PRB tests the QWGW library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "QWGW_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the QWGW library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QWGW_PRB:\n";
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
//    TEST01 tests QWGW for the Chebyshev Type 1 weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *aj;
  double b;
  double *bj;
  int j;
  double mu0;
  int n;
  const double pi = 3.141592653589793;
  double *w;
  double *x;
//
//  Set the quadrature interval and number of points.
//
  a = -1.0;
  b = +1.0;
  n = 5;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the Chebyshev Type 1 weight w(x) = 1/sqrt(1-x^2).\n";
  cout << "  Order N = " << n << "\n";
  cout << "  Interval = [" << a << "," << b << "]\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    aj[j] = 0.0;
  }

  bj[0] = 1.0 / 2.0;
  for ( j = 1; j < n - 1; j++ )
  {
    bj[j] = 1.0 / 4.0;
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = pi;
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests QWGW for the Chebyshev Type 2 weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *aj;
  double b;
  double *bj;
  int j;
  double mu0;
  int n;
  const double pi = 3.141592653589793;
  double *w;
  double *x;
//
//  Set the quadrature interval and number of points.
//
  a = -1.0;
  b = +1.0;
  n = 5;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the Chebyshev Type 2 weight w(x) = sqrt(1-x^2).\n";
  cout << "  Order N = " << n << "\n";
  cout << "  Interval = [" << a << "," << b << "]\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    aj[j] = 0.0;
  }

  for ( j = 0; j < n - 1; j++ )
  {
    bj[j] = 1.0 / 4.0;
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = pi / 2.0;
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests QWGW for the Gegenbauer weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *aj;
  double alpha;
  double b;
  double *bj;
  int j;
  double jr;
  double mu0;
  int n;
  const double pi = 3.141592653589793;
  double *w;
  double *x;
//
//  Set the quadrature interval and number of points.
//
  a = -1.0;
  b = +1.0;
  n = 5;
//
//  Set the weight function parameter.
//
  alpha = 0.25;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the Gegenbauer weight w(x) = (1-x^2)^alpha.\n";
  cout << "  Order N = " << n << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  Interval = [" << a << "," << b << "]\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    aj[j] = 0.0;
  }

  for ( j = 0; j < n - 1; j++ )
  {
    jr = ( double ) ( j + 1 );
    bj[j] = ( jr * ( 2.0 * alpha + jr ) ) 
      / ( 4.0 * pow ( alpha + jr, 2 ) - 1.0 );
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = r8_gamma ( alpha + 1.0 ) * r8_gamma ( 0.5 ) / r8_gamma ( alpha + 1.5 );
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests QWGW for the generalized Hermite weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *aj;
  double alpha;
  double *bj;
  int j;
  double jr;
  double mu0;
  int n;
  const double pi = 3.141592653589793;
  double *w;
  double *x;
//
//  The quadrature interval is (-oo,+oo).
//  Set the number of points.
//
  n = 5;
//
//  Set the weight function parameter.
//
  alpha = 2.0;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the generalized Hermite weight w(x) = |x|^alpha * exp(-x^2).\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  Order N = " << n << "\n";
  cout << "  Interval = (-oo,+oo)\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    aj[j] = 0.0;
  }

  for ( j = 0; j < n - 1; j++ )
  {
    jr = ( double ) ( j + 1 );
    if ( ( j % 2 ) == 0 )
    {
      bj[j] = ( jr + alpha ) / 2.0;
    }
    else
    {
      bj[j] = jr / 2.0;
    }
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = r8_gamma ( ( alpha + 1.0 ) / 2.0 );
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests QWGW for the generalized Laguerre weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *aj;
  double alpha;
  double *bj;
  int j;
  double jr;
  double mu0;
  int n;
  const double pi = 3.141592653589793;
  double *w;
  double *x;
//
//  The quadrature interval is [0,+oo).
//  Set the number of points.
//
  a = 0.0;
  n = 5;
//
//  Set the weight function parameter.
//
  alpha = 2.0;

  cout << "\n";
  cout << "TEST05:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the generalized Laguerre weight w(x) = x^alpha * exp(-x).\n";
  cout << "  Order N = " << n << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  Interval = [0,+oo)\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    jr = ( double ) ( j + 1 );
    aj[j] = alpha + 2.0 * jr - 1.0;
  }

  for ( j = 0; j < n - 1; j++ )
  {
    jr = ( double ) ( j + 1 );
    bj[j] = jr * ( alpha + jr );
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = r8_gamma ( alpha + 1.0 );
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests QWGW for the Hermite weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *aj;
  double *bj;
  int j;
  double jr;
  double mu0;
  int n;
  const double pi = 3.141592653589793;
  double *w;
  double *x;
//
//  The quadrature interval is (-oo,+oo).
//  Set the number of points.
//
  n = 5;

  cout << "\n";
  cout << "TEST06:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the Hermite weight w(x) = exp(-x^2).\n";
  cout << "  Order N = " << n << "\n";
  cout << "  Interval = (-oo,+oo)\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    aj[j] = 0.0;
  }

  for ( j = 0; j < n - 1; j++ )
  {
    jr = ( double ) ( j + 1 );
    bj[j] = jr / 2.0;
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = sqrt ( pi );
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests QWGW for the Jacobi weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *aj;
  double alpha;
  double b;
  double beta;
  double *bj;
  int j;
  double jr;
  double mu0;
  int n;
  double *w;
  double *x;
//
//  Set the quadrature interval and number of points.
//
  a = -1.0;
  b = +1.0;
  n = 5;
//
//  Set the weight function parameters.
//
  alpha = 0.25;
  beta = 0.75;

  cout << "\n";
  cout << "TEST07:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the Jacobi weight w(x) = (1-x^2)^alpha*(1+x)^beta\n";
  cout << "  Order N = " << n << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA =  " << beta << "\n";
  cout << "  Interval = [" << a << "," << b << "]\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    jr = ( double ) ( j + 1 );
    aj[j] = ( beta - alpha ) * ( beta + alpha ) 
      / ( alpha + beta + 2.0 * jr - 2.0 ) 
      / ( alpha + beta + 2.0 * jr );
  }

  for ( j = 0; j < n - 1; j++ )
  {
    jr = ( double ) ( j + 1 );
    bj[j] = 4.0 * jr * ( alpha + jr ) * ( beta + jr ) 
      * ( alpha + beta + jr ) 
      / ( pow ( alpha + beta + 2.0 * jr, 2 ) - 1.0 ) 
      / pow ( alpha + beta + 2.0 * jr, 2 );
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = pow ( 2.0, alpha + beta + 1.0 ) 
    * r8_gamma ( alpha + 1.0 ) * r8_gamma ( beta + 1.0 ) 
    / r8_gamma ( alpha + beta + 2.0 );
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests QWGW for the Laguerre weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *aj;
  double *bj;
  int j;
  double jr;
  double mu0;
  int n;
  const double pi = 3.141592653589793;
  double *w;
  double *x;
//
//  The quadrature interval is [a,+oo).
//  Set the number of points.
//
  a = 0.0;
  n = 5;

  cout << "\n";
  cout << "TEST08:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the Laguerre weight w(x) = exp(-x).\n";
  cout << "  Order N = " << n << "\n";
  cout << "  Interval = [0,+oo)\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    jr = ( double ) ( j + 1 );
    aj[j] = 2.0 * jr - 1.0;
  }

  for ( j = 0; j < n - 1; j++ )
  {
    jr = ( double ) ( j + 1 );
    bj[j] = jr * jr;
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = 1.0;
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests QWGW for the Legendre weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *aj;
  double b;
  double *bj;
  int j;
  double jr;
  double mu0;
  int n;
  const double pi = 3.141592653589793;
  double *w;
  double *x;
//
//  Set the quadrature interval and number of points.
//
  a = -1.0;
  b = +1.0;
  n = 5;

  cout << "\n";
  cout << "TEST09:\n";
  cout << "  Compute points and weights for Gauss quadrature\n";
  cout << "  with the Legendre weight w(x) = 1.\n";
  cout << "  Order N = " << n << "\n";
  cout << "  Interval = [" << a << "," << b << "]\n";
//
//  Set the recursion coefficients.
//
  aj = new double[n];
  bj = new double[n];

  for ( j = 0; j < n; j++ )
  {
    aj[j] = 0.0;
  }

  for ( j = 0; j < n - 1; j++ )
  {
    jr = ( double ) ( j + 1 );
    bj[j] = jr * jr / ( 4.0 * jr * jr - 1.0 );
  }
  bj[n-1] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    bj[j] = sqrt ( bj[j] );
  }

  mu0 = 2.0;
//
//  Compute the points and weights.
//
  x = new double[n];
  w = new double[n];

  sgqf ( n, aj, bj, mu0, x, w );

  r8vec_print ( n, x, "  Abscissas:" );
  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] aj;
  delete [] bj;
  delete [] w;
  delete [] x;

  return;
}
