# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "qwv.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QWV_PRB.
//
//  Discussion:
//
//    QWV_PRB tests the QWV library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "QWV_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the QWV library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QWV_PRB:\n";
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
//    TEST01 tests QWV for a Newton-Cotes rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int n;
  double pi = 3.141592653589793;
  double theta;
  double *w;
  double *x;

  a =  0.0;
  b = +1.0;
  n = 5;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Use the Vandermonde procedure to compute the\n";
  cout << "  quadrature weights for a Newton-Cotes rule.\n";
  cout << "  Order N = " << n << "\n";
  cout << "  Interval = [" << a << "," << b << "]\n";
//
//  Set the points.
//
  x = r8vec_even_new ( n, a, b );
  r8vec_print ( n, x, "  Abscissas:" );
//
//  Compute the weights.
//
  w = qwv ( n, a, b, x );

  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
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
//    TEST02 tests QWV for a Clenshaw-Curtis rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int n;
  double pi = 3.141592653589793;
  double theta;
  double *w;
  double *x;

  a = -1.0;
  b = +1.0;
  n = 5;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use the Vandermonde procedure to compute the\n";
  cout << "  quadrature weights for a Clenshaw-Curtis rule.\n";
  cout << "  Order N = " << n << "\n";
  cout << "  Interval is [" << a << "," << b << "]\n";
//
//  Set the points.
//
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    theta = ( double ) ( n - i - 1 ) * pi 
          / ( double ) ( n - 1 );

    x[i] = ( ( 1.0 - cos ( theta ) ) * a   
           + ( 1.0 + cos ( theta ) ) * b ) 
           /   2.0;
  }

  r8vec_print ( n, x, "  Abscissas:" );
//
//  Determine the corresponding weights.
//
  w = qwv ( n, a, b, x );

  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] w;
  delete [] x;

  return;
}
