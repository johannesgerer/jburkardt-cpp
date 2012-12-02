# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "chebyshev.hpp"

int main ( );
void test01 ( );
double f1 ( double x );
double f2 ( double x );
double f3 ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CHEBYSHEV_PRB.
//
//  Discussion:
//
//    CHEBYSHEV_PRB tests the CHEBYSHEV library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "CHEBYSHEV_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CHEBYSHEV library.\n";
 
  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CHEBYSHEV_PRB\n";
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
//    TEST01 tests CHEBYSHEV_COEFFICIENTS and CHEBYSHEV_INTERPOLANT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  double *fc;
  int i;
  int m;
  int n;
  double *x;
  
  cout << "\n";
  cout << "CHEBYSHEV_TEST01\n";
  cout << "  CHEBYSHEV_COEFFICIENTS computes the coefficients of the\n";
  cout << "  Chebyshev interpolant.\n";
  cout << "  CHEBYSHEV_INTERPOLANT evaluates the interpolant.\n";

  n = 5;
  a = -1.0;
  b = +1.0;

  c = chebyshev_coefficients ( a, b, n, f1 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  cout << "\n";
  cout << "  F(X) is a trig function:\n";
  cout << "\n";
  cout << "          X               C(I)            F(X)           C(F)(X)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << x[i]
         << "  " << setw(14) << c[i]
         << "  " << setw(14) << f1 ( x[i] )
         << "  " << setw(14) << fc[i] << "\n";
  }

  delete [] c;
  delete [] fc;
  delete [] x;
//
//  Try a variant interval.
//
  n = 5;
  a = 0.0;
  b = +3.0;

  c = chebyshev_coefficients ( a, b, n, f1 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  cout << "\n";
  cout << "  Consider the same F(X), but now over [0,3]:\n";
  cout << "\n";
  cout << "          X               C(I)            F(X)           C(F)(X)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << x[i]
         << "  " << setw(14) << c[i]
         << "  " << setw(14) << f1 ( x[i] )
         << "  " << setw(14) << fc[i] << "\n";
  }

  delete [] c;
  delete [] fc;
  delete [] x;
//
//  Try a higher order.
//
  n = 10;
  a = -1.0;
  b = +1.0;

  c = chebyshev_coefficients ( a, b, n, f1 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  cout << "\n";
  cout << "  Consider the same F(X), but now with higher order:\n";
  cout << "\n";
  cout << "          X               C(I)            F(X)           C(F)(X)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << x[i]
         << "  " << setw(14) << c[i]
         << "  " << setw(14) << f1 ( x[i] )
         << "  " << setw(14) << fc[i] << "\n";
  }

  delete [] c;
  delete [] fc;
  delete [] x;
//
//  Try a polynomial.
//
  n = 10;
  a = -1.0;
  b = +1.0;

  c = chebyshev_coefficients ( a, b, n, f3 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  cout << "\n";
  cout << "  F(X) is a degree 4 polynomial:\n";
  cout << "\n";
  cout << "          X               C(I)            F(X)           C(F)(X)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << x[i]
         << "  " << setw(14) << c[i]
         << "  " << setw(14) << f3 ( x[i] )
         << "  " << setw(14) << fc[i] << "\n";
  }

  delete [] c;
  delete [] fc;
  delete [] x;
//
//  Try a function with decaying behavior.
//
  n = 10;
  a = -1.0;
  b = +1.0;

  c = chebyshev_coefficients ( a, b, n, f2 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  cout << "\n";
  cout << "  The polynomial approximation to F(X) decays:\n";
  cout << "\n";
  cout << "          X               C(I)            F(X)           C(F)(X)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << x[i]
         << "  " << setw(14) << c[i]
         << "  " << setw(14) << f2 ( x[i] )
         << "  " << setw(14) << fc[i] << "\n";
  }

  delete [] c;
  delete [] fc;
  delete [] x;

  return;
}
//****************************************************************************80

double f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F1 evaluates a function that can be used for Chebyshev interpolation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, a point where the function is to be evaluated.
//
//    Output, double F1, the function value.
//
{
  double pi = 3.141592653589793;
  double value;

  value = cos ( 2.0 * pi * x ) * sin ( 3.0 * pi * x );

  return value;
}
//****************************************************************************80

double f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F2 evaluates a function that can be used for Chebyshev interpolation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input,  double X, a point where the function is to be evaluated.
//
//    Output, double F2, the function value.
//
{
  double value;

  value = exp ( x );

  return value;
}
//****************************************************************************80

double f3 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F3 evaluates a function that can be used for Chebyshev interpolation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, a point where the function is to be evaluated.
//
//    Output, double F3, the function values.
//
{
  double value;

  value = ( x - 3.0 ) * ( x - 1.0 ) * ( x - 1.0 ) * ( x + 2.0 );

  return value;
}
