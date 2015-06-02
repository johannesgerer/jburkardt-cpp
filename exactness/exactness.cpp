# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "exactness.hpp"

//****************************************************************************80

void chebyshev1_exactness ( int n, double x[], double w[], int p_max )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_EXACTNESS: monomial exactness for the Chebyshev1 integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points in the rule.
//
//    Input, double X[N], the quadrature points.
//
//    Input, double W[N], the quadrature weights.
//
//    Input, int P_MAX, the maximum exponent.
//    0 <= P_MAX.
//
{
  double e;
  int i;
  int p;
  double q;
  double s;
  double *v;

  cout << "\n";
  cout << "  Quadrature rule for the Chebyshev1 integral.\n";
  cout << "  Rule of order N = " << n << "\n";
  cout << "  Degree          Relative Error\n";
  cout << "\n";

  v = new double[n];

  for ( p = 0; p <= p_max; p++ )
  {
    s = chebyshev1_integral ( p );

    for ( i = 0; i < n; i++ )
    {
      v[i] = pow ( x[i], p );
    }

    q = r8vec_dot_product ( n, w, v );

    if ( s == 0.0 )
    {
      e = fabs ( q );
    }
    else
    {
      e = fabs ( q - s ) / fabs ( s );
    }

    cout << setw(6) << p << "  "
         << setw(24) << e << "\n";
  }

  delete [] v;

  return;
}
//****************************************************************************80

double chebyshev1_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_INTEGRAL evaluates a monomial Chebyshev type 1 integral.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x <= +1 ) x^n / sqrt ( 1 - x^2 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//
//    Output, double CHEBYSHEV1_INTEGRAL, the value of the exact integral.
//
{
  double bot;
  double exact;
  int i;
  const double r8_pi = 3.141592653589793;
  double top;
//
//  Get the exact value of the integral.
//
  if ( ( expon % 2 ) == 0 )
  {
    top = 1;
    bot = 1;
    for ( i = 2; i <= expon; i = i + 2 )
    {
      top = top * ( i - 1 );
      bot = bot *   i;
    }
	
    exact = r8_pi * ( double ) ( top ) / ( double ) ( bot );
  }
  else
  {
    exact = 0.0;	
  }

  return exact;
}
//****************************************************************************80

void chebyshev2_exactness ( int n, double x[], double w[], int p_max )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2_EXACTNESS: monomial exactness for the Chebyshev2 integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points in the rule.
//
//    Input, double X[N], the quadrature points.
//
//    Input, double W[N], the quadrature weights.
//
//    Input, int P_MAX, the maximum exponent.
//    0 <= P_MAX.
//
{
  double e;
  int i;
  int p;
  double q;
  double s;
  double *v;

  cout << "\n";
  cout << "  Quadrature rule for the Chebyshev2 integral.\n";
  cout << "  Rule of order N = " << n << "\n";
  cout << "  Degree          Relative Error\n";
  cout << "\n";

  v = new double[n];

  for ( p = 0; p <= p_max; p++ )
  {
    s = chebyshev2_integral ( p );

    for ( i = 0; i < n; i++ )
    {
      v[i] = pow ( x[i], p );
    }

    q = r8vec_dot_product ( n, w, v );

    if ( s == 0.0 )
    {
      e = fabs ( q );
    }
    else
    {
      e = fabs ( q - s ) / fabs ( s );
    }

    cout << setw(6) << p << "  "
         << setw(24) << e << "\n";
  }

  delete [] v;

  return;
}
//****************************************************************************80

double chebyshev2_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2_INTEGRAL evaluates a monomial Chebyshev type 2 integral.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x <= +1 ) x^n * sqrt ( 1 - x^2 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//
//    Output, double CHEBYSHEV2_INTEGRAL, the value of the exact integral.
//
{
  double bot;
  double exact;
  int i;
  const double r8_pi = 3.141592653589793;
  double top;
//
//  Get the exact value of the integral.
//
  if ( ( expon % 2 ) == 0 )
  {
    top = 1;
    bot = 1;
    for ( i = 2; i <= expon; i = i + 2 )
    {
      top = top * ( i - 1 );
      bot = bot *   i;
    }

    bot = bot * ( double ) ( expon + 2 );

    exact = r8_pi * ( double ) ( top ) / ( double ) ( bot );
  }
  else
  {
    exact = 0.0;
  }
  return exact;
}
//****************************************************************************80

void hermite_exactness ( int n, double x[], double w[], int p_max )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_EXACTNESS investigates exactness of Hermite quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points in the rule.
//
//    Input, double X[N], the quadrature points.
//
//    Input, double W[N], the quadrature weights.
//
//    Input, int P_MAX, the maximum exponent.
//    0 <= P_MAX.
//
{
  double e;
  int i;
  int p;
  double q;
  double s;
  double *v;

  cout << "\n";
  cout << "  Quadrature rule for the Hermite integral.\n";
  cout << "  Rule of order N = " << n << "\n";
  cout << "  Degree          Relative Error\n";
  cout << "\n";

  v = new double[n];

  for ( p = 0; p <= p_max; p++ )
  {
    s = chebyshev1_integral ( p );

    for ( i = 0; i < n; i++ )
    {
      v[i] = pow ( x[i], p );
    }

    q = r8vec_dot_product ( n, w, v );

    if ( s == 0.0 )
    {
      e = fabs ( q );
    }
    else
    {
      e = fabs ( q - s ) / fabs ( s );
    }

    cout << setw(6) << p << "  "
         << setw(24) << e << "\n";
  }

  delete [] v;

  return;
}
//****************************************************************************80

double hermite_integral ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_INTEGRAL evaluates a monomial Hermite integral.
//
//  Discussion:
//
//    Integral ( -oo < x < oo ) x^p exp(-x^2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the exponent of the monomial.  
//    0 <= P.
//
//    Output, double HERMITE_INTEGRAL, the value of the integral.
//
{
  const double r8_pi = 3.141592653589793;
  double value;

  if ( ( p % 2 ) == 0 )
  {
    value = r8_factorial2 ( p - 1 ) * sqrt ( r8_pi ) / pow ( 2.0, p / 2 );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

void laguerre_exactness ( int n, double x[], double w[], int p_max )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_EXACTNESS investigates exactness of Laguerre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points in the rule.
//
//    Input, double X[N], the quadrature points.
//
//    Input, double W[N], the quadrature weights.
//
//    Input, int P_MAX, the maximum exponent.
//    0 <= P_MAX.
//
{
  double e;
  int i;
  int p;
  double q;
  double s;
  double *v;

  cout << "\n";
  cout << "  Quadrature rule for the Laguerre integral.\n";
  cout << "  Rule of order N = " << n << "\n";
  cout << "  Degree          Relative Error\n";
  cout << "\n";

  v = new double[n];

  for ( p = 0; p <= p_max; p++ )
  {
    s = laguerre_integral ( p );

    for ( i = 0; i < n; i++ )
    {
      v[i] = pow ( x[i], p );
    }

    q = r8vec_dot_product ( n, w, v );

    if ( s == 0.0 )
    {
      e = fabs ( q );
    }
    else
    {
      e = fabs ( q - s ) / fabs ( s );
    }

    cout << setw(6) << p << "  "
         << setw(24) << e << "\n";
  }

  delete [] v;

  return;
}
//****************************************************************************80

double laguerre_integral ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.
//
//  Discussion:
//
//    The integral being computed is
//
//      integral ( 0 <= x < +oo ) x^p * exp ( -x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the exponent.
//    0 <= P.
//
//    Output, double LAGUERRE_INTEGRAL, the value of the integral.
//
{
  double s;

  s = r8_factorial ( p );

  return s;
}
//****************************************************************************80

void legendre_exactness ( int n, double x[], double w[], int p_max )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_EXACTNESS investigates exactness of Legendre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points in the rule.
//
//    Input, double X[N], the quadrature points.
//
//    Input, double W[N], the quadrature weights.
//
//    Input, int P_MAX, the maximum exponent.
//    0 <= P_MAX.
//
{
  double e;
  int i;
  int p;
  double q;
  double s;
  double *v;

  cout << "\n";
  cout << "  Quadrature rule for the Legendre integral.\n";
  cout << "  Rule of order N = " << n << "\n";
  cout << "  Degree          Relative Error\n";
  cout << "\n";

  v = new double[n];

  for ( p = 0; p <= p_max; p++ )
  {
    s = legendre_integral ( p );

    for ( i = 0; i < n; i++ )
    {
      v[i] = pow ( x[i], p );
    }

    q = r8vec_dot_product ( n, w, v );

    if ( s == 0.0 )
    {
      e = fabs ( q );
    }
    else
    {
      e = fabs ( q - s ) / fabs ( s );
    }

    cout << setw(6) << p << "  "
         << setw(24) << e << "\n";
  }

  delete [] v;

  return;
}
//****************************************************************************80

double legendre_integral ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.
//
//  Discussion:
//
//    Integral ( -1 <= x <= +1 ) x^p dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the exponent.
//    0 <= P.
//
//    Output, double LEGENDRE_INTEGRAL, the value of the exact integral.
//
{
  double s;

  if ( ( p % 2 ) == 0 )
  {
    s = 2.0 / ( double ) ( p + 1 );
  }
  else
  {
    s = 0.0;
  }

  return s;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
//****************************************************************************80

double r8_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL2 computes the double factorial function.
//
//  Discussion:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N Value
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial
//    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
//
//    Output, double R8_FACTORIAL2, the value of Factorial2(N).
//
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
