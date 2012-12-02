# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>

using namespace std;

# include "sine_transform.hpp"

//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *sine_transform_data ( int n, double d[] )

//****************************************************************************80
//
//  Purpose:
//
//    SINE_TRANSFORM_DATA does a sine transform on a vector of data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of data points.
//
//    Input, double D[N], the vector of data.
//
//    Output, double SINE_TRANSFORM_DATA[N], the sine transform coefficients.
//
{
  double angle;
  int i;
  int j;
  double pi = 3.141592653589793;
  double *s;

  s = new double[n];

  for ( i = 0; i < n; i++ )
  {
    s[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      angle = pi * ( double ) ( ( i + 1 ) * ( j + 1 ) ) / ( double ) ( n + 1 );
      s[i] = s[i] + sin ( angle ) * d[j];
    }
    s[i] = s[i] * sqrt ( 2.0 / ( double ) ( n + 1 ) );
  }
  return s;
}
//****************************************************************************80

double *sine_transform_function ( int n, double a, double b, 
  double f ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    SINE_TRANSFORM_FUNCTION does a sine transform on functional data.
//
//  Discussion:
//
//    The interval [A,B] is divided into N+1 intervals using N+2 points,
//    which are indexed by 0 through N+1.
//
//    The original function F(X) is regarded as the sum of a linear function 
//    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
//    which is 0 at A and B.
//
//    The sine transform coefficients for F2 are then computed.
//
//    To recover the interpolant of F(X), it is necessary to combine the
//    linear part F1 with the sine transform interpolant:
//
//      Interp(F)(X) = F1(X) + F2(X)
//
//    This can be done by calling SINE_TRANSFORM_INTERPOLANT().
//    
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data points.
//
//    Input, double A, B, the interval endpoints.
//
//    Input, double F ( double X ), a pointer to the function.
//
//    Output, SINE_TRANSFORM_FUNCTION[N], the sine transform coefficients.
//
{
  double angle;
  double *f2;
  double fa;
  double fb;;
  int i;
  int j;
  double pi = 3.141592653589793;
  double *s;
  double x;

  fa = f ( a );
  fb = f ( b );

  f2 = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x = ( ( double ) ( n - i     ) * a
        + ( double ) (     i + 1 ) * b )
        / ( double ) ( n     + 1 );

    f2[i] = f ( x )                
         - ( ( b - x     ) * fa   
         +   (     x - a ) * fb ) 
         /   ( b     - a );
  }

  s = new double[n];

  for ( i = 0; i < n; i++ )
  {
    s[i] = 0.0;

    for ( j = 0; j < n; j++ )
    {
      angle = pi * ( double ) ( ( i + 1 ) * ( j + 1 ) ) / ( double ) ( n + 1 );
      s[i] = s[i] + sin ( angle ) * f2[j];
    }
    s[i] = s[i] * sqrt ( 2.0 / ( double ) ( n + 1 ) );
  }

  delete [] f2;

  return s;
}
//****************************************************************************80

double *sine_transform_interpolant ( int n, double a, double b, double fa, 
  double fb, double s[], int nx, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SINE_TRANSFORM_INTERPOLANT evaluates the sine transform interpolant.
//
//  Discussion:
//
//    The interval [A,B] is divided into N+1 intervals using N+2 points,
//    which are indexed by 0 through N+1.
//
//    The original function F(X) is regarded as the sum of a linear function 
//    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
//    which is 0 at A and B.
//
//    The function F2 has been approximated using the sine transform,
//    and the interpolant is then evaluated as:
//
//      Interp(F)(X) = F1(X) + F2(X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of terms in the approximation.
//
//    Input, double A, B, the interval over which the approximant 
//    was defined.
//
//    Input, double FA, FB, the function values at A and B.
//
//    Input, double S[N], the approximant coefficients.
//
//    Input, int NX, the number of evaluation points.
//
//    Input, double X[NX], the evaluation points.
//
//    Output, double SINE_TRANSFORM_INTERPOLANT[NX], the value of the interpolant.
//
{
  double angle;
  double f1;
  double f2;
  int i;
  int j;
  double pi = 3.141592653589793;
  double *value;

  value = new double[nx];

  for ( i = 0; i < nx; i++ )
  {
    f1 = ( ( b - x[i]     ) * fa   
         + (     x[i] - a ) * fb ) 
         / ( b           - a );
    f2 = 0.0;
    for ( j = 0; j < n; j++ )
    {
      angle = ( double ) ( j + 1 ) * ( x[i] - a ) * pi / ( b - a );
      f2 = f2 + s[j] * sin ( angle );
    }
    f2 = f2 * sqrt ( 2.0 / ( double ) ( n + 1 ) );
    value[i] = f1 + f2;
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
