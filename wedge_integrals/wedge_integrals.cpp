# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "wedge_integrals.hpp"

//****************************************************************************80

double *monomial_value ( int m, int n, int e[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= i <= m ) x(i)^e(i)
//
//    The combination 0.0^0 is encountered is treated as 1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
//    Input, int E[M], the exponents.
//
//    Input, double X[M*N], the point coordinates.
//
//    Output, double MONOMIAL_VALUE[N], the monomial values.
//
{
  int i;
  int j;
  double *v;

  v = new double[n];
  for ( j = 0; j < n; j++)
  {
    v[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    if ( 0 != e[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        v[j] = v[j] * pow ( x[i+j*m], e[i] );
      }
    }
  }

  return v;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int &seed )

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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
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
//****************************************************************************80

double wedge01_integral ( int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE01_INTEGRAL returns the integral of a monomial in the unit wedge in 3D.
//
//  Discussion:
//
//    This routine returns the integral of
//
//      product ( 1 <= I <= 3 ) X(I)^E(I)
//
//    over the unit wedge.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1
//      -1 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int E[3], the exponents.
//
//    Output, double WEDGE01_INTEGRAL, the integral of the monomial.
//
{
  int i;
  int k;
  double value;

  value = 1.0;

  k = e[0];

  for ( i = 1; i <= e[1]; i++ )
  {
    k = k + 1;
    value = value * ( double ) ( i ) / ( double ) ( k );
  }

  k = k + 1;
  value = value / ( double ) ( k );

  k = k + 1;
  value = value / ( double ) ( k );
//
//  Now account for integration in Z.
//
  if ( e[2] == - 1 )
  {
    cerr << "\n";
    cerr << "WEDGE01_INTEGRAL - Fatal error!\n";
    cerr << "  E(3) = -1 is not a legal input.\n";
    exit ( 1 );
  }
  else if ( ( e[2] % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = value * 2.0 / ( double ) ( e[2] + 1 );
  }

  return value;
}
//****************************************************************************80

double *wedge01_sample ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//   WEDGE01_SAMPLE samples points uniformly from the unit wedge in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double WEDGE01_SAMPLE[3*N], the points.
//
{
  double *e;
  double e_sum;
  int i;
  int j;
  double *x;

  x = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    e = r8vec_uniform_01_new ( 4, seed );

    for ( i = 0; i < 3; i++ )
    {
      e[i] = - log ( e[i] );
    }

    e_sum = 0.0;
    for ( i = 0; i < 3; i++ )
    {
      e_sum = e_sum + e[i];
    }

    x[0+j*3] = e[0] / e_sum;
    x[1+j*3] = e[1] / e_sum;
    x[2+j*3] = 2.0 * e[3] - 1.0;

    delete [] e;
  }

  return x;
}
//****************************************************************************80

double wedge01_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE01_VOLUME returns the volume of the unit wedge in 3D.
//
//  Discussion:
//
//    The unit wedge is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1
//      -1 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double WEDGE01_VOLUME, the volume of the unit wedge.
//
{
  static double value = 1.0;

  return value;
}
