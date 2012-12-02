# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "walsh.hpp"

//****************************************************************************80

void ffwt ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FFWT performs an in-place fast Walsh transform.
//
//  Discussion:
//
//    This routine performs a fast Walsh transform on an input series X
//    leaving the transformed results in X. 
//    X is dimensioned N, which must be a power of 2.
//    The results of this Walsh transform are in sequency order.
//
//    The output sequence could be normalized by dividing by N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    Ken Beauchamp
//
//  Reference:
//
//    Ken Beauchamp,
//    Walsh functions and their applications,
//    Academic Press, 1975,
//    ISBN: 0-12-084050-2,
//    LC: QA404.5.B33.
//
//  Parameters:
//
//    Input, int N, the number of items in X.
//    N must be a power of 2.
//
//    Input/output, double X[N], the data to be transformed.
//
{
  double hold;
  int i;
  int ii;
  int j;
  int j2;
  int js;
  int k;
  int l;
  int m;
  int mw;
  int mw1;
  int nw;
  int nz;
  int nz2;
  int nzi;
  int nzn;
  int *two_power;
  double z;

  m = i4_log_2 ( n );

  two_power = new int[m];

  for ( i = 0; i < m; i++ )
  {
    two_power[i] = i4_power ( 2, m - 1 - i );
  }

  for ( l = 0; l < m; l++ )
  {
    nz = i4_power ( 2, l );
    nzi = 2 * nz;
    nzn = n / nzi;
    nz2 = nz / 2;
    if ( nz2 == 0 )
    {
      nz2 = 1;
    }

    for ( i = 0; i < nzn; i++ )
    {
      js = i * nzi;
      z = 1.0;
      for ( ii = 0; ii < 2; ii++ )
      {
        for ( j = 0; j < nz2; j++ )
        {
          js = js + 1;
          j2 = js + nz;
          hold = x[js-1] + z * x[j2-1];
          z = - z;
          x[j2-1] = x[js-1] + z * x[j2-1];
          x[js-1] = hold;
          z = - z;
        }
        if ( l == 0 )
        {
          break;
        }
        z = - 1.0;
      }
    }
  }
//
//  Bit reversal section.
//
  nw = 0;
  for ( k = 0; k < n; k++ )
  {
//
//  Choose correct index and switch elements if not already switched.
//
    if ( k < nw )
    {
      hold = x[nw];
      x[nw] = x[k];
      x[k] = hold;
    }
//
//  Bump up series by 1.
//
    for ( i = 0; i < m; i++ )
    {
      ii = i;
      if ( nw < two_power[i] )
      {
        break;
      }
      mw = nw / two_power[i];
      mw1 = mw / 2;
      if ( mw <= 2 * mw1 )
      {
        break;
      }
      nw = nw - two_power[i];
    }
    nw = nw + two_power[ii];
  }

  delete [] two_power;

  return;
}
//****************************************************************************80

void fwt ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FWT performs a fast Walsh transform.
//
//  Discussion:
//
//    This routine performs a fast Walsh transform on an input series X
//    leaving the transformed results in X. 
//    X is dimensioned N, which must be a power of 2.
//    The results of this Walsh transform are in sequency order.
//
//    The output sequence could be normalized by dividing by N.
//
//    Note that the program text in the reference included the line
//      y(jd) = abs ( x(j) - x(j2) )
//    which has been corrected to:
//      y(jd) = x(j) - x(j2)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    Ken Beauchamp
//
//  Reference:
//
//    Ken Beauchamp,
//    Walsh functions and their applications,
//    Academic Press, 1975,
//    ISBN: 0-12-084050-2,
//    LC: QA404.5.B33.
//
//  Parameters:
//
//    Input, int N, the number of items in X.
//    N must be a power of 2.
//
//    Input/output, double X[N], the data to be transformed.
//
{
  int i;
  int j;
  int j2;
  int jd;
  int js;
  int l;
  int m;
  int n2;
  int nx;
  int ny;
  int nz;
  int nzi;
  int nzn;
  double *y;

  y = new double[n];

  n2 = n / 2;
  m = i4_log_2 ( n );

  for ( l = 1; l <= m; l++ )
  {
    ny = 0;
    nz = i4_power ( 2, l - 1 );
    nzi = 2 * nz;
    nzn = n / nzi;
    for ( i = 1; i <= nzn; i++ )
    {
      nx = ny + 1;
      ny = ny + nz;
      js = ( i - 1 ) * nzi;
      jd = js + nzi + 1;
      for ( j = nx; j <= ny; j++ )
      {
        js = js + 1;
        j2 = j + n2;
        y[js-1] = x[j-1] + x[j2-1];
        jd = jd - 1;
        y[jd-1] = x[j-1] - x[j2-1];
      }
    }
    r8vec_copy ( n, y, x );
  }
  delete [] y;

  return;
}
//****************************************************************************80

void haar ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HAAR performs a Haar transform.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    Ken Beauchamp
//
//  Reference:
//
//    Ken Beauchamp,
//    Walsh functions and their applications,
//    Academic Press, 1975,
//    ISBN: 0-12-084050-2,
//    LC: QA404.5.B33.
//
//  Parameters:
//
//    Input, int N, the number of items in X.
//    N must be a power of 2.
//
//    Input/output, double X[N], the data to be transformed.
//
{
  int i;
  int i1;
  int j;
  int jj;
  int k;
  int l;
  int l2;
  int l3;
  double *y;

  y = new double[n];

  k = i4_log_2 ( n );

  for ( i = 1; i <= k; i++ )
  {
    l = k + 1 - i;
    l2 = i4_power ( 2, l - 1 );

    r8vec_copy ( 2 * l2, x, y );

    for ( j = 1; j <= l2; j++ )
    {
       l3 = l2 + j;
       jj = 2 * j - 1;
       x[j-1]  = y[jj-1] + y[jj];
       x[l3-1] = y[jj-1] - y[jj];
    }
  }
  delete [] y;

  return;
}
//****************************************************************************80

void haarin ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HAARIN inverts a Haar transform.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    Ken Beauchamp
//
//  Reference:
//
//    Ken Beauchamp,
//    Walsh functions and their applications,
//    Academic Press, 1975,
//    ISBN: 0-12-084050-2,
//    LC: QA404.5.B33.
//
//  Parameters:
//
//    Input, int N, the number of items in X.
//    N must be a power of 2.
//
//    Input/output, double X[N], the data to be transformed.
//
{
  int i;
  int i1;
  int j;
  int jj;
  int jj1;
  int k;
  int l;
  int lj;
  double *y;

  y = new double[n];

  k = i4_log_2 ( n );

  for ( i = 1; i <= k; i++ )
  {
    l = i4_power ( 2, i - 1 );
    r8vec_copy ( 2 * l, x, y );
    for ( j = 1; j <= l; j++ )
    {
      lj = l + j;
      jj = 2 * j;
      jj1 = jj - 1;
      x[jj-1]  = y[j-1] - y[lj-1];
      x[jj1-1] = y[j-1] + y[lj-1];
    }
  }

  delete [] y;

  return;
}
//****************************************************************************80

void hnorm ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HNORM computes normalization factors for a forward or inverse Haar transform.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    Ken Beauchamp
//
//  Reference:
//
//    Ken Beauchamp,
//    Walsh functions and their applications,
//    Academic Press, 1975,
//    ISBN: 0-12-084050-2,
//    LC: QA404.5.B33.
//
//  Parameters:
//
//    Input, int N, the number of items in X.
//    N must be a power of 2.
//
//    Input/output, double X[N], the data to be transformed.
//
{
  int i;
  int ii;
  int j;
  int jmax;
  int jmin;
  int k;
  double wlk;

  k = i4_log_2 ( n );

  x[0] = x[0] / pow ( 2.0, k );

  if ( 1 <= k )
  {
    x[1] = x[1] / pow ( 2.0, k );
  }

  for ( ii = 2; ii <= k; ii++ )
  {
    i = ii - 1;
    wlk = 1.0 / pow ( 2.0, k - i );
    jmin = i4_power ( 2, i );
    jmax = i4_power ( 2, ii ) - 1;
    for ( j = jmin; j <= jmax; j++ )
    {
      x[j] = x[j] * wlk;
    }
  }
  return;
}
//****************************************************************************80

int i4_log_2 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
//
//  Example:
//
//        I  I4_LOG_10
//    -----  --------
//        0    0
//        1    0
//        2    1
//        3    1
//        4    2
//        5    2
//        7    2
//        8    3
//        9    3
//     1000    9
//     1024   10
//
//  Discussion:
//
//    I4_LOG_2 ( I ) + 1 is the number of binary digits in I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number whose logarithm base 2 is desired.
//
//    Output, int I4_LOG_2, the integer part of the logarithm base 2 of
//    the absolute value of X.
//
{
  int i_abs;
  int two_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    two_pow = 2;

    i_abs = abs ( i );

    while ( two_pow <= i_abs )
    {
      value = value + 1;
      two_pow = two_pow * 2;
    }
  }

  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If 
//      NREM = I4_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
// 
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is 
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

double *r8vec_copy_new ( int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY_NEW copies an R8VEC to a "new" R8VEC.
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
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double R8VEC_COPY_NEW[N], the copy of A1.
//
{
  double *a2;
  int i;

  a2 = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
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

void walsh ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    WALSH performs a fast Walsh transform.
//
//  Discussion:
//
//    This routine performs a fast Wash transform on an input series X
//    leaving the transformed results in X.  The array Y is working space.
//    X and Y are dimensioned N, which must be a power of 2.
//    The results of this Walsh transform are in sequency order.
//
//    The output sequence could be normalized by dividing by N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    Ken Beauchamp
//
//  Reference:
//
//    Ken Beauchamp,
//    Walsh functions and their applications,
//    Academic Press, 1975,
//    ISBN: 0-12-084050-2,
//    LC: QA404.5.B33.
//
//  Parameters:
//
//    Input, int N, the number of items in X.
//    N must be a power of 2.
//
//    Input/output, double X[N], the data to be transformed.
//
{
  double a;
  int i;
  int i1;
  int is;
  int j;
  int j1;
  int l;
  int m;
  int n1;
  int n2;
  double w;
  double *y;
  double z;

  n2 = n / 2;
  y = new double[n2];
  m = i4_log_2 ( n );
  z = - 1.0;

  for ( j = 1; j <= m; j++ )
  {
    n1 = i4_power ( 2, m - j + 1 );
    j1 = i4_power ( 2, j - 1 );
    for ( l = 1; l <= j1; l++ )
    {
      is = ( l - 1 ) * n1 + 1;
      i1 = 0;
      w = z;
      for ( i = is; i <= is + n1 - 1; i = i + 2 )
      {
        a = x[i-1];
        x[is+i1-1] = a + x[i];
        i1 = i1 + 1;
        y[i1-1] = ( x[i] - a ) * w;
        w = w * z;
      }
      for ( i = 1; i <= n1 / 2; i++ )
      {
        x[n1/2+is+i-2] = y[i-1];
      }
    }
  }

  delete [] y;

  return;
}
