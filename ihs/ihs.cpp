# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>

using namespace std;

# include "ihs.hpp"

//****************************************************************************80

void covariance ( int dim_num, int n, int x[], double *average, double *std,
  double *covc )

//****************************************************************************80
//
//  Purpose:
//
//    COVARIANCE does a covariance calculation for IHS solutions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points to be generated.
//
//    Input, int X[M*N], the points.
//
//    Output, double *AVERAGE, the average minimum distance.
//
//    Output, double *STD, the standard deviation of the minimum distances.
//
//    Output, double *COVC, the covariance of the minimum distances.
//
{
  double dist;
  int i;
  int j;
  int k;
  double *mindist;
//
//  Find the minimum distance for each point.
//
  mindist = new double [n];

  for ( i = 0; i < n; i++ )
  {
    mindist[i] = r8_huge ( );
    for ( j = 0; j < n; j++ )
    {
      if ( i != j )
      {
        dist = 0.0;
        for ( k = 0; k < dim_num; k++ )
        {
          dist = dist + ( ( double )
            ( ( x[k+i*dim_num] - x[k+j*dim_num] )
            * ( x[k+i*dim_num] - x[k+j*dim_num] ) ) );
        }
        dist = sqrt ( dist );
        if ( dist < mindist[i] )
        {
          mindist[i] = dist;
        }
      }
    }
  }
//
//  Find the average minimum distance.
//
  *average = r8vec_average ( n, mindist );
//
//  Compute the standard deviation of the distances.
//
  *std = r8vec_std ( n, mindist );
//
//  Compute the covariance.
//
  *covc = *std / *average;

  delete [] mindist;

  return;
}
//****************************************************************************80

char digit_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.
//
//  Example:
//
//     I     C
//   -----  ---
//     0    '0'
//     1    '1'
//   ...    ...
//     9    '9'
//    10    '*'
//   -83    '*'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the digit, which should be between 0 and 9.
//
//    Output, char DIGIT_TO_CH, the appropriate character '0' through '9' or '*'.
//
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
//****************************************************************************80

int get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GET_SEED, a random seed value.
//
{
# define I_MAX 2147483647
  time_t clock;
  int i;
  int ihour;
  int imin;
  int isec;
  int seed;
  struct tm *lt;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  clock = time ( &tloc );
  lt = localtime ( &clock );
//
//  Hours is 1, 2, ..., 12.
//
  ihour = lt->tm_hour;

  if ( 12 < ihour )
  {
    ihour = ihour - 12;
  }
//
//  Move Hours to 0, 1, ..., 11
//
  ihour = ihour - 1;

  imin = lt->tm_min;

  isec = lt->tm_sec;

  seed = isec + 60 * ( imin + 60 * ihour );
//
//  We want values in [1,43200], not [0,43199].
//
  seed = seed + 1;
//
//  Remap ISEED from [1,43200] to [1,IMAX].
//
  seed = ( int )
    ( ( ( double ) seed )
    * ( ( double ) I_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;
# undef I_MAX
}
//****************************************************************************80

int i4_log_10 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_10 returns the whole part of the logarithm base 10 of an integer.
//
//  Discussion:
//
//    It should be the case that 10^I4_LOG_10(I) <= |I| < 10^(I4_LOG_10(I)+1).
//    (except for I = 0).
//
//    The number of decimal digits in I is I4_LOG_10(I) + 1.
//
//  Example:
//
//        I    I4_LOG_10(I)
//
//        0     0
//        1     0
//        2     0
//
//        9     0
//       10     1
//       11     1
//
//       99     1
//      100     2
//      101     2
//
//      999     2
//     1000     3
//     1001     3
//
//     9999     3
//    10000     4
//    10001     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer.
//
//    Output, int I4_LOG_10, the whole part of the logarithm of abs ( I ).
//
{
  int ten_pow;
  int value;

  i = abs ( i );

  ten_pow = 10;
  value = 0;

  while ( ten_pow <= i )
  {
    ten_pow = ten_pow * 10;
    value = value + 1;
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
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
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
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

char *i4_to_s ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_S converts an I4 to a string.
//
//  Example:
//
//    INTVAL  S
//
//         1  1
//        -1  -1
//         0  0
//      1952  1952
//    123456  123456
//   1234567  1234567
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer to be converted.
//
//    Output, char *I4_TO_S, the representation of the integer.
//
{
  int digit;
  int j;
  int length;
  int ten_power;
  char *s;

  length = i4_log_10 ( i );

  ten_power = ( int ) pow ( ( double ) 10, ( double ) length );

  if ( i < 0 )
  {
    length = length + 1;
  }
//
//  Add one position for the trailing null.
//
  length = length + 1;

  s = new char[length];

  if ( i == 0 )
  {
    s[0] = '0';
    s[1] = '\0';
    return s;
  }
//
//  Now take care of the sign.
//
  j = 0;
  if ( i < 0 )
  {
    s[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
//
//  Find the leading digit of I, strip it off, and stick it into the string.
//
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
//
//  Tack on the trailing NULL.
//
  s[j] = '\0';
  j = j + 1;

  return s;
}
//****************************************************************************80

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
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
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );

  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

void ihs ( int dim_num, int n, int d, int *seed, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    IHS implements the improved distributed hypercube sampling algorithm.
//
//  Discussion:
//
//    N Points in a DIM_NUM dimensional Latin hypercube are to be selected.
//
//    Each of the DIM_NUM coordinate dimensions is discretized to the values
//    1 through N.  The points are to be chosen in such a way that
//    no two points have any coordinate value in common.  This is
//    a standard Latin hypercube requirement, and there are many
//    solutions.
//
//    This algorithm differs in that it tries to pick a solution
//    which has the property that the points are "spread out"
//    as evenly as possible.  It does this by determining an optimal
//    even spacing, and using the duplication factor D to allow it
//    to choose the best of the various options available to it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Brian Beachkofski, Ramana Grandhi,
//    Improved Distributed Hypercube Sampling,
//    American Institute of Aeronautics and Astronautics Paper 2002-1274.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points to be generated.
//
//    Input, int D, the duplication factor.  This must
//    be at least 1.  A value of 5 is reasonable.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int X[DIM_NUM*N], the points.
//
{
  int *avail;
  int best;
  int count;
  double dist;
  int i;
  int j;
  int k;
  int *list;
  double min_all;
  double min_can;
  double opt;
  int *point;
  int point_index;

  avail = new int [ dim_num * n ];
  list = new int [ d * n ];
  point = new int [ dim_num * d * n ];

  opt = ( ( double ) n ) /
    pow ( ( double ) n, ( double ) ( 1.0 / ( double ) dim_num ) );
//
//  Pick the first point.
//
  for ( i = 0; i < dim_num; i++ )
  {
    x[i+(n-1)*dim_num] = i4_uniform ( 1, n, seed );
  }
//
//  Initialize AVAIL,
//  and set an entry in a random row of each column of AVAIL to N.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      avail[i+j*dim_num] = j+1;
    }
  }

  for ( i = 0; i < dim_num; i++ )
  {
    avail[i+(x[i+(n-1)*dim_num]-1)*dim_num] = n;
  }
//
//  Main loop:
//  Assign a value to X(1:M,COUNT) for COUNT = N-1 down to 2.
//
  for ( count = n-1; 2 <= count; count-- )
  {
//
//  Generate valid points.
//
    for ( i = 0; i < dim_num; i++ )
    {
      for ( k = 0; k < d; k++ )
      {
        for ( j = 0; j < count; j++ )
        {
          list[j+k*count] = avail[i+j*dim_num];
        }
      }

      for ( k = count*d - 1; 0 <= k; k-- )
      {
        point_index = i4_uniform ( 0, k, seed );
        point[i+k*dim_num] = list[point_index];
        list[point_index] = list[k];
      }
    }
//
//  For each candidate, determine the distance to all the
//  points that have already been selected, and save the minimum value.
//
    min_all = r8_huge ( );
    best = 0;

    for ( k = 0; k < d*count; k++ )
    {
      min_can = r8_huge ( );

      for ( j = count; j < n; j++ )
      {

        dist = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          dist = dist + ( point[i+k*dim_num] - x[i+j*dim_num] )
                      * ( point[i+k*dim_num] - x[i+j*dim_num] );
        }
        dist = sqrt ( dist );

        if ( dist < min_can )
        {
          min_can = dist;
        }
      }

      if ( fabs ( min_can - opt ) < min_all )
      {
        min_all = fabs ( min_can - opt );
        best = k;
      }

    }

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+(count-1)*dim_num] = point[i+best*dim_num];
    }
//
//  Having chosen X(*,COUNT), update AVAIL.
//
    for ( i = 0; i < dim_num; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        if ( avail[i+j*dim_num] == x[i+(count-1)*dim_num] )
        {
          avail[i+j*dim_num] = avail[i+(count-1)*dim_num];
        }
      }
    }
  }
//
//  For the last point, there's only one choice.
//
  for ( i = 0; i < dim_num; i++ )
  {
    x[i+0*dim_num] = avail[i+0*dim_num];
  }

  delete [] avail;
  delete [] list;
  delete [] point;

  return;
}
//****************************************************************************80

void ihs_write ( int dim_num, int n, int d, int seed_init, int seed,
  int r[], char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    IHS_WRITE writes an IHS dataset to a file.
//
//  Discussion:
//
//    The initial lines of the file are comments, which begin with a
//    "#" character.
//
//    Thereafter, each line of the file contains the DIM_NUM-dimensional
//    components of the next entry of the dataset.
//
//    Note that the actual values of the data are integers between 1
//    and N.  For our convenience, these are rescaled by the
//    mapping
//
//      I -> ( 2 * I - 1 )/ ( 2 * N ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int D, the duplication factor.
//
//    Input, int SEED_INIT, the initial random number seed.
//
//    Input, int SEED, the current random number seed.
//
//    Input, int R[DIM_NUM*N], the points.
//
//    Input, char *FILE_OUT_NAME, the name of
//    the output file.
//
{
  ofstream file_out;
  int i;
  int j;
  double x;

  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "IHS_WRITE - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  file_out << "#  " << file_out_name << "\n";
  file_out << "#  created by routine IHS_WRITE.C" << "\n";
  file_out << "#\n";

  file_out << "#  Spatial dimension DIM_NUM = "  << dim_num       << "\n";
  file_out << "#  Number of points N =        "  << n             << "\n";
  file_out << "#  EPSILON (unit roundoff) =   "  << r8_epsilon ( ) << "\n";
  file_out << "#  Duplication factor D =      "  << d             << "\n";
  file_out << "#  Initial SEED_INIT =         "  << seed_init     << "\n";
  file_out << "#  Current SEED =              "  << seed          << "\n";
  file_out << "#\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x = ( double ) ( 2 * r[i+j*dim_num] - 1 ) / ( double ) ( 2 * n );
      file_out << setw(10) << x << "  ";
    }
    file_out << "\n";
  }

  file_out.close ( );

  return;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
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
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 round off unit.
//
//  Discussion:
//
//    R8_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
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
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
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
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate,
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8vec_average ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_AVERAGE returns the average of an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2003
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
//    Output, double R8VEC_AVERAGE, the aveage of the vector.
//
{
  int i;
  double average;

  average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    average = average + a[i];
  }

  average = average / ( ( double ) n );

  return average;
}
//****************************************************************************80

double r8vec_std ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_STD returns the standard deviation of an R8VEC.
//
//  Discussion:
//
//    The standard deviation of a vector X of length N is defined as
//
//      mean ( X(1:n) ) = sum ( X(1:n) ) / n
//
//      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )**2 ) / ( n - 1 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//    N should be at least 2.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_STD, the standard deviation of the vector.
//
{
  double average;
  int i;
  double std;

  if ( n < 2 )
  {
    std = 0.0;
  }
  else
  {
    average = r8vec_average ( n, a );

    std = 0.0;
    for ( i = 0; i < n; i++ )
    {
      std = std + ( a[i] - average ) * ( a[i] - average );
    }
    std = sqrt ( std / ( ( double ) ( n - 1 ) ) );

  }
  return std;
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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
