# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>

using namespace std;

# include "ihs.hpp"

//****************************************************************************80

void covariance ( int dim_num, int n, int x[], double &average, double &std,
  double &covc )

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
//    Output, double &AVERAGE, the average minimum distance.
//
//    Output, double &STD, the standard deviation of the minimum distances.
//
//    Output, double &COVC, the covariance of the minimum distances.
//
{
  double dist;
  int i;
  int j;
  int k;
  double *mindist;
  const double r8_huge = 1.0E+30;
//
//  Find the minimum distance for each point.
//
  mindist = new double[n];

  for ( i = 0; i < n; i++ )
  {
    mindist[i] = r8_huge;
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
  average = r8vec_average ( n, mindist );
//
//  Compute the standard deviation of the distances.
//
  std = r8vec_std ( n, mindist );
//
//  Compute the covariance.
//
  covc = std / average;

  delete [] mindist;

  return;
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

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
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
//    02 October 2012
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
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void i4mat_transpose_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }
    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i - 1 << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = jlo;
    if ( j2lo < 1 )
    {
      j2lo = 1;
    }
    j2hi = jhi;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j - 1 << ":";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

int *ihs ( int dim_num, int n, int d, int &seed )

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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int IHS[DIM_NUM*N], the points.
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
  const double r8_huge = 1.0E+30;
  int *x;

  avail = new int [ dim_num * n ];
  list = new int [ d * n ];
  point = new int [ dim_num * d * n ];
  x = new int[dim_num*n];

  opt = ( ( double ) n ) /
    pow ( ( double ) n, ( double ) ( 1.0 / ( double ) dim_num ) );
//
//  Pick the first point.
//
  for ( i = 0; i < dim_num; i++ )
  {
    x[i+(n-1)*dim_num] = i4_uniform_ab ( 1, n, seed );
  }
//
//  Initialize AVAIL,
//  and set an entry in a random row of each column of AVAIL to N.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      avail[i+j*dim_num] = j + 1;
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
  for ( count = n - 1; 2 <= count; count-- )
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
        point_index = i4_uniform_ab ( 0, k, seed );
        point[i+k*dim_num] = list[point_index];
        list[point_index] = list[k];
      }
    }
//
//  For each candidate, determine the distance to all the
//  points that have already been selected, and save the minimum value.
//
    min_all = r8_huge;
    best = 0;

    for ( k = 0; k < d * count; k++ )
    {
      min_can = r8_huge;

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

  return x;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

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
//    09 April 2012
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  const int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int i2lo_hi;
  int i2lo_lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
//    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
      output << "  " << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
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
