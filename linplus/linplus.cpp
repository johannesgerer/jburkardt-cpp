# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

using namespace std;

# include "linplus.hpp"

//****************************************************************************80

void c8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_PRINT prints a double complex vector.
//
//  Discussion:
//
//    The complex vector of length N is stored as N pairs of real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[2*N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int k;

  k = 0;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  cout << "                  Real     Imaginary\n";
  cout << "                  Part     Part\n";
  cout << "\n";
  for ( i = 0; i <= n-1; i++ ) 
  {
    cout << setw(6)  << i + 1  << "  " 
         << setw(14) << a[k]   << "  "
         << setw(14) << a[k+1] << "\n";
    k = k + 2;
  }

  return;
}
//****************************************************************************80

void c8vec_sort_a2 ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_SORT_A2 ascending sorts a double complex array by L2 norm.
//
//  Discussion:
//
//    The double complex vector of length N is stored as N pairs of real values.
//
//    The L2 norm of A+Bi is sqrt ( A * A + B * B ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input/output, double X[2*N].
//    On input, an unsorted array.
//    On output, X has been sorted.
//
{
  int i;
  int indx;
  int isgn;
  int j;
  double normsq_i;
  double normsq_j;
  double temp;

  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;

  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );

    if ( 0 < indx )
    {
      temp     = x[0+(i-1)*2];
      x[0+(i-1)*2] = x[0+(j-1)*2];
      x[0+(j-1)*2] = temp;

      temp     = x[1+(i-1)*2];
      x[1+(i-1)*2] = x[1+(j-1)*2];
      x[1+(j-1)*2] = temp;
    }
    else if ( indx < 0 )
    {
      normsq_i = x[0+(i-1)*2] * x[0+(i-1)*2]
               + x[1+(i-1)*2] * x[1+(i-1)*2];

      normsq_j = x[0+(j-1)*2] * x[0+(j-1)*2]
               + x[1+(j-1)*2] * x[1+(j-1)*2];

      if ( normsq_i < normsq_j )
      {
        isgn = -1;
      }
      else
      {
        isgn = +1;
      }
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

double *c8vec_unity ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_UNITY returns the N roots of unity as a double complex vector.
//
//  Discussion:
//
//    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
//
//    X(1:N)^N = ( (1,0), (1,0), ..., (1,0) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, double C8VEC_UNITY[2*N], the N complex roots of unity.
//
{
  double *a;
  int i;
  double const pi = 3.141592653589793;
  double theta;

  a = new double[2*n];

  for ( i = 0; i < n; i++ )
  {
    theta = pi * ( double ) ( 2 * i ) / ( double ) ( n );
    a[0+i*2] = cos ( theta );
    a[1+i*2] = sin ( theta );
  }

  return a;
}
//****************************************************************************80

void daxpy ( int n, double sa, double x[], int incx, double y[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY adds a constant times one vector to another.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    Original FORTRAN77 version by Lawson, Hanson, Kincaid, Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double SA, the multiplier.
//
//    Input, double X[*], the vector to be scaled and added to Y.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], the vector to which a multiple of X is to
//    be added.
//
//    Input, int INCY, the increment between successive entries of Y.
//
{
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
  }
  else if ( sa == 0.0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] + sa * x[i];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy  )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      y[iy] = y[iy] + sa * x[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }

  return;
}
//****************************************************************************80

int file_delete ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_DELETE deletes a named file if it exists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file.
//
//    Output, int FILE_DELETE, is 0 if the file deletion was successful.
{
  int value;
//
//  Does the file exist?
//
  if ( !file_exist ( filename ) )
  {
    return 1;
  }
//
//  Try to remove it.
//
  value = remove ( filename.c_str ( ) );

  if ( value != 0 )
  {
    cerr << "\n";
    cerr << "FILE_DELETE: Warning!\n";
    cerr << "  Could not delete \"" << filename << "\".\n";
    return value;
  }

  cout << "\n";
  cout << "FILE_DELETE:\n";
  cout << "  Deleting old version of \"" << filename << "\".\n";

  return value;
}
//****************************************************************************80

bool file_exist ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_EXIST reports whether a file exists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file.
//
//    Output, bool FILE_EXIST, is TRUE if the file exists.
//
{
  ifstream file;

  file.open ( file_name.c_str ( ), ios::in );

  if ( !file )
  {
    return false;
  }
  else
  {
    return true;
  }
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
//    15 September 2003
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

  if ( ihour > 12 )
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

double *hilbert_inverse ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HILBERT_INVERSE returns the inverse of the Hilbert matrix.
//
//  Formula:
//
//    A(I,J) =  (-1)**(I+J) * (N+I-1)! * (N+J-1)! /
//           [ (I+J-1) * ((I-1)!*(J-1)!)**2 * (N-I)! * (N-J)! ]
//
//  Example:
//
//    N = 5
//
//       25    -300     1050    -1400     630
//     -300    4800   -18900    26880  -12600
//     1050  -18900    79380  -117600   56700
//    -1400   26880  -117600   179200  -88200
//      630  -12600    56700   -88200   44100
//
//  Properties:
//
//    A is symmetric.
//
//    Because A is symmetric, it is normal, so diagonalizable.
//
//    A is almost impossible to compute accurately by general routines
//    that compute the inverse.
//
//    A is integral.
//
//    The sum of the entries of A is N**2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Output, double HILBERT_INVERSE[N*N], the inverse Hilbert matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];
//
//  Set the (1,1) entry.
//
  a[0+0*n] = ( double ) ( n * n );
//
//  Define Row 1, Column J by recursion on Row 1 Column J-1
//
  i = 1;

  for ( j = 2; j <= n; j++ )
  {
    a[i-1+(j-1)*n] = -a[i-1+(j-2)*n] 
      * ( double ) ( ( n + j - 1 ) * ( i + j - 2 ) * ( n + 1 - j ) ) 
      / ( double ) ( ( i + j - 1 ) * ( j - 1 ) * ( j - 1 ) );
  }
//
//  Define Row I by recursion on row I-1
//
  for ( i = 2; i <= n; i++ )
  {
    for (  j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*n] = -a[i-2+(j-1)*n] 
        * ( double ) ( ( n + i - 1 ) * ( i + j - 2 ) * ( n + 1 - i ) ) 
        / ( double ) ( ( i + j - 1 ) * ( i - 1 ) * ( i - 1 ) );
    }
  }

  return a;
}
//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" integer value, usually the largest legal signed int.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" integer.
//
{
  return 2147483647;
}
//****************************************************************************80

int i4_log_10 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_10 returns the integer part of the logarithm base 10 of ABS(X).
//
//  Example:
//
//        I  I4_LOG_10
//    -----  --------
//        0    0
//        1    0
//        2    0
//        9    0
//       10    1
//       11    1
//       99    1
//      100    2
//      101    2
//      999    2
//     1000    3
//     1001    3
//     9999    3
//    10000    4
//
//  Discussion:
//
//    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
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
//    Input, int I, the number whose logarithm base 10 is desired.
//
//    Output, int I4_LOG_10, the integer part of the logarithm base 10 of
//    the absolute value of X.
//
{
  int i_abs;
  int ten_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    ten_pow = 10;

    i_abs = abs ( i );

    while ( ten_pow <= i_abs )
    {
      value = value + 1;
      ten_pow = ten_pow * 10;
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
//    I4_MAX returns the maximum of two integers.
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
//    I4_MIN returns the smaller of two integers.
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

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of integer division.
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
//  Example:
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
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
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
    *seed = *seed + 2147483647;
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

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an integer vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  for ( i = 0; i <= n-1; i++ ) 
  {
    cout << setw(6) << i + 1 << "  " 
         << setw(6) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

int i4vec_search_binary_a ( int n, int a[], int b )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SEARCH_BINARY_A searches an ascending sorted vector for a value.
//
//  Discussion:
//
//    Binary search is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Kreher, Douglas Simpson,
//    Algorithm 1.9,
//    Combinatorial Algorithms,
//    CRC Press, 1998, page 26.
//
//  Parameters:
//
//    Input, int N, the number of elements in the vector.
//
//    Input, int A[N], the array to be searched.  A must
//    be sorted in ascending order.
//
//    Input, int B, the value to be searched for.
//
//    Output, int I4VEC_SEARCH_BINARY_A, the result of the search.
//    -1, B does not occur in A.
//    I, A[I] = B.
//
{
  int high;
  int index;
  int low;
  int mid;
//
//  Check.
//
  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_SEARCH_BINARY_A - Fatal error!\n";
    cerr << "  The array dimension N is less than 1.\n";
    exit ( 1 );
  }

  index = -1;

  low = 1;
  high = n;

  while ( low <= high )
  {
    mid = ( low + high ) / 2;

    if ( a[mid-1] == b )
    {
      index = mid;
      break;
    }
    else if ( a[mid-1] < b )
    {
      low = mid + 1;
    }
    else if ( b < a[mid-1] )
    {
      high = mid - 1;
    }
  }
  return index;
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
//  Examples:
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

float r4_uniform ( float a, float b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM returns a scaled pseudorandom R4.
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
//    21 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, float R4_UNIFORM, a number strictly between A and B.
//
{
  int k;
  float value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  value = ( double ) ( *seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

float r4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01 returns a unit pseudorandom R4.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r4_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R4_UNIFORM_01
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
//    16 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  float r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( float ) ( *seed ) * 4.656612875E-10;

  return r;
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
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

bool r8_is_int ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    R8_IS_INT determines if a real number represents an integer value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the number to be checked.
//
//    Output, bool R8_IS_INT, is TRUE if R is an integer value.
//
{
  if ( ( ( double ) i4_huge ( ) ) < r )
  {
    return false;
  }
  else if ( r < -( ( double ) i4_huge ( ) ) )
  {
    return false;
  }
  else if ( r == ( ( double ) ( ( int ) r ) ) )
  {
    return true;
  }
  else
  {
    return false;
  }
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8s.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8s.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  } 
  else
  {
    return x;
  }
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

double r8_sign2 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN2 returns the first argument with the sign of the second.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the input arguments.
//
//    Output, double R8_SIGN2, is equal to the absolute value of X, and
//    has the sign of Y.
//
{
  if ( 0.0 <= y )
  {
    return x;
  } 
  else
  {
    return (-x);
  }
}
//****************************************************************************80

void r8_swap ( double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP switches two R8s.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  double z;

  z = *x;
  *x = *y;
  *y = z;
 
  return;
}
//****************************************************************************80

double r8_uniform ( double rlo, double rhi, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM randomizes a real number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RLO, RHI, the minimum and maximum values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_UNIFORM, the randomly chosen value.
//
{
  double t;

  t = r8_uniform_01 ( seed );

  t = ( 1.0 - t ) * rlo + t * rhi;

  return t;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
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
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double *r83_cr_fa ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_CR_FA decomposes a real tridiagonal matrix using cyclic reduction.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    Once R83_CR_FA has decomposed a matrix A, then R83_CR_SL may be used to solve
//    linear systems A * x = b.
//
//    R83_CR_FA does not employ pivoting.  Hence, the results can be more
//    sensitive to ill-conditioning than standard Gauss elimination.  In
//    particular, R83_CR_FA will fail if any diagonal element of the matrix
//    is zero.  Other matrices may also cause R83_CR_FA to fail.
//
//    R83_CR_FA can be guaranteed to work properly if the matrix is strictly
//    diagonally dominant, that is, if the absolute value of the diagonal
//    element is strictly greater than the sum of the absolute values of
//    the offdiagonal elements, for each equation.
//
//    The algorithm may be illustrated by the following figures:
//
//    The initial matrix is given by:
//
//          D1 U1
//          L1 D2 U2
//             L2 D3 U3
//                L3 D4 U4
//                   L4 D U5
//                      L5 D6
//
//    Rows and columns are permuted in an odd/even way to yield:
//
//          D1       U1
//             D3    L2 U3
//                D5    L4 U5
//          L1 U2    D2
//             L3 U4    D4
//                L5       D6
//
//    A block LU decomposition is performed to yield:
//
//          D1      |U1
//             D3   |L2 U3
//                D5|   L4 U5
//          --------+--------
//                  |D2'F3
//                  |F1 D4'F4
//                  |   F2 D6'
//
//    For large systems, this reduction is repeated on the lower right hand
//    tridiagonal subsystem until a completely upper triangular system
//    is obtained.  The system has now been factored into the product of a
//    lower triangular system and an upper triangular one, and the information
//    defining this factorization may be used by R83_CR_SL to solve linear
//    systems.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Roger Hockney,
//    A fast direct solution of Poisson's equation using Fourier Analysis,
//    Journal of the ACM,
//    Volume 12, Number 1, pages 95-113, January 1965.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Output, double R83_CR_FA[3*(2*N+1)], factorization information 
//    needed by R83_CR_SL.
//
{
  double *a_cr;
  int iful;
  int ifulp;
  int ihaf;
  int il;
  int ilp;
  int inc;
  int incr;
  int ipnt;
  int ipntp;
  int j;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "R83_CR_FA - Fatal error!\n";
    cerr << "  Nonpositive N = " << n << "\n";
    exit ( 1 );
  }

  a_cr = new double[3*(2*n+1)];

  if ( n == 1 )
  {
    a_cr[0+0*3] = 0.0;
    a_cr[0+1*3] = 0.0;
    a_cr[0+2*3] = 0.0;
    a_cr[1+0*3] = 0.0;
    a_cr[1+1*3] = 1.0 / a[1+0*3];
    a_cr[1+2*3] = 0.0;
    a_cr[2+0*3] = 0.0;
    a_cr[2+1*3] = 0.0;
    a_cr[2+2*3] = 0.0;

    return a_cr;
  }
//
//  Zero out the workspace entries.
//
  a_cr[0+0*3] = 0.0;
  for ( j = 1; j <= n-1; j++ )
  {
    a_cr[0+j*3] = a[0+j*3];
  }
  for ( j = n; j <= 2*n; j++ )
  {
    a_cr[0+j*3] = 0.0;
  }

  a_cr[1+0*3] = 0.0;
  for ( j = 1; j <= n; j++ )
  {
    a_cr[1+j*3] = a[1+(j-1)*3];
  }
  for ( j = n+1; j <= 2*n; j++ )
  {
    a_cr[1+j*3] = 0.0;
  }
  a_cr[2+0*3] = 0.0;
  for ( j = 1; j <= n-1; j++ )
  {
    a_cr[2+j*3] = a[2+(j-1)*3];
  }
  for ( j = n; j <= 2*n; j++ )
  {
    a_cr[2+j*3] = 0.0;
  }

  il = n;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    if ( ( il % 2 ) == 1 )
    {
      inc = il + 1;
    }
    else
    {
      inc = il;
    }

    incr = inc / 2;
    il = il / 2;
    ihaf = ipntp + incr + 1;
    ifulp = ipnt + inc + 2;

    for ( ilp = incr; 1 <= ilp; ilp-- )
    {
      ifulp = ifulp - 2;
      iful = ifulp - 1;
      ihaf = ihaf - 1;

      a_cr[1+iful*3] = 1.0 / a_cr[1+iful*3];
      a_cr[2+iful*3]  = a_cr[2+iful*3]  * a_cr[1+iful*3];
      a_cr[0+ifulp*3] = a_cr[0+ifulp*3] * a_cr[1+(ifulp+1)*3];
      a_cr[1+ihaf*3]  = a_cr[1+ifulp*3] 
        - a_cr[0+iful*3]  * a_cr[2+iful*3]
        - a_cr[0+ifulp*3] * a_cr[2+ifulp*3];
      a_cr[2+ihaf*3] = -a_cr[2+ifulp*3] * a_cr[2+(ifulp+1)*3];
      a_cr[0+ihaf*3] = -a_cr[0+ifulp*3] * a_cr[0+(ifulp+1)*3];
    }
  }

  a_cr[1+(ipntp+1)*3] = 1.0 / a_cr[1+(ipntp+1)*3];

  return a_cr;
}
//****************************************************************************80

double *r83_cr_sl ( int n, double a_cr[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_CR_SL solves a real linear system factored by R83_CR_FA.
//
//  Discussion:
//
//    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
//    LU factors of A.  It does so using a form of cyclic reduction.  If
//    the factors computed by R83_CR_FA are passed to R83_CR_SL, then one or 
//    many linear systems involving the matrix A may be solved.
//
//    Note that R83_CR_FA does not perform pivoting, and so the solution 
//    produced by R83_CR_SL may be less accurate than a solution produced 
//    by a standard Gauss algorithm.  However, such problems can be 
//    guaranteed not to occur if the matrix A is strictly diagonally 
//    dominant, that is, if the absolute value of the diagonal coefficient 
//    is greater than the sum of the absolute values of the two off diagonal 
//    coefficients, for each row of the matrix.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 January 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Roger Hockney,
//    A fast direct solution of Poisson's equation using Fourier Analysis,
//    Journal of the ACM,
//    Volume 12, Number 1, pages 95-113, January 1965.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_CR[3*(2*N+1)], factorization information computed by R83_CR_FA.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R83_CR_SL[N], the solution.
//
{
  int i;
  int iful;
  int ifulm;
  int ihaf;
  int il;
  int ipnt;
  int ipntp;
  int ndiv;
  double *rhs;
  double *x;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "R83_CR_SL - Fatal error!\n";
    cerr << "  Nonpositive N = " << n << "\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x = new double[1];
    x[0] = a_cr[1+1*3] * b[0];
    return x;
  }
//
//  Set up RHS.
//
  rhs = new double[2*n+1];

  rhs[0] = 0.0;
  for ( i = 1; i <= n; i++ )
  {
    rhs[i] = b[i-1];
  }
  for ( i = n+1; i <= 2*n; i++ )
  {
    rhs[i] = 0.0;
  }

  il = n;
  ndiv = 1;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    il = il / 2;
    ndiv = ndiv * 2;
    ihaf = ipntp;

    for ( iful = ipnt + 2; iful <= ipntp; iful = iful + 2 )
    {
      ihaf = ihaf + 1;
      rhs[ihaf] = rhs[iful] 
        - a_cr[2+(iful-1)*3] * rhs[iful-1]
        - a_cr[0+iful*3]     * rhs[iful+1];
    }
  }

  rhs[ihaf] = rhs[ihaf] * a_cr[1+ihaf*3];
  ipnt = ipntp;

  while ( 0 < ipnt )
  {
    ipntp = ipnt;
    ndiv = ndiv / 2;
    il = n / ndiv;
    ipnt = ipnt - il;
    ihaf = ipntp;

    for ( ifulm = ipnt + 1; ifulm <= ipntp; ifulm = ifulm + 2 )
    {
      iful = ifulm + 1;
      ihaf = ihaf + 1;
      rhs[iful] = rhs[ihaf];
      rhs[ifulm] = a_cr[1+ifulm*3] * ( 
                                rhs[ifulm] 
        - a_cr[2+(ifulm-1)*3] * rhs[ifulm-1] 
        - a_cr[0+ifulm*3]     * rhs[iful] );
    }
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = rhs[i+1];
  }

  delete [] rhs;

  return x;
}
//****************************************************************************80

double *r83_cr_sls ( int n, double a_cr[], int nb, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_CR_SLS solves several real linear systems factored by R83_CR_FA.
//
//  Discussion:
//
//    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
//    LU factors of A.  It does so using a form of cyclic reduction.  If
//    the factors computed by R83_CR_FA are passed to R83_CR_SLS, then one or 
//    many linear systems involving the matrix A may be solved.
//
//    Note that R83_CR_FA does not perform pivoting, and so the solutions
//    produced by R83_CR_SLS may be less accurate than a solution produced 
//    by a standard Gauss algorithm.  However, such problems can be 
//    guaranteed not to occur if the matrix A is strictly diagonally 
//    dominant, that is, if the absolute value of the diagonal coefficient 
//    is greater than the sum of the absolute values of the two off diagonal 
//    coefficients, for each row of the matrix.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Roger Hockney,
//    A fast direct solution of Poisson's equation using Fourier Analysis,
//    Journal of the ACM,
//    Volume 12, Number 1, pages 95-113, January 1965.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_CR[3*(2*N+1)], factorization information computed by R83_CR_FA.
//
//    Input, int NB, the number of systems.
//
//    Input, double B[N*NB], the right hand sides.
//
//    Output, double R83_CR_SL[N*NB], the solutions.
//
{
  int i;
  int iful;
  int ifulm;
  int ihaf;
  int il;
  int ipnt;
  int ipntp;
  int j;
  int ndiv;
  double *rhs;
  double *x;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "R83_CR_SLS - Fatal error!\n";
    cerr << "  Nonpositive N = " << n << "\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x = new double[1*nb];
    for ( j = 0; j < nb; j++ )
    {
      x[0+j*n] = a_cr[1+1*3] * b[0+j*n];
    }
    return x;
  }
//
//  Set up RHS.
//
  rhs = new double[(2*n+1)*nb];

  for ( j = 0; j < nb; j++ )
  {
    rhs[0+j*n] = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      rhs[i+j*n] = b[i-1+j*n];
    }
    for ( i = n + 1; i <= 2 * n; i++ )
    {
      rhs[i+j*n] = 0.0;
    }
  }

  il = n;
  ndiv = 1;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    il = il / 2;
    ndiv = ndiv * 2;

    for ( j = 0; j < nb; j++ )
    {
      ihaf = ipntp;
      for ( iful = ipnt + 2; iful <= ipntp; iful = iful + 2 )
      {
        ihaf = ihaf + 1;
        rhs[ihaf+j*n] = rhs[iful+j*n] 
          - a_cr[2+(iful-1)*3] * rhs[iful-1+j*n]
          - a_cr[0+iful*3]     * rhs[iful+1+j*n];
      }
    }
  }

  for ( j = 0; j < nb; j++ )
  {
    rhs[ihaf+j*n] = rhs[ihaf+j*n] * a_cr[1+ihaf*3];
  }

  ipnt = ipntp;

  while ( 0 < ipnt )
  {
    ipntp = ipnt;
    ndiv = ndiv / 2;
    il = n / ndiv;
    ipnt = ipnt - il;

    for ( j = 0; j < nb; j++ )
    {
      ihaf = ipntp;
      for ( ifulm = ipnt + 1; ifulm <= ipntp; ifulm = ifulm + 2 )
      {
        iful = ifulm + 1;
        ihaf = ihaf + 1;
        rhs[iful+j*n] = rhs[ihaf+j*n];
        rhs[ifulm+j*n] = a_cr[1+ifulm*3] * ( 
                                  rhs[ifulm+j*n] 
          - a_cr[2+(ifulm-1)*3] * rhs[ifulm-1+j*n] 
          - a_cr[0+ifulm*3]     * rhs[iful+j*n] );
      }
    }
  }

  x = new double[n*nb];

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = rhs[i+1+j*n];
    }
  }

  delete [] rhs;

  return x;
}
//****************************************************************************80

void r83_gs_sl ( int n, double a[], double b[], double x[], int it_max, 
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    R83_GS_SL solves a R83 system using Gauss-Seidel iteration.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    This routine simply applies a given number of steps of the
//    iteration to an input approximate solution.  On first call, you can
//    simply pass in the zero vector as an approximate solution.  If
//    the returned value is not acceptable, you may call again, using
//    it as the starting point for additional iterations.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Input/output, double X[N], an approximate solution to the system.
//
//    Input, int IT_MAX, the maximum number of iterations to take.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
{
  int i;
  int it_num;
//
//  No diagonal matrix entry can be zero.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      cerr << "\n";
      cerr << "R83_GS_SL - Fatal error!\n";
      cerr << "  Zero diagonal entry, index = " << i << "\n";
      exit ( 1 );
    }
  }

  if ( job == 0 )
  {
    for ( it_num = 1; it_num <= it_max; it_num++ )
    {
      x[0] =   ( b[0]                   - a[2+0*3] * x[1]     ) / a[1+0*3];
      for ( i = 1; i < n-1; i++ )
      {
        x[i] = ( b[i] - a[0+i*3] * x[i-1] - a[2+i*3] * x[i+1] ) / a[1+i*3];
      }
      x[n-1] =   ( b[n-1] - a[0+(n-1)*3] * x[n-2]             ) / a[1+(n-1)*3];
    }
  }
  else
  {
    for ( it_num = 1; it_num <= it_max; it_num++ )
    {
      x[0] =   ( b[0]                     - a[0+1*3] * x[1]     ) 
           / a[1+0*3];
      for ( i = 1; i < n-1; i++ )
      {
        x[i] = ( b[i] - a[2+(i-1)*3] * x[i-1] - a[0+(i+1)*3] * x[i+1] ) 
             / a[1+i*3];
      }
      x[n-1] =   ( b[n-1] - a[2+(n-2)*3] * x[n-2]                     ) 
             / a[1+(n-1)*3];
    }
  }
  return;
}
//****************************************************************************80

double *r83_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R83_INDICATOR sets up a R83 indicator matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//    Here are the values as stored in an indicator matrix:
//
//      00 12 23 34 45
//      11 22 33 44 55
//      21 32 43 54 00
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
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Output, double R83_INDICATOR[3*N], the R83 indicator matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[3*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  a[0+0*3] = 0.0;
  for ( j = 2; j <= n; j++ )
  {
    i = j - 1;
    a[0+(j-1)*3] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n; j++ )
  {
    i = j;
    a[1+(j-1)*3] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n-1; j++ )
  {
    i = j + 1;
    a[2+(j-1)*3] = ( double ) ( fac * i + j );
  }
  a[2+(n-1)*3] = 0.0;

  return a;
}
//****************************************************************************80

void r83_jac_sl ( int n, double a[], double b[], double x[], int it_max, 
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    R83_JAC_SL solves a R83 system using Jacobi iteration.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    This routine simply applies a given number of steps of the
//    iteration to an input approximate solution.  On first call, you can
//    simply pass in the zero vector as an approximate solution.  If
//    the returned value is not acceptable, you may call again, using
//    it as the starting point for additional iterations.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Input/output, double X[N], an approximate solution to the system.
//
//    Input, int IT_MAX, the maximum number of iterations to take.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
{
  int i;
  int it_num;
  double *xnew;

  xnew = new double[n];
//
//  No diagonal matrix entry can be zero.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      cerr << "\n";
      cerr << "R83_JAC_SL - Fatal error!\n";
      cerr << "  Zero diagonal entry, index = " << i << "\n";
      exit ( 1 );
    }
  }

  for ( it_num = 1; it_num <= it_max; it_num++ )
  {
//
//  Solve A*x=b:
//
    if ( job == 0 )
    {
      xnew[0] =   b[0]                           - a[2+0*3] * x[1];
      for ( i = 1; i < n-1; i++ )
      {
        xnew[i] = b[i]   - a[0+i*3]     * x[i-1] - a[2+i*3] * x[i+1];
      }
      xnew[n-1] = b[n-1] - a[0+(n-1)*3] * x[n-2];
    }
//
//  Solve A'*x=b:
//
    else
    {
      xnew[0] =   b[0]                     - a[0+1*3] * x[1];
      for ( i = 1; i < n-1; i++ )
      {
        xnew[i] = b[i] - a[2+(i-1)*3] * x[i-1] - a[0+(i+1)*3] * x[i+1];
      }
      xnew[n-1] =   b[n-1] - a[2+(n-2)*3] * x[n-2];
    }
//
//  Divide by the diagonal term, and overwrite X.
//
    for ( i = 0; i < n; i++ )
    {
      xnew[i] = xnew[i] / a[1+i*3];
    }

    for ( i = 0; i < n; i++ )
    {
      x[i] = xnew[i];
    }
  }

  delete [] xnew;

  return;
}
//****************************************************************************80

double *r83_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_MXV multiplies a R83 matrix times a vector.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R83_MXV[N], the product A * x.
//
{
  double *b;
  int i;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] =        a[1+i*3] * x[i];
  }
  for ( i = 0; i < n-1; i++ )
  {
    b[i] = b[i] + a[0+(i+1)*3] * x[i+1];
  }
  for ( i = 1; i < n; i++ )
  {
    b[i] = b[i] + a[2+(i-1)*3] * x[i-1];
  }

  return b;
}
//****************************************************************************80

double r83_np_det ( int n, double a_lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_NP_DET: determinant of a tridiagonal system factored by R83_NP_FA.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Input, double A_LU[3*N], the LU factors from R83_NP_FA.
//
//    Output, double R83_NP_DET, the determinant of the matrix.
//
{
  double det;
  int j;

  det = 1.0;
  for ( j = 0; j < n; j++ )
  {
    det = det * a_lu[1+j*3];
  }

  return det;
}
//****************************************************************************80

int r83_np_fa ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_NP_FA factors a R83 system without pivoting.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    Because this routine does not use pivoting, it can fail even when
//    the matrix is not singular, and it is liable to make larger
//    errors.
//
//    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
//    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
//    in one step, and does not save the factorization.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Input/output, double A[3*N].
//    On input, the tridiagonal matrix.  On output, factorization information.
//
//    Output, int R83_NP_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int i;

  for ( i = 1; i <= n-1; i++ )
  {
    if ( a[1+(i-1)*3] == 0.0 )
    {
      cerr << "\n";
      cerr << "R83_NP_FA - Fatal error!\n";
      cerr << "  Zero pivot on step " << i << "\n";
      exit ( 1 );
    }
//
//  Store the multiplier in L.
//
    a[2+(i-1)*3] = a[2+(i-1)*3] / a[1+(i-1)*3];
//
//  Modify the diagonal entry in the next column.
//
    a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3];
  }

  if ( a[1+(n-1)*3] == 0.0 )
  {
    cerr << "\n";
    cerr << "R83_NP_FA - Fatal error!\n";
    cerr << "  Zero pivot on step " << n << "\n";
    exit ( 1 );
  }

  return 0;
}
//****************************************************************************80

double *r83_np_fs ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_NP_FS factors and solves a R83 system.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    This algorithm requires that each diagonal entry be nonzero.
//    It does not use pivoting, and so can fail on systems that
//    are actually nonsingular.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, double A[3*N].
//    On input, the nonzero diagonals of the linear system.
//    On output, the data in these vectors has been overwritten
//    by factorization information.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R83_NP_FS[N], the solution of the linear system.
//    This is NULL if there was an error because one of the diagonal
//    entries was zero.
//
{
  int i;
  double *x;
  double xmult;
//
//  Check.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      return NULL;
    }
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( i = 1; i < n; i++ )
  {
    xmult = a[2+(i-1)*3] / a[1+(i-1)*3];
    a[1+i*3] = a[1+i*3] - xmult * a[0+i*3];
    x[i] = x[i] - xmult * x[i-1];
  }

  x[n-1] = x[n-1] / a[1+(n-1)*3];
  for ( i = n-2; 0 <= i; i-- )
  {
    x[i] = ( x[i] - a[0+(i+1)*3] * x[i+1] ) / a[1+i*3];
  }

  return x;
}
//****************************************************************************80

double *r83_np_ml ( int n, double a_lu[], double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R83_NP_ML computes Ax or xA, where A has been factored by R83_NP_FA.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Input, double A_LU[3*N], the LU factors from R83_FA.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double B[N], the product.
//
//    Input, int JOB, specifies the product to find.
//    0, compute A * x.
//    nonzero, compute A' * x.
//
{
  double *b;
  int i;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
//
//  Compute X := U * X
//
    for ( i = 1; i <= n; i++ )
    {
      b[i-1] = a_lu[1+(i-1)*3] * b[i-1];

      if ( i < n )
      {
        b[i-1] = b[i-1] + a_lu[0+i*3] * b[i];
      }
    }
//
//  Compute X: = L * X.
//
    for ( i = n; 2 <= i; i-- )
    {
      b[i-1] = b[i-1] + a_lu[2+(i-2)*3] * b[i-2];
    }
  }
  else
  {
//
//  Compute X: = L' * X.
//
    for ( i = 1; i <= n-1; i++ )
    {
      b[i-1] = b[i-1] + a_lu[2+(i-1)*3] * b[i];
    }
//
//  Compute X: = U' * X.
//
    for ( i = n; 1 <= i; i-- )
    {
      b[i-1] = a_lu[1+(i-1)*3] * b[i-1];
      if ( 1 < i )
      {
        b[i-1] = b[i-1] + a_lu[0+(i-1)*3] * b[i-2];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r83_np_sl ( int n, double a_lu[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R83_NP_SL solves a R83 system factored by R83_NP_FA.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Input, double A_LU[3*N], the LU factors from R83_NP_FA.
//
//    Input, double B[N], the right hand side of the linear system.
//    On output, B contains the solution of the linear system.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double R83_NP_SL[N], the solution of the linear system.
//
{
  int i;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
//
//  Solve L * Y = B.
//
    for ( i = 1; i < n; i++ )
    {
      x[i] = x[i] - a_lu[2+(i-1)*3] * x[i-1];
    }
//
//  Solve U * X = Y.
//
    for ( i = n; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( 1 < i )
      {
        x[i-2] = x[i-2] - a_lu[0+(i-1)*3] * x[i-1];
      }
    }
  }
  else
  {
//
//  Solve U' * Y = B
//
    for ( i = 1; i <= n; i++ )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( i < n )
      {
        x[i] = x[i] - a_lu[0+i*3] * x[i-1];
      }
    }
//
//  Solve L' * X = Y.
//
    for ( i = n-1; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] - a_lu[2+(i-1)*3] * x[i];
    }
  }

  return x;
}
//****************************************************************************80

void r83_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83_PRINT prints a R83 matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Input, string TITLE, a title.
//
{
  r83_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r83_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83_PRINT_SOME prints some of a R83 matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column, to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - 1 );

    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + 1 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(6) << i << ": ";

      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( 1 < i - j || 1 < j - i )
        {
          cout << "              ";
        }
        else if ( j == i + 1 )
        {
          cout << setw(12) << a[0+(j-1)*3] << "  ";
        }
        else if ( j == i )
        {
          cout << setw(12) << a[1+(j-1)*3] << "  ";
        }
        else if ( j == i - 1 )
        {
          cout << setw(12) << a[2+(j-1)*3] << "  ";
        }

      }
      cout << "\n";
    }
  }
  return;
# undef INCX
}
//****************************************************************************80

double *r83_random ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R83_RANDOM randomizes a R83 matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R83_RANDOM[3*N], the R83 matrix.
//
{
  double *a;
  int i;
  double *u;
  double *v;
  double *w;

  a = new double[3*n];

  u = r8vec_random ( n-1, 0.0, 1.0, seed );
  v = r8vec_random ( n,   0.0, 1.0, seed );
  w = r8vec_random ( n-1, 0.0, 1.0, seed );

  a[0+0*3] = 0.0;
  for ( i = 1; i < n; i++ )
  {
    a[0+i*3] = u[i-1];
  }
   for ( i = 0; i < n; i++ )
  {
    a[1+i*3] = v[i];
  }
  for ( i = 0; i < n-1; i++ )
  {
    a[2+i*3] = w[i];
  }
  a[2+(n-1)*3] = 0.0;

  delete [] u;
  delete [] v;
  delete [] w;

  return a;
}
//****************************************************************************80

double *r83_to_r8ge ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_TO_R8GE copies a R83 matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Output, double R83_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( j == i-1 )
      {
        b[i-1+(j-1)*n] = a[0+(i-1)*3];
      }
      else if ( i == j )
      {
        b[i-1+(j-1)*n] = a[1+(i-1)*3];
      }
      else if ( j == i+1 )
      {
        b[i-1+(j-1)*n] = a[2+(i-1)*3];
      }
      else
      {
        b[i-1+(j-1)*n] = 0.0;
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r83_vxm ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_VXM multiplies a vector times a R83 matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Input, double X[N], the vector to be multiplied by A'.
//
//    Output, double R83_VXM[N], the product A' * x.
//
{
  double *b;
  int i;

  b = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    b[i-1] = a[1+(i-1)*3] * x[i-1];
  }

  for ( i = 1; i <= n-1; i++ )
  {
    b[i-1] = b[i-1] + a[2+(i-1)*3] * x[i];
  }

  for ( i = 2; i <= n; i++ )
  {
    b[i-1] = b[i-1] + a[0+(i-1)*3] * x[i-2];
  }

  return b;
}
//****************************************************************************80

double *r83_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R83_ZERO zeros a R83 matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Output, double R83_ZERO[3*N], the R83 matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      a[i+j*3] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r83np_fs ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83NP_FS factors and solves an R83NP system.
//
//  Discussion:
//
//    The R83NP storage format is used for a tridiagonal matrix.
//    The subdiagonal   is in entries (0,1:N-1), 
//    the diagonal      is in entries (1,0:N-1), 
//    the superdiagonal is in entries (2,0:N-2). 
//
//    This algorithm requires that each diagonal entry be nonzero.
//    It does not use pivoting, and so can fail on systems that
//    are actually nonsingular.
//
//    The "R83NP" format used for this routine is different from the R83 format.
//    Here, we insist that the nonzero entries
//    for a given row now appear in the corresponding column of the
//    packed array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A21 A32 A43 A54
//      A11 A22 A33 A44 A55
//      A12 A23 A34 A45  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, double A[3*N].
//    On input, the nonzero diagonals of the linear system.
//    On output, the data in these vectors has been overwritten
//    by factorization information.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R83NP_FS[N], the solution of the linear system.
//
{
  int i;
  double *x;
//
//  Check.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      cerr << "\n";
      cerr << "R83NP_FS - Fatal error!\n";
      cerr << "  A[1+" << i << "*3] = 0.\n";
      exit ( 1 );
    }
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( i = 1; i < n; i++ )
  {
    a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3] / a[1+(i-1)*3];
    x[i]     = x[i]     - x[i-1]       * a[0+i*3] / a[1+(i-1)*3];
  }

  x[n-1] = x[n-1] / a[1+(n-1)*3];
  for ( i = n-2; 0 <= i; i-- )
  {
    x[i] = ( x[i] - a[2+i*3] * x[i+1] ) / a[1+i*3];
  }

  return x;
}
//****************************************************************************80

double r83p_det ( int n, double a_lu[], double work4 )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_DET computes the determinant of a matrix factored by R83P_FA.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored 
//    as a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input, double A_LU[3*N], the LU factors from R83P_FA.
//
//    Input, double WORK4, factorization information from R83P_FA.
//
//    Output, double R83P_DET, the determinant of the matrix.
//
{
  double det;
  int i;

  det = work4;
  for ( i = 0; i <= n-2; i++ )
  {
    det = det * a_lu[1+i*3];
  }

  return det;
}
//****************************************************************************80

int r83p_fa ( int n, double a[], double work2[], double work3[], double *work4 )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_FA factors a R83P matrix.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//    Once the matrix has been factored by R83P_FA, R83P_SL may be called
//    to solve linear systems involving the matrix.
//
//    The logical matrix has a form which is suggested by this diagram:
//
//      D1 U1          L1
//      L2 D2 U2
//         L3 R83 U3
//            L4 D4 U4
//               L5 R85 U5
//      U6          L6 D6
//
//    The algorithm treats the matrix as a border banded matrix:
//
//      ( A1  A2 )
//      ( A3  A4 )
//
//    where:
//
//      D1 U1          | L1
//      L2 D2 U2       |  0
//         L3 R83 U3    |  0
//            L4 D4 U4 |  0
//               L5 R85 | U5
//      ---------------+---
//      U6  0  0  0 L6 | D6
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Method:
//
//    The algorithm rewrites the system as:
//
//         X1 + inverse(A1) A2 X2 = inverse(A1) B1
//
//      A3 X1 +             A4 X2 = B2
//
//    The first equation can be "solved" for X1 in terms of X2:
//
//         X1 = - inverse(A1) A2 X2 + inverse(A1) B1
//
//    allowing us to rewrite the second equation for X2 explicitly:
//
//      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input/output, double A[3*N].
//    On input, the periodic tridiagonal matrix.  
//    On output, the arrays have been modified to hold information
//    defining the border-banded factorization of submatrices A1
//    and A3.
//
//    Output, int R83P_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
//    Output, double WORK2[N-1], WORK3[N-1], *WORK4, factorization information.
//
{
  int i;
  int info;
  int job;
  double *work1;

  work1 = new double[n-1];
//
//  Compute inverse(A1):
//
  info = r83_np_fa ( n-1, a );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "R83P_FA - Fatal error!\n";
    cerr << "  R83_NP_FA returned INFO = " << info << "\n";
    cerr << "  Factoring failed for column INFO.\n";
    cerr << "  The tridiagonal matrix A1 is singular.\n";
    cerr << "  This algorithm cannot continue!\n";
    exit ( 1 );
  }
//
//  WORK2 := inverse(A1) * A2.
//
  work2[0] = a[2+(n-1)*3];
  for ( i = 1; i < n-2; i++)
  {
    work2[i] = 0.0;
  }
  work2[n-2] = a[0+(n-1)*3];

  job = 0;
  work1 = r83_np_sl ( n-1, a, work2, job );
  for ( i = 0; i < n-1; i++ )
  {
    work2[i] = work1[i];
  }
//
//  WORK3 := inverse ( A1' ) * A3'.
//
  work3[0] = a[0+0*3];
  for ( i = 1; i < n-2; i++)
  {
    work3[i] = 0.0;
  }
  work3[n-2] = a[2+(n-2)*3];

  job = 1;
  work1 = r83_np_sl ( n-1, a, work3, job );
  for ( i = 0; i < n-1; i++ )
  {
    work3[i] = work1[i];
  }
//
//  A4 := ( A4 - A3 * inverse(A1) * A2 )
//
  *work4 = a[1+(n-1)*3] - a[0+0*3] * work2[0] - a[2+(n-2)*3] * work2[n-2];

  if ( *work4 == 0.0 )
  {
    cerr << "\n";
    cerr << "R83P_FA - Fatal error!\n";
    cerr << "  The factored A4 submatrix is zero.\n";
    cerr << "  This algorithm cannot continue!\n";
    exit ( 1 );
  }

  delete [] work1;

  return 0;
}
//****************************************************************************80

double *r83p_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_INDICATOR sets up a R83P indicator matrix.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored 
//    as a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//    Here are the values as stored in an indicator matrix:
//
//      51 12 23 34 45
//      11 22 33 44 55
//      21 32 43 54 15
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
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Output, double R83P_INDICATOR[3*N], the R83P indicator matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[3*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  i = n;
  j = 1;
  a[0+(j-1)*3] = ( double ) ( fac * i + j );
  for ( j = 2; j <= n; j++ )
  {
    i = j - 1;
    a[0+(j-1)*3] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n; j++ )
  {
    i = j;
    a[1+(j-1)*3] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n-1; j++ )
  {
    i = j + 1;
    a[2+(j-1)*3] = ( double ) ( fac * i + j );
  }
  i = 1;
  j = n;
  a[2+(j-1)*3] = ( double ) ( fac * i + j );

  return a;
}
//****************************************************************************80

double *r83p_ml ( int n, double a_lu[], double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_ML computes A * x or x * A, where A has been factored by R83P_FA.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored 
//    as a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input, double A_LU[3*N], the LU factors from R83P_FA.
//
//    Input, double X[N], the vector to be multiplied by the matrix.
//
//    Input, int JOB, indicates what product should be computed.
//    0, compute A * x.
//    nonzero, compute A' * x.
//
//    Output, double R83P_ML[N], the result of the multiplication.
//
{
  double *b;
  double *b_short;
  int i;
//
//  Multiply A(1:N-1,1:N-1) and X(1:N-1).
//
  b_short = r83_np_ml ( n-1, a_lu, x, job );

  b = new double[n];

  for ( i = 0; i < n-1; i++ )
  {
    b[i] = b_short[i];
  }
  b[n-1] = 0.0;

  delete [] b_short;
//
//  Add terms from the border.
//
  if ( job == 0 )
  {
    b[0] = b[0] + a_lu[2+(n-1)*3] * x[n-1];
    b[n-2] = b[n-2] + a_lu[0+(n-1)*3] * x[n-1];
    b[n-1] = a_lu[0+0*3] * x[0] + a_lu[2+(n-2)*3] * x[n-2] 
      + a_lu[1+(n-1)*3] * x[n-1];
  }
  else
  {
    b[0] = b[0] + a_lu[0+0*3] * x[n-1];
    b[n-2] = b[n-2] + a_lu[2+(n-2)*3] * x[n-1];
    b[n-1] = a_lu[2+(n-1)*3] * x[0] + a_lu[0+(n-1)*3] * x[n-2] 
           + a_lu[1+(n-1)*3] * x[n-1];
  }

  return b;
}
//****************************************************************************80

double *r83p_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_MXV multiplies a R83P matrix times a vector.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input, double A[3*N], the R83P matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R83P_MXV[N], the product A * x.
//
{
  double *b;
  int i;

  b = new double[n];

  b[0] =   a[2+(n-1)*3] * x[n-1] + a[1+0*3]     * x[0]   + a[0+1*3]     * x[1];

  for ( i = 1; i < n-1; i++ )
  {
    b[i] = a[2+(i-1)*3] * x[i-1] + a[1+i*3]     * x[i]   + a[0+(i+1)*3] * x[i+1];
  }

  b[n-1] = a[2+(n-2)*3] * x[n-2] + a[1+(n-1)*3] * x[n-1] + a[0+0*3]     * x[0];

  return b;
}
//****************************************************************************80

void r83p_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_PRINT prints a R83P matrix.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the R83P matrix.
//
//    Input, string TITLE, a title.
//
{
  r83p_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r83p_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_PRINT_SOME prints some of a R83P matrix.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the R83P matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column, to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );

    if ( 1 < i2lo || j2hi < n )
    {
      i2lo = i4_max ( i2lo, j2lo - 1 );
    }

    i2hi = i4_min ( ihi, n );

    if ( i2hi < n || 1 < j2lo )
    {
      i2hi = i4_min ( i2hi, j2hi + 1 );
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(4) << i << "  ";

      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( i == n && j == 1 )
        {
          cout << setw(12) << a[0+(j-1)*3] << "  ";
        }
        else if ( i == 1 && j == n )
        {
          cout << setw(12) << a[2+(j-1)*3] << "  ";
        }
        else if ( 1 < i-j || 1 < j-i )
        {
          cout << "              ";
        }
        else if ( j == i+1 )
        {
          cout << setw(12) << a[0+(j-1)*3] << "  ";
        }
        else if ( j == i )
        {
          cout << setw(12) << a[1+(j-1)*3] << "  ";
        }
        else if ( j == i-1 )
        {
          cout << setw(12) << a[2+(j-1)*3] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r83p_random ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_RANDOM randomizes a R83P matrix.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R83P_RANDOM[3*N], the R83P matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      a[i+j*3] = r8_uniform_01 ( seed );
    }
  }
  return a;
}
//****************************************************************************80

double *r83p_sl ( int n, double a_lu[], double b[], int job, double work2[], 
  double work3[], double work4 )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_SL solves a R83P system factored by R83P_FA.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input, double A_LU[3*N], the LU factors from R83P_FA.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Input, double WORK2(N-1), WORK3(N-1), WORK4, factor data from R83P_FA.
//
//    Output, double R83P_SL[N], the solution to the linear system.
//
{
  int i;
  double *x;
  double *xnm1;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
//
//  Solve A1 * X1 = B1.
//
    xnm1 = r83_np_sl ( n-1, a_lu, x, job );
//
//  X2 = B2 - A3 * X1
//
    for ( i = 0; i < n-1; i++ )
    {
      x[i] = xnm1[i];
    }
    delete [] xnm1;

    x[n-1] = x[n-1] - a_lu[0+0*3] * x[0] - a_lu[2+(n-2)*3] * x[n-2];
//
//  Solve A4 * X2 = X2
//
    x[n-1] = x[n-1] / work4;
//
//  X1 := X1 - inverse ( A1 ) * A2 * X2.
//
    for ( i = 0; i < n-1; i++ )
    {
      x[i] = x[i] - work2[i] * x[n-1];
    }
  }
  else
  {
//
//  Solve A1' * X1 = B1.
//
    xnm1 = r83_np_sl ( n-1, a_lu, x, job );
//
//  X2 := X2 - A2' * B1
//
    for ( i = 0; i < n-1; i++ )
    {
      x[i] = xnm1[i];
    }
    delete [] xnm1;

    x[n-1] = x[n-1] - a_lu[2+(n-1)*3] * x[0] - a_lu[0+(n-1)*3] * x[n-2];
//
//  Solve A4 * X2 = X2.
//
    x[n-1] = x[n-1] / work4;
//
//  X1 := X1 - transpose ( inverse ( A1 ) * A3 ) * X2.
//
    for ( i = 0; i < n-1; i++ )
    {
      x[i] = x[i] - work3[i] * x[n-1];
    }
  }
  return x;
}
//****************************************************************************80

double *r83p_to_r8ge ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_TO_R8GE copies a R83P matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input, double A[3*N], the R83P matrix.
//
//    Output, double R83P_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( i == j )
      {
        b[i-1+(j-1)*n] = a[1+(j-1)*3];
      }
      else if ( j == i-1 )
      {
        b[i-1+(j-1)*n] = a[2+(j-1)*3];
      }
      else if ( j == i+1 )
      {
        b[i-1+(j-1)*n] = a[0+(j-1)*3];
      }
      else if ( i == 1 && j == n )
      {
        b[i-1+(j-1)*n] = a[2+(j-1)*3];
      }
      else if ( i == n && j == 1 )
      {
        b[i-1+(j-1)*n] = a[0+(j-1)*3];
      }
      else
      {
        b[i-1+(j-1)*n] = 0.0;
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r83p_vxm ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_VXM multiplies a vector times a R83P matrix.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input, double A[3*N], the R83P matrix.
//
//    Input, double X, the vector to be multiplied by A.
//
//    Output, double R83P_VXM[N], the product X * A.
//
{
  double *b;
  int i;

  b = new double[n];

  b[0] = a[0+0*3] * x[n-1] + a[1+0*3] * x[0] + a[2+0*3] * x[1];

  for ( i = 2; i <= n-1; i++ )
  {
    b[i-1] = a[0+(i-1)*3] * x[i-2] + a[1+(i-1)*3] * x[i-1] + a[2+(i-1)*3] * x[i];
  }

  b[n-1] = a[0+(n-1)*3] * x[n-2] + a[1+(n-1)*3] * x[n-1] + a[2+(n-1)*3] * x[0];

  return b;
}
//****************************************************************************80

double *r83p_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R83P_ZERO zeros a R83P matrix.
//
//  Discussion:
//
//    The R83P storage format stores a periodic tridiagonal matrix is stored as 
//    a 3 by N array, in which each row corresponds to a diagonal, and 
//    column locations are preserved.  The matrix value 
//    A(1,N) is stored as the array entry A(1,1), and the matrix value
//    A(N,1) is stored as the array entry A(3,N).
//
//  Example:
//
//    Here is how a R83P matrix of order 5 would be stored:
//
//      A51 A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54 A15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Output, double S3P[3*N], the R83P matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      a[i+j*3] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r83s_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83S_MXV multiplies an R83S matrix times a vector.
//
//  Discussion:
//
//    The R83S storage format is used for a tridiagonal scalar matrix.
//    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
//    values that occur on every row.
//
//  Example:
//
//    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
//    be interpreted:
//
//      A2  A3   0   0   0
//      A1  A2  A3   0   0
//       0  A1  A2  A3   0 
//       0   0  A1  A2  A3
//       0   0   0  A1  A2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A[3], the R83S matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R83S_MXV[N], the product A * x.
//
{
  double *b;
  int i;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 1; i < n; i++ )
  {
    b[i] = b[i] + a[0] * x[i-1];
  }

  for ( i = 0; i < n; i++ )
  {
    b[i] = b[i] + a[1] * x[i];
  }

  for ( i = 0; i < n - 1; i++ )
  {
    b[i] = b[i] + a[2] * x[i+1];
  }

  return b;
}
//****************************************************************************80

void r83s_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83S_PRINT prints an R83S matrix.
//
//  Discussion:
//
//    The R83S storage format is used for a tridiagonal scalar matrix.
//    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
//    values that occur on every row.
//
//  Example:
//
//    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
//    be interpreted:
//
//      A2  A3   0   0   0
//      A1  A2  A3   0   0
//       0  A1  A2  A3   0 
//       0   0  A1  A2  A3
//       0   0   0  A1  A2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3], the R83S matrix.
//
//    Input, string TITLE, a title.
//
{
  r83s_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r83s_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83S_PRINT_SOME prints some of a R83S matrix.
//
//  Discussion:
//
//    The R83S storage format is used for a tridiagonal scalar matrix.
//    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
//    values that occur on every row.
//
//  Example:
//
//    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
//    be interpreted:
//
//      A2  A3   0   0   0
//      A1  A2  A3   0   0
//       0  A1  A2  A3   0 
//       0   0  A1  A2  A3
//       0   0   0  A1  A2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3], the R83S matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column, to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - 1 );

    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + 1 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(6) << i << ": ";

      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( 1 < i - j || 1 < j - i )
        {
          cout << "              ";
        }
        else if ( j == i - 1 )
        {
          cout << setw(12) << a[0] << "  ";
        }
        else if ( j == i )
        {
          cout << setw(12) << a[1] << "  ";
        }
        else if ( j == i + 1 )
        {
          cout << setw(12) << a[2] << "  ";
        }

      }
      cout << "\n";
    }
  }
  return;
# undef INCX
}
//****************************************************************************80

double *r83s_random ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R83S_RANDOM randomizes an R83S matrix.
//
//  Discussion:
//
//    The R83S storage format is used for a tridiagonal scalar matrix.
//    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
//    values that occur on every row.
//
//  Example:
//
//    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
//    be interpreted:
//
//      A2  A3   0   0   0
//      A1  A2  A3   0   0
//       0  A1  A2  A3   0 
//       0   0  A1  A2  A3
//       0   0   0  A1  A2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R83_RANDOM[3], the R83S matrix.
//
{
  double *a;

  a = r8vec_uniform_01_new ( 3, seed );

  return a;
}
//****************************************************************************80

double *r85_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R85_INDICATOR sets up a R85 indicator matrix.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//    Here are the values as stored in an indicator matrix:
//
//      00 00 13 24 35 46
//      00 12 23 34 45 56
//      11 22 33 44 55 66
//      21 32 43 54 65 00
//      31 42 53 64 00 00
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Output, double R85_INDICATOR[3*N], the R85 indicator matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[5*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  a[0+0*5] = 0.0;
  a[0+1*5] = 0.0;
  for ( j = 3; j <= n; j++ )
  {
    i = j - 2;
    a[0+(j-1)*5] = ( double ) ( fac * i + j );
  }

  a[1+0*5] = 0.0;
  for ( j = 2; j <= n; j++ )
  {
    i = j - 1;
    a[1+(j-1)*5] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n; j++ )
  {
    i = j;
    a[2+(j-1)*5] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n-1; j++ )
  {
    i = j + 1;
    a[3+(j-1)*5] = ( double ) ( fac * i + j );
  }
  a[3+(n-1)*5] = 0.0;

  for ( j = 1; j <= n-2; j++ )
  {
    i = j + 2;
    a[4+(j-1)*5] = ( double ) ( fac * i + j );
  }
  a[4+(n-2)*5] = 0.0;
  a[4+(n-1)*5] = 0.0;

  return a;
}
//****************************************************************************80

double *r85_np_fs ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R85_NP_FS factors and solves a R85 system.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//    This algorithm requires that each diagonal entry be nonzero.
//
//    No pivoting is performed, and therefore the algorithm may fail
//    in simple cases where the matrix is not singular.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Cheney, Kincaid.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Ward Cheney, David Kincaid,
//    Numerical Mathematics and Computing,
//    1985, pages 233-236.
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, double A[5*N],
//    On input, the pentadiagonal matrix.
//    On output, the data in these vectors has been overwritten
//    by factorization information.
//
//    Input/output, double B[N].
//    On input, B contains the right hand side of the linear system.
//    On output, B has been overwritten by factorization information.
//
//    Output, double R85_NP_FS[N], the solution of the linear system.
//
{
  int i;
  double *x;
  double xmult;

  for ( i = 0; i < n; i++ )
  {
    if ( a[2+i*5] == 0.0 )
    {
      return NULL;
    }
  }

  x = new double[n];

  for ( i = 2; i <= n-1; i++ )
  {
    xmult = a[1+(i-1)*5] / a[2+(i-2)*5];
    a[2+(i-1)*5] = a[2+(i-1)*5] - xmult * a[3+(i-2)*5];
    a[3+(i-1)*5] = a[3+(i-1)*5] - xmult * a[4+(i-2)*5];

    b[i-1] = b[i-1] - xmult * b[i-2];

    xmult = a[0+i*5] / a[2+(i-2)*5];
    a[1+i*5] = a[1+i*5] - xmult * a[3+(i-2)*5];
    a[2+i*5] = a[2+i*5] - xmult * a[4+(i-2)*5];

    b[i] = b[i] - xmult * b[i-2];
  }

  xmult = a[1+(n-1)*5] / a[2+(n-2)*5];
  a[2+(n-1)*5] = a[2+(n-1)*5] - xmult * a[3+(n-2)*5];

  x[n-1] = ( b[n-1] - xmult * b[n-2] ) / a[2+(n-1)*5];
  x[n-2] = ( b[n-2] - a[3+(n-2)*5] * x[n-1] ) / a[2+(n-2)*5];

  for ( i = n - 2; 1 <= i; i-- )
  {
    x[i-1] = ( b[i-1] - a[3+(i-1)*5] * x[i] - a[4+(i-1)*5] * x[i+1] ) 
      / a[2+(i-1)*5];
  }

  return x;
}
//****************************************************************************80

double *r85_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R85_MXV multiplies a R85 matrix times a vector.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A[5*N], the pentadiagonal matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R85_MXV[N], the product A * x.
//
{
  double *b;
  int i;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = a[2+i*5] * x[i];
  }
  for ( i = 2; i < n; i++ )
  {
    b[i] = b[i] + a[0+i*5] * x[i-2];
  }
  for ( i = 1; i < n; i++ )
  {
    b[i] = b[i] + a[1+i*5] * x[i-1];
  }

  for ( i = 0; i < n-1; i++ )
  {
    b[i] = b[i] + a[3+i*5] * x[i+1];
  }
  for ( i = 0; i < n-2; i++ )
  {
    b[i] = b[i] + a[4+i*5] * x[i+2];
  }

  return b;
}
//****************************************************************************80

void r85_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R85_PRINT prints a R85 matrix.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[5*N], the pentadiagonal matrix.
//
//    Input, string TITLE, a title.
//
{
  r85_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r85_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    R85_PRINT_SOME prints some of a R85 matrix.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[5*N], the pentadiagonal matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column, to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col:  ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - 2 );

    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + 2 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(6) << i << "  ";

      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( 2 < i-j || 2 < j-i )
        {
          cout << "            ";
        }
        else if ( j == i+2 )
        {
          cout << setw(10) << a[0+(j-1)*5] << "  ";
        }
        else if ( j == i+1 )
        {
          cout << setw(10) << a[1+(j-1)*5] << "  ";
        }
        else if ( j == i )
        {
          cout << setw(10) << a[2+(j-1)*5] << "  ";
        }
        else if ( j == i-1 )
        {
          cout << setw(10) << a[3+(j-1)*5] << "  ";
        }
        else if ( j == i-2 )
        {
          cout << setw(10) << a[4+(j-1)*5] << "  ";
        }
      }
      cout << "\n";
    }
    cout << "\n";
  }
  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r85_random ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R85_RANDOM randomizes a R85 matrix.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R85_RANDOM[5*N], the pentadiagonal matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[5*n];

  i = 0;
  a[0+0*5]     = 0.0;
  a[0+1*5]     = 0.0;
  for ( j = 2; j < n; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }

  i = 1;
  a[1+0*5]     = 0.0;
  for ( j = 1; j < n; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }

  i = 2;
  for ( j = 0; j < n; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }

  i = 3;
  for ( j = 0; j < n-1; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }
  a[3+(n-1)*5] = 0.0;

  i = 4;
  for ( j = 0; j < n-2; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }
  a[4+(n-2)*5] = 0.0;
  a[4+(n-1)*5] = 0.0;

  return a;
}
//****************************************************************************80

double *r85_to_r8ge ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R85_TO_R8GE copies a R85 matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 3.
//
//    Input, double A[5*N], the nonzero diagonals of the matrix.
//
//    Output, double R85_TO_R8GE[N*N], the pentadiagonal matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( j == i-2 )
      {
        b[i+j*5] = a[0+i*5];
      }
      else if ( j == i-1 )
      {
        b[i+j*5] = a[1+i*5];
      }
      else if ( i == j )
      {
        b[i+j*5] = a[2+i*5];
      }
      else if ( j == i+1 )
      {
        b[i+j*5] = a[3+i*5];
      }
      else if ( j == i+2 )
      {
        b[i+j*5] = a[4+i*5];
      }
      else
      {
        b[i+j*5] = 0.0;
      }
    }
  }
  return b;
}
//****************************************************************************80

double *r85_vxm ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R85_VXM multiplies a vector times a R85 matrix.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A[5*N], the pentadiagonal matrix.
//
//    Input, double X[N], the vector to be multiplied by A'.
//
//    Output, double R85_VXM[N], the product A' * x.
//
{
  double *b;
  int j;

  b = new double[n];

  for ( j = 0; j < n; j++ )
  {
    b[j] = a[2+j*5] * x[j];
  }

  for ( j = 1; j < n; j++ )
  {
    b[j] = b[j] + a[3+j*5] * x[j-1];
  }

  for ( j = 2; j < n; j++ )
  {
    b[j] = b[j] + a[4+j*5] * x[j-2];
  }

  for ( j = 0; j < n-1; j++ )
  {
    b[j] = b[j] + a[1+j*5] * x[j+1];
  }

  for ( j = 0; j < n-2; j++ )
  {
    b[j] = b[j] + a[0+j*5] * x[j+2];
  }

  return b;
}
//****************************************************************************80

double *r85_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R85_ZERO zeros a R85 matrix.
//
//  Discussion:
//
//    The R85 storage format represents a pentadiagonal matrix as a 5
//    by N array, in which each row corresponds to a diagonal, and
//    column locations are preserved.  Thus, the original matrix is
//    "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R85 matrix of order 6 would be stored:
//
//       *   *  A13 A24 A35 A46
//       *  A12 A23 A34 A45 A56
//      A11 A22 A33 A44 A55 A66
//      A21 A32 A43 A54 A65  *
//      A31 A42 A53 A64  *   *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Output, double R85_ZERO[5*N], the R85 matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[5*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 5; i++ )
    {
      a[i+j*5] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

void r8bb_add ( int n1, int n2, int ml, int mu, double a[], int i, int j, 
  double value )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_ADD adds a value to an entry in a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, 
//    and N2 by N2, respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input/output, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
//    Input, int I, J, the row and column of the entry to be incremented.
//    Some combinations of I and J are illegal.
//
//    Input, double VALUE, the value to be added to the (I,J)-th entry.
//
{
  int ij;

  if ( value == 0.0 )
  {
    return;
  }

  if ( i <= 0 || n1 + n2 < i )
  {
    cerr << "\n";
    cerr << "R8BB_ADD - Fatal error!\n";
    cerr << "  Illegal input value of row index I = " << i << "\n";
    exit ( 1 );
  }

  if ( j <= 0 || n1 + n2 < j )
  {
    cerr << "\n";
    cerr << "R8BB_ADD - Fatal error!\n";
    cerr << "  Illegal input value of column index J = " << j << "\n";
    exit ( 1 );
  }
//
//  The A1 block of the matrix.
//
//  Check for out of band problems.
//
//  Normally, we would check the condition MU < (J-I), but the storage
//  format requires extra entries be set aside in case of pivoting, which
//  means that the condition becomes MU+ML < (J-I).
//
  if ( i <= n1 && j <= n1 )
  {
    if ( (mu+ml) < (j-i) || ml < (i-j) )
    {
      cout << "\n";
      cout << "R8BB_ADD - Warning!\n";
      cout << "  Unable to add to entry (" << i << ", " << j << ").\n";
    }
    else
    {
      ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1);
    }
  }
//
//  The A2 block of the matrix.
//
  else if ( i <= n1 && n1 < j )
  {
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1 + i;
  }
//
//  The A3 and A4 blocks of the matrix.
//
  else if ( n1 < i )
  {
    ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1);
  }

  a[ij-1] = a[ij-1] + value;

  return;
}
//****************************************************************************80

int r8bb_fa ( int n1, int n2, int ml, int mu, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_FA factors a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Input/output, double A[(2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2 ].
//    On input, the border-banded matrix to be factored.
//    On output, information describing a partial factorization
//    of the original coefficient matrix.  This information is required
//    by R8BB_SL in order to solve linear systems associated with that
//    matrix.
//
//    Output, int PIVOT[N1+N2], contains pivoting information.
//
//    Output, int R8BB_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  double *b;
  int i;
  int ij;
  int ik;
  int info;
  int j;
  int jk;
  int job;
  int k;
  int nband;
  double *x;

  nband = (2*ml+mu+1) * n1;
//
//  Factor the A1 band matrix, overwriting A1 by its factors.
//
  if ( 0 < n1 )
  {
    info = r8gb_fa ( n1, ml, mu, a, pivot );

    if ( info != 0 )
    {
      return info;
    }
  }

  if ( 0 < n1 && 0 < n2 )
  {
//
//  Solve A1 * x = -A2 for x, and overwrite A2 by the results.
//
    for ( i = nband+1; i <= nband+n1*n2; i++ )
    {
      a[i-1] = -a[i-1];
    }

    b = new double[n1];
    x = new double[n1];

    job = 0;
    for ( j = 1; j <= n2; j++ )
    {
      for ( i = 0; i < n1; i++ )
      {
        b[i] = a[nband+(j-1)*n1+i];
      }
      x = r8gb_sl ( n1, ml, mu, a, pivot, b, job );
      for ( i = 0; i < n1; i++ )
      {
        a[nband+(j-1)*n1+i] = x[i];
      }
      delete [] x;
    }
    delete [] b;
//
//  A4 := A4 + A3 * A2.
//
    for ( i = 1; i <= n2; i++ )
    {
      for ( j = 1; j <= n1; j++ )
      {
        ij = nband + n1*n2 + (j-1)*n2 + i;
        for ( k = 1; k <= n2; k++ )
        {
          ik = nband + 2*n1*n2 + (k-1)*n2 + i;
          jk = nband + (k-1) * n1 + j;
          a[ik-1] = a[ik-1] + a[ij-1] * a[jk-1];
        }
      }
    }
  }
//
//  Factor A4.
//
  if ( 0 < n2 )
  {
    info = r8ge_fa ( n2, a+(nband+2*n1*n2), pivot+n1 );

    if ( info != 0 )
    {
      return info;
    }
  }

  return 0;
}
//****************************************************************************80

double r8bb_get ( int n1, int n2, int ml, int mu, double a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_GET gets a value of a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input/output, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
//    Input, int I, J, the row and column of the entry to be incremented.
//    Some combinations of I and J are illegal.
//
//    Output, double R8BB_GET, the value of the (I,J)-th entry.
//
{
  int ij;

  if ( i <= 0 || n1 + n2 < i )
  {
    cerr << "\n";
    cerr << "R8BB_GET - Fatal error!\n";
    cerr << "  Illegal input value of row index I = " << i << "\n";
    exit ( 1 );
  }

  if ( j <= 0 || n1 + n2 < j )
  {
    cerr << "\n";
    cerr << "R8BB_GET - Fatal error!\n";
    cerr << "  Illegal input value of column index J = " << j << "\n";
    exit ( 1 );
  }
//
//  The A1 block of the matrix.
//
//  Check for out of band problems.
//
//  Normally, we would check the condition MU < (J-I), but the storage
//  format requires extra entries be set aside in case of pivoting, which
//  means that the condition becomes MU+ML < (J-I).
//
  if ( i <= n1 && j <= n1 )
  {
    if ( (mu+ml) < (j-i) || ml < (i-j) )
    {
      return 0.0;
    }
    else
    {
      ij = (i-j+ml+mu+1) + (j-1)*(2*ml+mu+1);
    }
  }
//
//  The A2 block of the matrix.
//
  else if ( i <= n1 && n1 < j )
  {
    ij = (2*ml+mu+1)*n1 + (j-n1-1)*n1 + i;
  }
//
//  The A3 and A4 blocks of the matrix.
//
  else if ( n1 < i )
  {
    ij = (2*ml+mu+1)*n1 + n2*n1 + (j-1)*n2 + (i-n1);
  }

  return a[ij-1];
}
//****************************************************************************80

double *r8bb_indicator ( int n1, int n2, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_INDICATOR sets up a R8BB indicator matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//    The matrix is actually stored as a vector, and we will simply suggest
//    the structure and values of the indicator matrix as:
//
//      00 00 00 00 00
//      00 00 13 24 35     16 17     61 62 63 64 65     66 67
//      00 12 23 34 45  +  26 27  +  71 72 73 74 75  +  76 77
//      11 22 33 44 55     36 37     
//      21 32 43 54 00     46 47     
//                         56 57     
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Output, double R8BB_INDICATOR[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], 
//    the matrix.
//
{
  double *a;
  int base;
  int fac;
  int i;
  int j;
  int row;

  a = new double[(2*ml+mu+1)*n1+2*n1*n2+n2*n2];

  fac = i4_power ( 10, i4_log_10 ( n1 + n2 ) + 1 );
//
//  Set the banded matrix A1.
//
  for ( j = 1; j <= n1; j++ )
  {
    for ( row = 1; row <= 2 * ml + mu + 1; row++ )
    {
      i = row + j - ml - mu - 1;
      if ( ml < row && 1 <= i && i <= n1 )
      {
        a[row-1+(j-1)*(2*ml+mu+1)] = ( double ) ( fac * i + j );
      }
      else
      {
        a[row-1+(j-1)*(2*ml+mu+1)] = 0.0;
      }
    }
  }
//
//  Set the N1 by N2 rectangular strip A2.
//
  base = ( 2 * ml + mu + 1 ) * n1;

  for ( i = 1; i <= n1; i++ )
  {
    for ( j = n1 + 1; j <= n1 + n2; j++ )
    {
      a[base + i-1 + (j-n1-1)*n1 ] = ( double ) ( fac * i + j );
    }
  }
//
//  Set the N2 by N1 rectangular strip A3.
//
  base = ( 2 * ml + mu + 1 ) * n1 + n1 * n2;

  for ( i = n1 + 1; i <= n1 + n2; i++ )
  {
    for ( j = 1; j <= n1; j++ )
    {
      a[base + i-n1-1 + (j-1)*n2 ] = ( double ) ( fac * i + j );
    }
  }
//
//  Set the N2 by N2 square A4.
//
  base = ( 2 * ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1;

  for ( i = n1 + 1; i <= n1 + n2; i++ )
  {
    for ( j = n1 + 1; j <= n1 + n2; j++ )
    {
      a[base + i-n1-1 + (j-n1-1)*n2 ] = ( double ) ( fac * i + j );
    }
  }

  return a;
}
//****************************************************************************80

double *r8bb_mxv ( int n1, int n2, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_MXV multiplies a R8BB matrix times a vector.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Input, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
//    Input, double X[N1+N2], the vector to be multiplied by A.
//
//    Output, double R8BB_MXV[N1+N2], the result of multiplying A by X.
//
{
  double *b;
  int i;
  int ihi;
  int ij;
  int ilo;
  int j;
//
//  Initialize B.
//
  b = new double[n1+n2];
  for ( i = 0; i < n1 + n2; i++ )
  {
    b[i] = 0.0;
  }
//
//  Multiply by A1.
//
  for ( j = 1; j <= n1; j++ )
  {
    ilo = i4_max ( 1, j - mu - ml );
    ihi = i4_min ( n1, j + ml );
    ij = (j-1) * (2*ml+mu+1) - j + ml + mu + 1;
    for ( i = ilo; i <= ihi; i++ )
    {
      b[i-1] = b[i-1] + a[ij+i-1] * x[j-1];
    }
  }
//
//  Multiply by A2.
//
  for ( j = n1+1; j <= n1 + n2; j++ )
  {
    ij = (2*ml+mu+1)*n1 + (j-n1-1)*n1;
    for ( i = 1; i <= n1; i++ )
    {
      b[i-1] = b[i-1] + a[ij+i-1] * x[j-1];
    }
  }
//
//  Multiply by A3 and A4.
//
  for ( j = 1; j <= n1 + n2; j++ )
  {
    ij = (2*ml+mu+1)*n1 + n1*n2 + (j-1)*n2 - n1;
    for ( i = n1+1; i <= n1+n2; i++ )
    {
      b[i-1] = b[i-1] + a[ij+i-1] * x[j-1];
    }
  }

  return b;
}
//****************************************************************************80

void r8bb_print ( int n1, int n2, int ml, int mu, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_PRINT prints a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
//    Input, string TITLE, a title.
//
{
  r8bb_print_some ( n1, n2, ml, mu, a, 1, 1, n1+n2, n1+n2, title );

  return;
}
//****************************************************************************80

void r8bb_print_some ( int n1, int n2, int ml, int mu, double a[], int ilo, 
  int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_PRINT_SOME prints some of a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int ij;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n1+n2 );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n1+n2 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        aij = 0.0;

        if ( i <= n1 && j <= n1 )
        {
          if ( (j-i) <= mu+ml && (i-j) <= ml )
          {
            ij = (i-j+ml+mu+1) + (j-1)*(2*ml+mu+1);
            aij = a[ij-1];
          }
        }
        else if ( i <= n1 && n1 < j )
        {
          ij = (2*ml+mu+1)*n1 + (j-n1-1)*n1 + i;
          aij = a[ij-1];
        }
        else if ( n1 < i )
        {
          ij = (2*ml+mu+1)*n1 + n2*n1 + (j-1)*n2 + (i-n1);
          aij = a[ij-1];
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8bb_random ( int n1, int n2, int ml, int mu, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_RANDOM randomizes a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8BB_RANDOM[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
{
  double *a;
  int i;
  int j;
  double r;
  int row;

  a = new double[(2*ml+mu+1)*n1+2*n1*n2+n2*n2];
//
//  Randomize the banded matrix A1.
//  We still believe that the "junk" entries should be set to 0.
//
  for ( j = 1; j <= n1; j++ )
  {
    for ( row = 1; row <= 2 * ml + mu + 1; row++ )
    {
      i = row + j - ml - mu - 1;
      if ( ml < row && 1 <= i && i <= n1 )
      {
        r = r8_uniform_01 ( seed );
      }
      else
      {
        r = 0.0;
      }
      a[row-1+(j-1)*(2*ml+mu+1)] = r;
    }
  }
//
//  Randomize the rectangular strips A2+A3+A4.
//
  for ( i = (2*ml+mu+1)*n1; i < (2*ml+mu+1)*n1+2*n1*n2+n2*n2; i++ )
  {
    a[i] = r8_uniform_01 ( seed );
  }

  return a;
}
//****************************************************************************80

void r8bb_set ( int n1, int n2, int ml, int mu, double a[], int i, int j, 
  double value )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_SET sets a value of a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input/output, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
//    Input, int I, J, the row and column of the entry to be incremented.
//    Some combinations of I and J are illegal.
//
//    Input, double VALUE, the value to be assigned to the (I,J)-th entry.
//
{
  int ij;

  if ( i <= 0 || n1 + n2 < i )
  {
    cerr << "\n";
    cerr << "R8BB_SET - Fatal error!\n";
    cerr << "  Illegal input value of row index I = " << i << "\n";
    exit ( 1 );
  }

  if ( j <= 0 || n1 + n2 < j )
  {
    cerr << "\n";
    cerr << "R8BB_SET - Fatal error!\n";
    cerr << "  Illegal input value of column index J = " << j << "\n";
    exit ( 1 );
  }
//
//  The A1 block of the matrix.
//
//  Check for out of band problems.
//
//  Normally, we would check the condition MU < (J-I), but the storage
//  format requires extra entries be set aside in case of pivoting, which
//  means that the condition becomes MU+ML < (J-I).
//
  if ( i <= n1 && j <= n1 )
  {
    if ( (mu+ml) < (j-i) || ml < (i-j) )
    {
      cout << "\n";
      cout << "R8BB_SET - Warning!\n";
      cout << "  Unable to set entry (" << i << ", " << j << ").\n";
    }
    else
    {
      ij = (i-j+ml+mu+1) + (j-1)*(2*ml+mu+1);
    }
  }
//
//  The A2 block of the matrix.
//
  else if ( i <= n1 && n1 < j )
  {
    ij = (2*ml+mu+1)*n1 + (j-n1-1)*n1 + i;
  }
//
//  The A3 and A4 blocks of the matrix.
//
  else if ( n1 < i )
  {
    ij = (2*ml+mu+1)*n1 + n2*n1 + (j-1)*n2 + (i-n1);
  }

  a[ij-1] = value;

  return;
}
//****************************************************************************80

double *r8bb_sl ( int n1, int n2, int ml, int mu, double a_lu[], int pivot[], 
  double b[] )

//****************************************************************************80
//
//  Discussion:
//
//    R8BB_SL solves a R8BB system factored by SBB_FA.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Input, double A_LU[(2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the LU factors from R8BB_FA.
//
//    Input, int PIVOT[N1+N2], the pivoting information from R8BB_FA.
//
//    Input, double B[N1+N2], the right hand side.
//
//    Output, double R8BB_SL[N1+N2], the solution.
//
{
  double *b22;
  int i;
  int ij;
  int j;
  int job;
  int nband;
  double *x;
  double *x1;
  double *x2;

  nband = (2*ml+mu+1)*n1;
//
//  Set X1 := inverse(A1) * B1.
//
  if ( 0 < n1 )
  {
    job = 0;
    x1 = r8gb_sl ( n1, ml, mu, a_lu, pivot, b, job );
  }
//
//  Modify the right hand side of the second linear subsystem.
//  Set B22 := B2 - A3*X1.
//
  if ( 0 < n2 )
  {
    b22 = new double[n2];

    for ( i = 0; i < n2; i++ )
    {
      b22[i] = b[n1+i];
      for ( j = 0; j < n1; j++ )
      {
        ij = nband + n1*n2 + j*n2 + i;
        b22[i] = b22[i] - a_lu[ij] * x1[j];
      }
    }
  }
//
//  Set X2 := inverse(A4) * B22.
//
  if ( 0 < n2 )
  {
    job = 0;
    x2 = r8ge_sl_new ( n2, a_lu+(nband+2*n1*n2), pivot+n1, b22, job );
    delete [] b22;
  }
//
//  Modify the first subsolution.
//  Set X1 := X1 + A2*X2.
//
  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n2; j++ )
    {
      ij = nband + j*n1 + i;
      x1[i] = x1[i] + a_lu[ij] * x2[j];
    }
  }
//
//  Set X = [ X1 | X2 ].
//
  x = new double[n1+n2];

  if ( 0 < n1 )
  {
    for ( i = 0; i < n1; i++ )
    {
      x[i] = x1[i];
    }
    delete [] x1;
  }

  if ( 0 < n2 )
  {
    for ( i = 0; i < n2; i++ )
    {
      x[n1+i] = x2[i];
    }
    delete [] x2;
  }

  return x;
}
//****************************************************************************80

double *r8bb_to_r8ge ( int n1, int n2, int ml, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_TO_R8GE copies a R8BB matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
//    Output, double R8BB_TO_R8GE[(N1+N2)*(N1+N2)], the R8GE matrix.
//
{
  double *b;
  int i;
  int ij;
  int j;

  b = new double[(1+n2)*(n1+n2)];

  for ( i = 1; i <= n1; i++ )
  {
    for ( j = 1; j <= n1; j++ )
    {
      if ( mu+ml < (j-i) || ml < (i-j) )
      {
        b[i-1+(j-1)*(n1+n2)] = 0.0;
      }
      else
      {
        ij = (i-j+ml+mu+1) + (j-1)*(2*ml+mu+1);
        b[i-1+(j-1)*(n1+n2)]  = a[ij-1];
      }
    }
  }

  for ( i = 1; i <= n1; i++ )
  {
    for ( j = n1+1; j <= n2; j++ )
    {
      ij = (2*ml+mu+1)*n1 + (j-n1-1)*n1 + i;
      b[i-1+(j-1)*(n1+n2)]  = a[ij-1];
    }
  }

  for ( i = n1+1; i <= n2; i++ )
  {
    for ( j = 1; j <= n1+n2; j++ )
    {
      ij = (2*ml+mu+1)*n1 + n2*n1 + (j-1)*n2 + (i-n1);
      b[i-1+(j-1)*(n1+n2)]  = a[ij-1];
    }
  }

  return b;
}
//****************************************************************************80

double *r8bb_vxm ( int n1, int n2, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_VXM multiplies a vector by a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Input, double A[(2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8BB matrix.
//
//    Input, double X[N1+N2], the vector to multiply A.
//
//    Output, double R8BB_VXM[N1+N2], the product X times A.
//
{
  double *b;
  int i;
  int ihi;
  int ij;
  int ilo;
  int j;
//
//  Initialize B.
//
  b = new double[n1+n2];
  for ( i = 0; i < n1+n2; i++ )
  {
    b[i] = 0.0;
  }
//
//  Multiply by A1.
//
  for ( j = 1; j <= n1; j++ )
  {
    ilo = i4_max ( 1, j - mu - ml );
    ihi = i4_min ( n1, j + ml );
    ij = (j-1) * (2*ml+mu+1) - j + ml + mu + 1;
    for ( i = ilo; i <= ihi; i++ )
    {
      b[j-1] = b[j-1] + x[i-1] * a[ij+i-1];
    }
  }
//
//  Multiply by A2.
//
  for ( j = n1+1; j <= n1+n2; j++ )
  {
    ij = (2*ml+mu+1)*n1 + (j-n1-1)*n1;
    for ( i = 1; i <= n1; i++ )
    {
      b[j-1] = b[j-1] + x[i-1] * a[ij+i-1];
    }
  }
//
//  Multiply by A3 and A4.
//
  for ( j = 1; j <= n1+n2; j++ )
  {
    ij = (2*ml+mu+1)*n1 + n1*n2 + (j-1)*n2 - n1;
    for ( i = n1+1; i <= n1+n2; i++ )
    {
      b[j-1] = b[j-1] + x[i-1] * a[ij+i-1];
    }
  }

  return b;
}
//****************************************************************************80

double *r8bb_zero ( int n1, int n2, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8BB_ZERO zeros a R8BB matrix.
//
//  Discussion:
//
//    The R8BB storage format is for a border banded matrix.  Such a
//    matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
//    general band format.  The reason for the factor of 2 in front of
//    ML is to allocate space that may be required if pivoting occurs.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Example:
//
//    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
//
//       00
//       00  00
//       00  00  00 --- ---
//      A11 A12 A13  00 ---  A16 A17
//      A21 A22 A23 A24  00  A26 A27
//      --- A32 A33 A34 A35  A36 A37
//      --- --- A43 A44 A45  A46 A47
//      --- --- --- A54 A55  A56 A57
//                       00
//
//      A61 A62 A63 A64 A65  A66 A67
//      A71 A72 A73 A74 A75  A76 A77
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Output, double R8BB_ZERO[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
//
{
  double *a;
  int i;

  a = new double[(2*ml+mu+1)*n1+2*n1*n2+n2*n2];

  for ( i = 0; i < (2*ml+mu+1)*n1 + 2*n1*n2 + n2*n2; i++ )
  {
    a[i] = 0.0;
  }

  return a;
}
//****************************************************************************80

double r8blt_det ( int n, int ml, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_DET computes the determinant of a R8BLT matrix.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Input, double A[(ML+1)*N], the R8BLT matrix.
//
//    Output, double R8BLT_DET, the determinant of A.
//
{
  double det;
  int j;

  det = 1.0;
  for ( j = 0; j < n; j++ )
  {
    det = det * a[0+j*(ml+1)];
  }

  return det;
}
//****************************************************************************80

double *r8blt_indicator ( int n, int ml )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_INDICATOR sets up a R8BLT indicator matrix.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Example:
//
//    N = 5, ML = 2
//
//    A11   0   0   0   0
//    A21 A22   0   0   0
//    A31 A32 A33   0   0
//      0 A42 A43 A44   0
//      0   0 A53 A54 A55
//                --- ---
//                    ---
//
//    The indicator matrix is stored as:
//
//      11  22  33  44  55
//      21  32  43  54   0
//      31  42  53   0   0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Output, double R8BLT_INDICATOR[(ML+1)*N], the R8BLT matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;
  int jlo;

  a = new double[(ml+1)*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= n; i++ )
  {
    jlo = i4_max ( 1, i - ml );
    for ( j = jlo; j <= i; j++ )
    {
      a[i-j+(j-1)*(ml+1)] = ( double ) ( fac * i + j );
    }
  }
//
//  The junk entries can be thought of as corresponding to
//  elements of a phantom portion of the matrix.
//
  for ( i = n; i < n + ml; i++ )
  {
    for ( j = i - ml; j < n; j++ )
    {
      a[i-j+j*(ml+1)] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8blt_mxv ( int n, int ml, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_MXV multiplies a R8BLT matrix times a vector.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Input, double A[(ML+1)*N], the R8BLT matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8BLT_MXV[N], the product A * x.
//
{
  double *b;
  int i;
  int j;
  int jhi;
  int jlo;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    jlo = i4_max ( 0, i - ml );
    jhi = i;
    for ( j = jlo; j <= jhi; j++ )
    {
      b[i] = b[i] + a[i-j+j*(ml+1)] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void r8blt_print ( int n, int ml, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_PRINT prints a R8BLT matrix.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Input, double A[(ML+1)*N], the R8BLT matrix.
//
//    Input, string TITLE, a title.
//
{
  r8blt_print_some ( n, ml, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8blt_print_some ( int n, int ml, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_PRINT_SOME prints some of a R8BLT matrix.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Input, double A[(ML+1)*N], the R8BLT matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo );
    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + ml );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(5) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( ml < i-j || 0 < j-i )
        {
          cout << "              ";
        }
        else
        {
          cout << setw(12) << a[i-j+(j-1)*(ml+1)] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8blt_random ( int n, int ml, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_RANDOM randomizes a R8BLT matrix.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8BLT_RANDOM[(ML+1)*N], the R8BLT matrix.
//
{
  double *a;
  int i;
  int j;
  int jlo;

  a = new double[(ml+1)*n];

  for ( i = 0; i < n; i++ )
  {
    jlo = i4_max ( 0, i - ml );
    for ( j = jlo; j <= i; j++ )
    {
      a[i-j+j*(ml+1)] = r8_uniform_01 ( seed );
    }
  }
//
//  The junk entries can be thought of as corresponding to
//  elements of a phantom portion of the matrix.
//
  for ( i = n; i < n + ml; i++ )
  {
    for ( j = i - ml; j < n; j++ )
    {
      a[i-j+j*(ml+1)] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8blt_sl ( int n, int ml, double a[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_SL solves a R8BLT system.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//    No factorization of the lower triangular matrix is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Input, double A[(ML+1)*N], the R8BLT matrix.
//
//    Input, double B(N), the right hand side.
//
//    Input, int JOB, is 0 to solve the untransposed system,
//    nonzero to solve the transposed system.
//
//    Output, double R8BLT_SL[N], the solution vector.
//
{
  int i;
  int ihi;
  int ilo;
  int j;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = x[j] / a[0+j*(ml+1)];
      ihi = i4_min ( j + ml, n - 1 );
      for ( i = j+1; i <= ihi; i++ )
      {
        x[i] = x[i] - a[i-j+j*(ml+1)] * x[j];
      }
    }
  }
  else
  {
    for ( j = n-1; 0 <= j; j-- )
    {
      x[j] = x[j] / a[0+j*(ml+1)];
      ilo = i4_max ( j - ml, 0 );
      for ( i = ilo; i <= j-1; i++ )
      {
        x[i] = x[i] - a[j-i+i*(ml+1)] * x[j];
      }
    }
  }

  return x;
}
//****************************************************************************80

double *r8blt_to_r8ge ( int n, int ml, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_TO_R8GE copies a R8BLT matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8BLT storage format is used for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Example:
//
//    N = 5, ML = 2
//
//    A11   0   0   0   0
//    A21 A22   0   0   0
//    A31 A32 A33   0   0
//      0 A42 A43 A44   0
//      0   0 A53 A54 A55
//                --- ---
//                    ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, int ML, the lower bandwidth of A.
//    ML must be nonnegative, and no greater than N-1.
//
//    Input, double A[(ML+1)*N], the R8BLT matrix.
//
//    Output, double R8BLT_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( j <= i && i <= j + ml )
      {
        b[i-1+(j-1)*n] = a[i-j+(j-1)*(ml+1)];
      }
      else
      {
        b[i-1+(j-1)*n] = 0.0;
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8blt_vxm ( int n, int ml, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_VXM multiplies a vector by a R8BLT matrix.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Input, double A[(ML+1)*N], the R8BLT matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8BLT_VXM[N], the product X*A.
//
{
  double *b;
  int i;
  int j;
  int jhi;
  int jlo;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    jlo = i4_max ( 0, i - ml );
    jhi = i;
    for ( j = jlo; j <= jhi; j++ )
    {
      b[j] = b[j] + x[i] * a[i-j+j*(ml+1)];
    }
  }

  return b;
}
//****************************************************************************80

double *r8blt_zero ( int n, int ml )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLT_ZERO zeros a R8BLT matrix.
//
//  Discussion:
//
//    The R8BLT storage format is appropriate for a banded lower triangular matrix.
//    The matrix is assumed to be zero below the ML-th subdiagonal.
//    The matrix is stored in an ML+1 by N array, in which the diagonal
//    appears in the first row, followed by successive subdiagonals.
//    Columns are preserved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int ML, the lower bandwidth.
//
//    Output, double R8BLT_ZERO[(ML+1)*N], the R8BLT matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[(ml+1)*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < ml+1; i++ )
    {
      a[i+j*(ml+1)] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8bto_indicator ( int m, int l )

//****************************************************************************80
//
//  Purpose:
//
//    R8BTO_INDICATOR sets up a R8BTO indicator matrix.
//
//  Discussion:
//
//    The R8BTO storage format is for a block Toeplitz matrix. The matrix
//    can be regarded as an L by L array of blocks, each of size M by M.
//    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
//    that is, along its diagonal, the blocks repeat.
//
//    Storage for the matrix consists of the L blocks of the first row,
//    followed by the L-1 blocks of the first column (skipping the first row).
//    These items are stored in the natural way in an (M,M,2*L-1) array.
//
//  Example:
//
//    M = 2, L = 3
//
//    1 2 | 3 4 | 5 6
//    5 5 | 6 6 | 7 7
//    ----+-----+-----
//    7 8 | 1 2 | 3 4
//    8 8 | 5 5 | 6 6
//    ----+-----+-----
//    9 0 | 7 8 | 1 2
//    9 9 | 8 8 | 5 5
//
//    X = (/ 1, 2, 3, 4, 5, 6 /)
//
//    B = (/ 91, 134, 73, 125, 97, 129 /)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the order of the blocks of the matrix A.
//
//    Input, int L, the number of blocks in a row or column of A.
//
//    Output, double R8BTO_INDICATOR[M*M*(2*L-1)], the R8BTO matrix.
//
{
  double *a;
  int fac;
  int i;
  int i2;
  int j;
  int j2;
  int k;

  a = new double[m*m*(2*l-1)];

  fac = i4_power ( 10, i4_log_10 ( m * l ) + 1 );
//
//  Blocks 1 to L form the first row.
//
  j = 0;

  for ( k = 1; k <= l; k++ )
  {
    for ( j2 = 1; j2 <= m; j2++ )
    {
      j = j + 1;
      for ( i = 1; i <= m; i++ )
      {
        a[i-1+(j2-1)*m+(k-1)*m*m] = ( double ) ( fac * i + j );
      }
    }
  }
//
//  Blocks L+1 through 2*L-1 form the remainder of the first column.
//
  i = m;

  for ( k = l+1; k <= 2*l-1; k++ )
  {
    for ( i2 = 1; i2 <= m; i2++ )
    {
      i = i + 1;
      for ( j = 1; j <= m; j++ )
      {
        a[i2-1+(j-1)*m+(k-1)*m*m] = ( double ) ( fac * i + j );
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8bto_mxv ( int m, int l, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BTO_MXV multiplies a R8BTO matrix times a vector.
//
//  Discussion:
//
//    The R8BTO storage format is for a block Toeplitz matrix. The matrix
//    can be regarded as an L by L array of blocks, each of size M by M.
//    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
//    that is, along its diagonal, the blocks repeat.
//
//    Storage for the matrix consists of the L blocks of the first row,
//    followed by the L-1 blocks of the first column (skipping the first row).
//    These items are stored in the natural way in an (M,M,2*L-1) array.
//
//  Example:
//
//    M = 2, L = 3
//
//    1 2 | 3 4 | 5 6
//    5 5 | 6 6 | 7 7
//    ----+-----+-----
//    7 8 | 1 2 | 3 4
//    8 8 | 5 5 | 6 6
//    ----+-----+-----
//    9 0 | 7 8 | 1 2
//    9 9 | 8 8 | 5 5
//
//    X = (/ 1, 2, 3, 4, 5, 6 /)
//
//    B = (/ 91, 134, 73, 125, 79, 138 /)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the order of the blocks of the matrix A.
//
//    Input, int L, the number of blocks in a row or column of A.
//
//    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
//
//    Input, double X[M*L], the vector to be multiplied.
//
//    Output, double R8BTO_MXV[M*L], the product A * X.
//
{
  double *b;
  int i;
  int i2;
  int j;
  int k;

  b = new double[m*l];
//
//  Construct the right hand side by blocks.
//
  for ( j = 0; j < l; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+j*m] = 0.0;
    }

    for ( k = 0; k <= j-1; k++ )
    {
      for ( i = 0; i < m; i++ )
      {
        for ( i2 = 0; i2 < m; i2++ )
        {
          b[i+j*m] = b[i+j*m] + a[i+i2*m+(l+j-k-1)*m*m] * x[i2+k*m];
        }
      }
    }

    for ( k = j; k < l; k++ )
    {
      for ( i = 0; i < m; i++ )
      {
        for ( i2 = 0; i2 < m; i2++ )
        {
          b[i+j*m] = b[i+j*m] + a[i+i2*m+(k-j)*m*m] * x[i2+k*m];
        }
      }
    }
  }

  return b;
}
//****************************************************************************80

void r8bto_print ( int m, int l, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BTO_PRINT prints a R8BTO matrix.
//
//  Discussion:
//
//    The R8BTO storage format is for a block Toeplitz matrix. The matrix
//    can be regarded as an L by L array of blocks, each of size M by M.
//    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
//    that is, along its diagonal, the blocks repeat.
//
//    Storage for the matrix consists of the L blocks of the first row,
//    followed by the L-1 blocks of the first column (skipping the first row).
//    These items are stored in the natural way in an (M,M,2*L-1) array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the order of the blocks of the matrix A.
//
//    Input, int L, the number of blocks in a row or column of A.
//
//    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
//
//    Input, string TITLE, a title.
//
{
  r8bto_print_some ( m, l, a, 1, 1, m*l, m*l, title );

  return;
}
//****************************************************************************80

void r8bto_print_some ( int m, int l, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BTO_PRINT_SOME prints some of a R8BTO matrix.
//
//  Discussion:
//
//    The R8BTO storage format is for a block Toeplitz matrix. The matrix
//    can be regarded as an L by L array of blocks, each of size M by M.
//    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
//    that is, along its diagonal, the blocks repeat.
//
//    Storage for the matrix consists of the L blocks of the first row,
//    followed by the L-1 blocks of the first column (skipping the first row).
//    These items are stored in the natural way in an (M,M,2*L-1) array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the order of the blocks of the matrix A.
//
//    Input, int L, the number of blocks in a row or column of A.
//
//    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i1;
  int i2;
  int i3hi;
  int i3lo;
  int inc;
  int j;
  int j1;
  int j2;
  int j3hi;
  int j3lo;
  int n;

  n = m * l;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j3lo = jlo; j3lo <= jhi; j3lo = j3lo + INCX )
  {
    j3hi = j3lo + INCX - 1;
    j3hi = i4_min ( j3hi, n );
    j3hi = i4_min ( j3hi, jhi );

    inc = j3hi + 1 - j3lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j3lo; j <= j3hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i3lo = i4_max ( ilo, 1 );
    i3hi = i4_min ( ihi, n );

    for ( i = i3lo; i <= i3hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j3lo; j <= j3lo + inc - 1; j++ )
      {
//
//  i = M * ( i1 - 1 ) + i2
//  j = M * ( j1 - 1 ) + j2
//
        i1 = ( i - 1 ) / m + 1;
        i2 = i - m * ( i1 - 1 );
        j1 = ( j - 1 ) / m + 1;
        j2 = j - m * ( j1 - 1 );

        if ( i1 <= j1 )
        {
          cout << setw(12) << a[i2-1+(j2-1)*m+(j1-i1)*m*m] << "  ";
        }
        else
        {
          cout << setw(12) << a[i2-1+(j2-1)*m+(l-1+i1-j1)*m*m] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8bto_random ( int m, int l, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8BTO_RANDOM randomizes a R8BTO matrix.
//
//  Discussion:
//
//    The R8BTO storage format is for a block Toeplitz matrix. The matrix
//    can be regarded as an L by L array of blocks, each of size M by M.
//    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
//    that is, along its diagonal, the blocks repeat.
//
//    Storage for the matrix consists of the L blocks of the first row,
//    followed by the L-1 blocks of the first column (skipping the first row).
//    These items are stored in the natural way in an (M,M,2*L-1) array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the order of the blocks of the matrix A.
//
//    Input, int L, the number of blocks in a row or column of A.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8BTO_RANDOM[M*M*(2*L-1)], the R8BTO matrix.
//
{
  double *a;
  int i;
  int j;
  int k;

  a = new double[m*m*(2*l-1)];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < 2 * l - 1; k++ )
      {
        a[i+j*m+k*m*m] = r8_uniform_01 ( seed );
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8bto_to_r8ge ( int m, int l, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BTO_TO_R8GE copies a R8BTO matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8BTO storage format is for a block Toeplitz matrix. The matrix
//    can be regarded as an L by L array of blocks, each of size M by M.
//    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
//    that is, along its diagonal, the blocks repeat.
//
//    Storage for the matrix consists of the L blocks of the first row,
//    followed by the L-1 blocks of the first column (skipping the first row).
//    These items are stored in the natural way in an (M,M,2*L-1) array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the order of the blocks of the R8BTO matrix.
//
//    Input, int L, the number of blocks in a row or column of the
//    R8BTO matrix.
//
//    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
//
//    Output, double R8BTO_TO_R8GE[(M*L)*(M*L)], the R8GE matrix.
//
{
  double *b;
  int i;
  int i1;
  int i2;
  int j;
  int j1;
  int j2;
  int n;

  n = m * l;
  b = new double[n*n];

  for ( i = 1; i <= n; i++ )
  {
    i1 = ( i - 1 ) / m + 1;
    i2 = i - m * ( i1 - 1 );

    for ( j = 1; j <= n; j++ )
    {
      j1 = ( j - 1 ) / m + 1;
      j2 = j - m * ( j1 - 1 );

      if ( i1 <= j1 )
      {
        b[i-1+(j-1)*n] = a[i2-1+(j2-1)*m+(j1-i1)*m*m];
      }
      else
      {
        b[i-1+(j-1)*n] = a[i2-1+(j2-1)*m+(l+i1-j1-1)*m*m];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8bto_vxm ( int m, int l, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BTO_VXM multiplies a vector times a R8BTO matrix.
//
//  Discussion:
//
//    The R8BTO storage format is for a block Toeplitz matrix. The matrix
//    can be regarded as an L by L array of blocks, each of size M by M.
//    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
//    that is, along its diagonal, the blocks repeat.
//
//    Storage for the matrix consists of the L blocks of the first row,
//    followed by the L-1 blocks of the first column (skipping the first row).
//    These items are stored in the natural way in an (M,M,2*L-1) array.
//
//  Example:
//
//    M = 2, L = 3
//
//    1 2 | 3 4 | 5 6
//    5 5 | 6 6 | 7 7
//    ----+-----+-----
//    7 8 | 1 2 | 3 4
//    8 8 | 5 5 | 6 6
//    ----+-----+-----
//    9 0 | 7 8 | 1 2
//    9 9 | 8 8 | 5 5
//
//    X = (/ 1, 2, 3, 4, 5, 6 /)
//
//    B = (/ 163, 122, 121, 130, 87, 96 /)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the order of the blocks of the matrix A.
//
//    Input, int L, the number of blocks in a row or column of A.
//
//    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
//
//    Input, double X[M*L], the vector to be multiplied.
//
//    Output, double R8BTO_VXM[M*L], the product X * A.
//
{
  double *b;
  int i;
  int i2;
  int j;
  int k;

  b = new double[m*l];
//
//  Construct the right hand side by blocks.
//
  for ( j = 1; j <= l; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      b[i-1+(j-1)*m] = 0.0;
    }

    for ( k = 1; k <= j; k++ )
    {
      for ( i = 1; i <= m; i++ )
      {
        for ( i2 = 1; i2 <= m; i2++ )
        {
          b[i-1+(j-1)*m] = b[i-1+(j-1)*m] 
          + a[i2-1+(i-1)*m+(j-k)*m*m] * x[i2-1+(k-1)*m];
        }
      }
    }
    for ( k = j+1; k <= l; k++ )
    {
      for ( i = 1; i <= m; i++ )
      {
        for ( i2 = 1; i2 <= m; i2++ )
        {
          b[i-1+(j-1)*m] = b[i-1+(j-1)*m] 
          + a[i2-1+(i-1)*m+(l+k-j-1)*m*m] * x[i2-1+(k-1)*m];
        }
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8bto_zero ( int m, int l )

//****************************************************************************80
//
//  Purpose:
//
//    R8BTO_ZERO zeros a R8BTO matrix.
//
//  Discussion:
//
//    The R8BTO storage format is for a block Toeplitz matrix. The matrix
//    can be regarded as an L by L array of blocks, each of size M by M.
//    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
//    that is, along its diagonal, the blocks repeat.
//
//    Storage for the matrix consists of the L blocks of the first row,
//    followed by the L-1 blocks of the first column (skipping the first row).
//    These items are stored in the natural way in an (M,M,2*L-1) array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the order of the blocks of the matrix A.
//
//    Input, int L, the number of blocks in a row or column of A.
//
//    Output, double R8BTO_ZERO[M*M*(2*L-1)], the R8BTO matrix.
//
{
  double *a;
  int i;
  int j;
  int k;

  a = new double[m*m*(2*l-1)];

  for ( k = 0; k <= 2 * l - 1; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        a[i+m*(j+m*k)] = 0.0;
      }
    }
  }
  return a;
}
//****************************************************************************80

double r8but_det ( int n, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_DET computes the determinant of a R8BUT matrix.
//
//  Discussion:
//
//    The R8BUT storage format is used for a banded upper triangular matrix.
//    The matrix is assumed to be zero above the MU-th superdiagonal.
//    The matrix is stored in an MU+1 by N array.
//    Columns are preserved.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Example:
//
//    N = 5, MU = 2
//
//    A11 A12 A13   0   0
//      0 A22 A23 A24   0
//      0   0 A33 A34 A35
//      0   0   0 A44 A45
//      0   0   0   0 A55
//                --- ---
//                    ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int MU, the upper bandwidth.
//
//    Input, double A[(MU+1)*N], the R8BUT matrix.
//
//    Output, double R8BUT_DET, the determinant of A.
//
{
  double det;
  int j;

  det = 1.0;
  for ( j = 1; j <= n; j++ )
  {
    det = det * a[(mu+1-1)+(j-1)*(mu+1)];
  }

  return det;
}
//****************************************************************************80

double *r8but_indicator ( int n, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_INDICATOR sets up a R8BUT indicator matrix.
//
//  Discussion:
//
//    The R8BUT storage format is used for a banded upper triangular matrix.
//    The matrix is assumed to be zero above the MU-th superdiagonal.
//    The matrix is stored in an MU+1 by N array.
//    Columns are preserved.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Example:
//
//    N = 5, MU = 2
//
//    A11 A12 A13   0   0
//      0 A22 A23 A24   0
//      0   0 A33 A34 A35
//      0   0   0 A44 A45
//      0   0   0   0 A55
//                --- ---
//                    ---
//
//    The indicator matrix is stored as:
//
//       0   0  13  24  35
//       0  12  23  34  45
//      11  22  33  44  55
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int MU, the upper bandwidth.
//
//    Output, double A[(MU+1)*N], the R8BUT matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[(mu+1)*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= n; i++ )
  {
    for ( j = i; j <= i4_min ( n, i + mu ); j++ )
    {
      a[i-j+mu+1-1+(j-1)*(mu+1)] = ( double ) ( fac * i + j );
    }
  }

  for ( i = 1; i <= mu; i++ )
  {
    for ( j = 1; j <= mu+1-i; j++ )
    {
      a[i-1+(j-1)*(mu+1)] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8but_mxv ( int n, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_MXV multiplies a R8BUT matrix times a vector.
//
//  Discussion:
//
//    The R8BUT storage format is used for a banded upper triangular matrix.
//    The matrix is assumed to be zero above the MU-th superdiagonal.
//    The matrix is stored in an MU+1 by N array.
//    Columns are preserved.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Example:
//
//    N = 5, MU = 2
//
//    A11 A12 A13   0   0
//      0 A22 A23 A24   0
//      0   0 A33 A34 A35
//      0   0   0 A44 A45
//      0   0   0   0 A55
//                --- ---
//                    ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int MU, the upper bandwidth.
//
//    Input, double A[(MU+1)*N], the R8BUT matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8BUT_MXV[N], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    b[i-1] = 0.0;
    for ( j = i; j <= i4_min ( n, i + mu ); j++ )
    {
      b[i-1] = b[i-1] + a[i-j+mu+1-1+(j-1)*(mu+1)] * x[j-1];
    }
  }

  return b;
}
//****************************************************************************80

void r8but_print ( int n, int mu, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_PRINT prints a R8BUT matrix.
//
//  Discussion:
//
//    The R8BUT storage format is used for a banded upper triangular matrix.
//    The matrix is assumed to be zero above the MU-th superdiagonal.
//    The matrix is stored in an MU+1 by N array.
//    Columns are preserved.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Example:
//
//    N = 5, MU = 2
//
//    A11 A12 A13   0   0
//      0 A22 A23 A24   0
//      0   0 A33 A34 A35
//      0   0   0 A44 A45
//      0   0   0   0 A55
//                --- ---
//                    ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int MU, the upper bandwidth.
//
//    Input, double A[(MU+1)*N], the R8BUT matrix.
//
//    Input, string TITLE, a title.
//
{
  r8but_print_some ( n, mu, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8but_print_some ( int n, int mu, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_PRINT_SOME prints some of a R8BUT matrix.
//
//  Discussion:
//
//    The R8BUT storage format is used for a banded upper triangular matrix.
//    The matrix is assumed to be zero above the MU-th superdiagonal.
//    The matrix is stored in an MU+1 by N array.
//    Columns are preserved.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Example:
//
//    N = 5, MU = 2
//
//    A11 A12 A13   0   0
//      0 A22 A23 A24   0
//      0   0 A33 A34 A35
//      0   0   0 A44 A45
//      0   0   0   0 A55
//                --- ---
//                    ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int MU, the upper bandwidth.
//
//    Input, double A[(MU+1)*N], the R8BUT matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo );
    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + mu );

    for ( i = i2lo; i <= i2hi; i++ )
    {

      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( i <= j && j <= i + mu )
        {
          cout << setw(12) << a[i-j+mu+1-1+(j-1)*(mu+1)] << "  ";
        }
        else
        {
          cout << "              ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8but_random ( int n, int mu, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_RANDOM randomizes a R8BUT matrix.
//
//  Discussion:
//
//    The R8BUT storage format is used for a banded upper triangular matrix.
//    The matrix is assumed to be zero above the MU-th superdiagonal.
//    The matrix is stored in an MU+1 by N array.
//    Columns are preserved.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Example:
//
//    N = 5, MU = 2
//
//    A11 A12 A13   0   0
//      0 A22 A23 A24   0
//      0   0 A33 A34 A35
//      0   0   0 A44 A45
//      0   0   0   0 A55
//                --- ---
//                    ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int MU, the upper bandwidth.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8BUT_RANDOM[(MU+1)*N], the R8BUT matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[(mu+1)*n];

  for ( i = 1; i <= mu + 1; i++ )
  {
    for ( j = 1; j <= mu + 1 - i; j++ )
    {
      a[i-1+(j-1)*(mu+1)] = 0.0;
    }

    for ( j = i4_max ( 1, mu + 2 - i ); j <= n; j++ )
    {
      a[i-1+(j-1)*(mu+1)] = r8_uniform_01 ( seed );
    }

  }

  return a;
}
//****************************************************************************80

double *r8but_sl ( int n, int mu, double a[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_SL solves a R8BUT system.
//
//  Discussion:
//
//    The R8BUT storage format is used for a banded upper triangular matrix.
//    The matrix is assumed to be zero above the MU-th superdiagonal.
//    The matrix is stored in an MU+1 by N array.
//    Columns are preserved.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Example:
//
//    N = 5, MU = 2
//
//    A11 A12 A13   0   0
//      0 A22 A23 A24   0
//      0   0 A33 A34 A35
//      0   0   0 A44 A45
//      0   0   0   0 A55
//                --- ---
//                    ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int MU, the upper bandwidth.
//
//    Input, double A[(MU+1)*N], the R8BUT matrix.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB, is 0 to solve the untransposed system,
//    nonzero to solve the transposed system.
//
//    Output, double X[N], the solution vector.
//
{
  int i;
  int ihi;
  int j;
  int jlo;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = n; 1 <= j; j-- )
    {
      x[j-1] = x[j-1] / a[j-j+mu+(j-1)*(mu+1)];
      jlo = i4_max ( 1, j - mu );
      for ( i = jlo; i <= j-1; i++ )
      {
        x[i-1] = x[i-1] - a[i-j+mu+(j-1)*(mu+1)] * x[j-1];
      }
    }
  }
  else
  {
    for ( j = 1; j <= n; j++ )
    {
      x[j-1] = x[j-1] / a[j-j+mu+(j-1)*(mu+1)];
      ihi = i4_min ( n, j + mu );
      for ( i = j + 1; i <= ihi; i++ )
      {
        x[i-1] = x[i-1] - a[j-i+mu+(i-1)*(mu+1)] * x[j-1];
      }
    }

  }

  return x;
}
//****************************************************************************80

double *r8but_to_r8ge ( int n, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_TO_R8GE copies a R8BUT matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8BUT storage format is for a banded upper triangular matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, int MU, the upper bandwidth of A.
//    MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(MU+1)*N], the R8BUT matrix.
//
//    Output, double R8BUT_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( i <= j && j <= i+mu )
      {
        b[i-1+(j-1)*n] = a[mu+i-j+(j-1)*(mu+1)];
      }
      else
      {
        b[i-1+(j-1)*n] = 0.0;
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8but_vxm ( int n, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8BUT_VXM multiplies a vector by a R8BUT matrix.
//
//  Discussion:
//
//    The R8BUT storage format is used for a banded upper triangular matrix.
//    The matrix is assumed to be zero above the MU-th superdiagonal.
//    The matrix is stored in an MU+1 by N array.
//    Columns are preserved.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Example:
//
//    N = 5, MU = 2
//
//    A11 A12 A13   0   0
//      0 A22 A23 A24   0
//      0   0 A33 A34 A35
//      0   0   0 A44 A45
//      0   0   0   0 A55
//                --- ---
//                    ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int MU, the upper bandwidth.
//
//    Input, double A[(MU+1)*N], the R8BUT matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8BUT_VXM(N), the product X*A.
//
{
  double *b;
  int i;
  int ilo;
  int j;

  b = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    b[i-1] = 0.0;
    ilo = i4_max ( 1, i - mu );
    for ( j = ilo; j <= i; j++ )
    {
      b[i-1] = b[i-1] + x[j-1] * a[j-i+mu+(i-1)*(mu+1)];
    }
  }

  return b;
}
//****************************************************************************80

double r8cb_det ( int n, int ml, int mu, double a_lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_DET computes the determinant of a R8CB matrix factored by R8CB_NP_FA.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(ML+MU+1)*N], the LU factors from R8CB_FA.
//
//    Output, double R8CB_DET, the determinant of the matrix.
//
{
  double det;
  int j;

  det = 1.0;
  for ( j = 0; j < n; j++ )
  {
    det = det * a_lu[mu+j*(ml+mu+1)];
  }

  return det;
}
//****************************************************************************80

double *r8cb_indicator ( int m, int n, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_INDICATOR sets up a R8CB indicator matrix.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically ML+MU+1 by N.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2004
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Output, double R8CB_INDICATOR[(ML+MU+1)*N], the R8CB matrix.
//
{
  double *a;
  int col = ml + mu + 1;
  int diag;
  int fac;
  int i;
  int j;
  int k;

  a = new double[(ml+mu+1)*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );
  k = 0;

  for ( j = 1; j <= n; j++ )
  {
    for ( diag = 1; diag <= ml + mu + 1; diag++ )
    {
      i = diag + j - mu - 1;

      if ( 1 <= i && i <= m && i - ml <= j && j <= i + mu )
      {
        a[diag-1+(j-1)*col] = ( double ) ( fac * i + j );
      }
      else
      {
        k = k + 1;
        a[diag-1+(j-1)*col] = - ( double ) k;
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8cb_ml ( int n, int ml, int mu, double a_lu[], double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_ML computes A * x or A' * X, using R8CB_NP_FA factors.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//    It is assumed that R8CB_NP_FA has overwritten the original matrix
//    information by LU factors.  R8CB_ML is able to reconstruct the
//    original matrix from the LU factor data.
//
//    R8CB_ML allows the user to check that the solution of a linear
//    system is correct, without having to save an unfactored copy
//    of the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(ML+MU+1)*N], the LU factors from R8CB_NP_FA.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Input, int JOB, specifies the operation to be done:
//    JOB = 0, compute A * x.
//    JOB nonzero, compute A' * x.
//
//    Output, double R8CB_ML[N], the result of the multiplication.
//
{
  double *b;
  int i;
  int ihi;
  int ilo;
  int j;
  int jhi;
  int nrow = ml + mu + 1;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
//
//  Y = U * X.
//
    for ( j = 0; j < n; j++ )
    {
      ilo = i4_max ( 0, j - mu );
      for ( i = ilo; i < j; i++ )
      {
        b[i] = b[i] + a_lu[i-j+mu+j*nrow] * b[j];
      }
      b[j] = a_lu[j-j+mu+j*nrow] * b[j];
    }
//
//  B = PL * Y = PL * U * X = A * x.
//
    for ( j = n - 2; 0 <= j; j-- )
    {
      ihi = i4_min ( n - 1, j + ml );
      for ( i = j + 1; i <= ihi; i++ )
      {
        b[i] = b[i] - a_lu[i-j+mu+j*nrow] * b[j];
      }
    }
  }
  else
  {
//
//  Y = ( PL )' * X.
//
    for ( j = 0; j < n - 1; j++ )
    {
      jhi = i4_min ( n - 1, j + ml );
      for ( i = j + 1; i <= jhi; i++ )
      {
        b[j] = b[j] - b[i] * a_lu[i-j+mu+j*nrow];
      }
    }
//
//  B = U' * Y = ( PL * U )' * X = A' * X.
//
    for ( i = n - 1; 0 <= i; i-- )
    {
      jhi = i4_min ( n - 1, i + mu );
      for ( j = i+1; j <= jhi; j++ )
      {
        b[j] = b[j] + b[i] * a_lu[i-j+mu+j*nrow];
      }
      b[i] = b[i] * a_lu[i-i+mu+i*nrow];
    }
  }

  return b;
}
//****************************************************************************80

double *r8cb_mxv ( int n, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_MXV multiplies a R8CB matrix times a vector.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(ML+MU+1)*N], the R8CB matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8CB_MXV[N], the product A * x.
//
{
  double *b;
  int i;
  int j;
  int jhi;
  int jlo;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    jlo = i4_max ( 0, i - ml );
    jhi = i4_min ( n-1, i + mu );
    for ( j = jlo; j <= jhi; j++ )
    {
      b[i] = b[i] + a[i-j+mu+j*(ml+mu+1)] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

int r8cb_np_fa ( int n, int ml, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_NP_FA factors a R8CB matrix by Gaussian elimination.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//    R8CB_NP_FA is a version of the LINPACK routine SGBFA, modifed to use
//    no pivoting, and to be applied to the R8CB compressed band matrix storage
//    format.  It will fail if the matrix is singular, or if any zero
//    pivot is encountered.
//
//    If R8CB_NP_FA successfully factors the matrix, R8CB_NP_SL may be called
//    to solve linear systems involving the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input/output, double A[(ML+MU+1)*N], the compact band matrix.
//    On input, the coefficient matrix of the linear system.
//    On output, the LU factors of the matrix.
//
//    Output, int R8CB_NP_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int i;
  int j;
  int ju;
  int k;
  int lm;
  int m;
  int mm;
//
//  The value of M is MU + 1 rather than ML + MU + 1.
//
  m = mu + 1;
  ju = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  If our pivot entry A(MU+1,K) is zero, then we must give up.
//
    if ( a[m-1+(k-1)*(ml+mu+1)] == 0.0 )
    {
      cerr << "\n";
      cerr << "R8CB_FA - Fatal error!\n";
      cerr << "  Zero pivot on step " << k << "\n";
      exit ( 1 );
    }
//
//  LM counts the number of nonzero elements that lie below the current
//  diagonal entry, A(K,K).
//
//  Multiply the LM entries below the diagonal by -1/A(K,K), turning
//  them into the appropriate "multiplier" terms in the L matrix.
//
    lm = i4_min ( ml, n-k );
    for ( i = m+1; i <= m+lm; i++ )
    {
      a[i-1+(k-1)*(ml+mu+1)] = -a[i-1+(k-1)*(ml+mu+1)] / a[m-1+(k-1)*(ml+mu+1)];
    }
//
//  MM points to the row in which the next entry of the K-th row is, A(K,J).
//  We then add L(I,K)*A(K,J) to A(I,J) for rows I = K+1 to K+LM.
//
    ju = i4_max ( ju, mu + k );
    ju = i4_min ( ju, n );
    mm = m;

    for ( j = k+1; j <= ju; j++ )
    {
      mm = mm - 1;
      for ( i = 1; i <= lm; i++ )
      {
        a[mm+i-1+(j-1)*(ml+mu+1)] = a[mm+i-1+(j-1)*(ml+mu+1)] 
          + a[mm-1+(j-1)*(ml+mu+1)] * a[m+i-1+(k-1)*(ml+mu+1)];
      }
    }
  }

  if ( a[m-1+(n-1)*(ml+mu+1)] == 0.0 )
  {
    cerr << "\n";
    cerr << "R8CB_FA - Fatal error!\n";
    cerr << "  Zero pivot on step " << n << "\n";
    exit ( 1 );
  }

  return 0;
}
//****************************************************************************80

double *r8cb_np_sl ( int n, int ml, int mu, double a_lu[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_NP_SL solves a R8CB system factored by R8CB_NP_FA.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//    R8CB_NP_SL can also solve the related system A' * x = b.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(ML+MU+1)*N], the LU factors from R8CB_NP_FA.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Input, int JOB.
//    If JOB is zero, the routine will solve A * x = b.
//    If JOB is nonzero, the routine will solve A' * x = b.
//
//    Output, double R8CB_NP_SL[N], the solution of the linear system, X.
//
{
  int i;
  int k;
  int la;
  int lb;
  int lm;
  int m;
  double *x;

  x = new double[n];

  m = mu + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }
//
//  Solve A * x = b.
//
  if ( job == 0 )
  {
//
//  Solve PL * Y = B.
//
    if ( 0 < ml )
    {
      for ( k = 1; k <= n-1; k++ )
      {
        lm = i4_min ( ml, n-k );
        for ( i = 0; i < lm; i++ )
        {
          x[k+i] = x[k+i] + x[k-1] * a_lu[m+i+(k-1)*(ml+mu+1)];
        }
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a_lu[m-1+(k-1)*(ml+mu+1)];
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[lb+i-1] = x[lb+i-1] - x[k-1] * a_lu[la+i-1+(k-1)*(ml+mu+1)];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
//
//  Solve U' * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[k-1] = x[k-1] - a_lu[la+i-1+(k-1)*(ml+mu+1)] * x[lb+i-1];
      }
      x[k-1] = x[k-1] / a_lu[m-1+(k-1)*(ml+mu+1)];

    }
//
//  Solve ( PL )' * X = Y.
//
    if ( 0 < ml )
    {
      for ( k = n-1; 1 <= k; k-- )
      {
        lm = i4_min ( ml, n-k );
        for ( i = 0; i < lm; i++ )
        {
          x[k-1] = x[k-1] + a_lu[m+i+(k-1)*(ml+mu+1)] * x[k+i];
        }
      }
    }
  }

  return x;
}
//****************************************************************************80

void r8cb_print ( int m, int n, int ml, int mu, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_PRINT prints a R8CB matrix.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1..
//
//    Input, double A[(ML+MU+1)*N], the R8CB matrix.
//
//    Input, string TITLE, a title.
//
{
  r8cb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8cb_print_some ( int m, int n, int ml, int mu, double a[], int ilo, 
  int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_PRINT_SOME prints some of a R8CB matrix.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(ML+MU+1)*N], the R8CB matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - mu );
    i2hi = i4_min ( ihi, m );
    i2hi = i4_min ( i2hi, j2hi + ml );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( ml < i-j || mu < j-i )
        {
          cout << "              ";
        }
        else
        {
          cout << setw(12) << a[i-j+mu+(j-1)*(ml+mu+1)] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8cb_random ( int n, int ml, int mu, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_RANDOM randomizes a R8CB matrix.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8CB_RANDOM[(ML+MU+1)*N], the R8CB matrix.
//
{
  double *a;
  int i;
  int ihi;
  int ilo;
  int j;

  a = new double[(ml+mu+1)*n];
//
//  Set the entries that correspond to matrix elements.
//
  for ( j = 0; j < n; j++ )
  {
    ilo = i4_max ( 0, j - mu );
    ihi = i4_min ( n-1, j + ml );

    for ( i = j - mu; i < 0; i++ )
    {
      a[i-j+mu+j*(ml+mu+1)] = 0.0;
    }
    for ( i = ilo; i <= ihi; i++ )
    {
      a[i-j+mu+j*(ml+mu+1)] = r8_uniform_01 ( seed );
    }
    for ( i = n; i <= j+ml; i++ )
    {
      a[i-j+mu*j*(ml+mu+1)] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8cb_to_r8vec ( int m, int n, int ml, int mu, double *a )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_TO_R8VEC copies a R8CB matrix to a real vector.
//
//  Discussion:
//
//    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
//    a data item carries its dimensionality implicitly, and so cannot be
//    regarded sometimes as a vector and sometimes as an array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//
//    Input, double A[(ML+MU+1)*N], the array to be copied.
//
//    Output, double R8CB_TO_R8VEC[(ML+MU+1)*N], the vector.
//
{
  int i;
  int ihi;
  int ilo;
  int j;
  double *x;

  x = new double[(ml+mu+1)*n];

  for ( j = 1; j <= n; j++ )
  {
    ihi = i4_min ( mu, mu + 1 - j );
    for ( i = 1; i <= ihi; i++ )
    {
      x[i-1+(j-1)*(ml+mu+1)] = 0.0;
    }

    ilo = i4_max ( ihi + 1, 1 );
    ihi = i4_min ( ml+mu+1, mu+1+m-j );
    for ( i = ilo; i <= ihi; i++ )
    {
      x[i-1+(j-1)*(ml+mu+1)] = a[i-1+(j-1)*(ml+mu+1)];
    }

    ilo = ihi + 1;
    ihi = ml+mu+1;
    for ( i = ilo; i <= ihi; i++ )
    {
      x[i-1+(j-1)*(ml+mu+1)] = 0.0;
    }

  }

  return x;
}
//****************************************************************************80

double *r8cb_to_r8ge ( int n, int ml, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_TO_R8GE copies a R8CB matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths of A.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(ML+MU+1)*N], the R8CB matrix.
//
//    Output, double R8CB_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( j-mu <= i && i <= j+ml )
      {
        b[i+j*n] = a[mu+i-j+j*(ml+mu+1)];
      }
      else
      {
        b[i+j*n] = 0.0;
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8cb_vxm ( int n, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_VXM multiplies a vector by a R8CB matrix.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(ML+MU+1)*N], the R8CB matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8CB_VXM[N], the product X*A.
//
{
  double *b;
  int i;
  int j;
  int jhi;
  int jlo;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    jlo = i4_max ( 0, i - ml );
    jhi = i4_min ( n-1, i + mu );
    for ( j = jlo; j <= jhi; j++ )
    {
      b[j] = b[j] + x[i] * a[i-j+mu+j*(ml+mu+1)];
    }
  }

  return b;
}
//****************************************************************************80

double *r8cb_zero ( int n, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8CB_ZERO zeros a R8CB matrix.
//
//  Discussion:
//
//    The R8CB storage format is appropriate for a compact banded matrix.
//    It is assumed that the matrix has lower and upper bandwidths ML and MU,
//    respectively.  The matrix is stored in a way similar to that used
//    by LINPACK and LAPACK for a general banded matrix, except that in
//    this mode, no extra rows are set aside for possible fillin during pivoting.
//    Thus, this storage format is suitable if you do not intend to factor
//    the matrix, or if you can guarantee that the matrix can be factored
//    without pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be nonnegative.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N-1.
//
//    Output, double R8CB_ZERO[(ML+MU+1)*N), the R8CB matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[(ml+mu+1)*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < ml+mu+1; i++ )
    {
      a[i+j*(ml+mu+1)] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

void r8cbb_add ( int n1, int n2, int ml, int mu, double a[], int i, int j, 
  double value )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_ADD adds a value to an entry of a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input/output, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
//
//    Input, int I, J, the indices of the entry to be incremented.
//
//    Input, double VALUE, the value to be added to the (I,J) entry.
//
{
  int ij;

  if ( value == 0.0 )
  {
    return;
  }
//
//  Check for I or J out of bounds.
//
  if ( i <= 0 || n1+n2 < i )
  {
    cerr << "\n";
    cerr << "R8CBB_ADD - Fatal error!\n";
    cerr << "  Illegal input value of row index I = " << i << "\n";
    exit ( 1 );
  }

  if ( j <= 0 || n1+n2 < j )
  {
    cerr << "\n";
    cerr << "R8CBB_ADD - Fatal error!\n";
    cerr << "  Illegal input value of column index J = " << j << "\n";
    exit ( 1 );
  }
//
//  The A1 block of the matrix.
//
//  Check for out of band problems.
//
  if ( i <= n1 && j <= n1 )
  {
    if ( mu < (j-i) || ml < (i-j) )
    {
      cout << "\n";
      cout << "R8CBB_ADD - Warning!\n";
      cout << "  Unable to add to entry (" << i << ", " << j << ").\n";
      return;
    }
    else
    {
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1);
    }
  }
//
//  The A2 block of the matrix:
//
  else if ( i <= n1 && n1 < j )
  {
    ij = (ml+mu+1)*n1+(j-n1-1)*n1 + i;
  }
//
//  The A3 and A4 blocks of the matrix.
//
  else if ( n1 < i )
  {
    ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1);
  }

  a[ij-1] = a[ij-1] + value;

  return;
}
//****************************************************************************80

bool r8cbb_error ( int n1, int n2, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_ERROR checks the dimensions of a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1 - 1.
//
//    Output, bool R8CBB_ERROR, is TRUE if an error was detected.
//
{
  if ( ml < 0 ) 
  {
    cout << "\n";
    cout << "R8CBB_ERROR - Illegal ML = " <<  ml << "\n";
    cout << "  but ML must be greater than or equal to 0.\n";
    return true;
  }

  if ( i4_max ( n1 - 1, 0 ) < ml )
  {
    cout << "\n";
    cout << "R8CBB_ERROR - Illegal ML = " << ml << "\n";
    cout << "  but ML must be <= Max ( N1 - 1, 0 ).\n";
    return true;
  }

  if ( mu < 0  )
  {
    cout << "\n";
    cout << "R8CBB_ERROR - Illegal MU = " << mu << "\n";
    cout << "  but MU must be greater than or equal to 0.\n";
    return true;
  }

  if ( i4_max ( n1 - 1, 0 ) < ml )
  {
    cout << "\n";
    cout << "R8CBB_ERROR - Illegal MU = " << mu << "\n";
    cout << "  but MU must be <= Max ( N1 - 1, 0 ).\n";
    return true;
  }

  if ( n1 < 0 )
  {
    cout << "\n";
    cout << "R8CBB_ERROR - Illegal N1 = " << n1 << "\n";
    return true;
  }

  if ( n2 < 0 )
  {
    cout << "\n";
    cout << "R8CBB_ERROR - Illegal N2 = " << n2 << "\n";
    return true;
  }

  if ( n1 + n2 <= 0 )
  {
    cout << "\n";
    cout << "R8CBB_ERROR - Illegal N1+N2 = " << n1 + n2 << "\n";
    return true;
  }

  return false;
}
//****************************************************************************80

int r8cbb_fa ( int n1, int n2, int ml, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_FA factors a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//
//    Once the matrix has been factored by SCCB_FA, SCCB_SL may be called
//    to solve linear systems involving the matrix.
//
//    SCCB_FA uses special non-pivoting versions of LINPACK routines to
//    carry out the factorization.  The special version of the banded
//    LINPACK solver also results in a space saving, since no entries
//    need be set aside for fill in due to pivoting.
//
//    The linear system must be border banded, of the form:
//
//      ( A1 A2 ) (X1) = (B1)
//      ( A3 A4 ) (X2)   (B2)
//
//    where A1 is a (usually big) banded square matrix, A2 and A3 are
//    column and row strips which may be nonzero, and A4 is a dense
//    square matrix.
//
//    The algorithm rewrites the system as:
//
//         X1 + inverse(A1) A2 X2 = inverse(A1) B1
//
//      A3 X1 +             A4 X2 = B2
//
//    and then rewrites the second equation as
//
//      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
//
//    The algorithm will certainly fail if the matrix A1 is singular,
//    or requires pivoting.  The algorithm will also fail if the A4 matrix,
//    as modified during the process, is singular, or requires pivoting.
//    All these possibilities are in addition to the failure that will
//    if the total matrix A is singular.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input/output, double A[ (ML+MU+1)*N1 + 2*N1*N2 + N2*N2].
//    On input, A contains the compact border-banded coefficient matrix.
//    On output, A contains information describing a partial factorization
//    of the original coefficient matrix.  
//
//    Output, int R8CBB_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  double *b1;
  int i;
  int ij;
  int ik;
  int info;
  int j;
  int jk;
  int job;
  int k;
  int nband;
  double *x1;

  nband = (ml+mu+1)*n1;
//
//  Factor the A1 band matrix, overwriting A1 by its factors.
//
  if ( 0 < n1 )
  {
    info = r8cb_np_fa ( n1, ml, mu, a );
    if ( info != 0 )
    {
      cerr << "\n";
      cerr << "R8CBB_FA - Fatal error!\n";
      cerr << "  R8CB_NP_FA returned INFO = " << info << "\n";
      cerr << "  Factoring failed for column INFO.\n";
      cerr << "  The band matrix A1 is singular.\n";
      cerr << "  This algorithm cannot continue!\n";
      exit ( 1 );
    }
  }

  if ( 0 < n1 && 0 < n2 )
  {
//
//  Set A2 := -inverse(A1) * A2.
//
    for ( j = 0; j < n2; j++ )
    {
      for ( i = 0; i < n1; i++ )
      {
        a[nband+i+j*n1] = -a[nband+i+j*n1];
      }
    }

    b1 = new double[n1];
    x1 = new double[n1];
    job = 0;

    for ( j = 0; j < n2; j++ )
    {
      for ( i = 0; i < n1; i++ )
      {
        b1[i] = a[nband+i+j*n1];
      }
      x1 = r8cb_np_sl ( n1, ml, mu, a, b1, job );
      for ( i = 0; i < n1; i++ )
      {
        a[nband+i+j*n1] = x1[i];
      }
    }
    delete [] b1;
    delete [] x1;
//
//  Set A4 := A4 + A3*A2
//
    for ( i = 1; i <= n2; i++ )
    {
      for ( j = 1; j <= n1; j++ )
      {
        ij = nband + n1*n2 + (j-1)*n2 + i - 1;
        for ( k = 1; k <= n2; k++ )
        {
          ik = nband + 2*n1*n2 + (k-1)*n2 + i - 1;
          jk = nband + (k-1)*n1 + j - 1;
          a[ik] = a[ik] + a[ij] * a[jk];
        }
      }
    }
  }
//
//  Factor A4.
//
  if ( 0 < n2 )
  {
    info = r8ge_np_fa ( n2, a+(nband+2*n1*n2) );

    if ( info != 0 )
    {
      cerr << "\n";
      cerr << "R8CBB_FA - Fatal error!\n";
      cerr << "  R8GE_NP_FA returned INFO = " << info << "\n";
      cerr << "  This indicates singularity in column " << n1+info << ".\n";
      cerr << "  The dense matrix A4 is singular.\n";
      cerr << "  This algorithm cannot continue!\n";
      exit ( 1 );
    }
  }

  return 0;
}
//****************************************************************************80

double r8cbb_get ( int n1, int n2, int ml, int mu, double a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_GET gets the value of an entry of a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input/output, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
//
//    Input, int I, J, the indices of the entry to be incremented.
//
//    Output, double R8CBB_GET, the value of the (I,J) entry.
//
{
  int ij;
//
//  Check for I or J out of bounds.
//
  if ( i <= 0 || n1+n2 < i )
  {
    return 0.0;
  }

  if ( j <= 0 || n1+n2 < j )
  {
    return 0.0;
  }
//
//  The A1 block of the matrix.
//
//  Check for out of band problems.
//
  if ( i <= n1 && j <= n1 )
  {
    if ( mu < (j-i) || ml < (i-j) )
    {
      return 0.0;
    }
    else
    {
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1);
    }
  }
//
//  The A2 block of the matrix:
//
  else if ( i <= n1 && n1 < j )
  {
    ij = (ml+mu+1)*n1+(j-n1-1)*n1 + i;
  }
//
//  The A3 and A4 blocks of the matrix.
//
  else if ( n1 < i )
  {
    ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1);
  }

  return a[ij-1];
}
//****************************************************************************80

double *r8cbb_indicator ( int n1, int n2, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_INDICATOR sets up a R8CBB indicator matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Output, double R8CBB_INDICATOR[(ML+MU+1)*N1+2*N1*N2+N2*N2], the R8CBB indicator matrix.
//
{
  double *a;
  int base;
  int fac;
  int i;
  int j;
  int row;

  a = new double[(ml+mu+1)*n1+2*n1*n2+n2*n2];

  fac = i4_power ( 10, i4_log_10 ( n1 + n2 ) + 1 );
//
//  Set the banded matrix A1.
//
  for ( j = 1; j <= n1; j++ )
  {
    for ( row = 1; row <= ml + mu + 1; row++ )
    {
      i = row + j - mu - 1;
      if ( 1 <= i && i <= n1 )
      {
        a[row-1+(j-1)*(ml+mu+1)] = ( double ) ( fac * i + j );
      }
      else
      {
        a[row-1+(j-1)*(ml+mu+1)] = 0.0;
      }
    }
  }
//
//  Set the N1 by N2 rectangular strip A2.
//
  base = ( ml + mu + 1 ) * n1;

  for ( i = 1; i <= n1; i++ )
  {
    for ( j = n1 + 1; j <= n1 + n2; j++ )
    {
      a[base + i-1 + (j-n1-1)*n1 ] = ( double ) ( fac * i + j );
    }
  }
//
//  Set the N2 by N1 rectangular strip A3.
//
  base = ( ml + mu + 1 ) * n1 + n1 * n2;

  for ( i = n1 + 1; i <= n1 + n2; i++ )
  {
    for ( j = 1; j <= n1; j++ )
    {
      a[base + i-n1-1 + (j-1)*n2 ] = ( double ) ( fac * i + j );
    }
  }
//
//  Set the N2 by N2 square A4.
//
  base = ( ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1;

  for ( i = n1 + 1; i <= n1 + n2; i++ )
  {
    for ( j = n1 + 1; j <= n1 + n2; j++ )
    {
      a[base + i-n1-1 + (j-n1-1)*n2 ] = ( double ) ( fac * i + j );
    }
  }

  return a;
}
//****************************************************************************80

double *r8cbb_mxv ( int n1, int n2, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_MXV multiplies a R8CBB matrix times a vector.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
//
//    Input, double X[N1+N2], the vector to be multiplied by A.
//
//    Output, double R8CBB_MXV[N1+N2], the result of multiplying A by X.
//
{
  double *b;
  int i;
  int ihi;
  int ij;
  int ilo;
  int j;
//
//  Set B to zero.
//
  b = new double[n1+n2];

  for ( i = 0; i < n1+n2; i++ )
  {
    b[i] = 0.0;
  }
//
//  Multiply by A1.
//
  for ( j = 1; j <= n1; j++ )
  {
    ilo = i4_max ( 1, j-mu );
    ihi = i4_min ( n1, j+ml );
    ij = (j-1)*(ml+mu+1)-j+mu+1;
    for ( i = ilo; i <= ihi; i++ )
    {
      b[i-1] = b[i-1] + a[ij+i-1] * x[j-1];
    }
  }
//
//  Multiply by A2.
//
  for ( j = n1+1; j <= n1+n2; j++ )
  {
    ij = (ml+mu+1)*n1+(j-n1-1)*n1;
    for ( i = 1; i <= n1; i++ )
    {
      b[i-1] = b[i-1] + a[ij+i-1] * x[j-1];
    }
  }
//
//  Multiply by A3 and A4.
//
  for ( j = 1; j <= n1+n2; j++ )
  {
    ij = (ml+mu+1)*n1+n1*n2+(j-1)*n2-n1;
    for ( i = n1+1; i <= n1+n2; i++ )
    {
      b[i-1] = b[i-1] + a[ij+i-1] * x[j-1];
    }
  }

  return b;
}
//****************************************************************************80

void r8cbb_print ( int n1, int n2, int ml, int mu, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_PRINT prints a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, double A[(ML+MU+1)*N1+2*N1*N2+N2*N2], the R8CBB matrix.
//
//    Input, string TITLE, a title.
//
{
  r8cbb_print_some ( n1, n2, ml, mu, a, 1, 1, n1+n2, n1+n2, title );

  return;
}
//****************************************************************************80

void r8cbb_print_some ( int n1, int n2, int ml, int mu, double a[], int ilo, 
  int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_PRINT_SOME prints some of a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, double A[(ML+MU+1)*N1+2*N1*N2+N2*N2], the R8CBB matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int ij;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n1+n2 );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n1+n2 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
    cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        aij = 0.0;

        if ( i <= n1 && j <= n1 )
        {
          if ( j - i <= mu && i - j <= ml )
          {
            ij = (i-j+mu+1)+(j-1)*(ml+mu+1);
            aij = a[ij-1];
          }
        }
        else if ( i <= n1 && n1 < j )
        {
          ij = (ml+mu+1)*n1+(j-n1-1)*n1+i;
          aij = a[ij-1];
        }
        else if ( n1 < i )
        {
          ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1);
          aij = a[ij-1];
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8cbb_random ( int n1, int n2, int ml, int mu, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_RANDOM randomizes a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than N1-1.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8CBB_RANDOM[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
//
{
  double *a;
  int i;
  int j;
  double r;
  int row;

  a = new double[(ml+mu+1)*n1+2*n1*n2+n2*n2];
//
//  Randomize the banded matrix A1.
//  We still believe that the "junk" entries should be set to 0.
//
  for ( j = 1; j <= n1; j++ )
  {
    for ( row = 1; row <= ml+mu+1; row++ )
    {
      i = row + j - mu - 1;
      if ( 1 <= i && i <= n1 )
      {
        r = r8_uniform_01 ( seed );
      }
      else
      {
        r = 0.0;
      }
      a[row-1+(j-1)*(ml+mu+1)] = r;
    }
  }
//
//  Randomize the rectangular strips A2+A3+A4.
//
  for ( i = (ml+mu+1)*n1+1; i <= (ml+mu+1)*n1+2*n1*n2+n2*n2; i++ )
  {
    a[i-1] = r8_uniform_01 ( seed );
  }

  return a;
}
//****************************************************************************80

void r8cbb_set ( int n1, int n2, int ml, int mu, double a[], int i, int j, 
  double value )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_SET sets an entry of a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input/output, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
//
//    Input, int I, J, the indices of the entry to be incremented.
//
//    Input, double VALUE, the value to be assigned to the (I,J) entry.
//
{
  int ij;
//
//  Check for I or J out of bounds.
//
  if ( i <= 0 || n1+n2 < i )
  {
    cerr << "\n";
    cerr << "R8CBB_SET - Fatal error!\n";
    cerr << "  Illegal input value of row index I = " << i << "\n";
    exit ( 1 );
  }

  if ( j <= 0 || n1+n2 < j )
  {
    cerr << "\n";
    cerr << "R8CBB_SET - Fatal error!\n";
    cerr << "  Illegal input value of column index J = " << j << "\n";
    exit ( 1 );
  }
//
//  The A1 block of the matrix.
//
//  Check for out of band problems.
//
  if ( i <= n1 && j <= n1 )
  {
    if ( mu < (j-i) || ml < (i-j) )
    {
      cout << "\n";
      cout << "R8CBB_SET - Warning!\n";
      cout << "  Unable to set entry (" << i << ", " << j << ").\n";
      return;
    }
    else
    {
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1);
    }
  }
//
//  The A2 block of the matrix:
//
  else if ( i <= n1 && n1 < j )
  {
    ij = (ml+mu+1)*n1+(j-n1-1)*n1 + i;
  }
//
//  The A3 and A4 blocks of the matrix.
//
  else if ( n1 < i )
  {
    ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1);
  }

  a[ij-1] = value;

  return;
}
//****************************************************************************80

double *r8cbb_sl ( int n1, int n2, int ml, int mu, double a_lu[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_SL solves a R8CBB system factored by R8CBB_FA.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//
//    The linear system A * x = b is decomposable into the block system:
//
//      ( A1 A2 ) * (X1) = (B1)
//      ( A3 A4 )   (X2)   (B2)
//
//    where A1 is a (usually big) banded square matrix, A2 and A3 are
//    column and row strips which may be nonzero, and A4 is a dense
//    square matrix.
//
//    All the arguments except B are input quantities only, which are
//    not changed by the routine.  They should have exactly the same values
//    they had on exit from R8CBB_FA.
//
//    If more than one right hand side is to be solved, with the same
//    matrix, R8CBB_SL should be called repeatedly.  However, R8CBB_FA only
//    needs to be called once to create the factorization.
//
//    See the documentation of R8CBB_FA for details on the matrix storage.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, double A_LU[ (ML+MU+1)*N1 + 2*N1*N2 + N2*N2].
//    the LU factors from R8CBB_FA.
//
//    Input, double B[N1+N2], the right hand side of the linear system.
//
//    Output, double R8CBB_SL[N1+N2], the solution.
//
{
  double *b2;
  int i;
  int ij;
  int j;
  int job;
  int nband;
  double *x;
  double *x1;
  double *x2;

  nband = (ml+mu+1)*n1;
//
//  Set X1 := inverse(A1) * B1.
//
  if ( 0 < n1 )
  {
    job = 0;
    x1 = r8cb_np_sl ( n1, ml, mu, a_lu, b, job );
  }
//
//  Modify the right hand side of the second linear subsystem.
//  Set B2 = B2-A3*X1.
//
  b2 = new double[n2];

  for ( i = 0; i < n2; i++ )
  {
    ij = nband + n1*n2 + j*n2 + i;
    b2[i] = b[n1+i];
  }

  for ( j = 0; j < n1; j++ )
  {
    for ( i = 0; i < n2; i++ )
    {
      ij = nband + n1*n2 + j*n2 + i;
      b2[i] = b2[i] - a_lu[ij] * x1[j];
    }
  }
//
//  Solve A4*X2 = B2.
//
  if ( 0 < n2 )
  {
    job = 0;
    x2 = r8ge_np_sl ( n2, a_lu+(nband+2*n1*n2), b2, job );
  }
//
//  Modify the first subsolution.
//  Set X1 = X1+A2*X2.
//
  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n2; j++ )
    {
      ij = nband + j*n1 + i;
      x1[i] = x1[i] + a_lu[ij] * x2[j];
    }
  }
//
//  Collect X1 and X2 into X.
//
  x = new double[n1+n2];

  for ( i = 0; i < n1; i++ )
  {
    x[i] = x1[i];
  }
  for ( i = 0; i < n2; i++ )
  {
    x[n1+i] = x2[i];
  }

  delete [] b2;
  delete [] x1;
  delete [] x2;

  return x;
}
//****************************************************************************80

double *r8cbb_to_r8ge ( int n1, int n2, int ml, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_TO_R8GE copies a R8CBB matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, double A[(ML+MU+1)*N1+2*N1*N2+N2*N2], the R8CBB matrix.
//
//    Output, double R8CBB_TO_R8GE[(N1+N2)*(N1+N2)], the R8GE matrix.
//
{
  double *b;
  int i;
  int ij;
  int j;

  b = new double[(n1+n2)*(n1+n2)];

  for ( i = 1; i <= n1; i++ )
  {
    for ( j = 1; j <= n1; j++ )
    {
      if ( mu+ml < (j-i) || ml < (i-j) )
      {
        b[i-1+(j-1)*(n1+n2)] = 0.0;
      }
      else
      {
        ij = (i-j+mu+1)+(j-1)*(ml+mu+1);
        b[i-1+(j-1)*(n1+n2)] = a[ij-1];
      }
    }
  }

  for ( i = 1; i <= n1; i++ )
  {
    for ( j = n1+1; j <= n2; j++ )
    {
      ij = (ml+mu+1)*n1+(j-n1-1)*n1+i;
      b[i-1+(j-1)*(n1+n2)] = a[ij-1];
    }
  }

  for ( i = n1+1; i <= n2; i++ )
  {
    for ( j = 1; j <= n1+n2; j++)
    {
      ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1);
      b[i-1+(j-1)*(n1+n2)] = a[ij-1];
    }
  }

  return b;
}
//****************************************************************************80

double *r8cbb_vxm ( int n1, int n2, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_VXM multiplies a vector by a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
//
//    Input, double X[N1+N2], the vector to multiply the matrix.
//
//    Output, double R8CBB_VXM[N1+N2], the product X * A.
//
{
  double *b;
  int i;
  int ihi;
  int ij;
  int ilo;
  int j;
//
//  Set B to zero.
//
  b = new double[n1+n2];
  for ( i = 0; i < n1+n2; i++ )
  {
    b[i] = 0.0;
  }
//
//  Multiply by A1.
//
  for ( j = 1; j <= n1; j++ )
  {
    ilo = i4_max ( 1, j-mu );
    ihi = i4_min ( n1, j+ml );
    ij = (j-1)*(ml+mu+1)-j+mu+1;
    for ( i = ilo; i <= ihi; i++ )
    {
      b[j] = b[j] + x[i-1] * a[ij+i-1];
    }
  }
//
//  Multiply by A2.
//
  for ( j = n1+1; j <= n1+n2; j++ )
  {
    ij = (ml+mu+1)*n1+(j-n1-1)*n1;
    for ( i = 1; i <= n1; i++ )
    {
      b[j] = b[j] + x[i-1] * a[ij+i-1];
    }
  }
//
//  Multiply by A3 and A4.
//
  for ( j = 1; j <= n1+n2; j++ )
  {
    ij = (ml+mu+1)*n1+n1*n2+(j-1)*n2-n1;
    for ( i = n1+1; i <= n1+n2; i++ )
    {
      b[j-1] = b[j-1] + x[i-1] * a[ij+i-1];
    }
  }

  return b;
}
//****************************************************************************80

double *r8cbb_zero ( int n1, int n2, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8CBB_ZERO zeros a R8CBB matrix.
//
//  Discussion:
//
//    The R8CBB storage format is for a compressed border banded matrix.  
//    Such a matrix has the logical form:
//
//      A1 | A2
//      ---+---
//      A3 | A4
//
//    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
//    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
//    respectively.  
//
//    The R8CBB format is the same as the R8BB format, except that the banded
//    matrix A1 is stored in compressed band form rather than standard
//    banded form.  In other words, we do not include the extra room
//    set aside for fill in during pivoting.
//
//    A should be defined as a vector.  The user must then store
//    the entries of the four blocks of the matrix into the vector A.
//    Each block is stored by columns.
//
//    A1, the banded portion of the matrix, is stored in
//    the first (ML+MU+1)*N1 entries of A, using the obvious variant
//    of the LINPACK general band format.
//
//    The following formulas should be used to determine how to store
//    the entry corresponding to row I and column J in the original matrix:
//
//    Entries of A1:
//
//      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
//
//      Store the I, J entry into location
//      (I-J+MU+1)+(J-1)*(ML+MU+1).
//
//    Entries of A2:
//
//      1 <= I <= N1, N1+1 <= J <= N1+N2.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+(J-N1-1)*N1+I.
//
//    Entries of A3:
//
//      N1+1 <= I <= N1+N2, 1 <= J <= N1.
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//
//    Entries of A4:
//
//      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
//
//      Store the I, J entry into location
//      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
//      (same formula used for A3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the banded and dense blocks.
//    N1 and N2 must be nonnegative, and at least one must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N1-1.
//
//    Output, double R8CBB_ZERO[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
//
{
  double *a;
  int i;

  a = new double[(ml+mu+1)*n1+2*n1*n2+n2*n2];

  for ( i = 0; i < (ml+mu+1)*n1+2*n1*n2+n2*n2; i++ )
  {
    a[i] = 0.0;
  }

  return a;
}
//****************************************************************************80

double r8cc_get ( int m, int n, int nz_num, int col[], int row[], 
  double a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_GET gets a value of a R8CC matrix.
//
//  Discussion:
//
//    It is legal to request entries of the matrix for which no storage
//    was set aside.  In that case, a zero value will be returned.
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero entries.
//
//    Input, int COL[N+1], indicate where each column's data begins.
//
//    Input, int ROW[NZ_NUM], the row indices.
//
//    Input, double A[NZ_NUM], the nonzero entries.
//
//    Input, int I, J, the indices of the value to retrieve.
//
//    Output, double R8CC_GET, the value of A(I,J).
//
{
  double aij;
  int k;
//
//  Seek sparse index K corresponding to full index (I,J).
//
  k = r8cc_ijk ( m, n, nz_num, col, row, i, j );
//
//  If no K was found, then be merciful, and simply return 0.
//
  if ( k == -1 )
  {
    aij = 0.0;
  }
  else
  {
    aij = a[k-1];
  }

  return aij;
}
//****************************************************************************80

int r8cc_ijk ( int m, int n, int nz_num, int col[], int row[], int i, 
  int j )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_IJK seeks K, the sparse index of (I,J), the full index of a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero entries.
//
//    Input, int COL[N+1], indicate where each column's data begins.
//
//    Input, int ROW[NZ_NUM], the row indices.
//
//    Input, int I, J, the indices of the value to retrieve.
//
//    Output, int R8CC_IJK, the index of the sparse matrix in which entry
//    (I,J) is stored, or -1 if no such entry exists.
//
{
  int k;
  int k1;
  int k2;
//
//  Determine the part of ROW containing row indices of entries
//  in column J.
//
  k1 = col[j-1];
  k2 = col[j]-1;
//
//  Seek the location K for which ROW(K) = I.
//  
  k = i4vec_search_binary_a ( k2+1-k1, row+k1-1, i );

  if ( k != -1 )
  {
    k = k + k1 - 1;
  }

  return k;
}
//****************************************************************************80

void r8cc_inc ( int m, int n, int nz_num, int col[], int row[], double a[], 
  int i, int j, double aij )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_INC increments a value of a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero entries.
//
//    Input, int COL[N+1], indicate where each column's data begins.
//
//    Input, int ROW[NZ_NUM], the row indices.
//
//    Input/output, double A[NZ_NUM], the nonzero entries.
//    On output, entry (I,J) has been incremented.
//
//    Input, int I, J, the indices of the value to retrieve.
//
//    Input, double AIJ, the value to be added to A(I,J).
//
{
  int k;
//
//  Seek sparse index K corresponding to full index (I,J).
//
  k = r8cc_ijk ( m, n, nz_num, col, row, i, j );
//
//  If no K was found, we fail.
//
  if ( k == -1 )
  {
    cerr << "\n";
    cerr << "R8CC_INC - Fatal error!\n";
    cerr << "  R8CC_IJK could not find the entry.\n";
    cerr << "  Row I = " << i << "\n";
    cerr << "  Col J = " << j << "\n";
    exit ( 1 );
  }
  a[k-1] = a[k-1] + aij;

  return;
}
//****************************************************************************80

double *r8cc_indicator ( int m, int n, int nz_num, int col[], int row[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_INDICATOR sets up a R8CC indicator matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in A.
//
//    Input, int COL[N+1], points to the first element of each column.
//
//    Input, int ROW[NZ_NUM], contains the row indices of the elements.
//
//    Output, double R8CC_INDICATOR[NZ_NUM], the R8CC matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;
  int k;

  a = new double[nz_num];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( j = 1; j <= n; j++ )
  {
    for ( k = col[j-1]; k <= col[j] - 1; k++ )
    {
      i = row[k-1];
      a[k-1] = ( double ) ( fac * i + j );
    }
  }

  return a;
}
//****************************************************************************80

void r8cc_kij ( int m, int n, int nz_num, int col[], int row[], int k, 
  int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_KIJ seeks (I,J), the full index of K, the sparse index of a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero entries.
//
//    Input, int COL[N+1], indicate where each column's data begins.
//
//    Input, int ROW[NZ_NUM], the row indices.
//
//    Input, int K, the sparse index of an entry of the matrix.
//    1 <= K <= NZ_NUM.
//
//    Output, int *I, *J, the full indices corresponding to the sparse
//    index K.
//
{
  int jj;
  int k1;
  int k2;

  *i = -1;
  *j = -1;

  if ( k < 1 || nz_num < k )
  {
    return;
  }
//
//  The row index is easy.
//
  *i = row[k-1];
//
//  Determine the column by bracketing in COl.
//
  for ( jj = 1; jj <= n; jj++ )
  {
    k1 = col[jj-1];
    k2 = col[jj]-1;
    if ( k1 <= k && k <= k2 )
    {
      *j = jj;
      break;
    }
  }

  if ( *j == -1 )
  {
    return;
  }
  return;
}
//****************************************************************************80

double *r8cc_mxv ( int m, int n, int nz_num, int col[], int row[], 
  double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_MXV multiplies a R8CC matrix times a vector.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in A.
//
//    Input, int COL[N+1], points to the first element of each column.
//
//    Input, int ROW[NZ_NUM], contains the row indices of the elements.
//
//    Input, double A[NZ_NUM], the R8CC matrix.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Output, double R8CC_MXV[M], the product A * X.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = col[j]; k <= col[j+1] - 1; k++ )
    {
      i = row[k-1] - 1;
      b[i] = b[i] + a[k-1] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void r8cc_print ( int m, int n, int nz_num, int col[], int row[], 
  double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_PRINT prints a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in A.
//
//    Input, int COL[N+1], points to the first element of each column.
//
//    Input, int ROW[NZ_NUM], contains the row indices of the elements.
//
//    Input, double A[NZ_NUM], the R8CC matrix.
//
//    Input, string TITLE, a title.
//
{
  r8cc_print_some ( m, n, nz_num, col, row, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8cc_print_some ( int m, int n, int nz_num, int col[], int row[], 
  double a[], int ilo, int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_PRINT_SOME prints some of a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in A.
//
//    Input, int COL[N+1], points to the first element of each column.
//
//    Input, int ROW[NZ_NUM], contains the row indices of the elements.
//
//    Input, double A[NZ_NUM], the R8CC matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;
  int k;
  double value;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col:  ";

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";

    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
//
//  Now consider each column J in J2LO to J2HI,
//  and look at every nonzero, and check if it occurs in row I.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        value = 0.0;
        for ( k = col[j-1]; k <= col[j]-1; k++ )
        {
          if ( row[k-1] == i )
          {
            value = a[k-1];
          }
        }
        cout << setw(12) << value << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8cc_random ( int m, int n, int nz_num, int col[], int row[], 
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_RANDOM randomizes a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in A.
//
//    Input, int COL[N+1], points to the first element of each column.
//
//    Input, int ROW[NZ_NUM], contains the row indices of the elements.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8CC_RANDOM[NZ_NUM], the R8CC matrix.
//
{
  double *a;
  int j;
  int k;

  a = new double[nz_num];

  for ( j = 0; j < n; j++ )
  {
    for ( k = col[j]; k <= col[j+1] - 1; k++ )
    {
      a[k-1] = r8_uniform_01 ( seed );
    }
  }

  return a;
}
//****************************************************************************80

void r8cc_read ( string col_file, string row_file, string a_file, int m, 
  int n, int nz_num, int col[], int row[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_READ reads a R8CC matrix from three files.
//
//  Discussion:
//
//    This routine needs the values of M, N, and NZ_NUM, which can be 
//    determined by a call to R8CC_READ_SIZE.
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, string COL_FILE, ROW_FILE, A_FILE, the names of the 
//    files containing the column pointers, row indices, and matrix entries.
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Output, int COL[N+1], the column pointers.
//
//    Output, int ROW[NZ_NUM], the row indices.
//
//    Output, double A[NZ_NUM], the nonzero elements of the matrix.
//
{
  ifstream input;
  int k;

  input.open ( col_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8CC_READ - Fatal error!\n";
    cerr << "  Could not open the file \"" << col_file << "\".\n";
    exit ( 1 );
  }

  for ( k = 0; k < n+1; k++ )
  {
    input >> col[k];
  }

  input.close ( );
//
//  Read the row information.
//
  input.open ( row_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8CC_READ - Fatal error!\n";
    cerr << "  Could not open the file \"" << row_file << "\".\n";
    exit ( 1 );
  }

  for ( k = 0; k < nz_num; k++ )
  {
    input >> row[k];
  }

  input.close ( );
//
//  Read the value information.
//
  input.open ( a_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8CC_READ - Fatal error!\n";
    cerr << "  Could not open the file \"" << a_file << "\".\n";
    exit ( 1 );
  }

  for ( k = 0; k < nz_num; k++ )
  {
    input >> a[k];
  }

  input.close ( );

  return;
}
//****************************************************************************80

void r8cc_read_size ( string col_file, string row_file, int *m, int *n, 
  int *nz_num, int *base )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_READ_SIZE reads the sizes of a R8CC sparse matrix from a file.
//
//  Discussion:
//
//    The value of M is "guessed" to be the largest value that occurs in
//    the ROW file.  However, if a row index of 0 is encountered, then
//    the value of M is incremented by 1.
//
//    The value of N is the number of records in the COL file minus 1.
//
//    The value of NZ_NUM is simply the number of records in the ROW file.
//
//    The value of BASE is 0 or 1, depending on whether the program
//    "guesses" that the row and column indices are 0-based or 1-based.
//    Although the first entry of the COL array might be used as evidence,
//    this program makes its determination based on whether it encounters
//    a 0 index in the ROW file.
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, string COL_FILE, *ROW_FILE, the names of the 
//    column and row files that describe the structure of the matrix.
//
//    Output, int *M, *N, the inferred number of rows and columns 
//    in the sparse matrix.
//
//    Output, int *NZ_NUM, the number of nonzero entries in the
//    sparse matrix.
//
//    Output, int *BASE, is 0 if the row indexing is believed
//    to be 0-based, and 1 if the row-index is believed to be
//    1-based.  In uncertain cases, BASE = 1 is the default.
//
{
  int col;
  ifstream input;
  ifstream input2;
  int row;
//
//  Default values.
//
  *m = -1;
  *n = -1;
  *nz_num = -1;
  *base = -1;
//
//  Check the COL file first.
//
  input.open ( col_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8CC_READ_SIZE - Fatal error!\n";
    cerr << "  Could not open the file \"" << col_file << "\".\n";
    exit ( 1 );
  }

  *n = -1;

  for ( ; ; )
  {
    input >> col;

    if ( input.eof( ) )
    {
      break;
    }
    *n = *n + 1;
  }
  input.close ( );
//
//  Check the ROW file.
//
//  For unfathomable reasons, if I use "INPUT" for this file,
//  I can get a file open failure.  Rather than make right the
//  world, I gave up and accessed "INPUT2".
//
  input2.open ( row_file.c_str ( ) );

  if ( !input2 )
  {
    cerr << "\n";
    cerr << "R8CC_READ_SIZE - Fatal error!\n";
    cerr << "  Could not open the file \"" << row_file << "\".\n";
    exit ( 1 );
  }

  *base = 1;
  *m = 0;
  *nz_num = 0;
  
  for ( ; ; )
  {
    input2 >> row;

    if ( input2.eof ( ) )
    {
      break;
    }
    *nz_num = *nz_num + 1;
    *m = i4_max ( *m, row );
    if ( row == 0 )
    {
      *base = 0;
    }
  }
  input2.close ( );

  return;
}
//****************************************************************************80

void r8cc_set ( int m, int n, int nz_num, int col[], int row[], double a[], 
  int i, int j, double aij )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_SET sets a value of a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero entries.
//
//    Input, int COL[N+1], indicate where each column's data begins.
//
//    Input, int ROW[NZ_NUM], the row indices.
//
//    Input/output, double A[NZ_NUM], the nonzero entries.
//    On output, the entry of A corresponding to (I,J) has been reset.
//
//    Input, int I, J, the indices of the value to retrieve.
//
//    Input, double AIJ, the new value of A(I,J).
//
{
  int k;
//
//  Seek sparse index K corresponding to full index (I,J).
//
  k = r8cc_ijk ( m, n, nz_num, col, row, i, j );
//
//  If no K was found, we fail.
//
  if ( k == -1 )
  {
    cerr << "\n";
    cerr << "R8CC_SET - Fatal error!\n";
    cerr << "  R8CC_IJK could not find the entry.\n";
    cerr << "  Row I = " << i << "\n";
    cerr << "  Col J = " << j << "\n";
    exit ( 1 );
  }
  a[k-1] = aij;

  return;
}
//****************************************************************************80

double *r8cc_to_r8ge ( int m, int n, int nz_num, int col[], int row[], 
  double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_TO_R8GE converts a R8CC matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in A.
//
//    Input, int COL[N+1], points to the first element of each column.
//
//    Input, int ROW[NZ_NUM], contains the row indices of the elements.
//
//    Input, double A[NZ_NUM], the R8CC matrix.
//
//    Input, double R8CC_TO_R8GE[M*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+j*m] = 0.0;
    }
  }

  if ( col[0] < 0 || nz_num < col[0] )
  {
    cerr << "\n";
    cerr << "R8CC_TO_R8GE - Fatal error!\n";
    cerr << "  COL[" << j << "] = " << col[j] << "\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    if ( col[j+1] < 0 || nz_num < col[j+1] - 1 )
    {
      cerr << "\n";
      cerr << "R8CC_TO_R8GE - Fatal error!\n";
      cerr << "  COL[" << j+1 << "] = " << col[j+1] << "\n";
      exit ( 1 );
    }

    for ( k = col[j]; k <= col[j+1] - 1; k++ )
    {
      i = row[k-1] - 1;
      if ( i < 0 || m <= i )
      {
        cerr << "\n";
        cerr << "R8CC_TO_R8GE - Fatal error!\n";
        cerr << "  ROW[" << k-1 << "] = " << i << "\n";
        exit ( 1 );
      }
      b[i+j*m] = a[k-1];
    }
  }

  return b;
}
//****************************************************************************80

double *r8cc_vxm ( int m, int n, int nz_num, int col[], int row[], 
  double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_VXM multiplies a vector times a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in A.
//
//    Input, int COL[N+1], points to the first element of each column.
//
//    Input, int ROW[NZ_NUM], contains the row indices of the elements.
//
//    Input, double A[NZ_NUM], the R8CC matrix.
//
//    Input, double X[M], the vector to be multiplied.
//
//    Output, double R8CC_VXM[N], the product A' * X.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n];

  for ( j = 0; j < n; j++ )
  {
    b[j] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = col[j]; k <= col[j+1] - 1; k++ )
    {
      i = row[k-1] - 1;
      b[j] = b[j] + a[k-1] * x[i];
    }
  }

  return b;
}
//****************************************************************************80

void r8cc_write ( string col_file, string row_file, string a_file, int m, int n,
  int nz_num, int col[], int row[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_WRITE writes a R8CC matrix to three files.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, string COL_FILE, ROW_FILE, A_FILE, the names of the 
//    files containing the column pointers, row entries, and matrix entries.
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int COL[N+1], the column pointers.
//
//    Input, int ROW[NZ_NUM], the row indices.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
{
  ofstream output;
  int k;

  output.open ( col_file.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8CC_WRITE - Fatal error!\n";
    cerr << "  Could not open the file \"" << col_file << "\".\n";
    exit ( 1 );
  }

  for ( k = 0; k < n+1; k++ )
  {
    output << col[k] << "\n";
  }

  output.close ( );
//
//  Write the row information.
//
  output.open ( row_file.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8CC_WRITE - Fatal error!\n";
    cerr << "  Could not open the file \"" << row_file << "\".\n";
    exit ( 1 );
  }

  for ( k = 0; k < nz_num; k++ )
  {
    output << row[k] << "\n";
  }

  output.close ( );
//
//  Write the value information.
//
  output.open ( a_file.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8CC_WRITE - Fatal error!\n";
    cerr << "  Could not open the file \"" << a_file << "\".\n";
    exit ( 1 );
  }

  for ( k = 0; k < nz_num; k++ )
  {
    output << a[k] << "\n";
  }
  output.close ( );

  return;
}
//****************************************************************************80

double *r8cc_zero ( int m, int n, int nz_num, int col[], int row[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CC_ZERO zeros a R8CC matrix.
//
//  Discussion:
//
//    The R8CC format is the double precision sparse compressed column
//    format.  Associated with this format, we have an M by N matrix
//    with NZ_NUM nonzero entries.  We construct the column pointer
//    vector COL of length N+1, such that entries of column J will be
//    stored in positions COL(J) through COL(J+1)-1.  This indexing
//    refers to both the ROW and A vectors, which store the row indices
//    and the values of the nonzero entries.  The entries of the
//    ROW vector corresponding to each column are assumed to be
//    ascending sorted.
//
//    The R8CC format is equivalent to the MATLAB "sparse" format,
//    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in A.
//
//    Input, int COL[N+1], points to the first element of each column.
//
//    Input, int ROW[NZ_NUM], contains the row indices of the elements.
//
//    Output, double R8CC_ZERO[NZ_NUM], the R8CC matrix.
//
{
  double *a;
  int j;
  int k;

  a = new double[nz_num];

  for ( j = 0; j < n; j++ )
  {
    for ( k = col[j]; k <= col[j+1] - 1; k++ )
    {
      a[k-1] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8ci_eval ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_EVAL returns the eigenvalues of a R8CI matrix.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis,
//    Circulant Matrices,
//    Wiley, 1979.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N], the R8CI matrix.
//
//    Output, double R8CI_EVAL[2*N], the complex eigenvalues.
//
{
  int i;
  int j;
  double *lambda;
  double li;
  double lr;
  double *w;

  lambda = new double[2*n];

  w = c8vec_unity ( n );

  for ( i = 0; i < n; i++ )
  {
    lambda[0+i*2] = a[n-1];
    lambda[1+i*2] = 0.0;
  }

  for ( i = n-2; 0 <= i; i-- )
  {
    for ( j = 0; j < n; j++ )
    {
      lr = lambda[0+j*2] * w[0+j*2] 
         - lambda[1+j*2] * w[1+j*2] + a[i];

      li = lambda[0+j*2] * w[1+j*2]
         + lambda[1+j*2] * w[0+j*2];

      lambda[0+j*2] = lr;

      lambda[1+j*2] = li;
    }
  }

  c8vec_sort_a2 ( n, lambda );

  delete [] w;

  return lambda;
}
//****************************************************************************80

double *r8ci_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_INDICATOR sets up a R8CI indicator matrix.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//    The R8CI format simply records the first row of the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Output, double R8CI_INDICATOR[N], the R8CI matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  i = 1;

  for ( j = 1; j <= n; j++ )
  {
    a[j-1] = ( double ) ( fac * i + j );
  }

  return a;
}
//****************************************************************************80

double *r8ci_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_MXV multiplies a R8CI matrix times a vector.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N], the R8CI matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8CI_MXV[N], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j <= i-1; j++ )
    {
      b[i] = b[i] + a[j-i+n] * x[j];
    }
    for ( j = i; j < n; j++ )
    {
      b[i] = b[i] + a[j-i] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void r8ci_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_PRINT prints a R8CI matrix.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N], the R8CI matrix.
//
//    Input, string TITLE, a title.
//
{
  r8ci_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8ci_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_PRINT_SOME prints some of a R8CI matrix.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N], the R8CI matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(6) << i << "  ";

      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j )
        {
          cout << setw(12) << a[j-i] << "  ";
        }
        else
        {
          cout << setw(12) << a[n+j-i] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8ci_random ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_RANDOM randomizes a R8CI matrix.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8CI_RANDOM[N], the R8CI matrix.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = r8_uniform_01 ( seed );
  }

  return a;
}
//****************************************************************************80

double *r8ci_sl ( int n, double a[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_SL solves a R8CI system.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N], the R8CI matrix.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double R8CI_SL[N], the solution of the linear system.
//
{
  int i;
  int nsub;
  double r1;
  double r2;
  double r3;
  double r5;
  double r6;
  double *work;
  double *x;

  work = new double[2*n-2];
  x = new double[n];

  if ( job == 0 )
  {
//
//  Solve the system with the principal minor of order 1.
//
    r1 = a[0];
    x[0] = b[0] / r1;

    r2 = 0.0;
//
//  Recurrent process for solving the system.
//
    for ( nsub = 2; nsub <= n; nsub++ )
    {
//
//  Compute multiples of the first and last columns of
//  the inverse of the principal minor of order N.
//
      r5 = a[n+2-nsub-1];
      r6 = a[nsub-1];

      if ( 2 < nsub )
      {
        work[nsub-2] = r2;

        for ( i = 1; i <= nsub - 2; i++ )
        {
          r5 = r5 + a[n-i] * work[nsub-i-1];
          r6 = r6 + a[i] * work[n-2+i];
        }
      }

      r2 = - r5 / r1;
      r3 = - r6 / r1;
      r1 = r1 + r5 * r3;

      if ( 2 < nsub )
      {
        r6 = work[n-1];
        work[n+nsub-3] = 0.0;
        for ( i = 2; i <= nsub - 1; i++ )
        {
          r5 = work[n-2+i];
          work[n-2+i] = work[i-1] * r3 + r6;
          work[i-1] = work[i-1] + r6 * r2;
          r6 = r5;
        }
      }

      work[n-1] = r3;
//
//  Compute the solution of the system with the principal minor of order NSUB.
//
      r5 = 0.0;
      for ( i = 1; i <= nsub - 1; i++ )
      {
        r5 = r5 + a[n-i] * x[nsub-i-1];
      }

      r6 = ( b[nsub-1] - r5 ) / r1;
      for ( i = 1; i <= nsub-1; i++ )
      {
        x[i-1] = x[i-1] + work[n+i-2] * r6;
      }
      x[nsub-1] = r6;
    }
  }
  else
  {
//
//  Solve the system with the principal minor of order 1.
//
    r1 = a[0];
    x[0] = b[0] / r1;

    r2 = 0.0;
//
//  Recurrent process for solving the system.
//
    for ( nsub = 2; nsub <= n; nsub++ )
    {
//
//  Compute multiples of the first and last columns of
//  the inverse of the principal minor of order N.
//
      r5 = a[nsub-1];
      r6 = a[n+1-nsub];

      if ( 2 < nsub )
      {
        work[nsub-2] = r2;
        for ( i = 1; i <= nsub - 2; i++ )
        {
          r5 = r5 + a[i] * work[nsub-i-1];
          r6 = r6 + a[n-i] * work[n-2+i];
        }
      }

      r2 = - r5 / r1;
      r3 = - r6 / r1;
      r1 = r1 + r5 * r3;

      if ( 2 < nsub )
      {
        r6 = work[n-1];
        work[n+nsub-3] = 0.0;
        for ( i = 2; i <= nsub-1; i++ )
        {
          r5 = work[n-2+i];
          work[n-2+i] = work[i-1] * r3 + r6;
          work[i-1] = work[i-1] + r6 * r2;
          r6 = r5;
        }
      }

      work[n-1] = r3;
//
//  Compute the solution of the system with the principal minor of order NSUB.
//
      r5 = 0.0;
      for ( i = 1; i <= nsub - 1; i++ )
      {
        r5 = r5 + a[i] * x[nsub-i-1];
      }

      r6 = ( b[nsub-1] - r5 ) / r1;
      for ( i = 1; i <= nsub - 1; i++ )
      {
        x[i-1] = x[i-1] + work[n-2+i] * r6;
      }

      x[nsub-1] = r6;
    }
  }

  delete [] work;

  return x;
}
//****************************************************************************80

double *r8ci_to_r8ge ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_TO_R8GE copies a R8CI matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N], the R8CI matrix.
//
//    Output, double R8CI_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 1; i < n; i++ )
    {
      k = i4_modp (  j - i, n );
      b[i+j*n] = a[k];
    }
  }

  return b;
}
//****************************************************************************80

double *r8ci_vxm ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_VXM multiplies a vector times a R8CI matrix.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N], the R8CI matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8CI_VXM[N], the product A' * X.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j <= i; j++ )
    {
      b[i] = b[i] + a[i-j] * x[j];
    }
    for ( j = i+1; j < n; j++ )
    {
      b[i] = b[i] + a[n+i-j] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

double *r8ci_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8CI_ZERO zeros a R8CI matrix.
//
//  Discussion:
//
//    The R8CI storage format is used for an N by N circulant matrix.
//    An N by N circulant matrix A has the property that the entries on
//    row I appear again on row I+1, shifted one position to the right,
//    with the final entry of row I appearing as the first of row I+1.
//
//    A circulant matrix data structure simply records the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Output, double R8CI_ZERO[N], the R8CI matrix.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }

  return a;
}
//****************************************************************************80

double r8gb_det ( int n, int ml, int mu, double a_lu[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_DET computes the determinant of a matrix factored by R8GB_FA or R8GB_TRF.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_FA or R8GB_TRF.
//
//    Input, int PIVOT[N], the pivot vector, as computed by R8GB_FA
//    or R8GB_TRF.
//
//    Output, double R8GB_DET, the determinant of the matrix.
//
{
  int col = 2 * ml + mu + 1;
  double det;
  int i;

  det = 1.0;

  for ( i = 0; i < n; i++ )
  {
    det = det * a_lu[ml+mu+i*col];
  }

  for ( i = 0; i < n; i++ )
  {
    if ( pivot[i] != i+1 ) 
    {
      det = -det;
    }
  }

  return det;
}
//****************************************************************************80

int r8gb_fa ( int n, int ml, int mu, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_FA performs a LINPACK-style PLU factorization of a R8GB matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input/output, double A[(2*ML+MU+1)*N], the matrix in band storage.  
//    On output, A has been overwritten by the LU factors.
//
//    Output, int PIVOT[N], the pivot vector.
//
//    Output, int R8GB_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int col = 2 * ml + mu + 1;
  int i;
  int i0;
  int j;
  int j0;
  int j1;
  int ju;
  int jz;
  int k;
  int l;
  int lm;
  int m;
  int mm;
  double t;

  m = ml + mu + 1;
//
//  Zero out the initial fill-in columns.
//
  j0 = mu + 2;
  j1 = i4_min ( n, m ) - 1;

  for ( jz = j0; jz <= j1; jz++ )
  {
    i0 = m + 1 - jz;
    for ( i = i0; i <= ml; i++ )
    {
      a[i-1+(jz-1)*col] = 0.0;
    }
  }

  jz = j1;
  ju = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Zero out the next fill-in column.
//
    jz = jz + 1;
    if ( jz <= n ) 
    {
      for ( i = 1; i <= ml; i++ )
      {
        a[i-1+(jz-1)*col] = 0.0;
      }
    }
//
//  Find L = pivot index.
//
    lm = i4_min ( ml, n-k );
    l = m;

    for ( j = m+1; j <= m + lm; j++ )
    {
      if ( r8_abs ( a[l-1+(k-1)*col] ) < r8_abs ( a[j-1+(k-1)*col] ) )
      {
        l = j;
      }
    }

    pivot[k-1] = l + k - m;
//
//  Zero pivot implies this column already triangularized.
//
    if ( a[l-1+(k-1)*col] == 0.0 )
    {
      cerr << "\n";
      cerr << "R8GB_FA - Fatal error!\n";
      cerr << "  Zero pivot on step " << k << "\n";
      exit ( 1 );
    }
//
//  Interchange if necessary.
//
    t                = a[l-1+(k-1)*col];
    a[l-1+(k-1)*col] = a[m-1+(k-1)*col];
    a[m-1+(k-1)*col] = t;
//
//  Compute multipliers.
//
    for ( i = m+1; i <= m+lm; i++ )
    {
      a[i-1+(k-1)*col] = - a[i-1+(k-1)*col] / a[m-1+(k-1)*col];
    }
//
//  Row elimination with column indexing.
//
    ju = i4_max ( ju, mu + pivot[k-1] );
    ju = i4_min ( ju, n );
    mm = m;

    for ( j = k+1; j <= ju; j++ )
    {
      l = l - 1;
      mm = mm - 1;

      if ( l != mm )
      {
        t                 = a[l-1+(j-1)*col];
        a[l-1+(j-1)*col]  = a[mm-1+(j-1)*col];
        a[mm-1+(j-1)*col] = t;
      }
      for ( i = 1; i <= lm; i++ )
      {
        a[mm+i-1+(j-1)*col] = a[mm+i-1+(j-1)*col] 
          + a[mm-1+(j-1)*col] * a[m+i-1+(k-1)*col];
      }
    }
  }

  pivot[n-1] = n;

  if ( a[m-1+(n-1)*col] == 0.0 )
  {
    cerr << "\n";
    cerr << "R8GB_FA - Fatal error!\n";
    cerr << "  Zero pivot on step " << n << "\n";
    exit ( 1 );
  }

  return 0;
}
//****************************************************************************80

double *r8gb_indicator ( int m, int n, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_INDICATOR sets up a R8GB indicator matrix.
//
//  Discussion:
//
//    Note that the R8GB storage format includes extra room for
//    fillin entries that occur during Gauss elimination.  This routine
//    will leave those values at 0.
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2005
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Output, double R8GB_INDICATOR[(2*ML+MU+1)*N], the R8GB matrix.
//
{
  double *a;
  int col = 2 * ml + mu + 1;
  int diag;
  int fac;
  int i;
  int j;
  int k;

  a = new double[(2*ml+mu+1)*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );
  k = 0;

  for ( j = 1; j <= n; j++ )
  {
    for ( diag = 1; diag <= 2 * ml + mu + 1; diag++ )
    {
      i = diag + j - ml - mu - 1;

      if ( 1 <= i && i <= m && i - ml <= j && j <= i + mu )
      {
        a[diag-1+(j-1)*col] = ( double ) ( fac * i + j );
      }
      else if ( 1 <= i && i <= m && i - ml <= j && j <= i + mu + ml )
      {
        a[diag-1+(j-1)*col] = 0.0;
      }
      else
      {
        k = k + 1;
        a[diag-1+(j-1)*col] = - ( double ) k;
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8gb_ml ( int n, int ml, int mu, double a_lu[], int pivot[], double x[], 
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_ML computes A * x or A' * X, using R8GB_FA factors.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//    It is assumed that R8GB_FA has overwritten the original matrix
//    information by LU factors.  R8GB_ML is able to reconstruct the
//    original matrix from the LU factor data.
//
//    R8GB_ML allows the user to check that the solution of a linear
//    system is correct, without having to save an unfactored copy
//    of the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_FA.
//
//    Input, int PIVOT[N], the pivot vector computed by R8GB_FA.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Input, int JOB, specifies the operation to be done:
//    JOB = 0, compute A * x.
//    JOB nonzero, compute A' * X.
//
//    Output, double R8GB_ML[N], the result of the multiplication.
//
{
  double *b;
  int col = 2 * ml + mu + 1;
  int i;
  int ihi;
  int ilo;
  int j;
  int jhi;
  int k;
  double temp;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
//
//  Y = U * X.
//
    for ( j = 1; j <= n; j++ )
    {
      ilo = i4_max ( 1, j - ml - mu );
      for ( i = ilo; i <= j-1; i++ )
      {
        b[i-1] = b[i-1] + a_lu[i-j+ml+mu+(j-1)*col] * b[j-1];
      }
      b[j-1] = a_lu[j-j+ml+mu+(j-1)*col] * b[j-1];
    }
//
//  B = PL * Y = PL * U * X = A * x.
//
    for ( j = n-1; 1 <= j; j-- )
    {
      ihi = i4_min ( n, j + ml );
      for ( i = j+1; i <= ihi; i++ )
      {
        b[i-1] = b[i-1] - a_lu[i-j+ml+mu+(j-1)*col] * b[j-1];
      }

      k = pivot[j-1];

      if ( k != j )
      {
        temp   = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = temp;
      }
    }
  }
  else
  {
//
//  Y = ( PL )' * X.
//
    for ( j = 1; j <= n-1; j++ )
    {
      k = pivot[j-1];

      if ( k != j )
      {
        temp   = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = temp;
      }

      jhi = i4_min ( n, j + ml );
      for ( i = j+1; i <= jhi; i++ )
      {
        b[j-1] = b[j-1] - b[i-1] * a_lu[i-j+ml+mu+(j-1)*col];
      }

    }
//
//  B = U' * Y = ( PL * U )' * X = A' * X.
//
    for ( i = n; 1 <= i; i-- )
    {
      jhi = i4_min ( n, i + ml + mu );
      for ( j = i+1; j <= jhi; j++ )
      {
        b[j-1] = b[j-1] + b[i-1] * a_lu[i-j+ml+mu+(j-1)*col];
      }
      b[i-1] = b[i-1] * a_lu[i-i+ml+mu+(i-1)*col];
    }

  }

  return b;
}
//****************************************************************************80

double *r8gb_mu ( int n, int ml, int mu, double a_lu[], int pivot[], double x[], 
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_MU computes A * x or A' * X, using R8GB_TRF factors.
//
//  Warning:
//
//    This routine must be updated to allow for rectangular matrices.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    It is assumed that R8GB_TRF has overwritten the original matrix
//    information by LU factors.  R8GB_MU is able to reconstruct the
//    original matrix from the LU factor data.
//
//    R8GB_MU allows the user to check that the solution of a linear
//    system is correct, without having to save an unfactored copy
//    of the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_TRF.
//
//    Input, int PIVOT[N], the pivot vector computed by R8GB_TRF.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Input, int JOB, specifies the operation to be done:
//    JOB = 0, compute A * x.
//    JOB nonzero, compute A' * X.
//
//    Output, double R8GB_MU[N], the result of the multiplication.
//
{
  double *b;
  int i;
  int ihi;
  int ilo;
  int j;
  int jhi;
  int k;
  double t;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
//
//  Y = U * X.
//
    for ( j = 1; j <= n; j++ )
    {
      ilo = i4_max ( 1, j - ml - mu );
      for ( i = ilo; i <= j - 1; i++ )
      {
        b[i-1] = b[i-1] + a_lu[i-j+ml+mu+(j-1)*(2*ml+mu+1)] * b[j-1];
      }
      b[j-1] = a_lu[j-j+ml+mu+(j-1)*(2*ml+mu+1)] * b[j-1];
    }
//
//  B = PL * Y = PL * U * X = A * x.
//
    for ( j = n-1; 1 <= j; j-- )
    {
      ihi = i4_min ( n, j + ml );
      for ( i = j+1; i <= ihi; i++ )
      {
        b[i-1] = b[i-1] + a_lu[i-j+ml+mu+(j-1)*(2*ml+mu+1)] * b[j-1];
      }

      k = pivot[j-1];

      if ( k != j )
      {
        t      = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = t;
      }
    }
  }
  else
  {
//
//  Y = ( PL )' * X.
//
    for ( j = 1; j <= n-1; j++ )
    {
      k = pivot[j-1];

      if ( k != j )
      {
        t      = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = t;
      }

      jhi = i4_min ( n, j + ml );
      for ( i = j+1; i <= jhi; i++ )
      {
        b[j-1] = b[j-1] + b[i-1] * a_lu[i-j+ml+mu+(j-1)*(2*ml+mu+1)];
      }
    }
//
//  B = U' * Y = ( PL * U )' * X = A' * X.
//
    for ( i = n; 1 <= i; i-- )
    {
      jhi = i4_min ( n, i + ml + mu );
      for ( j = i+1; j <= jhi; j++ )
      {
        b[j-1] = b[j-1] + b[i-1] * a_lu[i-j+ml+mu+(j-1)*(2*ml+mu+1)];
      }
      b[i-1] = b[i-1] * a_lu[i-i+ml+mu+(i-1)*(2*ml+mu+1)];
    }
  }

  return b;
}
//****************************************************************************80

double *r8gb_mxv ( int m, int n, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_MXV multiplies a R8GB matrix times a vector.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//    LINPACK and LAPACK storage of general band matrices requires
//    an extra ML upper diagonals for possible fill in entries during
//    Gauss elimination.  This routine does not access any entries
//    in the fill in diagonals, because it assumes that the matrix
//    has NOT had Gauss elimination applied to it.  If the matrix
//    has been Gauss eliminated, then the routine R8GB_MU must be
//    used instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8GB_MXV[M], the product A * x.
//
{
  double *b;
  int col = 2 * ml + mu + 1;
  int i;
  int j;
  int jhi;
  int jlo;

  b = new double[m];

  for ( i = 1; i <= m; i++ )
  {
    b[i-1] = 0.0;
    jlo = i4_max ( 1, i - ml );
    jhi = i4_min ( n, i + mu );
    for ( j = jlo; j <= jhi; j++ )
    {
      b[i-1] = b[i-1] + a[i-j+ml+mu+(j-1)*col] * x[j-1];
    }
  }

  return b;
}
//****************************************************************************80

int r8gb_nz_num ( int m, int n, int ml, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_NZ_NUM counts the nonzeroes in a R8GB matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    LINPACK and LAPACK band storage requires that an extra ML
//    superdiagonals be supplied to allow for fillin during Gauss
//    elimination.  Even though a band matrix is described as
//    having an upper bandwidth of MU, it effectively has an
//    upper bandwidth of MU+ML.  This routine will examine
//    values it finds in these extra bands, so that both unfactored
//    and factored matrices can be handled.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2003
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
//
//    Output, int R8GB_NZ_NUM, the number of nonzero entries in A.
//
{
  int i;
  int j;
  int jhi;
  int jlo;
  int nz_num;

  nz_num = 0;

  for ( i = 0; i < m; i++ )
  {
    jlo = i4_max ( 0, i - ml );
    jhi = i4_min ( n-1, i + mu + ml );
    for ( j = jlo; j <= jhi; j++ )
    {
      if ( a[i-j+ml+mu+j*(2*ml+mu+1)] != 0.0 )
      {
        nz_num = nz_num + 1;
      }
    }
  }

  return nz_num;
}
//****************************************************************************80

void r8gb_print ( int m, int n, int ml, int mu, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_PRINT prints a R8GB matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1..
//
//    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
//
//    Input, string TITLE, a title.
//
{
  r8gb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8gb_print_some ( int m, int n, int ml, int mu, double a[], int ilo, 
  int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_PRINT_SOME prints some of a R8GB matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1..
//
//    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int col = 2 * ml + mu + 1;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - mu - ml );

    i2hi = i4_min ( ihi, m );
    i2hi = i4_min ( i2hi, j2hi + ml );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(6) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i < j - mu - ml || j + ml < i )
        {
          cout << "            ";
        }
        else
        {
          cout << setw(10) << a[i-j+ml+mu+(j-1)*col] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8gb_random ( int m, int n, int ml, int mu, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_RANDOM randomizes a R8GB matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//    LINPACK and LAPACK band storage requires that an extra ML
//    superdiagonals be supplied to allow for fillin during Gauss
//    elimination.  Even though a band matrix is described as
//    having an upper bandwidth of MU, it effectively has an
//    upper bandwidth of MU+ML.  This routine assumes it is setting
//    up an unfactored matrix, so it only uses the first MU upper bands,
//    and does not place nonzero values in the fillin bands.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8GB_RANDOM[(2*ML+MU+1)*N], the R8GB matrix.
//
{
  double *a;
  int col = 2 * ml + mu + 1;
  int i;
  int j;
  int row;

  a = new double[col*n];

  for ( j = 1; j <= n; j++ )
  {
    for ( row = 1; row <= col; row++ )
    {
      i = row + j - ml - mu - 1;
      if ( ml < row && 1 <= i && i <= m )
      {
        a[row-1+(j-1)*col] = r8_uniform_01 ( seed );
      }
      else
      {
        a[(row-1)+(j-1)*col] = 0.0;
      }
    }
  }
  return a;
}
//****************************************************************************80

double *r8gb_sl ( int n, int ml, int mu, double a_lu[], int pivot[], 
  double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_SL solves a system factored by R8GB_FA.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_FA.
//
//    Input, int PIVOT[N], the pivot vector from R8GB_FA.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, int JOB.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double R8GB_SL[N], the solution.
//
{
  int col = 2 * ml + mu + 1;
  int i;
  int k;
  int l;
  int la;
  int lb;
  int lm;
  int m;
  double t;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  m = mu + ml + 1;
//
//  Solve A * x = b.
//
  if ( job == 0 )
  {
//
//  Solve L * Y = B.
//
    if ( 1 <= ml )
    {
      for ( k = 1; k <= n-1; k++ )
      {
        lm = i4_min ( ml, n-k );
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
        for ( i = 1; i <= lm; i++ )
        {
          x[k+i-1] = x[k+i-1] + x[k-1] * a_lu[m+i-1+(k-1)*col];
        }
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a_lu[m-1+(k-1)*col];
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[lb+i-1] = x[lb+i-1] - x[k-1] * a_lu[la+i-1+(k-1)*col];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
//
//  Solve U' * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[k-1] = x[k-1] - x[lb+i-1] * a_lu[la+i-1+(k-1)*col];
      }
      x[k-1] = x[k-1] / a_lu[m-1+(k-1)*col];
    }
//
//  Solve L' * X = Y.
//
    if ( 1 <= ml )
    {
      for ( k = n-1; 1 <= k; k-- )
      {
        lm = i4_min ( ml, n-k );
        for ( i = 1; i <= lm; i++ )
        {
          x[k-1] = x[k-1] + x[k+i-1] * a_lu[m+i-1+(k-1)*col];
        }
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
      }
    }
  }

  return x;
}
//****************************************************************************80

void r8gb_to_r8s3 ( int m, int n, int ml, int mu, double a[], int nz_num, 
  int *isym, int row[], int col[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_TO_R8S3 copies a R8GB matrix to a R8S3 matrix.
//
//  Discussion:
//
//    The R8GB storage format is for an M by N banded matrix, with lower 
//    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
//    extra superdiagonals, which may be required to store nonzero entries 
//    generated during Gaussian elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    LINPACK and LAPACK band storage requires that an extra ML
//    superdiagonals be supplied to allow for fillin during Gauss
//    elimination.  Even though a band matrix is described as
//    having an upper bandwidth of MU, it effectively has an
//    upper bandwidth of MU+ML.  This routine will copy nonzero
//    values it finds in these extra bands, so that both unfactored
//    and factored matrices can be handled.
//
//    The R8S3 storage format corresponds to the SLAP Triad format.
//
//    The R8S3 storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.  The entries may be given in any order.  No
//    check is made for the erroneous case in which a given matrix entry is
//    specified more than once.
//
//    There is a symmetry option for square matrices.  If the symmetric storage
//    option is used, the format specifies that only nonzeroes on the diagonal
//    and lower triangle are stored.  However, this routine makes no attempt
//    to enforce this.  The only thing it does is to "reflect" any nonzero
//    offdiagonal value.  Moreover, no check is made for the erroneous case
//    in which both A(I,J) and A(J,I) are specified, but with different values.
//
//    This routine reorders the entries of A so that the first N entries
//    are exactly the diagonal entries of the matrix, in order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrices.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrices.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths of A1.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
//
//    Input, int NZ_NUM, the number of nonzero entries in A.
//    This number can be obtained by calling R8GB_NZ_NUM.
//
//    Output, int *ISYM, is 0 if the matrix is not symmetric, and 1
//    if the matrix is symmetric.  If the matrix is symmetric, then
//    only the nonzeroes on the diagonal and in the lower triangle are stored.
//    For this routine, ISYM is always output 0.
//
//    Output, int ROW[NZ_NUM]), the row indices.
//
//    Output, int COL[NZ_NUM], the column indices.
//
//    Output, double B[NZ_NUM], the R8S3 matrix.
//
{
  int i;
  int j;
  int nz;

  *isym = 0;
  nz = 0;

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( i - ml <= j && j <= i + mu + ml )
      {
        if ( a[ml+mu+1+i-j-1+(j-1)*(2*ml+mu+1)] != 0.0 )
        {
          if ( nz_num <= nz )
          {
            cerr << "\n";
            cerr << "R8GB_TO_R8S3 - Fatal error!\n";
            cerr << "  NZ_NUM = " << nz_num << "\n";
            cerr << "  But the matrix has more nonzeros than that!\n";
            exit ( 1 );
          }

          row[nz] = i;
          col[nz] = j;
          b[nz] = a[ml+mu+1+i-j-1+(j-1)*(2*ml+mu+1)];
          nz = nz + 1;
        }
      }
    }
  }

  if ( nz < nz_num )
  {
    cout << "\n";
    cout << "R8GB_TO_R8S3 - Warning!\n";
    cout << "  NZ_NUM = " << nz_num << "\n";
    cout << "  But the number of nonzeros is " << nz << "\n";
  }

  return;
}
//****************************************************************************80

void r8gb_to_r8sp ( int m, int n, int ml, int mu, double a[], int nz_num, 
  int row[], int col[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_TO_R8SP copies a R8GB matrix to a R8SP matrix.
//
//  Discussion:
//
//    The R8GB storage format is for an M by N banded matrix, with lower 
//    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
//    extra superdiagonals, which may be required to store nonzero entries 
//    generated during Gaussian elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    LINPACK and LAPACK band storage requires that an extra ML
//    superdiagonals be supplied to allow for fillin during Gauss
//    elimination.  Even though a band matrix is described as
//    having an upper bandwidth of MU, it effectively has an
//    upper bandwidth of MU+ML.  This routine will copy nonzero
//    values it finds in these extra bands, so that both unfactored
//    and factored matrices can be handled.
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrices.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrices.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths of A1.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
//
//    Input, int NZ_NUM, the number of nonzero entries in A.
//    This number can be obtained by calling R8GB_NZ_NUM.
//
//    Output, int ROW[NZ_NUM]), the row indices.
//
//    Output, int COL[NZ_NUM], the column indices.
//
//    Output, double B[NZ_NUM], the R8S3 matrix.
//
{
  int i;
  int j;
  int jhi;
  int jlo;
  int nz;

  nz = 0;

  for ( i = 1; i <= m; i++ )
  {
    jlo = i4_max ( 1, i - ml );
    jhi = i4_min ( n, i + ml + mu );

    for ( j = jlo; j <= jhi; j++ )
    {
      if ( a[ml+mu+1+i-j-1+(j-1)*(2*ml+mu+1)] == 0.0 )
      {
        continue;
      }
      if ( nz_num <= nz )
      {
        cerr << "\n";
        cerr << "R8GB_TO_R8SP - Fatal error!\n";
        cerr << "  NZ_NUM = " << nz_num << "\n";
        cerr << "  But the matrix has more nonzeros than that!\n";
        exit ( 1 );
      }

      row[nz] = i;
      col[nz] = j;
      b[nz] = a[ml+mu+1+i-j-1+(j-1)*(2*ml+mu+1)];
      nz = nz + 1;
    }
  }

  if ( nz < nz_num )
  {
    cout << "\n";
    cout << "R8GB_TO_R8SP - Warning!\n";
    cout << "  NZ_NUM = " << nz_num << "\n";
    cout << "  But the number of nonzeros is " << nz << "\n";
  }

  return;
}
//****************************************************************************80

double *r8gb_to_r8vec ( int m, int n, int ml, int mu, double *a )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_TO_R8VEC copies a R8GB matrix to a real vector.
//
//  Discussion:
//
//    In C++  and FORTRAN, this routine is not really needed.  In MATLAB,
//    a data item carries its dimensionality implicitly, and so cannot be
//    regarded sometimes as a vector and sometimes as an array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//
//    Input, double A[(2*ML+MU+1)*N], the array.
//
//    Output, double R8GB_TO_R8VEC[(2*ML+MU+1)*N], the vector.
//
{
  int i;
  int ihi;
  int ilo;
  int j;
  double *x;

  x = new double[(2*ml+mu+1)*n];

  for ( j = 1; j <= n; j++ )
  {
    ihi = i4_min ( ml + mu, ml + mu + 1 - j );
    for ( i = 1; i <= ihi; i++ )
    {
      x[i-1+(j-1)*(2*ml+mu+1)] = 0.0;
    }

    ilo = i4_max ( ihi + 1, 1 );
    ihi = i4_min ( 2*ml+mu+1, ml+mu+m+1-j );
    for ( i = ilo; i <= ihi; i++ )
    {
      x[i-1+(j-1)*(2*ml+mu+1)] = a[i-1+(j-1)*(2*ml+mu+1)];
    }

    ilo = ihi + 1;
    ihi = 2*ml+mu+1;
    for ( i = ilo; i <= ihi; i++ )
    {
      x[i-1+(j-1)*(2*ml+mu+1)] = 0.0;
    }

  }

  return x;
}
//****************************************************************************80

double *r8gb_to_r8ge ( int m, int n, int ml, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_TO_R8GE copies a R8GB matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//    LINPACK and LAPACK band storage requires that an extra ML
//    superdiagonals be supplied to allow for fillin during Gauss
//    elimination.  Even though a band matrix is described as
//    having an upper bandwidth of MU, it effectively has an
//    upper bandwidth of MU+ML.  This routine will copy nonzero
//    values it finds in these extra bands, so that both unfactored
//    and factored matrices can be handled.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrices.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrices.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths of A1.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
//
//    Output, double R8GB_TO_R8GE[M*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[m*n];

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( i - ml <= j && j <= i + mu )
      {
        b[i-1+(j-1)*m] = a[ml+mu+i-j+(j-1)*(2*ml+mu+1)];
      }
      else
      {
        b[i-1+(j-1)*m] = 0.0;
      }
    }
  }

  return b;
}
//****************************************************************************80

int r8gb_trf ( int m, int n, int ml, int mu, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_TRF performs a LAPACK-style PLU factorization of a R8GB matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//    This is a simplified, standalone version of the LAPACK
//    routine SGBTRF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 November 2003
//
//  Author:
//
//    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
//    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix A.  0 <= M.
//
//    Input, int N, the number of columns of the matrix A.  0 <= N.
//
//    Input, int ML, the number of subdiagonals within the band of A.
//    0 <= ML.
//
//    Input, int MU, the number of superdiagonals within the band of A.
//    0 <= MU.
//
//    Input/output, double A[(2*ML+MU+1)*N].  On input, the matrix A in band 
//    storage, and on output, information about the PLU factorization.
//
//    Output, int PIVOT(min(M,N)), the pivot indices;
//    for 1 <= i <= min(M,N), row i of the matrix was interchanged with
//    row IPIV(i).
//
//    Output, int R8GB_TRF, error flag.
//    = 0: successful exit;
//    < 0: an input argument was illegal;
//    > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
//         has been completed, but the factor U is exactly
//         singular, and division by zero will occur if it is used
//         to solve a system of equations.
//
{
  int i;
  int info;
  int j;
  int jp;
  int ju;
  int k;
  int km;
  int kv;
  double piv;
  double temp;

  info = 0;
//
//  KV is the number of superdiagonals in the factor U, allowing for fill-in.
//
  kv = mu + ml;
//
//  Set fill-in elements in columns MU+2 to KV to zero.
//
  for ( j = mu+2; j <= i4_min ( kv, n ); j++ )
  {
    for ( i = kv - j + 1; i <= ml; i++ )
    {
      a[i-1+(j-1)*(2*ml+mu+1)] = 0.0;
    }
  }
//
//  JU is the index of the last column affected by the current stage
//  of the factorization.
//
  ju = 1;

  for ( j = 1; j <= i4_min ( m, n ); j++ )
  {
//
//  Set the fill-in elements in column J+KV to zero.
//
    if ( j + kv <= n )
    {
      for ( i = 1; i <= ml; i++ )
      {
        a[i-1+(j+kv-1)*(2*ml+mu+1)] = 0.0;
      }
    }
//
//  Find the pivot and test for singularity.
//  KM is the number of subdiagonal elements in the current column.
//
    km = i4_min ( ml, m-j );

    piv = r8_abs ( a[kv+(j-1)*(2*ml+mu+1)] );
    jp = kv + 1;

    for ( i = kv + 2; i <= kv + km + 1; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(j-1)*(2*ml+mu+1)] ) )
      {
        piv = r8_abs ( a[i-1+(j-1)*(2*ml+mu+1)] );
        jp = i;
      }
    }

    jp = jp - kv;

    pivot[j-1] = jp + j - 1;

    if ( a[kv+jp-1+(j-1)*(2*ml+mu+1)] != 0.0 )
    {
      ju = i4_max ( ju, i4_min ( j+mu+jp-1, n ) );
//
//  Apply interchange to columns J to JU.
//
      if ( jp != 1 )
      {
        for ( i = 0; i <= ju - j; i++ )
        {
          temp = a[kv+jp-i-1+(j+i-1)*(2*ml+mu+1)];
          a[kv+jp-i-1+(j+i-1)*(2*ml+mu+1)] = a[kv+1-i-1+(j+i-1)*(2*ml+mu+1)];
          a[kv+1-i-1+(j+i-1)*(2*ml+mu+1)] = temp; 
        }
      }
//
//  Compute the multipliers.
//
      if ( 0 < km )
      {
        for ( i = kv+2; i <= kv+km+1; i++ )
        {
          a[i-1+(j-1)*(2*ml+mu+1)] = a[i-1+(j-1)*(2*ml+mu+1)] 
            / a[kv+(j-1)*(2*ml+mu+1)];
        }
//
//  Update the trailing submatrix within the band.
//
        if ( j < ju )
        {
          for ( k = 1; k <= ju - j; k++ )
          {
            if ( a[kv-k+(j+k-1)*(2*ml+mu+1)] != 0.0 )
            {
              for ( i = 1; i <= km; i++ )
              {
                a[kv+i-k+(j+k-1)*(2*ml+mu+1)] = a[kv+i-k+(j+k-1)*(2*ml+mu+1)] 
                  - a[kv+i+(j-1)*(2*ml+mu+1)] * a[kv-k+(j+k-1)*(2*ml+mu+1)];
              }
            }
          }
        }
      }
    }
    else
//
//  If pivot is zero, set INFO to the index of the pivot
//  unless a zero pivot has already been found.
//
    {
      if ( info == 0 )
      {
        info = j;
      }
    }
  }

  return info;
}
//****************************************************************************80

double *r8gb_trs ( int n, int ml, int mu, int nrhs, char trans, double a[], 
  int pivot[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_TRS solves a linear system factored by R8GB_TRF.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2004
//
//  Author:
//
//    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
//    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//    N must be positive.
//
//    Input, int ML, the number of subdiagonals within the band of A.
//    ML must be at least 0, and no greater than N - 1.
//
//    Input, int MU, the number of superdiagonals within the band of A.
//    MU must be at least 0, and no greater than N - 1.
//
//    Input, int NRHS, the number of right hand sides and the number of
//    columns of the matrix B.  NRHS must be positive.
//
//    Input, char TRANS, specifies the form of the system.
//    'N':  A * x = b  (No transpose)
//    'T':  A'* X = B  (Transpose)
//    'C':  A'* X = B  (Conjugate transpose = Transpose)
//
//    Input, double A[(2*ML+MU+1)*N], the LU factorization of the band matrix
//    A, computed by R8GB_TRF.  
//
//    Input, int PIVOT[N], the pivot indices; for 1 <= I <= N, row I
//    of the matrix was interchanged with row PIVOT(I).
//
//    Input, double B[N*NRHS], the right hand side vectors.
//
//    Output, double R8GB_TRS[N*NRHS], the solution vectors.
//
{
  int i;
  int j;
  int k;
  int kd;
  int l;
  int lm;
  double temp;
  double *x;
//
//  Test the input parameters.
//
  if ( trans != 'N' && trans != 'n' &&
       trans != 'T' && trans != 't' &&
       trans != 'C' && trans != 'c' )
  {
    return NULL;
  }
  else if ( n <= 0 )
  {
    return NULL;
  }
  else if ( ml < 0 )
  {
    return NULL;
  }
  else if ( mu < 0 )
  {
    return NULL;
  }
  else if ( nrhs <= 0 )
  {
    return NULL;
  }

  x = new double[n*nrhs];

  for ( k = 0; k < nrhs; k++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = b[i+k*n];
    }
  }

  kd = mu + ml + 1;
//
//  Solve A * x = b.
//
//  Solve L * x = b, overwriting b with x.
//
//  L is represented as a product of permutations and unit lower
//  triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
//  where each transformation L(i) is a rank-one modification of
//  the identity matrix.
//
  if ( trans == 'N' || trans == 'n' )
  {
    if ( 0 < ml )
    {
      for ( j = 1; j <= n-1; j++ )
      {
        lm = i4_min ( ml, n-j );
        l = pivot[j-1];

        for ( k = 0; k < nrhs; k++ )
        {
          temp       = x[l-1+k*n];
          x[l-1+k*n] = x[j-1+k*n];
          x[j-1+k*n] = temp;
        }

        for ( k = 0; k < nrhs; k++ )
        {
          if ( x[j-1+k*n] != 0.0 )
          {
            for ( i = 1; i <= lm; i++ )
            {
              x[j+i-1+k*n] = x[j+i-1+k*n] 
                - a[kd+i-1+(j-1)*(2*ml+mu+1)] * x[j-1+k*n];
            }
          }
        }
      }
    }
//
//  Solve U * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 1 <= j; j-- )
      {
        if ( x[j-1+k*n] != 0.0 )
        {
          l = ml + mu + 1 - j;
          x[j-1+k*n] = x[j-1+k*n] / a[ml+mu+(j-1)*(2*ml+mu+1)];
          for ( i = j - 1; i4_max ( 1, j - ml - mu ) <= i; i-- )
          {
            x[i-1+k*n] = x[i-1+k*n] 
              - a[l+i-1+(j-1)*(2*ml+mu+1)] * x[j-1+k*n];
          }
        }
      }
    }
  }
  else
  {
//
//  Solve A' * x = b.
//
//  Solve U' * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n; j++ )
      {
        temp = x[j-1+k*n];
        l = ml + mu + 1 - j;
        for ( i = i4_max ( 1, j - ml - mu ); i <= j - 1; i++ )
        {
          temp = temp - a[l+i-1+(j-1)*(2*ml+mu+1)] * x[i-1+k*n];
        }
        x[j-1+k*n] = temp / a[ml+mu+(j-1)*(2*ml+mu+1)];
      }
    }
//
//  Solve L' * x = b, overwriting b with x.
//
    if ( 0 < ml )
    {
      for ( j = n-1; 1 <= j; j-- )
      {
        lm = i4_min ( ml, n-j );

        for ( k = 0; k < nrhs; k++ )
        {
          for ( i = 1; i <= lm; i++ )
          {
            x[j-1+k*n] = x[j-1+k*n] 
              - x[j+i-1+k*n] * a[kd+i-1+(j-1)*(2*ml+mu+1)];
          }
        }

        l = pivot[j-1];

        for ( k = 0; k < nrhs; k++ )
        {
          temp       = x[l-1+k*n];
          x[l-1+k*n] = x[j-1+k*n];
          x[j-1+k*n] = temp;
        }
      }
    }
  }

  return x;
}
//****************************************************************************80

double *r8gb_vxm ( int m, int n, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_VXM multilies a vector times a R8GB matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//    For our purposes, X*A and A'*X mean the same thing.
//
//    LINPACK and LAPACK storage of general band matrices requires
//    an extra ML upper diagonals for possible fill in entries during
//    Gauss elimination.  This routine does not access any entries
//    in the fill in diagonals, because it assumes that the matrix
//    has NOT had Gauss elimination applied to it.  If the matrix
//    has been Gauss eliminated, then the routine R8GB_MU must be
//    used instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2003
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8GB_VXM[N], the product X*A or A'*X.
//
{
  double *b;
  int col = 2 * ml + mu + 1;
  int i;
  int ihi;
  int ilo;
  int j;

  b = new double[n];

  for ( j = 1; j <= n; j++ )
  {
    b[j-1] = 0.0;
    ilo = i4_max ( 1, j - mu );
    ihi = i4_min ( m, j + ml );
    for ( i = ilo; i <= ihi; i++ )
    {
      b[j-1] = b[j-1] + x[i-1] * a[i-j+ml+mu+(j-1)*col];
    }
  }

  return b;
}
//****************************************************************************80

double *r8gb_zero ( int m, int n, int ml, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8GB_ZERO zeros a R8GB matrix.
//
//  Discussion:
//
//    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
//    which may be required to store nonzero entries generated during Gaussian 
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.  
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be nonnegative.
//
//    Input, int N, the number of columns of the matrix.
//    N must be nonnegative.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative and no greater than min(M,N)-1.
//
//    Output, double R8GB_ZERO[(2*ML+MU+1)*N], the R8GB matrix.
//
{
  double *a;
  int col = 2 * ml + mu + 1;
  int j;
  int row;

  a = new double[col*n];

  for ( j = 0; j < n; j++ )
  {
    for ( row = 0; row < col; row++ )
    {
      a[row+j*col] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

bool r8gd_error ( int n, int ndiag )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_ERROR checks the dimensions of a R8GD matrix.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Output, bool R8GD_ERROR, is TRUE if an error was detected.
//
{
  if ( n < 1 )
  {
    cout << "\n";
    cout << "R8GD_ERROR - Illegal N = " << n << "\n";
    return true;
  }

  if ( ndiag < 1 || 2 * n - 1 < ndiag )
  {
    cout << "\n";
    cout << "R8GD_ERROR - Illegal NDIAG = " << ndiag << "\n";
    return true;
  }

  return false;
}
//****************************************************************************80

double *r8gd_indicator ( int n, int ndiag, int offset[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_INDICATOR sets up a R8GD indicator matrix.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Output, double R8GD_INDICATOR[N*NDIAG], the R8GD matrix.
//
{
  double *a;
  int diag;
  int fac;
  int i;
  int j;

  a = new double[n*ndiag];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= n; i++ )
  {
    for ( diag = 1; diag <= ndiag; diag++ )
    {
      j = i + offset[diag-1];
      if ( 1 <= j && j <= n )
      {
        a[i-1+(diag-1)*n] = ( double ) ( fac * i + j );
      }
      else
      {
        a[i-1+(diag-1)*n] = 0.0;
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8gd_mxv ( int n, int ndiag, int offset[], double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_MXV multiplies a R8GD matrix by a vector.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Example:
//
//    The "offset" value is printed near the first entry of each diagonal
//    in the original matrix, and above the columns in the new matrix.
//
//    Original matrix               New Matrix
//
//      0    1   2   3   4   5        -3  -2   0   1   3   5
//                             
//        11  12   0  14   0  16      --  --  11  12  14  16
//   -1 =  0  22  23   0  25   0      --  --  22  23  25  --
//   -2 = 31   0  33  34   0  36      --  31  33  34  36  --
//   -3 = 41  42   0  44  45   0      41  42  44  45  --  --
//   -4 =  0  52  53   0  55  56      52  53  55  56  --  --
//   -5 =  0   0  63  64  65  66      63  64  66  --  --  --
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8GD matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8GD_MXV[N], the product A * x.
//
{
  double *b;
  int diag;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( diag = 0; diag < ndiag; diag++ )
    {
      j = i + offset[diag];
      if ( 0 <= j && j < n )
      {
        b[i] = b[i] + a[i+diag*n] * x[j];
      }
    }
  }

  return b;
}
//****************************************************************************80

void r8gd_print ( int n, int ndiag, int offset[], double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_PRINT prints a R8GD matrix.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8GD matrix.
//
//    Input, string TITLE, a title.
//
{
  r8gd_print_some ( n, ndiag, offset, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8gd_print_some ( int n, int ndiag, int offset[], double a[], int ilo, 
  int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_PRINT_SOME prints some of a R8GD matrix.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8GD matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int diag;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col:  ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j2 = j2lo; j2 <= j2hi; j2++ )
      {
        aij = 0.0;
        for ( diag = 0; diag < ndiag; diag++ )
        {
          if ( j2 - i == offset[diag] )
          {
            aij = a[i-1+diag*n];
          }
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8gd_random ( int n, int ndiag, int offset[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_RANDOM randomizes a R8GD matrix.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8GD_RANDOM[N*NDIAG], the R8GD matrix.
//
{
  double *a;
  int diag;
  int i;
  int j;

  a = new double[n*ndiag];

  for ( i = 1; i <= n; i++ )
  {
    for ( diag = 0; diag < ndiag; diag++ )
    {
      j = i + offset[diag];
      if ( 1 <= j && j <= n )
      {
        a[i-1+diag*n] = r8_uniform_01 ( seed );
      }
      else
      {
        a[i-1+diag*n] = 0.0;
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8gd_to_r8ge ( int n, int ndiag, int offset[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_TO_R8GE copies a R8GD matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8GD matrix.
//
//    Output, double R8GD_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int diag;
  int i;
  int j;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    for ( diag = 0; diag < ndiag; diag++ )
    {
      j = i + offset[diag];
      if ( 0 <= j && j <= n-1 )
      {
        b[i+j*n] = a[i+diag*n];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8gd_vxm ( int n, int ndiag, int offset[], double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_VXM multiplies a vector by a R8GD matrix.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8GD matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8GD_VXM[N], the product X*A.
//
{
  double *b;
  int diag;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( diag = 0; diag < ndiag; diag++ )
    {
      j = i + offset[diag];
      if ( 0 <= j && j < n )
      {
        b[j] = b[j] + x[i] * a[i+diag*n];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8gd_zero ( int n, int ndiag )

//****************************************************************************80
//
//  Purpose:
//
//    R8GD_ZERO zeros a R8GD matrix.
//
//  Discussion:
//
//    The R8GD storage format is used for matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0.
//    Each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
//
//    Now, assuming that only a few of these diagonals contain nonzeros,
//    then for the I-th diagonal to be saved, we stored its offset in
//    OFFSET(I), and its entries in column I of the matrix.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than 2 * N - 1.
//
//    Output, double R8GD_ZERO[N*NDIAG], the R8GD matrix.
//
{
  double *a;
  int diag;
  int i;

  a = new double[n*ndiag];

  for ( diag = 0; diag < ndiag; diag++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+diag*n] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

double r8ge_co ( int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_CO factors a R8GE matrix and estimates its condition number.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    For the system A * X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    EPSILON/RCOND.
//
//    If RCOND is so small that the logical expression
//      1.0 + rcond == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate
//    underflows.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2003
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input/output, double A[N*N].  On input, a matrix to be factored.
//    On output, the LU factorization of the matrix.
//
//    Output, int PIVOT[N], the pivot indices.
//
//    Output, double R8GE_CO, an estimate of the reciprocal condition number of A.
//
{
  double anorm;
  double ek;
  int i;
  int info;
  int j;
  int k;
  int l;
  double rcond;
  double s;
  double sm;
  double t;
  double wk;
  double wkm;
  double ynorm;
  double *z;
//
//  Compute the L1 norm of A.
//
  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    s = 0.0;
    for ( i = 0; i < n; i++ )
    {
      s = s + r8_abs ( a[i+j*n] );
    }
    anorm = r8_max ( anorm, s );
  }
//
//  Compute the LU factorization.
//
  info = r8ge_fa ( n, a, pivot );

  if ( info != 0 ) 
  {
    rcond = 0.0;
    return rcond;
  }
//
//  RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )
//
//  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
//
//  where
//    A * Z = Y
//  and
//    A' * Y = E
//
//  The components of E are chosen to cause maximum local growth in the
//  elements of W, where U'*W = E.  The vectors are frequently rescaled
//  to avoid overflow.
//
//  Solve U' * W = E.
//
  ek = 1.0;
  z = new double[n];
  for ( i = 0; i < n; i++ )
  {
    z[i] = 0.0;
  }

  for ( k = 0; k < n; k++ )
  {
    if ( z[k] != 0.0 ) 
    {
      ek = - r8_sign2 ( ek, z[k] );
    }

    if ( r8_abs ( a[k+k*n] ) < r8_abs ( ek - z[k] ) )
    {
      s = r8_abs ( a[k+k*n] ) / r8_abs ( ek - z[k] );
      for ( i = 0; i < n; i++ )
      {
        z[i] = s * z[i];
      }
      ek = s * ek;
    }

    wk = ek - z[k];
    wkm = -ek - z[k];
    s = r8_abs ( wk );
    sm = r8_abs ( wkm );

    if ( a[k+k*n] != 0.0 )
    {
      wk = wk / a[k+k*n];
      wkm = wkm / a[k+k*n];
    }
    else
    {
      wk = 1.0;
      wkm = 1.0;
    }

    if ( k + 2 <= n )
    {
      for ( j = k+1; j < n; j++ )
      {
        sm = sm + r8_abs ( z[j] + wkm * a[k+j*n] );
        z[j] = z[j] + wk * a[k+j*n];
        s = s + r8_abs ( z[j] );
      }

      if ( s < sm )
      {
        t = wkm - wk;
        wk = wkm;
        for ( j = k+1; j < n; j++ )
        {
          z[j] = z[j] + t * a[k+j*n];
        }
      }
    }
    z[k] = wk;
  }

  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / s;
  }
//
//  Solve L' * Y = W
//
  for ( k = n-1; 0 <= k; k-- )
  {
    for ( i = k+1; i < n; i++ )
    {
      z[k] = z[k] + z[i] * a[i+k*n];
    }
    t = r8_abs ( z[k] );
    if ( 1.0 < t )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i] = z[i] / t;
      }
    }

    l = pivot[k] - 1;

    t    = z[l];
    z[l] = z[k];
    z[k] = t;
  }

  t = 0.0;
  for ( i = 0; i < n; i++ )
  {
    t = t + r8_abs ( z[i] );
  }
  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / t;
  }

  ynorm = 1.0;
//
//  Solve L * V = Y.
//
  for ( k = 0; k < n; k++ )
  {
    l = pivot[k] - 1;

    t    = z[l];
    z[l] = z[k];
    z[k] = t;

    for ( i = k+1; i < n; i++ )
    {
      z[i] = z[i] + t * a[i+k*n];
    }

    t = r8_abs ( z[k] );

    if ( 1.0 < t )
    {
      ynorm = ynorm / t;
      for ( i = 0; i < n; i++ )
      {
        z[i] = z[i] / t;
      }
    }
  }
  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / s;
  }
  ynorm = ynorm / s;
//
//  Solve U * Z = V.
//
  for ( k = n-1; 0 <= k; k-- )
  {
    if ( r8_abs ( a[k+k*n] ) < r8_abs ( z[k] ) )
    {
      s = r8_abs ( a[k+k*n] ) / r8_abs ( z[k] );
      for ( i = 0; i < n; i++ )
      {
        z[i] = s * z[i];
      }
      ynorm = s * ynorm;
    }

    if ( a[k+k*n] != 0.0 )
    {
      z[k] = z[k] / a[k+k*n];
    }
    else
    {
      z[k] = 1.0;
    }

    for ( i = 0; i < k; i++ )
    {
      z[i] = z[i] - a[i+k*n] * z[k];
    }
  }
//
//  Normalize Z in the L1 norm.
//
  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }
  s = 1.0 / s;

  for ( i = 0; i < n; i++ )
  {
    z[i] = s * z[i];
  }
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;

  return rcond;
}
//****************************************************************************80

double r8ge_det ( int n, double a_lu[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_DET computes the determinant of a matrix factored by R8GE_FA or R8GE_TRF.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA or R8GE_TRF.
//
//    Input, int PIVOT[N], as computed by R8GE_FA or R8GE_TRF.
//
//    Output, double R8GE_DET, the determinant of the matrix.
//
{
  double det;
  int i;

  det = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    det = det * a_lu[i-1+(i-1)*n];
    if ( pivot[i-1] != i )
    {
      det = -det;
    }
  }

  return det;
}
//****************************************************************************80

double *r8ge_dilu ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_DILU produces the diagonal incomplete LU factor of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2003
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
//    Input, double A[M*N], the R8GE matrix.
//
//    Output, double R8GE_DILU[M], the D-ILU factor.
//
{
  double *d;
  int i;
  int j;

  d = new double[m];

  for ( i = 0; i < m; i++ )
  {
    if ( i < n ) 
    {
      d[i] = a[i+i*m];
    }
    else
    {
      d[i] = 0.0;
    }
  }

  for ( i = 0; i < m && i < n; i++ )
  {
    d[i] = 1.0 / d[i];
    for ( j = i+1; j < m && j < n; j++ )
    {
      d[j] = d[j] - a[j+i*m] * d[i] * a[i+j*m];
    }
  }

  return d;
}
//****************************************************************************80

int r8ge_fa ( int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FA performs a LINPACK-style PLU factorization of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
//
//    The two dimensional array is stored by columns in a one dimensional
//    array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N], the matrix to be factored.
//    On output, A contains an upper triangular matrix and the multipliers
//    which were used to obtain it.  The factorization can be written
//    A = L * U, where L is a product of permutation and unit lower
//    triangular matrices and U is upper triangular.
//
//    Output, int PIVOT[N], a vector of pivot indices.
//
//    Output, int R8GE_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int i;
  int j;
  int k;
  int l;
  double t;
//
  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L, the index of the pivot row.
//
    l = k;

    for ( i = k+1; i <= n; i++ )
    {
      if ( r8_abs ( a[l-1+(k-1)*n] ) < r8_abs ( a[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;
//
//  If the pivot index is zero, the algorithm has failed.
//
    if ( a[l-1+(k-1)*n] == 0.0 )
    {
      cerr << "\n";
      cerr << "R8GE_FA - Fatal error!\n";
      cerr << "  Zero pivot on step " << k << "\n";
      exit ( 1 );
    }
//
//  Interchange rows L and K if necessary.
//
    if ( l != k )
    {
      t              = a[l-1+(k-1)*n];
      a[l-1+(k-1)*n] = a[k-1+(k-1)*n];
      a[k-1+(k-1)*n] = t;
    }
//
//  Normalize the values that lie below the pivot entry A(K,K).
//
    for ( i = k+1; i <= n; i++ )
    {
      a[i-1+(k-1)*n] = -a[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        t              = a[l-1+(j-1)*n];
        a[l-1+(j-1)*n] = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = t;
      }

      for ( i = k+1; i <= n; i++ )
      {
        a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + a[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }

    }

  }

  pivot[n-1] = n;

  if ( a[n-1+(n-1)*n] == 0.0 )
  {
    cerr << "\n";
    cerr << "R8GE_FA - Fatal error!\n";
    cerr << "  Zero pivot on step " << n << "\n";
    exit ( 1 );
  }

  return 0;
}
//****************************************************************************80

void r8ge_fs ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FS factors and solves a R8GE system.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    The function does not save the LU factors of the matrix, and hence cannot
//    be used to efficiently solve multiple linear systems, or even to
//    factor A at one time, and solve a single linear system at a later time.
//
//    The function uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input/output, double X[N], on input, the right hand side of the linear system.
//    On output, the solution of the linear system.
//
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8GE_FS - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }
//
//  Scale the pivot row.
//
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  return;
}
//****************************************************************************80

double *r8ge_fs_new ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FS_NEW factors and solves a R8GE system.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    The function does not save the LU factors of the matrix, and hence cannot
//    be used to efficiently solve multiple linear systems, or even to
//    factor A at one time, and solve a single linear system at a later time.
//
//    The function uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double R8GE_FS_NEW[N], the solution of the linear system.
//
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8GE_FS_NEW - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }
//
//  Scale the pivot row.
//
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  return x;
}
//****************************************************************************80

void r8ge_fss ( int n, double a[], int nb, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FSS factors and solves multiple R8GE systems.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    This routine does not save the LU factors of the matrix, and hence cannot
//    be used to efficiently solve multiple linear systems, or even to
//    factor A at one time, and solve a single linear system at a later time.
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, int NB, the number of right hand sides.
//
//    Input/output, double X[N*NB], on input, the right hand sides of the
//    linear systems.  On output, the solutions of the linear systems.
//
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8GE_FSS - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
//
//  Scale the pivot row.
//
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return;
}
//****************************************************************************80

double *r8ge_fss_new ( int n, double a[], int nb, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FSS_NEW factors and solves multiple R8GE systems.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    This routine does not save the LU factors of the matrix, and hence cannot
//    be used to efficiently solve multiple linear systems, or even to
//    factor A at one time, and solve a single linear system at a later time.
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, int NB, the number of right hand sides.
//
//    Input, double B[N*NB], the right hand sides of the linear systems.
//
//    Output, double R8GE_FSS_NEW[N*NB], the solutions of the linear systems.
//
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  x = new double[n*nb];

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }
  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8GE_FSS_NEW - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
//
//  Scale the pivot row.
//
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return x;
}
//****************************************************************************80

double *r8ge_identity ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_IDENTITY sets a R8GE matrix to the identity.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Output, double R8GE_IDENTITY[N*N], the N by N identity matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  return a;
}
//****************************************************************************80

void r8ge_ilu ( int m, int n, double a[], double l[], double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_ILU produces the incomplete LU factors of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    The incomplete LU factors of the M by N matrix A are:
//
//      L, an M by M unit lower triangular matrix,
//      U, an M by N upper triangular matrix
//
//    with the property that L and U are computed in the same way as
//    the usual LU factors, except that, whenever an off diagonal element
//    of the original matrix is zero, then the corresponding value of
//    U is forced to be zero.
//
//    This condition means that it is no longer the case that A = L*U.
//
//    On the other hand, L and U will have a simple sparsity structure
//    related to that of A.  The incomplete LU factorization is generally
//    used as a preconditioner in iterative schemes applied to sparse
//    matrices.  It is presented here merely for illustration.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 November 2003
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
//    Input, double A[M*N], the R8GE matrix.
//
//    Output, double L[M*M], the M by M unit lower triangular factor.
//
//    Output, double U[M*N], the M by N upper triangular factor.
//
{
  int i;
  int j;
  int jhi;
  int k;
//
//  Initialize:
//
//    L := M by M Identity
//    U := A
//
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i == j )
      {
        l[i+j*m] = 1.0;
      }
      else
      {
        l[i+j*m] = 0.0;
      }
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      u[i+j*m] = a[i+j*m];
    }
  }

  jhi = i4_min ( m - 1, n );

  for ( j = 0; j < jhi; j++ )
  {
//
//  Zero out the entries in column J, from row J+1 to M.
//
    for ( i = j+1; i < m; i++ )
    {
      if ( u[i+j*m] != 0.0 )
      {
        l[i+j*m] = u[i+j*m] / u[j+j*m];
        u[i+j*m] = 0.0;

        for ( k = j + 1; k < n; k++ )
        {
          if ( u[i+k*m] != 0.0 )
          {
            u[i+k*m] = u[i+k*m] - l[i+j*m] * u[j+k*m];
          }
        }
      }
    }
  }

  return;
}
//****************************************************************************80

double *r8ge_indicator ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_INDICATOR sets up a R8GE indicator matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2005
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
//    Output, double R8GE_INDICATOR[M*N], the R8GE matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[m*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = ( double ) ( fac * i + j );
    }
  }

  return a;
}
//****************************************************************************80

double *r8ge_inverse ( int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_INVERSE computes the inverse of a R8GE matrix factored by R8GE_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
//    SGEDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, double A[N*N], the factor information computed by R8GE_FA.
//
//    Input, int PIVOT(N), the pivot vector from R8GE_FA.
//
//    Output, double R8GE_INVERSE[N*N], the inverse matrix.
//
{
  double *b;
  int i;
  int j;
  int k;
  double temp;

  b = new double[n*n];
//
//  Compute Inverse(U).
//
  for ( k = 1; k <= n; k++ )
  {
    for ( i = 1; i <= k-1; i++ )
    {
      b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
    b[k-1+(k-1)*n] = 1.0 / a[k-1+(k-1)*n];

    for ( j = k+1; j <= n; j++ )
    {
      b[k-1+(j-1)*n] = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }
    }
  }
//
//  Multiply Inverse(U) by Inverse(L).
//
  for ( k = n-1; 1 <= k; k-- )
  {
    for ( i = k+1; i <= n; i++ )
    {
      b[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = b[i-1+(k-1)*n] + b[i-1+(j-1)*n] * a[j-1+(k-1)*n];
      }
    }

    if ( pivot[k-1] != k )
    {
      for ( i = 1; i <= n; i++ )
      {
        temp = b[i-1+(k-1)*n];
        b[i-1+(k-1)*n] = b[i-1+(pivot[k-1]-1)*n];
        b[i-1+(pivot[k-1]-1)*n] = temp;
      }

    }

  }

  return b;
}
//****************************************************************************80

double *r8ge_ml ( int n, double a_lu[], int pivot[], double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_ML computes A * x or A' * x, using R8GE_FA factors.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    It is assumed that R8GE_FA has overwritten the original matrix
//    information by LU factors.  R8GE_ML is able to reconstruct the
//    original matrix from the LU factor data.
//
//    R8GE_ML allows the user to check that the solution of a linear
//    system is correct, without having to save an unfactored copy
//    of the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA.
//
//    Input, int PIVOT[N], the pivot vector computed by R8GE_FA.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Input, int JOB, specifies the operation to be done:
//    JOB = 0, compute A * x.
//    JOB nonzero, compute A' * X.
//
//    Output, double R8GE_ML[N], the result of the multiplication.
//
{
  double *b;
  int i;
  int j;
  int k;
  double temp;
//
  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
//
//  Y = U * X.
//
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= j-1; i++ )
      {
        b[i-1] = b[i-1] + a_lu[i-1+(j-1)*n] * b[j-1];
      }
      b[j-1] = a_lu[j-1+(j-1)*n] * b[j-1];
    }
//
//  B = PL * Y = PL * U * X = A * x.
//
    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j+1; i <= n; i++ )
      {
        b[i-1] = b[i-1] - a_lu[i-1+(j-1)*n] * b[j-1];
      }

      k = pivot[j-1];

      if ( k != j )
      {
        temp   = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = temp;
      }
    }
  }
  else
  {
//
//  Y = (PL)' * X:
//
    for ( j = 1; j <= n-1; j++ )
    {
      k = pivot[j-1];

      if ( k != j )
      {
        temp   = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = temp;
      }

      temp = 0.0;
      for ( i = j+1; i <= n; i++ )
      {
        temp = temp + b[i-1] * a_lu[i-1+(j-1)*n];
      }
      b[j-1] = b[j-1] - temp;

    }
//
//  B = U' * Y = ( PL * U )' * X = A' * X.
//
    for ( i = n; 1 <= i; i-- )
    {
      for ( j = i+1; j <= n; j++ )
      {
        b[j-1] = b[j-1] + b[i-1] * a_lu[i-1+(j-1)*n];
      }
      b[i-1] = b[i-1] * a_lu[i-1+(i-1)*n];
    }

  }

  return b;
}
//****************************************************************************80

double *r8ge_mu ( int m, int n, double a_lu[], char trans, int pivot[], 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_MU computes A * x or A' * x, using R8GE_TRF factors.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    It is assumed that R8GE_TRF has overwritten the original matrix
//    information by PLU factors.  R8GE_MU is able to reconstruct the
//    original matrix from the PLU factor data.
//
//    R8GE_MU allows the user to check that the solution of a linear
//    system is correct, without having to save an unfactored copy
//    of the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//  Parameters:
//
//    Input, int M, the number of rows in the matrix.
//
//    Input, int N, the number of columns in the matrix.
//
//    Input, double A_LU[M*N], the LU factors from R8GE_TRF.
//
//    Input, char TRANS, specifies the form of the system of equations:
//    'N':  A * x = b  (No transpose)
//    'T':  A'* X = B  (Transpose)
//    'C':  A'* X = B  (Conjugate transpose = Transpose)
//
//    Input, int PIVOT[*], the pivot vector computed by R8GE_TRF.
//
//    Input, double X[*], the vector to be multiplied.
//    For the untransposed case, X should have N entries.
//    For the transposed case, X should have M entries.
//
//    Output, double R8GE_MU[*], the result of the multiplication.
//    For the untransposed case, the result should have M entries.
//    For the transposed case, the result should have N entries.
//
{
  double *b;
  int i;
  int j;
  int k;
  int mn_max;
  int npiv;
  double temp;
  double *y;

  npiv = i4_min ( m - 1, n );
  mn_max = i4_max ( m, n );
  y = new double[mn_max];

  if ( trans == 'n' || trans == 'N' )
  {
    b = new double[m];
//
//  Y[MN] = U[MNxN] * X[N].
//
    for ( i = 0; i < n; i++ )
    {
      y[i] = 0.0;
    }

    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= i4_min ( j, m ); i++ )
      {
        y[i-1] = y[i-1] + a_lu[i-1+(j-1)*m] * x[j-1];
      }
    }
//
//  Z[M] = L[MxMN] * Y[MN] = L[MxMN] * U[MNxN] * X[N].
//
    for ( i = 0; i < m; i++ )
    {
      if ( i < n ) 
      {
        b[i] = y[i];
      }
      else
      {
        b[i] = 0.0;
      }
    }

    for ( j = i4_min ( m-1, n ); 1 <= j; j-- )
    {
      for ( i = j+1; i <= m; i++ )
      {
        b[i-1] = b[i-1] + a_lu[i-1+(j-1)*m] * y[j-1];
      }
    }
//
//  B = P * Z = P * L * Y = P * L * U * X = A * x.
//
    for ( j = npiv; 1 <= j; j-- )
    {
      k = pivot[j-1];

      if ( k != j )
      {
        temp = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = temp;
      }
    }
  }
  else if ( trans == 't' || trans == 'T' ||
            trans == 'c' || trans == 'C' )
  {
    b = new double[n];
//
//  Y = P' * X:
//
    for ( i = 1; i <= npiv; i++ )
    {
      k = pivot[i-1];

      if ( k != i )
      {
        temp = x[k-1];
        x[k-1] = x[i-1];
        x[i-1] = temp;
      }
    }

    for ( i = 0; i < n; i++ )
    {
      if ( i < m ) 
      {
        b[i] = x[i];
      }
      else
      {
        b[i] = 0.0;
      }
    }
//
//  Z = L' * Y:
//
    for ( j = 1; j <= i4_min ( m - 1, n ); j++ )
    {
      for ( i = j+1; i <= m; i++ )
      {
        b[j-1] = b[j-1] + x[i-1] * a_lu[i-1+(j-1)*m];
      }
    }
//
//  B = U' * Z.
//
    for ( i = m; 1 <= i; i-- )
    {
      for ( j = i+1; j <= n; j++ )
      {
        b[j-1] = b[j-1] + b[i-1] * a_lu[i-1+(j-1)*m];
      }
      if ( i <= n )
      {
        b[i-1] = b[i-1] * a_lu[i-1+(i-1)*m];
      }
    }
//
//  Now restore X.
//
    for ( i = npiv; 1 <= i; i-- )
    {
      k = pivot[i-1];

      if ( k != i )
      {
        temp = x[k-1];
        x[k-1] = x[i-1];
        x[i-1] = temp;
      }
    }
  }
//
//  Illegal value of TRANS.
//
  else
  {
    cerr << "\n";
    cerr << "R8GE_MU - Fatal error!\n";
    cerr << "  Illegal value of TRANS = \"" << trans << "\"\n";
    exit ( 1 );
  }

  delete [] y;

  return b;
}
//****************************************************************************80

double *r8ge_mxm ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_MXM multiplies two R8GE matrices.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, double A[N*N], B[N*N], the R8GE factor matrices.
//
//    Output, double C[N*N], the R8GE product matrix.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        c[i+j*n] = c[i+j*n] + a[i+k*n] * b[k+j*n];
      }
    }
  }

  return c;
}
//****************************************************************************80

double *r8ge_mxv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_MXV multiplies a R8GE matrix times a vector.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2003
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
//    Input, double A[M*N], the R8GE matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8GE_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

double r8ge_np_det ( int n, double a_lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_NP_DET computes the determinant of a matrix factored by R8GE_NP_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.
//
//    Output, double R8GE_NP_DET, the determinant of the matrix.
//
{
  double det;
  int i;

  det = 1.0;
  for ( i = 0; i < n; i++ )
  {
    det = det * a_lu[i+i*n];
  }

  return det;
}
//****************************************************************************80

int r8ge_np_fa ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_NP_FA factors a R8GE matrix by nonpivoting Gaussian elimination.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_NP_FA is a version of the LINPACK routine SGEFA, but uses no
//    pivoting.  It will fail if the matrix is singular, or if any zero
//    pivot is encountered.
//
//    If R8GE_NP_FA successfully factors the matrix, R8GE_NP_SL may be called
//    to solve linear systems involving the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N].
//    On input, A contains the matrix to be factored.
//    On output, A contains information about the factorization,
//    which must be passed unchanged to R8GE_NP_SL for solutions.
//
//    Output, int R8GE_NP_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int i;
  int j;
  int k;

  for ( k = 1; k <= n-1; k++ )
  {
    if ( a[k-1+(k-1)*n] == 0.0 )
    {
      return k;
    }

    for ( i = k+1; i <= n; i++ )
    {
      a[i-1+(k-1)*n] = -a[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = k+1; i <= n; i++ )
      {
        a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + a[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }
    }
  }

  if ( a[n-1+(n-1)*n] == 0.0 )
  {
    return n;
  }

  return 0;
}
//****************************************************************************80

double *r8ge_np_inverse ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_NP_INVERSE computes the inverse of a matrix factored by R8GE_NP_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, double A[N*N], the factor information computed by R8GE_NP_FA.
//
//    Output, double R8GE_NP_INVERSE[N*N], the inverse matrix.
//
{
  double *b;
  int i;
  int j;
  int k;
  double temp;
  double *work;

  b = new double[n*n];
  work = new double[n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }
//
//  Compute Inverse(U).
//
  for ( k = 1; k <= n; k++ )
  {
    b[k-1+(k-1)*n] = 1.0 / b[k-1+(k-1)*n];
    for ( i = 1; i <= k-1; i++ )
    {
      b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] * b[k-1+(k-1)*n];
    }
    for ( j = k+1; j <= n; j++ )
    {
      temp = b[k-1+(j-1)*n];
      b[k-1+(j-1)*n] = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + temp * b[i-1+(k-1)*n];
      }
    }
  }
//
//  Form Inverse(U) * Inverse(L).
//
  for ( k = n-1; 1 <= k; k-- )
  {
    for ( i = k+1; i <= n; i++ )
    {
      work[i-1] = b[i-1+(k-1)*n];
      b[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = b[i-1+(k-1)*n] + b[i-1+(j-1)*n] * work[j-1];
      }
    }
  }

  delete [] work;

  return b;
}
//****************************************************************************80

double *r8ge_np_ml ( int n, double a_lu[], double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_NP_ML computes A * x or x * A, for a matrix factored by R8GE_NP_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    The matrix A is assumed to have been factored by R8GE_NP_FA.
//
//    R8GE_NP_ML allows the user to check that the solution of a linear
//    system is correct, without having to save an unfactored copy
//    of the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Input, int JOB, determines the multiplication to
//    be carried out:
//    JOB = 0, compute A * x.
//    JOB nonzero, compute A' * X.
//
//    Output, double R8GE_NP_ML[N], the result of the multiplication.
//
{
  double *b;
  int i;
  int j;
  double t;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
//
//  Compute U * X = Y:
//
    for ( i = 0; i < n; i++ )
    {
      t = 0.0;
      for ( j = i; j < n; j++ )
      {
        t = t + a_lu[i+j*n] * b[j];
      }
      b[i] = t;
    }
//
//  Compute L * Y = B:
//
    for ( j = n-2; 0 <= j; j-- )
    {
      for ( i = j+1; i < n; i++ )
      {
        b[i] = b[i] - a_lu[i+j*n] * b[j];
      }
    }
  }
  else
  {
//
//  Compute L' * X = Y:
//
    for ( j = 0; j < n-1; j++ )
    {
      for ( i = j+1; i < n; i++ )
      {
        b[j] = b[j] - b[i] * a_lu[i+j*n];
      }
    }
//
//  Compute U' * Y = B:
//
    for ( j = n-1; 0 <= j; j-- )
    {
      b[j] = b[j] * a_lu[j+j*n];
      for ( i = 0; i < j; i++ )
      {
        b[j] = b[j] + b[i] * a_lu[i+j*n];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8ge_np_sl ( int n, double a_lu[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_NP_SL solves a system factored by R8GE_NP_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB.
//    If JOB is zero, the routine will solve A * x = b.
//    If JOB is nonzero, the routine will solve A' * x = b.
//
//    Output, double R8GE_NP_SL[N], the solution.
//
{
  int i;
  int k;
  double *x;
//
//  Solve A * x = b.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( k = 0; k < n-1; k++ )
    {
      for ( i = k+1; i < n; i++ )
      {
        x[i] = x[i] + a_lu[i+k*n] * x[k];
      }
    }

    for ( k = n-1; 0 <= k; k-- )
    {
      x[k] = x[k] / a_lu[k+k*n];
      for ( i = 0; i <= k-1; i++ )
      {
        x[i] = x[i] - a_lu[i+k*n] * x[k];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
    for ( k = 0; k < n; k++ )
    {
      for ( i = 0; i <= k-1; i++ )
      {
        x[k] = x[k] - x[i] * a_lu[i+k*n];
      }
      x[k] = x[k] / a_lu[k+k*n];
    }

    for ( k = n-2; 0 <= k; k-- )
    {
      for ( i = k+1; i < n; i++ )
      {
        x[k] = x[k] + x[i] * a_lu[i+k*n];
      }
    }

  }

  return x;
}
//****************************************************************************80

int r8ge_np_trf ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_NP_TRF computes the LU factorization of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_NP_TRF is a nonpivoting version of R8GE_TRF, and will fail if
//    a zero element is encountered along the diagonal.
//
//    The factorization has the form
//      A = L * U
//    where L is lower triangular with unit diagonal elements (lower
//    trapezoidal if N < M), and U is upper triangular (upper trapezoidal
//    if M < N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix A.  0 <= M.
//
//    Input, int N, the number of columns of the matrix A.  0 <= N.
//
//    Input/output, double A[M*N].
//    On entry, the M by N matrix to be factored.
//    On exit, the factors L and U from the factorization
//    A = L*U; the unit diagonal elements of L are not stored.
//
//    Output, int R8GE_NP_TRF.
//    = 0: successful exit
//    = -K, the K-th argument had an illegal value
//    = K, U(K,K) is exactly zero. The factorization
//         has been completed, but the factor U is exactly
//         singular, and division by zero will occur if it is used
//         to solve a system of equations.
//
{
  int i;
  int ii;
  int info;
  int j;
//
//  Test the input parameters.
//
  info = 0;

  if ( m < 0 )
  {
    return (-1);
  }
  else if ( n < 0 )
  {
    return (-2);
  }

  if ( m == 0 || n == 0 )
  {
    return 0;
  }

  for ( j = 1; j <= i4_min ( m, n ); j++ )
  {
//
//  Compute elements J+1:M of the J-th column.
//
    if ( a[j-1+(j-1)*m] != 0.0 )
    {
      for ( i = j+1; i <= m; i++ )
      {
        a[i-1+(j-1)*m] = a[i-1+(j-1)*m] / a[j-1+(j-1)*m];
      }
    }
    else if ( info == 0 )
    {
      info = j;
    }
//
//  Update the trailing submatrix.
//
    if ( j < i4_min ( m, n ) )
    {
      for ( ii = j+1; ii <= m; ii++ )
      {
        for ( i = j+1; i <= n; i++ )
        {
          a[ii-1+(i-1)*m] = a[ii-1+(i-1)*m] - a[ii-1+(j-1)*m] * a[j-1+(i-1)*m];
        }
      }
    }
  }

  return info;
}
//****************************************************************************80

double *r8ge_np_trm ( int m, int n, double a[], double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_NP_TRM computes A * x or A' * x, for a matrix factored by R8GE_NP_TRF.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    The matrix A is assumed to have been factored by R8GE_NP_TRF.
//
//    R8GE_NP_TRM allows the user to check that the solution of a linear
//    system is correct, without having to save an unfactored copy
//    of the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//    M and N must be positive.
//
//    Input, double A[M*N], the M by N matrix factors computed by R8GE_NP_TRF.
//
//    Input, double X[*], the vector to be multiplied.
//    If JOB is 0, X must have dimension N.
//    If JOB is nonzero, X must have dimension M.
//
//    Input, int JOB, determines the multiplication to
//    be carried out:
//    JOB = 0, compute A * x.
//    JOB nonzero, compute A' * X.
//
//    Output, double R8GE_NP_TRM[*], the result of the multiplication.
//    If JOB is 0, the output has dimension M.
//    If JOB is nonzero, the output has dimension N.
//
{
  double *b;
  int i;
  int j;
  double temp;

  if ( job == 0 )
  {
    b = new double[m];
    for ( i = 0; i < m; i++ )
    {
      b[i] = 0.0;
    }
//
//  Compute U * X = Y:
//
    for ( i = 0; i < i4_min ( m, n ); i++ )
    {
      for ( j = i; j < n; j++ )
      {
        b[i] = b[i] + a[i+j*m] * x[j];
      }
    }
//
//  Compute L * Y = B:
//
    for ( i = i4_min ( m-1, n ); 1 <= i; i-- )
    {
      for ( j = 0; j < i; j++ )
      {
        b[i] = b[i] + a[i+j*m] * b[j];
      }
    }
  }
  else
  {
    b = new double[n];
    for ( i = 0; i < n; i++ )
    {
      b[i] = 0.0;
    }
//
//  Compute L' * X = Y:
//
    for ( i = 0; i < i4_min ( m, n ); i++ )
    {
      b[i] = x[i];
      for ( j = i+1; j < m; j++ )
      {
        b[i] = b[i] + a[j+i*m] * x[j];
      }
    }
//
//  Compute U' * Y = B:
//
    for ( i = i4_min ( m, n ) - 1; 0 <= i; i-- )
    {
      temp = 0.0;
      for ( j = 0; j <= i; j++ )
      {
        temp = temp + a[j+i*m] * b[j];
      }
      b[i] = temp;
    }

  }

  return b;
}
//****************************************************************************80

double *r8ge_np_trs ( int n, int nrhs, char trans, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_NP_TRS solves a system of linear equations factored by R8GE_NP_TRF.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_NP_TRS is a nonpivoting version of R8GE_TRS.
//
//    R8GE_TRS solves a system of linear equations
//      A * x = b  or  A' * X = B
//    with a general N by N matrix A using the LU factorization computed
//    by R8GE_NP_TRF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.  0 <= N.
//
//    Input, int NRHS, the number of right hand sides.  0 <= NRHS.
//
//    Input, char TRANS, pecifies the form of the system of equations:
//    'N':  A * x = b  (No transpose)
//    'T':  A'* X = B  (Transpose)
//    'C':  A'* X = B  (Conjugate transpose = Transpose)
//
//    Input, double A[N*N], the factors L and U from the factorization
//    A = L*U as computed by R8GE_NP_TRF.
//
//    Input, double B[N*NRHS], the right hand side matrix B.
//
//    Output, double R8GE_NP_TRS[N*NRHS], the solution matrix X.
//
{
  int i;
  int j;
  int k;
  double *x;

  if ( trans != 'n' && trans != 'N' && 
       trans != 't' && trans != 'T' && 
       trans != 'c' && trans != 'C' )
  {
    return NULL;
  }
  if ( n < 0 )
  {
    return NULL;
  }
  if ( nrhs < 0 )
  {
    return NULL;
  }

  if ( n == 0 || nrhs == 0 )
  {
    return NULL;
  }

  x = new double[n*nrhs];

  for ( j = 0; j < nrhs; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }

  if ( trans == 'n' || trans == 'N' )
  {
//
//  Solve L * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n-1; j++ )
      {
        for ( i = j+1; i <= n; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[i-1+(j-1)*n] * x[j-1+k*n];
        }
      }
    }
//
//  Solve U * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 1 <= j; j-- )
      {
        x[j-1+k*n] = x[j-1+k*n] / a[j-1+(j-1)*n];
        for ( i = 1; i <= j-1; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[i-1+(j-1)*n] * x[j-1+k*n];
        }
      }
    }
  }
  else
//
//  Solve U' * x = b, overwriting b with x.
//
  {
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n; j++ )
      {
        x[j-1+k*n] = x[j-1+k*n] / a[j-1+(j-1)*n];
        for ( i = j+1; i <= n; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[j-1+(i-1)*n] * x[j-1+k*n];
        }
      }
    }
//
//  Solve L' * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 2 <= j; j-- )
      {
        for ( i = 1; i <= j-1; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[j-1+(i-1)*n] * x[j-1+k*n];
        }
      }
    }
  }

  return x;
}
//****************************************************************************80

void r8ge_plu ( int m, int n, double a[], double p[], double l[], double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_PLU produces the PLU factors of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    The PLU factors of the M by N matrix A are:
//
//      P, an M by M permutation matrix P,
//      L, an M by M unit lower triangular matrix,
//      U, an M by N upper triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2003
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
//    Input, double A[M,N], the M by N matrix to be factored.
//
//    Output, double P[M*M], the M by M permutation factor.
//
//    Output, double L[M*M], the M by M unit lower triangular factor.
//
//    Output, double U[M*N], the M by N upper triangular factor.
//
{
  int i;
  int j;
  int k;
  int pivot_row;
  double pivot_value;
  double t;
//
//  Initialize:
//
//    P: = M by M Identity
//    L: = M by M Identity
//    U: = A
//

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      if ( i == j ) 
      {
        l[i+j*m] = 1.0;
        p[i+j*m] = 1.0;
      }
      else
      {
        l[i+j*m] = 0.0;
        p[i+j*m] = 0.0;
      }
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      u[i+j*m] = a[i+j*m];
    }
  }
//
//  On step J, find the pivot row and the pivot value.
//
  for ( j = 0; j < i4_min ( m-2, n-1 ); j++ )
  {
    pivot_value = 0.0;
    pivot_row = -1;

    for ( i = j; i < m; i++ )
    {
      if ( pivot_value < r8_abs ( u[i+j*m] ) )
      {
        pivot_value = r8_abs ( u[i+j*m] );
        pivot_row = i;
      }
    }
//
//  If the pivot row is nonzero swap:
//  * rows J and PIVOT_ROW of U;
//  * rows J and PIVOT_ROW of L and cols J and PIVOT_ROW of L;
//  * cols J and PIVOT_ROW of P.
//
    if ( pivot_row != -1 )
    {
      for ( k = 0; k < n; k++ )
      {
        t                = u[j+k*m];
        u[j+k*m]         = u[pivot_row+k*m];
        u[pivot_row+k*m] = t;
      }

      for ( k = 0; k < m; k++ )
      {
        t                = l[j+k*m];
        l[j+k*m]         = l[pivot_row+k*m];
        l[pivot_row+k*m] = t;
      }
      for ( k = 0; k < m; k++ )
      {
        t                = l[k+j*m];
        l[k+j*m]         = l[k+pivot_row*m];
        l[k+pivot_row*m] = t;
      }

      for ( k = 0; k < m; k++ )
      {
        t                = p[k+j*m];
        p[k+j*m]         = p[k+pivot_row*m];
        p[k+pivot_row*m] = t;
      }
//
//  Zero out the entries in column J, from row J+1 to M.
//
      for ( i = j+1; i < m; i++ )
      {
        if ( u[i+j*m] != 0.0 )
        {
          l[i+j*m] = u[i+j*m] / u[j+j*m];
          u[i+j*m] = 0.0;
          for ( k = j+1; k < n; k++ )
          {
            u[i+k*m] = u[i+k*m] - l[i+j*m] * u[j+k*m];
          }
        }
      }
    }
  }

  return;
}
//****************************************************************************80

double *r8ge_poly ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_POLY computes the characteristic polynomial of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N*N], the R8GE matrix.
//
//    Output, double R8GE_POLY[N+1], the coefficients of the characteristic
//    polynomial of A.  P(I) contains the coefficient of X**I.
//
{
  int i;
  int j;
  int k;
  int order;
  double *p;
  double trace;
  double *work1;
  double *work2;

  p = new double[n+1];
  work2 = new double[n*n];
//
//  Initialize WORK1 to the identity matrix.
//
  work1 = r8ge_identity ( n );

  p[n] = 1.0;

  for ( order = n-1; 0 <= order; order-- )
  {
//
//  Work2 = A * WORK1.
//
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        work2[i+j*n] = 0.0;
        for ( k = 0; k < n; k++ )
        {
          work2[i+j*n] = work2[i+j*n] + a[i+k*n] * work1[k+j*n];
        }
      }
    }
//
//  Take the trace.
//
    trace = 0.0;
    for ( i = 0; i < n; i++ )
    {
      trace = trace + work2[i+i*n];
    }
//
//  P(ORDER) = - Trace ( WORK2 ) / ( N - ORDER )
//
    p[order] = -trace / ( double ) ( n - order );
//
//  WORK1 := WORK2 + P(ORDER) * Identity.
//
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        work1[i+j*n] = work2[i+j*n];
      }
    }
    for ( j = 0; j < n; j++ )
    {
      work1[j+j*n] = work1[j+j*n] + p[order];
    }
  }

  delete [] work1;
  delete [] work2;

  return p;
}
//****************************************************************************80

void r8ge_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_PRINT prints a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[M*N], the R8GE matrix.
//
//    Input, string TITLE, a title.
//
{
  r8ge_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8ge_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_PRINT_SOME prints some of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[M*N], the R8GE matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8ge_random ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_RANDOM randomizes a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2004
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
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8GE_RANDOM[M*N], the randomized M by N matrix, 
//    with entries between 0 and 1.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = r8_uniform_01 ( seed );
    }
  }
  return a;
}
//****************************************************************************80

double *r8ge_res ( int m, int n, double a[], double x[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_RES computes the residual for a R8GE system.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the original, UNFACTORED matrix.
//
//    Input, double X[N], an estimate of the solution the linear system.
//
//    Input, double B[M], the right hand side vector.
//
//    Output, double R8GE_RES[M], the residual vector: b - A * x.
//
{
  double *r;
  int i;
  int j;

  r = new double[n];

  for ( i = 0; i < m; i++ )
  {
    r[i] = b[i];
    for ( j = 0; j < n; j++ )
    {
      r[i] = r[i] - a[i+j*m] * x[j];
    }
  }

  return r;
}
//****************************************************************************80

void r8ge_sl ( int n, double a_lu[], int pivot[], double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_SL solves a R8GE system factored by R8GE_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_SL is a simplified version of the LINPACK routine SGESL.
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
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA.
//
//    Input, int PIVOT[N], the pivot vector from R8GE_FA.
//
//    Input/output, double X[N], on input, the right hand side vector.
//    On output, the solution vector.
//
//    Input, int JOB, specifies the operation.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
{
  int i;
  int k;
  int l;
  double t;
//
//  Solve A * x = b.
//
  if ( job == 0 )
  {
//
//  Solve PL * Y = B.
//
    for ( k = 1; k <= n-1; k++ )
    {
      l = pivot[k-1];

      if ( l != k )
      {
        t      = x[l-1];
        x[l-1] = x[k-1];
        x[k-1] = t;
      }
      for ( i = k+1; i <= n; i++ )
      {
        x[i-1] = x[i-1] + a_lu[i-1+(k-1)*n] * x[k-1];
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a_lu[k-1+(k-1)*n];
      for ( i = 1; i <= k-1; i++ )
      {
        x[i-1] = x[i-1] - a_lu[i-1+(k-1)*n] * x[k-1];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
//
//  Solve U' * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      t = 0.0;
      for ( i = 1; i <= k-1; i++ )
      {
        t = t + x[i-1] * a_lu[i-1+(k-1)*n];
      }
      x[k-1] = ( x[k-1] - t ) / a_lu[k-1+(k-1)*n];
    }
//
//  Solve ( PL )' * X = Y.
//
    for ( k = n-1; 1 <= k; k-- )
    {
      t = 0.0;
      for ( i = k+1; i <= n; i++ )
      {
        t = t + x[i-1] * a_lu[i-1+(k-1)*n];
      }
      x[k-1] = x[k-1] + t;

      l = pivot[k-1];

      if ( l != k )
      {
        t      = x[l-1];
        x[l-1] = x[k-1];
        x[k-1] = t;
      }
    }
  }

  return;
}
//****************************************************************************80

double *r8ge_sl_it ( int n, double a[], double a_lu[], int pivot[], double b[], 
  int job, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_SL_IT applies one step of iterative refinement following R8GE_SL.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    It is assumed that:
//
//    * the original matrix A has been factored by R8GE_FA;
//    * the linear system A * x = b has been solved once by R8GE_SL.
//
//    (Actually, it is not necessary to solve the system once using R8GE_SL.
//    You may simply supply the initial estimated solution X = 0.)
//
//    Each time this routine is called, it will compute the residual in
//    the linear system, apply one step of iterative refinement, and
//    add the computed correction to the current solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N*N], the original, UNFACTORED R8GE matrix.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA.
//
//    Input, int PIVOT[N], the pivot vector from R8GE_FA.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, int JOB, specifies the operation.
//    0, solve A*X=B.
//    nonzero, solve A'*X=B.
//
//    Input, double X[N], an estimate of the solution of A * x = b.
//
//    Output, double R8GE_SL_IT[N], the solution after one step of 
//    iterative refinement.
//
{
  double *dx;
  int i;
  double *r;
  double *x_new;
//
//  Compute the residual vector.
//
  r = r8ge_res ( n, n, a, x, b );
//
//  Solve A * dx = r
//
  dx = r8ge_sl_new ( n, a_lu, pivot, r, job );
//
//  Add dx to x.
//
  x_new = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x_new[i] = x[i] + dx[i];
  }

  delete [] r;
  delete [] dx;

  return x_new;
}
//****************************************************************************80

double *r8ge_sl_new ( int n, double a_lu[], int pivot[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_SL_NEW solves a R8GE system factored by R8GE_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_SL is a simplified version of the LINPACK routine SGESL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA.
//
//    Input, int PIVOT[N], the pivot vector from R8GE_FA.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, int JOB, specifies the operation.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double R8GE_SL[N], the solution vector.
//
{
  int i;
  int k;
  int l;
  double t;
  double *x;
//
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }
//
//  Solve A * x = b.
//
  if ( job == 0 )
  {
//
//  Solve PL * Y = B.
//
    for ( k = 1; k <= n-1; k++ )
    {
      l = pivot[k-1];

      if ( l != k )
      {
        t      = x[l-1];
        x[l-1] = x[k-1];
        x[k-1] = t;
      }
      for ( i = k+1; i <= n; i++ )
      {
        x[i-1] = x[i-1] + a_lu[i-1+(k-1)*n] * x[k-1];
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a_lu[k-1+(k-1)*n];
      for ( i = 1; i <= k-1; i++ )
      {
        x[i-1] = x[i-1] - a_lu[i-1+(k-1)*n] * x[k-1];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
//
//  Solve U' * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      t = 0.0;
      for ( i = 1; i <= k-1; i++ )
      {
        t = t + x[i-1] * a_lu[i-1+(k-1)*n];
      }
      x[k-1] = ( x[k-1] - t ) / a_lu[k-1+(k-1)*n];
    }
//
//  Solve ( PL )' * X = Y.
//
    for ( k = n-1; 1 <= k; k-- )
    {
      t = 0.0;
      for ( i = k+1; i <= n; i++ )
      {
        t = t + x[i-1] * a_lu[i-1+(k-1)*n];
      }
      x[k-1] = x[k-1] + t;

      l = pivot[k-1];

      if ( l != k )
      {
        t      = x[l-1];
        x[l-1] = x[k-1];
        x[k-1] = t;
      }

    }

  }

  return x;
}
//****************************************************************************80

double *r8ge_to_r8gb ( int m, int n, int ml, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_TO_R8GB copies a R8GE matrix to a R8GB matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    The R8GB storage format is for an M by N banded matrix, with lower
//    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
//    extra superdiagonals, which may be required to store nonzero entries
//    generated during Gaussian elimination.
//
//    It usually doesn't make sense to try to store a general matrix
//    in a band matrix format.  You can always do it, but it will take
//    more space, unless the general matrix is actually banded.
//
//    The purpose of this routine is to allow a user to set up a
//    banded matrix in the easy-to-use general format, and have this
//    routine take care of the compression of the data into general
//    format.  All the user has to do is specify the bandwidths.
//
//    Note that this routine "believes" what the user says about the
//    bandwidth.  It will assume that all entries in the general matrix
//    outside of the bandwidth are zero.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    LINPACK and LAPACK band storage requires that an extra ML
//    superdiagonals be supplied to allow for fillin during Gauss
//    elimination.  Even though a band matrix is described as
//    having an upper bandwidth of MU, it effectively has an
//    upper bandwidth of MU+ML.  This routine will copy nonzero
//    values it finds in these extra bands, so that both unfactored
//    and factored matrices can be handled.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrices.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrices.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths of A1.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Output, double A[M*N], the R8GE matrix.
//
//    Input, double R8GE_TO_R8GB[(2*ML+MU+1)*N], the R8GB matrix.
//
{
  double *b;
  int i;
  int j;
  int jhi;
  int jlo;
  int k;

  b = new double[(2*ml+mu+1)*n];

  for ( k = 0; k < (2*ml+mu+1)*n; k++ )
  {
    b[k] = 0.0;
  }

  for ( i = 1; i <= m; i++ )
  {
    jlo = i4_max ( i - ml, 1 );
    jhi = i4_min ( i + mu, n );

    for ( j = jlo; j <= jhi; j++ )
    {
      b[ml+mu+i-j+(j-1)*(2*ml+mu+1)] = a[i-1+(j-1)*m];
    }
  }

  return b;
}
//****************************************************************************80

void r8ge_to_r8ri ( int n, double a[], int nz, int ija[], double sa[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_TO_R8RI converts a R8GE matrix to R8RI form.
//
//  Discussion:
//
//    A R8GE matrix is in general storage.
//
//    An R8RI matrix is in row indexed sparse storage form.
//
//    The size of the arrays IJA and SA can be determined by calling
//    R8GE_TO_R8RI_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix stored in GE 
//    or "general" format.
//
//    Input, int NZ, the size required for the RI
//    or "row indexed" sparse storage.
//
//    Output, int IJA[NZ], the index vector.
//
//    Output, double SA[NZ], the value vector.
//
{
  int i;
  int im;
  int j;
  int k;
  int l;

  for ( k = 0; k < n; k++ )
  {
    i = k;
    j = k;
    sa[k] = a[i+j*n];
  }

  k = n;
  sa[k] = 0.0;

  for ( i = 0; i <= n; i++ )
  {
    ija[i] = 0;
  }

  im = 0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i != j )
      {
        if ( a[i+j*n] != 0.0 )
        {
          k = k + 1;
          if ( ija[i] == 0 )
          {
            for ( l = im; l <= i; l++ )
            {
              ija[l] = k;
            }
            im = i + 1;
          }
          ija[k] = j;
          sa[k] = a[i+j*n];
        }
      }
    }
  }

  ija[n] = k + 1;

  return;
}
//****************************************************************************80

int r8ge_to_r8ri_size ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_TO_R8RI_SIZE determines the size of an R8RI matrix.
//
//  Discussion:
//
//    N spaces are always used for the diagonal entries, plus a dummy.
//    The remaining spaces store off-diagonal nonzeros.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix stored in GE or "general" format.
//
//    Output, int R8GE_TO_R8RI_SIZE, the size required for the RI
//    or "row indexed" sparse storage.
//
{
  int i;
  int j;
  int nz;

  nz = n + 1;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i != j )
      {
        if ( a[i+j*n] != 0.0 )
        {
          nz = nz + 1;
        }
      }
    }
  }

  return nz;
}
//****************************************************************************80

double *r8ge_to_r8vec ( int m, int n, double *a )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_TO_R8VEC copies a R8GE matrix to a real vector.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    In C++  and FORTRAN, this routine is not really needed.  In MATLAB,
//    a data item carries its dimensionality implicitly, and so cannot be
//    regarded sometimes as a vector and sometimes as an array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, double R8VEC_TO_R8GE[M*N], the array to be copied.
//
//    Output, double X[M*N], the vector.
//
{
  int i;
  int j;
  int k;
  double *x;

  x = new double[m*n];

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[k] = a[i+j*m];
      k = k + 1;
    }
  }

  return x;
}
//****************************************************************************80

int r8ge_trf ( int m, int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_TRF performs a LAPACK-style PLU factorization of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_TRF is a standalone version of the LAPACK routine SGETRF.
//
//    The factorization uses partial pivoting with row interchanges,
//    and has the form
//      A = P * L * U
//    where P is a permutation matrix, L is lower triangular with unit
//    diagonal elements (lower trapezoidal if N < M), and U is upper
//    triangular (upper trapezoidal if M < N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 November 2003
//
//  Author:
//
//    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
//    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix A.  0 <= M.
//
//    Input, int N, the number of columns of the matrix A.  0 <= N.
//
//    Input/output, double A[M*N].
//    On entry, the M by N matrix to be factored.
//    On exit, the factors L and U from the factorization
//    A = P*L*U; the unit diagonal elements of L are not stored.
//
//    Output, int PIVOT[min(M,N)], the pivot indices.
//
//    Output, int R8GE_TRF.
//    = 0: successful exit
//    = -K, the K-th argument had an illegal value
//    = K: U(K,K) is exactly zero. The factorization
//         has been completed, but the factor U is exactly
//         singular, and division by zero will occur if it is used
//         to solve a system of equations.
//
{
  int i;
  int ii;
  int info;
  int j;
  int jj;
  int jp;
  double temp;
//
//  Test the input parameters.
//
  info = 0;

  if ( m < 0 )
  {
    return (-1);
  }
  else if ( n < 0 )
  {
    return (-2);
  }

  if ( m == 0 || n == 0 )
  {
    return 0;
  }

  for ( j = 1; j <= i4_min ( m, n ); j++ )
  {
//
//  Find the pivot.
//
    temp = r8_abs ( a[j-1+(j-1)*m] );
    jp = j;
    for ( i = j+1; i <= m; i++ )
    {
      if ( temp < r8_abs ( a[i-1+(j-1)*m] ) )
      {
        temp = r8_abs ( a[i-1+(j-1)*m] );
        jp = i;
      }
    }

    pivot[j-1] = jp;
//
//  Apply the interchange to columns 1:N.
//  Compute elements J+1:M of the J-th column.
//
    if ( a[jp-1+(j-1)*m] != 0.0 )
    {
      if ( jp != j )
      {
        for ( jj = 1; jj <= n; jj++ )
        {
          temp             = a[j-1+(jj-1)*m];
          a[j-1+(jj-1)*m]  = a[jp-1+(jj-1)*m];
          a[jp-1+(jj-1)*m] = temp;
        }
      }

      if ( j < m )
      {
        for ( i = j+1; i <= m; i++ )
        {
          a[i-1+(j-1)*m] = a[i-1+(j-1)*m] / a[j-1+(j-1)*m];
        }
      }
    }
    else if ( info == 0 )
    {
      info = j;
    }
//
//  Update the trailing submatrix.
//
    if ( j < i4_min ( m, n ) )
    {
      for ( ii = j+1; ii <= m; ii++ )
      {
        for ( i = j+1; i <= n; i++ )
        {
          a[ii-1+(i-1)*m] = a[ii-1+(i-1)*m] - a[ii-1+(j-1)*m] * a[j-1+(i-1)*m];
        }
      }
    }
  }

  return info;
}
//****************************************************************************80

double *r8ge_trs ( int n, int nrhs, char trans, double a[], int pivot[], 
  double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_TRS solves a system of linear equations factored by R8GE_TRF.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_TRS is a standalone version of the LAPACK routine SGETRS.
//
//    R8GE_TRS solves a system of linear equations
//      A * x = b  or  A' * X = B
//    with a general N by N matrix A using the PLU factorization computed
//    by R8GE_TRF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 November 2003
//
//  Author:
//
//    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
//    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
//    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
//    Sven Hammarling, Alan McKenney, Danny Sorensen,
//    LAPACK User's Guide,
//    Second Edition,
//    SIAM, 1995.
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.  0 <= N.
//
//    Input, int NRHS, the number of right hand sides.  0 <= NRHS.
//
//    Input, char TRANS, specifies the form of the system of equations:
//    'N':  A * x = b  (No transpose)
//    'T':  A'* X = B  (Transpose)
//    'C':  A'* X = B  (Conjugate transpose = Transpose)
//
//    Input, double A[N*N], the factors L and U from the factorization
//    A = P*L*U as computed by R8GE_TRF.
//
//    Input, int PIVOT[N], the pivot indices from R8GE_TRF.
//
//    Input, double B[N*NRHS], the right hand side matrix.
//
//    Output, double R8GE_TRS[N*NRHS], the solution matrix X.
//
{
  int i;
  int j;
  int k;
  double temp;
  double *x;

  if ( trans != 'n' && trans != 'N' && 
       trans != 't' && trans != 'T' && 
       trans != 'c' && trans != 'C' )
  {
    return NULL;
  }

  if ( n < 0 )
  {
    return NULL;
  }

  if ( nrhs < 0 )
  {
    return NULL;
  }

  if ( n == 0 )
  {
    return NULL;
  }
  if ( nrhs == 0 )
  {
    return NULL;
  }

  x = new double[n*nrhs];
  for ( k = 0; k < nrhs; k++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = b[i+k*n];
    }
  }

  if ( trans == 'n' || trans == 'N' )
  {
//
//  Apply row interchanges to the right hand sides.
//
    for ( i = 1; i <= n; i++ )
    {
      if ( pivot[i-1] != i )
      {
        for ( k = 0; k < nrhs; k++ )
        {
          temp                = x[i-1+k*n];
          x[i-1+k*n]          = x[pivot[i-1]-1+k*n];
          x[pivot[i-1]-1+k*n] = temp;
        }
      }
    }
//
//  Solve L * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n-1; j++ )
      {
        for ( i = j+1; i <= n; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[i-1+(j-1)*n] * x[j-1+k*n];
        }
      }
    }
//
//  Solve U * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 1 <= j; j-- )
      {
        x[j-1+k*n] = x[j-1+k*n] / a[j-1+(j-1)*n];
        for ( i = 1; i < j; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[i-1+(j-1)*n] * x[j-1+k*n];
        }
      }
    }
  }
  else
  {
//
//  Solve U' * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n; j++ )
      {
        x[j-1+k*n] = x[j-1+k*n] / a[j-1+(j-1)*n];
        for ( i = j+1; i <= n; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[j-1+(i-1)*n] * x[j-1+k*n];
        }
      }
    }
//
//  Solve L' * x = b, overwriting b with x.
//
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 2 <= j; j-- )
      {
        for ( i = 1; i < j; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[j-1+(i-1)*n] * x[j-1+k*n];
        }
      }
    }
//
//  Apply row interchanges to the solution vectors.
//
    for ( i = n; 1 <= i; i-- )
    {
      if ( pivot[i-1] != i )
      {
        for ( k = 0; k < nrhs; k++ )
        {
          temp                = x[i-1+k*n];
          x[i-1+k*n]          = x[pivot[i-1]-1+k*n];
          x[pivot[i-1]-1+k*n] = temp;
        }
      }
    }
  }

  return x;
}
//****************************************************************************80

double *r8ge_vxm ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_VXM multiplies a vector times a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2003
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
//    Input, double A[M*N], the R8GE matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8GE_VXM[N], the product A' * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < m; j++ )
    {
      b[i] = b[i] + a[j+i*m] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

double *r8ge_zero ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_ZERO zeros a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2003
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
//    Output, double R8GE_ZERO[M*N], the M by N matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double r8lt_det ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_DET computes the determinant of a R8LT matrix.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N*N], the R8LT matrix.
//
//    Output, double R8LT_DET, the determinant of the matrix.
//
{
  double det;
  int i;

  det = 1.0;
  for ( i = 0; i < n; i++ )
  {
    det = det * a[i+i*n];
  }

  return det;
}
//****************************************************************************80

double *r8lt_indicator ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_INDICATOR sets up a R8LT indicator matrix.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//    M and N must be positive.
//
//    Output, double R8LT_INDICATOR[M*N], the R8LT matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[m*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= i4_min ( i, n ); j++ )
    {
      a[i-1+(j-1)*m] = ( double ) ( fac * i + j );
    }
    for ( j = i+1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8lt_inverse ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_INVERSE computes the inverse of a R8LT matrix.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Second edition,
//    Academic Press, 1978,
//    ISBN 0-12-519260-6
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the R8LT matrix.
//
//    Output, double R8LT_INVERSE[N*N], the inverse of the matrix.
//
{
  double *b;
  int i;
  int j;
  int k;
  double t;
//
//  Check.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[i+i*n] == 0.0 )
    {
      cerr << "\n";
      cerr << "R8LT_INVERSE - Fatal error!\n";
      cerr << "  Zero diagonal element.\n";
      exit ( 1 );
    }
  }

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i < j )
      {
        b[i+j*n] = 0.0;
      }
      else if ( i == j )
      {
        b[i+j*n] = 1.0 / b[i+j*n];
      }
      else if ( j < i )
      {
        t = 0.0;
        for ( k = j; k <= i-1; k++ )
        {
          t = t - b[i+k*n] * b[k+j*n];
        }
        b[i+j*n] = t / b[i+i*n];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8lt_mxm ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_MXM multiplies two R8LT matrices.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, double A[N*N], B[N*N], the R8LT factor matrices.
//
//    Output, double R8LT_MXM[N*N], the R8LT product matrix.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      c[i+j*n] = 0.0;
      for ( k = j; k <= i; k++ )
      {
        c[i+j*n] = c[i+j*n] + a[i+k*n] * b[k+j*n];
      }
    }
  }

  return c;
}
//****************************************************************************80

double *r8lt_mxv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_MXV multiplies a R8LT matrix times a vector.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2003
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
//    Input, double A[M*N], the R8LT matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8LT_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;
  int jmax;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    jmax = i4_min ( i, n-1 );
    for ( j = 0; j <= jmax; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void r8lt_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_PRINT prints a R8LT matrix.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[M*N], the R8LT matrix.
//
//    Input, string TITLE, a title.
//
{
  r8lt_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8lt_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_PRINT_SOME prints some of a R8LT matrix.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[M*N], the R8LT matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( ilo < jlo )
  {
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo );

    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(4) << i << "  ";

      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i < j )
        {
          cout << "              ";
        }
        else
        {
          cout << setw(12) << a[i-1+(j-1)*m] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8lt_random ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_RANDOM randomizes a R8LT matrix.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//    M and N must be positive.
//
//    Input/output, int SEED, a seed for the random number generator.
//
//    Output, double R8LT_RANDOM[M*N], the R8LT matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j <= i4_min ( i, n-1); j++ )
    {
      a[i+j*m] = r8_uniform_01 ( seed );
    }
    for ( j = i+1; j < n; j++ )
    {
      a[i+j*m] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8lt_sl ( int n, double a[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_SL solves a R8LT system.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//    No factorization of the lower triangular matrix is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the R8LT matrix.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB, is 0 to solve the untransposed system,
//    nonzero to solve the transposed system.
//
//    Output, double R8LT_SL[N], the solution vector.
//
{
  int i;
  int j;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = x[j] / a[j+j*n];
      for ( i = j+1; i < n; i++ )
      {
        x[i] = x[i] - a[i+j*n] * x[j];
      }
    }
  }
  else
  {
    for ( j = n-1; 0 <= j; j-- )
    {
      x[j] = x[j] / a[j+j*n];
      for ( i = 0; i < j; i++ )
      {
        x[i] = x[i] - a[j+i*n] * x[j];
      }
    }
  }

  return x;
}
//****************************************************************************80

double *r8lt_vxm ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_VXM multiplies a vector times a R8LT matrix.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2003
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
//    Input, double A[M*N], the R8LT matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8LT_VXM[N], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( j = 0; j < n; j++)
  {
    b[j] = 0.0;
    for ( i = j; i < m; i++ )
    {
      b[j] = b[j] + x[i] * a[i+j*m];
    }
  }

  return b;
}
//****************************************************************************80

double *r8lt_zero ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8LT_ZERO zeros a R8LT matrix.
//
//  Discussion:
//
//    The R8LT storage format is used for an M by N lower triangular 
//    matrix A, and allocates storage even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2004
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
//    Output, double R8LT_ZERO[M*N], the R8LT matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
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
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
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
//    20 August 2010
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
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

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
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
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
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

double *r8ncf_indicator ( int m, int n, int nz_num, int rowcol[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8NCF_INDICATOR sets up a R8NCF indicator matrix.
//
//  Discussion:
//
//    The R8NCF storage format stores NZ_NUM, the number of nonzeros,
//    a real array containing the nonzero values, a 2 by NZ_NUM integer
//    array storing the row and column of each nonzero entry.
//
//    The R8NCF format is used by NSPCG.  NSPCG requires that the information
//    for the diagonal entries of the matrix must come first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, int NZ_NUM, the number of nonzero entries.
//
//    Input, int ROWCOL[2*NZ_NUM], the coordinates of the nonzero entries.
//
//    Output, double A[NZ_NUM], the indicator matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;
  int k;

  a = new double[nz_num];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( k = 0; k < nz_num; k++ )
  {
    i = rowcol[0+k*2];
    j = rowcol[1+k*2];
    a[k] = ( double ) ( fac * i + j );
  }

  return a;
}
//****************************************************************************80

void r8ncf_print ( int m, int n, int nz_num, int rowcol[], 
  double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8NCF_PRINT prints a R8NCF matrix.
//
//  Discussion:
//
//    The R8NCF storage format stores NZ_NUM, the number of nonzeros,
//    a real array containing the nonzero values, a 2 by NZ_NUM integer
//    array storing the row and column of each nonzero entry.
//
//    The R8NCF format is used by NSPCG.  NSPCG requires that the information
//    for the diagonal entries of the matrix must come first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROWCOL[2*NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
//    Input, string TITLE, a title.
//
{
  r8ncf_print_some ( m, n, nz_num, rowcol, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8ncf_print_some ( int m, int n, int nz_num, int rowcol[], 
  double a[], int ilo, int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8NCF_PRINT_SOME prints some of a R8NCF matrix.
//
//  Discussion:
//
//    The R8NCF storage format stores NZ_NUM, the number of nonzeros,
//    a real array containing the nonzero values, a 2 by NZ_NUM integer
//    array storing the row and column of each nonzero entry.
//
//    The R8NCF format is used by NSPCG.  NSPCG requires that the information
//    for the diagonal entries of the matrix must come first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROWCOL[2*NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;
  int k;
  bool nonzero;
  double temp[INCX];

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";

    cout << "  Col:  ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      nonzero = false;
      for ( j2 = 1; j2 <= INCX; j2++ )
      {
        temp[j2-1] = 0.0;
      }

      for ( k = 1; k <= nz_num; k++ )
      {
        if ( i == rowcol[0+(k-1)*2] && 
          j2lo <= rowcol[1+(k-1)*2] && 
          rowcol[1+(k-1)*2] <= j2hi )
        {
          j2 = rowcol[1+(k-1)*2] - j2lo + 1;

          if ( a[k-1] == 0.0 )
          {
            continue;
          }

          nonzero = true;
          temp[j2-1] = a[k-1];
        }
      }

      if ( nonzero )
      {
        cout << setw(6) << i;
        for ( j2 = 1; j2 <= inc; j2++ )
        {
          cout << setw(12) << temp[j2-1] << "  ";
        }
        cout << "\n";
      }
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double r8pbl_det ( int n, int mu, double a_lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBL_DET computes the determinant of a matrix factored by R8PBL_FA.
//
//  Discussion:
//
//    The R8PBL storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and lower triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row 1 of the array.
//    The first subdiagonal in row 2, columns 1 through MU.
//    The second subdiagonal in row 3, columns 1 through MU-1.
//    The MU-th subdiagonal in row MU+1, columns 1 through 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the upper (and lower) bandwidth.
//    MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(MU+1)*N], the LU factors from R8PBL_FA.
//
//    Output, double R8PBL_DET, the determinant of the matrix.
//
{
  double det;
  int j;

  det = 1.0;

  for ( j = 0; j < n; j++ )
  {
    det = det * a_lu[0+j*(mu+1)] * a_lu[0+j*(mu+1)];
  }

  return det;
}
//****************************************************************************80

double *r8pbl_indicator ( int n, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBL_INDICATOR sets up a R8PBL indicator matrix.
//
//  Discussion:
//
//    The R8PBL storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and lower triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row 1 of the array.
//    The first subdiagonal in row 2, columns 1 through MU.
//    The second subdiagonal in row 3, columns 1 through MU-1.
//    The MU-th subdiagonal in row MU+1, columns 1 through 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of subdiagonals in the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Output, double R8PBL_INDICATOR[(MU+1)*N], the R8PBL matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[(mu+1)*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );
//
//  Zero out the "junk" entries.
//
  for ( j = n - mu + 1; j <= n; j++ )
  {
    for ( i = n + 1; i <= j + mu; i++ )
    {
      a[i-j+(j-1)*(mu+1)] = 0.0;
    }
  }
//
//  Set the meaningful values.
//
  for ( i = 0; i <= n; i++ )
  {
    for ( j = i4_max ( 1, i - mu ); j <= i; j++ )
    {
      a[i-j+(j-1)*(mu+1)] = ( double ) ( fac * i + j );
    }
  }
  return a;
}
//****************************************************************************80

void r8pbl_print ( int n, int mu, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBL_PRINT prints a R8PBL matrix.
//
//  Discussion:
//
//    The R8PBL storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and lower triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row 1 of the array.
//    The first subdiagonal in row 2, columns 1 through MU.
//    The second subdiagonal in row 3, columns 1 through MU-1.
//    The MU-th subdiagonal in row MU+1, columns 1 through 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the upper (and lower) bandwidth.
//    MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBL matrix.
//
//    Input, string TITLE, a title.
//
{
  r8pbl_print_some ( n, mu, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8pbl_print_some ( int n, int mu, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBL_PRINT_SOME prints some of a R8PBL matrix.
//
//  Discussion:
//
//    The R8PBL storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and lower triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row 1 of the array.
//    The first subdiagonal in row 2, columns 1 through MU.
//    The second subdiagonal in row 3, columns 1 through MU-1.
//    The MU-th subdiagonal in row MU+1, columns 1 through 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the upper (and lower) bandwidth.
//    MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBL matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - mu );

    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + mu );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j && j <= i + mu )
        {
          cout << setw(12) << a[j-i+(i-1)*(mu+1)] << "  ";
        }
        else if ( j <= i && i <= j + mu )
        {
          cout << setw(12) << a[i-j+(j-1)*(mu+1)] << "  ";
        }
        else
        {
          cout << "              ";
        }
      }
      cout << "\n";
    }
  }
  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8pbl_random ( int n, int mu, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBL_RANDOM randomizes a R8PBL matrix.
//
//  Discussion:
//
//    The R8PBL storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and lower triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row 1 of the array.
//    The first subdiagonal in row 2, columns 1 through MU.
//    The second subdiagonal in row 3, columns 1 through MU-1.
//    The MU-th subdiagonal in row MU+1, columns 1 through 1.
//
//    The matrix returned will be positive definite, but of limited
//    randomness.  The off diagonal elements are random values between
//    0 and 1, and the diagonal element of each row is selected to
//    ensure strict diagonal dominance.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of subdiagonals in the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8PBL_RANDOM[(MU+1)*N], the R8PBL matrix.
//
{
  double *a;
  int i;
  int j;
  double r;
  double sum2;

  a = new double[(mu+1)*n];
//
//  Zero out the "junk" entries.
//
  for ( j = n - mu; j <= n-1; j++ )
  {
    for ( i = n-j; i <= mu; i++ )
    {
      a[i+j*(mu+1)] = 0.0;
    }
  }
//
//  Set the off diagonal values.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = i4_max ( 0, i - mu ); j <= i - 1; j++ )
    {
      a[i-j+j*(mu+1)] = r8_uniform_01 ( seed );
    }
  }
//
//  Set the diagonal values.
//
  for ( i = 0; i < n; i++ )
  {
    sum2 = 0.0;

    for ( j = i4_max ( 0, i - mu ); j <= i-1; j++ )
    {
      sum2 = sum2 + r8_abs ( a[i-j+j*(mu+1)] );
    }

    for ( j = i+1; j <= i4_min ( i + mu, n -1 ); j++ )
    {
      sum2 = sum2 + r8_abs ( a[j-i+i*(mu+1)] );
    }

    r = r8_uniform_01 ( seed );

    a[0+i*(mu+1)] = ( 1.0 + r ) * ( sum2 + 0.01 );
  }

  return a;
}
//****************************************************************************80

double *r8pbl_to_r8ge ( int n, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBL_TO_R8GE copies a R8PBL matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8PBL storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and lower triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row 1 of the array.
//    The first subdiagonal in row 2, columns 1 through MU.
//    The second subdiagonal in row 3, columns 1 through MU-1.
//    The MU-th subdiagonal in row MU+1, columns 1 through 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, int MU, the upper bandwidth of A1.
//    MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBL matrix.
//
//    Output, double R8PBL_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i <= j && j <= i + mu )
      {
        b[i+j*n] = a[j-i+i*(mu+1)];
      }
      else if ( i - mu <= j && j < i )
      {
        b[i+j*n] = a[i-j+j*(mu+1)];
      }
      else
      {
        b[i+j*n] = 0.0;
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8pbl_zero ( int n, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBL_ZERO zeros a R8PBL matrix.
//
//  Discussion:
//
//    The R8PBL storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and lower triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row 1 of the array.
//    The first subdiagonal in row 2, columns 1 through MU.
//    The second subdiagonal in row 3, columns 1 through MU-1.
//    The MU-th subdiagonal in row MU+1, columns 1 through 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of subdiagonals in the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Output, double R8PBL_ZERO[(MU+1)*N], the R8PBL matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[(mu+1)*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < mu+1; i++ )
    {
      a[i+j*(mu+1)] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8pbu_cg ( int n, int mu, double a[], double b[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_CG uses the conjugate gradient method on a R8PBU system.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//    The matrix A must be a positive definite symmetric band matrix.
//
//    The method is designed to reach the solution after N computational
//    steps.  However, roundoff may introduce unacceptably large errors for
//    some problems.  In such a case, calling the routine again, using
//    the computed solution as the new starting estimate, should improve
//    the results.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Frank Beckman,
//    The Solution of Linear Equations by the Conjugate Gradient Method,
//    in Mathematical Methods for Digital Computers,
//    edited by John Ralston, Herbert Wilf,
//    Wiley, 1967,
//    ISBN: 0471706892,
//    LC: QA76.5.R3.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals.
//    MU must be at least 0, and no more than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBU matrix.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, double X[N], an estimate for the solution.
//
//    Output, double R8PBU_CG[N], the approximate solution vector.
//
{
  double alpha;
  double *ap;
  double beta;
  int i;
  int it;
  double *p;
  double pap;
  double pr;
  double *r;
  double rap;
  double *x2;

  ap = new double[n];
  p = new double[n];
  r = new double[n];
  x2 = new double[n];
//
//  Initialize
//    AP = A * x,
//    R  = b - A * x,
//    P  = b - A * x.
//
  for ( i = 0; i < n; i++ )
  {
    x2[i] = x[i];
  }

  ap = r8pbu_mxv ( n, mu, a, x2 );
  for ( i = 0; i < n; i++ )
  {
    r[i] = b[i] - ap[i];
  }

  for ( i = 0; i < n; i++ )
  {
    p[i] = b[i] - ap[i];
  }
//
//  Do the N steps of the conjugate gradient method.
//
  for ( it = 1; it <= n; it++ )
  {
//
//  Compute the matrix*vector product AP=A*P.
//
    delete [] ap;

    ap = r8pbu_mxv ( n, mu, a, p );
//
//  Compute the dot products
//    PAP = P*AP,
//    PR  = P*R
//  Set
//    ALPHA = PR / PAP.
//
    pap = 0.0;
    for ( i = 0; i < n; i++ )
    {
      pap = pap + p[i] * ap[i];
    }
    if ( pap == 0.0 )
    {
      break;
    }

    pr = 0.0;
    for ( i = 0; i < n; i++ )
    {
      pr = pr + p[i] * r[i];
    }
    alpha = pr / pap;
//
//  Set
//    X = X + ALPHA * P
//    R = R - ALPHA * AP.
//
    for ( i = 0; i < n; i++ )
    {
      x2[i] = x2[i] + alpha * p[i];
    }

    for ( i = 0; i < n; i++ )
    {
      r[i] = r[i] - alpha * ap[i];
    }
//
//  Compute the vector dot product
//    RAP = R*AP
//  Set
//    BETA = - RAP / PAP.
//
    rap = 0.0;
    for ( i = 0; i < n; i++ )
    {
      rap = rap + r[i] * ap[i];
    }
    beta = - rap / pap;
//
//  Update the perturbation vector
//    P = R + BETA * P.
//
    for ( i = 0; i < n; i++ )
    {
      p[i] = r[i] + beta * p[i];
    }
  }

  delete [] ap;
  delete [] p;
  delete [] r;

  return x2;
}
//****************************************************************************80

double r8pbu_det ( int n, int mu, double a_lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_DET computes the determinant of a matrix factored by R8PBU_FA.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals of the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Input, double A_LU[(MU+1)*N], the LU factors from R8PBU_FA.
//
//    Output, double R8PBU_DET, the determinant of the matrix.
//
{
  double det;
  int j;

  det = 1.0;

  for ( j = 0; j < n; j++ )
  {
    det = det * a_lu[mu+j*(mu+1)] * a_lu[mu+j*(mu+1)];
  }

  return det;
}
//****************************************************************************80

double *r8pbu_fa ( int n, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_FA factors a R8PBU matrix.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//    The matrix A must be a positive definite symmetric band matrix.
//
//    Once factored, linear systems A*x=b involving the matrix can be solved
//    by calling R8PBU_SL.  No pivoting is performed.  Pivoting is not necessary
//    for positive definite symmetric matrices.  If the matrix is not positive
//    definite, the algorithm may behave correctly, but it is also possible
//    that an illegal divide by zero will occur.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals of the matrix.
//    MU must be at least 0, and no more than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBU matrix.
//
//    Output, double R8PBU_FA[(MU+1)*N], information describing a factored
//    form of the matrix.
//
{
  double *b;
  int i;
  int ik;
  int j;
  int jk;
  int k;
  int mm;
  double s;
  double t;

  b = new double[(mu+1)*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < mu+1; i++ )
    {
      b[i+j*(mu+1)] = a[i+j*(mu+1)];
    }
  }

  for ( j = 1; j <= n; j++ )
  {

    ik = mu + 1;
    jk = i4_max ( j - mu, 1 );
    mm = i4_max ( mu + 2 - j, 1 );

    s = 0.0;

    for ( k = mm; k <= mu; k++ )
    {
      t = 0.0;
      for ( i = 0; i <= k-mm-1; i++ )
      {
        t = t + b[ik+i-1+(jk-1)*(mu+1)] * b[mm+i-1+(j-1)*(mu+1)];
      }
      b[k-1+(j-1)*(mu+1)] = ( b[k-1+(j-1)*(mu+1)] - t ) /
        b[mu+(jk-1)*(mu+1)];

      s = s + b[k-1+(j-1)*(mu+1)] * b[k-1+(j-1)*(mu+1)];
      ik = ik - 1;
      jk = jk + 1;
    }

    s = b[mu+(j-1)*(mu+1)] - s;

    if ( s <= 0.0 )
    {
      return NULL;
    }

    b[mu+(j-1)*(mu+1)] = sqrt ( s );
  }

  return b;
}
//****************************************************************************80

double *r8pbu_indicator ( int n, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_INDICATOR sets up a R8PBU indicator matrix.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals in the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Output, double R8PBU_INDICATOR[(MU+1)*N], the R8PBU matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[(mu+1)*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );
//
//  Zero out the "junk" entries.
//
  for ( j = 0; j < mu; j++ )
  {
    for ( i = 0; i <= mu - j; i++ )
    {
      a[i+j*(mu+1)] = 0.0;
    }
  }
//
//  Set the meaningful values.
//
  for ( i = 1; i <= n; i++ )
  {
    for ( j = i; j <= i4_min ( i + mu, n ); j++ )
    {
      a[mu+i-j+(j-1)*(mu+1)] = ( double ) ( fac * i + j );
    }
  }
  return a;
}
//****************************************************************************80

double *r8pbu_ml ( int n, int mu, double a_lu[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_ML multiplies a vector times a matrix that was factored by R8PBU_FA.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals of the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Input, double A_LU[(MU+1)*N], the LU factors from R8PBU_FA.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8PBU_ML[N], the product A * x.
//
{
  double *b;
  int i;
  int ilo;
  int j;
  int jhi;
  int k;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }
//
//  Multiply U * X = Y.
//
  for ( k = 1; k <= n; k++ )
  {
    ilo = i4_max ( 1, k - mu );
    for ( i = ilo; i <= k - 1; i++ )
    {
      b[i-1] = b[i-1] + a_lu[mu+i-k+(k-1)*(mu+1)] * b[k-1];
    }
    b[k-1] = a_lu[mu+(k-1)*(mu+1)] * b[k-1];
  }
//
//  Multiply L * Y = B.
//
  for ( k = n; 1 <= k; k-- )
  {
    jhi = i4_min ( k + mu, n );
    for ( j = k+1; j <= jhi; j++ )
    {
      b[j-1] = b[j-1] + a_lu[mu+k-j+(j-1)*(mu+1)] * b[k-1];
    }

    b[k-1] = a_lu[mu+(k-1)*(mu+1)] * b[k-1];
  }

  return b;
}
//****************************************************************************80

double *r8pbu_mxv ( int n, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_MXV multiplies a R8PBU matrix times a vector.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals in the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBU matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8PBU_MXV[N], the result vector A * x.
//
{
  double *b;
  int i;
  int ieqn;
  int j;

  b = new double[n];
//
//  Multiply X by the diagonal of the matrix.
//
  for ( j = 0; j < n; j++ )
  {
    b[j] = a[mu+j*(mu+1)] * x[j];
  }
//
//  Multiply X by the superdiagonals of the matrix.
//
  for ( i = mu; 1 <= i; i-- )
  {
    for ( j = mu+2-i; j <= n; j++ )
    {
      ieqn = i + j - mu - 1;
      b[ieqn-1] = b[ieqn-1] + a[i-1+(j-1)*(mu+1)] * x[j-1];
      b[j-1] = b[j-1] + a[i-1+(j-1)*(mu+1)] * x[ieqn-1];
    }
  }

  return b;
}
//****************************************************************************80

void r8pbu_print ( int n, int mu, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_PRINT prints a R8PBU matrix.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the upper (and lower) bandwidth.
//    MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBU matrix.
//
//    Input, string TITLE, a title.
//
{
  r8pbu_print_some ( n, mu, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8pbu_print_some ( int n, int mu, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_PRINT_SOME prints some of a R8PBU matrix.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the upper (and lower) bandwidth.
//    MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBU matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";

    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - mu );
    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + mu );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( mu < i-j || mu < j-i ) 
        {
          cout << "              ";
        }
        else if ( i <= j && j <= i + mu )
        {
          cout << setw(12) << a[mu+i-j+(j-1)*(mu+1)] << "  ";
        }
        else if ( i - mu <= j && j <= i )
        {
          cout << setw(12) << a[mu+j-i+(i-1)*(mu+1)] << "  ";
        }
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8pbu_random ( int n, int mu, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_RANDOM randomizes a R8PBU matrix.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//    The matrix returned will be positive definite, but of limited
//    randomness.  The off diagonal elements are random values between
//    0 and 1, and the diagonal element of each row is selected to
//    ensure strict diagonal dominance.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals in the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8PBU_RANDOM[(MU+1)*N], the R8PBU matrix.
//
{
  double *a;
  int i;
  int j;
  int jhi;
  int jlo;
  double r;
  double sum2;

  a = new double[(mu+1)*n];
//
//  Zero out the "junk" entries.
//
  for ( j = 0; j < mu; j++ )
  {
    for ( i = 0; i <= mu - j; i++ )
    {
      a[i+j*(mu+1)] = 0.0;
    }
  }
//
//  Set the off diagonal values.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = i+1; j <= i4_min ( i + mu, n-1 ); j++ )
    {
      a[mu+i-j+j*(mu+1)] = r8_uniform_01 ( seed );
    }
  }
//
//  Set the diagonal values.
//
  for ( i = 1; i <= n; i++ )
  {
    sum2 = 0.0;

    jlo = i4_max ( 1, i - mu );
    for ( j = jlo; j <= i-1; j++ )
    {
      sum2 = sum2 + r8_abs ( a[(mu+j-i)+(i-1)*(mu+1)] );
    }

    jhi = i4_min ( i + mu, n );
    for ( j = i+1; j <= jhi; j++ )
    {
      sum2 = sum2 + r8_abs ( a[mu+i-j+(j-1)*(mu+1)] );
    }

    r = r8_uniform_01 ( seed );

    a[mu+(i-1)*(mu+1)] = ( 1.0 + r ) * ( sum2 + 0.01 );

  }

  return a;
}
//****************************************************************************80

double *r8pbu_sl ( int n, int mu, double a_lu[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_SL solves a R8PBU system factored by R8PBU_FA.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 October 2003
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals of the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Input, double A_LU[(MU+1)*N], the LU factors from R8PBU_FA.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double R8PBU_SL[N], the solution vector.
//
{
  int i;
  int ilo;
  int k;
  double *x;

  x = new double[n];
//
//  Solve L * Y = B.
//
  for ( k = 0; k < n; k++ )
  {
    ilo = i4_max ( 0, k - mu );

    x[k] = b[k];
    for ( i = ilo; i <= k-1; i++ )
    {
      x[k] = x[k] - x[i] * a_lu[mu+i-k+k*(mu+1)];
    }
    x[k] = x[k] / a_lu[mu+k*(mu+1)];
  }
//
//  Solve U * X = Y.
//
  for ( k = n-1; 0 < k; k-- )
  {
    x[k] = x[k] / a_lu[mu+k*(mu+1)];

    ilo = i4_max ( 0, k - mu );
    for ( i = ilo; i <= k-1; i++ )
    {
      x[i] = x[i] - x[k] * a_lu[mu+i-k+k*(mu+1)];
    }
  }

  return x;
}
//****************************************************************************80

double *r8pbu_sor ( int n, int mu, double a[], double b[], double eps, int itchk, 
  int itmax, double omega, double x_init[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_SOR uses SOR iteration to solve a R8PBU linear system.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//    The matrix A must be a positive definite symmetric band matrix.
//
//    A relaxation factor OMEGA may be used.
//
//    The iteration will proceed until a convergence test is met,
//    or the iteration limit is reached.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals in the matrix.
//    MU must be at least 0, and no more than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBU matrix.
//
//    Input, double B[N], the right hand side of the system.
//
//    Input, double EPS, convergence tolerance for the system.  The vector
//    b - A * x is computed every ITCHK iterations, and if the maximum
//    entry of this vector is of norm less than EPS, the program
//    will return.
//
//    Input, int ITCHK, the interval between convergence checks.  ITCHK steps
//    will be taken before any check is made on whether the iteration
//    has converged.  ITCHK should be at least 1 and no greater
//    than ITMAX.
//
//    Input, int ITMAX, the maximum number of iterations allowed.  The
//    program will return to the user if this many iterations are taken
//    without convergence.
//
//    Input, double OMEGA, the relaxation factor.  OMEGA must be strictly between
//    0 and 2.  Use OMEGA = 1 for no relaxation, classical Jacobi iteration.
//
//    Input, double X_INIT[N], a starting vector for the iteration.
//
//    Output, double R8PBU_SOR[N], the approximation to the solution.
//
{
  double err;
  int i;
  int it;
  int itknt;
  double *x;
  double *xtemp;

  if ( itchk <= 0 || itmax < itchk )
  {
    cerr << "\n";
    cerr << "R8PBU_SOR - Fatal error!\n";
    cerr << "  Illegal ITCHK = " << itchk << "\n";
    exit ( 1 );
  }

  if ( itmax <= 0 )
  {
    cerr << "\n";
    cerr << "R8PBU_SOR - Fatal error!\n";
    cerr << "  Nonpositive ITMAX = " << itmax << "\n";
    exit ( 1 );
  }

  if ( omega <= 0.0 || 2.0 <= omega )
  {
    cerr << "\n";
    cerr << "R8PBU_SOR - Fatal error!\n";
    cerr << "  Illegal value of OMEGA = " << omega << "\n";
    exit ( 1 );
  }

  itknt = 0;

  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = x_init[i];
  }
//
//  Take ITCHK steps of the iteration before doing a convergence check.
//
  while ( itknt <= itmax )
  {
    for ( it = 1; it <= itchk; it++ )
    {
//
//  Compute XTEMP(I) = B(I) + A(I,I) * X(I) - SUM ( J=1 to N ) A(I,J) * X(J).
//
      xtemp = r8pbu_mxv ( n, mu, a, x );

      for ( i = 0; i < n; i++ )
      {
        xtemp[i] = x[i] + ( b[i] - xtemp[i] ) / a[mu+i*(mu+1)];
      }
//
//  Compute the next iterate as a weighted combination of the
//  old iterate and the just computed standard Jacobi iterate.
//
      if ( omega != 1.0 )
      {
        for ( i = 0; i < n; i++ )
        {
          xtemp[i] = ( 1.0 - omega ) * x[i] + omega * xtemp[i];
        }
      }
//
//  Copy the new result into the old result vector.
//
      for ( i = 0; i < n; i++ )
      {
        x[i] = xtemp[i];
      }
    }
    delete [] xtemp;
//
//  Compute the maximum residual, the greatest entry in the vector
//  RESID(I) = B(I) - A(I,J) * X(J).
//
    xtemp = r8pbu_mxv ( n, mu, a, x );

    err = 0.0;
    for ( i = 0; i < n; i++ )
    {
      err = r8_max ( err, r8_abs ( b[i] - xtemp[i] ) );
    }
    delete [] xtemp;
//
//  Test to see if we can quit because of convergence,
//
    if ( err <= eps )
    {
      return x;
    }

  }

  cout << "\n";
  cout << "R8PBU_SOR - Warning!\n";
  cout << "  The iteration did not converge.\n";

  return x;
}
//****************************************************************************80

double *r8pbu_to_r8ge ( int n, int mu, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_TO_R8GE copies a R8PBU matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, int MU, the upper bandwidth of A1.
//    MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(MU+1)*N], the R8PBU matrix.
//
//    Output, double R8PBU_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i <= j && j <= i+mu )
      {
        b[i+j*n] = a[mu+i-j+j*(mu+1)];
      }
      else if ( i-mu <= j && j < i )
      {
        b[i+j*n] = a[mu+j-i+i*(mu+1)];
      }
      else
      {
        b[i+j*n] = 0.0;
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8pbu_zero ( int n, int mu )

//****************************************************************************80
//
//  Purpose:
//
//    R8PBU_ZERO zeros a R8PBU matrix.
//
//  Discussion:
//
//    The R8PBU storage format is used for a symmetric positive definite band matrix.
//
//    To save storage, only the diagonal and upper triangle of A is stored,
//    in a compact diagonal format that preserves columns.
//
//    The diagonal is stored in row MU+1 of the array.
//    The first superdiagonal in row MU, columns 2 through N.
//    The second superdiagonal in row MU-1, columns 3 through N.
//    The MU-th superdiagonal in row 1, columns MU+1 through N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int MU, the number of superdiagonals in the matrix.
//    MU must be at least 0 and no more than N-1.
//
//    Output, double R8PBU_ZERO[(MU+1)*N], the R8PBU matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[(mu+1)*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < mu+1; i++ )
    {
      a[i+j*(mu+1)] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double r8po_det ( int n, double a_lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_DET computes the determinant of a matrix factored by R8PO_FA.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A_LU[N*N], the LU factors from R8PO_FA.
//
//    Output, double R8PO_DET, the determinant of A.
//
{
  double det;
  int i;

  det = 1.0;

  for ( i = 0; i < n; i++ )
  {
    det = det * a_lu[i+i*n];
  }

  return det;
}
//****************************************************************************80

double *r8po_fa ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_FA factors a R8PO matrix.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//    The positive definite symmetric matrix A has a Cholesky factorization
//    of the form:
//
//      A = R' * R
//
//    where R is an upper triangular matrix with positive elements on
//    its diagonal.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix in R8PO storage.
//
//    Output, double R8PO_FA[N*N], the Cholesky factor in SGE
//    storage, or NULL if there was an error.
//
{
  double *b;
  int i;
  int j;
  int k;
  double s;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = 0; k <= j-1; k++ )
    {
      for ( i = 0; i <= k-1; i++ )
      {
        b[k+j*n] = b[k+j*n] - b[i+k*n] * b[i+j*n];
      }
      b[k+j*n] = b[k+j*n] / b[k+k*n];
    }

    s = b[j+j*n];
    for ( i = 0; i <= j-1; i++ )
    {
      s = s - b[i+j*n] * b[i+j*n];
    }

    if ( s <= 0.0 )
    {
      delete [] b;
      return NULL;
    }

    b[j+j*n] = sqrt ( s );
  }
//
//  Since the Cholesky factor is in R8GE format, zero out the lower triangle.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      b[i+j*n] = 0.0;
    }
  }

  return b;
}
//****************************************************************************80

double *r8po_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_INDICATOR sets up a R8PO indicator matrix.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//    N must be positive.
//
//    Output, double R8PO_INDICATOR[N*N], the R8PO matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[n*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= i-1; j++ )
    {
      a[i-1+(j-1)*n] = 0.0;
    }
    for ( j = i; j <= n; j++ )
    {
      a[i-1+(j-1)*n] = ( double ) ( fac * i + j );
    }
  }

  return a;
}
//****************************************************************************80

double *r8po_inverse ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_INVERSE computes the inverse of a matrix factored by R8PO_FA.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the Cholesky factor, in R8GE storage, returned by R8PO_FA.
//
//    Output, double R8PO_INVERSE[N*N], the inverse, in R8PO storage.
//
{
  double *b;
  int i;
  int j;
  int k;
  double t;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }
//
//  Compute Inverse ( R ).
//
  for ( k = 0; k < n; k++ )
  {
    b[k+k*n] = 1.0 / b[k+k*n];
    for ( i = 0; i < k; i++ )
    {
      b[i+k*n] = -b[i+k*n] * b[k+k*n];
    }

    for ( j = k+1; j < n; j++ )
    {
      t = b[k+j*n];
      b[k+j*n] = 0.0;
      for ( i = 0; i <= k; i++ )
      {
        b[i+j*n] = b[i+j*n] + t * b[i+k*n];
      }
    }
  }
//
//  Compute Inverse ( R ) * Transpose ( Inverse ( R ) ).
//
  for ( j = 0; j < n; j++ )
  {
    for ( k = 0; k < j; k++ )
    {
      t = b[k+j*n];
      for ( i = 0; i <= k; i++ )
      {
        b[i+k*n] = b[i+k*n] + t * b[i+j*n];
      }
    }
    t = b[j+j*n];
    for ( i = 0; i <= j; i++ )
    {
      b[i+j*n] = b[i+j*n] * t;
    }
  }

  return b;
}
//****************************************************************************80

double *r8po_ml ( int n, double a_lu[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_ML computes A * x = b after A has been factored by R8PO_FA.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A_LU[N*N], the Cholesky factor from R8PO_FA.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8PO_ML[N], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];
//
//  Compute R * x = y.
//
  for ( i = 0; i < n; i++ )
  {
    b[i] = a_lu[i+i*n] * x[i];
    for ( j = i+1; j < n; j++ )
    {
      b[i] = b[i] + a_lu[i+j*n] * x[j];
    }
  }
//
//  Compute R' * y = b.
//
  for ( j = n-1; 0 <= j; j-- )
  {
    b[j] = a_lu[j+j*n] * b[j];
    for ( i = 0; i < j; i++ )
    {
      b[j] = b[j] + b[i] * a_lu[i+j*n];
    }
  }

  return b;
}
//****************************************************************************80

double *r8po_mxm ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_MXM multiplies two R8PO matrices.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, double A[N*N], B[N*N], the R8PO factor matrices.
//
//    Output, double R8PO_MXM[N*N], the R8PO product matrix.
//
{
  double aik;
  double bkj;
  double *c;
  int i;
  int j;
  int k;

  c = new double[n*n];

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c[i-1+(j-1)*n] = 0.0;
    }
  }

  for ( i = 1; i <= n; i++ )
  {
    for ( j = i; j <= n; j++ )
    {
      for ( k = 1; k <= n; k++ )
      {
        if ( i <= k )
        {
          aik = a[i-1+(k-1)*n];
        }
        else
        {
          aik = a[k-1+(i-1)*n];
        }

        if ( k <= j )
        {
          bkj = b[k-1+(j-1)*n];
        }
        else
        {
          bkj = b[j-1+(k-1)*n];
        }

        c[i-1+(j-1)*n] = c[i-1+(j-1)*n] + aik * bkj;

      }
    }

  }

  return c;
}
//****************************************************************************80

double *r8po_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_MXV multiplies a R8PO matrix times a vector.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the R8PO matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8PO_MXV(N), the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < i; j++ )
    {
      b[i] = b[i] + a[j+i*n] * x[j];
    }
    for ( j = i; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*n] * x[j];
    }
  }
  return b;
}
//****************************************************************************80

void r8po_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_PRINT prints a R8PO matrix.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[M*N], the R8PO matrix.
//
//    Input, string TITLE, a title.
//
{
  r8po_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8po_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_PRINT_SOME prints some of a R8PO matrix.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[M*N], the R8PO matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j )
        {
          cout << setw(12) << a[i-1+(j-1)*n] << "  ";
        }
        else
        {
          cout << setw(12) << a[j-1+(i-1)*n] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8po_random ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_RANDOM randomizes a R8PO matrix.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//    The matrix computed here is not simply a set of random numbers in
//    the nonzero slots of the R8PO array.  It is also a positive definite
//    matrix.  It is computed by setting a "random" upper triangular
//    Cholesky factor R, and then computing A = R'*R.
//    The randomness is limited by the fact that all the entries of
//    R will be between 0 and 1.  A truly random R is only required
//    to have positive entries on the diagonal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8PO_RANDOM[N*N], the R8PO matrix.
//
{
  double *a;
  int i;
  int j;
  int k;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }

  for ( i = n; 1 <= i; i-- )
  {
//
//  Set row I of R.
//
    for ( j = i; j <= n; j++ )
    {
      a[i-1+(j-1)*n] = r8_uniform_01 ( seed );
    }
//
//  Consider element J of row I, last to first.
//
    for ( j = n; i <= j; j-- )
    {
//
//  Add multiples of row I to lower elements of column J.
//
      for ( k = i+1; k <= j; k++ )
      {
        a[k-1+(j-1)*n] = a[k-1+(j-1)*n] + a[i-1+(k-1)*n] * a[i-1+(j-1)*n];
      }
//
//  Reset element J.
//
      a[i-1+(j-1)*n] = a[i-1+(i-1)*n] * a[i-1+(j-1)*n];
    }
  }

  return a;
}
//****************************************************************************80

double *r8po_sl ( int n, double a_lu[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_SL solves a linear system that has been factored by R8PO_FA.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A_LU[N*N], the Cholesky factor from R8PO_FA.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R8PO_SL[N], the solution vector.
//
{
  int i;
  int k;
  double *x;

  x = new double[n];

  for ( k = 0; k < n; k++ )
  {
    x[k] = b[k];
  }
//
//  Solve R' * y = b.
//
  for ( k = 0; k < n; k++ )
  {
    for ( i = 0; i < k; i++ )
    {
      x[k] = x[k] - x[i] * a_lu[i+k*n];
    }
    x[k] = x[k] / a_lu[k+k*n];
  }
//
//  Solve R * x = y.
//
  for ( k = n-1; 0 <= k; k-- )
  {
    x[k] = x[k] / a_lu[k+k*n];
    for ( i = 0; i < k; i++ )
    {
      x[i] = x[i] - a_lu[i+k*n] * x[k];
    }
  }

  return x;
}
//****************************************************************************80

double *r8po_to_r8ge ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_TO_R8GE copies a R8PO matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the R8PO matrix.
//
//    Output, double R8PO_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i <= j )
      {
        b[i+j*n] = a[i+j*n];
      }
      else
      {
        b[i+j*n] = a[j+i*n];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8po_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_ZERO zeros a R8PO matrix.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//    N must be positive.
//
//    Output, double R8PO_ZERO[N*N], the R8PO matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double r8pp_det ( int n, double a_lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_DET computes the determinant of a matrix factored by R8PP_FA.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A_LU[(N*(N+1))/2], the LU factors from R8PO_FA.
//
//    Output, double R8PP_DET, the determinant of A.
//
{
  double det;
  int i;
  int k;

  det = 1.0;

  k = 0;
  for ( i = 0; i < n; i++ )
  {
    k = k + i;
    det = det * a_lu[k];
  }

  return det;
}
//****************************************************************************80

double *r8pp_fa ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_FA factors a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Output, double R8PP_FA[(N*(N+1))/2], an upper triangular matrix R, stored 
//    in packed form, so that A = R'*R.
//
{
  double *b;
  int i;
  int j;
  int jj;
  int k;
  int kj;
  int kk;
  double s;
  double t;

  b = new double[(n*(n+1))/2];

  for ( i = 0; i < (n*(n+1))/2; i++ )
  {
    b[i] = a[i];
  }

  jj = 0;

  for ( j = 1; j <= n; j++ )
  {
    s = 0.0;
    kj = jj;
    kk = 0;

    for ( k = 1; k <= j-1; k++ )
    {
      kj = kj + 1;
      t = b[kj-1];
      for ( i = 1; i <= k-1; i++ )
      {
        t = t - b[kk+i-1] * b[jj+i-1];
      }
      kk = kk + k;
      t = t / b[kk-1];
      b[kj-1] = t;
      s = s + t * t;
    }

    jj = jj + j;
    s = b[jj-1] - s;

    if ( s <= 0.0 )
    {
      return NULL;
    }

    b[jj-1] = sqrt ( s );
  }

  return b;
}
//****************************************************************************80

double *r8pp_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_INDICATOR sets up a R8PP indicator matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 February 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Output, double R8PP_INDICATOR((N*(N+1))/2), the R8PP matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;
  int k;

  a = new double[(n*(n+1))/2];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  k = 0;
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= j; i++ )
    {
      a[k] = ( double ) ( fac * i + j );
      k = k + 1;
    }
  }

  return a;
}
//****************************************************************************80

double *r8pp_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_MXV multiplies a R8PP matrix times a vector.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8PP_MXV[N], the product A * x.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < i; j++ )
    {
      k = j + ( i * ( i + 1 ) ) / 2;
      b[i] = b[i] + a[k] * x[j];
    }
    for ( j = i; j < n; j++ )
    {
      k = i + ( j * ( j + 1 ) ) / 2;
      b[i] = b[i] + a[k] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void r8pp_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_PRINT prints a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Input, string TITLE, a title.
//
{
  r8pp_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8pp_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_PRINT_SOME prints some of a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j )
        {
          aij = a[i-1+(j*(j-1))/2];
        }
        else
        {
          aij = a[j-1+(i*(i-1))/2];
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8pp_random ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_RANDOM randomizes a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//    The matrix is computed by setting a "random" upper triangular
//    Cholesky factor R, and then computing A = R'*R.
//    The randomness is limited by the fact that all the entries of
//    R will be between 0 and 1.  A truly random R is only required
//    to have positive entries on the diagonal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8PP_RANDOM[(N*(N+1))/2], the R8PP matrix.
//
{
  double *a;
  int i;
  int ii;
  int ij;
  int ik;
  int j;
  int k;
  int kj;

  a = new double[(n*(n+1))/2];

  for ( i = 0; i < (n*(n+1))/2; i++ )
  {
    a[i] = 0.0;
  }

  for ( i = n; 1 <=i; i-- )
  {
//
//  Set row I of R.
//
    for ( j = i; j <= n; j++ )
    {
      ij = i + ( j * ( j - 1 ) ) / 2;
      a[ij-1] = r8_uniform_01 ( seed );
    }
//
//  Consider element J of row I, last to first.
//
    for ( j = n; i <= j; j-- )
    {
//
//  Add multiples of row I to lower elements of column J.
//
      ij = i + ( j * ( j - 1 ) ) / 2;

      for ( k = i+1; k <= j; k++ )
      {
        kj = k + (j*(j-1))/2;
        ik = i + (k*(k-1))/2;
        a[kj-1] = a[kj-1] + a[ik-1] * a[ij-1];
      }
//
//  Reset element J.
//
      ii = i + (i*(i-1))/2;
      a[ij-1] = a[ii-1] * a[ij-1];
    }
  }

  return a;
}
//****************************************************************************80

double *r8pp_sl ( int n, double a_lu[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_SL solves a R8PP system factored by R8PP_FA.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A_LU[(N*(N+1))/2], the LU factors from R8PP_FA.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R8PP_SL[N], the solution.
//
{
  int i;
  int k;
  int kk;
  double t;
  double *x;

  x = new double[n];

  kk = 0;

  for ( k = 1; k <= n; k++ )
  {
    t = 0.0;
    for ( i = 0; i < k-1; i++ )
    {
      t = t + a_lu[kk+i] * x[i];
    }
    kk = kk + k;
    x[k-1] = ( b[k-1] - t ) / a_lu[kk-1];
  }

  for ( k = n; 1 <= k; k-- )
  {
    x[k-1] = x[k-1] / a_lu[kk-1];
    kk = kk - k;
    t = -x[k-1];
    for ( i = 0; i < k-1; i++ )
    {
      x[i] = x[i] + t * a_lu[kk+i];
    }
  }

  return x;
}
//****************************************************************************80

double *r8pp_to_r8ge ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_TO_R8GE copies a R8PP matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Output, double R8PP_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i <= j )
      {
        b[i+j*n] = a[i+(j*(j+1))/2];
      }
      else
      {
        b[i+j*n] = a[j+(i*(i+1))/2];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8pp_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_ZERO zeros a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//    N must be positive.
//
//    Output, double R8PP_ZERO[(N*(N+1))/2], the R8PP matrix.
//
{
  double *a;
  int k;

  a = new double[(n*(n+1))/2];

  for ( k = 0; k < (n*(n+1))/2; k++ )
  {
    a[k] = 0.0;
  }

  return a;
}
//****************************************************************************80

double *r8ri_to_r8ge ( int nz, int ija[], double sa[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8RI_TO_R8GE converts an R8RI matrix to R8GE form.
//
//  Discussion:
//
//    An R8RI matrix is in row indexed sparse storage form.
//
//    A R8GE matrix is in general storage.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, int NZ, the size required for the RI
//    or "row indexed" sparse storage.
//
//    Input, int IJA[NZ], the index vector.
//
//    Input, double SA[NZ], the value vector.
//
//    Input, int N, the order of the matrix.
//
//    Output, double R8RI_TO_R8GE[N*N], the matrix stored in GE 
//    or "general" format.
//
{
  double *a;
  int i;
  int j;
  int k;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }

  for ( k = 0; k < n; k++ )
  {
    i = k;
    j = k;
    a[i+j*n] = sa[k];
  }

  for ( i = 0; i < n; i++ )
  {
    for ( k = ija[i]; k < ija[i+1]; k++ )
    {
      j = ija[k];
      a[i+j*n] = sa[k];
    }
  }

  return a;
}
//****************************************************************************80

void r8row_swap ( int m, int n, double a[], int row1, int row2 )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_SWAP swaps two rows of a real array.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based!  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, double A[M*N], an array of data.
//
//    Input, int ROW1, ROW2, the two rows to swap.
//    These indices should be between 1 and M.
//
{
# define OFFSET 1

  int j;
  double t;
//
//  Check.
//
  if ( row1 < 0+OFFSET || m-1+OFFSET < row1 )
  {
    cerr << "\n";
    cerr << "R8ROW_SWAP - Fatal error!\n";
    cerr << "  ROW1 is out of range.\n";
    exit ( 1 );
  }

  if ( row2 < 0+OFFSET || m-1+OFFSET < row2 )
  {
    cerr << "\n";
    cerr << "R8ROW_SWAP - Fatal error!\n";
    cerr << "  ROW2 is out of range.\n";
    exit ( 1 );
  }

  if ( row1 == row2 )
  {
    return;
  }
  for ( j = 0; j < n; j++ )
  {
    t                  = a[row1-OFFSET+j*m];
    a[row1-OFFSET+j*m] = a[row2-OFFSET+j*m];
    a[row2-OFFSET+j*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

double *r8s3_indicator ( int n, int nz_num, int isym, int row[], int col[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8S3_INDICATOR sets up a R8S3 indicator matrix.
//
//  Discussion:
//
//    The R8S3 storage format corresponds to the SLAP Triad format.
//
//    The R8S3 storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.  The entries may be given in any order.  No
//    check is made for the erroneous case in which a given matrix entry is
//    specified more than once.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NZ_NUM, the number of nonzero entries.
//
//    Input, int ISYM, is 0 if the matrix is not symmetric, and 1
//    if the matrix is symmetric.  If the matrix is symmetric, then
//    only the nonzeroes on the diagonal and in the lower triangle are stored.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Output, double R8S3_INDICATOR[NZ_NUM], the indicator matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;
  int k;

  a = new double[nz_num];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( k = 0; k < nz_num; k++ )
  {
    i = row[k];
    j = col[k];
    a[k] = ( double ) ( fac * i + j );
  }

  return a;
}
//****************************************************************************80

void r8s3_print ( int m, int n, int nz_num, int isym, int row[], int col[], 
  double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8S3_PRINT prints a R8S3 matrix.
//
//  Discussion:
//
//    The R8S3 storage format corresponds to the SLAP Triad format.
//
//    The R8S3 storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.  The entries may be given in any order.  No
//    check is made for the erroneous case in which a given matrix entry is
//    specified more than once.
//
//    There is a symmetry option for square matrices.  If the symmetric storage
//    option is used, the format specifies that only nonzeroes on the diagonal
//    and lower triangle are stored.  However, this routine makes no attempt
//    to enforce this.  The only thing it does is to "reflect" any nonzero
//    offdiagonal value.  Moreover, no check is made for the erroneous case
//    in which both A(I,J) and A(J,I) are specified, but with different values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ISYM, is 0 if the matrix is not symmetric, and 1
//    if the matrix is symmetric.  If the matrix is symmetric, then
//    only the nonzeroes on the diagonal and in the lower triangle are stored.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements 
//    of the matrix.
//
//    Input, string TITLE, a title.
//
{
  r8s3_print_some ( m, n, nz_num, isym, row, col, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8s3_print_some ( int m, int n, int nz_num, int isym, int row[], int col[], 
  double a[], int ilo, int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8S3_PRINT_SOME prints some of a R8S3 matrix.
//
//  Discussion:
//
//    The R8S3 storage format corresponds to the SLAP Triad format.
//
//    The R8S3 storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.  The entries may be given in any order.  No
//    check is made for the erroneous case in which a given matrix entry is
//    specified more than once.
//
//    There is a symmetry option for square matrices.  If the symmetric storage
//    option is used, the format specifies that only nonzeroes on the diagonal
//    and lower triangle are stored.  However, this routine makes no attempt
//    to enforce this.  The only thing it does is to "reflect" any nonzero
//    offdiagonal value.  Moreover, no check is made for the erroneous case
//    in which both A(I,J) and A(J,I) are specified, but with different values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ISYM, is 0 if the matrix is not symmetric, and 1
//    if the matrix is symmetric.  If the matrix is symmetric, then
//    only the nonzeroes on the diagonal and in the lower triangle are stored.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements 
//    of the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int index[INCX];
  int j;
  int j2;
  int j2hi;
  int j2lo;
  int k;
  bool nonzero;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";

    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      nonzero = false;
 
      for ( j2 = 0; j2 < inc; j2++ )
      {
        index[j2] = -1;
      }

      for ( k = 0; k < nz_num; k++ )
      {
        if ( i == row[k] && j2lo <= col[k] && col[k] <= j2hi )
        {
          j2 = col[k] - j2lo + 1;

          if ( a[k] != 0.0 )
          {
            index[j2-1] = k;
            nonzero = true;
          }
        }
        else if ( isym == 1 && m == n &&
          i == col[k] && j2lo <= row[k] && row[k] <= j2hi )
        {
          j2 = row[k] - j2lo + 1;

          if ( a[k] != 0.0 )
          {
            index[j2-1] = k;
            nonzero = true;
          }
        }
      }

      if ( nonzero )
      {
        cout << setw(5) << i << " ";
        for ( j2 = 0; j2 < inc; j2++ )
        {
          if ( 0 <= index[j2] )
          {
            aij = a[index[j2]];
          }
          else
          {
            aij = 0.0;
          }
          cout << setw(14) << aij;
        }
        cout << "\n";
      }
    }
  }
  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void r8s3_read ( string input_file, int n, int nz_num, int row[], int col[], 
  double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8S3_READ reads a square R8S3 matrix from a file.
//
//  Discussion:
//
//    This routine needs the value of NZ_NUM, which can be determined
//    by a call to R8S3_READ_SIZE.
//
//    The R8S3 storage format corresponds to the SLAP Triad format.
//
//    The R8S3 storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.  The entries may be given in any order.  No
//    check is made for the erroneous case in which a given matrix entry is
//    specified more than once.
//
//    There is a symmetry option for square matrices.  If the symmetric storage
//    option is used, the format specifies that only nonzeroes on the diagonal
//    and lower triangle are stored.  However, this routine makes no attempt
//    to enforce this.  The only thing it does is to "reflect" any nonzero
//    offdiagonal value.  Moreover, no check is made for the erroneous case
//    in which both A(I,J) and A(J,I) are specified, but with different values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILE, the name of the file to be read.
//
//    Unused, int N, the order of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Output, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Output, double A[NZ_NUM], the nonzero elements of the matrix.
//
{
  ifstream input;
  int k;

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8S3_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_file << "\"\n";
    exit ( 1 );
  }

  for ( k = 0; k < nz_num; k++ )
  {
    input >> row[k] >> col[k] >> a[k];
  }

  input.close ( );

  return;
}
//****************************************************************************80

void r8s3_read_size ( string input_file, int *n, int *nz_num )

//****************************************************************************80
//
//  Purpose:
//
//    R8S3_READ_SIZE reads the size of a square R8S3 matrix from a file.
//
//  Discussion:
//
//    The value of NZ_NUM is simply the number of records in the input file.
//
//    The value of N is determined as the maximum entry in the row and column
//    vectors.
//
//    The R8S3 storage format corresponds to the SLAP Triad format.
//
//    The R8S3 storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.  The entries may be given in any order.  No
//    check is made for the erroneous case in which a given matrix entry is
//    specified more than once.
//
//    There is a symmetry option for square matrices.  If the symmetric storage
//    option is used, the format specifies that only nonzeroes on the diagonal
//    and lower triangle are stored.  However, this routine makes no attempt
//    to enforce this.  The only thing it does is to "reflect" any nonzero
//    offdiagonal value.  Moreover, no check is made for the erroneous case
//    in which both A(I,J) and A(J,I) are specified, but with different values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILE, the name of the file to 
//    be read.
//
//    Output, int *N, the order of the matrix.
//
//    Output, int *NZ_NUM, the number of nonzero elements in the matrix.
//
{
  double a_k;
  int col_k;
  ifstream input;
  int row_k;

  *nz_num = 0;
  *n = 0;

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8S3_READ_SIZE - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_file << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    input >> row_k >> col_k >> a_k;

    if ( input.eof ( ) )
    {
      break;
    }

    *nz_num = *nz_num + 1;
    *n = i4_max ( *n, row_k );
    *n = i4_max ( *n, col_k );
  }

  input.close ( );

  return;
}
//****************************************************************************80

void r8s3_write ( int n, int nz_num, int isym, int row[], int col[], 
  double a[], string output_file )

//****************************************************************************80
//
//  Purpose:
//
//    R8S3_WRITE writes a square R8S3 matrix to a file.
//
//  Discussion:
//
//    The R8S3 storage format corresponds to the SLAP Triad format.
//
//    The R8S3 storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.  The entries may be given in any order.  No
//    check is made for the erroneous case in which a given matrix entry is
//    specified more than once.
//
//    There is a symmetry option for square matrices.  If the symmetric storage
//    option is used, the format specifies that only nonzeroes on the diagonal
//    and lower triangle are stored.  However, this routine makes no attempt
//    to enforce this.  The only thing it does is to "reflect" any nonzero
//    offdiagonal value.  Moreover, no check is made for the erroneous case
//    in which both A(I,J) and A(J,I) are specified, but with different values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ISYM, is 0 if the matrix is not symmetric, and 1
//    if the matrix is symmetric.  If the matrix is symmetric, then
//    only the nonzeroes on the diagonal and in the lower triangle are stored.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements 
//    of the matrix.
//
//    Input, string OUTPUT_FILE, the name of the file to which
//    the information is to be written.
//
{
  int k;
  ofstream output;

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8S3_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }

  for ( k = 0; k < nz_num; k++ )
  {
    output << "  " << setw(8) << row[k]
           << "  " << setw(8) << col[k]
           << "  " << setw(16) << a[k] << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80

double *r8sd_cg ( int n, int ndiag, int offset[], double a[], double b[], 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SD_CG uses the conjugate gradient method on a R8SD linear system.
//
//  Discussion:
//
//    The R8SD storage format is used for symmetric matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0, and 
//    each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//
//    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
//    we then create an array B that has N rows and NDIAG columns, and simply
//    "collapse" the matrix A to the left:
//
//    For the conjugate gradient method to be applicable, the matrix A must 
//    be a positive definite symmetric matrix.
//
//    The method is designed to reach the solution to the linear system
//      A * x = b
//    after N computational steps.  However, roundoff may introduce
//    unacceptably large errors for some problems.  In such a case,
//    calling the routine a second time, using the current solution estimate
//    as the new starting guess, should result in improved results.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Frank Beckman,
//    The Solution of Linear Equations by the Conjugate Gradient Method,
//    in Mathematical Methods for Digital Computers,
//    edited by John Ralston, Herbert Wilf,
//    Wiley, 1967,
//    ISBN: 0471706892,
//    LC: QA76.5.R3.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals that are stored.
//    NDIAG must be at least 1 and no more than N.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8SD matrix.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, double X[N], an estimate for the solution, which may be 0.
//
//    Output, double R8SD_CG[N], the approximate solution vector.  Note that 
//    repeated calls to this routine, using the approximate solution
//    output on the previous call, MAY improve the solution.
//
{
  double alpha;
  double *ap;
  double beta;
  int i;
  int it;
  double *p;
  double pap;
  double pr;
  double *r;
  double rap;
  double *x_new;
//
//  Initialize
//    AP = A * x,
//    R  = b - A * x,
//    P  = b - A * x.
//
  x_new = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x_new[i] = x[i];
  }

  ap = r8sd_mxv ( n, ndiag, offset, a, x_new );

  r = new double[n];
  for ( i = 0; i < n; i++ )
  {
    r[i] = b[i] - ap[i];
  }

  p = new double[n];
  for ( i = 0; i < n; i++ )
  {
    p[i] = b[i] - ap[i];
  }
//
//  Do the N steps of the conjugate gradient method.
//
  for ( it = 1; it < n; it++ )
  {
//
//  Compute the matrix*vector product AP = A*P.
//
    delete [] ap;
    ap = r8sd_mxv ( n, ndiag, offset, a, p );
//
//  Compute the dot products
//    PAP = P*AP,
//    PR  = P*R
//  Set
//    ALPHA = PR / PAP.
//
    pap = r8vec_dot_product ( n, p, ap );

    if ( pap == 0.0 )
    {
      break;
    }

    pr = r8vec_dot_product ( n, p, r );

    alpha = pr / pap;
//
//  Set
//    X = X + ALPHA * P
//    R = R - ALPHA * AP.
//
    for ( i = 0; i < n; i++ )
    {
      x_new[i] = x_new[i] + alpha * p[i];
    }
    for ( i = 0; i < n; i++ )
    {
      r[i] = r[i] - alpha * ap[i];
    }
//
//  Compute the vector dot product
//    RAP = R*AP
//  Set
//    BETA = - RAP / PAP.
//
    rap = r8vec_dot_product ( n, r, ap );

    beta = -rap / pap;
//
//  Update the perturbation vector
//    P = R + BETA * P.
//
    for ( i = 0; i < n; i++ )
    {
      p[i] = r[i] + beta * p[i];
    }
  }

  delete [] ap;
  delete [] p;
  delete [] r;

  return x_new;
}
//****************************************************************************80

double *r8sd_indicator ( int n, int ndiag, int offset[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SD_INDICATOR sets up a R8SD indicator matrix.
//
//  Discussion:
//
//    The R8SD storage format is used for symmetric matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0, and 
//    each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//
//    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
//    we then create an array B that has N rows and NDIAG columns, and simply
//    "collapse" the matrix A to the left.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals that are stored.
//    NDIAG must be at least 1 and no more than N.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Output, double R8SD_INDICATOR[N*NDIAG], the R8SD matrix.
//
{
  double *a;
  int diag;
  int fac;
  int i;
  int j;

  a = new double[n*ndiag];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= n; i++ )
  {
    for ( diag = 1; diag <= ndiag; diag++ )
    {
      j = i + offset[diag-1];
      if ( 1 <= j && j <= n )
      {
        a[i-1+(diag-1)*n] = ( double ) ( fac * i + j );
      }
      else
      {
        a[i-1+(diag-1)*n] = 0.0;
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8sd_mxv ( int n, int ndiag, int offset[], double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SD_MXV multiplies a R8SD matrix times a vector.
//
//  Discussion:
//
//    The R8SD storage format is used for symmetric matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0, and 
//    each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//
//    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
//    we then create an array B that has N rows and NDIAG columns, and simply
//    "collapse" the matrix A to the left.
//
//  Example:
//
//    The "offset" value is printed above each column.
//
//    Original matrix               New Matrix
//
//       0   1   2   3   4   5       0   1   3   5
//
//      11  12   0  14   0  16      11  12  14  16
//      21  22  23   0  25   0      22  23  25  --
//       0  32  33  34   0  36      33  34  36  --
//      41   0  43  44  45   0      44  45  --  --
//       0  52   0  54  55  56      55  56  --  --
//      61   0  63   0  65  66      66  --  --  --
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals that are stored.
//    NDIAG must be at least 1 and no more than N.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8SD matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8SD_MXV[N], the product A * x.
//
{
  double *b;
  int i;
  int j;
  int jdiag;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( jdiag = 0; jdiag < ndiag; jdiag++ )
    {
      j = i + offset[jdiag];
      if ( 0 <= j && j < n )
      {
        b[i] = b[i] + a[i+jdiag*n] * x[j];
        if ( offset[jdiag] != 0 )
        {
          b[j] = b[j] + a[i+jdiag*n] * x[i];
        }
      }
    }
  }

  return b;
}
//****************************************************************************80

void r8sd_print ( int n, int ndiag, int offset[], double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SD_PRINT prints a R8SD matrix.
//
//  Discussion:
//
//    The R8SD storage format is used for symmetric matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0, and 
//    each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//
//    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
//    we then create an array B that has N rows and NDIAG columns, and simply
//    "collapse" the matrix A to the left:
//
//  Example:
//
//    The "offset" value is printed above each column.
//
//    Original matrix               New Matrix
//
//       0   1   2   3   4   5       0   1   3   5
//
//      11  12   0  14   0  16      11  12  14  16
//      21  22  23   0  25   0      22  23  25  --
//       0  32  33  34   0  36      33  34  36  --
//      41   0  43  44  45   0      44  45  --  --
//       0  52   0  54  55  56      55  56  --  --
//      61   0  63   0  65  66      66  --  --  --
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than N.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8SD matrix.
//
//    Input, string TITLE, a title.
//
{
  r8sd_print_some ( n, ndiag, offset, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8sd_print_some ( int n, int ndiag, int offset[], double a[], int ilo, 
  int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SD_PRINT_SOME prints some of a R8SD matrix.
//
//  Discussion:
//
//    The R8SD storage format is used for symmetric matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0, and 
//    each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//
//    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
//    we then create an array B that has N rows and NDIAG columns, and simply
//    "collapse" the matrix A to the left:
//
//  Example:
//
//    The "offset" value is printed above each column.
//
//    Original matrix               New Matrix
//
//       0   1   2   3   4   5       0   1   3   5
//
//      11  12   0  14   0  16      11  12  14  16
//      21  22  23   0  25   0      22  23  25  --
//       0  32  33  34   0  36      33  34  36  --
//      41   0  43  44  45   0      44  45  --  --
//       0  52   0  54  55  56      55  56  --  --
//      61   0  63   0  65  66      66  --  --  --
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals of the matrix
//    that are stored in the array.
//    NDIAG must be at least 1, and no more than N.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8SD matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;
  int jdiag;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
    cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        aij = 0.0;

        for ( jdiag = 0; jdiag < ndiag; jdiag++ )
        {
          if ( j - i == offset[jdiag] )
          {
            aij = a[i-1+jdiag*n];
          }
          else if ( j - i == - offset[jdiag] )
          {
            aij = a[j-1+jdiag*n];
          }
        }
        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8sd_random ( int n, int ndiag, int offset[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8SD_RANDOM randomizes a R8SD matrix.
//
//  Discussion:
//
//    The R8SD storage format is used for symmetric matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0, and 
//    each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//
//    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
//    we then create an array B that has N rows and NDIAG columns, and simply
//    "collapse" the matrix A to the left.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals that are stored.
//    NDIAG must be at least 1 and no more than N.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8SD_RANDOM[N*NDIAG], the R8SD matrix.
//
{
  double *a;
  int i;
  int j;
  int jj;

  a = new double[n*ndiag];

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= ndiag; j++ )
    {
      jj = i + offset[j-1];
      if ( 1 <= jj && jj <= n )
      {
        a[i-1+(j-1)*n] = r8_uniform_01 ( seed );
      }
      else
      {
        a[i-1+(j-1)*n] = 0.0;
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8sd_to_r8ge ( int n, int ndiag, int offset[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SD_TO_R8GE copies a R8SD matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8SD storage format is used for symmetric matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0, and 
//    each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//
//    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
//    we then create an array B that has N rows and NDIAG columns, and simply
//    "collapse" the matrix A to the left:
//
//  Example:
//
//    The "offset" value is printed above each column.
//
//    Original matrix               New Matrix
//
//       0   1   2   3   4   5       0   1   3   5
//
//      11  12   0  14   0  16      11  12  14  16
//      21  22  23   0  25   0      22  23  25  --
//       0  32  33  34   0  36      33  34  36  --
//      41   0  43  44  45   0      44  45  --  --
//       0  52   0  54  55  56      55  56  --  --
//      61   0  63   0  65  66      66  --  --  --
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals that are stored.
//    NDIAG must be at least 1 and no more than N.
//
//    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
//
//    Input, double A[N*NDIAG], the R8SD matrix.
//
//    Output, double R8SD_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;
  int jj;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < ndiag; j++ )
    {
      jj = i + offset[j];
      if ( 0 <= jj && jj <= n-1 )
      {
        b[i+jj*n] = a[i+j*n];
        if ( i != jj )
        {
          b[jj+i*n] = a[i+j*n];
        }
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8sd_zero ( int n, int ndiag )

//****************************************************************************80
//
//  Purpose:
//
//    R8SD_ZERO zeros a R8SD matrix.
//
//  Discussion:
//
//    The R8SD storage format is used for symmetric matrices whose only nonzero entries
//    occur along a few diagonals, but for which these diagonals are not all
//    close enough to the main diagonal for band storage to be efficient.
//
//    In that case, we assign the main diagonal the offset value 0, and 
//    each successive superdiagonal gets an offset value 1 higher, until
//    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
//
//    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
//    we then create an array B that has N rows and NDIAG columns, and simply
//    "collapse" the matrix A to the left.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NDIAG, the number of diagonals that are stored.
//    NDIAG must be at least 1 and no more than N.
//
//    Output, double R8SD_ZERO[N*NDIAG], the R8SD matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*ndiag];

  for ( j = 0; j < ndiag; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

double *r8sm_ml ( int n, double a_lu[], double u[], double v[], int pivot[], 
  double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_ML multiplies a factored square R8SM matrix times a vector.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA.
//
//    Input, double U[N], V[N], the Sherman Morrison vectors.
//
//    Input, int PIVOT[N], the pivot vector computed by R8GE_FA.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Input, int JOB, specifies the operation to be done:
//    JOB = 0, compute (A-u*v') * x.
//    JOB nonzero, compute (A-u*v')' * x.
//
//    Output, double R8SM_ML[N], the result of the multiplication.
//
{
  double *b;
  int i;
  double ux;
  double vx;

  b = r8ge_ml ( n, a_lu, pivot, x, job );

  if ( job == 0 )
  {
    vx = 0.0;
    for ( i = 0; i < n; i++ )
    {
      vx = vx + v[i] * x[i];
    }
    for ( i = 0; i < n; i++ )
    {
      b[i] = b[i] - u[i] * vx;
    }
  }
  else
  {
    ux = 0.0;
    for ( i = 0; i < n; i++ )
    {
      ux = ux + u[i] * x[i];
    }
    for ( i = 0; i < n; i++ )
    {
      b[i] = b[i] - v[i] * ux;
    }
  }

  return b;
}
//****************************************************************************80

double *r8sm_mxv ( int m, int n, double a[], double u[], double v[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_MXV multiplies a R8SM matrix times a vector.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M*N], the R8SM matrix.
//
//    Input, double U[M], V[N], the R8SM vectors U and V.
//
//    Input, double X[N], the vector to be multiplied by (A-u*v').
//
//    Output, double R8SM_MXV[M], the product (A-u*v') * x.
//
{
  double *b;
  int i;
  int j;
  double vx;

  b = new double[m];

  vx = 0.0;
  for ( j = 0; j < n; j++ )
  {
    vx = vx + v[j] * x[j];
  }

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
    b[i] = b[i] - u[i] * vx;
  }

  return b;
}
//****************************************************************************80

void r8sm_print ( int m, int n, double a[], double u[], double v[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_PRINT prints a R8SM matrix.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M*N], the R8SM matrix.
//
//    Input, double U[M], V[N], the R8SM vectors.
//
//    Input, string TITLE, a title.
//
{
  r8sm_print_some ( m, n, a, u, v, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8sm_print_some ( int m, int n, double a[], double u[], double v[], int ilo, 
  int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_PRINT_SOME prints some of a R8SM matrix.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M*N], the R8SM matrix.
//
//    Input, double U[M], V[N], the R8SM vectors.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*n] - u[i-1] * v[j-1] << "  ";
      }
      cout << "\n";
    }
  }
  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void r8sm_random ( int m, int n, double a[], double u[], double v[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_RANDOM randomizes a R8SM matrix.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Output, double A[M*N], the R8SM matrix.
//
//    Output, double U[M], V[N], the R8SM vectors.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = r8_uniform_01 ( seed );
    }
  }
  for ( i = 0; i < m; i++ )
  {
    u[i] = r8_uniform_01 ( seed );
  }
  for ( j = 0; j < n; j++ )
  {
    v[j] = r8_uniform_01 ( seed );
  }

  return;
}
//****************************************************************************80

double *r8sm_sl ( int n, double a_lu[], double u[], double v[], double b[], 
  int pivot[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_SL solves a square R8SM system that has been factored.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//    It is assumed that A has been decomposed into its LU factors
//    by R8GE_FA.  The Sherman Morrison formula allows
//    us to solve linear systems involving (A-u*v') by solving linear
//    systems involving A and adjusting the results.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Stephen Nash
//    Numerical Methods and Software,
//    Prentice Hall, 1989
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA.
//
//    Input, double U[N], V[N], the R8SM vectors U and V.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, int PIVOT[N], the pivot vector produced by R8GE_FA.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve (A-u*v') * X = B.
//    nonzero, solve (A-u*v') * X = B.
//
//    Output, double R8SM_SL[N], the solution vector, or NULL if
//    an error occurred.
//
{
  double alpha;
  double beta;
  int i;
  int job_local;
  double *w;
  double *x;

  x = new double[n];

  if ( job == 0 )
  {
//
//  Solve A' * w = v.
//
    job_local = 1;
    w = r8ge_sl_new ( n, a_lu, pivot, v, job_local );
//
//  Set beta = w' * b.
//
    beta = 0.0;
    for ( i = 0; i < n; i++ )
    {
      beta = beta + w[i] * b[i];
    }
//
//  Solve A * x = b.
//
    job_local = 0;
    x = r8ge_sl_new ( n, a_lu, pivot, b, job_local );
//
//  Solve A * w = u.
//
    job_local = 0;
    delete [] w;
    w = r8ge_sl_new ( n, a_lu, pivot, u, job_local );
//
//  Set alpha = 1 / ( 1 - v' * w ).
//
    alpha = 1.0;
    for ( i = 0; i < n; i++ )
    {
      alpha = alpha - v[i] * w[i];
    }
  }
  else
  {
//
//  Solve A * w = u.
//
    job_local = 0;
    w = r8ge_sl_new ( n, a_lu, pivot, u, job_local );
//
//  Set beta = w' * b.
//
    beta = 0.0;
    for ( i = 0; i < n; i++ )
    {
      beta = beta + w[i] * b[i];
    }
//
//  Solve A' * x = b.
//
    job_local = 1;
    x = r8ge_sl_new ( n, a_lu, pivot, b, job_local );
//
//  Solve A' * w = v.
//
    job_local = 1;
    delete [] w;
    w = r8ge_sl_new ( n, a_lu, pivot, v, job_local );
//
//  Set alpha = 1 / ( 1 - u' * w ).
//
    alpha = 1.0;
    for ( i = 0; i < n; i++ ) 
    {
      alpha = alpha - u[i] * w[i];
    }
  }

  if ( alpha == 0.0 )
  {
    cerr << "\n";
    cerr << "R8SM_SL - Fatal error!\n";
    cerr << "  The divisor ALPHA is zero.\n";
    exit ( 1 );
  }

  alpha = 1.0 / alpha;
//
//  Set b = b + alpha * beta * w.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = x[i] + alpha * beta * w[i];
  }
  delete [] w;
  return x;
}
//****************************************************************************80

double *r8sm_to_r8ge ( int m, int n, double a[], double u[], double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_TO_R8GE copies a R8SM matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M*N], the R8SM matrix.
//
//    Input, double U[M], V[N], the R8SM vectors.
//
//    Output, double R8SM_TO_R8GE[M*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i+j*m] = a[i+j*m] - u[i] * v[j];
    }
  }

  return b;
}
//****************************************************************************80

double *r8sm_vxm ( int m, int n, double a[], double u[], double v[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_VXM multiplies a vector by a R8SM matrix.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M*N], the R8SM matrix.
//
//    Input, double U[M], V[N], the R8SM vectors.
//
//    Input, double X[M], the vector to be multiplied.
//
//    Output, double R8SM_VXM[N], the product (A-u*v')' * X.
//
{
  double *b;
  double dot;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < m; j++ )
    {
      b[i] = b[i] + x[j] * a[j+i*m];
    }
    dot = 0.0;
    for ( j = 0; j < m; j++ )
    {
      dot = dot + u[j] * x[j];
    }
    b[i] = b[i] - v[i] * dot;
  }

  return b;
}
//****************************************************************************80

void r8sm_zero ( int m, int n, double a[], double u[], double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SM_ZERO zeros a R8SM matrix.
//
//  Discussion:
//
//    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
//    which is defined by an M by N matrix A, an M vector U, and
//    an N vector V, by B = A - U * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Output, double A[M*N], the R8SM matrix.
//
//    Output, double U[M], V[N], the R8SM vectors.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  for ( i = 0; i < m; i++ )
  {
    u[i] = 0.0;
  }
  for ( i = 0; i < n; i++ )
  {
    v[i] = 0.0;
  }

  return;
}
//****************************************************************************80

bool r8sp_check ( int m, int n, int nz_num, int row[], int col[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_CHECK checks that a R8SP matrix data structure is properly sorted.
//
//  Discussion:
//
//    This routine assumes that the data structure has been sorted,
//    so that the entries of ROW are ascending sorted, and that the
//    entries of COL are ascending sorted, within the group of entries
//    that have a common value of ROW.
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of
//    the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in
//    the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and
//    column indices of the nonzero elements.
//
//    Output, bool R8SP_CHECK, is TRUE if the matrix is properly defined.
//
{
  bool check;
  int k;

  check = true;
//
//  Check 1 <= ROW(*) <= M.
//
  for ( k = 0; k < nz_num; k++ )
  {
    if ( row[k] < 1 || m < row[k] )
    {
      check = false;
      return check;
    }
  }
//
//  Check 1 <= COL(*) <= N.
//
  for ( k = 0; k < nz_num; k++ )
  {
    if ( col[k] < 1 || n < col[k] )
    {
      check = false;
      return check;
    }
  }
//
//  Check that ROW(K) <= ROW(K+1).
//
  for ( k = 0; k < nz_num - 1; k++ )
  {
    if ( row[k+1] < row[k] )
    {
      check = false;
      return check;
    }
  }
//
//  Check that, if ROW(K) == ROW(K+1), that COL(K) < COL(K+1).
//
  for ( k = 0; k < nz_num - 1; k++ )
  {
    if ( row[k] == row[k+1] )
    {
      if ( col[k+1] <= col[k] )
      {
        check = false;
        return check;
      }
    }
  }
  return check;
}
//****************************************************************************80

int r8sp_ij_to_k ( int nz_num, int row[], int col[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_IJ_TO_K seeks the compressed index of the (I,J) entry of A.
//
//  Discussion:
//
//    If A(I,J) is nonzero, then its value is stored in location K.
//
//    This routine searches the R8SP storage structure for the index K
//    corresponding to (I,J), returning -1 if no such entry was found.
//
//    This routine assumes that the data structure has been sorted,
//    so that the entries of ROW are ascending sorted, and that the
//    entries of COL are ascending sorted, within the group of entries
//    that have a common value of ROW.
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NZ_NUM, the number of nonzero elements in
//    the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and
//    column indices of the nonzero elements.
//
//    Input, int I, J, the row and column indices of the
//    matrix entry.
//
//    Output, int R8SP_IJ_TO_K, the R8SP index of the (I,J) entry.
//
{
  int hi;
  int k;
  int lo;
  int md;

  lo = 1;
  hi = nz_num;

  for ( ; ; )
  {
    if ( hi < lo )
    {
      k = -1;
      break;
    }

    md = ( lo + hi ) / 2;

    if ( row[md-1] < i || ( row[md-1] == i && col[md-1] < j ) )
    {
      lo = md + 1;
    }
    else if ( i < row[md-1] || ( row[md-1] == i && j < col[md-1] ) )
    {
      hi = md - 1;
    }
    else
    {
      k = md;
      break;
    }
  }

  return k;
}
//****************************************************************************80

double *r8sp_indicator ( int m, int n, int nz_num, int row[], int col[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_INDICATOR sets up a R8SP indicator matrix.
//
//  Discussion:
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Output, double R8SP_INDICATOR[NZ_NUM], the nonzero elements of the matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;
  int k;

  a = new double[nz_num];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( k = 0; k < nz_num; k++ )
  {
    i = row[k];
    j = col[k];
    a[k] = ( double ) ( fac * i + j );
  }

  return a;
}
//****************************************************************************80

double *r8sp_mxv ( int m, int n, int nz_num, int row[], int col[], 
  double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_MXV multiplies a R8SP matrix times a vector.
//
//  Discussion:
//
//    The R8SP storage format stores the row, column and value of each nonzero
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8SP_MXV[M], the product vector A*X.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( k = 0; k < nz_num; k++ )
  {
    i = row[k];
    j = col[k];
    b[i-1] = b[i-1] + a[k] * x[j-1];
  }

  return b;
}
//****************************************************************************80

void r8sp_print ( int m, int n, int nz_num, int row[], int col[], 
  double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_PRINT prints a R8SP matrix.
//
//  Discussion:
//
//    This version of R8SP_PRINT has been specifically modified to allow,
//    and correctly handle, the case in which a single matrix location
//    A(I,J) is referenced more than once by the sparse matrix structure.
//    In such cases, the routine prints out the sum of all the values.
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
//    Input, string TITLE, a title.
//
{
  r8sp_print_some ( m, n, nz_num, row, col, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8sp_print_some ( int m, int n, int nz_num, int row[], int col[], 
  double a[], int ilo, int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_PRINT_SOME prints some of a R8SP matrix.
//
//  Discussion:
//
//    This version of R8SP_PRINT_SOME has been specifically modified to allow,
//    and correctly handle, the case in which a single matrix location
//    A(I,J) is referenced more than once by the sparse matrix structure.
//    In such cases, the routine prints out the sum of all the values.
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij[INCX];
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;
  int k;
  bool nonzero;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";

    cout << "  Col:  ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      nonzero = false;
      for ( j2 = 0; j2 < INCX; j2++ )
      {
        aij[j2] = 0.0;
      }

      for ( k = 1; k <= nz_num; k++ )
      {
        if ( i == row[k-1] && j2lo <= col[k-1] && col[k-1] <= j2hi )
        {
          j2 = col[k-1] - j2lo;

          if ( a[k-1] == 0.0 )
          {
            continue;
          }

          nonzero = true;
          aij[j2] = aij[j2] + a[k-1];
        }
      }

      if ( nonzero )
      {
        cout << setw(6) << i;
        for ( j2 = 0; j2 < inc; j2++ )
        {
          cout << setw(12) << aij[j2] << "  ";
        }
        cout << "\n";
      }
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8sp_random ( int m, int n, int nz_num, int row[], int col[], 
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_RANDOM sets a random R8SP matrix.
//
//  Discussion:
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8SP_RANDOM[NZ_NUM], the nonzero elements of the matrix.
//
{
  return ( r8vec_random ( nz_num, 0.0, 1.0, seed ) );
}
//****************************************************************************80

void r8sp_read ( string input_file, int m, int n, int nz_num, int row[], 
  int col[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_READ reads a R8SP matrix from a file.
//
//  Discussion:
//
//    This routine needs the value of NZ_NUM, which can be determined
//    by a call to R8SP_READ_SIZE.
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILE, the name of the file to be read.
//
//    Unused, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Output, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Output, double A[NZ_NUM], the nonzero elements of the matrix.
//
{
  ifstream input;
  int k;

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8SP_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_file << "\"\n";
    exit ( 1 );
  }

  for ( k = 0; k < nz_num; k++ )
  {
    input >> row[k] >> col[k] >> a[k];
  }

  input.close ( );

  return;
}
//****************************************************************************80

void r8sp_read_size ( string input_file, int *m, int *n, int *nz_num )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_READ_SIZE reads the size of a R8SP matrix from a file.
//
//  Discussion:
//
//    The value of NZ_NUM is simply the number of records in the input file.
//
//    The values of M and N are determined as the maximum entry in the row 
//    and column vectors.
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILE, the name of the file to 
//    be read.
//
//    Output, int *M, *N, the number of rows and columns of the matrix.
//
//    Output, int *NZ_NUM, the number of nonzero elements in the matrix.
//
{
  double a_k;
  int col_k;
  ifstream input;
  int row_k;

  *m = 0;
  *n = 0;
  *nz_num = 0;

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8SP_READ_SIZE - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_file << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    input >> row_k >> col_k >> a_k;

    if ( input.eof ( ) )
    {
      break;
    }

    *m = i4_max ( *m, row_k );
    *n = i4_max ( *n, col_k );
    *nz_num = *nz_num + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

double *r8sp_to_r8ge ( int m, int n, int nz_num, int row[], int col[], 
  double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_TO_R8GE converts a R8SP matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
//    Output, double B[M*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+j*m] = 0.0;
    }
  }

  for ( k = 0; k < nz_num; k++ )
  {
    b[row[k]-1+(col[k]-1)*m] = a[k];
  }

  return b;
}
//****************************************************************************80

void r8sp_to_r8ncf ( int m, int n, int nz_num, int row[], int col[], 
  double a[], int rowcol[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_TO_R8NCF converts a R8SP matrix to a R8NCF matrix.
//
//  Discussion:
//
//    The R8SP and R8NCF formats are essentially identical, except that
//    R8SP keeps separate ROW and COLUMN vectors, while R8NCF uses a single
//    ROWCOL array.  Therefore, the input values NZ_NUM and A used in
//    the R8SP representation can be regarded as part of the output
//    values used for the R8NCF representation.
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//    The R8NCF storage format stores NZ_NUM, the number of nonzeros,
//    a real array containing the nonzero values, a 2 by NZ_NUM integer
//    array storing the row and column of each nonzero entry.
//
//    The R8NCF format is used by NSPCG.  NSPCG requires that the information
//    for the diagonal entries of the matrix must come first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
//    Output, int ROWCOL[2*NZ_NUM], the R8NCF row and column index vector.
//
{
  int i;

  for ( i = 0; i < nz_num; i++ )
  {
    rowcol[0+i*2] = row[i];
    rowcol[1+i*2] = col[i];
  }

  return;
}
//****************************************************************************80

double *r8sp_vxm ( int m, int n, int nz_num, int row[], int col[], 
  double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_VXM multiplies a vector times a R8SP matrix.
//
//  Discussion:
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements of the matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double B[N], the product vector A'*X.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n];

  for ( j = 0; j < n; j++ )
  {
    b[j] = 0.0;
  }

  for ( k = 0; k < nz_num; k++ )
  {
    i = row[k];
    if ( i < 0 || m < i )
    {
      cerr << "\n";
      cerr << "R8SP_VXM - Fatal error!\n";
      cerr << "  I < 0 or " << m << " M < I.\n";
      exit ( 1 );
    }
    j = col[k];
    if ( j < 0 || n < j )
    {
      cerr << "\n";
      cerr << "R8SP_VXM - Fatal error!\n";
      cerr << "  J < 0 or " << n << " N < J.\n";
      exit ( 1 );
    }
    b[j-1] = b[j-1] + a[k] * x[i-1];
  }

  return b;
}
//****************************************************************************80

void r8sp_write ( int m, int n, int nz_num, int row[], int col[], double a[], 
  string output_file )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_WRITE writes a R8SP matrix to a file.
//
//  Discussion:
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    It is possible that a pair of indices (I,J) may occur more than
//    once.  Presumably, in this case, the intent is that the actual value
//    of A(I,J) is the sum of all such entries.  This is not a good thing
//    to do, but I seem to have come across this in MATLAB.
//
//    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
//    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Input, double A[NZ_NUM], the nonzero elements 
//    of the matrix.
//
//    Input, string OUTPUT_FILE, the name of the file to which
//    the information is to be written.
//
{
  int k;
  ofstream output;

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8SP_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }

  for ( k = 0; k < nz_num; k++ )
  {
    output << "  " << setw(8) << row[k]
           << "  " << setw(8) << col[k]
           << "  " << setw(16) << a[k] << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80

double *r8sp_zero ( int m, int n, int nz_num, int row[], int col[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SP_ZERO returns a zero R8SP matrix.
//
//  Discussion:
//
//    The R8SP storage format stores the row, column and value of each nonzero
//    entry of a sparse matrix.
//
//    The R8SP format is used by DLAP/SLAP (nonsymmetric case), by MATLAB,
//    and by SPARSEKIT ("COO" format).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, int NZ_NUM, the number of nonzero elements in the matrix.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
//    of the nonzero elements.
//
//    Output, double R8SP_ZERO[NZ_NUM], the (potentially) nonzero elements 
//    of the matrix.
//
{
  double *a;
  int i;

  a = new double[nz_num];

  for ( i = 0; i < nz_num; i++ )
  {
    a[i] = 0.0;
  }

  return a;
}
//****************************************************************************80

void r8sr_indicator ( int n, int nz, int row[], int col[], double diag[], 
  double off[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SR_INDICATOR sets up a R8SR indicator matrix.
//
//  Discussion:
//
//    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
//    The off-diagonal entries of row I are stored in entries ROW(I)
//    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
//    the entry stored in OFF(J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ, the number of offdiagonal nonzero elements in A.
//
//    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
//    are contained in A(ROW(I)) through A(ROW(I+1)-1).
//
//    Input, int COL[NZ], contains the column index of the element
//    in the corresponding position in A.
//
//    Output, double DIAG[N], the diagonal elements of A.
//
//    Output, double OFF[NZ], the off-diagonal elements of A.
//
{
  int fac;
  int i;
  int j;
  int k;

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= n; i++ )
  {
    j = i;
    diag[i-1] = ( double ) ( fac * i + j );

    for ( k = row[i-1]; k <= row[i]-1; k++ )
    {
      j = col[k-1];
      off[k-1] = ( double ) ( fac * i + j );
    }
  }

  return;
}
//****************************************************************************80

double *r8sr_mxv ( int n, int nz, int row[], int col[], double diag[], 
  double off[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SR_MXV multiplies a R8SR matrix times a vector.
//
//  Discussion:
//
//    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
//    The off-diagonal entries of row I are stored in entries ROW(I)
//    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
//    the entry stored in OFF(J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ, the number of offdiagonal nonzero elements in A.
//
//    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
//    are contained in A(ROW(I)) through A(ROW(I+1)-1).
//
//    Input, int COL[NZ], contains the column index of the element
//    in the corresponding position in A.
//
//    Input, double DIAG[N], the diagonal elements of A.
//
//    Input, double OFF[NZ], the off-diagonal elements of A.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8SR_MXV[N], the product A * X.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = diag[i] * x[i];
  }

  for ( i = 0; i < n; i++ )
  {
    for ( k = row[i]; k <= row[i+1] - 1; k++ )
    {
      j = col[k-1];
      b[i] = b[i] + off[k-1] * x[j-1];
    }
  }

  return b;
}
//****************************************************************************80

void r8sr_print ( int n, int nz, int row[], int col[], double diag[], 
  double off[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SR_PRINT prints a R8SR matrix.
//
//  Discussion:
//
//    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
//    The off-diagonal entries of row I are stored in entries ROW(I)
//    through ROW(I+1)-1 of OFF.  COL(J) records the column index
//    of the entry in A(J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ, the number of offdiagonal nonzero elements in A.
//
//    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
//    are contained in A(ROW(I)) through A(ROW(I+1)-1).
//
//    Input, int COL[NZ], contains the column index of the element
//    in the corresponding position in A.
//
//    Input, double DIAG[N], the diagonal elements of A.
//
//    Input, double OFF[NZ], the off-diagonal elements of A.
//
//    Input, string TITLE, a title.
//
{
  r8sr_print_some ( n, nz, row, col, diag, off, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8sr_print_some ( int n, int nz, int row[], int col[], double diag[], 
  double off[], int ilo, int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SR_PRINT_SOME prints some of a R8SR matrix.
//
//  Discussion:
//
//    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
//    The off-diagonal entries of row I are stored in entries ROW(I)
//    through ROW(I+1)-1 of OFF.  COL(J) records the column index
//    of the entry in A(J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ, the number of offdiagonal nonzero elements in A.
//
//    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
//    are contained in A(ROW(I)) through A(ROW(I+1)-1).
//
//    Input, int COL[NZ], contains the column index of the element
//    in the corresponding position in A.
//
//    Input, double DIAG[N], the diagonal elements of A.
//
//    Input, double OFF[NZ], the off-diagonal elements of A.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;
  int k;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col:  ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";

    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        aij = 0.0;
        if ( j == i )
        {
          aij = diag[i-1];
        }
        else
        {
          for ( k = row[i-1]; k <= row[i]-1; k++ )
          {
            if ( j == col[k-1] )
            {
              aij = off[k-1];
            }
          }
        }
        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void r8sr_random ( int n, int nz, int row[], int col[], double diag[], 
  double off[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8SR_RANDOM randomizes a R8SR matrix.
//
//  Discussion:
//
//    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
//    The off-diagonal entries of row I are stored in entries ROW(I)
//    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
//    the entry stored in OFF(J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ, the number of offdiagonal nonzero elements in A.
//
//    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
//    are contained in A(ROW(I)) through A(ROW(I+1)-1).
//
//    Input, int COL[NZ], contains the column index of the element
//    in the corresponding position in A.
//
//    Output, double DIAG[N], the diagonal elements of A.
//
//    Output, double OFF[NZ], the off-diagonal elements of A.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
{
  int i;
  int j;

  for ( i = 0; i < n; i++ )
  {
    diag[i] = r8_uniform_01 ( seed );
    for ( j = row[i]-1; j <= row[i+1]-2; j++ )
    {
      off[j] = r8_uniform_01 ( seed );
    }
  }

  return;
}
//****************************************************************************80

double *r8sr_to_r8ge ( int n, int nz, int row[], int col[], double diag[], 
  double off[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SR_TO_R8GE converts a R8SR matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
//    The off-diagonal entries of row I are stored in entries ROW(I)
//    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
//    the entry stored in OFF(J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ, the number of offdiagonal nonzero elements in A.
//
//    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
//    are contained in A(ROW(I)) through A(ROW(I+1)-1).
//
//    Input, int COL[NZ], contains the column index of the element
//    in the corresponding position in A.
//
//    Input, double DIAG[N], the diagonal elements of A.
//
//    Input, double OFF[NZ], the off-diagonal elements of A.
//
//    Output, double R8SR_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    b[i+i*n] = diag[i];
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = row[i]-1; j <= row[i+1]-2; j++ )
    {
      b[i+(col[j]-1)*n] = off[j];
    }
  }

  return b;
}
//****************************************************************************80

double *r8sr_vxm ( int n, int nz, int row[], int col[], double diag[], 
  double off[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SR_VXM multiplies a vector times a R8SR matrix.
//
//  Discussion:
//
//    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
//    The off-diagonal entries of row I are stored in entries ROW(I)
//    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
//    the entry stored in OFF(J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ, the number of offdiagonal nonzero elements in A.
//
//    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
//    are contained in A(ROW(I)) through A(ROW(I+1)-1).
//
//    Input, int COL[NZ], contains the column index of the element
//    in the corresponding position in A.
//
//    Input, double DIAG[N], the diagonal elements of A.
//
//    Input, double OFF[NZ], the off-diagonal elements of A.
//
//    Input, double X[N], the vector to be multiplies by A.
//
//    Output, double R8SR_VXM[N], the product A' * X.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = diag[i] * x[i];
  }

  for ( i = 0; i < n; i++ )
  {
    for ( k = row[i]; k <= row[i+1] - 1; k++ )
    {
      j = col[k-1];
      b[j-1] = b[j-1] + off[k-1] * x[i];
    }
  }

  return b;
}
//****************************************************************************80

void r8sr_zero ( int n, int nz, int row[], int col[], double diag[], 
  double off[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SR_ZERO zeros a R8SR matrix.
//
//  Discussion:
//
//    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
//    The off-diagonal entries of row I are stored in entries ROW(I)
//    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
//    the entry stored in OFF(J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ, the number of offdiagonal nonzero elements in A.
//
//    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
//    are contained in A(ROW(I)) through A(ROW(I+1)-1).
//
//    Input, int COL[NZ], contains the column index of the element
//    in the corresponding position in A.
//
//    Output, double DIAG[N], the diagonal elements of A.
//
//    Output, double OFF[NZ], the off-diagonal elements of A.
//
{
  int i;
  int j;

  for ( i = 0; i < n; i++ )
  {
    diag[i] = 0.0;
    for ( j = row[i]-1; j <= row[i+1]-2; j++ )
    {
      off[j] = 0.0;
    }
  }

  return;
}
//****************************************************************************80

bool r8ss_error ( int diag[], int n, int na )

//****************************************************************************80
//
//  Purpose:
//
//    R8SS_ERROR checks dimensions for a R8SS matrix.
//
//  Discussion:
//
//    The R8SS storage format is used for real symmetric skyline matrices.
//    This storage is appropriate when the nonzero entries of the
//    matrix are generally close to the diagonal, but the number
//    of nonzeroes above each diagonal varies in an irregular fashion.
//
//    In this case, the strategy is essentially to assign column J
//    its own bandwidth, and store the strips of nonzeros one after
//    another.   Note that what's important is the location of the
//    furthest nonzero from the diagonal.  A slot will be set up for
//    every entry between that and the diagonal, whether or not
//    those entries are zero.
//
//    A skyline matrix can be Gauss-eliminated without disrupting
//    the storage format, as long as no pivoting is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIAG[N], the indices in A of the N diagonal elements.
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NA, the dimension of the array A.
//    NA must be at least N.
//
//    Output, bool R8SS_ERROR, is TRUE if an error was detected.
//
{
  int i;

  if ( n < 1 )
  {
    cout << "\n";
    cout << "R8SS_ERROR - Illegal N = " << n << "\n";
    return true;
  }

  if ( na < n )
  {
    cout << "\n";
    cout << "R8SS_ERROR - Illegal NA < N = " << n << "\n";
    return true;
  }

  if ( diag[0] != 1 )
  {
    cout << "\n";
    cout << "R8SS_ERROR - DIAG[0] != 1.\n";
    return true;
  }

  for ( i = 0; i < n-1; i++ )
  { 
    if ( diag[i+1] <= diag[i] ) 
    {
      cout << "\n";
      cout << "R8SS_ERROR - DIAG[I+1] <= DIAG[I].\n";
      return true;
    }
  }

  if ( na < diag[n-1] ) 
  {
    cout << "\n";
    cout << "R8SS_ERROR - NA < DIAG[N-1].\n";
    return true;
  }

  return false;
}
//****************************************************************************80

double *r8ss_indicator ( int n, int *na, int diag[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SS_INDICATOR sets up a R8SS indicator matrix.
//
//  Discussion:
//
//    The R8SS storage format is used for real symmetric skyline matrices.
//    This storage is appropriate when the nonzero entries of the
//    matrix are generally close to the diagonal, but the number
//    of nonzeroes above each diagonal varies in an irregular fashion.
//
//    In this case, the strategy is essentially to assign column J
//    its own bandwidth, and store the strips of nonzeros one after
//    another.   Note that what's important is the location of the
//    furthest nonzero from the diagonal.  A slot will be set up for
//    every entry between that and the diagonal, whether or not
//    those entries are zero.
//
//    A skyline matrix can be Gauss-eliminated without disrupting
//    the storage format, as long as no pivoting is required.
//
//    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
//    although the actual storage needed will generally be about half of
//    that.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Output, int *NA, the dimension of the array A, which for this
//    special case will be the maximum, ( N * ( N + 1 ) ) / 2.
//
//    Output, int DIAG[N], the indices in A of the N diagonal elements.
//
//    Output, double R8SS_INDICATOR[(N*(N+1))/2], the R8SS matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[(n*(n+1))/2];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  *na = 0;

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= j; i++ )
    {
      a[*na] = ( double ) ( fac * i + j );
      *na = *na + 1;
    }
    diag[j-1] = *na;
  }

  return a;
}
//****************************************************************************80

double *r8ss_mxv ( int n, int na, int diag[], double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SS_MXV multiplies a R8SS matrix times a vector.
//
//  Discussion:
//
//    The R8SS storage format is used for real symmetric skyline matrices.
//    This storage is appropriate when the nonzero entries of the
//    matrix are generally close to the diagonal, but the number
//    of nonzeroes above each diagonal varies in an irregular fashion.
//
//    In this case, the strategy is essentially to assign column J
//    its own bandwidth, and store the strips of nonzeros one after
//    another.   Note that what's important is the location of the
//    furthest nonzero from the diagonal.  A slot will be set up for
//    every entry between that and the diagonal, whether or not
//    those entries are zero.
//
//    A skyline matrix can be Gauss-eliminated without disrupting
//    the storage format, as long as no pivoting is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NA, the dimension of the array A.
//    NA must be at least N.
//
//    Input, int DIAG[N], the indices in A of the N diagonal elements.
//
//    Input, double A[NA], the R8SS matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8SS_MXV[N], the product vector A*x.
//
{
  double *b;
  int diagold;
  int i;
  int ilo;
  int j;
  int k;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  diagold = 0;
  k = 0;

  for ( j = 1; j <= n; j++ )
  {
    ilo = j + 1 + diagold - diag[j-1];

    for ( i = ilo; i < j; i++ )
    {
      b[i-1] = b[i-1] + a[k] * x[j-1];
      b[j-1] = b[j-1] + a[k] * x[i-1];
      k = k + 1;
    }

    b[j-1] = b[j-1] + a[k] * x[j-1];
    k = k + 1;
    diagold = diag[j-1];
  }

  return b;
}
//****************************************************************************80

void r8ss_print ( int n, int na, int diag[], double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SS_PRINT prints a R8SS matrix.
//
//  Discussion:
//
//    The R8SS storage format is used for real symmetric skyline matrices.
//    This storage is appropriate when the nonzero entries of the
//    matrix are generally close to the diagonal, but the number
//    of nonzeroes above each diagonal varies in an irregular fashion.
//
//    In this case, the strategy is essentially to assign column J
//    its own bandwidth, and store the strips of nonzeros one after
//    another.   Note that what's important is the location of the
//    furthest nonzero from the diagonal.  A slot will be set up for
//    every entry between that and the diagonal, whether or not
//    those entries are zero.
//
//    A skyline matrix can be Gauss-eliminated without disrupting
//    the storage format, as long as no pivoting is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NA, the dimension of the array A.
//
//    Input, int DIAG[N], the indices in A of the N diagonal elements.
//
//    Input, double A[NA], the R8SS matrix.
//
//    Input, string TITLE, a title.
//
{
  r8ss_print_some ( n, na, diag, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8ss_print_some ( int n, int na, int diag[], double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8SS_PRINT_SOME prints some of a R8SS matrix.
//
//  Discussion:
//
//    The R8SS storage format is used for real symmetric skyline matrices.
//    This storage is appropriate when the nonzero entries of the
//    matrix are generally close to the diagonal, but the number
//    of nonzeroes above each diagonal varies in an irregular fashion.
//
//    In this case, the strategy is essentially to assign column J
//    its own bandwidth, and store the strips of nonzeros one after
//    another.   Note that what's important is the location of the
//    furthest nonzero from the diagonal.  A slot will be set up for
//    every entry between that and the diagonal, whether or not
//    those entries are zero.
//
//    A skyline matrix can be Gauss-eliminated without disrupting
//    the storage format, as long as no pivoting is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NA, the dimension of the array A.
//
//    Input, int DIAG[N], the indices in A of the N diagonal elements.
//
//    Input, double A[NA], the R8SS matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int ij;
  int ijm1;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        aij = 0.0;

        if ( j < i )
        {
          if ( i == 1 )
          {
            ijm1 = 0;
          }
          else
          {
            ijm1 = diag[i-2];
          }
          ij = diag[i-1];
          if ( ijm1 < ij+j-i )
          {
            aij = a[ij+j-i-1];
          }
        }
        else if ( j == i )
        {
          ij = diag[j-1];
          aij = a[ij-1];
        }
        else if ( i < j )
        {
          if ( j == 1 )
          {
            ijm1 = 0;
          }
          else
          {
            ijm1 = diag[j-2];
          }
          ij = diag[j-1];
          if ( ijm1 < ij+i-j )
          {
            aij = a[ij+i-j-1];
          }
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void r8ss_random ( int n, int *na, int diag[], double a[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8SS_RANDOM randomizes a R8SS matrix.
//
//  Discussion:
//
//    The R8SS storage format is used for real symmetric skyline matrices.
//    This storage is appropriate when the nonzero entries of the
//    matrix are generally close to the diagonal, but the number
//    of nonzeroes above each diagonal varies in an irregular fashion.
//
//    In this case, the strategy is essentially to assign column J
//    its own bandwidth, and store the strips of nonzeros one after
//    another.   Note that what's important is the location of the
//    furthest nonzero from the diagonal.  A slot will be set up for
//    every entry between that and the diagonal, whether or not
//    those entries are zero.
//
//    A skyline matrix can be Gauss-eliminated without disrupting
//    the storage format, as long as no pivoting is required.
//
//    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
//    although the actual storage needed will generally be about half of
//    that.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Output, int *NA, the dimension of the array A.
//    NA will be at least N and no greater than ( N * ( N + 1 ) ) / 2.
//
//    Output, int DIAG[N], the indices in A of the N diagonal elements.
//
//    Output, double A[(N*(N+1))/2], the R8SS matrix.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
{
  int diagold;
  int i;
  int ilo;
  int j;
  int k;
//
//  Set the values of DIAG.
//
  diag[0] = 1;
  *na = 1;
  for ( i = 1; i < n; i++ )
  {
    k = i4_uniform ( 1, i, seed );
    diag[i] = diag[i-1] + k;
    *na = *na + k;
  }
//
//  Now set the values of A.
//
  diagold = 0;
  k = 0;

  for ( j = 0; j < n; j++ )
  {
    ilo = j + diagold + 1 - diag[j];

    for ( i = ilo; i <= j; i++ )
    {
      a[k] = r8_uniform_01 ( seed );
      k = k + 1;
    }
    diagold = diag[j];
  }

  return;
}
//****************************************************************************80

double *r8ss_to_r8ge ( int n, int na, int diag[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SS_TO_R8GE copies a R8SS matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8SS storage format is used for real symmetric skyline matrices.
//    This storage is appropriate when the nonzero entries of the
//    matrix are generally close to the diagonal, but the number
//    of nonzeroes above each diagonal varies in an irregular fashion.
//
//    In this case, the strategy is essentially to assign column J
//    its own bandwidth, and store the strips of nonzeros one after
//    another.   Note that what's important is the location of the
//    furthest nonzero from the diagonal.  A slot will be set up for
//    every entry between that and the diagonal, whether or not
//    those entries are zero.
//
//    A skyline matrix can be Gauss-eliminated without disrupting
//    the storage format, as long as no pivoting is required.
//
//  Example:
//
//    11   0  13  0 15
//     0  22  23  0  0
//    31  32  33 34  0
//     0   0  43 44  0
//    51   0   0  0 55
//
//    A = ( 11 | 22 | 13, 23, 33 | 34, 44 | 15, 0, 0, 0, 55 )
//    NA = 12
//    DIAG = ( 1, 2, 5, 7, 12 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NA, the dimension of the array A.
//    NA must be at least N.
//
//    Input, int DIAG[N], the indices in A of the N diagonal elements.
//
//    Input, double A[NA], the R8SS matrix.
//
//    Output, double R8SS_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int diagold;
  int i;
  int ilo;
  int j;
  int k;

  b = new double[n*n];

  diagold = 0;
  k = 0;

  for ( j = 0; j < n; j++ )
  {
    ilo = j + diagold + 1 - diag[j];

    for ( i = 0; i < ilo; i++ )
    {
      b[i+j*n] = 0.0;
      b[j+i*n] = 0.0;
    }

    for ( i = ilo; i < j; i++ )
    {
      b[i+j*n] = a[k];
      b[j+i*n] = a[k];
      k = k + 1;
    }

    b[j+j*n] = a[k];
    k = k + 1;

    diagold = diag[j];
  }

  return b;
}
//****************************************************************************80

double *r8ss_zero ( int n, int na, int diag[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8SS_ZERO zeros a R8SS matrix.
//
//  Discussion:
//
//    The R8SS storage format is used for real symmetric skyline matrices.
//    This storage is appropriate when the nonzero entries of the
//    matrix are generally close to the diagonal, but the number
//    of nonzeroes above each diagonal varies in an irregular fashion.
//
//    In this case, the strategy is essentially to assign column J
//    its own bandwidth, and store the strips of nonzeros one after
//    another.   Note that what's important is the location of the
//    furthest nonzero from the diagonal.  A slot will be set up for
//    every entry between that and the diagonal, whether or not
//    those entries are zero.
//
//    A skyline matrix can be Gauss-eliminated without disrupting
//    the storage format, as long as no pivoting is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int NA, the dimension of the array A.
//    NA must be at least N.
//
//    Input, int DIAG[N], the indices in A of the N diagonal elements.
//
//    Output, double R8SS_ZERO[NA], the R8SS matrix.
//
{
  double *a;
  int diagold;
  int i;
  int ihi;
  int ilo;
  int j;
  int k;

  a = new double[na];

  diagold = 0;
  k = 0;

  for ( j = 0; j < n; j++ )
  {
    ilo = j + diagold + 1 - diag[j];
    ihi = j;

    for ( i = ilo; i <= ihi; i++ )
    {
      a[k] = 0.0;
      k = k + 1;
    }

    diagold = diag[j];
  }

  return a;
}
//****************************************************************************80

double *r8sto_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_INDICATOR sets up a R8STO indicator matrix.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Output, double R8STO_INDICATOR[N], the R8STO matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;
  int k;

  a = new double[n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  i = 1;
  k = 0;
  for ( j = 1; j <= n; j++ )
  {
    a[k] = ( double ) ( fac * i + j );
    k = k + 1;
  }
  
  return a;
}
//****************************************************************************80

double *r8sto_inverse ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_INVERSE computes the inverse of a R8STO matrix.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//    For this routine, the matrix is also required to be positive definite.
//
//    The original implementation of the algorithm assumed that the
//    diagonal element was 1.  The algorithm has been modified so that
//    this is no longer necessary.
//
//    The inverse matrix is NOT guaranteed to be a Toeplitz matrix.  
//    It is guaranteed to be symmetric and persymmetric.
//    The inverse matrix is returned in general storage, that is,
//    as an "SGE" matrix.
//
//  Example:
//
//    To compute the inverse of
//
//     1.0 0.5 0.2
//     0.5 1.0 0.5
//     0.2 0.5 1.0
//
//    we input:
//
//      N = 3
//      A = { 1.0, 0.5, 0.2 }
//
//    with output:
//
//      B = ( 1/56) * [ 75, -40,   5,
//                     -40,  96, -40,
//                       5, -40,  75 ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, Charles Van Loan,
//    Section 4.7.3, "Computing the Inverse",
//    Matrix Computations,
//    Third Edition,
//    Johns Hopkins, 1996.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, double A[N], the R8STO matrix.
//
//    Output, double R8STO_INVERSE[N*N], the inverse of the matrix.
//
{
  double *a2;
  double *b;
  int i;
  int j;
  double t;
  double *v;
  double vn;

  a2 = new double[n-1];
  b = new double[n*n];

  for ( i = 0; i < n-1; i++ )
  {
    a2[i] = a[i+1] / a[0];
  }

  v = r8sto_yw_sl ( n-1, a2 );
//
//  Compute the N-th entry of V.
//
  t = 0.0;
  for ( i = 0; i < n-1; i++ )
  {
    t = t + a2[i] * v[i];
  }
  vn = 1.0 / ( 1.0 + t );
//
//  Reverse the first N-1 entries of V.
//
  for ( i = 0; i < (n-1)/2; i++ )
  {
    j = n - 2 - i;
    t    = v[i];
    v[i] = v[j];
    v[j] = t;
  }
//
//  Scale the entries.
//
  for ( i = 0; i < n-1; i++ )
  {
    v[i] = vn * v[i];
  }
//
//  Set the boundaries of B.
//
  b[0+0*n] = vn;
  for ( j = 1; j < n; j++ )
  {
  b[0+j*n] = v[n-j-1];
  }

  for ( j = 0; j < n-1; j++ )
  {
    b[n-1+j*n] = v[j];
  }
  b[n-1+(n-1)*n] = vn;

  for ( i = 1; i < n-1; i++ )
  {
    b[i+0*n]     = v[n-1-i];
    b[i+(n-1)*n] = v[i];
  }
//
//  Fill the interior.
//
  for ( i = 2; i <= 1+(n-1)/2; i++ )
  {
    for ( j = i; j <= n - i + 1; j++ )
    {
      t = b[i-2+(j-2)*n] + ( v[n-j] * v[n-i] - v[i-2] * v[j-2] ) / vn;
      b[i-1+(j-1)*n] = t;
      b[j-1+(i-1)*n] = t;
      b[n-i+(n-j)*n] = t;
      b[n-j+(n-i)*n] = t;
    }
  }
//
//  Scale B.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i+j*n] = b[i+j*n] / a[0];
    }
  }

  delete [] a2;
  delete [] v;

  return b;
}
//****************************************************************************80

double *r8sto_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_MXV multiplies a R8STO matrix times a vector.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N], the R8STO matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8STO_MXV[N], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j <= i-1; j++ )
    {
      b[i] = b[i] + a[i-j] * x[j];
    }
    for ( j = i; j < n; j++ )
    {
      b[i] = b[i] + a[j-i] * x[j];
    }

  }

  return b;
}
//****************************************************************************80

void r8sto_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_PRINT prints a R8STO matrix.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N], the R8STO matrix.
//
//    Input, string TITLE, a title.
//
{
  r8sto_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8sto_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_PRINT_SOME prints some of am R8STO matrix.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N], the R8STO matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(4) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j )
        {
          aij = a[j-i];
        }
        else
        {
          aij = a[i-j];
        }
        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8sto_random ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_RANDOM randomizes a R8STO matrix.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8STO_RANDOM[N], the R8STO matrix.
//
{
  return ( r8vec_random ( n, 0.0, 1.0, seed ) );
}
//****************************************************************************80

double *r8sto_sl ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_SL solves a R8STO system.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//    The matrix is also required to be positive definite.
//
//    This implementation of the algorithm assumes that the diagonal element
//    (the first element of A) is 1.
//
//    Note that there is a typographical error in the presentation
//    of this algorithm in the reference, and another in the presentation
//    of a sample problem.  Both involve sign errors.  A minor error
//    makes the algorithm incorrect for the case N = 1.
//
//  Example:
//
//    To solve
//
//     1.0 0.5 0.2    x1    4.0
//     0.5 1.0 0.5 *  x2 = -1.0
//     0.2 0.5 1.0    x3    3.0
//
//    we input:
//
//      N = 3
//      A = (/ 1.0, 0.5, 0.2 /)
//      B = (/ 4.0, -1.0, 3.0 /)
//
//    with output:
//
//      X = (/ 355, -376, 285 /) / 56
//        = (/ 6.339, -6.714, 5.089 /)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, Charles Van Loan,
//    Section 4.7.3, "The General Right Hand Side Problem",
//    Matrix Computations,
//    Third Edition,
//    Johns Hopkins, 1996.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, double A[N], the R8STO matrix, with the EXTRA CONDITION
//    that the first entry is 1.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double R8STO_SL[N], the solution of the linear system.
//
{
  double beta;
  int i;
  int k;
  double *x;
  double *y;

  x = new double[n];
  y = new double[n];

  k = 0;
  beta = 1.0;
  x[k] = b[k] / beta;

  if ( k < n-1 )
  {
    y[k] = -a[k+1] / beta;
  }

  for ( k = 1; k <= n-1; k++ )
  {
    beta = ( 1.0 - y[k-1] * y[k-1] ) * beta;

    x[k] = b[k];
    for ( i = 1; i <= k; i++ )
    {
      x[k] = x[k] - a[i] * x[k-i];
    }
    x[k] = x[k] / beta;

    for ( i = 1; i <= k; i++ )
    {
      x[i-1] = x[i-1] + x[k] * y[k-i]; 
    }

    if ( k < n - 1 )
    {
      y[k] = -a[k+1];
      for ( i = 1; i <= k; i++ )
      {
        y[k] = y[k] - a[i] * y[k-i];
      }
      y[k] = y[k] / beta;

      for ( i = 1; i <= k; i++ )
      {
        y[i-1] = y[i-1] + y[k] * y[k-i];
      }
    }
  }

  delete [] y;

  return x;
}
//****************************************************************************80

double *r8sto_to_r8ge ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_TO_R8GE copies a R8STO matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double R8STO_TO_R8GE[N], the R8STO matrix.
//
//    Output, double R8STO_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      b[i+j*n] = a[i-j];
    }
    for ( j = i; j < n; j++ )
    {
      b[i+j*n] = a[j-i];
    }
  }

  return b;
}
//****************************************************************************80

double *r8sto_yw_sl ( int n, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_YW_SL solves the Yule-Walker equations for a R8STO matrix.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//    The matrix is also required to be positive definite.
//
//    This implementation of the algorithm assumes that the diagonal element
//    is 1.
//
//    The real symmetric Toeplitz matrix can be described by N numbers, which,
//    for convenience, we will label B(0:N-1).  We assume there is one more
//    number, B(N).  If we let A be the symmetric Toeplitz matrix whose first
//    row is B(0:N-1), then the Yule-Walker equations are:
//
//      A * X = -B(1:N)
//
//  Example:
//
//    To solve
//
//     1.0 0.5 0.2    x1   0.5
//     0.5 1.0 0.5 *  x2 = 0.2
//     0.2 0.5 1.0    x3   0.1
//
//    we input:
//
//      N = 3
//      B = (/ 0.5, 0.2, 0.1 /)
//
//    with output:
//
//      X = (/ -75, 12, -5 /) / 140
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, Charles Van Loan,
//    Section 4.7.2, "Solving the Yule-Walker Equations",
//    Matrix Computations,
//    Third Edition,
//    Johns Hopkins, 1996.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, double B[N], defines the linear system.  The first entry of the
//    symmetric Toeplitz matrix is assumed to be a 1, which is NOT stored.  The N-1
//    remaining elements of the first row of are stored in B, followed by
//    the remaining scalar that defines the linear system.
//
//    Output, double R8STO_YW_SL[N], the solution of the linear system.
//
{
  double alpha;
  double beta;
  int i;
  int j;
  double *x;
  double *x2;

  x = new double[n];
  x2 = new double[n];

  x[0] = -b[0];
  beta = 1.0;
  alpha = -b[0];

  for ( i = 1; i <= n-1; i++ )
  {
    beta = ( 1.0 - alpha * alpha ) * beta;

    alpha = b[i];
    for ( j = 1; j <= i; j++ )
    {
      alpha = alpha + b[i-j] * x[j-1];
    }

    alpha = -alpha / beta;

    for ( j = 1; j <= i; j++ )
    {
      x2[j-1] = x[j-1];
    }

    for ( j = 1; j <= i; j++ )
    {
      x[j-1] = x[j-1] + alpha * x2[i-j];
    }

    x[i] = alpha;
  }

  delete [] x2;

  return x;
}
//****************************************************************************80

double *r8sto_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8STO_ZERO zeros a R8STO matrix.
//
//  Discussion:
//
//    The R8STO storage format is used for a symmetric Toeplitz matrix.
//    It stores the N elements of the first row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Output, double R8STO_ZERO[N], the R8STO matrix.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }

  return a;
}
//****************************************************************************80

double *r8to_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_INDICATOR sets up a R8TO indicator matrix.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Output, double R8TO_INDICATOR[2*N-1], the R8TO matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;
  int k;

  a = new double[2*n-1];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  i = 1;
  k = 0;
  for ( j = 1; j <= n; j++ )
  {
    a[k] = ( double ) ( fac * i + j );
    k = k + 1;
  }

  j = 1;
  for ( i = 2; i <= n; i++ )
  {
    a[k] = ( double ) ( fac * i + j );
    k = k + 1;
  }
  
  return a;
}
//****************************************************************************80

double *r8to_mxv ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_MXV multiplies a R8TO matrix times a vector.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[2*N-1], the R8TO matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8TO_MXV[N], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < i; j++ )
    {
      b[i] = b[i] + a[n+i-j-1] * x[j];
    }
    for ( j = i; j < n; j++ )
    {
      b[i] = b[i] + a[j-i] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void r8to_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_PRINT prints a R8TO matrix.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[2*N-1], the R8TO matrix.
//
//    Input, string TITLE, a title.
//
{
  r8to_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8to_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_PRINT_SOME prints some of a R8TO matrix.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[2*N-1], the R8TO matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "  Col: ";

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(4) << i << "  ";

      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j ) 
        {
          cout << setw(12) << a[j-i] << "  ";
        }
        else
        {
          cout << setw(12) << a[n+i-j-1] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8to_random ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_RANDOM randomizes a R8TO matrix.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8TO_RANDOM[2*N-1], the R8TO matrix.
//
{
  double *a;
  int i;

  a = new double[2*n-1];

  for ( i = 0; i < 2*n-1; i++ )
  {
    a[i] = r8_uniform_01 ( seed );
  }

  return a;
}
//****************************************************************************80

double *r8to_sl ( int n, double a[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_SL solves a R8TO system.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2003
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[2*N-1], the R8TO matrix.
//
//    Input, double B[N] the right hand side vector.
//
//    Input, int JOB,
//    0 to solve A*X=B,
//    nonzero to solve A'*X=B.
//
//    Output, double R8TO_SL[N], the solution vector.  X and B may share the
//    same storage.
//
{
  double *c1;
  double *c2;
  int i;
  int nsub;
  double r1;
  double r2;
  double r3;
  double r5;
  double r6;
  double *x;

  if ( n < 1 )
  {
    return NULL;
  }

  x = new double[n];
//
//  Solve the system with the principal minor of order 1.
//
  r1 = a[0];
  x[0] = b[0] / r1;

  if ( n == 1 )
  {
    return x;
  }

  c1 = new double[n-1];
  c2 = new double[n-1];
//
//  Recurrent process for solving the system with the Toeplitz matrix.
//
  for ( nsub = 2; nsub <= n; nsub++ )
  {
//
//  Compute multiples of the first and last columns of the inverse of
//  the principal minor of order NSUB.
//
    if ( job == 0 )
    {
      r5 = a[n+nsub-2];
      r6 = a[nsub-1];
    }
    else
    {
      r5 = a[nsub-1];
      r6 = a[n+nsub-2];
    }

    if ( 2 < nsub )
    {
      c1[nsub-2] = r2;

      for ( i = 1; i <= nsub-2; i++ )
      {
        if ( job == 0 )
        {
          r5 = r5 + a[n+i-1] * c1[nsub-i-1];
          r6 = r6 + a[i] * c2[i-1];
        }
        else
        {
          r5 = r5 + a[i] * c1[nsub-i-1];
          r6 = r6 + a[n+i-1] * c2[i-1];
        }
      }
    }

    r2 = - r5 / r1;
    r3 = - r6 / r1;
    r1 = r1 + r5 * r3;

    if ( 2 < nsub )
    {
      r6 = c2[0];
      c2[nsub-2] = 0.0;

      for ( i = 2; i <= nsub-1; i++ )
      {
        r5 = c2[i-1];
        c2[i-1] = c1[i-1] * r3 + r6;
        c1[i-1] = c1[i-1] + r6 * r2;
        r6 = r5;
      }
    }

    c2[0] = r3;
//
//  Compute the solution of the system with the principal minor of order NSUB.
//
    if ( job == 0 )
    {
      r5 = 0.0;
      for ( i = nsub-1; 1 <= i; i-- )
      {
        r5 = r5 + a[n+nsub-i-1] * x[i-1];
      }
    }
    else
    {
      r5 = 0.0;
      for ( i = nsub-1; 1 <= i; i-- )
      {
        r5 = r5 + a[nsub-i] * x[i-1];
      }
    }

    r6 = ( b[nsub-1] - r5 ) / r1;

    for ( i = 0; i < nsub-1; i++ )
    {
      x[i] = x[i] + c2[i] * r6;
    }
    x[nsub-1] = r6;
  }

  delete [] c1;
  delete [] c2;

  return x;
}
//****************************************************************************80

double *r8to_to_r8ge ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_TO_R8GE copies a R8TO matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[2*N-1], the R8TO matrix.
//
//    Output, double R8TO_TO_R8GE[N*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      b[i+j*n] = a[n+i-j-1];
    }
    for ( j = i; j < n; j++ )
    {
      b[i+j*n] = a[j-i];
    }
  }

  return b;
}
//****************************************************************************80

double *r8to_vxm ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_VXM multiplies a vector times a R8TO matrix.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[2*N-1], the R8TO matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8TO_VXM[N], the product A' * X.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j <= i; j++ )
    {
      b[i] = b[i] + a[i-j] * x[j];
    }
    for ( j = i+1; j < n; j++ )
    {
      b[i] = b[i] + a[n+j-i-1] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

double *r8to_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8TO_ZERO zeros a R8TO matrix.
//
//  Discussion:
//
//    The R8TO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//    N must be positive.
//
//    Output, double R8TO_ZERO[2*N-1], the R8TO matrix.
//
{
  double *a;
  int i;

  a = new double[2*n-1];

  for ( i = 0; i < 2*n-1; i++ )
  {
    a[i] = 0.0;
  }

  return a;
}
//****************************************************************************80

double r8ut_det ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_DET computes the determinant of a R8UT matrix.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N*N], the R8UT matrix.
//
//    Output, double R8UT_DET, the determinant of the matrix.
//
{
  double det;
  int i;

  det = 1.0;

  for ( i = 0; i < n; i++ )
  {
    det = det * a[i+i*n];
  }

  return det;
}
//****************************************************************************80

double *r8ut_indicator ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_INDICATOR sets up a R8UT indicator matrix.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//    M and N must be positive.
//
//    Output, double R8UT_INDICATOR[M*N], the R8UT matrix.
//
{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[m*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= i4_min ( i-1, n ); j++ )
    {
      a[i-1+(j-1)*m] = 0.0;
    }
    for ( j = i; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = ( double ) ( fac * i + j );
    }
  }

  return a;
}
//****************************************************************************80

double *r8ut_inverse ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_INVERSE computes the inverse of a R8UT matrix.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the R8UT matrix.
//
//    Output, double R8UT_INVERSE[N*N], the inverse of the upper
//    triangular matrix.
//
{
  double *b;
  int i;
  int j;
  int k;
//
//  Check.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[i+i*n] == 0.0 )
    {
      cerr << "\n";
      cerr << "R8UT_INVERSE - Fatal error!\n";
      cerr << "  Zero diagonal element.\n";
      exit ( 1 );
    }
  }

  b = new double[n*n];

  for ( j = n-1; 0 <= j; j-- )
  {
    for ( i = n-1; 0 <= i; i-- )
    {
      if ( j < i )
      {
        b[i+j*n] = 0.0;
      }
      else if ( i == j )
      {
        b[i+j*n] = 1.0 / a[i+j*n];
      }
      else if ( i < j )
      {
        b[i+j*n] = 0.0;

        for ( k = i+1; k <= j; k++ )
        {
          b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
        }
        b[i+j*n] = b[i+j*n] / a[i+i*n];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8ut_mxm ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_MXM multiplies two R8UT matrices.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//    N must be positive.
//
//    Input, double A[N*N], B[N*N], the R8UT factor matrices.
//
//    Output, double R8UT_MXM[N*N], the R8UT product matrix.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = i; j < n; j++ )
    {
      c[i+j*n] = 0.0;
      for ( k = i; k <= j; k++ )
      {
        c[i+j*n] = c[i+j*n] + a[i+k*n] * b[k+j*n];
      }
    }
  }

  return c;
}
//****************************************************************************80

double *r8ut_mxv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_MXV multiplies a R8UT matrix times a vector.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2003
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
//    Input, double A[M*N], the R8UT matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8UT_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    for ( j = i; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void r8ut_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_PRINT prints a R8UT matrix.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[M*N], the R8UT matrix.
//
//    Input, string TITLE, a title.
//
{
  r8ut_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8ut_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_PRINT_SOME prints some of a R8UT matrix.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[M*N], the R8UT matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( ilo < jlo )
  {
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo );

    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(4) << i << "  ";

      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( j < i )
        {
          cout << "              ";
        }
        else
        {
          cout << setw(12) << a[i-1+(j-1)*m] << "  ";
        }
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8ut_random ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_RANDOM randomizes a R8UT matrix.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//    M and N must be positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8UT_RANDOM[M*N], the R8UT matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < i4_min ( i, n ); j++ )
    {
      a[i+j*m] = 0.0;
    }
    for ( j = i; j < n; j++ )
    {
      a[i+j*m] = r8_uniform_01 ( seed );
    }
  }

  return a;
}
//****************************************************************************80

double *r8ut_sl ( int n, double a[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_SL solves a R8UT system.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//    No factorization of the upper triangular matrix is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the R8UT matrix.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB, is 0 to solve the untransposed system,
//    nonzero to solve the transposed system.
//
//    Output, double R8UT_SL[N], the solution vector.
//
{
  int i;
  int j;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = n-1; 0 <= j; j-- )
    {
      x[j] = x[j] / a[j+j*n];
      for ( i = 0; i < j; i++ )
      {
        x[i] = x[i] - a[i+j*n] * x[j];
      }
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = x[j] / a[j+j*n];
      for ( i = j+1; i < n; i++ )
      {
        x[i] = x[i] - a[j+i*n] * x[j];
      }
    }
  }

  return x;
}
//****************************************************************************80

double *r8ut_vxm ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_VXM multiplies a vector times a R8UT matrix.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2003
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
//    Input, double A[M*N], the R8UT matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8UT_VXM[N], the product A' * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( j = 0; j < n; j++ )
  {
    b[j] = 0.0;
    for ( i = 0; i <= j && i < m; i++ )
    {
      b[j] = b[j] + x[i] * a[i+j*m];
    }
  }

  return b;
}
//****************************************************************************80

double *r8ut_zero ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8UT_ZERO zeros a R8UT matrix.
//
//  Discussion:
//
//    The R8UT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2003
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
//    Output, double R8UT_ZERO[M*N], the R8UT matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

void r8utp_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8UTP_PRINT prints a R8UTP matrix.
//
//  Discussion:
//
//    The R8UTP storage format is appropriate for an upper triangular
//    matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the matrix.
//
//    Input, string TITLE, a title.
//
{
  r8utp_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8utp_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8UTP_PRINT_SOME prints some of an R8UTP matrix.
//
//  Discussion:
//
//    The R8UTP storage format is appropriate for an upper triangular
//    matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j )
        {
          aij = a[i-1+(j*(j-1))/2];
        }
        else
        {
          aij = 0.0;
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double r8vec_dot_product ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of two R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements in the vectors.
//
//    Input, double X[N], Y[N], the two vectors.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  double dot;
  int i;

  dot = 0.0;

  for ( i = 0; i < n; i++ ) 
  {
    dot = dot + x[i] * y[i];
  }

  return dot;
}
//****************************************************************************80

double *r8vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDICATOR_NEW sets an R8VEC to the indicator vector {1,2,3...}.
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
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, double R8VEC_INDICATOR_NEW[N], the array to be initialized.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i <= n-1; i++ ) 
  {
    a[i] = ( double ) ( i + 1 );
  }

  return a;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  for ( i = 0; i < n; i++ ) 
  {
    cout << setw(6)  << i + 1 << "  " 
         << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_SOME prints "some" of an R8VEC.
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
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, integer I_LO, I_HI, the first and last indices to print.
//    The routine expects 1 <= I_LO <= I_HI <= N.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i-1]  << "\n";
  }

  return;
}
//****************************************************************************80

double *r8vec_random ( int n, double alo, double ahi, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_RANDOM randomizes an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double ALO, AHI, the range allowed for the entries.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_RANDOM[N], the vector of randomly chosen integers.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = r8_uniform ( alo, ahi, seed );
  }

  return a;
}
//****************************************************************************80

double *r8vec_read ( string input_file, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_READ reads an R8VEC from a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILE, the name of the file to be read.
//
//    Input, int N, the size of the R8VEC.
//
//    Output, double R(N), the R8VEC.
//
{
  ifstream input;
  int k;
  double *r;

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8VEC_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_file << "\"\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( k = 0; k < n; k++ )
  {
    input >> r[k];
  }

  input.close ( );

  return r;
}
//****************************************************************************80

int r8vec_read_size ( string input_file )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_READ_SIZE reads the size of an R8VEC from a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILE, the name of the file to 
//    be read.
//
//    Output, int R8VEC_READ_SIZE, the size of the R8VEC.
//
{
  ifstream input;
  char line[100];
  int n;

  n = 0;

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8VEC_READ_SIZE - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_file << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    n = n + 1;
  }

  input.close ( );

  return n;
}
//****************************************************************************80

double *r8vec_to_r8cb ( int m, int n, int ml, int mu, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TO_R8CB copies an R8VEC into a R8CB matrix.
//
//  Discussion:
//
//    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
//    a data item carries its dimensionality implicitly, and so cannot be
//    regarded sometimes as a vector and sometimes as an array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//
//    Input, double X[(ML+MU+1)*N], the vector to be copied into the array.
//
//    Output, double R8VEC_TO_R8CB[(ML+MU+1)*N], the array.
//
{
  double *a;
  int i;
  int j;

  a = new double[(ml+mu+1)*n];

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= ml + mu + 1; i++ )
    {
      if ( 1 <= i + j - mu - 1 & i + j - mu - 1 <= m )
      {
        a[i-1+(j-1)*(ml+mu+1)] = x[i-1+(j-1)*(ml+mu+1)];
      }
      else
      {
        a[i-1+(j-1)*(ml+mu+1)] = 0.0;
      }
    }
  }

  return a;
}
//****************************************************************************80

double *r8vec_to_r8gb ( int m, int n, int ml, int mu, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TO_R8GB copies an R8VEC into a R8GB matrix.
//
//  Discussion:
//
//    In C++  and FORTRAN, this routine is not really needed.  In MATLAB,
//    a data item carries its dimensionality implicitly, and so cannot be
//    regarded sometimes as a vector and sometimes as an array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//
//    Input, double X[(2*ML+MU+1)*N], the vector to be copied into the array.
//
//    Output, double R8VEC_TO_R8GB[(2*ML+MU+1)*N], the array.
//
{
  double *a;
  int i;
  int j;

  a = new double[(2*ml+mu+1)*n];

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= 2 * ml + mu + 1; i++ )
    {
      if ( 1 <= i + j - ml - mu - 1 & i + j - ml - mu - 1 <= m )
      {
        a[i-1+(j-1)*(2*ml+mu+1)] = x[i-1+(j-1)*(2*ml+mu+1)];
      }
      else
      {
        a[i-1+(j-1)*(2*ml+mu+1)] = 0.0;
      }

    }
  }

  return a;
}
//****************************************************************************80

double *r8vec_to_r8ge ( int m, int n, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TO_R8GE copies an R8VEC into a R8GE matrix.
//
//  Discussion:
//
//    In C++  and FORTRAN, this routine is not really needed.  In MATLAB,
//    a data item carries its dimensionality implicitly, and so cannot be
//    regarded sometimes as a vector and sometimes as an array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, double X[M*N], the vector to be copied into the array.
//
//    Output, double R8VEC_TO_R8GE[M*N], the array.
//
{
  double *a;
  int i;
  int j;
  int k;

  a = new double[m*n];

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = x[k];
      k = k + 1;
    }
  }

  return a;
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

void r8vec_write ( string output_filename, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_WRITE writes an R8VEC file.
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
//    10 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the data.
//
{
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8VEC_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    output << "  " << setw(24) << setprecision(16) << x[j] << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

void r8vec2_print_some ( int n, double x1[], double x2[], int max_print, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_PRINT_SOME prints "some" of two real vectors.
//
//  Discussion:
//
//    The user specifies MAX_PRINT, the maximum number of lines to print.
//
//    If N, the size of the vectors, is no more than MAX_PRINT, then
//    the entire vectors are printed, one entry of each per line.
//
//    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
//    followed by a line of periods suggesting an omission,
//    and the last entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters: 
//
//    Input, int N, the number of entries of the vectors.
//
//    Input, double X1[N], X2[N], the vector to be printed.
//
//    Input, int MAX_PRINT, the maximum number of lines to print.
//
//    Input, string TITLE, a title.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print-2; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
    cout << "......  ..............  ..............\n";
    i = n - 1;
    cout << setw(6)  << i + 1 << "  "
         << setw(14) << x1[i] << "  "
         << setw(14) << x2[i] << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
    i = max_print - 1;
    cout << setw(6)  << i + 1 << "  "
         << setw(14) << x1[i] << "  "
         << setw(14) << x2[i] << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

double r8vm_det ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_DET computes the determinant of a R8VM matrix.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Output, double R8VM_DET, the determinant of the matrix.
//
{
  double det;
  int i;
  int j;

  det = 1.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = j+1; i < n; i++ )
    {
      det = det * ( a[i] - a[j] );
    }
  }

  return det;
}
//****************************************************************************80

double *r8vm_mxv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_MXV multiplies a R8VM matrix times a vector.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8VM_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      if ( i == 0 )
      {
        b[i] = b[i] + x[j];
      }
      else
      {
        b[i] = b[i] + pow ( a[j], i ) * x[j];
      }
    }
  }

  return b;
}
//****************************************************************************80

void r8vm_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_PRINT prints a R8VM matrix.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Input, string TITLE, a title.
//
{
  r8vm_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8vm_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_PRINT_SOME prints some of a R8VM matrix.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i == 1 )
        {
          aij = 1.0;
        }
        else
        {
          aij = pow ( a[j-1], i-1 );
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *r8vm_random ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_RANDOM randomizes a R8VM matrix.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//    The parameter M is not actually needed by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VM_RANDOM[N], the R8VM matrix.
//
{
  double *a;

  a = r8vec_uniform_01_new ( n, seed );

  return a;
}
//****************************************************************************80

void r8vm_sl ( int n, double a[], double b[], int job, double x[], int *info )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_SL solves a R8VM system.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//    Vandermonde systems are very close to singularity.  The singularity
//    gets worse as N increases, and as any pair of values defining
//    the matrix get close.  Even a system as small as N = 10 will
//    involve the 9th power of the defining values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    Original FORTRAN77 version by Golub, VanLoan.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Gene Golub, Charles Van Loan,
//    Matrix Computations,
//    Third Edition,
//    Johns Hopkins, 1996.
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double X[N], the solution of the linear system.
//
//    Output, int *INFO.
//    0, no error.
//    nonzero, at least two of the values in A are equal.
//
{
  int i;
  int j;
//
//  Check for explicit singularity.
//
  *info = 0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = j+1; i < n; i++ )
    {
      if ( a[i] == a[j] )
      {
        *info = 1;
        return;
      }
    }
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = 1; j <= n-1; j++ )
    {
      for ( i = n; j+1 <= i; i-- )
      {
        x[i-1] = x[i-1] - a[j-1] * x[i-2];
      }
    }

    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j+1; i <= n; i++ )
      {
        x[i-1] = x[i-1] / ( a[i-1] - a[i-j-1] );
      }

      for ( i = j; i <= n-1; i++ )
      {
        x[i-1] = x[i-1] - x[i];
      }
    }
  }
  else
  {
    for ( j = 1; j <= n-1; j++ )
    {
      for ( i = n; j+1 <= i; i-- )
      {
        x[i-1] = ( x[i-1] - x[i-2] ) / ( a[i-1] - a[i-j-1] );
      }
    }

    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j; i <= n-1; i++ )
      {
        x[i-1] = x[i-1] - x[i] * a[j-1];
      }
    }

  }

  return;
}
//****************************************************************************80

double *r8vm_sl_new ( int n, double a[], double b[], int job, int *info )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_SL_NEW solves a R8VM system.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//    Vandermonde systems are very close to singularity.  The singularity
//    gets worse as N increases, and as any pair of values defining
//    the matrix get close.  Even a system as small as N = 10 will
//    involve the 9th power of the defining values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    Original FORTRAN77 version by Golub, VanLoan.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Gene Golub, Charles Van Loan,
//    Matrix Computations,
//    Third Edition,
//    Johns Hopkins, 1996.
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, int *INFO.
//    0, no error.
//    nonzero, at least two of the values in A are equal.
//
//    Output, double R8VM_SL_NEW[N], the solution of the linear system.
//
{
  int i;
  int j;
  double *x;
//
//  Check for explicit singularity.
//
  *info = 0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = j+1; i < n; i++ )
    {
      if ( a[i] == a[j] )
      {
        *info = 1;
        return NULL;
      }
    }
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = 1; j <= n-1; j++ )
    {
      for ( i = n; j+1 <= i; i-- )
      {
        x[i-1] = x[i-1] - a[j-1] * x[i-2];
      }
    }

    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j+1; i <= n; i++ )
      {
        x[i-1] = x[i-1] / ( a[i-1] - a[i-j-1] );
      }

      for ( i = j; i <= n-1; i++ )
      {
        x[i-1] = x[i-1] - x[i];
      }
    }
  }
  else
  {
    for ( j = 1; j <= n-1; j++ )
    {
      for ( i = n; j+1 <= i; i-- )
      {
        x[i-1] = ( x[i-1] - x[i-2] ) / ( a[i-1] - a[i-j-1] );
      }
    }

    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j; i <= n-1; i++ )
      {
        x[i-1] = x[i-1] - x[i] * a[j-1];
      }
    }

  }

  return x;
}
//****************************************************************************80

double *r8vm_to_r8ge ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_TO_R8GE copies a R8VM matrix to a R8GE matrix.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Output, double R8VM_TO_R8GE[M*N], the R8GE matrix.
//
{
  double *b;
  int i;
  int j;

  b = new double[m*n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == 0 )
      {
        b[i+j*m] = 1.0;
      }
      else
      {
        b[i+j*m] = b[i-1+j*m] * a[j];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8vm_vxm ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_VXM multiplies a vector times a R8VM matrix.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8VM_VXM[N], the product A' * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[n];

  for ( j = 0; j < n; j++ )
  {
    b[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      if ( i == 0 )
      {
        b[j] = b[j] + x[i];
      }
      else
      {
        b[j] = b[j] + pow ( a[j], i ) * x[i];
      }
    }
  }

  return b;
}
//****************************************************************************80

double *r8vm_zero ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_ZERO zeros a R8VM matrix.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Output, double R8VM_ZERO[M*N], the zero R8VM matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }

  return a;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n ) 
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Nijenhuis and Wilf.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }

    k = k - 1;
    k1 = k;

  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ;; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
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
