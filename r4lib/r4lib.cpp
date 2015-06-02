# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>
# include <cstring>

using namespace std;

# include "r4lib.hpp"

//****************************************************************************80

float c4_argument ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ARGUMENT returns the argument of a C4.
//
//  Discussion:
//
//    A C4 is a single precision complex value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose argument is desired.
//
//    Output, float C4_ARGUMENT, the argument of X.
//
{
  float value;

  if ( imag ( x ) == 0.0 && real ( x ) == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = atan2 ( imag ( x ), real ( x ) );
  }

  return value;
}
//****************************************************************************80

float c4_magnitude ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_MAGNITUDE returns the magnitude of a C4.
//
//  Discussion:
//
//    A C4 is a single precision complex value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_MAGNITUDE, the magnitude of X.
//
{
  float magnitude;

  magnitude = sqrt ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) );

  return magnitude;
}
//****************************************************************************80

complex <float> c4_sqrt ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SQRT returns the principal square root of a C4.
//
//  Discussion:
//
//    A C4 is a single precision complex value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the number whose square root is desired.
//
//    Output, complex <float> C4_SQRT, the square root of X.
//
{
  float argument;
  float magnitude;
  complex <float> value;

  argument = c4_argument ( x );
  magnitude = c4_magnitude ( x );

  if ( magnitude == 0.0 )
  {
    value = complex <float> ( 0.0, 0.0 );
  }
  else
  {
    value = sqrt ( magnitude )
      * complex <float> ( cos ( argument / 2.0 ), sin ( argument / 2.0 ) );
  }

  return value;
}
//****************************************************************************80

int i4_log_10 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
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

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  int value;

  if ( i < 0 )
  {
    value = -1;
  }
  else
  {
    value = 1;
  }
  return value;
}
//****************************************************************************80

int i4_uniform ( int a, int b, int &seed )

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
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
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

float i4int_to_r4int ( int imin, int imax, int i, float rmin, float rmax )

//****************************************************************************80
//
//  Purpose:
//
//    I4INT_TO_R4INT maps an I4 interval to an R4 interval.
//
//  Discussion:
//
//    The formula is
//
//      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IMIN, IMAX, the range.
//
//    Input, int I, the integer to be converted.
//
//    Input, float RMIN, RMAX, the range.
//
//    Output, float R, the corresponding value in [RMIN,RMAX].
//
{
  float r;

  if ( imax == imin )
  {
    r = 0.5 * ( rmin + rmax );
  }
  else
  {
    r = ( ( float ) ( imax - i        ) * rmin
        + ( float ) (        i - imin ) * rmax )
        / ( float ) ( imax     - imin );
  }

  return r;
}
//****************************************************************************80

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY copies an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int A2[N], the copy of A1.
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

int *i4vec_indicator0_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR0_NEW sets an I4VEC to the indicator vector {0,1,2,...}.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int I4VEC_INDICATOR0_NEW[N], the array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i;
  }
  return a;
}
//****************************************************************************80

int *i4vec_indicator1_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR1_NEW sets an I4VEC to the indicator vector {1,2,3,...}.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int I4VEC_INDICATOR1_NEW[N], the array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return a;
}
//****************************************************************************80

void i4vec_permute ( int n, int p[], int base, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PERMUTE permutes an I4VEC in place.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    This routine permutes an array of integer "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   1,   3,   4,   0,   2 )
//      A = (   1,   2,   3,   4,   5 )
//
//    Output:
//
//      A    = (   2,   4,   5,   1,   3 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for
//    a 1-based permutation.
//
//    Input/output, int A[N], the array to be permuted.
//
{
  int a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "I4VEC_PERMUTE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
//  So temporarily add 1-BASE to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "I4VEC_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the base of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 + base;
  }

  return;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
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
//    Input, int A[N], the vector to be printed.
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
    cout << "  " << setw(8) << i
         << ": " << setw(8) << a[i]  << "\n";
  }
  return;
}
//****************************************************************************80

void i4vec_zero ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}
//****************************************************************************80

int *i4vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}
//****************************************************************************80

bool perm_check ( int n, int p[], int base )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from BASE to
//    to BASE+N-1 occurs among the N entries of the permutation.
//
//    Set the input quantity BASE to 0, if P is a 0-based permutation,
//    or to 1 if P is a 1-based permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Input, int BASE, the index base.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
  bool found;
  int i;
  int seek;

  for ( seek = base; seek < base + n; seek++ )
  {
    found = false;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = true;
        break;
      }
    }

    if ( !found )
    {
      return false;
    }

  }

  return true;
}
//****************************************************************************80

int *perm_uniform_new ( int n, int base, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM_NEW selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
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
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for
//    a 1-based permutation.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int PERM_UNIFORM_NEW[N], a permutation of
//    (BASE, BASE+1, ..., BASE+N-1).
//
{
  int i;
  int j;
  int k;
  int *p;

  p = new int[n];

  for ( i = 0; i < n; i++ )
  {
    p[i] = i + base;
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_uniform ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
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
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

float r4_add ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ADD adds two R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, Y, the numbers to be added.
//
//    Output, float R4_ADD, the sum of X and Y.
//
{
  float value;

  value = x + y;

  return value;
}
//****************************************************************************80

float r4_atan ( float y, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ATAN computes the inverse tangent of the ratio Y / X.
//
//  Discussion:
//
//    R4_ATAN returns an angle whose tangent is ( Y / X ), a job which
//    the built in functions ATAN and ATAN2 already do.
//
//    However:
//
//    * R4_ATAN always returns a positive angle, between 0 and 2 PI,
//      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
//      and [-PI,+PI] respectively;
//
//    * R4_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
//     function by contrast always returns an angle in the first or fourth
//     quadrants.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float Y, X, two quantities which represent the tangent of
//    an angle.  If Y is not zero, then the tangent is (Y/X).
//
//    Output, float R4_ATAN, an angle between 0 and 2 * PI, whose tangent is
//    (Y/X), and which lies in the appropriate quadrant so that the signs
//    of its cosine and sine match those of X and Y.
//
{
  float abs_x;
  float abs_y;
  const float r4_pi = 3.141592653589793;
  float theta;
  float theta_0;
//
//  Special cases:
//
  if ( x == 0.0 )
  {
    if ( 0.0 < y )
    {
      theta = r4_pi / 2.0;
    }
    else if ( y < 0.0 )
    {
      theta = 3.0 * r4_pi / 2.0;
    }
    else if ( y == 0.0 )
    {
      theta = 0.0;
    }
  }
  else if ( y == 0.0 )
  {
    if ( 0.0 < x )
    {
      theta = 0.0;
    }
    else if ( x < 0.0 )
    {
      theta = r4_pi;
    }
  }
//
//  We assume that ATAN2 is correct when both arguments are positive.
//
  else
  {
    abs_y = r4_abs ( y );
    abs_x = r4_abs ( x );

    theta_0 = atan2 ( abs_y, abs_x );

    if ( 0.0 < x && 0.0 < y )
    {
      theta = theta_0;
    }
    else if ( x < 0.0 && 0.0 < y )
    {
      theta = r4_pi - theta_0;
    }
    else if ( x < 0.0 && y < 0.0 )
    {
      theta = r4_pi + theta_0;
    }
    else if ( 0.0 < x && y < 0.0 )
    {
      theta = 2.0 * r4_pi - theta_0;
    }
  }
  return theta;
}
//****************************************************************************80

float r4_big ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BIG returns a "big" R4.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R4.  We simply want a "very large" but 
//    non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, float R4_BIG, a "big" R4 value.
//
{
  float value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

float r4_cas ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CAS returns the "casine" of an R4.
//
//  Discussion:
//
//    The "casine", used in the discrete Hartley transform, is abbreviated
//    CAS(X), and defined by:
//
//      CAS(X) = cos ( X ) + sin( X )
//             = sqrt ( 2 ) * sin ( X + pi/4 )
//             = sqrt ( 2 ) * cos ( X - pi/4 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose casine is desired.
//
//    Output, float R4_CAS, the casine of X, which will be between
//    plus or minus the square root of 2.
//
{
  float value;

  value = cos ( x ) + sin ( x );

  return value;
}
//****************************************************************************80

float r4_ceiling ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CEILING rounds an R4 "up" (towards +oo) to the next integer.
//
//  Example:
//
//    X        R4_CEILING(X)
//
//   -1.1      -1.0
//   -1.0      -1.0
//   -0.9       0.0
//   -0.1       0.0
//    0.0       0.0
//    0.1       1.0
//    0.9       1.0
//    1.0       1.0
//    1.1       2.0
//    2.9       3.0
//    3.0       3.0
//    3.14159   4.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose ceiling is desired.
//
//    Output, float R4_CEILING, the ceiling of X.
//
{
  float value;

  value = ( int ) x;

  if ( value < x )
  {
    value = value + 1.0;
  }

  return value;
}
//****************************************************************************80

float r4_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CHOOSE computes the binomial coefficient C(N,K) as an R4.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R4 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, float R4_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( float ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( float ) ( mx + i ) ) / ( float ) i;
    }
  }

  return value;
}
//****************************************************************************80

float r4_chop ( int place, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CHOP chops an R4 to a given number of binary places.
//
//  Example:
//
//    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.
//
//    The following values would be returned for the 'chopped' value of
//    3.875:
//
//    PLACE  Value
//
//       1      2
//       2      3     = 2 + 1
//       3      3.5   = 2 + 1 + 1/2
//       4      3.75  = 2 + 1 + 1/2 + 1/4
//       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PLACE, the number of binary places to preserve.
//    PLACE = 0 means return the integer part of X.
//    PLACE = 1 means return the value of X, correct to 1/2.
//    PLACE = 2 means return the value of X, correct to 1/4.
//    PLACE = -1 means return the value of X, correct to 2.
//
//    Input, float X, the number to be chopped.
//
//    Output, float R4_CHOP, the chopped number.
//
{
  float fac;
  int temp;
  float value;

  temp = ( int ) ( r4_log_2 ( x ) );
  fac = pow ( 2.0, ( temp - place + 1 ) );
  value = ( float ) ( ( int ) ( x / fac ) ) * fac;

  return value;
}
//****************************************************************************80

complex <float> r4_csqrt ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CSQRT returns the complex square root of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose square root is desired.
//
//    Output, complex <float> R4_CSQRT, the square root of X:
//
{
  float argument;
  float magnitude;
  const float r4_pi = 3.141592653589793;
  complex <float> value;

  if ( 0.0 < x )
  {
    magnitude = x;
    argument = 0.0;
  }
  else if ( 0.0 == x )
  {
    magnitude = 0.0;
    argument = 0.0;
  }
  else if ( x < 0.0 )
  {
    magnitude = -x;
    argument = r4_pi;
  }

  magnitude = sqrt ( magnitude );
  argument = argument / 2.0;

  value = magnitude * complex <float> ( cos ( argument ), sin ( argument ) );

  return value;
}
//****************************************************************************80

float r4_cube_root ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CUBE_ROOT returns the cube root of an R4.
//
//  Discussion:
//
//    This routine is designed to avoid the possible problems that can occur
//    when formulas like 0.0**(1/3) or (-1.0)**(1/3) are to be evaluated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input float X, the number whose cube root is desired.
//
//    Output, float R4_CUBE_ROOT, the cube root of X.
//
{
  float e;
  float value;

  if ( 0.0 < x )
  {
    e = 1.0 / 3.0;
    value = pow ( x, e );
  }
  else if ( x == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    e = 1.0 / 3.0;
    value = - pow ( r4_abs ( x ), e );
  }

  return value;
}
//****************************************************************************80

float r4_diff ( float x, float y, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4_DIFF computes (X-Y) to a specified accuracy.
//
//  Discussion:
//
//    The user controls how many binary digits of accuracy
//    are to be used.
//
//    N determines the accuracy of the value.  If N = 10,
//    for example, only 11 binary places will be used in the arithmetic.
//    In general, only N+1 binary places will be used.
//
//    N may be zero.  However, a negative value of N should
//    not be used, since this will cause both X and Y to look like 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, Y, the two values whose difference is desired.
//
//    Input, int N, the number of binary digits to use.
//
//    Output, float R4_DIFF, the value of X-Y.
//
{
  float cx;
  float cy;
  float pow2;
  float size;
  float value;

  if ( x == y )
  {
    value = 0.0;
    return value;
  }

  pow2 = pow ( 2.0, n );
//
//  Compute the magnitude of X and Y, and take the larger of the
//  two.  At least one of the two values is not zero//
//
  size = r4_max ( r4_abs ( x ), r4_abs ( y ) );
//
//  Make normalized copies of X and Y.  One of the two values will
//  actually be equal to 1.
//
  cx = x / size;
  cy = y / size;
//
//  Here's where rounding comes in.  We know that the larger of the
//  the two values equals 1.  We multiply both values by 2**N,
//  where N+1 is the number of binary digits of accuracy we want
//  to use, truncate the values, and divide back by 2**N.
//
  cx = ( float ) ( ( int ) ( cx * pow2 + 0.5 * r4_sign ( cx ) ) ) / pow2;
  cy = ( float ) ( ( int ) ( cy * pow2 + 0.5 * r4_sign ( cy ) ) ) / pow2;
//
//  Take the difference now.
//
  value = cx - cy;
//
//  Undo the scaling.
//
  value = value * size;

  return value;
}
//****************************************************************************80

int r4_digit ( float x, int idigit )

//****************************************************************************80
//
//  Purpose:
//
//    R4_DIGIT returns a particular decimal digit of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose IDIGIT-th decimal digit is desired.
//    Note that if X is zero, all digits will be returned as 0.
//
//    Input, int IDIGIT, the position of the desired decimal digit.
//    A value of 1 means the leading digit, a value of 2 the second digit
//    and so on.
//
//    Output, int R4_DIGIT, the value of the IDIGIT-th decimal digit of X.
//
{
  int digit;
  int i;
  int ival;

  if ( x == 0.0 )
  {
    digit = 0;
    return digit;
  }

  if ( idigit <= 0 )
  {
    digit = 0;
    return digit;
  }
//
//  Force X to lie between 1 and 10.
//
  x = r4_abs ( x );

  while ( x < 1.0 )
  {
    x = x * 10.0;
  }

  while ( 10.0 <= x )
  {
    x = x / 10.0;
  }

  for ( i = 1; i <= idigit; i++ )
  {
    ival = ( int ) ( x );
    x = ( x - ( float ) ival ) * 10.0;
  }

  digit = ival;

  return digit;
}
//****************************************************************************80

float r4_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_EPSILON returns the R4 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
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
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R4_EPSILON, the R4 round-off unit.
//
{
  const float value = 1.19209290E-07;

  return value;
}
//****************************************************************************80

float r4_epsilon_compute ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_EPSILON_COMPUTE computes the R4 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
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
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, float R4_EPSILON_COMPUTE, the R4 round-off unit.
//
{
  float one;
  float temp;
  float test;
  float value;

  one = ( float ) ( 1 );

  value = one;
  temp = value / 2.0;
  test = r4_add ( one, temp );

  while ( one < test )
  {
    value = temp;
    temp = value / 2.0;
    test = r4_add ( one, temp );
  }

  return value;
}
//****************************************************************************80

float r4_exp ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_EXP computes the exponential function, avoiding overflow and underflow.
//
//  Discussion:
//
//    For arguments of very large magnitude, the evaluation of the
//    exponential function can cause computational problems.  Some languages
//    and compilers may return an infinite value or a "Not-a-Number".  
//    An alternative, when dealing with a wide range of inputs, is simply
//    to truncate the calculation for arguments whose magnitude is too large.
//    Whether this is the right or convenient approach depends on the problem
//    you are dealing with, and whether or not you really need accurate
//    results for large magnitude inputs, or you just want your code to
//    stop crashing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the argument of the exponential function.
//
//    Output, float R4_EXP, the value of exp ( X ).
//
{
  const float r4_big = 1.0E+30;
  const float r4_log_max = +69.0776;
  const float r4_log_min = -69.0776;
  float value;

  if ( x <= r4_log_min )
  {
    value = 0.0;
  }
  else if ( x < r4_log_max )
  {
    value = exp ( x );
  }
  else
  {
    value = r4_big;
  }

  return value;
}
//****************************************************************************80

float r4_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4_FACTORIAL computes the factorial of N.
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
//    11 August 2010
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
//    Output, float R4_FACTORIAL, the factorial of N.
//
{
  int i;
  float value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( float ) ( i );
  }

  return value;
}
//****************************************************************************80

float r4_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4_FACTORIAL2 computes the float factorial function.
//
//  Discussion:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    Factorial2(N)
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
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the float factorial
//    function.  If N is less than 1, R4_FACTORIAL2 is returned as 1.0.
//
//    Output, float R4_FACTORIAL2, the value of Factorial2(N).
//
{
  int n_copy;
  float value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( float ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
//****************************************************************************80

float r4_floor ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_FLOOR rounds an R4 "down" (towards -oo) to the next integer.
//
//  Example:
//
//    X        R4_FLOOR(X)
//
//   -1.1      -2.0
//   -1.0      -1.0
//   -0.9      -1.0
//   -0.1      -1.0
//    0.0       0.0
//    0.1       0.0
//    0.9       0.0
//    1.0       1.0
//    1.1       1.0
//    2.9       2.0
//    3.0       3.0
//    3.14159   3.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose floor is desired.
//
//    Output, float R4_FLOOR, the floor of X.
//
{
  float value;

  value = ( int ) x;

  if ( x < value )
  {
    value = value - 1.0;
  }

  return value;
}
//****************************************************************************80

float r4_fraction ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    R4_FRACTION uses real arithmetic on an integer ratio.
//
//  Discussion:
//
//    Given integer variables I and J, both FORTRAN and C will evaluate
//    an expression such as "I/J" using what is called "integer division",
//    with the result being an integer.  It is often convenient to express
//    the parts of a fraction as integers but expect the result to be computed
//    using real arithmetic.  This function carries out that operation.
//
//  Example:
//
//       I     J   I/J  R4_FRACTION
//
//       1     2     0  0.5
//       7     4     1  1.75
//       8     4     2  2.00
//       9     4     2  2.25
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the arguments.
//
//    Output, float R4_FRACTION, the value of the ratio.
//
{
  float value;

  value = ( float ) ( i ) / ( float ) ( j );

  return value;
}
//****************************************************************************80

float r4_fractional ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_FRACTIONAL returns the fractional part of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_FRACTIONAL, the fractional part of X.
//
{
  float value;

  value = r4_abs ( x ) - ( float ) ( ( int ) r4_abs ( x ) );

  return value;
}
//****************************************************************************80

float r4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_HUGE returns a "huge" R4.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R4.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, float R4_HUGE, a "huge" R4 value.
//
{
  float value;

  value = 3.40282347E+38;

  return value;
}
//****************************************************************************80

bool r4_in_01 ( float a )

//****************************************************************************80
//
//  Purpose:
//
//    R4_IN_01 is TRUE if an R4 is in the range [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, the value.
//
//    Output, bool R4_IN_01, is TRUE if A is between 0 and 1.
//
{
  bool value;

  value = ( 0.0 <= a && a <= 1.0 );

  return value;
}
//****************************************************************************80

bool r4_is_int ( float r )

//****************************************************************************80
//
//  Purpose:
//
//    R4_IS_INT determines if an R4 represents an integer value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R, the number to be checked.
//
//    Output, bool R4_IS_INT, is TRUE if R is an integer value.
//
{
  const int i4_huge = 2147483647;
  bool value;

  if ( ( float ) ( i4_huge ) < r )
  {
    value = false;
  }
  else if ( r < - ( float ) ( i4_huge ) )
  {
    value = false;
  }
  else if ( r == ( float ) ( ( int ) ( r ) ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

float r4_log_10 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LOG_10 returns the logarithm base 10 of the absolute value of an R4.
//
//  Discussion:
//
//    value = Log10 ( |X| )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose base 2 logarithm is desired.
//    X should not be 0.
//
//    Output, float R4_LOG_10, the logarithm base 10 of the absolute
//    value of X.  It should be true that |X| = 10^R_LOG_10.
//
{
  float value;

  if ( x == 0.0 )
  {
    value = - r4_big ( );
  }
  else
  {
    value = log10 ( r4_abs ( x ) );
  }

  return value;
}
//****************************************************************************80

float r4_log_2 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LOG_2 returns the logarithm base 2 of the absolute value of an R4.
//
//  Discussion:
//
//    value = Log ( |X| ) / Log ( 2.0 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose base 2 logarithm is desired.
//    X should not be 0.
//
//    Output, float R4_LOG_2, the logarithm base 2 of the absolute
//    value of X.  It should be true that |X| = 2^R_LOG_2.
//
{
  float value;

  if ( x == 0.0 )
  {
    value = - r4_big ( );
  }
  else
  {
    value = log ( r4_abs ( x ) ) / log ( 2.0 );
  }

  return value;
}
//****************************************************************************80

float r4_log_b ( float x, float b )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LOG_B returns the logarithm base B of an R4.
//
//  Discussion:
//
//    value = log ( |X| ) / log ( |B| )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose base B logarithm is desired.
//    X should not be 0.
//
//    Input, float B, the base, which should not be 0, 1 or -1.
//
//    Output, float R4_LOG_B, the logarithm base B of the absolute
//    value of X.  It should be true that |X| = |B|^R_LOG_B.
//
{
  float value;

  if ( b == 0.0 || b == 1.0 || b == -1.0 )
  {
    value = - r4_big ( );
  }
  else if ( r4_abs ( x ) == 0.0 )
  {
    value = - r4_big ( );
  }
  else
  {
    value = log ( r4_abs ( x ) ) / log ( r4_abs ( b ) );
  }

  return value;
}
//****************************************************************************80

void r4_mant ( float x, int *s, float *r, int *l )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MANT computes the "mantissa" or "fraction part" of an R4.
//
//  Discussion:
//
//    X = S * R * 2^L
//
//    S is +1 or -1,
//    R is a real between 1.0 and 2.0,
//    L is an integer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the real number to be decomposed.
//
//    Output, int *S, the "sign" of the number.
//    S will be -1 if X is less than 0, and +1 if X is greater
//    than or equal to zero.
//
//    Output, float *R, the mantissa of X.  R will be greater
//    than or equal to 1, and strictly less than 2.  The one
//    exception occurs if X is zero, in which case R will also
//    be zero.
//
//    Output, int *L, the integer part of the logarithm (base 2) of X.
//
{
//
//  Determine the sign.
//
  if ( x < 0.0 )
  {
    *s = -1;
  }
  else
  {
    *s = 1;
  }
//
//  Set R to the absolute value of X, and L to zero.
//  Then force R to lie between 1 and 2.
//
  if ( x < 0.0 )
  {
    *r = -x;
  }
  else
  {
    *r = x;
  }

  *l = 0;
//
//  Time to bail out if X is zero.
//
  if ( x == 0.0 )
  {
    return;
  }

  while ( 2.0 <= *r )
  {
    *r = *r / 2.0;
    *l = *l + 1;
  }

  while ( *r < 1.0 )
  {
    *r = *r * 2.0;
    *l = *l - 1;
  }

  return;
}
//****************************************************************************80

float r4_max ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MAX returns the maximum of two R4's.
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
//    Input, float X, Y, the quantities to compare.
//
//    Output, float R4_MAX, the maximum of X and Y.
//
{
  float value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

float r4_min ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MIN returns the minimum of two R4's..
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
//    Input, float X, Y, the quantities to compare.
//
//    Output, float R_MIN, the minimum of X and Y.
//
{
  float value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

float r4_mod ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MOD returns the remainder of R4 division.
//
//  Discussion:
//
//    If
//      REM = R4_MOD ( X, Y )
//      RMULT = ( X - REM ) / Y
//    then
//      X = Y * RMULT + REM
//    where REM has the same sign as X, and abs ( REM ) < Y.
//
//  Example:
//
//        X         Y     R4_MOD   R4_MOD  Factorization
//
//      107        50       7     107 =  2 *  50 + 7
//      107       -50       7     107 = -2 * -50 + 7
//     -107        50      -7    -107 = -2 *  50 - 7
//     -107       -50      -7    -107 =  2 * -50 - 7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number to be divided.
//
//    Input, float Y, the number that divides X.
//
//    Output, float R4_MOD, the remainder when X is divided by Y.
//
{
  float value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R4_MOD - Fatal error!\n";
    cerr << "  R4_MOD ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( float ) ( ( int ) ( x / y ) ) ) * y;

  if ( x < 0.0 && 0.0 < value )
  {
    value = value - r4_abs ( y );
  }
  else if ( 0.0 < x && value < 0.0 )
  {
    value = value + r4_abs ( y );
  }

  return value;
}
//****************************************************************************80

float r4_modp ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MODP returns the nonnegative remainder of R4 division.
//
//  Discussion:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360.0) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, R4_MODP(A,360.0) is between 0 and 360, always.
//
//    If
//      REM = R4_MODP ( X, Y )
//      RMULT = ( X - REM ) / Y
//    then
//      X = Y * RMULT + REM
//    where REM is always nonnegative.
//
//  Example:
//
//        I         J     MOD  R4_MODP   R4_MODP Factorization
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
//    22 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number to be divided.
//
//    Input, float Y, the number that divides X.
//
//    Output, float R4_MODP, the nonnegative remainder when X is divided by Y.
//
{
  float value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R4_MODP - Fatal error!\n";
    cerr << "  R4_MODP ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( float ) ( ( int ) ( x / y ) ) ) * y;

  if ( value < 0.0 )
  {
    value = value + r4_abs ( y );
  }

  return value;
}
//****************************************************************************80

float r4_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MOP returns the I-th power of -1 as an R4 value.
//
//  Discussion:
//
//    An R4 is an float value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, float R4_MOP, the I-th power of -1.
//
{
  float value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
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
    value = - ( int ) ( - x + 0.5 );
  }
  else
  {
    value =   ( int ) (  x + 0.5 );
  }
  return value;
}
//****************************************************************************80

float r4_normal ( float a, float b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NORMAL returns a scaled pseudonormal R4.
//
//  Discussion:
//
//    The normal probability distribution function (PDF) is sampled,
//    with mean A and standard deviation B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, the mean of the PDF.
//
//    Input, float B, the standard deviation of the PDF.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, float R4_NORMAL, a sample of the normal PDF.
//
{
  float value;

  value = a + b * r4_normal_01 ( seed );

  return value;
}
//****************************************************************************80

float r4_normal_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NORMAL_01 returns a unit pseudonormal R4.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, float R4_NORMAL_01, a normally distributed random value.
//
{
  float r1;
  float r2;
  const float r4_pi = 3.141592653589793;
  static int used = -1;
  static int seed1 = 0;
  static int seed2 = 0;
  static int seed3 = 0;
  float x;
  static float y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
//
//  If we've used an even number of values so far, generate two more, return one,
//  and save one.
//
  if ( ( used % 2 ) == 0 )
  {
    seed1 = seed;
    r1 = r4_uniform_01 ( seed );

    if ( r1 == 0.0 )
    {
      cerr << "\n";
      cerr << "R4_NORMAL_01 - Fatal error!\n";
      cerr << "  R4_UNIFORM_01 returned a value of 0.\n";
      exit ( 1 );
    }

    seed2 = seed;
    r2 = r4_uniform_01 ( seed );
    seed3 = seed;
    seed = seed2;

    x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r4_pi * r2 );
    y = sqrt ( - 2.0 * log ( r1 ) ) * sin ( 2.0 * r4_pi * r2 );
  }
//
//  Otherwise, return the second, saved, value and the corresponding
//  value of SEED.
//
  else
  {
    x = y;
    seed = seed3;
  }

  used = used + 1;

  return x;
# undef PI
}
//****************************************************************************80

float r4_pi ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_PI returns the value of PI as an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, float R4_PI, the value of PI.
//
{
  float value = 3.141592653589793;

  return value;
}
//****************************************************************************80

float r4_power ( float r, int p )

//****************************************************************************80
//
//  Purpose:
//
//    R4_POWER computes an integer power of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R, the base.
//
//    Input, int P, the power, which may be negative.
//
//    Output, float R4_POWER, the value of R^P.
//
{
  float value;
//
//  Special case.  R^0 = 1.
//
  if ( p == 0 )
  {
    value = 1.0;
  }
//
//  Special case.  Positive powers of 0 are 0.
//  We go ahead and compute negative powers, relying on the software to complain.
//
  else if ( r == 0.0 )
  {
    if ( 0 < p )
    {
      value = 0.0;
    }
    else
    {
      value = pow ( r, p );
    }
  }
  else if ( 1 <= p )
  {
    value = pow ( r, p );
  }
  else
  {
    value = pow ( r, p );
  }

  return value;
}
//****************************************************************************80

float r4_power_fast ( float r, int p, int *mults )

//****************************************************************************80
//
//  Purpose:
//
//    R4_POWER_FAST computes the P-th power of R, for real R and integer P.
//
//  Discussion:
//
//    Obviously, R^P can be computed using P-1 multiplications.
//
//    However, R^P can also be computed using at most 2*LOG2(P) multiplications.
//    To do the calculation this way, let N = LOG2(P).
//    Compute A, A^2, A^4, ..., A^N by N-1 successive squarings.
//    Start the value of R^P at A, and each time that there is a 1 in
//    the binary expansion of P, multiply by the current result of the squarings.
//
//    This algorithm is not optimal.  For small exponents, and for special
//    cases, the result can be computed even more quickly.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R, the base.
//
//    Input, int P, the power, which may be negative.
//
//    Output, int *MULTS, the number of multiplications and divisions.
//
//    Output, float R4_POWER_FAST, the value of R^P.
//
{
  int p_mag;
  int p_sign;
  float r2;
  float value;

  *mults = 0;
//
//  Special bases.
//
  if ( r == 1.0 )
  {
    value = 1.0;
    return value;
  }

  if ( r == -1.0 )
  {
    if ( ( p % 2 ) == 1 )
    {
      value = -1.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  if ( r == 0.0 )
  {
    if ( p <= 0 )
    {
      cerr << "\n";
      cerr << "R4_POWER_FAST - Fatal error!\n";
      cerr << "  Base is zero, and exponent is negative.\n";
      exit ( 1 );
    }

    value = 0.0;
    return value;
  }
//
//  Special powers.
//
  if ( p == -1 )
  {
    value = 1.0 / r;
    *mults = *mults + 1;
    return value;
  }
  else if ( p == 0 )
  {
    value = 1.0;
    return value;
  }
  else if ( p == 1 )
  {
    value = r;
    return value;
  }
//
//  Some work to do.
//
  p_mag = abs ( p );
  p_sign = i4_sign ( p );

  value = 1.0;
  r2 = r;

  while ( 0 < p_mag )
  {
    if ( ( p_mag % 2 ) == 1 )
    {
      value = value * r2;
      *mults = *mults + 1;
    }

    p_mag = p_mag / 2;
    r2 = r2 * r2;
    *mults = *mults + 1;
  }

  if ( p_sign == -1 )
  {
    value = 1.0 / value;
    *mults = *mults + 1;
  }

  return value;
}
//****************************************************************************80

float r4_pythag ( float a, float b )

//****************************************************************************80
//
//  Purpose:
//
//    R4_PYTHAG computes sqrt ( A*A + B*B ), avoiding overflow and underflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, B, the values for which sqrt ( A*A + B*B ) is desired.
//
//    Output, float R4_PYTHAG, the value of sqrt ( A*A + B*B ).
//
{
  float a_abs;
  float b_abs;
  float result;

  a_abs = r4_abs ( a );
  b_abs = r4_abs ( b );

  if ( b_abs < a_abs )
  {
    result = a_abs * sqrt ( 1.0 + pow ( b_abs / a_abs, 2 ) );
  }
  else if ( b_abs == 0.0 )
  {
    result = 0.0;
  }
  else if ( a_abs <= b_abs )
  {
    result = b_abs * sqrt ( 1.0 + pow ( a_abs / b_abs, 2 ) );
  }

  return result;
}
//****************************************************************************80

float r4_reverse_bytes ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_REVERSE_BYTES reverses the four bytes in an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, a value whose bytes are to be reversed.
//
//    Output, R4_REVERSE_BYTES, a value with bytes in reverse order;
//
{
  char c;
  union
  {
    float yfloat;
    char ychar[4];
  } y;

  y.yfloat = x;

  c = y.ychar[0];
  y.ychar[0] = y.ychar[3];
  y.ychar[3] = c;

  c = y.ychar[1];
  y.ychar[1] = y.ychar[2];
  y.ychar[2] = c;

  return ( y.yfloat );
}
//****************************************************************************80

int r4_round_i4 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ROUND_I4 rounds an R4, returning an I4.
//
//  Example:
//
//        X         Value
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
//    03 April 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_ROUND_I4, the rounded value.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - floor ( - x + 0.5 );
  }
  else
  {
    value =   floor (   x + 0.5 );
  }

  return value;
}
//****************************************************************************80

float r4_round2 ( int nplace, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ROUND2 rounds an R4 in base 2.
//
//  Discussion:
//
//    Assume that the input quantity X has the form
//
//      X = S * J * 2**L
//
//    where S is plus or minus 1, L is an integer, and J is a binary
//    mantissa which is either exactly zero, or greater than or equal
//    to 0.5 and less than 1.0.
//
//    Then on return, XROUND = R4_ROUND2 ( NPLACE, X ) will satisfy
//
//      XROUND = S * K * 2**L
//
//    where S and L are unchanged, and K is a binary mantissa which
//    agrees with J in the first NPLACE binary digits and is zero
//    thereafter.
//
//    If NPLACE is 0, XROUND will always be zero.
//
//    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.
//
//    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
//    or 0.75.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NPLACE, the number of binary digits to
//    preserve.  NPLACE should be 0 or positive.
//
//    Input, float X, the real number to be decomposed.
//
//    Output, float R4_ROUND2, the rounded value of X.
//
{
  int iplace;
  int l;
  int s;
  float xmant;
  float xtemp;
  float value;

  value = 0.0;
//
//  1: Handle the special case of 0.
//
  if ( x == 0.0 )
  {
    return value;
  }

  if ( nplace <= 0 )
  {
    return value;
  }
//
//  2: Determine the sign S.
//
  if ( 0.0 < x )
  {
    s = 1;
    xtemp = x;
  }
  else
  {
    s = -1;
    xtemp = -x;
  }
//
//  3: Force XTEMP to lie between 1 and 2, and compute the
//  logarithm L.
//
  l = 0;

  while ( 2.0 <= xtemp )
  {
    xtemp = xtemp / 2.0;
    l = l + 1;
  }

  while ( xtemp < 1.0 )
  {
    xtemp = xtemp * 2.0;
    l = l - 1;
  }
//
//  4: Strip out the digits of the mantissa as XMANT, and decrease L.
//
  xmant = 0.0;
  iplace = 0;

  for ( ; ; )
  {
    xmant = 2.0 * xmant;

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + 1.0;
      xtemp = xtemp - 1.0;
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      value = s * xmant * pow ( 2.0, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * 2.0;
  }

  return value;
}
//****************************************************************************80

float r4_roundb ( int base, int nplace, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ROUNDB rounds an R4 in a given base.
//
//  Discussion:
//
//    The code does not seem to do a good job of rounding when
//    the base is negative.
//
//    Assume that the input quantity X has the form
//
//      X = S * J * BASE^L
//
//    where S is plus or minus 1, L is an integer, and J is a
//    mantissa base BASE which is either exactly zero, or greater
//    than or equal to (1/BASE) and less than 1.0.
//
//    Then on return, XROUND will satisfy
//
//      XROUND = S * K * BASE^L
//
//    where S and L are unchanged, and K is a mantissa base BASE
//    which agrees with J in the first NPLACE digits and is zero
//    thereafter.
//
//    Note that because of rounding, for most bases, most numbers
//    with a fractional quantities cannot be stored exactly in the
//    computer, and hence will have trailing "bogus" digits.
//
//    If NPLACE is 0, XROUND will always be zero.
//
//    If NPLACE is 1, the mantissa of XROUND will be 0,
//    1/BASE, 2/BASE, ..., (BASE-1)/BASE.
//
//    If NPLACE is 2, the mantissa of XROUND will be 0,
//    BASE/BASE^2, (BASE+1)/BASE^2, ...,
//    BASE^2-2/BASE^2, BASE^2-1/BASE^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int BASE, the base of the arithmetic.
//    BASE must not be zero.  Theoretically, BASE may be negative.
//
//    Input, int NPLACE, the number of digits base BASE to
//    preserve.  NPLACE should be 0 or positive.
//
//    Input, float X, the number to be decomposed.
//
//    Output, float R4_ROUNDB, the rounded value of X.
//
{
  int iplace;
  int is;
  int js;
  int l;
  float r4_base;
  float value;
  float xmant;
  float xtemp;

  value = 0.0;
  r4_base = ( float ) base;
//
//  0: Error checks.
//
  if ( base == 0 )
  {
    cerr << "\n";
    cerr << "R4_ROUNDB - Fatal error!\n";
    cerr << "  The base BASE cannot be zero.\n";
    exit ( 1 );
  }
//
//  1: Handle the special case of 0.
//
  if ( x == 0.0 )
  {
    return value;
  }

  if ( nplace <= 0 )
  {
    return value;
  }
//
//  2: Determine the sign IS.
//
  if ( 0.0 < x )
  {
    is = 1;
    xtemp = x;
  }
  else
  {
    is = -1;
    xtemp = -x;
  }
//
//  3: Force XTEMP to lie between 1 and ABS(BASE), and compute the
//  logarithm L.
//
  l = 0;

  while ( r4_abs ( r4_base ) <= r4_abs ( xtemp ) )
  {
    xtemp = xtemp / r4_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }
    l = l + 1;
  }

  while ( r4_abs ( xtemp ) < 1.0 )
  {
    xtemp = xtemp * r4_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }

    l = l - 1;
  }
//
//  4: Now strip out the digits of the mantissa as XMANT, and
//  decrease L.
//
  xmant = 0.0;
  iplace = 0;
  js = is;

  for ( ; ; )
  {
    xmant = r4_base * xmant;

    if ( xmant < 0.0 )
    {
      js = -js;
      xmant = -xmant;
    }

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + ( int ) ( xtemp );
      xtemp = xtemp - ( int ) ( xtemp );
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      value = ( float ) js * xmant * pow ( r4_base, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * r4_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }
  }

  return value;
}
//****************************************************************************80

float r4_roundx ( int nplace, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ROUNDX rounds an R4 in base 10.
//
//  Discussion:
//
//    Assume that the input quantity X has the form
//
//      X = S * J * 10^L
//
//    where S is plus or minus 1, L is an integer, and J is a decimal
//    mantissa which is either exactly zero, or greater than or equal
//    to 0.1 and less than 1.0.
//
//    Then on return, XROUND will satisfy
//
//      XROUND = S * K * 10^L
//
//    where S and L are unchanged, and K is a decimal mantissa which
//    agrees with J in the first NPLACE decimal digits and is zero
//    thereafter.
//
//    Note that because of rounding, most decimal fraction quantities
//    cannot be stored exactly in the computer, and hence will have
//    trailing "bogus" digits.
//
//    If NPLACE is 0, XROUND will always be zero.
//
//    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
//    0.2, ..., or 0.9.
//
//    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
//    0.03, ..., 0.98, 0.99.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NPLACE, the number of decimal digits to
//    preserve.  NPLACE should be 0 or positive.
//
//    Input, float X, the number to be decomposed.
//
//    Output, float R4_ROUNDX, the rounded value of X.
//
{
  int iplace;
  int is;
  int l;
  float xmant;
  float xround;
  float xtemp;

  xround = 0.0;
//
//  1: Handle the special case of 0.
//
  if ( x == 0.0 )
  {
    return xround;
  }

  if ( nplace <= 0 )
  {
    return xround;
  }
//
//  2: Determine the sign IS.
//
  if ( 0.0 < x )
  {
    is = 1;
    xtemp = x;
  }
  else
  {
    is = -1;
    xtemp = -x;
  }
//
//  3: Force XTEMP to lie between 1 and 10, and compute the
//  logarithm L.
//
  l = 0;

  while ( 10.0 <= x )
  {
    xtemp = xtemp / 10.0;
    l = l + 1;
  }

  while ( xtemp < 1.0 )
  {
    xtemp = xtemp * 10.0;
    l = l - 1;
  }
//
//  4: Now strip out the digits of the mantissa as XMANT, and
//  decrease L.
//
  xmant = 0.0;
  iplace = 0;

  for ( ; ; )
  {
    xmant = 10.0 * xmant;

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + ( int ) xtemp;
      xtemp = xtemp - ( int ) xtemp;
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      xround = is * xmant * pow ( 10.0, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * 10.0;
  }

  return xround;
}
//****************************************************************************80

float r4_sign ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIGN returns the sign of an R4.
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
//    Input, float X, the number whose sign is desired.
//
//    Output, float R4_SIGN, the sign of X.
//
{
  float value;

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

float r4_sign3 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIGN3 returns the threew-way sign of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose sign is desired.
//
//    Output, float R4_SIGN3, the sign of X.
//
{
  float value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else if ( x == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

bool r4_sign_opposite ( float r1, float r2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIGN_OPPOSITE is TRUE if two R4's are not of the same sign.
//
//  Discussion:
//
//    This test could be coded numerically as
//
//      if ( r1 * r2 <= 0.0 ) ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R1, R2, the values to check.
//
//    Output, bool R4_SIGN_OPPOSITE, is TRUE if ( R1 <= 0 and 0 <= R2 )
//    or ( R2 <= 0 and 0 <= R1 ).
//
{
  bool value;

  value = ( r1 <= 0.0 && 0.0 <= r2 ) || ( r2 <= 0.0 && 0.0 <= r1 );

  return value;
}
//****************************************************************************80

bool r4_sign_opposite_strict ( float r1, float r2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIGN_OPPOSITE_STRICT is TRUE if two R4's are strictly of opposite sign.
//
//  Discussion:
//
//    This test could be coded numerically as
//
//      if ( r1 * r2 < 0.0 ) ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R1, R2, the values to check.
//
//    Output, bool R4_SIGN_OPPOSITE_STRICT, is TRUE if ( R1 < 0 and 0 < R2 )
//    or ( R2 < 0 and 0 < R1 ).
//
{
  bool value;

  value = ( r1 < 0.0 && 0.0 < r2 ) || ( r2 < 0.0 && 0.0 < r1 );

  return value;
}
//****************************************************************************80

void r4_swap ( float &x, float &y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SWAP switches two R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, float &X, &Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  float z;

  z = x;
  x = y;
  y = z;

  return;
}
//****************************************************************************80

void r4_swap3 ( float &x, float &y, float &z )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SWAP3 swaps three R4's.
//
//  Example:
//
//    Input:
//
//      X = 1, Y = 2, Z = 3
//
//    Output:
//
//      X = 2, Y = 3, Z = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, float &X, &Y, &Z, three values to be swapped.
//
{
  double w;

  w = x;
  x = y;
  y = z;
  z = w;

  return;
}
//****************************************************************************80

float r4_tiny ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_TINY returns a "tiny" R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, float R4_TINY, a "tiny" R4 value.
//
{
  const float value = 0.1175494350822E-37;

  return value;
}
//****************************************************************************80

void r4_to_dhms ( float r, int *d, int *h, int *m, int *s )

//****************************************************************************80
//
//  Purpose:
//
//    R4_TO_DHMS converts an R4 day value into days, hours, minutes, seconds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R, a real number representing a time period measured in days.
//
//    Output, int D, H, M, S, the equivalent number of days, hours,
//    minutes and seconds.
//
{
  int sign;

  if ( 0.0 <= r )
  {
    sign = 1;
  }
  else if ( r < 0.0 )
  {
    sign = -1;
    r = -r;
  }

  *d = ( int ) r;

  r = r - ( float ) *d;
  r = 24.0 * r;
  *h = ( int ) r;

  r = r - ( float ) *h;
  r = 60.0 * r;
  *m = ( int ) r;

  r = r - ( float ) *m;
  r = 60.0 * r;
  *s = ( int ) r;

  if ( sign == -1 )
  {
    *d = -(*d);
    *h = -(*h);
    *m = -(*m);
    *s = -(*s);
  }

  return;
}
//****************************************************************************80

int r4_to_i4 ( float x, float xmin, float xmax, int ixmin, int ixmax )

//****************************************************************************80
//
//  Purpose:
//
//    R4_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
//
//  Discussion:
//
//    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
//    IX := min ( IX, max ( IXMIN, IXMAX ) )
//    IX := max ( IX, min ( IXMIN, IXMAX ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the real number to be converted.
//
//    Input, float XMIN, XMAX, the real range.  XMAX and XMIN must not be
//    equal.  It is not necessary that XMIN be less than XMAX.
//
//    Input, int IXMIN, IXMAX, the allowed range of the output
//    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
//    It is not necessary that IXMIN be less than IXMAX.
//
//    Output, int R4_TO_I4, the value in the range [IXMIN,IXMAX] that
//    corresponds to X.
//
{
  int ix;
  float temp;

  if ( xmax == xmin )
  {
    cerr << "\n";
    cerr << "R4_TO_I4 - Fatal error!\n";
    cerr << "  XMAX = XMIN, making a zero divisor.\n";
    cerr << "  XMAX = " << xmax << "\n";
    cerr << "  XMIN = " << xmin << "\n";
    exit ( 1 );
  }

  temp =
      ( ( xmax - x        ) * ( float ) ixmin
      + (        x - xmin ) * ( float ) ixmax )
      / ( xmax     - xmin );

  if ( 0.0 <= temp )
  {
    temp = temp + 0.5;
  }
  else
  {
    temp = temp - 0.5;
  }

  ix = ( int ) temp;

  return ix;
}
//****************************************************************************80

float r4_to_r4_discrete ( float r, float rmin, float rmax, int nr )

//****************************************************************************80
//
//  Purpose:
//
//    R4_TO_R4_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
//
//  Discussion:
//
//    if ( R < RMIN ) then
//      RD = RMIN
//    else if ( RMAX < R ) then
//      RD = RMAX
//    else
//      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
//      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )
//
//    In the special case where NR = 1, when
//
//      XD = 0.5 * ( RMAX + RMIN )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R, the number to be converted.
//
//    Input, float RMAX, RMIN, the maximum and minimum
//    values for RD.
//
//    Input, int NR, the number of allowed values for XD.
//    NR should be at least 1.
//
//    Output, float RD, the corresponding discrete value.
//
{
  int f;
  float rd;
//
//  Check for errors.
//
  if ( nr < 1 )
  {
    cerr << "\n";
    cerr << "R4_TO_R4_DISCRETE - Fatal error!\n";
    cerr << "  NR = " << nr << "\n";
    cerr << "  but NR must be at least 1.\n";
    exit ( 1 );
  }

  if ( nr == 1 )
  {
    rd = 0.5 * ( rmin + rmax );
    return rd;
  }

  if ( rmax == rmin )
  {
    rd = rmax;
    return rd;
  }

  f = r4_nint ( ( float ) ( nr ) * ( rmax - r ) / ( rmax - rmin ) );
  f = i4_max ( f, 0 );
  f = i4_min ( f, nr );

  rd = ( ( float ) (      f ) * rmin
       + ( float ) ( nr - f ) * rmax )
       / ( float ) ( nr     );

  return rd;
}
//****************************************************************************80

float r4_uniform_01 ( int &seed )

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
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      r4_uniform_01 = seed / ( 2^31 - 1 )
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  float value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  value = ( float ) ( seed ) * 4.656612875E-10;

  return value;
}
//****************************************************************************80

float r4_uniform_ab ( float a, float b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_AB returns a scaled pseudorandom R4.
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
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_AB, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  float value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( float ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

void r4_unswap3 ( float &x, float &y, float &z )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNSWAP3 unswaps three R4's.
//
//  Example:
//
//    Input:
//
//      X = 2, Y = 3, Z = 1
//
//    Output:
//
//      X = 1, Y = 2, Z = 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, float &X, &Y, &Z, three values to be swapped.
//
{
  float w;

  w = z;
  z = y;
  y = x;
  x = w;

  return;
}
//****************************************************************************80

float r4_walsh_1d ( float x, int digit )

//****************************************************************************80
//
//  Purpose:
//
//    R4_WALSH_1D evaluates the Walsh function of a real scalar argument.
//
//  Discussion:
//
//    Consider the binary representation of X, and number the digits
//    in descending order, from leading to lowest, with the units digit
//    being numbered 0.
//
//    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the argument of the Walsh function.
//
//    Input, int DIGIT, the index of the Walsh function.
//
//    Output, float R4_WALSH_1D, the value of the Walsh function.
//
{
  int n;
  float value;
//
//  Hide the effect of the sign of X.
//
  x = r4_abs ( x );
//
//  If DIGIT is positive, divide by 2 DIGIT times.
//  If DIGIT is negative, multiply by 2 (-DIGIT) times.
//
  x = x / pow ( 2.0, digit );
//
//  Make it an integer.
//  Because it's positive, and we're using INT, we don't change the
//  units digit.
//
  n = ( int ) x;
//
//  Is the units digit odd or even?
//
  if ( ( n % 2 ) == 0 )
  {
    value = 0.0;
  }
  else
  {
    value = 1.0;
  }

  return value;
}
//****************************************************************************80

float *r42_cheby ( int n, float alo, float ahi )

//****************************************************************************80
//
//  Purpose:
//
//    R42_CHEBY sets up the Chebyshev abscissas in an R4 interval.
//
//  Discussion:
//
//    The routine sets up a vector of X values spaced between the values
//    XLO and XHI in a similar way to the spacing of the Chebyshev
//    points of the same order in the interval [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points to compute.
//
//    Input, float ALO, AHI, the range.
//
//    Output, float R42_CHEBY[N], the computed X values.
//
{
  float *a;
  float arg;
  int i;
  const float r4_pi = 3.141592653589793;

  a = new float[n];

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else if ( 1 < n )
  {
    for ( i = 0; i < n; i++ )
    {
      arg = ( float ) ( 2 * i + 1 ) * r4_pi / ( float ) ( 2 * n );

      a[i] = 0.5 * ( ( 1.0 + cos ( arg ) ) * alo
                   + ( 1.0 - cos ( arg ) ) * ahi );

    }
  }

  return a;
}
//****************************************************************************80

void r42_print ( float a[2], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R42_PRINT prints an R42.
//
//  Discussion:
//
//    An R42 is an R4VEC with two entries.
//
//    A format is used which suggests a coordinate pair:
//
//  Example:
//
//    Center : ( 1.23, 7.45 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[2], the coordinates of the vector.
//
//    Input, string TITLE, a title.
//
{
  cout << "  " << title << " : ";
  cout << ": ( " << setw(12) << a[0]
       << ", "   << setw(12) << a[1] << " )\n";

  return;
}
//****************************************************************************80

void r42_uniform ( float b, float c, int &seed, float r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R42_UNIFORM returns a random R42 value in a given range.
//
//  Discussion:
//
//    An R42 is an R4VEC with two entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float B, C, the minimum and maximum values.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, float R[2], the randomly chosen value.
//
{
  int i;

  for ( i = 0; i < 2; i++ )
  {
    r[i] = r4_uniform_ab ( b, c, seed );
  }

  return;
}
//****************************************************************************80

void r42poly2_print ( float a, float b, float c, float d, float e,
  float f )

//****************************************************************************80
//
//  Purpose:
//
//    R42POLY2_PRINT prints a second order polynomial in two variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, B, C, D, E, F, the coefficients.
//
{
  cout << "  " << setw(8) << a
       << " * x^2 + " << setw(8) << b
       << " * y^2 + " << setw(8) << c
       << " * xy  + " << "\n";
  cout << "  " << setw(8) << d
       << " * x   + " << setw(8) << e
       << " * y   + " << setw(8) << f << "\n";

  return;
}
//****************************************************************************80

int r42poly2_type ( float a, float b, float c, float d, float e, float f )

//****************************************************************************80
//
//  Purpose:
//
//    R42POLY2_TYPE analyzes a second order polynomial in two variables.
//
//  Discussion:
//
//    The polynomial has the form
//
//      A x^2 + B y^2 + C xy + Dx + Ey + F = 0
//
//    The possible types of the solution set are:
//
//     1: a hyperbola;
//        9x^2 -  4y^2       -36x - 24y -  36 = 0
//     2: a parabola;
//        4x^2 +  1y^2 - 4xy + 3x -  4y +   1 = 0;
//     3: an ellipse;
//        9x^2 + 16y^2       +36x - 32y -  92 = 0;
//     4: an imaginary ellipse (no real solutions);
//         x^2 +   y^2       - 6x - 10y + 115 = 0;
//     5: a pair of intersecting lines;
//                        xy + 3x -   y -   3 = 0
//     6: one point;
//         x^2 +  2y^2       - 2x + 16y +  33 = 0;
//     7: a pair of distinct parallel lines;
//                 y^2            -  6y +   8 = 0
//     8: a pair of imaginary parallel lines (no real solutions);
//                 y^2            -  6y +  10 = 0
//     9: a pair of coincident lines.
//                 y^2            -  2y +   1 = 0
//    10: a single line;
//                             2x -   y +   1 = 0;
//    11; all space;
//                                          0 = 0;
//    12; no solutions;
//                                          1 = 0;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    CRC Press, 30th Edition, 1996, pages 282-284.
//
//  Parameters:
//
//    Input, float A, B, C, D, E, F, the coefficients.
//
//    Output, int TYPE, indicates the type of the solution set.
//
{
  float delta;
  float j;
  float k;
  int type;
//
//  Handle the degenerate case.
//
  if ( a == 0.0 && b == 0.0 && c == 0.0 )
  {
    if ( d == 0.0 && e == 0.0 )
    {
      if ( f == 0.0 )
      {
        type = 11;
      }
      else
      {
        type = 12;
      }
    }
    else
    {
      type = 10;
    }
    return type;
  }

  delta =
      8.0 * a * b * f
    + 2.0 * c * e * d
    - 2.0 * a * e * e
    - 2.0 * b * d * d
    - 2.0 * f * c * c;

  j = 4.0 * a * b - c * c;

  if ( delta != 0.0 )
  {
    if ( j < 0.0 )
    {
      type = 1;
    }
    else if ( j == 0.0 )
    {
      type = 2;
    }
    else if ( 0.0 < j )
    {
      if ( r4_sign ( delta ) != r4_sign ( a + b ) )
      {
        type = 3;
      }
      else if ( r4_sign ( delta ) == r4_sign ( a + b ) )
      {
        type = 4;
      }
    }
  }
  else if ( delta == 0.0 )
  {
    if ( j < 0.0 )
    {
      type = 5;
    }
    else if ( 0.0 < j )
    {
      type = 6;
    }
    else if ( j == 0.0 )
    {
      k = 4.0 * ( a + b ) * f - d * d - e * e;

      if ( k < 0.0 )
      {
        type = 7;
      }
      else if ( 0.0 < k )
      {
        type = 8;
      }
      else if ( k == 0.0 )
      {
        type = 9;
      }
    }
  }

  return type;
}
//****************************************************************************80

void r42poly2_type_print ( int type )

//****************************************************************************80
//
//  Purpose:
//
//    R42POLY2_TYPE_PRINT prints the meaning of the output from R42POLY2_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TYPE, the type index returned by R42POLY2_TYPE.
//
{
  if ( type == 1 )
  {
    cout << "  The set of solutions forms a hyperbola.\n";
  }
  else if ( type == 2 )
  {
    cout << "  The set of solutions forms a parabola.\n";
  }
  else if ( type == 3 )
  {
    cout << "  The set of solutions forms an ellipse.\n";
  }
  else if ( type == 4 )
  {
    cout << "  The set of solutions forms an imaginary ellipse.\n";
    cout << "  (There are no real solutions).\n";
  }
  else if ( type == 5 )
  {
    cout << "  The set of solutions forms a pair of intersecting lines.\n";
  }
  else if ( type == 6 )
  {
    cout << "  The set of solutions is a single point.\n";
  }
  else if ( type == 7 )
  {
    cout << "  The set of solutions form a pair of distinct parallel lines.\n";
  }
  else if ( type == 8 )
  {
    cout << "  The set of solutions forms a pair of imaginary parallel lines.\n";
    cout << "  (There are no real solutions).\n";
  }
  else if ( type == 9 )
  {
    cout << "  The set of solutions forms a pair of coincident lines.\n";
  }
  else if ( type == 10 )
  {
    cout << "  The set of solutions forms a single line.\n";
  }
  else if ( type == 11 )
  {
    cout << "  The set of solutions is all space.\n";
  }
  else if ( type == 12 )
  {
    cout << "  The set of solutions is empty.\n";
  }
  else
  {
    cout << "  This type index is unknown.\n";
  }
  return;
}
//****************************************************************************80

float *r42vec_max ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R42VEC_MAX returns the maximum value in an R42VEC.
//
//  Discussion:
//
//    An R42VEC is an array of pairs of float precision real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[2*N], the array.
//
//    Output, float R42VEC_MAX[2]; the largest entries in each row.
//
{
# define DIM_NUM 2

  float *amax = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amax = new float[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amax[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[0+j*DIM_NUM] )
      {
        amax[i] = a[0+j*DIM_NUM];
      }
    }
  }
  return amax;
# undef DIM_NUM
}
//****************************************************************************80

float *r42vec_min ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R42VEC_MIN returns the minimum value in an R42VEC.
//
//  Discussion:
//
//    An R42VEC is an array of pairs of float precision real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[2*N], the array.
//
//    Output, float R42VEC_MIN[2]; the smallest entries in each row.
//
{
# define DIM_NUM 2

  float *amin = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amin = new float[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amin[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( a[0+j*DIM_NUM] < amin[i] )
      {
        amin[i] = a[0+j*DIM_NUM];
      }
    }
  }
  return amin;
# undef DIM_NUM
}
//****************************************************************************80

int r42vec_order_type ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R42VEC_ORDER_TYPE finds if an R42VEC is (non)strictly ascending/descending.
//
//  Discussion:
//
//    An R42VEC is a vector whose entries are R42's.
//    An R42 is a vector of type float precision with two entries.
//    An R42VEC may be stored as a 2 by N array.
//
//    The dictionary or lexicographic ordering is used.
//
//    (X1,Y1) < (X2,Y2)  <=>  X1 < X2 or ( X1 = X2 and Y1 < Y2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the array.
//
//    Input, float A[2*N], the array to be checked.
//
//    Output, int R42VEC_ORDER_TYPE, order indicator:
//    -1, no discernable order;
//    0, all entries are equal;
//    1, ascending order;
//    2, strictly ascending order;
//    3, descending order;
//    4, strictly descending order.
//
{
  int i;
  int order;
//
//  Search for the first value not equal to A(1,1).
//
  i = 0;

  for ( ; ; )
  {
    i = i + 1;

    if ( n <= i )
    {
      order = 0;
      return order;
    }

    if ( a[0+0*2] < a[0+i*2] || ( a[0+0*2] == a[0+i*2] && a[1+0*2] < a[1+i*2] ) )
    {
      if ( i == 2 )
      {
        order = 2;
      }
      else
      {
        order = 1;
      }
      break;
    }
    else if ( a[0+i*2] < a[0+0*2] || ( a[0+i*2] == a[0+0*2] && a[1+i*2] < a[1+0*2] ) )
    {
      if ( i == 2 )
      {
        order = 4;
      }
      else
      {
        order = 3;
      }
      break;
    }
  }
//
//  Now we have a "direction".  Examine subsequent entries.
//
  for ( ; ; )
  {
    i = i + 1;
    if ( n <= i )
    {
      break;
    }

    if ( order == 1 )
    {
      if ( a[0+i*2] < a[0+(i-1)*2] ||
        ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] < a[1+(i-1)*2] ) )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 2 )
    {
      if ( a[0+i*2] < a[0+(i-1)*2] ||
        ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] < a[1+(i-1)*2] ) )
      {
        order = -1;
        break;
      }
      else if ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] == a[1+(i-1)*2] )
      {
        order = 1;
      }
    }
    else if ( order == 3 )
    {
      if ( a[0+(i-1)*2] < a[0+i*2] ||
        ( a[0+(i-1)*2] == a[0+i*2] && a[1+(i-1)*2] < a[1+i*2] ) )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 4 )
    {
      if ( a[0+(i-1)*2] < a[0+i*2] ||
        ( a[0+(i-1)*2] == a[0+i*2] && a[1+(i-1)*2] < a[1+i*2] ) )
      {
        order = -1;
        break;
      }
      else if ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] == a[1+(i-1)*2] )
      {
        order = 3;
      }
    }
  }
  return order;
}
//****************************************************************************80

void r42vec_part_quick_a ( int n, float a[], int *l, int *r )

//****************************************************************************80
//
//  Purpose:
//
//    R42VEC_PART_QUICK_A reorders an R42VEC as part of a quick sort.
//
//  Discussion:
//
//    An R42VEC is a vector whose entries are R42's.
//    An R42 is a vector of type float precision with two entries.
//    An R42VEC may be stored as a 2 by N array.
//
//    The routine reorders the entries of A.  Using A(1:2,1) as a
//    key, all entries of A that are less than or equal to the key will
//    precede the key, which precedes all entries that are greater than the key.
//
//  Example:
//
//    Input:
//
//      N = 8
//
//      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
//
//    Output:
//
//      L = 2, R = 4
//
//      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
//             -----------          ----------------------------------
//             LEFT          KEY    RIGHT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, float A[N*2].  On input, the array to be checked.
//    On output, A has been reordered as described above.
//
//    Output, int *L, *R, the indices of A that define the three segments.
//    Let KEY = the input value of A(1:2,1).  Then
//    I <= L                 A(1:2,I) < KEY;
//         L < I < R         A(1:2,I) = KEY;
//                 R <= I    A(1:2,I) > KEY.
//
{
  int i;
  int j;
  float key[2];
  int ll;
  int m;
  int rr;
//
  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R42VEC_PART_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key[0] = a[2*0+0];
  key[1] = a[2*0+1];
  m = 1;
//
//  The elements of unknown size have indices between L+1 and R-1.
//
  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( r4vec_gt ( 2, a+2*ll, key ) )
    {
      rr = rr - 1;
      r4vec_swap ( 2, a+2*(rr-1), a+2*ll );
    }
    else if ( r4vec_eq ( 2, a+2*ll, key ) )
    {
      m = m + 1;
      r4vec_swap ( 2, a+2*(m-1), a+2*ll );
      ll = ll + 1;
    }
    else if ( r4vec_lt ( 2, a+2*ll, key ) )
    {
      ll = ll + 1;
    }

  }
//
//  Now shift small elements to the left, and KEY elements to center.
//
  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = a[2*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = key[j];
    }
  }

  *l = ll;
  *r = rr;

  return;
}
//****************************************************************************80

void r42vec_permute ( int n, int p[], int base, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R42VEC_PERMUTE permutes an R42VEC in place.
//
//  Discussion:
//
//    An R42VEC is a vector whose entries are R42's.
//    An R42 is a vector of type float precision with two entries.
//    An R42VEC may be stored as a 2 by N array.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.
//
//    Input/output, float A[2*N], the array to be permuted.
//
{
  float a_temp[2];
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "R42VEC_PERMUTE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
//  So temporarily add 1-BASE to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp[0] = a[0+(istart-1)*2];
      a_temp[1] = a[1+(istart-1)*2];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "R42VEC_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[0+(iput-1)*2] = a_temp[0];
          a[1+(iput-1)*2] = a_temp[1];
          break;
        }
        a[0+(iput-1)*2] = a[0+(iget-1)*2];
        a[1+(iput-1)*2] = a[1+(iget-1)*2];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the base of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 + base;
  }
  return;
}
//****************************************************************************80

void r42vec_print ( int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R42VEC_PRINT prints an R42VEC.
//
//  Discussion:
//
//    An R42VEC is a vector whose entries are R42's.
//    An R42 is a vector of type float precision with two entries.
//    An R42VEC may be stored as a 2 by N array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, float A[2*N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int j;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout << "  " << setw(8)  << j
         << ": " << setw(14) << a[0+j*2]
         << "  " << setw(14) << a[1+j*2] << "\n";
  }

  return;
}
//****************************************************************************80

int *r42vec_sort_heap_index_a ( int n, int base, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R42VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R42VEC.
//
//  Discussion:
//
//    An R42VEC is a vector whose entries are R42's.
//    An R42 is a vector of type float precision with two entries.
//    An R42VEC may be stored as a 2 by N array.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(*,indx(*))
//
//    or explicitly, by the call
//
//      r42vec_permute ( n, indx, base, a )
//
//    after which a(*,*) is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int BASE, the desired indexing for the sort index:
//    0 for 0-based indexing,
//    1 for 1-based indexing.
//
//    Input, float A[2*N], an array to be index-sorted.
//
//    Output, int R42VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
//    I-th element of the sorted array is A(0:1,R42VEC_SORT_HEAP_INDEX_A(I)).
//
{
  float aval[2];
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0] + base;
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval[0] = a[0+indxt*2];
      aval[1] = a[1+indxt*2];
    }
    else
    {
      indxt = indx[ir-1];
      aval[0] = a[0+indxt*2];
      aval[1] = a[1+indxt*2];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }
    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if (   a[0+indx[j-1]*2] <  a[0+indx[j]*2] ||
             ( a[0+indx[j-1]*2] == a[0+indx[j]*2] &&
               a[1+indx[j-1]*2] <  a[1+indx[j]*2] ) )
        {
          j = j + 1;
        }
      }

      if (   aval[0] <  a[0+indx[j-1]*2] ||
           ( aval[0] == a[0+indx[j-1]*2] &&
             aval[1] <  a[1+indx[j-1]*2] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
//
//  Take care of the base.
//
  for ( i = 0; i < n; i++ )
  {
    indx[i] = indx[i] + base;
  }

  return indx;
}
//****************************************************************************80

void r42vec_sort_quick_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R42VEC_SORT_QUICK_A ascending sorts an R42VEC using quick sort.
//
//  Discussion:
//
//    An R42VEC is a vector whose entries are R42's.
//    An R42 is a vector of type float precision with two entries.
//    An R42VEC may be stored as a 2 by N array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, float A[N*2].
//    On input, the array to be sorted.
//    On output, the array has been sorted.
//
{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R42VEC_SORT_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
//
//  Partition the segment.
//
    r42vec_part_quick_a ( n_segment, a+2*(base-1)+0, &l_segment, &r_segment );
//
//  If the left segment has more than one element, we need to partition it.
//
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        cerr << "\n";
        cerr<< "R42VEC_SORT_QUICK_A - Fatal error!\n";
        cerr << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
//
//  The left segment and the middle segment are sorted.
//  Must the right segment be partitioned?
//
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
//
//  Otherwise, we back up a level if there is an earlier one.
//
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }
      }
    }
  }
  return;
# undef LEVEL_MAX
}
//****************************************************************************80

float *r43vec_max ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R43VEC_MAX returns the maximum value in an R43VEC.
//
//  Discussion:
//
//    An R43VEC is an array of triples of float precision real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[3*N], the array.
//
//    Output, float R43VEC_MAX[3]; the largest entries in each row.
//
{
# define DIM_NUM 3

  float *amax = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amax = new float[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amax[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[i+j*DIM_NUM] )
      {
        amax[i] = a[i+j*DIM_NUM];
      }
    }
  }
  return amax;
# undef DIM_NUM
}
//****************************************************************************80

float *r43vec_min ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R43VEC_MIN returns the minimum value in an R43VEC.
//
//  Discussion:
//
//    An R43VEC is an array of triples of float precision real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[3*N], the array.
//
//    Output, float R43VEC_MIN[3]; the smallest entries in each row.
//
{
# define DIM_NUM 3

  float *amin = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amin = new float[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amin[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( a[i+j*DIM_NUM] < amin[i] )
      {
        amin[i] = a[i+j*DIM_NUM];
      }
    }
  }
  return amin;
# undef DIM_NUM
}
//****************************************************************************80

void r43vec_part_quick_a ( int n, float a[], int *l, int *r )

//****************************************************************************80
//
//  Purpose:
//
//    R43VEC_PART_QUICK_A reorders an R43VEC as part of a quick sort.
//
//  Discussion:
//
//    An R43VEC is a vector whose entries are R43's.
//    An R43 is a vector of type float precision with three entries.
//    An R43VEC may be stored as a 3 by N array.
//
//    The routine reorders the entries of A.  Using A(1:3,1) as a
//    key, all entries of A that are less than or equal to the key will
//    precede the key, which precedes all entries that are greater than the key.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, float A[3*N].  On input, the array to be checked.
//    On output, A has been reordered as described above.
//
//    Output, int *L, *R, the indices of A that define the three segments.
//    Let KEY = the input value of A(1:3,1).  Then
//    I <= L                 A(1:3,I) < KEY;
//         L < I < R         A(1:3,I) = KEY;
//                 R <= I    A(1:3,I) > KEY.
//
{
  int i;
  int j;
  float key[3];
  int ll;
  int m;
  int rr;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R43VEC_PART_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key[0] = a[3*0+0];
  key[1] = a[3*0+1];
  key[2] = a[3*0+2];
  m = 1;
//
//  The elements of unknown size have indices between L+1 and R-1.
//
  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( r4vec_gt ( 3, a+3*ll, key ) )
    {
      rr = rr - 1;
      r4vec_swap ( 3, a+3*(rr-1), a+3*ll );
    }
    else if ( r4vec_eq ( 3, a+3*ll, key ) )
    {
      m = m + 1;
      r4vec_swap ( 3, a+3*(m-1), a+3*ll );
      ll = ll + 1;
    }
    else if ( r4vec_lt ( 3, a+3*ll, key ) )
    {
      ll = ll + 1;
    }
  }
//
//  Now shift small elements to the left, and KEY elements to center.
//
  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = a[3*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = key[j];
    }
  }

  *l = ll;
  *r = rr;

  return;
}
//****************************************************************************80

void r43vec_sort_quick_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R43VEC_SORT_QUICK_A ascending sorts an R43VEC using quick sort.
//
//  Discussion:
//
//    An R43VEC is a vector whose entries are R43's.
//    An R43 is a vector of type float precision with three entries.
//    An R43VEC may be stored as a 3 by N array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, float A[N*3].
//    On input, the array to be sorted.
//    On output, the array has been sorted.
//
{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R43VEC_SORT_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
//
//  Partition the segment.
//
    r43vec_part_quick_a ( n_segment, a+3*(base-1)+0, &l_segment, &r_segment );
//
//  If the left segment has more than one element, we need to partition it.
//
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        cerr << "\n";
        cerr << "R43VEC_SORT_QUICK_A - Fatal error!\n";
        cerr << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }
      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
//
//  The left segment and the middle segment are sorted.
//  Must the right segment be partitioned?
//
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
//
//  Otherwise, we back up a level if there is an earlier one.
//
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }
      }
    }
  }
  return;
# undef LEVEL_MAX
}
//****************************************************************************80

float *r4block_expand_linear ( int l, int m, int n, float x[], int lfat,
  int mfat, int nfat )

//****************************************************************************80
//
//  Purpose:
//
//    R4BLOCK_EXPAND_LINEAR linearly interpolates new data into a 3D block.
//
//  Discussion:
//
//    In this routine, the expansion is specified by giving the number
//    of intermediate values to generate between each pair of original
//    data rows and columns.
//
//    The interpolation is not actually linear.  It uses the functions
//
//      1, x, y, z, xy, xz, yz, xyz.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L, M, N, the dimensions of the input data.
//
//    Input, float X[L*M*N], the original data.
//
//    Input, int LFAT, MFAT, NFAT, the number of data values to interpolate
//    original data values in the first, second and third dimensions.
//
//    Output, float XFAT[L2*M2*N2], the fattened data, where
//    L2 = (L-1)*(LFAT+1)+1,
//    M2 = (M-1)*(MFAT+1)+1,
//    N2 = (N-1)*(NFAT+1)+1.
//
{
  int i;
  int ihi;
  int ii;
  int iii;
  int ip1;
  int j;
  int jhi;
  int jj;
  int jjj;
  int jp1;
  int k;
  int khi;
  int kk;
  int kkk;
  int kp1;
  int l2;
  int m2;
  int n2;
  float r;
  float s;
  float t;
  float x000;
  float x001;
  float x010;
  float x011;
  float x100;
  float x101;
  float x110;
  float x111;
  float *xfat;

  l2 = ( l - 1 ) * ( lfat + 1 ) + 1;
  m2 = ( m - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( n - 1 ) * ( nfat + 1 ) + 1;

  xfat = new float[l2*m2*n2];

  for ( i = 1; i <= l; i++ )
  {
    if ( i < l )
    {
      ihi = lfat;
    }
    else
    {
      ihi = 0;
    }

    for ( j = 1; j <= m; j++ )
    {
      if ( j < m )
      {
        jhi = mfat;
      }
      else
      {
        jhi = 0;
      }

      for ( k = 1; k <= n; k++ )
      {
        if ( k < n )
        {
          khi = nfat;
        }
        else
        {
          khi = 0;
        }

        if ( i < l )
        {
          ip1 = i + 1;
        }
        else
        {
          ip1 = i;
        }

        if ( j < m )
        {
          jp1 = j + 1;
        }
        else
        {
          jp1 = j;
        }

        if ( k < n )
        {
          kp1 = k + 1;
        }
        else
        {
          kp1 = k;
        }

        x000 = x[i-1+(j-1)*l+(k-1)*l*m];
        x001 = x[i-1+(j-1)*l+(kp1-1)*l*m];
        x100 = x[ip1-1+(j-1)*l+(k-1)*l*m];
        x101 = x[ip1-1+(j-1)*l+(kp1-1)*l*m];
        x010 = x[i-1+(jp1-1)*l+(k-1)*l*m];
        x011 = x[i-1+(jp1-1)*l+(kp1-1)*l*m];
        x110 = x[ip1-1+(jp1-1)*l+(k-1)*l*m];
        x111 = x[ip1-1+(jp1-1)*l+(kp1-1)*l*m];

        for ( ii = 0; ii <= ihi; ii++ )
        {
          r = ( float ) ( ii ) / ( float ) ( ihi + 1 );

          for ( jj = 0; jj <= jhi; jj++ )
          {
            s = ( float ) ( jj ) / ( float ) ( jhi + 1 );

            for ( kk = 0; kk <= khi; kk++ )
            {
              t = ( float ) ( kk ) / ( float ) ( khi + 1 );

              iii = 1 + ( i - 1 ) * ( lfat + 1 ) + ii;
              jjj = 1 + ( j - 1 ) * ( mfat + 1 ) + jj;
              kkk = 1 + ( k - 1 ) * ( nfat + 1 ) + kk;

              xfat[iii-1+(jjj-1)*l2+(kkk-1)*l2*m2] =
                  x000 * ( 1.0 - r ) * ( 1.0 - s ) * ( 1.0 - t )
                + x001 * ( 1.0 - r ) * ( 1.0 - s ) * (       t )
                + x010 * ( 1.0 - r ) * (       s ) * ( 1.0 - t )
                + x011 * ( 1.0 - r ) * (       s ) * (       t )
                + x100 * (       r ) * ( 1.0 - s ) * ( 1.0 - t )
                + x101 * (       r ) * ( 1.0 - s ) * (       t )
                + x110 * (       r ) * (       s ) * ( 1.0 - t )
                + x111 * (       r ) * (       s ) * (       t );
            }
          }
        }
      }
    }
  }

  return xfat;
}
//****************************************************************************80

void r4block_print ( int l, int m, int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4BLOCK_PRINT prints an R4BLOCK block (a 3D matrix).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L, M, N, the dimensions of the block.
//
//    Input, float A[L*M*N], the matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int j;
  int jhi;
  int jlo;
  int k;

  cout << "\n";
  cout << title << "\n";

  for ( k = 1; k <= n; k++ )
  {
    cout << "\n";
    cout << "  K = " << k << "\n";
    cout << "\n";
    for ( jlo = 1; jlo <= m; jlo = jlo + 5 )
    {
      jhi = i4_min ( jlo + 4, m );
      cout << "\n";
      cout << "      ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(7) << j << "       ";
      }
      cout << "\n";
      cout << "\n";
      for ( i = 1; i <= l; i++ )
      {
        cout << setw(5) << i << ":";
        for ( j = jlo; j <= jhi; j++ )
        {
          cout << "  " << setw(12) << a[i-1+(j-1)*l+(k-1)*l*m];
        }
        cout << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

float *r4block_zero_new ( int l, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4BLOCK_ZERO_NEW returns a new zeroed R4BLOCK.
//
//  Discussion:
//
//    An R4BLOCK is a triple dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L, M, N, the number of rows and columns.
//
//    Output, float R4BLOCK_ZERO_NEW[L*M*N], the new zeroed matrix.
//
{
  float *a;
  int i;
  int j;
  int k;

  a = new float[l*m*n];

  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        a[i+j*l+k*l*m] = 0.0;
      }
    }
  }
  return a;
}
//****************************************************************************80

int r4col_compare ( int m, int n, float a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_COMPARE compares two columns in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1.  2.  3.  4.
//        5.  6.  7.  8.
//        9. 10. 11. 12. )
//
//    Output:
//
//      R4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the M by N array.
//
//    Input, int I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, int R4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
  int k;
  int value;
//
//  Check.
//
  if ( i < 1 || n < i )
  {
    cerr << "\n";
    cerr << "R4COL_COMPARE - Fatal error!\n";
    cerr << "  Column index I is out of bounds.\n";
    cerr << "  I = " << i << "\n";
    exit ( 1 );
  }

  if ( j < 1 || n < j )
  {
    cerr << "\n";
    cerr << "R4COL_COMPARE - Fatal error!\n";
    cerr << "  Column index J is out of bounds.\n";
    cerr << "  J = " << j << "\n";
    exit ( 1 );
  }

  value = 0;

  if ( i == j )
  {
    return value;
  }

  k = 0;

  while ( k < m )
  {
    if ( a[k+(i-1)*m] < a[k+(j-1)*m] )
    {
      value = -1;
      return value;
    }
    else if ( a[k+(j-1)*m] < a[k+(i-1)*m] )
    {
      value = +1;
      return value;
    }
    k = k + 1;
  }

  return value;
}
//****************************************************************************80

float *r4col_duplicates ( int m, int n, int n_unique, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_DUPLICATES generates an R4COL with some duplicate columns.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    This routine generates a random R4COL with a specified number of
//    duplicate columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in each column of A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int N_UNIQUE, the number of unique columns in A.
//    1 <= N_UNIQUE <= N.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, float R4COL_DUPLICATES[M*N], the array.
//
{
  float *a;
  int i;
  int j1;
  int j2;
  float temp;

  if ( n_unique < 1 || n < n_unique )
  {
    cout << "\n";
    cout << "R4COL_DUPLICATES - Fatal error!\n";
    cout << "  1 <= N_UNIQUE <= N is required.\n";
    exit ( 1 );
  }

  a = new float[m*n];

  r4mat_uniform_01 ( m, n_unique, seed, a );
//
//  Randomly copy unique columns.
//
  for ( j1 = n_unique; j1 < n; j1++ )
  {
    j2 = i4_uniform ( 0, n_unique - 1, seed );
    for ( i = 0; i < m; i++ )
    {
      a[i+j1*m] = a[i+j2*m];
    }
  }
//
//  Permute the columns.
//
  for ( j1 = 0; j1 < n; j1 ++ )
  {
    j2 = i4_uniform ( j1, n - 1, seed );
    for ( i = 0; i < m; i++ )
    {
      temp      = a[i+j1*m];
      a[i+j1*m] = a[i+j2*m];
      a[i+j2*m] = temp;
    }
  }
  return a;
}
//****************************************************************************80

int r4col_find ( int m, int n, float a[], float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_FIND seeks a column value in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Example:
//
//    Input:
//
//      M = 3,
//      N = 4,
//
//      A = (
//        1.  2.  3.  4.
//        5.  6.  7.  8.
//        9. 10. 11. 12. )
//
//      x = ( 3.,
//            7.,
//           11. )
//
//    Output:
//
//      R4COL_FIND = 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], a table of numbers, regarded as
//    N columns of vectors of length M.
//
//    Input, float X[M], a vector to be matched with a column of A.
//
//    Output, int R4COL_FIND, the (one-based) index of the first column of A
//    which exactly matches every entry of X, or -1 if no match
//    could be found.
//
{
  int col;
  int i;
  int j;

  col = -1;

  for ( j = 1; j <= n; j++ )
  {
    col = j;

    for ( i = 1; i <= m; i++ )
    {
      if ( x[i-1] != a[i-1+(j-1)*m] )
      {
        col = -1;
        break;
      }
    }
    if ( col != -1 )
    {
      return col;
    }
  }
  return col;
}
//****************************************************************************80

int *r4col_first_index ( int m, int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_FIRST_INDEX indexes the first occurrence of values in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    For element A(1:M,J) of the matrix, FIRST_INDEX(J) is the index in A of
//    the first column whose entries are equal to A(1:M,J).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of A.
//    The length of an "element" of A, and the number of "elements".
//
//    Input, float A[M*N], the array.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int R4COL_FIRST_INDEX[N], the first occurrence index.
//
{
  float diff;
  int *first_index;
  int i;
  int j1;
  int j2;

  first_index = new int[n];

  for ( j1 = 0; j1 < n; j1++ )
  {
    first_index[j1] = -1;
  }
  for ( j1 = 0; j1 < n; j1++ )
  {
    if ( first_index[j1] == -1 )
    {
      first_index[j1] = j1;

      for ( j2 = j1 + 1; j2 < n; j2++ )
      {
        diff = 0.0;
        for ( i = 0; i < m; i++ )
        {
          diff = r4_max ( diff, r4_abs ( a[i+j1*m] - a[i+j2*m] ) );
        }
        if ( diff <= tol )
        {
          first_index[j2] = j1;
        }
      }
    }
  }
  return first_index;
}
//****************************************************************************80

int r4col_insert ( int n_max, int m, int n, float a[], float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_INSERT inserts a column into an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Example:
//
//    Input:
//
//      N_MAX = 10,
//      M = 3,
//      N = 4,
//
//      A = (
//        1.  2.  3.  4.
//        5.  6.  7.  8.
//        9. 10. 11. 12. )
//
//      X = ( 3., 4., 18. )
//
//    Output:
//
//      N = 5,
//
//      A = (
//        1.  2.  3.  3.  4.
//        5.  6.  4.  7.  8.
//        9. 10. 18. 11. 12. )
//
//      R4COL_INSERT = 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_MAX, the maximum number of columns in A.
//
//    Input, int M, the number of rows.
//
//    Input/output, int N, the number of columns.
//    If the new column is inserted into the table, then the output
//    value of N will be increased by 1.
//
//    Input/output, float A[M*N_MAX], a table of numbers, regarded
//    as an array of columns.  The columns must have been sorted
//    lexicographically.
//
//    Input, float X[M], a vector of data which will be inserted
//    into the table if it does not already occur.
//
//    Output, int R4COL_INSERT.
//    I, X was inserted into column I.
//    -I, column I was already equal to X.
//    0, N = N_MAX.
//
{
  int col;
  int high;
  int i;
  int isgn;
  int j;
  int low;
  int mid;
//
//  Refuse to work if N_MAX <= N.
//
  if ( n_max <= n )
  {
    col = 0;
    return col;
  }
//
//  Stick X temporarily in column N+1, just so it's easy to use R4COL_COMPARE.
//
  for ( i = 0; i < m; i++ )
  {
    a[i+n*m] = x[i];
  }
//
//  Do a binary search.
//
  low = 1;
  high = n;

  for ( ; ; )
  {
    if ( high < low )
    {
      col = low;
      break;
    }

    mid = ( low + high ) / 2;

    isgn = r4col_compare ( m, n+1, a, mid, n+1 );

    if ( isgn == 0 )
    {
      col = -mid;
      return col;
    }
    else if ( isgn == -1 )
    {
      low = mid + 1;
    }
    else if ( isgn == +1 )
    {
      high = mid - 1;
    }
  }
//
//  Shift part of the table up to make room.
//
  for ( j = n-1; col-1 <= j; j-- )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+(j+1)*m] = a[i+j*m];
    }
  }
//
//  Insert the new column.
//
  for ( i = 0; i < m; i++ )
  {
    a[i+(col-1)*m] = x[i];
  }

  n = n + 1;

  return col;
}
//****************************************************************************80

float *r4col_max ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_MAX returns the column maximums of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, float R4COL_MAX[N], the maximums of the columns.
//
{
  float *amax;
  int i;
  int j;

  amax = new float[n];

  for ( j = 0; j < n; j++ )
  {
    amax[j] = a[0+j*m];
    for ( i = 0; i < m; i++ )
    {
      amax[j] = r4_max ( amax[j], a[i+j*m] );
    }
  }

  return amax;
}
//****************************************************************************80

int *r4col_max_index ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_MAX_INDEX returns the indices of column maximums in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, int R4COL_MAX_INDEX[N]; entry I is the row of A in which
//    the maximum for column I occurs.
//
{
  float amax;
  int i;
  int *imax;
  int j;

  imax = new int[n];

  for ( j = 0; j < n; j++ )
  {
    imax[j] = 1;
    amax = a[0+j*m];

    for ( i = 1; i < m; i++ )
    {
      if ( amax < a[i+j*m] )
      {
        imax[j] = i+1;
        amax = a[i+j*m];
      }
    }
  }

  return imax;
}
//****************************************************************************80

void r4col_max_one ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_MAX_ONE rescales an R4COL so each column maximum is 1.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, float A[M*N], the array to be rescaled.
//
{
  int i;
  int i_big;
  int j;
  float temp;

  for ( j = 0; j < n; j++ )
  {
    i_big = 0;
    for ( i = 1; i < m; i++ )
    {
      if ( r4_abs ( a[i_big+j*m] ) < r4_abs ( a[i+j*m] ) )
      {
        i_big = i;
      }
    }
    temp = a[i_big+j*m];

    if ( temp != 0.0 )
    {
      for ( i = 0; i < m; i++ )
      {
        a[i+j*m] = a[i+j*m] / temp;
      }
    }
  }
  return;
}
//****************************************************************************80

float *r4col_mean ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_MEAN returns the column means of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    R4COL_MEAN =
//      1.5  4.0  5.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, float R4COL_MEAN[N], the means, or averages, of the columns.
//
{
  int i;
  int j;
  float *mean;

  mean = new float[n];

  for ( j = 0; j < n; j++ )
  {
    mean[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      mean[j] = mean[j] + a[i+j*m];
    }
    mean[j] = mean[j] / ( float ) ( m );
  }

  return mean;
}
//****************************************************************************80

float *r4col_min ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_MIN returns the column minimums of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, float R4COL_MIN[N], the minimums of the columns.
//
{
  float *amin;
  int i;
  int j;

  amin = new float[n];

  for ( j = 0; j < n; j++ )
  {
    amin[j] = a[0+j*m];
    for ( i = 0; i < m; i++ )
    {
      amin[j] = r4_min ( amin[j], a[i+j*m] );
    }
  }

  return amin;
}
//****************************************************************************80

int *r4col_min_index ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_MIN_INDEX returns the indices of column minimums in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, int R4COL_MIN_INDEX[N]; entry I is the row of A in which
//    the minimum for column I occurs.
//
{
  float amin;
  int i;
  int *imin;
  int j;

  imin = new int[n];

  for ( j = 0; j < n; j++ )
  {
    imin[j] = 1;
    amin = a[0+j*m];

    for ( i = 1; i < m; i++ )
    {
      if ( a[i+j*m] < amin )
      {
        imin[j] = i+1;
        amin = a[i+j*m];
      }
    }
  }

  return imin;
}
//****************************************************************************80

void r4col_part_quick_a ( int m, int n, float a[], int *l, int *r )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_PART_QUICK_A reorders the columns of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The routine reorders the columns of A.  Using A(1:M,1) as a
//    key, all entries of A that are less than or equal to the key will
//    precede the key, which precedes all entries that are greater than the key.
//
//  Example:
//
//    Input:
//
//      M = 2, N = 8
//      A = ( 2  8  6  0 10 10  0  5
//            4  8  2  2  6  0  6  8 )
//
//    Output:
//
//      L = 2, R = 4
//
//      A = (  0  0  2  8  6 10 10  4
//             2  6  4  8  2  6  0  8 )
//             ----     -------------
//             LEFT KEY     RIGHT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the row dimension of A, and the length of a column.
//
//    Input, int N, the column dimension of A.
//
//    Input/output, float A[M*N].  On input, the array to be checked.
//    On output, A has been reordered as described above.
//
//    Output, int *L, *R, the indices of A that define the three segments.
//    Let KEY = the input value of A(1:M,1).  Then
//    I <= L                 A(1:M,I) < KEY;
//         L < I < R         A(1:M,I) = KEY;
//                 R <= I    KEY < A(1:M,I).
//
{
  int i;
  int j;
  int k;
  float *key;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R4COL_PART_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key = new float[m];

  for ( i = 0; i < m; i++ )
  {
    key[i] = a[i+0*m];
  }
  k = 1;
//
//  The elements of unknown size have indices between L+1 and R-1.
//
  *l = 1;
  *r = n + 1;

  for ( j = 1; j < n; j++ )
  {
    if ( r4vec_gt ( m, a+(*l)*m, key ) )
    {
      *r = *r - 1;
      r4vec_swap ( m, a+(*r-1)*m, a+(*l)*m );
    }
    else if ( r4vec_eq ( m, a+(*l)*m, key ) )
    {
      k = k + 1;
      r4vec_swap ( m, a+(k-1)*m, a+(*l)*m );
      *l = *l + 1;
    }
    else if ( r4vec_lt ( m, a+(*l)*m, key ) )
    {
      *l = *l + 1;
    }
  }
//
//  Shift small elements to the left.
//
  for ( j = 0; j < *l - k; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+(j+k)*m];
    }
  }
//
//  Shift KEY elements to center.
//
  for ( j = *l-k; j < *l; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = key[i];
    }
  }
//
//  Update L.
//
  *l = *l - k;

  delete [] key;

  return;
}
//****************************************************************************80

void r4col_permute ( int m, int n, int p[], int base, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_PERMUTE permutes an R4COL in place.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      M = 2
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//      BASE = 1
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the length of objects.
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for a
//    1-based permutation.
//
//    Input/output, float A[M*N], the array to be permuted.
//
{
  float *a_temp;
  int i;
  int iget;
  int iput;
  int istart;
  int j;

  if ( !perm_check ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "R4COL_PERMUTE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
//  So temporarily add 1-BASE to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }

  a_temp = new float[m];
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      for ( i = 0; i < m; i++ )
      {
        a_temp[i] = a[i+(istart-1)*m];
      }
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "R4COL_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          for ( i = 0; i < m; i++ )
          {
            a[i+(iput-1)*m] = a_temp[i];
          }
          break;
        }
        for ( i = 0; i < m; i++ )
        {
          a[i+(iput-1)*m] = a[i+(iget-1)*m];
        }
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( j = 0; j < n; j++ )
  {
    p[j] = - p[j];
  }
//
//  Restore the base of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 +  base;
  }

  delete [] a_temp;

  return;
}
//****************************************************************************80

void r4col_sort_heap_a ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORT_HEAP_A ascending heapsorts an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    In lexicographic order, the statement "X < Y", applied to two real
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, the first time they differ, X is smaller.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, float A[M*N].
//    On input, the array of N columns of M-vectors.
//    On output, the columns of A have been sorted in lexicographic order.
//
{
  int i;
  int indx;
  int isgn;
  int j;

  if ( m <= 0 )
  {
    return;
  }

  if ( n <= 1 )
  {
    return;
  }
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      r4col_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = r4col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

int *r4col_sort_heap_index_a ( int m, int n, int base, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2) is negative.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(*,INDX(*)) is sorted,
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in each column of A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int BASE, the desired indexing for the sort index:
//    0 for 0-based indexing,
//    1 for 1-based indexing.
//
//    Input, float A[M*N], the array.
//
//    Output, int R4COL_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th column of the sorted array is A(*,INDX(I)).
//
{
  float *column;
  int i;
  int *indx;
  int indxt;
  int ir;
  int isgn;
  int j;
  int k;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0] + base;
    return indx;
  }

  column = new float[m];

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      for ( k = 0; k < m; k++ )
      {
        column[k] = a[k+indxt*m];
      }
    }
    else
    {
      indxt = indx[ir-1];
      for ( k = 0; k < m; k++ )
      {
        column[k] = a[k+indxt*m];
      }
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        isgn = r4vec_compare ( m, a+indx[j-1]*m, a+indx[j]*m );

        if ( isgn < 0 )
        {
          j = j + 1;
        }
      }

      isgn = r4vec_compare ( m, column, a+indx[j-1]*m );

      if ( isgn < 0 )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
  delete [] column;
//
//  Take care of the base.
//
  for ( i = 0; i < n; i++ )
  {
    indx[i] = indx[i] + base;
  }

  return indx;
}
//****************************************************************************80

void r4col_sort_quick_a ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORT_QUICK_A ascending quick sorts an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the row order of A, and the length of a column.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, float A[M*N].
//    On input, the array to be sorted.
//    On output, the array has been sorted.
//
{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( m <= 0 )
  {
    return;
  }

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R4COL_SORT_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  for ( ; ; )
  {
//
//  Partition the segment.
//
    r4col_part_quick_a ( m, n_segment, a+(base-1)*m, &l_segment, &r_segment );
//
//  If the left segment has more than one element, we need to partition it.
//
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        cerr << "\n";
        cerr << "R4COL_SORT_QUICK_A - Fatal error!\n";
        cerr << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
//
//  The left segment and the middle segment are sorted.
//  Must the right segment be partitioned?
//
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
//
//  Otherwise, we back up a level if there is an earlier one.
//
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          return;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }
      }
    }
  }
  return;
# undef LEVEL_MAX
}
//****************************************************************************80

void r4col_sorted_tol_undex ( int m, int n, float a[], int unique_num,
  float tol, int undx[], int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORTED_TOL_UNDEX: index tolerably unique entries of a sorted R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The goal of this routine is to determine a vector UNDX,
//    which points, to the unique elements of A, in sorted order,
//    and a vector XDNU, which identifies, for each entry of A, the index of
//    the unique sorted element of A.
//
//    This is all done with index vectors, so that the elements of
//    A are never moved.
//
//    Assuming A is already sorted, we examine the entries of A in order,
//    noting the unique entries, creating the entries of XDNU and
//    UNDX as we go.
//
//    Once this process has been completed, the vector A could be
//    replaced by a compressed vector XU, containing the unique entries
//    of X in sorted order, using the formula
//
//      XU(*) = A(UNDX(*)).
//
//    We could then, if we wished, reconstruct the entire vector A, or
//    any element of it, by index, as follows:
//
//      A(I) = XU(XDNU(I)).
//
//    We could then replace A by the combination of XU and XDNU.
//
//    Later, when we need the I-th entry of A, we can locate it as
//    the XDNU(I)-th entry of XU.
//
//    Here is an example of a vector A, the unique sort and inverse unique
//    sort vectors and the compressed unique sorted vector.
//
//      I      A      XU  Undx  Xdnu
//    ----+------+------+-----+-----+
//      0 | 11.0 |  11.0    0     0
//      1 | 11.0 |  22.0    4     0
//      2 | 11.0 |  33.0    7     0
//      3 | 11.0 |  55.0    8     0
//      4 | 22.0 |                1
//      5 | 22.0 |                1
//      6 | 22.0 |                1
//      7 | 33.0 |                2
//      8 | 55.0 |                3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the data values.
//
//    Input, int N, the number of data values,
//
//    Input, float A[M*N], the data values.
//
//    Input, int UNIQUE_NUM, the number of unique values in A.
//    This value is only required for languages in which the size of
//    UNDX must be known in advance.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
//
//    Output, int XDNU[N], the XDNU vector.
//
{
  float diff;
  int i;
  int i2;
  int i3;
  int j;
  int k;
  bool unique;
//
//  Consider entry I = 0.
//  It is unique, so set the number of unique items to K.
//  Set the K-th unique item to I.
//  Set the representative of item I to the K-th unique item.
//
  i = 0;
  k = 0;
  undx[k] = i;
  xdnu[i] = k;
//
//  Consider entry I.
//
//  If it is unique, increase the unique count K, set the
//  K-th unique item to I, and set the representative of I to K.
//
//  If it is not unique, set the representative of item I to a
//  previously determined unique item that is close to it.
//
  for ( i = 1; i < n; i++ )
  {
    unique = true;

    for ( j = 0; j <= k; j++ )
    {
      i2 = undx[j];
      diff = 0.0;
      for ( i3 = 0; i3 < m; i3++ )
      {
        diff = r4_max ( diff, r4_abs ( a[i3+i*m] - a[i3+i2*m] ) );
      }
      if ( diff <= tol )
      {
        unique = false;
        xdnu[i] = j;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = i;
      xdnu[i] = k;
    }
  }
  return;
}
//****************************************************************************80

int r4col_sorted_tol_unique ( int m, int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORTED_TOL_UNIQUE keeps tolerably unique elements in a sorted R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The columns of the array can be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, float A(M,N).
//    On input, the sorted array of N columns of M-vectors.
//    On output, a sorted array of columns of M-vectors.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int R4COL_SORTED_TOL_UNIQUE, the number of unique columns.
//
{
  float diff;
  int i;
  int j;
  int k;
  bool unique;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    unique = true;
    for ( j = 0; j < unique_num; j++ )
    {
      diff = 0.0;
      for ( k = 0; k < m; k++ )
      {
        diff = r4_max ( diff, r4_abs ( a[k+i*m] - a[k+j*m] ) );
      }
      if ( diff < tol )
      {
        unique = false;
        break;
      }
    }
    if ( unique )
    {
      for ( k = 0; k < m; k++ )
      {
        a[k+unique_num*m] = a[k+i*m];
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_num;
}
//****************************************************************************80

int r4col_sorted_tol_unique_count ( int m, int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORTED_TOL_UNIQUE_COUNT counts tolerably unique elements in a sorted R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The columns of the array may be ascending or descending sorted.
//
//    If the tolerance is large enough, then the concept of uniqueness
//    can become ambiguous.  If we have a tolerance of 1.5, then in the
//    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
//    one unique entry?  That would be because 1 may be regarded as unique,
//    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
//    be unique and so on.
//
//    This seems wrongheaded.  So I prefer the idea that an item is not
//    unique under a tolerance only if it is close to something that IS unique.
//    Thus, the unique items are guaranteed to cover the space if we include
//    a disk of radius TOL around each one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], a sorted array, containing
//    N columns of data.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int R4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
//
{
  float diff;
  int i;
  int i2;
  int i3;
  int j;
  int k;
  int *undx;
  bool unique;

  undx = new int[n];
//
//  Consider entry I = 0.
//  It is unique, so set the number of unique items to K.
//  Set the K-th unique item to I.
//  Set the representative of item I to the K-th unique item.
//
  i = 0;
  k = 0;
  undx[k] = i;
//
//  Consider entry I.
//
//  If it is unique, increase the unique count K, set the
//  K-th unique item to I, and set the representative of I to K.
//
//  If it is not unique, set the representative of item I to a
//  previously determined unique item that is close to it.
//
  for ( i = 1; i < n; i++ )
  {
    unique = true;

    for ( j = 0; j <= k; j++ )
    {
      i2 = undx[j];
      diff = 0.0;
      for ( i3 = 0; i3 < m; i3++ )
      {
        diff = r4_max ( diff, r4_abs ( a[i3+i*m] - a[i3+i2*m] ) );
      }
      if ( diff <= tol )
      {
        unique = false;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = i;
    }
  }

  delete [] undx;

  k = k + 1;

  return k;
}
//****************************************************************************80

void r4col_sorted_undex ( int m, int n, float a[], int unique_num,
  int undx[], int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORTED_UNDEX returns unique sorted indexes for a sorted R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The goal of this routine is to determine a vector UNDX,
//    which points, to the unique elements of A, in sorted order,
//    and a vector XDNU, which identifies, for each entry of A, the index of
//    the unique sorted element of A.
//
//    This is all done with index vectors, so that the elements of
//    A are never moved.
//
//    Assuming A is already sorted, we examine the entries of A in order,
//    noting the unique entries, creating the entries of XDNU and
//    UNDX as we go.
//
//    Once this process has been completed, the vector A could be
//    replaced by a compressed vector XU, containing the unique entries
//    of X in sorted order, using the formula
//
//      XU(*) = A(UNDX(*)).
//
//    We could then, if we wished, reconstruct the entire vector A, or
//    any element of it, by index, as follows:
//
//      A(I) = XU(XDNU(I)).
//
//    We could then replace A by the combination of XU and XDNU.
//
//    Later, when we need the I-th entry of A, we can locate it as
//    the XDNU(I)-th entry of XU.
//
//    Here is an example of a vector A, the unique sort and inverse unique
//    sort vectors and the compressed unique sorted vector.
//
//      I      A      XU  Undx  Xdnu
//    ----+------+------+-----+-----+
//      0 | 11.0 |  11.0    0     0
//      1 | 11.0 |  22.0    4     0
//      2 | 11.0 |  33.0    7     0
//      3 | 11.0 |  55.0    8     0
//      4 | 22.0 |                1
//      5 | 22.0 |                1
//      6 | 22.0 |                1
//      7 | 33.0 |                2
//      8 | 55.0 |                3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the data values.
//
//    Input, int N, the number of data values,
//
//    Input, float A[M*N], the data values.
//
//    Input, int UNIQUE_NUM, the number of unique values in A.
//    This value is only required for languages in which the size of
//    UNDX must be known in advance.
//
//    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
//
//    Output, int XDNU[N], the XDNU vector.
//
{
  float diff;
  int i;
  int j;
  int k;
//
//  Walk through the sorted array.
//
  i = 0;

  j = 0;
  undx[j] = i;

  xdnu[i] = j;

  for ( i = 1; i < n; i++ )
  {
    diff = 0.0;
    for ( k = 0; k < m; k++ )
    {
      diff = r4_max ( diff, r4_abs ( a[k+i*m] - a[k+undx[j]*m] ) );
    }
    if ( 0.0 < diff )
    {
      j = j + 1;
      undx[j] = i;
    }
    xdnu[i] = j;
  }
  return;
}
//****************************************************************************80

int r4col_sorted_unique ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORTED_UNIQUE keeps unique elements in a sorted R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The columns of the array can be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, float A(M,N).
//    On input, the sorted array of N columns of M-vectors.
//    On output, a sorted array of columns of M-vectors.
//
//    Output, int UNIQUE_NUM, the number of unique columns.
//
{
  bool equal;
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    equal = true;
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        equal = false;
        break;
      }
    }
    if ( !equal )
    {
      j1 = j1 + 1;
      for ( i = 0; i < m; i++ )
      {
        a[i+j1*m] = a[i+j2*m];
      }
    }
  }

  unique_num = j1 + 1;

  return unique_num;
}
//****************************************************************************80

int r4col_sorted_unique_count ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORTED_UNIQUE_COUNT counts unique elements in a sorted R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The columns of the array may be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], a sorted array, containing
//    N columns of data.
//
//    Output, int R4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
//
{
  bool equal;
  int i;
  int j1;
  int j2;
  int unique_num;

  unique_num = 0;

  if ( n <= 0 )
  {
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    equal = true;
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        equal = false;
        break;
      }
    }
    if ( !equal )
    {
      unique_num = unique_num + 1;
      j1 = j2;
    }
  }

  return unique_num;
}
//****************************************************************************80

void r4col_sortr_a ( int m, int n, float a[], int key )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SORTR_A ascending sorts one column of an R4COL, adjusting all entries.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, float A[M*N].
//    On input, an unsorted M by N array.
//    On output, rows of the array have been shifted in such
//    a way that column KEY of the array is in nondecreasing order.
//
//    Input, int KEY, the column in which the "key" value
//    is stored.  On output, column KEY of the array will be
//    in nondecreasing order.
//
{
  int i;
  int indx;
  int isgn;
  int j;

  if ( m <= 0 )
  {
    return;
  }

  if ( key < 1 || n < key )
  {
    cerr << "\n";
    cerr << "R4COL_SORTR_A - Fatal error!\n";
    cerr << "  The value of KEY is not a legal column index.\n";
    cerr << "  KEY = " << key << "\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( m, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      r4row_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      if ( a[i-1+(key-1)*m] < a[j-1+(key-1)*m] )
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

float *r4col_sum ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SUM sums the columns of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, float R4COL_SUM[N], the sums of the columns.
//
{
  float *colsum;
  int i;
  int j;

  colsum = new float[n];

  for ( j = 0; j < n; j++ )
  {
    colsum[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      colsum[j] = colsum[j] + a[i+j*m];
    }
  }
  return colsum;
}
//****************************************************************************80

void r4col_swap ( int m, int n, float a[], int j1, int j2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_SWAP swaps columns J1 and J2 of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, J1 = 2, J2 = 4
//
//      A = (
//        1.  2.  3.  4.
//        5.  6.  7.  8.
//        9. 10. 11. 12. )
//
//    Output:
//
//      A = (
//        1.  4.  3.  2.
//        5.  8.  7.  6.
//        9. 12. 11. 10. )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, float A[M*N], the M by N array.
//
//    Input, int J1, J2, the columns to be swapped.
//    These columns are 1-based.
//
{
  int i;
  float temp;

  if ( j1 < 1 || n < j1 || j2 < 1 || n < j2 )
  {
    cerr << "\n";
    cerr << "R4COL_SWAP - Fatal error!\n";
    cerr << "  J1 or J2 is out of bounds.\n";
    cerr << "  J1 =   " << j1 << "\n";
    cerr << "  J2 =   " << j2 << "\n";
    cerr << "  NCOL = " << n << "\n";
    exit ( 1 );
  }

  if ( j1 == j2 )
  {
    return;
  }

  for ( i = 0; i < m; i++ )
  {
    temp          = a[i+(j1-1)*m];
    a[i+(j1-1)*m] = a[i+(j2-1)*m];
    a[i+(j2-1)*m] = temp;
  }

  return;
}
//****************************************************************************80

float *r4col_to_r4vec ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_TO_R4VEC converts an R4COL to an R4VEC.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    This routine is not really useful in our C++ implementation, since
//    we actually store an M by N matrix exactly as a vector already.
//
//  Example:
//
//    M = 3, N = 4
//
//    A =
//      11 12 13 14
//      21 22 23 24
//      31 32 33 34
//
//    R4COL_TO_R4VEC = ( 11, 21, 31, 12, 22, 32, 13, 23, 33, 14, 24, 34 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the M by N array.
//
//    Output, float X[M*N], a vector containing the N columns of A.
//
{
  int i;
  int j;
  int k;
  float *x;

  x = new float[m*n];

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

void r4col_tol_undex ( int m, int n, float a[], int unique_num, float tol,
  int undx[], int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_TOL_UNDEX indexes tolerably unique entries of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The goal of this routine is to determine a vector UNDX,
//    which points to the unique elements of A, in sorted order,
//    and a vector XDNU, which identifies, for each entry of A, the index of
//    the unique sorted element of A.
//
//    This is all done with index vectors, so that the elements of
//    A are never moved.
//
//    The first step of the algorithm requires the indexed sorting
//    of A, which creates arrays INDX and XDNI.  (If all the entries
//    of A are unique, then these arrays are the same as UNDX and XDNU.)
//
//    We then use INDX to examine the entries of A in sorted order,
//    noting the unique entries, creating the entries of XDNU and
//    UNDX as we go.
//
//    Once this process has been completed, the vector A could be
//    replaced by a compressed vector XU, containing the unique entries
//    of A in sorted order, using the formula
//
//      XU(*) = A(UNDX(*)).
//
//    We could then, if we wished, reconstruct the entire vector A, or
//    any element of it, by index, as follows:
//
//      A(I) = XU(XDNU(I)).
//
//    We could then replace A by the combination of XU and XDNU.
//
//    Later, when we need the I-th entry of A, we can locate it as
//    the XDNU(I)-th entry of XU.
//
//    Here is an example of a vector A, the sort and inverse sort
//    index vectors, and the unique sort and inverse unique sort vectors
//    and the compressed unique sorted vector.
//
//      I     A  Indx  Xdni       XU  Undx  Xdnu
//    ----+-----+-----+-----+--------+-----+-----+
//      0 | 11.     0     0 |    11.     0     0
//      1 | 22.     2     4 |    22.     1     1
//      2 | 11.     5     1 |    33.     3     0
//      3 | 33.     8     7 |    55.     4     2
//      4 | 55.     1     8 |                  3
//      5 | 11.     6     2 |                  0
//      6 | 22.     7     5 |                  1
//      7 | 22.     3     6 |                  1
//      8 | 11.     4     3 |                  0
//
//    INDX(2) = 3 means that sorted item(2) is A(3).
//    XDNI(2) = 5 means that A(2) is sorted item(5).
//
//    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
//    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
//
//    XU(XDNU(I))) = X(I).
//    XU(I)        = X(UNDX(I)).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the data values.
//
//    Input, int N, the number of data values,
//
//    Input, float A[M*N], the data values.
//
//    Input, int UNIQUE_NUM, the number of unique values in A.
//    This value is only required for languages in which the size of
//    UNDX must be known in advance.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
//
//    Output, int XDNU[N], the XDNU vector.
//
{
  int base = 0;
  float diff;
  int i;
  int i2;
  int *indx;
  int j;
  int k;
  bool unique;
//
//  Implicitly sort the array.
//
  indx = r4col_sort_heap_index_a ( m, n, base, a );
//
//  Consider entry I = 0.
//  It is unique, so set the number of unique items to K.
//  Set the K-th unique item to I.
//  Set the representative of item I to the K-th unique item.
//
  i = 0;
  k = 0;
  undx[k] = indx[i];
  xdnu[indx[i]] = k;
//
//  Consider entry I.
//
//  If it is unique, increase the unique count K, set the
//  K-th unique item to I, and set the representative of I to K.
//
//  If it is not unique, set the representative of item I to a
//  previously determined unique item that is close to it.
//
  for ( i = 1; i < n; i++ )
  {
    unique = true;
    for ( j = 0; j <= k; j++ )
    {
      diff = 0.0;
      for ( i2 = 0; i2 < m; i2++ )
      {
        diff = r4_max ( diff, r4_abs ( a[i2+indx[i]*m] - a[i2+undx[j]*m] ) );
      }
      if ( diff <= tol )
      {
        unique = false;
        xdnu[indx[i]] = j;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = indx[i];
      xdnu[indx[i]] = k;
    }
  }
  delete [] indx;

  return;
}
//****************************************************************************80

int r4col_tol_unique_count ( int m, int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_TOL_UNIQUE_COUNT counts tolerably unique entries in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The columns of the array may be ascending or descending sorted.
//
//    If the tolerance is large enough, then the concept of uniqueness
//    can become ambiguous.  If we have a tolerance of 1.5, then in the
//    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
//    one unique entry?  That would be because 1 may be regarded as unique,
//    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
//    be unique and so on.
//
//    This seems wrongheaded.  So I prefer the idea that an item is not
//    unique under a tolerance only if it is close to something that IS unique.
//    Thus, the unique items are guaranteed to cover the space if we include
//    a disk of radius TOL around each one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the array of N columns of data.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int R4COL_TOL_UNIQUE_COUNT, the number of unique columns.
//
{
  int base = 0;
  float diff;
  int i;
  int i2;
  int *indx;
  int j;
  int k;
  bool unique;
  int *undx;

  undx = new int[n];
//
//  Implicitly sort the array.
//
  indx = r4col_sort_heap_index_a ( m, n, base, a );
//
//  Consider entry I = 0.
//  It is unique, so set the number of unique items to K.
//  Set the K-th unique item to I.
//  Set the representative of item I to the K-th unique item.
//
  i = 0;
  k = 0;
  undx[k] = indx[i];
//
//  Consider entry I.
//
//  If it is unique, increase the unique count K, set the
//  K-th unique item to I, and set the representative of I to K.
//
//  If it is not unique, set the representative of item I to a
//  previously determined unique item that is close to it.
//
  for ( i = 1; i < n; i++ )
  {
    unique = true;
    for ( j = 0; j <= k; j++ )
    {
      diff = 0.0;
      for ( i2 = 0; i2 < m; i2++ )
      {
        diff = r4_max ( diff, r4_abs ( a[i2+indx[i]*m] - a[i2+undx[j]*m] ) );
      }
      if ( diff <= tol )
      {
        unique = false;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = indx[i];
    }
  }
  delete [] indx;
  delete [] undx;

  k = k + 1;

  return k;
}
//****************************************************************************80

int *r4col_tol_unique_index ( int m, int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_TOL_UNIQUE_INDEX indexes tolerably unique entries in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
//    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
//    gathered in order, then
//
//      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of A.
//
//    Input, float A[M*N], the array.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int R4COL_TOL_UNIQUE_INDEX[N], the unique index.
//
{
  float diff;
  int i;
  int j1;
  int j2;
  int *unique_index;
  int unique_num;

  unique_index = new int[n];

  for ( j1 = 0; j1 < n; j1++ )
  {
    unique_index[j1] = -1;
  }
  unique_num = 0;

  for ( j1 = 0; j1 < n; j1++ )
  {
    if ( unique_index[j1] == -1 )
    {
      unique_index[j1] = unique_num;

      for ( j2 = j1 + 1; j2 < n; j2++ )
      {
        diff = 0.0;
        for ( i = 0; i < m; i++ )
        {
          diff = r4_max ( diff, r4_abs ( a[i+j1*m] - a[i+j2*m] ) );
        }
        if ( diff <= tol )
        {
          unique_index[j2] = unique_num;
        }
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_index;
}
//****************************************************************************80

void r4col_undex ( int m, int n, float a[], int unique_num, int undx[],
  int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_UNDEX indexes unique entries in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The goal of this routine is to determine a vector UNDX,
//    which points to the unique elements of A, in sorted order,
//    and a vector XDNU, which identifies, for each entry of A, the index of
//    the unique sorted element of A.
//
//    This is all done with index vectors, so that the elements of
//    A are never moved.
//
//    The first step of the algorithm requires the indexed sorting
//    of A, which creates arrays INDX and XDNI.  (If all the entries
//    of A are unique, then these arrays are the same as UNDX and XDNU.)
//
//    We then use INDX to examine the entries of A in sorted order,
//    noting the unique entries, creating the entries of XDNU and
//    UNDX as we go.
//
//    Once this process has been completed, the vector A could be
//    replaced by a compressed vector XU, containing the unique entries
//    of A in sorted order, using the formula
//
//      XU(*) = A(UNDX(*)).
//
//    We could then, if we wished, reconstruct the entire vector A, or
//    any element of it, by index, as follows:
//
//      A(I) = XU(XDNU(I)).
//
//    We could then replace A by the combination of XU and XDNU.
//
//    Later, when we need the I-th entry of A, we can locate it as
//    the XDNU(I)-th entry of XU.
//
//    Here is an example of a vector A, the sort and inverse sort
//    index vectors, and the unique sort and inverse unique sort vectors
//    and the compressed unique sorted vector.
//
//      I     A  Indx  Xdni       XU  Undx  Xdnu
//    ----+-----+-----+-----+--------+-----+-----+
//      0 | 11.     0     0 |    11.     0     0
//      1 | 22.     2     4 |    22.     1     1
//      2 | 11.     5     1 |    33.     3     0
//      3 | 33.     8     7 |    55.     4     2
//      4 | 55.     1     8 |                  3
//      5 | 11.     6     2 |                  0
//      6 | 22.     7     5 |                  1
//      7 | 22.     3     6 |                  1
//      8 | 11.     4     3 |                  0
//
//    INDX(2) = 3 means that sorted item(2) is A(3).
//    XDNI(2) = 5 means that A(2) is sorted item(5).
//
//    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
//    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
//
//    XU(XDNU(I))) = A(I).
//    XU(I)        = A(UNDX(I)).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the data values.
//
//    Input, int N, the number of data values,
//
//    Input, float A[M*N], the data values.
//
//    Input, int UNIQUE_NUM, the number of unique values in A.
//    This value is only required for languages in which the size of
//    UNDX must be known in advance.
//
//    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
//
//    Output, int XDNU[N], the XDNU vector.
//
{
  int base = 0;
  float diff;
  int i;
  int *indx;
  int j;
  int k;
//
//  Implicitly sort the array.
//
  indx = r4col_sort_heap_index_a ( m, n, base, a );
//
//  Walk through the implicitly sorted array.
//
  i = 0;

  j = 0;
  undx[j] = indx[i];

  xdnu[indx[i]] = j;

  for ( i = 1; i < n; i++ )
  {
    diff = 0.0;
    for ( k = 0; k < m; k++ )
    {
      diff = r4_max ( diff, r4_abs ( a[k+indx[i]*m] - a[k+undx[j]*m] ) );
    }
    if ( 0.0 < diff )
    {
      j = j + 1;
      undx[j] = indx[i];
    }
    xdnu[indx[i]] = j;
  }
  delete [] indx;

  return;
}
//****************************************************************************80

int r4col_unique_count ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_UNIQUE_COUNT counts unique entries in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    The columns of the array may be ascending or descending sorted.
//
//    If the tolerance is large enough, then the concept of uniqueness
//    can become ambiguous.  If we have a tolerance of 1.5, then in the
//    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
//    one unique entry?  That would be because 1 may be regarded as unique,
//    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
//    be unique and so on.
//
//    This seems wrongheaded.  So I prefer the idea that an item is not
//    unique under a tolerance only if it is close to something that IS unique.
//    Thus, the unique items are guaranteed to cover the space if we include
//    a disk of radius TOL around each one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the array of N columns of data.
//
//    Output, int R4COL_UNIQUE_COUNT, the number of unique columns.
//
{
  float diff;
  int i;
  int j1;
  int j2;
  bool *unique;
  int unique_num;

  unique_num = 0;

  unique = new bool[n];

  for ( j1 = 0; j1 < n; j1++ )
  {
    unique_num = unique_num + 1;
    unique[j1] = true;

    for ( j2 = 0; j2 < j1; j2++ )
    {
      diff = 0.0;
      for ( i = 0; i < m; i++ )
      {
        diff = r4_max ( diff, r4_abs ( a[i+j1*m] - a[i+j2*m] ) );
      }
      if ( diff == 0.0 )
      {
        unique_num = unique_num - 1;
        unique[j1] = false;
        break;
      }
    }
  }

  delete [] unique;

  return unique_num;
}
//****************************************************************************80

int *r4col_unique_index ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_UNIQUE_INDEX indexes unique entries in an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
//    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
//    gathered in order, then
//
//      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of A.
//    The length of an "element" of A, and the number of "elements".
//
//    Input, float A[M*N], the array.
//
//    Output, int R4COL_UNIQUE_INDEX[N], the unique index.
//
{
  float diff;
  int i;
  int j1;
  int j2;
  int *unique_index;
  int unique_num;

  unique_index = new int[n];

  for ( j1 = 0; j1 < n; j1++ )
  {
    unique_index[j1] = -1;
  }
  unique_num = 0;

  for ( j1 = 0; j1 < n; j1++ )
  {
    if ( unique_index[j1] == -1 )
    {
      unique_index[j1] = unique_num;

      for ( j2 = j1 + 1; j2 < n; j2++ )
      {
        diff = 0.0;
        for ( i = 0; i < m; i++ )
        {
          diff = r4_max ( diff, r4_abs ( a[i+j1*m] - a[i+j2*m] ) );
        }
        if ( diff == 0.0 )
        {
          unique_index[j2] = unique_num;
        }
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_index;
}
//****************************************************************************80

float *r4col_variance ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4COL_VARIANCE returns the variances of an R4COL.
//
//  Discussion:
//
//    An R4COL is an M by N array of R4's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, float A[M*N], the array whose variances are desired.
//
//    Output, float R4COL_VARIANCE[N], the variances of the rows.
//
{
  int i;
  int j;
  float mean;
  float *variance;

  variance = new float[n];

  for ( j = 0; j < n; j++ )
  {
    mean = 0.0;
    for ( i = 0; i < m; i++ )
    {
      mean = mean + a[i+j*m];
    }
    mean = mean / ( float ) ( m );

    variance[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      variance[j] = variance[j] + pow ( a[i+j*m] - mean, 2 );
    }

    if ( 1 < m )
    {
      variance[j] = variance[j] / ( float ) ( m - 1 );
    }
    else
    {
      variance[j] = 0.0;
    }
  }

  return variance;
}
//****************************************************************************80

float r4int_to_r4int ( float rmin, float rmax, float r, float r2min,
  float r2max )

//****************************************************************************80
//
//  Purpose:
//
//    R4INT_TO_R4INT maps one R4 interval to another.
//
//  Discussion:
//
//    The formula used is
//
//      R2 := R2MIN + ( R2MAX - R2MIN ) * ( R - RMIN ) / ( RMAX - RMIN )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float RMIN, RMAX, the first range.
//
//    Input, float R, the number to be converted.
//
//    Input, float R2MAX, R2MIN, the second range.
//
//    Output, float R4INT_TO_R4INT, the corresponding value in
//    the range [R2MIN,R2MAX].
//
{
  float  r2;

  if ( rmax == rmin )
  {
    r2 = ( r2max + r2min ) / 2.0;
  }
  else
  {
    r2 = ( ( ( rmax - r        ) * r2min
           + (        r - rmin ) * r2max )
           / ( rmax     - rmin ) );
  }

  return r2;
}
//****************************************************************************80

int r4int_to_i4int ( float rmin, float rmax, float r, int imin, int imax )

//****************************************************************************80
//
//  Purpose:
//
//    R4INT_TO_I4INT maps an R4 interval to an integer interval.
//
//  Discussion:
//
//    The formula used is
//
//      I := IMIN + ( IMAX - IMIN ) * ( R - RMIN ) / ( RMAX - RMIN )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float RMIN, RMAX, the range.
//
//    Input, float R, the number to be converted.
//
//    Input, int IMAX, IMIN, the integer range.
//
//    Output, int I, the corresponding value in the range [IMIN,IMAX].
//
{
  int i;

  if ( rmax == rmin )
  {
    i = ( imax + imin ) / 2;
  }
  else
  {
    i = r4_nint (
      ( ( rmax - r        ) * ( float ) ( imin )
      + (        r - rmin ) * ( float ) ( imax ) )
      / ( rmax     - rmin ) );
  }

  return i;
}
//****************************************************************************80

float *r4mat_border_add ( int m, int n, float table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_BORDER_ADD adds a "border" to an R4MAT.
//
//  Discussion:
//
//    We suppose the input data gives values of a quantity on nodes
//    in the interior of a 2D grid, and we wish to create a new table
//    with additional positions for the nodes that would be on the
//    border of the 2D grid.
//
//                  0 0 0 0 0 0
//      * * * *     0 * * * * 0
//      * * * * --> 0 * * * * 0
//      * * * *     0 * * * * 0
//                  0 0 0 0 0 0
//
//    The illustration suggests the situation in which a 3 by 4 array
//    is input, and a 5 by 6 array is to be output.
//
//    The old data is shifted to its correct positions in the new array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, float TABLE[M*N], the table data.
//
//    Output, float TABLE2[(M+2)*(N+2)], the augmented table data.
//
{
  int i;
  int j;
  float *table2;

  table2 = new float[(m+2)*(n+2)];

  for ( j = 0; j < n+2; j++ )
  {
    for ( i = 0; i < m+2; i++ )
    {
      if ( i == 0 || i == m+1 || j == 0 || j == n+1 )
      {
        table2[i+j*(m+2)] = 0.0;
      }
      else
      {
        table2[i+j*(m+2)] = table[(i-1)+(j-1)*m];
      }
    }
  }
  return table2;
}
//****************************************************************************80

float *r4mat_border_cut ( int m, int n, float table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_BORDER_CUT cuts the "border" of an R4MAT.
//
//  Discussion:
//
//    We suppose the input data gives values of a quantity on nodes
//    on a 2D grid, and we wish to create a new table corresponding only
//    to those nodes in the interior of the 2D grid.
//
//      0 0 0 0 0 0
//      0 * * * * 0    * * * *
//      0 * * * * 0 -> * * * *
//      0 * * * * 0    * * * *
//      0 0 0 0 0 0
//
//    The illustration suggests the situation in which a 5 by 6 array
//    is input, and a 3 by 4 array is to be output.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, float TABLE[M*N], the table data.
//
//    Output, float TABLE2[(M-2)*(N-2)], the "interior" table data.
//
{
  int i;
  int j;
  float *table2;

  if ( m <= 2 || n <= 2 )
  {
    return NULL;
  }

  table2 = new float[(m-2)*(n-2)];

  for ( j = 0; j < n-2; j++ )
  {
    for ( i = 0; i < m-2; i++ )
    {
      table2[i+j*(m-2)] = table[(i+1)+(j+1)*m];
    }
  }

  return table2;
}
//****************************************************************************80

float *r4mat_cholesky_factor ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The matrix must be symmetric and positive semidefinite.
//
//    For a positive semidefinite symmetric matrix A, the Cholesky factorization
//    is a lower triangular matrix L such that:
//
//      A = L * L'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix A.
//
//    Input, float A[N*N], the N by N matrix.
//
//    Output, float R4MAT_CHOLESKY_FACTOR[N*N], the N by N lower triangular
//    Cholesky factor.
//
{
  float *c;
  int i;
  int j;
  int k;
  float sum2;

  c = r4mat_copy_new ( n, n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      c[i+j*n] = 0.0;
    }
    for ( i = j; i < n; i++ )
    {
      sum2 = c[j+i*n];
      for ( k = 0; k < j; k++ )
      {
        sum2 = sum2 - c[j+k*n] * c[i+k*n];
      }
      if ( i == j )
      {
        if ( sum2 <= 0.0 )
        {
          cerr << "\n";
          cerr << "R4MAT_CHOLESKY_FACTOR - Fatal error!\n";
          cerr << "  Matrix is not positive definite.\n";
          exit ( 1 );
        }
        c[i+j*n] = sqrt ( sum2 );
      }
      else
      {
        if ( c[j+j*n] != 0.0 )
        {
          c[i+j*n] = sum2 / c[j+j*n];
        }
        else
        {
          c[i+j*n] = 0.0;
        }
      }
    }
  }

  return c;
}
//****************************************************************************80

float *r4mat_cholesky_solve ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix A.
//
//    Input, float A[N*N], the N by N Cholesky factor of the
//    system matrix.
//
//    Input, float B[N], the right hand side of the linear system.
//
//    Output, float R4MAT_CHOLESKY_SOLVE[N], the solution of the linear system.
//
{
  float *x;
  float *y;
//
//  Solve L * y = b.
//
  y = r4mat_l_solve ( n, a, b );
//
//  Solve L' * x = y.
//
  x = r4mat_lt_solve ( n, a, y );

  delete [] y;

  return x;
}
//****************************************************************************80

float *r4mat_choresky_factor ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_CHORESKY_FACTOR computes the "Choresky" factor of a symmetric R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The matrix must be symmetric and positive semidefinite.
//
//    For a positive semidefinite symmetric matrix A, the Cholesky factorization
//    is an upper triangular matrix R such that:
//
//      A = R * R'
//
//    Note that the usual Cholesky factor is a LOWER triangular matrix L
//    such that
//
//      A = L * L'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix A.
//
//    Input, float A[N*N], the N by N matrix.
//
//    Output, float R4MAT_CHORESKY_FACTOR[N*N], the N by N upper triangular
//    "Choresky" factor.
//
{
  float *c;
  int i;
  int j;
  int k;
  float sum2;

  c = r4mat_copy_new ( n, n, a );

  r4mat_flip_rows ( n, n, c );
  r4mat_flip_cols ( n, n, c );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      c[i+j*n] = 0.0;
    }
    for ( i = j; i < n; i++ )
    {
      sum2 = c[j+i*n];
      for ( k = 0; k < j; k++ )
      {
        sum2 = sum2 - c[j+k*n] * c[i+k*n];
      }
      if ( i == j )
      {
        if ( sum2 <= 0.0 )
        {
          cerr << "\n";
          cerr << "R4MAT_CHORESKY_FACTOR - Fatal error!\n";
          cerr << "  Matrix is not positive definite.\n";
          exit ( 1 );
        }
        c[i+j*n] = sqrt ( sum2 );
      }
      else
      {
        if ( c[j+j*n] != 0.0 )
        {
          c[i+j*n] = sum2 / c[j+j*n];
        }
        else
        {
          c[i+j*n] = 0.0;
        }
      }
    }
  }

  r4mat_flip_cols ( n, n, c );
  r4mat_flip_rows ( n, n, c );

  return c;
}
//****************************************************************************80

void r4mat_copy ( int m, int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_COPY copies one R4MAT to another.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A1[M*N], the matrix to be copied.
//
//    Output, float A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
}
//****************************************************************************80

float *r4mat_copy_new ( int m, int n, float a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_COPY_NEW copies one R4MAT to a "new" R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A1[M*N], the matrix to be copied.
//
//    Output, float R4MAT_COPY_NEW[M*N], the copy of A1.
//
{
  float *a2;
  int i;
  int j;

  a2 = new float[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}
//****************************************************************************80

void r4mat_delete ( float **a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DELETE frees memory associated with an R4MAT.
//
//  Discussion:
//
//    This function releases the memory associated with an array that was 
//    created by a command like
//      float **a;
//      a = r4mat_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float **A, the pointer to the array.
//
//    Input, int M, N, the number of rows and columns in the array.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}
//****************************************************************************80

float r4mat_det ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DET computes the determinant of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    Original FORTRAN77 version by Helmut Spaeth
//    C++ version by John Burkardt
//
//  Reference:
//
//    Helmut Spaeth,
//    Cluster Analysis Algorithms
//    for Data Reduction and Classification of Objects,
//    Ellis Horwood, 1980, page 125-127.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, float A[N*N], the matrix whose determinant is desired.
//
//    Output, float R4MAT_DET, the determinant of the matrix.
//
{
  float *b;
  float det;
  int i;
  int j;
  int k;
  int kk;
  int m;
  float temp;

  b = new float[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  det = 1.0;

  for ( k = 1; k <= n; k++ )
  {
    m = k;
    for ( kk = k+1; kk <= n; kk++ )
    {
      if ( r4_abs ( b[m-1+(k-1)*n] ) < r4_abs ( b[kk-1+(k-1)*n] ) )
      {
        m = kk;
      }
    }

    if ( m != k )
    {
      det = -det;

      temp = b[m-1+(k-1)*n];
      b[m-1+(k-1)*n] = b[k-1+(k-1)*n];
      b[k-1+(k-1)*n] = temp;
    }

    det = det * b[k-1+(k-1)*n];

    if ( b[k-1+(k-1)*n] != 0.0 )
    {
      for ( i = k+1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / b[k-1+(k-1)*n];
      }

      for ( j = k+1; j <= n; j++ )
      {
        if ( m != k )
        {
          temp = b[m-1+(j-1)*n];
          b[m-1+(j-1)*n] = b[k-1+(j-1)*n];
          b[k-1+(j-1)*n] = temp;
        }
        for ( i = k+1; i <= n; i++ )
        {
          b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * b[k-1+(j-1)*n];
        }
      }
    }
  }

  delete [] b;

  return det;
}
//****************************************************************************80

float r4mat_det_2d ( float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DET_2D computes the determinant of a 2 by 2 R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Discussion:
//
//    The determinant of a 2 by 2 matrix is
//
//      a11 * a22 - a12 * a21.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[2*2], the matrix whose determinant is desired.
//
//    Output, float R4MAT_DET_2D, the determinant of the matrix.
//
{
  float det;

  det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];

  return det;
}
//****************************************************************************80

float r4mat_det_3d ( float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DET_3D computes the determinant of a 3 by 3 R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The determinant of a 3 by 3 matrix is
//
//        a11 * a22 * a33 - a11 * a23 * a32
//      + a12 * a23 * a31 - a12 * a21 * a33
//      + a13 * a21 * a32 - a13 * a22 * a31
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[3*3], the matrix whose determinant is desired.
//
//    Output, float R4MAT_DET_3D, the determinant of the matrix.
//
{
  float det;

  det =
      a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
    + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
    + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );

  return det;
}
//****************************************************************************80

float r4mat_det_4d ( float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DET_4D computes the determinant of a 4 by 4 R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[4*4], the matrix whose determinant is desired.
//
//    Output, float R4MAT_DET_4D, the determinant of the matrix.
//
{
  float det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}
//****************************************************************************80

float r4mat_det_5d ( float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DET_5D computes the determinant of a 5 by 5 R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[5*5], the matrix whose determinant is desired.
//
//    Output, float R4MAT_DET_5D, the determinant of the matrix.
//
{
  float b[4*4];
  float det;
  int i;
  int inc;
  int j;
  int k;
  float sign;
//
//  Expand the determinant into the sum of the determinants of the
//  five 4 by 4 matrices created by dropping row 1, and column k.
//
  det = 0.0;
  sign = 1.0;

  for ( k = 0; k < 5; k++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      for ( j = 0; j < 4; j++ )
      {
        if ( j < k )
        {
          inc = 0;
        }
        else
        {
          inc = 1;
        }
        b[i+j*4] = a[i+1+(j+inc)*5];
      }
    }

    det = det + sign * a[0+k*5] * r4mat_det_4d ( b );

    sign = - sign;
  }

  return det;
}
//****************************************************************************80

void r4mat_diag_add_scalar ( int n, float a[], float s )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input/output, float A[N*N], the N by N matrix to be modified.
//
//    Input, float S, the value to be added to the diagonal
//    of the matrix.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = a[i+i*n] + s;
  }

  return;
}
//****************************************************************************80

void r4mat_diag_add_vector ( int n, float a[], float v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DIAG_ADD_VECTOR adds a vector to the diagonal of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input/output, float A[N*N], the N by N matrix.
//
//    Input, float V[N], the vector to be added to the diagonal of A.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = a[i+i*n] + v[i];
  }

  return;
}
//****************************************************************************80

float *r4mat_diag_get_vector ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, float A[N*N], the N by N matrix.
//
//    Output, float R4MAT_DIAG_GET_VECTOR[N], the diagonal entries
//    of the matrix.
//
{
  int i;
  float *v;

  v = new float[n];

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i+i*n];
  }

  return v;
}
//****************************************************************************80

void r4mat_diag_set_scalar ( int n, float a[], float s )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DIAG_SET_SCALAR sets the diagonal of an R4MAT to a scalar value.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input/output, float A[N*N], the N by N matrix to be modified.
//
//    Input, float S, the value to be assigned to the diagonal
//    of the matrix.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = s;
  }

  return;
}
//****************************************************************************80

void r4mat_diag_set_vector ( int n, float a[], float v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DIAG_SET_VECTOR sets the diagonal of an R4MAT to a vector.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input/output, float A[N*N], the N by N matrix.
//
//    Input, float V[N], the vector to be assigned to the
//    diagonal of A.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = v[i];
  }

  return;
}
//****************************************************************************80

float *r4mat_expand_linear ( int m, int n, float x[], int mfat, int nfat )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_EXPAND_LINEAR linearly interpolates new data into an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    In this routine, the expansion is specified by giving the number
//    of intermediate values to generate between each pair of original
//    data rows and columns.
//
//    The interpolation is not actually linear.  It uses the functions
//
//      1, x, y, and xy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of input data.
//
//    Input, float X[M*N], the original data.
//
//    Input, int MFAT, NFAT, the number of data values to interpolate
//    between each row, and each column, of original data values.
//
//    Output, float XFAT[M2*N2], the fattened data, where
//    M2 = (M-1)*(MFAT+1)+1,
//    N2 = (N-1)*(NFAT+1)+1.
//
{
  int i;
  int ihi;
  int ii;
  int iii;
  int ip1;
  int j;
  int jhi;
  int jj;
  int jjj;
  int jp1;
  int m2;
  int n2;
  float s;
  float t;
  float x00;
  float x01;
  float x10;
  float x11;
  float *xfat;

  m2 = ( m - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( n - 1 ) * ( nfat + 1 ) + 1;

  xfat = new float[m2*n2];

  for ( i = 1; i <= m; i++ )
  {
    if ( i < m )
    {
      ihi = mfat;
    }
    else
    {
      ihi = 0;
    }

    for ( j = 1; j <= n; j++ )
    {
      if ( j < n )
      {
        jhi = nfat;
      }
      else
      {
        jhi = 0;
      }

      if ( i < m )
      {
        ip1 = i + 1;
      }
      else
      {
        ip1 = i;
      }

      if ( j < n )
      {
        jp1 = j + 1;
      }
      else
      {
        jp1 = j;
      }

      x00 = x[i-1+(j-1)*m];
      x10 = x[ip1-1+(j-1)*m];
      x01 = x[i-1+(jp1-1)*m];
      x11 = x[ip1-1+(jp1-1)*m];

      for ( ii = 0; ii <= ihi; ii++ )
      {
        s = ( float ) ( ii ) / ( float ) ( ihi + 1 );

        for ( jj = 0; jj <= jhi; jj++ )
        {
          t = ( float ) ( jj ) / ( float ) ( jhi + 1 );

          iii = 1 + ( i - 1 ) * ( mfat + 1 ) + ii;
          jjj = 1 + ( j - 1 ) * ( nfat + 1 ) + jj;

          xfat[iii-1+(jjj-1)*m2] =
                                            x00
              + s     * (       x10       - x00 )
              + t     * (             x01 - x00 )
              + s * t * ( x11 - x10 - x01 + x00 );
        }
      }
    }
  }

  return xfat;
}
//****************************************************************************80

float *r4mat_expand_linear2 ( int m, int n, float a[], int m2, int n2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_EXPAND_LINEAR2 expands an R4MAT by linear interpolation.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    In this version of the routine, the expansion is indicated
//    by specifying the dimensions of the expanded array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in A.
//
//    Input, float A(M,N), a "small" M by N array.
//
//    Input, int M2, N2, the number of rows and columns in A2.
//
//    Output, float R4MAT_EXPAND_LINEAR2[M2*N2], the expanded array,
//    which contains an interpolated version of the data in A.
//
{
  float *a2;
  int i;
  int i1;
  int i2;
  int j;
  int j1;
  int j2;
  float r;
  float r1;
  float r2;
  float s;
  float s1;
  float s2;

  a2 = new float[m2*n2];

  for ( i = 1; i <= m2; i++ )
  {
    if ( m2 == 1 )
    {
      r = 0.5;
    }
    else
    {
      r = ( float ) ( i - 1 ) / ( float ) ( m2 - 1 );
    }

    i1 = 1 + ( int ) ( r * ( float ) ( m - 1 ) );
    i2 = i1 + 1;

    if ( m < i2 )
    {
      i1 = m - 1;
      i2 = m;
    }

    r1 = ( float ) ( i1 - 1 ) / ( float ) ( m - 1 );
    r2 = ( float ) ( i2 - 1 ) / ( float ) ( m - 1 );

    for ( j = 1; j <= n2; j++ )
    {
      if ( n2 == 1 )
      {
        s = 0.5;
      }
      else
      {
        s = ( float ) ( j - 1 ) / ( float ) ( n2 - 1 );
      }

      j1 = 1 + ( int ) ( s * ( float ) ( n - 1 ) );
      j2 = j1 + 1;

      if ( n < j2 )
      {
        j1 = n - 1;
        j2 = n;
      }

      s1 = ( float ) ( j1 - 1 ) / ( float ) ( n - 1 );
      s2 = ( float ) ( j2 - 1 ) / ( float ) ( n - 1 );

      a2[i-1+(j-1)*m2] =
        ( ( r2 - r ) * ( s2 - s ) * a[i1-1+(j1-1)*m]
        + ( r - r1 ) * ( s2 - s ) * a[i2-1+(j1-1)*m]
        + ( r2 - r ) * ( s - s1 ) * a[i1-1+(j2-1)*m]
        + ( r - r1 ) * ( s - s1 ) * a[i2-1+(j2-1)*m] )
        / ( ( r2 - r1 ) * ( s2 - s1 ) );
    }
  }

  return a2;
}
//****************************************************************************80

void r4mat_flip_cols ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_FLIP_COLS swaps the columns of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
//
//    To "flip" the columns of an R4MAT is to start with something like
//
//      11 12 13 14 15
//      21 22 23 24 25
//      31 32 33 34 35
//      41 42 43 44 45
//      51 52 53 54 55
//
//    and return
//
//      15 14 13 12 11
//      25 24 23 22 21
//      35 34 33 32 31
//      45 44 43 42 41
//      55 54 53 52 51
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, float A[M*N], the matrix whose columns are to be flipped.
//
{
  int i;
  int j;
  float t;

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < ( n / 2 ); j++ )
    {
      t              = a[i+     j *m];
      a[i+     j *m] = a[i+(n-1-j)*m];
      a[i+(n-1-j)*m] = t;
    }
  }
  return;
}
//****************************************************************************80

void r4mat_flip_rows ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_FLIP_ROWS swaps the rows of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
//
//    To "flip" the rows of an R4MAT is to start with something like
//
//      11 12 13 14 15
//      21 22 23 24 25
//      31 32 33 34 35
//      41 42 43 44 45
//      51 52 53 54 55
//
//    and return
//
//      51 52 53 54 55
//      41 42 43 44 45
//      31 32 33 34 35
//      21 22 23 24 25
//      11 12 13 14 15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, float A[M*N], the matrix whose rows are to be flipped.
//
{
  int i;
  int j;
  float t;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < ( m / 2 ); i++ )
    {
      t            = a[    i+j*m];
      a[    i+j*m] = a[m-1-i+j*m];
      a[m-1-i+j*m] = t;
    }
  }
  return;
}
//****************************************************************************80

float *r4mat_givens_post ( int n, float a[], int row, int col )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_GIVENS_POST computes the Givens postmultiplier rotation matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Givens post-multiplier matrix G(ROW,COL) has the property that
//    the (ROW,COL)-th entry of A*G is zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices A and G.
//
//    Input, float A[N*N], the matrix to be operated upon.
//
//    Input, int ROW, COL, the row and column of the
//    entry of A*G which is to be zeroed out.
//
//    Output, float R4MAT_GIVENS_POST[N*N], the Givens rotation matrix.
//    G is an orthogonal matrix, that is, the inverse of
//    G is the transpose of G.
//
{
  float *g;
  float theta;

  g = r4mat_identity ( n );

  theta = atan2 ( a[row-1+(col-1)*n], a[row-1+(row-1)*n] );

  g[row-1+(row-1)*n] =  cos ( theta );
  g[row-1+(col-1)*n] = -sin ( theta );
  g[col-1+(row-1)*n] =  sin ( theta );
  g[col-1+(col-1)*n] =  cos ( theta );

  return g;
}
//****************************************************************************80

float *r4mat_givens_pre ( int n, float a[], int row, int col )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_GIVENS_PRE computes the Givens premultiplier rotation matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Givens premultiplier rotation matrix G(ROW,COL) has the
//    property that the (ROW,COL)-th entry of G*A is zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices A and G.
//
//    Input, float A[N*N], the matrix to be operated upon.
//
//    Input, int ROW, COL, the row and column of the
//    entry of the G*A which is to be zeroed out.
//
//    Output, float R4MAT_GIVENS_PRE[N*N], the Givens rotation matrix.
//    G is an orthogonal matrix, that is, the inverse of
//    G is the transpose of G.
//
{
  float *g;
  float theta;

  g = r4mat_identity ( n );

  theta = atan2 ( a[row-1+(col-1)*n], a[col-1+(col-1)*n] );

  g[row-1+(row-1)*n] =  cos ( theta );
  g[row-1+(col-1)*n] = -sin ( theta );
  g[col-1+(row-1)*n] =  sin ( theta );
  g[col-1+(col-1)*n] =  cos ( theta );

  return g;
}
//****************************************************************************80

float *r4mat_hess ( float (*fx) ( int n, float x[] ), int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_HESS approximates a Hessian matrix via finite differences.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    H(I,J) = d2 F / d X(I) d X(J)
//
//    The values returned by this routine will be only approximate.
//    In some cases, they will be so poor that they are useless.
//    However, one of the best applications of this routine is for
//    checking your own Hessian calculations, since as Heraclitus
//    said, you'll never get the same result twice when you differentiate
//    a complicated expression by hand.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float *FX ( int N, float X[] ), the name of the user
//    function routine.
//
//    Input, int N, the number of variables.
//
//    Input, float X[N], the values of the variables.
//
//    Output, float H[N*N], the approximated N by N Hessian matrix.
//
{
  float eps;
  float f00;
  float fmm;
  float fmp;
  float fpm;
  float fpp;
  float *h;
  int i;
  int j;
  float *s;
  float xi;
  float xj;
//
//  Choose the stepsizes.
//
  s = new float[n];

  eps = pow ( ( double ) r4_epsilon ( ), 0.33 );

  for ( i = 0; i < n; i++ )
  {
    s[i] = eps * r4_max ( r4_abs ( x[i] ), 1.0 );
  }
//
//  Calculate the diagonal elements.
//
  h = new float[n*n];

  for ( i = 0; i < n; i++ )
  {
    xi = x[i];

    f00 = fx ( n, x );

    x[i] = xi + s[i];
    fpp = fx ( n, x );

    x[i] = xi - s[i];
    fmm = fx ( n, x );

    h[i+i*n] = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s[i] / s[i];

    x[i] = xi;
  }
//
//  Calculate the off diagonal elements.
//
  for ( i = 0; i < n; i++ )
  {
    xi = x[i];

    for ( j = i+1; j < n; j++ )
    {
      xj = x[j];

      x[i] = xi + s[i];
      x[j] = xj + s[j];
      fpp = fx ( n, x );

      x[i] = xi + s[i];
      x[j] = xj - s[j];
      fpm = fx ( n, x );

      x[i] = xi - s[i];
      x[j] = xj + s[j];
      fmp = fx ( n, x );

      x[i] = xi - s[i];
      x[j] = xj - s[j];
      fmm = fx ( n, x );

      h[j+i*n] = ( ( fpp - fpm ) + ( fmm - fmp ) ) / ( 4.0 * s[i] * s[j] );

      h[i+j*n] = h[j+i*n];

      x[j] = xj;
    }
    x[i] = xi;
  }

  return h;
}
//****************************************************************************80

void r4mat_house_axh ( int n, float a[], float v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Householder matrix H(V) is defined by
//
//      H(V) = I - 2 * v * v' / ( v' * v )
//
//    This routine is not particularly efficient.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Input/output, float A[N*N], on input, the matrix to be postmultiplied.
//    On output, A has been replaced by A*H.
//
//    Input, float V[N], a vector defining a Householder matrix.
//
{
  int i;
  int j;
  int k;
  float v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }
//
//  Compute A*H' = A*H
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        a[i+j*n] = a[i+j*n] - 2.0 * a[i+k*n] * v[k] * v[j] / v_normsq;
      }
    }
  }

  return;
}
//****************************************************************************80

float *r4mat_house_axh_new ( int n, float a[], float v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_HOUSE_AXH_NEW computes A*H where H is a compact Householder matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Householder matrix H(V) is defined by
//
//      H(V) = I - 2 * v * v' / ( v' * v )
//
//    This routine is not particularly efficient.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Input, float A[N*N], the matrix to be postmultiplied.
//
//    Input, float V[N], a vector defining a Householder matrix.
//
//    Output, float R4MAT_HOUSE_AXH[N*N], the product A*H.
//
{
  float *ah;
  int i;
  int j;
  int k;
  float v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }
//
//  Compute A*H' = A*H
//
  ah = new float[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      ah[i+j*n] = a[i+j*n];
      for ( k = 0; k < n; k++ )
      {
        ah[i+j*n] = ah[i+j*n] - 2.0 * a[i+k*n] * v[k] * v[j] / v_normsq;
      }
    }
  }

  return ah;
}
//****************************************************************************80

float *r4mat_house_form ( int n, float v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_HOUSE_FORM constructs a Householder matrix from its compact form.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    H(v) = I - 2 * v * v' / ( v' * v )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, float V[N], the vector defining the Householder matrix.
//
//    Output, float R4MAT_HOUSE_FORM[N*N], the Householder matrix.
//
{
  float beta;
  float *h;
  int i;
  int j;
//
//  Compute the L2 norm of V.
//
  beta = 0.0;
  for ( i = 0; i < n; i++ )
  {
    beta = beta + v[i] * v[i];
  }
//
//  Form the matrix H.
//
  h = r4mat_identity ( n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      h[i+j*n] = h[i+j*n] - 2.0 * v[i] * v[j] / beta;
    }
  }

  return h;
}
//****************************************************************************80

float *r4mat_house_hxa ( int n, float a[], float v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Householder matrix H(V) is defined by
//
//      H(V) = I - 2 * v * v' / ( v' * v )
//
//    This routine is not particularly efficient.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Input, float A[N*N], the matrix to be premultiplied.
//
//    Input, float V[N], a vector defining a Householder matrix.
//
//    Output, float R4MAT_HOUSE_HXA[N*N], the product H*A.
//
{
  float *ha;
  int i;
  int j;
  int k;
  float v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }
//
//  Compute A*H' = A*H
//
  ha = new float[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      ha[i+j*n] = a[i+j*n];
      for ( k = 0; k < n; k++ )
      {
        ha[i+j*n] = ha[i+j*n] - 2.0 * v[i] * v[k] * a[k+j*n] / v_normsq;
      }
    }
  }

  return ha;
}
//****************************************************************************80

float *r4mat_house_post ( int n, float a[], int row, int col )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_HOUSE_POST computes a Householder post-multiplier matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    H(ROW,COL) has the property that the ROW-th column of
//    A*H(ROW,COL) is zero from entry COL+1 to the end.
//
//    In the most common case, where a QR factorization is being computed,
//    ROW = COL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//
//    Input, float A[N*N], the matrix whose Householder matrix
//    is to be computed.
//
//    Input, int ROW, COL, specify the location of the
//    entry of the matrix A which is to be preserved.  The entries in
//    the same row, but higher column, will be zeroed out if
//    A is postmultiplied by H.
//
//    Output, float R4MAT_HOUSE_POST[N*N], the Householder matrix.
//
{
  float *h;
  int j;
  float *v;
  float *w;
//
//  Set up the vector V.
//
  w = new float[n];

  for ( j = 0; j < col-1; j++ )
  {
    w[j] = 0.0;
  }
  for ( j = col-1; j < n; j++ )
  {
    w[j] = a[row+j*n];
  }

  v = r4vec_house_column ( n, w, col );
//
//  Form the matrix H(V).
//
  h = r4mat_house_form ( n, v );

  delete [] v;
  delete [] w;

  return h;
}
//****************************************************************************80

float *r4mat_house_pre ( int n, float a[], int row, int col )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_HOUSE_PRE computes a Householder pre-multiplier matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    H(ROW,COL) has the property that the COL-th column of
//    H(ROW,COL)*A is zero from entry ROW+1 to the end.
//
//    In the most common case, where a QR factorization is being computed,
//    ROW = COL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//
//    Input, float A[N*N], the matrix whose Householder matrix
//    is to be computed.
//
//    Input, int ROW, COL, specify the location of the
//    entry of the matrix A which is to be preserved.  The entries in
//    the same column, but higher rows, will be zeroed out if A is
//    premultiplied by H.
//
//    Output, float R4MAT_HOUSE_PRE[N*N], the Householder matrix.
//
{
  float *h;
  int i;
  float *v;
  float *w;

  w = new float[n];
//
//  Set up the vector V.
//
  for ( i = 0; i < row-1; i++ )
  {
    w[i] = 0.0;
  }
  for ( i = row-1; i < n; i++ )
  {
    w[i] = a[i+col*n];
  }

  v = r4vec_house_column ( n, w, row );
//
//  Form the matrix H(V).
//
  h = r4mat_house_form ( n, v );

  delete [] v;
  delete [] w;

  return h;
}
//****************************************************************************80

float *r4mat_identity ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_IDENTITY sets the square matrix A to the identity.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Output, float A[N*N], the N by N identity matrix.
//
{
  float *a;
  int i;
  int j;
  int k;

  a = new float[n*n];

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return a;
}
//****************************************************************************80

bool r4mat_in_01 ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_IN_01 is TRUE if the entries of an R4MAT are in the range [0,1].
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the matrix.
//
//    Output, bool R4MAT_IN_01, is TRUE if every entry of A is
//    between 0 and 1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < 0.0 || 1.0 < a[i+j*m] )
      {
        return false;
      }
    }
  }

  return true;
}
//****************************************************************************80

float *r4mat_indicator_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_INDICATOR_NEW sets up an "indicator" R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The value of each entry suggests its location, as in:
//
//      11  12  13  14
//      21  22  23  24
//      31  32  33  34
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Output, float R4MAT_INDICATOR_NEW[M*N], the table.
//
{
  float *a;
  int fac;
  int i;
  int j;

  a = new float[m*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = ( float ) ( fac * i + j );
    }
  }
  return a;
}
//****************************************************************************80

float *r4mat_inverse_2d ( float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_INVERSE_2D inverts a 2 by 2 matrix using Cramer's rule.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[2*2], the matrix to be inverted.
//
//    Output, float R4MAT_INVERSE_2D[2*2], the inverse of the matrix A.
//
{
  float *b;
  float det;
//
//  Compute the determinant of A.
//
  det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];
//
//  If the determinant is zero, bail out.
//
  if ( det == 0.0 )
  {
    return NULL;
  }
//
//  Compute the entries of the inverse matrix using an explicit formula.
//
  b = new float[2*2];

  b[0+0*2] = + a[1+1*2] / det;
  b[0+1*2] = - a[0+1*2] / det;
  b[1+0*2] = - a[1+0*2] / det;
  b[1+1*2] = + a[0+0*2] / det;

  return b;
}
//****************************************************************************80

float *r4mat_inverse_3d ( float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_INVERSE_3D inverts a 3 by 3 matrix using Cramer's rule.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    If the determinant is zero, A is singular, and does not have an
//    inverse.  In that case, the output is set to NULL.
//
//    If the determinant is nonzero, its value is an estimate
//    of how nonsingular the matrix A is.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[3*3], the matrix to be inverted.
//
//    Output, float R4MAT_INVERSE_3D[3*3], the inverse of the matrix A.
//
{
  float *b;
  float det;
//
//  Compute the determinant of A.
//
  det =
     a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
   + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
   + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );

  if ( det == 0.0 )
  {
    return NULL;
  }

  b = new float[3*3];

  b[0+0*3] =   ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) / det;
  b[0+1*3] = - ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] ) / det;
  b[0+2*3] =   ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] ) / det;

  b[1+0*3] = - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] ) / det;
  b[1+1*3] =   ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] ) / det;
  b[1+2*3] = - ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] ) / det;

  b[2+0*3] =   ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] ) / det;
  b[2+1*3] = - ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] ) / det;
  b[2+2*3] =   ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] ) / det;

  return b;
}
//****************************************************************************80

float *r4mat_inverse_4d ( float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_INVERSE_4D inverts a 4 by 4 matrix using Cramer's rule.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[4][4], the matrix to be inverted.
//
//    Output, float R4MAT_INVERSE_4D[4][4], the inverse of the matrix A.
//
{
  float *b;
  float det;
//
//  Compute the determinant of A.
//
  det = r4mat_det_4d ( a );
//
//  If the determinant is zero, bail out.
//
  if ( det == 0.0 )
  {
    return NULL;
  }
//
//  Compute the entries of the inverse matrix using an explicit formula.
//
  b = new float[4*4];

  b[0+0*4] =
    +(
    + a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+3*4] * a[3+1*4] - a[2+1*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    ) / det;

  b[1+0*4] =
    -(
    + a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
    ) / det;

  b[2+0*4] =
    +(
    + a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
    + a[1+1*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[3+0*4] =
    -(
    + a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    + a[1+1*4] * ( a[2+2*4] * a[3+0*4] - a[2+0*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[0+1*4] =
    -(
    + a[0+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+3*4] * a[3+1*4] - a[2+1*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    ) / det;

  b[1+1*4] =
    +(
    + a[0+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
    ) / det;

  b[2+1*4] =
    -(
    + a[0+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
    + a[0+1*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[3+1*4] =
    +(
    + a[0+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    + a[0+1*4] * ( a[2+2*4] * a[3+0*4] - a[2+0*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[0+2*4] =
    +(
    + a[0+1*4] * ( a[1+2*4] * a[3+3*4] - a[1+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[3+1*4] - a[1+1*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+1*4] * a[3+2*4] - a[1+2*4] * a[3+1*4] )
    ) / det;

  b[1+2*4] =
    -(
    + a[0+0*4] * ( a[1+2*4] * a[3+3*4] - a[1+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[3+0*4] - a[1+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[3+2*4] - a[1+2*4] * a[3+0*4] )
    ) / det;

  b[2+2*4] =
    +(
    + a[0+0*4] * ( a[1+1*4] * a[3+3*4] - a[1+3*4] * a[3+1*4] )
    + a[0+1*4] * ( a[1+3*4] * a[3+0*4] - a[1+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[3+1*4] - a[1+1*4] * a[3+0*4] )
    ) / det;

  b[3+2*4] =
    -(
    + a[0+0*4] * ( a[1+1*4] * a[3+2*4] - a[1+2*4] * a[3+1*4] )
    + a[0+1*4] * ( a[1+2*4] * a[3+0*4] - a[1+0*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+0*4] * a[3+1*4] - a[1+1*4] * a[3+0*4] )
    ) / det;

  b[0+3*4] =
    -(
    + a[0+1*4] * ( a[1+2*4] * a[2+3*4] - a[1+3*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[2+1*4] - a[1+1*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+1*4] * a[2+2*4] - a[1+2*4] * a[2+1*4] )
    ) / det;

  b[1+3*4] =
    +(
    + a[0+0*4] * ( a[1+2*4] * a[2+3*4] - a[1+3*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[2+0*4] - a[1+0*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[2+2*4] - a[1+2*4] * a[2+0*4] )
    ) / det;

  b[2+3*4] =
    -(
    + a[0+0*4] * ( a[1+1*4] * a[2+3*4] - a[1+3*4] * a[2+1*4] )
    + a[0+1*4] * ( a[1+3*4] * a[2+0*4] - a[1+0*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[2+1*4] - a[1+1*4] * a[2+0*4] )
    ) / det;

  b[3+3*4] =
    +(
    + a[0+0*4] * ( a[1+1*4] * a[2+2*4] - a[1+2*4] * a[2+1*4] )
    + a[0+1*4] * ( a[1+2*4] * a[2+0*4] - a[1+0*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+0*4] * a[2+1*4] - a[1+1*4] * a[2+0*4] )
    ) / det;

  return b;
}
//****************************************************************************80

float *r4mat_jac ( int m, int n, float eps,
  float *(*fx) ( int m, int n, float x[] ), float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_JAC estimates a dense jacobian matrix of the function FX.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    FPRIME(I,J) = d F(I) / d X(J).
//
//    The jacobian is assumed to be dense, and the LINPACK/LAPACK
//    float precision general matrix storage mode ("DGE") is used.
//
//    Forward differences are used, requiring N+1 function evaluations.
//
//    Values of EPS have typically been chosen between
//    sqrt ( EPSMCH ) and sqrt ( sqrt ( EPSMCH ) ) where EPSMCH is the
//    machine tolerance.
//
//    If EPS is too small, then F(X+EPS) will be the same as
//    F(X), and the jacobian will be full of zero entries.
//
//    If EPS is too large, the finite difference estimate will
//    be inaccurate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of functions.
//
//    Input, int N, the number of variables.
//
//    Input, float EPS, a tolerance to be used for shifting the
//    X values during the finite differencing.  No single value
//    of EPS will be reliable for all vectors X and functions FX.
//
//    Input, float *(*FX) ( int m, int n, float x[] ), the name of
//    the user written routine which evaluates the M-dimensional
//    function at a given N-dimensional point X.
//
//    Input, float X[N], the point where the jacobian
//    is to be estimated.
//
//    Output, float R4MAT_JAC[M*N], the estimated jacobian matrix.
//
{
  float del;
  float *fprime;
  int i;
  int j;
  float xsave;
  float *work1;
  float *work2;

  fprime = new float[m*n];
//
//  Evaluate the function at the base point, X.
//
  work2 = fx ( m, n, x );
//
//  Now, one by one, vary each component J of the base point X, and
//  estimate DF(I)/DX(J) = ( F(X+) - F(X) )/ DEL.
//
  for ( j = 0; j < n; j++ )
  {
    xsave = x[j];
    del = eps * ( 1.0 + r4_abs ( x[j] ) );
    x[j] = x[j] + del;
    work1 = fx ( m, n, x );
    x[j] = xsave;
    for ( i = 0; i < m; i++ )
    {
      fprime[i+j*m] = ( work1[i] - work2[i] ) / del;
    }
    delete [] work1;
  }
  delete [] work2;

  return fprime;
}
//****************************************************************************80

float *r4mat_l_inverse ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_L_INVERSE inverts a lower triangular R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    A lower triangular matrix is a matrix whose only nonzero entries
//    occur on or below the diagonal.
//
//    The inverse of a lower triangular matrix is a lower triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, number of rows and columns in the matrix.
//
//    Input, float A[N*N], the lower triangular matrix.
//
//    Output, float R4MAT_L_INVERSE[N*N], the inverse matrix.
//
{
  float *b;
  int i;
  int j;
  int k;
  float temp;

  b = new float[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i < j )
      {
        b[i+j*n] = 0.0;
      }
      else if ( j == i )
      {
        b[i+j*n] = 1.0 / a[i+j*n];
      }
      else
      {
        temp = 0.0;
        for ( k = 0; k < i; k++ )
        {
          temp = temp + a[i+k*n] * b[k+j*n];
        }
        b[i+j*n] = -temp / a[i+i*n];
      }
    }
  }

  return b;
}
//****************************************************************************80

void r4mat_l_print ( int m, int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_L_PRINT prints a lower triangular R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Example:
//
//    M = 5, N = 5
//    A = (/ 11, 21, 31, 41, 51, 22, 32, 42, 52, 33, 43, 53, 44, 54, 55 /)
//
//    11
//    21 22
//    31 32 33
//    41 42 43 44
//    51 52 53 54 55
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[*], the M by N matrix.  Only the lower
//    triangular elements are stored, in column major order.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int indx[10];
  int j;
  int jhi;
  int jlo;
  int jmax;
  int nn;
  int size;

  cout << "\n";
  cout << title << "\n";

  jmax = i4_min ( n, m );

  if ( m <= n )
  {
    size = ( m * ( m + 1 ) ) / 2;
  }
  else if ( n < m )
  {
    size = ( n * ( n + 1 ) ) / 2 + ( m - n ) * n;
  }

  if ( r4vec_is_int ( size, a ) )
  {
    nn = 10;
    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m, jmax ) );
      cout << "\n";
      cout << "  Col   ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(6) << j;
      }
      cout << "\n";
      cout << "  Row  \n";
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        cout << "  " << setw(6) << i;
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          cout << setw(6) << a[indx[j]-1];
        }
        cout << "\n";
      }
    }
  }
  else if ( r4vec_amax ( size, a ) < 1000000.0 )
  {
    nn = 5;
    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn );
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
      cout << "\n";
      cout << "  Col ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(14) << j;
      }
      cout << "\n";
      cout << "  Row  \n";
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        cout << "  " << setw(6) << i;
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          cout << setw(14) << a[indx[j]-1];
        }
        cout << "\n";
      }
    }
  }
  else
  {
    nn = 5;

    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
      cout << "\n";
      cout << "  Col ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(7) << j << "       ";
      }
      cout << "\n";
      cout << "  Row \n";
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        cout << setw(6) << i;
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          cout << setw(14) << a[indx[j]-1];
        }
      }
    }
  }

  return;
}
//****************************************************************************80

float *r4mat_l_solve ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_L_SOLVE solves a lower triangular linear system.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix A.
//
//    Input, float A[N*N], the N by N lower triangular matrix.
//
//    Input, float B[N], the right hand side of the linear system.
//
//    Output, float R4MAT_L_SOLVE[N], the solution of the linear system.
//
{
  int i;
  int j;
  float temp;
  float *x;

  x = new float[n];
//
//  Solve L * x = b.
//
  for ( i = 0; i < n; i++ )
  {
    temp = 0.0;
    for ( j = 0; j < i; j++ )
    {
      temp = temp + a[i+j*n] * x[j];
    }
    x[i] = ( b[i] - temp ) / a[i+i*n];
  }

  return x;
}
//****************************************************************************80

float *r4mat_l1_inverse ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_L1_INVERSE inverts a unit lower triangular R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    A unit lower triangular matrix is a matrix with only 1's on the main
//    diagonal, and only 0's above the main diagonal.
//
//    The inverse of a unit lower triangular matrix is also
//    a unit lower triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, number of rows and columns in the matrix.
//
//    Input, float A[N*N], the unit lower triangular matrix.
//
//    Output, float R4MAT_L1_INVERSE[N*N], the inverse matrix.
//
{
  float *b;
  int i;
  int j;
  int k;

  b = new float[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i < j )
      {
        b[i+j*n] = 0.0;
      }
      else if ( j == i )
      {
        b[i+j*n] = 1.0;
      }
      else
      {
        b[i+j*n] = 0.0;
        for ( k = 0; k < i; k++ )
        {
          b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
        }
      }
    }
  }

  return b;
}
//****************************************************************************80

float *r4mat_lt_solve ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_LT_SOLVE solves a transposed lower triangular linear system.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Given the lower triangular matrix A, the linear system to be solved is:
//
//      A' * x = b
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix A.
//
//    Input, float A[N*N], the N by N lower triangular matrix.
//
//    Input, float B[N], the right hand side of the linear system.
//
//    Output, float R4MAT_LT_SOLVE[N], the solution of the linear system.
//
{
  int i;
  int j;
  float *x;

  x = new float[n];

  for ( j = n-1; 0 <= j; j-- )
  {
    x[j] = b[j];
    for ( i = j+1; i < n; i++ )
    {
      x[j] = x[j] - x[i] * a[i+j*n];
    }
    x[j] = x[j] / a[j+j*n];
  }

  return x;
}
//****************************************************************************80

void r4mat_lu ( int m, int n, float a[], float l[], float p[], float u[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_LU computes the LU factorization of a rectangular R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The routine is given an M by N matrix A, and produces
//
//      L, an M by M unit lower triangular matrix,
//      U, an M by N upper triangular matrix, and
//      P, an M by M permutation matrix P,
//
//    so that
//
//      A = P' * L * U.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the M by N matrix to be factored.
//
//    Output, float L[M*M], the M by M unit lower triangular factor.
//
//    Output, float P[M*M], the M by M permutation matrix.
//
//    Output, float U[M*N], the M by N upper triangular factor.
//
{
  int i;
  int ipiv;
  int j;
  int k;
  float pivot;
//
//  Initialize:
//
//    U:=A
//    L:=Identity
//    P:=Identity
//
  r4mat_copy ( m, n, a, u );

  r4mat_zero ( m, m, l );
  r4mat_zero ( m, m, p );
  for ( i = 0; i < m; i++ )
  {
    l[i+i*m] = 1.0;
    p[i+i*m] = 1.0;
  }
//
//  On step J, find the pivot row, IPIV, and the pivot value PIVOT.
//
  for ( j = 0; j < i4_min ( m - 1, n ); j++ )
  {
    pivot = 0.0;
    ipiv = -1;

    for ( i = j; i < m; i++ )
    {
      if ( pivot < r4_abs ( u[i+j*m] ) )
      {
        pivot = r4_abs ( u[i+j*m] );
        ipiv = i;
      }
    }
//
//  Unless IPIV is zero, swap rows J and IPIV.
//
    if ( ipiv != -1 )
    {
      r4row_swap ( m, n, u, j+1, ipiv+1 );

      r4row_swap ( m, m, l, j+1, ipiv+1 );

      r4row_swap ( m, m, p, j+1, ipiv+1 );
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

float r4mat_max ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MAX returns the maximum entry of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
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
//    Input, float A[M*N], the M by N matrix.
//
//    Output, float R4MAT_MAX, the maximum entry of A.
//
{
  int i;
  int j;
  float value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

void r4mat_max_index ( int m, int n, float a[], int *i_max, int *j_max )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MAX_INDEX returns the location of the maximum entry of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the M by N matrix.
//
//    Output, int *I_MAX, *J_MAX, the indices of the maximum entry of A.
//
{
  int i;
  int i2;
  int j;
  int j2;

  i2 = -1;
  j2 = -1;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i2 == -1 && j2 == -1 )
      {
        i2 = i;
        j2 = j;
      }
      else if ( a[i2+j2*m] < a[i+j*m] )
      {
        i2 = i;
        j2 = j;
      }
    }
  }

  *i_max = i2 + 1;
  *j_max = j2 + 1;

  return;
}
//****************************************************************************80

float r4mat_maxcol_minrow ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MAXCOL_MINROW gets the maximum column minimum row of an M by N matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    R4MAT_MAXCOL_MINROW = max ( 1 <= I <= N ) ( min ( 1 <= J <= M ) A(I,J) )
//
//    For a given matrix, R4MAT_MAXCOL_MINROW <= R4MAT_MINROW_MAXCOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the matrix.
//
//    Output, float R4MAT_MAXCOL_MINROW, the maximum column
//    minimum row entry of A.
//
{
  int i;
  int j;
  float minrow;
  float value;

  value = - r4_huge ( );

  for ( i = 0; i < m; i++ )
  {
    minrow = r4_huge ( );

    for ( j = 0; j < n; j++ )
    {
      minrow = r4_min ( minrow, a[i+j*m] );
    }
    value = r4_max ( value, minrow );
  }

  return value;
}
//****************************************************************************80

float r4mat_maxrow_mincol ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MAXROW_MINCOL gets the maximum row minimum column of an M by N matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    R4MAT_MAXROW_MINCOL = max ( 1 <= J <= N ) ( min ( 1 <= I <= M ) A(I,J) )
//
//    For a given matrix, R4MAT_MAXROW_MINCOL <= R4MAT_MINCOL_MAXROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the matrix.
//
//    Output, float R4MAT_MAXROW_MINCOL, the maximum row
//    minimum column entry of A.
//
{
  int i;
  int j;
  float mincol;
  float value;

  value = - r4_huge ( );

  for ( j = 0; j < n; j++ )
  {
    mincol = r4_huge ( );
    for ( i = 0; i < m; i++ )
    {
      mincol = r4_min ( mincol, a[i+j*m] );
    }
    value = r4_max ( value, mincol );
  }
  return value;
}
//****************************************************************************80

float r4mat_min ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MIN returns the minimum entry of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
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
//    Input, float A[M*N], the M by N matrix.
//
//    Output, float DMIN_MAX, the minimum entry of A.
//
{
  int i;
  int j;
  float value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < value )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

void r4mat_min_index ( int m, int n, float a[], int *i_min, int *j_min )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MIN_INDEX returns the location of the minimum entry of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the M by N matrix.
//
//    Output, int *I_MIN, *J_MIN, the indices of the minimum entry of A.
//
{
  int i;
  int i2;
  int j;
  int j2;

  i2 = -1;
  j2 = -1;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i2 == -1 && j2 == -1 )
      {
        i2 = i;
        j2 = j;
      }
      else if ( a[i+j*m] < a[i2+j2*m] )
      {
        i2 = i;
        j2 = j;
      }
    }
  }

  *i_min = i2 + 1;
  *j_min = j2 + 1;

  return;
}
//****************************************************************************80

float r4mat_mincol_maxrow ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MINCOL_MAXROW gets the minimum column maximum row of an M by N matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    R4MAT_MINCOL_MAXROW = min ( 1 <= I <= N ) ( max ( 1 <= J <= M ) A(I,J) )
//
//    For a given matrix, R4MAT_MAXROW_MINCOL <= R4MAT_MINCOL_MAXROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A(M,N), the matrix.
//
//    Output, float R4MAT_MINCOL_MAXROW, the minimum column
//    maximum row entry of A.
//
{
  int i;
  int j;
  float maxrow;
  float value;

  value = r4_huge ( );

  for ( i = 0; i < m; i++ )
  {
    maxrow = - r4_huge ( );
    for ( j = 0; j < n; j++ )
    {
      maxrow = r4_max ( maxrow, a[i+j*m] );
    }
    value = r4_min ( value, maxrow );
  }

  return value;
}
//****************************************************************************80

float r4mat_minrow_maxcol ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MINROW_MAXCOL gets the minimum row maximum column of an M by N matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    R4MAT_MINROW_MAXCOL = min ( 1 <= J <= N ) ( max ( 1 <= I <= M ) A(I,J) )
//
//    For a given matrix, R4MAT_MAXCOL_MINROW <= R4MAT_MINROW_MAXCOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the matrix.
//
//    Output, float R4MAT_MINROW_MAXCOL, the minimum row
//    maximum column entry of A.
//
{
  int i;
  int j;
  float maxcol;
  float value;;

  value = r4_huge ( );

  for ( j = 0; j < n; j++ )
  {
    maxcol = - r4_huge ( );
    for ( i = 0; i < m; i++ )
    {
      maxcol = r4_max ( maxcol, a[i+j*m] );
    }
    value = r4_min ( value, maxcol );
  }

  return value;
}
//****************************************************************************80

void r4mat_mm ( int n1, int n2, int n3, float a[], float b[], float c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MM multiplies two matrices.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, float A[N1*N2], float B[N2*N3], the matrices to multiply.
//
//    Output, float C[N1*N3], the product matrix C = A * B.
//
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return;
}
//****************************************************************************80

float *r4mat_mm_new ( int n1, int n2, int n3, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, float A[N1*N2], float B[N2*N3], the matrices to multiply.
//
//    Output, float R4MAT_MM[N1*N3], the product matrix C = A * B.
//
{
  float *c;
  int i;
  int j;
  int k;

  c = new float[n1*n3];

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
//****************************************************************************80

float *r4mat_mtv ( int m, int n, float a[], float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MTV multiplies a transposed matrix times a vector.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, float A[M,N], the M by N matrix.
//
//    Input, float X[M], the vector to be multiplied by A.
//
//    Output, float R4MAT_MTV[N], the product A'*X.
//
{
  int i;
  int j;
  float *y;

  y = new float[n];

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  return y;
}
//****************************************************************************80

void r4mat_mtxv ( int m, int n, float a[], float x[], float atx[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MTXV multiplies a transposed matrix times a vector.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as an argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, float A[M,N], the M by N matrix.
//
//    Input, float X[M], the vector to be multiplied by A.
//
//    Output, float ATX[N], the product A'*X.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    atx[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      atx[j] = atx[j] + a[i+j*m] * x[i];
    }
  }

  return;
}
//****************************************************************************80

float *r4mat_mv ( int m, int n, float a[], float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MV multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, float A[M,N], the M by N matrix.
//
//    Input, float X[N], the vector to be multiplied by A.
//
//    Output, float R4MAT_MV[M], the product A*X.
//
{
  int i;
  int j;
  float *y;

  y = new float[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
//****************************************************************************80

void r4mat_mxm ( int n1, int n2, int n3, float a[], float b[], float c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MXM multiplies two matrices.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as an argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, float A[N1*N2], float B[N2*N3], the matrices to multiply.
//
//    Output, float C[N1*N3], the product matrix C = A * B.
//
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return;
}
//****************************************************************************80

float *r4mat_mxm_new ( int n1, int n2, int n3, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MXM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, float A[N1*N2], float B[N2*N3], the matrices to multiply.
//
//    Output, float R4MAT_MXM_NEW[N1*N3], the product matrix C = A * B.
//
{
  float *c;
  int i;
  int j;
  int k;

  c = new float[n1*n3];

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
//****************************************************************************80

void r4mat_mxv ( int m, int n, float a[], float x[], float ax[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MXV multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as an argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, float A[M,N], the M by N matrix.
//
//    Input, float X[N], the vector to be multiplied by A.
//
//    Output, float AX[M], the product A*X.
//
{
  int i;
  int j;

  for ( i = 0; i < m; i++ )
  {
    ax[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      ax[i] = ax[i] + a[i+j*m] * x[j];
    }
  }

  return;
}
//****************************************************************************80

float **r4mat_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_NEW allocates a new R4MAT.
//
//  Discussion:
//
//    A declaration of the form
//      float **a;
//    is necesary.  Then an assignment of the form:
//      a = r4mat_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17.0;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, float R4MAT_NEW[M][N], a new matrix.
//
{
  float **a;
  int i;

  a = new float *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R4MAT_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new float[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R4MAT_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}
//****************************************************************************80

void r4mat_nint ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_NINT rounds the entries of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of A.
//
//    Input/output, float A[M*N], the matrix to be NINT'ed.
//
{
  int i;
  int j;
  int s;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*m] < 0.0 )
      {
        s = -1;
      }
      else
      {
        s = 1;
      }
      a[i+j*m] = s * ( int ) ( r4_abs ( a[i+j*m] ) + 0.5 );
    }
  }

  return;
}
//****************************************************************************80

float r4mat_norm_eis ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_NORM_EIS returns the EISPACK norm of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The EISPACK norm is defined as:
//
//      R4MAT_NORM_EIS =
//        sum ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the matrix whose EISPACK norm is desired.
//
//    Output, float R4MAT_NORM_EIS, the EISPACK norm of A.
//
{
  int i;
  int j;
  float value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + r4_abs ( a[i+j*m] );
    }
  }

  return value;
}
//****************************************************************************80

float r4mat_norm_fro ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_NORM_FRO returns the Frobenius norm of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R4MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)**2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r4vec_norm_l2 ( A * x ) <= r4mat_norm_fro ( A ) * r4vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the matrix whose Frobenius
//    norm is desired.
//
//    Output, float R4MAT_NORM_FRO, the Frobenius norm of A.
//
{
  int i;
  int j;
  float value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

float *r4mat_nullspace ( int m, int n, float a[], int nullspace_size )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_NULLSPACE computes the nullspace of a matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Let A be an MxN matrix.
//
//    If X is an N-vector, and A*X = 0, then X is a null vector of A.
//
//    The set of all null vectors of A is called the nullspace of A.
//
//    The 0 vector is always in the null space.
//
//    If the 0 vector is the only vector in the nullspace of A, then A
//    is said to have maximum column rank.  (Because A*X=0 can be regarded
//    as a linear combination of the columns of A).  In particular, if A
//    is square, and has maximum column rank, it is nonsingular.
//
//    The dimension of the nullspace is the number of linearly independent
//    vectors that span the nullspace.  If A has maximum column rank,
//    its nullspace has dimension 0.
//
//    This routine uses the reduced row echelon form of A to determine
//    a set of NULLSPACE_SIZE independent null vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of
//    the matrix A.
//
//    Input, float A[M*N], the matrix to be analyzed.
//
//    Input, int NULLSPACE_SIZE, the size of the nullspace.
//
//    Output, float R4MAT_NULLSPACE[N*NULLSPACE_SIZE], vectors that
//    span the nullspace.
//
{
  int *col;
  int i;
  int i2;
  int j;
  int j2;
  float *nullspace;
  int *row;
  float *rref;
//
//  Make a copy of A.
//
  rref = r4mat_copy_new ( m, n, a );
//
//  Get the reduced row echelon form of A.
//
  r4mat_rref ( m, n, rref );
//
//  Note in ROW the columns of the leading nonzeros.
//  COL(J) = +J if there is a leading 1 in that column, and -J otherwise.
//
  row = new int[m];
  for ( i = 0; i < m; i++ )
  {
    row[i] = 0;
  }

  col = new int[n];
  for ( j = 0; j < n; j++ )
  {
    col[j] = - ( j + 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( rref[i+j*m] == 1.0 )
      {
        row[i] = ( j + 1 );
        col[j] = ( j + 1 );
        break;
      }
    }
  }

  nullspace = r4mat_zero_new ( n, nullspace_size );

  j2 = 0;
//
//  If column J does not contain a leading 1, then it contains
//  information about a null vector.
//
  for ( j = 0; j < n; j++ )
  {
    if ( col[j] < 0 )
    {
      for ( i = 0; i < m; i++ )
      {
        if ( rref[i+j*m] != 0.0 )
        {
          i2 = row[i] - 1;
          nullspace[i2+j2*n] = - rref[i+j*m];
        }
      }
      nullspace[j+j2*n] = 1.0;
      j2 = j2 + 1;
    }
  }
  delete [] col;
  delete [] row;
  delete [] rref;

  return nullspace;
}
//****************************************************************************80

int r4mat_nullspace_size ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Let A be an MxN matrix.
//
//    If X is an N-vector, and A*X = 0, then X is a null vector of A.
//
//    The set of all null vectors of A is called the nullspace of A.
//
//    The 0 vector is always in the null space.
//
//    If the 0 vector is the only vector in the nullspace of A, then A
//    is said to have maximum column rank.  (Because A*X=0 can be regarded
//    as a linear combination of the columns of A).  In particular, if A
//    is square, and has maximum column rank, it is nonsingular.
//
//    The dimension of the nullspace is the number of linearly independent
//    vectors that span the nullspace.  If A has maximum column rank,
//    its nullspace has dimension 0.
//
//    This routine ESTIMATES the dimension of the nullspace.  Cases of
//    singularity that depend on exact arithmetic will probably be missed.
//
//    The nullspace will be estimated by counting the leading 1's in the
//    reduced row echelon form of A, and subtracting this from N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of
//    the matrix A.
//
//    Input, float A[M*N], the matrix to be analyzed.
//
//    Output, int R4MAT_NULLSPACE_SIZE, the estimated size
//    of the nullspace.
//
{
  int i;
  int j;
  int leading;
  int nullspace_size;
  float *rref;
//
//  Make a copy of A.
//
  rref = r4mat_copy_new ( m, n, a );
//
//  Get the reduced row echelon form of A.
//
  r4mat_rref ( m, n, rref );
//
//  Count the leading 1's in A.
//
  leading = 0;
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( rref[i+j*m] == 1.0 )
      {
        leading = leading + 1;
        break;
      }
    }
  }
  nullspace_size = n - leading;

  return nullspace_size;
}
//****************************************************************************80

float *r4mat_orth_uniform_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_ORTH_UNIFORM_NEW returns a random orthogonal matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The inverse of A is equal to A'.
//
//    A * A'  = A' * A = I.
//
//    Columns and rows of A have unit Euclidean norm.
//
//    Distinct pairs of columns of A are orthogonal.
//
//    Distinct pairs of rows of A are orthogonal.
//
//    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
//
//    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
//
//    The determinant of A is +1 or -1.
//
//    All the eigenvalues of A have modulus 1.
//
//    All singular values of A are 1.
//
//    All entries of A are between -1 and 1.
//
//  Discussion:
//
//    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
//    National Academy of Sciences of Belarus, for convincingly
//    pointing out the severe deficiencies of an earlier version of
//    this routine.
//
//    Essentially, the computation involves saving the Q factor of the
//    QR factorization of a matrix whose entries are normally distributed.
//    However, it is only necessary to generate this matrix a column at
//    a time, since it can be shown that when it comes time to annihilate
//    the subdiagonal elements of column K, these (transformed) elements of
//    column K are still normally distributed random values.  Hence, there
//    is no need to generate them at the beginning of the process and
//    transform them K-1 times.
//
//    For computational efficiency, the individual Householder transformations
//    could be saved, as recommended in the reference, instead of being
//    accumulated into an explicit matrix format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pete Stewart,
//    Efficient Generation of Random Orthogonal Matrices With an Application
//    to Condition Estimators,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 3, June 1980, pages 403-409.
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, float R4MAT_ORTH_UNIFORM_NEW[N*N], the orthogonal matrix.
//
{
  float *a;
  float *a2;
  int i;
  int j;
  float *v;
  float *x;
//
//  Start with A = the identity matrix.
//
  a = r4mat_identity ( n );
//
//  Now behave as though we were computing the QR factorization of
//  some other random matrix.  Generate the N elements of the first column,
//  compute the Householder matrix H1 that annihilates the subdiagonal elements,
//  and set A := A * H1' = A * H.
//
//  On the second step, generate the lower N-1 elements of the second column,
//  compute the Householder matrix H2 that annihilates them,
//  and set A := A * H2' = A * H2 = H1 * H2.
//
//  On the N-1 step, generate the lower 2 elements of column N-1,
//  compute the Householder matrix HN-1 that annihilates them, and
//  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
//  This is our random orthogonal matrix.
//
  x = new float[n];

  for ( j = 1; j < n; j++ )
  {
//
//  Set the vector that represents the J-th column to be annihilated.
//
    for ( i = 1; i < j; i++ )
    {
      x[i-1] = 0.0;
    }
    for ( i = j; i <= n; i++ )
    {
      x[i-1] = r4_normal_01 ( seed );
    }
//
//  Compute the vector V that defines a Householder transformation matrix
//  H(V) that annihilates the subdiagonal elements of X.
//
    v = r4vec_house_column ( n, x, j );
//
//  Postmultiply the matrix A by H'(V) = H(V).
//
    a2 = r4mat_house_axh_new ( n, a, v );

    delete [] v;

    r4mat_copy ( n, n, a2, a );

    delete [] a2;
  }

  delete [] x;

  return a;
}
//****************************************************************************80

void r4mat_plot ( int m, int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_PLOT "plots" an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
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
//    Input, float A[M*N], the matrix.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int j;
  int jhi;
  int jlo;

  cout << "\n";
  cout << title << "\n";

  for ( jlo = 1; jlo <= n; jlo = jlo + 70 )
  {
    jhi = i4_min ( jlo + 70-1, n );
    cout << "\n";
    cout << "          ";
    for ( j = jlo; j <= jhi; j++ )
    {
      cout <<  ( j % 10 );
    }
    cout << "\n";
    cout << "\n";

    for ( i = 1; i <= m; i++ )
    {
      cout << setw(6) << i << "    ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << r4mat_plot_symbol ( a[i-1+(j-1)*m] );
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

char r4mat_plot_symbol ( float r )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_PLOT_SYMBOL returns a symbol for entries of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R, a value whose symbol is desired.
//
//    Output, char R4MAT_PLOT_SYMBOL, is
//    '-' if R is negative,
//    '0' if R is zero,
//    '+' if R is positive.
//
{
  char c;

  if ( r < 0.0 )
  {
    c = '-';
  }
  else if ( r == 0.0 )
  {
    c = '0';
  }
  else if ( 0.0 < r )
  {
    c = '+';
  }

  return c;
}
//****************************************************************************80

float *r4mat_poly_char ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_POLY_CHAR computes the characteristic polynomial of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, float A[N*N], the N by N matrix.
//
//    Output, float R4MAT_POLY_CHAR[N+1], the coefficients of the characteristic
//    polynomial of A.  P(N) contains the coefficient of X**N
//    (which will be 1), P(I) contains the coefficient of X**I,
//    and P(0) contains the constant term.
//
{
  int i;
  int order;
  float *p;
  float trace;
  float *work1;
  float *work2;

  p = new float[n+1];
//
//  Initialize WORK1 to the identity matrix.
//
  work1 = r4mat_identity ( n );

  p[n] = 1.0;

  for ( order = n-1; 0 <= order; order-- )
  {
//
//  Work2 = A * WORK1.
//
    work2 = r4mat_mm_new ( n, n, n, a, work1 );
//
//  Take the trace.
//
    trace = r4mat_trace ( n, work2 );
//
//  P(ORDER) = -Trace ( WORK2 ) / ( N - ORDER )
//
    p[order] = -trace / ( float ) ( n - order );
//
//  WORK1 := WORK2 + P(IORDER) * Identity.
//
    delete [] work1;

    r4mat_copy ( n, n, work2, work1 );

    delete [] work2;

    for ( i = 0; i < n; i++ )
    {
      work1[i+i*n] = work1[i+i*n] + p[order];
    }
  }

  delete [] work1;

  return p;
}
//****************************************************************************80

float *r4mat_power ( int n, float a[], int npow )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_POWER computes a nonnegative power of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The algorithm is:
//
//      B = I
//      do NPOW times:
//        B = A * B
//      end
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Input, float A[N*N], the matrix to be raised to a power.
//
//    Input, int NPOW, the power to which A is to be raised.
//    NPOW must be nonnegative.
//
//    Output, float B[N*N], the value of A**NPOW.
//
{
  float *b;
  float *c;
  int ipow;

  if ( npow < 0 )
  {
    cerr << "\n";
    cerr << "R4MAT_POWER - Fatal error!\n";
    cerr << "  Input value of NPOW < 0.\n";
    cerr << "  NPOW = " << npow << "\n";
    exit ( 1 );
  }

  b = r4mat_identity ( n );

  for ( ipow = 1; ipow <= npow; ipow++ )
  {
    c = r4mat_mm_new ( n, n, n, a, b );
    r4mat_copy ( n, n, c, b );
    delete [] c;
  }

  return b;
}
//****************************************************************************80

void r4mat_power_method ( int n, float a[], float *r, float v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_POWER_METHOD applies the power method to a matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    If the power method has not converged, then calling the routine
//    again immediately with the output from the previous call will
//    continue the iteration.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Input, float A[N*N], the matrix.
//
//    Output, float *R, the estimated eigenvalue.
//
//    Input/output, float V[N], on input, an estimate
//    for the eigenvector.  On output, an improved estimate for the
//    eigenvector.
//
{
  float *av;
  float eps;
  int i;
  int it;
  float it_eps = 0.0001;
  int it_max = 100;
  int it_min = 10;
  int j;
  float r2;
  float r_old;

  eps = sqrt ( r4_epsilon ( ) );

  *r = r4vec_norm ( n, v );

  if ( *r == 0.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v[i] = 1.0;
    }
    *r = sqrt ( ( float ) n );
  }

  for ( i = 0; i < n; i++ )
  {
    v[i] = v[i] / *r;
  }

  for ( it = 1; it <= it_max; it++ )
  {
    av = r4mat_mv ( n, n, a, v );

    r_old = *r;
    *r = r4vec_norm ( n, av );

    if ( it_min < it )
    {
      if ( r4_abs ( *r - r_old ) <= it_eps * ( 1.0 + r4_abs ( *r ) ) )
      {
        break;
      }
    }

    r4vec_copy ( n, av, v );

    delete [] av;

    if ( *r != 0.0 )
    {
      for ( i = 0; i < n; i++ )
      {
        v[i] = v[i] / *r;
      }
    }
//
//  Perturb V a bit, to avoid cases where the initial guess is exactly
//  the eigenvector of a smaller eigenvalue.
//
    if ( it < it_max / 2 )
    {
      j = ( ( it - 1 ) % n );
      v[j] = v[j] + eps * ( 1.0 + r4_abs ( v[j] ) );
      r2 = r4vec_norm ( n, v );
      for ( i = 0; i < n; i++ )
      {
        v[i] = v[i] / r2;
      }
    }
  }
  return;
}
//****************************************************************************80

void r4mat_print ( int m, int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_PRINT prints an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
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
//    Input, float A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_PRINT_SOME prints some of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
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
//    Input, float A[M*N], the matrix.
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
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

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
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }

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

void r4mat_ref ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_REF computes the row echelon form of a matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    A matrix is in row echelon form if:
//
//    * The first nonzero entry in each row is 1.
//
//    * The leading 1 in a given row occurs in a column to
//      the right of the leading 1 in the previous row.
//
//    * Rows which are entirely zero must occur last.
//
//  Example:
//
//    Input matrix:
//
//     1.0  3.0  0.0  2.0  6.0  3.0  1.0
//    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
//     3.0  9.0  0.0  0.0  6.0  6.0  2.0
//    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
//
//    Output matrix:
//
//     1.0  3.0  0.0  2.0  6.0  3.0  1.0
//     0.0  0.0  0.0  1.0  2.0  4.5  1.5
//     0.0  0.0  0.0  0.0  0.0  1.0  0.3
//     0.0  0.0  0.0  0.0  0.0  0.0  0.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of
//    the matrix A.
//
//    Input/output, float A[M*N].  On input, the matrix to be
//    analyzed.  On output, the REF form of the matrix.
//
{
  int i;
  int j;
  int lead;
  int r;
  float temp;

  lead = 0;

  for ( r = 0; r < m; r++ )
  {
    if ( n - 1 < lead )
    {
      break;
    }

    i = r;

    while ( a[i+lead*m] == 0.0 )
    {
      i = i + 1;

      if ( m - 1 < i )
      {
        i = r;
        lead = lead + 1;
        if ( n - 1 < lead )
        {
          lead = -1;
          break;
         }
      }
    }

    if ( lead < 0 )
    {
      break;
    }

    for ( j = 0; j < n; j++ )
    {
      temp     = a[i+j*m];
      a[i+j*m] = a[r+j*m];
      a[r+j*m] = temp;
    }

    temp = a[r+lead*m];

    for ( j = 0; j < n; j++ )
    {
      a[r+j*m] = a[r+j*m] / temp;
    }

    for ( i = r + 1; i < m; i++ )
    {
      temp = a[i+lead*m];
      for ( j = 0; j < n; j++ )
      {
        a[i+j*m] = a[i+j*m] - temp * a[r+j*m];
      }
    }
    lead = lead + 1;
  }
  return;
}
//****************************************************************************80

void r4mat_rref ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_RREF computes the reduced row echelon form of a matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    A matrix is in row echelon form if:
//
//    * The first nonzero entry in each row is 1.
//
//    * The leading 1 in a given row occurs in a column to
//      the right of the leading 1 in the previous row.
//
//    * Rows which are entirely zero must occur last.
//
//    The matrix is in reduced row echelon form if, in addition to
//    the first three conditions, it also satisfies:
//
//    * Each column containing a leading 1 has no other nonzero entries.
//
//  Example:
//
//    Input matrix:
//
//     1.0  3.0  0.0  2.0  6.0  3.0  1.0
//    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
//     3.0  9.0  0.0  0.0  6.0  6.0  2.0
//    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
//
//    Output matrix:
//
//     1.0  3.0  0.0  0.0  2.0  0.0  0.0
//     0.0  0.0  0.0  1.0  2.0  0.0  0.0
//     0.0  0.0  0.0  0.0  0.0  1.0  0.3
//     0.0  0.0  0.0  0.0  0.0  0.0  0.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of
//    the matrix A.
//
//    Input/output, float A[M*N].  On input, the matrix to be
//    analyzed.  On output, the RREF form of the matrix.
//
{
  int i;
  int j;
  int lead;
  int r;
  float temp;

  lead = 0;

  for ( r = 0; r < m; r++ )
  {
    if ( n - 1 < lead )
    {
      break;
    }

    i = r;

    while ( a[i+lead*m] == 0.0 )
    {
      i = i + 1;

      if ( m - 1 < i )
      {
        i = r;
        lead = lead + 1;
        if ( n - 1 < lead )
        {
          lead = -1;
          break;
         }
      }
    }

    if ( lead < 0 )
    {
      break;
    }

    for ( j = 0; j < n; j++ )
    {
      temp     = a[i+j*m];
      a[i+j*m] = a[r+j*m];
      a[r+j*m] = temp;
    }

    temp = a[r+lead*m];

    for ( j = 0; j < n; j++ )
    {
      a[r+j*m] = a[r+j*m] / temp;
    }

    for ( i = 0; i < m; i++ )
    {
      if ( i != r )
      {
        temp = a[i+lead*m];
        for ( j = 0; j < n; j++ )
        {
          a[i+j*m] = a[i+j*m] - temp * a[r+j*m];
        }
      }
    }
    lead = lead + 1;

  }
  return;
}
//****************************************************************************80

int r4mat_solve ( int n, int rhs_num, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*N]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
//    must be at least 0.
//
//    Input/output, float A[N*(N+RHS_NUM)], contains in rows and columns 1
//    to N the coefficient matrix, and in columns N+1 through
//    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
//    area has been destroyed, while the right hand sides have
//    been overwritten with the corresponding solutions.
//
//    Output, int R4MAT_SOLVE, singularity flag.
//    0, the matrix was not singular, the solutions were computed;
//    J, factorization failed on step J, and the solutions could not
//    be computed.
//
{
  float apivot;
  float factor;
  int i;
  int ipivot;
  int j;
  int k;
  float temp;

  for ( j = 0; j < n; j++ )
  {
//
//  Choose a pivot row.
//
    ipivot = j;
    apivot = a[j+j*n];

    for ( i = j; i < n; i++ )
    {
      if ( r4_abs ( apivot ) < r4_abs ( a[i+j*n] ) )
      {
        apivot = a[i+j*n];
        ipivot = i;
      }
    }

    if ( apivot == 0.0 )
    {
      return j;
    }
//
//  Interchange.
//
    for ( i = 0; i < n + rhs_num; i++ )
    {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }
//
//  A(J,J) becomes 1.
//
    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
    {
      a[j+k*n] = a[j+k*n] / apivot;
    }
//
//  A(I,J) becomes 0.
//
    for ( i = 0; i < n; i++ )
    {
      if ( i != j )
      {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
        {
          a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
        }
      }
    }
  }

  return 0;
}
//****************************************************************************80

float *r4mat_solve_2d ( float a[], float b[], float *det )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    If the determinant DET is returned as zero, then the matrix A is
//    singular, and does not have an inverse.  In that case, X is
//    returned as the NULL vector.
//
//    If DET is nonzero, then its value is roughly an estimate
//    of how nonsingular the matrix A is.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[2*2], the matrix.
//
//    Input, float B[2], the right hand side.
//
//    Output, float *DET, the determinant of the system.
//
//    Output, float R4MAT_SOLVE_2D[2], the solution of the system,
//    if DET is nonzero.  Otherwise, the NULL vector.
//
{
  float *x;
//
//  Compute the determinant.
//
  *det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];
//
//  If the determinant is zero, bail out.
//
  if ( *det == 0.0 )
  {
    return NULL;
  }
//
//  Compute the solution.
//
  x = new float[2];

  x[0] = (  a[1+1*2] * b[0] - a[0+1*2] * b[1] ) / ( *det );
  x[1] = ( -a[1+0*2] * b[0] + a[0+0*2] * b[1] ) / ( *det );

  return x;
}
//****************************************************************************80

float *r4mat_solve_3d ( float a[], float b[], float *det )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_SOLVE_3D solves a 3 by 3 linear system using Cramer's rule.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    If the determinant DET is returned as zero, then the matrix A is
//    singular, and does not have an inverse.  In that case, X is
//    returned as the NULL vector.
//
//    If DET is nonzero, then its value is roughly an estimate
//    of how nonsingular the matrix A is.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A[3*3], the matrix.
//
//    Input, float B[3], the right hand side.
//
//    Output, float *DET, the determinant of the system.
//
//    Output, float R4MAT_SOLVE_3D[3], the solution of the system,
//    if DET is nonzero.  Otherwise, the NULL vector.
//
{
  float *x;
//
//  Compute the determinant.
//
  *det =  a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
        + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
        + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );
//
//  If the determinant is zero, bail out.
//
  if ( *det == 0.0 )
  {
    return NULL;
  }
//
//  Compute the solution.
//
  x = new float[3];

  x[0] = (   ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) * b[0]
           - ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] ) * b[1]
           + ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] ) * b[2] ) / ( *det );

  x[1] = ( - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] ) * b[0]
           + ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] ) * b[1]
           - ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] ) * b[2] ) / ( *det );

  x[2] = (   ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] ) * b[0]
           - ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] ) * b[1]
           + ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] ) * b[2] ) / ( *det );

  return x;
}
//****************************************************************************80

float *r4mat_solve2 ( int n, float a[], float b[], int *ierror )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_SOLVE2 computes the solution of an N by N linear system.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    The linear system may be represented as
//
//      A*X = B
//
//    If the linear system is singular, but consistent, then the routine will
//    still produce a solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of equations.
//
//    Input/output, float A[N*N].
//    On input, A is the coefficient matrix to be inverted.
//    On output, A has been overwritten.
//
//    Input/output, float B[N].
//    On input, B is the right hand side of the system.
//    On output, B has been overwritten.
//
//    Output, float R4MAT_SOLVE2[N], the solution of the linear system.
//
//    Output, int *IERROR.
//    0, no error detected.
//    1, consistent singularity.
//    2, inconsistent singularity.
//
{
  float amax;
  int i;
  int imax;
  int j;
  int k;
  int *piv;
  float *x;

  *ierror = 0;

  piv = i4vec_zero_new ( n );
  x = r4vec_zero_new ( n );
//
//  Process the matrix.
//
  for ( k = 1; k <= n; k++ )
  {
//
//  In column K:
//    Seek the row IMAX with the properties that:
//      IMAX has not already been used as a pivot;
//      A(IMAX,K) is larger in magnitude than any other candidate.
//
    amax = 0.0;
    imax = 0;
    for ( i = 1; i <= n; i++ )
    {
      if ( piv[i-1] == 0 )
      {
        if ( amax < r4_abs ( a[i-1+(k-1)*n] ) )
        {
          imax = i;
          amax = r4_abs ( a[i-1+(k-1)*n] );
        }
      }
    }
//
//  If you found a pivot row IMAX, then,
//    eliminate the K-th entry in all rows that have not been used for pivoting.
//
    if ( imax != 0 )
    {
      piv[imax-1] = k;
      for ( j = k+1; j <= n; j++ )
      {
        a[imax-1+(j-1)*n] = a[imax-1+(j-1)*n] / a[imax-1+(k-1)*n];
      }
      b[imax-1] = b[imax-1] / a[imax-1+(k-1)*n];
      a[imax-1+(k-1)*n] = 1.0;

      for ( i = 1; i <= n; i++ )
      {
        if ( piv[i-1] == 0 )
        {
          for ( j = k+1; j <= n; j++ )
          {
            a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - a[i-1+(k-1)*n] * a[imax-1+(j-1)*n];
          }
          b[i-1] = b[i-1] - a[i-1+(k-1)*n] * b[imax-1];
          a[i-1+(k-1)*n] = 0.0;
        }
      }
    }
  }
//
//  Now, every row with nonzero PIV begins with a 1, and
//  all other rows are all zero.  Begin solution.
//
  for ( j = n; 1 <= j; j-- )
  {
    imax = 0;
    for ( k = 1; k <= n; k++ )
    {
      if ( piv[k-1] == j )
      {
        imax = k;
      }
    }

    if ( imax == 0 )
    {
      x[j-1] = 0.0;

      if ( b[j-1] == 0.0 )
      {
        *ierror = 1;
        cout << "\n";
        cout << "R4MAT_SOLVE2 - Warning:\n";
        cout << "  Consistent singularity, equation = " << j << "\n";
      }
      else
      {
        *ierror = 2;
        cout << "\n";
        cout << "R4MAT_SOLVE2 - Warning:\n";
        cout << "  Inconsistent singularity, equation = " << j << "\n";
      }
    }
    else
    {
      x[j-1] = b[imax-1];

      for ( i = 1; i <= n; i++ )
      {
        if ( i != imax )
        {
          b[i-1] = b[i-1] - a[i-1+(j-1)*n] * x[j-1];
        }
      }
    }
  }

  delete [] piv;

  return x;
}
//****************************************************************************80

float *r4mat_symm_eigen ( int n, float x[], float q[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_SYMM_EIGEN returns a symmetric matrix with given eigensystem.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    The user must supply the desired eigenvalue vector, and the desired
//    eigenvector matrix.  The eigenvector matrix must be orthogonal.  A
//    suitable random orthogonal matrix can be generated by
//    R4MAT_ORTH_UNIFORM_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Input, float X[N], the desired eigenvalues for the matrix.
//
//    Input, float Q[N*N], the eigenvector matrix of A.
//
//    Output, float R4MAT_SYMM_EIGEN[N*N], a symmetric N by N matrix with
//    eigenvalues X and eigenvectors the columns of Q.
//
{
  float *a;
  int i;
  int j;
  int k;
//
//  Set A = Q * Lambda * Q'.
//
  a = new float[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        a[i+j*n] = a[i+j*n] + q[i+k*n] * x[k] * q[j+k*n];
      }
    }
  }

  return a;
}
//****************************************************************************80

void r4mat_symm_jacobi ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_SYMM_JACOBI applies Jacobi eigenvalue iteration to a symmetric matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    This code was modified so that it treats as zero the off-diagonal
//    elements that are sufficiently close to, but not exactly, zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Input/output, float A[N*N], a symmetric N by N matrix.
//    On output, the matrix has been overwritten by an approximately
//    diagonal matrix, with the eigenvalues on the diagonal.
//
{
  float c;
  float eps = 0.00001;
  int i;
  int it;
  int it_max = 100;
  int j;
  int k;
  float norm_fro;
  float s;
  float sum2;
  float t;
  float t1;
  float t2;
  float u;

  norm_fro = r4mat_norm_fro ( n, n, a );

  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < i; j++ )
      {
        if ( eps * norm_fro < r4_abs ( a[i+j*n] ) + r4_abs ( a[j+i*n] ) )
        {
          u = ( a[j+j*n] - a[i+i*n] ) / ( a[i+j*n] + a[j+i*n] );

          t = r4_sign ( u ) / ( r4_abs ( u ) + sqrt ( u * u + 1.0 ) );
          c = 1.0 / sqrt ( t * t + 1.0 );
          s = t * c;
//
//  A -> A * Q.
//
          for ( k = 0; k < n; k++ )
          {
            t1 = a[i+k*n];
            t2 = a[j+k*n];
            a[i+k*n] = t1 * c - t2 * s;
            a[j+k*n] = t1 * s + t2 * c;
          }
//
//  A -> QT * A
//
          for ( k = 0; k < n; k++ )
          {
            t1 = a[k+i*n];
            t2 = a[k+j*n];
            a[k+i*n] = c * t1 - s * t2;
            a[k+j*n] = s * t1 + c * t2;
          }
        }
      }
    }
//
//  Test the size of the off-diagonal elements.
//
    sum2 = 0.0;
    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < i; j++ )
      {
        sum2 = sum2 + r4_abs ( a[i+j*n] );
      }
    }

    if ( sum2 <= eps * ( norm_fro + 1.0 ) )
    {
      break;
    }

    if ( it_max <= it )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

int r4mat_to_r4plu ( int n, float a[], int pivot[], float lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_TO_R4PLU factors a general matrix.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    This routine is a simplified version of the LINPACK routine DGEFA.
//    Fortran conventions are used to index doubly-dimensioned arrays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, float A[N*N], the matrix to be factored.
//
//    Output, int PIVOT[N], a vector of pivot indices.
//
//    Output, float LU[N*N], an upper triangular matrix U and the multipliers
//    L which were used to obtain it.  The factorization can be written
//    A = L * U, where L is a product of permutation and unit lower
//    triangular matrices and U is upper triangular.
//
//    Output, int R4MAT_TO_R4PLU, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the R4MAT_TO_R4PLU-th step.
//
{
  int i;
  int info;
  int j;
  int k;
  int l;
  float temp;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      lu[i+j*n] = a[i+j*n];
    }
  }
  info = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L, the index of the pivot row.
//
    l = k;
    for ( i = k+1; i <= n; i++ )
    {
      if ( r4_abs ( lu[l-1+(k-1)*n] ) < r4_abs ( lu[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;
//
//  If the pivot index is zero, the algorithm has failed.
//
    if ( lu[l-1+(k-1)*n] == 0.0 )
    {
      info = k;
      return info;
    }
//
//  Interchange rows L and K if necessary.
//
    if ( l != k )
    {
      temp            = lu[l-1+(k-1)*n];
      lu[l-1+(k-1)*n] = lu[k-1+(k-1)*n];
      lu[k-1+(k-1)*n] = temp;
    }
//
//  Normalize the values that lie below the pivot entry A(K,K).
//
    for ( i = k+1; i <= n; i++ )
    {
      lu[i-1+(k-1)*n] = -lu[i-1+(k-1)*n] / lu[k-1+(k-1)*n];
    }
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        temp            = lu[l-1+(j-1)*n];
        lu[l-1+(j-1)*n] = lu[k-1+(j-1)*n];
        lu[k-1+(j-1)*n] = temp;
      }

      for ( i = k+1; i <= n; i++ )
      {
        lu[i-1+(j-1)*n] = lu[i-1+(j-1)*n] + lu[i-1+(k-1)*n] * lu[k-1+(j-1)*n];
      }
    }
  }

  pivot[n-1] = n;

  if ( lu[n-1+(n-1)*n] == 0.0 )
  {
    info = n;
  }

  return info;
}
//****************************************************************************80

float r4mat_trace ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_TRACE computes the trace of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    The trace of a square matrix is the sum of the diagonal elements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, float A[N*N], the matrix whose trace is desired.
//
//    Output, float R4MAT_TRACE, the trace of the matrix.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i+i*n];
  }

  return value;
}
//****************************************************************************80

float *r4mat_transpose ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_TRANSPOSE returns the transpose of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix A.
//
//    Input, float A[M*N], the matrix whose transpose is desired.
//
//    Output, float R4MAT_TRANSPOSE[N*M], the transposed matrix.
//
{
  float *b;
  int i;
  int j;

  b = new float[n*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[j+i*n] = a[i+j*m];
    }
  }
  return b;
}
//****************************************************************************80

void r4mat_transpose_in_place ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_TRANSPOSE_IN_PLACE transposes a square R4MAT in place.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix A.
//
//    Input/output, float A[N*N], the matrix to be transposed.
//
{
  int i;
  int j;
  float t;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      t        = a[i+j*n];
      a[i+j*n] = a[j+i*n];
      a[j+i*n] = t;
    }
  }
  return;
}
//****************************************************************************80

void r4mat_transpose_print ( int m, int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_TRANSPOSE_PRINT prints an R4MAT, transposed.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
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
//    Input, float A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r4mat_transpose_print_some ( int m, int n, float a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_TRANSPOSE_PRINT_SOME prints some of an R4MAT, transposed.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
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
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], an M by N matrix to be printed.
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

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

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

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

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

float *r4mat_u_inverse ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_U_INVERSE inverts an upper triangular R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    An upper triangular matrix is a matrix whose only nonzero entries
//    occur on or above the diagonal.
//
//    The inverse of an upper triangular matrix is an upper triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, number of rows and columns in the matrix.
//
//    Input, float A[N*N], the upper triangular matrix.
//
//    Output, float R4MAT_U_INVERSE[N*N], the inverse matrix.
//
{
  float *b;
  int i;
  int j;
  int k;

  b = new float[n*n];

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
      else
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

float *r4mat_u1_inverse ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_U1_INVERSE inverts a unit upper triangular R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    A unit upper triangular matrix is a matrix with only 1's on the main
//    diagonal, and only 0's below the main diagonal.
//
//    The inverse of a unit upper triangular matrix is also
//    a unit upper triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    C++ translation by John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, number of rows and columns in the matrix.
//
//    Input, float A[N*N], the unit upper triangular matrix.
//
//    Output, float R4MAT_U1_INVERSE[N*N), the inverse matrix.
//
{
  float *b;
  int i;
  int j;
  int k;

  b = new float[n*n];

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
        b[i+j*n] = 1.0;
      }
      else
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

void r4mat_uniform_01 ( int m, int n, int &seed, float r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_UNIFORM_01 returns a unit pseudorandom R4MAT.
//
//  Discussion:
//
//    An R4MAT is an array of R4's.
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
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has
//    been updated.
//
//    Output, float R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = ( float ) ( seed ) * 4.656612875E-10;
    }
  }
  return;
}
//****************************************************************************80

float *r4mat_uniform_new ( int m, int n, float a, float b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_UNIFORM_NEW returns a scaled pseudorandom R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4's,  stored as a vector
//    in column-major order.
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
//    Input, float A, B, the limits of the pseudorandom values.
//
//    Input/output, int &SEED, the "seed" value.  On output, SEED has
//    been updated.
//
//    Output, float R4MAT_UNIFORM_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  float *r;

  r = new float[m*n];

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
      r[i+j*m] = a + ( b - a ) * ( float ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

float *r4mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_UNIFORM_01_NEW returns a unit pseudorandom R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4's,  stored as a vector
//    in column-major order.
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
//    Input/output, int &SEED, the "seed" value.  On output, SEED has
//    been updated.
//
//    Output, float R4MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  float *r;

  r = new float[m*n];

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
      r[i+j*m] = ( float ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

float *r4mat_vand2 ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_VAND2 returns the N by N row Vandermonde matrix A.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    The row Vandermonde matrix returned by this routine reads "across"
//    rather than down.  In particular, each row begins with a 1, followed by
//    some value X, followed by successive powers of X.
//
//    The formula for the matrix entries is:
//
//      A(I,J) = X(I)**(J-1)
//
//  Properties:
//
//    A is nonsingular if, and only if, the X values are distinct.
//
//    The determinant of A is
//
//      det(A) = product ( 2 <= I <= N ) (
//        product ( 1 <= J <= I-1 ) ( ( X(I) - X(J) ) ) ).
//
//    The matrix A is generally ill-conditioned.
//
//  Example:
//
//    N = 5, X = (2, 3, 4, 5, 6)
//
//    1 2  4   8   16
//    1 3  9  27   81
//    1 4 16  64  256
//    1 5 25 125  625
//    1 6 36 216 1296
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix desired.
//
//    Input, float X[N], the values that define A.
//
//    Output, float R4MAT_VAND2[N*N], the N by N row Vandermonde matrix.
//
{
  float *a;
  int i;
  int j;

  a = new float[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( j == 0 && x[i] == 0.0 )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = pow ( x[i], j );
      }
    }
  }

  return a;
}
//****************************************************************************80

void r4mat_zero ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_ZERO zeroes an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, float A[M*N], a matrix of zeroes.
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
  return;
}
//****************************************************************************80

float *r4mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_ZERO_NEW returns a new zeroed R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, float R4MAT_ZERO_NEW[M*N], the new zeroed matrix.
//
{
  float *a;
  int i;
  int j;

  a = new float[m*n];

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

float r4plu_det ( int n, int pivot[], float lu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4PLU_DET computes the determinant of a real PLU matrix.
//
//  Discussion:
//
//    The matrix should have been factored by R4MAT_TO_R4PLU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int PIVOT[N], the pivot vector computed by R4MAT_TO_R4PLU.
//
//    Input, float LU[N*N], the LU factors computed by R4MAT_TO_R4PLU.
//
//    Output, float R4PLU_DET, the determinant of the matrix.
//
{
  float det;
  int i;

  det = 1.0;

  for ( i = 0; i < n; i++ )
  {
    det = det * lu[i+i*n];
    if ( pivot[i] != i+1 )
    {
      det = -det;
    }
  }

  return det;
}
//****************************************************************************80

void r4plu_inverse ( int n, int pivot[], float lu[], float a_inverse[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4PLU_INVERSE computes the inverse of a real PLU matrix.
//
//  Discussion:
//
//    The matrix should have been factored by R4MAT_TO_R4PLU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, int PIVOT[N], the pivot vector from R4MAT_TO_R4PLU.
//
//    Input, float LU[N*N], the LU factors computed by R4MAT_TO_R4PLU.
//
//    Output, float A_INVERSE[N*N], the inverse of the original matrix
//    A that was factored by R4MAT_TO_R4PLU.
//
{
  int i;
  int j;
  int k;
  float temp;
  float *work;
//
  work = new float[n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a_inverse[i+j*n] = lu[i+j*n];
    }
  }
//
//  Compute Inverse(U).
//
  for ( k = 1; k <= n; k++ )
  {
    a_inverse[k-1+(k-1)*n]     = 1.0 / a_inverse[k-1+(k-1)*n];
    for ( i = 1; i <= k-1; i++ )
    {
      a_inverse[i-1+(k-1)*n] = -a_inverse[i-1+(k-1)*n] * a_inverse[k-1+(k-1)*n];
    }

    for ( j = k+1; j <= n; j++ )
    {
      temp                     = a_inverse[k-1+(j-1)*n];
      a_inverse[k-1+(j-1)*n]   = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        a_inverse[i-1+(j-1)*n] = a_inverse[i-1+(j-1)*n]
          + temp * a_inverse[i-1+(k-1)*n];
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
      work[i-1] = a_inverse[i-1+(k-1)*n];
      a_inverse[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        a_inverse[i-1+(k-1)*n] = a_inverse[i-1+(k-1)*n]
          + a_inverse[i-1+(j-1)*n] * work[j-1];
      }
    }

    if ( pivot[k-1] != k )
    {
      for ( i = 1; i <= n; i++ )
      {
        temp                            = a_inverse[i-1+(k-1)*n];
        a_inverse[i-1+(k-1)*n]          = a_inverse[i-1+(pivot[k-1]-1)*n];
        a_inverse[i-1+(pivot[k-1]-1)*n] = temp;
      }
    }
  }

  delete [] work;

  return;
}
//****************************************************************************80

void r4plu_mul ( int n, int pivot[], float lu[], float x[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4PLU_MUL computes A * x using the PLU factors of A.
//
//  Discussion:
//
//    It is assumed that R4MAT_TO_R4PLU has computed the PLU factors of
//    the matrix A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2003
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
//    Input, int PIVOT[N], the pivot vector computed by R4MAT_TO_R4PLU.
//
//    Input, float LU[N*N], the matrix factors computed by R4MAT_TO_R4PLU.
//
//    Input, float X[N], the vector to be multiplied.
//
//    Output, float B[N], the result of the multiplication.
//
{
  int i;
  int j;
  int k;
  float temp;
//
  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }
//
//  Y = U * X.
//
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 0; i < j-1; i++ )
    {
      b[i] = b[i] + lu[i+(j-1)*n] * b[j-1];
    }
    b[j-1] = lu[j-1+(j-1)*n] * b[j-1];
  }
//
//  B = PL * Y = PL * U * X = A * x.
//
  for ( j = n-1; 1 <= j; j-- )
  {
    for ( i = j; i < n; i++ )
    {
      b[i] = b[i] - lu[i+(j-1)*n] * b[j-1];
    }

    k = pivot[j-1];

    if ( k != j )
    {
      temp = b[k-1];
      b[k-1] = b[j-1];
      b[j-1] = temp;
    }
  }

  return;
}
//****************************************************************************80

void r4plu_sol ( int n, int pivot[], float lu[], float b[], float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4PLU_SOL solves a linear system A*x=b from the PLU factors.
//
//  Discussion:
//
//    The PLU factors should have been computed by R4MAT_TO_R4PLU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int PIVOT[N], the pivot vector from R4MAT_TO_R4PLU.
//
//    Input, float LU[N*N], the LU factors from R4MAT_TO_R4PLU.
//
//    Input, float B[N], the right hand side vector.
//
//    Output, float X[N], the solution vector.
//
{
  int i;
  int j;
  int k;
  float temp;
//
//  Solve PL * Y = B.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( k = 1; k <= n-1; k++ )
  {
    j = pivot[k-1];

    if ( j != k )
    {
      temp   = x[j-1];
      x[j-1] = x[k-1];
      x[k-1] = temp;
    }

    for ( i = k+1; i <= n; i++ )
    {
      x[i-1] = x[i-1] + lu[i-1+(k-1)*n] * x[k-1];
    }
  }
//
//  Solve U * X = Y.
//
  for ( k = n; 1 <= k; k-- )
  {
    x[k-1] = x[k-1] / lu[k-1+(k-1)*n];
    for ( i = 1; i <= k-1; i++ )
    {
      x[i-1] = x[i-1] - lu[i-1+(k-1)*n] * x[k-1];
    }
  }

  return;
}
//****************************************************************************80

void r4plu_to_r4mat ( int n, int pivot[], float lu[], float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4PLU_TO_R4MAT recovers the matrix A that was factored by R4MAT_TO_R4PLU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2003
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
//    Input, int PIVOT[N], the pivot vector computed by R4MAT_TO_R4PLU.
//
//    Input, float LU[N*N], the matrix factors computed by R4MAT_TO_R4PLU.
//
//    Output, float A[N*N], the matrix whose factors are represented by
//    LU and PIVOT.
//
{
  int i;
  int j;
  int k;
  float temp;

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

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n; i++ )
    {
      for ( k = 1; k <= i-1; k++ )
      {
        a[k-1+(j-1)*n] = a[k-1+(j-1)*n] + lu[k-1+(i-1)*n] * a[i-1+(j-1)*n];
      }
      a[i-1+(j-1)*n] = lu[i-1+(i-1)*n] * a[i-1+(j-1)*n];
    }
//
//  B = PL * Y = PL * U * X = A * x.
//
    for ( i = n-1; 1 <= i; i-- )
    {
      for ( k = i+1; k <= n; k++ )
      {
        a[k-1+(j-1)*n] = a[k-1+(j-1)*n] - lu[k-1+(i-1)*n] * a[i-1+(j-1)*n];
      }

      k = pivot[i-1];

      if ( k != i )
      {
        temp           = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = a[i-1+(j-1)*n];
        a[i-1+(j-1)*n] = temp;
      }
    }
  }

  return;
}
//****************************************************************************80

int r4poly_degree ( int na, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_DEGREE returns the degree of a polynomial.
//
//  Discussion:
//
//    The degree of a polynomial is the index of the highest power
//    of X with a nonzero coefficient.
//
//    The degree of a constant polynomial is 0.  The degree of the
//    zero polynomial is debatable, but this routine returns the
//    degree as 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, float A[NA+1], the coefficients of the polynomials.
//
//    Output, int R4POLY_DEGREE, the degree of A.
//
{
  int degree;

  degree = na;

  while ( 0 < degree )
  {
    if ( a[degree] != 0.0 )
    {
      return degree;
    }
    degree = degree - 1;
  }

  return degree;
}
//****************************************************************************80

float *r4poly_deriv ( int n, float c[], int p )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_DERIV returns the derivative of a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the degree of the polynomial.
//
//    Input, float C[N+1], the polynomial coefficients.
//    C[I] is the coefficient of X^I.
//
//    Input, int P, the order of the derivative.
//    0 means no derivative is taken.
//    1 means first derivative,
//    2 means second derivative and so on.
//    Values of P less than 0 are meaningless.  Values of P greater
//    than N are meaningful, but the code will behave as though the
//    value of P was N+1.
//
//    Output, float R4POLY_DERIV CP[N-P+1], the polynomial coefficients of
//    the derivative.
//
{
  float *cp;
  float *cp_temp;
  int d;
  int i;

  if ( n < p )
  {
    return NULL;
  }
  cp_temp = r4vec_copy_new ( n+1, c );

  for ( d = 1; d <= p; d++ )
  {
    for ( i = 0; i <= n-d; i++ )
    {
      cp_temp[i] = ( float ) ( i + 1 ) * cp_temp[i+1];
    }
    cp_temp[n-d+1] = 0.0;
  }

  cp = r4vec_copy_new ( n - p + 1, cp_temp );

  delete [] cp_temp;

  return cp;
}
//****************************************************************************80

float r4poly_lagrange_0 ( int npol, float xpol[], float xval )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_LAGRANGE_0 evaluates the Lagrange factor at a point.
//
//  Discussion:
//
//    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
//
//  Discussion:
//
//    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
//    polynomial L(IPOL)(X), has the property:
//
//      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
//
//    and may be expressed as:
//
//      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
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
//    Input, int NPOL, the number of abscissas.
//    NPOL must be at least 1.
//
//    Input, float XPOL[NPOL], the abscissas, which should be distinct.
//
//    Input, float XVAL, the point at which the Lagrange factor is to be
//    evaluated.
//
//    Output, float R4POLY_LAGRANGE_0, the value of the Lagrange factor at XVAL.
//
{
  int i;
  float wval;

  wval = 1.0;
  for ( i = 0; i < npol; i++ )
  {
    wval = wval * ( xval - xpol[i] );
  }

  return wval;
}
//****************************************************************************80

float r4poly_lagrange_1 ( int npol, float xpol[], float xval )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_LAGRANGE_1 evaluates the first derivative of the Lagrange factor.
//
//  Discussion:
//
//    W(XPOL(1:NPOL))(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
//
//    W'(XPOL(1:NPOL))(X)
//      = Sum ( 1 <= J <= NPOL ) Product ( I /= J ) ( X - XPOL(I) )
//
//    We also have the recursion:
//
//      W'(XPOL(1:NPOL))(X) = d/dX ( ( X - XPOL(NPOL) ) * W(XPOL(1:NPOL-1))(X) )
//                    = W(XPOL(1:NPOL-1))(X)
//                    + ( X - XPOL(NPOL) ) * W'(XPOL(1:NPOL-1))(X)
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
//    Input, int NPOL, the number of abscissas.
//
//    Input, float XPOL[NPOL], the abscissas, which should be distinct.
//
//    Input, float XVAL, the point at which the Lagrange factor is to be
//    evaluated.
//
//    Output, float R4POLY_LAGRANGE_1, the derivative of W with respect to XVAL.
//
{
  float dwdx;
  int i;
  float w;

  dwdx = 0.0;
  w = 1.0;

  for ( i = 0; i < npol; i++ )
  {
    dwdx = w + ( xval - xpol[i] ) * dwdx;
    w = w * ( xval - xpol[i] );
  }

  return dwdx;
}
//****************************************************************************80

float r4poly_lagrange_2 ( int npol, float xpol[], float xval )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_LAGRANGE_2 evaluates the second derivative of the Lagrange factor.
//
//  Discussion:
//
//    W(X)  = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
//
//    W'(X) = Sum ( 1 <= J <= NPOL )
//            Product ( I /= J ) ( X - XPOL(I) )
//
//    W"(X) = Sum ( 1 <= K <= NPOL )
//            Sum ( J =/ K )
//            Product ( I /= K, J ) ( X - XPOL(I) )
//
//    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
//    polynomial L(IPOL)(X), has the property:
//
//      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
//
//    and may be expressed as:
//
//      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
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
//    Input, int NPOL, the number of abscissas.
//    NPOL must be at least 1.
//
//    Input, float XPOL[NPOL], the abscissas, which should be distinct.
//
//    Input, float XVAL, the point at which the Lagrange factor is to be
//    evaluated.
//
//    Output, float R4POLY_LAGRANGE_2, the second derivative of W with respect to XVAL.
//
{
  float dw2dx2;
  int i;
  int j;
  int k;
  float term;

  dw2dx2 = 0.0;

  for ( k = 0; k < npol; k++ )
  {
    for ( j = 0; j < npol; j++ )
    {
      if ( j != k )
      {
        term = 1.0;
        for ( i = 0; i < npol; i++ )
        {
          if ( i != j && i != k )
          {
            term = term * ( xval - xpol[i] );
          }
        }
        dw2dx2 = dw2dx2 + term;
      }
    }
  }

  return dw2dx2;
}
//****************************************************************************80

float *r4poly_lagrange_coef ( int npol, int ipol, float xpol[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_LAGRANGE_COEF returns the coefficients of a Lagrange polynomial.
//
//  Discussion:
//
//    Given NPOL distinct abscissas, XPOL(*), the IPOL-th Lagrange
//    polynomial P(IPOL)(X) is defined as the polynomial of degree
//    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
//    abscissas.
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
//    Input, int NPOL, the number of abscissas.
//    NPOL must be at least 1.
//
//    Input, int IPOL, the index of the polynomial to evaluate.
//    IPOL must be between 1 and NPOL.
//
//    Input, float XPOL[NPOL], the abscissas of the Lagrange polynomials.
//    The entries in XPOL must be distinct.
//
//    Output, float R4POLY_LAGRANGE_COEF[NPOL], the polynomial coefficients
//    of the IPOL-th Lagrange polynomial.
//
{
  int i;
  int index;
  int j;
  float *pcof;
//
//  Make sure IPOL is legal.
//
  if ( ipol < 1 || npol < ipol )
  {
    cerr << "\n";
    cerr << "R4POLY_LAGRANGE_COEF - Fatal error!\n";
    cerr << "  1 <= IPOL <= NPOL is required.\n";
    cerr << "  but IPOL = " << ipol << "\n";
    cerr << "  and NPOL = " << npol << "\n";
    exit ( 1 );
  }
//
//  Check that the abscissas are distinct.
//
  if ( !r4vec_distinct ( npol, xpol ) )
  {
    cerr << "\n";
    cerr << "R4POLY_LAGRANGE_COEF - Fatal error!\n";
    cerr << "  Two entries of XPOL are equal:\n";
    exit ( 1 );
  }

  pcof = new float[npol];

  pcof[0] = 1.0;
  for ( i = 1; i < npol; i++ )
  {
    pcof[i] = 0.0;
  }

  index = 0;

  for ( i = 1; i <= npol; i++ )
  {
    if ( i != ipol )
    {
      index = index + 1;

      for ( j = index; 0 <= j; j-- )
      {
        pcof[j] = - xpol[i-1] * pcof[j] / ( xpol[ipol-1] - xpol[i-1] );

        if ( 0 < j )
        {
          pcof[j] = pcof[j] + pcof[j-1] / ( xpol[ipol-1] - xpol[i-1] );
        }
      }
    }
  }

  return pcof;
}
//****************************************************************************80

void r4poly_lagrange_factor ( int npol, float xpol[], float xval,
  float *wval, float *dwdx )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_LAGRANGE_FACTOR evaluates the polynomial Lagrange factor at a point.
//
//  Discussion:
//
//    Suppose F(X) is at least N times continuously differentiable in the
//    interval [A,B].  Pick NPOL distinct points XPOL(I) in [A,B] and compute
//    the interpolating polynomial P(X) of order NPOL ( and degree NPOL-1)
//    which passes through all the points ( XPOL(I), F(XPOL(I)) ).
//    Then in the interval [A,B], the maximum error
//
//      abs ( F(X) - P(X) )
//
//    is bounded by:
//
//      C * FNMAX * W(X)
//
//    where
//
//      C is a constant,
//      FNMAX is the maximum value of the NPOL-th derivative of F in [A,B],
//      W(X) is the Lagrange factor.
//
//    Thus, the value of W(X) is useful as part of an estimated bound
//    for the interpolation error.
//
//    The formula is:
//
//      W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
//
//    Note that the Chebyshev abscissas have the property that they minimize
//    the value of W(X) over the interval [A,B].  Hence, if the abscissas may
//    be chosen arbitrarily, the Chebyshev abscissas have this advantage over
//    other choices.
//
//    For a set of points XPOL[I], 0 <= I <= NPOL-1, the IPOL-th Lagrange basis
//    polynomial L(IPOL)(X), has the property:
//
//      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
//
//    and may be expressed as:
//
//      L(IPOL)(X) = W(X) / ( ( X - XPOL[IPOL] ) * W'(XPOL[IPOL]) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NPOL, the number of abscissas.
//    NPOL must be at least 1.
//
//    Input, float XPOL[NPOL], the abscissas, which should be distinct.
//
//    Input, float XVAL, the point at which the Lagrange factor is to be evaluated.
//
//    Output, float *WVAL, the value of the Lagrange factor at XVAL.
//
//    Output, float *DWDX, the derivative of W with respect to XVAL.
//
{
  int i;
  int j;
  float term;

  *wval = 1.0;
  for ( i = 0; i < npol; i++ )
  {
    *wval = *wval * ( xval - xpol[i] );
  }

  *dwdx = 0.0;

  for ( i = 0; i < npol; i++ )
  {
    term = 1.0;

    for ( j = 0; j < npol; j++ )
    {
      if ( i != j )
      {
        term = term * ( xval - xpol[j] );
      }
    }
    *dwdx = *dwdx + term;
  }

  return;
}
//****************************************************************************80

int r4poly_lagrange_val ( int npol, int ipol, float xpol[], float xval,
  float *pval, float *dpdx )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_LAGRANGE_VAL evaluates the IPOL-th Lagrange polynomial.
//
//  Discussion:
//
//    Given NPOL distinct abscissas, XPOL[*], the IPOL-th Lagrange
//    polynomial P(IPOL)(X) is defined as the polynomial of degree
//    NPOL - 1 which is 1 at XPOL[IPOL] and 0 at the NPOL - 1 other
//    abscissas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NPOL, the number of abscissas.
//    NPOL must be at least 1.
//
//    Input, int IPOL, the index of the polynomial to evaluate.
//    IPOL must be between 0 and NPOL-1.
//
//    Input, float XPOL[NPOL], the abscissas of the Lagrange polynomials.
//    The entries in XPOL must be distinct.
//
//    Input, float XVAL, the point at which the IPOL-th Lagrange polynomial
//    is to be evaluated.
//
//    Output, float *PVAL, the value of the IPOL-th Lagrange polynomial at XVAL.
//
//    Output, float *DPDX, the derivative of the IPOL-th Lagrange polynomial at XVAL.
//
//    Output, int R4POLY_LAGRANGE_VAL, 0 if no error.
{
  int i;
  int j;
  float p2;
//
//  Make sure IPOL is legal.
//
  if ( ipol < 0 || npol-1 < ipol )
  {
    cerr << "\n";
    cerr << "R4POLY_LAGRANGE_VAL - Fatal error!\n";
    cerr << "  0 <= IPOL <= NPOL-1 is required.\n";
    exit ( 1 );
  }
//
//  Check that the abscissas are distinct.
//
  for ( i = 1; i < npol; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      if ( xpol[i] == xpol[j] )
      {
        cerr << "\n";
        cerr << "R4POLY_LAGRANGE_VAL - Fatal error!\n";
        cerr << "  Two entries of XPOL are equal:\n";
        cerr << "  XPOL(" << i << ") = " << xpol[i] << ".\n";
        cerr << "  XPOL(" << j << ") = " << xpol[j] << ".\n";
        exit ( 1 );
      }
    }
  }
//
//  Evaluate the polynomial.
//
  *pval = 1.0;

  for ( i = 0; i < npol; i++ )
  {
    if ( i != ipol )
    {
      *pval = *pval * ( xval - xpol[i] ) / ( xpol[ipol] - xpol[i] );
    }
  }
//
//  Evaluate the derivative, which can be found by summing up the result
//  of differentiating one factor at a time, successively.
//
  *dpdx = 0.0;

  for ( i = 0; i < npol; i++ )
  {
    if ( i != ipol )
    {
      p2 = 1.0;

      for ( j = 0; j < npol; j++ )
      {
        if ( j == i )
        {
          p2 = p2                      / ( xpol[ipol] - xpol[j] );
        }
        else if ( j != ipol )
        {
          p2 = p2 * ( xval - xpol[j] ) / ( xpol[ipol] - xpol[j] );
        }
      }
      *dpdx = *dpdx + p2;
    }
  }

  return 0;
}
//****************************************************************************80

int r4poly_order ( int na, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_ORDER returns the order of a polynomial.
//
//  Discussion:
//
//    The order of a polynomial is one more than the degree.
//
//    The order of a constant polynomial is 1.  The order of the
//    zero polynomial is debatable, but this routine returns the
//    order as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, float A[NA+1], the coefficients of the polynomials.
//
//    Output, int R4POLY_ORDER, the order of A.
//
{
  int order;

  order = na + 1;

  while ( 1 < order )
  {
    if ( a[order-1] != 0.0 )
    {
      return order;
    }
    order = order - 1;
  }

  return order;
}
//****************************************************************************80

void r4poly_print ( int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_PRINT prints out a polynomial.
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
//    Input, int N, the dimension of A.
//
//    Input, float A[N+1], the polynomial coefficients.
//    A(0) is the constant term and
//    A(N) is the coefficient of X**N.
//
//    Input, string TITLE, a title.
//
{
  int i;
  float mag;
  int n2;
  char plus_minus;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  n2 = r4poly_degree ( n, a );

  if ( n2 <= 0 )
  {
    cout << "  p(x) = 0\n";
    return;
  }

  if ( a[n2] < 0.0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = r4_abs ( a[n2] );

  if ( 2 <= n2 )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << " * x ^ " << n2 << "\n";
  }
  else if ( n2 == 1 )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << " * x\n";
  }
  else if ( n2 == 0 )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << "\n";
  }

  for ( i = n2-1; 0 <= i; i-- )
  {
    if ( a[i] < 0.0 )
    {
      plus_minus = '-';
    }
    else
    {
      plus_minus = '+';
    }

    mag = r4_abs ( a[i] );

    if ( mag != 0.0 )
    {
      if ( 2 <= i )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << " * x ^ " << i << "\n";
      }
      else if ( i == 1 )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << " * x\n";
      }
      else if ( i == 0 )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void r4poly_shift ( float scale, float shift, int n, float poly_cof[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
//
//  Discussion:
//
//    Assuming P(X) is a polynomial in the argument X, of the form:
//
//      P(X) =
//          C(N) * X**N
//        + ...
//        + C(1) * X
//        + C(0),
//
//    and that Z is related to X by the formula:
//
//      Z = SCALE * X + SHIFT
//
//    then this routine computes coefficients C for the polynomial Q(Z):
//
//      Q(Z) =
//          C(N) * Z**N
//        + ...
//        + C(1) * Z
//        + C(0)
//
//    so that:
//
//      Q(Z(X)) = P(X)
//
//  Example:
//
//    P(X) = 2 * X**2 - X + 6
//
//    Z = 2.0D+00 * X + 3.0D+00
//
//    Q(Z) = 0.5 *         Z**2 -  3.5 * Z + 12
//
//    Q(Z(X)) = 0.5 * ( 4.0D+00 * X**2 + 12.0D+00 * X +  9 )
//            - 3.5 * (               2.0D+00 * X +  3 )
//                                            + 12
//
//            = 2.0D+00         * X**2 -  1.0D+00 * X +  6
//
//            = P(X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2005
//
//  Reference:
//
//    Press, Flannery, Teukolsky, Vetterling,
//    Numerical Recipes: The Art of Scientific Computing,
//    Cambridge University Press.
//
//  Parameters:
//
//    Input, float SHIFT, SCALE, the shift and scale applied to X,
//    so that Z = SCALE * X + SHIFT.
//
//    Input, int N, the number of coefficients.
//
//    Input/output, float POLY_COF[N+1].
//    On input, the coefficient array in terms of the X variable.
//    On output, the coefficient array in terms of the Z variable.
//
{
  int i;
  int j;

  for ( i = 1; i <= n; i++ )
  {
    for ( j = i; j <= n; j++ )
    {
      poly_cof[j] = poly_cof[j] / scale;
    }
  }

  for ( i = 0; i <= n-1; i++ )
  {
    for ( j = n-1; i <= j; j-- )
    {
      poly_cof[j] = poly_cof[j] - shift * poly_cof[j+1];
    }
  }

  return;
}
//****************************************************************************80

float r4poly_val_horner ( int n, float c[], float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_VAL_HORNER evaluates a polynomial using Horner's method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of C.
//
//    Input, float C[N+1], the polynomial coefficients.
//    C(I) is the coefficient of X**I.
//
//    Input, float X, the point at which the polynomial is
//    to be evaluated.
//
//    Output, float R4POLY_VAL_HORNER, the value of the polynomial at X.
//
{
  int i;
  float value;

  value = c[n];
  for ( i = n-1; 0 <= i; i-- )
  {
    value = value * x + c[i];
  }

  return value;
}
//****************************************************************************80

float r4poly_value ( int n, float a[], float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY_VALUE evaluates an R4POLY.
//
//  Discussion:
//
//    For sanity's sake, the value of N indicates the NUMBER of
//    coefficients, or more precisely, the ORDER of the polynomial,
//    rather than the DEGREE of the polynomial.  The two quantities
//    differ by 1, but cause a great deal of confusion.
//
//    Given N and A, the form of the polynomial is:
//
//      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, float A[N], the coefficients of the polynomial.
//    A[0] is the constant term.
//
//    Input, float X, the point at which the polynomial is to be evaluated.
//
//    Output, float R4POLY_VALUE, the value of the polynomial at X.
//
{
  int i;
  float value;

  value = 0.0;

  for ( i = n-1; 0 <= i; i-- )
  {
    value = value * x + a[i];
  }

  return value;
}
//****************************************************************************80

int r4poly2_ex ( float x1, float y1, float x2, float y2, float x3,
  float y3, float *x, float *y )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY2_EX finds the extremal point of a parabola determined by three points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
//    on the parabola.  X1, X2 and X3 must be distinct.
//
//    Output, float *X, *Y, the X coordinate of the extremal point of the
//    parabola, and the value of the parabola at that point.
//
//    Output, int R4POLY2_EX, error flag.
//    0, no error.
//    1, two of the X values are equal.
//    2, the data lies on a straight line; there is no finite extremal
//    point.
//    3, the data lies on a horizontal line; every point is "extremal".
//
{
  float bot;

  *x = 0.0;
  *y = 0.0;

  if ( x1 == x2 || x2 == x3 || x3 == x1 )
  {
    return 1;
  }

  if ( y1 == y2 && y2 == y3 && y3 == y1 )
  {
    *x = x1;
    *y = y1;
    return 3;
  }

  bot = ( x2 - x3 ) * y1 + ( x3 - x1 ) * y2 + ( x1 - x2 ) * y3;

  if ( bot == 0.0 )
  {
    return 2;
  }

  *x = 0.5 * (
      x1 * x1 * ( y3 - y2 )
    + x2 * x2 * ( y1 - y3 )
    + x3 * x3 * ( y2 - y1 ) ) /
    ( ( x2 - x3 ) * y1 + ( x3 - x1 ) * y2 + ( x1 - x2 ) * y3 );

  *y = - (
      ( *x - x2 ) * ( *x - x3 ) * ( x2 - x3 ) * y1
    + ( *x - x1 ) * ( *x - x3 ) * ( x3 - x1 ) * y2
    + ( *x - x1 ) * ( *x - x2 ) * ( x1 - x2 ) * y3 ) /
    ( ( x1 - x2 ) * ( x2 - x3 ) * ( x3 - x1 ) );

  return 0;
}
//****************************************************************************80

int r4poly2_ex2 ( float x1, float y1, float x2, float y2, float x3,
  float y3, float *x, float *y, float *a, float *b, float *c )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY2_EX2 finds the extremal point of a parabola determined by three points.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, float X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
//    on the parabola.  X1, X2 and X3 must be distinct.
//
//    Output, float *X, *Y, the X coordinate of the extremal point of the
//    parabola, and the value of the parabola at that point.
//
//    Output, float *A, *B, *C, the coefficients that define the parabola:
//    P(X) = A * X**2 + B * X + C.
//
//    Output, int R4POLY2_EX2, error flag.
//    0, no error.
//    1, two of the X values are equal.
//    2, the data lies on a straight line; there is no finite extremal
//    point.
//    3, the data lies on a horizontal line; any point is an "extremal point".
//
{
  float v[3*3];
  float *w;

  *a = 0.0;
  *b = 0.0;
  *c = 0.0;
  *x = 0.0;
  *y = 0.0;

  if ( x1 == x2 || x2 == x3 || x3 == x1 )
  {
    return 1;
  }

  if ( y1 == y2 && y2 == y3 && y3 == y1 )
  {
    *x = x1;
    *y = y1;
    return 3;
  }
//
//  Set up the Vandermonde matrix.
//
  v[0+0*3] = 1.0;
  v[0+1*3] = x1;
  v[0+2*3] = x1 * x1;

  v[1+0*3] = 1.0;
  v[1+1*3] = x2;
  v[1+2*3] = x2 * x2;

  v[2+0*3] = 1.0;
  v[2+1*3] = x3;
  v[2+2*3] = x3 * x3;
//
//  Get the inverse.
//
  w = r4mat_inverse_3d ( v );
//
//  Compute the parabolic coefficients.
//
  *c = w[0+0*3] * y1 + w[0+1*3] * y2 + w[0+2*3] * y3;
  *b = w[1+0*3] * y1 + w[1+1*3] * y2 + w[1+2*3] * y3;
  *a = w[2+0*3] * y1 + w[2+1*3] * y2 + w[2+2*3] * y3;
//
//  Determine the extremal point.
//
  if ( *a == 0.0 )
  {
    return 2;
  }

  *x = - *b / ( 2.0 * *a );
  *y = *a * *x * *x + *b * *x + *c;

  return 0;
}
//****************************************************************************80

void r4poly2_root ( float a, float b, float c, complex <float> *r1,
  complex <float> *r2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY2_ROOT returns the two roots of a quadratic polynomial.
//
//  Discussion:
//
//    The polynomial has the form:
//
//      A * X**2 + B * X + C = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2005
//
//  Parameters:
//
//    Input, float A, B, C, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <float> *R1, *R2, the roots of the polynomial, which
//    might be real and distinct, real and equal, or complex conjugates.
//
{
  float disc;
  complex <float> q;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R4POLY2_ROOT - Fatal error!\n";
    cerr << "  The coefficient A is zero.\n";
    exit ( 1 );
  }

  disc = b * b - 4.0 * a * c;
//
//  For some frap-dapping reason, C++ is unable to combine these quantities!
//
//q = - 0.5 * ( r4_sign ( b ) * r4_csqrt ( disc ) );

  q = - complex <float> ( 0.5, 0.0 ) * ( r4_sign ( b ) * r4_csqrt ( disc ) );

  *r1 = q / a;
  *r2 = c / q;

  return;
}
//****************************************************************************80

void r4poly2_rroot ( float a, float b, float c, float *r1, float *r2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
//
//  Example:
//
//    A    B    C       roots              R1   R2
//   --   --   --     ------------------   --   --
//    1   -4    3     1          3          1    3
//    1    0    4     2*i      - 2*i        0    0
//    2   -6    5     3 +   i    3 -   i    3    3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, B, C, the coefficients of the quadratic
//    polynomial A * X**2 + B * X + C = 0 whose roots are desired.
//    A must not be zero.
//
//    Output, float *R1, *R2, the real parts of the roots
//    of the polynomial.
//
{
  float disc;
  float q;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R4POLY2_RROOT - Fatal error!\n";
    cerr << "  The coefficient A is zero.\n";
    exit ( 1 );
  }

  disc = b * b - 4.0 * a * c;
  disc = r4_max ( disc, 0.0 );

  q = ( b + r4_sign ( b ) * sqrt ( disc ) );
  *r1 = -0.5 * q / a;
  *r2 = -2.0 * c / q;

  return;
}
//****************************************************************************80

void r4poly2_val ( float x1, float y1, float x2, float y2,
  float x3, float y3, float x, float *y, float *yp, float *ypp )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY2_VAL evaluates a parabola defined by three data values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X1, Y1, X2, Y2, X3, Y3, three pairs of data values.
//    If the X values are distinct, then all the Y values represent
//    actual values of the parabola.
//
//    Three special cases are allowed:
//
//      X1 = X2 =/= X3: Y2 is the derivative at X1;
//      X1 =/= X2 = X3: Y3 is the derivative at X3;
//      X1 = X2 = X3:   Y2 is the derivative at X1, and
//                      Y3 is the second derivative at X1.
//
//    Input, float X, an abscissa at which the parabola is to be
//    evaluated.
//
//    Output, float *Y, *YP, *YPP, the values of the parabola and
//    its first and second derivatives at X.
//
{
  int distinct;
  float dif1;
  float dif2;
  float temp;
//
//  If any X's are equal, put them and the Y data first.
//
  if ( x1 == x2 && x2 == x3 )
  {
    distinct = 1;
  }
  else if ( x1 == x2 )
  {
    distinct = 2;
  }
  else if ( x1 == x3 )
  {
    cerr << "\n";
    cerr << "R4POLY2_VAL - Fatal error!\n";
    cerr << "  X1 = X3 =/= X2.\n";
    return;
  }
  else if ( x2 == x3 )
  {
    distinct = 2;
    temp = x1;
    x1 = x3;
    x3 = temp;
    temp = y1;
    y1 = y2;
    y2 = y3;
    y3 = y1;
  }
  else
  {
    distinct = 3;
  }
//
//  Set up the coefficients.
//
  if ( distinct == 1 )
  {
    dif1 = y2;
    dif2 = 0.5 * y3;
  }
  else if ( distinct == 2 )
  {
    dif1 = y2;
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 )
             - y2 ) / ( x3 - x2 );
  }
  else if ( distinct == 3 )
  {
    dif1 = ( y2 - y1 ) / ( x2 - x1 );
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 )
           - ( y2 - y1 ) / ( x2 - x1 ) ) / ( x3 - x2 );
  }
//
//  Evaluate.
//
  *y = y1 + ( x - x1 ) * dif1 + ( x - x1 ) * ( x - x2 ) * dif2;
  *yp = dif1 + ( 2.0 * x - x1 - x2 ) * dif2;
  *ypp = 2.0 * dif2;

  return;
}
//****************************************************************************80

void r4poly2_val2 ( int ndata, float tdata[],
  float ydata[], int left, float tval, float *yval )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY2_VAL2 evaluates a parabolic function through 3 points in a table.
//
//  Discussion:
//
//    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
//    It constructs the parabolic interpolant through the data in
//    3 consecutive entries of a table and evaluates this interpolant
//    at a given abscissa value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points.
//    NDATA must be at least 3.
//
//    Input, float TDATA[NDATA], the abscissas of the data points.  The
//    values in TDATA must be in strictly ascending order.
//
//    Input, float YDATA[NDATA], the data points corresponding to
//    the abscissas.
//
//    Input, int LEFT, the location of the first of the three
//    consecutive data points through which the parabolic interpolant
//    must pass.  0 <= LEFT <= NDATA - 3.
//
//    Input, float TVAL, the value of T at which the parabolic interpolant
//    is to be evaluated.  Normally, TDATA[0] <= TVAL <= T[NDATA-1], and
//    the data will be interpolated.  For TVAL outside this range,
//    extrapolation will be used.
//
//    Output, float *YVAL, the value of the parabolic interpolant
//    at TVAL.
//
{
  float dif1;
  float dif2;
  float t1;
  float t2;
  float t3;
  float y1;
  float y2;
  float y3;
//
//  Check.
//
  if ( left < 0 || ndata-3 < left )
  {
    cerr << "\n";
    cerr << "RPOLY2_VAL2 - Fatal error!\n";
    cerr << "  LEFT < 0 or NDATA-3 < LEFT.\n";
    exit ( 1 );
  }
//
//  Copy out the three abscissas.
//
  t1 = tdata[left];
  t2 = tdata[left+1];
  t3 = tdata[left+2];

  if ( t2 <= t1 || t3 <= t2 )
  {
    cerr << "\n";
    cerr << "RPOLY2_VAL2 - Fatal error!\n";
    cerr << "  T2 <= T1 or T3 <= T2.\n";
    cerr << "  T1 = " << t1 << "\n";
    cerr << "  T2 = " << t2 << "\n";
    cerr << "  T3 = " << t3 << "\n";
    exit ( 1 );
  }
//
//  Construct and evaluate a parabolic interpolant for the data.
//
  y1 = ydata[left];
  y2 = ydata[left+1];
  y3 = ydata[left+2];

  dif1 = ( y2 - y1 ) / ( t2 - t1 );
  dif2 =
    ( ( y3 - y1 ) / ( t3 - t1 )
    - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 );

  *yval = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 );

  return;
}
//****************************************************************************80

void r4poly3_root ( float a, float b, float c, float d,
  complex <float> *r1, complex <float> *r2, complex <float> *r3 )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY3_ROOT returns the three roots of a cubic polynomial.
//
//  Discussion:
//
//    The polynomial has the form
//
//      A * X**3 + B * X**2 + C * X + D = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2005
//
//  Parameters:
//
//    Input, float A, B, C, D, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <float> *R1, *R2, *R3, the roots of the polynomial, which
//    will include at least one real root.
//
{
  complex <float> i;
  float q;
  float r;
  const float r4_pi = 3.141592653589793;
  float s1;
  float s2;
  float temp;
  float theta;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R4POLY3_ROOT - Fatal error!\n";
    cerr << "  A must not be zero.\n";
    exit ( 1 );
  }

  i = complex <float> ( 0.0, 1.0 );

  q = ( pow ( b / a, 2 ) - 3.0 * ( c / a ) ) / 9.0;

  r = ( 2.0 * pow ( b / a, 3 ) - 9.0 * ( b / a ) * ( c / a )
      + 27.0 * ( d / a ) ) / 54.0;

  if ( r * r < q * q * q )
  {
    theta = acos ( r / sqrt ( pow ( q, 3 ) ) );
    *r1 = -2.0 * sqrt ( q ) * cos (   theta                 / 3.0 );
    *r2 = -2.0 * sqrt ( q ) * cos ( ( theta + 2.0 * r4_pi ) / 3.0 );
    *r3 = -2.0 * sqrt ( q ) * cos ( ( theta + 4.0 * r4_pi ) / 3.0 );
  }
  else if ( q * q * q <= r * r )
  {
    temp = -r + sqrt ( r * r - q * q * q );
    s1 = r4_sign ( temp ) * pow ( ( double ) r4_abs ( temp ), 1.0 / 3.0 );

    temp = -r - sqrt ( r * r - q * q * q );
    s2 = r4_sign ( temp ) * pow ( ( double ) r4_abs ( temp ), 1.0 / 3.0 );

    *r1 = s1 + s2;
    *r2 = complex <float> ( - 0.5 * ( s1 + s2 ),
                              0.5 * sqrt ( 3.0 ) * ( s1 - s2 ) );

    *r3 = complex <float> ( - 0.5 * ( s1 + s2 ), 
                            - 0.5 * sqrt ( 3.0 ) * ( s1 - s2 ) );
  }

  *r1 = *r1 - complex <float> ( b / ( 3.0 * a ) );
  *r2 = *r2 - complex <float> ( b / ( 3.0 * a ) );
  *r3 = *r3 - complex <float> ( b / ( 3.0 * a ) );

  return;
}
//****************************************************************************80

void r4poly4_root ( float a, float b, float c, float d, float e,
  complex <float> *r1, complex <float> *r2, complex <float> *r3,
  complex <float> *r4 )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY4_ROOT returns the four roots of a quartic polynomial.
//
//  Discussion:
//
//    The polynomial has the form:
//
//      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2005
//
//  Parameters:
//
//    Input, float A, B, C, D, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <float> *R1, *R2, *R3, *R4, the roots of the polynomial.
//
{
  float a3;
  float a4;
  float b3;
  float b4;
  float c3;
  float c4;
  float d3;
  float d4;
  complex <float> p;
  complex <float> q;
  complex <float> r;
  complex <float> zero;

  zero = complex <float> ( 0.0, 0.0 );

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R4POLY4_ROOT - Fatal error!\n";
    cerr << "  A must not be zero.\n";
    exit ( 1 );
  }

  a4 = b / a;
  b4 = c / a;
  c4 = d / a;
  d4 = e / a;
//
//  Set the coefficients of the resolvent cubic equation.
//
  a3 = 1.0;
  b3 = -b4;
  c3 = a4 * c4 - 4.0 * d4;
  d3 = -a4 * a4 * d4 + 4.0 * b4 * d4 - c4 * c4;
//
//  Find the roots of the resolvent cubic.
//
  r4poly3_root ( a3, b3, c3, d3, r1, r2, r3 );
//
//  Choose one root of the cubic, here R1.
//
//  Set R = sqrt ( 0.25 * A4**2 - B4 + R1 )
//
  r = c4_sqrt ( complex <float> ( 0.25, 0.0 ) * a4 * a4 - b4 + (*r1) );

  if ( real ( r ) != 0.0 || imag ( r ) != 0.0 )
  {
    p = c4_sqrt ( complex <float> ( 0.75, 0.0 ) * a4 * a4 - r * r 
      - complex <float> ( 2.0, 0.0 ) * b4
      + complex <float> ( 0.25, 0.0 ) 
      * ( complex <float> ( 4.0, 0.0 ) * a4 * b4 
      - complex <float> ( 8.0, 0.0 ) * c4 - a4 * a4 * a4 ) / r );

    q = c4_sqrt ( complex <float> ( 0.75, 0.0 ) * a4 * a4 - r * r 
      - complex <float> ( 2.0, 0.0 ) * b4
      - complex <float> ( 0.25, 0.0 ) 
      * ( complex <float> ( 4.0, 0.0 ) * a4 * b4 
      - complex <float> ( 8.0, 0.0 ) * c4 - a4 * a4 * a4 ) / r );
  }
  else
  {
    p = c4_sqrt ( complex <float> ( 0.75, 0.0 ) * a4 * a4 
      - complex <float> ( 2.0, 0.0 ) * b4
      + complex <float> ( 2.0, 0.0 ) * c4_sqrt ( (*r1) * (*r1) 
      - complex <float> ( 4.0, 0.0 ) * d4 ) );

    q = c4_sqrt ( complex <float> ( 0.75, 0.0 ) * a4 * a4 
      - complex <float> ( 2.0, 0.0 ) * b4
      - complex <float> ( 2.0, 0.0 ) * c4_sqrt ( (*r1) * (*r1) 
      - complex <float> ( 4.0, 0.0 ) * d4 ) );
  }
//
//  Set the roots.
//
  *r1 = - complex <float> ( 0.25, 0.0 ) * a4 + complex <float> ( 0.5, 0.0 ) * r 
    + complex <float> ( 0.5, 0.0 ) * p;
  *r2 = - complex <float> ( 0.25, 0.0 ) * a4 + complex <float> ( 0.5, 0.0 ) * r 
    - complex <float> ( 0.5, 0.0 ) * p;
  *r3 = - complex <float> ( 0.25, 0.0 ) * a4 - complex <float> ( 0.5, 0.0 ) * r 
    + complex <float> ( 0.5, 0.0 ) * q;
  *r4 = - complex <float> ( 0.25, 0.0 ) * a4 - complex <float> ( 0.5, 0.0 ) * r 
    - complex <float> ( 0.5, 0.0 ) * q;

  return;
}
//****************************************************************************80

void r4pp_delete ( float **a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4PP_DELETE frees the memory set aside by R4PP_NEW.
//
//  Discussion:
//
//    An R4PP is a pointer to pointers to R4's, and is a sort of
//    variably-dimensioned matrix.
//
//    This function releases the memory associated with an array that was 
//    created by a command like:
//
//      float **a;
//      a = r4pp_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float **A, the pointer to the pointers.
//
//    Input, int M, N, the number of rows and columns in the array.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}
//****************************************************************************80

float **r4pp_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4PP_NEW allocates a new R4PP.
//
//  Discussion:
//
//    An R4PP is a pointer to pointers to R4's, and is a sort of
//    variably-dimensioned matrix.
//
//    A declaration of the form
//      float **a;
//    is necesary.  Then an assignment of the form:
//      a = r4pp_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, float **R4PP_NEW, a pointer to the pointers to the M by N array.
//
{
  float **a;
  int i;

  a = new float *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R4PP_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new float[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R4PP_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}

//****************************************************************************80

int r4r4_compare ( float x1, float y1, float x2, float y2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4R4_COMPARE compares two R4R4's.
//
//  Discussion:
//
//    An R4R4 is simply a pair of R4 values, stored separately.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X1, Y1, the first vector.
//
//    Input, float X2, Y2, the second vector.
//
//    Output, int R4R4_COMPARE:
//    -1, (X1,Y1) < (X2,Y2);
//     0, (X1,Y1) = (X2,Y2);
//    +1, (X1,Y1) > (X2,Y2).
//
{
  int value;

  if ( x1 < x2 )
  {
    value = -1;
  }
  else if ( x2 < x1 )
  {
    value = +1;
  }
  else if ( y1 < y2 )
  {
    value = -1;
  }
  else if ( y2 < y1 )
  {
    value = +1;
  }
  else
  {
    value = 0;
  }

  return value;
}
//****************************************************************************80

void r4r4_print ( float a1, float a2, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4R4_PRINT prints an R4R4.
//
//  Discussion:
//
//    An R4R4 is a pair of R4 values, regarded as a single item.
//
//    A format is used which suggests a coordinate pair:
//
//  Example:
//
//    Center : ( 1.23, 7.45 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A1, A2, the coordinates of the vector.
//
//    Input, string TITLE, a title.
//
{
  cout << "  " << title << " : ";
  cout << "  ( " << setw(12) << a1
       << ", "   << setw(12) << a2 << " )\n";

  return;
}
//****************************************************************************80

int r4r4r4_compare ( float x1, float y1, float z1, float x2, float y2,
  float z2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4R4R4_COMPARE compares two R4R4R4's.
//
//  Discussion:
//
//    An R4R4R4 is simply 3 R4 values, stored as scalars.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X1, Y1, Z1, the first vector.
//
//    Input, float X2, Y2, Z2, the second vector.
//
//    Output, int R4R4R4_COMPARE:
//    -1, (X1,Y1,Z1) < (X2,Y2,Z2);
//     0, (X1,Y1,Z1) = (X2,Y2,Z2);
//    +1, (X1,Y1,Z1) > (X2,Y2,Z2).
//
{
  int value;

  if ( x1 < x2 )
  {
    value = -1;
  }
  else if ( x2 < x1 )
  {
    value = +1;
  }
  else if ( y1 < y2 )
  {
    value = -1;
  }
  else if ( y2 < y1 )
  {
    value = +1;
  }
  else if ( z1 < z2 )
  {
    value = -1;
  }
  else if ( z2 < z1 )
  {
    value = +1;
  }
  else
  {
    value = 0;
  }

  return value;
}
//****************************************************************************80

void r4r4r4vec_index_insert_unique ( int maxn, int *n, float x[], float y[],
  float z[], int indx[], float xval, float yval, float zval, int *ival,
  int *ierror )

//****************************************************************************80
//
//  Purpose:
//
//    R4R4R4VEC_INDEX_INSERT_UNIQUE inserts a unique R4R4R4 value in an indexed sorted list.
//
//  Discussion:
//
//    If the input value does not occur in the current list, it is added,
//    and N, X, Y, Z and INDX are updated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MAXN, the maximum size of the list.
//
//    Input/output, int *N, the size of the list.
//
//    Input/output, float X[N], Y[N], Z[N], the R4R4R4 vector.
//
//    Input/output, int INDX[N], the sort index of the list.
//
//    Input, float XVAL, YVAL, ZVAL, the value to be inserted
//    if it is not already in the list.
//
//    Output, int *IVAL, the index in X, Y, Z corresponding to the
//    value XVAL, YVAL, ZVAL.
//
//    Output, int *IERROR, 0 for no error, 1 if an error occurred.
//
{
  int equal;
  int i;
  int less;
  int more;

  *ierror = 0;

  if ( *n <= 0 )
  {
    if ( maxn <= 0 )
    {
      *ierror = 1;
      cerr << "\n";
      cerr << "R4R4R4VEC_INDEX_INSERT_UNIQUE - Fatal error!\n";
      cerr << "  Not enough space to store new data.\n";
      return;
    }
    *n = 1;
    x[0] = xval;
    y[0] = yval;
    z[0] = zval;
    indx[0] = 1;
    *ival = 1;
    return;
  }
//
//  Does ( XVAL, YVAL, ZVAL ) already occur in ( X, Y, Z)?
//
  r4r4r4vec_index_search ( *n, x, y, z, indx, xval, yval, zval,
    &less, &equal, &more );

  if ( equal == 0 )
  {
    if ( maxn <= *n )
    {
      *ierror = 1;
      cerr << "\n";
      cerr << "R4R4R4VEC_INDEX_INSERT_UNIQUE - Fatal error!\n";
      cerr << "  Not enough space to store new data.\n";
      return;
    }

    x[*n] = xval;
    y[*n] = yval;
    z[*n] = zval;
    *ival = *n + 1;
    for ( i = *n-1; more-1 <= i; i-- )
    {
      indx[i+1] = indx[i];
    }
    indx[more-1] = *n + 1;
    *n = *n + 1;
  }
  else
  {
    *ival = indx[equal-1];
  }

  return;
}
//****************************************************************************80

void r4r4r4vec_index_search ( int n, float x[], float y[], float z[],
  int indx[], float xval, float yval, float zval, int *less, int *equal,
  int *more )

//****************************************************************************80
//
//  Purpose:
//
//    R4R4R4VEC_INDEX_SEARCH searches for an R4R4R4 value in an indexed sorted list.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the list.
//
//    Input, float X[N], Y[N], Z[N], the list.
//
//    Input, int INDX[N], the sort index of the list.
//
//    Input, float XVAL, YVAL, ZVAL, the value to be sought.
//
//    Output, int *LESS, *EQUAL, *MORE, the indexes in INDX of the
//    entries of X that are just less than, equal to, and just greater
//    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
//    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
//    is the greatest entry of X, then MORE is N+1.
//
{
  int compare;
  int hi;
  int lo;
  int mid;
  float xhi;
  float xlo;
  float xmid;
  float yhi;
  float ylo;
  float ymid;
  float zhi;
  float zlo;
  float zmid;

  if ( n <= 0 )
  {
    *less = 0;
    *equal = 0;
    *more = 0;
    return;
  }

  lo = 1;
  hi = n;

  xlo = x[indx[lo-1]-1];
  ylo = y[indx[lo-1]-1];
  zlo = z[indx[lo-1]-1];

  xhi = x[indx[hi-1]-1];
  yhi = y[indx[hi-1]-1];
  zhi = z[indx[hi-1]-1];

  compare = r4r4r4_compare ( xval, yval, zval, xlo, ylo, zlo );

  if ( compare == -1 )
  {
    *less = 0;
    *equal = 0;
    *more = 1;
    return;
  }
  else if ( compare == 0 )
  {
    *less = 0;
    *equal = 1;
    *more = 2;
    return;
  }

  compare = r4r4r4_compare ( xval, yval, zval, xhi, yhi, zhi );

  if ( compare == 1 )
  {
    *less = n;
    *equal = 0;
    *more = n + 1;
    return;
  }
  else if ( compare == 0 )
  {
    *less = n - 1;
    *equal = n;
    *more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      *less = lo;
      *equal = 0;
      *more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];
    ymid = y[indx[mid-1]-1];
    zmid = z[indx[mid-1]-1];

    compare = r4r4r4_compare ( xval, yval, zval, xmid, ymid, zmid );

    if ( compare == 0 )
    {
      *equal = mid;
      *less = mid - 1;
      *more = mid + 1;
      return;
    }
    else if ( compare == -1 )
    {
      hi = mid;
    }
    else if ( compare == +1 )
    {
      lo = mid;
    }
  }

  return;
}
//****************************************************************************80

void r4r4vec_index_insert_unique ( int maxn, int *n, float x[], float y[],
  int indx[], float xval, float yval, int *ival, int *ierror )

//****************************************************************************80
//
//  Purpose:
//
//    R4R4VEC_INDEX_INSERT_UNIQUE inserts a unique R4R4 value in an indexed sorted list.
//
//  Discussion:
//
//    If the input value does not occur in the current list, it is added,
//    and N, X, Y and INDX are updated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MAXN, the maximum size of the list.
//
//    Input/output, int *N, the size of the list.
//
//    Input/output, float X[N], Y[N], the list of R4R4 vectors.
//
//    Input/output, int INDX[N], the sort index of the list.
//
//    Input, float XVAL, YVAL, the value to be inserted if it is
//    not already in the list.
//
//    Output, int *IVAL, the index in X, Y corresponding to the
//    value XVAL, YVAL.
//
//    Output, int *IERROR, 0 for no error, 1 if an error occurred.
//
{
  int equal;
  int i;
  int less;
  int more;

  *ierror = 0;

  if ( *n <= 0 )
  {
    if ( maxn <= 0 )
    {
      cerr << "\n";
      cerr << "R4R4VEC_INDEX_INSERT_UNIQUE - Fatal error!\n";
      cerr << "  Not enough space to store new data.\n";
      exit ( 1 );
    }

    *n = 1;
    x[0] = xval;
    y[0] = yval;
    indx[0] = 1;
    *ival = 1;
    return;
  }
//
//  Does ( XVAL, YVAL ) already occur in ( X, Y )?
//
  r4r4vec_index_search ( *n, x, y, indx, xval, yval, &less, &equal, &more );

  if ( equal == 0 )
  {
    if ( maxn <= *n )
    {
      cerr << "\n";
      cerr << "R4R4VEC_INDEX_INSERT_UNIQUE - Fatal error!\n";
      cerr << "  Not enough space to store new data.\n";
      exit ( 1 );
    }

    x[*n] = xval;
    y[*n] = yval;
    *ival = *n + 1;
    for ( i = *n-1; more-1 <= i; i-- )
    {
      indx[i+1] = indx[i];
    }
    indx[more-1] = *n + 1;
    *n = *n + 1;
  }
  else
  {
    *ival = indx[equal-1];
  }

  return;
}
//****************************************************************************80

void r4r4vec_index_search ( int n, float x[], float y[], int indx[],
  float xval, float yval, int *less, int *equal, int *more )

//****************************************************************************80
//
//  Purpose:
//
//    R4R4VEC_INDEX_SEARCH searches for an R4R4 value in an indexed sorted list.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the current list.
//
//    Input, float X[N], Y[N], the list.
//
//    Input, int INDX[N], the sort index of the list.
//
//    Input, float XVAL, YVAL, the value to be sought.
//
//    Output, int *LESS, *EQUAL, *MORE, the indexes in INDX of the
//    entries of X that are just less than, equal to, and just greater
//    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
//    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
//    is the greatest entry of X, then MORE is N+1.
//
{
  int compare;
  int hi;
  int lo;
  int mid;
  float xhi;
  float xlo;
  float xmid;
  float yhi;
  float ylo;
  float ymid;

  if ( n <= 0 )
  {
    *less = 0;
    *equal = 0;
    *more = 0;
    return;
  }

  lo = 1;
  hi = n;

  xlo = x[indx[lo-1]-1];
  ylo = y[indx[lo-1]-1];

  xhi = x[indx[hi-1]-1];
  yhi = y[indx[hi-1]-1];

  compare = r4r4_compare ( xval, yval, xlo, ylo );

  if ( compare == -1 )
  {
    *less = 0;
    *equal = 0;
    *more = 1;
    return;
  }
  else if ( compare == 0 )
  {
    *less = 0;
    *equal = 1;
    *more = 2;
    return;
  }

  compare = r4r4_compare ( xval, yval, xhi, yhi );

  if ( compare == 1 )
  {
    *less = n;
    *equal = 0;
    *more = n + 1;
    return;
  }
  else if ( compare == 0 )
  {
    *less = n - 1;
    *equal = n;
    *more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      *less = lo;
      *equal = 0;
      *more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];
    ymid = y[indx[mid-1]-1];

    compare = r4r4_compare ( xval, yval, xmid, ymid );

    if ( compare == 0 )
    {
      *equal = mid;
      *less = mid - 1;
      *more = mid + 1;
      return;
    }
    else if ( compare == -1 )
    {
      hi = mid;
    }
    else if ( compare == +1 )
    {
      lo = mid;
    }
  }

  return;
}
//****************************************************************************80

float *r4row_max ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4ROW_MAX returns the row maximums of an R4ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MAX =
//      3
//      7
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
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, float R4ROW_MAX[M], the maximums of the rows.
//
{
  int i;
  int j;
  float *amax;

  amax = new float[m];

  for ( i = 0; i < m; i++ )
  {
    amax[i] = a[i+0*m];

    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[i+j*m] )
      {
        amax[i] = a[i+j*m];
      }
    }
  }

  return amax;
}
//****************************************************************************80

float *r4row_mean ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4ROW_MEAN returns the row means of an R4ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MEAN =
//      2
//      5
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
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, float R4ROW_MEAN[M], the means, or averages, of the rows.
//
{
  int i;
  int j;
  float *mean;

  mean = new float[m];

  for ( i = 0; i < m; i++ )
  {
    mean[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      mean[i] = mean[i] + a[i+j*m];
    }
    mean[i] = mean[i] / ( float ) ( n );
  }

  return mean;
}
//****************************************************************************80

float *r4row_min ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4ROW_MIN returns the row minimums of an R4ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MIN =
//      1
//      2
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
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, float A[M*N], the array to be examined.
//
//    Output, float R4ROW_MIN[M], the minimums of the rows.
//
{
  int i;
  int j;
  float *amin;

  amin = new float[m];

  for ( i = 0; i < m; i++ )
  {
    amin[i] = a[i+0*m];
    for ( j = 1; j < n; j++ )
    {
      if ( a[i+j*m] < amin[i] )
      {
        amin[i] = a[i+j*m];
      }
    }
  }

  return amin;
}
//****************************************************************************80

float *r4row_sum ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4ROW_SUM returns the sums of the rows of an R4ROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the M by N array.
//
//    Output, float ROWSUM[M], the sum of the entries of
//    each row.
//
{
  int i;
  int j;
  float *rowsum;

  rowsum = new float[m];

  for ( i = 0; i < m; i++ )
  {
    rowsum[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      rowsum[i] = rowsum[i] + a[i+j*m];
    }
  }

  return rowsum;
}
//****************************************************************************80

void r4row_swap ( int m, int n, float a[], int irow1, int irow2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4ROW_SWAP swaps two rows of an R4ROW.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based.  However, a preprocessor
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
//    Input/output, float A[M*N], an array of data.
//
//    Input, int IROW1, IROW2, the two rows to swap.
//    These indices should be between 1 and M.
//
{
# define OFFSET 1

  int j;
  float t;
//
//  Check.
//
  if ( irow1 < 1 || m < irow1 )
  {
    cerr << "\n";
    cerr << "R4ROW_SWAP - Fatal error!\n";
    cerr << "  IROW1 is out of range.\n";
    exit ( 1 );
  }

  if ( irow2 < 1 || m < irow2 )
  {
    cerr << "\n";
    cerr << "R4ROW_SWAP - Fatal error!\n";
    cerr << "  IROW2 is out of range.\n";
    exit ( 1 );
  }

  if ( irow1 == irow2 )
  {
    return;
  }

  for ( j = 0; j < n; j++ )
  {
    t              = a[irow1-1+j*m];
    a[irow1-1+j*m] = a[irow2-1+j*m];
    a[irow2-1+j*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

float *r4row_to_r4vec ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4ROW_TO_R4VEC converts an R4ROW into an R4VEC.
//
//  Example:
//
//    M = 3, N = 4
//
//    A =
//      11 12 13 14
//      21 22 23 24
//      31 32 33 34
//
//    R4ROW_TO_R4VEC = ( 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float A[M*N], the M by N array.
//
//    Output, float R4ROW_TO_R4VEC[M*N], a vector containing the M rows of A.
//
{
  int i;
  int j;
  int k;
  float *x;

  x = new float[m*n];

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

float *r4row_variance ( int m, int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4ROW_VARIANCE returns the variances of an R4ROW.
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
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, float A[M*N], the array whose variances are desired.
//
//    Output, float R4ROW_VARIANCE[M], the variances of the rows.
//
{
  int i;
  int j;
  float mean;
  float *variance;

  variance = new float[m];

  for ( i = 0; i < m; i++ )
  {
    mean = 0.0;
    for ( j = 0; j < n; j++ )
    {
      mean = mean + a[i+j*m];
    }
    mean = mean / ( float ) ( n );

    variance[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      variance[i] = variance[i] + pow ( ( a[i+j*m] - mean ), 2 );
    }

    if ( 1 < n )
    {
      variance[i] = variance[i] / ( float ) ( n - 1 );
    }
    else
    {
      variance[i] = 0.0;
    }

  }

  return variance;
}
//****************************************************************************80

void r4slmat_print ( int m, int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4SLMAT_PRINT prints a strict lower triangular R4MAT.
//
//  Example:
//
//    M = 5, N = 5
//    A = (/ 21, 31, 41, 51, 32, 42, 52, 43, 53, 54 /)
//
//    21
//    31 32
//    41 42 43
//    51 52 53 54
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2005
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
//    Input, float A[*], the M by N matrix.  Only the strict
//    lower triangular elements are stored, in column major order.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int indx;
  int j;
  int jhi;
  int jlo;
  int jmax;
  int nn;

  cout << "\n";
  cout << title << "\n";

  jmax = i4_min ( n, m - 1 );

  nn = 5;

  for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
  {
    jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
    cout << "\n";
    cout << "  Col   ";
    for ( j = jlo; j <= jhi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    for ( i = jlo + 1; i <= m; i++ )
    {
      cout << setw(5) << i << ":";
      jhi = i4_min ( jlo + nn - 1, i4_min ( i - 1, jmax ) );
      for ( j = jlo; j <= jhi; j++ )
      {
        indx = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2;
        cout << " " << setw(12) << a[indx-1];
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

void r4vec_01_to_ab ( int n, float a[], float amax, float amin )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_01_TO_AB shifts and rescales data to lie within given bounds.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    On input, A contains the original data, which is presumed to lie
//    between 0 and 1.  However, it is not necessary that this be so.
//
//    On output, A has been shifted and rescaled so that all entries which
//    on input lay in [0,1] now lie between AMIN and AMAX.  Other entries will
//    be mapped in a corresponding way.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input/output, float A[N], the vector to be rescaled.
//
//    Input, float AMAX, AMIN, the maximum and minimum values
//    allowed for A.
//
{
  float amax2;
  float amax3;
  float amin2;
  float amin3;
  int i;

  if ( amax == amin )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = amin;
    }
    return;
  }

  amax2 = r4_max ( amax, amin );
  amin2 = r4_min ( amax, amin );

  amin3 = r4vec_min ( n, a );
  amax3 = r4vec_max ( n, a );

  if ( amax3 != amin3 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( amax3 - a[i]         ) * amin2
             + (         a[i] - amin3 ) * amax2 )
             / ( amax3          - amin3 );
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0.5 * ( amax2 + amin2 );
    }
  }

  return;
}
//****************************************************************************80

void r4vec_ab_to_01 ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_AB_TO_01 shifts and rescales data to lie within [0,1].
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    On input, A contains the original data.  On output, A has been shifted
//    and scaled so that all entries lie between 0 and 1.
//
//    The formula is:
//
//      A(I) := ( A(I) - AMIN ) / ( AMAX - AMIN )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input/output, float A[N], the data to be rescaled.
//
{
  float amax;
  float amin;
  int i;

  amax = r4vec_max ( n, a );
  amin = r4vec_min ( n, a );

  if ( amin == amax )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0.5;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( a[i] - amin ) / ( amax - amin );
    }
  }

  return;
}
//****************************************************************************80

float *r4vec_ab_to_cd ( int n, float a[], float bmin, float bmax )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_AB_TO_CD shifts and rescales data to lie within a given pair of bounds.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The mininum entry of A is mapped to BMIN, the maximum entry
//    to BMAX, and values in between are mapped linearly.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, float A[N], the data to be remapped.
//
//    Input, float BMIN, BMAX, the values to which min(A) and max(A)
//    are to be assigned.
//
//    Output, float R4VEC_AB_TO_CD[N], the remapped data.
//
{
  float amax;
  float amin;
  float *b;
  int i;

  b = new float[n];

  if ( bmax == bmin )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = bmin;
    }
    return b;
  }

  amax = r4vec_max ( n, a );
  amin = r4vec_min ( n, a );

  if ( amin == amax )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = 0.5 * ( bmax + bmin );
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = ( ( amax - a[i]        ) * bmin
             + (        a[i] - amin ) * bmax )
             / ( amax        - amin );
    }
  }

  return b;
}
//****************************************************************************80

float r4vec_amax ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_AMAX returns the maximum absolute value in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], the array.
//
//    Output, float AMAX, the value of the entry
//    of largest magnitude.
//
{
  float amax;
  int i;

  amax = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( amax < r4_abs ( a[i] ) )
    {
      amax = r4_abs ( a[i] );
    }
  }

  return amax;
}
//****************************************************************************80

int r4vec_amax_index ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_AMAX_INDEX returns the index of the maximum absolute value in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], the array.
//
//    Output, int R4VEC_AMAX_INDEX, the index of the entry of largest magnitude.
//
{
  float amax;
  int amax_index;
  int i;

  if ( n <= 0 )
  {
    amax_index = -1;
  }
  else
  {
    amax_index = 1;
    amax = r4_abs ( a[0] );

    for ( i = 2; i <= n; i++ )
    {
      if ( amax < r4_abs ( a[i-1] ) )
      {
        amax_index = i;
        amax = r4_abs ( a[i-1] );
      }
    }
  }

  return amax_index;
}
//****************************************************************************80

float r4vec_amin ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_AMIN returns the minimum absolute value in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], the array.
//
//    Output, float R4VEC_AMIN, the value of the entry
//    of smallest magnitude.
//
{
  float amin;
  int i;

  amin = r4_huge ( );
  for ( i = 0; i < n; i++ )
  {
    if ( r4_abs ( a[i] ) < amin )
    {
      amin = r4_abs ( a[i] );
    }
  }

  return amin;
}
//****************************************************************************80

int r4vec_amin_index ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_AMIN_INDEX returns the index of the minimum absolute value in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], the array.
//
//    Output, int R4VEC_AMIN_INDEX, the index of the entry of smallest magnitude.
//
{
  float amin;
  int amin_index;
  int i;

  if ( n <= 0 )
  {
    amin_index = -1;
  }
  else
  {
    amin_index = 1;
    amin = r4_abs ( a[0] );

    for ( i = 2; i <= n; i++ )
    {
      if ( r4_abs ( a[i-1] ) < amin )
      {
        amin_index = i;
        amin = r4_abs ( a[i-1] );
      }
    }
  }

  return amin_index;
}
//****************************************************************************80

float *r4vec_any_normal ( int dim_num, float v1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_ANY_NORMAL returns some normal vector to V1.
//
//  Discussion:
//
//    If DIM_NUM < 2, then no normal vector can be returned.
//
//    If V1 is the zero vector, then any unit vector will do.
//
//    No doubt, there are better, more robust algorithms.  But I will take
//    just about ANY reasonable unit vector that is normal to V1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, float V1[DIM_NUM], the vector.
//
//    Output, float R4VEC_ANY_NORMAL[DIM_NUM], a vector that is
//    normal to V2, and has unit Euclidean length.
//
{
  int i;
  int j;
  int k;
  float *v2;
  float vj;
  float vk;

  if ( dim_num < 2 )
  {
    cerr << "\n";
    cerr << "R4VEC_ANY_NORMAL - Fatal error!\n";
    cerr << "  Called with DIM_NUM < 2.\n";
    exit ( 1 );
  }

  v2 = new float[dim_num];

  if ( r4vec_norm ( dim_num, v1 ) == 0.0 )
  {
    r4vec_zero ( dim_num, v2 );
    v2[0] = 1.0;
    return v2;
  }
//
//  Seek the largest entry in V1, VJ = V1(J), and the
//  second largest, VK = V1(K).
//
//  Since V1 does not have zero norm, we are guaranteed that
//  VJ, at least, is not zero.
//
  j = -1;
  vj = 0.0;

  k = -1;
  vk = 0.0;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( r4_abs ( vk ) < r4_abs ( v1[i] ) || k == -1 )
    {
      if ( r4_abs ( vj ) < r4_abs ( v1[i] ) || j == -1 )
      {
        k = j;
        vk = vj;
        j = i;
        vj = v1[i];
      }
      else
      {
        k = i;
        vk = v1[i];
      }
    }
  }
//
//  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
//  will just about do the trick.
//
  r4vec_zero ( dim_num, v2 );

  v2[j] = -vk / sqrt ( vk * vk + vj * vj );
  v2[k] =  vj / sqrt ( vk * vk + vj * vj );

  return v2;
}
//****************************************************************************80

void r4vec_bracket ( int n, float x[], float xval, int *left,
  int *right )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_BRACKET searches a sorted array for successive brackets of a value.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    If the values in the vector are thought of as defining intervals
//    on the real line, then this routine searches for the interval
//    nearest to or containing the given value.
//
//    It is always true that RIGHT = LEFT+1.
//
//    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
//      XVAL   < X[0] < X[1];
//    If X(1) <= XVAL < X[N-1], then
//      X[LEFT-1] <= XVAL < X[RIGHT-1];
//    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
//      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
//
//    For consistency, this routine computes indices RIGHT and LEFT
//    that are 1-based, although it would be more natural in C and
//    C++ to use 0-based values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input, float X[N], an array that has been sorted into ascending order.
//
//    Input, float XVAL, a value to be bracketed.
//
//    Output, int *LEFT, *RIGHT, the results of the search.
//
{
  int i;

  for ( i = 2; i <= n - 1; i++ )
  {
    if ( xval < x[i-1] )
    {
      *left = i - 1;
      *right = i;
      return;
    }

   }

  *left = n - 1;
  *right = n;

  return;
}
//****************************************************************************80

void r4vec_bracket2 ( int n, float x[], float xval, int start, int *left,
  int *right )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_BRACKET2 searches a sorted array for successive brackets of a value.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    If the values in the vector are thought of as defining intervals
//    on the real line, then this routine searches for the interval
//    containing the given value.
//
//    R4VEC_BRACKET2 is a variation on R4VEC_BRACKET.  It seeks to reduce
//    the search time by allowing the user to suggest an interval that
//    probably contains the value.  The routine will look in that interval
//    and the intervals to the immediate left and right.  If this does
//    not locate the point, a binary search will be carried out on
//    appropriate subportion of the sorted array.
//
//    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
//    and X(LEFT) <= XVAL <= X(RIGHT).
//
//    Special cases:
//      Value is less than all data values:
//    LEFT = -1, RIGHT = 1, and XVAL < X(RIGHT).
//      Value is greater than all data values:
//    LEFT = N, RIGHT = -1, and X(LEFT) < XVAL.
//      Value is equal to a data value:
//    LEFT = RIGHT, and X(LEFT) = X(RIGHT) = XVAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of the input array.
//
//    Input, float X[N], an array that has been sorted into
//    ascending order.
//
//    Input, float XVAL, a value to be bracketed by entries of X.
//
//    Input, int START, between 1 and N, specifies that XVAL
//    is likely to be in the interval:
//      [ X(START), X(START+1) ]
//    or, if not in that interval, then either
//      [ X(START+1), X(START+2) ]
//    or
//      [ X(START-1), X(START) ].
//
//    Output, int *LEFT, *RIGHT, the results of the search.
//
{
  int high;
  int low;
//
//  Check.
//
  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R4VEC_BRACKET2 - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( start < 1 || n < start )
  {
    start = ( n + 1 ) / 2;
  }
//
//  XVAL = X(START)?
//
  if ( x[start-1] == xval )
  {
    *left = start;
    *right = start;
    return;
  }
//
//  X(START) < XVAL?
//
  else if ( x[start-1] < xval )
  {
//
//  X(START) = X(N) < XVAL < oo?
//
    if ( n < start + 1 )
    {
      *left = start;
      *right = -1;
      return;
    }
//
//  XVAL = X(START+1)?
//
    else if ( xval == x[start] )
    {
      *left = start + 1;
      *right = start + 1;
      return;
    }
//
//  X(START) < XVAL < X(START+1)?
//
    else if ( xval < x[start] )
    {
      *left = start;
      *right = start + 1;
      return;
    }
//
//  X(START+1) = X(N) < XVAL < oo?
//
    else if ( n < start + 2 )
    {
      *left = start + 1;
      *right = -1;
      return;
    }
//
//  XVAL = X(START+2)?
//
    else if ( xval == x[start+1] )
    {
      *left = start + 2;
      *right = start + 2;
      return;
    }
//
//  X(START+1) < XVAL < X(START+2)?
//
    else if ( xval < x[start+1] )
    {
      *left = start + 1;
      *right = start + 2;
      return;
    }
//
//  Binary search for XVAL in [ X(START+2), X(N) ],
//  where XVAL is guaranteed to be greater than X(START+2).
//
    else
    {
      low = start + 2;
      high = n;

      r4vec_bracket ( high + 1 - low, x+low-1, xval, left, right );

      *left = *left + low - 1;
      *right = *right + low - 1;
    }
  }
//
//  -oo < XVAL < X(START) = X(1).
//
  else if ( start == 1 )
  {
    *left = -1;
    *right = start;
    return;
  }
//
//  XVAL = X(START-1)?
//
  else if ( xval == x[start-2] )
  {
    *left = start - 1;
    *right = start - 1;
    return;
  }
//
//  X(START-1) < XVAL < X(START)?
//
  else if ( x[start-2] <= xval )
  {
    *left = start - 1;
    *right = start;
    return;
  }
//
//  Binary search for XVAL in [ X(1), X(START-1) ],
//  where XVAL is guaranteed to be less than X(START-1).
//
  else
  {
    low = 1;
    high = start - 1;
    r4vec_bracket ( high + 1 - low, x, xval, left, right );
  }

  return;
}
//****************************************************************************80

void r4vec_bracket3 ( int n, float t[], float tval, int *left )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_BRACKET3 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[N-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    This version of the function has been revised so that the value of
//    LEFT that is returned uses the 0-based indexing natural to C++.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of the input array.
//
//    Input, float T[N], an array that has been sorted into ascending order.
//
//    Input, float TVAL, a value to be bracketed by entries of T.
//
//    Input/output, int *LEFT.
//    On input, if 0 <= LEFT <= N-2, LEFT is taken as a suggestion for the
//    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
//    is searched first, followed by the appropriate interval to the left
//    or right.  After that, a binary search is used.
//    On output, LEFT is set so that the interval [ T[LEFT], T[LEFT+1] ]
//    is the closest to TVAL; it either contains TVAL, or else TVAL
//    lies outside the interval [ T[0], T[N-1] ].
//
{
  int high;
  int low;
  int mid;
//
//  Check the input data.
//
  if ( n < 2 )
  {
    cerr << "\n";
    cerr << "R4VEC_BRACKET3 - Fatal error!\n";
    cerr << "  N must be at least 2.\n";
    exit ( 1 );
  }
//
//  If *LEFT is not between 0 and N-2, set it to the middle value.
//
  if ( *left < 0 || n - 2 < *left )
  {
    *left = ( n - 1 ) / 2;
  }
//
//  CASE 1: TVAL < T[*LEFT]:
//  Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
//
  if ( tval < t[*left] )
  {
    if ( *left == 0 )
    {
      return;
    }
    else if ( *left == 1 )
    {
      *left = 0;
      return;
    }
    else if ( t[*left-1] <= tval )
    {
      *left = *left - 1;
      return;
    }
    else if ( tval <= t[1] )
    {
      *left = 0;
      return;
    }
//
//  ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
//
    low = 1;
    high = *left - 2;

    for ( ; ; )
    {
      if ( low == high )
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval )
      {
        low = mid;
      }
      else
      {
        high = mid - 1;
      }
    }
  }
//
//  CASE 2: T[*LEFT+1] < TVAL:
//  Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
//
  else if ( t[*left+1] < tval )
  {
    if ( *left == n - 2 )
    {
      return;
    }
    else if ( *left == n - 3 )
    {
      *left = *left + 1;
      return;
    }
    else if ( tval <= t[*left+2] )
    {
      *left = *left + 1;
      return;
    }
    else if ( t[n-2] <= tval )
    {
      *left = n - 2;
      return;
    }
//
//  ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
//
    low = *left + 2;
    high = n - 3;

    for ( ; ; )
    {

      if ( low == high )
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval )
      {
        low = mid;
      }
      else
      {
        high = mid - 1;
      }
    }
  }
//
//  CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
//  T is just where the user said it might be.
//
  else
  {
  }

  return;
}
//****************************************************************************80

void r4vec_bracket4 ( int nt, float t[], int ns, float s[], int left[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_BRACKET4 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[NT-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    This version of the function has been revised so that the value of
//    LEFT that is returned uses the 0-based indexing natural to C++.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, length of the input array.
//
//    Input, float T[NT], an array that has been sorted
//    into ascending order.
//
//    Input, int NS, the number of points to be bracketed.
//
//    Input, float S[NS], values to be bracketed by entries of T.
//
//    Output, int LEFT[NS].
//    LEFT[I] is set so that the interval [ T[LEFT[I]], T[LEFT[I]+1] ]
//    is the closest to S[I]; it either contains S[I], or else S[I]
//    lies outside the interval [ T[0], T[NT-1] ].
//
{
  int high;
  int i;
  int low;
  int mid;
//
//  Check the input data.
//
  if ( nt < 2 )
  {
    cerr << "\n";
    cerr << "R4VEC_BRACKET4 - Fatal error!\n";
    cerr << "  NT must be at least 2.\n";
    exit ( 1 );
  }

  for ( i = 0; i < ns; i++ )
  {
    left[i] = ( nt - 1 ) / 2;
//
//  CASE 1: S[I] < T[LEFT]:
//  Search for S[I] in (T[I],T[I+1]), for I = 0 to LEFT-1.
//
    if ( s[i] < t[left[i]] )
    {
      if ( left[i] == 0 )
      {
        continue;
      }
      else if ( left[i] == 1 )
      {
        left[i] = 0;
        continue;
      }
      else if ( t[left[i]-1] <= s[i] )
      {
        left[i] = left[i] - 1;
        continue;
      }
      else if ( s[i] <= t[1] )
      {
        left[i] = 0;
        continue;
      }
//
//  ...Binary search for S[I] in (T[I],T[I+1]), for I = 1 to *LEFT-2.
//
      low = 1;
      high = left[i] - 2;

      for ( ; ; )
      {
        if ( low == high )
        {
          left[i] = low;
          break;
        }

        mid = ( low + high + 1 ) / 2;

        if ( t[mid] <= s[i] )
        {
          low = mid;
        }
        else
        {
          high = mid - 1;
        }
      }
    }
//
//  CASE 2: T[LEFT+1] < S[I]:
//  Search for S[I] in (T[I],T[I+1]) for intervals I = LEFT+1 to NT-2.
//
    else if ( t[left[i]+1] < s[i] )
    {
      if ( left[i] == nt - 2 )
      {
        continue;
      }
      else if ( left[i] == nt - 3 )
      {
        left[i] = left[i] + 1;
        continue;
      }
      else if ( s[i] <= t[left[i]+2] )
      {
        left[i] = left[i] + 1;
        continue;
      }
      else if ( t[nt-2] <= s[i] )
      {
        left[i] = nt - 2;
        continue;
      }
//
//  ...Binary search for S[I] in (T[I],T[I+1]) for intervals I = LEFT+2 to NT-3.
//
      low = left[i] + 2;
      high = nt - 3;

      for ( ; ; )
      {

        if ( low == high )
        {
          left[i] = low;
          break;
        }

        mid = ( low + high + 1 ) / 2;

        if ( t[mid] <= s[i] )
        {
          low = mid;
        }
        else
        {
          high = mid - 1;
        }
      }
    }
//
//  CASE 3: T[LEFT] <= S[I] <= T[LEFT+1]:
//
    else
    {
    }
  }
  return;
}
//****************************************************************************80

float r4vec_circular_variance ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_CIRCULAR_VARIANCE returns the circular variance of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[N], the vector whose variance is desired.
//
//    Output, float R4VEC_CIRCULAR VARIANCE, the circular variance
//    of the vector entries.
//
{
  int i;
  float mean;
  float sum_c;
  float sum_s;
  float value;

  mean = r4vec_mean ( n, x );

  sum_c = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum_c = sum_c + cos ( x[i] - mean );
  }

  sum_s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum_s = sum_s + sin ( x[i] - mean );
  }

  value = sqrt ( sum_c * sum_c + sum_s * sum_s ) / ( float ) n;

  value = 1.0 - value;

  return value;
}
//****************************************************************************80

int r4vec_compare ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_COMPARE compares two R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The lexicographic ordering is used.
//
//  Example:
//
//    Input:
//
//      A1 = ( 2.0, 6.0, 2.0 )
//      A2 = ( 2.0, 8.0, 12.0 )
//
//    Output:
//
//      ISGN = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, float A[N], B[N], the vectors to be compared.
//
//    Output, int R4VEC_COMPARE, the results of the comparison:
//    -1, A is lexicographically less than B,
//     0, A is equal to B,
//    +1, A is lexicographically greater than B.
//
{
  int isgn;
  int k;

  isgn = 0;

  for ( k = 0; k < n; k++ )
  {
    if ( a[k] < b[k] )
    {
      isgn = -1;
      return isgn;
    }
    else if ( b[k] < a[k] )
    {
      isgn = +1;
      return isgn;
    }
  }
  return isgn;
}
//****************************************************************************80

float *r4vec_convolve_circ ( int n, float x[], float y[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_CONVOLVE_CIRC returns the discrete circular convolution of two R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    z(1+m) = xCCy(m) = sum ( 0 <= k <= n-1 ) x(1+k) * y(1+m-k)
//
//    Here, if the index of Y becomes nonpositive, it is "wrapped around"
//    by having N added to it.
//
//    The circular convolution is equivalent to multiplication of Y by a
//    circulant matrix formed from the vector X.
//
//  Example:
//
//    Input:
//
//      X = (/ 1, 2, 3, 4 /)
//      Y = (/ 1, 2, 4, 8 /)
//
//    Output:
//
//      Circulant form:
//
//      Z = ( 1 4 3 2 )   ( 1 )
//          ( 2 1 4 3 )   ( 2 )
//          ( 3 2 1 4 ) * ( 4 )
//          ( 4 3 2 1 )   ( 8 )
//
//      The formula:
//
//      Z = (/ 1*1 + 2*8 + 3*4 + 4*2,
//             1*2 + 2*1 + 3*8 + 4*4,
//             1*4 + 2*2 + 3*1 + 4*8,
//             1*8 + 2*4 + 3*2 + 4*1 /)
//
//      Result:
//
//      Z = (/ 37, 44, 43, 26 /)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, float X[N], Y[N], the vectors to be convolved.
//
//    Output, float R4VEC_CONVOLVE_CIRC[N], the circular convolution of X and Y.
//
{
  int i;
  int m;
  float *z;

  z = new float[n];

  for ( m = 1; m <= n; m++ )
  {
    z[m-1] = 0.0;
    for ( i = 1; i <= m; i++ )
    {
      z[m-1] = z[m-1] + x[i-1] * y[m-i];
    }
    for ( i = m+1; i <= n; i++ )
    {
      z[m-1] = z[m-1] + x[i-1] * y[n+m-i];
    }
  }

  return z;
}
//****************************************************************************80

void r4vec_copy ( int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_COPY copies an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, float A1[N], the vector to be copied.
//
//    Output, float A2[N], the copy of A1.
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

float *r4vec_copy_new ( int n, float a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_COPY_NEW copies an R4VEC to a "new" R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, float A1[N], the vector to be copied.
//
//    Output, float R4VEC_COPY_NEW[N], the copy of A1.
//
{
  float *a2;
  int i;

  a2 = new float[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

float r4vec_correlation ( int n, float x[], float y[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_CORRELATION returns the correlation of two R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    If X and Y are two nonzero vectors of length N, then
//
//      correlation = (x/||x||)' (y/||y||)
//
//    It is the cosine of the angle between the two vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, float X[N], Y[N], the vectors to be convolved.
//
//    Output, float R4VEC_CORRELATION, the correlation of X and Y.
//
{
  float correlation;
  float x_norm;
  float xy_dot;
  float y_norm;

  x_norm = r4vec_norm ( n, x );
  y_norm = r4vec_norm ( n, y );
  xy_dot = r4vec_dot_product ( n, x, y );

  if ( x_norm == 0.0 || y_norm == 0.0 )
  {
    correlation = 0.0;
  }
  else
  {
    correlation = xy_dot / x_norm / y_norm;
  }

  return correlation;
}
//****************************************************************************80

float r4vec_covar ( int n, float x[], float y[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_COVAR computes the covariance of two vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 April 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the dimension of the two vectors.
//
//    Input, float X[N], Y[N], the two vectors.
//
//    Output, float R4VEC_COVAR, the covariance of the two vectors.
//
{
  int i;
  float value;
  float x_average;
  float y_average;

  x_average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    x_average = x_average + x[i];
  }
  x_average = x_average / ( float ) ( n );

  y_average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    y_average = y_average + x[i];
  }
  y_average = y_average / ( float ) ( n );

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + ( x[i] - x_average ) * ( y[i] - y_average );
  }

  value = value / ( float ) ( n - 1 );

  return value;
}
//****************************************************************************80

float r4vec_cross_product_2d ( float v1[2], float v2[2] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_CROSS_PRODUCT_2D finds the cross product of a pair of R4VEC's in 2D.
//
//  Discussion:
//
//    Strictly speaking, the vectors lie in the (X,Y) plane, and
//    the cross product here is a vector in the Z direction.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float V1[2], V2[2], the vectors.
//
//    Output, float R4VEC_CROSS_PRODUCT_2D, the Z component of the cross product
//    of V1 and V2.
//
{
  float value;

  value = v1[0] * v2[1] - v1[1] * v2[0];

  return value;
}
//****************************************************************************80

float r4vec_cross_product_affine_2d ( float v0[2], float v1[2],
  float v2[2] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_CROSS_PRODUCT_AFFINE_2D finds the affine cross product in 2D.
//
//  Discussion:
//
//    Strictly speaking, the vectors lie in the (X,Y) plane, and
//    the cross product here is a vector in the Z direction.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float V0[2], the base vector.
//
//    Input, float V1[2], V2[2], the vectors.
//
//    Output, float R4VEC_CROSS_PRODUCT_AFFINE_2D, the Z component of the
//    cross product of V1 and V2.
//
{
  float value;

  value =
      ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
    - ( v2[0] - v0[0] ) * ( v1[1] - v0[1] );

  return value;
}
//****************************************************************************80

float *r4vec_cross_product_3d ( float v1[3], float v2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_CROSS_PRODUCT_3D computes the cross product of two R4VEC's in 3D.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float V1[3], V2[3], the coordinates of the vectors.
//
//    Output, float R4VEC_CROSS_PRODUCT_3D[3], the cross product vector.
//
{
  float *v3;

  v3 = new float[3];

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}
//****************************************************************************80

float *r4vec_cross_product_affine_3d ( float v0[3], float v1[3],
  float v2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_CROSS_PRODUCT_AFFINE_3D computes the affine cross product in 3D.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float V0[3], the base vector.
//
//    Input, float V1[3], V2[3], the coordinates of the vectors.
//
//    Output, float R4VEC_CROSS_PRODUCT_AFFINE_3D[3], the cross product vector.
//
{
  float *v3;

  v3 = ( float * ) malloc ( 3 * sizeof ( float ) );

  v3[0] =
      ( v1[1] - v0[1] ) * ( v2[2] - v0[2] )
    - ( v2[1] - v0[1] ) * ( v1[2] - v0[2] );

  v3[1] =
      ( v1[2] - v0[2] ) * ( v2[0] - v0[0] )
    - ( v2[2] - v0[2] ) * ( v1[0] - v0[0] );

  v3[2] =
      ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
    - ( v2[0] - v0[0] ) * ( v1[1] - v0[1] );

  return v3;
}
//****************************************************************************80

float *r4vec_dif ( int n, float h )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIF computes coefficients for estimating the N-th derivative.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The routine computes the N+1 coefficients for a centered finite difference
//    estimate of the N-th derivative of a function.
//
//    The estimate has the form
//
//      FDIF(N,X) = Sum (I = 0 to N) COF(I) * F ( X(I) )
//
//    To understand the computation of the coefficients, it is enough
//    to realize that the first difference approximation is
//
//      FDIF(1,X) = F(X+DX) - F(X-DX) ) / (2*DX)
//
//    and that the second difference approximation can be regarded as
//    the first difference approximation repeated:
//
//      FDIF(2,X) = FDIF(1,X+DX) - FDIF(1,X-DX) / (2*DX)
//         = F(X+2*DX) - 2 F(X) + F(X-2*DX) / (4*DX)
//
//    and so on for higher order differences.
//
//    Thus, the next thing to consider is the integer coefficients of
//    the sampled values of F, which are clearly the Pascal coefficients,
//    but with an alternating negative sign.  In particular, if we
//    consider row I of Pascal's triangle to have entries j = 0 through I,
//    then P(I,J) = P(I-1,J-1) - P(I-1,J), where P(*,-1) is taken to be 0,
//    and P(0,0) = 1.
//
//       1
//      -1  1
//       1 -2   1
//      -1  3  -3   1
//       1 -4   6  -4   1
//      -1  5 -10  10  -5  1
//       1 -6  15 -20  15 -6 1
//
//    Next, note that the denominator of the approximation for the
//    N-th derivative will be (2*DX)**N.
//
//    And finally, consider the location of the N+1 sampling
//    points for F:
//
//      X-N*DX, X-(N-2)*DX, X-(N-4)*DX, ..., X+(N-4)*DX, X+(N-2*DX), X+N*DX.
//
//    Thus, a formula for evaluating FDIF(N,X) is
//
//      fdif = 0.0
//      do i = 0, n
//        xi = x + (2*i-n) * h
//        fdif = fdif + cof(i) * f(xi)
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the derivative to be approximated.
//    N must be 0 or greater.
//
//    Input, float H, the half spacing between points.
//    H must be positive.
//
//    Output, float R4VEC_DIF[N+1], the coefficients needed to approximate
//    the N-th derivative of a function F.
//
{
  float *cof;
  int i;
  int j;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_DIF - Fatal error!\n";
    cerr << "  Derivative order N = " << n << "\n";
    cerr << "  but N must be at least 0.\n";
    exit ( 1 );
  }

  if ( h <= 0.0 )
  {
    cerr << "\n";
    cerr << "R4VEC_DIF - Fatal error!\n";
    cerr << "  The half sampling spacing is H = " << h << "\n";
    cerr << "  but H must be positive.\n";
    exit ( 1 );
  }

  cof = new float[n+1];

  for ( i = 0; i <= n; i++ )
  {
    cof[i] = 1.0;

    for ( j = i - 1; 1 <= j; j-- )
    {
      cof[j] = -cof[j] + cof[j-1];
    }

    if ( 0 < i )
    {
      cof[0] = - cof[0];
    }
  }

  for ( i = 0; i <= n; i++ )
  {
    cof[i] = cof[i] / pow ( 2.0 * h, n );
  }

  return cof;
}
//****************************************************************************80

float r4vec_diff_norm ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIFF_NORM returns the L2 norm of the difference of R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L2 norm is defined as:
//
//      R4VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], B[N], the vectors.
//
//    Output, float R4VEC_DIFF_NORM, the L2 norm of A - B.
//
{
  int i;
  float value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

float r4vec_diff_norm_l1 ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIFF_NORM_L1 returns the L1 norm of the difference of R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L1 norm is defined as:
//
//      R4VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], B[N], the vectors.
//
//    Output, float R4VEC_DIFF_NORM_L1, the L1 norm of A - B.
//
{
  int i;
  float value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + r4_abs ( a[i] - b[i] );
  }
  return value;
}
//****************************************************************************80

float r4vec_diff_norm_l2 ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIFF_NORM_L2 returns the L2 norm of the difference of R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L2 norm is defined as:
//
//      R4VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], B[N], the vectors.
//
//    Output, float R4VEC_DIFF_NORM_L2, the L2 norm of A - B.
//
{
  int i;
  float value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

float r4vec_diff_norm_li ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIFF_NORM_LI returns the L-oo norm of the difference of R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L-oo norm is defined as:
//
//      R4VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], B[N], the vectors.
//
//    Output, float R4VEC_DIFF_NORM_LI, the L-oo norm of A - B.
//
{
  int i;
  float value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = r4_max ( value, r4_abs ( a[i] - b[i] ) );
  }
  return value;
}
//****************************************************************************80

float r4vec_diff_norm_squared ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIFF_NORM_L2 returns the square of the L2 norm of the difference of R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L2 norm is defined as:
//
//      R4VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], B[N], the vectors.
//
//    Output, float R4VEC_DIFF_NORM_SQUARED, the L2 norm of A - B.
//
{
  int i;
  float value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }

  return value;
}
//****************************************************************************80

void r4vec_direct_product ( int factor_index, int factor_order,
  float factor_value[], int factor_num, int point_num, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIRECT_PRODUCT creates a direct product of R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    To explain what is going on here, suppose we had to construct
//    a multidimensional quadrature rule as the product of K rules
//    for 1D quadrature.
//
//    The product rule will be represented as a list of points and weights.
//
//    The J-th item in the product rule will be associated with
//      item J1 of 1D rule 1,
//      item J2 of 1D rule 2,
//      ...,
//      item JK of 1D rule K.
//
//    In particular,
//      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
//    and
//      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
//
//    So we can construct the quadrature rule if we can properly
//    distribute the information in the 1D quadrature rules.
//
//    This routine carries out that task.
//
//    Another way to do this would be to compute, one by one, the
//    set of all possible indices (J1,J2,...,JK), and then index
//    the appropriate information.  An advantage of the method shown
//    here is that you can process the K-th set of information and
//    then discard it.
//
//  Example:
//
//    Rule 1:
//      Order = 4
//      X(1:4) = ( 1, 2, 3, 4 )
//
//    Rule 2:
//      Order = 3
//      X(1:3) = ( 10, 20, 30 )
//
//    Rule 3:
//      Order = 2
//      X(1:2) = ( 100, 200 )
//
//    Product Rule:
//      Order = 24
//      X(1:24) =
//        ( 1, 10, 100 )
//        ( 2, 10, 100 )
//        ( 3, 10, 100 )
//        ( 4, 10, 100 )
//        ( 1, 20, 100 )
//        ( 2, 20, 100 )
//        ( 3, 20, 100 )
//        ( 4, 20, 100 )
//        ( 1, 30, 100 )
//        ( 2, 30, 100 )
//        ( 3, 30, 100 )
//        ( 4, 30, 100 )
//        ( 1, 10, 200 )
//        ( 2, 10, 200 )
//        ( 3, 10, 200 )
//        ( 4, 10, 200 )
//        ( 1, 20, 200 )
//        ( 2, 20, 200 )
//        ( 3, 20, 200 )
//        ( 4, 20, 200 )
//        ( 1, 30, 200 )
//        ( 2, 30, 200 )
//        ( 3, 30, 200 )
//        ( 4, 30, 200 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR_INDEX, the index of the factor being processed.
//    The first factor processed must be factor 0.
//
//    Input, int FACTOR_ORDER, the order of the factor.
//
//    Input, float FACTOR_VALUE[FACTOR_ORDER], the factor values
//    for factor FACTOR_INDEX.
//
//    Input, int FACTOR_NUM, the number of factors.
//
//    Input, int POINT_NUM, the number of elements in the direct product.
//
//    Input/output, float X[FACTOR_NUM*POINT_NUM], the elements of the
//    direct product, which are built up gradually.
//
//  Local Parameters:
//
//    Local, int START, the first location of a block of values to set.
//
//    Local, int CONTIG, the number of consecutive values to set.
//
//    Local, int SKIP, the distance from the current value of START
//    to the next location of a block of values to set.
//
//    Local, int REP, the number of blocks of values to set.
//
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( j = 0; j < point_num; j++ )
    {
      for ( i = 0; i < factor_num; i++ )
      {
        x[i+j*factor_num] = 0.0;
      }
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( i = 0; i < factor_order; i++ )
  {
    start = 0 + i * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( j = start; j < start + contig; j++ )
      {
        x[factor_index+j*factor_num] = factor_value[i];
      }
      start = start + skip;
    }
  }
  contig = contig * factor_order;

  return;
}
//****************************************************************************80

void r4vec_direct_product2 ( int factor_index, int factor_order,
  float factor_value[], int factor_num, int point_num, float w[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIRECT_PRODUCT2 creates a direct product of R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    To explain what is going on here, suppose we had to construct
//    a multidimensional quadrature rule as the product of K rules
//    for 1D quadrature.
//
//    The product rule will be represented as a list of points and weights.
//
//    The J-th item in the product rule will be associated with
//      item J1 of 1D rule 1,
//      item J2 of 1D rule 2,
//      ...,
//      item JK of 1D rule K.
//
//    In particular,
//      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
//    and
//      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
//
//    So we can construct the quadrature rule if we can properly
//    distribute the information in the 1D quadrature rules.
//
//    This routine carries out that task for the weights W.
//
//    Another way to do this would be to compute, one by one, the
//    set of all possible indices (J1,J2,...,JK), and then index
//    the appropriate information.  An advantage of the method shown
//    here is that you can process the K-th set of information and
//    then discard it.
//
//  Example:
//
//    Rule 1:
//      Order = 4
//      W(1:4) = ( 2, 3, 5, 7 )
//
//    Rule 2:
//      Order = 3
//      W(1:3) = ( 11, 13, 17 )
//
//    Rule 3:
//      Order = 2
//      W(1:2) = ( 19, 23 )
//
//    Product Rule:
//      Order = 24
//      W(1:24) =
//        ( 2 * 11 * 19 )
//        ( 3 * 11 * 19 )
//        ( 4 * 11 * 19 )
//        ( 7 * 11 * 19 )
//        ( 2 * 13 * 19 )
//        ( 3 * 13 * 19 )
//        ( 5 * 13 * 19 )
//        ( 7 * 13 * 19 )
//        ( 2 * 17 * 19 )
//        ( 3 * 17 * 19 )
//        ( 5 * 17 * 19 )
//        ( 7 * 17 * 19 )
//        ( 2 * 11 * 23 )
//        ( 3 * 11 * 23 )
//        ( 5 * 11 * 23 )
//        ( 7 * 11 * 23 )
//        ( 2 * 13 * 23 )
//        ( 3 * 13 * 23 )
//        ( 5 * 13 * 23 )
//        ( 7 * 13 * 23 )
//        ( 2 * 17 * 23 )
//        ( 3 * 17 * 23 )
//        ( 5 * 17 * 23 )
//        ( 7 * 17 * 23 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR_INDEX, the index of the factor being processed.
//    The first factor processed must be factor 0.
//
//    Input, int FACTOR_ORDER, the order of the factor.
//
//    Input, float FACTOR_VALUE[FACTOR_ORDER], the factor values for
//    factor FACTOR_INDEX.
//
//    Input, int FACTOR_NUM, the number of factors.
//
//    Input, int POINT_NUM, the number of elements in the direct product.
//
//    Input/output, float W[POINT_NUM], the elements of the
//    direct product, which are built up gradually.
//
//  Local Parameters:
//
//    Local, integer START, the first location of a block of values to set.
//
//    Local, integer CONTIG, the number of consecutive values to set.
//
//    Local, integer SKIP, the distance from the current value of START
//    to the next location of a block of values to set.
//
//    Local, integer REP, the number of blocks of values to set.
//
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 1.0;
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        w[i] = w[i] * factor_value[j];
      }
      start = start + skip;
    }
  }

  contig = contig * factor_order;

  return;
}
//****************************************************************************80

float r4vec_distance ( int dim_num, float v1[], float v2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DISTANCE returns the Euclidean distance between two R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, float V1[DIM_NUM], V2[DIM_NUM], the vectors.
//
//    Output, float R4VEC_DISTANCE, the Euclidean distance
//    between the vectors.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = pow ( v1[i] - v2[i], 2 );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

bool r4vec_distinct ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DISTINCT is true if the entries in an R4VEC are distinct.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[N], the vector to be checked.
//
//    Output, bool R4VEC_DISTINCT is true if all N elements of X
//    are distinct.
//
{
  int i;
  int j;

  for ( i = 1; i <= n-1; i++ )
  {
    for ( j = 1; j <= i - 1; j++ )
    {
      if ( x[i] == x[j] )
      {
        return false;
      }
    }
  }
  return true;
}
//****************************************************************************80

void r4vec_divide ( int n, float a[], float s )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DIVIDE divides an R4VEC by a nonzero scalar.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, float A[N].  On input, the vector to be scaled.
//    On output, each entry has been divided by S.
//
//    Input, float S, the divisor.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / s;
  }
  return;
}
//****************************************************************************80

float r4vec_dot_product ( int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DOT_PRODUCT computes the dot product of a pair of R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, float A1[N], A2[N], the two vectors to be considered.
//
//    Output, float R4VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

float r4vec_dot_product_affine ( int n, float v0[], float v1[], float v2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_DOT_PRODUCT_AFFINE computes the affine dot product.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, float V0[N], the base vector.
//
//    Input, float V1[N], V2[N], the two vectors to be considered.
//
//    Output, float R4VEC_DOT_PRODUCT_AFFINE, the dot product of the vectors.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v2[i] - v0[i] );
  }
  return value;
}
//****************************************************************************80

bool r4vec_eq ( int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_EQ is true if every pair of entries in two R4VEC's is equal.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, float A1[N], A2[N], two vectors to compare.
//
//    Output, bool R4VEC_EQ, is TRUE if every pair of elements A1(I)
//    and A2(I) are equal, and FALSE otherwise.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

void r4vec_even ( int n, float alo, float ahi, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_EVEN returns an R4VEC of values evenly spaced between ALO and AHI.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, float ALO, AHI, the low and high values.
//
//    Output, float A[N], N evenly spaced values.
//    Normally, A[0] = ALO and A[N-1] = AHI.
//    However, if N = 1, then A[0] = 0.5*(ALO+AHI).
//
{
  int i;

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( float ) ( n - i - 1 ) * alo
             + ( float ) (     i     ) * ahi )
             / ( float ) ( n     - 1 );
    }
  }

  return;
}
//****************************************************************************80

float *r4vec_even_new ( int n, float alo, float ahi )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_EVEN_NEW returns an R4VEC of values evenly spaced between ALO and AHI.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, float ALO, AHI, the low and high values.
//
//    Output, float R4VEC_EVEN_NEW[N], N evenly spaced values.
//    Normally, A[0] = ALO and A[N-1] = AHI.
//    However, if N = 1, then A[0] = 0.5*(ALO+AHI).
//
{
  float *a;
  int i;

  a = new float[n];

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( float ) ( n - i - 1 ) * alo
             + ( float ) (     i     ) * ahi )
             / ( float ) ( n     - 1 );
    }
  }

  return a;
}
//****************************************************************************80

float r4vec_even_select ( int n, float xlo, float xhi, int ival )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / ( N - 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, float XLO, XHI, the low and high values.
//
//    Input, int IVAL, the index of the desired point.
//    IVAL is normally between 1 and N, but may be any integer value.
//
//    Output, float R4VEC_EVEN_SELECT, the IVAL-th of N evenly spaced values
//    between XLO and XHI.
//    Unless N = 1, X(1) = XLO and X(N) = XHI.
//    If N = 1, then X(1) = 0.5*(XLO+XHI).
//
{
  float xval;

  if ( n == 1 )
  {
    xval = 0.5 * ( xlo + xhi );
  }
  else
  {
    xval = ( ( float ) ( n - ival     ) * xlo
           + ( float ) (     ival - 1 ) * xhi )
           / ( float ) ( n        - 1 );
  }

  return xval;
}
//****************************************************************************80

void r4vec_even2 ( int maxval, int nfill[], int nold, float xold[],
  int *nval, float xval[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_EVEN2 linearly interpolates new numbers into an R4VECa.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The number of values created between two old values can vary from
//    one pair of values to the next.
//
//    The interpolated values are evenly spaced.
//
//    This routine is a generalization of R4VEC_EVEN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MAXVAL, the size of the XVAL array, as declared by the
//    user.  MAXVAL must be large enough to hold the NVAL values computed by
//    this routine.  In other words, MAXVAL must be at least equal to
//    NOLD + SUM (1 <= I <= NOLD-1) NFILL(I).
//
//    Input, int NFILL[NOLD-1], the number of values
//    to be interpolated between XOLD(I) and XOLD(I+1).
//    NFILL(I) does not count the endpoints.  Thus, if
//    NFILL(I) is 1, there will be one new point generated
//    between XOLD(I) and XOLD(I+1).
//    NFILL(I) must be nonnegative.
//
//    Input, int NOLD, the number of values XOLD,
//    between which extra values are to be interpolated.
//
//    Input, float XOLD[NOLD], the original vector of numbers
//    between which new values are to be interpolated.
//
//    Output, int *NVAL, the number of values computed
//    in the XVAL array.
//    NVAL = NOLD + SUM ( 1 <= I <= NOLD-1 ) NFILL(I)
//
//    Output, float XVAL[MAXVAL].  On output, XVAL contains the
//    NOLD values of XOLD, as well as the interpolated
//    values, making a total of NVAL values.
//
{
  int i;
  int j;
  int nadd;

  *nval = 1;

  for ( i = 1; i <= nold - 1; i++ )
  {

    if ( nfill[i-1] < 0 )
    {
      cerr << "\n";
      cerr << "R4VEC_EVEN2 - Fatal error!\n";
      cerr << "  NFILL[I-1] is negative for I = " << i << "\n";
      cerr << "  NFILL[I-1] = " << nfill[i-1] << "\n";
      exit ( 1 );
    }

    if ( maxval < *nval + nfill[i-1] + 1 )
    {
      cerr << "\n";
      cerr << "R4VEC_EVEN2 - Fatal error!\n";
      cerr << "  MAXVAL = " << maxval << " is not large enough.\n";
      cerr << "  for the storage for interval I = " << i << "\n";
      exit ( 1 );
    }

    nadd = nfill[i-1] + 2;

    for ( j = 1; j <= nadd; j++ )
    {
      xval[*nval+j-2] = ( ( float ) ( nadd - j     ) * xold[i-1]
                        + ( float ) (        j - 1 ) * xold[i] )
                        / ( float ) ( nadd     - 1 );
    }

    *nval = *nval + nfill[i-1] + 1;
  }

  return;
}
//****************************************************************************80

void r4vec_even3 ( int nold, int nval, float xold[], float xval[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_EVEN3 evenly interpolates new data into an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    This routine accepts a short vector of numbers, and returns a longer
//    vector of numbers, created by interpolating new values between
//    the given values.
//
//    Between any two original values, new values are evenly interpolated.
//
//    Over the whole vector, the new numbers are interpolated in
//    such a way as to try to minimize the largest distance interval size.
//
//    The algorithm employed is not "perfect".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NOLD, the number of values XOLD, between which extra
//    values are to be interpolated.
//
//    Input, int NVAL, the number of values to be computed
//    in the XVAL array.  NVAL should be at least NOLD.
//
//    Input, float XOLD[NOLD], the original vector of numbers
//    between which new values are to be interpolated.
//
//    Output, float XVAL[NVAL].  On output, XVAL contains the
//    NOLD values of XOLD, as well as interpolated
//    values, making a total of NVAL values.
//
{
  float density;
  int i;
  int ival;
  int j;
  int nmaybe;
  int npts;
  int ntemp;
  int ntot;
  float xlen;
  float xleni;
  float xlentot;

  xlen = 0.0;
  for ( i = 1; i <= nold - 1; i++ )
  {
    xlen = xlen + r4_abs ( xold[i] - xold[i-1] );
  }

  ntemp = nval - nold;

  density = ( float ) ( ntemp ) / xlen;

  ival = 1;
  ntot = 0;
  xlentot = 0.0;

  for ( i = 1; i <= nold - 1; i++ )
  {
    xleni = r4_abs ( xold[i] - xold[i-1] );
    npts = ( int ) ( density * xleni );
    ntot = ntot + npts;
//
//  Determine if we have enough left-over density that it should
//  be changed into a point.  A better algorithm would agonize
//  more over where that point should go.
//
    xlentot = xlentot + xleni;
    nmaybe = r4_nint ( xlentot * density );

    if ( ntot < nmaybe )
    {
      npts = npts + nmaybe - ntot;
      ntot = nmaybe;
    }
    for ( j = 1; j <= npts + 2; j++ )
    {
      xval[ival+j-2] = ( ( float ) ( npts+2 - j     ) * xold[i-1]
                       + ( float ) (          j - 1 ) * xold[i] )
                       / ( float ) ( npts+2     - 1 );
    }
    ival = ival + npts + 1;
  }

  return;
}
//****************************************************************************80

float *r4vec_expand_linear ( int n, float x[], int fat )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_EXPAND_LINEAR linearly interpolates new data into an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of input data values.
//
//    Input, float X[N], the original data.
//
//    Input, int FAT, the number of data values to interpolate
//    between each pair of original data values.
//
//    Output, float R4VEC_EXPAND_LINEAR[(N-1)*(FAT+1)+1], the "fattened" data.
//
{
  int i;
  int j;
  int k;
  float *xfat;

  xfat = new float[(n-1)*(fat+1)+1];

  k = 0;

  for ( i = 0; i < n-1; i++ )
  {
    xfat[k] = x[i];
    k = k + 1;

    for ( j = 1; j <= fat; j++ )
    {
      xfat[k] = ( ( float ) ( fat - j + 1 ) * x[i]
                + ( float ) (       j     ) * x[i+1] )
                / ( float ) ( fat     + 1 );
      k = k + 1;
    }
  }

  xfat[k] = x[n-1];
  k = k + 1;

  return xfat;
}
//****************************************************************************80

int *r4vec_first_index ( int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_FIRST_INDEX indexes the first occurrence of values in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    For element A(I) of the vector, FIRST_INDEX(I) is the index in A of
//    the first occurrence of the value A(I).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, float A[N], the unsorted array to examine.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int R4VEC_FIRST_INDEX[N], the first occurrence index.
//
{
  int *first_index;
  int i;
  int j;

  first_index = new int[n];

  for ( i = 0; i < n; i++ )
  {
    first_index[i] = -1;
  }
  for ( i = 0; i < n; i++ )
  {
    if ( first_index[i] == -1 )
    {
      first_index[i] = i;
      for ( j = i + 1; j < n; j++ )
      {
        if ( r4_abs ( a[i] - a[j] ) <= tol )
        {
          first_index[j] = i;
        }
      }
    }
  }
  return first_index;
}
//****************************************************************************80

float r4vec_frac ( int n, float a[], int k )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_FRAC searches for the K-th smallest entry in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Hoare's algorithm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input/output, float A[N].
//    On input, A is the array to search.
//    On output, the elements of A have been somewhat rearranged.
//
//    Input, int K, the fractile to be sought.  If K = 1, the minimum
//    entry is sought.  If K = N, the maximum is sought.  Other values
//    of K search for the entry which is K-th in size.  K must be at
//    least 1, and no greater than N.
//
//    Output, float R4VEC_FRAC, the value of the K-th fractile of A.
//
{
  float frac;
  int i;
  int iryt;
  int j;
  int left;
  float temp;
  float x;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( k <= 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of K = " << k << "\n";
    exit ( 1 );
  }

  if ( n < k )
  {
    cerr << "\n";
    cerr << "R4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal N < K, K = " << k << "\n";
    exit ( 1 );
  }

  left = 1;
  iryt = n;

  for ( ; ; )
  {
    if ( iryt <= left )
    {
      frac = a[k-1];
      break;
    }

    x = a[k-1];
    i = left;
    j = iryt;

    for ( ; ; )
    {
      if ( j < i )
      {
        if ( j < k )
        {
          left = i;
        }
        if ( k < i )
        {
          iryt = j;
        }
        break;
      }
//
//  Find I so that X <= A(I).
//
      while ( a[i-1] < x )
      {
        i = i + 1;
      }
//
//  Find J so that A(J) <= X.
//
      while ( x < a[j-1] )
      {
        j = j - 1;
      }

      if ( i <= j )
      {
        temp   = a[i-1];
        a[i-1] = a[j-1];
        a[j-1] = temp;
        i = i + 1;
        j = j - 1;
      }
    }
  }

  return frac;
}
//****************************************************************************80

float *r4vec_fraction ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_FRACTION returns the fraction parts of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    If we regard a real number as
//
//      R4 = SIGN * ( WHOLE + FRACTION )
//
//    where
//
//      SIGN is +1 or -1,
//      WHOLE is a nonnegative integer
//      FRACTION is a nonnegative real number strictly less than 1,
//
//    then this routine returns the value of FRACTION.
//
//  Example:
//
//     R4    R4_FRACTION
//
//    0.00      0.00
//    1.01      0.01
//    2.02      0.02
//   19.73      0.73
//   -4.34      0.34
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, float X[N], the arguments.
//
//    Output, float R4_FRACTION[N], the fraction parts.
//
{
  float *fraction;
  int i;

  fraction = new float[n];

  for ( i = 0; i < n; i++ )
  {
    fraction[i] = r4_abs ( x[i] ) - ( float ) ( ( int ) ( r4_abs ( x[i] ) ) );
  }

  return fraction;
}
//****************************************************************************80

bool r4vec_gt ( int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_GT == ( A1 > A2 ) for two R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The comparison is lexicographic.
//
//    A1 > A2  <=>                              A1(1) > A2(1) or
//                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
//                 ...
//                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, float A1[N], A2[N], the vectors to be compared.
//
//    Output, bool R4VEC_GT, is TRUE if and only if A1 > A2.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {

    if ( a2[i] < a1[i] )
    {
       return true;
    }
    else if ( a1[i] < a2[i] )
    {
      return false;
    }

  }

  return false;
}
//****************************************************************************80

void r4vec_heap_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_HEAP_A reorders an R4VEC into a ascending heap.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    An ascending heap is an array A with the property that, for every index J,
//    A[J] <= A[2*J+1] and A[J] <= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//
//            A(1)         A(2)
//
//      A(3)       A(4)  A(5) A(6)
//
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
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
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, float A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//
{
  int i;
  int ifree;
  float key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m+1] < a[m] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( a[m] <= key )
        {
          break;
        }
        a[ifree] = a[m];
        ifree = m;
      }
    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;
  }

  return;
}
//****************************************************************************80

void r4vec_heap_d ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_HEAP_D reorders an R4VEC into a descending heap.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//
//            A(1)         A(2)
//
//      A(3)       A(4)  A(5) A(6)
//
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
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
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, float A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//
{
  int i;
  int ifree;
  float key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }
      }
    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;
  }

  return;
}
//****************************************************************************80

int *r4vec_histogram ( int n, float a[], float a_lo, float a_hi,
  int histo_num )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_HISTOGRAM histograms an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Values between A_LO and A_HI will be histogrammed into the bins
//    1 through HISTO_NUM.  Values below A_LO are counted in bin 0,
//    and values greater than A_HI are counted in bin HISTO_NUM+1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, float A[N], the array to examine.
//
//    Input, float A_LO, A_HI, the lowest and highest
//    values to be histogrammed.  These values will also define the bins.
//
//    Input, int HISTO_NUM, the number of bins to use.
//
//    Output, int HISTO_GRAM[HISTO_NUM+2], contains the number of
//    entries of A in each bin.
//
{
  float delta;
  int *histo_gram;
  int i;
  int j;

  histo_gram = new int[histo_num+2];

  i4vec_zero ( histo_num+2, histo_gram );

  delta = ( a_hi - a_lo ) / ( float ) ( 2 * histo_num );

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < a_lo )
    {
      histo_gram[0] = histo_gram[0] + 1;
    }
    else if ( a[i] <= a_hi )
    {
      j = r4_nint (
        ( ( a_hi -       delta - a[i]        ) * ( float ) ( 1         )
        + (      -       delta + a[i] - a_lo ) * ( float ) ( histo_num ) )
        / ( a_hi - 2.0 * delta        - a_lo ) );

      histo_gram[j] = histo_gram[j] + 1;
    }
    else if ( a_hi < a[i] )
    {
      histo_gram[histo_num+1] = histo_gram[histo_num+1] + 1;
    }
  }

  return histo_gram;
}
//****************************************************************************80

float *r4vec_house_column ( int n, float a[], int k )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The routine returns a vector V that defines a Householder
//    premultiplier matrix H(V) that zeros out the subdiagonal entries of
//    column K of the matrix A.
//
//       H(V) = I - 2 * v * v'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, float A[N], column K of the matrix A.
//
//    Input, int K, the column of the matrix to be modified.
//
//    Output, float R4VEC_HOUSE_COLUMN[N], a vector of unit L2 norm which
//    defines an orthogonal Householder premultiplier matrix H with the property
//    that the K-th column of H*A is zero below the diagonal.
//
{
  int i;
  float s;
  float *v;

  v = r4vec_zero_new ( n );

  if ( k < 1 || n <= k )
  {
    return v;
  }

  s = r4vec_norm_l2 ( n+1-k, a+k-1 );

  if ( s == 0.0 )
  {
    return v;
  }

  v[k-1] = a[k-1] + r4_abs ( s ) * r4_sign ( a[k-1] );

  r4vec_copy ( n-k, a+k, v+k );

  s = r4vec_norm_l2 ( n-k+1, v+k-1 );

  for ( i = k-1; i < n; i++ )
  {
    v[i] = v[i] / s;
  }

  return v;
}
//****************************************************************************80

float r4vec_i4vec_dot_product ( int n, float r4vec[], int i4vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_I4VEC_DOT_PRODUCT computes the dot product of an R4VEC and an I4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, float R4VEC[N], the first vector.
//
//    Input, int I4VEC[N], the second vector.
//
//    Output, float R4VEC_I4VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + r4vec[i] * ( float ) ( i4vec[i] );
  }
  return value;
}
//****************************************************************************80

bool r4vec_in_01 ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_IN_01 is TRUE if the entries of an R4VEC are in the range [0,1].
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], the vector
//
//    Output, bool R4VEC_IN_01, is TRUE if every entry of A is
//    between 0 and 1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 0.0 || 1.0 < a[i] )
    {
      return false;
    }
  }

  return true;
}
//****************************************************************************80

void r4vec_index_delete_all ( int n, float x[], int indx[], float xval,
  int *n2, float x2[], int indx2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_DELETE_ALL deletes all occurrences of a value from an indexed sorted list.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Note that the value of N is adjusted because of the deletions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the current list.
//
//    Input, float X[N], the list.
//
//    Input, int INDX[N], the sort index of the list.
//
//    Input, float XVAL, the value to be sought.
//
//    Output, int *N2, the size of the current list.
//
//    Output, float X2[N2], the list.
//
//    Output, int INDX2[N2], the sort index of the list.
//
{
  int equal;
  int equal1;
  int equal2;
  int get;
  int i;
  int less;
  int more;
  int put;

  if ( n < 1 )
  {
    *n2 = 0;
    return;
  }

  i4vec_copy ( n, indx, indx2 );
  r4vec_copy ( n, x, x2 );
  *n2 = n;

  r4vec_index_search ( *n2, x2, indx2, xval, &less, &equal, &more );

  if ( equal == 0 )
  {
    return;
  }

  equal1 = equal;

  for ( ; ; )
  {
    if ( equal1 <= 1 )
    {
      break;
    }

    if ( x2[indx2[equal1-2]-1] != xval )
    {
      break;
    }
    equal1 = equal1 - 1;
  }

  equal2 = equal;

  for ( ; ; )
  {
    if ( *n2 <= equal2 )
    {
      break;
    }

    if ( x2[indx2[equal2]-1] != xval )
    {
      break;
    }
    equal2 = equal2 + 1;
  }
//
//  Discard certain X values.
//
  put = 0;

  for ( get = 1; get <= *n2; get++ )
  {
    if ( x2[get-1] != xval )
    {
      put = put + 1;
      x2[put-1] = x2[get-1];
    }
  }
//
//  Adjust the INDX values.
//
  for ( equal = equal1; equal <= equal2; equal++ )
  {
    for ( i = 1; i <= *n2; i++ )
    {
      if ( indx2[equal-1] < indx2[i-1] )
      {
        indx2[i-1] = indx2[i-1] - 1;
      }
    }
  }
//
//  Discard certain INDX values.
//
  for ( i = 0; i <= *n2 - equal2 - 1; i++ )
  {
    indx2[equal1+i-1] = indx2[equal2+i];
  }
  for ( i = *n2 + equal1 - equal2; i <= *n2; i++ )
  {
    indx2[i-1] = 0;
  }
//
//  Adjust N.
//
  *n2 = put;

  return;
}
//****************************************************************************80

void r4vec_index_delete_dupes ( int n, float x[], int indx[],
  int *n2, float x2[], int indx2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted list.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The output quantities N2, X2, and INDX2 are computed from the
//    input quantities by sorting, and eliminating duplicates.
//
//    The output arrays should be dimensioned of size N, unless the user
//    knows in advance what the value of N2 will be.
//
//    The output arrays may be identified with the input arrays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the input list.
//
//    Input, float X[N], the list.
//
//    Input, int INDX[N], the sort index of the list.
//
//    Output, int *N2, the number of unique entries in X.
//
//    Output, float X2[N2], a copy of the list which has
//    been sorted, and made unique.
//
//    Output, int INDX2[N2], the sort index of the new list.
//
{
  int i;
  int n3;
  float *x3;

  i = 0;
  n3 = 0;
  x3 = new float[n];

  for ( ; ; )
  {
    i = i + 1;

    if ( n < i )
    {
      break;
    }

    if ( 1 < i )
    {
      if ( x[indx[i-1]-1] == x3[n3-1] )
      {
        continue;
      }
    }
    n3 = n3 + 1;
    x3[n3-1] = x[indx[i-1]-1];
  }
//
//  Set the output data.
//
  *n2 = n3;
  r4vec_copy ( n3, x3, x2 );
  for ( i = 0; i < n3; i++ )
  {
    indx2[i] = i + 1;
  }

  delete [] x3;

  return;
}
//****************************************************************************80

void r4vec_index_delete_one ( int n, float x[], int indx[], float xval,
  int *n2, float x2[], int indx2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_DELETE_ONE deletes one copy of a value from an indexed sorted list.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    If the value occurs in the list more than once, only one copy is deleted.
//
//    Note that the value of N is adjusted because of the deletions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the current list.
//
//    Input, float X[N], the list.
//
//    Input, int INDX[N], the sort index of the list.
//
//    Input, float XVAL, the value to be sought.
//
//    Output, int *N2, the size of the current list.
//
//    Output, float X2[N2], the list.
//
//    Output, int INDX2[N2], the sort index of the list.
//
{
  int equal;
  int i;
  int j;
  int less;
  int more;

  if ( n < 1 )
  {
    *n2 = 0;
    return;
  }

  *n2 = n;
  i4vec_copy ( *n2, indx, indx2 );
  r4vec_copy ( *n2, x, x2 );

  r4vec_index_search ( *n2, x2, indx2, xval, &less, &equal, &more );

  if ( equal != 0 )
  {
    j = indx2[equal-1];
    for ( i = j; i <= *n2-1; i++ )
    {
      x2[i-1] = x[i];
    }
    for ( i = equal; i <= *n2-1; i++ )
    {
      indx2[i-1] = indx2[i];
    }
    for ( i = 1; i <= *n2 - 1; i++ )
    {
      if ( j < indx2[i-1] )
      {
        indx2[i-1] = indx2[i-1] - 1;
      }
    }
    *n2 = *n2 - 1;
  }

  return;
}
//****************************************************************************80

void r4vec_index_insert ( int *n, float x[], int indx[], float xval )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_INSERT inserts a value in an indexed sorted list.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *N, the size of the current list.
//
//    Input, float X[N], the list.
//
//    Input, int INDX[N], the sort index of the list.
//
//    Input, float XVAL, the value to be sought.
//
{
  int equal;
  int i;
  int less;
  int more;

  if ( *n <= 0 )
  {
    *n = 1;
    x[0] = xval;
    indx[0] = 1;
    return;
  }

  r4vec_index_search ( *n, x, indx, xval, &less, &equal, &more );

  x[*n] = xval;
  for ( i = *n; more <= i; i-- )
  {
    indx[i] = indx[i-1];
  }
  indx[more-1] = *n + 1;
  *n = *n + 1;

  return;
}
//****************************************************************************80

void r4vec_index_insert_unique ( int *n, float x[], int indx[], float xval )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_INSERT_UNIQUE inserts a unique value in an indexed sorted list.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *N, the size of the current list.
//    If the input value XVAL does not already occur in X, then N is increased.
//
//    Input/output, float X[N], the list.
//    If the input value XVAL does not already occur in X, then it is added
//    to X.
//
//    Input/output, int INDX[N], the sort index of the list.
//    If the input value XVAL does not already occur in X, then INDX is updated.
//
//    Input, float XVAL, the value which will be inserted into the X
//    vector if it is not there already.
//
{
  int equal;
  int i;
  int less;
  int more;

  if ( *n <= 0 )
  {
    *n = 1;
    x[0] = xval;
    indx[0] = 1;
    return;
  }
//
//  Does XVAL already occur in X?
//
  r4vec_index_search ( *n, x, indx, xval, &less, &equal, &more );

  if ( equal == 0 )
  {
    x[*n] = xval;
    for ( i = *n; more <= i; i-- )
    {
      indx[i] = indx[i-1];
    }
    indx[more-1] = *n + 1;
    *n = *n + 1;
  }

  return;
}
//****************************************************************************80

void r4vec_index_order ( int n, float x[], int indx[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_ORDER sorts an R4VEC using an index vector.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The index vector itself is not modified.  Therefore, the pair
//    (X,INDX) no longer represents an index sorted vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the current list.
//
//    Input/output, float X[N], the list.  On output, the list
//    has been sorted.
//
//    Input, int INDX[N], the sort index of the list.
//
{
  int i;
  float *y;

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[indx[i]-1];
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = y[i];
  }
  delete [] y;

  return;
}
//****************************************************************************80

void r4vec_index_search ( int n, float x[], int indx[], float xval, int *less,
  int *equal, int *more )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_SEARCH searches for a value in an indexed sorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the current list.
//
//    Input, float X[N], the list.
//
//    Input, int INDX[N], the sort index of the list.
//
//    Input, float XVAL, the value to be sought.
//
//    Output, int *LESS, *EQUAL, *MORE, the indexes in INDX of the
//    entries of X that are just less than, equal to, and just greater
//    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
//    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
//    is the greatest entry of X, then MORE is N+1.
//
{
  int hi;
  int lo;
  int mid;
  float xhi;
  float xlo;
  float xmid;

  if ( n <= 0 )
  {
    *less = 0;
    *equal = 0;
    *more = 0;
    return;
  }

  lo = 1;
  hi = n;
  xlo = x[indx[lo-1]-1];
  xhi = x[indx[hi-1]-1];

  if ( xval < xlo )
  {
    *less = 0;
    *equal = 0;
    *more = 1;
    return;
  }
  else if ( xval == xlo )
  {
    *less = 0;
    *equal = 1;
    *more = 2;
    return;
  }

  if ( xhi < xval )
  {
    *less = n;
    *equal = 0;
    *more = n + 1;
    return;
  }
  else if ( xval == xhi )
  {
    *less = n - 1;
    *equal = n;
    *more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      *less = lo;
      *equal = 0;
      *more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];

    if ( xval == xmid )
    {
      *equal = mid;
      *less = mid - 1;
      *more = mid + 1;
      return;
    }
    else if ( xval < xmid )
    {
      hi = mid;
    }
    else if ( xmid < xval )
    {
      lo = mid;
    }
  }
  return;
}
//****************************************************************************80

void r4vec_index_sort_unique ( int n, float x[], int *n2, float x2[],
  int indx2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_SORT_UNIQUE creates a sort index for an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the current list.
//
//    Input, float X[N], the list.
//
//    Output, int *N2, the number of unique elements in X.
//
//    Output, float X2[N2], a list of the unique elements of X.
//
//    Output, int INDX2[N2], the sort index of the list.
//
{
  int i;

  *n2 = 0;

  for ( i = 0; i < n; i++ )
  {
    r4vec_index_insert_unique ( n2, x2, indx2, x[i] );
  }

  for ( i = *n2; i < n; i++ )
  {
    x2[i] = -1;
  }
  for ( i = *n2; i < n; i++ )
  {
    indx2[i] = -1;
  }

  return;
}
//****************************************************************************80

void r4vec_index_sorted_range ( int n, float r[], int indx[], float r_lo,
  float r_hi, int *i_lo, int *i_hi )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEX_SORTED_RANGE: search index sorted vector for elements in a range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items in the vector.
//
//    Input, float R[N], the index sorted vector.
//
//    Input, int INDX[N], the vector used to sort R.
//    The vector R[INDX[*]] is sorted.
//
//    Input, float R_LO, R_HI, the limits of the range.
//
//    Output, int *I_LO, *I_HI, the range of indices
//    so that I_LO <= I <= I_HI => R_LO <= R[INDX[I]] <= R_HI.  If no
//    values in R lie in the range, then I_HI < I_LO will be returned.
//
{
  int i1;
  int i2;
  int j1;
  int j2;
//
//  Cases we can handle immediately.
//
  if ( r[indx[n-1]] < r_lo )
  {
    *i_lo = n;
    *i_hi = n - 1;
    return;
  }

  if ( r_hi < r[indx[0]] )
  {
    *i_lo = 0;
    *i_hi = -1;
    return;
  }
//
//  Are there are least two intervals?
//
  if ( n == 1 )
  {
    if ( r_lo <= r[indx[0]] && r[indx[0]] <= r_hi )
    {
      *i_lo = 0;
      *i_hi = 0;
    }
    else
    {
      *i_lo = -1;
      *i_hi = -2;
    }
    return;
  }
//
//  Bracket R_LO.
//
  if ( r_lo <= r[indx[0]] )
  {
    *i_lo = 0;
  }
  else
  {
//
//  R_LO is in one of the intervals spanned by R(INDX(J1)) to R(INDX(J2)).
//  Examine the intermediate interval [R(INDX(I1)), R(INDX(I1+1))].
//  Does R_LO lie here, or below or above?
//
    j1 = 0;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_lo < r[indx[i1]] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[indx[i2]] < r_lo )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_lo = i1;
        break;
      }
    }
  }
//
//  Bracket R_HI.
//
  if ( r[indx[n-1]] <= r_hi )
  {
    *i_hi = n - 1;
  }
  else
  {
    j1 = *i_lo;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_hi < r[indx[i1]] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[indx[i2]] < r_hi )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_hi = i2;
        break;
      }
    }
  }
//
//  We expect to have computed the largest I_LO and smallest I_HI such that
//    R(INDX(I_LO)) <= R_LO <= R_HI <= R(INDX(I_HI))
//  but what we want is actually
//    R_LO <= R(INDX(I_LO)) <= R(INDX(I_HI)) <= R_HI
//  which we can usually get simply by incrementing I_LO and decrementing I_HI.
//
  if ( r[indx[*i_lo]] < r_lo )
  {
    *i_lo = *i_lo + 1;
    if ( n - 1 < *i_lo )
    {
      *i_hi = *i_lo - 1;
    }
  }

  if ( r_hi < r[indx[*i_hi]] )
  {
    *i_hi = *i_hi - 1;
    if ( i_hi < 0 )
    {
      *i_lo = *i_hi + 1;
    }
  }

  return;
}
//****************************************************************************80

void r4vec_indexed_heap_d ( int n, float a[], int indx[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEXED_HEAP_D creates a descending heap from an indexed R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    An indexed R4VEC is an R4VEC of data values, and an R4VEC of N indices,
//    each referencing an entry of the data vector.
//
//    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
//    we have:
//      A[INDX[2*J+1]]   <= A[INDX[J]]
//    and
//      A[INDX[2*J+2]] <= A[INDX[J]]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the index array.
//
//    Input, float A[*], the data vector.
//
//    Input/output, int INDX[N], the index array.
//    Each entry of INDX must be a valid index for the array A.
//    On output, the indices have been reordered into a descending heap.
//
{
  int i;
  int ifree;
  int key;
  int m;
//
//  Only nodes N/2 - 1 down to 0 can be "parent" nodes.
//
  for ( i = ( n / 2 ) - 1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = indx[i];
    ifree = i;

    for ( ; ; )
    {
//
//  Positions 2*IFREE+1 and 2*IFREE+2 are the descendants of position
//  IFREE.  (One or both may not exist because they exceed N-1.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n - 1 < m )
      {
        break;
      }
//
//  Does the second position exist?
//
      if ( m + 1 <= n - 1 )
      {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
        if ( a[indx[m]] < a[indx[m+1]] )
        {
          m = m + 1;
        }
      }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
      if ( a[indx[m]] <= a[key] )
      {
        break;
      }

      indx[ifree] = indx[m];
      ifree = m;
    }
//
//  Once there is no more shifting to do, KEY moves into the free spot IFREE.
//
    indx[ifree] = key;
  }

  return;
}
//****************************************************************************80

int r4vec_indexed_heap_d_extract ( int *n, float a[], int indx[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    An indexed R4VEC is an R4VEC of data values, and an R4VEC of N indices,
//    each referencing an entry of the data vector.
//
//    The routine finds the maximum value in the heap, returns that value to the
//    user, deletes that value from the heap, and restores the heap to its
//    proper form.
//
//    Note that the argument N must be a variable, which will be decremented
//    before return, and that INDX will hold one less value on output than it
//    held on input.
//
//    This is one of three functions needed to model a priority queue.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Thomas Cormen, Charles Leiserson, Ronald Rivest,
//    Introduction to Algorithms,
//    MIT Press, 2001,
//    ISBN: 0262032937,
//    LC: QA76.C662.
//
//  Parameters:
//
//    Input/output, int *N, the number of items in the index vector.
//
//    Input, float A[*], the data vector.
//
//    Input/output, int INDX[N], the index vector.
//
//    Output, int R4VEC_INDEXED_HEAP_D_EXTRACT, the index in A of the item of
//    maximum value, which has now been removed from the heap.
//
{
  int indx_extract;

  if ( *n < 1 )
  {
    cerr << "\n";
    cerr << "R4VEC_INDEXED_HEAP_D_EXTRACT - Fatal error!\n";
    cerr << "  The heap is empty.\n";
    exit ( 1 );
  }
//
//  Get the index of the maximum value.
//
  indx_extract = indx[0];

  if ( *n == 1 )
  {
    *n = 0;
    return indx_extract;
  }
//
//  Shift the last index down.
//
  indx[0] = indx[*n-1];
//
//  Restore the heap structure.
//
  *n = *n - 1;
  r4vec_indexed_heap_d ( *n, a, indx );

  return indx_extract;
}
//****************************************************************************80

void r4vec_indexed_heap_d_insert ( int *n, float a[], int indx[],
  int indx_insert )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    An indexed R4VEC is an R4VEC of data values, and an R4VEC of N indices,
//    each referencing an entry of the data vector.
//
//    Note that the argument N must be a variable, and will be incremented before
//    return, and that INDX must be able to hold one more entry on output than
//    it held on input.
//
//    This is one of three functions needed to model a priority queue.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Thomas Cormen, Charles Leiserson, Ronald Rivest,
//    Introduction to Algorithms,
//    MIT Press, 2001,
//    ISBN: 0262032937,
//    LC: QA76.C662.
//
//  Parameters:
//
//    Input/output, int *N, the number of items in the index vector.
//
//    Input, float A[*], the data vector.
//
//    Input/output, int INDX[N], the index vector.
//
//    Input, int INDX_INSERT, the index in A of the value
//    to be inserted into the heap.
//
{
  int i;
  int parent;

  *n = *n + 1;
  i = *n - 1;

  while ( 0 < i )
  {
    parent = ( i - 1 ) / 2;

    if ( a[indx_insert] <= a[indx[parent]] )
    {
      break;
    }

    indx[i] = indx[parent];
    i = parent;
  }

  indx[i] = indx_insert;

  return;
}
//****************************************************************************80

int r4vec_indexed_heap_d_max ( int n, float a[], int indx[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    An indexed R4VEC is an R4VEC of data values, and an R4VEC of N indices,
//    each referencing an entry of the data vector.
//
//    This is one of three functions needed to model a priority queue.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Thomas Cormen, Charles Leiserson, Ronald Rivest,
//    Introduction to Algorithms,
//    MIT Press, 2001,
//    ISBN: 0262032937,
//    LC: QA76.C662.
//
//  Parameters:
//
//    Input, int N, the number of items in the index vector.
//
//    Input, float A[*], the data vector.
//
//    Input, int INDX[N], the index vector.
//
//    Output, int R4VEC_INDEXED_HEAP_D_MAX, the index in A of the maximum value
//    in the heap.
//
{
  int indx_max;

  indx_max = indx[0];

  return indx_max;
}
//****************************************************************************80

void r4vec_indicator0 ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDICATOR0 sets an R4VEC to the indicator vector {0,1,2,...}.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, float A[N], the array to be initialized.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( float ) i;
  }

  return;
}
//****************************************************************************80

float *r4vec_indicator0_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDICATOR0_NEW sets an R4VEC to the indicator vector {0,1,2,...}.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, float R4VEC_INDICATOR0_NEW[N], the indicator array.
//
{
  float *a;
  int i;

  a = new float[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( float ) i;
  }

  return a;
}
//****************************************************************************80

void r4vec_indicator1 ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDICATOR1 sets an R4VEC to the indicator vector {1,2,3,...}.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, float A[N], the array to be initialized.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( float ) ( i + 1 );
  }

  return;
}
//****************************************************************************80

float *r4vec_indicator1_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INDICATOR1_NEW sets an R4VEC to the indicator vector {1,2,3,...}.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, float R4VEC_INDICATOR1_NEW[N], the indicator array.
//
{
  float *a;
  int i;

  a = new float[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( float ) ( i + 1 );
  }

  return a;
}
//****************************************************************************80

void r4vec_insert ( int n, float a[], int pos, float value )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_INSERT inserts a value into an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the array on input.
//
//    Input/output, float A[N+1], the array.  On input, A is
//    assumed to contain only N entries, while on output, A actually
//    contains N+1 entries.
//
//    Input, int POS, the position to be assigned the new entry.
//    1 <= POS <= N+1.
//
//    Input, float VALUE, the value to be inserted.
//
{
  int i;

  if ( pos < 1 || n + 1 < pos )
  {
    cerr << "\n";
    cerr << "R4VEC_INSERT - Fatal error!\n";
    cerr << "  Illegal insertion position = " << pos << "\n";;
    exit ( 1 );
  }
  else
  {
    for ( i = n+1; pos+1 <= i; i-- )
    {
      a[i-1] = a[i-2];
    }

    a[pos-1] = value;
  }

  return;
}
//****************************************************************************80

bool r4vec_is_int ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_IS_INT is TRUE if an R4VEC is integral.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], the vector
//
//    Output, bool R4VEC_IS_INT, is TRUE if every entry of A is an integer.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] != ( float ) ( int ) a[i] )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

bool r4vec_is_nonnegative ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_IS_NONNEGATIVE is true if all entries in an R4VEC are nonnegative.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[N], the vector to be checked.
//
//    Output, bool R4VEC_IS_NONNEGATIVE is true if all elements of X
//    are nonnegative.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < 0.0 )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

bool r4vec_is_zero ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_IS_ZERO is true if the entries in an R4VEC are all zero.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[N], the vector to be checked.
//
//    Output, bool R4VEC_IS_ZERO is true if all N elements of X
//    are zero.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] != 0.0 )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

float *r4vec_linspace_new ( int n, float a_first, float a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float A_FIRST, A_LAST, the first and last entries.
//
//    Output, float R4VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  float *a;
  int i;

  a = new float[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( float ) ( n - 1 - i ) * a_first 
             + ( float ) (         i ) * a_last ) 
             / ( float ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

bool r4vec_lt ( int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_LT == ( A1 < A2 ) for two R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The comparison is lexicographic.
//
//    A1 < A2  <=>                              A1(1) < A2(1) or
//                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
//                 ...
//                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, float A1[N], A2[N], the vectors to be compared.
//
//    Output, bool R4VEC_LT, is TRUE if and only if A1 < A2.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] < a2[i] )
    {
      return true;
    }
    else if ( a2[i] < a1[i] )
    {
      return false;
    }

  }

  return false;
}
//****************************************************************************80

void r4vec_mask_print ( int n, float a[], int mask_num, int mask[],
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MASK_PRINT prints a masked R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, float A[N], the vector to be printed.
//
//    Input, int MASK_NUM, the number of masked elements.
//
//    Input, int MASK[MASK_NUM], the indices of the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << "  Masked vector printout:\n";

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < mask_num; i++ )
  {
    cout << "  " << setw(6)  << i
         << ": " << setw(6)  << mask[i]
         << "  " << setw(12) << a[mask[i]-1] << "\n";
  }

  return;
}
//****************************************************************************80

float r4vec_max ( int n, float r4vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MAX returns the value of the maximum element in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float R4VEC[N], a pointer to the first entry of the array.
//
//    Output, float R4VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  float value;

  value = - r4_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( value < r4vec[i] )
    {
      value = r4vec[i];
    }
  }
  return value;
}
//****************************************************************************80

int r4vec_max_index ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MAX_INDEX returns the index of the maximum value in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], the array.
//
//    Output, int R4VEC_MAX_INDEX, the index of the largest entry.
//
{
  int i;
  int max_index;

  if ( n <= 0 )
  {
    max_index = -1;
  }
  else
  {
    max_index = 0;

    for ( i = 1; i < n; i++ )
    {
      if ( a[max_index] < a[i] )
      {
        max_index = i;
      }
    }
  }

  return max_index;
}
//****************************************************************************80

float r4vec_mean ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MEAN returns the mean of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[N], the vector whose mean is desired.
//
//    Output, float R4VEC_MEAN, the mean, or average, of the vector entries.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + x[i];
  }

  value = value / ( float ) n;

  return value;
}
//****************************************************************************80

float r4vec_median ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MEDIAN returns the median of an unsorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Hoare's algorithm is used.  The values of the vector are
//    rearranged by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input/output, float A[N], the array to search.  On output,
//    the order of the elements of A has been somewhat changed.
//
//    Output, float R4VEC_MEDIAN, the value of the median of A.
//
{
  int k;
  float median;

  k = ( n + 1 ) / 2;

  median = r4vec_frac ( n, a, k );

  return median;
}
//****************************************************************************80

float r4vec_min ( int n, float r4vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MIN returns the value of the minimum element in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float R4VEC[N], the array to be checked.
//
//    Output, float R4VEC_MIN, the value of the minimum element.
//
{
  int i;
  float value;

  value = r4_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( r4vec[i] < value )
    {
      value = r4vec[i];
    }
  }
  return value;
}
//****************************************************************************80

int r4vec_min_index ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MIN_INDEX returns the index of the minimum value in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], the array.
//
//    Output, int R4VEC_MIN_INDEX, the index of the smallest entry.
//
{
  int i;
  int min_index;

  if ( n <= 0 )
  {
    min_index = -1;
  }
  else
  {
    min_index = 0;

    for ( i = 1; i < n; i++ )
    {
      if ( a[i] < a[min_index] )
      {
        min_index = i;
      }
    }
  }

  return min_index;
}
//****************************************************************************80

float r4vec_min_pos ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MIN_POS returns the minimum positive value of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, float A[N], the array.
//
//    Output, float R4VEC_MIN_POS, the smallest positive entry.
//
{
  int i;
  float value;

  value = r4_huge ( );

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      if ( a[i] < value )
      {
        value = a[i];
      }
    }
  }
  return value;
}
//****************************************************************************80

bool r4vec_mirror_next ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MIRROR_NEXT steps through all sign variations of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    In normal use, the user would set every element of A to be positive.
//    The routine will take the input value of A, and output a copy in
//    which the signs of one or more entries have been changed.  Repeatedly
//    calling the routine with the output from the previous call will generate
//    every distinct "variation" of A; that is, all possible sign variations.
//
//    When the output variable DONE is TRUE (or equal to 1), then the
//    output value of A_NEW is the last in the series.
//
//    Note that A may have some zero values.  The routine will essentially
//    ignore such entries; more exactly, it will not stupidly assume that -0
//    is a proper "variation" of 0.
//
//    Also, it is possible to call this routine with the signs of A set
//    in any way you like.  The routine will operate properly, but it
//    will nonethess terminate when it reaches the value of A in which
//    every nonzero entry has negative sign.
//
//
//    More efficient algorithms using the Gray code seem to require internal
//    memory in the routine, which is not one of MATLAB's strong points,
//    or the passing back and forth of a "memory array", or the use of
//    global variables, or unnatural demands on the user.  This form of
//    the routine is about as clean as I can make it.
//
//  Example:
//
//      Input         Output
//    ---------    --------------
//    A            A         DONE
//    ---------    --------  ----
//     1  2  3     -1  2  3  false
//    -1  2  3      1 -2  3  false
//     1 -2  3     -1 -2  3  false
//    -1 -2  3      1  2 -3  false
//     1  2 -3     -1  2 -3  false
//    -1  2 -3      1 -2 -3  false
//     1 -2 -3     -1 -2 -3  false
//    -1 -2 -3      1  2  3  true
//
//     1  0  3     -1  0  3  false
//    -1  0  3      1  0 -3  false
//     1  0 -3     -1  0 -3  false
//    -1  0 -3      1  0  3  true
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2008
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
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, float A[N], a vector of real numbers.  On
//    output, some signs have been changed.
//
//    Output, bool R4VEC_MIRROR_NEXT, is TRUE if the input vector A was
//    the last element
//    in the series (every entry was nonpositive); the output vector is reset
//    so that all entries are nonnegative, but presumably the ride is over.
//
{
  bool done;
  int i;
  int positive;
//
//  Seek the first strictly positive entry of A.
//
  positive = -1;
  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      positive = i;
      break;
    }
  }
//
//  If there is no strictly positive entry of A, there is no successor.
//
  if ( positive == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = - a[i];
    }
    done = true;
    return done;
  }
//
//  Otherwise, negate A up to the positive entry.
//
  for ( i = 0; i <= positive; i++ )
  {
    a[i] = - a[i];
  }
  done = false;

  return done;
}
//****************************************************************************80

bool r4vec_negative_strict ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NEGATIVE_STRICT: all entries of R4VEC are strictly negative.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.
//
//    Input, float A[N], the vector.
//
//    Output, bool R4VEC_NEGATIVE_STRICT, is TRUE if every entry of
//    A is strictly negative.
//
{
  int i;
  bool value;

  for ( i = 0; i < n; i++ )
  {
    if ( 0 <= a[i] )
    {
      value = false;
      return value;
    }
  }
  value = true;
  return value;
}
//****************************************************************************80

float *r4vec_nint ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NINT rounds the entries of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], the vector to be rounded.
//
//    Output, float B[N], the rounded values.
//
{
  float *b;
  int i;
  int s;

  b = new float[n];

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 0.0 )
    {
      s = -1;
    }
    else
    {
      s = 1;
    }
    b[i] = ( float ) ( s * ( int ) ( r4_abs ( a[i] ) + 0.5 ) );
  }

  return b;
}
//****************************************************************************80

float r4vec_norm ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORM returns the L2 norm of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L2 norm is defined as:
//
//      R4VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], the vector whose L2 norm is desired.
//
//    Output, float R4VEC_NORM, the L2 norm of A.
//
{
  int i;
  float v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
//****************************************************************************80

float r4vec_norm_affine ( int n, float v0[], float v1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORM_AFFINE returns the affine L2 norm of an R4VEC.
//
//  Discussion:
//
//    The affine vector L2 norm is defined as:
//
//      R4VEC_NORM_AFFINE(V0,V1)
//        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, float V0[N], the base vector.
//
//    Input, float V1[N], the vector.
//
//    Output, float R4VEC_NORM_AFFINE, the affine L2 norm.
//
{
  int i;
  float value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

float r4vec_norm_l1 ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORM_L1 returns the L1 norm of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L1 norm is defined as:
//
//      R4VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], the vector whose L1 norm is desired.
//
//    Output, float R4VEC_NORM_L1, the L1 norm of A.
//
{
  int i;
  float v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + r4_abs ( a[i] );
  }

  return v;
}
//****************************************************************************80

float r4vec_norm_l2 ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORM_L2 returns the L2 norm of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L2 norm is defined as:
//
//      R4VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], the vector whose L2 norm is desired.
//
//    Output, float R4VEC_NORM_L2, the L2 norm of A.
//
{
  int i;
  float v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
//****************************************************************************80

float r4vec_norm_li ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORM_LI returns the L-oo norm of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector L-oo norm is defined as:
//
//      R4VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], the vector whose L-oo norm is desired.
//
//    Output, float R4VEC_NORM_LI, the L-oo norm of A.
//
{
  int i;
  float v1;
  float v2;

  v1 = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v2 = r4_abs ( a[i] );

    if ( v1 < v2 )
    {
      v1 = v2;
    }
  }

  return v1;
}
//****************************************************************************80

float r4vec_norm_lp ( int n, float a[], float p )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORM_LP returns the LP norm of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The vector LP norm is defined as:
//
//      R4VEC_NORM_LP = ( sum ( 1 <= I <= N ) ( abs ( A(I) ) )^P )^(1/P).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], the vector whose LP norm is desired.
//
//    Input, float P, the index of the norm.
//
//    Output, float R4VEC_NORML_LP, the LP norm of A.
//
{
  int i;
  float v;

  v = 0.0;

  if ( p == 1.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + r4_abs ( a[i] );
    }
  }
  else if ( p == 2.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + a[i] * a[i];
    }
    v = sqrt ( v );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + pow ( r4_abs ( a[i] ), p );
    }
    v = pow (  ( double ) v, 1.0 / p );
  }

  return v;
}
//****************************************************************************80

float *r4vec_normal_01 ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORMAL_01 returns a unit pseudonormal R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, float X[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, float R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, float Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
  int i;
  int m;
  static int made = 0;
  float *r;
  const float r4_pi = 3.141592653589793;
  static int saved = 0;
  float *x;
  int x_hi;
  int x_lo;
  static float y = 0.0;
//
//  I'd like to allow the user to reset the internal data.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }

  x = new float[n];
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
//
//  Maybe we don't need any more values.
//
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r4vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * r4_pi * r[1] );
    y =         sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * r4_pi * r[1] );

    saved = 1;

    made = made + 2;

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r4vec_uniform_01_new ( 2 * m, seed );

    for ( i = 0; i <= 2 * m - 2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r4_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r4_pi * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r4vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r4_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r4_pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r4_pi * r[i+1] );
    y           = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r4_pi * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
}
//****************************************************************************80

void r4vec_normalize ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORMALIZE normalizes an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, float A[N], the vector to be normalized.
//    On output, A should have unit Euclidean norm.
//
{
  int i;
  float norm;

  norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    norm = norm + a[i] * a[i];
  }
  norm = sqrt ( norm );

  if ( norm == 0.0 )
  {
    cerr << "\n";
    cerr << "R4VEC_NORMALIZE - Fatal error!\n";
    cerr << "  The vector norm is 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / norm;
  }

  return;
}
//****************************************************************************80

void r4vec_normalize_l1 ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORMALIZE_L1 normalizes an R4VEC to have unit sum.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, float A[N], the vector to be normalized.
//    On output, the entries of A should have unit sum.  However, if
//    the input vector has zero sum, the routine halts.
//
{
  float a_sum;
  int i;

  a_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    a_sum = a_sum + a[i];
  }

  if ( a_sum == 0.0 )
  {
    cerr << "\n";
    cerr << "R4VEC_NORMALIZE_L1 - Fatal error!\n";
    cerr << "  The vector entries sum to 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / a_sum;
  }

  return;
}
//****************************************************************************80

float r4vec_normsq ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORMSQ returns the squared L2 norm of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The squared vector L2 norm is defined as:
//
//      R4VEC_NORMSQ = sum ( 1 <= I <= N ) A(I)^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the vector dimension.
//
//    Input, float A[N], the vector.
//
//    Output, float R4VEC_NORMSQ, the squared L2 norm.
//
{
  int i;
  float v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  return v;
}
//****************************************************************************80

float r4vec_normsq_affine ( int n, float v0[], float v1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_NORMSQ_AFFINE returns the squared affine L2 norm of an R4VEC.
//
//  Discussion:
//
//    The squared affine vector L2 norm is defined as:
//
//      R4VEC_NORMSQ_AFFINE(V0,V1)
//        = sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, float V0[N], the base vector.
//
//    Input, float V1[N], the vector whose squared affine L2 norm is desired.
//
//    Output, float R4VEC_NORMSQ_AFFINE, the squared affine L2 norm.
//
{
  int i;
  float value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  return value;
}
//****************************************************************************80

float *r4vec_ones_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_ONES_NEW creates a vector of 1's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, float R4VEC_ONES_NEW[N], a vector of 1's.
//
{
  float *a;
  int i;

  a = new float[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 1.0;
  }
  return a;
}
//****************************************************************************80

int r4vec_order_type ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_ORDER_TYPE determines if an R4VEC is (non)strictly ascending/descending.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the array.
//
//    Input, float X[N], the array to be checked.
//
//    Output, int R4VEC_ORDER_TYPE, order indicator:
//    -1, no discernable order;
//    0, all entries are equal;
//    1, ascending order;
//    2, strictly ascending order;
//    3, descending order;
//    4, strictly descending order.
//
{
  int i;
  int order;
//
//  Search for the first value not equal to X(0).
//
  i = 0;

  for ( ; ; )
  {
    i = i + 1;
    if ( n-1 < i )
    {
      order = 0;
      return order;
    }

    if ( x[0] < x[i] )
    {
      if ( i == 1 )
      {
        order = 2;
        break;
      }
      else
      {
        order = 1;
        break;
      }
    }
    else if ( x[i] < x[0] )
    {
      if ( i == 1 )
      {
        order = 4;
        break;
      }
      else
      {
        order = 3;
        break;
      }
    }
  }
//
//  Now we have a "direction".  Examine subsequent entries.
//
  for ( ; ; )
  {
    i = i + 1;
    if ( n - 1 < i )
    {
      break;
    }

    if ( order == 1 )
    {
      if ( x[i] < x[i-1] )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 2 )
    {
      if ( x[i] < x[i-1] )
      {
        order = -1;
        break;
      }
      else if ( x[i] == x[i-1] )
      {
        order = 1;
      }
    }
    else if ( order == 3 )
    {
      if ( x[i-1] < x[i] )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 4 )
    {
      if ( x[i-1] < x[i] )
      {
        order = -1;
        break;
      }
      else if ( x[i] == x[i-1] )
      {
        order = 3;
      }
    }
  }
  return order;
}
//****************************************************************************80

void r4vec_part_quick_a ( int n, float a[], int *l, int *r )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PART_QUICK_A reorders an R4VEC as part of a quick sort.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The routine reorders the entries of A.  Using A[0] as a
//    key, all entries of A that are less than or equal to A[0] will
//    precede A[0] which precedes all entries that are greater than A[0].
//
//  Example:
//
//    Input:
//
//  N = 8
//
//  A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
//
//    Output:
//
//  L = 3, R = 6
//
//  A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
//        -------        -------
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, float A[N].  On input, the array to be checked.
//    On output, A has been reordered as described above.
//
//    Output, int L, R, the indices of A that define the three segments.
//    Let KEY = the input value of A[0].  Then
//    I <= L             A(I) < KEY;
//     L < I < R         A(I) = KEY;
//             R <= I    A(I) > KEY.
//
{
  int i;
  float key;
  int m;
  float temp;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R4VEC_PART_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }
  else if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key = a[0];
  m = 1;
//
//  The elements of unknown size have indices between L+1 and R-1.
//
  *l = 1;
  *r = n + 1;

  for ( i = 2; i <= n; i++ )
  {

    if ( key < a[*l] )
    {
      *r = *r - 1;
      temp = a[*r-1];
      a[*r-1] = a[*l];
      a[*l] = temp;
    }
    else if ( a[*l] == key )
    {
      m = m + 1;
      temp = a[m-1];
      a[m-1] = a[*l];
      a[*l] = temp;
      *l = *l + 1;
    }
    else if ( a[*l] < key )
    {
      *l = *l + 1;
    }

  }
//
//  Now shift small elements to the left, and KEY elements to center.
//
  for ( i = 1; i <= *l -m; i++ )
  {
    a[i-1] = a[i+m-1];
  }

  *l = *l - m;

  for ( i = *l+1; i <= *l+m; i++ )
  {
    a[i-1] = key;
  }

  return;
}
//****************************************************************************80

void r4vec_permute ( int n, int p[], int base, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PERMUTE permutes an R4VEC in place.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   2,   4,   5,   1,   3 )
//      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
//      BASE = 1
//
//    Output:
//
//      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.
//
//    Input/output, float A[N], the array to be permuted.
//
{
  float a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "R4VEC_PERMUTE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
//  So temporarily add 1-BASE to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "R4VEC_PERMUTE - Fatal error!\n";
          cerr << "  A permutation index is out of range.\n";
          cerr << "  P(" << iput << ") = " << iget << "\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the base of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 +  base;
  }
  return;
}
//****************************************************************************80

void r4vec_permute_cyclic ( int n, int k, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PERMUTE_CYCLIC performs a cyclic permutation of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    For 0 <= K < N, this function cyclically permutes the input vector
//    to have the form
//
//     ( A[K], A[K+1], ..., A[N-1], A[0], ..., A[K-1] )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int K, the increment used.
//
//    Input/output, float A[N], the array to be permuted.
//
{
  float *b;
  int i;
  int ipk;

  b = new float[n];

  for ( i = 0; i < n; i++ )
  {
    ipk = i4_wrap ( i + k, 0, n - 1 );
    b[i] = a[ipk];
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = b[i];
  }

  delete [] b;

  return;
}
//****************************************************************************80

void r4vec_permute_uniform ( int n, float a[], int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PERMUTE_UNIFORM randomly permutes an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input/output, float A[N], the array to be permuted.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  int base = 0;
  int *p;

  p = perm_uniform_new ( n, base, seed );

  r4vec_permute ( n, p, base, a );

  delete [] p;

  return;
}
//****************************************************************************80

void r4vec_polarize ( int n, float a[], float p[], float a_normal[],
  float a_parallel[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_POLARIZE decomposes an R4VEC into normal and parallel components.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The (nonzero) vector P defines a direction.
//
//    The vector A can be written as the sum
//
//      A = A_normal + A_parallel
//
//    where A_parallel is a linear multiple of P, and A_normal
//    is perpendicular to P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], the vector to be polarized.
//
//    Input, float P[N], the polarizing direction.
//
//    Output, float A_NORMAL[N], A_PARALLEL[N], the normal
//    and parallel components of A.
//
{
  float a_dot_p;
  int i;
  float p_norm;

  p_norm = r4vec_norm ( n, p );

  if ( p_norm == 0.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      a_normal[i] = a[i];
    }
    for ( i = 0; i < n; i++ )
    {
      a_parallel[i] = 0.0;
    }
    return;
  }
  a_dot_p = r4vec_dot_product ( n, a, p ) / p_norm;

  for ( i = 0; i < n; i++ )
  {
    a_parallel[i] = a_dot_p * p[i] / p_norm;
  }

  for ( i = 0; i < n; i++ )
  {
    a_normal[i] = a[i] - a_parallel[i];
  }

  return;
}
//****************************************************************************80

bool r4vec_positive_strict ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_POSITIVE_STRICT: all entries of R4VEC are strictly positive.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.
//
//    Input, float A[N], the vector.
//
//    Output, bool R4VEC_POSITIVE_STRICT, is TRUE if every entry of
//    A is strictly positive.
//
{
  int i;
  bool value;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] <= 0.0 )
    {
      value = false;
      return value;
    }
  }
  value = true;
  return value;
}
//****************************************************************************80

void r4vec_print ( int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PRINT prints an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, float A[N], the vector to be printed.
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
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void r4vec_print_part ( int n, float a[], int max_print, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PRINT_PART prints "part" of an R4VEC.
//
//  Discussion:
//
//    The user specifies MAX_PRINT, the maximum number of lines to print.
//
//    If N, the size of the vector, is no more than MAX_PRINT, then
//    the entire vector is printed, one entry per line.
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
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, float A[N], the vector to be printed.
//
//    Input, int MAX_PRINT, the maximum number of lines
//    to print.
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
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i] << "\n";
    }
    cout << "  ........  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i] << "\n";
  }
  else
  {
    for ( i= 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i] << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i]
         << "  " << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

float r4vec_product ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PRODUCT returns the product of the entries of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    Input, int N, the number of entries in the vector.
//
//    Input, float A[N], the vector.
//
//    Output, float R4VEC_PRODUCT, the product of the vector.
//
{
  int i;
  float product;

  product = 1.0;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
//****************************************************************************80

void r4vec_range ( int n, float x[], float xmin, float xmax, float y[],
  float *ymin, float *ymax )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_RANGE finds the range of Y's within a restricted X range.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The routine is given a set of pairs of points (X,Y), and a range
//    XMIN to XMAX of valid X values.  Over this range, it seeks
//    YMIN and YMAX, the minimum and maximum values of Y for
//    valid X's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float X[N], the X array.
//
//    Input, float XMIN, XMAX, the range of X values to check.
//
//    Input, float Y[N], the Y array.
//
//    Output, float *YMIN, *YMAX, the range of Y values whose
//    X value is within the X range.
//
{
  int i;

  *ymin =   r4_huge ( );
  *ymax = - r4_huge ( );

  for ( i = 0; i < n; i++ )
  {
    if ( xmin <= x[i] && x[i] <= xmax )
    {
      *ymin = r4_min ( *ymin, y[i] );
      *ymax = r4_max ( *ymax, y[i] );
    }
  }

  return;
}
//****************************************************************************80

void r4vec_range_2 ( int n, float a[], float *amin, float *amax )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_RANGE_2 updates a range to include a new R4VEC
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Given a range AMIN to AMAX, and an array A, the routine will
//    decrease AMIN if necessary, or increase AMAX if necessary, so that
//    every entry of A is between AMIN and AMAX.
//
//    However, AMIN will not be increased, nor AMAX decreased.
//
//    This routine may be used to compute the maximum and minimum of a
//    collection of arrays one at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], the array.
//
//    Input/output, float *AMIN, *AMAX.  On input, the
//    current legal range of values for A.  On output, AMIN and AMAX
//    are either unchanged, or else "widened" so that all entries
//    of A are within the range.
//
{
  *amax = r4_max ( *amax, r4vec_max ( n, a ) );
  *amin = r4_min ( *amin, r4vec_min ( n, a ) );

  return;
}
//****************************************************************************80

void r4vec_reverse ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_REVERSE reverses the elements of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Example:
//
//    Input:
//
//      N = 5, A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
//
//    Output:
//
//      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, float A[N], the array to be reversed.
//
{
  int i;
  float temp;

  for ( i = 1; i <= n/2; i++ )
  {
    temp   = a[i-1];
    a[i-1] = a[n-i];
    a[n-i] = temp;
  }

  return;
}
//****************************************************************************80

void r4vec_rotate ( int n, float a[], int m )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_ROTATE "rotates" the entries of an R4VEC in place.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    This routine rotates an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5, M = 2
//      A    = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
//
//    Output:
//
//      A    = ( 4.0, 5.0, 1.0, 2.0, 3.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int M, the number of positions to the right that
//    each element should be moved.  Elements that shift pass position
//    N "wrap around" to the beginning of the array.
//
//    Input/output, float A[N], the array to be rotated.
//
{
  int iget;
  int iput;
  int istart;
  int mcopy;
  int nset;
  float temp;
//
//  Force M to be positive, between 0 and N-1.
//
  mcopy = i4_modp ( m, n );

  if ( mcopy == 0 )
  {
    return;
  }

  istart = 0;
  nset = 0;

  for ( ; ; )
  {
    istart = istart + 1;

    if ( n < istart )
    {
      break;
    }

    temp = a[istart-1];
    iget = istart;
//
//  Copy the new value into the vacated entry.
//
    for ( ; ; )
    {
      iput = iget;

      iget = iget - mcopy;
      if ( iget < 1 )
      {
        iget = iget + n;
      }

      if ( iget == istart )
      {
        break;
      }

      a[iput-1] = a[iget-1];
      nset = nset + 1;
    }

    a[iput-1] = temp;
    nset = nset + 1;

    if ( n <= nset )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

float r4vec_scalar_triple_product ( float v1[3], float v2[3], float v3[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SCALAR_TRIPLE_PRODUCT computes the scalar triple product.
//
//  Discussion:
//
//    STRIPLE = V1 dot ( V2 x V3 ).
//
//    STRIPLE is the volume of the parallelogram whose sides are
//    formed by V1, V2 and V3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float V1[3], V2[3], V3[3], the three vectors.
//
//    Output, float R4VEC_SCALAR_TRIPLE_PRODUCT, the scalar
//    triple product.
//
{
  float value;

  value =
      v1[0] * ( v2[1] * v3[2] - v2[2] * v3[1] )
    + v1[1] * ( v2[2] * v3[0] - v2[0] * v3[2] )
    + v1[2] * ( v2[0] * v3[1] - v2[1] * v3[0] );

  return value;
}
//****************************************************************************80

int r4vec_search_binary_a ( int n, float a[], float aval )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SEARCH_BINARY_A searches an ascending sorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Binary search is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
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
//    Input, int N, the number of elements in the array.
//
//    Input, float A[N], the array to be searched.  The array must
//    be sorted in ascending order.
//
//    Input, float AVAL, the value to be searched for.
//
//    Output, int R4VEC_SEARCH_BINARY_A, the result of the search.
//    -1, AVAL does not occur in the array.
//    I, A(I) = AVAL.
//
{
  int high;
  int indx;
  int low;
  int mid;

  indx = -1;

  low = 1;
  high = n;

  while ( low <= high )
  {
    mid = ( low + high ) / 2;

    if ( a[mid-1] == aval )
    {
      indx = mid;
      break;
    }
    else if ( a[mid-1] < aval )
    {
      low = mid + 1;
    }
    else if ( aval < a[mid-1] )
    {
      high = mid - 1;
    }
  }

  return indx;
}
//****************************************************************************80

void r4vec_shift ( int shift, int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SHIFT performs a shift on an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SHIFT, the amount by which each entry is to
//    be shifted.
//
//    Input, int N, the length of the vector.
//
//    Input/output, float X[N], the vector to be shifted.
//
{
  int i;
  int ihi;
  int ilo;
  float *y;

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[i];
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  ilo = i4_max ( 0, shift );
  ihi = i4_min ( n, n + shift );

  for ( i = ilo; i < ihi; i++ )
  {
    x[i] = y[i-shift];
  }

  delete [] y;

  return;
}
//****************************************************************************80

void r4vec_shift_circular ( int shift, int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SHIFT_CIRCULAR performs a circular shift on an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SHIFT, the amount by which each entry is to
//    be shifted.
//
//    Input, int N, the length of the vector.
//
//    Input/output, float X[N], the vector to be shifted.
//
{
  int i;
  int j;
  float *y;

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[i];
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_wrap ( i - shift, 0, n - 1 );
    x[i] = y[j];
  }
  delete [] y;
  return;
}
//****************************************************************************80

void r4vec_sort_bubble_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_BUBBLE_A ascending sorts an R4VEC using bubble sort.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input/output, float A[N].
//    On input, an unsorted array of floats.
//    On output, A has been sorted.
//
{
  int i;
  int j;
  float temp;

  for ( i = 0; i < n-1; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( a[j] < a[i] )
      {
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
      }
    }
  }
  return;
}
//****************************************************************************80

void r4vec_sort_bubble_d ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_BUBBLE_D descending sorts an R4VEC using bubble sort.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input/output, float A[N].
//    On input, an unsorted array of floats.
//    On output, A has been sorted.
//
{
  int i;
  int j;
  float temp;

  for ( i = 0; i < n-1; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( a[i] < a[j] )
      {
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
      }
    }
  }
  return;
}
//****************************************************************************80

void r4vec_sort_heap_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_HEAP_A ascending sorts an R4VEC using heap sort.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
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
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, float A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  float temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  r4vec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    r4vec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}
//****************************************************************************80

void r4vec_sort_heap_d ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_HEAP_D descending sorts an R4VEC using heap sort.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
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
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, float A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  float temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into ascending heap form.
//
  r4vec_heap_a ( n, a );
//
//  2: Sort A.
//
//  The smallest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    r4vec_heap_a ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}
//****************************************************************************80

void r4vec_sort_heap_index_a ( int n, float a[], int indx[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R4VEC
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(indx(*))
//
//    or explicitly, by the call
//
//      r4vec_permute ( n, indx, 0, a )
//
//    after which a(*) is sorted.
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], an array to be index-sorted.
//
//    Output, int INDX[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  float aval;
  int i;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return;
}
//****************************************************************************80

int *r4vec_sort_heap_index_a_new ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_HEAP_INDEX_A_NEW does an indexed heap ascending sort of an R4VEC
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(indx(*))
//
//    or explicitly, by the call
//
//      r4vec_permute ( n, indx, 0, a )
//
//    after which a(*) is sorted.
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], an array to be index-sorted.
//
//    Output, int R4VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  float aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
//****************************************************************************80

void r4vec_sort_heap_index_d ( int n, float a[], int indx[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_HEAP_INDEX_D_NEW does an indexed heap descending sort of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(indx(*))
//
//    or explicitly, by the call
//
//      r4vec_permute ( n, indx, 0, a )
//
//    after which a(*) is sorted.
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], an array to be index-sorted.
//
//    Output, int INDX[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  float aval;
  int i;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j]] < a[indx[j-1]] )
        {
          j = j + 1;
        }
      }

      if ( a[indx[j-1]] < aval )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }

    indx[i-1] = indxt;
  }
  return;
}
//****************************************************************************80

int *r4vec_sort_heap_index_d_new ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_HEAP_INDEX_D_NEW does an indexed heap descending sort of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(indx(*))
//
//    or explicitly, by the call
//
//      r4vec_permute ( n, indx, 0, a )
//
//    after which a(*) is sorted.
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], an array to be index-sorted.
//
//    Output, int R4VEC_SORT_HEAP_INDEX_D[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  float aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j]] < a[indx[j-1]] )
        {
          j = j + 1;
        }
      }

      if ( a[indx[j-1]] < aval )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }

    indx[i-1] = indxt;
  }
  return indx;
}
//****************************************************************************80

int *r4vec_sort_heap_mask_a ( int n, float a[], int mask_num, int mask[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_HEAP_MASK_A: indexed heap ascending sort of a masked R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    An array A is given.  An array MASK of indices into A is given.
//    The routine produces a vector INDX, which is a permutation of the
//    entries of MASK, so that:
//
//      A(MASK(INDX(I)) <= A(MASK(INDX(J))
//
//    whenever
//
//      I <= J
//
//    In other words, only the elements of A that are indexed by MASK
//    are to be considered, and the only thing that happens is that
//    a rearrangment of the indices in MASK is returned that orders the
//    masked elements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], an array to be index-sorted.
//
//    Input, int MASK_NUM, the number of mask elements.
//
//    Input, int MASK[MASK_NUM], the mask array.  This is
//    simply a list of indices of A.  The entries of MASK should
//    be unique, and each one should be between 1 and N.
//
//    Output, int INDX[MASK_NUM], the sort index.  There are MASK_NUM
//    elements of A selected by MASK.  If we want to list those elements
//    in order, then the I-th element is A(MASK(INDX(I))).
//
{
  float aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  if ( mask_num < 1 )
  {
    return NULL;
  }

  if ( mask_num == 1 )
  {
    indx = new int[1];
    indx[0] = 1;
    return indx;
  }

  indx = i4vec_indicator1_new ( mask_num );

  l = mask_num / 2 + 1;
  ir = mask_num;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[mask[indxt-1]-1];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[mask[indxt-1]-1];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[mask[indx[j-1]-1]-1] < a[mask[indx[j]-1]-1] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[mask[indx[j-1]-1]-1] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
//****************************************************************************80

void r4vec_sort_insert_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_INSERT_A ascending sorts an R4VEC using an insertion sort.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Kreher, Douglas Simpson,
//    Algorithm 1.1,
//    Combinatorial Algorithms,
//    CRC Press, 1998, page 11.
//
//  Parameters:
//
//    Input, int N, the number of items in the vector.
//    N must be positive.
//
//    Input/output, float A[N].
//
//    On input, A contains data to be sorted.
//    On output, the entries of A have been sorted in ascending order.
//
{
  int i;
  int j;
  float x;

  for ( i = 1; i < n; i++ )
  {
    x = a[i];

    j = i;

    while ( 1 <= j && x < a[j-1] )
    {
      a[j] = a[j-1];
      j = j - 1;
    }
    a[j] = x;
  }

  return;
}
//****************************************************************************80

int *r4vec_sort_insert_index_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_INSERT_INDEX_A ascending index sorts an R4VEC using insertion.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Kreher, Douglas Simpson,
//    Combinatorial Algorithms,
//    CRC Press, 1998, page 11.
//
//  Parameters:
//
//    Input, int N, the number of items in the vector.
//    N must be positive.
//
//    Input, float A[N], the array to be sorted.
//
//    Output, int R4VEC_SORT_INSET_INDEX_A[N], the sorted indices.  The array
//    is sorted when listed from A(INDX(1)) through A(INDX(N)).
//
{
  int i;
  int *indx;
  int j;
  float x;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = i4vec_indicator1_new ( n );

  for ( i = 2; i <= n; i++ )
  {
    x = a[i-1];

    j = i - 1;

    while ( 1 <= j )
    {
      if ( a[indx[j-1]-1] <= x )
      {
        break;
      }

      indx[j] = indx[j-1];
      j = j - 1;
    }
    indx[j] = i;
  }

  return indx;
}
//****************************************************************************80

void r4vec_sort_quick_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_QUICK_A ascending sorts an R4VEC using quick sort.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Example:
//
//    Input:
//
//      N = 7
//
//      A = ( 6, 7, 3, 2, 9, 1, 8 )
//
//    Output:
//
//      A = ( 1, 2, 3, 6, 7, 8, 9 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, float A[N].  On input, the array to be sorted.
//    On output, A has been reordered into ascending order.
//
{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R4VEC_SORT_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }
  else if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[0] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
//
//  Partition the segment.
//
    r4vec_part_quick_a ( n_segment, a+base-1, &l_segment, &r_segment );
//
//  If the left segment has more than one element, we need to partition it.
//
    if ( 1 < l_segment )
    {

      if ( LEVEL_MAX < level )
      {
        cerr << "\n";
        cerr << "R4VEC_SORT_QUICK_A - Fatal error!\n";
        cerr << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
//
//  The left segment and the middle segment are sorted.
//  Must the right segment be partitioned?
//
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
//
//  Otherwise, we back up a level if there is an earlier one.
//
    else
    {
      for ( ; ; )
      {
        if ( 1 < level )
        {
          base = rsave[level-1];
          n_segment = rsave[level-2] - rsave[level-1];
          level = level - 1;
          if ( 0 < n_segment )
          {
            break;
          }
        }
        else
        {
          n_segment = 0;
          break;
        }
      }
    }
  }

  return;
# undef LEVEL_MAX
}
//****************************************************************************80

void r4vec_sort_shell_a ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORT_SHELL_A ascending sorts an R4VEC using Shell's sort.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, float A[N].
//    On input, an array to be sorted.
//    On output, the sorted array.
//
{
  float asave;
  int i;
  int ifree;
  int inc;
  int ipow;
  int j;
  int k;
  int maxpow;
  int test;

  if ( n <= 1 )
  {
    return;
  }
//
//  Determine the smallest MAXPOW so that
//    N <= ( 3**MAXPOW - 1 ) / 2
//
  maxpow = 1;
  test = 3;

  while ( test < 2 * n + 1 )
  {
    maxpow = maxpow + 1;
    test = test * 3;
  }

  if ( 1 < maxpow )
  {
    maxpow = maxpow - 1;
    test = test / 3;
  }
//
//  Now sort groups of size ( 3**IPOW - 1 ) / 2.
//
  for ( ipow = maxpow; 1 <= ipow; ipow-- )
  {
    inc = ( test - 1 ) / 2;
    test = test / 3;
//
//  Sort the values with indices equal to K mod INC.
//
    for ( k = 1; k <= inc; k++ )
    {
//
//  Insertion sort of the items with index
//  INC+K, 2*INC+K, 3*INC+K, ...
//
      for ( i = inc+k; i <= n; i = i + inc )
      {
        asave = a[i-1];
        ifree = i;
        j = i - inc;

        for ( ; ; )
        {
          if ( j < 1 )
          {
            break;
          }

          if ( a[j-1] <= asave )
          {
            break;
          }

          ifree = j;
          a[j+inc-1] = a[j-1];
          j = j - inc;
        }
        a[ifree-1] = asave;
      }
    }
  }

  return;
}
//****************************************************************************80

float *r4vec_sorted_merge_a ( int na, float a[], int nb, float b[], int *nc )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORTED_MERGE_A merges two ascending sorted R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The elements of A and B should be sorted in ascending order.
//
//    The elements in the output array C will also be in ascending order,
//    and unique.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, float A[NA], the first sorted array.
//
//    Input, int NB, the dimension of B.
//
//    Input, float B[NB], the second sorted array.
//
//    Output, int *NC, the number of entries in the merged vector.
//
//    Output, float R4VEC_SORTED_MERGE_A[NC], the merged unique sorted array.
//
{
  float *c;
  float *d;
  int j;
  int ja;
  int jb;
  int na2;
  int nb2;
  int nd;
  int order;

  na2 = na;
  nb2 = nb;

  ja = 0;
  jb = 0;
  *nc = 0;
  nd = 0;
  d = new float[na+nb];

  order = r4vec_order_type ( na2, a );

  if ( order < 0 || 2 < order )
  {
    cerr << "\n";
    cerr << "R4VEC_SORTED_MERGE_A - Fatal error!\n";
    cerr << "  The input array A is not ascending sorted.\n";
    return NULL;
  }

  order = r4vec_order_type ( nb2, b );

  if ( order < 0 || 2 < order )
  {
    cerr << "\n";
    cerr << "R4VEC_SORTED_MERGE_A - Fatal error!\n";
    cerr << "  The input array B is not ascending sorted.\n";
    return NULL;
  }

  for ( ; ; )
  {
//
//  If we've used up all the entries of A, stick the rest of B on the end.
//
    if ( na2 <= ja )
    {
      for ( j = 1; j <= nb2 - jb; j++ )
      {
        jb = jb + 1;
        if ( nd == 0 )
        {
          nd = nd + 1;
          d[nd-1] = b[jb-1];
        }
        else if ( d[nd-1] < b[jb-1] )
        {
          nd = nd + 1;
          d[nd-1] = b[jb-1];
        }
      }
      break;
    }
//
//  If we've used up all the entries of B, stick the rest of A on the end.
//
    else if ( nb2 <= jb )
    {
      for ( j = 1; j <= na2 - ja; j++ )
      {
        ja = ja + 1;
        if ( nd == 0 )
        {
          nd = nd + 1;
          d[nd-1] = a[ja-1];
        }
        else if ( d[nd-1] < a[ja-1] )
        {
          nd = nd + 1;
          d[nd-1] = a[ja-1];
        }
      }
      break;
    }
//
//  Otherwise, if the next entry of A is smaller, that's our candidate.
//
    else if ( a[ja] <= b[jb] )
    {
      ja = ja + 1;
      if ( nd == 0 )
      {
        nd = nd + 1;
        d[nd-1] = a[ja-1];
      }
      else if ( d[nd-1] < a[ja-1] )
      {
        nd = nd + 1;
        d[nd-1] = a[ja-1];
      }
    }
//
//  ...or if the next entry of B is the smaller, consider that.
//
    else
    {
      jb = jb + 1;
      if ( nd == 0 )
      {
        nd = nd + 1;
        d[nd-1] = b[jb-1];
      }
      else if ( d[nd-1] < b[jb-1] )
      {
        nd = nd + 1;
        d[nd-1] = b[jb-1];
      }
    }
  }

  *nc = nd;

  c = r4vec_copy_new ( nd, d );

  delete [] d;

  return c;
}
//****************************************************************************80

int r4vec_sorted_nearest ( int n, float a[], float value )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORTED_NEAREST returns the nearest element in a sorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, float A[N], a sorted vector.
//
//    Input, float VALUE, the value whose nearest vector entry is sought.
//
//    Output, int R4VEC_SORTED_NEAREST, the index of the nearest
//    entry in the vector.
//
{
  int hi;
  int lo;
  int mid;

  if ( n < 1 )
  {
    return (-1);
  }

  if ( n == 1 )
  {
    return 1;
  }

  if ( a[0] < a[n-1] )
  {
    if ( value < a[0] )
    {
      return 1;
    }
    else if ( a[n-1] < value )
    {
      return n;
    }
//
//  Seek an interval containing the value.
//
    lo = 1;
    hi = n;

    while ( lo < hi - 1 )
    {
      mid = ( lo + hi ) / 2;

      if ( value == a[mid-1] )
      {
        return mid;
      }
      else if ( value < a[mid-1] )
      {
        hi = mid;
      }
      else
      {
        lo = mid;
      }
    }
//
//  Take the nearest.
//
    if ( r4_abs ( value - a[lo-1] ) < r4_abs ( value - a[hi-1] ) )
    {
      return lo;
    }
    else
    {
      return hi;
    }
  }
//
//  A descending sorted vector A.
//
  else
  {
    if ( value < a[n-1] )
    {
      return n;
    }
    else if ( a[0] < value )
    {
      return 1;
    }
//
//  Seek an interval containing the value.
//
    lo = n;
    hi = 1;

    while ( lo < hi - 1 )
    {
      mid = ( lo + hi ) / 2;

      if ( value == a[mid-1] )
      {
        return mid;
      }
      else if ( value < a[mid-1] )
      {
        hi = mid;
      }
      else
      {
        lo = mid;
      }
    }
//
//  Take the nearest.
//
    if ( r4_abs ( value - a[lo-1] ) < r4_abs ( value - a[hi-1] ) )
    {
      return lo;
    }
    else
    {
      return hi;
    }
  }
}
//****************************************************************************80

void r4vec_sorted_range ( int n, float r[], float r_lo, float r_hi,
  int *i_lo, int *i_hi )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORTED_RANGE searches a sorted vector for elements in a range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items in the vector.
//
//    Input, float R[N], the sorted vector.
//
//    Input, float R_LO, R_HI, the limits of the range.
//
//    Output, int *I_LO, *I_HI, the range of indices
//    so that I_LO <= I <= I_HI => R_LO <= R(I) <= R_HI.  If no
//    values in R lie in the range, then I_HI < I_LO will be returned.
//
{
  int i1;
  int i2;
  int j1;
  int j2;
//
//  Cases we can handle immediately.
//
  if ( r[n-1] < r_lo )
  {
    *i_lo = - 1;
    *i_hi = - 2;
    return;
  }

  if ( r_hi < r[0] )
  {
    *i_lo = - 1;
    *i_hi = - 2;
    return;
  }
//
//  Are there are least two intervals?
//
  if ( n == 1 )
  {
    if ( r_lo <= r[0] && r[0] <= r_hi )
    {
      *i_lo = 1;
      *i_hi = 1;
    }
    else
    {
      *i_lo = - 1;
      *i_hi = - 2;
    }
    return;
  }
//
//  Bracket R_LO.
//
  if ( r_lo <= r[0] )
  {
    *i_lo = 0;
  }
  else
  {
//
//  R_LO is in one of the intervals spanned by R(J1) to R(J2).
//  Examine the intermediate interval [R(I1), R(I1+1)].
//  Does R_LO lie here, or below or above?
//
    j1 = 0;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_lo < r[i1] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[i2] < r_lo )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_lo = i1;
        break;
      }
    }
  }
//
//  Bracket R_HI
//
  if ( r[n-1] <= r_hi )
  {
    *i_hi = n - 1;
  }
  else
  {
    j1 = *i_lo;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_hi < r[i1] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[i2] < r_hi )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_hi = i2;
        break;
      }
    }
  }
//
//  We expect to have computed the largest I_LO and smallest I_HI such that
//    R(I_LO) <= R_LO <= R_HI <= R(I_HI)
//  but what we want is actually
//    R_LO <= R(I_LO) <= R(I_HI) <= R_HI
//  which we can usually get simply by incrementing I_LO and decrementing I_HI.
//
  if ( r[*i_lo] < r_lo )
  {
    *i_lo = *i_lo + 1;
    if ( n - 1 < *i_lo )
    {
      *i_hi = *i_lo - 1;
    }
  }

  if ( r_hi < r[*i_hi] )
  {
    *i_hi = *i_hi - 1;
    if ( *i_hi < 0 )
    {
      *i_lo = *i_hi + 1;
    }
  }

  return;
}
//****************************************************************************80

void r4vec_sorted_split ( int n, float a[], float split, int *i_lt,
  int *i_gt )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORTED_SPLIT "splits" a sorted R4VEC, given a splitting value.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Given a splitting value SPLIT, the routine seeks indices
//    I_LT and I_GT so that
//
//      A(I_LT) < SPLIT < A(I_GT),
//
//    and if there are intermediate index values between I_LT and
//    I_GT, then those entries of A are exactly equal to SPLIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int N, the number of entries in A.
//
//    Input, float A[N], a sorted array.
//
//    Input, float SPLIT, a value to which the entries in A are
//    to be compared.
//
//    Output, int *I_LT:
//    0 if no entries are less than SPLIT;
//    N if all entries are less than SPLIT;
//    otherwise, the index of the last entry in A less than SPLIT.
//
//    Output, int *I_GT:
//    1 if all entries are greater than SPLIT;
//    N+1 if no entries are greater than SPLIT;
//    otherwise the index of the first entry in A greater than SPLIT.
//
{
  int hi;
  int i;
  int lo;
  int mid;

  if ( n < 1 )
  {
    *i_lt = -1;
    *i_gt = -1;
    return;
  }

  if ( split < a[0] )
  {
    *i_lt = 0;
    *i_gt = 1;
    return;
  }

  if ( a[n-1] < split )
  {
    *i_lt = n;
    *i_gt = n + 1;
    return;
  }

  lo = 1;
  hi = n;

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      *i_lt = lo;
      break;
    }

    mid = ( lo + hi ) / 2;

    if ( split <= a[mid-1] )
    {
      hi = mid;
    }
    else
    {
      lo = mid;
    }
  }

  for ( i = *i_lt + 1; i <= n; i++ )
  {
    if ( split < a[i-1] )
    {
      *i_gt = i;
      return;
    }
  }

  *i_gt = n + 1;

  return;
}
//****************************************************************************80

void r4vec_sorted_undex ( int x_num, float x_val[], int x_unique_num,
  float tol, int undx[], int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORTED_UNDEX returns unique sorted indexes for a sorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The goal of this routine is to determine a vector UNDX,
//    which points, to the unique elements of X, in sorted order,
//    and a vector XDNU, which identifies, for each entry of X, the index of
//    the unique sorted element of X.
//
//    This is all done with index vectors, so that the elements of
//    X are never moved.
//
//    Assuming X is already sorted, we examine the entries of X in order,
//    noting the unique entries, creating the entries of XDNU and
//    UNDX as we go.
//
//    Once this process has been completed, the vector X could be
//    replaced by a compressed vector XU, containing the unique entries
//    of X in sorted order, using the formula
//
//      XU(I) = X(UNDX(I)).
//
//    We could then, if we wished, reconstruct the entire vector X, or
//    any element of it, by index, as follows:
//
//      X(I) = XU(XDNU(I)).
//
//    We could then replace X by the combination of XU and XDNU.
//
//    Later, when we need the I-th entry of X, we can locate it as
//    the XDNU(I)-th entry of XU.
//
//    Here is an example of a vector X, the sort and inverse sort
//    index vectors, and the unique sort and inverse unique sort vectors
//    and the compressed unique sorted vector.
//
//      I      X      XU  Undx  Xdnu
//    ----+------+------+-----+-----+
//      0 | 11.0 |  11.0    0     0
//      1 | 11.0 |  22.0    4     0
//      2 | 11.0 |  33.0    7     0
//      3 | 11.0 |  55.0    8     0
//      4 | 22.0 |                1
//      5 | 22.0 |                1
//      6 | 22.0 |                1
//      7 | 33.0 |                2
//      8 | 55.0 |                3
//
//    INDX(2) = 3 means that sorted item(2) is X(3).
//    XDNI(2) = 5 means that X(2) is sorted item(5).
//
//    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
//    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
//
//    XU(XDNU(I))) = X(I).
//    XU(I)        = X(UNDX(I)).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of data values.
//
//    Input, float X_VAL[X_NUM], the data values.
//
//    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
//    This value is only required for languages in which the size of
//    UNDX must be known in advance.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.
//
//    Output, int XDNU[X_NUM], the XDNU vector.
//
{
  int i;
  int j;
//
//  Walk through the sorted array X.
//
  i = 0;

  j = 0;
  undx[j] = i;

  xdnu[i] = j;

  for ( i = 1; i < x_num; i++ )
  {
    if ( tol < r4_abs ( x_val[i] - x_val[undx[j]] ) )
    {
      j = j + 1;
      undx[j] = i;
    }
    xdnu[i] = j;
  }

  return;
}
//****************************************************************************80

float *r4vec_sorted_unique ( int n, float a[], float tol, int *unique_num )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORTED_UNIQUE finds the unique elements in a sorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    If the data is not sorted, the results of the routine will
//    be garbage.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, float A[N], the sorted array of N elements;
//
//    Input, float TOL, a tolerance for checking equality.
//
//    Output, int *UNIQUE_NUM, the number of unique elements of A.
//
//    Output, float R4VEC_SORTED_UNIQUE[UNIQUE_NUM], the unique elements of A.
//
{
  float *a_unique;
  int i;
  int iuniq;

  *unique_num = 0;

  if ( n <= 0 )
  {
    return NULL;
  }
//
//  Determine the number of unique elements.
//
  iuniq = 0;
  *unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < r4_abs ( a[i] - a[iuniq] ) )
    {
       iuniq = i;
      *unique_num = *unique_num + 1;
    }
  }
//
//  Set aside space for the unique elements.
//
  a_unique = new float[*unique_num];
//
//  Repeat the search, but now store the unique elements.
//
  *unique_num = 0;

  a_unique[*unique_num] = a[0];
  *unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < r4_abs ( a[i] - a_unique[*unique_num-1] ) )
    {
      a_unique[*unique_num] = a[i];
      *unique_num = *unique_num + 1;
    }
  }

  return a_unique;
}
//****************************************************************************80

int r4vec_sorted_unique_count ( int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Because the array is sorted, this algorithm is O(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, float A[N], the sorted array to examine.
//
//    Input, float TOL, a tolerance for checking equality.
//
//    Output, int R4VEC_SORTED_UNIQUE_COUNT, the number of unique elements of A.
//
{
  int i;
  int unique_num;

  unique_num = 0;

  if ( n < 1 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < r4_abs ( a[i-1] - a[i] ) )
    {
      unique_num = unique_num + 1;
    }
  }

  return unique_num;
}
//****************************************************************************80

void r4vec_sorted_unique_hist ( int n, float a[], float tol, int maxuniq,
  int *unique_num, float auniq[], int acount[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SORTED_UNIQUE_HIST histograms unique elements of a sorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, float A[N], the array to examine, which must have been
//    sorted.
//
//    Input, float TOL, a tolerance for checking equality.
//
//    Input, int MAXUNIQ, the maximum number of unique elements
//    that can be handled.  If there are more than MAXUNIQ unique
//    elements in A, the excess will be ignored.
//
//    Output, int *UNIQUE_NUM, the number of unique elements of A.
//
//    Output, float AUNIQ[UNIQUE_NUM], the unique elements of A.
//
//    Output, int ACOUNT[UNIQUE_NUM], the number of times each element
//    of AUNIQ occurs in A.
//
{
  int i;
  int index;
//
//  Start taking statistics.
//
  index = -1;

  for ( i = 0; i < n; i++ )
  {

    if ( i == 0 )
    {
      index = 0;
      auniq[index] = a[0];
      acount[index] = 1;
    }
    else if ( r4_abs ( a[i] - auniq[index] ) <= tol )
    {
      acount[index] = acount[index] + 1;
    }
    else if ( index + 1 < maxuniq )
    {
      index = index + 1;
      auniq[index] = a[i];
      acount[index] = 1;
    }
  }

  *unique_num = index + 1;

  return;
}
//****************************************************************************80

int r4vec_split ( int n, float a[], float split )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SPLIT "splits" an unsorted R4VEC based on a splitting value.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    If the vector is already sorted, it is simpler to do a binary search
//    on the data than to call this routine.
//
//    The vector is not assumed to be sorted before input, and is not
//    sorted during processing.  If sorting is not needed, then it is
//    more efficient to use this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input/output, float A[N], the array to split.  On output,
//    all the entries of A that are less than or equal to SPLIT
//    are in A(1:ISPLIT).
//
//    Input, float SPLIT, the value used to split the vector.
//    It is not necessary that any value of A actually equal SPLIT.
//
//    Output, int R4VEC_SPLIT, indicates the position of the last
//    entry of the split vector that is less than or equal to SPLIT.
//
{
  int i;
  int i1;
  int i2;
  int i3;
  int isplit;
  int j1;
  int j2;
  int j3;
  float temp;
//
//  Partition the vector into A1, A2, A3, where
//    A1 = A(I1:J1) holds values <= SPLIT,
//    A2 = A(I2:J2) holds untested values,
//    A3 = A(I3:J3) holds values > SPLIT.
//
  i1 = 1;
  j1 = 0;

  i2 = 1;
  j2 = n;

  i3 = n+1;
  j3 = n;
//
//  Pick the next item from A2, and move it into A1 or A3.
//  Adjust indices appropriately.
//
  for ( i = 1; i <= n; i++ )
  {
    if ( a[i2-1] <= split )
    {
      i2 = i2 + 1;
      j1 = j1 + 1;
    }
    else
    {
      temp = a[i2-1];
      a[i2-1] = a[i3-2];
      a[i3-2] = temp;
      i3 = i3 - 1;
      j2 = j2 - 1;
    }
  }

  isplit = j1;

  return isplit;
}
//****************************************************************************80

float r4vec_std ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_STD returns the standard deviation of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    14 April 2003
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
//    Input, float A[N], the vector.
//
//    Output, float R4VEC_STD, the standard deviation of the vector.
//
{
  int i;
  float mean;
  float std;

  if ( n < 2 )
  {
    std = 0.0;
  }
  else
  {
    mean = 0.0;
    for ( i = 0; i < n; i++ )
    {
      mean = mean + a[i];
    }
    mean = mean / ( ( float ) n );

    std = 0.0;
    for ( i = 0; i < n; i++ )
    {
      std = std + ( a[i] - mean ) * ( a[i] - mean );
    }
    std = sqrt ( std / ( ( float ) ( n - 1 ) ) );
  }

  return std;
}
//****************************************************************************80

void r4vec_stutter ( int n, float a[], int m, float am[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_STUTTER makes a "stuttering" copy of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
//    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the input vector.
//
//    Input, float A[N], the vector.
//
//    Input, int M, the "stuttering factor".
//
//    Output, float AM[M*N], the stuttering vector.
//
{
  int i;
  int j;
  int k;

  k = 0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      am[k] = a[i];
      k = k + 1;
    }
  }
  return;
}
//****************************************************************************80

float *r4vec_stutter_new ( int n, float a[], int m )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_STUTTER_NEW makes a "stuttering" copy of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
//    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the input vector.
//
//    Input, float A[N], the vector.
//
//    Input, int M, the "stuttering factor".
//
//    Output, float R4VEC_STUTTER_NEW[M*N], the stuttering vector.
//
{
  float *am;
  int i;
  int j;
  int k;

  am = new float[m*n];

  k = 0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      am[k] = a[i];
      k = k + 1;
    }
  }
  return am;
}
//****************************************************************************80

float r4vec_sum ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SUM returns the sum of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    Input, float A[N], the vector.
//
//    Output, float R4VEC_SUM, the sum of the vector.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
//****************************************************************************80

void r4vec_swap ( int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SWAP swaps the entries of two R4VEC's.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the arrays.
//
//    Input/output, float A1[N], A2[N], the vectors to swap.
//
{
  int i;
  float temp;

  for ( i = 0; i < n; i++ )
  {
    temp  = a1[i];
    a1[i] = a2[i];
    a2[i] = temp;
  }

  return;
}
//****************************************************************************80

void r4vec_transpose_print ( int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_TRANSPOSE_PRINT prints an R4VEC "transposed".
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Example:
//
//    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
//    TITLE = 'My vector:  '
//
//    My vector:
//
//        1.0    2.1    3.2    4.3    5.4
//        6.5    7.6    8.7    9.8   10.9
//       11.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, float A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int ihi;
  int ilo;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= 0 )
  {
    cout << "  (Empty)\n";
    return;
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      cout << "  " << setw(12) << a[i];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void r4vec_undex ( int x_num, float x_val[], int x_unique_num, float tol,
  int undx[], int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNDEX returns unique sorted indexes for an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    The goal of this routine is to determine a vector UNDX,
//    which points, to the unique elements of X, in sorted order,
//    and a vector XDNU, which identifies, for each entry of X, the index of
//    the unique sorted element of X.
//
//    This is all done with index vectors, so that the elements of
//    X are never moved.
//
//    The first step of the algorithm requires the indexed sorting
//    of X, which creates arrays INDX and XDNI.  (If all the entries
//    of X are unique, then these arrays are the same as UNDX and XDNU.)
//
//    We then use INDX to examine the entries of X in sorted order,
//    noting the unique entries, creating the entries of XDNU and
//    UNDX as we go.
//
//    Once this process has been completed, the vector X could be
//    replaced by a compressed vector XU, containing the unique entries
//    of X in sorted order, using the formula
//
//      XU(*) = X(UNDX(*)).
//
//    We could then, if we wished, reconstruct the entire vector X, or
//    any element of it, by index, as follows:
//
//      X(I) = XU(XDNU(I)).
//
//    We could then replace X by the combination of XU and XDNU.
//
//    Later, when we need the I-th entry of X, we can locate it as
//    the XDNU(I)-th entry of XU.
//
//    Here is an example of a vector X, the sort and inverse sort
//    index vectors, and the unique sort and inverse unique sort vectors
//    and the compressed unique sorted vector.
//
//      I     X  Indx  Xdni       XU  Undx  Xdnu
//    ----+-----+-----+-----+--------+-----+-----+
//      0 | 11.     0     0 |    11.     0     0
//      1 | 22.     2     4 |    22.     1     1
//      2 | 11.     5     1 |    33.     3     0
//      3 | 33.     8     7 |    55.     4     2
//      4 | 55.     1     8 |                  3
//      5 | 11.     6     2 |                  0
//      6 | 22.     7     5 |                  1
//      7 | 22.     3     6 |                  1
//      8 | 11.     4     3 |                  0
//
//    INDX(2) = 3 means that sorted item(2) is X(3).
//    XDNI(2) = 5 means that X(2) is sorted item(5).
//
//    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
//    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
//
//    XU(XDNU(I))) = X(I).
//    XU(I)        = X(UNDX(I)).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of data values.
//
//    Input, float X_VAL[X_NUM], the data values.
//
//    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
//    This value is only required for languages in which the size of
//    UNDX must be known in advance.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.
//
//    Output, int XDNU[X_NUM], the XDNU vector.
//
{
  int i;
  int *indx;
  int j;
//
//  Implicitly sort the array.
//
  indx = r4vec_sort_heap_index_a_new ( x_num, x_val );
//
//  Walk through the implicitly sorted array X.
//
  i = 0;

  j = 0;
  undx[j] = indx[i];

  xdnu[indx[i]] = j;

  for ( i = 1; i < x_num; i++ )
  {
    if ( tol < r4_abs ( x_val[indx[i]] - x_val[undx[j]] ) )
    {
      j = j + 1;
      undx[j] = indx[i];
    }
    xdnu[indx[i]] = j;
  }
  delete [] indx;

  return;
}
//****************************************************************************80

void r4vec_uniform ( int n, float b, float c, int &seed, float r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIFORM returns a scaled pseudorandom R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    23 April 2008
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
//    Input, float B, C, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, float R[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = b + ( c - b ) * ( float ) ( seed ) * 4.656612875E-10;
  }

  return;
}
//****************************************************************************80

float *r4vec_uniform_new ( int n, float b, float c, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIFORM_NEW returns a scaled pseudorandom R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    23 April 2008
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
//    Input, float B, C, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, float R4VEC_UNIFORM_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  float *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_UNIFORM_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new float[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = b + ( c - b ) * ( float ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void r4vec_uniform_01 ( int n, int &seed, float r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIFORM_01 returns a unit unit pseudorandom R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    22 May 2006
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
//    Output, float R[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( float ) ( seed ) * 4.656612875E-10;
  }

  return;
}
//****************************************************************************80

float *r4vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIFORM_01_NEW returns a unit unit pseudorandom R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    22 May 2006
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
//    Output, float R4VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  float *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new float[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( float ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

int r4vec_unique_count ( int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIQUE_COUNT counts the unique elements in an unsorted R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    Because the array is unsorted, this algorithm is O(N^2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, float A[N], the array to examine, which does NOT have to
//    be sorted.
//
//    Input, float TOL, a tolerance for checking equality.
//
//    Output, int R4VEC_UNIQUE_COUNT, the number of unique elements of A.
//
{
  int i;
  int j;
  int unique_num;

  unique_num = 0;

  for ( i = 0; i < n; i++ )
  {
    unique_num = unique_num + 1;

    for ( j = 0; j < i; j++ )
    {
      if ( r4_abs ( a[i] - a[j] ) <= tol )
      {
        unique_num = unique_num - 1;
        break;
      }
    }
  }
  return unique_num;
}
//****************************************************************************80

int *r4vec_unique_index ( int n, float a[], float tol )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIQUE_INDEX indexes the unique occurrence of values in an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//    For element A(I) of the vector, UNIQUE_INDEX(I) is the uniqueness index
//    of A(I).  That is, if A_UNIQUE contains the unique elements of A,
//    gathered in order, then
//
//      A_UNIQUE ( UNIQUE_INDEX(I) ) = A(I)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, float A[N], the unsorted array to examine.
//
//    Input, float TOL, a tolerance for equality.
//
//    Output, int R4VEC_UNIQUE_INDEX[N], the unique index.
//
{
  int i;
  int j;
  int *unique_index;
  int unique_num;

  unique_index = new int[n];

  for ( i = 0; i < n; i++ )
  {
    unique_index[i] = -1;
  }
  unique_num = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( unique_index[i] == -1 )
    {
      unique_index[i] = unique_num;
      for ( j = i + 1; j < n; j++ )
      {
        if ( r4_abs ( a[i] - a[j] ) <= tol )
        {
          unique_index[j] = unique_num;
        }
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_index;
}
//****************************************************************************80

float r4vec_variance ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_VARIANCE returns the variance of an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[N], the vector whose variance is desired.
//
//    Output, float R4VEC_VARIANCE, the variance of the vector entries.
//
{
  int i;
  float mean;
  float variance;

  mean = r4vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ )
  {
    variance = variance + ( x[i] - mean ) * ( x[i] - mean );
  }

  if ( 1 < n )
  {
    variance = variance / ( float ) ( n - 1 );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
}
//****************************************************************************80

float *r4vec_vector_triple_product ( float v1[3], float v2[3], float v3[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_VECTOR_TRIPLE_PRODUCT computes the vector triple product.
//
//  Discussion:
//
//    VTRIPLE = V1 x (V2 x V3)
//
//    VTRIPLE is a vector perpendicular to V1, lying in the plane
//    spanned by V2 and V3.  The norm of VTRIPLE is the product
//    of the norms of V1, V2 and V3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float V1[3], V2[3], V3[3], the coordinates
//    of the three vectors.
//
//    Output, float R4VEC_VECTOR_TRIPLE_PRODUCT[3], the vector triple product.
//
{
  float *v123;
  float *v23;

  v23 = r4vec_cross_product_3d ( v2, v3 );

  v123 = r4vec_cross_product_3d ( v1, v23 );

  delete [] v23;

  return v123;
}
//****************************************************************************80

void r4vec_write ( int n, float r[], string output_file )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_WRITE writes an R4VEC to a file.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, float R[N], the vector to be written.
//
//    Input, string OUTPUT_FILE, the name of the file to which
//    the information is to be written.
//
{
  int i;
  ofstream output;

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R4VEC_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    output << "  " << setw(16) << r[i] << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80

void r4vec_zero ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_ZERO zeroes an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    Input, int N, the number of entries in the vector.
//
//    Output, float A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
//****************************************************************************80

float *r4vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_ZERO_NEW creates and zeroes an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, float R4VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  float *a;
  int i;

  a = new float[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
}
//****************************************************************************80

int r4vec2_compare ( int n, float a1[], float a2[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_COMPARE compares two elements of an R4VEC2.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data items.
//
//    Input, float A1[N], A2[N], contain the two components of each item.
//
//    Input, int I, J, the items to be compared.  These values will be
//    1-based indices for the arrays A1 and A2.
//
//    Output, int R4VEC2_COMPARE, the results of the comparison:
//    -1, item I < item J,
//     0, item I = item J,
//    +1, item J < item I.
//
{
  int isgn;

  isgn = 0;

  if ( a1[i-1] < a1[j-1] )
  {
    isgn = -1;
  }
  else if ( a1[i-1] == a1[j-1] )
  {
    if ( a2[i-1] < a2[j-1] )
    {
      isgn = -1;
    }
    else if ( a2[i-1] < a2[j-1] )
    {
      isgn = 0;
    }
    else if ( a2[j-1] < a2[i-1] )
    {
      isgn = +1;
    }
  }
  else if ( a1[j-1] < a1[i-1] )
  {
    isgn = +1;
  }

  return isgn;
}
//****************************************************************************80

void r4vec2_print ( int n, float a1[], float a2[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_PRINT prints an R4VEC2.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, float A1[N], float A2[N], the vectors to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(6)  << i
         << ": " << setw(14) << a1[i]
         << "  " << setw(14) << a2[i] << "\n";
  }

  return;
}
//****************************************************************************80

void r4vec2_print_some ( int n, float x1[], float x2[], int max_print,
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_PRINT_SOME prints "some" of an R4VEC2.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
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
//    Input, float X1[N], X2[N], the vector to be printed.
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
      cout << setw(6)  << i << ": "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print-2; i++ )
    {
      cout << setw(6)  << i     << ": "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
    cout << "......  ..............  ..............\n";
    i = n - 1;
    cout << setw(6)  << i     << ": "
         << setw(14) << x1[i] << "  "
         << setw(14) << x2[i] << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << setw(6)  << i     << ": "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
    i = max_print - 1;
    cout << setw(6)  << i     << ": "
         << setw(14) << x1[i] << "  "
         << setw(14) << x2[i] << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

void r4vec2_sort_a ( int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_SORT_A ascending sorts an R4VEC2.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//    Each item to be sorted is a pair of reals (X,Y), with the X
//    and Y values stored in separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items of data.
//
//    Input/output, float A1[N], A2[N], the data to be sorted.
//
{
  int i;
  int indx;
  int isgn;
  int j;
  float temp;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      temp    = a1[i-1];
      a1[i-1] = a1[j-1];
      a1[j-1] = temp;

      temp    = a2[i-1];
      a2[i-1] = a2[j-1];
      a2[j-1] = temp;
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = r4vec2_compare ( n, a1, a2, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void r4vec2_sort_d ( int n, float a1[], float a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_SORT_D descending sorts an R4VEC2.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//    Each item to be sorted is a pair of reals (X,Y), with the X
//    and Y values stored in separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items of data.
//
//    Input/output, float A1[N], A2[N], the data to be sorted.
//
{
  int i;
  int indx;
  int isgn;
  int j;
  float temp;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      temp    = a1[i-1];
      a1[i-1] = a1[j-1];
      a1[j-1] = temp;

      temp    = a2[i-1];
      a2[i-1] = a2[j-1];
      a2[j-1] = temp;
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = - r4vec2_compare ( n, a1, a2, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

int *r4vec2_sort_heap_index_a ( int n, int base, float x[], float y[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R4VEC2.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    ( X(I), Y(I) ) < ( X(J), Y(J) ) if:
//
//    * X(I) < X(J), or
//
//    * X(I) = X(J), and Y(I) < Y(J).
//
//    Once the index array is computed, the sorting can be carried out
//    implicitly:
//
//      ( x(indx(*)), y(indx(*) )
//
//    or explicitly, by the calls
//
//      r4vec_permute ( n, indx, 0, x )
//      r4vec_permute ( n, indx, 0, y )
//
//    after which ( x(*), y(*) ), is sorted.
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int BASE, the desired indexing for the sort index:
//    0 for 0-based indexing,
//    1 for 1-based indexing.
//
//    Input, float X[N], Y[N], pairs of X, Y coordinates of points.
//
//    Output, int INDX[N], the sort index.  The
//    I-th element of the sorted array has coordinates
//    ( X(INDX(I)), Y(INDX(I) ).
//
{
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;
  float xval;
  float yval;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0] + base;
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      xval = x[indxt];
      yval = y[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      xval = x[indxt];
      yval = y[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( x[indx[j-1]] < x[indx[j]] ||
          ( x[indx[j-1]] == x[indx[j]] && y[indx[j-1]] < y[indx[j]] ) )
        {
          j = j + 1;
        }
      }

      if ( xval < x[indx[j-1]] ||
         ( xval == x[indx[j-1]] && yval < y[indx[j-1]] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
//
//  Take care of the base.
//
  for ( i = 0; i < n; i++ )
  {
    indx[i] = indx[i] + base;
  }

  return indx;
}
//****************************************************************************80

void r4vec2_sorted_unique ( int n, float a1[], float a2[], int *unique_num )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_SORTED_UNIQUE keeps the unique elements in an R4VEC2.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//    Item I is stored as the pair A1(I), A2(I).
//
//    The items must have been sorted, or at least it must be the
//    case that equal items are stored in adjacent vector locations.
//
//    If the items were not sorted, then this routine will only
//    replace a string of equal values by a single representative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items.
//
//    Input/output, float A1[N], A2[N].
//    On input, the array of N items.
//    On output, an array of UNIQUE_NUM unique items.
//
//    Output, int *UNIQUE_NUM, the number of unique items.
//
{
  int itest;

  *unique_num = 0;

  if ( n <= 0 )
  {
    return;
  }

  *unique_num = 1;

  for ( itest = 1; itest < n; itest++ )
  {
    if ( a1[itest] != a1[*unique_num-1] ||
         a2[itest] != a2[*unique_num-1] )
    {
      a1[*unique_num] = a1[itest];
      a2[*unique_num] = a2[itest];
      *unique_num = *unique_num + 1;
    }
  }

  return;
}
//****************************************************************************80

void r4vec2_sorted_unique_index ( int n, float a1[], float a2[],
  int *unique_num, int indx[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_SORTED_UNIQUE_INDEX indexes unique elements in a sorted R4VEC2.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//    Item I is stored as the pair A1(I), A2(I).
//
//    The items must have been sorted, or at least it should be the
//    case that equal items are stored in adjacent vector locations.
//
//    If the items are not sorted, then this routine will only
//    replace a string of equal values by a single representative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items.
//
//    Input/output, float A1[N], A2[N].
//    On input, the array of N items.
//    On output, an array of unique items.
//
//    Output, int *UNIQUE_NUM, the number of unique items.
//
//    Output, int INDX[N], contains in entries 1 through UNIQUE_NUM an index
//    array of the unique items.  To build new arrays with no repeated elements:
//      B1(*) = A1(INDX(*))
//
{
  int itest;

  if ( n <= 0 )
  {
    *unique_num = 0;
    return;
  }
  i4vec_zero ( n, indx );

  *unique_num = 1;
  indx[0] = 1;

  for ( itest = 2; itest <= n; itest++ )
  {
    if ( a1[itest-2] != a1[itest-1] || a2[itest-2] != a2[itest-1] )
    {
      *unique_num = *unique_num + 1;
      indx[*unique_num-1] = itest;
    }
  }

  return;
}
//****************************************************************************80

int r4vec2_sum_max_index ( int n, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC2_SUM_MAX_INDEX returns the index of the maximum sum of two R4VEC's.
//
//  Discussion:
//
//    An R4VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, float A[N], B[N], two arrays whose sum
//    is to be examined.
//
//    Output, int R4VEC2_SUM_MAX_INDEX, the index of the largest entry in A+B.
//
{
  int i;
  float sum_max;
  int sum_max_index;

  if ( n <= 0 )
  {
    sum_max_index = -1;
  }
  else
  {
    sum_max_index = 1;
    sum_max = a[0] + b[0];

    for ( i = 2; i <= n; i++ )
    {
      if ( sum_max < a[i-1] + b[i-1] )
      {
        sum_max = a[i-1] + b[i-1];
        sum_max_index = i;
      }
    }
  }
  return sum_max_index;
}
//****************************************************************************80

void r4vec3_print ( int n, float a1[], float a2[], float a3[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC3_PRINT prints a triple of real vectors.
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
//    Input, int N, the number of components of the vector.
//
//    Input, float A1[N], float A2[N], float A3[N], the vectors
//    to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(4)  << i     << ": "
         << setw(10) << a1[i] << "  "
         << setw(10) << a2[i] << "  "
         << setw(10) << a3[i] << "\n";
  }

  return;
}
//****************************************************************************80

float *roots_to_r4poly ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    ROOTS_TO_R4POLY converts polynomial roots to polynomial coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of roots specified.
//
//    Input, float X[N], the roots.
//
//    Output, float ROOTS_TO_R4POLY[N+1], the coefficients of the polynomial.
//
{
  float *c;
  int i;
  int j;

  c = r4vec_zero_new ( n + 1 );
//
//  Initialize C to (0, 0, ..., 0, 1).
//  Essentially, we are setting up a divided difference table.
//
  c[n] = 1.0;
//
//  Convert to standard polynomial form by shifting the abscissas
//  of the divided difference table to 0.
//
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n+1-j; i++ )
    {
      c[n-i] = c[n-i] - x[n+1-i-j] * c[n-i+1];
    }
  }
  return c;
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
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
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

  for ( ; ; )
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
