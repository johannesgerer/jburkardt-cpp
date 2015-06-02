# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "monomial.hpp"

//****************************************************************************80

int i4_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
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
//    06 January 2013
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
//    Output, int I4_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = k;
  if ( n - k < mn )
  {
    mn = n - k;
  }

  if ( mn < 0 )
  {
    value = 0;
  }
  else if ( mn == 0 )
  {
    value = 1;
  }
  else
  {
    mx = k;
    if ( mx < n - k )
    {
      mx = n - k;
    }
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
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

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
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
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

int *i4vec_uniform_ab_new ( int n, int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIFORM_AB_NEW returns a scaled pseudorandom I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    The pseudorandom numbers should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 May 2012
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
//    Input, int N, the dimension of the vector.
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int IVEC_UNIFORM_AB_NEW[N], a vector of random values 
//    between A and B.
//
{
  int c;
  int i;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;
  int *x;
  
  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_UNIFORM_AB_NEW - Fatal error!\n";
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

  x = new int[n];

  for ( i = 0; i < n; i++ )
  {
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

    x[i] = value;
  }

  return x;
}
//****************************************************************************80

int mono_between_enum ( int m, int n1, int n2 )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_ENUM enumerates monomials in M dimensions of degrees in a range.
//
//  Discussion:
//
//    For D = 3, we have the following table:
//
//     N2 0  1  2  3  4  5  6   7   8
//   N1 +----------------------------
//    0 | 1  4 10 20 35 56 84 120 165
//    1 | 0  3  9 19 34 55 83 119 164
//    2 | 0  0  6 16 31 52 80 116 161
//    3 | 0  0  0 10 25 46 74 110 155
//    4 | 0  0  0  0 15 36 64 100 145
//    5 | 0  0  0  0  0 21 49  85 130
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N1, N2, the minimum and maximum degrees.
//    0 <= N1 <= N2.
//
//    Output, int MONO_BETWEEN_ENUM, the number of monomials 
//    in D variables, of total degree between N1 and N2 inclusive.
//
{
  int n0;
  int n1_copy;
  int value;

  n1_copy = i4_max ( n1, 0 );

  if ( n2 < n1_copy )
  {
    value = 0;
    return value;
  }

  if ( n1_copy == 0 )
  {
    value = i4_choose ( n2 + m, n2 );
  }
  else if ( n1_copy == n2 )
  {
    value = i4_choose ( n2 + m - 1, n2 );
  }
  else
  {
    n0 = n1_copy - 1;
    value = i4_choose ( n2 + m, n2 ) - i4_choose ( n0 + m, n0 );
  }

  return value;
}
//****************************************************************************80

void mono_between_next_grevlex ( int m, int n1, int n2, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_NEXT_GREVLEX: grevlex next monomial, degree between N1 and N2.
//
//  Discussion:
//
//    We consider all monomials in an M dimensional space, with total
//    degree N between N1 and N2, inclusive.
//
//    For example:
//
//    M = 3
//    N1 = 2
//    N2 = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     2        2
//    2 |  0     1     1        2
//    3 |  1     0     1        2
//    4 |  0     2     0        2
//    5 |  1     1     0        2
//    6 |  2     0     0        2
//      |
//    7 |  0     0     3        3
//    8 |  0     1     2        3
//    9 |  1     0     2        3
//   10 |  0     2     1        3
//   11 |  1     1     1        3
//   12 |  2     0     1        3
//   13 |  0     3     0        3
//   14 |  1     2     0        3
//   15 |  2     1     0        3
//   16 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N1, N2, the minimum and maximum degrees.
//    0 <= N1 <= N2.
//
//    Input, int X[M], the current monomial.
//    To start the sequence, set X = [ N1, 0, 0, ... ].
//
//    Output, int X[M], the next monomial.
//    The last value in the sequence is X = [ 0, 0, ..., 0, N2 ].
//
{
  int i;
  int j;
  int t;

  if ( n1 < 0 )
  {
    cout << "\n";
    cout << "MONO_BETWEEN_NEXT_GREVLEX - Fatal error!\n";
    cout << "  N1 < 0.\n";
    exit ( 1 );
  }

  if ( n2 < n1 ) 
  {
    cout << "\n";
    cout << "MONO_BETWEEN_NEXT_GREVLEX - Fatal error!\n";
    cout << "  N2 < N1.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( m, x ) < n1 )
  {
    cout << "\n";
    cout << "MONO_BETWEEN_NEXT_GREVLEX - Fatal error!\n";
    cout << "  Input X sums to less than N1.\n";
    exit ( 1 );
  }

  if ( n2 < i4vec_sum ( m, x ) )
  {
    cout << "\n";
    cout << "MONO_BETWEEN_NEXT_GREVLEX - Fatal error!\n";
    cout << "  Input X sums to more than N2.\n";
    exit ( 1 );
  }

  if ( n2 == 0 )
  {
    return;
  }

  if ( x[0] == n2 )
  {
    x[0] = 0;
    x[m-1] = n1;
  }
  else
  {
    mono_next_grevlex ( m, x );
  }

  return;
}
//****************************************************************************80

void mono_between_next_grlex ( int m, int n1, int n2, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_NEXT_GRLEX: grlex next monomial, degree between N1 and N2.
//
//  Discussion:
//
//    We consider all monomials in an M dimensional space, with total
//    degree N between N1 and N2, inclusive.
//
//    For example:
//
//    M = 3
//    N1 = 2
//    N2 = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     2        2
//    2 |  0     1     1        2
//    3 |  0     2     0        2
//    4 |  1     0     1        2
//    5 |  1     1     0        2
//    6 |  2     0     0        2
//      |
//    7 |  0     0     3        3
//    8 |  0     1     2        3
//    9 |  0     2     1        3
//   10 |  0     3     0        3
//   11 |  1     0     2        3
//   12 |  1     1     1        3
//   13 |  1     2     0        3
//   14 |  2     0     1        3
//   15 |  2     1     0        3
//   16 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N1, N2, the minimum and maximum degrees.
//    0 <= N1 <= N2.
//
//    Input/output, int X[M], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, N1 ].
//    The last value in the sequence is X = [ N2, 0, ..., 0, 0 ].
//
{
  int i;
  int im1;
  int j;
  int t;

  if ( n1 < 0 )
  {
    cerr << "\n";
    cerr << "MONO_BETWEEN_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N1 < 0.\n";
    exit ( 1 );
  }

  if ( n2 < n1 )
  {
    cerr << "\n";
    cerr << "MONO_BETWEEN_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N2 < N1.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( m, x ) < n1 )
  {
    cerr << "\n";
    cerr << "MONO_BETWEEN_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X sums to less than N1.\n";
    exit ( 1 );
  }

  if ( n2 < i4vec_sum ( m, x ) )
  {
    cerr << "\n";
    cerr << "MONO_BETWEEN_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X sums to more than N2.\n";
    exit ( 1 );
  }

  if ( n2 == 0 )
  {
    return;
  }

  if ( x[0] == n2 )
  {
    x[0] = 0;
    x[m-1] = n1;
  }
  else
  {
    mono_next_grlex ( m, x );
  }

  return;
}
//****************************************************************************80

int *mono_between_random ( int m, int n1, int n2, int &seed, int &rank )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_RANDOM: random monomial with total degree between N1 and N2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N1, N2, the minimum and maximum degrees.
//    0 <= N1 <= N2.
//
//    Input/output, int *SEED, the random number seed.
//
//    Output int *RANK, the rank of the monomial.
//
//    Output int MONO_BETWEEN_RANDOM[M], the random monomial.
//
{
  int n1_copy;
  int rank_max;
  int rank_min;
  int *x;

  n1_copy = i4_max ( n1, 0 );
  rank_min = mono_upto_enum ( m, n1_copy - 1 ) + 1;
  rank_max = mono_upto_enum ( m, n2 );
  rank = i4_uniform_ab ( rank_min, rank_max, seed );
  x = mono_unrank_grlex ( m, rank );

  return x;
}
//****************************************************************************80

void mono_next_grevlex ( int m, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_NEXT_GREVLEX: grevlex next monomial.
//
//  Discussion:
//
//    Example:
//
//    M = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     0        0
//      |
//    2 |  0     0     1        1
//    3 |  0     1     0        1
//    4 |  1     0     0        1
//      |
//    5 |  0     0     2        2
//    6 |  0     1     1        2
//    7 |  1     0     1        2
//    8 |  0     2     0        2
//    9 |  1     1     0        2
//   10 |  2     0     0        2
//      |
//   11 |  0     0     3        3
//   12 |  0     1     2        3
//   13 |  1     0     2        3
//   14 |  0     2     1        3
//   15 |  1     1     1        3
//   16 |  2     0     1        3
//   17 |  0     3     0        3
//   18 |  1     2     0        3
//   19 |  2     1     0        3
//   20 |  3     0     0        3
//
//    Thanks to Stefan Klus for pointing out a discrepancy in a previous
//    version of this code, 05 February 2015.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int X[M], the current monomial.
//    The first item is X = [ 0, 0, ..., 0, 0 ].
//
{
  int i;
  int j;
  int t;

  if ( i4vec_sum ( m, x ) < 0 )
  {
    cout << "\n";
    cout << "MONO_UPTO_NEXT_GREVLEX - Fatal error!\n";
    cout << "  Input X sums to less than 0.\n";
    exit ( 1 );
  }
//
//  Seeking the first index 1 < I for which 0 < X(I).
//
  j = 0;

  for ( i = 1; i < m; i++ )
  {
    if ( 0 < x[i] )
    {
      j = i;
      break;
    }
  }

  if ( j == 0 )
  {
    t = x[0];
    x[0] = 0;
    x[m-1] = t + 1;
  }
  else if ( j < m - 1 )
  {
    x[j] = x[j] - 1;
    t = x[0] + 1;
    x[0] = 0;
    x[j-1] = x[j-1] + t;
  }
  else if ( j == m - 1 )
  {
    t = x[0];
    x[0] = 0;
    x[j-1] = t + 1;
    x[j] = x[j] - 1;
  }

  return;
}
//****************************************************************************80

void mono_next_grlex ( int m, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_NEXT_GRLEX returns the next monomial in grlex order.
//
//  Discussion:
//
//    Example:
//
//    M = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     0        0
//      |
//    2 |  0     0     1        1
//    3 |  0     1     0        1
//    4 |  1     0     0        1
//      |
//    5 |  0     0     2        2
//    6 |  0     1     1        2
//    7 |  0     2     0        2
//    8 |  1     0     1        2
//    9 |  1     1     0        2
//   10 |  2     0     0        2
//      |
//   11 |  0     0     3        3
//   12 |  0     1     2        3
//   13 |  0     2     1        3
//   14 |  0     3     0        3
//   15 |  1     0     2        3
//   16 |  1     1     1        3
//   17 |  1     2     0        3
//   18 |  2     0     1        3
//   19 |  2     1     0        3
//   20 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int X[M], the current monomial.
//    The first element is X = [ 0, 0, ..., 0, 0 ].
//
{
  int i;
  int im1;
  int j;
  int t;
//
//  Ensure that 1 <= M.
//
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= X(I).
//
  for ( i = 0; i < m; i++ )
  {
    if ( x[i] < 0 )
    {
      cerr << "\n";
      cerr << "MONO_NEXT_GRLEX - Fatal error!\n";
      cerr << "  X[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  Find I, the index of the rightmost nonzero entry of X.
//
  i = 0;
  for ( j = m; 1 <= j; j-- )
  {
    if ( 0 < x[j-1] )
    {
      i = j;
      break;
    }
  }
//
//  set T = X(I)
//  set X(I) to zero,
//  increase X(I-1) by 1,
//  increment X(D) by T-1.
//
  if ( i == 0 )
  {
    x[m-1] = 1;
    return;
  }
  else if ( i == 1 )
  {
    t = x[0] + 1;
    im1 = m;
  }
  else if ( 1 < i )
  {
    t = x[i-1];
    im1 = i - 1;
  }

  x[i-1] = 0;
  x[im1-1] = x[im1-1] + 1;
  x[m-1] = x[m-1] + t - 1;

  return;
}
//****************************************************************************80

void mono_print ( int m, int f[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_PRINT prints a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int F[M], the exponents.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << title;
  cout << "x^(";
  for ( i = 0; i < m; i++ )
  {
    cout << f[i];
    if ( i < m - 1 )
    {
      cout << ",";
    }
    else
    {
      cout << ").\n";
    }
  }

  return;
}
//****************************************************************************80

int mono_rank_grlex ( int m, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_RANK_GRLEX computes the graded lexicographic rank of a monomial.
//
//  Discussion:
//
//    The graded lexicographic ordering is used, over all monomials of
//    dimension M, for degree NM = 0, 1, 2, ...
//
//    For example, if M = 3, the ranking begins:
//
//    Rank  Sum    1  2  3
//    ----  ---   -- -- --
//       1    0    0  0  0
//
//       2    1    0  0  1
//       3    1    0  1  0
//       4    1    1  0  1
//
//       5    2    0  0  2
//       6    2    0  1  1
//       7    2    0  2  0
//       8    2    1  0  1
//       9    2    1  1  0
//      10    2    2  0  0
//
//      11    3    0  0  3
//      12    3    0  1  2
//      13    3    0  2  1
//      14    3    0  3  0
//      15    3    1  0  2
//      16    3    1  1  1
//      17    3    1  2  0
//      18    3    2  0  1
//      19    3    2  1  0
//      20    3    3  0  0
//
//      21    4    0  0  4
//      ..   ..   .. .. ..
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    1 <= D.
//
//    Input, int XC[M], the monomial.
//    For each 1 <= I <= M, we have 0 <= XC(I).
//
//    Output, int MONO_RANK_GRLEX, the rank.
//
{
  int i;
  int j;
  int ks;
  int n;
  int nm;
  int ns;
  int rank;
  int tim1;
  int *xs;
//
//  Ensure that 1 <= M.
//
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_RANK_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= X(I).
//
  for ( i = 0; i < m; i++ )
  {
    if ( x[i] < 0 )
    {
      cerr << "\n";
      cerr << "MONO_RANK_GRLEX - Fatal error!\n";
      cerr << "  X[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  NM = sum ( X )
//
  nm = i4vec_sum ( m, x );
//
//  Convert to KSUBSET format.
//
  ns = nm + m - 1;
  ks = m - 1;
  xs = new int[ks];
  xs[0] = x[0] + 1;
  for ( i = 2; i < m; i++ )
  {
    xs[i-1] = xs[i-2] + x[i-1] + 1;
  }
//
//  Compute the rank.
//
  rank = 1;

  for ( i = 1; i <= ks; i++ )
  {
    if ( i == 1 )
    {
      tim1 = 0;
    }
    else
    {
      tim1 = xs[i-2];
    }

    if ( tim1 + 1 <= xs[i-1] - 1 )
    {
      for ( j = tim1 + 1; j <= xs[i-1] - 1; j++ )
      {
        rank = rank + i4_choose ( ns - j, ks - i );
      }
    }
  }

  for ( n = 0; n < nm; n++ )
  {
    rank = rank + i4_choose ( n + m - 1, n );
  }

  delete [] xs;

  return rank;
}
//****************************************************************************80

int mono_total_enum ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_ENUM enumerates monomials in M dimensions of degree equal to N.
//
//  Discussion:
//
//    For M = 3, we have the following values:
//
//    N  VALUE
//
//    0    1
//    1    3
//    2    6
//    3   10
//    4   15
//    5   21
//
//    In particular, VALUE(3,3) = 10 because we have the 10 monomials:
//
//      x^3, x^2y, x^2z, xy^2, xyz, xz^3, y^3, y^2z, yz^2, z^3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the maximum degree.
//
//    Output, int MONO_TOTAL_ENUM, the number of monomials in D variables,
//    of total degree N.
//
{
  int value;

  value = i4_choose ( n + m - 1, n );

  return value;
}
//****************************************************************************80

void mono_total_next_grevlex ( int m, int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_NEXT_GREVLEX: grevlex next monomial with total degree equal to N.
//
//  Discussion:
//
//    We consider all monomials in a M dimensional space, with total degree N.
//
//    For example:
//
//    M = 3
//    N = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     3        3
//    2 |  0     1     2        3
//    3 |  1     0     2        3
//    4 |  0     2     1        3
//    5 |  1     1     1        3
//    6 |  2     0     1        3
//    7 |  0     3     0        3
//    8 |  1     2     0        3
//    9 |  2     1     0        3
//   10 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the degree.
//    0 <= N.
//
//    Input, int X[M], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, N ].
//
//    Output, int X[M], the next monomial.
//    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
//
{
  int i;
  int j;
  int t;

  if ( n < 0 )
  {
    cout << "\n";
    cout << "MONO_TOTAL_NEXT_GREVLEX - Fatal error!\n";
    cout << "  N < 0.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( m, x ) != n )
  {
    cout << "\n";
    cout << "MONO_TOTAL_NEXT_GREVLEX - Fatal error!\n";
    cout << "  Input X does not sum to N.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    return;
  }
  if ( x[0] == n )
  {
    x[0] = 0;
    x[m-1] = n;
  }
  else
  {
    mono_next_grevlex ( m, x );
  }

  return;
}
//****************************************************************************80

void mono_total_next_grlex ( int m, int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_NEXT_GRLEX: grlex next monomial with total degree equal to N.
//
//  Discussion:
//
//    We consider all monomials in a M dimensional space, with total degree N.
//
//    For example:
//
//    M = 3
//    N = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     3        3
//    2 |  0     1     2        3
//    3 |  0     2     1        3
//    4 |  0     3     0        3
//    5 |  1     0     2        3
//    6 |  1     1     1        3
//    7 |  1     2     0        3
//    8 |  2     0     1        3
//    9 |  2     1     0        3
//   10 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the degree.
//    0 <= N.
//
//    Input/output, int X[M], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, N ].
//    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
//
{
  int i;
  int im1;
  int j;
  int t;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "MONO_TOTAL_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( m, x ) != n )
  {
    cerr << "\n";
    cerr << "MONO_TOTAL_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X does not sum to N.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    return;
  }

  if ( x[0] == n )
  {
    x[0] = 0;
    x[m-1] = n;
  }
  else
  {
    mono_next_grlex ( m, x );
  }

  return;
}
//****************************************************************************80

int *mono_total_random ( int m, int n, int &seed, int &rank )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_RANDOM: random monomial with total degree equal to N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the degree.
//    0 <= N.
//
//    Input/output, int &SEED, the random number seed.
//
//    Output, int &RANK, the rank of the monomial.
//
//    Output, int MONO_TOTAL_RANDOM[M], the random monomial.
//
{
  int rank_max;
  int rank_min;
  int *x;

  rank_min = mono_upto_enum ( m, n - 1 ) + 1;
  rank_max = mono_upto_enum ( m, n );
  rank = i4_uniform_ab ( rank_min, rank_max, seed );
  x = mono_unrank_grlex ( m, rank );

  return x;
}
//****************************************************************************80

int *mono_unrank_grlex ( int m, int rank )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UNRANK_GRLEX computes the composition of given grlex rank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    1 <= M.
//
//    Input, int RANK, the rank.
//    1 <= RANK.
//
//    Output, int MONO_UNRANK_GRLEX[M], the monomial X of the given rank.
//    For each I, 0 <= XC[I] <= NM, and 
//    sum ( 1 <= I <= M ) XC[I] = NM.
//
{
  int i;
  int j;
  int ks;
  int nksub;
  int nm;
  int ns;
  int r;
  int rank1;
  int rank2;
  int *x;
  int *xs;
//
//  Ensure that 1 <= M.
//
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 1 <= RANK.
//
  if ( rank < 1 )
  {
    cerr << "\n";
    cerr << "MONO_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  RANK < 1\n";
    exit ( 1 );
  }
//
//  Special case M == 1.
//
  if ( m == 1 )
  {
    x = new int[m];
    x[0] = rank - 1;
    return x;
  }
//
//  Determine the appropriate value of NM.
//  Do this by adding up the number of compositions of sum 0, 1, 2, 
//  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
//  gives you the rank of the composition within the set of compositions
//  of sum NM.  And that's the number you need in order to do the
//  unranking.
//
  rank1 = 1;
  nm = -1;
  for ( ; ; )
  {
    nm = nm + 1;
    r = i4_choose ( nm + m - 1, nm );
    if ( rank < rank1 + r )
    {
      break;
    }
    rank1 = rank1 + r;
  }

  rank2 = rank - rank1;
//
//  Convert to KSUBSET format.
//  Apology: an unranking algorithm was available for KSUBSETS,
//  but not immediately for compositions.  One day we will come back
//  and simplify all this.
//
  ks = m - 1;
  ns = nm + m - 1;
  xs = new int[ks];

  nksub = i4_choose ( ns, ks );

  j = 1;

  for ( i = 1; i <= ks; i++ )
  {
    r = i4_choose ( ns - j, ks - i );

    while ( r <= rank2 && 0 < r )
    {
      rank2 = rank2 - r;
      j = j + 1;
      r = i4_choose ( ns - j, ks - i );
    }
    xs[i-1] = j;
    j = j + 1;
  }
//
//  Convert from KSUBSET format to COMP format.
//
  x = new int[m];
  x[0] = xs[0] - 1;
  for ( i = 2; i < m; i++ )
  {
    x[i-1] = xs[i-1] - xs[i-2] - 1;
  }
  x[m-1] = ns - xs[ks-1];

  delete [] xs;

  return x;
}
//****************************************************************************80

int mono_upto_enum ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_ENUM enumerates monomials in M dimensions of degree up to N.
//
//  Discussion:
//
//    For M = 2, we have the following values:
//
//    N  VALUE
//
//    0    1
//    1    3
//    2    6
//    3   10
//    4   15
//    5   21
//
//    In particular, VALUE(2,3) = 10 because we have the 10 monomials:
//
//      1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the maximum degree.
//
//    Output, int MONO_UPTO_ENUM, the number of monomials in
//    M variables, of total degree N or less.
//
{
  int value;

  value = i4_choose ( n + m, n );

  return value;
}
//****************************************************************************80

void mono_upto_next_grevlex ( int m, int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_NEXT_GREVLEX: grevlex next monomial with total degree up to N.
//
//  Discussion:
//
//    We consider all monomials in a M dimensional space, with total
//    degree up to N.
//
//    For example:
//
//    M = 3
//    N = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     0        0
//      |
//    2 |  0     0     1        1
//    3 |  0     1     0        1
//    4 |  1     0     0        1
//      |
//    5 |  0     0     2        2
//    6 |  0     1     1        2
//    7 |  1     0     1        2
//    8 |  0     2     0        2
//    9 |  1     1     0        2
//   10 |  2     0     0        2
//      |
//   11 |  0     0     3        3
//   12 |  0     1     2        3
//   13 |  1     0     2        3
//   14 |  0     2     1        3
//   15 |  1     1     1        3
//   16 |  2     0     1        3
//   17 |  0     3     0        3
//   18 |  1     2     0        3
//   19 |  2     1     0        3
//   20 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the maximum degree.
//    0 <= N.
//
//    Input, int X[M], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, 0 ].
//
//    Output, int X[M], the next monomial.
//    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
//
{
  int i;
  int j;
  int t;

  if ( n < 0 )
  {
    cout << "\n";
    cout << "MONO_UPTO_NEXT_GREVLEX - Fatal error!\n";
    cout << "  N < 0.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( m, x ) < 0 )
  {
    cout << "\n";
    cout << "MONO_UPTO_NEXT_GREVLEX - Fatal error!\n";
    cout << "  Input X sums to less than 0.\n";
    exit ( 1 );
  }

  if ( n < i4vec_sum ( m, x ) )
  {
    cout << "\n";
    cout << "MONO_UPTO_NEXT_GREVLEX - Fatal error!\n";
    cout << "  Input X sums to more than N.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    return;
  }

  if ( x[0] == n )
  {
    x[0] = 0;
  }
  else
  {
    mono_next_grevlex ( m, x );
  }

  return;
}
//****************************************************************************80

void mono_upto_next_grlex ( int m, int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_NEXT_GRLEX: grlex next monomial with total degree up to N.
//
//  Discussion:
//
//    We consider all monomials in a M dimensional space, with total
//    degree up to N.
//
//    For example:
//
//    M = 3
//    N = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     0        0
//      |
//    2 |  0     0     1        1
//    3 |  0     1     0        1
//    4 |  1     0     0        1
//      |
//    5 |  0     0     2        2
//    6 |  0     1     1        2
//    7 |  0     2     0        2
//    8 |  1     0     1        2
//    9 |  1     1     0        2
//   10 |  2     0     0        2
//      |
//   11 |  0     0     3        3
//   12 |  0     1     2        3
//   13 |  0     2     1        3
//   14 |  0     3     0        3
//   15 |  1     0     2        3
//   16 |  1     1     1        3
//   17 |  1     2     0        3
//   18 |  2     0     1        3
//   19 |  2     1     0        3
//   20 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the maximum degree.
//    0 <= N.
//
//    Input/output, int X[M], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, 0 ].
//    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
//
{
  int i;
  int im1;
  int j;
  int t;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "MONO_UPTO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( m, x ) < 0 )
  {
    cerr << "\n";
    cerr << "MONO_UPTO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X sums to less than 0.\n";
    exit ( 1 );
  }

  if ( n < i4vec_sum ( m, x ) )
  {
    cerr << "\n";
    cerr << "MONO_UPTO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X sums to more than N.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    return;
  }
  if ( x[0] == n )
  {
    x[0] = 0;
    x[m-1] = n;
  }
  else
  {
    mono_next_grlex ( m, x );
  }

  return;
}
//****************************************************************************80

int *mono_upto_random ( int m, int n, int &seed, int &rank )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_RANDOM: random monomial with total degree less than or equal to N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the degree.
//    0 <= N.
//
//    Input/output, int &SEED, the random number seed.
//
//    Output, int &RANK, the rank of the monomial.
//
//    Output, int MONO_UPTO_RANDOM[M], the random monomial.
//
{
  int rank_max;
  int rank_min;
  int *x;

  rank_min = 1;
  rank_max = mono_upto_enum ( m, n );
  rank = i4_uniform_ab ( rank_min, rank_max, seed );
  x = mono_unrank_grlex ( m, rank );

  return x;
}
//****************************************************************************80

double *mono_value ( int m, int n, int f[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_VALUE evaluates a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2013
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
//    Input, int F[M], the exponents of the monomial.
//
//    Input, double X[M*N], the coordinates of the evaluation points.
//
//    Output, double MONO_VALUE[N], the value of the monomial at X.
//
{
  int i;
  int j;
  double *v;

  v = new double[n];

  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
    for ( i = 0; i < m; i++ )
    {
      v[j] = v[j] * pow ( x[i+j*m], f[i] );
    }
  }

  return v;
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
