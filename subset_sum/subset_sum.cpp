# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <string>

using namespace std;

# include "subset_sum.hpp"

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

int *i4_to_digits_binary ( int i, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_DIGITS_BINARY produces the binary digits of an I4.
//
//  Example:
//
//     I    N     C               Binary
//    --  ---   ---         ------------
//     0    1   0                      0
//     0    2   0, 0                  00
//     1    3   1, 0, 0              100
//     2    3   0, 1, 0              010
//     3    3   1, 1, 0              011
//     4    3   0, 0, 1              100
//     8    3   0, 0, 0           (1)000
//     8    5   0, 0, 0, 1, 0      01000
//    -8    5   0, 0, 0, 1, 0  (-) 01000
//
//     0    3   0, 0, 0
//     1    3   1, 0, 0
//     2    3   0, 1, 0
//     3    3   1, 1, 0
//     4    3   0, 0, 1
//     5    3   1, 0, 1
//     6    3   0, 1, 1
//     7    3   1, 1, 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer to be analyzed.
//
//    Input, int N, the number of digits to determine.
//
//    Output, int I4_TO_DIGITS_BINARY[N], the first N binary digits of I.
//    Entry 0 is the units digit.
//
{
  int *c;
  int j;

  c = new int[n];

  i = abs ( i );

  for ( j = 0; j < n; j++ )
  {
    c[j] = i % 2;
    i = ( i - c[j] ) / 2;
  }

  return c;
}
//****************************************************************************80

int *i4vec_copy_new ( int n, int a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY_NEW copies an I4VEC to a "new" I4VEC.
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
//    04 July 2008
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
//    Output, int I4VEC_COPY_NEW[N], the copy of A1.
//
{
  int *a2;
  int i;

  a2 = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

int i4vec_dot_product ( int n, int x[], int y[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_DOT_PRODUCT computes the dot product of two I4VEC's.
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
//    19 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the array.
//
//    Input, int X[N], Y[N], the arrays.
//
//    Output, int I4VEC_DOT_PRODUCT, the dot product of X and Y.
//
{
  int i;
  int value;

  value = 0;
  for ( i = 0; i < n; i++ )
  {
    value = value + x[i] * y[i];
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

int subset_sum_count ( int n, int w[], int t, int ind_min, int ind_max )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_SUM_COUNT counts solutions to the subset sum problem in a range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the subset.
//
//    Input, int W[N], a set of weights.  The length of this
//    array must be no more than 31.
//
//    Input, int T, the target value.
//
//    Input, int IND_MIN, IND_MAX, the lower and upper
//    limits to be searched.  0 <= IND_MIN <= IND_MAX <= (2^N)-1.
//
//    Output, int SUBSET_SUM_COUNT, the number of distinct
//    solutions of the subset sum problem found within the given range.
//
{
  int *c;
  int ind;
  int ind_max2;
  int ind_min2;
  int solution_num;
//
//  Check the data.
//
  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "SUBSET_SUM_COUNT - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( 31 < n )
  {
    cerr << "\n";
    cerr << "SUBSET_SUM_COUNT - Fatal error!\n";
    cerr << "  31 < N.\n";
    exit ( 1 );
  }

  ind_min2 = i4_max ( ind_min, 0 );
  ind_max2 = i4_min ( ind_max, i4_power ( 2, n ) - 1 );
//
//  Run through the range.
//
  cout << "\n";
  cout << "  Searching from IND_MIN = " << ind_min2;
  cout << "  through IND_MAX = " << ind_max2 << "\n";

  solution_num = 0;

  for ( ind = ind_min2; ind <= ind_max2; ind++ )
  {
//
//  Convert INDEX into vector of indices in W.
//
    c = i4_to_digits_binary ( ind, n );
//
//  If the sum of those weights matches the target, return combination.
//
    if ( i4vec_dot_product ( n, c, w ) == t )
    {
      solution_num = solution_num + 1;
    }

    delete [] c;
  }

  return solution_num;
}
//****************************************************************************80

int *subset_sum_find ( int n, int w[], int t, int ind_min, int ind_max, 
  int &ind )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_SUM seeks a subset of a set that has a given sum.
//
//  Discussion:
//
//    This function tries to compute a target value as the sum of
//    a selected subset of a given set of weights.
//
//    This function works by brute force, that is, it tries every
//    possible subset to see if it sums to the desired value.
//
//    Given N weights, every possible selection can be described by
//    one of the N-digit binary numbers from 0 to 2^N-1.
//
//    This function includes a range, which allows the user to
//    control which subsets are to be checked.  Thus, if there are
//    N weights, specifying a range of [ 0, 2^N-1] indicates that
//    all subsets should be checked.  On the other hand, this full
//    range could be broken down into smaller subranges, each of
//    which could be checked independently.
//
//    It is possible that, in the given range, there may be multiple
//    solutions of the problem.  This function will only return
//    one such solution, if found.  However, the function may be called
//    again, with an appropriate restriction of the range, to continue
//    the search for other solutions.
//
//  Example:
//
//    w = [ 1, 2, 4, 8, 16, 32 ];
//    t = 22;
//    r = [ 0, 2^6 - 1 ];
//
//    call subset_sum ( w, t, r, c, ind )
//
//    c = [ 2, 3, 5 ]
//    index = 22
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the subset.
//
//    Input, int W[N], a set of weights.  The length of this
//    array must be no more than 31.
//
//    Input, int T, the target value.
//
//    Input, int IND_MIN, IND_MAX, the lower and upper
//    limits to be searched.  0 <= IND_MIN <= IND_MAX <= (2^N)-1.
//
//    Output, int &IND, the index of the solution.
//    If IND is -1, no solution was found in the range.
//
//    Output, int SUBSET_SUM_FIND[N], indicates the solution, assuming
//    that IND is not -1.  In that case, the sum T is made by selecting
//    those weights W(I) for which C(I) is 1.  In fact,
//    T = sum ( 1 <= I <= N ) C(I) * W(I).
//
{
  int *c;
  int ind_max2;
  int ind_min2;
//
//  Check the data.
//
  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "SUBSET_SUM_COUNT - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( 31 < n )
  {
    cerr << "\n";
    cerr << "SUBSET_SUM_COUNT - Fatal error!\n";
    cerr << "  31 < N.\n";
    exit ( 1 );
  }

  ind_min2 = i4_max ( ind_min, 0 );
  ind_max2 = i4_min ( ind_max, i4_power ( 2, n ) - 1 );
//
//  Run through the range.
//
  cout << "\n";
  cout << "  Searching from IND_MIN = " << ind_min2;
  cout << "  through IND_MAX = " << ind_max2 << "\n";

  for ( ind = ind_min2; ind <= ind_max2; ind++ )
  {
//
//  Convert INDEX into vector of indices in W.
//
    c = i4_to_digits_binary ( ind, n );
//
//  If the sum of those weights matches the target, return combination.
//
    if ( i4vec_dot_product ( n, c, w ) == t )
    {
      return c;
    }

    delete [] c;
  }

  ind = - 1;

  return NULL;
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
