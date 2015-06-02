# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "subset_sum_serial.hpp"

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
/******************************************************************************/

int *i4_to_digits_binary ( int i, int n )

/******************************************************************************/
/*
  Purpose:

    I4_TO_DIGITS_BINARY produces the binary digits of an I4.

  Example:

     I    N     C               Binary
    --  ---   ---         ------------
     0    1   0                      0
     0    2   0, 0                  00
     1    3   1, 0, 0              100
     2    3   0, 1, 0              010
     3    3   1, 1, 0              011
     4    3   0, 0, 1              100
     8    3   0, 0, 0           (1)000
     8    5   0, 0, 0, 1, 0      01000
    -8    5   0, 0, 0, 1, 0  (-) 01000

     0    3   0, 0, 0
     1    3   1, 0, 0
     2    3   0, 1, 0
     3    3   1, 1, 0
     4    3   0, 0, 1
     5    3   1, 0, 1
     6    3   0, 1, 1
     7    3   1, 1, 1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int I, the integer to be analyzed.

    Input, int N, the number of digits to determine.

    Output, int I4_TO_DIGITS_BINARY[N], the first N binary digits of I.
    Entry 0 is the units digit.
*/
{
  int *c;
  int j;

  c = ( int * ) malloc ( n * sizeof ( int ) );

  i = abs ( i );

  for ( j = 0; j < n; j++ )
  {
    c[j] = i % 2;
    i = ( i - c[j] ) / 2;
  }

  return c;
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

int *subset_sum_serial ( int n, int weight[], int target )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_SUM_SERIAL seeks a subset of a set that has a given sum.
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
//    It is possible that there may be multiple solutions of the problem.  
//    This function will only return the first solution found.
//
//  Example:
//
//    n = 6
//    target = 22
//    w = (/ 1, 2, 4, 8, 16, 32 /)
//
//    choice = (/ 0, 1, 1, 0, 1, 0 /)
//    w(choice) = 2 + 4 + 16 = 22
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of weights.
//
//    Input, int WEIGHT[N], the weights.
//
//    Input, int TARGET, the target value.
//
//    Output, int SUBSET_SUM_SERIAL[N], contains a 1 for each
//    weight that is chosen.  If no solution was found, all entries
//    are returned as -1.
//
{
  int *choice;
  int i;
  int i_max;
  int w_sum;

  i_max = i4_power ( 2, n );

  for ( i = 0; i < i_max; i++ )
  {
//
//  Convert I to a string of binary digits.
//
    choice = i4_to_digits_binary ( i, n );
//
//  Combine the weights whose binary digit is 1.
//
    w_sum = i4vec_dot_product ( n, choice, weight );
//
//  Return if we matched our target sum.
//
    if ( w_sum == target )
    {
      return choice;
    }
    delete [] choice;
  }

  choice = new int[n];
  for ( i = 0; i < n; i++ )
  {
    choice[i] = -1;
  }

  return choice;
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
