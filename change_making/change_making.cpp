# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "change_making.hpp"

//****************************************************************************80

int *change_making_list ( int coin_num, int coin_value[], int target )

//****************************************************************************80
//
//  Purpose:
//
//    CHANGE_MAKING_LIST solves the change making problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int COIN_NUM, the number of coin denomiations.
//
//    Input, int COIN_VALUE[COIN_NUM], the value of each coin.
//    These values should be positive integers.
//
//    Input, int TARGET, the desired sum.
//
//    Output, int CHANGE_MAKING_LIST[TARGET+1], A(T) lists 
//    the smallest number of coins needed to form the sum T, or "Inf" if 
//    it is not possible to form this sum.
//
{
  int *a;
  int i;
  int i4_huge = 2147483647;
  int j;

  a = new int[target+1];
  a[0] = 0;
  for ( j = 1; j <= target; j++ )
  {
    a[j] = i4_huge;
  }
//
//   If T is the value of a coin, then A(T) is 1.
//
  for ( i = 0; i < coin_num; i++ )
  {
    a[coin_value[i]] = 1;
  }
//
//  To compute A(T) in general, consider getting there by adding
//  one coin of value V, and looking at A(T-V).
//
  for ( j = 1; j <= target; j++ )
  {
    for ( i = 0; i < coin_num; i++ )
    {
      if ( 0 <= j - coin_value[i] )
      {
        a[j] = i4_min ( a[j] - 1, a[j-coin_value[i]] ) + 1;
      }
    }
  }

  return a;
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
//    29 August 2006
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
