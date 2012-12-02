# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "cycle_floyd.hpp"

//****************************************************************************80

void cycle_floyd ( int f ( int i ), int x0, int &lam, int &mu )

//****************************************************************************80
//
//  Purpose:
//
//    CYCLE_FLOYD finds a cycle in an iterated mapping using Floyd's method.
//
//  Discussion:
//
//    Suppose we a repeatedly apply a function f(), starting with the argument
//    x0, then f(x0), f(f(x0)) and so on.  Suppose that the range of f is finite.
//    Then eventually the iteration must reach a cycle.  Once the cycle is reached,
//    succeeding values stay within that cycle.
//
//    Starting at x0, there is a "nearest element" of the cycle, which is
//    reached after MU applications of f.
//
//    Once the cycle is entered, the cycle has a length LAM, which is the number
//    of steps required to first return to a given value.
//
//    This function uses Floyd's method to determine the values of MU and LAM,
//    given F and X0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 2, Seminumerical Algorithms,
//    Third Edition,
//    Addison Wesley, 1997,
//    ISBN: 0201896842,
//    LC: QA76.6.K64.
//
//  Parameters:
//
//    Input, int F ( int i ), the name of the function 
//    to be analyzed.
//
//    Input, int X0, the starting point.
//
//    Output, int &LAM, the length of the cycle.
//
//    Output, int &MU, the index in the sequence starting
//    at X0, of the first appearance of an element of the cycle.
//
{
  int hare;
  int i;
  int tortoise;

  tortoise = f ( x0 );
  hare = f ( tortoise );

  while ( tortoise != hare )
  {
    tortoise = f ( tortoise );
    hare = f ( f ( hare ) );
  }

  mu = 0;
  tortoise = x0;

  while ( tortoise != hare )
  {
    tortoise = f ( tortoise );
    hare = f ( hare );
    mu = mu + 1;
  }

  lam = 1;
  hare = f ( tortoise );
  while ( tortoise != hare )
  {
    hare = f ( hare );
    lam = lam + 1;
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
