# include <cstdlib>
# include <ctime>

# include "wtime.hpp"

using namespace std;

//****************************************************************************80

double wtime ( )

//****************************************************************************80
//
//  Purpose:
// 
//    WTIME reports the elapsed wallclock time.
//
//  Discussion:
//
//    The reliability of this function depends in part on the value of
//    CLOCKS_PER_SECOND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double WTIME, the a reading of the wall clock timer,
//    in seconds.
//
{
  double value;

  value = ( double ) clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

  return value;
}
