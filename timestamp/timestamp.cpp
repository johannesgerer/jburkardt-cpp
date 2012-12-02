# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "timestamp.hpp"

//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Discussion:
//
//    The data available to this routine through "CLOCK" is not very reliable,
//    and hence the values of CPU_TIME returned should not be taken too 
//    seriously, especially when short intervals are being timed.
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
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) std::clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

int *time_numbers ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIME_NUMBERS returns the data as a string of integers.
//
//  Example:
//
//    2001  Year
//    5     Month
//    31    Day
//    9     Hour (0-23)
//    45    Minute
//    12    Second
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int TIME_NUMBERS[6], the year, month, day, hour, minute and second.
//
{
  const struct std::tm *tm_ptr;
  std::time_t now;
  int *value;

  now = std::time ( 0 );
  tm_ptr = std::localtime ( &now );

  value = new int[6];

  value[0] = 1900 + tm_ptr->tm_year;
  value[1] = 1 + tm_ptr->tm_mon;
  value[2] = tm_ptr->tm_mday;
  value[3] = tm_ptr->tm_hour;
  value[4] = tm_ptr->tm_min;
  value[5] = tm_ptr->tm_sec;

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
//****************************************************************************80

char *timestring ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTRING returns the current YMDHMS date as a string.
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
//    Output, char *TIMESTRING, a string containing the current YMDHMS date.
//
{
# define TIME_SIZE 40

  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;
  char *s;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  s = new char[TIME_SIZE];

  len = std::strftime ( s, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  return s;
# undef TIME_SIZE
}
