# include <cstdlib>
# include <iostream>
# include <ctime>
# include <cstring>
# include <cmath>

using namespace std;

# include "file_name_sequence.hpp"

//****************************************************************************80

char digit_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.
//
//  Example:
//
//     I     C
//   -----  ---
//     0    '0'
//     1    '1'
//   ...    ...
//     9    '9'
//    10    '*'
//   -83    '*'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the digit, which should be between 0 and 9.
//
//    Output, char DIGIT_TO_CH, the appropriate character '0'
//    through '9' or '*'.
//
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
//****************************************************************************80

void filename_inc ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the name is empty, then the routine stops.
//
//    If the name contains no digits, the empty string is returned.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a9to99.txt"     "a0to00.txt"  (wrap around)
//      "cat.txt"        " "           (no digits to increment)
//      " "              STOP!         (error)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_INC - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;

      if ( c == '9' )
      {
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}
//****************************************************************************80

int i4_log_10 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_10 returns the whole part of the logarithm base 10 of an I4.
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

string i4_to_s ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_S converts an I4 to a string.
//
//  Example:
//
//    INTVAL  S
//
//         1  1
//        -1  -1
//         0  0
//      1952  1952
//    123456  123456
//   1234567  1234567
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer to be converted.
//
//    Output, string I4_TO_S, the representation of the integer.
//
{
  int digit;
  int j;
  int length;
  int ten_power;
  string s;
  char s_char[80];
  static double ten = 10.0;

  length = i4_log_10 ( i );

  ten_power = ( int ) ( pow ( ten, length ) );

  if ( i < 0 )
  {
    length = length + 1;
  }
//
//  Add one position for the trailing null.
//
  length = length + 1;

  if ( i == 0 )
  {
    s_char[0] = '0';
    s_char[1] = '\0';
    s = string ( s_char );
    return s;
  }
//
//  Now take care of the sign.
//
  j = 0;
  if ( i < 0 )
  {
    s_char[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
//
//  Find the leading digit of I, strip it off, and stick it into the string.
//
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s_char[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
//
//  Tack on the trailing NULL.
//
  s_char[j] = '\0';
  j = j + 1;

  s = string ( s_char );

  return s;
}
//****************************************************************************80

string i4_to_s0 ( int i, int digits )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_S0 converts an I4 to a string with leading zeros.
//
//  Example:
//
//         I  DIGITS       S
//
//         1       3     001
//        -1       3     -01
//         0       3     000
//      1952       4    1952
//    123456       6  123456
//   1234567       6  ******
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer to be converted.
//
//    Input, int DIGITS, the number of positions available.
//
//    Output, string I4_TO_S0, the representation of the integer.
//
{
  int digit;
  int j;
  int length;
  int length_max = 80;
  string s;
  char s_char[81];
//
//  Determine LEGNTH, the natural "length" of the number.
//
  length = i4_log_10 ( i ) + 1;

  if ( i < 0 )
  {
    length = length + 1;
  }
//
//  If DIGITS is too large, reset it.
//
  digits = i4_min ( digits, length_max );
//
//  If I is too large, return stars.
//
  if ( digits < length )
  {
    for ( j = 0; j < digits; j++ )
    {
      s_char[j] = '*';
    }
    s_char[digits] = '\0';
    s = string ( s_char );
    return s;
  }
//
//  Set the digits to zero.
//
  for ( j = 0; j < digits; j++ )
  {
    s_char[j] = '0';
  }
  s_char[digits] = '\0';
//
//  If I = 0, we can return now.
//
  if ( i == 0 )
  {
    s = string ( s_char );
    return s;
  }
//
//  If negative, take care of the sign.
//
  if ( i < 0 )
  {
    s_char[0] = '-';
    i = abs ( i );
  }
//
//  Mod away.
//
  j = digits - 1;
  while ( i != 0 )
  {
    s_char[j] = digit_to_ch ( i % 10 );
    i = i / 10;
    j = j - 1;
  }
  s = string ( s_char );

  return s;
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
