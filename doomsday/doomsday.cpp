# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "doomsday.hpp"

//****************************************************************************80

int doomsday_gregorian ( int y, int m, int d )

//****************************************************************************80
//
//  Purpose:
//
//    DOOMSDAY_GREGORIAN: weekday given any date in Gregorian calendar.
//
//  Discussion:
//
//    This procedure does not include any procedure to switch to the Julian
//    calendar for dates early enough that that calendar was used instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Conway,
//    Tomorrow is the Day After Doomsday,
//    Eureka,
//    Volume 36, October 1973, pages 28-31.
//
//  Parameters:
//
//    Input, int Y, M, D, the year, month and day of the date.
//    Note that the year must be positive.
//
//    Output, int DOOMSDAY_GREGORIAN, the weekday of the date.
//
{
  int anchor[4] = { 1, 6, 4, 3 };
  int c;
  int drd;
  int drdr;
  int l;
  int mdoom[12] = { 3, 28, 0, 4, 9, 6, 11, 8, 5, 10, 7, 12 };
  int w;
  int ydoom;
  int yy;
  int yy12d;
  int yy12r;
  int yy12r4d;
//
//  Refuse to handle Y <= 0.
//
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "DOOMSDAY_GREGORIAN - Fatal error!\n";
    cerr << "  Y <= 0.\n";
    exit ( 1 );
  }
//
//  Determine the century C.
//
  c = y / 100;
//
//  Determine the last two digits of the year, YY
//
  yy = y % 100;
//
//  Divide the last two digits of the year by 12.
//
  yy = y % 100;
  yy12d = yy / 12;
  yy12r = yy % 12; 
  yy12r4d = yy12r / 4;
  drd = yy12d + yy12r + yy12r4d;
  drdr = drd % 7;
  ydoom = anchor[ ( c - 1 ) % 4 ] + drdr;
  ydoom = i4_wrap ( ydoom, 1, 7 );
//
//  If M = 1 or 2, and leap year, add 1.
//
  if ( ( m == 1 || m == 2 ) && year_is_leap_gregorian ( y ) )
  {
    l = 1;
  }
  else
  {
    l = 0;
  }

  w = ydoom + ( d -  mdoom[m-1] - l );
  w = i4_wrap ( w, 1, 7 );

  return w;
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

string weekday_to_name_common ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_COMMON returns the name of a Common weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string WEEKDAY_TO_NAME_COMMON, the weekday name.
//
{
  string s;
  static string weekday_list[7] = { 
    "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", 
    "Friday", "Saturday" };

  if ( 1 <= w && w <= 7 )
  {
    s = weekday_list[w-1];
    return s;
  }
  else
  {
    cerr << "\n";
    cerr << "WEEKDAY_TO_NAME_COMMON\n";
    cerr << "  Illegal weekday.\n";
    exit ( 1 );
  }
}
//****************************************************************************80

void weekday_values ( int &n_data, int &y, int &m, int &d, int &w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_VALUES returns the day of the week for various dates.
//
//  Discussion:
//
//    The CE or Common Era calendar is used, under the
//    hybrid Julian/Gregorian Calendar, with a transition from Julian
//    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
//
//    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
//    years BC/BCE are indicated by a negative year value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz,
//    Calendrical Calculations: The Millennium Edition,
//    Cambridge University Press, 2001,
//    ISBN: 0 521 77752 6
//    LC: CE12.R45.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA by 1,
//    and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &Y, &M, &D, the Common Era date.
//
//    Output, int &W, the day of the week.  Sunday = 1.
//
{
# define N_MAX 34

  static int d_vec[N_MAX] = {
    30,
     8,
    26,
     3,
     7,
    18,
     7,
    19,
    14,
    18,
    16,
     3,
    26,
    20,
     4,
    25,
    31,
     9,
    24,
    10,
    30,
    24,
    19,
     2,
    27,
    19,
    25,
    29,
    19,
     7,
    17,
    25,
    10,
    18 };
  static int m_vec[N_MAX] = {
     7,
    12,
     9,
    10,
     1,
     5,
    11,
     4,
    10,
     5,
     3,
     3,
     3,
     4,
     6,
     1,
     3,
     9,
     2,
     6,
     6,
     7,
     6,
     8,
     3,
     4,
     8,
     9,
     4,
    10,
     3,
     2,
    11,
     7 };
  static int w_vec[N_MAX] = {
    1,
    4,
    4,
    1,
    4,
    2,
    7,
    1,
    7,
    1,
    6,
    7,
    6,
    1,
    1,
    4,
    7,
    7,
    7,
    4,
    1,
    6,
    1,
    2,
    4,
    1,
    1,
    2,
    2,
    5,
    3,
    1,
    4,
    1 };
  static int y_vec[N_MAX] = {
    - 587,
    - 169,
       70,
      135,
      470,
      576,
      694,
     1013,
     1066,
     1096,
     1190,
     1240,
     1288,
     1298,
     1391,
     1436,
     1492,
     1553,
     1560,
     1648,
     1680,
     1716,
     1768,
     1819,
     1839,
     1903,
     1929,
     1941,
     1943,
     1943,
     1992,
     1996,
     2038,
     2094 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    y = 0;
    m = 0;
    d = 0;
    w = 0;
  }
  else
  {
    y = y_vec[n_data-1];
    m = m_vec[n_data-1];
    d = d_vec[n_data-1];
    w = w_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool year_is_leap_gregorian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_GREGORIAN returns TRUE if the Gregorian year was a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_GREGORIAN, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "YEAR_IS_LEAP_GREGORIAN - Fatal error!\n";
    cerr << "  This function will not accept nonpositive years.\n";
    exit ( 1 );
  }

  if ( ( y % 400 ) == 0 )
  {
    value = true;
  }
  else if ( ( y % 100 ) == 0 )
  {
    value = false;
  }
  else if ( ( y % 4 ) == 0 )
  {
    value = true;
  }
  else
  {
    value = false;
  }

  return value;
}

