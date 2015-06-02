# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "calpak.hpp"

//****************************************************************************80

double cws_to_jed_gps ( int c, int w, double s )

//****************************************************************************80
//
//  Purpose:
//
//    CWS_TO_JED_GPS converts a GPS CWS date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int C, int W, double S, 
//    the GPS cycle/week/second date.
//
//    Output, double CWS_TO_JED_GPS, the corresponding Julian Ephemeris Date.
//
{
  double d;
  double jed;
  double jed_epoch;

  jed_epoch = epoch_to_jed_gps ( );

  d = ( double ) ( 7 * ( 1024 * c + w ) ) + s / ( 24.0 * 60.0 * 60.0 );

  jed = jed_epoch + d;

  return jed;
}
//****************************************************************************80

double datenum_to_jed ( double dn )

//****************************************************************************80
//
//  Purpose:
//
//    DATENUM_TO_JED converts a MATLAB date number to a JED.
//
//  Discussion:
//
//    The MATLAB "datenum" function accepts a string defining
//    a date and returns a MATLAB date number:
//
//      dn = datenum ( 'Aug 17 1939' )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double DN, a MATLAB date number.
//
//    Output, double JED, the Julian Ephemeris Date.
//
{
  double jed;

  jed = dn + epoch_to_jed_matlab ( );

  return jed;
}
//****************************************************************************80

void day_borrow_alexandrian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_ALEXANDRIAN borrows days from months in an Alexandrian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_alexandrian ( y, m );

    days = month_length_alexandrian ( y, m );

    d = d + days;
  }
  return;
}
//****************************************************************************80

void day_borrow_common ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_COMMON borrows days from months in a Common date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_common ( y, m );

    days = month_length_common ( y, m );

    d = d + days;
  }

  return;
}
//****************************************************************************80

void day_borrow_eg_civil ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_EG_CIVIL borrows days from months in an Egyptian Civil date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_eg_civil ( y, m );

    days = month_length_eg_civil ( y, m );

    d = d + days;
  }

  return;
}
//****************************************************************************80

void day_borrow_english ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_ENGLISH borrows days from months in an English date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_english ( y, m );

    days = month_length_english ( y, m );

    d = d + days;
  }

  return;
}
//****************************************************************************80

void day_borrow_gregorian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_GREGORIAN borrows days from months in a Gregorian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_gregorian ( y, m );

    days = month_length_gregorian ( y, m );

    d = d + days;
  }
  return;
}
//****************************************************************************80

void day_borrow_hebrew ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_HEBREW borrows days from months in a Hebrew date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_hebrew ( y, m );

    days = month_length_hebrew ( y, m );

    d = d + days;
  }

  return;
}
//****************************************************************************80

void day_borrow_islamic ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_ISLAMIC borrows days from months in an Islamic date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_islamic ( y, m );

    days = month_length_islamic ( y, m );

    d = d + days;
  }

  return;
}
//****************************************************************************80

void day_borrow_julian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_JULIAN borrows days from months in a Julian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_julian ( y, m );

    days = month_length_julian ( y, m );

    d = d + days;
  }
  return;
}
//****************************************************************************80

void day_borrow_republican ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_REPUBLICAN borrows days from months in a Republican date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_republican ( y, m );

    days = month_length_republican ( y, m );

    d = d + days;
  }
  return;
}
//****************************************************************************80

void day_borrow_roman ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_BORROW_ROMAN borrows days from months in a Roman date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, a year, month, and day
//    representing a date.  On input, D might be negative.  On output,
//    M should have decreased by one month, and D gone up by the
//    number of days in the month we "cashed in".  Y may be affected
//    if the input value of M was 1.
//
{
  int days;

  while ( d <= 0 )
  {
    m = m - 1;

    month_borrow_roman ( y, m );

    days = month_length_roman ( y, m );

    d = d + days;
  }
  return;
}
//****************************************************************************80

void day_carry_alexandrian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_ALEXANDRIAN carries days to months in an Alexandrian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;

  days = month_length_alexandrian ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_alexandrian ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_alexandrian ( y, m );
  }
  return;
}
//****************************************************************************80

void day_carry_common ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_COMMON carries days to months in a Common date.
//
//  Discussion:
//
//    While ( number of days in M ) < D:
//      decrease the day D by the number of days in the month M;
//      increase M by 1;
//      if necessary, adjust Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, int &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;
//
//  If the date is in the transition month, deflate it,
//  so we can perform ordinary arithmetic.
//
  deflate_common ( y, m, d );

  days = month_length_common ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_common ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_common ( y, m );
  }
//
//  If the date is in the transition month, inflate it.
//
  inflate_common ( y, m, d );

  return;
}
//****************************************************************************80

void day_carry_eg_civil ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_EG_CIVIL carries days to months in an Egyptian Civil date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, int &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;

  days = month_length_eg_civil ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_eg_civil ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_eg_civil ( y, m );
  }
  return;
}
//****************************************************************************80

void day_carry_english ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_ENGLISH carries days to months in an English date.
//
//  Discussion:
//
//    While ( number of days in M ) < D:
//      decrease the day D by the number of days in the month M;
//      increase M by 1;
//      if necessary, adjust Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;
//
//  If the date is in the transition month, deflate it,
//  so we can perform ordinary arithmetic.
//
  deflate_english ( y, m, d );

  days = month_length_english ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_english ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_english ( y, m );
  }
//
//  If the date is in the transition month, inflate it.
//
  inflate_english ( y, m, d );

  return;
}
//****************************************************************************80

void day_carry_gregorian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_GREGORIAN carries days to months in a Gregorian date.
//
//  Discussion:
//
//    While ( number of days in M ) < D:
//      decrease the day D by the number of days in the month M;
//      increase M by 1;
//      if necessary, adjust Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;

  days = month_length_gregorian ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_gregorian ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_gregorian ( y, m );
  }
  return;
}
//****************************************************************************80

void day_carry_hebrew ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_HEBREW carries days to months in a Hebrew date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;

  days = month_length_hebrew ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_hebrew ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_hebrew ( y, m );
  }
  return;
}
//****************************************************************************80

void day_carry_islamic ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_ISLAMIC carries days to months in an Islamic date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;

  days = month_length_islamic ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_islamic ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_islamic ( y, m );
  }
  return;
}
//****************************************************************************80

void day_carry_julian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_JULIAN carries days to months in a Julian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int Y, &M, &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;

  days = month_length_julian ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_julian ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_julian ( y, m );
  }
  return;
}
//****************************************************************************80

void day_carry_republican ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_REPUBLICAN carries days to months in a Republican date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int Y, &M, &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;

  days = month_length_republican ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_republican ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_republican ( y, m );
  }
  return;
}
//****************************************************************************80

void day_carry_roman ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_CARRY_ROMAN carries days to months in a Roman date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//    On output, D is between 1 and the number of days in M.
//
{
  int days;

  days = month_length_roman ( y, m );

  while ( days < d )
  {
    d = d - days;
    m = m + 1;
    days = month_length_roman ( y, m );
//
//  Make sure the month is not too big.
//
    month_carry_roman ( y, m );
  }
  return;
}
//****************************************************************************80

void day_list_common ( int y1, int m1, int d1, int y2, int m2, int d2 )

//****************************************************************************80
//
//  Purpose:
//
//    DAY_LIST_COMMON prints a list of days between two dates.
//
//  Discussion:
//
//    Given the dates of September 25, 2005 and October 2, 2005,
//    the routine should print out:
//
//    Sun, 25 Sep 2005 
//    Mon, 26 Sep 2005 
//    Tue, 27 Sep 2005 
//    Wed, 28 Sep 2005 
//    Thu, 29 Sep 2005 
//    Fri, 30 Sep 2005 
//    Sat,  1 Oct 2005 
//    Sun,  2 Oct 2005 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, the first date.
//
//    Input, int Y2, M2, D2, the second date.
//
{
  char cmp;
  int d;
  double f;
  int m;
  string m_name;
  int w;
  string w_name;
  int y;

  y = y1;
  m = m1;
  d = d1;
  f = 0.0;

  cmp = '<';
  
  while ( cmp != '>' )
  {
    w = ymdf_to_weekday_common ( y, m, d, f );

    w_name = weekday_to_name_common3 ( w );

    m_name = month_to_month_name_common3 ( m );

    cout << w_name << ", "
         << setw(2) << d << " " << m_name << " " << y << "\n";

    ymdf_next_common ( y, m, d, f, y, m, d, f );

    cmp = ymdf_compare ( y, m, d, f, y2, m2, d2, f );
  }

  return;
}
//****************************************************************************80

int days_before_month_common ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    DAYS_BEFORE_MONTH_COMMON returns the number of days before a Common month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int DAYS_BEFORE_MONTH_COMMON, the number of 
//    days in the year before the first day of the given month.
//
{
  int days;
  int m2;
  int mdays[12] = {
     0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Check the input.
//
  ym_check_common ( y2, m2 );

  days = mdays[m2-1];

  if ( 2 < m2 && year_is_leap_common ( y2 ) )
  {
    days = days + 1;
  }

  return days;
}
//****************************************************************************80

int days_before_month_gregorian ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    DAYS_BEFORE_MONTH_GREGORIAN returns the number of days before a Gregorian month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int DAYS_BEFORE_MONTH_GREGORIAN, the number of 
//    days in the year before the first day of the given month.
//
{
  int days;
  int m2;
  int mdays[12] = {
     0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Check the input.
//
  ym_check_gregorian ( y2, m2 );

  days = mdays[m2-1];

  if ( 2 < m2 && year_is_leap_gregorian ( y2 ) )
  {
    days = days + 1;
  }

  return days;
}
//****************************************************************************80

int days_before_month_julian ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    DAYS_BEFORE_MONTH_JULIAN returns the number of days before a Julian month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int DAYS_BEFORE_MONTH_JULIAN, the number of 
//    days in the year before the first day of the given month.
//
{
  int days;
  int m2;
  int mdays[12] = {
     0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Check the input.
//
  ym_check_julian ( y2, m2 );

  days = mdays[m2-1];

  if ( 2 < m2 && year_is_leap_julian ( y2 ) )
  {
    days = days + 1;
  }

  return days;
}
//****************************************************************************80

void deflate_common ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DEFLATE_COMMON "deflates" dates in the Common Calendar transition month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//
{
  if ( y == 1582 )
  {
    if ( m == 10 )
    {
      if ( 15 <= d )
      {
        d = d - 10;
      }
    }
  }
  return;
}
//****************************************************************************80

void deflate_english ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DEFLATE_ENGLISH "deflates" dates in the English Calendar transition month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//
{
  if ( y == 1752 )
  {
    if ( m == 9 )
    {
      if ( 14 <= d )
      {
        d = d - 11;
      }
    }
  }
  return;
}
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

void easter_ds ( int y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    EASTER_DS computes the month and day of Easter for a Common year.
//
//  Example:
//
//    Input:
//
//      Y = 2000
//
//    Output:
//
//      M = 4
//      D = 23
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Duffett-Smith,
//    Practical Astronomy With Your Calculator,
//    Third Edition,
//    Cambridge University Press, 1996,
//    ISBN: 0-521-35699-7,
//    LC: QB62.5.D83.
//
//  Parameters:
//
//    Input, int Y, the year, which must be 1583 or greater.
//    (The formula is only valid for years after the Gregorian calendar
//    was adopted.)
//
//    Output, int &M, &D, the month and day of Easter.
//
{
  int a;
  int b;
  int c;
  int dd;
  int e;
  int f;
  int g;
  int h;
  int i;
  int k;
  int l;
  int mm;

  if ( y <= 0 )
  {
    m = -1;
    d = -1;
    return;
  }

  a = year_to_golden_number ( y );

  a = a - 1;

  b = y / 100;
  c = y % 100;

  dd = b / 4;
  e = b % 4;

  f = ( b + 8 ) / 25;
  g = ( b - f + 1 ) / 3;
  h = ( 19 * a + b - dd - g + 15 ) % 30;

  i = c / 4;
  k = c % 4;

  l = ( 32 + 2 * e + 2 * i - h - k ) % 7;
  mm = ( a + 11 * h + 22 * l ) / 451;

  m = ( h + l - 7 * mm + 114 ) / 31;
  d = ( ( h + l - 7 * mm + 114 ) % 31 ) + 1;

  return;
}
//****************************************************************************80

void easter_egr ( int y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    EASTER_EGR computes the month and day of Easter for a Common year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm O,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 375.
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &M, &D, the month and day of Easter.
//
{
  int c;
  int e;
  int g;
  int h;
  int n;
  int p;
  int q;
  int r;
  int s;
  int u;
  int vp;

  if ( y <= 0 )
  {
    m = -1;
    d = -1;
    return;
  }

  p = y + ( y / 4 ) - ( y / 100 ) + ( y / 400 ) - 1;
  n = 7 - ( p % 7 );
  h = y / 100;
  q = h - h / 4;
  g = 1 + ( y % 19 );
  e = ( 57 + 11 * g - q + ( h - ( h - 17 ) / 25 ) / 3 ) % 30;
  u = ( 53 - e ) % 30;
  vp = ( g - 1 + 11 * u ) / 319;
  r = 22 + u - vp;
  c = i4_wrap ( r + 3, 1, 7 );
  s = r + ( ( 7 + n - c ) % 7 );

  m = 3 + ( s / 32 );
  d = i4_wrap ( s, 1, 31 );

  return;
}
//****************************************************************************80

void easter_egr2 ( int y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    EASTER_EGR2 computes the month and day of Easter for a Common year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm P,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 376.
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &M, &D, the month and day of Easter.
//
{
  int a;
  int b;
  int c;
  int dd;
  int e;
  int s;

  if ( y <= 0 )
  {
    m = -1;
    d = -1;
    return;
  }

  a = y / 100;
  b = a - ( a / 4 );
  c = ( y % 19 );
  dd = ( 15 + 19 * c + b - ( a - ( a - 17 ) / 25 ) / 3 ) % 30;
  e = dd - ( c + 11 * dd ) / 319;
  s = 22 + e + ( 140004 - y - ( y / 4 ) + b - e ) % 7;

  m = 3 + ( s / 32 );
  d = i4_wrap ( s, 1, 31 );

  return;
}
//****************************************************************************80

void easter_julian ( int y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    EASTER_JULIAN computes the date of Easter in the Julian calendar.
//
//  Discussion:
//
//    This computation for the date of Easter uses the Dionysian
//    canon that applied to the Julian calendar.  The determination
//    of the date of Easter changed at the same time that the calendar
//    was modified to use the Gregorian system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm M,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 365.
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &M, &D, the month and day of the Julian 
//    calendar on which Easter occurs.
//
{
  int c;
  int e;
  int n;
  int p;
  int r;
  int s;

  if ( y <= 0 )
  {
    m = -1;
    d = -1;
    return;
  }

  p = y + ( y / 4 ) + 4;
  n = 7 - ( p % 7 );

  e = year_to_epact_julian ( y );

  r = 22 + ( ( 53 - e ) % 30 );

  c = i4_wrap ( r + 3, 1, 7 );

  s = r + ( ( 7 + n - c ) % 7 );

  m = 3 + ( s / 32 );
//
//  Use wrapping so that 1 <= D <= 31.
//
  d = i4_wrap ( s, 1, 31 );

  return;
}
//****************************************************************************80

void easter_julian2 ( int y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    EASTER_JULIAN2 computes the date of Easter in the Julian calendar.
//
//  Discussion:
//
//    This computation for the date of Easter uses the Dionysian
//    canon that applied to the Julian calendar.  The determination
//    of the date of Easter changed at the same time that the calendar
//    was modified to use the Gregorian system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm N,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 365.
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &M, &D, the month and day of the Julian calendar
//    on which Easter occurs.
//
{
  int a;
  int b;
  int s;

  if ( y <= 0 )
  {
    m = -1;
    d = -1;
    return;
  }

  a = year_to_golden_number ( y );
  a = a - 1;

  b = 22 + ( ( 225 - 11 * a ) % 30 );
  s = b + ( ( 56 + 6 * y - ( y / 4 ) - b ) % 7 );

  m = 3 + ( s / 32 );
//
//  Use wrapping to ensure that 1 <= D <= 31.
//
  d = i4_wrap ( s, 1, 31 );

  return;
}
//****************************************************************************80

double epoch_to_jed_akbar ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_AKBAR returns the epoch of the Akbar calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_AKBAR, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2289425.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_alexandrian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_ALEXANDRIAN: epoch of the Alexandrian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_ALEXANDRIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1713262.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_armenian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_ARMENIAN returns the epoch of the Armenian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_ARMENIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1922867.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_bahai ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_BAHAI returns the epoch of the Bahai calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_BAHAI, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2394646.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_bessel ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_BESSEL returns the epoch of the Bessel calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_BESSEL, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2415020.31352;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_chinese ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_CHINESE returns the epoch of the Chinese calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_CHINESE, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 758325.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_common ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_COMMON returns the epoch of the Common calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_COMMON, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1721423.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_coptic ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_COPTIC returns the epoch of the Coptic calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_COPTIC, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1825029.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_deccan ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_DECCAN returns the epoch of the Fasli Deccan calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_DECCAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1936747.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_eg_civil ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_EG_CIVIL: epoch of the Egyptian Civil calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_EG_CIVIL, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1448637.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_eg_lunar ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_EG_LUNAR: epoch of the Egyptian Lunar calendar as a JED.
//
//  Discussion:
//
//    This is just a fake value, making the Egyptian Lunar calendar start
//    at the same data as the Egyptian Civil calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_EG_LUNAR, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1448637.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_english ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_ENGLISH returns the epoch of the English calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_ENGLISH, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1721423.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_ethiopian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_ETHIOPIAN returns the epoch of the Ethiopian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_ETHIOPIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1724220.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_gps ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_GPS returns the epoch of the GPS calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_GPS, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2444244.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_greek ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_GREEK returns the epoch of the Greek calendar as a JED.
//
//  Discussion:
//
//    The Greek Olympiad calendar began on 9 July 776 BC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_GREEK, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1438178.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_gregorian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_GREGORIAN returns the epoch of the Gregorian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_GREGORIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1721425.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_hebrew ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_HEBREW returns the epoch of the Hebrew calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_HEBREW, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 347998.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_hindu_lunar ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_HINDU_LUNAR: epoch of the Hindu lunar calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_HINDU_LUNAR, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1741959.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_hindu_solar ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_HINDU_SOLAR: epoch of the Hindu solar calendar as a JED.
//
//  Discussion:
//
//    This is the beginning of the Kali Yuga era.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_HINDU_SOLAR, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 588465.75;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_islamic_a ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_ISLAMIC_A returns the epoch of the Islamic A calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_ISLAMIC_A, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1948438.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_islamic_b ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_ISLAMIC_B returns the epoch of the Islamic B calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_ISLAMIC_B, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1948439.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_jed ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_JED returns the epoch of the JED as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_JED, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 0.0;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_jelali ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_JELALI returns the epoch of the Jelali calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_JELALI, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2114872.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_julian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_JULIAN returns the epoch of the Julian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_JULIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1721423.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_khwarizmian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_KHWARIZMIAN: epoch of the Khwarizmian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_KHWARIZMIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1952067.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_macedonian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_MACEDONIAN: epoch of the Macedonian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_MACEDONIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1607708.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_matlab ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_MATLAB: epoch of the "MATLAB calendar" as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_MATLAB, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1721058.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_mayan_long ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_MAYAN_LONG: epoch of the Mayan long count calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_MAYAN_LONG, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 584282.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_mjd ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_MJD returns the epoch of the MJD calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_MJD, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2400000.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_nyt ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_NYT returns the epoch of the NYT calendar as a JED.
//
//  Discussion:
//
//    The "epoch" of the NYT calendar is the mythical date when issue "0"
//    would have been printed, namely, a tad past midnight, 17 September 1851.
//
//    Volume #1, Issue #1 was printed on 18 September 1851.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_NYT, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2397382.5;
//
//  The following value is effectively the JED we are using for an
//  epoch set to the nominal issue number 50,000.
//
// jed = 2449790.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_persian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_PERSIAN returns the epoch of the Persian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_PERSIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1952062.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_persian_solar ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_PERSIAN_SOLAR: epoch of the Persian solar calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_PERSIAN_SOLAR, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1948320.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_rd ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_RD returns the epoch of the RD calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_RD, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1721425.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_republican ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_REPUBLICAN: epoch of the Republican calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_REPUBLICAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2375839.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_roman ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_ROMAN returns the epoch of the Roman calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//     12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_ROMAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1446389.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_saka ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_SAKA returns the epoch of the Saka calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_SAKA, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1749994.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_soor_san ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_SOOR_SAN: epoch of the Fasli Soor San calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_SOOR_SAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1940351.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_syrian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_SYRIAN returns the epoch of the Syrian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_SYRIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1607738.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_unix ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_UNIX returns the epoch of the UNIX calendar as a JED.
//
//  Discussion:
//
//    The UNIX Epoch is taken to be the first second of 1 January 1970.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_UNIX, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2440587.50;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_y2k ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_Y2K returns the epoch of the Y2K calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_Y2K, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 2451544.5;

  return jed;
}
//****************************************************************************80

double epoch_to_jed_zoroastrian ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPOCH_TO_JED_ZOROASTRIAN: epoch of the Zoroastrian calendar as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EPOCH_TO_JED_ZOROASTRIAN, the Julian Ephemeris Date of the epoch.
//
{
  double jed;

  jed = 1862836.5;

  return jed;
}
//****************************************************************************80

void frac_borrow_common ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_BORROW_COMMON borrows fractions from days in a Common YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  while ( f < 0.0 )
  {
    f = f + 1.0;
    d = d - 1;
  }
  day_borrow_common ( y, m, d );

  return;
}
//****************************************************************************80

void frac_borrow_english ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_BORROW_ENGLISH borrows fractions from days in an English YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, int &D, double &F, a YMDF date.
//
{
  while ( f < 0.0 )
  {
    f = f + 1.0;
    d = d - 1;
  }
  day_borrow_english ( y, m, d );

  return;
}
//****************************************************************************80

void frac_borrow_gregorian ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_BORROW_GREGORIAN borrows fractions from days in a Gregorian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  while ( f < 0.0 )
  {
    f = f + 1.0;
    d = d - 1;
  }
  day_borrow_gregorian ( y, m, d );

  return;
}
//****************************************************************************80

void frac_borrow_hebrew ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_BORROW_HEBREW borrows fractions from days in a Hebrew YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  while ( f < 0.0 )
  {
    f = f + 1.0;
    d = d - 1;
  }
  day_borrow_hebrew ( y, m, d );

  return;
}
//****************************************************************************80

void frac_borrow_islamic ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_BORROW_ISLAMIC borrows fractions from days in an Islamic YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  while ( f < 0.0 )
  {
    f = f + 1.0;
    d = d - 1;
  }
  day_borrow_islamic ( y, m, d );

  return;
}
//****************************************************************************80

void frac_borrow_julian ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_BORROW_JULIAN borrows fractions from days in a Julian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  while ( f < 0.0 )
  {
    f = f + 1.0;
    d = d - 1;
  }
  day_borrow_julian ( y, m, d );

  return;
}
//****************************************************************************80

void frac_borrow_republican ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_BORROW_REPUBLICAN borrows fractions from days in a Republican YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  while ( f < 0.0 )
  {
    f = f + 1.0;
    d = d - 1;
  }
  day_borrow_republican ( y, m, d );

  return;
}
//****************************************************************************80

void frac_borrow_roman ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_BORROW_ROMAN borrows fractions from days in a Roman YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  while ( f < 0.0 )
  {
    f = f + 1.0;
    d = d - 1;
  }
  day_borrow_roman ( y, m, d );

  return;
}
//****************************************************************************80

void frac_carry_common ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_CARRY_COMMON carries fractions from days in a Common YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  int days;

  if ( f < 1.0 )
  {
    return;
  }

  days = ( int ) ( f );

  f = f - ( double ) ( days );
  d = d + days;

  day_carry_common ( y, m, d );

  return;
}
//****************************************************************************80

void frac_carry_english ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_CARRY_ENGLISH carries fractions from days in an English YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  int days;

  if ( f < 1.0 )
  {
    return;
  }

  days = ( int ) ( f );

  f = f - ( double ) ( days );
  d = d + days;

  day_carry_english ( y, m, d );

  return;
}
//****************************************************************************80

void frac_carry_gregorian ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_CARRY_GREGORIAN carries fractions from days in a Gregorian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  int days;

  if ( f < 1.0 )
  {
    return;
  }

  days = ( int ) ( f );

  f = f - ( double ) ( days );
  d = d + days;

  day_carry_gregorian ( y, m, d );

  return;
}
//****************************************************************************80

void frac_carry_hebrew ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_CARRY_HEBREW carries fractions from days in a Hebrew YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  int days;

  if ( f < 1.0 )
  {
    return;
  }

  days = ( int ) ( f );

  f = f - ( double ) ( days );
  d = d + days;

  day_carry_hebrew ( y, m, d );

  return;
}
//****************************************************************************80

void frac_carry_islamic ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_CARRY_ISLAMIC carries fractions from days in an Islamic YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  int days;

  if ( f < 1.0 )
  {
    return;
  }

  days = ( int ) ( f );

  f = f - ( double ) ( days );
  d = d + days;

  day_carry_islamic ( y, m, d );

  return;
}
//****************************************************************************80

void frac_carry_julian ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_CARRY_JULIAN carries fractions from days in a Julian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  int days;

  if ( f < 1.0 )
  {
    return;
  }

  days = ( int ) ( f );

  f = f - ( double ) ( days );
  d = d + days;

  day_carry_julian ( y, m, d );

  return;
}
//****************************************************************************80

void frac_carry_republican ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_CARRY_REPUBLICAN carries fractions from days in a Republican YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  int days;

  if ( f < 1.0 )
  {
    return;
  }

  days = ( int ) ( f );

  f = f - ( double ) ( days );
  d = d + days;

  day_carry_republican ( y, m, d );

  return;
}
//****************************************************************************80

void frac_carry_roman ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_CARRY_ROMAN carries fractions from days in a Roman YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, a YMDF date.
//
{
  int days;

  if ( f < 1.0 )
  {
    return;
  }

  days = ( int ) ( f );

  f = f - ( double ) ( days );
  d = d + days;

  day_carry_roman ( y, m, d );

  return;
}
//****************************************************************************80

void frac_to_hms ( double f, int &h, int &m, int &s )

//****************************************************************************80
//
//  Purpose:
//
//    FRAC_TO_HMS converts a fractional day into hours, minutes, seconds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double F, a day fraction between 0.0 and 1.0.
//
//    Output, int &H, &M, &S, the equivalent hours, minutes
//    and seconds.
//
{
  double f2;

  f2 = f;

  f2 = 24.0 * f2;
  h = ( int ) ( f2 );
  f2 = f2 - ( double ) ( h );

  f2 = 60.0 * f2;
  m = ( int ) ( f2 );
  f2 = f2 - ( double ) ( m );

  f2 = 60.0 * f2;
  s = ( int ) ( f2 );
  f2 = f2 - ( double ) ( s );

  return;
}
//****************************************************************************80

void hour_borrow_common ( int &y, int &m, int &d, int &h )

//****************************************************************************80
//
//  Purpose:
//
//    HOUR_BORROW_COMMON "borrows" a day of hours.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, &H, the year, month, day
//    and hour of the date.  The value of H is presumably negative, and
//    so hours will be "borrowed" to make H positive.
//
{
  while ( h <= 0 )
  {
    h = h + 24;
    d = d - 1;

    day_borrow_common ( y, m, d );
  }

  return;
}
//****************************************************************************80

void hour_carry_common ( int &y, int &m, int &d, int &h )

//****************************************************************************80
//
//  Purpose:
//
//    HOUR_CARRY_COMMON is given a YMDH date, and carries hours to days.
//
//  Algorithm:
//
//    While 24 < H:
//
//      decrease H by the number of hours in a day;
//      increase D by 1;
//      if necessary, adjust M and Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, &H, the year, month, day
//    and hour of the date.  On input, H is presumably 24 or greater.
//
{
  while ( 24 < h )
  {
    h = h - 24;
    d = d + 1;

    day_carry_common ( y, m, d );
  }

  return;
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

void i4_swap ( int &i, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &I, &J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = i;
  i = j;
  j = k;
 
  return;
}
//****************************************************************************80

char i4_to_a ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_A returns the I-th alphabetic character.
//
//  Example:
//
//    I  I4_TO_A
//
//   -8  ' '
//    0  ' '
//    1  'A'
//    2  'B'
//   ..
//   26  'Z'
//   27  'a'
//   52  'z'
//   53  ' '
//   99  ' '
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the index of the letter to be returned.
//    0 is a space;
//    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
//    27 through 52 requests 'a' through 'z', (ASCII 97:122);
//
//    Output, char I4_TO_A, the requested alphabetic letter.
//
{
  char value;

  if ( i <= 0 )
  {
    value = ' ';
  }
  else if ( 1 <= i && i <= 26 )
  {
    value = 'A' + i - 1;
  }
  else if ( 27 <= i && i <= 52 )
  {
    value = 'a' + i - 27;
  }
  else
  {
    value = ' ';
  }
  return value;
}
//****************************************************************************80

string i4_to_string ( int i4, string format )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  char i4_char[80];
  string i4_string;

  sprintf ( i4_char, format.c_str ( ), i4 );

  i4_string = string ( i4_char );

  return i4_string;
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

void inflate_common ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    INFLATE_COMMON "inflates" dates in the Common Calendar transition month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//
{
  if ( y == 1582 )
  {
    if ( m == 10 )
    {
      if ( 5 <= d )
      {
        d = d + 10;
      }
    }
  }

  return;
}
//****************************************************************************80

void inflate_english ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    INFLATE_ENGLISH "inflates" dates in the English Calendar transition month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date.
//
{
  if ( y == 1752 )
  {
    if ( m == 9 )
    {
      if ( 3 <= d )
      {
        d = d + 11;
      }
    }
  }
  return;
}
//****************************************************************************80

void j_borrow_common ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_BORROW_COMMON borrows year-days from years in a Common date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  while ( j <= 0 )
  {
    y = y - 1;

    days = year_length_common ( y );

    j = j + days;
  }
  return;
}
//****************************************************************************80

void j_borrow_english ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_BORROW_ENGLISH borrows year-days from years in an English date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  while ( j <= 0 )
  {
    y = y - 1;

    days = year_length_english ( y );

    j = j + days;
  }
  return;
}
//****************************************************************************80

void j_borrow_gregorian ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_BORROW_GREGORIAN borrows year-days from years in a Gregorian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  while ( j <= 0 )
  {
    y = y - 1;

    days = year_length_gregorian ( y );

    j = j + days;
  }
  return;
}
//****************************************************************************80

void j_borrow_hebrew ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_BORROW_HEBREW borrows year-days from years in a Hebrew date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  while ( j <= 0 )
  {
    y = y - 1;

    days = year_length_hebrew ( y );

    j = j + days;
  }
  return;
}
//****************************************************************************80

void j_borrow_islamic ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_BORROW_ISLAMIC borrows year-days from years in an Islamic date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  while ( j <= 0 )
  {
    y = y - 1;

    days = year_length_islamic ( y );

    j = j + days;
  }
  return;
}
//****************************************************************************80

void j_borrow_julian ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_BORROW_JULIAN borrows year-days from years in a Julian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  while ( j <= 0 )
  {
    y = y - 1;

    days = year_length_julian ( y );

    j = j + days;
  }
  return;
}
//****************************************************************************80

void j_borrow_republican ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_BORROW_REPUBLICAN borrows year-days from years in a Republican date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  while ( j <= 0 )
  {
    y = y - 1;

    days = year_length_republican ( y );

    j = j + days;
  }
  return;
}
//****************************************************************************80

void j_borrow_roman ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_BORROW_ROMAN borrows year-days from years in a Roman date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  while ( j <= 0 )
  {
    y = y - 1;

    days = year_length_roman ( y );

    j = j + days;
  }
  return;
}
//****************************************************************************80

void j_carry_common ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_CARRY_COMMON carries year-days to years in a Common date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  for ( ; ; )
  {
    days = year_length_common ( y );

    if ( j < days )
    {
      break;
    }
    j = j - days;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void j_carry_english ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_CARRY_ENGLISH carries year-days to years in an English date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  for ( ; ; )
  {
    days = year_length_english ( y );

    if ( j < days )
    {
      break;
    }
    j = j - days;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void j_carry_gregorian ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_CARRY_GREGORIAN carries year-days to years in a Gregorian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  for ( ; ; )
  {
    days = year_length_gregorian ( y );

    if ( j < days )
    {
      break;
    }
    j = j - days;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void j_carry_hebrew ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_CARRY_HEBREW carries year-days to years in a Hebrew date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  for ( ; ; )
  {
    days = year_length_hebrew ( y );

    if ( j < days )
    {
      break;
    }
    j = j - days;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void j_carry_islamic ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_CARRY_ISLAMIC carries year-days to years in an Islamic date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  for ( ; ; )
  {
    days = year_length_islamic ( y );

    if ( j < days )
    {
      break;
    }
    j = j - days;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void j_carry_julian ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_CARRY_JULIAN carries year-days to years in a Julian date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  for ( ; ; )
  {
    days = year_length_julian ( y );

    if ( j < days )
    {
      break;
    }
    j = j - days;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void j_carry_republican ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_CARRY_REPUBLICAN carries year-days to years in a Republican date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  for ( ; ; )
  {
    days = year_length_republican ( y );

    if ( j < days )
    {
      break;
    }
    j = j - days;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void j_carry_roman ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    J_CARRY_ROMAN carries year-days to years in a Roman date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, a YJ date.
//
{
  int days;

  for ( ; ; )
  {
    days = year_length_roman ( y );

    if ( j < days )
    {
      break;
    }
    j = j - days;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void jed_check ( double jed )

//****************************************************************************80
//
//  Purpose:
//
//    JED_CHECK checks a Julian Ephemeris Date.
//
//  Discussion:
//
//    The routine returns an error if JED < 0, although there is no
//    reason why such dates are invalid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
{
  if ( jed < 0.0 )
  {
    cerr << "\n";
    cerr << "JED_CHECK - Fatal error!\n";
    cerr << "  Input JED < 0.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double jed_test ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TEST returns some "interesting" JED's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Bonnie Blackburn, Leofranc Holford-Stevens,
//    The Oxford Companion to the Year,
//    Oxford, 1999.
//
//    Frank Parise, editor,
//    The Book of Calendars,
//    Facts on File, Inc, 1982,
//    CE11.K4 / 529.3.
//
//    Edward Reingold, Nachum Dershowitz,
//    Calendrical Calculations, the Millennium Edition,
//    Cambridge, 2002,
//    CE12.R45 / 529.3-dc21
//
//    Edward Richards,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999.
//
//  Parameters:
//
//    Input, int I, the test date requested.
//
//    Output, double JED_TEST, the Julian Ephemeris Date.
//    If I is less than 1, or greater than the number of test dates
//    available, JED is returned as -1.0.
//
{
  double jed;
  double jed_epoch_50000;
//
//  JED Epoch:
//  Beginning of current Scaliger cycle.
//  Monday, Noon, 1 January 4713 BCE/Julian
//
  if ( i == 1 )
  {
    jed = epoch_to_jed_jed ( );
  }
//
//  The day after the JED Epoch.
//  Tuesday, Noon, 2 January 4713 BCE/Julian
//
  else if ( i == 2 )
  {
    jed = epoch_to_jed_jed ( );
    jed = jed + 1.0;
  }
//
//  Archbishop James Ussher's estimate of the date of Creation,
//  (Noon), 23 October 4004 BCE/Julian
//
  else if ( i == 3 )
  {
    jed = 259258.000;
  }
//
//  Hebrew Epoch.
//  7 October 3761 BCE/Julian
//
  else if ( i == 4 )
  {
    jed = epoch_to_jed_hebrew ( );
  }
//
//  Mayan Long Count Epoch.
//  6 September 3114 BCE/Julian
//  (Reingold and Dershowitz)
//
  else if ( i == 5 )
  {
    jed = epoch_to_jed_mayan_long ( );
  }
//
//  Hindu Solar Epoch.
//  Beginning of the Kali Yuga age.
//  18 February 3102 BCE/Julian
//
  else if ( i == 6 )
  {
    jed = epoch_to_jed_hindu_solar ( );
  }
//
//  Chinese Epoch.
//  8 March 2637 BCE/Julian
//
  else if ( i == 7 )
  {
    jed = epoch_to_jed_chinese ( );
  }
//
//  Greek Olympiad Epoch
//  9 July 776 BCE/Julian
//
  else if ( i == 8 )
  {
    jed = epoch_to_jed_greek ( );
  }
//
//  Roman Epoch
//  Ab Urbe Condita
//  1 January 753 BCE/Julian
//
  else if ( i == 9 )
  {
    jed = epoch_to_jed_roman ( );
  }
//
//  Egyptian Civil Calendar Epoch.
//  Ascension of Nabonassar to throne of Babylon.
//  26 February 747 BCE/Julian
//
  else if ( i == 10 )
  {
    jed = epoch_to_jed_eg_civil ( );
  }
//
//  Egyptian Lunar Calendar Epoch.
//  (Don't really know where to set this...)
//  Ascension of Nabonassar to throne of Babylon.
//  26 February 747 BCE/Julian
//
  else if ( i == 11 )
  {
    jed = epoch_to_jed_eg_lunar ( );
  }
//
//  Macedonian Epoch
//  1 September 312 BCE/Julian
//
  else if ( i == 12 )
  {
    jed = epoch_to_jed_macedonian ( );
  }
//
//  Syrian Epoch
//  1 October 312 BCE/Julian
//
  else if ( i == 13 )
  {
    jed = epoch_to_jed_syrian ( );
  }
//
//  Alexandrian Epoch
//  29 August 23 BCE/Julian
//
  else if ( i == 14 )
  {
    jed = epoch_to_jed_alexandrian ( );
  }
//
//  "1 January, 0 BC"?  MATLAB epoch?
//
  else if ( i == 15 )

    jed = epoch_to_jed_matlab ( );
//
//  Julian Epoch MINUS ONE DAY
//  Friday, 31 December 1 BCE/Julian
//
  else if ( i == 16 )
  {
    jed = epoch_to_jed_julian ( );
    jed = jed - 1.0;
  }
//
//  Julian Epoch
//  Saturday, 1 January 1 CE/Julian
//
  else if ( i == 17 )
  {
    jed = epoch_to_jed_julian ( );
  }
//
//  Gregorian Epoch
//  Monday, 3 January 1 CE/Julian
//  Monday, 1 January 1 Gregorian
//
  else if ( i == 18 )
  {
    jed = epoch_to_jed_gregorian ( );
  }
//
//  RD: Reingold and Dershowitz Epoch
//  Monday, 3 January 1 CE/Julian
//  Monday, 1 January 1 Gregorian
//
  else if ( i == 19 )
  {
    jed = epoch_to_jed_rd ( );
  }
//
//  Ethiopian Epoch
//  29 August 8 CE/Julian
//  (Reingold and Dershowitz)
//
  else if ( i == 20 )
  {
    jed = epoch_to_jed_ethiopian ( );
  }
//
//  Hindu Lunar Epoch, the Vikrama
//  24 March 57 CE/Julian
//  (The actual day and month are not specified by RD)
//  (Reingold and Dershowitz)
//
  else if ( i == 21 )
  {
    jed = epoch_to_jed_hindu_lunar ( );
  }
//
//  Saka Epoch
//  4 March 79 CE/Julian
//
  else if ( i == 22 )
  {
    jed = epoch_to_jed_saka ( );
  }
//
//  Coptic Epoch
//  29 August 284 CE/Julian
//
  else if ( i == 23 )
  {
    jed = epoch_to_jed_coptic ( );
  }
//
//  Zoroastrian Epoch.
//  3 March 388 CE/Julian
//
   else if ( i == 24 )
  {
     jed = epoch_to_jed_zoroastrian ( );
  }
//
//  Armenian Epoch
//  11 July 552 CE/Julian
//
  else if ( i == 25 )
  {
    jed = epoch_to_jed_armenian ( );
  }
//
//  Fasli Deccan Epoch
//  12 July 590 CE/Julian
//
   else if ( i == 26 )
  {
     jed = epoch_to_jed_deccan ( );
  }
//
//  Fasli Soor San Epoch
//  24 May 600 CE/Julian
//
   else if ( i == 27 )
  {
     jed = epoch_to_jed_soor_san ( );
  }
//
//  Persian Solar Epoch
//  19 March 622 CE/Julian
//
  else if ( i == 28 )
  {
    jed = epoch_to_jed_persian_solar ( );
  }
//
//  Islamic A Epoch
//  Thursday, 15 July 622 CE/Julian
//
  else if ( i == 29 )
  {
    jed = epoch_to_jed_islamic_a ( );
  }
//
//  Islamic B Epoch
//  Friday, 16 July 622 CE/Julian
//
  else if ( i == 30 )
  {
    jed = epoch_to_jed_islamic_b ( );
  }
//
//  Yazdegerd Epoch
//  16 June 632 CE
//
  else if ( i == 31 )
  {
    jed = epoch_to_jed_persian ( );
  }
//
//  Khwarizmian Epoch
//  21 June 632 CE/Julian
//
  else if ( i == 32 )
  {
    jed = epoch_to_jed_khwarizmian ( );
  }
//
//  Battle of Hastings.
//  Saturday, 14 October 1066 CE/Julian.
//           (20 October 1066 Gregorian.)
//
  else if ( i == 33 )
  {
    jed = 2110700.5;
  }
//
//  Jelali Epoch
//  17 March 1078 CE/Julian
//
  else if ( i == 34 )
  {
    jed = epoch_to_jed_jelali ( );
  }
//
//  Akbar Epoch
//  9 February 1556 CE/Julian
//  19 February 1556 Gregorian
//
  else if ( i == 35 )
  {
    jed = epoch_to_jed_akbar ( );
  }
//
//  Common Era calendar transition:
//  Noon of the last day of Julian calendar usage.
//  Thursday, 04 October 1582 CE/English/Julian
//  Thursday, 14 October 1582 Gregorian
//
  else if ( i == 36 )
  {
    jed = transition_to_jed_common ( );
    jed = jed - 0.5;
  }
//
//  Common Era calendar transition:
//  Noon of the first day of Gregorian calendar usage.
//  Friday, 05 October 1582 English/Julian
//  Friday, 15 October 1582 CE/Gregorian
//
  else if ( i == 37 )
  {
    jed = transition_to_jed_common ( );
    jed = jed + 0.5;
  }
//
//  A day chosen by Lewis Carroll to test his day-of-the-week algorithm,
//  Wednesday, 4 March 1676 CE/Gregorian
//  Wednesday, 23 February 1676 English/Julian
//
  else if ( i == 38 )
  {
    jed = 2333269.5;
  }
//
//  English calendar
//  noon of the last day of Julian calendar usage.
//  02 September 1752 English/Julian
//  13 September 1752 CE/Gregorian
//
  else if ( i == 39 )
  {
    jed = transition_to_jed_english ( );
    jed = jed - 0.5;
  }
//
//  English calendar,
//  noon of the first day of Gregorian calendar usage.
//  03 September 1752 Julian
//  14 September 1752 CE/English/Gregorian
//
  else if ( i == 40 )
  {
    jed = transition_to_jed_english ( );
    jed = jed + 0.5;
  }
//
//  A day chosen by Lewis Carroll to test his day-of-the-week algorithm,
//  Thursday, 18 September 1783 CE/Gregorian
//
  else if ( i == 41 )
  {
    jed = 2372547.5;
  }
//
//  French Republican Epoch
//  Saturday, 11 September 1792 Julian
//  Saturday, 22 September 1792 CE/Gregorian
//
  else if ( i == 42 )
  {
    jed = epoch_to_jed_republican ( );
  }
//
//  Bahai Epoch.
//  9 March 1844 Julian
//  21 March 1844 CE/Gregorian
//
  else if ( i == 43 )
  {
    jed = epoch_to_jed_bahai ( );
  }
//
//  Clive James Lucas test date.
//
  else if ( i == 44 )
  {
    jed = 2394710.50;
  }
//
//  New York Times "epoch" date,
//  fictitious Volume 1, issue #0,
//  17 September 1851
//  (issue #1 was on 18 September 1851):
//
  else if ( i == 45 )
  {
    jed = 2397383.50;
  }
//
//  Modified Julian Date Epoch.
//  17 November 1858 CE/Gregorian
//
  else if ( i == 46 )
  {
    jed = epoch_to_jed_mjd ( );
  }
//
//  NYT issue 10,000
//  24 September 1883
//
  else if ( i == 47 )
  {
    jed_epoch_50000 = 2449790.5;
    jed = jed_epoch_50000 - 40000.0 - 88.0;
  }
//
//  Bessel Year Count Epoch.
//  1 January 1900 CE/Gregorian
//
  else if ( i == 48 )
  {
    jed = epoch_to_jed_bessel ( );
  }
//
//  NYT issue 30,000
//  14 March 1940
//
  else if ( i == 49 )
  {
    jed_epoch_50000 = 2449790.5;
    jed = jed_epoch_50000 - 20000.0 - 88.0;
  }
//
//  NYT issue 40,000
//  ???
//
  else if ( i == 50 )
  {
    jed_epoch_50000 = 2449790.5;
    jed = jed_epoch_50000 - 10000.0 - 88.0;
  }
//
//  UNIX epoch.
//  1 January 1970 CE/Gregorian.
//
  else if ( i == 51 )
  {
    jed = epoch_to_jed_unix ( );
  }
//
//  NYT issue 44027
//  ???
//
  else if ( i == 52 )
  {
    jed_epoch_50000 = 2449790.5;
    jed = jed_epoch_50000 - 5973;
  }
//
//  NYT issue 44028
//  ???
//
  else if ( i == 53 )
  {
    jed_epoch_50000 = 2449790.5;
    jed = jed_epoch_50000 - 5972;
  }
//
//  GPS epoch.
//  6 January 1980 CE/Gregorian
//
  else if ( i == 54 )
  {
    jed = epoch_to_jed_gps ( );
  }
//
//  NYT issue 50,000
//  14 March 1995
//
  else if ( i == 55 )
  {
    jed_epoch_50000 = 2449790.5;
    jed = jed_epoch_50000;
  }
//
//  25 February 1996
//  A Reingold/Dershowitz test date.
//
  else if ( i == 56 )
  {
    jed = 2450138.5;
  }
//
//  Y2K day
//  1 January 2000 CE/Gregorian
//
  else if ( i == 57 )
  {
    jed = epoch_to_jed_y2k ( );
  }
//
//  Today
//
  else if ( i == 58 )
  {
    jed = now_to_jed ( );
  }
//
//  End of Current Mayan Great Cycle
//  21 December 2012 CE/Gregorian
//
  else if ( i == 59 )
  {
    jed = transition_to_jed_mayan_long ( );
  }
//
//  Scaliger cycle repeats.
//  1 January 3266 CE/Gregorian
//
  else if ( i == 60 )
  {
    jed = transition_to_jed_jed ( );
  }
  else
  {
    jed = -1.0;
  }

  return jed;
}
//****************************************************************************80

void jed_to_cws_gps ( double jed, int &c, int &w, double &s )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_CWS_GPS converts a JED to a GPS CWS date.
//
//  Discussion:
//
//    The GPS time keeping is in terms of seconds, weeks, and cycles
//    of 1024 weeks.  The weeks and cycles begin numbering at 0.
//
//    The computation is only valid for dates after the GPS epoch,
//    that is, after 6 January 1980.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &C, &W, double &S, the 
//    corresponding GPS cycles/weeks/seconds date.
//
{
  double d;
  double jed_epoch;

  jed_epoch = epoch_to_jed_gps ( );

  d = jed - jed_epoch;

  if ( d < 0.0 )
  {
    s = -1.0;
    w = -1;
    c = -1;
    return;
  }

  w = ( int ) ( d ) / 7;
  d = d - ( double ) ( 7 * w );

  c = w / 1024;
  w = w - 1024 * c;

  s = d * ( double ) ( 24.0 * 60.0 * 60.0 );

  return;
}
//****************************************************************************80

double jed_to_datenum ( double jed )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_DATENUM converts a JED to a MATLAB date number.
//
//  Discussion:
//
//    The MATLAB "datenum" function accepts a string defining
//    a date and returns a datenumber:
//
//      dn = datenum ( 'Aug 17 1939' )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, double JED_TO_DATENUM, a MATLAB date number.
//
{
  double dn;

  dn = jed - 1721058.5;

  return dn;
}
//****************************************************************************80

void jed_to_mayan_long ( double jed, int &pictun, int &baktun, int &katun, 
  int &tun, int &uinal, int &kin, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_MAYAN_LONG converts a JED to a Mayan long count date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, chapter 27.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &PICTUN, &BAKTUN, &KATUN, &TUN, &UINAL, &KIN, values
//    defining the Mayan long date.
//
//    Output, double &F, the fractional part of the date.
//
{
  int days;
  int j;
  double jed_epoch;

  jed_epoch = epoch_to_jed_mayan_long ( );

  j = ( int ) ( jed - jed_epoch );
  f = ( jed - jed_epoch ) - ( double ) ( j );

  days = j;

  if ( 0 <= days )
  {
    pictun = days / 2880000;
    days = days - pictun * 2880000;
  }
  else
  {
    pictun = 0;
    while ( days < 0 )
    {
      pictun = pictun - 1;
      days = days + 2880000;
    }
  }

  baktun = days / 144000;
  days = days - baktun * 144000;
  katun = days / 7200;
  days = days - katun * 7200;
  tun = days / 360;
  days = days - tun * 360;
  uinal = days / 20;
  days = days - uinal * 20;
  kin = days;

  return;
}
//****************************************************************************80

void jed_to_mayan_round ( double jed, int &y, int &a, int &b, int &c, int &d, 
  double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_MAYAN_ROUND converts a JED to a Mayan round date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm K,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 340.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &A, &B, &C, &D, values defining the Mayan 
//    round date.
//
//    Output, double &F, the fractional part of the date.
//
{
  int days;
  int j;
  double jed_epoch;
  int n;

  jed_epoch = epoch_to_jed_mayan_long ( );

  j = ( int ) ( jed - jed_epoch );
  f = ( jed - jed_epoch ) - ( double ) ( j );

  days = j;

  y = 0;

  while ( days < 0 )
  {
    days = days + 18980;
    y = y - 1;
  }

  y = y + days / 18980;

  days = ( days % 18980 );

  a = i4_wrap ( days + 4, 1, 13 );
  b = i4_wrap ( days, 1, 20 );

  n = ( days + 348 ) % 365;
  c = ( n % 20 );
  d = n / 20;

  return;
}
//****************************************************************************80

double jed_to_mjd ( double jed )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_MJD converts a JED to a modified JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, double JED_TO_MJD, the modified Julian Ephemeris Date.
//
{
  double jed_epoch;
  double mjd;

  jed_epoch = epoch_to_jed_mjd ( );

  mjd = jed - jed_epoch;

  return mjd;
}
//****************************************************************************80

double jed_to_nearest_noon ( double jed1 )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_NEAREST_NOON converts a JED to the JED of the nearest noon.
//
//  Discussion:
//
//    This is primarily to make a fair test of the weekday routines,
//    which have trouble when the JED is at midnight.
//
//    Note that noon corresponds to an integral JED value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED1, the Julian Ephemeris Date.
//
//    Output, double JED_TO_NEAREST_NOON, the Julian Ephemeris Date
//    of the nearest noon.  If JED1 was at midnight, JED2 is
//    advanced to the NEXT noon, not the previous one.
//
{
  double jed2;

  jed2 = ( double ) r8_nint ( jed1 );

  return jed2;
}
//****************************************************************************80

double jed_to_next_noon ( double jed1 )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_NEXT_NOON converts a JED to the JED of the next noon.
//
//  Discussion:
//
//    This is primarily to make a fair test of the weekday routines,
//    which have trouble when the JED is at midnight.
//
//    Note that noon corresponds to an integral JED value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED1, the Julian Ephemeris Date.
//
//    Output, double JED_TO_NEXT_NOON, the Julian Ephemeris Date
//    of the next noon.
//
{
  double jed2;

  jed2 = ( double ) ( ( int ) ( jed1 ) );
//
//  The integer part of JED1 is one of the two integers that
//  bracket JED1.  If it's the smaller one (which it should
//  be as long as JED1 is positive), make it the bigger one.
//
//  This correctly leaves undisturbed cases where JED1 is
//  already an integer, and where JED1 is negative (which
//  is not a case we expect to occur often).
//
  if ( jed2 < jed1 )
  {
    jed2 = jed2 + 1.0;
  }

  return jed2;
}
//****************************************************************************80

double jed_to_rd ( double jed )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_RD converts a JED to an RD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz,
//    Calendrical Calculations, the Millennium Edition,
//    Cambridge, 2002.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, double JED_TO_RD, the RD date.
//
{
  double rd;
  double rd_epoch;

  rd_epoch = epoch_to_jed_rd ( );

  rd = jed - rd_epoch;

  return rd;
}
//*****************************************************************************80

double jed_to_ss_unix ( double jed )

//*****************************************************************************80
//
//  Purpose:
//
//    JED_TO_SS_UNIX converts a JED to a UNIX SS date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, double JED_TO_SS_UNIX, the corresponding UNIX SS date.
//
{
  double d;
  double jed_epoch;
  double s;

  jed_epoch = epoch_to_jed_unix ( );

  d = jed - jed_epoch;

  s = d * 24.0 * 60.0 * 60.0;

  return s;
}
//****************************************************************************80

void jed_to_weekday ( double jed, int &w, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_WEEKDAY computes the day of the week from a JED.
//
//  Discussion:
//
//    BC 4713/01/01 => JED = 0.0 was noon on a Monday.
//
//    jedmod = mod ( 0.0, 7.0 ) = 0.0
//    j = mod ( nint ( 0 ), 7 ) = 0
//    f = ( 0.0 + 0.5 ) - real ( j ) = 0.5
//    w = i4_wrap ( 0 + 2, 1, 7 ) = 2 = MONDAY
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int W, the day of the week of the date.
//    The days are numbered from Sunday through Saturday, 1 through 7.
//
//    Output, double F, the fractional part of the day.
//
{
  int j;
  double jedmod;

  jedmod = r8_mod ( jed, 7.0 );

  j = ( r8_nint ( jedmod ) % 7 );

  f = ( jedmod + 0.5 ) - ( double ) ( j );

  w = i4_wrap ( j + 2, 1, 7 );

  return;
}
//****************************************************************************80

int jed_to_year_hebrew ( double jed )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YEAR_HEBREW: the year in the Hebrew calendar when a JED occurred.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm H,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 331.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int JED_TO_YEAR_HEBREW, the year in the Hebrew calendar that 
//    included the JED.  If the input JED is less than the epoch of the Hebrew 
//    calendar, then Y is always returned as -1.
//
{
  double jed2;
  double jed_epoch;
  int m;
  int y;

  jed_epoch = epoch_to_jed_hebrew ( );

  if ( jed < jed_epoch )
  {
    y = - 1;
    return y;
  }
//
//  Using integer arithmetic in this computation may cause overflow.
//
//  Compute the number of months elapsed up to the date.
//
  m = 1 + ( int ) ( ( 25920.0 * ( jed - jed_epoch + 2.5 ) ) / 765433.0 );
//
//  Estimate the number of years represented by these months.
//
  y = 19 * ( m / 235 ) + ( 19 * ( i4_modp ( m, 235 ) - 2 ) ) / 235 + 1;
//
//  Determine the JED of the first day of that year.
//
  jed2 = new_year_to_jed_hebrew ( y );
//
//  We might have been off by 1 year.
//
  if ( jed < jed2 )
  {
    y = y - 1;
  }

  return y;
}
//****************************************************************************80

double jed_to_yearcount_bessel ( double jed )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YEARCOUNT_BESSEL converts a JED to Bessel year count.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, double BESSEL, the Bessel year.
//
{
  double bessel;
  double jed_epoch;
  double year_length = 365.242198781;

  jed_epoch = epoch_to_jed_bessel ( );
  bessel = 1900.0 + ( jed - jed_epoch ) / year_length;

  return bessel;
}
//****************************************************************************80

double jed_to_yearcount_julian ( double jed )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YEARCOUNT_JULIAN converts a JED to a Julian year count.
//
//  Discussion:
//
//    An average year in the Julian calendar is exactly 365.25 days long.
//    This calculation counts the number of average Julian years from
//    the beginning of the year 2000.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, double JED_TO_YEARCOUNT_JULIAN, the Julian year.
//
{
  double jed_epoch;
  double julian;
  double year_length = 365.25;

  jed_epoch = epoch_to_jed_y2k ( );

  julian = 2000.0 + ( jed - jed_epoch ) / year_length;

  return julian;
}
//****************************************************************************80

void jed_to_yjf_common ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_COMMON converts a JED to a Common YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_common ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_common ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_yjf_english ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_ENGLISH converts a JED to an English YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_english ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_english ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_yjf_gregorian ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_GREGORIAN converts a JED to a Gregorian YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_gregorian ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_gregorian ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_yjf_hebrew ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_HEBREW converts a JED to a Hebrew YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_hebrew ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_hebrew ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_yjf_islamic_a ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_ISLAMIC_A converts a JED to an Islamic-A YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_islamic_a ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_islamic ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_yjf_islamic_b ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_ISLAMIC_B converts a JED to an Islamic-B YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_islamic_b ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_islamic ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_yjf_julian ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_JULIAN converts a JED to a Julian YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_julian ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_julian ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_yjf_republican ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_REPUBLICAN converts a JED to a Republican YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_republican ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_republican ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_yjf_roman ( double jed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YJF_ROMAN converts a JED to a Roman YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  int d1;
  double f1;
  int m1;
  int y1;

  jed_to_ymdf_roman ( jed, y1, m1, d1, f1 );

  ymdf_to_yjf_roman ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void jed_to_ymdf_alexandrian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_ALEXANDRIAN converts a JED to an Alexandrian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, 
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 124;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime = t_prime / 30;
  d_prime = t_prime % 30;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 13 ) + 1;
  y = y_prime - 4690 + ( 13 - m ) / 13;

  return;
}
//****************************************************************************80

void jed_to_ymdf_armenian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_ARMENIAN converts a JED to an Armenian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 317;

  y_prime =   j_prime / 365;
  t_prime = ( j_prime % 365 );
  m_prime =   t_prime / 30;
  d_prime = ( t_prime % 30 );
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 13 ) + 1;
  y = y_prime - 5268 + ( 13 - m ) / 13;

  return;
}
//****************************************************************************80

void jed_to_ymdf_bahai ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_BAHAI converts a JED to a Bahai YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int g;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  g = ( 3 * ( ( 4 * j + 274273 ) / 146097 ) / 4 ) - 50;
  j_prime = j + 1412 + g;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime =       t_prime / 19;
  d_prime =       t_prime % 19;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( ( m_prime + 19 ) % 20 ) + 1;
  y = y_prime - 6560 + ( 39 - m ) / 20;

  return;
}
//****************************************************************************80

void jed_to_ymdf_common ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_COMMON converts a JED to a Common YMDF date.
//
//  Discussion:
//
//    The "common" calendar is meant to be the calendar which is Julian up to
//    JED = 2299160.5, and Gregorian thereafter.
//
//    There is no year 0.  BC years are specified using a negative value.
//
//  Example:
//
//        JED            Y    M   D
//    -------    ------------------
//          0    BCE  4713  Jan   1
//    2440000    CE   1968  May  23
//    2446065    CE   1984  Dec  31
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  double jed_transition;

  jed_transition = transition_to_jed_common ( );

  if ( jed <= jed_transition )
  {
    jed_to_ymdf_julian ( jed, y, m, d, f );
  }
  else
  {
    jed_to_ymdf_gregorian ( jed, y, m, d, f );
  }

  return;
}
//****************************************************************************80

void jed_to_ymdf_coptic ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_COPTIC converts a JED to a Coptic YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 124;

  y_prime =     ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime = t_prime / 30;
  d_prime = t_prime % 30;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 13 ) + 1;
  y = y_prime - 4996 + ( 13 - m ) / 13;

  return;
}
//****************************************************************************80

void jed_to_ymdf_eg_civil ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_EG_CIVIL converts a JED to an Egyptian Civil YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 47 - 1;

  y_prime = j_prime / 365;
  t_prime = j_prime % 365;
  m_prime = t_prime / 30;
  d_prime = t_prime % 30;
//
//  Convert the computational date to calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 13 ) + 1;
  y = y_prime - 3968 + ( 13 - m ) / 13;

  return;
}
//****************************************************************************80

void jed_to_ymdf_eg_lunar ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_EG_LUNAR converts a JED to an Egyptian Lunar YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int j;
  double jed_epoch;
  int ncycle;
 
  jed_epoch = epoch_to_jed_eg_lunar ( );

  j = ( int ) ( jed - jed_epoch );
  f = ( jed - jed_epoch ) - ( double ) ( j );

  d = 1 + j;
  m = 1;
  y = 1;
//
//  Account for the number of 25 year cycles of 9125 days.
//
  ncycle = d / 9125;
  y = y + 25 * ncycle;
  d = d - ncycle * 9125;

  while ( year_length_eg_lunar ( y ) < d )
  {
    d = d - year_length_eg_lunar ( y );
    y = y + 1;
  }

  while ( month_length_eg_lunar ( y, m ) < d )
  {
    d = d - month_length_eg_lunar ( y, m );
    m = m + 1;
  }
  return;
}
//****************************************************************************80

void jed_to_ymdf_english ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_ENGLISH converts a JED to an English YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  double jed_transition;

  jed_transition = transition_to_jed_english ( );

  if ( jed <= jed_transition )
  {
    jed_to_ymdf_julian ( jed, y, m, d, f );
  }
  else
  {
    jed_to_ymdf_gregorian ( jed, y, m, d, f );
  }
  return;
}
//****************************************************************************80

void jed_to_ymdf_ethiopian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_ETHIOPIAN converts a JED to an Ethiopian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 124;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime = t_prime / 30;
  d_prime = t_prime % 30;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 13 ) + 1;
  y = y_prime - 4720 + ( 13 - m ) / 13;

  return;
}
//****************************************************************************80

void jed_to_ymdf_gregorian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  PURPOSE:
//
//    JED_TO_YMDF_GREGORIAN converts a JED to a Gregorian YMDF date.
//
//  Discussion:
//
//    This Gregorian calendar is extended backwards in time before
//    its actual adoption.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int g;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  g = ( 3 * ( ( 4 * j + 274277 ) / 146097 ) / 4 ) - 38;
  j_prime = j + 1401 + g;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime =   ( 5 * t_prime + 2 ) / 153;
  d_prime = ( ( 5 * t_prime + 2 ) % 153 ) / 5;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( ( m_prime + 2 ) % 12 ) + 1;
  y = y_prime - 4716 + ( 14 - m ) / 12;
//
//  Any year before 1 AD must be moved one year further back, since
//  this calendar does not include a year 0.
//
  y = y_astronomical_to_common ( y );

  return;
}
//****************************************************************************80

void jed_to_ymdf_gregorian2 ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_GREGORIAN2 converts a JED to a Gregorian YMDF date.
//
//  Discussion:
//
//    The theory behind this routine is very clean.  The Gregorian
//    calendar has cycles of 1, 4, 100 and 400 years, and we can
//    analyze a date by determining where it lies within these cycles.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz, Stewart Clamen,
//    Calendrical Calculations, II: Three Historical Calendars,
//    Software - Practice and Experience,
//    Volume 23, Number 4, pages 383-404, April 1993.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  int d0;
  int d1;
  int d2;
  int d3;
  int d4;
  double f1;
  int g1 = 365;
  int g4 = 1461;
  int g100 = 36524;
  int g400 = 146097;
  int j;
  double jed_epoch;
  int j1;
  int n1;
  int n4;
  int n100;
  int n400;
  int y1;

  jed_epoch = epoch_to_jed_gregorian ( );

  j = ( int ) ( jed - jed_epoch );
  f1 = ( jed - jed_epoch ) - ( double ) ( j );

  d0 = j;
  n400 = 0;

  while ( d0 < 0 )
  {
    d0 = d0 + g400;
    n400 = n400 - 1;
  }

  n400 = n400 + d0 / g400;
  d1 = i4_modp ( d0, g400 );

  n100 = d1 / g100;
  d2 = i4_modp ( d1, g100 );

  n4 = d2 / g4;
  d3 = i4_modp ( d2, g4 );

  n1 = d3 / g1;
  d4 = i4_modp ( d3, g1 );

  if ( n100 == 4 || n1 == 4 )
  {
    j1 = 366;
    y1 = 400 * n400 + 100 * n100 + 4 * n4 + n1;
  }
  else
  {
    j1 = d4 + 1;
    y1 = 400 * n400 + 100 * n100 + 4 * n4 + n1 + 1;
  }
//
//  Any year before 1 AD must be moved one year further back, since
//  this calendar does not include a year 0.
//
  y1 = y_astronomical_to_common ( y1 );

  yjf_to_ymdf_gregorian ( y1, j1, f1, y, m, d, f );

  return;
}
//****************************************************************************80

void jed_to_ymdf_hebrew ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_HEBREW converts a JED to a Hebrew YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm I,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 334.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  double f1;
  int j1;
  double jed2;
  int type;
  int y1;

  y1 = jed_to_year_hebrew ( jed );

  jed2 = new_year_to_jed_hebrew ( y1 );

  type = year_to_type_hebrew ( y1 );

  j1 = ( int ) ( jed - jed2 );
  f1 = ( jed - jed2 ) - ( double ) ( j1 );

  j1 = j1 + 1;

  yjf_to_ymdf_hebrew ( y1, j1, f1, y, m, d, f );

  return;
}
//****************************************************************************80

void jed_to_ymdf_hindu_solar ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_HINDU_SOLAR converts a JED to a Hindu solar YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz, Stewart Clamen,
//    Calendrical Calculations, II: Three Historical Calendars,
//    Software - Practice and Experience,
//    Volume 23, Number 4, pages 383-404, April 1993.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  double jed_epoch;
  double jf;

//  Date = JF.
//
  jf = jed;
//
//  Date = JED_EPOCH + JF
//
  jed_epoch = epoch_to_jed_hindu_solar ( );
  jf = jf - jed_epoch;
//
//  Date = JED_EPOCH + Y years + JF
//
  y = ( int ) ( jf / year_length_hindu_solar ( ) );
  jf = jf - y * year_length_hindu_solar ( );
//
//  Date = JED_EPOCH + Y years + ( M - 1 ) months + JF
//
  m = 1 + ( int ) ( jf / month_length_hindu_solar ( ) );
  jf = jf - ( ( double ) ( m - 1 ) * month_length_hindu_solar ( ) );
//
//  Date = JED_EPOCH + Y years + ( M - 1 ) months + ( D - 1 ) days + f
//
  d = int ( jf ) + 1;
  f = jf - ( d - 1 );

  return;
}
//****************************************************************************80

void jed_to_ymdf_islamic_a ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_ISLAMIC_A converts a JED to an Islamic A YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 7665;

  y_prime =   ( 30 * j_prime + 15 ) / 10631;
  t_prime = ( ( 30 * j_prime + 15 ) % 10631 ) / 30;
  m_prime =   ( 100 * t_prime + 10 ) / 2951;
  d_prime = ( ( 100 * t_prime + 10 ) % 2951 ) / 100;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 12 ) + 1;
  y = y_prime - 5519 + ( 12 - m ) / 12;

  return;
}
//****************************************************************************80

void jed_to_ymdf_islamic_b ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_ISLAMIC_B converts a JED to an Islamic B YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 7664;

  y_prime =   ( 30 * j_prime + 15 ) / 10631;
  t_prime = ( ( 30 * j_prime + 15 ) % 10631 ) / 30;
  m_prime =   ( 100 * t_prime + 10 ) / 2951;
  d_prime = ( ( 100 * t_prime + 10 ) % 2951 ) / 100;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 12 ) + 1;
  y = y_prime - 5519 + ( 12 - m ) / 12;

  return;
}
//****************************************************************************80

void jed_to_ymdf_jelali ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_JELALI converts a JED to a Jelali YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  int day;
  int j;
  double jed_epoch;
  int n;

  jed_epoch = epoch_to_jed_jelali ( );

  j = ( int ) ( jed - jed_epoch );
  f = ( jed - jed_epoch ) - ( double ) ( j );

  d = 1 + j;
  m = 1;
  y = 1;
//
//  Account for the number of completed 4 year cycles of 1461 days.
//
  n = ( d - 1 ) / 1461;
  y = y + 4 * n;
  d = d - n * 1461;
//
//  Account for the number of completed 365 day years.
//
  n = ( d - 1 ) / 365;
  y = y + n;
  d = d - n * 365;
//
//  Account for the number of completed 30 day months.
//
  n = ( d - 1 ) / 30;
  m = m + n;
  d = d - n * 30;

  return;
}
//****************************************************************************80

void jed_to_ymdf_julian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_JULIAN converts a JED to a Julian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date (Y'/M'/D').
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 1401;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime =   ( 5 * t_prime + 2 ) / 153;
  d_prime = ( ( 5 * t_prime + 2 ) % 153 ) / 5;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( ( m_prime + 2 ) % 12 ) + 1;
  y = y_prime - 4716 + ( 14 - m ) / 12;
//
//  Any year before 1 AD must be moved one year further back, since
//  this calendar does not include a year 0.
//
  y = y_astronomical_to_common ( y );

  return;
}
//****************************************************************************80

void jed_to_ymdf_julian2 ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_JULIAN2 converts a JED to a Julian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  int j;
  int jd;
  int je;
  int jg;
//
//  Check the input.
//
  jed_check ( jed );

  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  jd = ( int ) ( ( ( double ) ( j + 1524 ) - 122.1 ) / 365.25 );
  je = ( int ) ( 365.25 * ( double ) ( jd ) );
  jg = ( int ) ( ( double ) ( j + 1524 - je ) / 30.6001 );
//
//  Now compute D, M and Y.
//
  d = j + 1524 - je - ( int ) ( 30.6001 * jg );

  if ( jg <= 13 )
  {
    m = jg - 1;
  }
  else
  {
    m = jg - 13;
  }

  if ( 2 < m )
  {
    y = jd - 4716;
  }
  else
  {
    y = jd - 4715;
  }
//
//  Any year before 1 AD must be moved one year further back, since
//  this calendar does not include a year 0.
//
  y = y_astronomical_to_common ( y );

  return;
}
//****************************************************************************80

void jed_to_ymdf_julian3 ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_JULIAN3 converts a JED to a Julian YMDF date.
//
//  Discussion:
//
//    The theory behind this routine is very clean.  The Julian
//    calendar has cycles of 1 and 4 years, and we can analyze a date
//    by determining where it lies within these cycles.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz, Stewart Clamen,
//    Calendrical Calculations, II: Three Historical Calendars,
//    Software - Practice and Experience,
//    Volume 23, Number 4, pages 383-404, April 1993.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int cycle_length = 1461;
  int d0;
  int d1;
  int d2;
  double f1;
  int j;
  int j1;
  double jed_epoch;
  int n1;
  int n4;
  int y1;
  int year_length = 365;

  jed_epoch = epoch_to_jed_julian ( );

  j = ( int ) ( jed - jed_epoch );
  f1 = ( jed - jed_epoch ) - ( double ) ( j );

  if ( f1 < 0.0 )
  {
    f1 = f1 + 1.0;
    j = j - 1;
  }

  d0 = j;
  n4 = 0;

  while ( d0 <= 0 )
  {
    d0 = d0 + cycle_length;
    n4 = n4 - 1;
  }

  n4 = n4 + d0 / cycle_length;
  d1 = i4_modp ( d0, cycle_length );
  n1 = d1 / year_length;
  d2 = i4_modp ( d1, year_length );

  if ( n1 == 4 )
  {
    j1 = 366;
    y1 = 4 * n4 + n1;
  }
  else
  {
    j1 = d2 + 1;
    y1 = 4 * n4 + n1 + 1;
  }
//
//  Any year before 1 AD must be moved one year further back, since
//  this calendar does not include a year 0.
//
  y1 = y_astronomical_to_common ( y1 );

  yjf_to_ymdf_julian ( y1, j1, f1, y, m, d, f );

  return;
}
//****************************************************************************80

void jed_to_ymdf_khwarizmian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_KHWARIZMIAN converts a JED to a Khwarizmian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date.
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 317;

  y_prime =   j_prime / 365;
  t_prime = ( j_prime % 365 );
  m_prime =   t_prime / 30;
  d_prime = ( t_prime % 30 );
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 13 ) + 1;
  y = y_prime - 5348 + ( 13 - m ) / 13;

  return;
}
//****************************************************************************80

void jed_to_ymdf_macedonian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_MACEDONIAN converts a JED to a Macedonian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &M, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date.
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 1401;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime =   ( 5 * t_prime + 2 ) / 153;
  d_prime = ( ( 5 * t_prime + 2 ) % 153 ) / 5;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( ( m_prime + 6 ) % 12 ) + 1;
  y = y_prime - 4405 + ( 18 - m ) / 12;

  return;
}
//****************************************************************************80

void jed_to_ymdf_persian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_PERSIAN converts a JED to a Persian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date.
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 77;

  y_prime =   j_prime / 365;
  t_prime = ( j_prime % 365 );
  m_prime =   t_prime / 30;
  d_prime = ( t_prime % 30 );
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( ( m_prime + 9 ) % 13 ) + 1;
  y = y_prime - 5348 + ( 22 - m ) / 13;

  return;
}
//****************************************************************************80

void jed_to_ymdf_republican ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_REPUBLICAN converts a JED to a Republican YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int g;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date.
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  g = ( 3 * ( ( 4 * j + 578797 ) / 146097 ) / 4 ) - 51;
  j_prime = j + 111 + g;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime =   t_prime / 30;
  d_prime = ( t_prime % 30 );
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( m_prime % 13 ) + 1;
  y = y_prime - 6504 + ( 13 - m ) / 13;

  return;
}
//****************************************************************************80

void jed_to_ymdf_roman ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_ROMAN converts a JED to a Roman YMDF date.
//
//  Discussion:
//
//    The Roman calendar used here is artificial.  It is assumed to begin
//    on the Julian calendar date 1 January 753 BC, and to be simply a
//    copy of the Julian calendar, shifted by 753 years.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int yj;

  jed_to_ymdf_julian ( jed, yj, m, d, f );

  y = y_julian_to_roman ( yj );

  return;
}
//****************************************************************************80

void jed_to_ymdf_saka ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_SAKA converts a JED to a Saka YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int g;
  int j;
  int j_prime;
  int m_prime;
  int s;
  int t_prime;
  int x;
  int y_prime;
  int z;
//
//  Determine the computational date.
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  g = ( 3 * ( ( 4 * j + 274073 ) / 146097 ) / 4 ) - 36;

  j_prime = j + 1348 + g;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;

  x = t_prime / 365;
  z = t_prime / 185 - x;
  s = 31 - z;

  m_prime =           ( t_prime - 5 * z ) / s;
  d_prime = 6 * x + ( ( t_prime - 5 * z ) % s );
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( ( m_prime + 1 ) % 12 ) + 1;
  y = y_prime - 4794 + ( 13 - m ) / 12;

  return;
}
//****************************************************************************80

void jed_to_ymdf_soor_san ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_SOOR_SAN converts a JED to a Soor San YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int j;
  double jed_epoch;
  int n;

  jed_epoch = epoch_to_jed_soor_san ( );

  j = ( int ) ( jed - jed_epoch );
  f = ( jed - jed_epoch ) - ( double ) ( j );

  d = 1 + j;
  m = 1;
  y = 1;
//
//  Account for the number of completed 4 year cycles of 1461 days.
//
  n = ( d - 1 ) / 1461;
  y = y + 4 * n;
  d = d - n * 1461;
//
//  Account for the number of completed 365 day years.
//
  n = ( d - 1 ) / 365;
  y = y + n;
  d = d - n * 365;
//
//  Account for the number of completed 30 day months.
//
  n = ( d - 1 ) / 30;
  m = m + n;
  d = d - n * 30;

  return;
}
//****************************************************************************80

void jed_to_ymdf_syrian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_SYRIAN converts a JED to a Syrian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm F,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 324-325.
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int d_prime;
  int j;
  int j_prime;
  int m_prime;
  int t_prime;
  int y_prime;
//
//  Determine the computational date.
//
  j = ( int ) ( jed + 0.5 );
  f = ( jed + 0.5 ) - ( double ) ( j );

  j_prime = j + 1401;

  y_prime =   ( 4 * j_prime + 3 ) / 1461;
  t_prime = ( ( 4 * j_prime + 3 ) % 1461 ) / 4;
  m_prime =   ( 5 * t_prime + 2 ) / 153;
  d_prime = ( ( 5 * t_prime + 2 ) % 153 ) / 5;
//
//  Convert the computational date to a calendar date.
//
  d = d_prime + 1;
  m = ( ( m_prime + 5 ) % 12 ) + 1;
  y = y_prime - 4405 + ( 17 - m ) / 12;

  return;
}
//****************************************************************************80

void jed_to_ymdf_zoroastrian ( double jed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_TO_YMDF_ZOROASTRIAN converts a JED to a Zoroastrian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, double &F,
//    the YMDF date.
//
{
  int j;
  double jed_epoch;
  int months;
  int years;

  jed_epoch = epoch_to_jed_zoroastrian ( );

  j = ( int ) ( jed - jed_epoch );
  f = ( jed - jed_epoch ) - ( double ) ( j );

  d = 1 + j;
  m = 1;
  y = 1;

  years = ( d - 1 ) / 365;
  y = y + years;
  d = d - years * 365;

  months = ( d - 1 ) / 30;
  m = m + months;
  d = d - months * 30;

  return;
}
//****************************************************************************80

double mayan_long_to_jed ( int pictun, int baktun, int katun, int tun, 
  int uinal, int kin, double f )

//****************************************************************************80
//
//  Purpose:
//
//    MAYAN_LONG_TO_JED converts a Mayan long count date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm L,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 341.
//
//  Parameters:
//
//    Input, int PICTUN, BAKTUN, KATUN, TUN, UINAL, KIN, values
//    defining the Mayan long date.
//
//    Input, double F, the fractional part of the date.
//
//    Output, double MAYAN_LONG_TO_JED, the Julian Ephemeris Date.
//
{
  int days;
  double jed;
  double jed_epoch;

  days = (((((   pictun   * 20 + baktun ) * 20 + katun  ) * 20 
    + tun    ) * 18 + uinal  ) * 20 + kin );

  jed_epoch = epoch_to_jed_mayan_long ( );

  jed = jed_epoch + ( double ) ( days ) + f;

  return jed;
}
//****************************************************************************80

double mayan_round_to_jed ( int y, int a, int b, int c, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    MAYAN_ROUND_TO_JED converts a Mayan round date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm L,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 341.
//
//  Parameters:
//
//    Input, int Y, A, B, C, D, values defining the Mayan 
//    round date.
//
//    Input, double F, the fractional part of the date.
//
//    Output, double MAYAN_ROUND_TO_JED, the Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch;
  int m;
  int n;
  int r;

  m = 13 * i4_modp ( 60 + 3 * ( a - b ), 20 ) + a - 1;
  m = i4_modp ( m + 101, 260 );
  n = 20 * d + c;
  n = i4_modp ( n + 17, 365 );
  r = 365 * i4_modp ( 364 + m - n, 52 ) + n;

  jed_epoch = epoch_to_jed_mayan_long ( );
  jed = jed_epoch + ( double ) ( 18980 * y + r ) + f;

  return jed;
}
//****************************************************************************80

void minute_borrow_common ( int &y, int &m, int &d, int &h, int &n )

//****************************************************************************80
//
//  Purpose:
//
//    MINUTE_BORROW_COMMON "borrows" an hour of minutes in a Common date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, &H, &N, the year,
//    month, day, hour and minute representing a date.  On input, N
//    might be negative.
//    On output, H should have decreased by one, and N gone up by 60.
//
{
  while ( n < 0 )
  {
    n = n + 60;
    h = h - 1;

    hour_borrow_common ( y, m, d, h );
  }
  return;
}
//****************************************************************************80

void minute_carry_common ( int &y, int &m, int &d, int &h, int &n )

//****************************************************************************80
//
//  Purpose:
//
//    MINUTE_CARRY_COMMON: given a Common YMDHMS date, carries minutes to hours.
//
//  Algorithm:
//
//    While 60 <= N:
//
//      decrease N by the number of minutes in an hour;
//      increase H by 1;
//      if necessary, adjust Y, M and D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, &H, &N, the date.
//    On output, N is between 0 and 59.
//
{
  while ( 60 <= n )
  {
    n = n - 60;
    h = h + 1;

    hour_carry_common ( y, m, d, h );
  }

  return;
}
//****************************************************************************80

double mjd_to_jed ( double mjd )

//****************************************************************************80
//
//  Purpose:
//
//    MJD_TO_JED converts a modified JED to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MJD, the modified Julian Ephemeris Date.
//
//    Output, double MJD_TO_JED, the Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch;

  jed_epoch = epoch_to_jed_mjd ( );
  jed = mjd + jed_epoch;

  return jed;
}
//****************************************************************************80

void month_borrow_alexandrian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_ALEXANDRIAN borrows a year of months on Alexandrian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//    On input, M might be negative.  On output, Y should have decreased by
//    one, and M gone up by the number of months in the year that we
//    "cashed in".  The routine knows there was no year 0.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_alexandrian ( y );

    m = m + months;
    y = y - 1;
  }
  return;
}
//****************************************************************************80

void month_borrow_bahai ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_BAHAI borrows a year of months on the Bahai calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//    On input, M might be negative.  On output, Y should have decreased by
//    one, and M gone up by the number of months in the year that we
//    "cashed in".  The routine knows there was no year 0.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_bahai ( y );

    m = m + months;
    y = y - 1;
  }
  return;
}
//****************************************************************************80

void month_borrow_common ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_COMMON borrows a year of months on the Common calendar.
//
//  Discussion:
//
//    If the month index is legal, nothing is done.  If the month index
//    is too small, then one or more years are "cashed in" to bring the
//    month index up to a legal value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_common ( y );

    m = m + months;
    y = y - 1;

    if ( y == 0 )
    {
      y = - 1;
    }
  }
  return;
}
//****************************************************************************80

void month_borrow_eg_civil ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_EG_CIVIL borrows a year of months on Egyptian Civil calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_eg_civil ( y );

    m = m + months;
    y = y - 1;
  }
  return;
}
//****************************************************************************80

void month_borrow_english ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_ENGLISH borrows a year of months on the English calendar.
//
//  Discussion:
//
//    If the month index is legal, nothing is done.  If the month index
//    is too small, then one or more years are "cashed in" to bring the
//    month index up to a legal value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_english ( y );

    m = m + months;
    y = y - 1;

    if ( y == 0 )
    {
      y = - 1;
    }
  }
  return;
}
//****************************************************************************80

void month_borrow_gregorian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_GREGORIAN borrows a year of months on the Gregorian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_gregorian ( y );

    m = m + months;
    y = y - 1;

    if ( y == 0 )
    {
      y = - 1;
    }
  }
  return;
}
//****************************************************************************80

void month_borrow_hebrew ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_HEBREW borrows a year of months on the Hebrew calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_hebrew ( y );

    m = m + months;
    y = y - 1;
  }
  return;
}
//****************************************************************************80

void month_borrow_islamic ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_ISLAMIC borrows a year of months on the Islamic calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_islamic ( y );

    m = m + months;
    y = y - 1;
  }
  return;
}
//****************************************************************************80

void month_borrow_julian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_JULIAN borrows a year of months on the Julian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_julian ( y );

    m = m + months;
    y = y - 1;

    if ( y == 0 )
    {
      y = - 1;
    }
  }
  return;
}
//****************************************************************************80

void month_borrow_republican ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_REPUBLICAN borrows a year of months on the Republican calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_republican ( y );

    m = m + months;
    y = y - 1;
  }
  return;
}
//****************************************************************************80

void month_borrow_roman ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_BORROW_ROMAN borrows a year of months on the Roman calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  int months;

  while ( m <= 0 )
  {
    months = year_length_months_roman ( y );

    m = m + months;
    y = y - 1;

    if ( y == 0 )
    {
      y = - 1;
    }
  }
  return;
}
//****************************************************************************80

void month_carry_alexandrian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_ALEXANDRIAN carries a year of months on the Alexandrian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the year and month.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_alexandrian ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_bahai ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_BAHAI carries a year of months on the Bahai calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the year and month.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_bahai ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_common ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_COMMON carries a year of months on the Common calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the year and month.
//    On output, M is no greater than 12.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_common ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_eg_civil ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_EG_CIVIL carries a year of months on the Egyptian Civil calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the year and month.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_eg_civil ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_english ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_ENGLISH carries a year of months on the English calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the year and month.
//    On output, M is no greater than 12.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_english ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_gregorian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_GREGORIAN carries a year of months on the Gregorian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the year and month.
//    On output, M is no greater than 12.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_gregorian ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_hebrew ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_HEBREW carries a year of months on the Hebrew calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the year and month.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_hebrew ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_islamic ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_ISLAMIC carries a year of months on the Islamic calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the year and month.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_islamic ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_julian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_JULIAN carries a year of months on the Julian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the year and month.
//    On output, M is no greater than 12.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_julian ( y );

    if ( m <= months )
    {
      break;
    }

    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_republican ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_REPUBLICAN carries a year of months on the Republican calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the year and month.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_republican ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

void month_carry_roman ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_CARRY_ROMAN carries a year of months on the Roman calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the year and month.
//
{
  int months;

  for ( ; ; )
  {
    months = year_length_months_roman ( y );

    if ( m <= months )
    {
      break;
    }
    m = m - months;
    y = y + 1;
  }
  return;
}
//****************************************************************************80

int month_length_alexandrian ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_ALEXANDRIAN returns the number of days in an Alexandrian month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_ALEXANDRIAN, the number of 
//    days in the month.
//
{
  int days;
  int m2;
  int mdays[13] = {
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;

  if ( m2 < 1 || 13 < m2 )
  {
    days = 0;
  }
  else
  {
    days = mdays[m2-1];
  }

  if ( m2 == 13 && year_is_leap_alexandrian ( y2 ) )
  {
    days = days + 1;
  }

  return days;
}
//****************************************************************************80

int month_length_bahai ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_BAHAI returns the number of days in a Bahai month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_BAHAI, the number of 
//    days in the month.
//
{
  int days;
  int m2;
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Check the input.
//
  ym_check_bahai ( y2, m2 );

  if ( m2 <= 18 || m2 == 20 )
  {
    days = 19;
  }
  else if ( year_is_leap_bahai ( y2 ) )
  {
    days = 5;
  }
  else
  {
    days = 4;
  }

  return days;
}
//****************************************************************************80

int month_length_common ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_COMMON returns the number of days in a Common month.
//
//  Discussion:
//
//    The "common" calendar is meant to be the calendar which is Julian up to
//    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
//
//    The routine knows that February has 28 days, except in leap years,
//    when it has 29.
//
//    In the Common calendar, October 1582 had only 21 days
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_COMMON, the number of days 
//    in the month.
//
{
  int days;
  int m2;
  int mdays[12] = {
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Check the input.
//
  ym_check_common ( y2, m2 );
//
//  Take care of the special case.
//
  if ( y2 == 1582 )
  {
    if ( m2 == 10 )
    {
      days = 21;
      return days;
    }
  }
//
//  Get the number of days in the month.
//
  days = mdays[m2-1];
//
//  If necessary, add 1 day for February 29.
//
  if ( m2 == 2 && year_is_leap_common ( y2 ) )
  {
    days = days + 1;
  }
  return days;
}
//****************************************************************************80

int month_length_coptic ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_COPTIC returns the number of days in a Coptic month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_COPTIC, the number of days 
//    in the month.
//
{
  int days;
  int m2;
  int mdays[13] = {
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;

  if ( m2 < 1 || 13 < m2 )
  {
    days = 0;
  }
  else
  {
    days = mdays[m2-1];
  }

  if ( m2 == 13 && year_is_leap_coptic ( y2 ) )
  {
    days = days + 1;
  }
  return days;
}
//****************************************************************************80

int month_length_eg_civil ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_EG_CIVIL returns the number of days in an Egyptian Civil month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_EG_CIVIL, the number of days 
//    in the month.
//
{
  int days;

  if ( m < 1 )
  {
    days = 0;
  }
  else if ( m <= 12 )
  {
    days = 30;
  }
  else if ( m == 13 )
  {
    days = 5;
  }
  else
  {
    days = 0;
  }
  return days;
}
//****************************************************************************80

int month_length_eg_lunar ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_EG_LUNAR returns the number of days in an Egyptian Lunar month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_EG_LUNAR, the number of days 
//    in the month.
//
{
  int days;
  int last;
  int m2;
  int mdays[13] = {
    29, 30, 29, 30, 29, 30, 29, 30, 29, 30, 29, 30, 30 };

  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;

  last = year_length_months_eg_lunar ( y2 );

  if ( m2 < 1 )
  {
    days = 0;
  }
  else if ( m2 <= last )
  {
    days = mdays[m2-1];
  }
  else
  {
    days = 0;
  }

  if ( m2 == last )
  {
    if ( year_is_leap_eg_lunar ( y ) )
    {
      days = days + 1;
    }
  }

  return days;
}
//****************************************************************************80

int month_length_english ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_ENGLISH returns the number of days in an English month.
//
//  Discussion:
//
//    In the English calendar, September 1752 had only 19 days.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_ENGLISH, the number of days in the month.
//
{
  int days;
  int m2;
  int mdays[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int y2;
//
//  Check the input.
//
  m2 = m;
  y2 = y;

  ym_check_english ( y2, m2 );
//
//  Take care of special cases:
//
  if ( y2 == 1752 )
  {
    if ( m2 == 9 )
    {
      days = 19;
      return days;
    }
  }
//
//  Get the number of days in the month.
//
  days = mdays[m2-1];
//
//  If necessary, add 1 day for February 29.
//
  if ( m2 == 2 && year_is_leap_english ( y2 ) )
  {
    days = days + 1;
  }

  return days;
}
//****************************************************************************80

int month_length_ethiopian ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_ETHIOPIAN returns the number of days in an Ethiopian month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_ETHIOPIAN, the number of days 
//    in the month.
//
{
  int days;
  int m2;
  int mdays[13] = {
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;

  if ( m2 < 1 || 13 < m2 )
  {
    days = 0;
  }
  else
  {
    days = mdays[m2-1];
  }

  if ( m2 == 13 && year_is_leap_ethiopian ( y2 ) )
  {
    days = days + 1;
  }

  return days;
}
//****************************************************************************80

int month_length_greek ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_GREEK returns the number of days in a Greek month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_GREEK, the number of days 
//    in the month.
//
{
  int days;
  int m2;
  int mdays[13] = {
    30, 29, 30, 29, 30, 29, 29, 30, 29, 30, 29, 30, 29 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;

  if ( m2 < 1 )
  {
    days = 0;
    return days;
  }
//
//  A 13-month year.
//
  if ( year_is_embolismic_greek ( y2 ) )
  {
    if ( 13 < m2 )
    {
      days = 0;
      return days;
    }

    days = mdays[m2-1];

    if ( m2 == 7 && year_is_leap_greek ( y2 ) )
    {
      days = days + 1;
    }
  }
//
//  A 12 month year.
//
  else
  {
    if ( m2 <= 6 )
    {
      days = mdays[m2-1];
    }
    else if ( m2 <= 12 )
    {
      days = mdays[m2];
    }
    else
    {
      days = 0;
    }
  }

  return days;
}
//****************************************************************************80

int month_length_gregorian ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_GREGORIAN returns the number of days in a Gregorian month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_GREGORIAN, the number of days in the month.
//
{
  int days;
  int mdays[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
//
//  Check the input.
//
  ym_check_gregorian ( y, m );
//
//  Get the number of days in the month.
//
  days = mdays[m-1];
//
//  If necessary, add 1 day for February 29.
//
  if ( m == 2 )
  {
    if ( year_is_leap_gregorian ( y ) )
    {
      days = days + 1;
    }
  }
  return days;
}
//****************************************************************************80

int month_length_hebrew ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_HEBREW returns the number of days in a Hebrew month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 333.
//
//  Parameters:
//
//    Input, int Y, M, the year and month number.  Note that
//    some years only had 12 months in them, while others have 13.  If
//    Y only has 12 months, then the length of the 13th month is
//    returned as 0 days.
//
//    Output, int MONTH_LENGTH_HEBREW, the number of days 
//    in the month.
//
{
  int a[6*13] = {
    30,  30,  30,  30,  30,  30, 
    29,  29,  30,  29,  29,  30, 
    29,  30,  30,  29,  30,  30, 
    29,  29,  29,  29,  29,  29, 
    30,  30,  30,  30,  30,  30, 
    29,  29,  29,  30,  30,  30, 
    30,  30,  30,  29,  29,  30, 
    29,  29,  29,  30,  30,  29, 
    30,  30,  30,  29,  29,  29, 
    29,  29,  29,  30,  30,  30, 
    30,  30,  30,  29,  29,  29, 
    29,  29,  29,  30,  30,  30, 
     0,   0,   0,  29,  29,  29 };
  int days;
  int m2;
  int type;
  int y2;
//
//  Copy the input
//
  y2 = y;
  m2 = m;
//
//  Check the input.
//
  ym_check_hebrew ( y2, m2 );

  type = year_to_type_hebrew ( y2 );

  if ( type < 1 || 6 < type )
  {
    cerr << "\n";
    cerr << "MONTH_LENGTH_HEBREW - Fatal error!\n";
    cerr << "  Illegal year TYPE = " << type << "\n";
    cerr << "  Y = " << y2 << "\n";
    exit ( 1 );
  }

  if ( m2 < 1 || 13 < m2 )
  {
    cerr << "\n";
    cerr << "MONTH_LENGTH_HEBREW - Fatal error!\n";
    cerr << "  Illegal MONTH = " << m2 << "\n";
    exit ( 1 );
  }

  days = a[(type-1)+(m2-1)*6];

  return days;
}
//****************************************************************************80

double month_length_hindu_solar ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_HINDU_SOLAR returns the number of days in a Hindu solar month.
//
//  Discussion:
//
//    Warning: this is a DOUBLE PRECISION quantity, with a fractional part!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double MONTH_LENGTH_HINDU_SOLAR, the number of
//    days in the month.
//
{
  double days;

  days = year_length_hindu_solar ( ) / 12.0;

  return days;
}
//****************************************************************************80

int month_length_iranian ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_IRANIAN returns the number of days in an Iranian month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_IRANIAN, the number of 
//    days in the month.
//
{
  int days;
  int m2;
  int mdays[12] = {
    31, 31, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Get the number of days in the month.
//
  days = mdays[m2-1];
//
//  If necessary, add 1 day for a leap year..
//
  if ( m2 == 12 && year_is_leap_iranian ( y2 ) )
  {
    days = days + 1;
  }

  return days;
}
//****************************************************************************80

int month_length_islamic ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_ISLAMIC returns the number of days in an Islamic month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_ISLAMIC, the number of days 
//    in the month.
//
{
  int days;
  int m2;
  int mdays[12] = {
    30, 29, 30, 29, 30, 29, 30, 29, 30, 29, 30, 29 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Check the input.
//
  ym_check_islamic ( y2, m2 );
//
//  Get the number of days in the month.
//
  days = mdays[m2-1];
//
//  If necessary, add 1 day for a leap year.
//
  if ( m2 == 12 && year_is_leap_islamic ( y2 ) )
  {
    days = days + 1;
  }
  return days;
}
//****************************************************************************80

int month_length_julian ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_JULIAN returns the number of days in a Julian month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_JULIAN, the number of days in the month.
//
{
  int days;
  int m2;
  int mdays[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int y2;
//
//  Check the input.
//
  m2 = m;
  y2 = y;

  ym_check_julian ( y2, m2 );
//
//  Get the number of days in the month.
//
  days = mdays[m2-1];
//
//  If necessary, add 1 day for February 29.
//
  if ( m2 == 2 && year_is_leap_julian ( y2 ) )
  {
    days = days + 1;
  }
  return days;
}
//****************************************************************************80

double month_length_lunar ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_LUNAR returns the number of days in a lunar month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, double MONTH_LENGTH_LUNAR, the number of days in 
//    the month.
//
{
  int days;

  if ( m < 1 || 12 < m )
  {
    days = 0;
  }
  else
  {
    days = 29.53058;
  }

  return days;
}
//****************************************************************************80

int month_length_persian ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_PERSIAN returns the number of days in a Persian month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_PERSIAN, the number of days 
//    in the month.
//
{
  int days;
  int m2;
  int mdays[12] = {
    31, 31, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Get the number of days in the month.
//
  days = mdays[m2-1];
//
//  If necessary, add 1 day for a leap year.
//
  if ( m2 == 12 && year_is_leap_persian ( y2 ) )
  {
    days = days + 1;
  }

  return days;
}
//****************************************************************************80

int month_length_republican ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_REPUBLICAN returns the number of days in a Republican month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_REPUBLICAN, the number of days 
//    in the month.
//
{
  int days;
  int m2;
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Check the input.
//
  ym_check_republican ( y2, m2 );
//
//  Get the number of days in the month.
//
  if ( 1 <= m2 && m2 <= 12 )
  {
    days = 30;
  }
  else if ( m2 == 13 )
  {
    if ( year_is_leap_republican ( y2 ) )
    {
      days = 6;
    }
    else
    {
      days = 5;
    }
  }
  return days;
}
//****************************************************************************80

int month_length_roman ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_ROMAN returns the number of days in a Roman month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year in which the month occurred.
//
//    Input, int M, the number of the month.
//
//    Output, int MONTH_LENGTH_ROMAN, the number of days 
//    in the month.
//
{
  int days;
  int m2;
  int mdays[12] ={
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int y2;
//
//  Copy the input.
//
  m2 = m;
  y2 = y;
//
//  Check the input.
//
  ym_check_roman ( y2, m2 );

  days = mdays[m2-1];

  if ( m2 == 2 && year_is_leap_roman ( y2 ) )
  {
    days = days + 1;
  }
  return days;
}
//****************************************************************************80

double month_length_synodic ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH_SYNODIC returns the mean synodic month length.
//
//  Discussion:
//
//    The synodic month is the time from one new moon to the next, that is,
//    when the moon and Sun are in conjunction.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double MONTH_LENGTH_SYNODIC, the length of the 
//    mean synodic month,
//    in days.
//
{
  double days;

  days = 29.53058885;

  return days;
}
//****************************************************************************80

string month_to_month_name_common ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_TO_MONTH_NAME_COMMON returns the name of a Common month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the month index.
//
//    Output, string MONTH_TO_MONTH_NAME_COMMON, the month name.
//
{
  string month_name;
  string name[12] = {
    "January", 
    "February", 
    "March", 
    "April", 
    "May", 
    "June", 
    "July", 
    "August", 
    "September", 
    "October", 
    "November", 
    "December" };

  if ( m < 1 || 12 < m )
  {
    cerr << "\n";
    cerr << "MONTH_TO_MONTH_NAME_COMMON - Fatal error!\n";
    cerr << "  Illegal month index.\n";
    exit ( 1 );
  }
  month_name = name[m-1];

  return month_name;
}
//****************************************************************************80

string month_to_month_name_common3 ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_TO_MONTH_NAME_COMMON3 returns the abbreviated name of a Common month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the month index.
//
//    Output, string MONTH_TO_MONTH_NAME_COMMON3, the month name.
//
{
  string month_name;
  string name[12] = {
    "Jan", 
    "Feb", 
    "Mar", 
    "Apr", 
    "May", 
    "Jun", 
    "Jul", 
    "Aug", 
    "Sep", 
    "Oct", 
    "Nov", 
    "Dec" };

  if ( m < 1 || 12 < m )
  {
    cerr << "\n";
    cerr << "MONTH_TO_MONTH_NAME_COMMON3 - Fatal error!\n";
    cerr << "  Illegal month index.\n";
    exit ( 1 );
  }
  month_name = name[m-1];

  return month_name;
}
//****************************************************************************80

int month_to_nones_roman ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_TO_NONES_ROMAN returns the day of the nones of a Roman month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the month index.
//
//    Output, int MONTH_TO_NONES_ROMAN, the day of the nones of the month.
//
{
  int d;
  int nones[12] = {
    5, 5, 7, 5, 7, 5, 7, 5, 5, 7, 5, 5 };

  if ( m < 1 || 12 < m )
  {
    d = -1;
  }
  else
  {
    d = nones[m-1];
  }
  return d;
}
//****************************************************************************80

void mothers_day ( int y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    MOTHERS_DAY computes the date of Mother's Day (US) for a Common year.
//
//  Discussion:
//
//    Mother's Day occurs on the second Sunday in May.
//
//  Example:
//
//    Input:
//
//      Y = 2003
//
//    Output:
//
//      M = 5
//      D = 11
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &M, &D, the month and day of Mother's Day.
//
{
  double f;
  int w;
//
//  Determine the day of the week for 8 May, the earliest day
//  that Mother's day can occur.  
//
  m = 5;
  d = 8;
  f = 0.0;

  w = ymdf_to_weekday_common ( y, m, d, f );
//
//  W = 1 means this day is Sunday, and day D is Mother's day.
//  Otherwise, figure out how to increment W to 8 (Sunday again);
//  The same increment makes D the correct day number.
//
  if ( w != 1 )
  {
    d = d + 8 - w;
  }
  return;
}
//****************************************************************************80

double new_year_to_jed_hebrew ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    NEW_YEAR_TO_JED_HEBREW returns the JED of the beginning of a Hebrew year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm G,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 330.
//
//  Parameters:
//
//    Input, int Y, the Hebrew year.
//
//    Output, double NEW_YEAR_TO_JED_HEBREW, the Julian Ephemeris Date.
//
{
  int d;
  int e;
  int e_prime;
  int f;
  int g;
  double jed;
  double jed_epoch;
  int mu;
  int t_prime;
  int tc;
  int th;
  int w;

  mu = ( 235 * y - 234 ) / 19;
  tc = 204 + 793 * mu;
  th = 5 + 12 * mu + tc / 1080;
  d = 1 + 29 * mu + th / 24;
  t_prime = ( tc % 1080 ) + 1080 * ( th % 24 );

  w = i4_wrap ( d + 1, 1, 7 );

  e = ( ( 7 * y + 13 ) % 19 ) / 12;
  e_prime = ( ( 7 * y + 6 ) % 19 ) / 12;

  if ( 19940 <= t_prime ||
     (  9924 <= t_prime && w == 3 && e == 0 ) ||
     ( 16788 <= t_prime && w == 2 && e == 0 && e_prime == 1 ) )
  {
    d = d + 1;
  }

  jed_epoch = epoch_to_jed_hebrew ( );

  f = ( d + 5 ) % 7;
  g = f % 2;

  jed = jed_epoch - 1.0 + ( double ) ( d + g );

  return jed;
}
//****************************************************************************80

double now_to_jed ( )

//****************************************************************************80
//
//  Purpose:
//
//    NOW_TO_JED expresses the current date as JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double JED, the Julian Ephemeris Date.
//
{
  time_t clock;
  int d;
  double f;
  int h;
  double jed;
  struct tm *lt;
  int m;
  int mu;
  int n;
  int s;
  time_t tloc;
  int y;

  clock = time ( &tloc );
  lt = localtime ( &clock );

  y = lt->tm_year + 1900;
  m = lt->tm_mon + 1;
  d = lt->tm_mday + 1;
  h = lt->tm_hour;
  n = lt->tm_min;
  s = lt->tm_sec;
  mu = 0;

  f = ( double ) mu;
  f = ( double ) ( s ) + f / 1000.0;
  f = ( double ) ( n ) + f / 60.0;
  f = ( double ) ( h ) + f / 60.0;
  f = f / 24.0;

  jed = ymdf_to_jed_common ( y, m, d, f );

  return jed;
}
//****************************************************************************80

void now_to_yjf_common ( int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    NOW_TO_YJF_COMMON expresses the current date as a Common YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int &Y, &J, double &F, the YJF date.
//
{
  time_t clock;
  int d1;
  double f1;
  int h1;
  struct tm *lt;
  int m1;
  int mu1;
  int n1;
  int s1;
  time_t tloc;
  int y1;

  clock = time ( &tloc );
  lt = localtime ( &clock );

  y1 = lt->tm_year + 1900;
  m1 = lt->tm_mon + 1;
  d1 = lt->tm_mday + 1;
  h1 = lt->tm_hour;
  n1 = lt->tm_min;
  s1 = lt->tm_sec;
  mu1 = 0;

  f1 = ( double ) mu1;
  f1 = ( double ) ( s1 ) + f1 / 1000.0;
  f1 = ( double ) ( n1 ) + f1 / 60.0;
  f1 = ( double ) ( h1 ) + f1 / 60.0;
  f1 = f1 / 24.0;

  ymdf_to_yjf_common ( y1, m1, d1, f1, y, j, f );

  return;
}
//****************************************************************************80

void now_to_ymdf_common ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    NOW_TO_YMDF_COMMON expresses the current date as a Common YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int &Y, &M, &D, double &F, the YMDF date.
//
{
  time_t clock;
  int h;
  struct tm *lt;
  int mu;
  int n;
  int s;
  time_t tloc;

  clock = time ( &tloc );
  lt = localtime ( &clock );

  y = lt->tm_year + 1900;
  m = lt->tm_mon + 1;
  d = lt->tm_mday + 1;
  h = lt->tm_hour;
  n = lt->tm_min;
  s = lt->tm_sec;
  mu = 0;

  f = ( double ) mu;
  f = ( double ) ( s ) + f / 1000.0;
  f = ( double ) ( n ) + f / 60.0;
  f = ( double ) ( h ) + f / 60.0;
  f = f / 24.0;

  return;
}
//****************************************************************************80

void now_to_ymdhms_common ( int &y, int &m, int &d, int &h, int &n, int &s )

//****************************************************************************80
//
//  Purpose:
//
//    NOW_TO_YMDHMS_COMMON expresses the current date as a Common YMDHMS date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int &Y, &M, &D, &H, &N, &S, the YMDHMS date.
//
{
  time_t clock;
  struct tm *lt;
  time_t tloc;

  clock = time ( &tloc );
  lt = localtime ( &clock );

  y = lt->tm_year + 1900;
  m = lt->tm_mon + 1;
  d = lt->tm_mday + 1;
  h = lt->tm_hour;
  n = lt->tm_min;
  s = lt->tm_sec;

  return;
}
//****************************************************************************80

double nyt_to_jed ( int volume, int issue )

//****************************************************************************80
//
//  Purpose:
//
//    NYT_TO_JED converts an NYT date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int VOLUME, ISSUE, the New York Times 
//    volume and issue.
//
//    Output, double NYT_TO_JED, the Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch_50000 = 2449790.5;

  if ( 149 < volume )
  {
    jed = jed_epoch_50000 + ( double ) ( issue - 50000 + 500 );
  }
//
//  Take care of the bizarre case of the second half of Volume 149,
//  Jan 1 2000 to Sep 17 2000, issues 51254 through ?, which were also
//  lowered by 500.
//
  else if ( volume == 149 && issue < 51600 )
  {
    jed = jed_epoch_50000 + ( double ) ( issue - 50000 + 500 );
  }
  else if ( 44028 <= issue )
  {
    jed = jed_epoch_50000 + ( double ) ( issue - 50000 );
  }
//
//  Factor in the strike of 1978.
//
  else
  {
    jed = jed_epoch_50000 + ( double ) ( issue - 50000 - 88 );
  }
  return jed;
}
//****************************************************************************80

void nyt_to_ymd ( int volume, int issue, int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    NYT_TO_YMD converts an NYT date to a YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int VOLUME, ISSUE, the New York Times 
//    volume and issue.
//
//    Output, int &Y, &M, &D, the year, month and day.
//
{
  double f;
  double jed;

  jed = nyt_to_jed ( volume, issue );

  jed_to_ymdf_common ( jed, y, m, d, f );

  return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
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
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_mod ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOD returns the remainder of R8 division.
//
//  Discussion:
//
//    If
//      REM = R8_MOD ( X, Y )
//      RMULT = ( X - REM ) / Y
//    then
//      X = Y * RMULT + REM
//    where REM has the same sign as X, and abs ( REM ) < Y.
//
//  Example:
//
//        X         Y     R8_MOD   R8_MOD  Factorization
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
//    Input, double X, the number to be divided.
//
//    Input, double Y, the number that divides X.
//
//    Output, double R8_MOD, the remainder when X is divided by Y.
//
{
  double value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_MOD - Fatal error!\n";
    cerr << "  R8_MOD ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( x < 0.0 && 0.0 < value )
  {
    value = value - r8_abs ( y );
  }
  else if ( 0.0 < x && value < 0.0 )
  {
    value = value + r8_abs ( y );
  }

  return value;
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
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
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r8_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r8_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

double r8_round ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ROUND rounds an R8 to the nearest integral value.
//
//  Example:
//
//        X         Value
//
//      1.3         1.0
//      1.4         1.0
//      1.5         1.0 or 2.0
//      1.6         2.0
//      0.0         0.0
//     -0.7        -1.0
//     -1.1        -1.0
//     -1.6        -2.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, double R8_ROUND, the rounded value.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( double ) floor ( - x + 0.5 );
  }
  else
  {
    value =   ( double ) floor (   x + 0.5 );
  }

  return value;
}
//****************************************************************************80

void r8_swap ( double &x, double &y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP switches two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double &X, &Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  double z;

  z = x;
  x = y;
  y = z;
 
  return;
}
//****************************************************************************80

double r8_uniform_ab ( double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
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
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  double value;

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed <= 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

double rd_to_jed ( double rd )

//****************************************************************************80
//
//  Purpose:
//
//    RD_TO_JED converts an RD to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RD, the RD Date.
//
//    Output, double RD_TO_JED, the Julian Ephemeris Date.
//
{
  double jed;
  double rd_epoch;

  rd_epoch = epoch_to_jed_rd ( );
  jed = rd_epoch + rd;

  return jed;
}
//****************************************************************************80

void second_borrow_common ( int &y, int &m, int &d, int &h, int &n, int &s )

//****************************************************************************80
//
//  Purpose:
//
//    SECOND_BORROW_COMMON "borrows" a minute of seconds in a common date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, &H, &N, &S, the YMDHMS date.
//
{
  while ( s < 0 )
  {
    s = s + 60;
    n = n - 1;
    minute_borrow_common ( y, m, d, h, n );
  }
  return;
}
//****************************************************************************80

void second_carry_common ( int &y, int &m, int &d, int &h, int &n, int &s )

//****************************************************************************80
//
//  Purpose:
//
//    SECOND_CARRY_COMMON: given a Common YMDHMS date, carries seconds to minutes.
//
//  Algorithm:
//
//    While 60 <= S:
//
//      decrease S by 60;
//      increase N by 1;
//      if necessary, adjust H, D, M and Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, &H, &N, &S,
//    the year, month, day, hours, minutes, seconds,
//    On output, S is between 0 and 59.
//
{
  while ( 60 <= s )
  {
    s = s - 60;
    n = n + 1;
    minute_carry_common ( y, m, d, h, n );
  }
  return;
}
//****************************************************************************80

double ss_to_jed_unix ( double s )

//****************************************************************************80
//
//  Purpose:
//
//    SS_TO_JED_UNIX converts a UNIX SS date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the UNIX date.
//
//    Output, double JED, the corresponding Julian Ephemeris Date.
//
{
  double d;
  double jed;
  double jed_epoch;

  jed_epoch = epoch_to_jed_unix ( );

  d = s / ( 24.0 * 60.0 * 60.0 );

  jed = jed_epoch + d;

  return jed;
}
//****************************************************************************80

void thanksgiving_canada ( int y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    THANKSGIVING_CANADA computes Canadian Thanksgiving for a Common year.
//
//  Discussion:
//
//    Canadian Thanksgiving occurs on the second Monday in October.
//
//  Example:
//
//    Input:
//
//      Y = 2002
//
//    Output:
//
//      M = 11
//      D = 28
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &M, &D, the month and day of Thanksgiving.
//
{
  double f;
  int w;
//
//  Determine the day of the week for 8 October, the earliest day
//  that Thanksgiving can occur.
//
  m = 10;
  d = 8;
  f = 0.0;

  w = ymdf_to_weekday_common ( y, m, d, f );
//
//  If W = 2 means this day is Monday, and day D is Thanksgiving.
//  Otherwise, figure out how to increment W to 2;
//  The same increment makes D the correct day number.
//
  if ( w < 2 )
  {
    d = d + 2 - w;
  }
  else if ( 2 < w )
  {
    d = d + 2 + 7 - w;
  }
  return;
}
//****************************************************************************80

void thanksgiving_us ( int y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    THANKSGIVING_US computes the date of Thanksgiving (US) for a Common year.
//
//  Discussion:
//
//    Thanksgiving (US) occurs on the fourth Thursday in November.
//
//  Example:
//
//    Input:
//
//      Y = 2002
//
//    Output:
//
//      M = 11
//      D = 28
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &M, &D, the month and day of Thanksgiving.
//
{
  double f;
  int w;
//
//  Determine the day of the week for 22 November, the earliest day
//  that Thanksgiving can occur.
//
  m = 11;
  d = 22;
  f = 0.0;

  w = ymdf_to_weekday_common ( y, m, d, f );
//
//  W = 5 means this day is Thursday, and day D is Thanksgiving.
//  Otherwise, figure out how to increment W to 5;
//  The same increment makes D the correct day number.
//
  if ( w < 5 )
  {
    d = d + 5 - w;
  }
  else if ( 5 < w )
  {
    d = d + 12 - w;
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
//****************************************************************************80

double transition_to_jed_common ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRANSITION_TO_JED_COMMON returns the Common calendar transition as a JED.
//
//  Discussion:
//
//    In the Common calendar, the last moment of the Julian calendar was
//      11:59 pm, 4 October 1582 Julian/CE,
//      11:59 pm, 14 October 1582 Gregorian.
//    The first minute of the Gregorian calendar ended at
//      12:01 am, 5 October 1582 Julian,
//      12:01 am, 15 October 1582 Gregorian/CE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TRANSITION_TO_JED_COMMON, the Julian Ephemeris Date 
//    of the date.
//
{
  double jed;

  jed = 2299160.5;

  return jed;
}
//****************************************************************************80

double transition_to_jed_english ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRANSITION_TO_JED_ENGLISH returns the English calendar transition as a JED.
//
//  Discussion:
//
//    In the English calendar, the last moment of the Julian calendar was
//      11:59 pm, 2 September 1752 Julian/English,
//      11:59 pm, 13 September 1752 Gregorian/CE.
//    The first minute of the Gregorian calendar ended at
//      12:01 am, 3 September 1752 Julian,
//      12:01 am, 15 September 1752 Gregorian/CE/English.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TRANSITION_TO_JED_ENGLISH, the Julian Ephemeris Date 
//    of the date.
//
{
  double jed;

  jed = 2361221.5;

  return jed;
}
//****************************************************************************80

double transition_to_jed_jed ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRANSITION_TO_JED_JED returns the JED calendar transition as a JED.
//
//  Discussion:
//
//    In Scaliger's design of the JED, three cycles with different periods
//    began on JED = 0.  These three cycles coincide once more on the
//    transition day.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TRANSITION_TO_JED_JED, the Julian Ephemeris Date 
//    of the date.
//
{
  double jed;

  jed = 2913943.0;

  return jed;
}
//****************************************************************************80

double transition_to_jed_mayan_long ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRANSITION_TO_JED_MAYAN_LONG: Mayan long count calendar transition as a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TRANSITION_TO_JED_MAYAN_LONG, the Julian Ephemeris Date 
//    of the date.
//
{
  double jed;

  jed = 2456282.5;

  return jed;
}
//****************************************************************************80

int weekday_check_common ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_CHECK_COMMON makes sure the Common weekday number is between 1 and 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, int WEEKDAY_CHECK_COMMON, the corrected weekday index.
//
{
  int w2;

  w2 = i4_wrap ( w, 1, 7 );

  return w2;
}
//****************************************************************************80

string weekday_to_name_bahai ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_BAHAI returns the name of a Bahai weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string S, the weekday name.
//
{
  string s;
  int w2;
  string weekday_name[7] = {
    "Jalal", 
    "Jamal", 
    "Kamal", 
    "Fidal",
    "Idal", 
    "Istijlal", 
    "Istiqlal" };
//
//  Check the weekday number.
//
  w2 = weekday_check_common ( w );
//
//  Return the weekday name.
//
  s = weekday_name[w2-1];

  return s;
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
//    11 May 2010
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
  string weekday_name[7] = {
    "Sunday", 
    "Monday",
    "Tuesday",
    "Wednesday",
    "Thursday",
    "Friday",
    "Saturday" };
//
//  Check the weekday number.
//
  w = weekday_check_common ( w );
//
//  Return the value.
//
  s = weekday_name[w-1];

  return s;
}
//****************************************************************************80

string weekday_to_name_common2 ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_COMMON2 returns the abbreviated name of a Common weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string WEEKDAY_TO_NAME_COMMON2, the weekday name.
//
{
  string s;
  string weekday_name[7] = {
    "Su", 
    "M",
    "Tu",
    "W",
    "Th",
    "F",
    "Sa" };
//
//  Check the weekday number.
//
  w = weekday_check_common ( w );
//
//  Return the value.
//
  s = weekday_name[w-1];

  return s;
}
//****************************************************************************80

string weekday_to_name_common3 ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_COMMON3 returns the abbreviated name of a Common weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string WEEKDAY_TO_NAME_COMMON3, the weekday name.
//
{
  string s;
  string weekday_name[7] = {
    "Sun", 
    "Mon",
    "Tue",
    "Wed",
    "Thu",
    "Fri",
    "Sat" };
//
//  Check the weekday number.
//
  w = weekday_check_common ( w );
//
//  Return the value.
//
  s = weekday_name[w-1];

  return s;
}
//****************************************************************************80

string weekday_to_name_french ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_FRENCH returns the name of a French weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string S, the weekday name.
//
{
  string s;
  int w2;
  string weekday_name[7] = {
    "Dimanche", 
    "Lundi", 
    "Mardi", 
    "Mercredi", 
    "Jeudi",
    "Vendredi", 
    "Samedi" };
//
//  Check the weekday number.
//
  w2 = weekday_check_common ( w );
//
//  Return the weekday name.
//
  s = weekday_name[w2-1];

  return s;
}
//****************************************************************************80

string weekday_to_name_german ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_GERMAN returns the name of a German weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string S, the weekday name.
//
{
  string s;
  int w2;
  string weekday_name[7] = {
    "Sonntag",
    "Montag",
    "Dienstag",
    "Mittwoch", 
    "Donnerstag",
    "Freitag",
    "Samstag" };
//
//  Check the weekday number.
//
  w2 = weekday_check_common ( w );
//
//  Return the weekday name.
//
  s = weekday_name[w2-1];

  return s;
}
//****************************************************************************80

string weekday_to_name_hebrew ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_HEBREW returns the name of a Hebrew weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string S, the weekday name.
//
{
  string s;
  int w2;
  string weekday_name[7] = {
    "Yom rishon", 
    "Yom sheni", 
    "Yom shelishi", 
    "Yom revii", 
    "Yom hamishi", 
    "Yom shishi", 
    "Sabbath" };
//
//  Check the weekday number.
//
  w2 = weekday_check_common ( w );
//
//  Return the weekday name.
//
  s = weekday_name[w2-1];

  return s;
}
//****************************************************************************80

string weekday_to_name_islamic ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_ISLAMIC returns the name of an Islamic weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string S, the weekday name.
//
{
  string s;
  int w2;
  string weekday_name[7] = {
    "Yom ilHadd", 
    "Yom litneen", 
    "Yom ittalat", 
    "Yom larba", 
    "Yom ilkhamiis", 
    "Yom ilguma", 
    "Yom issabt" };
//
//  Check the weekday number.
//
  w2 = weekday_check_common ( w );
//
//  Return the weekday name.
//
  s = weekday_name[w2-1];

  return s;
}
//****************************************************************************80

string weekday_to_name_italian ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_ITALIAN returns the name of an Italian weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string S, the weekday name.
//
{
  string s;
  int w2;
  string weekday_name[7] = {
    "domenica", 
    "lunedi", 
    "martedi", 
    "mercoledi", 
    "giovedi", 
    "venerdi", 
    "sabato" };
//
//  Check the weekday number.
//
  w2 = weekday_check_common ( w );
//
//  Return the weekday name.
//
  s = weekday_name[w2-1];

  return s;
}
//****************************************************************************80

string weekday_to_name_republican ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_REPUBLICAN returns the name of a Republican weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string S, the weekday name.
//
{
  string s;
  int w2;
  string weekday_name[10] = {
    "Primedi", 
    "Duodi", 
    "Tridi", 
    "Quartidi", 
    "Quintidi", 
    "Sextidi", 
    "Septidi", 
    "Octidi", 
    "Nonidi", 
    "Decadi" };

  if ( w < 1 || 10 < 2 )
  {
    s = "????";
  }
  else
  {
    s = weekday_name[w2-1];
  }

  return s;
}
//****************************************************************************80

string weekday_to_name_roman ( int w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_TO_NAME_ROMAN returns the name of a Roman weekday.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int W, the weekday index.
//
//    Output, string S, the weekday name.
//
{
  string s;
  int w2;
  string weekday_name[7] = {
    "Dies Solis",
    "Dies Lunae",
    "Dies Martis",
    "Dies Mercurii",
    "Dies Iovis",
    "Dies Veneris",
    "Dies Saturni" };
//
//  Check the weekday number.
//
  w2 = weekday_check_common ( w );
//
//  Return the weekday name.
//
  s = weekday_name[w2-1];

  return s;
}
//****************************************************************************80

int y_astronomical_to_common ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_ASTRONOMICAL_TO_COMMON converts an Astronomical year to a Common year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the astronomical year.
//
//    Output, int Y_ASTRONOMICAL_TO_COMMON, the Common year.
//
{
  int y2;

  if ( y <= 0 )
  {
    y2 = y - 1;
  }
  else
  {
    y2 = y;
  }
  return y2;
}
//****************************************************************************80

void y_check_alexandrian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_ALEXANDRIAN checks an Alexandrian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must be strictly positive.
{
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_ALEXANDRIAN - Fatal error!\n";
    cerr << "  Year <= 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_bahai ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_BAHAI checks a Bahai year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must be strictly positive.
{
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_BAHAI - Fatal error!\n";
    cerr << "  Year <= 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_common ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_COMMON checks a Common year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must not be 0.
{
  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_COMMON - Fatal error!\n";
    cerr << "  Year 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_eg_civil ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_EG_CIVIL checks an Egyptian Civil year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must be strictly positive.
{
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_EG_CIVIL - Fatal error!\n";
    cerr << "  Year <= 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_english ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_ENGLISH checks an English year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must not be 0.
//
{
  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_ENGLISH - Fatal error!\n";
    cerr << "  Year 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_greek ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_GREEK checks a Greek year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must be strictly positive.
{
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_GREEK - Fatal error!\n";
    cerr << "  Year <= 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_gregorian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_GREGORIAN checks a Gregorian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must not be 0.
//
{
  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_GREGORIAN - Fatal error!\n";
    cerr << "  Year 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_hebrew ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_HEBREW checks a Hebrew year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must be strictly positive.
{
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_HEBREW - Fatal error!\n";
    cerr << "  Year <= 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_islamic ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_ISLAMIC checks an Islamic year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must be strictly positive.
{
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_ISLAMIC - Fatal error!\n";
    cerr << "  Year <= 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_julian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_JULIAN checks a Julian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must not be 0.
//
{
  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_JULIAN - Fatal error!\n";
    cerr << "  Year 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_republican ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_REPUBLICAN checks a Republican year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must be positive.
//
{
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_REPUBLICAN - Fatal error!\n";
    cerr << "  Year <= 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void y_check_roman ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_CHECK_ROMAN checks a Roman year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year, which must be positive.
//
{
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "Y_CHECK_ROMAN - Fatal error!\n";
    cerr << "  Year <= 0 is illegal.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int y_common_to_astronomical ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_COMMON_TO_ASTRONOMICAL converts a Common year to an Astronomical year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the Common year.
//
//    Output, int Y_COMMON_TO_ASTRONOMICAL, the Astronomical year.
//
{
  int y2;

  if ( y < 0 )
  {
    y2 = y + 1;
  }
  else if ( y == 0 )
  {
    cerr << "\n";
    cerr << "Y_COMMON_TO_ASTRONOMICAL - Fatal error!\n";
    cerr << "  Common calendar does not have a year 0.\n";
    exit ( 1 );
  }
  else
  {
    y2 = y;
  }

  return y2;
}
//****************************************************************************80

int y_julian_to_roman ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_JULIAN_TO_ROMAN converts a Julian year to a Roman year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the Julian year.
//
//    Output, int Y_JULIAN_TO_ROMAN, the corresponding Roman year.
//
{
  int y2;

  y_check_julian ( y );

  if ( y < 0 )
  {
    y = y + 1;
  }

  y2 = y + 753;

  return y2;
}
//****************************************************************************80

int y_roman_to_julian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_ROMAN_TO_JULIAN converts a Roman year to a Julian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the Roman year.
//
//    Output, int Y_ROMAN_TO_JULIAN, the corresponding Julian year.
//
{
  int y2;

  y2 = y - 753;

  if ( y2 <= 0 )
  {
    y2 = y2 - 1;
  }

  return y2;
}
//****************************************************************************80

string y_to_s_common ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    Y_TO_S_COMMON writes a Common year into a string.
//
//  Format:
//
//    YearNumber BCE
//    YearNumber CE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, string S, a representation of the year.
//
{
  char ch_vec[10];
  string s;

  y_check_common ( y );

  if ( y < 0 )
  {
    sprintf ( ch_vec, "BCE %d", - y );
  }
  else
  {
    sprintf ( ch_vec, "CE %d", y );
  }
  s = string ( ch_vec );

  return s;
}
//****************************************************************************80

bool year_is_embolismic_eg_lunar ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_EMBOLISMIC_EG_LUNAR: TRUE if the Egyptian Lunar year was embolismic.
//
//  Discussion:
//
//    This is just a "fake" function, which does repeat every 25 years,
//    and has 9 embolismic and 16 common years in that cycle, but with
//    a pattern I just made up for now.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_EMBOLISMIC_EG_LUNAR, TRUE if the year
//    was embolismic.
//
{
  bool value;
  int y2;

  y2 = ( y - 1 ) % 25;

  if ( ( y2 % 3 ) == 0 )
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

bool year_is_embolismic_greek ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_EMBOLISMIC_GREEK returns TRUE if the Greek year was embolismic.
//
//  Discussion:
//
//    Apparently, the Greek calendar was emended haphazardly.  This
//    routine does not attempt to follow that historical pattern, and
//    just uses the Hebrew calendar pattern for now.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_EMBOLISMIC_GREEK, TRUE if the year was embolismic.
//
{
  bool value;

  if ( 12 <= i4_modp ( 7 * y + 13, 19 ) )
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

bool year_is_embolismic_hebrew ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_EMBOLISMIC_HEBREW returns TRUE if the Hebrew year was embolismic.
//
//  Discussion:
//
//    In a 19 year cycle, there are 7 embolismic years.  During these years,
//    an extra month, "Adar II", (sometimes called "Veadar") is inserted after
//    the month of Adar.  Nonembolismic years are called "common" years.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_EMBOLISMIC_HEBREW, TRUE if the year was embolismic.
//
{
  bool value;

  if ( 12 <= i4_modp ( 7 * y + 13, 19 ) )
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

bool year_is_leap_alexandrian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_ALEXANDRIAN: TRUE if the Alexandrian year was a leap year.
//
//  Discussion:
//
//    The Alexandrian year, which started on the 29th of August of the Julian
//    year, was a leap year if it included the bissextile day of the Julian
//    calendar.  In other words, if the Alexandrian year BEGAN in a Julian year
//    that preceded a Julian leap year, then the Alexandrian year was a leap
//    year.
//
//    We deem year AX 1 to have begun in Julian 23 BC.  Julian 21 BC was
//    theoretically a leap year, so AX 2 was a leap year, as was AX 6, AX 10,
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_ALEXANDRIAN, TRUE if the year was a leap year.
//
{
  bool value;

  if ( ( y % 4 ) == 2 )
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

bool year_is_leap_bahai ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_BAHAI returns TRUE if the Bahai year was a leap year.
//
//  Discussion:
//
//    The leap year rules are the same as those used in the Gregorian
//    calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_BAHAI, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  if ( y <= 0 )
  {
    value = false;
    return value;
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
//****************************************************************************80

bool year_is_leap_common ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_COMMON returns TRUE if the Common year was a leap year.
//
//  Discussion:
//
//    The "common" calendar is meant to be the calendar which is Julian up to
//    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
//
//  Algorithm:
//
//    If ( the year is less than 0 ) then
//
//      if the year+1 is divisible by 4 then
//        the year is a leap year.
//
//    else if ( the year is 0 ) then
//
//      the year is not a leap year ( in fact, it's illegal )
//
//    else if ( the year is no greater than 1582 ) then
//
//      if the year is divisible by 4 then
//        the year is a leap year.
//
//    else if (
//      the year is divisible by 4 and
//      ( the year is not divisible by 100
//      or
//      the year is divisible by 400 )
//      ) then
//        the year is a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_COMMON, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;
  int y2;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_LEAP_COMMON - Fatal error!\n";
    cerr << "  Year 0 is illegal.\n";
    exit ( 1 );
  }
//
//  BC years have to have 1 added to them to make a proper leap year evaluation.
//
  y2 = y_common_to_astronomical ( y );

  if ( y2 <= 1582 )
  {
    if ( i4_modp ( y2, 4 ) == 0 )
    {
      value = true;
    }
    else
    {
      value = false;
    }
  }
  else
  {
    if ( i4_modp ( y2, 400 ) == 0 )
    {
      value = true;
    }
    else if ( i4_modp ( y2, 100 ) == 0 )
    {
      value = false;
    }
    else if ( i4_modp ( y2, 4 ) == 0 )
    {
      value = true;
    }
    else
    {
      value = false;
    }
  }
  return value;
}
//****************************************************************************80

bool year_is_leap_coptic ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_COPTIC returns TRUE if the Coptic year was a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Nachum Dershowitz, Edward Reingold,
//    Calendrical Calculations,
//    Cambridge, 1997, page 58.
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_COPTIC, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  if ( y <= 0 )
  {
    value = false;
    return value;
  }

  if ( ( y % 4 ) == 3 )
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

bool year_is_leap_eg_lunar ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_EG_LUNAR: TRUE if the Egyptian Lunar year was a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_EG_LUNAR, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  if ( y <= 0 )
  {
    value = false;
    return value;
  }

  if ( ( y % 5 ) == 0 )
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

bool year_is_leap_english ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_ENGLISH returns TRUE if the English year was a leap year.
//
//  Algorithm:
//
//    If ( the year is less than 0 ) then
//
//      if the year+1 is divisible by 4 then
//        the year is a leap year.
//
//    else if ( the year is 0 ) then
//
//      the year is not a leap year ( in fact, it is illegal )
//
//    else if ( the year is no greater than 1752 ) then
//
//      if the year is divisible by 4 then
//        the year is a leap year.
//
//    else if (
//      the year is divisible by 4 and
//      ( the year is not divisible by 100
//      or
//      the year is divisible by 400 )
//      ) then
//        the year is a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_ENGLISH, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  if ( y == 0 )
  {
    value = false;
    return value;
  }
//
//  BC years have to have 1 added to them to make a proper leap year evaluation.
//
  y = y_common_to_astronomical ( y );

  if ( y <= 1752 )
  {
    if ( i4_modp ( y, 4 ) == 0 )
    {
      value = true;
    }
    else
    {
      value = false;
    }
  }
  else
  {
    if ( i4_modp ( y, 400 ) == 0 )
    {
      value = true;
    }
    else if ( i4_modp ( y, 100 ) == 0 )
    {
      value = false;
    }
    else if ( i4_modp ( y, 4 ) == 0 )
    {
      value = true;
    }
    else
    {
      value = false;
    }
  }
  return value;
}
//****************************************************************************80

bool year_is_leap_ethiopian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_ETHIOPIAN returns TRUE if the Ethiopian year was a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Nachum Dershowitz, Edward Reingold,
//    Calendrical Calculations,
//    Cambridge, 1997, page 58.
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_ETHIOPIAN, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  if ( y <= 0 )
  {
    value = false;
    return value;
  }

  if ( ( y % 4 ) == 3 )
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

bool year_is_leap_greek ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_GREEK returns TRUE if the Greek year was a leap year.
//
//  Discussion:
//
//    The actual practice of adding the extra day to the Greek calendar
//    seems to have been unmethodical.  Here, we simply make up a rule
//    as a placeholder for now.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_GREEK, TRUE if the year was a leap year.
//
{
  bool value;

  if ( year_is_embolismic_greek ( y ) && ( ( y % 3 ) == 0 ) )
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
//    23 March 2010
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

  if ( y == 0 ) 
  {
    value = false;
    return value;
  }
//
//  BC years have to have 1 added to them to make a proper leap year evaluation.
//
  y = y_common_to_astronomical ( y );

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
//****************************************************************************80

bool year_is_leap_iranian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_IRANIAN returns TRUE if the Iranian year was a leap year.
//
//  Discussion:
//
//    I don't know the rule for this, so I'm just setting it FALSE for now.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_IRANIAN, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  value = false;

  return value;
}
//****************************************************************************80

bool year_is_leap_islamic ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_ISLAMIC returns TRUE if the Islamic year was a leap year.
//
//  Discussion:
//
//    In a 30 year cycle, there are 11 leap years, years 2, 5, 7, 10, 13,
//    16, 18, 21, 24, 26 and 29.  During these years, the 12th month has
//    30 days instead of 29.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_ISLAMIC, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  if ( i4_modp ( 11 * y + 14, 30 ) < 11 )
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

bool year_is_leap_julian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_JULIAN returns TRUE if the Julian year was a leap year.
//
//  Algorithm:
//
//    If ( Y < 0 and Y+1 is divisible by 4 ) then
//      the year is a leap year.
//    else if ( Y == 0 ) then
//      the year is illegal
//    else if ( 0 < Y and Y is divisible by 4 ) then
//      the year is a leap year.
//    else
//      the year is NOT a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_JULIAN, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;

  if ( y == 0 )
  {
    value = false;
    return value;
  }

  y = y_common_to_astronomical ( y );

  if ( i4_modp ( y, 4 ) == 0 )
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

bool year_is_leap_persian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_PERSIAN returns TRUE if the Persian year was a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Nachum Dershowitz, Edward Reingold,
//    Calendrical Calculations,
//    Cambridge, 1997, page 58.
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_PERSIAN, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;
  int y2;
  int y3;

  if ( y <= 0 )
  {
    y2 = y - 473;
  }
  else
  {
    y2 = y - 474;
  }

  y3 = 474 + ( y2 % 2820 );

  if ( ( 682 * ( y3 + 38 ) ) % 2816 < 682 )
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

bool year_is_leap_republican ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_REPUBLICAN returns TRUE if the Republican year was a leap year.
//
//  Discussion:
//
//    The French Republican calendar was in use for 14 years.
//    In that time, years 3, 7 and 11 were designated as leap years.
//    The easiest way to harmonize the rules and history is to apply
//    the leap year rules to Y+1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_REPUBLICAN, TRUE if the year was a leap year,
//    FALSE otherwise.
//
{
  bool value;
  int y2;

  y2 = y;

  y_check_republican ( y2 );

  value = false;

  if ( ( y2 + 1 ) % 4 == 0 )
  {
    value = true;
    if ( ( y2 + 1 ) % 100 == 0 )
    {
      value = false;
      if ( ( y2 + 1 ) % 400 == 0 )
      {
        value = true;
        if ( ( y2 + 1 ) % 4000 == 0 )
        {
          value = false;
        }
      }
    }
  }
  return value;
}
//****************************************************************************80

bool year_is_leap_roman ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_IS_LEAP_ROMAN returns TRUE if the Roman year was a leap year.
//
//  Discussion:
//
//    For our unrealistic and idealized Roman calendar, we are going to
//    take a year to have been a leap year if the corresponding year in
//    the idealized Julian calendar was a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, bool YEAR_IS_LEAP_ROMAN, TRUE if the year was a leap year.
//
{
  bool value;
  int y2;

  y_check_roman ( y ); 

  y2 = y_roman_to_julian ( y );

  value = year_is_leap_julian ( y2 );

  return value;
}
//****************************************************************************80

int year_length_alexandrian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_ALEXANDRIAN returns the number of days in an Alexandrian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_ALEXANDRIAN, the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_alexandrian ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_bahai ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_BAHAI returns the number of days in a Bahai year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_BAHAI, the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_bahai ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_common ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_COMMON returns the number of days in a Common year.
//
//  Discussion:
//
//    The "common" calendar is meant to be the calendar which is Julian up to
//    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
//
//    If Y is 0, then the routine returns 0, reflecting the fact that
//    there was officially no year 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_COMMON, the number of 
//    days in the year.
//
{
  int days;

  if ( y == 0 ) 
  {
    days = 0;
  }
  else if ( y == 1582 )
  {
    days = 355;
  }
  else if ( year_is_leap_common ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_coptic ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_COPTIC returns the number of days in a Coptic year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_COPTIC, the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_coptic ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_eg_civil ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_EG_CIVIL returns the number of days in an Egyptian Civil year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_EG_CIVIL, the number of 
//    days in the year.
//
{
  int days;

  days = 365;

  return days;
}
//****************************************************************************80

int year_length_eg_lunar ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_EG_LUNAR returns the number of days in an Egyptian lunar year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_EG_LUNAR, the number of 
//    days in the year.
//
{
  int days;

  if ( ! year_is_embolismic_eg_lunar ( y ) )
  {
    days = 354;
  }
  else 
  {
    days = 384;
  }

  if ( year_is_leap_eg_lunar ( y ) )
  {
    days = days + 1;
  }
  return days;
}
//****************************************************************************80

int year_length_english ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_ENGLISH returns the number of days in an English year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_ENGLISH, the number of 
//    days in the year.
//
{
  int days;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_LENGTH_ENGLISH - Fatal error!\n";
    cerr << "  Illegal Y = 0.\n";
    exit ( 1 );
  }

  if ( y == 1752 )
  {
    days = 355;
  }
  else if ( year_is_leap_english ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_ethiopian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_ETHIOPIAN returns the number of days in an Ethiopian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_ETHIOPIAN, the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_ethiopian ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_greek ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_GREEK returns the number of days in a Greek year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_GREEK, the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_embolismic_greek ( y ) )
  {
    days = 386;
    if ( year_is_leap_greek ( y ) )
    {
      days = days + 1;
    }
  }
  else
  {
    days = 357;
  }
  return days;
}
//****************************************************************************80

int year_length_gregorian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_GREGORIAN returns the number of days in a Greegorian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_GREGORIAN, the number of 
//    days in the year.
//
{
  int days;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_LENGTH_GREGORIAN - Fatal error!\n";
    cerr << "  Illegal Y = 0.\n";
    exit ( 1 );
  }

  if ( year_is_leap_gregorian ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_hebrew ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_HEBREW returns the number of days in a Hebrew year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_HEBREW, the number of 
//    days in the year.
//
{
  int days;
  double jed;
  double jed2;
  int y2;

  jed = new_year_to_jed_hebrew ( y );

  y2 = y + 1;
  jed2 = new_year_to_jed_hebrew ( y2 );

  days = r8_nint ( jed2 - jed );

  return days;
}
//****************************************************************************80

double year_length_hindu_solar ( )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_HINDU_SOLAR returns the number of days in a Hindu solar year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double YEAR_LENGTH_HINDU_SOLAR, the number of 
//    days in the year.
//
{
  double days;

  days = 1577917828.0 / 4320000.0;

  return days;
}
//****************************************************************************80

int year_length_islamic ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_ISLAMIC returns the number of days in an Islamic year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_ISLAMIC the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_islamic ( y ) )
  {
    days = 355;
  }
  else
  {
    days = 354;
  }
  return days;
}
//****************************************************************************80

int year_length_julian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_JULIAN returns the number of days in a Julian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_JULIAN, the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_julian ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

double year_length_lunar ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_LUNAR returns the number of days in a "lunar year".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, double YEAR_LENGTH_LUNAR, the number of 
//    days in the year.
//
{
  double days;

  days = 354.3671;

  return days;
}
//****************************************************************************80

int year_length_persian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_PERSIAN returns the number of days in a Persian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_PERSIAN the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_persian ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_republican ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_REPUBLICAN returns the number of days in a Republican year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_REPUBLICAN, the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_republican ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

int year_length_roman ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_ROMAN returns the number of days in a Roman year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_ROMAN, the number of 
//    days in the year.
//
{
  int days;

  if ( year_is_leap_roman ( y ) )
  {
    days = 366;
  }
  else
  {
    days = 365;
  }
  return days;
}
//****************************************************************************80

double year_length_solar ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_SOLAR returns the number of days in a "solar year".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, double YEAR_LENGTH_SOLAR, the number of 
//    days in the year.
//
{
  double days;

  if ( y < 1 )
  {
    y = y + 1;
  }

  if ( y < - 4000 )
  {
    days = 365.2424992;
  }
  else if ( y <= 2000 )
  {
    days = ( ( double ) ( 2000 - y ) * 365.2424992 
           + ( double ) ( 4000 + y ) * 365.2421897 )
           /              6000.0;
  }
  else
  {
    days = 365.2421897;
  }
  return days;
}
//****************************************************************************80

int year_length_months_alexandrian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_ALEXANDRIAN: number of months in an Alexandrian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_ALEXANDRIAN, the 
//    number of months in the year.
//
{
  int months;

  months = 13;

  return months;
}
//****************************************************************************80

int year_length_months_bahai ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_BAHAI: number of months in a Bahai year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_BAHAI, the 
//    number of months in the year.
//
{
  int months;

  months = 20;

  return months;
}
//****************************************************************************80

int year_length_months_common ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_COMMON returns the number of months in a Common year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_COMMON, the number of months
//    in the year.
//
{
  int value;

  value = 12;

  return value;
}
//****************************************************************************80

int year_length_months_coptic ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_COPTIC: number of months in a Coptic year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_COPTIC, the 
//    number of months in the year.
//
{
  int months;

  months = 13;

  return months;
}
//****************************************************************************80

int year_length_months_eg_civil ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_EG_CIVIL: number of months in an Egyptian Civil year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_EG_CIVIL, the 
//    number of months in the year.
//
{
  int months;

  months = 13;

  return months;
}
//****************************************************************************80

int year_length_months_eg_lunar ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_EG_LUNAR: number of months in an Egyptian lunar year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_EG_LUNAR, the 
//    number of months in the year.
//
{
  int months;

  if ( year_is_embolismic_eg_lunar ( y ) )
  {
    months = 13;
  }
  else
  {
    months = 12;
  }
  return months;
}
//****************************************************************************80

int year_length_months_english ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_ENGLISH returns the number of months in an English year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_ENGLISH, the number of months
//    in the year.
//
{
  int value;

  value = 12;

  return value;
}
//****************************************************************************80

int year_length_months_ethiopian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_ETHIOPIAN returns the number of months in an Ethiopian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_ETHIOPIAN, the number of months
//    in the year.
//
{
  int value;

  value = 13;

  return value;
}
//****************************************************************************80

int year_length_months_greek ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_GREEK returns the number of months in a Greek year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_GREEK, the number of months
//    in the year.
//
{
  int value;

  if ( year_is_embolismic_greek ( y ) )
  {
    value = 13;
  }
  else
  {
    value = 12;
  }
  return value;
}
//****************************************************************************80

int year_length_months_gregorian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_GREGORIAN returns the number of months in a Gregorian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_GREGORIAN, the number of months
//    in the year.
//
{
  int value;

  value = 12;

  return value;
}
//****************************************************************************80

int year_length_months_hebrew ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_HEBREW returns the number of months in a Hebrew year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_HEBREW, the number of months
//    in the year.
//
{
  int value;

  if ( year_is_embolismic_hebrew ( y ) )
  {
    value = 13;
  }
  else
  {
    value = 12;
  }
  return value;
}
//****************************************************************************80

int year_length_months_hindu_lunar ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_HINDU_LUNAR returns the number of months in a Hindu lunar year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_HINDU_LUNAR, the number of months
//    in the year.
//
{
  int value;

  value = 12;

  return value;
}
//****************************************************************************80

int year_length_months_hindu_solar ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_HINDU_SOLAR returns the number of months in a Hindu solar year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_HINDU_SOLAR, the number of months
//    in the year.
//
{
  int value;

  value = 12;

  return value;
}
//****************************************************************************80

int year_length_months_islamic ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_ISLAMIC returns the number of months in an Islamic year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_ISLAMIC, the number of months
//    in the year.
//
{
  int value;

  value = 12;

  return value;
}
//****************************************************************************80

int year_length_months_julian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_JULIAN returns the number of months in a Julian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_JULIAN, the number of months
//    in the year.
//
{
  int value;

  value = 12;

  return value;
}
//****************************************************************************80

int year_length_months_persian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_PERSIAN returns the number of months in a Persian year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_PERSIAN, the number of months
//    in the year.
//
{
  int value;

  value = 12;

  return value;
}
//****************************************************************************80

int year_length_months_republican ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_REPUBLICAN returns the number of months in a French Republican year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_REPUBLICAN, the number of months
//    in the year.
//
{
  int value;

  value = 13;

  return value;
}
//****************************************************************************80

int year_length_months_roman ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_LENGTH_MONTHS_ROMAN returns the number of months in a Roman year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year to be checked.
//
//    Output, int YEAR_LENGTH_MONTHS_ROMAN, the number of months
//    in the year.
//
{
  int value;

  value = 13;

  return value;
}
//****************************************************************************80

void year_to_dominical_common ( int y, int &n1, int &n2 )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_TO_DOMINICAL_COMMON: dominical numbers, Common calendar.
//
//  Discussion:
//
//    The Julian calendar calculations are used through the year 1582,
//    and the Gregorian thereafter.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &N1, &N2, the dominical numbers for the year.
//    If Y is a leap year, then N1 applies before March 1, and N2 after.
//    If Y is not a leap year, then N1 applies throughout the year,
//    and N2 is returned as N1.
//
{
  if ( y <= 1582 )
  {
    year_to_dominical_julian ( y, n1, n2 );
  }
  else
  {
    year_to_dominical_gregorian ( y, n1, n2 );
  }
  return;
}
//****************************************************************************80

void year_to_dominical_gregorian ( int y, int &n1, int &n2 )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_TO_DOMINICAL_GREGORIAN: dominical numbers, Gregorian calendar.
//
//  Discussion:
//
//    The days of each year are numbered with "calendar letters", with
//    January 1 having letter 'A', January 7 having letter 'G', and
//    the cycle then repeating with January 8 having letter 'A'.
//
//    This cycle is independent of the weekday cycle.  If a year is
//    not a leap year, then all Sundays have the same calendar letter.
//    This is called the dominical letter of the year.  If a year is
//    a leap year, then all Sundays before March 1 have one calendar
//    letter, and all Sundays after have another (namely, the calendar
//    letter one position earlier in the cycle).
//
//    Using the correspondence A = 1, B = 2, ..., we may speak of
//    the dominical number of a year, or dominical numbers for a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int &N1, &N2, the dominical numbers for the year.
//    If Y is a leap year, then N1 applies before March 1, and N2 after.
//    If Y is not a leap year, then N1 applies throughout the year,
//    and N2 is returned as N1.
//
{
  int p1;
  int y2;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_DOMINICAL_GREGORIAN - Fatal error!\n";
    cerr << "  Illegal input Y = 0.\n";
    exit ( 1 );
  }

  y2 = y_common_to_astronomical ( y );

  p1 = y2 + ( y2 / 4 ) - ( y2 / 100 ) + ( y2 / 400 ) - 1;
  n1 = 7 - i4_modp ( p1, 7 );

  if ( year_is_leap_gregorian ( y2 ) )
  {
    n2 = n1;
    p1 = p1 - 1;
    n1 = 7 - i4_modp ( p1, 7 );
  }
  else
  {
    n2 = n1;
  }
  return;
}
//****************************************************************************80

void year_to_dominical_julian ( int y, int &n1, int &n2 )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_TO_DOMINICAL_JULIAN: dominical numbers, Julian calendar.
//
//  Discussion:
//
//    The days of each year are numbered with "calendar letters", with
//    January 1 having letter 'A', January 7 having letter 'G', and
//    the cycle then repeating with January 8 having letter 'A'.
//
//    This cycle is independent of the weekday cycle.  If a year is
//    not a leap year, then all Sundays have the same calendar letter.
//    This is called the dominical letter of the year.  If a year is
//    a leap year, then all Sundays before March 1 have one calendar
//    letter, and all Sundays after have another (namely, the calendar
//    letter one position earlier in the cycle).
//
//    Using the correspondence A = 1, B = 2, ..., we may speak of
//    the dominical number of a year, or dominical numbers for a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.  The year 0 is illegal input.
//
//    Output, int &N1, &N2, the dominical numbers for the year.
//    If Y is a leap year, then N1 applies before March 1, and N2 after.
//    If Y is not a leap year, then N1 applies throughout the year,
//    and N2 is returned as N1.
//
{
  int p1;
  int y2;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_DOMINICAL_JULIAN - Fatal error!\n";
    cerr << "  Illegal input Y = 0.\n";
    exit ( 1 );
  }

  y2 = y_common_to_astronomical ( y );

  p1 = y2 + ( y2 / 4 ) + 4;
  n1 = 7 - i4_modp ( p1, 7 );

  if ( year_is_leap_julian ( y2 ) )
  {
    n2 = n1;
    p1 = p1 - 1;
    n1 = 7 - i4_modp ( p1, 7 );
  }
  else
  {
    n2 = n1;
  }
  return;
}
//****************************************************************************80

int year_to_epact_gregorian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_TO_EPACT_GREGORIAN returns the epact of a Gregorian year.
//
//  Discussion:
//
//    The epact of a year is the age in days of the notional moon on
//    the first day of the year.  If the year begins with a new moon,
//    the epact is zero.  If the new moon occurred the day before,
//    the epact is 1.  There is a unique epact for every golden number.
//
//    The Gregorian epact calculation is an adjustment to the Julian
//    calculation that takes into account the shift of the calendar
//    to restore the vernal equinox to March 21, and the adjustment to
//    the average length of the year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999.
//
//  Parameters:
//
//    Input, int Y, the year.  The year 0 is illegal input.
//
//    Output, int YEAR_TO_EPACT_GREGORIAN, the epact, between 0 and 28.
//
{
  int e;
  int g;
  int h;
  int q;
  int y2;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_EPACT_GREGORIAN - Fatal error!\n";
    cerr << "  Illegal input Y = 0.\n";
    exit ( 1 );
  }

  y2 = y_common_to_astronomical ( y );

  g = year_to_golden_number ( y );

  h = ( y2 / 100 );

  q = h - ( h / 4 );

  e = ( 57 + 11 * g - q + ( h - ( h - 17 ) / 25 ) / 3 ) % 30;

  if ( e == 24 || ( e == 25 && 12 <= g ) )
  {
    e = e + 1;
  }

  return e;
}
//****************************************************************************80

int year_to_epact_julian ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_TO_EPACT_JULIAN returns the epact of a Julian year.
//
//  Discussion:
//
//    The epact of a year is the age in days of the notional moon on
//    the first day of the year.  If the year begins with a new moon,
//    the epact is zero.  If the new moon occurred the day before,
//    the epact is 1.  There is a unique epact for every golden number.
//
//    Bear in mind that the notional moon is not the one in the sky,
//    but a theoretical one that satisfactorily approximates the behavior
//    of the real one, but which is tame enough to be described by a formula.
//
//  Example:
//
//    Year  Golden Number  Epact
//
//      1 BC     1           8
//      1 AD     2          19
//      2 AD     3           0
//      3 AD     4          11
//      4 AD     5          22
//      5 AD     6           3
//      6 AD     7          14
//      7 AD     8          25
//      8 AD     9           6
//      9 AD    10          17
//     10 AD    11          28
//     11 AD    12           9
//     12 AD    13          20
//     13 AD    14           1
//     14 AD    15          12
//     15 AD    16          23
//     16 AD    17           4
//     17 AD    18          15
//     18 AD    19          26
//     19 AD     1           8
//     20 AD     2          19
//   1066 AD     3           0
//   1900 AD     1           8
//   1919 AD     1           8
//   1938 AD     1           8
//   1957 AD     1           8
//   1976 AD     1           8
//   1995 AD     1           8
//   2014 AD     1           8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999.
//
//  Parameters:
//
//    Input, int Y, the year.  The year 0 is illegal input.
//
//    Output, int YEAR_TO_EPACT_JULIAN, the epact, between 0 and 28.
//
{
  int e;
  int g;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_EPACT_JULIAN - Fatal error!\n";
    cerr << "  Illegal input Y = 0.\n";
    exit ( 1 );
  }

  g = year_to_golden_number ( y );

  e = i4_wrap ( 11 * g - 3, 0, 29 );

  return e;
}
//****************************************************************************80

int year_to_golden_number ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_TO_GOLDEN_NUMBER returns the golden number of a Common year.
//
//  Discussion:
//
//    Nineteen solar years are very close to 235 lunations.  Calendars
//    that try to keep track of both the sun and moon often make use of
//    this fact, ascribed to the Greek astronomer Meton.
//
//    While trying to determine a formula for Easter, Dionysus Exiguus
//    symbolized the place of each year in its Metonic cycle by a
//    "golden number" between 1 and 19.  The numbering began with the
//    year 1 BC, assigned the golden number of 1.  The following year,
//    1 AD, got the golden number of 2, and after that it gets easier.
//
//    The same golden year calculation is done for years in the Julian
//    or Gregorian calendar.
//
//  Example:
//
//    Year  Golden Number
//
//      1 BC     1
//      1 AD     2
//      2 AD     3
//     18 AD    19
//     19 AD     1
//     20 AD     2
//   1066 AD     3
//   1900 AD     1
//   1919 AD     1
//   1938 AD     1
//   1957 AD     1
//   1976 AD     1
//   1995 AD     1
//   2014 AD     1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int YEAR_TO_GOLDEN_NUMBER the golden number, between 1 and 19.  
//    This records the position of the year in the 19 year Metonic cycle.
//
{
  int g;
  int y2;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_GOLDEN_NUMBER - Fatal error!\n";
    cerr << "  Illegal input Y = 0.\n";
    exit ( 1 );
  }
//
//  We assume that BC years come in as negative numbers, and that
//  the year before 1 AD is 1 BC.  So add 1 to any negative value
//  so that the arithmetic works.
//
  y2 = y_common_to_astronomical ( y );

  g = i4_wrap ( y2 + 1, 1, 19 );

  return g;
}
//****************************************************************************80

int year_to_indiction_common ( int y )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_TO_INDICTION_COMMON returns the indiction number of a Common year.
//
//  Discussion:
//
//    The Roman empire had a taxation cycle that, at one time, comprised
//    15 years.  As is typical in calendrical matters, the actual length
//    of this cycle and the time that the cycle began varied from place
//    to place and time to time, and historians even disagree about the
//    indiction cycle given a specific place and time.  Nonetheless,
//    it is customary to retrospectively impose a uniform and regular
//    indiction cycle on the ancient world.  (The 15 year indiction cycle,
//    in fact, was factored into Scaliger's determination of an appropriate
//    starting point for the Julian Ephemeris Date.)
//
//  Example:
//
//    Year  Indiction Number
//
//      3 BC     1
//      2 BC     2
//      1 BC     3
//      1 AD     4
//     10 AD    13
//     11 AD    14
//     12 AD    15
//     13 AD     1
//     14 AD     2
//     15 AD     3
//     26 AD    14
//     27 AD    15
//     28 AD     1
//   1900 AD    13
//   2000 AD     8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the year.
//
//    Output, int YEAR_TO_INDICTION_COMMON, the indiction number, between 1 and 15.
//
{
  int  i;
  int y2;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_INDICTION_COMMON - Fatal error!\n";
    cerr << "  Illegal input Y = 0.\n";
    exit ( 1 );
  }
//
//  We assume that BC years come in as negative numbers, and that
//  the year before 1 AD is 1 BC.  So add 1 to any negative value
//  so that the arithmetic works.
//
  y2 = y_common_to_astronomical ( y );

  i = i4_wrap ( y2 + 3, 1, 15 );

  return i;
}
//****************************************************************************80

void year_to_scaliger_common ( int y, int &c1, int &c2, int &c3, int &r1, 
  int &r2, int &r3 )

//****************************************************************************80
//
//  Purpose:
//
//    YEAR_TO_SCALIGER_COMMON converts a Common year to its Scaliger indices.
//
//  Discussion:
//
//    The year 4713 BCE was chosen by Joseph Scaliger for the start of
//    his Julian Ephemeris Date system, because three cycles coincided
//    in that year, the 28 year Julian calendar cycle, the 19 year Metonic
//    cycle, and the 15 year Roman Indiction cycle.  Thus, the year
//    4713 BCE has Scaliger index (1,1,1).  Each subsequent year has a distinct
//    set of Scaliger indices until 7980 years later, when the year
//    3266 CE will again have the Scaliger index (1,1,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, the Common year.
//
//    Output, int &C1, &C2, &C3, the number of completed 
//    Julian, Metonic and Indiction cycles.
//
//    Output, int &R1, &R2, &R3, the Julian, Metonic and 
//    Indiction cycle numbers that make up the Scaliger index.
//
{
  int y2;

  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_SCALIGER_COMMON - Fatal error!\n";
    cerr << "  Illegal input Y = 0.\n";
    exit ( 1 );
  }
//
//  Adjust for missing year 0.
//
  if ( y < 0 )
  {
    y2 = y + 1;
  }
  else
  {
    y2 = y;
  }
//
//  Now shift so 4713 BC becomes the year 1.
//
  y2 = y2 + 4713;

  c1 = ( y2 - 1 ) / 28;
  c2 = ( y2 - 1 ) / 19;
  c3 = ( y2 - 1 ) / 15;

  r1 = i4_wrap ( y2, 1, 28 );
  r2 = i4_wrap ( y2, 1, 19 );
  r3 = i4_wrap ( y2, 1, 15 );

  return;
}
//****************************************************************************80

int year_to_type_hebrew ( int y )

//****************************************************************************80
//
//// YEAR_TO_TYPE_HEBREW returns the type of a Hebrew year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2000
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 332.
//
//  Parameters:
//
//    Input, int Y, the Hebrew year.  
//    Nonpositive years are illegal input.
//
//    Output, int YEAR_TO_TYPE_HEBREW, the year type.
//    1, Common, Deficient, 12 months, 353 days;
//    2, Common, Regular, 12 months, 354 days;
//    3, Common, Abundant, 12 months, 355 days;
//    4, Embolismic, Deficient, 13 months, 383 days;
//    5, Embolismic, Regular, 13 months, 384 days;
//    6, Embolismic, Abundant, 13 months, 385 days.
//
{
  double jed;
  double jed2;
  int type;
  int year_length;

  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "YEAR_TO_TYPE_HEBREW - Fatal error!\n";
    cerr << "  Illegal input Y = 0.\n";
    exit ( 1 );
  }

  jed = new_year_to_jed_hebrew ( y );

  jed2 = new_year_to_jed_hebrew ( y + 1 );

  year_length = ( r8_nint ) ( jed2 - jed );

       if ( year_length == 353 )
  {
    type = 1;
  }
  else if ( year_length == 354 )
  {
    type = 2;
  }
  else if ( year_length == 355 )
  {
    type = 3;
  }
  else if ( year_length == 383 )
  {
    type = 4;
  }
  else if ( year_length == 384 )
  {
    type = 5;
  }
  else if ( year_length == 385 )
  {
    type = 6;
  }
  else
  {
    type = 0;
    cerr << "\n";
    cerr << "YEAR_TO_TYPE_HEBREW - Fatal error!\n";
    cerr << "  Computed an illegal type = " << type << "\n";
    exit ( 1 );
  }

  return type;
}
//****************************************************************************80

void yj_check_common ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    YJ_CHECK_COMMON checks a Common YJ date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, the YJ date.
//
{
//
//  Check the year.
//
  y_check_common ( y );
//
//  Make sure J is not too small or too big.
//
  j_borrow_common ( y, j );

  j_carry_common ( y, j );

  return;
}
//****************************************************************************80

void yj_check_english ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    YJ_CHECK_ENGLISH checks an English YJ date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, the YJ date.
//
{
//
//  Check the year.
//
  y_check_english ( y );
//
//  Make sure J is not too small or too big.
//
  j_borrow_english ( y, j );

  j_carry_english ( y, j );

  return;
}
//****************************************************************************80

void yj_check_gregorian ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    YJ_CHECK_GREGORIAN checks a Gregorian YJ date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, the YJ date.
//
{
//
//  Check the year.
//
  y_check_gregorian ( y );
//
//  Make sure J is not too small or too big.
//
  j_borrow_gregorian ( y, j );

  j_carry_gregorian ( y, j );

  return;
}
//****************************************************************************80

void yj_check_hebrew ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    YJ_CHECK_HEBREW checks a Hebrew YJ date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, the YJ date.
//
{
//
//  Check the year.
//
  y_check_hebrew ( y );
//
//  Make sure J is not too small or too big.
//
  j_borrow_hebrew ( y, j );

  j_carry_hebrew ( y, j );

  return;
}
//****************************************************************************80

void yj_check_islamic ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    YJ_CHECK_ISLAMIC checks an Islamic YJ date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, the YJ date.
//
{
//
//  Check the year.
//
  y_check_islamic ( y );
//
//  Make sure J is not too small or too big.
//
  j_borrow_islamic ( y, j );

  j_carry_islamic ( y, j );

  return;
}
//****************************************************************************80

void yj_check_julian ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    YJ_CHECK_JULIAN checks a Julian YJ date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, the YJ date.
//
{
//
//  Check the year.
//
  y_check_julian ( y );
//
//  Make sure J is not too small or too big.
//
  j_borrow_julian ( y, j );

  j_carry_julian ( y, j );

  return;
}
//****************************************************************************80

void yj_check_republican ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    YJ_CHECK_REPUBLICAN checks a Republican YJ date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, the YJ date.
//
{
//
//  Check the year.
//
  y_check_republican ( y );
//
//  Make sure J is not too small or too big.
//
  j_borrow_republican ( y, j );

  j_carry_republican ( y, j );

  return;
}
//****************************************************************************80

void yj_check_roman ( int &y, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    YJ_CHECK_ROMAN checks a Roman YJ date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, the YJ date.
//
{
//
//  Check the year.
//
  y_check_roman ( y );
//
//  Make sure J is not too small or too big.
//
  j_borrow_roman ( y, j );

  j_carry_roman ( y, j );

  return;
}
//****************************************************************************80

void yjf_check_common ( int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_CHECK_COMMON normalizes a Common YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, double &F, the YJF date.
//
{
  yj_check_common ( y, j );
//
//  Force the fraction to lie between 0 and 1.
//
  while ( f < 0.0 )
  {
    f = f + 1.0;
    j = j - 1;
  }

  while ( 1.0 <= f )
  {
    f = f - 1.0;
    j = j + 1;
  }

  return;
}
//****************************************************************************80

void yjf_check_english ( int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_CHECK_ENGLISH normalizes an English YJF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &J, double &F, the YJF date.
//
{
  yj_check_english ( y, j );
//
//  Force the fraction to lie between 0 and 1.
//
  while ( f < 0.0 )
  {
    f = f + 1.0;
    j = j - 1;
  }

  while ( 1.0 <= f )
  {
    f = f - 1.0;
    j = j + 1;
  }

  return;
}
//****************************************************************************80

char yjf_compare ( int y1, int j1, double f1, int y2, int j2, double f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_COMPARE compares two YJF dates.
//
//  Discussion:
//
//    The routine is "generic" and does not assume a particular calendar.
//    However, it does assume that the calendar dates are "normalized".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, 
//    the first YJF date.
//
//    Input, int Y2, J2, double F2, 
//    the second YJF date.
//
//    Output, char YJF_COMPARE:
//    '<' if date 1 precedes date 2;
//    '=' if date 1 equals date 2;
//    '>' if date 1 follows date 2;
//
{
  char cmp;

  if ( y1 < y2 )
  {
    cmp = '<';
  }
  else if ( y1 > y2 )
  {
    cmp = '>';
  }
  else
  {
    if ( j1 < j2 )
    {
      cmp = '<';
    }
    else if ( j1 > j2 )
    {
      cmp = '>';
    }
    else
    {
      if ( f1 < f2 )
      {
        cmp = '<';
      }
      else if ( f1 > f2 )
      {
        cmp = '>';
      }
      else
      {
        cmp = '=';
      }
    }
  }
  return cmp;
}
//****************************************************************************80

double yjf_dif_common ( int y1, int j1, double f1, int y2, int j2, double f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_DIF_COMMON computes day difference between two Common YJF dates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//   John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, real ( kind = 8 ) F1, 
//    the first YJF date.
//
//    Input, int Y2, J2, real ( kind = 8 ) F2, 
//    the second YJF date.
//
//    Output, double YJF_DIF_COMMON, the day difference between the two dates.
//
{
  double days;
  double jed1;
  double jed2;
//
//  Check the dates.
//
  yjf_check_common ( y1, j1, f1 );

  yjf_check_common ( y2, j2, f2 );

  jed1 = yjf_to_jed_common ( y1, j1, f1 );

  jed2 = yjf_to_jed_common ( y2, j2, f2 );

  days = jed2 - jed1;

  return days;
}
//****************************************************************************80

void yjf_swap ( int &y1, int &j1, double &f1, int &y2, int &j2, double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_SWAP swaps the data defining two YJF dates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 June 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y1, &J1, double &F1,
//    the first date.
//
//    Input/output, int &Y2, &J2, double &F2, 
//    the second date.
//
{
  double d;
  int i;

  i  = y1;
  y1 = y2;
  y2 = i;

  i  = j1;
  j1 = j2;
  j2 = i;

  d  = f1;
  f1 = f2;
  f2 = d;

  return;
}
//****************************************************************************80

double yjf_to_jed_common ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_COMMON converts a Common YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_COMMON, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  int j1;
  double jed;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yjf_check_common ( y1, j1, f1 );
//
//  Convert the input.
//
  yjf_to_ymdf_common ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_common ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

double yjf_to_jed_english ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_ENGLISH converts an English YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_ENGLISH, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  int j1;
  double jed;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yj_check_english ( y1, j1 );
//
//  Convert the input.
//
  yjf_to_ymdf_english ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_english ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

double yjf_to_jed_gregorian ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_GREGORIAN converts a Gregorian YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_GREGORIAN, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  int j1;
  double jed;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yj_check_gregorian ( y1, j1 );
//
//  Convert the input.
//
  yjf_to_ymdf_gregorian ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_gregorian ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

double yjf_to_jed_hebrew ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_HEBREW converts a Hebrew YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_HEBREW, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  int j1;
  double jed;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yj_check_hebrew ( y1, j1 );
//
//  Convert the input.
//
  yjf_to_ymdf_hebrew ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_hebrew ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

double yjf_to_jed_islamic_a ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_ISLAMIC_A converts an Islamic-A YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_ISLAMIC_A, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  int j1;
  double jed;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yj_check_islamic ( y1, j1 );
//
//  Convert the input.
//
  yjf_to_ymdf_islamic ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_islamic_a ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

double yjf_to_jed_islamic_b ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_ISLAMIC_B converts an Islamic-B YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_ISLAMIC_B, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  int j1;
  double jed;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yj_check_islamic ( y1, j1 );
//
//  Convert the input.
//
  yjf_to_ymdf_islamic ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_islamic_b ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

double yjf_to_jed_julian ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_JULIAN converts a Julian YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_JULIAN, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  int j1;
  double jed;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yj_check_julian ( y1, j1 );
//
//  Convert the input.
//
  yjf_to_ymdf_julian ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_julian ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

double yjf_to_jed_republican ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_REPUBLICAN converts a Republican YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_REPUBLICAN, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  int j1;
  double jed;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yj_check_republican ( y1, j1 );
//
//  Convert the input.
//
  yjf_to_ymdf_republican ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_republican ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

double yjf_to_jed_roman ( int y, int j, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_JED_ROMAN converts a Roman YJF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, J, double F, the YJF date.
//
//    Output, double YJF_TO_JED_ROMAN, the Julian Ephemeris Date.
//
{
  int d2;
  double f1;
  double f2;
  double jed;
  int j1;
  int m2;
  int y1;
  int y2;
//
//  Copy the input.
//
  y1 = y;
  j1 = j;
  f1 = f;
//
//  Check the input.
//
  yj_check_roman ( y1, j1 );
//
//  Convert the input.
//
  yjf_to_ymdf_roman ( y1, j1, f1, y2, m2, d2, f2 );

  jed = ymdf_to_jed_roman ( y2, m2, d2, f2 );

  return jed;
}
//****************************************************************************80

void yjf_to_ymdf_common ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_YMDF_COMMON converts a Common date from YJF to YMDF format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, the YJF date.
//
//    Output, int &Y2, &M2, &D2, double &F2, the YMDF date.
//
{
  int j2;
//
//  Copy the input.
//
  y2 = y1;
  j2 = j1;
  f2 = f1;
//
//  Check the input.
//
  yjf_check_common ( y2, j2, f2 );
//
//  Convert the input.
//
  d2 = j2;
  m2 = 1;

  day_borrow_common ( y2, m2, d2 );

  day_carry_common ( y2, m2, d2 );

  return;
}
//****************************************************************************80

void yjf_to_ymdf_english ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_YMDF_ENGLISH converts an English date from YJF to YMDF format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, the YJF date.
//
//    Output, int &Y2, &M2, &D2, double &F2, the YMDF date.
//
{
  int j2;
//
//  Copy the input.
//
  y2 = y1;
  j2 = j1;
  f2 = f1;
//
//  Check the input.
//
  yj_check_english ( y2, j2 );
//
//  Convert the input.
//
  m2 = 1;
  d2 = j2;

  day_borrow_english ( y2, m2, d2 );

  day_carry_english ( y2, m2, d2 );

  return;
}
//****************************************************************************80

void yjf_to_ymdf_gregorian ( int y1, int j1, double f1, int &y2, int &m2, 
  int &d2, double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_YMDF_GREGORIAN converts a Gregorian date from YJF to YMDF format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, the YJF date.
//
//    Output, int &Y2, &M2, &D2, double &F2, the YMDF date.
//
{
  int j2;
//
//  Copy the input.
//
  y2 = y1;
  j2 = j1;
  f2 = f1;
//
//  Check the input.
//
  yj_check_gregorian ( y2, j2 );
//
//  Convert the input.
//
  m2 = 1;
  d2 = j2;

  day_borrow_gregorian ( y2, m2, d2 );

  day_carry_gregorian ( y2, m2, d2 );

  return;
}
//****************************************************************************80

void yjf_to_ymdf_hebrew ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_YMDF_HEBREW converts a YJF to YMDF date, both in the Hebrew calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, the YJF date.
//
//    Output, int &Y2, &M2, &D2, double &F2, the YMDF date.
//
{
  int j2;
//
//  Copy the input.
//
  y2 = y1;
  j2 = j1;
  f2 = f1;
//
//  Check the input.
//
  yj_check_hebrew ( y2, j2 );
//
//  Convert the input.
//
  m2 = 1;
  d2 = j2;

  day_borrow_hebrew ( y2, m2, d2 );

  day_carry_hebrew ( y2, m2, d2 );

  return;
}
//****************************************************************************80

void yjf_to_ymdf_islamic ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_YMDF_ISLAMIC: YJF to YMDF date, both in the Islamic calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, the YJF date.
//
//    Output, int &Y2, &M2, &D2, double &F2, the YMDF date..
//
{
  int j2;
//
//  Copy the input.
//
  y2 = y1;
  j2 = j1;
  f2 = f1;
//
//  Check the input.
//
  yj_check_islamic ( y2, j2 );
//
//  Convert the input.
//
  m2 = 1;
  d2 = j2;

  day_borrow_islamic ( y2, m2, d2 );

  day_carry_islamic ( y2, m2, d2 );

  return;
}
//****************************************************************************80

void yjf_to_ymdf_julian ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_YMDF_JULIAN converts a YJF to YMDF date, both in the Julian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//   23 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, the YJF date.
//
//    Output, int &Y2, &M2, &D2, double &F2, the YMDF date.
//
{
  int j2;
//
//  Copy the input.
//
  y2 = y1;
  j2 = j1;
  f2 = f1;
//
//  Check the input.
//
  yj_check_julian ( y2, j2 );
//
//  Convert the input.
//
  m2 = 1;
  d2 = j2;

  day_borrow_julian ( y2, m2, d2 );

  day_carry_julian ( y2, m2, d2 );

  return;
}
//****************************************************************************80

void yjf_to_ymdf_republican ( int y1, int j1, double f1, int &y2, int &m2, 
  int &d2, double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_YMDF_REPUBLICAN: YJF to YMDF date in the Republican calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, the YJF date.
//
//    Output, int &Y2, &M2, &D2, double &F2, the YMDF date.
//
{
  int j2;
//
//  Copy the input.
//
  y2 = y1;
  j2 = j1;
  f2 = f1;
//
//  Check the input.
//
  yj_check_republican ( y2, j2 );
//
//  Convert the input.
//
  m2 = 1;
  d2 = j2;

  day_borrow_republican ( y2, m2, d2 );

  day_carry_republican ( y2, m2, d2 );

  return;
}
//****************************************************************************80

void yjf_to_ymdf_roman ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_TO_YMDF_ROMAN converts a YJF to YMDF date in the Roman calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1, the YJF date.
//
//    Output, int &Y2, &M2, &D2, double &F2, the YMDF date.
//
{
  int j2;
//
//  Copy the input.
//
  y2 = y1;
  j2 = j1;
  f2 = f1;
//
//  Check the input.
//
  yj_check_roman ( y2, j2 );
//
//  Convert the input.
//
  m2 = 1;
  d2 = j2;

  day_borrow_roman ( y2, m2, d2 );

  day_carry_roman ( y2, m2, d2 );

  return;
}
//****************************************************************************80

void yjf_uniform_common ( int y1, int j1, double f1, int y2, int j2, double f2, 
  int &seed, int &y, int &j, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    YJF_UNIFORM_COMMON picks a random Common YJF date between two given dates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, J1, double F1,
//    the first YJF date.
//
//    Input, int Y2, J2, double F2,
//    the second YJF date.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, int &Y, &J, double &F, the random 
//    YJF date.
//
{
  double jed;
  double jed1;
  double jed2;

  jed1 = yjf_to_jed_common ( y1, j1, f1 );
  jed2 = yjf_to_jed_common ( y2, j2, f2 );

  jed = r8_uniform_ab ( jed1, jed2, seed );

  jed_to_yjf_common ( jed, y, j, f );

  return;
}
//****************************************************************************80

void ym_check_alexandrian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_ALEXANDRIAN checks an Alexandrian YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by the number of months in the PREVIOUS year, and the year is
//    decremented by 1.
//
//    If the month is greater than the number of months in the CURRENT year,
//    then the month is decremented by the number of months in the CURRENT year,
//    and the year incremented by 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//
{
//
//  Check the year.
//
  y_check_alexandrian ( y );
//
//  Make sure the month isn't too small or too big.
//
  month_borrow_alexandrian ( y, m );

  month_carry_alexandrian ( y, m );

  return;
}
//****************************************************************************80

void ym_check_bahai ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_BAHAI checks a Bahia YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by the number of months in the PREVIOUS year, and the year is
//    decremented by 1.
//
//    If the month is greater than the number of months in the CURRENT year,
//    then the month is decremented by the number of months in the CURRENT year,
//    and the year incremented by 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//
{
//
//  Check the year.
//
  y_check_bahai ( y );
//
//  Make sure the month isn't too small or too big.
//
  month_borrow_bahai ( y, m );

  month_carry_bahai ( y, m );

  return;
}
//****************************************************************************80

void ym_check_common ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_COMMON checks a Common YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by 12, and the year decremented by 1, repeatedly, until
//    the month is greater than or equal to 1.
//
//    If the month is greater than 12, then the month is decremented
//    by 12, and the year incremented by 1, repeatedly, until the
//    month is less than or equal to 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//
{
//
//  Check the year.
//
  y_check_common ( y );
//
//  Make sure the month isn't too small or too big.
//
  month_borrow_common ( y, m );

  month_carry_common ( y, m );

  return;
}
//****************************************************************************80

void ym_check_eg_civil ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_EG_CIVIL checks an Egyptian Civil YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by 12, and the year decremented by 1, repeatedly, until
//    the month is greater than or equal to 1.
//
//    If the month is greater than 12, then the month is decremented
//    by 12, and the year incremented by 1, repeatedly, until the
//    month is less than or equal to 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//
{
//
//  Check the year.
//
  y_check_eg_civil ( y );
//
//  Make sure the month isn't too small or too big.
//
  month_borrow_eg_civil ( y, m );

  month_carry_eg_civil ( y, m );

  return;
}
//****************************************************************************80

void ym_check_english ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_ENGLISH checks an English YM date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  y_check_english ( y );

  month_borrow_english ( y, m );

  month_carry_english ( y, m );

  return;
}
//****************************************************************************80

void ym_check_gregorian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_GREGORIAN checks a Gregorian YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by 12, and the year decremented by 1, repeatedly, until
//    the month is greater than or equal to 1.
//
//    If the month is greater than 12, then the month is decremented
//    by 12, and the year incremented by 1, repeatedly, until the
//    month is less than or equal to 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  y_check_gregorian ( y );

  month_borrow_gregorian ( y, m );

  month_carry_gregorian ( y, m );

  return;
}
//****************************************************************************80

void ym_check_hebrew ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_HEBREW checks a Hebrew YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by 12, and the year decremented by 1, repeatedly, until
//    the month is greater than or equal to 1.
//
//    If the month is greater than 12, then the month is decremented
//    by 12, and the year incremented by 1, repeatedly, until the
//    month is less than or equal to 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//
{
//
//  Check the year.
//
  y_check_hebrew ( y );
//
//  Make sure the month isn't too small or too big.
//
  month_borrow_hebrew ( y, m );

  month_carry_hebrew ( y, m );

  return;
}
//****************************************************************************80

void ym_check_islamic ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_ISLAMIC checks an Islamic YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by 12, and the year decremented by 1, repeatedly, until
//    the month is greater than or equal to 1.
//
//    If the month is greater than 12, then the month is decremented
//    by 12, and the year incremented by 1, repeatedly, until the
//    month is less than or equal to 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//
{
//
//  Check the year.
//
  y_check_islamic ( y );
//
//  Make sure the month isn't too small or too big.
//
  month_borrow_islamic ( y, m );

  month_carry_islamic ( y, m );

  return;
}
//****************************************************************************80

void ym_check_julian ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_JULIAN checks a Julian YM date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, the YM date.
//
{
  y_check_julian ( y );

  month_borrow_julian ( y, m );

  month_carry_julian ( y, m );

  return;
}
//****************************************************************************80

void ym_check_republican ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_REPUBLICAN checks a Republican YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by 12, and the year decremented by 1, repeatedly, until
//    the month is greater than or equal to 1.
//
//    If the month is greater than 12, then the month is decremented
//    by 12, and the year incremented by 1, repeatedly, until the
//    month is less than or equal to 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//
{
//
//  Check the year.
//
  y_check_republican ( y );
//
//  Make sure the month isn't too small or too big.
//
  month_borrow_republican ( y, m );

  month_carry_republican ( y, m );

  return;
}
//****************************************************************************80

void ym_check_roman ( int &y, int &m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_CHECK_ROMAN checks a Roman YM date.
//
//  Discussion:
//
//    If the month is less than 1, then the month is incremented
//    by 12, and the year decremented by 1, repeatedly, until
//    the month is greater than or equal to 1.
//
//    If the month is greater than 12, then the month is decremented
//    by 12, and the year incremented by 1, repeatedly, until the
//    month is less than or equal to 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, the YM date.
//
{
//
//  Check the year.
//
  y_check_roman ( y );
//
//  Make sure the month isn't too small or too big.
//
  month_borrow_roman ( y, m );

  month_carry_roman ( y, m );

  return;
}
//****************************************************************************80

double ym_to_decimal ( int y, int m )

//****************************************************************************80
//
//  Purpose:
//
//    YM_TO_DECIMAL converts a Y/M date to a Decimal YM date.
//
//  Discussion:
//
//    Each month is take to be 1/12 of a year long, and the decimal value
//    is returned for the middle of the month.
//
//    1980 January  => 1980.04
//    1980 February => 1980.12
//    1980 March    => 1980.21
//    1980 December => 1980.96
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, the YM date.
//
//    Output, double YM_TO_DECIMAL, the Decimal date.
//
{
  double ym;

  ym =  ( double ) ( y ) + ( double ) ( 2 * m - 1 ) / 24.0;

  return ym;
}
//****************************************************************************80

void ymd_check_alexandrian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_ALEXANDRIAN checks an Alexandrian YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 November 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the year.
//
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "YMD_CHECK_ALEXANDRIAN - Fatal error!\n";
    cerr << "  Y <= 0 is illegal.\n";
    exit ( 1 );
  }
//
//  Check the month.
//
  ym_check_alexandrian ( y, m );
//
//  Check the day.
//
  day_borrow_alexandrian ( y, m, d );

  day_carry_alexandrian ( y, m, d );

  return;
}
//****************************************************************************80

void ymd_check_common ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_COMMON checks a Common YMD date.
//
//  Discussion:
//
//    Certain simple errors in dates will be corrected, such as
//      "31 September 1996"
//    which will become
//      "1 October 1996".
//
//    The routine also knows that in the Common calendar, the dates
//    5 October 1582 through 14 October 1582 are illegal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date,
//    which may be corrected if necessary and possible.
//
{
//
//  Check the year.
//
  if ( y == 0 )
  {
    cerr << "\n";
    cerr << "YMD_CHECK_COMMON - Fatal error!\n";
    cerr << "  Y = 0 is illegal.\n";
    exit ( 1 );
  }
//
//  Check the month.
//
  month_borrow_common ( y, m );

  month_carry_common ( y, m );
//
//  Check the day.
//
  day_borrow_common ( y, m, d );

  day_carry_common ( y, m, d );
//
//  Now make sure that the date does not fall in the
//  Julian-to-Gregorian calendar switchover limbo.
//
  if ( y == 1582 )
  {
    if ( m == 10 )
    {
      if ( 5 <= d && d <= 14 )
      {
        cerr << "\n";
        cerr << "YMD_CHECK_COMMON - Fatal error!\n";
        cerr << "  Illegal date in Julian-to-Gregorian transition.\n";
        exit ( 1 );
      }
    }
  }
  return;
}
//****************************************************************************80

void ymd_check_eg_civil ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_EG_CIVIL checks an Egyptian Civil YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the year.
//
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "YMD_CHECK_EG_CIVIL - Fatal error!\n";
    cerr << "  Y <= 0 is illegal.\n";
    exit ( 1 );
  }
//
//  Check the month.
//
  ym_check_eg_civil ( y, m );
//
//  Check the day.
//
  day_borrow_eg_civil ( y, m, d );

  day_carry_eg_civil ( y, m, d );

  return;
}
//****************************************************************************80

void ymd_check_english ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_ENGLISH checks an English YMD date.
//
//  Discussion:
//
//    Certain simple errors in dates will be corrected, such as
//      "31 September 1996"
//    which will become
//      "1 October 1996".
//
//    The routine also knows that in the English calendar, the dates
//    3 September 1752 through 13 September 1752 are illegal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the month.
//
  ym_check_english ( y, m );
//
//  Check the day.
//
  day_borrow_english ( y, m, d );

  day_carry_english ( y, m, d );
//
//  Now make sure that the date does not fall in the
//  Julian-to-Gregorian calendar switchover limbo.
//
  if ( y == 1752 )
  {
    if ( m == 9 )
    {
      if ( 3 <= d && d <= 13 )
      {
        cerr << "\n";
        cerr << "YMD_CHECK_ENGLISH - Fatal error!\n";
        cerr << "  Illegal date: " << y << "  " << m << "  " << d << "\n";
        exit ( 1 ); 
      }
    }
  }

  return;
}
//****************************************************************************80

void ymd_check_gregorian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_GREGORIAN checks a Gregorian YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the month.
//
  ym_check_gregorian ( y, m );
//
//  Check the day.
//
  day_borrow_gregorian ( y, m, d );

  day_carry_gregorian ( y, m, d );

  return;
}
//****************************************************************************80

void ymd_check_hebrew ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_HEBREW checks a Hebrew YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the year.
//
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "YMD_CHECK_HEBREW - Fatal error!\n";
    cerr << "  Y <= 0 is illegal.\n";
    exit ( 1 );
  }
//
//  Check the month.
//
  ym_check_hebrew ( y, m );
//
//  Check the day.
//
  day_borrow_hebrew ( y, m, d );

  day_carry_hebrew ( y, m, d );

  return;
}
//****************************************************************************80

void ymd_check_islamic ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_ISLAMIC checks an Islamic YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the year.
//
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "YMD_CHECK_ISLAMIC - Fatal error!\n";
    cerr << "  Y <= 0 is illegal.\n";
    exit ( 1 );
  }
//
//  Check the month.
//
  ym_check_islamic ( y, m );
//
//  Check the day.
//
  day_borrow_islamic ( y, m, d );

  day_carry_islamic ( y, m, d );

  return;
}
//****************************************************************************80

void ymd_check_julian ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_JULIAN checks a Julian YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the month.
//
  ym_check_julian ( y, m );
//
//  Check the day.
//
  day_borrow_julian ( y, m, d );

  day_carry_julian ( y, m, d );

  return;
}
//****************************************************************************80

void ymd_check_republican ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_REPUBLICAN checks a Republican YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the year.
//
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "YMD_CHECK_REPUBLICAN - Fatal error!\n";
    cerr << "  Y <= 0 is illegal.\n";
    exit ( 1 );
  }
//
//  Check the month.
//
  ym_check_republican ( y, m );
//
//  Check the day.
//
  day_borrow_republican ( y, m, d );

  day_carry_republican ( y, m, d );

  return;
}
//****************************************************************************80

void ymd_check_roman ( int &y, int &m, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_CHECK_ROMAN checks a Roman YMD date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, the YMD date, which may
//    be corrected if necessary and possible.
//
{
//
//  Check the year.
//
  if ( y <= 0 )
  {
    cerr << "\n";
    cerr << "YMD_CHECK_ROMAN - Fatal error!\n";
    cerr << "  Y <= 0 is illegal.\n";
    exit ( 1 );
  }
//
//  Check the month.
//
  ym_check_roman ( y, m );
//
//  Check the day.
//
  day_borrow_roman ( y, m, d );

  day_carry_roman ( y, m, d );

  return;
}
//****************************************************************************80

char ymd_compare ( int y1, int m1, int d1, int y2, int m2, int d2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_COMPARE compares two YMD dates.
//
//  Discussion:
//
//    The comparison should work for a pair of dates in any calendar.
//
//    No check is made that the dates are actually legitimate.  It is
//    assumed that the calling routine has already ensured that the
//    dates are properly "normalized".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, int M1, int D1, the first YMD date.
//
//    Input, int Y2, int M2, int D2, the second YMD date.
//
//    Output, char YMD_COMPARE:
//    '<' if date 1 precedes date 2;
//    '=' if date 1 equals date 2;
//    '>' if date 1 follows date 2;
//
{
  char cmp;

  cmp = '?';

  if ( y1 < y2 )
  {
    cmp = '<';
  }
  else if ( y1 > y2 )
  {
    cmp = '>';
  }
  else
  {
    if ( m1 < m2 )
    {
      cmp = '<';
    }
    else if ( m1 > m2 )
    {
      cmp = '>';
    }
    else
    {
      if ( d1 < d2 )
      {
        cmp = '<';
      }
      else if ( d1 > d2 )
      {
        cmp = '>';
      }
      else
      {
        cmp = '=';
      }
    }
  }
  return cmp;
}
//****************************************************************************80

int ymd_dif_common ( int y1, int m1, int d1, int y2, int m2, int d2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_DIF_COMMON gets the day difference between two Common YMD dates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, the first YMD date.
//
//    Input, int Y2, M2, D2, the second YMD date.
//
//    Output, int YMD_DIF_COMMON, the number of days between the dates.
//
{
  int days;
  double jed1;
  double jed2;

  days = 0;
//
//  Check the dates.
//
  ymd_check_common ( y1, m1, d1 );

  ymd_check_common ( y2, m2, d2 );

  jed1 = ymd_to_jed_common ( y1, m1, d1 );

  jed2 = ymd_to_jed_common ( y2, m2, d2 );

  days = r8_round ( jed2 - jed1 );

  return days;
}
//****************************************************************************80

void ymd_inc_ymd_common ( int y1, int m1, int d1, int yn, int mn, int dn, 
  int &y2, int &m2, int &d2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_INC_YMD_COMMON increments a Common YMD date by a YMD increment.
//
//  Discussion:
//
//    You often see on old gravestones statements like
//
//      "Joe Blow died on May 8 1784 aged 38 Years, 7 Months and 5 Days."
//
//    It's not exactly clear how to interpret such a statement, since
//    we can't actually convert 38 Years, 7 Months and 5 Days to a number
//    of days.  (Years and months vary in their day length).  However,
//    we can assume that what was meant was, if you take the year, month
//    and day of Joe Blow's birthday, and you:
//
//      add 38 to the year,
//      add 7 to the month, and if you go past December, subtract 12 and
//        increment the year,
//      add 5 to the day, and if you go past the length of the month,
//        increment the month and decrement the day appropriately.
//
//    Notice, in particular, that if you do the operations in the reverse
//    order, you may get a different answer, since what you do with a large
//    day value depends on the month you assume you are working in.
//
//    Just warning you that this is a poorly posed problem.
//
//    Thanks to Charlie Cullen for pointing out this little problem to me.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, the YMD date.
//
//    Input, int YN, MN, DN, the increment to the YMD date.
//
//    Output, int &Y2, &M2, &D2, the incremented YMD date.
//
{
  double f1;
  double f2;
  double fn;

  f1 = 0.0;
  fn = 0.0;

  y2 = y1 + yn;
  m2 = m1 + mn;
  d2 = d1 + dn;
  f2 = f1 + fn;

  ymdf_check_common ( y2, m2, d2, f2 );

  return;
}
//****************************************************************************80

double ymd_to_decimal ( int y, int m, int d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_TO_DECIMAL converts a Y/M/D date to a Decimal Y.F date.
//
//  Discussion:
//
//    The day is assumed to be at noon.  In other words, 1983 January 1st has
//    a decimal value of 1983 + 0.5 / 365.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, the YMD date.
//
//    Output, double YMD_TO_DECIMAL, the Decimal date.
//
{
  int d2;
  int day_max;
  int days;
  double f;
  int ierror;
  int m2;
  double yf;
//
//  How many days between January 1st and day D?
//
  m2 = 1;
  d2 = 1;
  days = ymd_dif_common ( y, m2, d2, y, m, d );
//
//  How many days in this year total?
//
  day_max = year_length_common ( y );
//
//  The decimal part of the year is ( D + 0.5 ) / DMAX.
//
  f = ( ( double ) ( days ) + 0.5 ) 
      / ( double ) ( day_max );

  yf = ( double ) ( y ) + f;

  return yf;
}
//****************************************************************************80

double ymd_to_jed_common ( int y, int m, int d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_TO_JED_COMMON converts a Common YMD date to a JED.
//
//  Discussion:
//
//    The "common" calendar is meant to be the calendar which is Julian up to
//    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
//
//    The Julian Ephemeris Date is essentially a count of the number
//    of days that have elapsed since noon, 1 January 4713 BC, at
//    Greenwich, England.  Strictly speaking, the Julian Ephemeris Date
//    is counted from noon, and thus day "0" began at noon on 1 January 4713 BC,
//    and ended at noon on 2 January 4713 BC.
//
//    The Julian Ephemeris Date was devised by Joseph Scaliger in 1583.
//
//    The Julian Ephemeris Date has been adopted by astronomers as
//    a convenient reference for dates.
//
//  Example:
//
//       Y   M     D         JED
//    --------------     -------
//    BC 4713 Jan  1           0
//    AD 1968 May 23     2440000
//    AD 1984 Dec 31     2446065
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, the YMD date.
//
//    Output, double YMD_TO_JED, the Julian Ephemeris Date.
//
{
  char cmp;
  int d1;
  int d2;
  double f;
  double jed;
  int m1;
  int m2;
  int y1;
  int y2;
//
//  Copy the month and year.
//
  y1 = y;
  m1 = m;
  d1 = d;

  ymd_check_common ( y1, m1, d1 );

  y2 = 1582;
  m2 = 10;
  d2 = 4+1;

  cmp = ymd_compare ( y1, m1, d1, y2, m2, d2 );

  if ( cmp == '<' )
  {
    jed = ymd_to_jed_julian ( y1, m1, d1 );
    return jed;
  }
//
//  Use the Gregorian calendar for dates strictly after 1752/9/13.
//
  y2 = 1582;
  m2 = 10;
  d2 = 15-1;

  cmp = ymd_compare ( y1, m1, d1, y2, m2, d2 );

  if ( cmp == '>' )
  {
    jed = ymd_to_jed_gregorian ( y1, m1, d1 );
    return jed;
  }

  jed = -1.0;
  cerr << "\n";
  cerr << "YMD_TO_JED_COMMON - Error!\n";
  cerr << "  Illegal date!\n";
  exit ( 1 );
}
//****************************************************************************80

double ymd_to_jed_gregorian ( int y, int m, int d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_TO_JED_GREGORIAN converts a Gregorian YMD date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, the YMD date.
//
//    Output, double YMD_TO_JED_GREGORIAN, the corresponding JED.
//
{
  int d_prime;
  int g;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y2;
  int y_prime;
//
//  Check the date.
//
  ymd_check_gregorian ( y, m, d );
//
//  Account for the missing year 0 by moving negative years up one.
//
  y2 = y_common_to_astronomical ( y );
//
//  Convert the calendar date to a computational date.
//
  y_prime = y2 + 4716 - ( 14 - m ) / 12;
  m_prime = ( ( m + 9 ) % 12 );
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = ( 153 * m_prime + 2 ) / 5;

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 38;

  jed = ( double ) ( j1 + j2 + d_prime - 1401 - g ) - 0.5;

  return jed;
}
//****************************************************************************80

double ymd_to_jed_julian ( int y, int m, int d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_TO_JED_JULIAN converts a Julian YMD date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2001
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, the YMD date.
//
//    Output, double YMD_TO_JED_JULIAN, the Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y2;
  int y_prime;
//
//  Check the date.
//
  ymd_check_julian ( y, m, d );
//
//  Account for the missing year 0 by moving negative years up one.
//
  y2 = y_common_to_astronomical ( y );
//
//  Convert the calendar date to a computational date.
//
  y_prime = y2 + 4716 - ( 14 - m ) / 12;
  m_prime = ( ( m + 9 ) % 12 );
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = ( 153 * m_prime + 2 ) / 5;

  jed = ( double ) ( j1 + j2 + d_prime - 1401 ) - 0.5;

  return jed;
}
//****************************************************************************80

string ymd_to_s_common ( int y, int m, int d )

//****************************************************************************80
//
//  Purpose:
//
//    YMD_TO_S_COMMON writes a Common YMD date into a string.
//
//  Format:
//
//    CE YYYY/MM/DD
//    BCE YYYY/MM/DD
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, the YMD date.
//
//    Output, string S, a representation of the date.
//
{
  string s;
  char s_char[80];

  if ( 0 <= y )
  {
    sprintf ( s_char, "CE %d/%02d/%02d", y, m, d );
  }
  else
  {
    sprintf ( s_char, "BCE %d/%02d/%02d", -y, m, d );
  }

  s = string ( s_char );

  return s;
}
//****************************************************************************80

void ymdf_check_common ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_CHECK_COMMON checks a Common YMDF date.
//
//  Discussion:
//
//    Certain simple errors in dates will be corrected, such as
//      "31 September 1996"
//    which will become
//      "1 October 1996".
//
//    The routine also knows that in the Common calendar, the dates
//    5 October 1582 through 14 October 1582 are illegal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, &M, &D, double &F, the 
//    YMDF date, which may be corrected if necessary and possible.
//
{
  ymd_check_common ( y, m, d );

  frac_borrow_common ( y, m, d, f );

  frac_carry_common ( y, m, d, f );

  return;
}
//****************************************************************************80

void ymdf_check_julian ( int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_CHECK_JULIAN checks a Julian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &Y, int &M, int &D, double &F, the 
//    YMDF date, which may be corrected if necessary and possible.
//
{
  ymd_check_julian ( y, m, d );

  frac_borrow_julian ( y, m, d, f );

  frac_carry_julian ( y, m, d, f );

  return;
}
//****************************************************************************80

char ymdf_compare ( int y1, int m1, int d1, double f1, int y2, int m2, int d2, 
  double f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_COMPARE compares two YMDF dates.
//
//  Discussion:
//
//    The comparison should work for a pair of dates in any calendar.
//
//    No check is made that the dates are actually legitimate.  It is
//    assumed that the calling routine has already ensured that the
//    dates are properly "normalized".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, int M1, int D1, double F1, the 
//    first YMDF date.
//
//    Input, int Y2, int M2, int D2, double F2, the 
//    second YMDF date.
//
//    Output, char YMDF_COMPARE:
//    '<' if date 1 precedes date 2;
//    '=' if date 1 equals date 2;
//    '>' if date 1 follows date 2;
//
{
  char cmp;

  cmp = '?';

  if ( y1 < y2 )
  {
    cmp = '<';
  }
  else if ( y1 > y2 )
  {
    cmp = '>';
  }
  else
  {
    if ( m1 < m2 )
    {
      cmp = '<';
    }
    else if ( m1 > m2 )
    {
      cmp = '>';
    }
    else
    {
      if ( d1 < d2 )
      {
        cmp = '<';
      }
      else if ( d1 > d2 )
      {
        cmp = '>';
      }
      else
      {
        if ( f1 < f2 )
        {
          cmp = '<';
        }
        else if ( f1 > f2 )
        {
          cmp = '>';
        }
        else
        {
          cmp = '=';
        }
      }
    }
  }
  return cmp;
}
//****************************************************************************80

void ymdf_next_common ( int y1, int m1, int d1, double f1, int &y2, int &m2, 
  int &d2, double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_NEXT_COMMON returns the Common YMDF date of the next day.
//
//  Discussion:
//
//    The routine knows that in the Common calendar, the day after
//    4 October 1582 was 15 October 1582.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1, 
//    the YMDF date.
//
//    Output, int &Y2, &M2, &D2, double &F2,
//    tomorrow's YMDF date.
//
{
  y2 = y1;
  m2 = m1;
  d2 = d1 + 1;
  f2 = f1;

  day_carry_common ( y2, m2, d2 );

  return;
}
//****************************************************************************80

double ymdf_to_jed_alexandrian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ALEXANDRIAN converts an Alexandrian YMDF date to a JED.
//
//  Discussion:
//
//    This code needs to be adjusted to fit the Alexandrian model.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ALEXANDRIAN, the corresponding 
//    Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 4690 - ( 13 - m ) / 13;
  m_prime = ( m + 12 ) % 13;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  jed = ( double ) 
    ( ( 1461 * y_prime ) / 4 + 30 * m_prime + d_prime - 124 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_armenian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ARMENIAN converts an Armenian YMDF date to a JED.
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
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ARMENIAN, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 5268 - ( 13 - m ) / 13;
  m_prime = ( m + 12 ) % 13;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  jed = ( double ) ( 365 * y_prime + 30 * m_prime + d_prime - 317 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_bahai ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_BAHAI converts a Bahai YMDF date to a JED.
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
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_BAHAI, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  int g;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 6560 - ( 39 - m ) / 20;
  m_prime = m % 20;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = 19 * m_prime;

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 50;
  jed = ( double ) ( j1 + j2 + d_prime - 1412 - g ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_common ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_COMMON converts a Common YMDF date to a JED.
//
//  Discussion:
//
//    The "common" calendar is meant to be the calendar which is Julian up to
//    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
//
//    The Julian Ephemeris Date is essentially a count of the number
//    of days that have elapsed since noon, 1 January 4713 BC, at
//    Greenwich, England.  Strictly speaking, the Julian Ephemeris Date
//    is counted from noon, and thus day "0" began at noon on 1 January 4713 BC,
//    and ended at noon on 2 January 4713 BC.
//
//    The Julian Ephemeris Date was devised by Joseph Scaliger in 1583.
//
//    The Julian Ephemeris Date has been adopted by astronomers as
//    a convenient reference for dates.
//
//  Example:
//
//     Y   M     D         JED
//    --------------     -------
//    BC 4713 Jan  1           0
//    AD 1968 May 23     2440000
//    AD 1984 Dec 31     2446065
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_COMMON, the Julian Ephemeris Date.
//
{
  char cmp;
  int d1;
  int d2;
  double f1;
  double f2;
  int ierror;
  double jed;
  int m1;
  int m2;
  int y1;
  int y2;
//
//  Copy the month and year.
//
  y1 = y;
  m1 = m;
  d1 = d;
  f1 = f;

  ymdf_check_common ( y1, m1, d1, f1 );

  y2 = 1582;
  m2 = 10;
  d2 = 4+1;
  f2 = 0.0;

  cmp = ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2 );

  if ( cmp == '<' )
  {
    jed = ymdf_to_jed_julian ( y1, m1, d1, f1 );
    return jed;
  }
//
//  Use the Gregorian calendar for dates strictly after 1752/9/13.
//
  y2 = 1582;
  m2 = 10;
  d2 = 15-1;
  f2 = 0.0;

  cmp = ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2 );

  if ( cmp == '>' )
  {
    jed = ymdf_to_jed_gregorian ( y1, m1, d1, f1 );
    return jed;
  }

  jed = -1.0;
  cout << "\n";
  cout << "YMDF_TO_JED_COMMON - Error!\n";
  cout << "  Illegal date!\n";

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_coptic ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_COPTIC converts a Coptic YMDF date to a JED.
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
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_COPTIC, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int m_prime;
  int y_prime;

//  Convert the calendar date to a computational date.
//
  y_prime = y + 4996 - ( 13 - m ) / 13;
  m_prime = ( m + 12 ) % 13;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  jed = ( double ) ( 
    ( 1461 * y_prime ) / 4 + 30 * m_prime + d_prime - 124 ) - 0.5;

  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_eg_civil ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_EG_CIVIL converts an Egyptian Civil YMDF date to a JED.
//
//  Discussion:
//
//    The Egyptian Civil calendar used a year of 365 days.  The year comprised
//    12 months of 30 days, with 5 epagomenal days occurring at the end of
//    the year.  Since the observed year is about 365.25 days long, and no
//    attempt was made to adjust the Egyptian Civil year to the observed year,
//    the calendar dates gradually drifted with respect to the observed dates.
//
//    The epoch or first day of the Egyptian Civil calendar is taken as
//    JED = 1448638.5.
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
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_EG_CIVIL, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 3968 - ( 13 - m ) / 13;
  m_prime = ( m + 12 ) % 13;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  jed = ( double ) ( 365 * y_prime + 30 * m_prime + d_prime - 47 + 1 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_eg_lunar ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_EG_LUNAR converts an Egyptian Lunar YMDF date to a JED.
//
//  Discussion:
//
//    Count
//      the days up to the day before the start of the calendar,
//      the days in the current month,
//      the 29 days guaranteed in the previous months of this year,
//      the (months/2) 30th days in the previous months of this year,
//      the 354 days guaranteed in each of the previous years,
//      the extra leap days in the preceding years,
//      the extra 30 days in the leap months in the preceding years.
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
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_EG_LUNAR, the corresponding Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch;

  jed_epoch = epoch_to_jed_eg_lunar ( );

  jed = jed_epoch + ( double ) ( - 1 + d + 29 * ( m - 1 ) + ( m - 1 ) / 2 
    + 354 * ( y - 1 ) + ( y - 1 ) / 5 
    + 30 * ( ( ( y - 1 ) / 25 ) * 9 + ( ( ( y - 1 ) % 25 ) + 2 ) / 3 ) );

  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_english ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ENGLISH converts an English YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, int M, int D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ENGLISH, the Julian Ephemeris Date.
//
{
  char cmp;
  int d1;
  int d2;
  double f1;
  double f2;
  double jed;
  int m1;
  int m2;
  int y1;
  int y2;
//
//  Check the date.
//
  ymd_check_english ( y, m, d );
//
//  Use the Julian Calendar for dates strictly before 1752/9/3.
//
  y1 = 1752;
  m1 = 9;
  d1 = 3;
  f1 = 0.0;

  cmp = ymdf_compare ( y, m, d, f, y1, m1, d1, f1 );

  if ( cmp == '<' )
  {
    jed = ymdf_to_jed_julian ( y, m, d, f );
    return jed;
  }
//
//  Use the Gregorian calendar for dates strictly after 1752/9/13.
//
  y2 = 1752;
  m2 = 9;
  d2 = 13;
  f2 = 0.0;

  cmp = ymdf_compare ( y, m, d, f, y2, m2, d2, f2 );

  if ( cmp == '>' )
  {
    jed = ymdf_to_jed_gregorian ( y, m, d, f );
    return jed;
  }
//
//  Error return if the date falls between the transition dates.
//
  cerr << "\n";
  cerr << "YMDF_TO_JED_ENGLISH - Fatal error!\n";
  cerr << "  Date is illegal.\n";
  exit ( 1 );
}
//****************************************************************************80

double ymdf_to_jed_ethiopian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ETHIOPIAN converts an Ethiopian YMDF date to a JED.
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
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ETHIOPIAN, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 4720 - ( 13 - m ) / 13;
  m_prime = ( m + 12 ) % 13;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = 30 * m_prime;

  jed = ( double ) ( j1 + j2 + d_prime - 124 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_gregorian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_GREGORIAN converts a Gregorian YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, int M, int D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_GREGORIAN, the corresponding JED.
//
{
  int d_prime;
  int g;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Account for the missing year 0 by moving negative years up one.
//
  y = y_common_to_astronomical ( y );
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 4716 - ( 14 - m ) / 12;
  m_prime = ( ( m + 9 ) % 12 );
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = ( 153 * m_prime + 2 ) / 5;

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 38;

  jed = ( double ) ( j1 + j2 + d_prime - 1401 - g ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_hebrew ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_HEBREW converts a Hebrew YMDF date to a JED.
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
//  Reference:
//
//    Edward Richards,
//    Algorithm J,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, page 334.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_HEBREW, the corresponding JED.
//
{
  double jed;
  int m2;
//
//  Check the date.
//
  ymd_check_hebrew ( y, m, d );
//
//  Determine the JED of the beginning of the year.
//
  jed = new_year_to_jed_hebrew ( y );
//
//  Work through the preceding months.
//
  for ( m2 = 1; m2 < m; m2++ )
  {
    jed = jed + ( double ) ( month_length_hebrew ( y, m2 ) );
  }
//
//  Add on the days.
//
  jed = jed + ( double ) ( d - 1 );
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_hindu_solar ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_HINDU_SOLAR converts a Hindu solar YMDF date to a JED.
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
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz, Stewart Clamen,
//    Calendrical Calculations, II: Three Historical Calendars,
//    Software - Practice and Experience,
//    Volume 23, Number 4, pages 383-404, April 1993.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_HINDU_SOLAR, the Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch;

  jed_epoch = epoch_to_jed_hindu_solar ( );

  jed = jed_epoch +
      ( double ) ( d - 1 ) 
    + ( double ) ( m - 1 ) * month_length_hindu_solar ( ) 
    + ( double ) ( y ) * year_length_hindu_solar ( );

  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_islamic_a ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ISLAMIC_A converts an Islamic A YMDF date to a JED.
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
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ISLAMIC_A, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Check the date.
//
  ymd_check_islamic ( y, m, d );
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 5519 - ( 12 - m ) / 12;
  m_prime = ( m + 11 ) % 12;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 10631 * y_prime + 14 ) / 30;

  j2 = ( 2951 * m_prime + 51 ) / 100;

  jed = ( double ) ( j1 + j2 + d_prime - 7665 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_islamic_a2 ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ISLAMIC_A2 converts an Islamic A YMDF date to a JED.
//
//  Discussion:
//
//    The algorithm has the beauty of being comprehensible//
//
//    Count the days up to the day before the start of the calendar,
//    plus the days in the current month, the 29 days guaranteed
//    in the previous months of this year, the (months/2) 30th days,
//    the 354 days in each of the previous years, plus the total number
//    of leap days in the preceding years.
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
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz,
//    Calendrical Calculations I,
//    Software - Practice and Experience,
//    Volume 20, Number 9, September 1990, pages 899-928.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ISLAMIC_A2, the corresponding Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch;
//
//  Check the date.
//
  ymd_check_islamic ( y, m, d );

  jed_epoch = epoch_to_jed_islamic_a ( );

  jed = jed_epoch + ( double ) ( - 1 + d + 29 * ( m - 1 ) + ( m / 2 ) 
    + 354 * ( y - 1 ) + ( 11 * y + 3 ) / 30 );

  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_islamic_b ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ISLAMIC_B converts an Islamic B YMDF date to a JED.
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
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ISLAMIC_B, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Check the date.
//
  ymd_check_islamic ( y, m, d );
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 5519 - ( 12 - m ) / 12;
  m_prime = ( m + 11 ) % 12;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 10631 * y_prime + 14 ) / 30;

  j2 = ( 2951 * m_prime + 51 ) / 100;

  jed = ( double ) ( j1 + j2 + d_prime - 7664 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_jelali ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_JELALI converts a Jelali YMDF date to a JED.
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
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_JELALI, the corresponding Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch;

  jed_epoch = epoch_to_jed_jelali ( );

  jed = jed_epoch + ( double ) ( ( d - 1 ) + 30 * ( m - 1 ) + 365 * ( y - 1 ) 
    + ( y - 1 ) / 4 );

  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_julian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_JULIAN converts a Julian YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, int M, int D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_JULIAN, the Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Account for the missing year 0 by moving negative years up one.
//
  y = y_common_to_astronomical ( y );
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 4716 - ( 14 - m ) / 12;
  m_prime = ( ( m + 9 ) % 12 );
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = ( 153 * m_prime + 2 ) / 5;

  jed = ( double ) ( j1 + j2 + d_prime - 1401 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_julian2 ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_JULIAN2 converts a Julian YMDF date to a JED.
//
//  Example:
//
//          Y  M  D          JED
//    --------------     -------
//    BC 4713  1  1            0
//    AD    1  1  1      1721424
//    AD 1844  5 11      2394710
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_JULIAN2, the Julian Ephemeris Date.
//
 {
  int d2;
  double jed;
  int m2;
  int y2;
  int y3;

  y2 = y;
  m2 = m;
  d2 = d;

  ymd_check_julian ( y2, m2, d2 );
//
//  Account for the missing year 0 by moving negative years up one.
//
  y3 = y_common_to_astronomical ( y2 );
//
//  The JED is the number of days in years past, plus the number of days in
//  the previous months this year, plus the number of days.
//
  jed = ( double ) ( ( ( 1461 * ( y3 + 4715 ) ) / 4 ) - 1095 
    + days_before_month_julian ( y2, m2 ) + d2 - 1 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_khwarizmian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_KHWARIZMIAN converts a Khwarizmian YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_KHWARIZMIAN, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 5348 - ( 13 - m ) / 13;
  m_prime = ( m + 12 ) % 13;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  jed = ( double ) ( 365 * y_prime + 30 * m_prime + d_prime - 317 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_macedonian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_MACEDONIAN converts a Macedonian YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_MACEDONIAN, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 4405 - ( 18 - m ) / 12;
  m_prime = ( m + 5 ) % 12;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = ( 153 * m_prime + 2 ) / 5;

  jed = ( double ) ( j1 + j2 + d_prime - 1401 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_persian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_PERSIAN converts a Persian YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_PERSIAN, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 5348 - ( 22 - m ) / 13;
  m_prime = ( m + 3 ) % 13;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  jed = ( double ) ( 365 * y_prime + 30 * m_prime + d_prime - 77 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_republican ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_REPUBLICAN converts a Republican YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_REPUBLICAN, the corresponding JED.
//
{
  int d_prime;
  int g;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Check the date.
//
  ymd_check_republican ( y, m, d );
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 6504 - ( 13 - m ) / 13;
  m_prime = ( m + 12 ) % 13;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = 30 * m_prime;

  g = ( 3 * ( ( y_prime + 396 ) / 100 ) / 4 ) - 51;
  jed = ( double ) ( j1 + j2 + d_prime - 111 - g ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_roman ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ROMAN converts a Roman YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ROMAN, the Julian Ephemeris Date.
//
{
  double jed;
  int y2;
//
//  Check the date.
//
  ymd_check_roman ( y, m, d );

  y2 = y_roman_to_julian ( y );

  jed = ymdf_to_jed_julian ( y2, m, d, f );

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_saka ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_SAKA converts a Saka YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_SAKA, the corresponding JED.
//
{
  int d_prime;
  int g;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
  int z;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 4794 - ( 13 - m ) / 12;
  m_prime = ( m + 10 ) % 12;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  z = m_prime / 6;

  j2 = ( 31 - z ) * m_prime + 5 * z;

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 36;

  jed = ( double ) ( j1 + j2 + d_prime - 1348 - g ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_soor_san ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_SOOR_SAN converts a Soor San YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_SOOR_SAN, the corresponding Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch;
 
  jed_epoch = epoch_to_jed_soor_san ( );

  jed = jed_epoch + ( double ) ( ( d - 1 ) + 30 * ( m - 1 ) + 365 * ( y - 1 ) 
    + ( y - 1 ) / 4 );
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_syrian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_SYRIAN converts a Syrian YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Richards,
//    Algorithm E,
//    Mapping Time, The Calendar and Its History,
//    Oxford, 1999, pages 323-324.
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_SYRIAN, the corresponding Julian Ephemeris Date.
//
{
  int d_prime;
  double jed;
  int j1;
  int j2;
  int m_prime;
  int y_prime;
//
//  Convert the calendar date to a computational date.
//
  y_prime = y + 4405 - ( 17 - m ) / 12;
  m_prime = ( m + 6 ) % 12;
  d_prime = d - 1;
//
//  Convert the computational date to a JED.
//
  j1 = ( 1461 * y_prime ) / 4;

  j2 = ( 153 * m_prime + 2 ) / 5;

  jed = ( double ) ( j1 + j2 + d_prime - 1401 ) - 0.5;
  jed = jed + f;

  return jed;
}
//****************************************************************************80

double ymdf_to_jed_zoroastrian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_JED_ZOROASTRIAN converts a Zoroastrian YMDF date to a JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, M, D, double F, the YMDF date.
//
//    Output, double YMDF_TO_JED_ZOROASTRIAN, the corresponding Julian Ephemeris Date.
//
{
  double jed;
  double jed_epoch;

  jed_epoch = epoch_to_jed_zoroastrian ( );

  jed = jed_epoch + ( double ) ( 
    ( d - 1 ) + 30 * ( m - 1 ) + 365 * ( y - 1 ) );
  jed = jed + f;

  return jed;
}
//****************************************************************************80

int ymdf_to_weekday_common ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_WEEKDAY_COMMON returns the weekday of a Common YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, int M, int D, double F, the YMDF date.
//
//    Output, int W, is the week day number of the date, with
//    1 for Sunday, through 7 for Saturday.
//
{
  double f2;
  double jed;
  int w;

  jed = ymdf_to_jed_common ( y, m, d, f );

  jed_to_weekday ( jed, w, f2 );

  return w;
}
//****************************************************************************80

int ymdf_to_weekday_english ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_WEEKDAY_ENGLISH returns the weekday of an English YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, int M, int D, double F, the YMDF date.
//
//    Output, int W, is the week day number of the date, with
//    1 for Sunday, through 7 for Saturday.
//
{
  double f2;
  double jed;
  int w;

  jed = ymdf_to_jed_english ( y, m, d, f );

  jed_to_weekday ( jed, w, f2 );

  return w;
}
//****************************************************************************80

int ymdf_to_weekday_gregorian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_WEEKDAY_GREGORIAN returns the weekday of a Gregorian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, int M, int D, double F, the YMDF date.
//
//    Output, int W, is the week day number of the date, with
//    1 for Sunday, through 7 for Saturday.
//
{
  double f2;
  double jed;
  int w;

  jed = ymdf_to_jed_gregorian ( y, m, d, f );

  jed_to_weekday ( jed, w, f2 );

  return w;
}
//****************************************************************************80

int ymdf_to_weekday_julian ( int y, int m, int d, double f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_WEEKDAY_JULIAN returns the weekday of a Julian YMDF date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y, int M, int D, double F, the YMDF date.
//
//    Output, int W, is the week day number of the date, with
//    1 for Sunday, through 7 for Saturday.
//
{
  double f2;
  double jed;
  int w;

  jed = ymdf_to_jed_julian ( y, m, d, f );

  jed_to_weekday ( jed, w, f2 );

  return w;
}
//****************************************************************************80

void ymdf_to_yjf_common ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_YJF_COMMON converts from YMDF to YJF form in the Common calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1, 
//    the YMDF date.
//
//    Output, int &Y2, &J2, double &F2, YJF date.
//
{
  int m;

  y2 = y1;
  j2 = d1;
  f2 = f1;
//
//  Add in the days of the elapsed months.
//
  for ( m = 1; m < m1; m++ )
  {
    j2 = j2 + month_length_common ( y2, m );
  }

  return;
}
//****************************************************************************80

void ymdf_to_yjf_english ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_YJF_ENGLISH converts from YMDF to YJF form in the English calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1,
//    the YMDF date.
//
//    Output, int &Y2, &J2, double &F2, the YJF date.
//
{
  int m;

  y2 = y1;
  j2 = d1;
  f2 = f1;
//
//  Add in the days of the elapsed months.
//
  for ( m = 1; m < m1; m++ )
  {
    j2 = j2 + month_length_english ( y2, m );
  }

  return;
}
//****************************************************************************80

void ymdf_to_yjf_gregorian ( int y1, int m1, int d1, double f1, int &y2, 
  int &j2, double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_YJF_GREGORIAN: YMDF to YJF form in the Gregorian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1, 
//    the YMDF date.
//
//    Output, int &Y2, &J2, double &F2, the YJF date.
//
{
  int m;

  y2 = y1;
  j2 = d1;
  f2 = f1;
//
//  Add in the days of the elapsed months.
//
  for ( m = 1; m < m1; m++ )
  {
    j2 = j2 + month_length_gregorian ( y2, m );
  }

  return;
}
//****************************************************************************80

void ymdf_to_yjf_hebrew ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_YJF_HEBREW converts from YMDF to YJF form in the Hebrew calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1, 
//    the YMDF date.
//
//    Output, int &Y2, &J2, double &F2, the YJF date.
//
{
  int m;

  y2 = y1;
  j2 = d1;
  f2 = f1;

  for ( m = 1; m < m1; m++ )
  {
    j2 = j2 + month_length_hebrew ( y2, m );
  }

  return;
}
//****************************************************************************80

void ymdf_to_yjf_islamic ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_YJF_ISLAMIC converts from YMDF to YJF form in the Islamic calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1, 
//    the YMDF date.
//
//    Output, int &Y2, &J2, double &F2, the YJF date.
//
{
  int m;

  y2 = y1;
  j2 = d1;
  f2 = f1;
//
//  Add in the days of the elapsed months.
//
  for ( m = 1; m < m1; m++ )
  {
    j2 = j2 + month_length_islamic ( y2, m );
  }

  return;
}
//****************************************************************************80

void ymdf_to_yjf_julian ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_YJF_JULIAN converts from YMDF to YJF form in the Julian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1, the 
//    YMDF date.
//
//    Output, int &Y2, &J2, double &F2, the YJF date.
//
{
  int m;

  y2 = y1;
  j2 = d1;
  f2 = f1;
//
//  Add in the days of the elapsed months.
//
  for ( m = 1; m < m1; m++ )
  {
    j2 = j2 + month_length_julian ( y2, m );
  }

  return;
}
//****************************************************************************80

void ymdf_to_yjf_republican ( int y1, int m1, int d1, double f1, int &y2, 
  int &j2, double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_YJF_REPUBLICAN: YMDF to YJF form in the Republican calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1, the 
//    YMDF date.
//
//    Output, int &Y2, &J2, double &F2, the YJF date.
//
{
  int m;

  y2 = y1;
  j2 = d1;
  f2 = f1;
//
//  Add in the days of the elapsed months.
//
  for ( m = 1; m < m1; m++ )
  {
    j2 = j2 + month_length_republican ( y2, m );
  }

  return;
}
//****************************************************************************80

void ymdf_to_yjf_roman ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_TO_YJF_ROMAN converts from YMDF to YJF form in the Roman calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double F1, the YMDF date.
//
//    Output, int &Y2, &J2, double &F2, the YJF date.
//
{
  int m;

  y2 = y1;
  j2 = d1;
  f2 = f1;
//
//  Add in the days of the elapsed months.
//
  for ( m = 1; m < m1; m++ )
  {
    j2 = j2 + month_length_roman ( y2, m );
  }

  return;
}
//****************************************************************************80

void ymdf_uniform_common ( int y1, int m1, int d1, double f1, int y2, int m2, 
  int d2, double f2, int &seed, int &y, int &m, int &d, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    YMDF_UNIFORM_COMMON picks a random Common YMDF date between two given dates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int Y1, M1, D1, double &F1,
//    the first YMDF date.
//
//    Input, int Y2, M2, D2, double &F2,
//    the second YMDF date.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, int &Y, &M, &D, double &F,
//    the random YMDF date.
//
{
  double jed;
  double jed1;
  double jed2;

  jed1 = ymdf_to_jed_common ( y1, m1, d1, f1 );
  jed2 = ymdf_to_jed_common ( y2, m2, d2, f2 );

  jed = r8_uniform_ab ( jed1, jed2, seed );

  jed_to_ymdf_common ( jed, y, m, d, f );

  return;
}
//****************************************************************************80

double ymdhms_to_decimal ( int y, int m, int d, int h, int n, int s )

//****************************************************************************80
//
//  Purpose:
//
//    YMDHMS_TO_DECIMAL converts a Y/M/D/H/Mn/S date to a Decimal Y.F date.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer ( kind = 4 ) Y, M, D, H, N, S, the YMDHMS date.
//
//    Output, real ( kind = 8 ) YF, the Decimal date.
//
{
  int day_min;
  int day_max;
  int days;
  double f;
  int ierror;
  double yf;
//
//  How many days between January 1st and day D?
//
  day_min = 1;
  days = ymd_dif_common ( y, m, day_min, y, m, d );
//
//  How many days in this year total?
//
  day_max = year_length_common ( y );
//
//  The decimal part of the year is ( D + H/24 + N/24*60 + S/24*60*60 ) / DMAX.
//
  f =       ( double ) ( s )      / 60.0;
  f = ( f + ( double ) ( n )    ) / 60.0;
  f = ( f + ( double ) ( h )    ) / 24.0;
  f = ( f + ( double ) ( days ) ) / ( double ) ( day_max );

  yf = ( double ) ( y ) + f;

  return yf;
}
