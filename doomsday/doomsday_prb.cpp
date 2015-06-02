# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "doomsday.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DOOMSDAY_PRB.
//
//  Discussion:
//
//    DOOMSDAY_PRB tests the DOOMSDAY library.
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
{
  timestamp ( );
  cout << "\n";
  cout << "DOOMSDAY_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the DOOMSDAY library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DOOMSDAY_PRB:\n";
  cout << "  Test the DOOMSDAY library.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests DOOMSDAY against a couple of test dates.
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
{
  int d;
  int m;
  int n_data;
  int w;
  string s1;
  string s2;
  int y;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Try a couple selected dates.\n";
  cout << "\n";
  cout << "  YYYY  MM  DD  Weekday    Weekday\n";
  cout << "                Tabulated  Computed\n";
  cout << "\n";

  y = 1989;
  m = 7;
  d = 13;
  w = doomsday_gregorian ( y, m, d );
  s1 = weekday_to_name_common ( w );
  s2 = "Thursday";
  cout << "  " << setw(4) << y
       << "  " << setw(2) << m
       << "  " << setw(2) << d
       << "  " << setw(10) << s1
       << "  " << setw(10) << s2 << "\n";

  y = 2012;
  m = 5;
  d = 26;
  w = doomsday_gregorian ( y, m, d );
  s1 = weekday_to_name_common ( w );
  s2 = "Saturday";
  cout << "  " << setw(4) << y
       << "  " << setw(2) << m
       << "  " << setw(2) << d
       << "  " << setw(10) << s1
       << "  " << setw(10) << s2 << "\n";


  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests DOOMSDAY against a number of known values.
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
{
  int d;
  int m;
  int n_data;
  int w1;
  int w2;
  string s1;
  string s2;
  int y;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  WEEKDAY_VALUES supplies a list of dates and weekdays.\n";
  cout << "\n";
  cout << "  YYYY  MM  DD  Weekday    Weekday\n";
  cout << "                Tabulated  Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    weekday_values ( n_data, y, m, d, w1 );

    if ( n_data <= 0 )
    {
      break;
    }
//
//  The transition from Julian to Gregorian calendars occurred in 1582
//  (for some people).  The data in "WEEKDAY_VALUES" before the transition
//  is stored in Julian format, which DOOMSDAY_GREGORIAN can't handle.
//  So let's just refuse to handle 1582 or earlier//
//
    if ( y <= 1582 )
    {
      continue;
    }

    w2 = doomsday_gregorian ( y, m, d );

    s1 = weekday_to_name_common ( w1 );
    s2 = weekday_to_name_common ( w2 );

    cout << "  " << setw(4) << y
         << "  " << setw(2) << m
         << "  " << setw(2) << d
         << "  " << setw(10) << s1
         << "  " << setw(10) << s2 << "\n";
  }

  return;
}
