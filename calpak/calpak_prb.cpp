# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <string>

using namespace std;

# include "calpak.hpp"

int main ( );
void test005 ( );
void test39 ( );
void test689 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CALPAK_PRB.
//
//  Discussion:
//
//    CALPAK_PRB tests the CALPAK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CALPAK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CALPAK library.\n";

  test005 ( );
  test39 ( );
  test689 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CALPAK_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests JED_TO_WEEKDAY.
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
{
  int i;
  double f2;
  double jed1;
  double jed2;
  string s2;
  int w2;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  JED_TO_WEEKDAY reports the day of the week\n";
  cout << "  for a Julian Ephemeris Date.\n";
  cout << "\n";
  cout << "          JED     W  Name\n";
  cout << "\n";

  i = 0;

  for ( ; ; )
  {
    i = i + 1;
    jed1 = jed_test ( i );

    if ( jed1 < 0.0 )
    {
      break;
    }

    jed2 = jed_to_next_noon ( jed1 );

    jed_to_weekday ( jed2, w2, f2 );
 
    s2 = weekday_to_name_common ( w2 );

    cout << "  " << setw(11) << ( int ) jed2
         << "  " << setw(4)  << w2
         << "  " << s2 << "\n";
  }

  return;
}
//****************************************************************************80

void test39 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST39 tests MONTH_TO_MONTH_NAME_COMMON.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  string month_name;
  int months;
  int y;

  cout << "\n";
  cout << "TEST39\n";
  cout << "  For the Common calendar,\n";
  cout << "  MONTH_TO_MONTH_NAME_COMMON names the months:\n";
  cout << "\n";

  y = 1;
  months = year_length_months_common ( y );

  for ( m = 1; m <= months; m++ )
  {
    month_name = month_to_month_name_common ( m );
    cout << "  " << setw(2) << m
         << "  " << month_name << "\n";
  }

  return;
}
//****************************************************************************80

void test689 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST689 tests YMD_TO_DECIMAL
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
{
  int d;
  int dhi = 1;
  int dlo = 1;
  double f;
  double fhi = 0.0;
  double flo = 0.0;
  int i;
  int m;
  int mhi = 1;
  int mlo = 1;
  string s;
  int seed = 123456789;
  int y;
  double yf;
  int yhi = 1970;
  int ylo = 1960;

  cout << "\n";
  cout << "TEST689\n";
  cout << "  YMD_TO_DECIMAL converts a date to a year and decimal.\n";
  cout << "\n";
  cout << "  YMD                         Y.F\n";
  cout << "\n";
 
  for ( i = 1; i <= 10; i++ )
  {
    ymdf_uniform_common ( ylo, mlo, dlo, flo, yhi, mhi, dhi, fhi, 
      seed, y, m, d, f );

    s = ymd_to_s_common ( y, m, d );

    yf = ymd_to_decimal ( y, m, d );
//
//  Apparently, precision 8 means how many digits total, including
//  those BEFORE the decimal place.  Because the default for reals
//  is not "f" but "e" or "g"...  Yet another clumsy and misleading
//  feature of IOMANIP.
//
    cout << "  " << setw(13) << s << "   "
         << "  " << setw(14) << setprecision(8) << yf << "\n";
  }

  return;
}
