# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "weekday.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WEEKDAY_PRB.
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
{
  timestamp ( );

  cout << "\n";
  cout << "WEEKDAY_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the WEEKDAY library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WEEKDAY_PRB:\n";
  cout << "  Noraml end of execution.\n";

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
//    TEST01 tests YMD_TO_WEEKDAY_COMMON.
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
  string s1;
  string s2;
  string s3;
  int w1;
  int w2;
  int y;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For dates in the Common calendar:\n";
  cout << "  YMD_TO_WEEKDAY_COMMON returns the day of the week.\n";
  cout << "\n";
  cout << "  YMD                   Weekday    Weekday\n";
  cout << "                        Tabulated  Computed\n";
  cout << "\n";

  for ( ; ; )
  {
    weekday_values ( n_data, y, m, d, w1 );

    if ( n_data == 0 )
    {
      break;
    }

    s3 = ymd_to_s_common ( y, m, d );
    w2 = ymd_to_weekday_common ( y, m, d );
    s1 = weekday_to_name_common ( w1 );
    s2 = weekday_to_name_common ( w2 );

    cout << "  " << setw(20) << s3
         << "  " << setw(9)  << s1
         << "  " << setw(9)  << s2 << "\n";
  }
  return;
}
