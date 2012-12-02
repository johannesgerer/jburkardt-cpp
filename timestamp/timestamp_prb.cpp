# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "timestamp.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP_PRB demonstrates the use of TIMESTAMP.
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
{
  timestamp ( );

  cout << "\n";
  cout << "TIMESTAMP_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TIMESTAMP library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TIMESTAMP_PRB\n";
  cout << "  Normal end of execution.\n";
 
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
//    TEST01 demonstrates the use of TIMESTAMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST01\n";
  cout << "  TIMESTAMP prints out the current wallclock time,\n";
  cout << "  including the year, month, day, hours, minutes,\n";
  cout << "  seconds, thousandths of a second, and AM/PM.\n";
  cout << "\n";
  cout << "  This can be useful in keeping track of the date\n";
  cout << "  of execution of a particular program\n";
  cout << "  or to give a rough idea of the length of time\n";
  cout << "  required to run a program.\n";

  cout << "\n";
  timestamp ( );

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates the use of TIMESTRING.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  char *s;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  TIMESTRING returns the current wallclock time,\n";
  cout << "  including the year, month, day, hours, minutes,\n";
  cout << "  seconds, thousandths of a second, and AM/PM\n";
  cout << "  in a string, which the user may print or manipulate.\n";

  s = timestring ( );

  cout << "\n";
  cout << "  TIMESTRING returned the value \"" << s << "\"\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 demonstrates the use of TIME_NUMBERS
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
{
  int *time_vec;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  TIME_NUMBERS returns the date as a string of integers.\n";

  time_vec = time_numbers ( );

  cout << "\n";
  cout << "  Year =        " << time_vec[0] << "\n";
  cout << "  Month =       " << time_vec[1] << "\n";
  cout << "  Day =         " << time_vec[2] << "\n";
  cout << "  Hour =        " << time_vec[3] << "\n";
  cout << "  Minute =      " << time_vec[4] << "\n";
  cout << "  Second =      " << time_vec[5] << "\n";

  delete [] time_vec;

  return;
}
