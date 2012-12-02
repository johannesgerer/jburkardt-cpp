# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );
void test03 ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for G++_INTRINSICS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "G++_QUADMATH:\n";
  cout << "  C++ version\n";
  cout << "  Test the G++ quadmath library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
///
//  Terminate.
//
  cout << "\n";
  cout << "GCC_QUADMATH:\n";
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
//    TEST01 uses flaot arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int divs;
  float x;
  float x_old;
  float y;
  float z;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Using FLOAT arithmetic:\n";
  cout << "  Compute smallest 1/2^DIV that can be added to 1.\n";

  x = 1.0;
  z = 1.0;
  divs = 0;

  for ( ; ; )
  {
    x_old = x;
    x = x / 2.0;
    y = 1.0 + x;
    if ( y <= z )
    {
      break;
    }
    divs = divs + 1;
  }

  cout << "  Number of divisions DIV = " << divs << "\n";
  cout << "  1/2^DIV =         " << x_old << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses double arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int divs;
  double x;
  double x_old;
  double y;
  double z;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Using DOUBLE arithmetic:\n";
  cout << "  Compute smallest 1/2^DIV that can be added to 1.\n";

  x = 1.0;
  z = 1.0;
  divs = 0;

  for ( ; ; )
  {
    x_old = x;
    x = x / 2.0;
    y = 1.0 + x;
    if ( y <= z )
    {
      break;
    }
    divs = divs + 1;
  }

  cout << "  Number of divisions DIV = " << divs << "\n";
  cout << "  1/2^DIV =         " << x_old << "\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses __float80 arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int divs;
  __float80 x;
  __float80 x_old;
  __float80 y;
  __float80 z;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Using __FLOAT80 arithmetic:\n";
  cout << "  Compute smallest 1/2^DIV that can be added to 1.\n";

  x = 1.0;
  z = 1.0;
  divs = 0;

  for ( ; ; )
  {
    x_old = x;
    x = x / 2.0;
    y = 1.0 + x;
    if ( y <= z )
    {
      break;
    }
    divs = divs + 1;
  }

  cout << "  Number of divisions DIV = " << divs << "\n";
  cout << "  1/2^DIV =         " << x_old << "\n";

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
