# include <cstdlib>
# include <iostream>
# include <ctime>
# include <limits>

using namespace std;

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    LIMITS returns some information about numeric datatypes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "LIMITS\n";
  cout << "  C++ version\n";
  cout << "  Return the limits on various data types.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << "\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LIMITS\n";
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
//    TEST01 returns information about the FLOAT data type.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Information about the FLOAT data type:\n";
  cout << "\n";

  cout << "  min =            " << numeric_limits<float>::min() << "\n";
  cout << "  max =            " << numeric_limits<float>::max() << "\n";
  cout << "  epsilon =        " << numeric_limits<float>::epsilon() << "\n";
  cout << "  max_exponent10 = " 
    << numeric_limits<float>::max_exponent10 << "\n";
  cout << "  min_exponent10 = " 
    << numeric_limits<float>::min_exponent10 << "\n";
  cout << "  max_exponent =   " 
    << numeric_limits<float>::max_exponent << "\n";
  cout << "  min_exponent =   " 
    << numeric_limits<float>::min_exponent << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 returns information about the DOUBLE data type.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Information about the DOUBLE data type:\n";
  cout << "\n";

  cout << "  min =            " << numeric_limits<double>::min() << "\n";
  cout << "  max =            " << numeric_limits<double>::max() << "\n";
  cout << "  epsilon =        " << numeric_limits<double>::epsilon() << "\n";
  cout << "  max_exponent10 = " 
    << numeric_limits<double>::max_exponent10 << "\n";
  cout << "  min_exponent10 = " 
    << numeric_limits<double>::min_exponent10 << "\n";
  cout << "  max_exponent =   " 
    << numeric_limits<double>::max_exponent << "\n";
  cout << "  min_exponent =   " 
    << numeric_limits<double>::min_exponent << "\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 returns information about the LONG DOUBLE data type.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Information about the LONG DOUBLE data type:\n";
  cout << "\n";

  cout << "  min =            " << numeric_limits<long double>::min() << "\n";
  cout << "  max =            " << numeric_limits<long double>::max() << "\n";
  cout << "  epsilon =        " 
    << numeric_limits<long double>::epsilon() << "\n";
  cout << "  max_exponent10 = " 
    << numeric_limits<long double>::max_exponent10 << "\n";
  cout << "  min_exponent10 = " 
    << numeric_limits<long double>::min_exponent10 << "\n";
  cout << "  max_exponent =   " 
    << numeric_limits<long double>::max_exponent << "\n";
  cout << "  min_exponent =   " 
    << numeric_limits<long double>::min_exponent << "\n";

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 returns information about the SHORT data type.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST04:\n";
  cout << "  Information about the SHORT data type:\n";
  cout << "\n";

  cout << "  min =            " << numeric_limits<short>::min() << "\n";
  cout << "  max =            " << numeric_limits<short>::max() << "\n";

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 returns information about the INT data type.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST05:\n";
  cout << "  Information about the INT data type:\n";
  cout << "\n";

  cout << "  min =            " << numeric_limits<int>::min() << "\n";
  cout << "  max =            " << numeric_limits<int>::max() << "\n";

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 returns information about the LONG data type.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST06:\n";
  cout << "  Information about the LONG data type:\n";
  cout << "\n";

  cout << "  min =            " << numeric_limits<long>::min() << "\n";
  cout << "  max =            " << numeric_limits<long>::max() << "\n";

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
