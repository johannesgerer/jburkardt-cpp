# include <cstdlib>
# include <iostream>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIZES returns the size of various data types.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "SIZES\n";
  cout << "  C++ version\n";
  cout << "  Return the size of various data types.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << "\n";

  cout << "\n";
  cout << "         char " << sizeof ( char )        << " bytes.\n";
  cout << "\n";
  cout << "    short int " << sizeof ( short int )    << " bytes.\n";
  cout << "          int " << sizeof ( int )         << " bytes.\n";
  cout << "     long int " << sizeof ( long int )     << " bytes.\n";
  cout << "long long int " << sizeof ( long long int ) << " bytes.\n";
  cout << "\n";
  cout << "        float " << sizeof ( float )        << " bytes.\n";
  cout << "       double " << sizeof ( double )       << " bytes.\n";
  cout << "  long double " << sizeof ( long double )   << " bytes.\n";
  cout << "\n";
  cout << "         bool " << sizeof ( bool )         << " bytes.\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "SIZES\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
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
