# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void test01 ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    DRAND48_TEST generates random numbers using DRAND48().
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "DRAND48_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Generate random numbers using\n";
  cout << "  SRAND48 to set the seed, and\n";
  cout << "  DRAND48 to return the random values.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << "\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DRAND48_TEST:\n";
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
//    TEST01 simply calls DRAND48 a few times.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  long long int seed = 123456789LL;
  double x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Call SRAND48 to initialize the seed.\n";
  cout << "  Call DRAND48 to generate some values.\n";
  cout << "\n";
  cout << "  Initial SEED = " << seed << "\n";

  srand48 ( seed );
  cout << "\n";
  cout << "      Step    DRAND48()\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    x = drand48 ( );
    cout << "  " << setw(8) << i
         << "  " << setw(12) << x << "\n";
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
