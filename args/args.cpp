# include <cstdlib>
# include <iostream>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void timestamp ( void );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ARGS.
//
//  Discussion:
//
//    ARGS reports on the command line arguments.
//
//    On input, ARGC is the number of command line arguments (including the
//    program name itself), and ARGV is a pointer to an array of pointers, 
//    the null-terminated strings that constitute the command line arguments.
//
//    Thus, if we have invoked a program by:
//
//      fred 1 alpha 3.7
//
//    then on input
//
//      ARGC is 4, and 
//
//      ARGV --> ARGV(0) --> "fred"
//               ARGV(1) --> "1"
//               ARGV(2) --> "alpha"
//               ARGV(3) --> "3.7"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Deitel, Deitel,
//    C++: How to Program,
//    Third Edition,
//    Prentice Hall, 2001.
//
{
  int i;
  bool VERBOSE = true;

  if ( VERBOSE )
  {
    timestamp ( );

    cout << "\n";
    cout << "ARGS\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Print the command line arguments of a C++ program.\n";
    cout << "\n";
    cout << "  ARGC reports the number of arguments as " << argc << ".\n";
    cout << "\n";
  }

  for ( i = 0; i < argc; i++ ) 
  {
    cout << i << "  " << *argv << "\n";
    argv++;
  }

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "ARGS:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

void timestamp ( void )

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
//    02 October 2003
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
