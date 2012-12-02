# include <cstdlib>
# include <iostream>
# include <iomanip>
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
//    Here is why it's good to initialize array pointers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2006
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int *b = NULL;

  timestamp ( );

  cout << "\n";
  cout << "NOT_ALLOCATED_ARRAYS\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  When an array starts out as a pointer, you have to use the\n";
  cout << "  NEW command to allocate memory.  You should always initialize\n";
  cout << "  such array pointers to NULL, so you can tell if they've been\n";
  cout << "  allocated or not!\n";
  cout << "\n";
  cout << "  Unfortunately, when you DELETE an array, you also have to\n";
  cout << "  reset the pointer to NULL; that does not happen automatically\n";
  cout << "  either!\n";

  cout << "\n";
  cout << "  The pointer A is not preset to NULL.\n";
  cout << "  Before allocation, we check the value:\n";
  cout << "    a = " << a << "\n";
  cout << "  The test 'if ( !a )' is not guaranteed to return 1\n";
  cout << "  because we did not initialize A properly.\n";
  cout << "    !a = " << !a << "\n";
  cout << "  Now we allocate A.\n" << flush;

  a = new int[10];

  cout << "    a = " << a << "\n";
  cout << "    !a = " << !a << "\n";
  cout << "  Now we DELETE A.\n";
  delete [] a;
  cout << "    a = " << a << "\n";
  cout << "    !a = " << !a << "\n";
  cout << "  Now we RESET A to NULL!\n";
  a = NULL;
  cout << "    a = " << a << "\n";
  cout << "    !a = " << !a << "\n";

  cout << "\n";
  cout << "  The pointer B is preset to NULL.\n";
  cout << "  Before allocation, we check the value:\n";
  cout << "    b = " << b << "\n";
  cout << "  The test 'if ( !b )' is guaranteed to return 1\n";
  cout << "  because we initialized B properly.\n";
  cout << "    !b = " << !b << "\n";
  cout << "  Now we allocate B.\n" << flush;

  b = new int[10];

  cout << "    b = " << b << "\n";
  cout << "    !b = " << !b << "\n";
  cout << "  Now we DELETE B.\n";
  delete [] b;
  cout << "    b = " << b << "\n";
  cout << "    !b = " << !b << "\n";
  cout << "  Now we RESET B to NULL!\n";
  b = NULL;
  cout << "    b = " << b << "\n";
  cout << "    !b = " << !b << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "NOT_ALLOCATED_ARRAYS:\n";
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
