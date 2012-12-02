# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( );
int search ( int a, int b, int c );
int f ( int i );
void timestamp ( );
double cpu_time ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SEARCH_SERIAL.
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
  int a;
  int b;
  int c;
  int fj;
  int i4_huge = 2147483647;
  int j;
  double wtime;

  a = 1;
  b = i4_huge;
  c = 45;

  timestamp ( );
  cout << "\n";
  cout << "SEARCH_SERIAL:\n";
  cout << "  C++ version\n";
  cout << "  Search the integers from A to B\n";
  cout << "  for a value J such that F(J) = C.\n";
  cout << "\n";
  cout << "  A           = " << a << "\n";
  cout << "  B           = " << b << "\n";
  cout << "  C           = " << c << "\n";

  wtime = cpu_time ( );

  j = search ( a, b, c );

  wtime = cpu_time ( ) - wtime;

  if ( j == -1 )
  {
    cout << "\n";
    cout << "  No solution was found.\n";
  }
  else
  {
    cout << "\n";
    cout << "  Found     J = " << j << "\n";
    cout << "  Verify F(J) = " << f ( j ) << "\n";
  }

  cout << "  Elapsed CPU time is " << wtime << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "SEARCH_SERIAL:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int search ( int a, int b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    SEARCH searches integers in [A,B] for a J so that F(J) = C.
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
//  Parameters:
//
//    Input, int A, B, the search range.
//
//    Input, int C, the desired function value.
//
//    Output, int SEARCH, the computed solution, or -1
//    if no solution was found.
//
{
  int fi;
  int i;
  int j;

  j = -1;

  for ( i = a; i <= b; i++ )
  {
    fi = f ( i );

    if ( fi == c )
    {
      j = i;
      break;
    }
  }

  return j;
}
//****************************************************************************80

int f ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    F is the function we are analyzing.
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
//  Parameters:
//
//    Input, int I, the argument.
//
//    Input, int F, the value.
//
{
  int i4_huge = 2147483647;
  int j;
  int k;
  int value;

  value = i;

  for ( j = 1; j <= 5; j++ )
  {
    k = value / 127773;

    value = 16807 * ( value - k * 127773 ) - k * 2836;

    if ( value <= 0 )
    {
      value = value + i4_huge;
    }
  }

  return value;
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
//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Discussion:
//
//    The data available to this routine through "CLOCK" is not very reliable,
//    and hence the values of CPU_TIME returned should not be taken too 
//    seriously, especially when short intervals are being timed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

  return value;
}
