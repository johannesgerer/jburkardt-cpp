# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <ctime>
# include <float.h>
# include <limits.h>

using namespace std;

int main ( );
void test01 ( );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BIG_INTS_REAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BIG_INTS_REAL:\n";
  cout << "  C++ version\n";
  cout << "  Examine the transfer of integer values into and out of real variables.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BIG_INTS_REAL:\n";
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
//    TEST01 stores huge integers as reals.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i4;
  int i4r4i4;
  int i4r8i4;
  long long int i8;
  long long int i8r4i8;
  long long int i8r8i8;
  float r4;
  float r4i4;
  float r4i8;
  double r8;
  double r8i4;
  double r8i8;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Compute the largest possible integers.\n";
  cout << "  Try to store them as real values.\n";
  cout << "  Then copy them back.\n";

  cout << "\n";
  cout << "  'Huge' integers and huge reals:\n";
  cout << "\n";

  i4 = INT_MAX;
  i8 = LLONG_MAX;
  r4 = FLT_MAX;
  r8 = DBL_MAX;

  cout << "  i4 = INT_MAX   = " << i4 << "\n";
  cout << "  i8 = LLONG_MAX = " << i8 << "\n";
  cout << "  r4 = FLT_MAX   = " << r4 << "\n";
  cout << "  r8 = DBL_MAX   = " << r8 << "\n";

  cout << "\n";
  cout << "  Convert huge integers to real values:\n";
  cout << "\n";

  r4i4 = ( float ) ( i4 );
  r4i8 = ( float ) ( i8 );
  r8i4 = ( double ) ( i4 );
  r8i8 = ( double ) ( i8 );

  cout << "  r4i4 = ( float ) ( i4 )  = " << r4i4 << "\n";
  cout << "  r4i8 = ( float ) ( i8 )  = " << r4i8 << "\n";
  cout << "  r8i4 = ( double ) ( i4 ) = " << r8i4 << "\n";
  cout << "  r8i8 = ( double ) ( i8 ) = " << r8i8 << "\n";

  cout << "\n";
  cout << "  Convert real values of integers back to integers:\n";
  cout << "\n";

  i4r4i4 = ( int ) ( r4i4 );
  i4r8i4 = ( int ) ( r8i4 );
  i8r4i8 = ( long long int ) ( r4i8 );
  i8r8i8 = ( long long int ) ( r8i8 );

  cout << "  i4r4i4 = ( int ) ( r4i4 )           = " << i4r4i4 << "\n";
  cout << "  i4r8i4 = ( int ) ( r8i4 )           = " << i4r8i4 << "\n";
  cout << "  i8r4i8 = ( long long int ) ( r4i8 ) = " << i8r4i8 << "\n";
  cout << "  i8r8i8 = ( long long int ) ( r8i8 ) = " << i8r8i8 << "\n";

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
