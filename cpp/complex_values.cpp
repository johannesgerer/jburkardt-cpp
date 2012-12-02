# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <complex>

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
//    COMPLEX_VALUES demonstrates properties of the ANSI COMPLEX class.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "COMPLEX_VALUES\n";
  cout << "  C++ version\n";
  cout << "  Demonstrate the use of the ANSI COMPLEX class.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "COMPLEX_VALUES\n";
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
//    TEST01 looks at assignments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 April 2006
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
  complex <float> a = ( 1.0, 2.0 );
  complex <float> b ( 3.0, 4.0 );
  float c;
  float d;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Look at assignment operators.\n";
  cout << "  Note that it might seem natural to try the assignment\n";
  cout << "  statement in the form:\n";
  cout << "\n";
  cout << "    variable = ( float, float )\n";
  cout << "\n";
  cout << "  but as the examples will show, this is NOT the correct way!\n";
  cout << "\n";
  cout << "\n";
  cout << "  complex <float> a = ( 1.0, 2.0 ) in initialization.\n";
  cout << "  a = " << a << "\n";

  cout << "\n";
  cout << "  complex <float> b ( 3.0, 4.0 ) in initialization.\n";
  cout << "  b = " << b << "\n";

  cout << "\n";
  cout << "  a = ( 5.0, 6.0 ) in executable statement.\n";
  a = ( 5.0, 6.0 );
  cout << "  a = " << a << "\n";

  cout << "\n";
  cout << "  a = ( complex <float> ) ( 7.0, 8.0 ) in executable statement.\n";
  a = ( complex <float> ) ( 7.0, 8.0 );
  cout << "  a = " << a << "\n";

  cout << "\n";
  cout << "  a = complex <float> ( 9.0, 10.0 ) in executable statement.\n";
  a = complex <float> ( 9.0, 10.0 );
  cout << "  a = " << a << "\n";

  cout << "\n";
  cout << "  a = complex <float> ( 11.0 ) in executable statement.\n";
  a = complex <float> ( 11.0 );
  cout << "  a = " << a << "\n";

  c = 12.0;
  d = 13.0;

  cout << "\n";
  cout << "  Use two float variables for assignment:\n";
  cout << "  c = " << c << "\n";
  cout << "  d = " << d << "\n";
  cout << "  a = complex <float> ( c, d ) in executable statement.\n";
  a = complex <float> ( c, d );
  cout << "  a = " << a << "\n";

  cout << "\n";
  cout << "  Use one float variable for assignment:\n";
  cout << "  c = " << c << "\n";
  cout << "  a = complex <float> ( c ) in executable statement.\n";
  a = complex <float> ( c );
  cout << "  a = " << a << "\n";

  cout << "\n";
  cout << "CONCLUSION:\n";
  cout << "  To initialize a complex number in a declaration:\n";
  cout << "\n";
  cout << "    complex <float> a ( 1.0, 2.0 )\n";
  cout << "\n";
  cout << "  To assign a complex number a complex value:\n";
  cout << "\n";
  cout << "    a = complex <float> ( 3.0, 4.0 )\n";
  cout << "    a = complex <float> ( b, c )\n";
  cout << "\n";
  cout << "  To assign a complex number a real value:\n";
  cout << "\n";
  cout << "    a = complex <float> ( 5.0 )\n";
  cout << "    a = complex <float> ( d )\n";
  
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 looks at array initialization.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2006
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
  complex <float> a[3] = {
    ( 1.0, 2.0 ),
    ( 3.0, 4.0 ),
    ( 5.0, 6.0 ) };
  complex <float> b[3] = {
    complex <float> ( 1.0, 2.0 ),
    complex <float> ( 3.0, 4.0 ),
    complex <float> ( 5.0, 6.0 ) };
  int i;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Look at complex array initialization.\n";
  cout << "\n";
  cout << "    complex <float> a[3] = {\n";
  cout << "      ( 1.0, 2.0 ),\n";
  cout << "      ( 3.0, 4.0 ),\n";
  cout << "      ( 5.0, 6.0 ) }\n";
  cout << "\n";
  cout << "  A:\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << setw(12) << a[i] << "\n";
  }
  cout << "\n";
  cout << "    complex <float> b[3] = {\n";
  cout << "      complex <float> ( 1.0, 2.0 ),\n";
  cout << "      complex <float> ( 3.0, 4.0 ),\n";
  cout << "      complex <float> ( 5.0, 6.0 ) }\n";
  cout << "\n";
  cout << "  B:\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << setw(12) << b[i] << "\n";
  }
  cout << "\n";
  cout << "CONCLUSION:\n";
  cout << "  To initialize a complex vector in a declaration:\n";
  cout << "\n";
  cout << "    complex <float> b[3] = {\n";
  cout << "      complex <float> ( 1.0, 2.0 ),\n";
  cout << "      complex <float> ( 3.0, 4.0 ),\n";
  cout << "      complex <float> ( 5.0, 6.0 ) }\n";
  cout << "\n";
  
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 looks at function calls.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2007
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
  complex <float> a = complex <float> ( 3.0, 4.0 );
  float pi = 3.14159265;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Look at new function calls.\n";
  cout << "\n";
  cout << "\n";
  cout << "  A =          " << setw(12) << a << "\n";
  cout << "  real ( A ) = " << setw(12) << real ( a ) << "\n";
  cout << "  imag ( A ) = " << setw(12) << imag ( a ) << "\n";
  cout << "  arg ( A ) =  " << setw(12) << arg ( a ) << "\n";
  cout << "  norm ( A ) = " << setw(12) << norm ( a ) << "\n";
  cout << "  abs ( A ) =  " << setw(12) << abs ( a ) << "\n";
  cout << "  conj ( A ) = " << setw(12) << conj ( a ) << "\n";

  cout << "  polar ( abs ( a ), arg ( a ) ) = " << setw(12) 
       << polar ( abs ( a ), arg ( a ) ) << "\n";

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
//    24 September 2003
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
