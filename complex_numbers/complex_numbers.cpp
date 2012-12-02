# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <ctime>

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMPLEX_NUMBERS is a program which demonstrates the use of complex numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "COMPLEX_NUMBERS\n";
  cout << "  C++ version\n";
  cout << "  Demonstrate complex number usage.\n";
//
//  Single precision complex numbers: "complex <float>".
//
  test01 ( );
  test02 ( );
  test03 ( );
//
//  Double precision complex numbers: "complex <double>".
//
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "COMPLEX_NUMBERS\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 demonstrate declaration and assignment for complex <float> variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
{
//
//  Declare a complex <float> number A.
//  Declare a complex <float> vector B.
//  Declare a complex <float> array C.
//
  complex <float> a;
  complex <float> b[3];
  complex <float> c[2][2];
  int i;
  int j;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Declare a COMPLEX <FLOAT> variable.\n";
  cout << "  Assign value with an = statement.\n";
//
//  Assign values to A, B, and C.
//
  a = complex <float> ( 1.0, 2.0 );
  b[0] = complex <float> ( 1.0, 2.0 );
  b[1] = complex <float> ( 3.0, 4.0 );
  b[2] = complex <float> ( 5.0, 6.0 );

  c[0][0] = complex <float> ( 1.0, 0.1 );
  c[0][1] = complex <float> ( 1.0, 0.2 );
  c[1][0] = complex <float> ( 2.0, 0.1 );
  c[1][1] = complex <float> ( 2.0, 0.2 );
//
//  Print them.
//
  cout << "\n";
  cout << "  Scalar A:\n";
  cout << "\n";

  cout << "  " << a << "\n";

  cout << "\n";
  cout << "  Vector B:\n";
  cout << "\n";

  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << b[i] << "\n";
  }

  cout << "\n";
  cout << "  Array C:\n";
  cout << "\n";

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      cout << "  " << c[i][j];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02: declaration and iniitialization for complex <float> variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
{
//
//  Declare and initialize a complex <float> number A.
//  Declare and initialize a complex <float> vector B.
//  Declare and initialize a complex <float> array C.
//
  complex <float> a = complex <float> ( 1.0, 2.0 );
  complex <float> b[3] = {
    complex <float> ( 1.0, 2.0 ), 
    complex <float> ( 3.0, 4.0 ), 
    complex <float> ( 5.0, 6.0 ) };
  complex <float> c[2][2] = {
    { complex <float> ( 1.0, 0.1 ), complex <float> ( 1.0, 0.2 ) },
    { complex <float> ( 2.0, 0.1 ), complex <float> ( 2.0, 0.2 ) } };
  int i;
  int j;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Declare a COMPLEX <FLOAT> variable.\n";
  cout << "  Initialize as part of the declaration.\n";
//
//  Print them.
//
  cout << "\n";
  cout << "  Scalar A:\n";
  cout << "\n";

  cout << "  " << a << "\n";

  cout << "\n";
  cout << "  Vector B:\n";
  cout << "\n";

  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << b[i] << "\n";
  }

  cout << "\n";
  cout << "  Array C:\n";
  cout << "\n";

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      cout << "  " << c[i][j];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03: intrinsic functions for complex <float> variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> a = complex <float> ( 1.0, 2.0 );
  float pi = 3.141592653589793;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Apply operations and intrinsic functions \n";
  cout << "  to COMPLEX <FLOAT> variable.\n";

  cout << "\n";
  cout << "  a =                        " << a << "\n";
  cout << "  - a =                      " << - a << "\n";
//cout << "  a + 3 =                    " << a + 3 << "\n";
  cout << "  a + complex <float>(0,5) = " << a + complex <float>(0,5) << "\n";
//cout << "  4 * a =                    " << 4 * a << "\n";
//cout << "  a / 8 =                    " << a / 8 << "\n";
  cout << "  pow ( a, 2 ) =             " << pow ( a, 2 ) << "\n";
//cout << "  pow ( 2, a ) =             " << pow ( 2, a ) << "\n";
//cout << "  pow ( 2.1, a ) =             " << pow ( 2, a ) << "\n";
  cout << "  pow ( a, a ) =             " << pow ( a, a ) << "\n";
//cout << "  1 / a =                    " << 1 / a << "\n";
  cout << "\n";
  cout << "  abs ( a ) =                " << abs ( a ) << "\n";
  cout << "  arg ( a ) =                " << arg ( a ) << "\n";
  cout << "  conj ( a ) =               " << conj ( a ) << "\n";
  cout << "  cos ( a ) =                " << cos ( a ) << "\n";
  cout << "  cosh ( a ) =               " << cosh ( a ) << "\n";
  cout << "  exp ( a ) =                " << exp ( a ) << "\n";
  cout << "  imag ( a ) =               " << imag ( a ) << "\n";
  cout << "  log ( a ) =                " << log ( a ) << "\n";
  cout << "  log10 ( a ) =              " << log10 ( a ) << "\n";
  cout << "  norm ( a ) =               " << norm ( a ) << "\n";
  cout << "  polar (10.0,pi/4) =        " << polar ( 10.0, pi/4.0 ) << "\n";
  cout << "  real ( a ) =               " << real ( a ) << "\n";
  cout << "  sin ( a ) =                " << sin ( a ) << "\n";
  cout << "  sinh ( a ) =               " << sinh ( a ) << "\n";
  cout << "  sqrt ( a ) =               " << sqrt ( a ) << "\n";
  cout << "  tan ( a ) =                " << tan ( a ) << "\n";
  cout << "  tanh ( a ) =               " << tanh ( a ) << "\n";

  return;
}
//****************************************************************************80

void test04 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 demonstrate declaration and assignment for complex <double> variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
{
//
//  Declare a complex <double> number A.
//  Declare a complex <double> vector B.
//  Declare a complex <double> array C.
//
  complex <double> a;
  complex <double> b[3];
  complex <double> c[2][2];
  int i;
  int j;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Declare a COMPLEX <DOUBLE> variable.\n";
  cout << "  Assign value with an = statement.\n";
//
//  Assign values to A, B, and C.
//
  a = complex <double> ( 1.0, 2.0 );

  b[0] = complex <double> ( 1.0, 2.0 );
  b[1] = complex <double> ( 3.0, 4.0 );
  b[2] = complex <double> ( 5.0, 6.0 );

  c[0][0] = complex <double> ( 1.0, 0.1 );
  c[0][1] = complex <double> ( 1.0, 0.2 );
  c[1][0] = complex <double> ( 2.0, 0.1 );
  c[1][1] = complex <double> ( 2.0, 0.2 );
//
//  Print them.
//
  cout << "\n";
  cout << "  Scalar A:\n";
  cout << "\n";

  cout << "  " << a << "\n";

  cout << "\n";
  cout << "  Vector B:\n";
  cout << "\n";

  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << b[i] << "\n";
  }

  cout << "\n";
  cout << "  Array C:\n";
  cout << "\n";

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      cout << "  " << c[i][j];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test05 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05: declaration and iniitialization for complex <double> variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
{
//
//  Declare and initialize a complex <double> number A.
//  Declare and initialize a complex <double> vector B.
//  Declare and initialize a complex <double> array C.
//
  complex <double> a = complex <double> ( 1.0, 2.0 );
  complex <double> b[3] = {
    complex <double> ( 1.0, 2.0 ), 
    complex <double> ( 3.0, 4.0 ), 
    complex <double> ( 5.0, 6.0 ) };
  complex <double> c[2][2] = {
    { complex <double> ( 1.0, 0.1 ), complex <double> ( 1.0, 0.2 ) },
    { complex <double> ( 2.0, 0.1 ), complex <double> ( 2.0, 0.2 ) } };
  int i;
  int j;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Declare a COMPLEX <DOUBLE> variable.\n";
  cout << "  Initialize as part of the declaration.\n";
//
//  Print them.
//
  cout << "\n";
  cout << "  Scalar A:\n";
  cout << "\n";

  cout << "  " << a << "\n";

  cout << "\n";
  cout << "  Vector B:\n";
  cout << "\n";

  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << b[i] << "\n";
  }

  cout << "\n";
  cout << "  Array C:\n";
  cout << "\n";

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      cout << "  " << c[i][j];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06: intrinsic functions for complex <double> variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> a = complex <double> ( 1.0, 2.0 );
  double pi = 3.141592653589793;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Apply operations and intrinsic functions \n";
  cout << "  to COMPLEX <DOUBLE> variable.\n";

  cout << "\n";
  cout << "  a =                         " << a << "\n";
  cout << "  - a =                       " << - a << "\n";
//cout << "  a + 3 =                     " << a + 3 << "\n";
  cout << "  a + complex <double>(0,5) = " << a + complex <double>(0,5) << "\n";
//cout << "  4 * a =                     " << 4 * a << "\n";
//cout << "  a / 8 =                     " << a / 8 << "\n";
  cout << "  pow ( a, 2 ) =              " << pow ( a, 2 ) << "\n";
//cout << "  pow ( 2, a ) =              " << pow ( 2, a ) << "\n";
//cout << "  pow ( 2.1, a ) =            " << pow ( 2, a ) << "\n";
  cout << "  pow ( a, a ) =              " << pow ( a, a ) << "\n";
//cout << "  1 / a =                     " << 1 / a << "\n";
  cout << "\n";
  cout << "  abs ( a ) =                 " << abs ( a ) << "\n";
  cout << "  arg ( a ) =                 " << arg ( a ) << "\n";
  cout << "  conj ( a ) =                " << conj ( a ) << "\n";
  cout << "  cos ( a ) =                 " << cos ( a ) << "\n";
  cout << "  cosh ( a ) =                " << cosh ( a ) << "\n";
  cout << "  exp ( a ) =                 " << exp ( a ) << "\n";
  cout << "  imag ( a ) =                " << imag ( a ) << "\n";
  cout << "  log ( a ) =                 " << log ( a ) << "\n";
  cout << "  log10 ( a ) =               " << log10 ( a ) << "\n";
  cout << "  norm ( a ) =                " << norm ( a ) << "\n";
  cout << "  polar (10.0,pi/4) =         " << polar ( 10.0, pi/4.0 ) << "\n";
  cout << "  real ( a ) =                " << real ( a ) << "\n";
  cout << "  sin ( a ) =                 " << sin ( a ) << "\n";
  cout << "  sinh ( a ) =                " << sinh ( a ) << "\n";
  cout << "  sqrt ( a ) =                " << sqrt ( a ) << "\n";
  cout << "  tan ( a ) =                 " << tan ( a ) << "\n";
  cout << "  tanh ( a ) =                " << tanh ( a ) << "\n";

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
