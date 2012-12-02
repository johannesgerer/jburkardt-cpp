# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "blas1_c.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
complex <float> c4_uniform_01 ( int *seed );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLAS1_C_PRB.
//
//  Discussion:
//
//    BLAS1_C_PRB tests the BLAS1 single precision complex routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "BLAS1_C_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the BLAS1_C library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BLAS1_C_PRB:\n";
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
//    TEST01 tests CABS1.
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
{
  complex <float> c;
  float c_norm;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  CABS1 returns the L1 norm of a complex number.\n";
  cout << "\n";
  cout << "      Real      Imaginary\n";
  cout << "      Part      Part           CABS1(Z)\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
//
//  The compiler seems to be unhappy with the statement
//
//    c = 5.0 * c4_uniform_01 ( &seed );
//
//  Poor compiler.
//
    c = c4_uniform_01 ( &seed );
    c = complex <float> ( 5.0 ) * c;

    c_norm = cabs1 ( c );

    cout << "  " << setw(10) << real ( c )
         << "  " << setw(10) << imag ( c )
         << "  " << setw(10) << c_norm << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CABS2.
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
{
  complex <float> c;
  float c_norm;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  CABS2 returns the L2 norm of a complex number.\n";
  cout << "\n";
  cout << "      Real      Imaginary\n";
  cout << "      Part      Part           CABS2(Z)\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    c = c4_uniform_01 ( &seed );
    c = complex <float> ( 5.0 ) * c;

    c_norm = cabs2 ( c );

    cout << "  " << setw(10) << real ( c )
         << "  " << setw(10) << imag ( c )
         << "  " << setw(10) << c_norm << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests CAXPY.
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
{
# define N 5

  int i;
  complex <float> s;
  complex <float> x[N] = {
    complex <float> (  2.0, -1.0 ),
    complex <float> ( -4.0, -2.0 ),
    complex <float> (  3.0,  1.0 ),
    complex <float> (  2.0,  2.0 ),
    complex <float> ( -1.0, -1.0 ) };
  complex <float> y[N] = {
    complex <float> ( -1.0,  0.0 ),
    complex <float> (  0.0, -3.0 ),
    complex <float> (  4.0,  0.0 ),
    complex <float> ( -3.0,  4.0 ),
    complex <float> ( -2.0,  0.0 ) };

  cout << "\n";
  cout << "TEST03\n";
  cout << "  CAXPY adds a multiple of one complex vector to another.\n";

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }

  s = complex <float> ( 0.50, -1.00 );

  cout << "\n";
  cout << "  The scalar multiplier is: " << s << "\n";

  caxpy ( N, s, x, 1, y, 1 );

  cout << "\n";
  cout << "  A * X + Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests CCOPY.
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
{
# define N1 5
# define N2 5
# define N 10

  complex <float> a[N1*N2];
  int i;
  int j;
  complex <float> x[N];
  complex <float> y[N];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  CCOPY copies one complex vector into another.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 10 * ( i + 1 ), i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <float> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1] = complex <float> ( 10 * ( i + 1 ),  ( j + 1 ) );
    }
  }

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }
  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }
  cout << "\n";
  cout << "  A =\n";
  cout << "\n";
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      cout << "  " << setw(5) << real ( a[i+j*N1] )
           << "  " << setw(5) << imag ( a[i+j*N1] ) << "\n";
    }
    cout << "\n";
  }

  ccopy ( 5, x, 1, y, 1 );
  cout << "\n";
  cout << "  CCOPY ( 5, X, 1, Y, 1 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <float> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  ccopy ( 3, x, 2, y, 3 );

  cout << "\n";
  cout << "  CCOPY ( 3, X, 2, Y, 3 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }

  ccopy ( 5, x, 1, a, 1 );

  cout << "\n";
  cout << "  CCOPY ( 5, X, 1, A, 1 )\n";
  cout << "\n";
  cout << "\n";
  cout << "  A =\n";
  cout << "\n";
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      cout << "  " << setw(5) << real ( a[i+j*N1] )
           << "  " << setw(5) << imag ( a[i+j*N1] ) << "\n";
    }
    cout << "\n";
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1] = complex <float> ( 10 * ( i + 1 ),  ( j + 1 ) );
    }
  }

  ccopy ( 5, x, 2, a, 5 );

  cout << "\n";
  cout << "  CCOPY ( 5, X, 2, A, 5 )\n";
  cout << "\n";
  cout << "  A =\n";
  cout << "\n";
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      cout << "  " << setw(5) << real ( a[i+j*N1] )
           << "  " << setw(5) << imag ( a[i+j*N1] ) << "\n";
    }
    cout << "\n";
  }

  return;
# undef N
# undef N1
# undef N2
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests CDOTC.
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
{
# define N 5

  int i;
  complex <float> x_norm;
  complex <float> xy_dot;
  complex <float> x[N] = {
    complex <float> (  2.0, -1.0 ),
    complex <float> ( -4.0, -2.0 ),
    complex <float> (  3.0,  1.0 ),
    complex <float> (  2.0,  2.0 ),
    complex <float> ( -1.0, -1.0 ) };
  complex <float> y[N] = {
    complex <float> ( -1.0,  0.0 ),
    complex <float> (  0.0, -3.0 ),
    complex <float> (  4.0,  0.0 ),
    complex <float> ( -3.0,  4.0 ),
    complex <float> ( -2.0,  0.0 ) };

  cout << "\n";
  cout << "TEST05\n";
  cout << "  CDOTC computes the conjugated dot product of\n";
  cout << "  two complex vectors.\n";

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  x_norm = cdotc ( N, x, 1, x, 1 );

  cout << "\n";
  cout << "  The square of the norm of X, computed as\n";
  cout << "  CDOTC(X,X) = " << x_norm << "\n";

  xy_dot = cdotc ( N, x, 1, y, 1 );

  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }
  cout << "\n";
  cout << "  The dot product X.Y* is " << xy_dot << "\n";

  return;
# undef N
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests CDOTU.
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
{
# define N 5

  int i;
  complex <float> x_norm;
  complex <float> xy_dot;
  complex <float> x[N] = {
    complex <float> (  2.0, -1.0 ),
    complex <float> ( -4.0, -2.0 ),
    complex <float> (  3.0,  1.0 ),
    complex <float> (  2.0,  2.0 ),
    complex <float> ( -1.0, -1.0 ) };
  complex <float> y[N] = {
    complex <float> ( -1.0,  0.0 ),
    complex <float> (  0.0, -3.0 ),
    complex <float> (  4.0,  0.0 ),
    complex <float> ( -3.0,  4.0 ),
    complex <float> ( -2.0,  0.0 ) };

  cout << "\n";
  cout << "TEST06\n";
  cout << "  CDOTU computes the unconjugated dot product of\n";
  cout << "  two complex vectors.\n";

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  x_norm = cdotu ( N, x, 1, x, 1 );

  cout << "\n";
  cout << "  The unconjugated dot product ( X dot X )\n";
  cout << "  (which is NOT the square of the norm of X!):\n";
  cout << "  CDOTU(X,X) = " << x_norm << "\n";

  xy_dot = cdotu ( N, x, 1, y, 1 );

  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }

  cout << "\n";
  cout << "  The dot product ( X dot Y ) is " << xy_dot << "\n";

  return;
# undef N
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests CMACH.
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
{
  cout << "\n";
  cout << "TEST07\n";
  cout << "  CMACH computes several machine-dependent\n";
  cout << "  complex arithmetic parameters.\n";

  cout << "\n";
  cout << "  CMACH(1)  = machine epsilon = " << cmach ( 1 ) << "\n";
  cout << "  CMACH(2)  = a tiny value    = " << cmach ( 2 ) << "\n";
  cout << "  CMACH(3)  = a huge value    = " << cmach ( 3 ) << "\n";

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests CROTG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> a;
  complex <float> b;
  float c;
  complex <float> r;
  complex <float> s;
  complex <float> sa;
  complex <float> sb;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  CROTG generates a complex Givens rotation\n";
  cout << "    (  C  S ) * ( A ) = ( R )\n";
  cout << "    ( -S  C )   ( B )   ( 0 )\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = c4_uniform_01 ( &seed );
    b = c4_uniform_01 ( &seed );

    sa = a;
    sb = b;

    crotg ( &sa, sb, &c, &s );

    r = sa;

    cout << "\n";
    cout << "  A =  " << a << "\n";
    cout << "  B =  " << b << "\n";
    cout << "  C =  " << c << "\n";
    cout << "  S =  " << s << "\n";
    cout << "  R =  " << r << "\n";
    cout << "         C *A+S*B = " <<         c   * a + s * b << "\n";
    cout << "  -conjg(S)*A+C*B = " << -conj ( s ) * a + c * b << "\n";
  }

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests CSCAL.
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
{
# define N 6

  complex <float> da;
  int i;
  complex <float> x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 10 * ( i + 1 ), i + 1 );
  }

  cout << "\n";
  cout << "TEST09\n";
  cout << "  CSCAL multiplies a complex scalar times a vector.\n";

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  da = complex <float> ( 5.0, 0.0 );
  cscal ( N, da, x, 1 );
  cout << "\n";
  cout << "  CSCAL ( N, (" << da << "), X, 1 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 10 * ( i + 1 ), i + 1 );
  }

  da = complex <float> ( -2.0, 1.0 );
  cscal ( 3, da, x, 2 );
  cout << "\n";
  cout << "  CSCAL ( 3, (" << da << "), X, 2 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests CSIGN1.
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
{
  complex <float> c1;
  complex <float> c2;
  complex <float> c3;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  CSIGN1 ( C1, C2 ) transfers the sign of complex C2\n";
  cout << "  to the CABS1 magnitude of C1.\n";
  cout << "\n";
  cout <<
    "           C1                    C2                    C3\n";
  cout <<
    "  --------------------  --------------------  --------------------\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    c1 = complex <float> ( 5.0 ) * c4_uniform_01 ( &seed );
    c2 = complex <float> ( 5.0 ) * c4_uniform_01 ( &seed );
    c3 = csign1 ( c1, c2 );

    cout << "  " << setw(20) << c1
         << "  " << setw(20) << c2
         << "  " << setw(20) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests CSIGN2.
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
{
  complex <float> c1;
  complex <float> c2;
  complex <float> c3;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  CSIGN2 ( C1, C2 ) transfers the sign of complex C2\n";
  cout << "  to the CABS2 magnitude of C1.\n";
  cout << "\n";
  cout <<
    "           C1                    C2                    C3\n";
  cout <<
    "  --------------------  --------------------  --------------------\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    c1 = complex <float> ( 5.0 ) * c4_uniform_01 ( &seed );
    c2 = complex <float> ( 5.0 ) * c4_uniform_01 ( &seed );
    c3 = csign2 ( c1, c2 );

    cout << "  " << setw(20) << c1
         << "  " << setw(20) << c2
         << "  " << setw(20) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests CSROT.
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
{
# define N 6

  float c;
  int i;
  float s;
  complex <float> x[N];
  complex <float> y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 10 * ( i + 1 ), i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <float> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  cout << "\n";
  cout << "TEST12\n";
  cout << "  CSROT carries out a Givens rotation\n";
  cout << "  on a complex vector.\n";
  cout << "\n";
  cout << "  X and Y\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(20) << x[i]
         << "  " << setw(20) << y[i] << "\n";
  }

  c = 0.5;
  s = sqrt ( 1.0 - c * c );
  csrot ( N, x, 1, y, 1, c, s );
  cout << "\n";
  cout << "  CSROT ( N, X, 1, Y, 1, " << c << "," << s << " )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(20) << x[i]
         << "  " << setw(20) << y[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests CSSCAL.
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
{
# define N 6

  float da;
  int i;
  complex <float> x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 10 * ( i + 1 ), i + 1 );
  }

  cout << "\n";
  cout << "TEST13\n";
  cout << "  CSSCAL multiplies a real scalar times a complex vector.\n";

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  da = 5.0;
  csscal ( N, da, x, 1 );
  cout << "\n";
  cout << "  CSSCAL ( N, " << da << ", X, 1 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 10 * ( i + 1 ), i + 1 );
  }

  da = -2.0;
  csscal ( 3, da, x, 2 );
  cout << "\n";
  cout << "  CSSCAL ( 3, " << da << ", X, 2 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests CSWAP.
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
{
# define N 5

  int i;
  complex <float> x[N];
  complex <float> y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 10 * ( i + 1 ), i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <float> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  cout << "\n";
  cout << "TEST14\n";
  cout << "  CSWAP swaps two complex vectors.\n";
  cout << "\n";
  cout << "  X and Y\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(20) << x[i]
         << "  " << setw(20) << y[i] << "\n";
  }

  cswap ( N, x, 1, y, 1 );
  cout << "\n";
  cout << "  CSWAP ( N, X, 1, Y, 1 )\n";
  cout << "\n";
  cout << "  X and Y\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(20) << x[i]
         << "  " << setw(20) << y[i] << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 10 * ( i + 1 ), i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <float> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  cswap ( 3, x, 2, y, 1 );
  cout << "\n";
  cout << "  CSWAP ( 3, X, 2, Y, 1 )\n";
  cout << "\n";
  cout << "  X and Y\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(20) << x[i]
         << "  " << setw(20) << y[i] << "\n";
  }

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests ICAMAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  int incx;
  complex <float> x[N] = {
    complex <float>(  2.0, -1.0 ),
    complex <float>( -4.0, -2.0 ),
    complex <float>(  3.0,  1.0 ),
    complex <float>(  2.0,  2.0 ),
    complex <float>( -1.0, -1.0 ) };

  cout << "\n";
  cout << "TEST15\n";
  cout << "  ICAMAX returns the index of the entry of\n";
  cout << "  maximum magnitude in a complex vector.\n";

  cout << "\n";
  cout << "  The entries and CABS1 magnitudes:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(16) << x[i]
         << "  " << setw(8)  << cabs1 ( x[i] ) << "\n";
  }

  incx = 1;

  i = icamax ( N, x, incx );

  cout << "\n";
  cout << "  The index of maximum magnitude = " << i << "\n";
  cout << "\n";
  cout << "  Note that this is a 1-based index.\n";
  cout << "  Note that the L1 norm is used.\n";

  return;
# undef N
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests SCASUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2006
//
//  Author:
//
//    John Burkardt
//
{
# define MA 5
# define NA 4
# define NX 8

  complex <float> a[MA*NA] = {
    complex <float> ( -3.0,  4.0 ),
    complex <float> (  2.0,  0.0 ),
    complex <float> (  3.0, -4.0 ),
    complex <float> (  2.0,  0.0 ),
    complex <float> (  2.0, -1.0 ),
    complex <float> ( -1.0,  1.0 ),
    complex <float> (  0.0,  5.0 ),
    complex <float> ( -4.0, -2.0 ),
    complex <float> ( -4.0,  1.0 ),
    complex <float> ( -4.0, -3.0 ),
    complex <float> (  0.0, -2.0 ),
    complex <float> (  1.0,  3.0 ),
    complex <float> ( -3.0,  3.0 ),
    complex <float> ( -3.0,  3.0 ),
    complex <float> ( -1.0, -2.0 ),
    complex <float> ( -1.0,  2.0 ),
    complex <float> (  2.0, -4.0 ),
    complex <float> (  0.0, -1.0 ),
    complex <float> (  0.0, -1.0 ),
    complex <float> ( -2.0,  4.0 ) };
  int i;
  int j;
  complex <float> x[NX] = {
    complex <float> (  2.0, -1.0 ),
    complex <float> ( -4.0, -2.0 ),
    complex <float> (  3.0,  1.0 ),
    complex <float> (  2.0,  2.0 ),
    complex <float> ( -1.0, -1.0 ),
    complex <float> ( -1.0,  0.0 ),
    complex <float> (  0.0, -3.0 ),
    complex <float> (  4.0,  0.0 ) };

  cout << "\n";
  cout << "TEST16\n";
  cout << "  SCASUM adds the absolute values of elements\n";
  cout << "  of a complex vector.\n";
  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < NX; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(12) << x[i] << "\n";
  }

  cout << "\n";
  cout << "  SCASUM ( NX,   X, 1    ) = "
    << scasum ( NX,   x, 1 ) << "\n";
  cout << "  SCASUM ( NX/2, X, 2    ) = "
    << scasum ( NX/2, x, 2 ) << "\n";
  cout << "  SCASUM ( 2,    X, NX/2 ) = "
    << scasum ( 2,    x, NX/2 ) << "\n";

  cout << "\n";
  cout << "  Demonstrate with a matrix A:\n";
  cout << "\n";
  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      cout << "  " << setw(12) << a[i+j*MA];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  SCASUM ( MA, A[1,2], 1 )   = "
    << scasum ( MA, a+0+1*MA, 1 ) << "\n";
  cout << "  SCASUM ( NA, A[2,1], MA ) = "
    << scasum ( NA, a+1+0*MA, MA ) << "\n";

  return;
# undef MA
# undef NA
# undef NX
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests SCNRM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  int incx;
  float norm;
  complex <float> x[N] = {
    complex <float> (  2.0, -1.0 ),
    complex <float> ( -4.0, -2.0 ),
    complex <float> (  3.0,  1.0 ),
    complex <float> (  2.0,  2.0 ),
    complex <float> ( -1.0, -1.0 ) };

  cout << "\n";
  cout << "TEST17\n";
  cout << "  SCNRM2 returns the Euclidean norm of a complex vector.\n";

  cout << "\n";
  cout << "  The vector X:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(16) << x[i] << "\n";
  }

  incx = 1;
  norm = scnrm2 ( N, x, incx );

  cout << "\n";
  cout << "  The L2 norm of X is " << norm << "\n";

  return;
# undef N
}
//****************************************************************************80

complex <float> c4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4_UNIFORM_01 returns a unit complex pseudorandom number.
//
//  Discussion:
//
//    The angle should be uniformly distributed between 0 and 2 * PI,
//    the square root of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
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
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4_UNIFORM_01, a pseudorandom complex value.
//
{
  float r;
  int k;
  float pi = 3.1415926;
  float theta;
  complex <float> value;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = sqrt ( ( float ) ( ( double ) ( *seed ) * 4.656612875E-10 ) );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  theta = 2.0 * pi * ( float )
    ( ( double ) ( *seed ) * 4.656612875E-10 );

  value = complex <float> ( r * cos ( theta ), r * sin ( theta ) );

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

