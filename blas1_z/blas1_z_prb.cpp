# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "blas0.hpp"
# include "blas1_z.hpp"

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

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLAS1_Z_PRB.
//
//  Discussion:
//
//    BLAS1_Z_PRB tests the BLAS1_Z library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BLAS1_Z_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the BLAS1_Z library,\n";
 
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
  cout << "BLAS1_Z_PRB:\n";
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
//    TEST01 tests DZASUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MA 5
# define NA 4
# define NX 8

  complex <double> a[MA*NA] = {
    complex <double> ( -3.0,  4.0 ), 
    complex <double> (  2.0,  0.0 ), 
    complex <double> (  3.0, -4.0 ), 
    complex <double> (  2.0,  0.0 ), 
    complex <double> (  2.0, -1.0 ), 
    complex <double> ( -1.0,  1.0 ), 
    complex <double> (  0.0,  5.0 ), 
    complex <double> ( -4.0, -2.0 ), 
    complex <double> ( -4.0,  1.0 ), 
    complex <double> ( -4.0, -3.0 ), 
    complex <double> (  0.0, -2.0 ), 
    complex <double> (  1.0,  3.0 ), 
    complex <double> ( -3.0,  3.0 ), 
    complex <double> ( -3.0,  3.0 ), 
    complex <double> ( -1.0, -2.0 ), 
    complex <double> ( -1.0,  2.0 ), 
    complex <double> (  2.0, -4.0 ), 
    complex <double> (  0.0, -1.0 ), 
    complex <double> (  0.0, -1.0 ), 
    complex <double> ( -2.0,  4.0 ) };
  int i;
  int j;
  complex <double> x[NX] = {
    complex <double> (  2.0, -1.0 ), 
    complex <double> ( -4.0, -2.0 ), 
    complex <double> (  3.0,  1.0 ), 
    complex <double> (  2.0,  2.0 ), 
    complex <double> ( -1.0, -1.0 ), 
    complex <double> ( -1.0,  0.0 ), 
    complex <double> (  0.0, -3.0 ), 
    complex <double> (  4.0,  0.0 ) };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  DZASUM adds the absolute values of elements\n";
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
  cout << "  DZASUM ( NX,   X, 1    ) = "
    << dzasum ( NX,   x, 1 ) << "\n";
  cout << "  DZASUM ( NX/2, X, 2    ) = "
    << dzasum ( NX/2, x, 2 ) << "\n";
  cout << "  DZASUM ( 2,    X, NX/2 ) = "
    << dzasum ( 2,    x, NX/2 ) << "\n";

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
  cout << "  DZASUM ( MA, A[1,2], 1 )   = "
    << dzasum ( MA, a+0+1*MA, 1 ) << "\n";
  cout << "  DZASUM ( NA, A[2,1], MA ) = "
    << dzasum ( NA, a+1+0*MA, MA ) << "\n";

  return;
# undef MA
# undef NA
# undef NX
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests DZNRM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  int incx;
  double norm;
  complex <double> x[N] = {
    complex <double> (  2.0, -1.0 ), 
    complex <double> ( -4.0, -2.0 ), 
    complex <double> (  3.0,  1.0 ), 
    complex <double> (  2.0,  2.0 ), 
    complex <double> ( -1.0, -1.0 ) };

  cout << "\n";
  cout << "TEST02\n";
  cout << "  DZNRM2 returns the Euclidean norm of a complex vector.\n";
 
  cout << "\n";
  cout << "  The vector X:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(16) << x[i] << "\n";
  }

  incx = 1;
  norm = dznrm2 ( N, x, incx );

  cout << "\n"; 
  cout << "  The L2 norm of X is " << norm << "\n";

  return;
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests IZAMAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  int incx;
  complex <double> x[N] = {
    complex <double>(  2.0, -1.0 ), 
    complex <double>( -4.0, -2.0 ), 
    complex <double>(  3.0,  1.0 ), 
    complex <double>(  2.0,  2.0 ), 
    complex <double>( -1.0, -1.0 ) };

  cout << "\n";
  cout << "TEST03\n";
  cout << "  IZAMAX returns the index of the entry of\n";
  cout << "  maximum magnitude in a complex vector.\n";
 
  cout << "\n";
  cout << "  The entries and ZABS1 magnitudes:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(16) << x[i] 
         << "  " << setw(8)  << zabs1 ( x[i] ) << "\n";
  }

  incx = 1;

  i = izamax ( N, x, incx );

  cout << "\n";
  cout << "  The index of maximum magnitude = " << i << "\n";
  cout << "\n";
  cout << "  Note that this is a 1-based index.\n";
  cout << "  Note that the L1 norm is used.\n";

  return;
# undef N
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests ZABS1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c;
  double c_norm;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  ZABS1 returns the L1 norm of a complex number.\n";
  cout << "\n";
  cout << "      Real      Imaginary\n";
  cout << "      Part      Part           ZABS1(Z)\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
//
//  The compiler seems to be unhappy with the statement
//
//    c = 5.0 * c8_uniform_01 ( seed );
//
//  Poor compiler.
//
    c = c8_uniform_01 ( seed );
    c = complex <double> ( 5.0 ) * c;

    c_norm = zabs1 ( c );

    cout << "  " << setw(10) << real ( c )
         << "  " << setw(10) << imag ( c )
         << "  " << setw(10) << c_norm << "\n";
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests ZABS2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c;
  double c_norm;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  ZABS2 returns the L2 norm of a complex number.\n";
  cout << "\n";
  cout << "      Real      Imaginary\n";
  cout << "      Part      Part           ZABS2(Z\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    c = c8_uniform_01 ( seed );
    c = complex <double> ( 5.0 ) * c;

    c_norm = zabs2 ( c );

    cout << "  " << setw(10) << real ( c )
         << "  " << setw(10) << imag ( c )
         << "  " << setw(10) << c_norm << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests ZAXPY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  complex <double> s;
  complex <double> x[N] = {
    complex <double> (  2.0, -1.0 ), 
    complex <double> ( -4.0, -2.0 ), 
    complex <double> (  3.0,  1.0 ), 
    complex <double> (  2.0,  2.0 ), 
    complex <double> ( -1.0, -1.0 ) };
  complex <double> y[N] = {
    complex <double> ( -1.0,  0.0 ), 
    complex <double> (  0.0, -3.0 ),
    complex <double> (  4.0,  0.0 ), 
    complex <double> ( -3.0,  4.0 ), 
    complex <double> ( -2.0,  0.0 ) };

  cout << "\n";
  cout << "TEST06\n";
  cout << "  ZAXPY adds a multiple of one complex vector to another.\n";
 
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

  s = complex <double> ( 0.50, -1.00 );

  cout << "\n";
  cout << "  The scalar multiplier is: " << s << "\n";

  zaxpy ( N, s, x, 1, y, 1 );

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

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests ZCOPY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 5
# define N2 5
# define N 10

  complex <double> a[N1*N2];
  int i;
  int j;
  complex <double> x[N];
  complex <double> y[N];

  cout << "\n";
  cout << "TEST07\n";
  cout << "  ZCOPY copies one complex vector into another.\n";
 
  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( 10 * ( i + 1 ), i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <double> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1] = complex <double> ( 10 * ( i + 1 ),  ( j + 1 ) );
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

  zcopy ( 5, x, 1, y, 1 );
  cout << "\n";
  cout << "  ZCOPY ( 5, X, 1, Y, 1 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <double> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  zcopy ( 3, x, 2, y, 3 );

  cout << "\n";
  cout << "  ZCOPY ( 3, X, 2, Y, 3 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( y[i] )
         << "  " << setw(6) << imag ( y[i] ) << "\n";
  }

  zcopy ( 5, x, 1, a, 1 );

  cout << "\n";
  cout << "  ZCOPY ( 5, X, 1, A, 1 )\n";
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
      a[i+j*N1] = complex <double> ( 10 * ( i + 1 ),  ( j + 1 ) );
    }
  }

  zcopy ( 5, x, 2, a, 5 );

  cout << "\n";
  cout << "  ZCOPY ( 5, X, 2, A, 5 )\n";
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

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests ZDOTC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  complex <double> x_norm;
  complex <double> xy_dot;
  complex <double> x[N] = {
    complex <double> (  2.0, -1.0 ), 
    complex <double> ( -4.0, -2.0 ), 
    complex <double> (  3.0,  1.0 ), 
    complex <double> (  2.0,  2.0 ), 
    complex <double> ( -1.0, -1.0 ) };
  complex <double> y[N] = {
    complex <double> ( -1.0,  0.0 ), 
    complex <double> (  0.0, -3.0 ), 
    complex <double> (  4.0,  0.0 ), 
    complex <double> ( -3.0,  4.0 ), 
    complex <double> ( -2.0,  0.0 ) };

  cout << "\n";
  cout << "TEST08\n";
  cout << "  ZDOTC computes the conjugated dot product of\n";
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

  x_norm = zdotc ( N, x, 1, x, 1 );

  cout << "\n";
  cout << "  The square of the norm of X, computed as\n";
  cout << "  ZDOTC(X,X) = " << x_norm << "\n";

  xy_dot = zdotc ( N, x, 1, y, 1 );

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

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests ZDOTU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  complex <double> x_norm;
  complex <double> xy_dot;
  complex <double> x[N] = {
    complex <double> (  2.0, -1.0 ), 
    complex <double> ( -4.0, -2.0 ), 
    complex <double> (  3.0,  1.0 ), 
    complex <double> (  2.0,  2.0 ), 
    complex <double> ( -1.0, -1.0 ) };
  complex <double> y[N] = {
    complex <double> ( -1.0,  0.0 ), 
    complex <double> (  0.0, -3.0 ), 
    complex <double> (  4.0,  0.0 ), 
    complex <double> ( -3.0,  4.0 ), 
    complex <double> ( -2.0,  0.0 ) };

  cout << "\n";
  cout << "TEST09\n";
  cout << "  ZDOTU computes the unconjugated dot product of\n";
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

  x_norm = zdotu ( N, x, 1, x, 1 );

  cout << "\n";
  cout << "  The unconjugated dot product ( X dot X )\n";
  cout << "  (which is NOT the square of the norm of X!):\n";
  cout << "  ZDOTU(X,X) = " << x_norm << "\n";

  xy_dot = zdotu ( N, x, 1, y, 1 );

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

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests ZDROT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  double c;
  int i;
  double s;
  complex <double> x[N];
  complex <double> y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( 10 * ( i + 1 ), i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <double> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  cout << "\n";
  cout << "TEST10\n";
  cout << "  ZDROT carries out a Givens rotation\n";
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
  zdrot ( N, x, 1, y, 1, c, s );
  cout << "\n";
  cout << "  ZDROT ( N, X, 1, Y, 1, " << c << "," << s << " )\n";
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

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests ZDSCAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  double da;
  int i;
  complex <double> x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( 10 * ( i + 1 ), i + 1 );
  }

  cout << "\n";
  cout << "TEST11\n";
  cout << "  ZDSCAL multiplies a real scalar times a complex vector.\n";

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
  zdscal ( N, da, x, 1 );
  cout << "\n";
  cout << "  ZDSCAL ( N, " << da << ", X, 1 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( 10 * ( i + 1 ), i + 1 );
  }

  da = -2.0;
  zdscal ( 3, da, x, 2 );
  cout << "\n";
  cout << "  ZDSCAL ( 3, " << da << ", X, 2 )\n";
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

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests ZMACH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST12\n";
  cout << "  ZMACH computes several machine-dependent\n";
  cout << "  complex arithmetic parameters.\n";

  cout << "\n";
  cout << "  ZMACH(1)  = machine epsilon = " << zmach ( 1 ) << "\n";
  cout << "  ZMACH(2)  = a tiny value    = " << zmach ( 2 ) << "\n";
  cout << "  ZMACH(3)  = a huge value    = " << zmach ( 3 ) << "\n";

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests ZROTG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> a;
  complex <double> b;
  double c;
  complex <double> r;
  complex <double> s;
  complex <double> sa;
  complex <double> sb;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  ZROTG generates a complex Givens rotation\n";
  cout << "    (  C  S ) * ( A ) = ( R )\n";
  cout << "    ( -S  C )   ( B )   ( 0 )\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = c8_uniform_01 ( seed );
    b = c8_uniform_01 ( seed );

    sa = a;
    sb = b;

    zrotg ( &sa, sb, &c, &s );

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

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests ZSCAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  complex <double> da;
  int i;
  complex <double> x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( 10 * ( i + 1 ), i + 1 );
  }

  cout << "\n";
  cout << "TEST14\n";
  cout << "  ZSCAL multiplies a complex scalar times a vector.\n";

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  da = complex <double> ( 5.0, 0.0 );
  zscal ( N, da, x, 1 );
  cout << "\n";
  cout << "  ZSCAL ( N, (" << da << "), X, 1 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << real ( x[i] )
         << "  " << setw(6) << imag ( x[i] ) << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( 10 * ( i + 1 ), i + 1 );
  }

  da = complex <double> ( -2.0, 1.0 );
  zscal ( 3, da, x, 2 );
  cout << "\n";
  cout << "  ZSCAL ( 3, (" << da << "), X, 2 )\n";
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

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests ZSIGN1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  ZSIGN1 ( C1, C2 ) transfers the sign of complex C2\n";
  cout << "  to the ZABS1 magnitude of C1.\n";
  cout << "\n";
  cout << 
    "           C1                    C2                    C3\n";
  cout << 
    "  --------------------  --------------------  --------------------\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    c1 = complex <double> ( 5.0 ) * c8_uniform_01 ( seed );
    c2 = complex <double> ( 5.0 ) * c8_uniform_01 ( seed );
    c3 = zsign1 ( c1, c2 );

    cout << "  " << setw(20) << c1
         << "  " << setw(20) << c2
         << "  " << setw(20) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests ZSIGN2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST16\n";
  cout << "  ZSIGN2 ( C1, C2 ) transfers the sign of complex C2\n";
  cout << "  to the ZABS2 magnitude of C1.\n";
  cout << "\n";
  cout << 
    "           C1                    C2                    C3\n";
  cout << 
    "  --------------------  --------------------  --------------------\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    c1 = complex <double> ( 5.0 ) * c8_uniform_01 ( seed );
    c2 = complex <double> ( 5.0 ) * c8_uniform_01 ( seed );
    c3 = zsign2 ( c1, c2 );

    cout << "  " << setw(20) << c1
         << "  " << setw(20) << c2
         << "  " << setw(20) << c3 << "\n";
  }

  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests ZSWAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  complex <double> x[N];
  complex <double> y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( 10 * ( i + 1 ), i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <double> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  cout << "\n";
  cout << "TEST17\n";
  cout << "  ZSWAP swaps two complex vectors.\n";
  cout << "\n";
  cout << "  X and Y\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(20) << x[i]
         << "  " << setw(20) << y[i] << "\n";
  }

  zswap ( N, x, 1, y, 1 );
  cout << "\n";
  cout << "  ZSWAP ( N, X, 1, Y, 1 )\n";
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
    x[i] = complex <double> ( 10 * ( i + 1 ), i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = complex <double> ( 20 * ( i + 1 ), 2 * ( i + 1 ) );
  }

  zswap ( 3, x, 2, y, 1 );
  cout << "\n";
  cout << "  ZSWAP ( 3, X, 2, Y, 1 )\n";
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

