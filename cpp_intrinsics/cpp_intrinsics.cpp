# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void test_abs ( );
void test_acos ( );
void test_asin ( );
void test_atan ( );
void test_atan2 ( );
void test_ceil ( );
void test_cos ( );
void test_cosh ( );
void test_exp ( );
void test_fabs ( );
void test_floor ( );
void test_fmod ( );
void test_frexp ( );
void test_ldexp ( );
void test_log ( );
void test_log10 ( );
void test_modf ( );
void test_pow ( );
void test_sin ( );
void test_sinh ( );
void test_sqrt ( );
void test_tan ( );
void test_tanh ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform ( int a, int b, int *seed );
int r4_nint ( float x );
float r4_uniform ( float b, float c, int *seed );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CPP_INTRINSICS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CPP_INTRINSICS:\n";
  cout << "  Test the C++ intrinsic library.\n";

  test_abs ( );
  test_acos ( );
  test_asin ( );
  test_atan ( );
  test_atan2 ( );
  test_ceil ( );
  test_cos ( );
  test_cosh ( );
  test_exp ( );
  test_fabs ( );
  test_floor ( );
  test_fmod ( );
  test_frexp ( );
  test_ldexp ( );
  test_log ( );
  test_log10 ( );
  test_modf ( );
  test_pow ( );
  test_sin ( );
  test_sinh ( );
  test_sqrt ( );
  test_tan ( );
  test_tanh ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CPP_INTRINSICS:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test_abs ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ABS tests ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int arg;
  int result;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST_ABS:\n";
  cout << "  Test ABS, which evaluates the absolute value of an int.\n";
  cout << "\n";
  cout << "         I         ABS(I)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    arg = i4_uniform ( -10000, +10000, &seed );
    result = abs ( arg );

    cout << "  " << setw(10) << arg
         << "  " << setw(10) << result << "\n";
  }
  return;
}
//****************************************************************************80

void test_acos ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ACOS tests ACOS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_ACOS:\n";
  cout << "  Test ACOS, which evaluates the arc-cosine function.\n";
  cout << "\n";
  cout << "         X          Y           Z\n";
  cout << "                 ACOS(X)      COS(Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -1.0, +1.0, &seed );
    y = acos ( x );
    z = cos ( y );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y 
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_asin ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ASIN tests ASIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_ASIN:\n";
  cout << "  Test ASIN, which evaluates the arc-sine function.\n";
  cout << "\n";
  cout << "         X          Y           Z\n";
  cout << "                 ASIN(X)      SIN(Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -1.0, +1.0, &seed );
    y = asin ( x );
    z = sin ( y );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y 
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_atan ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ATAN tests ATAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_ATAN:\n";
  cout << "  Test ATAN, which evaluates the arc-tangent function.\n";
  cout << "\n";
  cout << "         X          Y           Z\n";
  cout << "                 ATAN(X)      TAN(Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -1.0, +1.0, &seed );
    y = atan ( x );
    z = tan ( y );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y 
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_atan2 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ATAN2 tests ATAN2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x1;
  float x2;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_ATAN2:\n";
  cout << "  Test ATAN2, which evaluates the arc-tangent function.\n";
  cout << "\n";
  cout << "         X1        X2          X1/X2        Y           Z\n";
  cout << "                                         ATAN(X1,X2)  TAN(Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x1 = r4_uniform ( -1.0, +1.0, &seed );
    x2 = r4_uniform ( -1.0, +1.0, &seed );
    y = atan2 ( x1, x2 );
    z = tan ( y );
    cout << "  " << setw(10) << x1
         << "  " << setw(10) << x2
         << "  " << setw(10) << x1 / x2
         << "  " << setw(10) << y 
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_ceil ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_CEIL tests CEIL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_CEIL:\n";
  cout << "  Test CEIL, which evaluates the ceiling function.\n";
  cout << "\n";
  cout << "         X           Y\n";
  cout << "                    CEIL(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -20.0, +20.0, &seed );
    y = ceil ( x );

    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_cos ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_COS tests COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_COS:\n";
  cout << "  Test COS, which evaluates the cosine function.\n";
  cout << "\n";
  cout << "         X          Y\n";
  cout << "                 COS(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -pi, +pi, &seed );
    y = cos ( x );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_cosh ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_COSH tests COSH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_COSH:\n";
  cout << "  Test COSH, which evaluates the hyperbolic cosine function.\n";
  cout << "\n";
  cout << "         X          Y\n";
  cout << "                 COSH(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -5.0, +5.0, &seed );
    y = cosh ( x );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_exp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_EXP tests EXP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_EXP:\n";
  cout << "  Test EXP, which evaluates the exponential function.\n";
  cout << "\n";
  cout << "         X          Y           Z\n";
  cout << "                  EXP(X)      LOG(Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -5.0, +10.0, &seed );
    y = exp ( x );
    z = log ( y );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y 
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_fabs ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_FABS tests FABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_FABS:\n";
  cout << "  Test FABS, which evaluates the absolute value of a real quantity.\n";
  cout << "\n";
  cout << "         X           Y\n";
  cout << "                   FABS(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -10000.0, +10000.0, &seed );
    y = fabs ( x );

    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_floor ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_FLOOR tests FLOOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_FLOOR:\n";
  cout << "  Test FLOOR, which evaluates the floor function.\n";
  cout << "\n";
  cout << "         X           Y\n";
  cout << "                    FLOOR(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -20.0, +20.0, &seed );
    y = floor ( x );

    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_fmod ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_FMOD tests FMOD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x1;
  float x2;
  float y;
  float z;
  float z1;
  double z2;

  cout << "\n";
  cout << "TEST_FMOD:\n";
  cout << "  Test FMOD, which returns the remainder of X1 / X2.\n";
  cout << "\n";
  cout << "         X1        X2          X1/X2        Y           Z\n";
  cout << "                                         FMOD(X1,X2)  X2*MODF(X1/X2,*)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x1 = r4_uniform ( -5.0, +5.0, &seed );
    x2 = r4_uniform ( -5.0, +5.0, &seed );
    y = fmod ( x1, x2 );
    z1 = modf ( x1 / x2, &z2 ); 
    z = z1 * x2;
    cout << "  " << setw(10) << x1
         << "  " << setw(10) << x2
         << "  " << setw(10) << x1 / x2
         << "  " << setw(10) << y
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_frexp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_FREXP tests FREXP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_FREXP:\n";
  cout << "  Test FREXP, which splits X into a normalized fraction\n";
  cout << "  and a power of 2.\n";
  cout << "\n";
  cout << "        X           Y              N         Z\n";
  cout << "                                           Y*2^N\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -20.0, +20.0, &seed );
    y = frexp ( x, &n );
    z = y * pow ( 2.0, n );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y
         << "  " << setw(10) << n
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_ldexp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LDEXP tests LDEXP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int result;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_LDEXP:\n";
  cout << "  Test LDEXP, which evaluates X*2^N.\n";
  cout << "\n";
  cout << "         X           N           Y\n";
  cout << "                             LDEXP(X,N)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -1.0, +1.0, &seed );
    n = i4_uniform ( -10, +10, &seed );
    y = ldexp ( x, n );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << n
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_log ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LOG tests LOG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_LOG:\n";
  cout << "  Test LOG, which evaluates the logarithm function.\n";
  cout << "\n";
  cout << "         X          Y           Z\n";
  cout << "                  LOG(X)      EXP(Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( 0, +10000.0, &seed );
    y = log ( x );
    z = exp ( y );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y 
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_log10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LOG10 tests LOG10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_LOG10:\n";
  cout << "  Test LOG10, which evaluates the logarithm base 10 function.\n";
  cout << "\n";
  cout << "         X          Y           Z\n";
  cout << "                 LOG10(X)    POW(10,Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( 0, +10000.0, &seed );
    y = log10 ( x );
    z = pow ( 10.0, y );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y 
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_modf ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_MODF tests MODF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y1;
  double y2;
  float z;

  cout << "\n";
  cout << "TEST_MODF:\n";
  cout << "  Test MODF, which splits X into integer and fractional\n";
  cout << "  parts.\n";
  cout << "\n";
  cout << "        X           Y1          Y2         Z\n";
  cout << "                                         Y1+Y2\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -20.0, +20.0, &seed );
    y1 = modf ( x, &y2 );
    z = y1 + y2;
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y1
         << "  " << setw(10) << y2
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_pow ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_POW tests POW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x1;
  float x2;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_POW:\n";
  cout << "  Test POW, which evaluates the power function X1^X2.\n";
  cout << "\n";
  cout << "         X1        X2          Y\n";
  cout << "                           POW(X1,X2)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x1 = r4_uniform ( 0.0, +10.0, &seed );
    x2 = r4_uniform ( -2.0, +10.0, &seed );
    y = pow ( x1, x2 );
    cout << "  " << setw(10) << x1
         << "  " << setw(10) << x2 
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_sin ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_SIN tests SIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_SIN:\n";
  cout << "  Test SIN, which evaluates the sine function.\n";
  cout << "\n";
  cout << "         X          Y\n";
  cout << "                 SIN(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -pi, +pi, &seed );
    y = sin ( x );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_sinh ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_SINH tests SINH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_SINH:\n";
  cout << "  Test SINH, which evaluates the hyperbolic sine function.\n";
  cout << "\n";
  cout << "         X          Y\n";
  cout << "                 SINH(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -5.0, +5.0, &seed );
    y = sinh ( x );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_sqrt ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_SQRT tests SQRT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  cout << "\n";
  cout << "TEST_SQRT:\n";
  cout << "  Test SQRT, which evaluates the square root function.\n";
  cout << "\n";
  cout << "         X          Y           Z\n";
  cout << "                 SQRT(X)       Y*Y\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( 0.0, +100.0, &seed );
    y = sqrt ( x );
    z = y * y;
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y 
         << "  " << setw(10) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test_tan ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_TAN tests TAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_TAN:\n";
  cout << "  Test TAN, which evaluates the tangent function.\n";
  cout << "\n";
  cout << "         X          Y\n";
  cout << "                 TAN(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -pi/2.0, +pi/2.0, &seed );
    y = tan ( x );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

void test_tanh ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_TANH tests TANH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  cout << "\n";
  cout << "TEST_TANH:\n";
  cout << "  Test TANH, which evaluates the hyperbolic tangent function.\n";
  cout << "\n";
  cout << "         X          Y\n";
  cout << "                 TANH(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -5.0, +5.0, &seed );
    y = tanh ( x );
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y << "\n";
  }
  return;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( - x + 0.5 );
  }
  else
  {
    value =   ( int ) (  x + 0.5 );
  }
  return value;
}
//****************************************************************************80

float r4_uniform ( float a, float b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM returns a scaled pseudorandom R4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, float R4_UNIFORM, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  float value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  value = ( float ) ( *seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

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
