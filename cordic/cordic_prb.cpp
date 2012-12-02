# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "cordic.hpp"

int main ( );

void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test008 ( );
void test009 ( );
void test010 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CORDIC_PRB.
//
//  Discussion:
//
//    CORDIC_PRB calls the CORDIC routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "CORDIC_PRB:\n";
  cout << "  C++ version,\n";
  cout << "  Test the CORDIC library.\n";

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test008 ( );
  test009 ( );
  test010 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CORDIC_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test001 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST001 demonstrates the use of COSSIN_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double c1;
  double c2;
  double d;
  int n;
  int n_data;
  double s2;

  cout << "\n";
  cout << "TEST001:\n";
  cout << "  COSSIN_CORDIC computes the cosine and sine\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "          A        N      Cos(A)           Cos(A)           Difference\n";
  cout << "                          Tabulated        CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cos_values ( &n_data, &a, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      cossin_cordic ( a, n, &c2, &s2 );

      d = c1 - c2;

      cout
           << setw(12) << setprecision ( 4 ) << a  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << c1  << "  "
           << setw(16) << setprecision ( 8 ) << c2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST002 demonstrates the use of COSSIN_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double c2;
  double d;
  int n;
  int n_data;
  double s1;
  double s2;

  cout << "\n";
  cout << "TEST002:\n";
  cout << "  COSSIN_CORDIC computes the cosine and sine\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "          A        N      Sin(A)           Sin(A)           Difference\n";
  cout << "                          Tabulated        CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sin_values ( &n_data, &a, &s1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      cossin_cordic ( a, n, &c2, &s2 );

      d = s1 - s2;

      cout
           << setw(12) << setprecision ( 4 ) << a  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << s1  << "  "
           << setw(16) << setprecision ( 8 ) << s2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 demonstrates the use of ARCTAN_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a1;
  double a2;
  double d;
  int n;
  int n_data;
  double r;
  double s;
  int seed;
  double x;
  double y;
  double z;

  seed = 123456789;

  cout << "\n";
  cout << "TEST003:\n";
  cout << "  ARCTAN_CORDIC computes the arctangent of Y/X\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "      X      Y    N       ArcTan(Y/X) ArcTan(Y/X)      Difference\n";
  cout <<
    "                           Tabulated   CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arctan_values ( &n_data, &z, &a1 );

    if ( n_data == 0 )
    {
      break;
    }

    r = r8_uniform_01 ( &seed );

    x = r;
    y = r * z;

    s = r8_uniform_01 ( &seed );

    if ( s < 0.5 )
    {
      x = -x;
      y = -y;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      a2 = arctan_cordic ( x, y, n );

      d = a1 - a2;

      cout
           << setw(12) << setprecision ( 4 ) << x  << "  "
           << setw(12) << setprecision ( 4 ) << y  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << a1  << "  "
           << setw(16) << setprecision ( 8 ) << a2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test004 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST004 demonstrates the use of ARCCOS_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a1;
  double a2;
  double d;
  int n;
  int n_data;
  double t;

  cout << "\n";
  cout << "TEST004:\n";
  cout << "  ARCCOS_CORDIC computes the arccosine of T\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "      T    N        ArcCos(T)  ArcCos(T)      Difference\n";
  cout <<
    "                   Tabulated   CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arccos_values ( &n_data, &t, &a1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      a2 = arccos_cordic ( t, n );

      d = a1 - a2;

      cout
           << setw(12) << setprecision ( 4 ) << t  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << a1  << "  "
           << setw(16) << setprecision ( 8 ) << a2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 demonstrates the use of ARCSIN_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a1;
  double a2;
  double d;
  int n;
  int n_data;
  double t;

  cout << "\n";
  cout << "TEST005:\n";
  cout << "  ARCSIN_CORDIC computes the arcsine of T\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "      T    N        ArcSin(T)  ArcSin(T)      Difference\n";
  cout <<
    "                   Tabulated   CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arcsin_values ( &n_data, &t, &a1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      a2 = arcsin_cordic ( t, n );

      d = a1 - a2;

      cout
           << setw(12) << setprecision ( 4 ) << t  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << a1  << "  "
           << setw(16) << setprecision ( 8 ) << a2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 demonstrates the use of TAN_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double d;
  int n;
  int n_data;
  double t1;
  double t2;
  double theta;

  cout << "\n";
  cout << "TEST006:\n";
  cout << "  TAN_CORDIC computes the tangent of THETA\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "  THETA    N         Tan(THETA)  Tan(THETA)      Difference\n";
  cout <<
    "                     Tabulated   CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tan_values ( &n_data, &theta, &t1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      t2 = tan_cordic ( theta, n );

      d = t1 - t2;

      cout
           << setw(12) << setprecision ( 4 ) << theta  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << t1  << "  "
           << setw(16) << setprecision ( 8 ) << t2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test007 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 demonstrates the use of EXP_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double d;
  double fx1;
  double fx2;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST007:\n";
  cout << "  EXP_CORDIC computes the exponential function\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "    X      N           Exp(X)      Exp(X)        Difference\n";
  cout <<
    "                     Tabulated   CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    exp_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      fx2 = exp_cordic ( x, n );

      d = fx1 - fx2;

      cout
           << setw(12) << setprecision ( 4 ) << x  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << fx1  << "  "
           << setw(16) << setprecision ( 8 ) << fx2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test008 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST008 demonstrates the use of LN_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double d;
  double fx1;
  double fx2;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST008:\n";
  cout << "  LN_CORDIC computes the natural logarithm\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "    X      N            Ln(X)       Ln(X)        Difference\n";
  cout <<
    "                     Tabulated   CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    ln_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      fx2 = ln_cordic ( x, n );

      d = fx1 - fx2;

      cout
           << setw(12) << setprecision ( 4 ) << x  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << fx1  << "  "
           << setw(16) << setprecision ( 8 ) << fx2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test009 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST009 demonstrates the use of SQRT_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double d;
  double fx1;
  double fx2;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST009:\n";
  cout << "  SQRT_CORDIC computes the square root\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "    X      N          Sqrt(X)     Sqrt(X)        Difference\n";
  cout <<
    "                     Tabulated   CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sqrt_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      fx2 = sqrt_cordic ( x, n );

      d = fx1 - fx2;

      cout
           << setw(12) << setprecision ( 4 ) << x  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << fx1  << "  "
           << setw(16) << setprecision ( 8 ) << fx2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test010 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST010 demonstrates the use of CBRT_CORDIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double d;
  double fx1;
  double fx2;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST010:\n";
  cout << "  CBRT_CORDIC computes the cube root\n";
  cout << "  using the CORDIC algorithm.\n";
  cout << "\n";
  cout <<
    "    X      N          Cbrt(X)     Cbrt(X)        Difference\n";
  cout <<
    "                     Tabulated   CORDIC\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cbrt_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "\n";
    for ( n = 0; n <= 25; n = n + 5 )
    {
      fx2 = cbrt_cordic ( x, n );

      d = fx1 - fx2;

      cout
           << setw(12) << setprecision ( 4 ) << x  << "  "
           << setw(4)                        << n  << "  "                                 << "  "
           << setw(16) << setprecision ( 8 ) << fx1  << "  "
           << setw(16) << setprecision ( 8 ) << fx2  << "  "
           << setw(12) << setprecision ( 4 ) << d << "\n";
    }
  }
  return;
}
