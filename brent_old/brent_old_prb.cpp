# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <string>

using namespace std;

# include "brent_old.hpp"

int main ( );

void test_zero_all ( );
void test_zero_rc_all ( );
void test_local_min_all ( );
void test_local_min_rc_all ( );
void test_glomin_all ( );

void test_zero_one ( double a, double b, double t, double f ( double x ), 
  string title );
void test_zero_rc_one ( double a, double b, double t, double f ( double x ), 
  string title );
void test_local_min_one ( double a, double b, double t, double f ( double x ), 
  string title );
void test_local_min_rc_one ( double a, double b, double t, 
  double f ( double x ), string title );
void test_glomin_one ( double a, double b, double c, double m, double e, 
  double t, double f ( double x ), string title );

double f_01 ( double x );
double f_02 ( double x );
double f_03 ( double x );
double f_04 ( double x );
double f_05 ( double x );
double g_01 ( double x );
double g_02 ( double x );
double g_03 ( double x );
double g_04 ( double x );
double g_05 ( double x );
double h_01 ( double x );
double h_02 ( double x );
double h_03 ( double x );
double h_04 ( double x );
double h_05 ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BRENT_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BRENT_OLD_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BRENT_OLD library.\n";

  test_zero_all ( );
  test_zero_rc_all ( );
  test_local_min_all ( );
  test_local_min_rc_all ( );
  test_glomin_all ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BRENT_OLD_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test_zero_all ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ZERO_ALL tests Brent's zero finding routine on all test functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double t;
  double z;

  cout << "\n";
  cout << "TEST_ZERO_ALL\n";
  cout << "  Test the Brent ZERO routine, which seeks\n";
  cout << "  a root of a function F(X)\n";
  cout << "  in an interval [A,B].\n";

  t = r8_epsilon ( );

  a = 1.0;
  b = 2.0;

  test_zero_one ( a, b, t, f_01, 
    "f_01(x) = sin ( x ) - x / 2" );

  a = 0.0;
  b = 1.0;

  test_zero_one ( a, b, t, f_02, 
    "f_02(x) = 2 * x - exp ( - x )" );

  a = -1.0;
  b =  0.5;

  test_zero_one ( a, b, t, f_03, 
    "f_03(x) = x * exp ( - x )" );

  a =  0.0001;
  b =  20.0;

  test_zero_one ( a, b, t, f_04, 
    "f_04(x) = exp ( x ) - 1 / ( 100 * x * x )" );

  a = -5.0;
  b =  2.0;

  test_zero_one ( a, b, t, f_05, 
    "f_05(x) = (x+3) * (x-1) * (x-1)" );

  return;
}
//****************************************************************************80

void test_zero_rc_all ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ZERO_RC_ALL tests ZERO_RC on all test functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double t;
  double z;

  cout << "\n";
  cout << "TEST_ZERO_RC_ALL\n";
  cout << "  Test the ZERO_RC routine, which seeks\n";
  cout << "  a root of a function F(X)\n";
  cout << "  in an interval [A,B].\n";

  t = r8_epsilon ( );

  a = 1.0;
  b = 2.0;

  test_zero_rc_one ( a, b, t, f_01, 
    "f_01(x) = sin ( x ) - x / 2" );

  a = 0.0;
  b = 1.0;

  test_zero_rc_one ( a, b, t, f_02, 
    "f_02(x) = 2 * x - exp ( - x )" );

  a = -1.0;
  b =  0.5;

  test_zero_rc_one ( a, b, t, f_03, 
    "f_03(x) = x * exp ( - x )" );

  a =  0.0001;
  b =  20.0;

  test_zero_rc_one ( a, b, t, f_04, 
    "f_04(x) = exp ( x ) - 1 / ( 100 * x * x )" );

  a = -5.0;
  b =  2.0;

  test_zero_rc_one ( a, b, t, f_05, 
    "f_05(x) = (x+3) * (x-1) * (x-1)" );

  return;
}
//****************************************************************************80

void test_local_min_all ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LOCAL_MIN_ALL tests Brent"s local minimizer on all test functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double t;
  double z;

  cout << "\n";
  cout << "TEST_LOCAL_MIN_ALL\n";
  cout << "  Test the Brent LOCAL_MIN routine, which seeks\n";
  cout << "  a local minimizer of a function F(X)\n";
  cout << "  in an interval [A,B].\n";

  t = r8_epsilon ( );

  a = 0.0;
  b = 3.141592653589793;

  test_local_min_one ( a, b, t, g_01, 
    "g_01(x) = ( x - 2 ) * ( x - 2 ) + 1" );

  a = 0.0;
  b = 1.0;

  test_local_min_one ( a, b, t, g_02, 
    "g_02(x) = x * x + exp ( - x )" );

  a = -2.0;
  b =  2.0;

  test_local_min_one ( a, b, t, g_03, 
    "g_03(x) = x^4 + 2x^2 + x + 3" );

  a =  0.0001;
  b =  1.0;

  test_local_min_one ( a, b, t, g_04, 
    "g_04(x) = exp ( x ) + 1 / ( 100 x )" );

  a =  0.0002;
  b = 2.0;

  test_local_min_one ( a, b, t, g_05, 
    "g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)" );

  return;
}
//****************************************************************************80

void test_local_min_rc_all ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LOCAL_MIN_RC_ALL tests LOCAL_MIN_RC on all test functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double t;
  double z;

  cout << "\n";
  cout << "TEST_LOCAL_MIN_RC_ALL\n";
  cout << "  Test the reverse communication version of\n";
  cout << "  the Brent LOCAL_MIN routine, which seeks\n";
  cout << "  a local minimizer of a function F(X)\n";
  cout << "  in an interval [A,B].\n";

  t = 10.0 * sqrt ( r8_epsilon ( ) );

  a = 0.0;
  b = 3.141592653589793;

  test_local_min_rc_one ( a, b, t, g_01, 
    "g_01(x) = ( x - 2 ) * ( x - 2 ) + 1" );

  a = 0.0;
  b = 1.0;

  test_local_min_rc_one ( a, b, t, g_02, 
    "g_02(x) = x * x + exp ( - x )" );

  a = -2.0;
  b =  2.0;

  test_local_min_rc_one ( a, b, t, g_03, 
    "g_03(x) = x^4 + 2x^2 + x + 3" );

  a =  0.0001;
  b =  1.0;

  test_local_min_rc_one ( a, b, t, g_04, 
    "g_04(x) = exp ( x ) + 1 / ( 100 x )" );

  a =  0.0002;
  b = 2.0;

  test_local_min_rc_one ( a, b, t, g_05, 
    "g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)" );

  return;
}
//****************************************************************************80

void test_glomin_all ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GLOMIN_ALL tests the Brent global minimizer on all test functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double e;
  double m;
  double t;
  double z;

  cout << "\n";
  cout << "TEST_GLOMIN_ALL\n";
  cout << "  Test the Brent GLOMIN routine, which seeks\n";
  cout << "  a global minimizer of a function F(X)\n";
  cout << "  in an interval [A,B],\n";
  cout << "  given some upper bound M \n";
  cout << "  for the second derivative of F.\n";

  e = sqrt ( r8_epsilon ( ) );
  t = sqrt ( r8_epsilon ( ) );

  a = 7.0;
  b = 9.0;
  c = ( a + b ) / 2.0;
  m = 0.0;

  test_glomin_one ( a, b, c, m, e, t, h_01, 
    "h_01(x) = 2 - x, M = 0" );

  a = 7.0;
  b = 9.0;
  c = ( a + b ) / 2.0;
  m = 100.0;

  test_glomin_one ( a, b, c, m, e, t, h_01, 
    "h_01(x) = 2 - x, M = 100" );

  a = -1.0;
  b = +2.0;
  c = ( a + b ) / 2.0;
  m = 2.0;

  test_glomin_one ( a, b, c, m, e, t, h_02, 
    "h_02(x) = x * x, M = 2" );

  a = -1.0;
  b = +2.0;
  c = ( a + b ) / 2.0;
  m = 2.1;

  test_glomin_one ( a, b, c, m, e, t, h_02, 
    "h_02(x) = x * x, M = 2.1" );

  a = -0.5;
  b =  +2.0;
  c = ( a + b ) / 2.0;
  m = 14.0;

  test_glomin_one ( a, b, c, m, e, t, h_03, 
    "h_03(x) = x^3 + x^2, M = 14" );

  a = -0.5;
  b =  +2.0;
  c = ( a + b ) / 2.0;
  m = 28.0;

  test_glomin_one ( a, b, c, m, e, t, h_03, 
    "h_03(x) = x^3 + x^2, M = 28" );

  a = -10.0;
  b = +10.0;
  c = ( a + b ) / 2.0;
  m = 72.0;

  test_glomin_one ( a, b, c, m, e, t, h_04, 
    "h_04(x) = ( x + sin(x) ) * exp(-x*x), M = 72" );

  a = -10.0;
  b = +10.0;
  c = ( a + b ) / 2.0;
  m = 72.0;

  test_glomin_one ( a, b, c, m, e, t, h_05, 
    "h_05(x) = ( x - sin(x) ) * exp(-x*x), M = 72" );

  return;
}
//****************************************************************************80

void test_zero_one ( double a, double b, double t, double f ( double x ), 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ZERO_ONE tests Brent's zero finding routine on one test function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the two endpoints of the change of sign
//    interval.
//
//    Input, double T, a positive error tolerance.
//
//    Input, double F ( double x ), the name of a user-supplied
//    function which evaluates the function whose zero is being sought.
//
//    Input, string TITLE, a title for the problem.
//
{
  double fa;
  double fb;
  double fz;
  double z;
  
  z = zero ( a, b, t, f );
  fz = f ( z );
  fa = f ( a );
  fb = f ( b );

  cout << "\n";
  cout << "  " << title << "\n";
  cout << "\n";
  cout << "           A                 Z             B\n";
  cout << "         F(A)              F(Z)          F(B)\n";
  cout << "  " << setw(14) << a
       << "  " << setw(14) << z
       << "  " << setw(14) << b << "\n";
  cout << "  " << setw(14) << fa
       << "  " << setw(14) << fz
       << "  " << setw(14) << fb << "\n";

  return;
}
//****************************************************************************80

void test_zero_rc_one ( double a, double b, double t, double f ( double x ), 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ZERO_RC_ONE tests ZERO_RC on one test function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the two endpoints of the change of sign
//    interval.
//
//    Input, double MACHEP, an estimate for the relative machine
//    precision.
//
//    Input, double T, a positive error tolerance.
//
//    Input, double F ( double x ), the name of a user-supplied
//    function which evaluates the function whose zero is being sought.
//
//    Input, string TITLE, a title for the problem.
//
{
  double arg;
  int status;
  double value;
  
  cout << "\n";
  cout << "  " << title << "\n";
  cout << "\n";
  cout << "    STATUS      X               F(X)\n";
  cout << "\n";

  status = 0;

  for ( ; ; ) 
  {
    zero_rc ( a, b, t, arg, status, value );

    if ( status < 0 )
    {
      cout << "\n";
      cout << "  ZERO_RC returned an error flag!\n";
      break;
    }

    value = f ( arg );

    cout << "  " << setw(8) << status
         << "  " << setw(14) << arg
         << "  " << setw(14) << value << "\n";

    if ( status == 0 )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void test_local_min_one ( double a, double b, double t, double f ( double x ), 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LOCAL_MIN_ONE tests Brent's local minimizer on one test function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, double T, a positive absolute error tolerance.
//
//    Input, double F ( double x ), the name of a user-supplied
//    function, whose local minimum is being sought.
//
//    Input, string TITLE, a title for the problem.
//
{
  double fa;
  double fb;
  double fx;
  double x;

  fx = local_min ( a, b, t, f, x );
  fa = f ( a );
  fb = f ( b );

  cout << "\n";
  cout << "  " << title << "\n";
  cout << "\n";
  cout << "           A                 X             B\n";
  cout << "         F(A)              F(X)          F(B)\n";
  cout << "  " << setw(14) << a
       << "  " << setw(14) << x
       << "  " << setw(14) << b << "\n";
  cout << "  " << setw(14) << fa
       << "  " << setw(14) << fx
       << "  " << setw(14) << fb << "\n";

  return;
}
//****************************************************************************80

void test_local_min_rc_one ( double a, double b, double t, 
  double f ( double x ), string title )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LOCAL_MIN_RC_ONE tests LOCAL_MIN_RC on one test function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, double T, a positive absolute error tolerance.
//
//    Input, double F ( double x ), the name of a user-supplied
//    function, whose local minimum is being sought.
//
//    Input, string TITLE, a title for the problem.
//
{
  double a2;
  double arg;
  double b2;
  int status;
  int step;
  double value;
  double x;

  cout << "\n";
  cout << "  " << title << "\n";
  cout << "\n";
  cout << "  Step      X                          F(X)\n";
  cout << "\n";
  step = 0;

  arg = a;
  value = f ( arg );
  cout << "  " << setw(4) << step
       << "  " << setprecision(16) << setw(24) << arg
       << "  " << setprecision(16) << setw(24) << value << "\n";

  arg = b;
  value = f ( arg );
  cout << "  " << setw(4) << step
       << "  " << setprecision(16) << setw(24) << arg
       << "  " << setprecision(16) << setw(24) << value << "\n";

  a2 = a;
  b2 = b;
  status = 0;

  for ( ; ; )
  {
    arg = local_min_rc ( a2, b2, status, value );
 
    if ( status < 0 )
    {
      cout << "\n";
      cout << "TEST_LOCAL_MIN_RC_ONE - Fatal error!\n";
      cout << "  LOCAL_MIN_RC returned negative status.\n";
      break;
    }

    value = f ( arg );

    step = step + 1;
    cout << "  " << setw(4) << step
         << "  " << setprecision(16) << setw(24) << arg
         << "  " << setprecision(16) << setw(24) << value << "\n";

    if ( 50 < step )
    {
      cout << "\n";
      cout << "TEST_LOCAL_MIN_RC_ONE - Fatal error!\n";
      cout << "  Too many steps!\n";
      break;
     }
    if ( status == 0 )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void test_glomin_one ( double a, double b, double c, double m, 
  double e, double t, double f ( double x ), string title )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GLOMIN_ONE tests the Brent global minimizer on one test function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, double C, an initial guess for the global
//    minimizer.  If no good guess is known, C = A or B is acceptable.
//
//    Input, double M, the bound on the second derivative.
//
//    Input, double E, a positive tolerance, a bound for the
//    absolute error in the evaluation of F(X) for any X in [A,B].
//
//    Input, double T, a positive absolute error tolerance.
//
//    Input, double F ( double x ), the name of a user-supplied
//    function whose global minimum is being sought.
//
//    Input, string TITLE, a title for the problem.
//
{
  double fa;
  double fb;
  double fx;
  double x;

  fx = glomin ( a, b, c, m, e, t, f, x );
  fa = f ( a );
  fb = f ( b );

  cout << "\n";
  cout << "  " << title << "\n";
  cout << "\n";
  cout << "           A                 X             B\n";
  cout << "         F(A)              F(X)          F(B)\n";
  cout << "  " << setprecision(6) << setw(14) << a
       << "  " << setw(14) << x
       << "  " << setw(14) << b << "\n";
  cout << "  " << setw(14) << fa
       << "  " << setw(14) << fx
       << "  " << setw(14) << fb << "\n";

  return;
}
//****************************************************************************80

double f_01 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F_01 evaluates sin ( x ) - x / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double F_01, the value of the function at X.
//
{
  double value;

  value = sin ( x ) - 0.5 * x;

  return value;
}
//****************************************************************************80

double f_02 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F_02 evaluates 2*x-exp(-x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double F_02, the value of the function at X.
//
{
  double value;

  value = 2.0 * x - exp ( - x );

  return value;
}
//****************************************************************************80

double f_03 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F_03 evaluates x*exp(-x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double F_03, the value of the function at X.
//
{
  double value;

  value = x * exp ( - x );

  return value;
}
//****************************************************************************80

double f_04 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F_04 evaluates exp(x) - 1 / (100*x*x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double F_04, the value of the function at X.
//
{
  double value;

  value = exp ( x ) - 1.0 / 100.0 / x / x;

  return value;
}
//****************************************************************************80

double f_05 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F_05 evaluates (x+3)*(x-1)*(x-1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double F_05, the value of the function at X.
//
{
  double value;

  value = ( x + 3.0 ) * ( x - 1.0 ) * ( x - 1.0 );

  return value;
}
//****************************************************************************80

double g_01 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    G_01 evaluates (x-2)^2 + 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double G_01, the value of the function at X.
//
{
  double value;

  value = ( x - 2.0 ) * ( x - 2.0 ) + 1.0;

  return value;
}
//****************************************************************************80

double g_02 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    G_02 evaluates x^2 + exp ( - x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double G_02, the value of the function at X.
//
{
  double value;

  value = x * x + exp ( - x );

  return value;
}
//****************************************************************************80

double g_03 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    G_03 evaluates x^4+2x^2+x+3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double G_03, the value of the function at X.
//
{
  double value;

  value = ( ( x * x + 2.0 ) * x + 1.0 ) * x + 3.0;

  return value;
}
//****************************************************************************80

double g_04 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    G_04 evaluates exp(x)+1/(100X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double G_04, the value of the function at X.
//
{
  double value;

  value = exp ( x ) + 0.01 / x;

  return value;
}
//****************************************************************************80

double g_05 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    G_05 evaluates exp(x) - 2x + 1/(100x) - 1/(1000000x^2)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double G_05, the value of the function at X.
//
{
  double value;

  value = exp ( x ) - 2.0 * x + 0.01 / x - 0.000001 / x / x;

  return value;
}
//****************************************************************************80

double h_01 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    H_01 evaluates 2 - x.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double H_01, the value of the function at X.
//
{
  double value;

  value = 2.0 - x;

  return value;
}
//****************************************************************************80

double h_02 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    H_02 evaluates x^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double H_02, the value of the function at X.
//
{
  double value;

  value = x * x;

  return value;
}
//****************************************************************************80

double h_03 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    H_03 evaluates x^3+x^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double H_03, the value of the function at X.
//
{
  double value;

  value = x * x * ( x + 1.0 );

  return value;
}
//****************************************************************************80

double h_04 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    H_04 evaluates ( x + sin ( x ) ) * exp ( - x * x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double H_04, the value of the function at X.
//
{
  double value;

  value = ( x + sin ( x ) ) * exp ( - x * x );

  return value;
}
//****************************************************************************80

double h_05 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    H_05 evaluates ( x - sin ( x ) ) * exp ( - x * x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which F is to be evaluated.
//
//    Output, double H_05, the value of the function at X.
//
{
  double value;

  value = ( x - sin ( x ) ) * exp ( - x * x );

  return value;
}
