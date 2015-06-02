# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>

using namespace std;

int main ( int argc, char *argv[] );
void erf_values ( int *n_data, double *x, double *fx );
void erfc_values ( int *n_data, double *x, double *fx );
void gamma_values ( int *n_data, double *x, double *fx );
double r8_abs ( double x );
void test_erf ( );
void test_erfc ( );
void test_isnan ( );
void test_lgamma ( );
void test_tgamma ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for G++_INTRINSICS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "G++_INTRINSICS:\n";
  cout << "  Demonstrate the use of some G++ intrinsic functions.\n";

  test_erf ( );
  test_erfc ( );
  test_isnan ( );
  test_lgamma ( );
  test_tgamma ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "G++_INTRINSICS:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test_erf ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ERF tests ERF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST_ERF:\n";
  cout << "  Test ERF, which evaluates the error function ERF(X).\n";
  cout << "\n";
  cout << "      X               ERF(X)          ERF(X)          DIFF\n";
  cout << "                     (tabulated)     (computed)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = erf ( x );

    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2
         << "  " << setw(14) << r8_abs ( fx - fx2 ) << "\n";
  }
  return;
}
//****************************************************************************80

void test_erfc ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ERFC tests ERFC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST_ERFC:\n";
  cout << "  Test ERFC, which evaluates the complementary error function.\n";
  cout << "\n";
  cout << "      X               ERFC(X)         ERFC(X)         DIFF\n";
  cout << "                     (tabulated)     (computed)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erfc_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = erfc ( x );

    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2
         << "  " << setw(14) << r8_abs ( fx - fx2 ) << "\n";
  }
  return;
}
//****************************************************************************80

void test_isnan ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_ISNAN tests ISNAN.
//
//  Discussion:
//
//    ISNAN ( X ) is a logical function which is TRUE if the number X
//    is a "NaN", or "Not a Number".  This is a special value set to
//    indicate that an illegal argument was input to a function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 4;
  float x[4] = { - 1.0, 0.0, + 1.0, + 2.0 };

  cout << "\n";
  cout << "TEST_ISNAN\n";
  cout << "  ISNAN(X) is TRUE if X is \"Not a Number\".\n";
  cout << "\n";
  cout << "   Function  -1  0 +1 +2\n";
  cout << "\n";
  
  cout << "  ACOS(X)   ";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << isnan ( acos ( x[i] ) );
  }
  cout << "\n";

  cout << "  ASIN(X)   ";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << isnan ( asin ( x[i] ) );
  }
  cout << "\n";

  cout << "  ATAN(X)   ";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << isnan ( atan ( x[i] ) );
  }
  cout << "\n";

  cout << "   LOG(X)   ";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << isnan ( log ( x[i] ) );
  }
  cout << "\n";

  cout << "  SQRT(X)   ";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << isnan ( sqrt ( x[i] ) );
  }
  cout << "\n";

  cout << "    1 / X   ";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << isnan ( 1.0 / x[i] );
  }
  cout << "\n";

  return;
}
//****************************************************************************80

void test_lgamma ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_LGAMMA tests LGAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << " TEST_LGAMMA:\n";
  cout << "   Test LGAMMA, which evaluates the log(Gamma) function.\n";
  cout << "\n";
  cout << "      X            Log(GAMMA(X))      LGAMMA(X)        DIFF\n";
  cout << "                    (tabulated)      (computed)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    if ( x <= 0.0 )
    {
      continue;
    }

    fx = log ( fx );

    fx2 = lgamma ( x );

    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2
         << "  " << setw(14) << r8_abs ( fx - fx2 ) << "\n";
  }
  return;
}
//****************************************************************************80

void test_tgamma ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_TGAMMA tests TGAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << " TEST_TGAMMA:\n";
  cout << "   Test TGAMMA, which evaluates the Gamma function.\n";
  cout << "\n";
  cout << "      X              GAMMA(X)         TGAMMA(X)        DIFF\n";
  cout << "                    (tabulated)      (computed)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = tgamma ( x );

    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2
         << "  " << setw(14) << r8_abs ( fx - fx2 ) << "\n";
  }
  return;
}
//****************************************************************************80

void erf_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    ERF_VALUES returns some values of the ERF or "error" function.
//
//  Discussion:
//
//    The error function is defined by:
//
//      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
//
//    In Mathematica, the function can be evaluated by:
//
//      Erf[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  double fx_vec[N_MAX] = { 
     0.0000000000000000E+00,  
     0.1124629160182849E+00,  
     0.2227025892104785E+00,  
     0.3286267594591274E+00,  
     0.4283923550466685E+00,  
     0.5204998778130465E+00,  
     0.6038560908479259E+00,  
     0.6778011938374185E+00,  
     0.7421009647076605E+00,  
     0.7969082124228321E+00,  
     0.8427007929497149E+00,  
     0.8802050695740817E+00,  
     0.9103139782296354E+00,  
     0.9340079449406524E+00,  
     0.9522851197626488E+00,  
     0.9661051464753107E+00,  
     0.9763483833446440E+00,  
     0.9837904585907746E+00,  
     0.9890905016357307E+00,  
     0.9927904292352575E+00,  
     0.9953222650189527E+00 }; 

  double x_vec[N_MAX] = { 
     0.0E+00,   
     0.1E+00,   
     0.2E+00,   
     0.3E+00,   
     0.4E+00,   
     0.5E+00,   
     0.6E+00,   
     0.7E+00,   
     0.8E+00,   
     0.9E+00,   
     1.0E+00,   
     1.1E+00,   
     1.2E+00,   
     1.3E+00,   
     1.4E+00,   
     1.5E+00,   
     1.6E+00,   
     1.7E+00,   
     1.8E+00,   
     1.9E+00,   
     2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void erfc_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    ERFC_VALUES returns some values of the ERFC or "complementary error" function.
//
//  Discussion:
//
//    The complementary error function is defined by:
//
//      ERFC(X) = 1 - ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
//
//    In Mathematica, the function can be evaluated by:
//
//      Erfc[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  double fx_vec[N_MAX] = { 
    1.000000000000000E+00, 
    0.7772974107895215E+00, 
    0.5716076449533315E+00, 
    0.3961439091520741E+00, 
    0.2578990352923395E+00, 
    0.1572992070502851E+00, 
    0.08968602177036462E+00, 
    0.04771488023735119E+00, 
    0.02365161665535599E+00, 
    0.01090949836426929E+00, 
    0.004677734981047266E+00,
    0.001862846297981891E+00, 
    0.0006885138966450786E+00, 
    0.0002360344165293492E+00, 
    0.00007501319466545902E+00, 
    0.00002209049699858544E+00, 
    6.025761151762095E-06, 
    1.521993362862285E-06, 
    3.558629930076853E-07, 
    7.700392745696413E-08, 
    1.541725790028002E-08 }; 

  double x_vec[N_MAX] = { 
    0.0E+00,  
    0.2E+00, 
    0.4E+00,  
    0.6E+00, 
    0.8E+00,  
    1.0E+00,  
    1.2E+00, 
    1.4E+00, 
    1.6E+00, 
    1.8E+00, 
    2.0E+00,
    2.2E+00, 
    2.4E+00, 
    2.6E+00, 
    2.8E+00, 
    3.0E+00, 
    3.2E+00, 
    3.4E+00, 
    3.6E+00, 
    3.8E+00, 
    4.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_VALUES returns some values of the Gamma function.
//
//  Discussion:
//
//    The Gamma function is defined as:
//
//      Gamma(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) exp(-T) dT
//
//    It satisfies the recursion:
//
//      Gamma(X+1) = X * Gamma(X)
//
//    Gamma is undefined for nonpositive integral X.
//    Gamma(0.5) = sqrt(PI)
//    For N a positive integer, Gamma(N+1) = N!, the standard factorial.
//
//    In Mathematica, the function can be evaluated by:
//
//      Gamma[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 25

  double fx_vec[N_MAX] = { 
     -0.3544907701811032E+01,  
     -0.1005871979644108E+03,  
      0.9943258511915060E+02,  
      0.9513507698668732E+01,  
      0.4590843711998803E+01,  
      0.2218159543757688E+01,  
      0.1772453850905516E+01,  
      0.1489192248812817E+01,  
      0.1164229713725303E+01,  
      0.1000000000000000E+01,  
      0.9513507698668732E+00,  
      0.9181687423997606E+00,  
      0.8974706963062772E+00,  
      0.8872638175030753E+00,  
      0.8862269254527580E+00,  
      0.8935153492876903E+00,  
      0.9086387328532904E+00,  
      0.9313837709802427E+00,  
      0.9617658319073874E+00,  
      0.1000000000000000E+01,  
      0.2000000000000000E+01,  
      0.6000000000000000E+01,  
      0.3628800000000000E+06,  
      0.1216451004088320E+18,  
      0.8841761993739702E+31 };

  double x_vec[N_MAX] = { 
     -0.50E+00,  
     -0.01E+00,  
      0.01E+00,  
      0.10E+00,  
      0.20E+00,  
      0.40E+00,  
      0.50E+00,  
      0.60E+00,  
      0.80E+00,  
      1.00E+00,  
      1.10E+00,  
      1.20E+00,  
      1.30E+00,  
      1.40E+00,  
      1.50E+00,  
      1.60E+00,  
      1.70E+00,  
      1.80E+00,  
      1.90E+00,  
      2.00E+00,  
      3.00E+00,  
      4.00E+00,  
     10.00E+00,  
     20.00E+00,  
     30.00E+00 }; 

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
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
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

void timestamp ( void )

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
