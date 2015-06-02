# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "filon.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
double *zero_integrand ( int n, double x[] );
double *one_integrand ( int n, double x[] );
double *two_integrand ( int n, double x[] );
double *log_integrand ( int n, double x[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FILON_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FILON_PRB\n";
  cout << "  C version\n";
  cout << "  Test the FILON library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FILON_PRB\n";
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
//    TEST01 tests FILON_TAB_COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  double exact;
  double *ftab;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
  double *x;

  a = 0.0;
  b = 2.0 * r8_pi;

  n = 11;
//
//  Set the X values.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  cout << "\n";
  cout << "TEST01\n";
  cout << "  FILON_TAB_COS estimates the integral of.\n";
  cout << "  F(X) * COS ( T * X )\n";
  cout << "  Use integrands F(X)=1, X, X^2.\n";
  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "       T                      Approximate             Exact\n";

  for ( k = 1; k <= 3; k++ )
  {
    if ( k == 1 )
    {
      t = 1.0;
    }
    else if ( k == 2 )
    {
      t = 2.0;
    }
    else if ( k == 3 )
    {
      t = 10.0;
    }

    cout << "\n";

    for ( i = 1; i <= 3; i++ )
    {
      if ( i == 1 )
      {
        ftab = zero_integrand ( n, x );
      }
      else if ( i == 2 )
      {
        ftab = one_integrand ( n, x );
      }
      else if ( i == 3 )
      {
        ftab = two_integrand ( n, x );
      }

      result = filon_tab_cos ( n, ftab, a, b, t );

      if ( i == 1 )
      {
        exact = ( sin ( t * b ) - sin ( t * a ) ) / t;
      }
      else if ( i == 2 )
      {
        exact = ( ( cos ( t * b ) + t * b * sin ( t * b ) ) 
                - ( cos ( t * a ) + t * a * sin ( t * a ) ) ) / t / t;
      }
      else if ( i == 3 )
      {
        exact = ( ( 2.0 * t * b * cos ( t * b ) 
              + ( t * t * b * b - 2.0 ) * sin ( t * b ) ) 
                - ( 2.0 * t * a * cos ( t * a ) 
              + ( t * t * a * a - 2.0 ) * sin ( t * a ) ) ) / t / t / t;
      }
      cout << setw(24) << t << "  "
           << setw(24) << result << "  "
           << setw(24) << exact << "\n";
    }
  }

  delete [] ftab;
  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests FILON_TAB_COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  double exact;
  double *ftab;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
  double *x;
//
//  Example suggested by James Roedder.
//
  cout << "\n";
  cout << "TEST02\n";
  cout << "  Integrate F(X) = log(1+X)*cos(T*X):\n";
  cout << "  Supply integrand as a table.\n";
  cout << "  T = 10, and N increases\n";
  cout << "\n";
  cout << "       N    Approximate             Exact                   Error\n";
  cout << "\n";

  a = 0.0;
  b = 2.0 * r8_pi;

  for ( j = 1; j <= 6; j++ )
  {
    n = i4_power ( 2, j ) * 10 + 1;
//
//  Set the X values.
//
    x = new double[n];
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i - 1 ) * a   
             + ( double ) (     i     ) * b ) 
             / ( double ) ( n     - 1 );
    }

    ftab = log_integrand ( n, x );
 
    t = 10.0;

    result = filon_tab_cos ( n, ftab, a, b, t );

    exact = -0.008446594405;
    error = result - exact;

    cout << setw(6) << n << "  "
         << setw(24) << result << "  "
         << setw(24) << exact << "  "
         << setw(24) << error << "\n";

    delete [] ftab;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests FILON_FUN_COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
//
//  Example suggested by James Roedder.
//
  cout << "\n";
  cout << "TEST03\n";
  cout << "  Integrate F(X)=log(1+X)*cos(T*X):\n";
  cout << "  Supply integrand as a function.\n";
  cout << "  T = 10, and N increases\n";
  cout << "\n";
  cout << "       N    Approximate             Exact                   Error\n";
  cout << "\n";

  a = 0.0;
  b = 2.0 * r8_pi;

  for ( j = 1; j <= 6; j++ )
  {
    n = i4_power ( 2, j ) * 10 + 1;
 
    t = 10.0;

    result = filon_fun_cos ( n, log_integrand, a, b, t );

    exact = -0.008446594405;
    error = result - exact;

    cout << setw(6) << n << "  "
         << setw(24) << result << "  "
         << setw(24) << exact << "  "
         << setw(24) << error << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests FILON_TAB_SIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  double exact;
  double *ftab;
  int i;
  int j;
  int k;
  int n;
  double r8_pi = 3.141592653589793;
  double result;
  double t;
  double *x;

  a = 0.0;
  b = 2.0 * r8_pi;
  n = 11;
//
//  Set the X values.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  cout << "\n";
  cout << "TEST04\n";
  cout << "  FILON_TAB_SIN estimates the integral of.\n";
  cout << "  F(X) * SIN ( T * X )\n";
  cout << "  Use integrands 1, X, X^2.\n";
  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "       T                      Approximate             Exact\n";
  cout << "\n";

  for ( k = 1; k <= 3; k++ )
  {
    if ( k == 1 )
    {
      t = 1.0;
    }
    else if ( k == 2 )
    {
      t = 2.0;
    }
    else if ( k == 3 )
    {
      t = 10.0;
    }

    cout << "\n";

    for ( i = 1; i <= 3; i++ )
    {
      if ( i == 1 )
      {
        ftab = zero_integrand ( n, x );
      }
      else if ( i == 2 )
      {
        ftab = one_integrand ( n, x );
      }
      else if ( i == 3 )
      {
        ftab = two_integrand ( n, x );
      }

      result = filon_tab_sin ( n, ftab, a, b, t );

      if ( i == 1 )
      {
        exact = ( - cos ( t * b ) + cos ( t * a ) ) / t;
      }
      else if ( i == 2 )
      {
        exact = ( ( sin ( t * b ) - t * b * cos ( t * b ) ) 
                - ( sin ( t * a ) - t * a * cos ( t * a ) ) ) / t / t;
      }
      else if ( i == 3 )
      {
        exact = ( ( 2.0 * t * b * sin ( t * b ) 
                + ( 2.0 - t * t * b * b ) * cos ( t * b ) ) 
                - ( 2.0 * t * a * sin ( t * a ) 
                + ( 2.0 - t * t * a * a ) * cos ( t * a ) ) ) / t / t / t;
      }
      cout << setw(24) << t << "  "
           << setw(24) << result << "  "
           << setw(24) << exact << "\n";
    }
  }

  delete [] ftab;
  delete [] x;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests FILON_TAB_COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  double exact;
  double *ftab;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
  double *x;
//
//  Example suggested by James Roedder.
//
  cout << "\n";
  cout << "TEST05\n";
  cout << "  Integrate F(X)=log(1+X)*sin(T*X):\n";
  cout << "  Supply integrand as a table.\n";
  cout << "  T = 10, and N increases\n";
  cout << "\n";
  cout << "       N    Approximate             Exact                   Error\n";
  cout << "\n";

  a = 0.0;
  b = 2.0 * r8_pi;

  for ( j = 1; j <= 6; j++ )
  {
    n = i4_power ( 2, j ) * 10 + 1;
//
//  Set the X values.
//
    x = new double[n];
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i - 1 ) * a   
             + ( double ) (     i     ) * b ) 
             / ( double ) ( n     - 1 );
    }

    ftab = log_integrand ( n, x );

    t = 10.0;

    result = filon_tab_sin ( n, ftab, a, b, t );

    exact = -0.19762680771872;
    error = result - exact;

    cout << setw(6) << n << "  "
         << setw(24) << result << "  "
         << setw(24) << exact << "  "
         << setw(24) << error << "\n";

    delete [] ftab;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests FILON_FUN_COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
//
//  Example suggested by James Roedder.
//
  cout << "\n";
  cout << "TEST06\n";
  cout << "  Integrate F(X)=log(1+X)*sin(T*X):\n";
  cout << "  Supply integrand as a function.\n";
  cout << "  T = 10, and N increases\n";
  cout << "\n";
  cout << "       N    Approximate             Exact                   Error\n";
  cout << "\n";

  a = 0.0;
  b = 2.0 * r8_pi;

  for ( j = 1; j <= 6; j++ )
  {
    n = i4_power ( 2, j ) * 10 + 1;

    t = 10.0;

    result = filon_fun_sin ( n, log_integrand, a, b, t );

    exact = -0.19762680771872;
    error = result - exact;

    cout << setw(6) << n << "  "
         << setw(24) << result << "  "
         << setw(24) << exact << "  "
         << setw(24) << error << "\n";
  }

  return;
}
//****************************************************************************80

double *zero_integrand ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZERO_INTEGRAND evaluates the integrand x^0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double ZERO_INTEGRAND[N], the function values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0;
  }
  return fx;
}
//****************************************************************************80

double *one_integrand ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    ONE_INTEGRAND evaluates the integrand X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double ONE_INTEGRAND[N], the function values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = x[i];
  }
  return fx;
}
//****************************************************************************80

double *two_integrand ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TWO_INTEGRAND evaluates the integrand X^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double TWO_INTEGRAND[N], the function values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = x[i] * x[i];
  }
  return fx;
}
//****************************************************************************80

double *log_integrand ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_INTEGRAND evaluates the logarithmic integrand.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double LOG_INTEGRAND[N], the function values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = log ( 1.0 + x[i] );
  }
  return fx;
}
