# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "legendre_polynomial.hpp"

int main ( );
void legendre_polynomial_test01 ( );
void legendre_polynomial_test015 ( );
void legendre_polynomial_test016 ( );
void legendre_polynomial_test02 ( );
void legendre_polynomial_test03 ( );
void legendre_polynomial_test04 ( );
void legendre_polynomial_test05 ( int p, double b );
void legendre_polynomial_test06 ( int p, int e );
void legendre_polynomial_test07 ( );
void legendre_polynomial_test08 ( );
void legendre_polynomial_test09 ( );
void legendre_polynomial_test095 ( );
void legendre_polynomial_test10 ( int p );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LEGENDRE_POLYNOMIAL_PRB.
//
//  Discussion:
//
//    LEGENDRE_POLYNOMIAL_PRB tests the LEGENDRE_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int e;
  int p;

  timestamp ( );
  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_PRB:\n";
  cout << "  C++ version.\n";
  cout << "  Test the LEGENDRE_POLYNOMIAL library.\n";

  legendre_polynomial_test01 ( );
  legendre_polynomial_test015 ( );
  legendre_polynomial_test016 ( );
  legendre_polynomial_test02 ( );
  legendre_polynomial_test03 ( );
  legendre_polynomial_test04 ( );

  p = 5;
  b = 0.0;
  legendre_polynomial_test05 ( p, b );

  p = 5;
  b = 1.0;
  legendre_polynomial_test05 ( p, b );

  p = 5;
  e = 0;
  legendre_polynomial_test06 ( p, e );

  p = 5;
  e = 1;
  legendre_polynomial_test06 ( p, e );

  legendre_polynomial_test07 ( );
  legendre_polynomial_test08 ( );
  legendre_polynomial_test09 ( );
  legendre_polynomial_test095 ( );

  p = 5;
  legendre_polynomial_test10 ( p );
//
//  Terminate.
//
  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void legendre_polynomial_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST01 tests P_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m = 1;
  int n;
  int prec;
  double *v;
  double x;
  double x_vec[1];

  prec = cout.precision ( );

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST01:\n";
  cout << "  P_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the Legendre polynomial P(n,x).\n";
  cout << "  P_POLYNOMIAL_VALUE evaluates the polynomial.\n";
  cout << "\n";
  cout << "                        Tabulated                 Computed\n";
  cout << "     N        X           P(N,X)                    P(N,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    p_polynomial_values ( n_data, n, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = p_polynomial_value ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    delete [] fx2_vec;

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setw(8) << e << "\n";
  }
//
//  Restore default precision.
//
  cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_polynomial_test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST015 tests P_POLYNOMIAL_PRIME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int j;
  int m;
  int n;
  int prec;
  double *vp;
  double *x;

  prec = cout.precision ( );

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST015:\n";
  cout << "  P_POLYNOMIAL_PRIME evaluates the derivative of the\n";
  cout << "  Legendre polynomial P(n,x).\n";
  cout << "\n";
  cout << "                        Computed\n";
  cout << "     N        X           P'(N,X)\n";
  cout << "\n";

  m = 11;
  a = - 1.0;
  b = + 1.0;
  x = r8vec_linspace_new ( m, a, b );

  n = 5;
  vp = p_polynomial_prime ( m, n, x );

  for ( i = 0; i < m; i++ )
  {
    cout << "\n";
    for ( j = 0; j <= n; j++ )
    {
      cout << "  " << setw(4) << j
           << "  " << setw(12) << x[i]
           << "  " << setprecision(16) << setw(24) << vp[i+j*m] << "\n";
    }
  }

  delete [] vp;
  delete [] x;
//
//  Restore default precision.
//
  cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_polynomial_test016 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST016 tests P_POLYNOMIAL_PRIME2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int j;
  int m;
  int n;
  int prec;
  double *vpp;
  double *x;

  prec = cout.precision ( );

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST016:\n";
  cout << "  P_POLYNOMIAL_PRIME2 evaluates the second derivative of the\n";
  cout << "  Legendre polynomial P(n,x).\n";
  cout << "\n";
  cout << "                        Computed\n";
  cout << "     N        X           P\"(N,X)\n";
  cout << "\n";

  m = 11;
  a = - 1.0;
  b = + 1.0;
  x = r8vec_linspace_new ( m, a, b );

  n = 5;
  vpp = p_polynomial_prime2 ( m, n, x );

  for ( i = 0; i < m; i++ )
  {
    cout << "\n";
    for ( j = 0; j <= n; j++ )
    {
      cout << "  " << setw(4) << j
           << "  " << setw(12) << x[i]
           << "  " << setprecision(16) << setw(24) << vpp[i+j*m] << "\n";
    }
  }

  delete [] vpp;
  delete [] x;
//
//  Restore default precision.
//
  cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_polynomial_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST02 tests P_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int j;
  int n = 10;

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST02\n";
  cout << "  P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).\n";

  c = p_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    cout << "\n";
    cout << "  P(" << i << ",x) =\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(n+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(n+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(n+1)] << " * x^" << j << "\n";
      }
    }
  }
  delete [] c;

  return;
}
//****************************************************************************80

void legendre_polynomial_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST03 tests P_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  double *hz;
  string title;
  double *z;

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST03:\n";
  cout << "  P_POLYNOMIAL_ZEROS computes the zeros of P(n,x)\n";
  cout << "  Check by calling P_POLYNOMIAL_VALUE there.\n";

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = p_polynomial_zeros ( degree );
    title = "  Computed zeros for P(" + i4_to_s ( degree ) + ",z):";
    r8vec_print ( degree, z, title );

    hz = p_polynomial_value ( degree, degree, z );
    title = "  Evaluate P(" + i4_to_s ( degree ) + ",z):";
    r8vec_print ( degree, hz+degree*degree, title );

    delete [] hz;
    delete [] z;
  }
  return;
}
//****************************************************************************80

void legendre_polynomial_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST04 tests P_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double *f;
  int i;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST04:\n";
  cout << "  P_QUADRATURE_RULE computes the quadrature rule\n";
  cout << "  associated with P(n,x)\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  p_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  cout << "\n";
  cout << "  Use the quadrature rule to estimate:\n";
  cout << "\n";
  cout << "    Q = Integral ( -1 <= X <= +1.0 ) X^E dx\n";
  cout << "\n";
  cout << "   E       Q_Estimate      Q_Exact\n";
  cout << "\n";

  f = new double[n];

  for ( e = 0; e <= 2 * n - 1; e++ )
  {
    if ( e == 0 )
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = pow ( x[i], e );
      }
    }
    q = r8vec_dot_product ( n, w, f );
    q_exact = p_integral ( e );
    cout << "  " << setw(2) << e
         << "  " << setw(14) << q
         << "  " << setw(14) << q_exact << "\n";
  }

  delete [] f;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void legendre_polynomial_test05 ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST05 tests P_EXPONENTIAL_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polynomial 
//    factors.
//
//    Input, double B, the coefficient of X in the exponential factor.
//
{
  double *table;

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST05\n";
  cout << "  Compute a Legendre exponential product table.\n";
  cout << "\n";
  cout << "  Tij = integral ( -1.0 <= X <= +1.0 ) exp(B*X) P(I,X) P(J,X) dx\n";
  cout << "\n";
  cout << "  where P(I,X) = Legendre polynomial of degree I.\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponential argument coefficient B = " << b << "\n";

  table = p_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void legendre_polynomial_test06 ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST06 tests P_POWER_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polynomial 
//    factors.
//
//    Input, int E, the exponent of X.
//
{
  double *table;

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST06\n";
  cout << "  Compute a Legendre power product table.\n";
  cout << "\n";
  cout << "  Tij = integral ( -1.0 <= X <= +1.0 ) X^E P(I,X) P(J,X) dx\n";
  cout << "\n";
  cout << "  where P(I,X) = Legendre polynomial of degree I.\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponent of X, E = " << e << "\n";

  table = p_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void legendre_polynomial_test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST07 tests PM_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  int prec;
  double *v;
  double x;
  double x_vec[1];

  prec = cout.precision ( );

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST07:\n";
  cout << "  PM_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the Legendre polynomial Pm(n,m,x).\n";
  cout << "  PM_POLYNOMIAL_VALUE evaluates the polynomial.\n";
  cout << "\n";
  cout << "                             Tabulated                 Computed\n";
  cout << "     N     M        X        Pm(N,M,X)                 Pm(N,M,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    pm_polynomial_values ( n_data, n, m, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = pm_polynomial_value ( 1, n, m, x_vec );
    fx2 = fx2_vec[n];
    delete [] fx2_vec;

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(4) << m
         << "  " << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setw(8) << e << "\n";
  }
//
//  Restore default precision.
//
  cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_polynomial_test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST08 tests PMN_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  int prec;
  double *v;
  double x;
  double x_vec[1];

  prec = cout.precision ( );

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST08:\n";
  cout << "  PMN_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the normalized Legendre polynomial Pmn(n,m,x).\n";
  cout << "  PMN_POLYNOMIAL_VALUE evaluates the polynomial.\n";
  cout << "\n";
  cout << "                             Tabulated                 Computed\n";
  cout << "     N     M        X       Pmn(N,M,X)                Pmn(N,M,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    pmn_polynomial_values ( n_data, n, m, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = pmn_polynomial_value ( 1, n, m, x_vec );
    fx2 = fx2_vec[n];
    delete [] fx2_vec;

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(4) << m
         << "  " << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setw(8) << e << "\n";
  }
//
//  Restore default precision.
//
  cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_polynomial_test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST09 tests PMNS_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  int prec;
  double *v;
  double x;
  double x_vec[1];

  prec = cout.precision ( );

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST09:\n";
  cout << "  PMNS_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the sphere-normalized Legendre polynomial Pmns(n,m,x).\n";
  cout << "  PMNS_POLYNOMIAL_VALUE evaluates the polynomial.\n";
  cout << "\n";
  cout << "                             Tabulated                 Computed\n";
  cout << "     N     M        X       Pmns(N,M,X)                Pmns(N,M,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    pmns_polynomial_values ( n_data, n, m, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = pmns_polynomial_value ( 1, n, m, x_vec );
    fx2 = fx2_vec[n];
    delete [] fx2_vec;

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(4) << m
         << "  " << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setw(8) << e << "\n";
  }
//
//  Restore default precision.
//
  cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_polynomial_test095 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST095 tests PN_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int j;
  int n = 10;

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST095\n";
  cout << "  PN_POLYNOMIAL_COEFFICIENTS: coefficients of normalized Legendre polynomial P(n,x).\n";

  c = pn_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    cout << "\n";
    cout << "  P(" << i << ",x) =\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(n+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(n+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(n+1)] << " * x^" << j << "\n";
      }
    }
  }
  delete [] c;

  return;
}
//****************************************************************************80

void legendre_polynomial_test10 ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLYNOMIAL_TEST10 tests PN_PAIR_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polynomial 
//    factors.
//
{
  double *table;

  cout << "\n";
  cout << "LEGENDRE_POLYNOMIAL_TEST10\n";
  cout << "  Compute a pair product table for Pn(n,x).\n";
  cout << "\n";
  cout << "  Tij = integral ( -1.0 <= X <= +1.0 ) Pn(I,X) Pn(J,X) dx\n";
  cout << "\n";
  cout << "  where Pn(I,X) = normalized Legendre polynomial of degree I.\n";
  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";

  table = pn_pair_product ( p );

  r8mat_print ( p + 1, p + 1, table, "  Pair product table:" );

  delete [] table;

  return;
}
