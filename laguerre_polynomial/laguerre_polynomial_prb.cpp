# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "laguerre_polynomial.hpp"

int main ( );
void laguerre_polynomial_test01 ( );
void laguerre_polynomial_test02 ( );
void laguerre_polynomial_test03 ( );
void laguerre_polynomial_test04 ( );
void laguerre_polynomial_test05 ( );
void laguerre_polynomial_test06 ( );
void laguerre_polynomial_test07 ( int p, double b );
void laguerre_polynomial_test08 ( int p, int e );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LAGUERRE_POLYNOMIAL_PRB.
//
//  Discussion:
//
//    LAGUERRE_POLYNOMIAL_PRB tests the LAGUERRE_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2012
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
  cout << "LAGUERRE_POLYNOMIAL_PRB:\n";
  cout << "  C++ version.\n";
  cout << "  Test the LAGUERRE_POLYNOMIAL library.\n";

  laguerre_polynomial_test01 ( );
  laguerre_polynomial_test02 ( );
  laguerre_polynomial_test03 ( );
  laguerre_polynomial_test04 ( );
  laguerre_polynomial_test05 ( );
  laguerre_polynomial_test06 ( );

  p = 5;
  b = 0.0;
  laguerre_polynomial_test07 ( p, b );

  p = 5;
  b = 1.0;
  laguerre_polynomial_test07 ( p, b );

  p = 5;
  e = 0;
  laguerre_polynomial_test08 ( p, e );

  p = 5;
  e = 1;
  laguerre_polynomial_test08 ( p, e );
//
//  Terminate.
//
  cout << "\n";
  cout << "LAGUERRE_POLYNOMIAL_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void laguerre_polynomial_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_TEST01 tests L_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2012
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
  int n;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "LAGUERRE_POLYNOMIAL_TEST01:\n";
  cout << "  L_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the Laguerre polynomials.\n";
  cout << "  L_POLYNOMIAL evaluates the polynomial.\n";
  cout << "\n";
  cout << "                        Tabulated                 Computed\n";
  cout << "     N        X           L(N,X)                    L(N,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    l_polynomial_values ( n_data, n, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = l_polynomial ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    delete [] fx2_vec;

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setw(8) << e << "\n";
  }
  return;
}
//****************************************************************************80

void laguerre_polynomial_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_TEST02 tests L_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 March 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *c;
  int i;
  int j;

  cout << "\n";
  cout << "LAGUERRE_POLYNOMIAL_TEST02\n";
  cout << "  L_POLYNOMIAL_COEFFICIENTS determines Laguerre polynomial coefficients.\n";

  c = l_polynomial_coefficients ( N );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  L(" << i << ",x) =\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(N+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }

  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void laguerre_polynomial_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_TEST03 tests L_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  double *hz;
  string title;
  char title_char[80];
  double *z;

  cout << "\n";
  cout << "LAGUERRE_POLYNOMIAL_TEST03:\n";
  cout << "  L_POLYNOMIAL_ZEROS computes the zeros of L(n,x)\n";
  cout << "  Check by calling L_POLYNOMIAL there.\n";

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = l_polynomial_zeros ( degree );
    sprintf ( title_char, "  Computed zeros for L(%d,z):", degree );
    title = string ( title_char );
    r8vec_print ( degree, z, title );

    hz = l_polynomial ( degree, degree, z );
    sprintf ( title_char, "  Evaluate L(%d,z):", degree );
    title = string ( title_char );
    r8vec_print ( degree, hz+degree*degree, title );

    delete [] hz;
    delete [] z;
  }
  return;
}
//****************************************************************************80

void laguerre_polynomial_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_TEST04 tests L_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2012
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
  cout << "LAGUERRE_POLYNOMIAL_TEST04:\n";
  cout << "  L_QUADRATURE_RULE computes the quadrature rule\n";
  cout << "  associated with L(n,x)\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  l_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  cout << "\n";
  cout << "  Use the quadrature rule to estimate:\n";
  cout << "\n";
  cout << "    Q = Integral ( 0 <= X < +00 ) X^E exp(-X) dx\n";
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
    q_exact = l_integral ( e );
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

void laguerre_polynomial_test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_TEST05 tests LM_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2012
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
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "LAGUERRE_POLYNOMIAL_TEST05:\n";
  cout << "  LM_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the Laguerre polynomials.\n";
  cout << "  LM_POLYNOMIAL evaluates the polynomial.\n";
  cout << "\n";
  cout << "                              Tabulated                 Computed\n";
  cout << "     N     M        X         Lm(N,M,X)                  Lm(N,M,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lm_polynomial_values ( n_data, n, m, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = lm_polynomial ( 1, n, m, x_vec );
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
  return;
}
//****************************************************************************80

void laguerre_polynomial_test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_TEST06 tests LM_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 March 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *c;
  int i;
  int j;
  int m;

  cout << "\n";
  cout << "LAGUERRE_POLYNOMIAL_TEST06\n";
  cout << "  LM_POLYNOMIAL_COEFFICIENTS determines Laguerre polynomial coefficients.\n";

  for ( m = 0; m <= 4; m++ )
  {
    c = lm_polynomial_coefficients ( N, m );
 
    for ( i = 0; i <= N; i++ )
    {
      cout << "\n";
      cout << "  Lm(" << i << "," << m << ",x) =\n";
      cout << "\n";
      for ( j = i; 0 <= j; j-- )
      {
        if ( c[i+j*(N+1)] == 0.0 )
        {
        }
        else if ( j == 0 )
        {
          cout << setw(14) << c[i+j*(N+1)] << "\n";;
        }
        else if ( j == 1 )
        {
          cout << setw(14) << c[i+j*(N+1)] << " * x\n";
        }
        else
        {
          cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
        }
      }
    }
    delete [] c;
  }

  return;
# undef N
}
//****************************************************************************80

void laguerre_polynomial_test07 ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_TEST07 tests L_EXPONENTIAL_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2012
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
  cout << "LAGUERREE_POLYNOMIAL_TEST07\n";
  cout << "  Compute an exponential product table for L(n,x):\n";
  cout << "\n";
  cout << "  Tij = integral ( 0 <= x < +oo ) exp(b*x) Ln(i,x) Ln(j,x) exp(-x) dx\n";
  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponential argument coefficient B = " << b << "\n";

  table = l_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void laguerre_polynomial_test08 ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_TEST08 tests L_POWER_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2012
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
  cout << "LAGUERRE_POLYNOMIAL_TEST08\n";
  cout << "  Compute a power product table for L(n,x).\n";
  cout << "\n";
  cout << "  Tij = integral ( 0 <= x < +oo ) x^e L(i,x) L(j,x) exp(-x) dx\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponent of X, E = " << e << "\n";

  table = l_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  delete [] table;

  return;
}
