# include <cstdlib>
# include <cmath>
# include <iostream>
# include <sstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "hermite_polynomial.hpp"

int main ( );
void hermite_polynomial_test01 ( );
void hermite_polynomial_test02 ( );
void hermite_polynomial_test03 ( );
void hermite_polynomial_test04 ( );
void hermite_polynomial_test05 ( );
void hermite_polynomial_test06 ( );
void hermite_polynomial_test07 ( );
void hermite_polynomial_test08 ( int p, double b );
void hermite_polynomial_test09 ( int p, int e );
void hermite_polynomial_test10 ( int p, double b );
void hermite_polynomial_test11 ( int p, int e );
void hermite_polynomial_test12 ( int p, double b );
void hermite_polynomial_test13 ( int p, int e );
void hermite_polynomial_test14 ( );
void hermite_polynomial_test15 ( );
string i4_to_string ( int i4 );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HERMITE_POLYNOMIAL_PRB.
//
//  Discussion:
//
//    HERMITE_POLYNOMIAL_PRB tests the HERMITE_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2012
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
  cout << "HERMITE_POLYNOMIAL_PRB:\n";
  cout << "  C++ version.\n";
  cout << "  Test the HERMITE_POLYNOMIAL library.\n";

  hermite_polynomial_test01 ( );
  hermite_polynomial_test02 ( );
  hermite_polynomial_test03 ( );
  hermite_polynomial_test04 ( );
  hermite_polynomial_test05 ( );
  hermite_polynomial_test06 ( );
  hermite_polynomial_test07 ( );

  p = 5;
  b = 0.0;
  hermite_polynomial_test08 ( p, b );

  p = 5;
  b = 1.0;
  hermite_polynomial_test08 ( p, b );

  p = 5;
  e = 0;
  hermite_polynomial_test09 ( p, e );

  p = 5;
  e = 1;
  hermite_polynomial_test09 ( p, e );

  p = 5;
  b = 0.0;
  hermite_polynomial_test10 ( p, b );

  p = 5;
  b = 1.0;
  hermite_polynomial_test10 ( p, b );

  p = 5;
  e = 0;
  hermite_polynomial_test11 ( p, e );

  p = 5;
  e = 1;
  hermite_polynomial_test11 ( p, e );

  p = 5;
  b = 0.0;
  hermite_polynomial_test12 ( p, b );

  p = 5;
  b = 1.0;
  hermite_polynomial_test12 ( p, b );

  p = 5;
  e = 0;
  hermite_polynomial_test13 ( p, e );

  p = 5;
  e = 1;
  hermite_polynomial_test13 ( p, e );

  hermite_polynomial_test14 ( );

  hermite_polynomial_test15 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HERMITE_POLYNOMIAL_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void hermite_polynomial_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST01 tests H_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST01:\n";
  cout << "  H_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the physicist's Hermite polynomials.\n";
  cout << "  H_POLYNOMIAL_VALUE evaluates the polynomial.\n";
  cout << "\n";
  cout << "                        Tabulated                 Computed\n";
  cout << "     N        X           H(N,X)                    H(N,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    h_polynomial_values ( n_data, n, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = h_polynomial_value ( 1, n, x_vec );
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

void hermite_polynomial_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST02 tests HE_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST02:\n";
  cout << "  HE_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the probabilist's Hermite polynomials.\n";
  cout << "  HE_POLYNOMIAL_VALUE evaluates the polynomial.\n";
  cout << "\n";
  cout << "                        Tabulated                 Computed\n";
  cout << "     N        X          He(N,X)                   He(N,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    he_polynomial_values ( n_data, n, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = he_polynomial_value ( 1, n, x_vec );
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

void hermite_polynomial_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST03 tests HF_FUNCTION_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST03:\n";
  cout << "  HF_FUNCTION_VALUES stores values of\n";
  cout << "  the Hermite function Hf(n,x).\n";
  cout << "  HF_FUNCTION_VALUE evaluates the function.\n";
  cout << "\n";
  cout << "                        Tabulated                 Computed\n";
  cout << "     N        X          Hf(N,X)                   Hf(N,X)                   Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hf_function_values ( n_data, n, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = hf_function_value ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    delete [] fx2_vec;

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(6)  << setw(14) << e << "\n";
  }
  return;
}
//****************************************************************************80

void hermite_polynomial_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST04 tests H_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2014
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
  cout << "HERMITE_POLYNOMIAL_TEST04:\n";
  cout << "  H_POLYNOMIAL_ZEROS computes the zeros of H(n,x)\n";
  cout << "  Check by calling H_POLYNOMIAL there.\n";

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = h_polynomial_zeros ( degree );
    title = "  Computed zeros for H(" + i4_to_string ( degree ) + ",z):";
    r8vec_print ( degree, z, title );

    hz = h_polynomial_value ( degree, degree, z );
    title = "  Evaluate H(" + i4_to_string ( degree ) + ",z):";
    r8vec_print ( degree, hz+degree*degree, title );

    delete [] hz;
    delete [] z;
  }
  return;
}
//****************************************************************************80

void hermite_polynomial_test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST05 tests HE_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST05:\n";
  cout << "  HE_POLYNOMIAL_ZEROS computes the zeros of He(n,x)\n";
  cout << "  Check by calling HE_POLYNOMIAL there.\n";

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = he_polynomial_zeros ( degree );
    title = "  Computed zeros for He(" + i4_to_string ( degree ) + ",z):";
    r8vec_print ( degree, z, title );

    hz = he_polynomial_value ( degree, degree, z );
    title = "  Evaluate He(" + i4_to_string ( degree ) + ",z):";
    r8vec_print ( degree, hz+degree*degree, title );

    delete [] hz;
    delete [] z;
  }
  return;
}
//****************************************************************************80

void hermite_polynomial_test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST06 tests H_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 March 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST06:\n";
  cout << "  H_QUADRATURE_RULE computes the quadrature rule\n";
  cout << "  associated with H(n,x)\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  h_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  cout << "\n";
  cout << "  Use the quadrature rule to estimate:\n";
  cout << "\n";
  cout << "    Q = Integral ( -oo < X < +00 ) X^E exp(-X^2) dx\n";
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
    q_exact = h_integral ( e );
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

void hermite_polynomial_test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST07 tests HE_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 March 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST07:\n";
  cout << "  HE_QUADRATURE_RULE computes the quadrature rule\n";
  cout << "  associated with He(n,x)\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  he_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  cout << "\n";
  cout << "  Use the quadrature rule to estimate:\n";
  cout << "\n";
  cout << "    Q = Integral ( -oo < X < +00 ) X^E exp(-X^2) dx\n";
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
    q_exact = he_integral ( e );
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

void hermite_polynomial_test08 ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST08 tests HN_EXPONENTIAL_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST08\n";
  cout << "  Compute a normalized physicist''s Hermite exponential product table.\n";
  cout << "\n";
  cout << "  Tij = integral ( -oo < X < +oo ) exp(B*X) Hn(I,X) Hn(J,X) exp(-X*X) dx\n";
  cout << "\n";
  cout << "  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponential argument coefficient B = " << b << "\n";

  table = hn_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void hermite_polynomial_test09 ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST09 tests HN_POWER_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST09\n";
  cout << "  Compute a normalized physicist''s Hermite power product table.\n";
  cout << "\n";
  cout << "  Tij = integral ( -oo < X < +oo ) X^E Hn(I,X) Hn(J,X) exp(-X*X) dx\n";
  cout << "\n";
  cout << "  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponent of X, E = " << e << "\n";

  table = hn_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void hermite_polynomial_test10 ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST10 tests HEN_EXPONENTIAL_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST10\n";
  cout << "  Compute a normalized probabilist''s Hermite exponential product table.\n";
  cout << "\n";
  cout << "  Tij = integral ( -oo < X < +oo ) exp(B*X) Hen(I,X) Hen(J,X) exp(-0.5*X*X) dx\n";
  cout << "\n";
  cout << "  where Hen(I,X) = normalized probabilist''s Hermite polynomial of degree I.\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponential argument coefficient B = " << b << "\n";

  table = hen_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void hermite_polynomial_test11 ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST11 tests HEN_POWER_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST11\n";
  cout << "  Compute a normalized probabilist''s Hermite power product table.\n";
  cout << "\n";
  cout << "  Tij = integral ( -oo < X < +oo ) X^E Hen(I,X) Hen(J,X) exp(-X*X) dx\n";
  cout << "\n";
  cout << "  where Hen(I,X) = normalized probabilist''s Hermite polynomial of degree I.\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponent of X, E = " << e << "\n";

  table = hen_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void hermite_polynomial_test12 ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST12 tests HF_EXPONENTIAL_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST12\n";
  cout << "  Compute a Hermite function exponential product table.\n";
  cout << "\n";
  cout << "  Tij = integral ( -oo < X < +oo ) exp(B*X) Hf(I,X) Hf(J,X) dx\n";
  cout << "\n";
  cout << "  where Hf(I,X) = Hermite function of \"degree\" I.\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponential argument coefficient B = " << b << "\n";

  table = hf_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void hermite_polynomial_test13 ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST13 tests HF_POWER_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
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
  cout << "HERMITE_POLYNOMIAL_TEST13\n";
  cout << "  Compute a Hermite function product table.\n";
  cout << "\n";
  cout << "  Tij = integral ( -oo < X < +oo ) X^E Hf(I,X) Hf(J,X) exp(-X*X) dx\n";
  cout << "\n";
  cout << "  where Hf(I,X) = Hermite function of \"degree\" I.\n";

  cout << "\n";
  cout << "  Maximum degree P = " << p << "\n";
  cout << "  Exponent of X, E = " << e << "\n";

  table = hf_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  delete [] table;

  return;
}
//****************************************************************************80

void hermite_polynomial_test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST14 tests H_POLYNOMIAL_COEFFICIENTS.
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
  double *c;
  int i;
  int j;
  int n = 10;

  cout << "\n";
  cout << "HERMITE_POLYNOMIAL_TEST14\n";
  cout << "  H_POLYNOMIAL_COEFFICIENTS determines physicist's Hermite polynomial coefficients.\n";

  c = h_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    cout << "\n";
    cout << "  H(" << i << ",x) =\n";
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

void hermite_polynomial_test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST15 tests HE_POLYNOMIAL_COEFFICIENTS.
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
  double *c;
  int i;
  int j;
  int n = 10;

  cout << "\n";
  cout << "HERMITE_POLYNOMIAL_TEST15\n";
  cout << "  HE_POLYNOMIAL_COEFFICIENTS determines probabilist's Hermite polynomial coefficients.\n";

  c = he_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    cout << "\n";
    cout << "  He(" << i << ") =\n";
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

string i4_to_string ( int i4 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  ostringstream fred;
  string value;

  fred << i4;

  value = fred.str ( );

  return value;
}

