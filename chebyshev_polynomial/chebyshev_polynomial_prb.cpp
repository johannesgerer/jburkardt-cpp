# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "chebyshev_polynomial.hpp"

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
//    MAIN is the main program for CHEBYSHEV_POLYNOMIAL_PRB.
//
//  Discussion:
//
//    CHEBYSHEV_POLYNOMIAL_PRB tests the CHEBYSHEV_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CHEBYSHEV_POLYNOMIAL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CHEBYSHEV_POLYNOMIAL library.\n";
 
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
  cout << "CHEBYSHEV_POLYNOMIAL_PRB\n";
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
//    TEST01 tests T_PROJECT_COEFFICIENTS_DATA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  double *d;
  double *d2;
  int i;
  int m;
  int n;
  int seed;
  double *x;

  cout << "\n";
  cout << "CHEBYSHEV_POLYNOMIAL_TEST01:\n";
  cout << "  T_PROJECT_COEFFICIENTS_DATA estimates the Chebyshev polynomial\n";
  cout << "  coefficients for a function given as data (x,fx).\n";
  cout << "\n";
  cout << "  Here, we use fx = f(x) = x^2 for the data.\n";
  cout << "\n";
  cout << "  Since T(0,x) = 1 and T(2,x) = 2*x^2 - 1, the correct expansion is\n";
  cout << "  f(x) = 1/2 T(0,x) + 0 T(1,x) + 1/2 T(2,x) + 0 * all other polys.\n";
//
//  Data in [0,1];
//
  a = 0.0;
  b = 1.0;
  m = 20;
  seed = 123456789;
  x = r8vec_uniform_01_new ( m, &seed );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = x[i] * x[i];
  }

  r8vec2_print ( m, x, d, "  Data ( X, D ):" );

  n = 4;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  
  r8vec_print ( n, c, "  Coefficients of Chebyshev expansion of degree 4." );
//
//  Compare Chebyshev expansion and original function.
//
  d2 = t_project_value ( m, n, x, c );

  cout << "\n";
  cout << "   I      X(I)     Data(I)      Chebyshev(X(I))\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << d[i]
         << "  " << setw(12) << d2[i] << "\n";
  }
  delete [] c;
  delete [] d;
  delete [] d2;
  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests T_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int j;
  int n;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  T_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients \n";
  cout << "  ot T(n,x).\n";

  n = 5;

  c = t_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    cout << "\n";
    cout << "  T(" << i << ",x)\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] != 0.0 )
      {
        if ( j == 0 )
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
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests T_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  T_POLYNOMIAL evaluates the Chebyshev polynomial T(n,x).\n";
  cout << "\n";
  cout << "                   Tabulated      Computed\n";
  cout << "     N      X        T(n,x)        T(n,x)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    t_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = t_polynomial ( 1, n, x_vec );

    cout << "  " << setw(8)  << n
         << "  " << setw(8)  << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests T_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *fx;
  int i;
  int n;
  int n_max = 5;
  double *z;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  T_POLYNOMIAL_ZEROS returns zeroes of T(n,x).\n";
  cout << "\n";
  cout << "       N      X        T(n,x)\n";
  cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    z = t_polynomial_zeros ( n );
    fx = t_polynomial ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << n
           << "  " << setw(8) << z[i]
           << "  " << setw(14) << fx[i+n*n] << "\n";
    }
    cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests T_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
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
  cout << "TEST05:\n";
  cout << "  T_QUADRATURE_RULE computes the quadrature rule\n";
  cout << "  associated with T(n,x);\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  t_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "    N      X            W" );

  cout << "\n";
  cout << "  Use the quadrature rule to estimate:\n";
  cout << "\n";
  cout << "    Q = Integral ( -1 <= X <= +1 ) X^E / sqrt ( 1-x^2) dx\n";
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
    q_exact = t_integral ( e );
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

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_POLYNOMIAL_TEST06 tests the projection of T(i,x) and T(j,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int j;
  int k;
  int n;
  double *phi;
  string title;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST06:\n";
  cout << "  As a sanity check, make sure that the projection of:\n";
  cout << "  T(i,x) onto T(j,x) is:\n";
  cout << "  0 if i is not equal to j;\n";
  cout << "  pi if i = j = 0;\n";
  cout << "  pi/2 if i = j =/= 0.\n";

  n = 3;

  x = new double[n+1];
  w = new double[n+1];

  t_quadrature_rule ( n + 1, x, w );

  c = new double[n+1];

  phi = t_polynomial ( n + 1, n, x );
 
  for ( j = 0; j <= n; j++ )
  {

    for ( i = 0; i <= n; i++ )
    {
      c[i] = 0.0;
      for ( k = 0; k <= n; k++ )
      {
        c[i] = c[i] + w[k] * phi[k+i*(n+1)] * phi[k+j*(n+1)];
      }
    }

    title = "  Chebyshev expansion coefficients for T(" + i4_to_string ( j, "%d" ) + ",x)";
    r8vec_print ( n + 1, c, title );
  }

  delete [] c;
  delete [] phi;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_POLYNOMIAL_TEST07 tests T_PROJECT_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  int n;

  cout << "\n";
  cout << "CHEBYSHEV_POLYNOMIAL_TEST07:\n";
  cout << "  T_PROJECT_COEFFICIENTS computes the Chebyshev coefficients\n";
  cout << "  of a function defined over [-1,+1].\n";
  cout << "  T_PROJECT_COEFFICIENTS_AB works in [A,B].\n";

  n = 3;
  c = new double[n+1];
  c = t_project_coefficients ( n, exp );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) in [-1,+1]" );
  delete [] c;

  n = 5;
  c = new double[n+1];
  c = t_project_coefficients ( n, exp );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) in [-1,+1]" );
  delete [] c;

  n = 5;
  c = new double[n+1];
  c = t_project_coefficients ( n, sin );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) in [-1,+1]" );
  delete [] c;
//
//  Repeat calculation with T_PROJECT_COEFFICIENTS_AB.
//
  n = 5;
  c = new double[n+1];
  a = -1.0;
  b = +1.0;
  c = t_project_coefficients_ab ( n, sin, a, b );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) in [-1,+1]" );
  delete [] c;
//
//  Now try a different interval.
//
  n = 5;
  c = new double[n+1];
  a = 0.0;
  b = 1.0;
  c = t_project_coefficients_ab ( n, sqrt, a, b );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sqrt(x) in [0,+1]" );
  delete [] c;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests T_PROJECT_COEFFICIENTS_DATA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  double *d;
  int i;
  int m;
  int n;
  int seed;
  double *x;

  cout << "\n";
  cout << "TEST08:\n";
  cout << "  T_PROJECT_COEFFICIENTS_DATA computes the Chebyshev\n";
  cout << "  coefficients of a function defined by data.\n";
  cout << "\n";
  cout << "  We are looking for an approximation that is good in [-1,+1].\n";
  cout << "\n";
  cout << "  Begin by using equally spaced points in [-1,+1].\n";

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = exp ( x[i] );
  }
  n = 3;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) on [-1,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = exp ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) on [-1,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [-1,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  cout << "\n";
  cout << "  Now sample equally spaced points in [0,+1].\n";
  cout << "  The approximation still applies to the interval [-1,+1].\n";

  a = 0.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [0,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  a = 0.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = sqrt ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sqrt(x) on [0,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  cout << "\n";
  cout << "  Now random points in [-1,+1].\n";

  a = -1.0;
  b = +1.0;
  m = 10;
  seed = 123456789;
  x = r8vec_uniform_new ( m, a, b, &seed );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [-1,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 compares a function and projection over [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  int i;
  int m;
  int n;
  double r;
  double *v;
  double *x;

  cout << "\n";
  cout << "TEST09:\n";
  cout << "  T_PROJECT_COEFFICIENTS computes the Chebyshev interpolant C(F)(n,x)\n";
  cout << "  of a function F(x) defined over [-1,+1].\n";
  cout << "  T_PROJECT_VALUE evaluates that projection.\n";

  cout << "\n";
  cout << "  Compute projections of order N to exp(x) over [-1,+1],\n";
  cout << "\n";
  cout << "   N   Max||F(x)-C(F)(n,x)||\n";
  cout << "\n";

  a = -1.0;
  b = +1.0;

  for ( n = 0; n <= 10; n++ )
  {
    c = t_project_coefficients ( n, exp );
    m = 101;
    x = r8vec_linspace_new ( m, a, b );
    v = t_project_value ( m, n, x, c );
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r8_max ( r, fabs ( v[i] - exp ( x[i] ) ) );
    }
    cout << "  " << setw(2) << n
         << "  " << setw(14) << r << "\n";
    delete [] c;
    delete [] v;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 compares a function and projection over [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  int i;
  int m;
  int n;
  double r;
  double *v;
  double *x;

  cout << "\n";
  cout << "TEST10:\n";
  cout << "  T_PROJECT_COEFFICIENTS_AB computes the Chebyshev interpolant C(F)(n,x)\n";
  cout << "  of a function F(x) defined over [A,B].\n";
  cout << "  T_PROJECT_VALUE_AB evaluates that projection.\n";

  a = 0.0;
  b = 1.5;

  cout << "\n";
  cout << "  Compute projections of order N to exp(x) over [" << a << "," << b << "]\n";
  cout << "\n";
  cout << "   N   Max||F(x)-C(F)(n,x)||\n";
  cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    c = t_project_coefficients_ab ( n, exp, a, b );
    m = 101;
    x = r8vec_linspace_new ( m, a, b );
    v = t_project_value_ab ( m, n, x, c, a, b );
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r8_max ( r, fabs ( v[i] - exp ( x[i] ) ) );
    }
    cout << "  " << setw(2) << n
         << "  " << setw(14) << r << "\n";
    delete [] c;
    delete [] v;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests U_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int j;
  int n;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  U_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients \n";
  cout << "  of U(n,x).\n";

  n = 5;

  c = u_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    cout << "\n";
    cout << "  U(" << i << ",x)\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] != 0.0 )
      {
        if ( j == 0 )
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
  }
  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests U_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "TEST12:\n";
  cout << "  U_POLYNOMIAL evaluates the Chebyshev polynomial U(n,x).\n";
  cout << "\n";
  cout << "                   Tabulated      Computed\n";
  cout << "     N      X        U(n,x)        U(n,x)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    u_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = u_polynomial ( 1, n, x_vec );

    cout << "  " << setw(8)  << n
         << "  " << setw(8)  << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests U_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *fx;
  int i;
  int n;
  int n_max = 5;
  double *z;

  cout << "\n";
  cout << "TEST13:\n";
  cout << "  U_POLYNOMIAL_ZEROS returns zeroes of U(n,x).\n";
  cout << "\n";
  cout << "       N      X        U(n,x)\n";
  cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    z = u_polynomial_zeros ( n );
    fx = u_polynomial ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << n
           << "  " << setw(8) << z[i]
           << "  " << setw(14) << fx[i+n*n] << "\n";
    }
    cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests U_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 April 2012
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
  cout << "TEST14:\n";
  cout << "  U_QUADRATURE_RULE computes the quadrature rule\n";
  cout << "  associated with U(n,x);\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  u_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "    N      X            W" );

  cout << "\n";
  cout << "  Use the quadrature rule to estimate:\n";
  cout << "\n";
  cout << "    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt ( 1-x^2) dx\n";
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
    q_exact = u_integral ( e );
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

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests V_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "TEST15:\n";
  cout << "  V_POLYNOMIAL evaluates the Chebyshev polynomial V(n,x).\n";
  cout << "\n";
  cout << "                   Tabulated      Computed\n";
  cout << "     N      X        V(n,x)        V(n,x)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    v_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = v_polynomial ( 1, n, x_vec );

    cout << "  " << setw(8)  << n
         << "  " << setw(8)  << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests W_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "TEST16:\n";
  cout << "  W_POLYNOMIAL evaluates the Chebyshev polynomial W(n,x).\n";
  cout << "\n";
  cout << "                   Tabulated      Computed\n";
  cout << "     N      X        W(n,x)        W(n,x)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    w_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = w_polynomial ( 1, n, x_vec );

    cout << "  " << setw(8)  << n
         << "  " << setw(8)  << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests T_TRIPLE_PRODUCT_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int i;
  int j;
  int k;
  int l;
  int n;
  int seed;
  int test;
  int test_num = 20;
  double ti;
  double tj;
  double tk;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST17:\n";
  cout << "  T_TRIPLE_PRODUCT_INTEGRAL computes the triple integral\n";
  cout << "  Tijk = integral ( -1 <= x <= 1 ) T(i,x) T(j,x) T(k,x) / sqrt ( 1-x^2) dx\n";
  cout << "\n";
  cout << "   I   J   K     Tijk           Tijk\n";
  cout << "                 computed       exact\n";
  cout << "\n";

  n = 15;
  x = new double[n];
  w = new double[n];

  t_quadrature_rule ( n, x, w );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform ( 2, 6, seed );
    j = i4_uniform ( 1, 3, seed );
    k = i4_uniform ( 0, 4, seed );
    fx1 = t_triple_product_integral ( i, j, k );
    fx2 = 0.0;
    for ( l = 0; l < n; l++ )
    {
      ti = t_polynomial_value ( i, x[l] );
      tj = t_polynomial_value ( j, x[l] );
      tk = t_polynomial_value ( k, x[l] );
      fx2 = fx2 + w[l] * ti * tj * tk;
    }
    cout << "  " << setw(2) << i
         << "  " << setw(2) << j
         << "  " << setw(2) << k
         << "  " << setw(14) << fx1
         << "  " << setw(14) << fx2 << "\n";
  }

  delete [] x;
  delete [] w;

  return;
}
