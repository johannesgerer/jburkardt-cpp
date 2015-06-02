# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "quadrule.hpp"

int main ( );
void bashforth_set_test ( );
void test02 ( );
void test03 ( );
void test04 ( );
void chebyshev1_compute_test ( );
void test06 ( );
void test065 ( int n );
void chebyshev3_compute_test ( );
void clenshaw_curtis_compute_test ( );
void clenshaw_curtis_set_test ( );
void fejer1_compute_test ( );
void fejer1_set_test ( );
void fejer2_compute_test ( );
void fejer2_set_test ( );
void gegenbauer_compute_test ( int order, double alpha );
void test08 ( );
void test085 ( );
void test087 ( );
void test089 ( );
void test09 ( );
void test095 ( );
void test096 ( );
void test10 ( );
void test105 ( );
void test108 ( int n );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test165 ( int order, double alpha );
void test17 ( );
void test18 ( int n );
void test185 ( int n );
void test19 ( );
void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );
void test28 ( );
void test29 ( );
void test30 ( );
void test31 ( );
void test32 ( );
void test33 ( );
void test34 ( );
void lobatto_compute_test ( );
void lobatto_set_test ( );
void moulton_set_test ( );
void ncc_set_test ( );
void test38 ( );
void test39 ( );
void test40 ( );
void test401 ( );
void test402 ( );
void test403 ( );
void test404 ( );
void test41 ( );
double f1sd1 ( double x );
char *function_name ( int function_index );
void function_set ( string action, int *i );
double function_value ( double x );
double fx1sd1 ( double x );
double fx2sd1 ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QUADRULE_PRB.
//
//  Discussion:
//
//    QUADRULE_PRB tests the QUADRULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  double alpha;
  int n;
  int order;

  cout << "\n";
  cout << "QUADRULE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the QUADRULE library.\n";

  bashforth_set_test ( );
  test02 ( );
  test03 ( );
  test04 ( );
  chebyshev1_compute_test ( );
  test06 ( );
  n = 10;
  test065 ( n );
  chebyshev3_compute_test ( );
  clenshaw_curtis_compute_test ( );
  clenshaw_curtis_set_test ( );
  fejer1_compute_test ( );
  fejer1_set_test ( );
  fejer2_compute_test ( );
  fejer2_set_test ( );

  n = 5;
  alpha = 0.5;
  gegenbauer_compute_test ( n, alpha );

  n = 10;
  alpha = - 0.5;
  gegenbauer_compute_test ( n, alpha );

  test08 ( );
  test085 ( );
  test087 ( );
  test089 ( );
  test09 ( );
  test095 ( );
  test096 ( );

  test10 ( );
  test105 ( );
  n = 10;
  test108 ( n );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );

  n = 11;
  alpha = 0.0;
  test165 ( n, alpha );

  n = 11;
  alpha = 0.5;
  test165 ( n, alpha );

  n = 11;
  alpha = 2.0;
  test165 ( n, alpha );

  test17 ( );
//
//  Compare computed and lookup versions of Gauss-Legendre rules.
//
  n = 31;
  test18 ( n );
  n = 64;
  test18 ( n );
  n = 129;
  test18 ( n );
  n = 255;
  test18 ( n );

  n = 31;
  test185 ( n );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  lobatto_compute_test ( );
  lobatto_set_test ( );
  moulton_set_test ( );
  ncc_set_test ( );
  test38 ( );
  test39 ( );

  test40 ( );
  test401 ( );
  test402 ( );
  test403 ( );
  test404 ( );
  test41 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QUADRULE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );
 
  return 0;
}
//****************************************************************************80

void bashforth_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASHFORTH_SET_TEST tests BASHFORTH_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  cout << "\n";
  cout << "BASHFORTH_SET_TEST\n";
  cout << "  BASHFORTH_SET sets up an Adams-Bashforth rule;\n";
  cout << "\n";
  cout << "  Index             X                   W\n";
  cout << "\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    bashforth_set ( n, x, w );

    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i 
           << "  " << setw(24) << x[i]
           << "  " << setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
 
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BDFC_SET and BDF_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int n;
  int n_max = 10;
  double *result;
  double *w;
  double *x;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  cout << "\n";
  cout << "TEST02\n";
  cout << "  BDFC_SET sets up a Backward Difference Corrector rule;\n";
  cout << "  BDF_SUM carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [0,1].\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( n = 1; n <= n_max; n++ )
    {
      x = new double[n];
      w = new double[n];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        bdfc_set ( n, x, w );
 
        result[i] = bdf_sum ( function_value, n, x, w );
 
      }
      cout << setw(2) << n << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] x;
      delete [] w;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests BDFP_SET and BDF_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 10;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  cout << "\n";
  cout << "TEST03\n";
  cout << "  BDFP_SET sets up a Backward Difference Predictor rule;\n";
  cout << "  BDF_SUM carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [0,1].\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        bdfp_set ( order, xtab, weight );
 
        result[i] = bdf_sum ( function_value, order, xtab, weight );
 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests CHEBYSHEV_SET and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int n;
  int n_max = 9;
  int nsub;
  double *result;
  double *w;
  double *x;
  double xhi;
  double xlo;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  CHEBYSHEV_SET sets up a Chebyshev rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( n = 1; n <= n_max; n++ )
    {
      if ( n == 8 )
      {
        continue;
      }

      x = new double[n];
      w = new double[n];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        chebyshev_set ( n, x, w );
 
        result[i] = sum_sub ( function_value, a, b, nsub, n, 
          xlo, xhi, x, w ); 
      }
      cout << setw(2) << n << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] x;
      delete [] w;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void chebyshev1_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_COMPUTE_TEST tests CHEBYSHEV1_COMPUTE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  cout << "\n";
  cout << "CHEBYSHEV1_COMPUTE_TEST\n";
  cout << "  CHEBYSHEV1_COMPUTE computes\n";
  cout << "  a Chebyshev Type 1 quadrature rule over [-1,1]\n";
  cout << "  of given order.\n";

  cout << "\n";
  cout << "  Index             X                   W\n";
  cout << "\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    chebyshev1_compute ( n, x, w );

    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i 
           << "  " << setw(24) << x[i]
           << "  " << setw(24) << w[i] << "\n";
    }
    delete [] w;
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
//    TEST06 tests CHEBYSHEV2_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int j;
  int nsub;
  int order;
  int order_max = 4;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = -1.0;
  b =  1.0;
  nsub = 1;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  CHEBYSHEV2_COMPUTE sets up a Gauss-Chebyshev type 2 rule,\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "  The weight function is 1 / sqrt ( 1 - X**2 )\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      if ( order == 8 )
      {
        continue;
      }

      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        chebyshev2_compute ( order, xtab, weight );
 
        result[i] = 0.0;
        for ( j = 0; j < order; j++ )
        {
          result[i] = result[i] + weight[j] * function_value ( xtab[j] );
        }  
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(6) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;

  return;
}
//****************************************************************************80

void test065 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST065 uses CHEBSHEV2_COMPUTE to integral over the semicircle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  double *f;
  int i;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST065\n";
  cout << "  Approximate the integral of f(x,y) over the semicircle\n";
  cout << "    -1 <= x <= 1, y = sqrt ( 1 - x^2 )\n";
  cout << "  using N Chebyshev points.\n";
  cout << "  If p(x,y) involves any term of odd degree in y,\n";
  cout << "  the estimate will only be approximate.\n";
  cout << "\n";
  cout << "  Polynomial    N    Integral        Estimate       Error\n";
  cout << "\n";

  f = new double[n];
  x = new double[n];
  w = new double[n];

  chebyshev2_compute ( n, x, w );
//
//  f(x,y) = 1
//
  exact = 1.5707963267948966192;
  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0;
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  cout << "  1            "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x
//
  exact = 0.0;
  for ( i = 0; i < n; i++ )
  {
    f[i] = x[i];
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  cout << "  x            "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = y = sqrt ( 1 - x^2 )
//
  exact = 0.66666666666666666667;
  for ( i = 0; i < n; i++ )
  {
    f[i] = sqrt ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 2.0;
  error = r8_abs ( q - exact );
  cout << "     y         "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^2
//
  exact = 0.39269908169872415481;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 );
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  cout << "  x^2          "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = xy = x * sqrt ( 1 - x^2 )
//
  exact = 0.0;
  for ( i = 0; i < n; i++ )
  {
    f[i] = x[i] * sqrt ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 2.0;
  error = r8_abs ( q - exact );
  cout << "  x  y         "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = y^2 -> ( 1 - x^2 )
//
  exact = 0.39269908169872415481;
  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0 - pow ( x[i], 2 );
  }
  q = r8vec_dot_product ( n, w, f ) / 3.0;
  error = r8_abs ( q - exact );
  cout << "     y^2       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^3
//
  exact = 0.0;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 3 );
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  cout << "  x^3          "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^2 y = x^2 sqrt ( 1 - x^2 )
//
  exact = 0.13333333333333333333;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 ) * sqrt ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 2.0;
  error = r8_abs ( q - exact );
  cout << "  x^2y         "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x y^2 = x * ( 1 - x^2 )
//
  exact = 0.0;
  for ( i = 0; i < n; i++ )
  {
    f[i] = x[i] * ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 3.0;
  error = r8_abs ( q - exact );
  cout << "  x  y^2       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = y^3
//
  exact = 0.26666666666666666667;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( 1.0 - pow ( x[i], 2 ), 1.5 );
  }
  q = r8vec_dot_product ( n, w, f ) / 4.0;
  error = r8_abs ( q - exact );
  cout << "     y^3       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^4
//
  exact = 0.19634954084936207740;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 4 );
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  cout << "  x^4          "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^2y^2 -> x^2( 1 - x^2 )
//
  exact = 0.065449846949787359135;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 ) * ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 3.0;
  error = r8_abs ( q - exact );
  cout << "  x^2y^2       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = y^4 -> ( 1 - x^2 )^2
//
  exact = 0.19634954084936207740;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( 1.0 - pow ( x[i], 2 ), 2 );
  }
  q = r8vec_dot_product ( n, w, f ) / 5.0;
  error = r8_abs ( q - exact );
  cout << "     y^4       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^4y = x^4 sqrt ( 1 - x^2 )
//
  exact = 0.057142857142857142857;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 4 ) * sqrt ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 2.0;
  error = r8_abs ( q - exact );
  cout << "  x^4y         "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  x^2y^3 = x^2 ( 1 - x^2 )^(3/2)
//
  exact = 0.038095238095238095238;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 ) * pow ( 1.0 - pow ( x[i], 2 ), 1.5 );
  }
  q = r8vec_dot_product ( n, w, f ) / 4.0;
  error = r8_abs ( q - exact );
  cout << "  x^2y^3       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = y^5
//
  exact = 0.15238095238095238095;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( 1.0 - pow ( x[i], 2 ), 2.5 );
  }
  q = r8vec_dot_product ( n, w, f ) / 6.0;
  error = r8_abs ( q - exact );
  cout << "     y^5       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^6
//
  exact = 0.12271846303085129838;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 6 );
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  cout << "  x^6          "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^4y^2 -> x^4( 1 - x^2 )
//
  exact = 0.024543692606170259675;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 4 ) * ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 3.0;
  error = r8_abs ( q - exact );
  cout << "  x^4y^2       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = x^2y^4 -> x^2( 1 - x^2 )^2
//
  exact = 0.024543692606170259675;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 ) * pow ( 1.0 - pow ( x[i], 2 ), 2 );
  }
  q = r8vec_dot_product ( n, w, f ) / 5.0;
  error = r8_abs ( q - exact );
  cout << "  x^2y^4       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";
//
//  f(x,y) = y^6 -> ( 1 - x^2 )^3
//
  exact = 0.12271846303085129838;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( 1.0 - pow ( x[i], 2 ), 3 );
  }
  q = r8vec_dot_product ( n, w, f ) / 7.0;
  error = r8_abs ( q - exact );
  cout << "     y^6       "
       << "  " << setw(2) << n
       << "  " << setw(14) << exact
       << "  " << setw(14) << q
       << "  " << setw(14) << error << "\n";

  delete [] f;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void chebyshev3_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV3_COMPUTE_TEST tests CHEBYSHEV3_COMPUTE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  cout << "\n";
  cout << "CHEBYSHEV3_COMPUTE_TEST\n";
  cout << "  CHEBYSHEV3_COMPUTE computes\n";
  cout << "  a Chebyshev Type 3 quadrature rule over [-1,1]\n";
  cout << "  of given order.\n";

  cout << "\n";
  cout << "  Index             X                   W\n";
  cout << "\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    chebyshev3_compute ( n, x, w );

    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i 
           << "  " << setw(24) << x[i]
           << "  " << setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void clenshaw_curtis_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE_TEST tests CLENSHAW_CURTIS_COMPUTE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  cout << "\n";
  cout << "CLENSHAW_CURTIS_COMPUTE_TEST\n";
  cout << "  CLENSHAW_CURTIS_COMPUTE computes\n";
  cout << "  a Clenshaw-Curtis quadrature rule over [-1,1]\n";
  cout << "  of given order.\n";

  cout << "\n";
  cout << "  Index             X                   W\n";
  cout << "\n";

  for ( order = 1; order <= order_max; order++ )
  {
    w = new double[order];
    x = new double[order];

    clenshaw_curtis_compute ( order, x, w );

    cout << "\n";

    for ( i = 0; i < order; i++ )
    {
      cout << "  " << setw(2) << i 
           << "  " << setw(24) << x[i]
           << "  " << setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void clenshaw_curtis_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_SET_TEST tests CLENSHAW_CURTIS_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double e;
  double exact;
  int i;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "CLENSHAW_CURTIS_SET_TEST\n";
  cout << "  CLENSHAW_CURTIS_SET sets up a Clenshaw-Curtis rule;\n";
  cout << "\n";
  cout << "  Estimate the integral of sqrt(abs(x)) over [-1,+1].\n";
  cout << "\n";
  cout << "   N           Estimate             Error\n";
  cout << "\n";

  exact = 4.0 / 3.0;

  for ( n = 1; n <= 10; n++ )
  {
    x = new double[n];
    w = new double[n];

    clenshaw_curtis_set ( n, x, w );

    q = 0.0;
    for ( i = 0; i < n; i++ )
    {
      q = q + w[i] * sqrt ( fabs ( x[i] ) );
    }

    e = fabs ( q - exact );

    cout << "  " << setw(2) << n
         << "  " << setw(24) << q
         << "  " << setw(14) << e << "\n";

    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void fejer1_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER1_COMPUTE_TEST tests FEJER1_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_max = 10;
  double *w;
  double *x;

  cout << "\n";
  cout << "FEJER1_COMPUTE_TEST\n";
  cout << "  FEJER1_COMPUTE computes a Fejer type 1 quadrature rule;\n";
  cout << "\n";
  cout << "     Order        W               X\n";
  cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    w = new double[n];
    x = new double[n];

    fejer1_compute ( n, x, w );

    cout << "\n";
    cout << "  " << setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "          "
           << "  " << setw(14) << w[i]
           << "  " << setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void fejer1_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER1_SET_TEST tests FEJER1_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_max = 10;
  double *w;
  double *x;

  cout << "\n";
  cout << "FEJER1_SET_TEST\n";
  cout << "  FEJER1_SET sets a Fejer type 1 quadrature rule;\n";
  cout << "\n";
  cout << "     Order        W               X\n";
  cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    w = new double[n];
    x = new double[n];

    fejer1_set ( n, x, w );

    cout << "\n";
    cout << "  " << setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "          "
           << "  " << setw(14) << w[i]
           << "  " << setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void fejer2_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER2_COMPUTE_TEST tests FEJER2_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_max = 10;
  double *w;
  double *x;

  cout << "\n";
  cout << "FEJER2_COMPUTE_TEST\n";
  cout << "  FEJER2_COMPUTE computes a Fejer type 2 quadrature rule;\n";
  cout << "\n";
  cout << "     Order        W               X\n";
  cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    w = new double[n];
    x = new double[n];

    fejer2_compute ( n, x, w );

    cout << "\n";
    cout << "  " << setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "          "
           << "  " << setw(14) << w[i]
           << "  " << setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void fejer2_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER2_SET_TEST tests FEJER2_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_max = 10;
  double *w;
  double *x;

  cout << "\n";
  cout << "FEJER2_SET_TEST\n";
  cout << "  FEJER2_SET sets a Fejer type 2 quadrature rule;\n";
  cout << "\n";
  cout << "     Order        W               X\n";
  cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    w = new double[n];
    x = new double[n];

    fejer2_set ( n, x, w );

    cout << "\n";
    cout << "  " << setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "          "
           << "  " << setw(14) << w[i]
           << "  " << setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void gegenbauer_compute_test ( int order, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_COMPUTE_TEST tests GEGENBAUER_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double ALPHA, the parameter.
//
{
  int i;
  double *w;
  double *x;

  cout << "\n";
  cout << "GEGENBAUER_COMPUTE_TEST\n";
  cout << "  GEGENBAUER_COMPUTE computes a Gauss-Gegenbauer rule;\n";
  cout << "\n";
  cout << "  The printed output of this routine can be inserted into\n";
  cout << "  a C++ program.\n";

  w = new double[order];
  x = new double[order];

  gegenbauer_compute ( order, alpha, x, w );

  cout << "//\n";
  cout << "//  Abscissas X and weights W for a Gauss Gegenbauer rule\n";
  cout << "//  of ORDER   = " << order << "\n";
  cout << "//  with ALPHA = " << alpha << "\n";
  cout << "//\n";

  for ( i = 0; i < order; i++ )
  {
    cout << "    x[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << x[i] << ";\n";
  }
  cout << "\n";
  for ( i = 0; i < order; i++ )
  {
    cout << "    w[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << w[i] << ";\n";
  }
//
//  Free memory.
//
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests HERMITE_EK_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int j;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  cout << "\n";
  cout << "TEST08\n";
  cout << "  HERMITE_EK_COMPUTE computes a Gauss-Hermite rule;\n";
  cout << "\n";
  cout << "  The integration interval is ( -oo, +oo ).\n";
  cout << "  The weight function is exp ( - x * x )\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        hermite_ek_compute ( order, xtab, weight );
 
        result[i] = 0.0;
        for ( j = 0; j < order; j++ )
        {
          result[i] = result[i] + weight[j] * function_value ( xtab[j] );
        } 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests HERMITE_EK_COMPUTE against HERMITE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  double *f_vec;
  int i;
  int n;
  int order;
  int order_max = 10;
  double *weight;
  double *xtab;

  cout << "\n";
  cout << "TEST085\n";
  cout << "  HERMITE_EK_COMPUTE computes a Gauss-Hermite rule\n";
  cout << "  which is appropriate for integrands of the form\n";
  cout << "    f(x) * exp(-x*x) from -oo to +oo.\n";
  cout << "\n";
  cout << "  HERMITE_INTEGRAL determines the exact value of\n";
  cout << "  this integal when f(x) = x^n.\n";
  cout << "\n";
  cout << "         N     Order       Estimate       Exact            Error\n";

  for ( n = 0; n <= 10; n = n + 2 )
  {
    exact = hermite_integral ( n );

    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      f_vec = new double[order];
      weight = new double[order];
      xtab = new double[order];

      hermite_ek_compute ( order, xtab, weight );
 
      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f_vec[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f_vec[i] = pow ( xtab[i], n );
        }
      }
      estimate = r8vec_dot_product ( order, weight, f_vec );

      error = r8_abs ( exact - estimate );
  
      cout << "  " << setw(8)  << n
           << "  " << setw(8)  << order
           << "  " << setw(14) << estimate
           << "  " << setw(14) << exact
           << "  " << setw(14) << error << "\n";

      delete [] f_vec;
      delete [] weight;
      delete [] xtab;
    }
  }
 
  return;
}
//****************************************************************************80

void test087 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST087 tests HERMITE_EK_COMPUTE.
//
//  Discussion:
//
//    I used this test to generate tabular values of weights and
//    abscissas for Gauss-Hermite quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 31;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST087\n";
  cout << "  HERMITE_EK_COMPUTE computes a Gauss-Hermite rule;\n";
  cout << "\n";
  cout << "  Compute the data for N = " << n << "\n";

  w = new double[n];
  x = new double[n];

  hermite_ek_compute ( n, x, w );
 
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "    x[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << x[i] << ";\n";
  }
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "    w[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << w[i] << ";\n";
  }
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test089 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST089 tests HERMITE_PROBABILIST_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int j;
  int order;
  int order_max = 10;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  cout << "\n";
  cout << "TEST089\n";
  cout << "  HERMITE_PROBABILIST_SET sets up a Hermite probabilist rule;\n";
  cout << "\n";
  cout << "  The integration interval is ( -oo, +oo ).\n";
  cout << "  The weight function is exp ( - x * x / 2 ) / sqrt ( 2 * pi ).\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        hermite_probabilist_set ( order, xtab, weight );
 
        result[i] = 0.0;
        for ( j = 0; j < order; j++ )
        {
          result[i] = result[i] + weight[j] * function_value ( xtab[j] );
        }
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests HERMITE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int j;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  cout << "\n";
  cout << "TEST09\n";
  cout << "  HERMITE_SET sets up a Gauss-Hermite rule;\n";
  cout << "\n";
  cout << "  The integration interval is ( -oo, +oo ).\n";
  cout << "  The weight function is exp ( - x * x )\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        hermite_set ( order, xtab, weight );
 
        result[i] = 0.0;
        for ( j = 0; j < order; j++ )
        {
          result[i] = result[i] + weight[j] * function_value ( xtab[j] );
        }
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test095 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST095 tests HERMITE_GK16_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
{
# define L_MAX 8

  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int l;
  int m;
  int n;
  int n_list[L_MAX+1] = { 1, 3, 7, 9, 17, 19, 31, 33, 35 };
  int p;
  int p_list[L_MAX+1] = { 1, 5, 7, 15, 17, 29, 31, 33, 51 };
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST095\n";
  cout << "  HERMITE_GK16_SET sets up a nested rule\n";
  cout << "  for the Hermite integration problem.\n";
  cout << "\n";
  cout << "  The integration interval is ( -oo, +oo ).\n";
  cout << "  The weight function is exp ( - x * x )\n";
  cout << "\n";
  cout << "  HERMITE_INTEGRAL determines the exact value of\n";
  cout << "  the integal when f(x) = x^m.\n";
  cout << "\n";
  cout << "         M         N       Estimate       Exact            Error\n";

  for ( l = 0; l <= L_MAX; l++ )
  {
    cout << "\n";
    n = n_list[l];

    f = new double[n];
    x = new double[n];
    w = new double[n];

    p = p_list[l];

    hermite_gk16_set ( n, x, w );

    for ( m = 0; m <= i4_min ( p + 2, 20 ); m = m + 2 )
    {
      exact = hermite_integral ( m );

      if ( m == 0 )
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
          f[i] = pow ( x[i], m );
        }
      }

      estimate = r8vec_dot_product ( n, w, f );

      error = r8_abs ( exact - estimate );
  
      cout << "  " << setw(8) << m
           << "  " << setw(8) << n
           << "  " << setw(14) << estimate
           << "  " << setw(14) << exact
           << "  " << setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }
  return;
# undef L_MAX
}
//****************************************************************************80

void test096 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST096 tests HERMITE_GK**_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int m;
  int n;
  int p;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST096\n";
  cout << "  HERMITE_GK**_SET sets up a nested rule\n";
  cout << "  for the Hermite integration problem.\n";
  cout << "\n";
  cout << "  The integration interval is ( -oo, +oo ).\n";
  cout << "  The weight function is exp ( - x * x )\n";
  cout << "\n";
  cout << "  HERMITE_INTEGRAL determines the exact value of\n";
  cout << "  the integal when f(x) = x^m.\n";
  cout << "\n";
  cout << "  Here, we just test the highest order rule.\n";

  cout << "\n";
  cout << "  HERMITE_GK16_SET:\n";
  cout << "\n";
  cout << "         M         N       Estimate       Exact            Error\n";
  cout << "\n";

  n = 35;
  p = 51;

  f = new double[n];
  w = new double[n];
  x = new double[n];

  hermite_gk16_set ( n, x, w );

  for ( m = 0; m <= i4_min ( p + 2, 20 ); m = m + 2 )
  {
    exact = hermite_integral ( m );

    if ( m == 0 )
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
        f[i] = pow ( x[i], m );
      }
    }

    estimate = r8vec_dot_product ( n, w, f );

    error = r8_abs ( exact - estimate );
  
    cout << "  " << setw(8) << m
         << "  " << setw(8) << n
         << "  " << setw(14) << estimate
         << "  " << setw(14) << exact
         << "  " << setw(14) << error << "\n";
  }
  delete [] f;
  delete [] w;
  delete [] x;

  cout << "\n";
  cout << "  HERMITE_GK18_SET:\n";
  cout << "\n";
  cout << "         M         N       Estimate       Exact            Error\n";
  cout << "\n";
  
  n = 37;
  p = 55;

  f = new double[n];
  w = new double[n];
  x = new double[n];

  hermite_gk18_set ( n, x, w );

  for ( m = 0; m <= i4_min ( p + 2, 20 ); m = m + 2 )
  {
    exact = hermite_integral ( m );

    if ( m == 0 )
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
        f[i] = pow ( x[i], m );
      }
    }

    estimate = r8vec_dot_product ( n, w, f );

    error = r8_abs ( exact - estimate );
  
    cout << "  " << setw(8) << m
         << "  " << setw(8) << n
         << "  " << setw(14) << estimate
         << "  " << setw(14) << exact
         << "  " << setw(14) << error << "\n";
  }
  delete [] f;
  delete [] w;
  delete [] x;

  cout << "\n";
  cout << "  HERMITE_GK22_SET:\n";
  cout << "\n";
  cout << "         M         N       Estimate       Exact            Error\n";
  cout << "\n";

  n = 41;
  p = 63;

  f = new double[n];
  w = new double[n];
  x = new double[n];

  hermite_gk22_set ( n, x, w );

  for ( m = 0; m <= i4_min ( p + 2, 20 ); m = m + 2 )
  {
    exact = hermite_integral ( m );

    if ( m == 0 )
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
        f[i] = pow ( x[i], m );
      }
    }

    estimate = r8vec_dot_product ( n, w, f );

    error = r8_abs ( exact - estimate );
  
    cout << "  " << setw(8) << m
         << "  " << setw(8) << n
         << "  " << setw(14) << estimate
         << "  " << setw(14) << exact
         << "  " << setw(14) << error << "\n";
  }
  delete [] f;
  delete [] w;
  delete [] x;

  cout << "\n";
  cout << "  HERMITE_GK24_SET:\n";
  cout << "\n";
  cout << "         M         N       Estimate       Exact            Error\n";
  cout << "\n";

  n = 43;
  p = 67;

  f = new double[n];
  w = new double[n];
  x = new double[n];

  hermite_gk24_set ( n, x, w );

  for ( m = 0; m <= i4_min ( p + 2, 20 ); m = m + 2 )
  {
    exact = hermite_integral ( m );

    if ( m == 0 )
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
        f[i] = pow ( x[i], m );
      }
    }

    estimate = r8vec_dot_product ( n, w, f );

    error = r8_abs ( exact - estimate );
  
    cout << "  " << setw(8) << m
         << "  " << setw(8) << n
         << "  " << setw(14) << estimate
         << "  " << setw(14) << exact
         << "  " << setw(14) << error << "\n";
  }
  delete [] f;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests JACOBI_EK_COMPUTE and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double alpha;
  double b;
  double beta;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int k;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  a = -1.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  JACOBI_EK_COMPUTE computes a Gauss-Jacobi rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( k = 1; k <= 2; k++ )
  {
    result = new double[function_num];

    if ( k == 1 )
    {
      alpha = 0.0;
      beta = 0.0;
    }
    else
    {
      alpha = 1.0;
      beta = 0.0;
    }
    cout << "\n";
    cout << "  ALPHA = " << alpha << "\n";
    cout << "  BETA =  " << beta << "\n";

    for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
    {
      ihi = i4_min ( ilo + 4, function_num - 1 );

      cout << "\n";
      cout << "Order  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(10) << function_name ( i ) << "    ";
      }
      cout << "\n";
      cout << "\n";

      for ( order = 1; order <= order_max; order++ )
      {
        xtab = new double[order];
        weight = new double[order];

        for ( i = ilo; i <= ihi; i++ )
        {
          function_set ( "SET", &i );

          jacobi_ek_compute ( order, alpha, beta, xtab, weight );
 
          result[i] = sum_sub ( function_value, a, b, nsub, order, 
            xlo, xhi, xtab, weight ); 
        }
        cout << setw(2) << order << "  ";
        for ( i = ilo; i <= ihi; i++ )
        {
          cout << "  " << setw(12) << setprecision(8) << result[i];
        }
        cout << "\n";

        delete [] xtab;
        delete [] weight;
      }
    }
    delete [] result;
  }
 
  return;
}
//****************************************************************************80

void test105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST105 tests JACOBI_EK_COMPUTE and JACOBI_SS_COMPUTE.
//
//  Discussion:
//
//    Compare with tabular values on page 178 of Stroud and Secrest.
//
//     In particular,
//
//             X              W
//
//     1  -0.9833999115   0.4615276287E-03
//     2  -0.9447138932   0.2732249104E-02
//     3  -0.8849310847   0.8045830455E-02
//    ..  .............   ................
//    19   0.9656375637   0.7613987785E-01
//    20   0.9934477866   0.3348337670E-01
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int order = 20;
  double *weight;
  double *xtab;

  a = -1.0;
  b =  1.0;

  cout << "\n";
  cout << "TEST105\n";
  cout << "  JACOBI_EK_COMPUTE computes a Gauss-Jacobi rule;\n";
  cout << "  JACOBI_SS_COMPUTE computes a Gauss-Jacobi rule;\n";
  cout << "  Here, we simply compute a single rule and\n";
  cout << "  print it, to check for accuracy.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";

  alpha = 1.5;
  beta = 0.5;

  cout << "\n";
  cout << "  N = " << order << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA =  " << beta << "\n";

  xtab = new double[order];
  weight = new double[order];
  
  jacobi_ek_compute ( order, alpha, beta, xtab, weight );
 
  cout << "\n";
  cout << "  JACOBI_EK_COMPUTE:\n";
  cout << "\n";
  cout << "     I        X(I)            W(I)\n";
  cout << "\n";

  for ( i = 0; i < order; i++ )
  {
    cout << "  " << setw(4)  << i+1
         << "  " << setw(14) << setprecision(8) << xtab[i]
         << "  " << setw(14) << setprecision(8) << weight[i] << "\n";
  }

  jacobi_ss_compute ( order, alpha, beta, xtab, weight );
 
  cout << "\n";
  cout << "  JACOBI_SS_COMPUTE:\n";
  cout << "\n";
  cout << "     I        X(I)            W(I)\n";
  cout << "\n";

  for ( i = 0; i < order; i++ )
  {
    cout << "  " << setw(4)  << i+1
         << "  " << setw(14) << setprecision(8) << xtab[i]
         << "  " << setw(14) << setprecision(8) << weight[i] << "\n";
  }


  delete [] xtab;
  delete [] weight;
 
  return;
}
//****************************************************************************80

void test108 ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    TEST108 tests JACOBI_EK_COMPUTE.
//
//  Discussion:
//
//    I used this test to generate tabular values of weights and
//    abscissas for Gauss-Jacobi quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  double alpha = 0.5;
  double beta  = 2.0;
  int i;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST108\n";
  cout << "  JACOBI_EK_COMPUTE computes a Gauss-Jacobi rule;\n";
  cout << "\n";
  cout << "  The printed output of this test can be inserted into\n";
  cout << "  a C++ program.\n";

  w = new double[order];
  x = new double[order];

  jacobi_ek_compute ( order, alpha, beta, x, w );

  cout << "//\n";
  cout << "//  Abscissas X and weights W for a Gauss Jacobi rule\n";
  cout << "//  of ORDER   = " << order << "\n";
  cout << "//  with ALPHA = " << alpha << "\n";
  cout << "//  and  BETA  = " << beta << "\n";
  cout << "//\n";

  cout << "\n";
  for ( i = 0; i < order; i++ )
  {
    cout << "    x[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << x[i] << ";\n";
  }
  cout << "\n";
  for ( i = 0; i < order; i++ )
  {
    cout << "    w[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << w[i] << ";\n";
  }

  delete [] w;
  delete [] x;

  return;
# undef ORDER
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests KRONROD_SET, LEGENDRE_SET and SUMMER_GK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ORDERG 10
# define ORDERK 2 * ORDERG + 1

  double resultg;
  double resultk;
  double weightg[ORDERG];
  double weightk[ORDERK];
  double xtabg[ORDERG];
  double xtabk[ORDERK];

  cout << "\n";
  cout << "TEST11\n";
  cout << "  KRONROD_SET sets up a Kronrod rule;\n";
  cout << "  LEGENDRE_SET sets up Gauss-Legendre rule;\n";
  cout << "  SUMMER_GK carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [-1, 1].\n";
  cout << "  Integrand is X**2 / SQRT ( 1.1 - X**2 ).\n";
  cout << "\n";

  legendre_set ( ORDERG, xtabg, weightg );

  kronrod_set ( ORDERK, xtabk, weightk );

  summer_gk ( fx2sd1, ORDERG, weightg, &resultg, 
    ORDERK, xtabk, weightk, &resultk );

  cout << "  " << setw(2) << ORDERG 
       << "  " << setw(16) << setprecision(10) << resultg << "\n";
  cout << "  " << setw(2) << ORDERK 
       << "  " << setw(16) << setprecision(10) << resultk << "\n";
  cout << "    "
       << "  " << setw(16) << setprecision(10) << resultg - resultk << "\n";

  return;
# undef ORDERG
# undef ORDERK
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests KRONROD_SET, LEGENDRE_SET and SUM_SUB_GK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ORDERG 7
# define ORDERK 2 * ORDERG + 1

  double a;
  double b;
  double error;
  int nsub;
  double resultg;
  double resultk;
  double weightg[ORDERG];
  double weightk[ORDERK];
  double xtabg[ORDERG];
  double xtabk[ORDERK];

  a = -1.0;
  b =   1.0;
  nsub = 5;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  KRONROD_SET sets up a Kronrod rule;\n";
  cout << "  LEGENDRE_SET sets up Gauss-Legendre rule;\n";
  cout << "  SUM_SUB_GK carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "]\n";
  cout << "  Number of subintervals is " << nsub << "\n";
  cout << "  Integrand is X**2 / SQRT ( 1.1 - X**2 ).\n";
  cout << "\n";

  legendre_set ( ORDERG, xtabg, weightg );

  kronrod_set ( ORDERK, xtabk, weightk );

  sum_sub_gk ( fx2sd1, a, b, nsub, ORDERG, weightg, &resultg, 
    ORDERK, xtabk, weightk, &resultk, &error );

  cout << "  " << setw(2) << ORDERG 
       << "  " << setw(16) << setprecision(10) << resultg << "\n";
  cout << "  " << setw(2) << ORDERK 
       << "  " << setw(16) << setprecision(10) << resultk << "\n";
  cout << "    "
       << "  " << setw(16) << setprecision(10) << error << "\n";

  return;
# undef ORDERG
# undef ORDERK
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests LAGUERRE_EK_COMPUTE and LAGUERRE_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 1.0;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  LAGUERRE_EK_COMPUTE computes a Gauss-Laguerre rule;\n";
  cout << "  LAGUERRE_SUM carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << ", +oo).\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "  The weight function is EXP ( - X ).\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        laguerre_ek_compute ( order, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests LAGUERRE_EK_COMPUTE and LAGUERRE_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  LAGUERRE_EK_COMPUTE sets up a Gauss-Laguerre rule;\n";
  cout << "  LAGUERRE_SUM carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << ", +oo).\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "  The weight function is EXP ( - X ).\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        laguerre_ek_compute ( order, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests GEN_LAGUERRE_EK_COMPUTE and LAGUERRE_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double alpha;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  alpha = 1.0;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  GEN_LAGUERRE_EK_COMPUTE computes a generalized Gauss-Laguerre rule;\n";
  cout << "  LAGUERRE_SUM carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << ", +oo).\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "  The weight function is EXP ( - X ) * X^ " << alpha << ".\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        gen_laguerre_ek_compute ( order, alpha, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests GEN_LAGUERRE_EK_COMPUTE and LAGUERRE_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double alpha;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  alpha = 2.0;

  cout << "\n";
  cout << "TEST16\n";
  cout << "  GEN_LAGUERRE_EK_COMPUTE computes a generalized Gauss-Laguerre rule;\n";
  cout << "  LAGUERRE_SUM carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << ", +oo).\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "  The weight function is EXP ( - X ) * X^ " << alpha << ".\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        gen_laguerre_ek_compute ( order, alpha, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test165 ( int order, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    TEST165 tests GEN_LAGUERRE_EK_COMPUTE.
//
//  Discussion:
//
//    I used this test to generate tabular values of weights and
//    abscissas for generalized Gauss-Laguerre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double ALPHA, the parameter.
//
{
  int i;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST165\n";
  cout << "  GEN_LAGUERRE_EK_COMPUTE computes a generalized Gauss-Laguerre rule;\n";
  cout << "\n";
  cout << "  The printed output of this routine can be inserted into\n";
  cout << "  a C++ program.\n";

  w = new double[order];
  x = new double[order];

  gen_laguerre_ek_compute ( order, alpha, x, w );

  cout << "//\n";
  cout << "//  Abscissas X and weights W for a Gauss Laguerre rule\n";
  cout << "//  of ORDER   = " << order << "\n";
  cout << "//  with ALPHA = " << alpha << "\n";
  cout << "//\n";

  for ( i = 0; i < order; i++ )
  {
    cout << "    x[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << x[i] << ";\n";
  }
  cout << "\n";
  for ( i = 0; i < order; i++ )
  {
    cout << "    w[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << w[i] << ";\n";
  }

  delete [] w;
  delete [] x;

  return;
# undef ORDER
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests LAGUERRE_SET and LAGUERRE_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  1.0;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  LAGUERRE_SET sets up a Gauss-Laguerre rule;\n";
  cout << "  LAGUERRE_SUM carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << ", +oo).\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "  The weight function is EXP ( - X ).\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        laguerre_set ( order, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test18 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 compares LEGENDRE_EK_COMPUTE and LEGENDRE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int iwdifmax;
  int ixdifmax;
  double *w1;
  double *w2;
  double wdifmax;
  double *x1;
  double *x2;
  double xdifmax;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  LEGENDRE_EK_COMPUTE computes a Gauss-Legendre rule;\n";
  cout << "  LEGENDRE_SET looks up the same data.\n";
  cout << "\n";
  cout << "  Compare the data for N = " << n << "\n";

  x1 = new double[n];
  w1 = new double[n];

  legendre_ek_compute ( n, x1, w1 );

  x2 = new double[n];
  w2 = new double[n];

  legendre_set ( n, x2, w2 );

  xdifmax = 0.0;
  ixdifmax = -1;

  wdifmax = 0.0;
  iwdifmax = -1;

  for ( i = 0; i < n; i++ )
  {
    if ( xdifmax < r8_abs ( x1[i] - x2[i] ) )
    {
      xdifmax = r8_abs ( x1[i] - x2[i] );
      ixdifmax = i;
    }

    if ( wdifmax < r8_abs ( w1[i] - w2[i] ) )
    {
      wdifmax = r8_abs ( w1[i] - w2[i] );
      iwdifmax = i;
    }

  }

  if ( -1 < ixdifmax )
  {
    cout << "\n";
    cout << "  Maximum abscissa difference is " << xdifmax << "\n";
    cout << "  for index I = " << ixdifmax << "\n";
    cout << "  Computed:" << x1[ixdifmax] << "\n";
    cout << "  Stored:  " << x2[ixdifmax] << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  The computed and stored abscissas are identical.\n";
  }

  if ( -1 < iwdifmax )
  {
    cout << "\n";
    cout << "  Maximum weight difference is   " << wdifmax << "\n";
    cout << "  for index I = " << iwdifmax << "\n";
    cout << "  Computed:" << w1[iwdifmax] << "\n";
    cout << "  Stored:  " << w2[iwdifmax] << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  The computed and stored weights are identical.\n";
  }

  delete [] w1;
  delete [] w2;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void test185 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST185 tests LEGENDRE_EK_COMPUTE.
//
//  Discussion:
//
//    I used this test to generate tabular values of weights and
//    abscissas for Gauss-Legendre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
{
  int i;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST185\n";
  cout << "  LEGENDRE_EK_COMPUTE computes a Gauss-Legendre rule;\n";
  cout << "\n";
  cout << "  Compute the data for N = " << n << "\n";

  x = new double[n];
  w = new double[n];

  legendre_ek_compute ( n, x, w );
 
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "    x[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << x[i] << ";\n";
  }
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "    w[" << setw(2) << i
         << "] = " << setw(24) << setprecision(16) << w[i] << ";\n";
  }

  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests LEGENDRE_EK_COMPUTE and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int iexp;
  int n = 2;
  int nsub;
  double result;
  double *w;
  double *x;
  double xhi;
  double xlo;

  a = 0.0;
  b = 1.0;

  xlo = -1.0;
  xhi = +1.0;

  x = new double[n];
  w = new double[n];

  cout << "\n";
  cout << "TEST19\n";
  cout << "  LEGENDRE_EK_COMPUTE computes a Gauss-Legendre rule;\n";
  cout << "  SUM_SUB carries it out over subintervals.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "]\n";
  cout << "  Here, we use a fixed order N = " << n << "\n";
  cout << "  and use more and more subintervals.\n";
  cout << "\n";
  cout << "  NSUB     Integral\n";
  cout << "\n";
 
  legendre_ek_compute ( n, x, w );
 
  for ( iexp = 0; iexp <= 9; iexp++ )
  {
    nsub = i4_power ( 2, iexp );

    result = sum_sub ( fx2sd1, a, b, nsub, n, xlo, xhi, x, w );

    cout << "  " << setw(4)  << nsub
         << "  " << setw(16) << setprecision(8) << result << "\n";
  } 

  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests LEGENDRE_EK_COMPUTE and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 10;
  double *result;
  double *weight;
  double *x;
  double xhi;
  double xlo;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = -1.0;
  xhi = +1.0;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  LEGENDRE_EK_COMPUTE sets up a Gauss-Legendre rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [ " << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << ".\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      x = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_ek_compute ( order, x, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          x, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] x;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests LEGENDRE_SET and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = -1.0;
  xhi = +1.0;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  LEGENDRE_SET sets up a Gauss-Legendre rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [ " << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << ".\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests LEGENDRE_SET, LEGENDRE_SET_X0_01 and RULE_ADJUST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ORDER 5

  double a;
  double b;
  double c;
  double d;
  int i;
  double weight1[ORDER];
  double weight2[ORDER];
  double weight3[ORDER];
  double xtab1[ORDER];
  double xtab2[ORDER];
  double xtab3[ORDER];

  a = -1.0;
  b = +1.0;
  c =  0.0;
  d =  1.0;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  LEGENDRE_SET sets up a Gauss-Legendre rule\n";
  cout << "    for integrating F(X) over [-1,1];\n";
  cout << "  RULE_ADJUST adjusts a rule for a new interval.\n";
  cout << "  LEGENDRE_SET_X0_01 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating F(X) over [0,1];\n";
  cout << "\n";
  cout << "  We will use LEGENDRE_SET to get a rule for [-1,1],\n";
  cout << "  adjust it to [0,1] using RULE_ADJUST,\n";
  cout << "  and compare the results of LEGENDRE_SET_X0_01.\n";
  cout << "\n";

  legendre_set ( ORDER, xtab1, weight1 );

  r8vec_copy ( ORDER, xtab1, xtab2 );
  r8vec_copy ( ORDER, weight1, weight2 );

  rule_adjust ( a, b, c, d, ORDER, xtab2, weight2 );

  legendre_set_x0_01 ( ORDER, xtab3, weight3 );

  cout << "\n";
  cout << "  Abscissas:\n";
  cout << "\n";
  cout << "          Original          Adjusted            Stored\n";
  cout << "\n";

  for ( i = 0; i < ORDER; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(16) << setprecision(12) << xtab1[i]
         << "  " << setw(16) << setprecision(12) << xtab2[i]
         << "  " << setw(16) << setprecision(12) << xtab3[i] << "\n";
  }

  cout << "\n";
  cout << "  Weights:\n";
  cout << "\n";
  cout << "          Original          Adjusted            Stored\n";
  cout << "\n";

  for ( i = 0; i < ORDER; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(16) << setprecision(12) << weight1[i]
         << "  " << setw(16) << setprecision(12) << weight2[i]
         << "  " << setw(16) << setprecision(12) << weight3[i] << "\n";
  }

  return;
# undef ORDER
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests LEGENDRE_SET_COS and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int iexp;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 20;
  double pi = 3.141592653589793;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = -0.5 * pi;
  b = +0.5 * pi;

  nsub = 1;

  xlo = -0.5 * pi;
  xhi = +0.5 * pi;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  LEGENDRE_SET_COS sets up a Gauss-Legendre rule\n";
  cout << "    over [-PI/2,PI/2] with weight function COS(X);\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( iexp = 0; iexp <= 4; iexp++ )
    {
      order = i4_power ( 2, iexp );

      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_cos ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests LEGENDRE_SET_SQRTX_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int iexp;
  int ihi;
  int ilo;
  int j;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  LEGENDRE_SET_SQRTX_01 sets up a Gauss-Legendre rule\n";
  cout << "  over [0,1] with weight function SQRT(X);\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( iexp = 0; iexp <= 3; iexp++ )
    {
      order = i4_power ( 2, iexp );

      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_sqrtx_01 ( order, xtab, weight );
 
        result[i] = 0.0;
        for ( j = 0; j < order; j++ )
        {
          result[i] = result[i] + weight[j] * function_value ( xtab[j] );
        }
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests LEGENDRE_SET_SQRTX2_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int iexp;
  int ihi;
  int ilo;
  int j;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  cout << "\n";
  cout << "TEST25\n";
  cout << "  LEGENDRE_SET_SQRTX2_01 sets up a Gauss-Legendre rule\n";
  cout << "  over [0,1] with weight function 1/SQRT(X);\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( iexp = 0; iexp <= 3; iexp++ )
    {
      order = i4_power ( 2, iexp );

      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_sqrtx2_01 ( order, xtab, weight );
 
        result[i] = 0.0;
        for ( j = 0; j < order; j++ )
        {
          result[i] = result[i] + weight[j] * function_value ( xtab[j] );
        }
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests LEGENDRE_SET_COS2 and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int iexp;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 20;
  double pi = 3.141592653589793;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = -0.5 * pi;
  b = +0.5 * pi;

  nsub = 1;

  xlo = -0.5 * pi;
  xhi = +0.5 * pi;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  LEGENDRE_SET_COS2 sets up a Gauss-Legendre rule\n";
  cout << "    over [0,PI/2] with weight function COS(X);\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( iexp = 1; iexp <= 4; iexp++ )
    {
      order = i4_power ( 2, iexp );

      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_cos2 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests LEGENDRE_SET_LOG and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 9

  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 20;
  int order_test[TEST_NUM] = { 1, 2, 3, 4, 5, 6, 7, 8, 16 };
  double *result;
  int test;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = 0.0;
  xhi = +1.0;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  LEGENDRE_SET_LOG sets up a Gauss-Legendre rule\n";
  cout << "    for integrating -LOG(X) * F(X) over [0,1];\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [ " << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << ".\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( test = 0; test < TEST_NUM; test++ )
    {
      order = order_test[test];
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_log ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests LEGENDRE_SET_X0_01 and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 8;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = 0.0;
  xhi = +1.0;

  cout << "\n";
  cout << "TEST28\n";
  cout << "  LEGENDRE_SET_X0_01 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating F(X) over [0,1];\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [ " << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << ".\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x0_01 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests LEGENDRE_SET_X1 and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = -1.0;
  xhi = +1.0;

  cout << "\n";
  cout << "TEST29\n";
  cout << "  LEGENDRE_SET_X1 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating ( 1 + X ) * F(X) over [-1,1];\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [ " << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << ".\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x1 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests LEGENDRE_SET_X1, LEGENDRE_SET_X1_01 and RULE_ADJUST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ORDER 5

  double a;
  double b;
  double c;
  double d;
  int i;
  double weight1[ORDER];
  double weight2[ORDER];
  double weight3[ORDER];
  double xtab1[ORDER];
  double xtab2[ORDER];
  double xtab3[ORDER];

  a = -1.0;
  b =  1.0;
  c =  0.0;
  d =  1.0;

  cout << "\n";
  cout << "TEST30:\n";
  cout << "  LEGENDRE_SET_X1 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating ( 1 + X ) * F(X) over [-1,1];\n";
  cout << "  RULE_ADJUST adjusts a rule for a new interval.\n";
  cout << "  LEGENDRE_SET_X1_01 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating X * F(X) over [0,1];\n";
  cout << "\n";
  cout << "  We will use LEGENDRE_SET_X1 to get a rule for [-1,1],\n";
  cout << "  adjust it to [0,1] using RULE_ADJUST,\n";
  cout << "  make further adjustments because the weight function\n";
  cout << "  is not 1,\n";
  cout << "  and compare the results of LEGENDRE_SET_X1_01.\n";
  cout << "\n";

  legendre_set_x1 ( ORDER, xtab1, weight1 );

  for ( i = 0; i < ORDER; i++ )
  {
    xtab2[i] = xtab1[i];
  }
  for ( i = 0; i < ORDER; i++ )
  {
    weight2[i] = weight1[i];
  }

  rule_adjust ( a, b, c, d, ORDER, xtab2, weight2 );
//
//  Because the weight function W(X) is not 1, we need to do more
//  adjustments to the weight vector.
//
  for ( i = 0; i < ORDER; i++ )
  {
    weight2[i] = weight2[i] / 2.0;
  }
  legendre_set_x1_01 ( ORDER, xtab3, weight3 );

  cout << "\n";
  cout << "  Abscissas:\n";
  cout << "\n";
  cout << "  Original  Adjusted Stored\n";
  cout << "\n";

  for ( i = 0; i < ORDER; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(16) << setprecision(12) << xtab1[i]
         << "  " << setw(16) << setprecision(12) << xtab2[i]
         << "  " << setw(16) << setprecision(12) << xtab3[i] << "\n";
  }

  cout << "\n";
  cout << "  Weights:\n";
  cout << "\n";
  cout << "  Original  Adjusted Stored\n";
  cout << "\n";

  for ( i = 0; i < ORDER; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(16) << setprecision(12) << weight1[i]
         << "  " << setw(16) << setprecision(12) << weight2[i]
         << "  " << setw(16) << setprecision(12) << weight3[i] << "\n";
  }

  return;
# undef ORDER
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests LEGENDRE_SET_X1_01 and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 8;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = 0.0;
  xhi = +1.0;

  cout << "\n";
  cout << "TEST31\n";
  cout << "  LEGENDRE_SET_X1_01 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating X * F(X) over [0,1];\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [ " << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << ".\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x1_01 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests LEGENDRE_SET_Xs and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = -1.0;
  xhi = +1.0;

  cout << "\n";
  cout << "TEST32\n";
  cout << "  LEGENDRE_SET_X2 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating (1 + X)**2 * F(X) over [-1,1];\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [ " << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << ".\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x2 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests LEGENDRE_SET_X2, LEGENDRE_SET_X2_01 and RULE_ADJUST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ORDER 5

  double a;
  double b;
  double c;
  double d;
  int i;
  double weight1[ORDER];
  double weight2[ORDER];
  double weight3[ORDER];
  double xtab1[ORDER];
  double xtab2[ORDER];
  double xtab3[ORDER];

  a = -1.0;
  b =  1.0;
  c =  0.0;
  d =  1.0;

  cout << "\n";
  cout << "TEST33:\n";
  cout << "  LEGENDRE_SET_X2 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating ( 1 + X )^2 * F(X) over [-1,1];\n";
  cout << "  RULE_ADJUST adjusts a rule for a new interval.\n";
  cout << "  LEGENDRE_SET_X2_01 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating X^2 * F(X) over [0,1];\n";
  cout << "\n";
  cout << "  We will use LEGENDRE_SET_X2 to get a rule for [-1,1],\n";
  cout << "  adjust it to [0,1] using RULE_ADJUST,\n";
  cout << "  make further adjustments because the weight function\n";
  cout << "  is not 1,\n";
  cout << "  and compare the results of LEGENDRE_SET_X2_01.\n";
  cout << "\n";

  legendre_set_x2 ( ORDER, xtab1, weight1 );

  for ( i = 0; i < ORDER; i++ )
  {
    xtab2[i] = xtab1[i];
  }
  for ( i = 0; i < ORDER; i++ )
  {
    weight2[i] = weight1[i];
  }

  rule_adjust ( a, b, c, d, ORDER, xtab2, weight2 );
//
//  Because the weight function W(X) is not 1, we need to do more
//  adjustments to the weight vector.
//
  for ( i = 0; i < ORDER; i++ )
  {
    weight2[i] = weight2[i] / 4.0;
  }
  legendre_set_x2_01 ( ORDER, xtab3, weight3 );

  cout << "\n";
  cout << "  Abscissas:\n";
  cout << "\n";
  cout << "  Original  Adjusted Stored\n";
  cout << "\n";

  for ( i = 0; i < ORDER; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(16) << setprecision(12) << xtab1[i]
         << "  " << setw(16) << setprecision(12) << xtab2[i]
         << "  " << setw(16) << setprecision(12) << xtab3[i] << "\n";
  }

  cout << "\n";
  cout << "  Weights:\n";
  cout << "\n";
  cout << "  Original  Adjusted Stored\n";
  cout << "\n";

  for ( i = 0; i < ORDER; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(16) << setprecision(12) << weight1[i]
         << "  " << setw(16) << setprecision(12) << weight2[i]
         << "  " << setw(16) << setprecision(12) << weight3[i] << "\n";
  }

  return;
# undef ORDER
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 tests LEGENDRE_SET_X2_01 and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 8;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = 0.0;
  xhi = +1.0;

  cout << "\n";
  cout << "TEST34\n";
  cout << "  LEGENDRE_SET_X2_01 sets up a Gauss-Legendre rule\n";
  cout << "    for integrating X*X * F(X) over [0,1];\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [ " << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << ".\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x2_01 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }
  delete [] result;

  return;
}
//****************************************************************************80

void lobatto_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_COMPUTE_TEST tests LOBATTO_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  cout << "\n";
  cout << "LOBATTO_COMPUTE_TEST\n";
  cout << "  LOBATTO_COMPUTE computes a Lobatto rule;\n";
  cout << "\n";
  cout << "         I      X             W\n";

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w = new double[n];
    x = new double[n];

    lobatto_compute ( n, x, w );

    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8)  << i
           << "  " << setw(12) << x[i]
           << "  " << setw(12) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void lobatto_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_SET_TEST tests LOBATTO_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  cout << "\n";
  cout << "LOBATTO_SET_TEST\n";
  cout << "  LOBATTO_SET sets a Lobatto rule;\n";
  cout << "\n";
  cout << "         I      X             W\n";

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w = new double[n];
    x = new double[n];

    lobatto_set ( n, x, w );

    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8)  << i
           << "  " << setw(12) << x[i]
           << "  " << setw(12) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void moulton_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MOULTON_SET_TEST tests MOULTON_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  cout << "\n";
  cout << "MOULTON_SET_TEST\n";
  cout << "  MOULTON_SET sets up an Adams-Moulton rule;\n";
  cout << "\n";
  cout << "  Index             X                   W\n";
  cout << "\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    moulton_set ( n, x, w );

    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i 
           << "  " << setw(24) << x[i]
           << "  " << setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
 
  return;
}
//****************************************************************************80

void ncc_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_SET_TEST tests NCC_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  cout << "\n";
  cout << "NCC_SET_TEST\n";
  cout << "  NCC_SET sets up a Newton-Cotes Closed rule;\n";
  cout << "\n";
  cout << "  Index             X                   W\n";
  cout << "\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    ncc_set ( n, x, w );

    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i 
           << "  " << setw(24) << x[i]
           << "  " << setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
 
  return;
}
//****************************************************************************80

void test38 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST38 tests NCC_COMPUTE and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 21;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST38\n";
  cout << "  NCC_COMPUTE computes a closed Newton-Cotes rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        ncc_compute ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test39 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST39 tests NCO_SET and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST39\n";
  cout << "  NCO_SET sets up an open Newton-Cotes rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      if ( order == 8 )
      {
        continue;
      }
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        nco_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test40 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST40 tests NCO_COMPUTE and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST40\n";
  cout << "  NCO_COMPUTE computes an open Newton-Cotes rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        nco_compute ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test401 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST401 tests NCOH_SET and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST401\n";
  cout << "  NCOH_SET sets up an open half Newton-Cotes rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        ncoh_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test402 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST402 tests NCOH_COMPUTE and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST402\n";
  cout << "  NCOH_COMPUTE computes an open half Newton-Cotes rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        ncoh_compute ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test403 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST403 tests PATTERSON_SET and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int level;
  int level_max = 7;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST403\n";
  cout << "  PATTERSON_SET sets up a Patterson rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order   ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( level = 1; level <= level_max; level++ )
    {
      order = i4_power ( 2, level ) - 1;

      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        patterson_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(3) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

void test404 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST404 tests RADAU_COMPUTE and RADAU_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST404\n";
  cout << "  RADAU_COMPUTE computes a Radau rule;\n";
  cout << "  RADAU_SET sets a rule from a table.\n";
  cout << "\n";
  cout << "         I      X1            X2            W1            W2\n";

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w1 = new double[n];
    w2 = new double[n];
    x1 = new double[n];
    x2 = new double[n];

    radau_compute ( n, x1, w1 );
    radau_set ( n, x2, w2 );

    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8)  << i+1
           << "  " << setw(12) << x1[i]
           << "  " << setw(12) << x2[i]
           << "  " << setw(12) << w1[i]
           << "  " << setw(12) << w2[i] << "\n";
    }
    delete [] w1;
    delete [] w2;
    delete [] x1;
    delete [] x2;
  }
  return;
}
//****************************************************************************80

void test41 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST41 tests RADAU_SET and SUM_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 15;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = new double[function_num];

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  cout << "\n";
  cout << "TEST41\n";
  cout << "  RADAU_SET sets up a Radau rule;\n";
  cout << "  SUM_SUB carries it out.\n";
  cout << "\n";
  cout << "  The integration interval is [" << a << "," << b << "].\n";
  cout << "  The number of subintervals is " << nsub << "\n";
  cout << "  Quadrature order will vary.\n";
  cout << "  Integrand will vary.\n";
  cout << "\n";

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    cout << "\n";
    cout << "Order  ";
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(10) << function_name ( i ) << "    ";
    }
    cout << "\n";
    cout << "\n";

    for ( order = 1; order <= order_max; order++ )
    {
      if ( order == 8 )
      {
        continue;
      }

      xtab = new double[order];
      weight = new double[order];

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        radau_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      cout << setw(2) << order << "  ";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << "  " << setw(12) << setprecision(8) << result[i];
      }
      cout << "\n";

      delete [] xtab;
      delete [] weight;
    }
  }

  delete [] result;
 
  return;
}
//****************************************************************************80

double f1sd1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F1SD1 evaluates the function 1.0D+00/ sqrt ( 1.1 - x**2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double F1SD1, the value of the function.
//
{
  double value;

  value = 1.0 / sqrt ( 1.1 - x * x );
 
  return value;
}
//****************************************************************************80

char *function_name ( int function_index )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_NAME returns the name of the function evaluated in FUNCTION_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer FUNCTION_INDEX, the index of the function.
//
//    Output, char *FUNCTION_NAME, the name of the function.
//
{
  char *value;

  value = new char[10];

  if ( function_index == 0 )
  {
    strcpy ( value, "         1" );
  }
  else if ( function_index == 1 )
  {
    strcpy ( value, "         X" );
  }
  else if ( function_index == 2 )
  {
    strcpy ( value, "       X^2" );
  }
  else if ( function_index == 3 )
  {
    strcpy ( value, "       X^3" );
  }
  else if ( function_index == 4 )
  {
    strcpy ( value, "       X^4" );
  }
  else if ( function_index == 5 )
  {
    strcpy ( value, "       X^5" );
  }
  else if ( function_index == 6 )
  {
    strcpy ( value, "       X^6" );
  }
  else if ( function_index == 7 )
  {
    strcpy ( value, "       X^7" );
  }
  else if ( function_index == 8 )
  {
    strcpy ( value, "    SIN(X)" );
  }
  else if ( function_index == 9 )
  {
    strcpy ( value, "    EXP(X)" );
  }
  else if ( function_index == 10 )
  {
    strcpy ( value, " SQRT(|X|)" );
  }
  else
  {
    strcpy ( value, "??????????" );
  }
  return value;
}
//****************************************************************************80

void function_set ( string action, int *i )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_SET sets the function to be returned by FUNCTION_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION, the action to be carried out.
//    "COUNT" means the call is made to count the number of functions available.
//    "GET" means the call is made to find out the current function index.
//    "SET" means the call is made to set the current function index.
//
//    Input/output, int *I.
//    For "COUNT", I is output as the number of functions available;
//    For "GET", I is output as the currently chosen function;
//    For "SET", I is input as the user's new choice for the function.
//
{
  static int function_index = -1;

  if ( s_eqi ( action, "COUNT" ) )
  {
    *i = 11;
  }
  else if ( s_eqi ( action, "GET" ) )
  {
    *i = function_index;
  }
  else if ( s_eqi ( action, "SET" ) )
  {
    function_index = *i;
  }
  else
  {
    cout << "\n";
    cout << "FUNCTION_SET - Warning!\n";
    cout << "  Unrecognized action = \"" << action << "\".\n";
  }

  return;
}
//****************************************************************************80

double function_value ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_VALUE evaluates a function of X, as chosen by the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double FUNCTION_VALUE, the value of the function.
//
{
  int function_index;
  double value;

  function_set ( "GET", &function_index );

  if ( function_index == 0 )
  {
    value = 1.0;
  }
  else if ( function_index == 1 )
  {
    value = x;
  }
  else if ( function_index == 2 )
  {
    value = pow ( x, 2 );
  }
  else if ( function_index == 3 )
  {
    value = pow ( x, 3 );
  }
  else if ( function_index == 4 )
  {
    value = pow ( x, 4 );
  }
  else if ( function_index == 5 )
  {
    value = pow ( x, 5 );
  }
  else if ( function_index == 6 )
  {
    value = pow ( x, 6 );
  }
  else if ( function_index == 7 )
  {
    value = pow ( x, 7 );
  }
  else if ( function_index == 8 )
  {
    value = sin ( x );
  }
  else if ( function_index == 9 )
  {
    value = exp ( x );
  }
  else if ( function_index == 10 )
  {
    value = sqrt ( fabs ( x ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double fxsd1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FXSD1 evaluates the function x / sqrt ( 1.1 - x**2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double FXSD1, the value of the function.
//
{
  double value;

  value = x / sqrt ( 1.1 - x * x );
 
  return value;
}
//****************************************************************************80

double fx2sd1 ( double x )

//****************************************************************************80
//
//  Purpose;
//
//    FX2SD1 evaluates the function x**2 / sqrt ( 1.1 - x**2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double FX2SD1, the value of the function.
//
{
  double value;

  value = x * x / sqrt ( 1.1 - x * x );
 
  return value;
}
