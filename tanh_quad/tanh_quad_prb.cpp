# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "tanh_quad.hpp"

int main ( );
void test01 ( );
void test012 ( );
void test02 ( );
void test025 ( );
void test03 ( );
void test032 ( int p );
void test04 ( int p );
void test05 ( int p );
void test06 ( int p );
void p00_ab ( int p, double *a, double *b );
void p00_e ( int p, double *exact );
void p00_f ( int p, int n, double x[], double fx[] );
void p01_ab ( double *a, double *b );
void p01_e ( double *exact );
void p01_f ( int n, double x[], double fx[] );
void p02_ab ( double *a, double *b );
void p02_e ( double *exact );
void p02_f ( int n, double x[], double fx[] );
void p03_ab ( double *a, double *b );
void p03_e ( double *exact );
void p03_f ( int n, double x[], double fx[] );
void p04_ab ( double *a, double *b );
void p04_e ( double *exact );
void p04_f ( int n, double x[], double fx[] );
void p05_ab ( double *a, double *b );
void p05_e ( double *exact );
void p05_f ( int n, double x[], double fx[] );
void p06_ab ( double *a, double *b );
void p06_e ( double *exact );
void p06_f ( int n, double x[], double fx[] );
void p07_ab ( double *a, double *b );
void p07_e ( double *exact );
void p07_f ( int n, double x[], double fx[] );
void p08_ab ( double *a, double *b );
void p08_e ( double *exact );
void p08_f ( int n, double x[], double fx[] );
void p09_ab ( double *a, double *b );
void p09_e ( double *exact );
void p09_f ( int n, double x[], double fx[] );
void p10_ab ( double *a, double *b );
void p10_e ( double *exact );
void p10_f ( int n, double x[], double fx[] );
void p11_ab ( double *a, double *b );
void p11_e ( double *exact );
void p11_f ( int n, double x[], double fx[] );
void p12_ab ( double *a, double *b );
void p12_e ( double *exact );
void p12_f ( int n, double x[], double fx[] );
void p13_ab ( double *a, double *b );
void p13_e ( double *exact );
void p13_f ( int n, double x[], double fx[] );
void p14_ab ( double *a, double *b );
void p14_e ( double *exact );
void p14_f ( int n, double x[], double fx[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TANH_QUAD_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int p;

  timestamp ( );

  cout << "\n";
  cout << "TANH_QUAD_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TANH_QUAD library.\n";

  test01 ( );
  test012 ( );
  test02 ( );
  test025 ( );
  test03 ( );

  for ( p = 1; p <= 14; p++ )
  {
    test032 ( p );
  }

  for ( p = 1; p <= 14; p++ )
  {
    test04 ( p );
  }

  for ( p = 1; p <= 14; p++ )
  {
    test05 ( p );
  }

  for ( p = 1; p <= 14; p++ )
  {
    test06 ( p );
  }
  cout << "\n";
  cout << "TANH_QUAD_PRB\n";
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
//    TEST01 demonstrates TANH_M_TO_H and TANH_H_TO_N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  double h;
  int m;
  int n;
  double tol;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  TANH_M_TO_H determines the spacing H from level M\n";
  cout << "  TANH_H_TO_N determines the quadrature order N from\n";
  cout << "  the spacing H and tolerance TOL.\n";

  tol = r8_epsilon ( );

  cout << "\n";
  cout << "  All tests use TOL = " << tol << "\n";
  cout << "\n";
  cout << "         M        H                  N\n";
  cout << "\n";

  for ( m = -2; m <= 8; m++ )
  {
    h = tanh_m_to_h ( m );

    n = tanh_h_to_n ( h, tol );

    cout << "  " << setw(8)  << m
         << "  " << setw(16) << h
         << "  " << setw(8)  << n << "\n";
  }

  return;
}
//****************************************************************************80

void test012 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST012 computes a midpoint quadrature rule W, X with N = 0, 1, 3, 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  double h;
  int i;
  int m;
  int n;
  int offset;
  int order;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST012\n";
  cout << "  Determine nested midpoint quadrature rules W, X\n";
  cout << "  by choosing N = 0, 1, 3, 7.\n";

  for ( m = 1; m <= 4; m++ )
  {
    order = i4_power ( 2, m ) - 1;
    n = ( ( order + 1 ) / 2 ) - 1;
    offset = n;
    h = 2.0 / ( double ) ( order + 1 );

    cout << "\n";
    cout << "  M = " << m
         << "  ORDER = " << order
         << "  N = " << n
         << "  H = " << h << "\n";

    w = new double[order];
    x = new double[order];

    midpoint_rule ( n, x, w );

    cout << "\n";
    cout << "         I      Wi                Xi\n";
    cout << "\n";

    for ( i = -n; i <= n; i++ )
    {
      cout << "  " << setw(8)  << i
           << "  " << setw(16) << w[i+offset]
           << "  " << setw(16) << x[i+offset] << "\n";
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
//    TEST02 computes a tanh quadrature rule W, X with N = 5.
//
//  Discussion:
//
//    These results should match the following values reported in 
//    Kahaner, Moler and Nash:
//
//       I    Wi         Xi
//      --  --------  --------
//      -5  0.000471  -0.999737
//      -4  0.002807  -0.998428
//      -3  0.016631  -0.990649
//      -2  0.094844  -0.945434
//      -1  0.439127  -0.713098
//       0  0.893459   0.000000
//       1  0.439127   0.713098
//       2  0.094844   0.945434
//       3  0.016631   0.990649
//       4  0.002807   0.998428
//       5  0.000471   0.999737
//
//    Note that these values do not sum to 2, although they come close//
//    Thus, a fundamental feature of most quadrature rules is ignored here.
//    This rule will not integrate f(x) = 1 exactly.  But it is not a
//    family of rules based on polynomial accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
{
  double h;
  int i;
  int j;
  int m;
  int n;
  int offset;
  int order;
  double tol;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Determine tanh quadrature rules W, X\n";
  cout << "  by choosing N = 5, 10, 20.\n";
  cout << "\n";
  cout << "  The rule for N = 5 appears in the reference\n";
  cout << "  Kahaner, Moler and Nash.\n";
  cout << "\n";
  cout << "  Also note that the rule for H(N) means that\n";
  cout << "  rules for doubled N do not nest.\n";

  tol = r8_epsilon ( );

  n = 5;

  for ( j = 1; j <= 3; j++ )
  {
    cout << "\n";
    cout << "  Quadrature order N = " << n << "\n";

    h = tanh_n_to_h ( n );

    offset = n;
    order = 2 * n + 1;

    w = new double[order];
    x = new double[order];

    tanh_rule ( n, h, x, w );

    cout << "\n";
    cout << "         I      Wi                Xi\n";
    cout << "\n";

    for ( i = -n; i <= n; i++ )
    {
      cout << "  " << setw(8)  << i
           << "  " << setw(16) << w[i+offset]
           << "  " << setw(16) << x[i+offset] << "\n";
    }

    delete [] w;
    delete [] x;

    n = 2 * n;
  }

  cout << "\n";
  cout << "  Note that the weights do not sum to 2!\n";
  cout << "\n";
  cout << "         N        sum(W)\n";
  cout << "\n";

  n = 5;

  for ( j = 1; j <= 5; j++ )
  {
    h = tanh_n_to_h ( n );

    offset = n;
    order = 2 * n + 1;

    w = new double[order];
    x = new double[order];

    tanh_rule ( n, h, x, w );

    cout << "  " << setw(8) << n
         << "  " << setw(16) << r8vec_sum ( order, w ) << "\n";

    delete [] w;
    delete [] x;

    n = 2 * n;
  }

  return;
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 computes tanh-sinh quadrature rules.
//
//  Discussion:
//
//    We are seeking a family of nested rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
{
  double h;
  int i;
  int j;
  int m;
  int n;
  int offset;
  int order;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST025\n";
  cout << "  Determine tanh-sinh quadrature rules W, X\n";
  cout << "  for N = 1, 3, 7, 15, 31, 63.\n";

  for ( m = -3; m <= 3; m++ )
  {
    order = i4_power ( 2, m + 4 ) - 1;
    n = ( order + 1 ) / 2 - 1;
    h = 4.0 / ( double ) ( order + 1 );
    offset = n;

    cout << "\n";
    cout << "  M = " << m
         << "  ORDER = " << order
         << "  N = " << n
         << "  H = " << h << "\n";

    w = new double[order];
    x = new double[order];

    tanh_sinh_rule ( n, h, x, w );

    cout << "\n";
    cout << "         I      Wi                Xi\n";
    cout << "\n";

    for ( i = -n; i <= n; i++ )
    {
      cout << "  " << setw(8)  << i
           << "  " << setw(16) << w[i+offset]
           << "  " << setw(16) << x[i+offset] << "\n";
    }

    delete [] w;
    delete [] x;
  }

  cout << "\n";
  cout << "  Note that, especially for low N, the weights need not sum to 2!\n";
  cout << "\n";
  cout << "         N     H            sum(W)\n";
  cout << "\n";

  for ( m = -3; m <= 10; m++ )
  {
    order = i4_power ( 2, m + 4 ) - 1;
    n = ( order + 1 ) / 2 - 1;
    h = 1.0 / pow ( 2.0, m );

    w = new double[order];
    x = new double[order];

    tanh_sinh_rule ( n, h, x, w );

    cout << "  " << setw(8) << n
         << "  " << setw(16) << r8vec_sum ( order, w ) << "\n";

    delete [] w;
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
//    TEST03 computes a quadrature rule W, X based on a tolerance.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  double h;
  int i;
  int m;
  int n;
  int offset;
  int order;
  double tol;
  double *w;
  double *x;
  double pi = 3.141592653589793;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Determine a quadrature rule W, X by specifying\n";
  cout << "  a tolerance.\n";

  tol = r8_epsilon ( );

  cout << "  Tolerance TOL = " << tol << "\n";

  for ( m = -1; m <= 2; m++ )
  {
    cout << "\n";
    cout << "  Level M = " << m << "\n";
    h = tanh_m_to_h ( m );
    cout << "  Spacing H = " << h << "\n";
    n = tanh_h_to_n ( h, tol );
    cout << "  Quadrature order N = " << n << "\n";

    order = 2 * n + 1;
    offset = n;

    w = new double[order];
    x = new double[order];

    tanh_rule ( n, h, x, w );

    cout << "\n";
    cout << "         I      Wi                Xi\n";
    cout << "\n";

    for ( i = -n; i <= n; i++ )
    {
      cout << "  " << setw(8)  << i
           << "  " << setw(16) << w[i+offset]
           << "  " << setw(16) << x[i+offset] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test032 ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    TEST032 applies a sequence of midpoint rules to a test integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the problem index.
//
{
  double a;
  double b;
  double c;
  double d;
  double error;
  double exact;
  double *fx;
  double h;
  double h_midpoint;
  int i;
  int m;
  int n;
  int offset;
  int order;
  double quad;
  double tol;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST032\n";
  cout << "  Apply a sequence of midpoint rules to a test integral.\n";

  tol = r8_epsilon ( );

  cout << "\n";
  cout << "  Problem index P = " << p << "\n";
  cout << "  Tolerance TOL = " << tol << "\n";

  cout << "\n";
  cout << "         M      H_MIDPOINT           N          Exact           Quad            Error\n";
  cout << "\n";

  for ( m = 1; m <= 8; m++ )
  {
//
//  Choose N based on the value of H that would be used by the TANH rule.
//
    h = tanh_m_to_h ( m );
    n = tanh_h_to_n ( h, tol );
    offset = n;
    order = 2 * n + 1;
    h_midpoint = 2.0 / ( double ) ( 2 * ( n + 1 ) );

    w = new double[order];
    x = new double[order];

    midpoint_rule ( n, x, w );
//
//  Adjust the rule from [-1,+1] to the actual integration limits.
//
    a = -1.0;
    b = +1.0;

    p00_ab ( p, &c, &d );

    rule_adjust ( a, b, c, d, order, x, w );
//
//  Evaluate the integrand.
//
    fx = new double[order];

    p00_f ( p, order, x, fx );
//
//  Form the quadrature estimate.
//
    quad = r8vec_dot ( order, w, fx );

    p00_e ( p, &exact );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8)  << m
         << "  " << setw(16) << h_midpoint
         << "  " << setw(8)  << n
         << "  " << setw(16) << exact
         << "  " << setw(16) << quad
         << "  " << setw(16) << error << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test04 ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 applies a sequence of trapezoid rules to a test integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the problem index.
//
{
  double a;
  double b;
  double c;
  double d;
  double error;
  double exact;
  double *fx;
  double h;
  int i;
  int m;
  int n;
  int offset;
  int order;
  double quad;
  double tol;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Apply a sequence of trapezoid rules to a test integral.\n";

  tol = r8_epsilon ( );

  cout << "\n";
  cout << "  Problem index P = " << p << "\n";
  cout << "  Tolerance TOL = " << tol << "\n";

  cout << "\n";
  cout << "         M      H                    N          Exact           Quad            Error\n";
  cout << "\n";

  for ( m = 1; m <= 8; m++ )
  {
    h = tanh_m_to_h ( m );

    n = tanh_h_to_n ( h, tol );

    order = 2 * n + 1;

    w = new double[order];
    x = new double[order];

    trap_rule ( n, x, w );
//
//  Adjust the rule from [-1,+1] to the actual integration limits.
//
    a = -1.0;
    b = +1.0;

    p00_ab ( p, &c, &d );

    rule_adjust ( a, b, c, d, order, x, w );
//
//  Evaluate the integrand.
//
    fx = new double[order];

    p00_f ( p, order, x, fx );
//
//  Form the quadrature estimate.
//
    quad = r8vec_dot ( order, w, fx );

    p00_e ( p, &exact );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8)  << m
         << "  " << setw(16) << h
         << "  " << setw(8)  << n
         << "  " << setw(16) << exact
         << "  " << setw(16) << quad
         << "  " << setw(16) << error << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test05 ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 applies a sequence of tanh rules to a test integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the problem index.
//
{
  double a;
  double b;
  double c;
  double d;
  double error;
  double exact;
  double *fx;
  double h;
  int i;
  int m;
  int n;
  int offset;
  int order;
  double quad;
  double tol;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Apply a sequence of tanh rules to a test integral.\n";

  tol = r8_epsilon ( );

  cout << "\n";
  cout << "  Problem index P = " << p << "\n";
  cout << "  Tolerance TOL = " << tol << "\n";

  cout << "\n";
  cout << "         M      H                    N        Exact           Quad            Error\n";
  cout << "\n";

  for ( m = -2; m <= 8; m++ )
  {
    h = tanh_m_to_h ( m );

    n = tanh_h_to_n ( h, tol );

    order = 2 * n + 1;

    w = new double[order];
    x = new double[order];

    tanh_rule ( n, h, x, w );
//
//  Adjust the rule from [-1,+1] to the actual integration limits.
//
    a = -1.0;
    b = +1.0;

    p00_ab ( p, &c, &d );

    rule_adjust ( a, b, c, d, order, x, w );
//
//  Evaluate the integrand.
//
    fx = new double[order];

    p00_f ( p, order, x, fx );
//
//  Form the quadrature estimate.
//
    quad = r8vec_dot ( order, w, fx );

    p00_e ( p, &exact );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8)  << m
         << "  " << setw(16) << h
         << "  " << setw(8)  << n
         << "  " << setw(16) << exact
         << "  " << setw(16) << quad
         << "  " << setw(16) << error << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test06 ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 applies a sequence of tanh-sinh rules to a test integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the problem index.
//
{
  double a;
  double b;
  double c;
  double d;
  double error;
  double exact;
  double *fx;
  double h;
  int i;
  int m;
  int n;
  int offset;
  int order;
  double quad;
  double tol;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Apply a sequence of tanh-sinh rules to a test integral.\n";

  tol = r8_epsilon ( );

  cout << "\n";
  cout << "  Problem index P = " << p << "\n";
  cout << "  Tolerance TOL = " << tol << "\n";

  cout << "\n";
  cout << "         M      H                    N        Exact           Quad            Error\n";
  cout << "\n";

  for ( m = -2; m <= 8; m++ )
  {
    h = tanh_m_to_h ( m );

    n = tanh_sinh_h_to_n ( h, tol );

    order = 2 * n + 1;

    w = new double[order];
    x = new double[order];

    tanh_sinh_rule ( n, h, x, w );
//
//  Adjust the rule from [-1,+1] to the actual integration limits.
//
    a = -1.0;
    b = +1.0;

    p00_ab ( p, &c, &d );

    rule_adjust ( a, b, c, d, order, x, w );
//
//  Evaluate the integrand.
//
    fx = new double[order];

    p00_f ( p, order, x, fx );
//
//  Form the quadrature estimate.
//
    quad = r8vec_dot ( order, w, fx );

    p00_e ( p, &exact );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8)  << m
         << "  " << setw(16) << h
         << "  " << setw(8)  << n
         << "  " << setw(16) << exact
         << "  " << setw(16) << quad
         << "  " << setw(16) << error << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void p00_ab ( int p, double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P00_AB returns the integration limits for a given problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the problem index.
//
//    Output, double *A, *B, the integration limits.
//
{
  if ( p == 1 )
  {
    p01_ab ( a, b );
  }
  else if ( p == 2 )
  {
    p02_ab ( a, b );
  }
  else if ( p == 3 )
  {
    p03_ab ( a, b );
  }
  else if ( p == 4 )
  {
    p04_ab ( a, b );
  }
  else if ( p == 5 )
  {
    p05_ab ( a, b );
  }
  else if ( p == 6 )
  {
    p06_ab ( a, b );
  }
  else if ( p == 7 )
  {
    p07_ab ( a, b );
  }
  else if ( p == 8 )
  {
    p08_ab ( a, b );
  }
  else if ( p == 9 )
  {
    p09_ab ( a, b );
  }
  else if ( p == 10 )
  {
    p10_ab ( a, b );
  }
  else if ( p == 11 )
  {
    p11_ab ( a, b );
  }
  else if ( p == 12 )
  {
    p12_ab ( a, b );
  }
  else if ( p == 13 )
  {
    p13_ab ( a, b );
  }
  else if ( p == 14 )
  {
    p14_ab ( a, b );
  }
  else
  {
    cout << "\n";
    cout << " 'P00_AB - Fatal error!\n";
    cout << " '  Illegal value of P = " << p << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void p00_e ( int p, double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P00_E returns the exact value of the integral for a given problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the problem index.
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  if ( p == 1 )
  {
    p01_e ( exact );
  }
  else if ( p == 2 )
  {
    p02_e ( exact );
  }
  else if ( p == 3 )
  {
    p03_e ( exact );
  }
  else if ( p == 4 )
  {
    p04_e ( exact );
  }
  else if ( p == 5 )
  {
    p05_e ( exact );
  }
  else if ( p == 6 )
  {
    p06_e ( exact );
  }
  else if ( p == 7 )
  {
    p07_e ( exact );
  }
  else if ( p == 8 )
  {
    p08_e ( exact );
  }
  else if ( p == 9 )
  {
    p09_e ( exact );
  }
  else if ( p == 10 )
  {
    p10_e ( exact );
  }
  else if ( p == 11 )
  {
    p11_e ( exact );
  }
  else if ( p == 12 )
  {
    p12_e ( exact );
  }
  else if ( p == 13 )
  {
    p13_e ( exact );
  }
  else if ( p == 14 )
  {
    p14_e ( exact );
  }
  else
  {
    cout << "\n";
    cout << " 'P00_E - Fatal error!\n";
    cout << "  Illegal value of P = " << p << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void p00_f ( int p, int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_F evaluates the integrand for a given problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the problem index.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  if ( p == 1 )
  {
    p01_f ( n, x, fx );
  }
  else if ( p == 2 )
  {
    p02_f ( n, x, fx );
  }
  else if ( p == 3 )
  {
    p03_f ( n, x, fx );
  }
  else if ( p == 4 )
  {
    p04_f ( n, x, fx );
  }
  else if ( p == 5 )
  {
    p05_f ( n, x, fx );
  }
  else if ( p == 6 )
  {
    p06_f ( n, x, fx );
  }
  else if ( p == 7 )
  {
    p07_f ( n, x, fx );
  }
  else if ( p == 8 )
  {
    p08_f ( n, x, fx );
  }
  else if ( p == 9 )
  {
    p09_f ( n, x, fx );
  }
  else if ( p == 10 )
  {
    p10_f ( n, x, fx );
  }
  else if ( p == 11 )
  {
    p11_f ( n, x, fx );
  }
  else if ( p == 12 )
  {
    p12_f ( n, x, fx );
  }
  else if ( p == 13 )
  {
    p13_f ( n, x, fx );
  }
  else if ( p == 14 )
  {
    p14_f ( n, x, fx );
  }
  else
  {
    cout << "\n";
    cout << "P00_F - Fatal error!\n";
    cout << "  Illegal value of P = " << p << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void p01_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//     P01_AB returns the integration limits for problem 01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p01_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P01_E returns the exact value of the integral for problem 01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  *exact = 1.0;

  return;
}
//****************************************************************************80

void p01_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_F evaluates the integrand for problem 01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0;
  }

  return;
}
//****************************************************************************80

void p02_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P02_AB returns the integration limits for problem 02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p02_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P02_E returns the exact value of the integral for problem 02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  *exact = 1.0;

  return;
}
//****************************************************************************80

void p02_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_F evaluates the integrand for problem 02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 2.0 * x[i];
  }

  return;
}
//****************************************************************************80

void p03_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P03_AB returns the integration limits for problem 03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p03_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P03_E returns the exact value of the integral for problem 03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  *exact = 1.0;

  return;
}
//****************************************************************************80

void p03_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_F evaluates the integrand for problem 03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 3.0 * x[i] * x[i];
  }

  return;
}
//****************************************************************************80

void p04_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P04_AB returns the integration limits for problem 04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = -1.0;
  *b = +1.0;

  return;
}
//****************************************************************************80

void p04_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P04_E returns the exact value of the integral for problem 04.
//
//  Discussion:
//
//    The 20 digit estimate for the exact value comes from Mathematica.
//
//    N [ Integrate [ Exp[-x*x]*Log[1 + x], { x, -1, 1 } ], 20 ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  *exact = -0.31671419631358172053;

  return;
}
//****************************************************************************80

void p04_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_F evaluates the integrand for problem 04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = exp ( - x[i] * x[i] ) * log ( 1.0 - x[i] );
  }

  return;
}
//****************************************************************************80

void p05_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P05_AB returns the integration limits for problem 05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p05_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P05_E returns the exact value of the integral for problem 05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  *exact = 0.25;

  return;
}
//****************************************************************************80

void p05_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_F evaluates the integrand for problem 05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = x[i] * log ( 1.0 + x[i] );
  }

  return;
}
//****************************************************************************80

void p06_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P06_AB  returns the integration limits for problem 06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p06_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P06_E returns the exact value of the integral for problem 06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  double pi = 3.141592653589793;

  *exact = ( pi - 2.0 + 2.0 * log ( 2.0 ) ) / 12.0;

  return;
}
//****************************************************************************80

void p06_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_F evaluates the integrand for problem 06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = x[i] * x[i] * atan ( x[i] );
  }

  return;
}
//****************************************************************************80

void p07_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P07_AB  returns the integration limits for problem 07.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  double pi = 3.141592653589793;

  *a = 0.0;
  *b = pi / 2.0;

  return;
}
//****************************************************************************80

void p07_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P07_E returns the exact value of the integral for problem 07.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  double pi = 3.141592653589793;

  *exact = ( exp ( pi / 2.0 ) - 1.0 ) / 2.0;

  return;
}
//****************************************************************************80

void p07_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_F evaluates the integrand for problem 07.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = exp ( x[i] ) * cos ( x[i] );
  }

  return;
}
//****************************************************************************80

void p08_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P08_AB  returns the integration limits for problem 08.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p08_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P08_E returns the exact value of the integral for problem 08.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  double pi = 3.141592653589793;

  *exact = ( 5.0 * pi * pi ) / 96.0;

  return;
}
//****************************************************************************80

void p08_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_F evaluates the integrand for problem 08.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = atan ( sqrt ( 2.0 + x[i] * x[i] ) ) 
    / ( 1.0 + x[i] * x[i] ) 
    / sqrt ( 2.0 + x[i] * x[i] );
  }

  return;
}
//****************************************************************************80

void p09_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P09_AB  returns the integration limits for problem 09.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p09_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P09_E returns the exact value of the integral for problem 09.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  *exact = - 4.0 / 9.0;

  return;
}
//****************************************************************************80

void p09_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P09_F evaluates the integrand for problem 09.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] =  sqrt ( x[i] ) * log ( x[i] );
  }

  return;
}
//****************************************************************************80

void p10_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P10_AB  returns the integration limits for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p10_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P10_E returns the exact value of the integral for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  double pi = 3.141592653589793;

  *exact = pi / 4.0;

  return;
}
//****************************************************************************80

void p10_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P10_F evaluates the integrand for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = sqrt ( 1.0 - x[i] * x[i] );
  }

  return;
}
//****************************************************************************80

void p11_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P11_AB  returns the integration limits for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p11_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P11_E returns the exact value of the integral for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  double pi = 3.141592653589793;

  *exact = 2.0 * sqrt ( pi ) * gamma ( 0.75 ) / gamma ( 0.25 );

  return;
}
//****************************************************************************80

void p11_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P11_F evaluates the integrand for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = sqrt ( x[i] ) / sqrt ( 1.0 - x[i] * x[i] );
  }

  return;
}
//****************************************************************************80

void p12_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P12_AB  returns the integration limits for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  *a = 0.0;
  *b = 1.0;

  return;
}
//****************************************************************************80

void p12_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P12_E returns the exact value of the integral for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  *exact = 2.0;

  return;
}
//****************************************************************************80

void p12_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P12_F evaluates the integrand for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = pow ( log ( x[i] ), 2 );
  }

  return;
}
//****************************************************************************80

void p13_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P13_AB  returns the integration limits for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  double pi = 3.141592653589793;

  *a = 0.0;
  *b = pi / 2.0;

  return;
}
//****************************************************************************80

void p13_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P13_E returns the exact value of the integral for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  double pi = 3.141592653589793;

  *exact = - pi * log ( 2.0 ) / 2.0;

  return;
}
//****************************************************************************80

void p13_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P13_F evaluates the integrand for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = log ( cos ( x[i] ) );
  }

  return;
}
//****************************************************************************80

void p14_ab ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P14_AB  returns the integration limits for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, the integration limits.
//
{
  double pi = 3.141592653589793;

  *a = 0.0;
  *b = pi / 2.0;

  return;
}
//****************************************************************************80

void p14_e ( double *exact )

//****************************************************************************80
//
//  Purpose:
//
//    P14_E returns the exact value of the integral for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *EXACT, the exact value of the integral.
//
{
  double pi = 3.141592653589793;

  *exact = pi * sqrt ( 2.0 ) / 2.0;

  return;
}
//****************************************************************************80

void p14_f ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P14_F evaluates the integrand for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Karthik Jeyabalan, Xiaoye Li,
//    A Comparison of Three High-Precision Quadrature Schemes,
//    Experimental Mathematics,
//    Volume 14, Number 3, pages 317-329.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fx[i] = sqrt ( tan ( x[i] ) );
  }

  return;
}

