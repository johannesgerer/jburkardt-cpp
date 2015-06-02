# include "sandia_rules.hpp"
# include "sandia_cubature.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

int main ( );
void cn_geg_tests ( );
void cn_geg_test ( int n, double alpha, int expon[] );
void cn_jac_tests ( );
void cn_jac_test ( int n, double alpha, double beta, int expon[] );
void cn_leg_tests ( );
void cn_leg_test ( int n, int expon[] );
void en_her_tests ( );
void en_her_test ( int n, int expon[] );
void epn_glg_tests ( );
void epn_glg_test ( int n, int expon[], double alpha );
void epn_lag_tests ( );
void epn_lag_test ( int n, int expon[] );
void gw_tests ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SANDIA_CUBATURE_PRB.
//
//  Discussion:
//
//    SANDIA_CUBATURE_PRB tests the SANDIA_CUBATURE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SANDIA_CUBATURE_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SANDIA_CUBATURE library.\n";

  cn_geg_tests ( );
  cn_jac_tests ( );
  cn_leg_tests ( );
  en_her_tests ( );
  epn_glg_tests ( );
  epn_lag_tests ( );
  gw_tests ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SANDIA_CUBATURE_PRB\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  webbur::timestamp ( );

  return 0;
}
//****************************************************************************80

void cn_geg_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_TESTS tests the rules for CN with Gegenbauer weight on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  double alpha;
  double alpha_test[TEST_NUM] = { -0.5, 0.0, 0.5, 1.0, 1.5 };
  int *expon;
  int i;
  int n;
  int test;

  std::cout << "\n";
  std::cout << "CN_GEG_TESTS\n";
  std::cout << "  Demonstrate the use of quadrature rules for the region\n";
  std::cout << "  CN_GEG, that is, the hypercube [-1,+1]^N, with the\n";
  std::cout << "  weight W(ALPHA;X) = product ( 1 <= I <= N )\n";
  std::cout << "    (1-X(I)^2)^ALPHA\n";
  std::cout << "\n";
  std::cout << "  We use the formulas to integrate various monomials of\n";
  std::cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  std::cout << "  and compare to the exact integral.\n";
  std::cout << "\n";
  std::cout << "  The precision of each formula is known, and we only use\n";
  std::cout << "  a formula if its precision indicates it should be able to\n";
  std::cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      webbur::i4vec_zero ( n, expon );
      cn_geg_test ( n, alpha, expon );
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      webbur::i4vec_zero ( n, expon );
      expon[n-1] = 1;
      cn_geg_test ( n, alpha, expon );
    }

    if ( 2 <= n )
    {
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];

        webbur::i4vec_zero ( n, expon );
        expon[0] = 1;
        expon[1] = 1;
        cn_geg_test ( n, alpha, expon );
      }
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      webbur::i4vec_zero ( n, expon );
      expon[0] = 2;
      cn_geg_test ( n, alpha, expon );

    }
    delete [] expon;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void cn_geg_test ( int n, double alpha, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_TEST tests the rules for CN with Gegenbauer weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double pi = 3.141592653589793;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "  N = " << n << "\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    std::cout << std::setw(4) << expon[i];
  }
  std::cout << "\n";
  d = webbur::i4vec_sum ( n, expon );
  std::cout << "  Degree = " << d << "\n";
  std::cout << "\n";

  exact = webbur::cn_geg_monomial_integral ( n, alpha, expon );

  p = 1;

  if ( d <= p )
  {
    o = webbur::cn_geg_01_1_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_geg_01_1 ( n, alpha, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_GEG_01_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = webbur::cn_geg_02_xiu_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_geg_02_xiu ( n, alpha, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_GEG_02_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = webbur::gw_02_xiu_size ( n );
    gamma0 = 1.0;
    delta0 = 0.0;
    c1 = 1.0 / ( 2.0 * alpha + 3.0 );
    volume_1d = std::sqrt ( pi ) * webbur::r8_gamma ( alpha + 1.0 ) 
      / webbur::r8_gamma ( alpha + 1.5 );
    x = new double[n*o];
    w = new double[o];
    webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:     "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 3;

  if ( d <= p )
  {
    o = webbur::cn_geg_03_xiu_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_geg_03_xiu ( n, alpha, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_GEG_03_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  std::cout << "  EXACT                  "
       << "  " << std::setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void cn_jac_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_TESTS tests the rules for CN with Jacobi weight on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double alpha;
  double alpha_test[TEST_NUM] = { 0.0, 1.0, 0.0, 0.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.0, 0.0, 2.0, 1.5 };
  int *expon;
  int i;
  int n;
  int test;

  std::cout << "\n";
  std::cout << "CN_JAC_TESTS\n";
  std::cout << "  Demonstrate the use of quadrature rules for the region\n";
  std::cout << "  CN_JAC, that is, the hypercube [-1,+1]^N, with the\n";
  std::cout << "  weight W(ALPHA,BETA;X) = product ( 1 <= I <= N )\n";
  std::cout << "    (1-X(I))^ALPHA (1+X(I))^BETA\n";
  std::cout << "\n";
  std::cout << "  We use the formulas to integrate various monomials of\n";
  std::cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  std::cout << "  and compare to the exact integral.\n";
  std::cout << "\n";
  std::cout << "  The precision of each formula is known, and we only use\n";
  std::cout << "  a formula if its precision indicates it should be able to\n";
  std::cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      webbur::i4vec_zero ( n, expon );
      cn_jac_test ( n, alpha, beta, expon );
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      webbur::i4vec_zero ( n, expon );
      expon[n-1] = 1;
      cn_jac_test ( n, alpha, beta, expon );
    }

    if ( 2 <= n )
    {
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];
        beta  = beta_test[test];

        webbur::i4vec_zero ( n, expon );
        expon[0] = 1;
        expon[1] = 1;
        cn_jac_test ( n, alpha, beta, expon );
      }
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      webbur::i4vec_zero ( n, expon );
      expon[0] = 2;
      cn_jac_test ( n, alpha, beta, expon );

    }
    delete [] expon;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void cn_jac_test ( int n, double alpha, double beta, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_TEST tests the rules for CN with Jacobi weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  int j;

  std::cout << "\n";
  std::cout << "  N = " << n << "\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "  BETA =  " << beta << "\n";
  std::cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    std::cout << std::setw(4) << expon[i];
  }
  std::cout << "\n";
  d = webbur::i4vec_sum ( n, expon );
  std::cout << "  Degree = " << d << "\n";
  std::cout << "\n";

  exact = webbur::cn_jac_monomial_integral ( n, alpha, beta, expon );

  p = 1;

  if ( d <= p )
  {
    o = webbur::cn_jac_01_1_size ( n, alpha, beta );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_jac_01_1 ( n, alpha, beta, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_JAC_01_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = webbur::cn_jac_02_xiu_size ( n, alpha, beta );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_jac_02_xiu ( n, alpha, beta, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_JAC_02_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = webbur::gw_02_xiu_size ( n );
    gamma0 = ( alpha + beta + 2.0 ) / 2.0;
    delta0 = ( alpha - beta ) / 2.0;
    c1 = 2.0 * ( alpha + 1.0 ) * ( beta + 1.0 ) / ( alpha + beta + 3.0 ) 
      / ( alpha + beta + 2.0 );
    volume_1d = std::pow ( 2.0, alpha + beta + 1.0 ) * webbur::r8_gamma ( alpha + 1.0 ) 
      * webbur::r8_gamma ( beta + 1.0 ) / ( alpha + beta + 1.0 ) / webbur::r8_gamma ( alpha + beta + 1.0 );
    x = new double[n*o];
    w = new double[o];
    webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:     "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  std::cout << "  EXACT                  "
       << "  " << std::setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void cn_leg_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_TESTS tests the rules for CN with Legendre weight on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *expon;
  int i;
  int n;
  int test;

  std::cout << "\n";
  std::cout << "CN_LEG_TESTS\n";
  std::cout << "  Demonstrate the use of quadrature rules for the region\n";
  std::cout << "  CN_LEG, that is, the hypercube [-1,+1]^N, with the\n";
  std::cout << "  Legendre weight W(X) = 1.\n";
  std::cout << "\n";
  std::cout << "  We use the formulas to integrate various monomials of\n";
  std::cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  std::cout << "  and compare to the exact integral.\n";
  std::cout << "\n";
  std::cout << "  The precision of each formula is known, and we only use\n";
  std::cout << "  a formula if its precision indicates it should be able to\n";
  std::cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    webbur::i4vec_zero ( n, expon );
    cn_leg_test ( n, expon );

    webbur::i4vec_zero ( n, expon );
    expon[n-1] = 1;
    cn_leg_test ( n, expon );

    if ( 2 <= n )
    {
      webbur::i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      cn_leg_test ( n, expon );
    }

    webbur::i4vec_zero ( n, expon );
    expon[0] = 2;
    cn_leg_test ( n, expon );

    webbur::i4vec_zero ( n, expon );
    expon[0] = 3;
    cn_leg_test ( n, expon );

    webbur::i4vec_zero ( n, expon );
    expon[n-1] = 4;
    cn_leg_test ( n, expon );

    if ( 2 <= n )
    {
      webbur::i4vec_zero ( n, expon );
      expon[0] = 3;
      expon[1] = 2;
      cn_leg_test ( n, expon );
    }
    delete [] expon;
  }

  return;
}
//****************************************************************************80

void cn_leg_test ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_TEST tests the rules for CN with Legendre weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "  N = " << n << "\n";
  std::cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    std::cout << std::setw(4) << expon[i];
  }
  std::cout << "\n";
  d = webbur::i4vec_sum ( n, expon );
  std::cout << "  Degree = " << d << "\n";
  std::cout << "\n";

  exact = webbur::cn_leg_monomial_integral ( n, expon );

  p = 1;

  if ( d <= p )
  {
    o = webbur::cn_leg_01_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_leg_01_1 ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_LEG_01_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = webbur::cn_leg_02_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_leg_02_xiu ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_LEG_02_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = webbur::gw_02_xiu_size ( n );
    gamma0 = 1.0;
    delta0 = 0.0;
    c1 = 1.0 / 3.0;
    volume_1d = 2.0;
    x = new double[n*o];
    w = new double[o];
    webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:     "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 3;

  if ( d <= p )
  {
    o = webbur::cn_leg_03_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_leg_03_1 ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_LEG_03_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = webbur::cn_leg_03_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::cn_leg_03_xiu ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  CN_LEG_03_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 5;

  if ( d <= p )
  {
    if ( 4 <= n && n <= 6 )
    {
      o = webbur::cn_leg_05_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      option = 1;
      webbur::cn_leg_05_1 ( n, option, o, x, w );
      v = webbur::monomial_value ( n, o, x, expon );
      quad = webbur::r8vec_dot_product ( o, w, v );
      err = webbur::r8_abs ( quad - exact );
      std::cout << "  CN_LEG_05_1(1):"
           << "  " << std::setw(6) << o
           << "  " << std::setw(14) << quad
           << "  " << std::setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
    if ( 4 <= n && n <= 5 )
    {
      o = webbur::cn_leg_05_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      option = 2;
      webbur::cn_leg_05_1 ( n, option, o, x, w );
      v = webbur::monomial_value ( n, o, x, expon );
      quad = webbur::r8vec_dot_product ( o, w, v );
      err = webbur::r8_abs ( quad - exact );
      std::cout << "  CN_LEG_05_1(2):"
           << "  " << std::setw(6) << o
           << "  " << std::setw(14) << quad
           << "  " << std::setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
    if ( 2 <= n )
    {
      o = webbur::cn_leg_05_2_size ( n );
      x = new double[n*o];
      w = new double[o];
      webbur::cn_leg_05_2 ( n, o, x, w );
      v = webbur::monomial_value ( n, o, x, expon );
      quad = webbur::r8vec_dot_product ( o, w, v );
      err = webbur::r8_abs ( quad - exact );
      std::cout << "  CN_LEG_05_2:   "
           << "  " << std::setw(6) << o
           << "  " << std::setw(14) << quad
           << "  " << std::setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
  }

  std::cout << "  EXACT                  "
       << "  " << std::setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void en_her_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_TESTS tests the Stroud EN_HER rules on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int *expon;
  int i;
  int n;

  std::cout << "\n";
  std::cout << "EN_HER_TESTS\n";
  std::cout << "  Demonstrate the use of Stroud rules for the region\n";
  std::cout << "  EN_HER, that is, all of N-dimensional space, with the\n";
  std::cout << "  weight function W(X) = exp ( - X1^2 - X2^2 ... -XN^2 )\n";
  std::cout << "\n";
  std::cout << "  We use the formulas to integrate various monomials of\n";
  std::cout << "  the form X1^ALPHA1 * X2^ALPHA2 * ... XN^ALPHAN\n";
  std::cout << "  and compare to the exact integral.\n";
  std::cout << "\n";
  std::cout << "  The precision of each formula is known, and we only use\n";
  std::cout << "  a formula if its precision indicates it should be able to\n";
  std::cout << "  produce an exact result.\n";

  for ( n = 1; n <= 7; n++ )
  {
    expon = new int[n];

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    en_her_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    expon[0] = 2;
    en_her_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    expon[1] = 4;
    en_her_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = i + 1;
    }
    d = webbur::i4vec_sum ( n, expon );
    if ( d <= 5 )
    {
      en_her_test ( n, expon );
    }

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 2;
    }
    d = webbur::i4vec_sum ( n, expon );
    if ( d <= 5 )
    {
      en_her_test ( n, expon );
    }

    delete [] expon;
  }

  return;
}
//****************************************************************************80

void en_her_test ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_TEST tests the Stroud EN_HER rules on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double pi = 3.141592653589793;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "  N = " << n << "\n";
  std::cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    std::cout << "  " << std::setw(4) << expon[i];
  }
  std::cout << "\n";
  d = webbur::i4vec_sum ( n, expon );
  std::cout << "  Degree = " << d << "\n";
  std::cout << "\n";

  exact = webbur::en_her_monomial_integral ( n, expon );

  p = 1;

  if ( d <= p )
  {
    o = webbur::en_her_01_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::en_her_01_1 ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EN_HER_01_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = webbur::en_her_02_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::en_her_02_xiu ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EN_HER_02_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = webbur::gw_02_xiu_size ( n );
    gamma0 = 2.0;
    delta0 = 0.0;
    c1 = 1.0;
    volume_1d = std::sqrt ( pi );
    x = new double[n*o];
    w = new double[o];
    webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:     "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 3;

  if ( d <= p )
  {
    o = webbur::en_her_03_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::en_her_03_1 ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EN_HER_03_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = webbur::en_her_03_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::en_her_03_xiu ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EN_HER_03_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 5;

  if ( d <= p )
  {
    if ( 2 <= n && n <= 7) 
    {
      option = 1;
      o = webbur::en_her_05_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      webbur::en_her_05_1 ( n, option, o, x, w );
      v = webbur::monomial_value ( n, o, x, expon );
      quad = webbur::r8vec_dot_product ( o, w, v );
      err = webbur::r8_abs ( quad - exact );
      std::cout << "  EN_HER_05_1(1):"
           << "  " << std::setw(6) << o
           << "  " << std::setw(14) << quad
           << "  " << std::setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }

    if ( n == 3 || n == 5 || n == 6 )
    {
      option = 2;
      o = webbur::en_her_05_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      webbur::en_her_05_1 ( n, option, o, x, w );
      v = webbur::monomial_value ( n, o, x, expon );
      quad = webbur::r8vec_dot_product ( o, w, v );
      err = webbur::r8_abs ( quad - exact );
      std::cout << "  EN_HER_05_1(2):"
           << "  " << std::setw(6) << o
           << "  " << std::setw(14) << quad
           << "  " << std::setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }

    o = webbur::en_her_05_2_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::en_her_05_2 ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EN_HER_05_2:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }
  std::cout << "  EXACT                  "
       << "  " << std::setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void epn_glg_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_TESTS tests the rules for EPN with GLG on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  double alpha;
  double alpha_test[TEST_NUM] = { -0.5, 0.0, 0.5, 1.0, 2.0 };
  int *expon;
  int i;
  int n;
  int test;

  std::cout << "\n";
  std::cout << "EPN_GLG_TESTS\n";
  std::cout << "  Demonstrate the use of quadrature rules for the region\n";
  std::cout << "  EPN_GLG, that is, the positive half space [0,+oo)^N, with the\n";
  std::cout << "  weight W(ALPHA;X) = product ( 1 <= I <= N ) X(I)^ALPHA exp ( -X(I) )\n";
  std::cout << "\n";
  std::cout << "  We use the formulas to integrate various monomials of\n";
  std::cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  std::cout << "  and compare to the exact integral.\n";
  std::cout << "\n";
  std::cout << "  The precision of each formula is known, and we only use\n";
  std::cout << "  a formula if its precision indicates it should be able to\n";
  std::cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    webbur::i4vec_zero ( n, expon );
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }

    webbur::i4vec_zero ( n, expon );
    expon[n-1] = 1;
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }

    if ( 2 <= n )
    {
      webbur::i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];
        epn_glg_test ( n, expon, alpha );
      }
    }

    webbur::i4vec_zero ( n, expon );
    expon[0] = 2;
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }
    delete [] expon;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void epn_glg_test ( int n, int expon[], double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_TEST tests the rules for EPN with GLG weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "  N = " << n << "\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    std::cout << std::setw(4) << expon[i];
  }
  std::cout << "\n";
  d = webbur::i4vec_sum ( n, expon );
  std::cout << "  Degree = " << d << "\n";
  std::cout << "\n";

  exact = webbur::epn_glg_monomial_integral ( n, expon, alpha );

  p = 1;

  if ( d <= p )
  {
    o = webbur::epn_glg_01_1_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    webbur::epn_glg_01_1 ( n, alpha, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EPN_GLG_01_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = webbur::epn_glg_02_xiu_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    webbur::epn_glg_02_xiu ( n, alpha, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EPN_GLG_02_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = webbur::gw_02_xiu_size ( n );
    gamma0 = - 1.0;
    delta0 = alpha + 1.0;
    c1 = - alpha - 1.0;
    volume_1d = webbur::r8_gamma ( 1.0 + alpha );
    x = new double[n*o];
    w = new double[o];
    webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:      "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  std::cout << "  EXACT                   "
       << "  " << std::setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void epn_lag_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_TESTS tests the rules for EPN with Laguerre weight on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *expon;
  int i;
  int n;
  int test;

  std::cout << "\n";
  std::cout << "EPN_LAG_TESTS\n";
  std::cout << "  Demonstrate the use of quadrature rules for the region\n";
  std::cout << "  EPN_LAG, that is, the positive half space [0,+oo)^N, with the\n";
  std::cout << "  weight W(X) = product ( 1 <= I <= N ) exp ( -X(I) )\n";
  std::cout << "\n";
  std::cout << "  We use the formulas to integrate various monomials of\n";
  std::cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  std::cout << "  and compare to the exact integral.\n";
  std::cout << "\n";
  std::cout << "  The precision of each formula is known, and we only use\n";
  std::cout << "  a formula if its precision indicates it should be able to\n";
  std::cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    webbur::i4vec_zero ( n, expon );
    epn_lag_test ( n, expon );

    webbur::i4vec_zero ( n, expon );
    expon[n-1] = 1;
    epn_lag_test ( n, expon );

    if ( 2 <= n )
    {
      webbur::i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      epn_lag_test ( n, expon );
    }

    webbur::i4vec_zero ( n, expon );
    expon[0] = 2;
    epn_lag_test ( n, expon );

    delete [] expon;
  }

  return;
}
//****************************************************************************80

void epn_lag_test ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_TEST tests the rules for EPN with Laguerre weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "  N = " << n << "\n";
  std::cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    std::cout << std::setw(4) << expon[i];
  }
  std::cout << "\n";
  d = webbur::i4vec_sum ( n, expon );
  std::cout << "  Degree = " << d << "\n";
  std::cout << "\n";

  exact = webbur::epn_lag_monomial_integral ( n, expon );

  p = 1;

  if ( d <= p )
  {
    o = webbur::epn_lag_01_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::epn_lag_01_1 ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EPN_LAG_01_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = webbur::epn_lag_02_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    webbur::epn_lag_02_xiu ( n, o, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  EPN_LAG_02_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = webbur::gw_02_xiu_size ( n );
    gamma0 = - 1.0;
    delta0 = 1.0;
    c1 = - 1.0;
    volume_1d = 1.0;
    x = new double[n*o];
    w = new double[o];
    webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = webbur::monomial_value ( n, o, x, expon );
    quad = webbur::r8vec_dot_product ( o, w, v );
    err = webbur::r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:      "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  std::cout << "  EXACT                   "
       << "  " << std::setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void gw_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    GW_TESTS tests the rules for GW on monomials.
//
//  Discussion:
//
//    Right now, this test simply calls the GW rule for each type of
//    weight function for which the orthogonal polynomials are known.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  double c1;
  double delta0;
  double gamma0;
  int *expon;
  int i;
  int j;
  int n;
  int o;
  double pi = 3.141592653589793;
  int test;
  double volume_1d;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "GW_TESTS\n";
  std::cout << "  Demonstrate the use of quadrature rules for a Golub Welsch rule\n";
  std::cout << "  defined over some interval and some weight function for which\n";
  std::cout << "  the three term recursion of the orthogonal polynomials is known.\n";
//
//  For a given dimension N, the rule is always the same size.
//
  n = 2;
  o = webbur::gw_02_xiu_size ( n );
  x = new double[n*o];
  w = new double[o];
//
//  Chebyshev Type 1.
//
  gamma0 = 1.0;
  delta0 = 0.0;
  c1 = 0.5;
  volume_1d = pi;
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Chebyshev1:\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Chebyshev Type 2.
//
  gamma0 = 2.0;
  delta0 = 0.0;
  c1 = 0.5;
  volume_1d = pi / 2.0;
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Chebyshev2:\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Gegenbauer
//
  alpha = 1.0;
  gamma0 = 1.0;
  delta0 = 0.0;
  c1 = 1.0 / ( 2.0 * alpha + 3.0 );
  volume_1d = sqrt ( pi ) * webbur::r8_gamma ( alpha + 1.0 ) / webbur::r8_gamma ( alpha + 1.5 );
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Gegenbauer:\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Generalized Hermite.
//
  alpha = 1.0;
  gamma0 = 2.0;
  delta0 = 0.0;
  c1 = 2.0 + 2.0 * alpha;
  volume_1d = webbur::r8_gamma ( ( alpha + 1.0 ) / 2.0 );
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Generalized Hermite:\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Generalized Laguerre.
//
  alpha = 1.0;
  gamma0 = - 1.0;
  delta0 = alpha + 1.0;
  c1 = - alpha - 1.0;
  volume_1d = webbur::r8_gamma ( alpha + 1.0 );
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Generalized Laguerre:\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Hermite (physicist)
//
  gamma0 = 2.0;
  delta0 = 0.0;
  c1 = 1.0;
  volume_1d = std::sqrt ( pi );
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Hermite (physicist):\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Hermite (probabilist)
//
  gamma0 = 1.0;
  delta0 = 0.0;
  c1 = 1.0;
  volume_1d = std::sqrt ( 2.0 * pi );
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Hermite ( probabilist):\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Jacobi.
//
  alpha = 0.5;
  beta = 1.5;
  gamma0 = ( alpha + beta + 2.0 ) / 2.0;
  delta0 = ( alpha - beta ) / 2.0;
  c1 = 2.0 * ( alpha + 1.0 ) * ( beta + 1.0 ) / ( alpha + beta + 3.0 )
    / ( alpha + beta + 2.0 );
  volume_1d = std::pow ( 2.0, alpha + beta + 1.0 ) * webbur::r8_gamma ( alpha + 1 )
    * webbur::r8_gamma ( beta + 1.0 ) / ( alpha + beta + 1.0 ) 
    / webbur::r8_gamma ( alpha + beta + 1.0 );
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Jacobi:\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "  BETA  = " << beta << "\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Laguerre.
//
  gamma0 = - 1.0;
  delta0 = 1.0;
  c1 = - 1.0;
  volume_1d = 1.0;
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Laguerre:\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }
//
//  Legendre.
//
  gamma0 = 1.0;
  delta0 = 0.0;
  c1 = 1.0 / 3.0;
  volume_1d = 2.0;
  webbur::gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
  std::cout << "\n";
  std::cout << "  Legendre:\n";
  std::cout << "\n";
  for ( j = 0; j < o; j++ )
  {
    std::cout << std::setw(14) << w[j];
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(14) << x[i+j*n];
    }
    std::cout << "\n";
  }

  delete [] w;
  delete [] x;

  return;
}
