# include "sandia_rules.hpp"

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>

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
void test18 ( );
void test19 ( );
void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( int r );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );
void test28 ( );
void test285 ( );
void test29 ( );
void test30 ( );
void test31 ( );
void test32 ( );
void test33 ( );
void test34 ( );
void test35 ( );
void test36 ( );
void test37 ( );
void test38 ( );
void test39 ( );

void test01_np ( );
void test02_np ( );
void test03_np ( );
void test04_np ( );
void test05_np ( );
void test06_np ( );
void test07_np ( );
void test08_np ( );
void test09_np ( );
void test10_np ( );
void test11_np ( );
void test12_np ( );
void test13_np ( );
void test14_np ( );
void test15_np ( );
void test16_np ( );
void test17_np ( );
void test18_np ( );
void test19_np ( );
void test20_np ( );
void test21_np ( );
void test22_np ( );

double cubic_antiderivative ( double x );
double cubic_integrate ( double a, double b );
void cubic_value ( double x, double *f, double *d, double *s, double *t );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN tests the SANDIA_RULES library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  int r;

  webbur::timestamp ( );

  std::cout << "\n";
  std::cout << "SANDIA_RULES_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SANDIA_RULES library.\n";

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
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
//
//  Repeat tests, but now call "NP" versions of routines.
//
  test01_np ( );
  test02_np ( );
  test03_np ( );
  test04_np ( );
  test05_np ( );
  test06_np ( );
  test07_np ( );
  test08_np ( );
  test09_np ( );

  test10_np ( );
  test11_np ( );
  test12_np ( );
  test13_np ( );
  test14_np ( );
  test15_np ( );
  test16_np ( );
  test17_np ( );
  test18_np ( );
  test19_np ( );

  test20_np ( );
  test21_np ( );
  test22_np ( );
//
//  TEST23 takes an input argument of R, a rule index.
//
  r = 1;
  test23 ( r );
  r = 3;
  test23 ( r );
  r = 4;
  test23 ( r );
  r = 11;
  test23 ( r );

  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test285 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test39 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SANDIA_RULES_PRB\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CHEBYSHEV1_COMPUTE against CHEBYSHEV1_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
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
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST01\n";
  std::cout << "  CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) / sqrt ( 1 - x^2 ) dx.\n";
  std::cout << "\n";
  std::cout << "  CHEBYSHEV1_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::chebyshev1_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::chebyshev1_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
                << "  " << std::setw(8)  << n
                << "  " << std::setprecision(6) << std::setw(14) << estimate
                << "  " << std::setprecision(6) << std::setw(14) << exact
                << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
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
//    TEST02 tests CHEBYSHEV1_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST02\n";
  std::cout << "  CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) / sqrt(1-x^2) dx.\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::chebyshev1_compute ( order, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
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
//    TEST03 tests CHEBYSHEV2_COMPUTE against CHEBYSHEV2_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
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
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST03\n";
  std::cout << "  CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) * sqrt ( 1 - x^2 ) dx.\n";
  std::cout << "\n";
  std::cout << "  CHEBYSHEV2_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::chebyshev2_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::chebyshev2_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests CHEBYSHEV2_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST04\n";
  std::cout << "  CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) * sqrt(1-x^2) dx.\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::chebyshev2_compute ( order, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests CLENSHAW_CURTIS_COMPUTE against LEGENDRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
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
  int n;
  int n_hi;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST05\n";
  std::cout << "  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n";
  std::cout << "\n";
  std::cout << "  LEGENDRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N up to\n";
  std::cout << "    N = ORDER+1 if ORDER is odd, or\n";
  std::cout << "    N = ORDER   if ORDER is even\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::clenshaw_curtis_compute ( order, x, w );

    if ( ( order % 2 ) == 0 )
    {
      n_hi = order + 2;
    }
    else
    {
      n_hi = order + 3;
    }

    for ( n = 0; n <= n_hi; n = n + 1 )
    {
      exact = webbur::legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
                << "  " << std::setw(8)  << n
                << "  " << std::setprecision(6) << std::setw(14) << estimate
                << "  " << std::setprecision(6) << std::setw(14) << exact
                << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
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
//    TEST06 tests CLENSHAW_CURTIS_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST06\n";
  std::cout << "  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) dx.\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::clenshaw_curtis_compute ( order, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests FEJER2_COMPUTE against LEGENDRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
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
  int n;
  int n_hi;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST07\n";
  std::cout << "  FEJER2_COMPUTE computes a Fejer Type 2 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n";
  std::cout << "\n";
  std::cout << "  LEGENDRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N up to\n";
  std::cout << "    N = ORDER+1 if ORDER is odd, or\n";
  std::cout << "    N = ORDER   if ORDER is even\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::fejer2_compute ( order, x, w );

    if ( ( order % 2 ) == 0 )
    {
      n_hi = order + 2;
    }
    else
    {
      n_hi = order + 3;
    }

    for ( n = 0; n <= n_hi; n = n + 1 )
    {
      exact = webbur::legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests FEJER2_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST08\n";
  std::cout << "  FEJER2_COMPUTE computes a Fejer Type 2 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) dx.\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::fejer2_compute ( order, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests GEGENBAUER_COMPUTE against GEGENBAUER_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST09\n";
  std::cout << "  GEGENBAUER_COMPUTE computes a generalized Gauss-Gegenbauer rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) (1-x^2)^alpha dx.\n";
  std::cout << "\n";
  std::cout << "  GEGENBAUER_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Alpha           Estimate       Exact            Error\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";

      f = new double[order];
      w = new double[order];
      x = new double[order];

      webbur::gegenbauer_compute ( order, alpha, x, w );

      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = webbur::gegenbauer_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = std::pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = webbur::r8_abs ( exact - estimate );

        std::cout << "  " << std::setw(8)  << order
             << "  " << std::setw(8)  << n
             << "  " << std::setw(14) << alpha
             << "  " << std::setprecision(6) << std::setw(14) << estimate
             << "  " << std::setprecision(6) << std::setw(14) << exact
             << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
      }
      delete [] f;
      delete [] w;
      delete [] x;
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests GEGENBAUER_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST10\n";
  std::cout << "  GEGENBAUER_COMPUTE computes a generalized Gauss-Gegenbauer rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) (1-x^2)^alpha dx.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";
      std::cout << "  Order = " << order << "\n";
      std::cout << "  ALPHA = " << alpha << "\n";

      w = new double[order];
      x = new double[order];

      webbur::gegenbauer_compute ( order, alpha, x, w );

      for ( i = 0; i < order; i =i + 1 )
      {
        std::cout << "  " << std::setw(8) << i
             << "  " << std::setprecision(16) << std::setw(24) << x[i]
             << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
      }
      delete [] w;
      delete [] x;
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests GEN_HERMITE_COMPUTE against GEN_HERMITE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST11\n";
  std::cout << "  GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.\n";
  std::cout << "\n";
  std::cout << "  GEN_HERMITE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Alpha           Estimate       Exact            Error\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";

      f = new double[order];
      w = new double[order];
      x = new double[order];

      webbur::gen_hermite_compute ( order, alpha, x, w );

      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = webbur::gen_hermite_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = std::pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = webbur::r8_abs ( exact - estimate );

        std::cout << "  " << std::setw(8)  << order
             << "  " << std::setw(8)  << n
             << "  " << std::setprecision(6) << std::setw(14) << alpha
             << "  " << std::setprecision(6) << std::setw(14) << estimate
             << "  " << std::setprecision(6) << std::setw(14) << exact
             << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
      }
      delete [] f;
      delete [] w;
      delete [] x;
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests GEN_HERMITE_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST12\n";
  std::cout << "  GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";
      std::cout << "  Order = " << order << "\n";
      std::cout << "  ALPHA = " << alpha << "\n";

      w = new double[order];
      x = new double[order];

      webbur::gen_hermite_compute ( order, alpha, x, w );

      for ( i = 0; i < order; i =i + 1 )
      {
        std::cout << "  " << std::setw(8) << i
             << "  " << std::setprecision(16) << std::setw(24) << x[i]
             << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
      }
      delete [] w;
      delete [] x;
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests GEN_LAGUERRE_COMPUTE against GEN_LAGUERRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST13\n";
  std::cout << "  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.\n";
  std::cout << "\n";
  std::cout << "  GEN_LAGUERRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Alpha           Estimate       Exact            Error\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";

      f = new double[order];
      w = new double[order];
      x = new double[order];

      webbur::gen_laguerre_compute ( order, alpha, x, w );

      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = webbur::gen_laguerre_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = std::pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = webbur::r8_abs ( exact - estimate );

        std::cout << "  " << std::setw(8)  << order
             << "  " << std::setw(8)  << n
             << "  " << std::setprecision(6) << std::setw(14) << alpha
             << "  " << std::setprecision(6) << std::setw(14) << estimate
             << "  " << std::setprecision(6) << std::setw(14) << exact
             << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
      }
      delete [] f;
      delete [] w;
      delete [] x;
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests GEN_LAGUERRE_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST14\n";
  std::cout << "  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";
      std::cout << "  Order = " << order << "\n";
      std::cout << "  ALPHA = " << alpha << "\n";

      w = new double[order];
      x = new double[order];

      webbur::gen_laguerre_compute ( order, alpha, x, w );

      for ( i = 0; i < order; i =i + 1 )
      {
        std::cout << "  " << std::setw(8) << i
             << "  " << std::setprecision(16) << std::setw(24) << x[i]
             << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
      }
      delete [] w;
      delete [] x;
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests HERMITE_COMPUTE against HERMITE_INTEGRAL.
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
  int d;
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int n_max = 10;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST15\n";
  std::cout << "  HERMITE_COMPUTE computes a Gauss-Hermite rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.\n";
  std::cout << "\n";
  std::cout << "  HERMITE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^d.\n";
  std::cout << "\n";
  std::cout << "  A rule of order N should be exact for monomials x^d\n";
  std::cout << "  up to degree D = 2*N-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "         N         D       Estimate       Exact            Error\n";

  for ( n = 1; n <= n_max; n++ )
  {
    std::cout << "\n";

    f = new double[n];
    w = new double[n];
    x = new double[n];

    webbur::hermite_compute ( n, x, w );

    for ( d = 0; d <= 2 * n + 2; d = d + 1 )
    {
      exact = webbur::hermite_integral ( d );

      if ( d == 0 )
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
          f[i] = std::pow ( x[i], d );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < n; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << n
           << "  " << std::setw(8)  << d
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests HERMITE_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST16\n";
  std::cout << "  HERMITE_COMPUTE computes a Gauss-Hermite rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::hermite_compute ( order, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests JACOBI_COMPUTE against JACOBI_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  int test1;
  int test2;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST17\n";
  std::cout << "  JACOBI_COMPUTE computes a Gauss-Jacobi rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n";
  std::cout << "\n";
  std::cout << "  JACOBI_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Alpha           Beta            Estimate       Exact            Error\n";

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];

      for ( order = 1; order <= order_max; order++ )
      {
        std::cout << "\n";

        f = new double[order];
        w = new double[order];
        x = new double[order];

        webbur::jacobi_compute ( order, alpha, beta, x, w );

        for ( n = 0; n <= 2 * order + 2; n = n + 1 )
        {
          exact = webbur::jacobi_integral ( n, alpha, beta );

          if ( n == 0 )
          {
            for ( i = 0; i < order; i++ )
            {
              f[i] = 1.0;
            }
          }
          else
          {
            for ( i = 0; i < order; i++ )
            {
              f[i] = std::pow ( x[i], n );
            }
          }
          estimate = 0.0;
          for ( i = 0; i < order; i++ )
          {
            estimate = estimate + w[i] * f[i];
          }

          error = webbur::r8_abs ( exact - estimate );

          std::cout << "  " << std::setw(8)  << order
               << "  " << std::setw(8)  << n
               << "  " << std::setprecision(6) << std::setw(14) << alpha
               << "  " << std::setprecision(6) << std::setw(14) << beta
               << "  " << std::setprecision(6) << std::setw(14) << estimate
               << "  " << std::setprecision(6) << std::setw(14) << exact
               << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
        }
        delete [] f;
        delete [] w;
        delete [] x;
      }
    }
  }
  return;
# undef TEST_NUM
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests JACOBI_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int n;
  int n_max = 10;
  int test1;
  int test2;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST18\n";
  std::cout << "  JACOBI_COMPUTE computes a Gauss-Jacobi rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n";

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];

      for ( n = 1; n <= n_max; n++ )
      {
        std::cout << "\n";
        std::cout << "  N = " << n << "\n";
        std::cout << "  ALPHA = " << alpha << "\n";
        std::cout << "  BETA = "  << beta << "\n";

        w = new double[n];
        x = new double[n];

        webbur::jacobi_compute ( n, alpha, beta, x, w );

        for ( i = 0; i < n; i++ )
        {
          std::cout << "  " << std::setw(8) << i
                    << "  " << std::setprecision(16) << std::setw(24) << x[i]
                    << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
        }
        delete [] w;
        delete [] x;
      }
    }
  }
  return;
# undef TEST_NUM
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests LAGUERRE_COMPUTE against LAGUERRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
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
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST19\n";
  std::cout << "  LAGUERRE_COMPUTE computes a Gauss-Laguerre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.\n";
  std::cout << "\n";
  std::cout << "  LAGUERRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::laguerre_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::laguerre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests LAGUERRE_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
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

  std::cout << "\n";
  std::cout << "TEST20\n";
  std::cout << "  LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::laguerre_compute ( order, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests LEGENDRE_COMPUTE against LEGENDRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
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
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST21\n";
  std::cout << "  LEGENDRE_COMPUTE computes a Gauss-Legendre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n";
  std::cout << "\n";
  std::cout << "  LEGENDRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::legendre_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests LEGENDRE_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST22\n";
  std::cout << "  LEGENDRE_COMPUTE computes a Gauss-Legendre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) dx.\n";

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::legendre_compute ( order, x, w );

    for ( i = 0; i < order; i++ )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test01_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01_NP tests CHEBYSHEV1_COMPUTE_NP against CHEBYSHEV1_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
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
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST01\n";
  std::cout << "  CHEBYSHEV1_COMPUTE_NP computes a Gauss-Chebyshev type 1 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) / sqrt ( 1 - x^2 ) dx.\n";
  std::cout << "\n";
  std::cout << "  CHEBYSHEV1_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::chebyshev1_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::chebyshev1_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test02_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02_NP tests CHEBYSHEV1_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST02\n";
  std::cout << "  CHEBYSHEV1_COMPUTE_NP computes a Gauss-Chebyshev type 1 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) / sqrt(1-x^2) dx.\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::chebyshev1_compute_np ( order, np, p, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test03_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03_NP tests CHEBYSHEV2_COMPUTE_NP against CHEBYSHEV2_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
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
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST03\n";
  std::cout << "  CHEBYSHEV2_COMPUTE_NP computes a Gauss-Chebyshev type 2 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) * sqrt ( 1 - x^2 ) dx.\n";
  std::cout << "\n";
  std::cout << "  CHEBYSHEV2_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::chebyshev2_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::chebyshev2_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test04_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04_NP tests CHEBYSHEV2_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST04\n";
  std::cout << "  CHEBYSHEV2_COMPUTE_NP computes a Gauss-Chebyshev type 2 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) * sqrt(1-x^2) dx.\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::chebyshev2_compute_np ( order, np, p, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test05_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05_NP tests CLENSHAW_CURTIS_COMPUTE_NP against LEGENDRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
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
  int n;
  int n_hi;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST05\n";
  std::cout << "  CLENSHAW_CURTIS_COMPUTE_NP computes a Clenshaw Curtis rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n";
  std::cout << "\n";
  std::cout << "  LEGENDRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N up to\n";
  std::cout << "    N = ORDER+1 if ORDER is odd, or\n";
  std::cout << "    N = ORDER   if ORDER is even\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::clenshaw_curtis_compute_np ( order, np, p, x, w );

    if ( ( order % 2 ) == 0 )
    {
      n_hi = order + 2;
    }
    else
    {
      n_hi = order + 3;
    }

    for ( n = 0; n <= n_hi; n = n + 1 )
    {
      exact = webbur::legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test06_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06_NP tests CLENSHAW_CURTIS_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST06\n";
  std::cout << "  CLENSHAW_CURTIS_COMPUTE_NP computes a Clenshaw Curtis rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) dx.\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::clenshaw_curtis_compute_np ( order, np, p, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test07_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07_NP tests FEJER2_COMPUTE_NP against LEGENDRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
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
  int n;
  int n_hi;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST07\n";
  std::cout << "  FEJER2_COMPUTE_NP computes a Fejer Type 2 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n";
  std::cout << "\n";
  std::cout << "  LEGENDRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N up to\n";
  std::cout << "    N = ORDER+1 if ORDER is odd, or\n";
  std::cout << "    N = ORDER   if ORDER is even\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::fejer2_compute_np ( order, np, p, x, w );

    if ( ( order % 2 ) == 0 )
    {
      n_hi = order + 2;
    }
    else
    {
      n_hi = order + 3;
    }

    for ( n = 0; n <= n_hi; n = n + 1 )
    {
      exact = webbur::legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test08_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08_NP tests FEJER2_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST08\n";
  std::cout << "  FEJER2_COMPUTE_NP computes a Fejer Type 2 rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) dx.\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::fejer2_compute_np ( order, np, p, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test09_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09_NP tests GEGENBAUER_COMPUTE_NP against GEGENBAUER_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST09\n";
  std::cout << "  GEGENBAUER_COMPUTE_NP computes a generalized Gauss-Gegenbauer rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) (1-x^2)^alpha dx.\n";
  std::cout << "\n";
  std::cout << "  GEGENBAUER_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Alpha           Estimate       Exact            Error\n";

  p = new double[np];

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";

      f = new double[order];
      w = new double[order];
      x = new double[order];

      webbur::gegenbauer_compute_np ( order, np, p, x, w );

      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = webbur::gegenbauer_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = std::pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = webbur::r8_abs ( exact - estimate );

        std::cout << "  " << std::setw(8)  << order
             << "  " << std::setw(8)  << n
             << "  " << std::setw(14) << alpha
             << "  " << std::setprecision(6) << std::setw(14) << estimate
             << "  " << std::setprecision(6) << std::setw(14) << exact
             << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
      }
      delete [] f;
      delete [] w;
      delete [] x;
    }
  }

  delete [] p;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test10_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10_NP tests GEGENBAUER_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST10\n";
  std::cout << "  GEGENBAUER_COMPUTE_NP computes a generalized Gauss-Gegenbauer rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) (1-x^2)^alpha dx.\n";

  p = new double[np];

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";
      std::cout << "  Order = " << order << "\n";
      std::cout << "  ALPHA = " << alpha << "\n";

      w = new double[order];
      x = new double[order];

      webbur::gegenbauer_compute_np ( order, np, p, x, w );

      for ( i = 0; i < order; i =i + 1 )
      {
        std::cout << "  " << std::setw(8) << i
             << "  " << std::setprecision(16) << std::setw(24) << x[i]
             << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
      }
      delete [] w;
      delete [] x;
    }
  }

  delete [] p;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test11_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11_NP tests GEN_HERMITE_COMPUTE_NP against GEN_HERMITE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST11\n";
  std::cout << "  GEN_HERMITE_COMPUTE_NP computes a generalized Gauss-Hermite rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.\n";
  std::cout << "\n";
  std::cout << "  GEN_HERMITE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Alpha           Estimate       Exact            Error\n";

  p = new double[np];

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";

      f = new double[order];
      w = new double[order];
      x = new double[order];

      webbur::gen_hermite_compute_np ( order, np, p, x, w );

      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = webbur::gen_hermite_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = std::pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = webbur::r8_abs ( exact - estimate );

        std::cout << "  " << std::setw(8)  << order
             << "  " << std::setw(8)  << n
             << "  " << std::setprecision(6) << std::setw(14) << alpha
             << "  " << std::setprecision(6) << std::setw(14) << estimate
             << "  " << std::setprecision(6) << std::setw(14) << exact
             << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
      }
      delete [] f;
      delete [] w;
      delete [] x;
    }
  }

  delete [] p;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test12_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12_NP tests GEN_HERMITE_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST12\n";
  std::cout << "  GEN_HERMITE_COMPUTE_NP computes a generalized Gauss-Hermite rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.\n";

  p = new double[np];

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";
      std::cout << "  Order = " << order << "\n";
      std::cout << "  ALPHA = " << alpha << "\n";

      w = new double[order];
      x = new double[order];

      webbur::gen_hermite_compute_np ( order, np, p, x, w );

      for ( i = 0; i < order; i =i + 1 )
      {
        std::cout << "  " << std::setw(8) << i
             << "  " << std::setprecision(16) << std::setw(24) << x[i]
             << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
      }
      delete [] w;
      delete [] x;
    }
  }

  delete [] p;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test13_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13_NP tests GEN_LAGUERRE_COMPUTE_NP against GEN_LAGUERRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST13\n";
  std::cout << "  GEN_LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.\n";
  std::cout << "\n";
  std::cout << "  GEN_LAGUERRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Alpha           Estimate       Exact            Error\n";

  p = new double[np];

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";

      f = new double[order];
      w = new double[order];
      x = new double[order];

      webbur::gen_laguerre_compute_np ( order, np, p, x, w );

      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = webbur::gen_laguerre_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = std::pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = webbur::r8_abs ( exact - estimate );

        std::cout << "  " << std::setw(8)  << order
             << "  " << std::setw(8)  << n
             << "  " << std::setprecision(6) << std::setw(14) << alpha
             << "  " << std::setprecision(6) << std::setw(14) << estimate
             << "  " << std::setprecision(6) << std::setw(14) << exact
             << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
      }
      delete [] f;
      delete [] w;
      delete [] x;
    }
  }

  delete [] p;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test14_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14_NP tests GEN_LAGUERRE_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST14\n";
  std::cout << "  GEN_LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.\n";

  p = new double[np];

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      std::cout << "\n";
      std::cout << "  Order = " << order << "\n";
      std::cout << "  ALPHA = " << alpha << "\n";

      w = new double[order];
      x = new double[order];

      webbur::gen_laguerre_compute_np ( order, np, p, x, w );

      for ( i = 0; i < order; i =i + 1 )
      {
        std::cout << "  " << std::setw(8) << i
             << "  " << std::setprecision(16) << std::setw(24) << x[i]
             << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
      }
      delete [] w;
      delete [] x;
    }
  }

  delete [] p;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test15_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15_NP tests HERMITE_COMPUTE_NP against HERMITE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
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
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST15\n";
  std::cout << "  HERMITE_COMPUTE_NP computes a Gauss-Hermite rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.\n";
  std::cout << "\n";
  std::cout << "  HERMITE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::hermite_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::hermite_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test16_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16_NP tests HERMITE_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST16\n";
  std::cout << "  HERMITE_COMPUTE_NP computes a Gauss-Hermite rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::hermite_compute_np ( order, np, p, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test17_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17_NP tests JACOBI_COMPUTE_NP against JACOBI_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 2;
  int order;
  int order_max = 10;
  double *p;
  int test1;
  int test2;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST17\n";
  std::cout << "  JACOBI_COMPUTE_NP computes a Gauss-Jacobi rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n";
  std::cout << "\n";
  std::cout << "  JACOBI_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Alpha           Beta            Estimate       Exact            Error\n";

  p = new double[np];

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];
    p[0] = alpha;

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];
      p[1] = beta;

      for ( order = 1; order <= order_max; order++ )
      {
        std::cout << "\n";

        f = new double[order];
        w = new double[order];
        x = new double[order];

        webbur::jacobi_compute_np ( order, np, p, x, w );

        for ( n = 0; n <= 2 * order + 2; n = n + 1 )
        {
          exact = webbur::jacobi_integral ( n, alpha, beta );

          if ( n == 0 )
          {
            for ( i = 0; i < order; i++ )
            {
              f[i] = 1.0;
            }
          }
          else
          {
            for ( i = 0; i < order; i++ )
            {
              f[i] = std::pow ( x[i], n );
            }
          }
          estimate = 0.0;
          for ( i = 0; i < order; i++ )
          {
            estimate = estimate + w[i] * f[i];
          }

          error = webbur::r8_abs ( exact - estimate );

          std::cout << "  " << std::setw(8)  << order
               << "  " << std::setw(8)  << n
               << "  " << std::setprecision(6) << std::setw(14) << alpha
               << "  " << std::setprecision(6) << std::setw(14) << beta
               << "  " << std::setprecision(6) << std::setw(14) << estimate
               << "  " << std::setprecision(6) << std::setw(14) << exact
               << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
        }
        delete [] f;
        delete [] w;
        delete [] x;
      }
    }
  }

  delete [] p;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test18_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18_NP tests JACOBI_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int np = 2;
  int order;
  int order_max = 10;
  double *p;
  int test1;
  int test2;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST18\n";
  std::cout << "  JACOBI_COMPUTE_NP computes a Gauss-Jacobi rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n";

  p = new double[np];

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];
    p[0] = alpha;

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];
      p[1] = beta;

      for ( order = 1; order <= order_max; order++ )
      {
        std::cout << "\n";
        std::cout << "  Order = " << order << "\n";
        std::cout << "  ALPHA = " << alpha << "\n";
        std::cout << "  BETA = "  << beta << "\n";

        w = new double[order];
        x = new double[order];

        webbur::jacobi_compute_np ( order, np, p, x, w );

        for ( i = 0; i < order; i =i + 1 )
        {
          std::cout << "  " << std::setw(8) << i
               << "  " << std::setprecision(16) << std::setw(24) << x[i]
               << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
        }
        delete [] w;
        delete [] x;
      }
    }
  }

  delete [] p;

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test19_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19_NP tests LAGUERRE_COMPUTE_NP against LAGUERRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
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
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST19\n";
  std::cout << "  LAGUERRE_COMPUTE_NP computes a Gauss-Laguerre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.\n";
  std::cout << "\n";
  std::cout << "  LAGUERRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::laguerre_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::laguerre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test20_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20_NP tests LAGUERRE_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST20\n";
  std::cout << "  LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::laguerre_compute_np ( order, np, p, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test21_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21_NP tests LEGENDRE_COMPUTE_NP against LEGENDRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
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
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST21\n";
  std::cout << "  LEGENDRE_COMPUTE_NP computes a Gauss-Legendre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n";
  std::cout << "\n";
  std::cout << "  LEGENDRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = 2*ORDER-1\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::legendre_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = webbur::legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test22_np ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22_NP tests LEGENDRE_COMPUTE_NP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST22\n";
  std::cout << "  LEGENDRE_COMPUTE_NP computes a Gauss-Legendre rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) dx.\n";

  p = new double[np];

  for ( order = 1; order <= order_max; order++ )
  {
    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";

    w = new double[order];
    x = new double[order];

    webbur::legendre_compute_np ( order, np, p, x, w );

    for ( i = 0; i < order; i =i + 1 )
    {
      std::cout << "  " << std::setw(8) << i
           << "  " << std::setprecision(16) << std::setw(24) << x[i]
           << "  " << std::setprecision(16) << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  delete [] p;

  return;
}
//****************************************************************************80

void test23 ( int r )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests LEVEL_GROWTH_TO_ORDER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int R, the index of the rule to be examined.
//
{
  int dim;
  int dim_num = 11;
  int g;
  int *growth;
  int *level;
  int *order;
  int *rule;

  std::cout << "\n";
  std::cout << "TEST23\n";
  std::cout << "  LEVEL_GROWTH_TO_ORDER uses Level, Growth and Rule\n";
  std::cout << "  to determine the orders of each entry of a vector of 1D rules.\n";
  std::cout << "\n";
  std::cout << "  Here we examine rule " << r << ".\n";
  std::cout << "\n";
  std::cout <<
    "       LEVEL:0     1     2     3     4     5     6     7     8     9    10\n";
  std::cout << "GROWTH\n";

  growth = new int[dim_num];
  level = new int[dim_num];
  order = new int[dim_num];
  rule = new int[dim_num];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    rule[dim] = r;
  }

  for ( g = 0; g <= 6; g++ )
  {
    if ( r == 3 || r == 10 )
    {
      if ( g == 1 || g == 2 || g == 3 )
      {
        continue;
      }
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      growth[dim] = g;
    }
    webbur::level_growth_to_order ( dim_num, level, rule, growth, order );

    std::cout << "  " << std::setw(4) << g << "  ";

    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(4) << order[dim];
    }
    std::cout << "\n";
  }

  delete [] growth;
  delete [] level;
  delete [] order;
  delete [] rule;

  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests LEVEL_TO_ORDER_DEFAULT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_num = 11;
  int *level;
  int *order;
  int r;
  int *rule;

  std::cout << "\n";
  std::cout << "TEST24\n";
  std::cout << "  LEVEL_TO_ORDER_DEFAULT uses a default rule to\n";
  std::cout << "  determine the order of a rule from its level.\n";
  std::cout << "\n";
  std::cout <<
    "RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10\n";
  std::cout << "\n";

  level = new int[dim_num];
  order = new int[dim_num];
  rule = new int[dim_num];

  for ( r = 1; r <= 16; r++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      rule[dim] = r;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    webbur::level_to_order_default ( dim_num, level, rule, order );

    std::cout << "  " << std::setw(4) << r << "  ";

    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(4) << order[dim];
    }
    std::cout << "\n";
  }

  delete [] level;
  delete [] order;
  delete [] rule;

  return;
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests LEVEL_TO_ORDER_EXPONENTIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 December 2009
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_num = 11;
  int *level;
  int *order;
  int r;
  int *rule;

  std::cout << "\n";
  std::cout << "TEST25\n";
  std::cout << "  LEVEL_TO_ORDER_EXPONENTIAL uses an exponential rule to\n";
  std::cout << "  determine the order of a rule from its level.\n";
  std::cout << "\n";
  std::cout <<
    "RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10\n";
  std::cout << "\n";

  level = new int[dim_num];
  order = new int[dim_num];
  rule = new int[dim_num];

  for ( r = 1; r <= 10; r++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      rule[dim] = r;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    webbur::level_to_order_exponential ( dim_num, level, rule, order );

    std::cout << "  " << std::setw(4) << r << "  ";

    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(4) << order[dim];
    }
    std::cout << "\n";
  }

  delete [] level;
  delete [] order;
  delete [] rule;

  return;
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests LEVEL_TO_ORDER_EXPONENTIAL_SLOW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 December 2009
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_num = 11;
  int *level;
  int *order;
  int r;
  int *rule;

  std::cout << "\n";
  std::cout << "TEST26\n";
  std::cout << "  LEVEL_TO_ORDER_EXPONENTIAL_SLOW uses a slow exponential rule to\n";
  std::cout << "  determine the order of a rule from its level.\n";
  std::cout << "\n";
  std::cout << "  Since it is really only useful for fully nested rules,\n";
  std::cout << "  we only consider rules 11, 12 and 13.\n";
  std::cout << "\n";
  std::cout <<
    "RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10\n";
  std::cout << "\n";

  level = new int[dim_num];
  order = new int[dim_num];
  rule = new int[dim_num];

  for ( r = 11; r <= 13; r++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      rule[dim] = r;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    webbur::level_to_order_exponential_slow ( dim_num, level, rule, order );

    std::cout << "  " << std::setw(4) << r << "  ";

    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(4) << order[dim];
    }
    std::cout << "\n";
  }

  delete [] level;
  delete [] order;
  delete [] rule;

  return;
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests LEVEL_TO_ORDER_LINEAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 December 2009
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_num = 11;
  int *level;
  int *order;
  int r;
  int *rule;

  std::cout << "\n";
  std::cout << "TEST27\n";
  std::cout << "  LEVEL_TO_ORDER_LINEAR uses a linear rule to\n";
  std::cout << "  determine the order of a rule from its level.\n";
  std::cout << "\n";
  std::cout <<
    "RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10\n";
  std::cout << "\n";

  level = new int[dim_num];
  order = new int[dim_num];
  rule = new int[dim_num];

  for ( r = 1; r <= 10; r++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      rule[dim] = r;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    webbur::level_to_order_linear ( dim_num, level, rule, order );

    std::cout << "  " << std::setw(4) << r << "  ";

    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(4) << order[dim];
    }
    std::cout << "\n";
  }

  delete [] level;
  delete [] order;
  delete [] rule;

  return;
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests PATTERSON_LOOKUP against LEGENDRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2011
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
  int level;
  int level_max = 8;
  int n;
  int order;
  int p;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST28\n";
  std::cout << "  PATTERSON_LOOKUP computes a Gauss-Patterson rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n";
  std::cout << "\n";
  std::cout << "  LEGENDRE_INTEGRAL determines the exact value of\n";
  std::cout << "  this integal when f(x) = x^n.\n";
  std::cout << "\n";
  std::cout << "  A rule of order ORDER should be exact for monomials X^N\n";
  std::cout << "  up to N = (3*ORDER+1)/2\n";
  std::cout << "\n";
  std::cout << "  In the following table, for each order, the LAST THREE estimates\n";
  std::cout << "  are made on monomials that exceed the exactness limit for the rule.\n";
  std::cout << "\n";
  std::cout << "     Order         N       Estimate       Exact            Error\n";

  for ( level = 0; level <= level_max; level++ )
  {
    order = webbur::i4_power ( 2, level + 1 ) - 1;

    std::cout << "\n";

    f = new double[order];
    w = new double[order];
    x = new double[order];

    webbur::patterson_lookup ( order, x, w );

    if ( order == 1 )
    {
      p = 1;
    }
    else
    {
      p = ( 3 * order + 1 ) / 2;
    }
//
//  Truncate at 50.
//
    if ( 50 < p )
    {
      p = 50;
    }

    for ( n = 0; n <= p + 3; n = n + 1 )
    {
      exact = webbur::legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = std::pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = webbur::r8_abs ( exact - estimate );

      std::cout << "  " << std::setw(8)  << order
           << "  " << std::setw(8)  << n
           << "  " << std::setprecision(6) << std::setw(14) << estimate
           << "  " << std::setprecision(6) << std::setw(14) << exact
           << "  " << std::setprecision(6) << std::setw(14) << error << "\n";
    }
    delete [] f;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void test285 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST285 tests POINT_RADIAL_TOL_UNIQUE_INDEX_INC1, INC2 and INC3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 2
# define N1 11
# define N2 8

  int m = M;
  int n1 = N1;
  int n2 = N2;

  double a1[M*N1] = {
    0.0, 0.0,
    0.5, 0.5,
    1.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    0.500000001, 0.5,
    0.0, 0.0,
    0.0, 0.5,
    0.5, 0.0,
    1.0, 0.5,
    0.5, 1.0 };
  double a2[M*N2] = {
    0.4999999999, 0.5,
    0.75,        0.25,
    0.500000001, 0.9999999999,
    0.500000001, 0.0000000001,
    0.25,        0.75,
    0.75,        0.25,
    0.250000001, 0.7499999999,
    0.75,        0.75 };
  double *a3;
  int i1;
  int i2;
  int i3;
  int *indx1;
  int *indx2;
  int *indx3;
  int j;
  int n3;
  double *r1;
  double *r2;
  double *r3;
  int seed;
  double tol;
  int undx_value;
  int *undx1;
  int *undx2;
  int *undx3;
  int unique_num1;
  int unique_num2;
  int unique_num3;
  bool *unique1;
  bool *unique2;
  bool *unique3;
  int *xdnu1;
  int *xdnu2;
  int *xdnu3;
  double *z;

  seed = 123456789;

  std::cout << "\n";
  std::cout << "TEST285\n";
  std::cout << "  POINT_RADIAL_TOL_UNIQUE_INDEX_INC1 can index unique\n";
  std::cout << "  points in a \"permanent\" point set;\n";
  std::cout << "  POINT_RADIAL_TOL_UNIQUE_INDEX_INC2 can incremented by\n";
  std::cout << "  \"temporary\" points.\n";
  std::cout << "  POINT_RADIAL_TOL_UNIQUE_INDEX_INC3 can merge permanent\n";
  std::cout << "  and temporary points.\n";

  tol = std::sqrt ( webbur::r8_epsilon ( ) );
  std::cout << "  Using tolerance TOL = " << tol << "\n";
//
//  Step 1
//
  indx1 = new int[n1];
  r1 = new double[n1];
  undx1 = new int[n1];
  unique1 = new bool[n1];
  xdnu1 = new int[n1];
  z = new double[m];

  webbur::point_radial_tol_unique_index_inc1 ( m, n1, a1, tol, &seed, z, r1,
    indx1, unique1, &unique_num1, undx1, xdnu1 );

  std::cout << "\n";
  std::cout << "  UNIQUE_NUM1 = " << unique_num1 << "\n";
  std::cout << "  Expected   =  " << 9 << "\n";

  std::cout << "\n";
  std::cout << "  Item I1, unique index XDNU1[I1], representative location UNDX1[XDNU1[I1]]:\n";
  std::cout << "\n";
  std::cout << "           I1  XDNU1  UNDX1\n";
  std::cout << "\n";
  for ( i1 = 0; i1 < n1; i1++ )
  {
    std::cout << "  " << "    "
              << "  " << std::setw(4) << i1
              << "  " << std::setw(4) << xdnu1[i1]
              << "  " << std::setw(4) << undx1[xdnu1[i1]] << "\n";
  }

  std::cout << "\n";
  std::cout << "  Unique item I1, location UNDX1(I1), value A1(:,UNDX1(I1)):\n";
  std::cout << "\n";
  std::cout << "          I1 UNDX1  --A1(1,*)---  --A1(2,*)---\n";
  std::cout << "\n";
  for ( i1 = 0; i1 < unique_num1; i1++ )
  {
    std::cout << "  " << "    "
              << "  " << std::setw(4) << i1
              << "  " << std::setw(4) << undx1[i1]
              << "  " << std::setw(12) << a1[0+undx1[i1]*m]
              << "  " << std::setw(12) << a1[1+undx1[i1]*m] << "\n";
  }
  std::cout << "\n";
  std::cout << "  Unique item I1, location UNDX1(I1), value A1(:,UNDX1(I1)):\n";
  std::cout << "\n";
  std::cout << "          I1 UNIQUE1       R1     --A1(1,I1)--  --A1(2,I1)--\n";
  std::cout << "\n";
  for ( i1 = 0; i1 < n1; i1++ )
  {
    std::cout << "  " << "    "
              << "  " << std::setw(4) << i1
              << "  " << std::setw(4) << unique1[i1]
              << "  " << std::setw(12) << r1[i1]
              << "  " << std::setw(12) << a1[0+i1*m]
              << "  " << std::setw(12) << a1[1+i1*m] << "\n";
  }
  std::cout << "\n";
  std::cout << "          I1   INDX1   R1(INDX1)\n";
  std::cout << "\n";
  for ( i1 = 0; i1 < n1; i1++ )
  {
    std::cout << "  " << "    "
              << "  " << std::setw(4) << i1
              << "  " << std::setw(4) << indx1[i1]
              << "  " << std::setw(12) << r1[indx1[i1]] << "\n";
  }
//
//  Step 2
//
  indx2 = new int[n2];
  r2 = new double[n2];
  undx2 = new int[n2];
  unique2 = new bool[n2];
  xdnu2 = new int[n2];

  webbur::point_radial_tol_unique_index_inc2 ( m, n1, a1, n2, a2, tol, z,
    r1, indx1, unique1,  unique_num1, undx1, xdnu1,
    r2, indx2, unique2, &unique_num2, undx2, xdnu2 );

  std::cout << "\n";
  std::cout << "  UNIQUE_NUM2 = " << unique_num2 << "\n";
  std::cout << "  Expected   =  " << 3 << "\n";

  std::cout << "\n";
  std::cout << "  Item I2, unique index XDNU2[I2], representative location UNDX2[XDNU2[I2]]:\n";
  std::cout << "\n";
  std::cout << "    I2 XDNU2 UNDX2\n";
  std::cout << "\n";
  std::cout << "  (Temporary data)\n";
  std::cout << "\n";
  for ( i2 = 0; i2 < n2; i2++ )
  {
    if ( xdnu2[i2] < unique_num1 )
    {
      undx_value = undx1[xdnu2[i2]];
    }
    else
    {
      undx_value = undx2[ xdnu2[i2] - unique_num1 ];
    }
    std::cout << "  " << std::setw(4) << i2
              << "  " << std::setw(4) << i2 + n1
              << "  " << std::setw(4) << xdnu2[i2]
              << "  " << std::setw(4) << undx_value << "\n";
  }
  std::cout << "\n";
  std::cout << "  Unique item I2, location UNDX2(I2), value A2(:,UNDX2(I2)):\n";
  std::cout << "\n";
  std::cout << "    I2 UNDX2  --A2(1,*)---  --A2(2,*)---\n";
  std::cout << "\n";
  std::cout << "  (Temporary data)\n";
  std::cout << "\n";
  for ( i2 = 0; i2 < unique_num2; i2++ )
  {
    std::cout << "  " << std::setw(4) << i2
              << "  " << std::setw(4) << i2 + unique_num1
              << "  " << std::setw(4) << undx2[i2]
              << "  " << std::setw(12) << a2[0+(undx2[i2]-n1)*m]
              << "  " << std::setw(12) << a2[1+(undx2[i2]-n1)*m] << "\n";
  }
//
//  Step 3.
//
  std::cout << "\n";
  std::cout << "  Merge the temporary data with the permanent data.\n";

  a3 = new double[m*(n1+n2)];
  indx3 = new int[n1+n2];
  r3 = new double[n1+n2];
  undx3 = new int[n1+n2];
  unique3 = new bool[n1+n2];
  xdnu3 = new int[n1+n2];

  webbur::point_radial_tol_unique_index_inc3 ( m,
     n1, a1, r1, indx1, unique1,  unique_num1, undx1, xdnu1,
     n2, a2, r2, indx2, unique2,  unique_num2, undx2, xdnu2,
    &n3, a3, r3, indx3, unique3, &unique_num3, undx3, xdnu3 );

  std::cout << "\n";
  std::cout << "  UNIQUE_NUM3 = " << unique_num3 << "\n";
  std::cout << "  Expected   =  12\n";

  std::cout << "\n";
  std::cout << "  Item I3, unique index XDNU3[I3], representative location UNDX3[XDNU3[I3]]:\n";
  std::cout << "\n";
  std::cout << "          I3 XDNU3 UNDX3\n";
  std::cout << "\n";
  for ( i3 = 0; i3 < n3; i3++ )
  {
    std::cout << "  " << "    "
              << "  " << std::setw(4) << i3
              << "  " << std::setw(4) << xdnu3[i3]
              << "  " << std::setw(4) << undx3[xdnu3[i3]] << "\n";
  }
  std::cout << "\n";
  std::cout << "  Unique item I3, location UNDX3(I3), value A3(:,UNDX3(I3)):\n";
  std::cout << "\n";
  std::cout << "          I3 UNDX3  --A3(1,*)---  --A3(2,*)---\n";
  std::cout << "\n";
  for ( i3 = 0; i3 < unique_num3; i3++ )
  {
    std::cout << "  " << "    "
              << "  " << std::setw(4) << i3
              << "  " << std::setw(4) << undx3[i3]
              << "  " << std::setw(12) << a3[0+undx3[i3]*m]
              << "  " << std::setw(12) << a3[1+undx3[i3]*m] << "\n";
  }
  std::cout << "\n";
  std::cout << "  Unique item I3, location UNDX3(I3), value A3(:,UNDX3(I3)):\n";
  std::cout << "\n";
  std::cout << "          I3 UNIQUE3       R3     --A3(1,I3)--  --A3(2,I3)--\n";
  std::cout << "\n";
  for ( i3 = 0; i3 < n3; i3++ )
  {
    std::cout << "  " << "    "
              << "  " << std::setw(4) << i3
              << "  " << std::setw(4) << unique3[i3]
              << "  " << std::setw(12) << r3[i3]
              << "  " << std::setw(12) << a3[0+i3*m]
              << "  " << std::setw(12) << a3[1+i3*m] << "\n";
  }
  std::cout << "\n";
  std::cout << "          I3   INDX3   R3(INDX3)\n";
  std::cout << "\n";
  for ( i3 = 0; i3 < n3; i3++ )
  {
    std::cout << "  " << "    "
              << "  " << std::setw(4) << i3
              << "  " << std::setw(4) << indx3[i3]
              << "  " << std::setw(12) << r3[indx3[i3]] << "\n";
  }

  delete [] a3;
  delete [] indx1;
  delete [] indx2;
  delete [] indx3;
  delete [] r1;
  delete [] r2;
  delete [] r3;
  delete [] undx1;
  delete [] undx2;
  delete [] undx3;
  delete [] unique1;
  delete [] unique2;
  delete [] unique3;
  delete [] xdnu1;
  delete [] xdnu2;
  delete [] xdnu3;
  delete [] z;

  return;
# undef M
# undef N1
# undef N2
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests R8COL_TOL_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  double tol;
  int *undx;
  int unique_num;
  int *xdnu;

  std::cout << "\n";
  std::cout << "TEST29\n";
  std::cout << "  R8COL_TOL_UNDEX produces index vectors which create a sorted\n";
  std::cout << "  list of the tolerably unique columns of an R8COL,\n";
  std::cout << "  and a map from the original R8COL to the (implicit)\n";
  std::cout << "  R8COL of sorted tolerably unique elements.\n";

  webbur::r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  tol = 0.25;

  std::cout << "\n";
  std::cout << "  Using tolerance = " << tol << "\n";

  n_unique = webbur::r8col_tol_unique_count ( m, n, a, tol );

  std::cout << "\n";
  std::cout << "  Number of tolerably unique columns is " << n_unique << "\n";

  au = new double[m*n_unique];
  undx = new int[n_unique];
  xdnu = new int[n];

  webbur::r8col_tol_undex ( m, n, a, n_unique, tol, undx, xdnu );

  std::cout << "\n";
  std::cout << "  XDNU points to the representative for each item.\n";
  std::cout << "  UNDX selects the representatives.\n";
  std::cout << "\n";
  std::cout << "     I  XDNU  UNDX\n";
  std::cout << "\n";
  for ( i = 0; i < n_unique; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(4) << xdnu[i]
              << "  " << std::setw(4) << undx[i] << "\n";
  }
  for ( i = n_unique; i < n; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(4) << xdnu[i] << "\n";
  }

  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  webbur::r8mat_transpose_print ( m, n_unique, au,
    "  The tolerably unique R8COL (transposed):" );

  delete [] au;
  delete [] undx;
  delete [] xdnu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests R8VEC_SORT_HEAP_INDEX_A_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int i;
  int *indx;
  int n = 20;
  int seed;

  std::cout << "\n";
  std::cout << "TEST30\n";
  std::cout << "  R8VEC_SORT_HEAP_INDEX_A_NEW creates an ascending\n";
  std::cout << "  sort index for a R8VEC.\n";

  seed = 123456789;

  a = webbur::r8vec_uniform_01_new ( n, &seed );

  webbur::r8vec_print ( n, a, "  The unsorted array:" );

  indx = webbur::r8vec_sort_heap_index_a_new ( n, a );

  webbur::i4vec_print ( n, indx, "  The index vector:" );

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = a[indx[i]];
  }

  webbur::r8vec_print ( n, b, "  The sorted array A(INDX(:)):" );

  delete [] a;
  delete [] b;
  delete [] indx;

  return;
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests HCE_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 March 2011
//
//  Author:
//
//    John Burkardt
//
{
# define M 11
# define N ( 2 * M )

  int i;
  int m = M;
  int n = N;
  double r[N];
  double w[N];
  double x[N];

  std::cout << "\n";
  std::cout << "TEST31:\n";
  std::cout << "  HCE_COMPUTE returns a quadrature rule\n";
  std::cout << "  for piecewise Hermite cubic splines which are based\n";
  std::cout << "  on equally spaced function and derivative data.\n";
  std::cout << "\n";
  std::cout << "  Here we compute a rule of order N = " << n << "\n";

  webbur::hce_compute ( n, x, w );

  std::cout << "\n";
  std::cout << "     I        X(I)        W(I)\n";
  std::cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(10) << x[i]
              << "  " << std::setw(10) << w[i] << "\n";
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests HCC_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2011
//
//  Author:
//
//    John Burkardt
//
{
# define M 11
# define N ( 2 * M )

  int i;
  int m = M;
  int n = N;
  double r[N];
  double w[N];
  double x[N];

  std::cout << "\n";
  std::cout << "TEST32:\n";
  std::cout << "  HCC_COMPUTE returns a quadrature rule\n";
  std::cout << "  for piecewise Hermite cubic splines which are based\n";
  std::cout << "  on Chebyshev-spaced function and derivative data.\n";
  std::cout << "\n";
  std::cout << "  Here we compute a rule of order N = " << n << "\n";

  webbur::hcc_compute ( n, x, w );

  std::cout << "\n";
  std::cout << "     I        X(I)        W(I)\n";
  std::cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(10) << x[i]
              << "  " << std::setw(10) << w[i] << "\n";
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests HC_COMPUTE_WEIGHTS_FROM_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  double dn;
  double fn;
  int j;
  int n = N;
  double q;
  double q_exact;
  double *r;
  double s;
  int seed;
  double t;
  int test;
  double *w;
  double x[N];

  std::cout << "\n";
  std::cout << "TEST33:\n";
  std::cout << "  HC_COMPUTE_WEIGHTS_FROM_POINTS returns quadrature weights\n";
  std::cout << "  given the points.\n";

  seed = 123456789;

  for ( test = 1; test <= 3; test++ )
  {
    r = webbur::r8vec_uniform_01_new ( n, &seed );

    x[0] = r[0];
    for ( j = 1; j < n; j++ )
    {
      x[j] = x[j-1] + r[j];
    }
    delete [] r;

    std::cout << "\n";
    std::cout << "  Trial #" << test << ":\n";
    std::cout << "  Random spacing\n";
    std::cout << "  Number of points N = " << n << "\n";
    std::cout << "  Interval = [" << x[0] << ", " << x[n-1] << "]\n";

    w = new double [ 2 * n ];

    webbur::hc_compute_weights_from_points ( n, x, w );

    q = 0.0;

    for ( j = 0; j < n; j++ )
    {
      cubic_value ( x[j], &fn, &dn, &s, &t );
      q = q + w[0+j*2] * fn + w[1+j*2] * dn;
    }

    q_exact = cubic_integrate ( x[0], x[n-1] );

    std::cout << "\n";
    std::cout << "  Q         = " << q << "\n";
    std::cout << "  Q (exact) = " << q_exact << "\n";

    delete [] w;
  }
  return;
# undef N
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 uses HERMITE_INTERPOLANT on the Runge function using equally spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double max_dif;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double xhi;
  double xlo;
  double xt;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double *ysp;
  double yt;

  std::cout << "\n";
  std::cout << "TEST34\n";
  std::cout << "  HERMITE_INTERPOLANT computes the Hermite interpolant to data.\n";
  std::cout << "  Here, f(x) is the Runge function\n";
  std::cout << "  and the data is evaluated at equally spaced points.\n";
  std::cout << "  As N increases, the maximum error grows.\n";
  std::cout << "\n";
  std::cout << "     N     Max | F(X) - H(F(X)) |\n";
  std::cout << "\n";

  for ( n = 3; n <= 15; n = n + 2 )
  {
    y = new double[n];
    yp = new double[n];

    nd = 2 * n;
    xd = new double[nd];
    yd = new double[nd];

    ndp = 2 * n - 1;
    xdp = new double[ndp];
    ydp = new double[ndp];

    ns = 10 * ( n - 1 ) + 1;

    xlo = -5.0;
    xhi = +5.0;
    x = webbur::r8vec_linspace_new ( n, xlo, xhi );

    for ( i = 0; i < n; i++ )
    {
      y[i] = 1.0 / ( 1.0 + x[i] * x[i] );
      yp[i] = - 2.0 * x[i] / ( 1.0 + x[i] * x[i] ) / ( 1.0 + x[i] * x[i] );
    }

    webbur::hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
//
//  Compare exact and interpolant at sample points.
//
    xs = webbur::r8vec_linspace_new ( ns, xlo, xhi );

    ys = new double[ns];
    ysp = new double[ns];

    webbur::hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp );

    max_dif = 0.0;
    for ( i = 0; i < ns; i++ )
    {
      xt = xs[i];
      yt = 1.0 / ( 1.0 + xt * xt );
      max_dif = webbur::r8_max ( max_dif, webbur::r8_abs ( ys[i] - yt ) );
    }

    std::cout << "  " << std::setw(4) << n
              << "  " << std::setw(14) << max_dif << "\n";

    delete [] x;
    delete [] xd;
    delete [] xdp;
    delete [] xs;
    delete [] y;
    delete [] yd;
    delete [] ydp;
    delete [] yp;
    delete [] ys;
    delete [] ysp;
  }

  return;
}
//****************************************************************************80

void test35 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST35 uses HERMITE_INTERPOLANT on the Runge function using Chebyshev spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double max_dif;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double xhi;
  double xlo;
  double xt;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double *ysp;
  double yt;

  std::cout << "\n";
  std::cout << "TEST35\n";
  std::cout << "  HERMITE_INTERPOLANT computes the Hermite interpolant to data.\n";
  std::cout << "  Here, f(x) is the Runge function\n";
  std::cout << "  and the data is evaluated at Chebyshev spaced points.\n";
  std::cout << "  As N increases, the maximum error decreases.\n";
  std::cout << "\n";
  std::cout << "     N     Max | F(X) - H(F(X)) |\n";
  std::cout << "\n";

  for ( n = 3; n <= 15; n = n + 2 )
  {
    y = new double[n];
    yp = new double[n];

    nd = 2 * n;
    xd = new double[nd];
    yd = new double[nd];

    ndp = 2 * n - 1;
    xdp = new double[ndp];
    ydp = new double[ndp];

    ns = 10 * ( n - 1 ) + 1;

    xlo = -5.0;
    xhi = +5.0;
    x = webbur::r8vec_chebyshev_new ( n, xlo, xhi );

    for ( i = 0; i < n; i++ )
    {
      y[i] = 1.0 / ( 1.0 + x[i] * x[i] );
      yp[i] = - 2.0 * x[i] / ( 1.0 + x[i] * x[i] ) / ( 1.0 + x[i] * x[i] );
    }

    webbur::hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
//
//  Compare exact and interpolant at sample points.
//
    xs = webbur::r8vec_linspace_new ( ns, xlo, xhi );

    ys = new double[ns];
    ysp = new double[ns];

    webbur::hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp );

    max_dif = 0.0;
    for ( i = 0; i < ns; i++ )
    {
      xt = xs[i];
      yt = 1.0 / ( 1.0 + xt * xt );
      max_dif = webbur::r8_max ( max_dif, webbur::r8_abs ( ys[i] - yt ) );
    }

    std::cout << "  " << std::setw(4) << n
              << "  " << std::setw(14) << max_dif << "\n";

    delete [] x;
    delete [] xd;
    delete [] xdp;
    delete [] xs;
    delete [] y;
    delete [] yd;
    delete [] ydp;
    delete [] yp;
    delete [] ys;
    delete [] ysp;
  }

  return;
}
//****************************************************************************80

void test36 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST36 tests HERMITE_INTERPOLANT_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int e;
  double error;
  double exact;
  int i;
  int k;
  int n;
  double q;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST36:\n";
  std::cout << "  HERMITE_INTERPOLANT_RULE\n";
  std::cout << "  is given a set of N abscissas for a Hermite interpolant\n";
  std::cout << "  and returns N pairs of quadrature weights\n";
  std::cout << "  for function and derivative values at the abscissas.\n";
//
//  1: Behavior with increasing N.
//
  a = 0.0;
  b = 1.0;

  std::cout << "\n";
  std::cout << "  Observe behavior of quadrature weights for increasing N\n";
  std::cout << "  We are working in " << a << " <= X <= " << b << "\n";

  for ( n = 3; n <= 11; n = n + 2 )
  {
    x = webbur::r8vec_linspace_new ( n, a, b );
    w = new double[2*n];
    webbur::hermite_interpolant_rule ( n, a, b, x, w );

    std::cout << "\n";
    std::cout << "     I       X               W(F(X))        W(F'(X))\n";
    std::cout << "\n";
    k = 0;
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(4) << i
                << "  " << std::setw(14) << x[i]
                << "  " << std::setw(14) << w[k]
                << "  " << std::setw(14) << w[k+1] << "\n";
      k = k + 2;
    }
    delete [] x;
    delete [] w;
  }
//
//  2: Integral estimates with equally spaced points.
//
  a = - 5.0;
  b = + 5.0;
  n = 11;

  std::cout << "\n";
  std::cout << "  Use the rule with N = " << n << " to estimate integrals.\n";
  std::cout << "  Points are equally spaced.\n";
  std::cout << "  We are working in " << a << " <= X <= " << b << "\n";

  x = webbur::r8vec_linspace_new ( n, a, b );
  w = new double[2*n];
  webbur::hermite_interpolant_rule ( n, a, b, x, w );

  std::cout << "\n";
  std::cout << "     I       X               W(F(X))        W(F'(X))\n";
  std::cout << "\n";
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(14) << x[i]
              << "  " << std::setw(14) << w[k]
              << "  " << std::setw(14) << w[k+1] << "\n";
    k = k + 2;
  }

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  std::cout << "\n";
  std::cout << "  Estimate integral of 1 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  std::cout << "  Estimate integral of X = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  std::cout << "  Estimate integral of X^2 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
          - w[k+1] * 2.0 * x[i] / std::pow ( 1.0 + x[i] * x[i], 2 );
    k = k + 2;
  }
  std::cout << "  Estimate integral of 1/(1+x^2) = " << q << "\n";

  delete [] w;
  delete [] x;
//
//  3: Integral estimates with Chebyshev spaced points.
//
  a = - 5.0;
  b = + 5.0;
  n = 11;

  std::cout << "\n";
  std::cout << "  Use the rule with N = " << n << " to estimate integrals.\n";
  std::cout << "  Points are Chebyshev spaced.\n";
  std::cout << "  We are working in " << a << " <= X <= " << b << "\n";

  x = webbur::r8vec_chebyshev_new ( n, a, b );
  w = new double[2*n];
  webbur::hermite_interpolant_rule ( n, a, b, x, w );

  std::cout << "\n";
  std::cout << "     I       X               W(F(X))        W(F'(X))\n";
  std::cout << "\n";
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(14) << x[i]
              << "  " << std::setw(14) << w[k]
              << "  " << std::setw(14) << w[k+1] << "\n";
    k = k + 2;
  }

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  std::cout << "\n";
  std::cout << "  Estimate integral of 1 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  std::cout << "  Estimate integral of X = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  std::cout << "  Estimate integral of X^2 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
          - w[k+1] * 2.0 * x[i] / std::pow ( 1.0 + x[i] * x[i], 2 );
    k = k + 2;
  }
  std::cout << "  Estimate integral of 1/(1+x^2) = " << q << "\n";

  delete [] w;
  delete [] x;
//
//  4: Integral estimates with Legendre spaced points.
//
  a = - 5.0;
  b = + 5.0;
  n = 11;

  std::cout << "\n";
  std::cout << "  Use the rule with N = " << n << " to estimate integrals.\n";
  std::cout << "  Points are Legendre spaced.\n";
  std::cout << "  We are working in " << a << " <= X <= " << b << "\n";

  x = webbur::r8vec_legendre_new ( n, a, b );
  w = new double[2*n];
  webbur::hermite_interpolant_rule ( n, a, b, x, w );

  std::cout << "\n";
  std::cout << "     I       X               W(F(X))        W(F'(X))\n";
  std::cout << "\n";
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(14) << x[i]
              << "  " << std::setw(14) << w[k]
              << "  " << std::setw(14) << w[k+1] << "\n";
    k = k + 2;
  }

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  std::cout << "\n";
  std::cout << "  Estimate integral of 1 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  std::cout << "  Estimate integral of X = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  std::cout << "  Estimate integral of X^2 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
          - w[k+1] * 2.0 * x[i] / std::pow ( 1.0 + x[i] * x[i], 2 );
    k = k + 2;
  }
  std::cout << "  Estimate integral of 1/(1+x^2) = " << q << "\n";

  delete [] w;
  delete [] x;
//
//  5: Integral estimates with Chebyshev spaced points on 1/(1+x^2), increasing N.
//
  a = - 5.0;
  b = + 5.0;

  std::cout << "\n";
  std::cout << "  Approximate integral of 1/(1+x^2) with increasing N.\n";
  std::cout << "  Points are Chebyshev spaced.\n";
  std::cout << "  We are working in " << a << " <= X <= " << b << "\n";
  std::cout << "\n";
  std::cout << "     N     Estimate         Error\n";
  std::cout << "\n";

  for ( n = 2; n <= 11; n++ )
  {
    x = webbur::r8vec_chebyshev_new ( n, a, b );
    w = new double[2*n];
    webbur::hermite_interpolant_rule ( n, a, b, x, w );

    q = 0.0;
    k = 0;
    for ( i = 0; i < n; i++ )
    {
      q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
            - w[k+1] * 2.0 * x[i] / std::pow ( 1.0 + x[i] * x[i], 2 );
      k = k + 2;
    }
    exact = std::atan ( b ) - std::atan ( a );
    error = webbur::r8_abs ( q - exact );
    std::cout << "  " << std::setw(4) << n
              << "  " << std::setw(14) << q
              << "  " << std::setw(14) << error << "\n";

    delete [] w;
    delete [] x;
  }
//
//  6: Integral estimates, with Chebyshev spaced points, for monomials, using N = 11.
//
  a = -1.0;
  b = 1.0;

  std::cout << "\n";
  std::cout << "  Compute integral estimates for monomials X^0 through X^15.\n";
  std::cout << "  Use N = 5, 9, 13, 17, 21 point rules.\n";
  std::cout << "  Points are Chebyshev spaced.\n";
  std::cout << "  We are working in " << a << " <= X <= " << b << "\n";

  for ( n = 5; n <= 21; n = n + 4 )
  {
    x = webbur::r8vec_chebyshev_new ( n, a, b );
    w = new double[2*n];
    webbur::hermite_interpolant_rule ( n, a, b, x, w );

    std::cout << "\n";
    std::cout << "  Estimates are made using N = " << n << "\n";
    std::cout << "  F(X)         Integral        Estimate           Error\n";
    std::cout << "\n";
    for ( e = 0; e <= 15; e++ )
    {
      q = 0.0;
      k = 0;
      for ( i = 0; i < n; i++ )
      {
        if ( e == 0 )
        {
          q = q + w[k];
        }
        else
        {
          q = q + w[k] * std::pow ( x[i], e ) + w[k+1] * e * std::pow ( x[i], e - 1 );
        }
        k = k + 2;
      }
      exact = ( std::pow ( b, e + 1 ) - std::pow ( a, e + 1 ) ) / ( double ) ( e + 1 );
      std::cout << "  X^" << std::setw(7) << e
                << "  " << std::setw(14) << exact
                << "  " << std::setw(14) << q
                << "  " << std::setw(14) << webbur::r8_abs ( exact - q ) << "\n";
    }

    q = 0.0;
    k = 0;
    for ( i = 0; i < n; i++ )
    {
      q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
            - w[k+1] * 2.0 * x[i] / std::pow ( 1.0 + x[i] * x[i], 2 );
      k = k + 2;
    }

    exact = std::atan ( b ) - std::atan ( a );
    std::cout << "  1/(1+x^2)"
              << "  " << std::setw(14) << exact
              << "  " << std::setw(14) << q
              << "  " << std::setw(14) << webbur::r8_abs ( exact - q ) << "\n";

    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test37 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST37 checks that the HGK weights are correctly scaled.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  int o;
  int order[8] = { 1, 3, 9, 19, 35, 37, 41, 43 };
  static double sqrtpi = 1.7724538509055159;
  int rule;
  double s;
  double *w;

  std::cout << "\n";
  std::cout << "TEST37\n";
  std::cout << "  HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS looks up weights\n";
  std::cout << "  for Genz-Keister quadrature rules for the Hermite weight function.\n";
  std::cout << "\n";
  std::cout << "  This test simply checks that, for each rule, the quadrature\n";
  std::cout << "  weights correctly sum to sqrt(pi).\n";
  std::cout << "\n";
  std::cout << " Index     Order         Sum of Weights\n";
  std::cout << "\n";

  for ( rule = 0; rule < 8; rule++ )
  {
    o = order[rule];

    w = new double[o];

    webbur::hermite_genz_keister_lookup_weights ( o, w );

    s = webbur::r8vec_sum ( o, w );

    std::cout << "  " << std::setw(4) << rule
              << "  " << std::setw(8) << o
              << "  " << std::setprecision(6) << std::setw(14) << s << "\n";

    delete [] w;
  }
  std::cout << "\n";
  std::cout << " Correct sum:            " << sqrtpi << "\n";

  return;
}
//****************************************************************************80

void test38 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST38 tabulates the Hermite interpolant and its derivative. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double *ysp;

  std::cout << "\n";
  std::cout << "TEST38\n";
  std::cout << "  HERMITE_INTERPOLANT sets up the Hermite interpolant.\n";
  std::cout << "  HERMITE_INTERPOLANT_VALUE evaluates it.\n";
  std::cout << "  Consider data for y=sin(x) at x=0,1,2,3,4.\n";

  n = 5;
  y = new double[n];
  yp = new double[n];

  nd = 2 * n;
  xd = new double[nd];
  yd = new double[nd];

  ndp = 2 * n - 1;
  xdp = new double[ndp];
  ydp = new double[ndp];

  x = webbur::r8vec_linspace_new ( n, 0.0, 4.0 );
  for ( i = 0; i < n; i++ )
  {
    y[i] = std::sin ( x[i] );
    yp[i] = std::cos ( x[i] );
  }

  webbur::hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
/*
  Now sample the interpolant at NS points, which include data values.
*/
  ns = 4 * ( n - 1 ) + 1;
  ys = new double[ns];
  ysp = new double[ns];

  xs = webbur::r8vec_linspace_new ( ns, 0.0, 4.0 );

  webbur::hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp );

  std::cout << "\n";
  std::cout << "  In the following table, there should be perfect\n";
  std::cout << "  agreement between F and H, and F' and H'\n";
  std::cout << "  at the data points X = 0, 1, 2, 3, and 4.\n";
  std::cout << "\n";
  std::cout << "  In between, H and H' approximate F and F'.\n";
  std::cout << "\n";
  std::cout << "     I       X(I)          F(X(I))         H(X(I)) ";
  std::cout << "        F'(X(I))        H'(X(I))\n";
  std::cout << "\n";
  for ( i = 0; i < ns; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(14) << xs[i]
              << "  " << std::setw(14) << std::sin ( xs[i] )
              << "  " << std::setw(14) << ys[i]
              << "  " << std::setw(14) << std::cos ( xs[i] )
              << "  " << std::setw(14) << ysp[i] << "\n";
  }

  delete [] x;
  delete [] xd;
  delete [] xdp;
  delete [] xs;
  delete [] y;
  delete [] yd;
  delete [] ydp;
  delete [] yp;
  delete [] ys;
  delete [] ysp;

  return;
}
//****************************************************************************80

void test39 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST39 tests the LEVEL_TO_ORDER_** functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  int growth;
  int level;
  int order;

  std::cout << "\n";
  std::cout << "TEST39:\n";
  std::cout << "  Test the various LEVEL_TO_ORDER_** functions,\n";
  std::cout << "  which, for a given type of quadrature rule, accept\n";
  std::cout << "  LEVEL, the index of the rule in the family, and\n";
  std::cout << "  GROWTH (0=slow, 1=moderate, 2=full), a growth rate, and\n";
  std::cout << "  return the appropriate corresponding order for the rule.\n";
  std::cout << "\n";
  std::cout << "  LEVEL_TO_ORDER_EXP_CC:\n";
  std::cout << "  Slow/moderate/full exponential growth typical of\n";
  std::cout << "  the Clenshaw Curtis rule.\n";
  std::cout << "\n";
  std::cout << "  LEVEL  ORDER  ORDER  ORDER\n";
  std::cout << "         G = 0  G = 1  G = 2\n";
  std::cout << "\n";

  for ( level = 0; level <= 8; level++ )
  {
    std::cout << "  " << std::setw(5) << level;
    for ( growth = 0; growth <= 2; growth++ )
    {
      order = webbur::level_to_order_exp_cc ( level, growth );
      std::cout << "  " << std::setw(5) << order;
    }
    std::cout << "\n";
  }

  std::cout << "\n";
  std::cout << "  LEVEL_TO_ORDER_EXP_F2:\n";
  std::cout << "  Slow/moderate/full exponential growth typical of\n";
  std::cout << "  the Fejer Type 2 rule.\n";
  std::cout << "\n";
  std::cout << "  LEVEL  ORDER  ORDER  ORDER\n";
  std::cout << "         G = 0  G = 1  G = 2\n";
  std::cout << "\n";

  for ( level = 0; level <= 8; level++ )
  {
    std::cout << "  " << std::setw(5) << level;
    for ( growth = 0; growth <= 2; growth++ )
    {
      order = webbur::level_to_order_exp_f2 ( level, growth );
      std::cout << "  " << std::setw(5) << order;
    }
    std::cout << "\n";
  }

  std::cout << "\n";
  std::cout << "  LEVEL_TO_ORDER_EXP_GAUSS:\n";
  std::cout << "  Slow/moderate/full exponential growth typical of\n";
  std::cout << "  a Gauss rule.\n";
  std::cout << "\n";
  std::cout << "  LEVEL  ORDER  ORDER  ORDER\n";
  std::cout << "         G = 0  G = 1  G = 2\n";
  std::cout << "\n";

  for ( level = 0; level <= 8; level++ )
  {
    std::cout << "  " << std::setw(5) << level;
    for ( growth = 0; growth <= 2; growth++ )
    {
      order = webbur::level_to_order_exp_gauss ( level, growth );
      std::cout << "  " << std::setw(5) << order;
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::cout << "  LEVEL_TO_ORDER_EXP_GP:\n";
  std::cout << "  Slow/moderate/full exponential growth typical of\n";
  std::cout << "  a Gauss-Patterson rule.\n";
  std::cout << "\n";
  std::cout << "  LEVEL  ORDER  ORDER  ORDER\n";
  std::cout << "         G = 0  G = 1  G = 2\n";
  std::cout << "\n";

  for ( level = 0; level <= 8; level++ )
  {
    std::cout << "  " << std::setw(5) << level;
    for ( growth = 0; growth <= 2; growth++ )
    {
      order = webbur::level_to_order_exp_gp ( level, growth );
      std::cout << "  " << std::setw(5) << order;
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::cout << "  LEVEL_TO_ORDER_EXP_HGK:\n";
  std::cout << "  Slow/moderate/full exponential growth typical of\n";
  std::cout << "  a Hermite Genz-Keister rule.\n";
  std::cout << "\n";
  std::cout << "  LEVEL  ORDER  ORDER  ORDER\n";
  std::cout << "         G = 0  G = 1  G = 2\n";
  std::cout << "\n";

  for ( level = 0; level <= 8; level++ )
  {
    std::cout << "  " << std::setw(5) << level;
    for ( growth = 0; growth <= 2; growth++ )
    {
      if ( growth == 2 && 5 < level )
      {
      }
      else
      {
        order = webbur::level_to_order_exp_hgk ( level, growth );
        std::cout << "  " << std::setw(5) << order;
      }
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::cout << "  LEVEL_TO_ORDER_LINEAR_NN:\n";
  std::cout << "  Slow/moderate linear growth typical of\n";
  std::cout << "  a non-nested Gauss rule.\n";
  std::cout << "\n";
  std::cout << "  LEVEL  ORDER  ORDER\n";
  std::cout << "         G = 0  G = 1\n";
  std::cout << "\n";

  for ( level = 0; level <= 8; level++ )
  {
    std::cout << "  " << std::setw(5) << level;
    for ( growth = 0; growth <= 1; growth++ )
    {
      order = webbur::level_to_order_linear_nn ( level, growth );
      std::cout << "  " << std::setw(5) << order;
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::cout << "  LEVEL_TO_ORDER_LINEAR_WN:\n";
  std::cout << "  Slow/moderate linear growth typical of\n";
  std::cout << "  a weakly-nested Gauss rule.\n";
  std::cout << "\n";
  std::cout << "  LEVEL  ORDER  ORDER\n";
  std::cout << "         G = 0  G = 1\n";
  std::cout << "\n";

  for ( level = 0; level <= 8; level++ )
  {
    std::cout << "  " << std::setw(5) << level;
    for ( growth = 0; growth <= 1; growth++ )
    {
      order = webbur::level_to_order_linear_wn ( level, growth );
      std::cout << "  " << std::setw(5) << order;
    }
    std::cout << "\n";
  }
  return;
}
//****************************************************************************80

double cubic_antiderivative ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    CUBIC_ANTIDERIVATIVE evaluates the antiderivative function of a cubic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double CUBIC_ANTIDERIVATIVE, the value.
//
{
  double value;

  value = x * x * ( 5.0 + x * ( - 7.0 / 3.0 + x * 1.0 / 4.0 ) );

  return value;
}
//****************************************************************************80

double cubic_integrate ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CUBIC_INTEGRATE integrates the cubic from A to B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the integration interval.
//
//    Output, double Q, the integral from A to B.
//
{
  double q;

  q = cubic_antiderivative ( b ) - cubic_antiderivative ( a );

  return q;
}
//****************************************************************************80

void cubic_value ( double x, double *f, double *d, double *s, double *t )

//****************************************************************************80
//
//  Purpose:
//
//    CUBIC_VALUE evaluates a cubic function.
//
//  Discussion:
//
//    f(x) =   x^3 -  7 x^2 + 10 x
//    d(x) = 3 x^2 - 14 x   + 10
//    s(x) = 6 x   - 14
//    t(x) = 6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double F, D, S, T, the value and first three
//    derivatives of the cubic function.
//
{
  *f = 0.0 + x * ( 10.0 + x * ( -  7.0 + x * 1.0 ) );
  *d =             10.0 + x * ( - 14.0 + x * 3.0 );
  *s =                          - 14.0 + x * 6.0;
  *t =                                       6.0;

  return;
}
