# include "sandia_rules.hpp"
# include "sandia_cubature.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

namespace webbur 
{
//****************************************************************************80

double c1_geg_monomial_integral ( double alpha, int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on C1.
//
//  Discussion:
//
//    C1_GEG is the interval [-1,+1] with the Gegenbauer weight function
//
//      w(alpha;x) = (1-x^2)^alpha
//
//    with -1.0 < alpha.
//
//    value = integral ( -1 <= x <= +1 ) x^expon (1-x^2)^alpha dx
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
//  Parameters:
//
//    Input, double ALPHA, the exponent of (1-X^2).
//    - 1.0 < ALPHA.
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double C1_GEG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double value;
  double value1;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "C1_GEG_MONOMIAL_INTEGRAL - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
    return value;
  }

  c = ( double ) ( expon );

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value1 = webbur::r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = 2.0 * webbur::r8_gamma ( 1.0 + c ) * webbur::r8_gamma ( 1.0 + alpha ) 
    * value1 / webbur::r8_gamma ( 2.0 + alpha + c );

  return value;
}
//****************************************************************************80

double c1_jac_monomial_integral ( double alpha, double beta, int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over C1.
//
//  Discussion:
//
//    value = integral ( -1 <= x <= +1 ) x^expon (1-x)^alpha (1+x)^beta dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the exponent of (1-X) in the weight factor.
//
//    Input, double BETA, the exponent of (1+X) in the weight factor.
//
//    Input, int EXPON, the exponent.
//
//    Output, double C1_JAC_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double s;
  double value;
  double value1;
  double value2;

  c = ( double ) ( expon );

  if ( ( expon % 2 ) == 0 )
  {
    s = +1.0;
  }
  else
  {
    s = -1.0;
  }

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + beta + c;
  arg4 = - 1.0;

  value1 = webbur::r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  arg1 = - beta;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value2 = webbur::r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = webbur::r8_gamma ( 1.0 + c ) * ( 
      s * webbur::r8_gamma ( 1.0 + beta  ) * value1 
    / webbur::r8_gamma ( 2.0 + beta  + c ) 
    +     webbur::r8_gamma ( 1.0 + alpha ) * value2 
    / webbur::r8_gamma ( 2.0 + alpha + c ) );

  return value;
}
//****************************************************************************80

double c1_leg_monomial_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on C1.
//
//  Discussion:
//
//    C1_LEG is the interval [-1,+1] with the Legendre weight function
//
//      w(x) = 1.
//
//    value = integral ( -1 <= x <= +1 ) x^expon dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double C1_LEG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double value;

  if ( expon < 0 )
  {
    std::cerr << "\n";
    std::cerr << "C1_LEG_MONOMIAL_INTEGRAL - Fatal error!\n";
    std::cerr << "  EXPON < 0.\n";
    std::exit ( 1 );
  }

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
    return value;
  }

  value = 2.0 / ( double) ( expon + 1 );

  return value;
}
//****************************************************************************80

void cn_geg_01_1 ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_01_1 implements a precision 1 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[[O], the weights.
//
{
  int expon;
  int i;
  int j;
  int k;
  double value2;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_GEG_01_1 - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  expon = 0;
  volume_1d = webbur::c1_geg_monomial_integral ( alpha, expon );
  volume = std::pow ( volume_1d, n );

  expon = 1;
  value2 = webbur::c1_geg_monomial_integral ( alpha, expon );

  webbur::r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = value2 / volume_1d;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int cn_geg_01_1_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_01_1_SIZE sizes a precision 1 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, int CN_GEG_01_1_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_GEG_01_1_SIZE - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void cn_geg_02_xiu ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_02_XIU implements the Xiu precision 2 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_GEG_02_XIU - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = webbur::r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = 1.0;
  delta0 = 0.0;
  c1 = 1.0 / ( 2.0 * alpha + 3.0 );

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( std::sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = webbur::c1_geg_monomial_integral ( alpha, expon );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int cn_geg_02_xiu_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_02_XIU_SIZE sizes the Xiu precision 2 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
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
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, int CN_GEG_02_XIU_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_GEG_02_XIU_SIZE - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  o = n + 1;

  return o;
}
//****************************************************************************80

void cn_geg_03_xiu ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_03_XIU implements the Xiu precision 3 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
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
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int expon;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_GEG_03_XIU - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( ( 2 * r - 1 ) * ( j + 1 ) ) * pi / ( double ) ( n );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg ) 
        / std::sqrt ( 2.0 * alpha + 3.0 );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg ) 
        / std::sqrt ( 2.0 * alpha + 3.0 );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = std::sqrt ( 2.0 ) * webbur::r8_mop ( j + 1 ) 
        / std::sqrt ( 2.0 * alpha + 3.0 );
      if ( n == 1 )
      {
        x[i+j*n] = x[i+j*n] / std::sqrt ( 2.0 );
      }
      i = i + 1;
    }
  }

  expon = 0;
  volume_1d = webbur::c1_geg_monomial_integral ( alpha, expon );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int cn_geg_03_xiu_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
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
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, int CN_GEG_03_XIU_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_GEG_03_XIU_SIZE - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  o = 2 * n;

  return o;
}
//****************************************************************************80

double cn_geg_monomial_integral ( int n, double alpha, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on CN.
//
//  Discussion:
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//    value = integral ( CN ) 
//      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i)^2)^alpha dx(i)
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of (1-X).
//    -1.0 < ALPHA.
//
//    Input, int EXPON[N], the exponents.
//
//    Output, double CN_GEG_MONOMIAL_INTEGRA, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_GEG_MONOMIAL_INTEGRAL - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = c1_geg_monomial_integral ( alpha, expon[i] );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void cn_jac_01_1 ( int n, double alpha, double beta, int o, double x[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_01_1 implements a precision 1 rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha. 
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int i;
  int k;
  double value2;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_JAC_01_1 - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_JAC_01_1 - Fatal error!\n";
    std::cerr << "  BETA <= -1.0\n";
    std::exit ( 1 );
  }

  expon = 0;
  volume_1d = webbur::c1_jac_monomial_integral ( alpha, beta, expon );
  volume = std::pow ( volume_1d, n );

  expon = 1;
  value2 = c1_jac_monomial_integral ( alpha, beta, expon );

  webbur::r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = value2 / volume_1d;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int cn_jac_01_1_size ( int n, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_01_1_SIZE sizes a precision 1 rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha. 
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Output, int CN_JAC_01_1_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_JAC_01_1_SIZE - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_JAC_01_1_SIZE - Fatal error!\n";
    std::cerr << "  BETA <= -1.0\n";
    std::exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void cn_jac_02_xiu ( int n, double alpha, double beta, int o, double x[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_02_XIU implements the Xiu precision 2 rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_JAC_02_XIU - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_JAC_02_XIU - Fatal error!\n";
    std::cerr << "  BETA <= -1.0\n";
    std::exit ( 1 );
  }

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = webbur::r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = ( alpha + beta + 2.0 ) / 2.0;
  delta0 = ( alpha - beta ) / 2.0;
  c1 = 2.0 * ( alpha + 1.0 ) * ( beta + 1.0 ) / ( alpha + beta + 3.0 )
    / ( alpha + beta + 2.0 );

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( std::sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = webbur::c1_jac_monomial_integral ( alpha, beta, expon );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int cn_jac_02_xiu_size ( int n, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_02_XIU_SIZE sizes the Xiu precision 2 rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Output, int CN_JAC_02_XIU_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_JAC_02_XIU_SIZE - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "CN_JAC_02_XIU_SIZE - Fatal error!\n";
    std::cerr << "  BETA <= -1.0\n";
    std::exit ( 1 );
  }

  o = n + 1;

  return o;
}
//****************************************************************************80

double cn_jac_monomial_integral ( int n, double alpha, double beta, 
  int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over CN.
//
//  Discussion:
//
//    value = integral ( CN ) 
//      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i))^alpha (1+x(i))^beta dx(i)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of (1-X) in the weight factor.
//
//    Input, double BETA, the exponent of (1+X) in the weight factor.
//
//    Input, int EXPON[N], the exponents.
//
//    Output, double CN_JAC_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = c1_jac_monomial_integral ( alpha, beta, expon[i] );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void cn_leg_01_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_01_1 implements the midpoint rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int k;
  double volume;
  double volume_1d;

  expon = 0;
  volume_1d = webbur::c1_leg_monomial_integral ( expon );
  volume = std::pow ( volume_1d, n );

  webbur::r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = volume;

  return;
}
//****************************************************************************80

int cn_leg_01_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_01_1_SIZE sizes the midpoint rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_01_1_SIZE, the order.
//
{
  int o;

  o = 1;

  return o;
}
//****************************************************************************80

void cn_leg_02_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_02_XIU implements the Xiu precision 2 rule for region CN_LEG.
//
//  Discussion:
//
//    This is the same as the Stroud rule CN:2-1.
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = webbur::r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = 1.0;
  delta0 = 0.0;
  c1 = 1.0 / 3.0;

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( std::sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = webbur::c1_leg_monomial_integral ( expon );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int cn_leg_02_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_02_XIU_SIZE sizes the Xiu precision 2 rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_02_XIU_SIZE, the order.
//
{
  int o;

  o = n + 1;

  return o;
}
//****************************************************************************80

void cn_leg_03_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_03_1 implements the Stroud rule CN:3-1 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//    The necessary treatment of the final coordinate of points when
//    N is odd seems to vary from what Stroud declares! 
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int expon;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  expon = 0;
  volume_1d = webbur::c1_leg_monomial_integral ( expon );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( ( 2 * r - 1 ) * ( j + 1 ) ) * pi / ( double ) ( n );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg ) / std::sqrt ( 3.0 );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg ) / std::sqrt ( 3.0 );
      i = i + 1;
    }
//
//  The following code does not correspond to what Stroud declares.
//
    if ( i < n )
    {
      if ( n == 1 )
      {
        x[i+j*n] =                     webbur::r8_mop ( j + 1 ) / std::sqrt ( 3.0 );
      }
      else
      {
        x[i+j*n] = std::sqrt ( 2.0 ) * webbur::r8_mop ( j + 1 ) / std::sqrt ( 3.0 );
      }
      i = i + 1;
    }
  }

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int cn_leg_03_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_03_1_SIZE sizes the Stroud rule CN:3-1 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_03_1_SIZE, the order.
//
{
  int o;

  o = 2 * n;

  return o;
}
//****************************************************************************80

void cn_leg_03_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_03_XIU implements the Xiu precision 3 rule for region CN_LEG.
//
//  Discussion:
//
//    This rule is identical to the Stroud rule CN:3-1.
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int expon;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( ( 2 * r - 1 ) * ( j + 1 ) ) * pi / ( double ) ( n );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg ) / std::sqrt ( 3.0 );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg ) / std::sqrt ( 3.0 );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = std::sqrt ( 2.0 ) * webbur::r8_mop ( j + 1 ) / std::sqrt ( 3.0 );
      if ( n == 1 )
      {
        x[i+j*n] = x[i+j*n] / std::sqrt ( 2.0 );
      }
      i = i + 1;
    }
  }

  expon = 0;
  volume_1d = webbur::c1_leg_monomial_integral ( expon );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int cn_leg_03_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_03_XIU_SIZE, the order.
//
{
  int o;

  o = 2 * n;

  return o;
}
//****************************************************************************80

void cn_leg_05_1 ( int n, int option, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_05_1 implements the Stroud rule CN:5-1 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N^2 + N + 2.
//
//    The rule has precision P = 5.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    N must be 4, 5, or 6.
//
//    Input, int OPTION, is only used in case N = 4 or 5.
//    In that case, OPTION should be 1 or 2 to select the
//    two available variants of the rule.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double arg;
  double b;
  double c;
  double eta;
  int expon;
  double gamma;
  int i;
  int i1;
  int i2;
  int k;
  double lambda;
  double mu;
  double volume;
  double volume_1d;
  double xsi;

  if ( n < 4 || 6 < n )
  {
    std::cerr << "\n";
    std::cerr << "CN_LEG_05_1 - Fatal error!\n";
    std::cerr << "  The value of N must be 4, 5, or 6.\n";
    std::exit ( 1 );
  }

  if ( n == 4 || n == 5 )
  {
    if ( option < 1 || 2 < option )
    {
      std::cerr << "\n";
      std::cerr << "CN_LEG_05_1 - Fatal error!\n";
      std::cerr << "  When N = 4 or 5, the value of OPTION must be 1 or 2.\n";
      std::exit ( 1 );
    }
  }

  expon = 0;
  volume_1d = webbur::c1_leg_monomial_integral ( expon );
  volume = std::pow ( volume_1d, n );

  if ( n == 4 && option == 1 )
  {
    eta    =   0.778984505799815E+00;
    lambda =   1.284565137874656E+00;
    xsi =    - 0.713647298819253E+00;
    mu =     - 0.715669761974162E+00;
    gamma =    0.217089151000943E+00;
    a =        0.206186096875899E-01 * volume;
    b =        0.975705820221664E-02 * volume;
    c =        0.733921929172573E-01 * volume;
  }
  else if ( n == 4 && option == 2 )
  {
    eta    =   0.546190755827425E+00;
    lambda =   0.745069130115661E+00;
    xsi =    - 0.413927294508700E+00;
    mu =     - 0.343989637454535E+00;
    gamma =    1.134017894600344E+00;
    a =        0.853094758323323E-01 * volume;
    b =        0.862099000096395E-01 * volume;
    c =        0.116418206881849E-01 * volume;
  }
  else if ( n == 5 && option == 1 )
  {
    eta    =   0.522478547481276E+00;
    lambda =   0.936135175985774E+00;
    xsi =    - 0.246351362101519E+00;
    mu =     - 0.496308106093758E+00;
    gamma =    0.827180176822930E+00;
    a =        0.631976901960153E-01 * volume;
    b =        0.511464127430166E-01 * volume;
    c =        0.181070246088902E-01 * volume;
  }
  else if ( n == 5 && option == 2 )
  {
    eta    =   0.798317301388741E+00;
    lambda =   0.637344273885728E+00;
    xsi =    - 0.455245909918377E+00;
    mu =     - 1.063446229997311E+00;
    gamma =    0.354482076665770E+00;
    a =        0.116952384292206E-01 * volume;
    b =        0.701731258612708E-01 * volume;
    c =        0.137439132264426E-01 * volume;
  }
  else if ( n == 6 )
  {
    eta    =   0.660225291773525E+00;
    lambda =   1.064581294844754E+00;
    xsi =      0.000000000000000E+00;
    mu =     - 0.660225291773525E+00;
    gamma =    0.660225291773525E+00;
    a =        0.182742214532872E-01 * volume;
    b =        0.346020761245675E-01 * volume;
    c =        0.182742214532872E-01 * volume;
  }
  k = -1;

  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = eta;
  }
  w[k] = a;

  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = - eta;
  }
  w[k] = a;

  for ( i1 = 0; i1 < n; i1++ )
  {
    k = k + 1;
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = xsi;
    }
    x[i1+k*n] = lambda;
    w[k] = b;
  }

  for ( i1 = 0; i1 < n; i1++ )
  {
    k = k + 1;
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = - xsi;
    }
    x[i1+k*n] = - lambda;
    w[k] = b;
  }

  for ( i1 = 0; i1 < n - 1; i1++ )
  {
    for ( i2 = i1 + 1; i2 < n; i2++ )
    {
      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = gamma;
      }
      x[i1+k*n] = mu;
      x[i2+k*n] = mu;
      w[k] = c;
    }
  }

  for ( i1 = 0; i1 < n - 1; i1++ )
  {
    for ( i2 = i1 + 1; i2 < n; i2++ )
    {
      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = - gamma;
      }
      x[i1+k*n] = - mu;
      x[i2+k*n] = - mu;
      w[k] = c;
    }
  }

  return;
}
//****************************************************************************80

int cn_leg_05_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_05_1_SIZE sizes the Stroud rule CN:5-1 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N^2 + N + 2.
//
//    The rule has precision P = 5.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_05_1_SIZE, the order.
//
{
  int o;

  o = n * n + n + 2;

  return o;
}
//****************************************************************************80

void cn_leg_05_2 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_05_2 implements the Stroud rule CN:5-2 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 N^2 + 1.
//
//    The rule has precision P = 5.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    N must be at least 2.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double b0;
  double b1;
  double b2;
  int expon;
  int i;
  int i1;
  int i2;
  int k;
  double r;
  double volume;
  double volume_1d;

  if ( n < 2 )
  {
    std::cerr << "\n";
    std::cerr << "CN_LEG_05_2 - Fatal error!\n";
    std::cerr << "  N must be at least 2.\n";
    std::exit ( 1 );
  }

  expon = 0;
  volume_1d = webbur::c1_leg_monomial_integral ( expon );
  volume = std::pow ( volume_1d, n );

  b0 = ( double ) ( 25 * n * n - 115 * n + 162 ) * volume / 162.0;
  b1 = ( double ) ( 70 - 25 * n ) * volume / 162.0;
  b2 = 25.0 * volume / 324.0;

  r = std::sqrt ( 3.0 / 5.0 );

  k = - 1;

  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = 0.0;
  }
  w[k] = b0;

  for ( i1 = 0; i1 < n; i1++ )
  {
    k = k + 1;
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = 0.0;
    }
    x[i1+k*n] = + r;
    w[k] = b1;

    k = k + 1;
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = 0.0;
    }
    x[i1+k*n] = - r;
    w[k] = b1;
  }

  for ( i1 = 0; i1 < n - 1; i1++ )
  {
    for ( i2 = i1 + 1; i2 < n; i2++ )
    {
      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = 0.0;
      }
      x[i1+k*n] = + r;
      x[i2+k*n] = + r;
      w[k] = b2;

      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = 0.0;
      }
      x[i1+k*n] = + r;
      x[i2+k*n] = - r;
      w[k] = b2;

      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = 0.0;
      }
      x[i1+k*n] = - r;
      x[i2+k*n] = + r;
      w[k] = b2;

      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = 0.0;
      }
      x[i1+k*n] = - r;
      x[i2+k*n] = - r;
      w[k] = b2;
    }
  }
  return;
}
//****************************************************************************80

int cn_leg_05_2_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_05_2_SIZE sizes the Stroud rule CN:5-2 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 N^2 + 1.
//
//    The rule has precision P = 5.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_05_2_SIZE, the order.
//
{
  int o;

  o = 2 * n * n + 1;

  return o;
}
//****************************************************************************80

double cn_leg_monomial_integral ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on CN.
//
//  Discussion:
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//    value = integral ( CN ) product ( 1 <= i <= n ) x(I)^expon(i) dx(i)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int EXPON(N), the exponents.
//
//    Output, double CN_LEG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = c1_leg_monomial_integral ( expon[i] );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void en_her_01_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_01_1 implements the Stroud rule 1.1 for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EN_HER is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int k;
  double pi = 3.141592653589793E+00;
  double volume;
  double volume_1d;

  volume_1d = std::sqrt ( pi );
  volume = std::pow ( volume_1d, n );

  webbur::r8vec_zero ( n * o, x );

  k = 0;
//
//  1 point.
//
  w[k] = volume;

  return;
}
//****************************************************************************80

int en_her_01_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_01_1_SIZE sizes the Stroud rule 1.1 for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = 1.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_HER_01_1_SIZE, the order.
//
{
  int o;

  o = 1;

  return o;
}
//****************************************************************************80

void en_her_02_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_02_XIU implements the Xiu precision 2 rule for region EN_HER.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EN_HER is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double delta0;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = webbur::r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = 2.0;
  delta0 = 0.0;
  c1 = 1.0;

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( std::sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  volume_1d = std::sqrt ( pi );
  volume = pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int en_her_02_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_02_XIU_SIZE sizes the Xiu precision 2 rule for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = N + 1;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_HER_01_1_SIZE, the order.
//
{
  int o;

  o = n + 1;

  return o;
}
//****************************************************************************80

void en_her_03_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_03_1 implements the Stroud rule 3.1 for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = 2 * N.
//
//    The rule has precision P = 3.
//
//    EN_HER is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  int i;
  int k;
  double pi = 3.141592653589793;
  double r;
  double volume;
  double volume_1d;

  volume_1d = std::sqrt ( pi );
  volume = std::pow ( volume_1d, n );

  a = volume / ( double ) ( o );
  r = std::sqrt ( ( double ) ( n ) / 2.0 );

  webbur::r8vec_zero ( n * o, x );

  k = - 1;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - r;
    w[k] = a;
    k = k + 1;
    x[i+k*n] = + r;
    w[k] = a;
  }

  return;
}
//****************************************************************************80

int en_her_03_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_03_1_SIZE sizes the Stroud rule 3.1 for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = 2 * N.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_HER_03_1_SIZE, the order.
//
{
  int o;

  o = 2 * n;

  return o;
}
//****************************************************************************80

void en_her_03_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_03_XIU implements the Xiu precision 3 rule for region EN_HER.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    EN_HER is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
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
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( ( 2 * r - 1 ) * ( j + 1 ) ) * pi / ( double ) ( n );
      x[i+j*n] = std::cos ( arg );
      i = i + 1;
      x[i+j*n] = std::sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = webbur::r8_mop ( j + 1 );
      if ( n == 1 )
      {
        x[i+j*n] = x[i+j*n] / std::sqrt ( 2.0 );
      }
      i = i + 1;
    }
  }

  volume_1d = std::sqrt ( pi );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int en_her_03_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_03_XIU_SIZE sizes the Xiu precision 3 rule for region EN_HER.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    EN_HER is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
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
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_HER_XIU_SIZE, the order.
//
{
  int o;

  o = 2 * n;

  return o;
}
//****************************************************************************80

void en_her_05_1 ( int n, int option, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_05_1 implements the Stroud rule 5.1 for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = N^2 + N + 2.
//
//    The rule has precision P = 5.
//
//    EN_HER is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
//    the OPTION variable to 1 or 2.
//
//    Versions of this rule are only available for N = 2 through 7.
//
//    There is a typographical error in the reference.
//    For the second version of the rule for N = 2, the line
//      gamma =    0.313300683022281E+00
//    should read
//      gamma =    0.312200683022281E+00
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    2 <= N <= 7.
//
//    Input, int OPTION, selects option 1 or 2.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double b;
  double c;
  double eta;
  double gamma;
  int i;
  int i1;
  int j;
  int k;
  double lambda;
  double mu;
  double pi = 3.141592653589793E+00;
  double volume;
  double volume_1d;
  double xsi;

  if ( n < 2 || 7 < n )
  {
    std::cerr << "\n";
    std::cerr << "EN_HER_05_1 - Fatal error!\n";
    std::cerr << "  2 <= N <= 7 required.\n";
    std::exit ( 1 );
  }

  if ( option < 1 || 2 < option )
  {
    std::cerr << "\n";
    std::cerr << "EN_HER_05_1 - Fatal error!\n";
    std::cerr << "  1 <= OPTION <= 2 required.\n";
    std::exit ( 1 );
  }

  if ( option == 2 )
  {
    if ( n != 3 && n != 5 && n != 6 )
    {
      std::cerr << "\n";
      std::cerr << "EN_HER_05_1 - Fatal error!\n";
      std::cerr << "  OPTION = 2 requires N = 3, 5 or 6.\n";
      std::exit ( 1 );
    }
  }

  volume_1d = std::sqrt ( pi );
  volume = std::pow ( volume_1d, n );

  if ( n == 2 )
  {
    eta =      0.446103183094540E+00;
    lambda =   0.136602540378444E+01;
    xsi =    - 0.366025403784439E+00;
    mu =       0.198167882945871E+01;
    gamma =    0.000000000000000E+00;
    a =        0.328774019778636E+00 * volume;
    b =        0.833333333333333E-01 * volume;
    c =        0.455931355469736E-02 * volume;
  }
  else if ( n == 3 && option == 1 )
  {
    eta =      0.476731294622796E+00;
    lambda =   0.935429018879534E+00;
    xsi =    - 0.731237647787132E+00;
    mu =       0.433155309477649E+00;
    gamma =    0.266922328697744E+01;
    a =        0.242000000000000E+00 * volume;
    b =        0.810000000000000E-01 * volume;
    c =        0.500000000000000E-02 * volume;
  }
//
//  The value of gamma that follows corrects an error in the reference.
//
  else if ( n == 3 && option == 2 )
  {
    eta =      0.476731294622796E+00;
    lambda =   0.128679320334269E+01;
    xsi =    - 0.379873463323979E+00;
    mu =     - 0.192386729447751E+01;
    gamma =    0.312200683022281E+00;
    a =        0.242000000000000E+00 * volume;
    b =        0.810000000000000E-01 * volume;
    c =        0.500000000000000E-02 * volume;
  }
  else if ( n == 4 ) 
  {
    eta =      0.523945658287507E+00;
    lambda =   0.119433782552719E+01;
    xsi =    - 0.398112608509063E+00;
    mu =     - 0.318569372920112E+00;
    gamma =    0.185675837424096E+01;
    a =        0.155502116982037E+00 * volume;
    b =        0.777510584910183E-01 * volume;
    c =        0.558227484231506E-02 * volume;
  }
  else if ( n == 5 && option == 1 )
  {
    eta =      0.214972564378798E+01;
    lambda =   0.464252986016289E+01;
    xsi =    - 0.623201054093728E+00;
    mu =     - 0.447108700673434E+00;
    gamma =    0.812171426076311E+00;
    a =        0.487749259189752E-03 * volume;
    b =        0.487749259189752E-03 * volume;
    c =        0.497073504444862E-01 * volume;
  }
  else if ( n == 5 && option == 2 )
  {
    eta =      0.615369528365158E+00;
    lambda =   0.132894698387445E+01;
    xsi =    - 0.178394363877324E+00;
    mu =     - 0.745963266507289E+00;
    gamma =    0.135503972310817E+01;
    a =        0.726415024414905E-01 * volume;
    b =        0.726415024414905E-01 * volume;
    c =        0.641509853510569E-02 * volume;
  }
  else if ( n == 6 && option == 1 )
  {
    eta =      0.100000000000000E+01;
    lambda =   0.141421356237309E+01;
    xsi =      0.000000000000000E+00;
    mu =     - 0.100000000000000E+01;
    gamma =    0.100000000000000E+01;
    a =        0.781250000000000E-02 * volume;
    b =        0.625000000000000E-01 * volume;
    c =        0.781250000000000E-02 * volume;
  }
  else if ( n == 6 && option == 2 )
  {
    eta =      0.100000000000000E+01;
    lambda =   0.942809041582063E+00;
    xsi =    - 0.471404520791032E+00;
    mu =     - 0.166666666666667E+01;
    gamma =    0.333333333333333E+00;
    a =        0.781250000000000E-02 * volume;
    b =        0.625000000000000E-01 * volume;
    c =        0.781250000000000E-02 * volume;
  }
  else if ( n == 7 )
  {
    eta =      0.000000000000000E+00;
    lambda =   0.959724318748357E+00;
    xsi =    - 0.772326488820521E+00;
    mu =     - 0.141214270131942E+01;
    gamma =    0.319908106249452E+00;
    a =        0.111111111111111E+00 * volume;
    b =        0.138888888888889E-01 * volume;
    c =        0.138888888888889E-01 * volume;
  }

  webbur::r8vec_zero ( n * o, x );

  k = - 1;
//
//  2 points.
//
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = - eta;
  }
  w[k] = a;
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = + eta;
  }
  w[k] = a;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    for ( i1 = 0; i1 < n; i1++ )
    {
      x[i1+k*n] = - xsi;
    }
    x[i+k*n] = - lambda;
    w[k] = b;
    k = k + 1;
    for ( i1 = 0; i1 < n; i1++ )
    {
      x[i1+k*n] = + xsi;
    }
    x[i+k*n] = + lambda;
    w[k] = b;
  }
//
//  2 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      for ( i1 = 0; i1 < n; i1++ )
      {
        x[i1+k*n] = - gamma;
      }
      x[i+k*n] = - mu;
      x[j+k*n] = - mu;
      w[k] = c;
      k = k + 1;
      for ( i1 = 0; i1 < n; i1++ )
      {
        x[i1+k*n] = + gamma;
      }
      x[i+k*n] = + mu;
      x[j+k*n] = + mu;
      w[k] = c;
    }
  }
  return;
}
//****************************************************************************80

int en_her_05_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_05_1_SIZE sizes the Stroud rule 5.1 for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = N^2 + N + 2.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_HER_05_1_SIZE, the order.
//
{
  int o;

  o = n * n + n + 2;

  return o;
}
//****************************************************************************80

void en_her_05_2 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_05_2 implements the Stroud rule 5.2 for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = 2 * N^2 + 1.
//
//    The rule has precision P = 5.
//
//    EN_HER is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double b;
  double c;
  int i;
  int j;
  int k;
  double pi = 3.141592653589793;
  double r;
  double s;
  double volume;
  double volume_1d;

  volume_1d = std::sqrt ( pi );
  volume = std::pow ( volume_1d, n );

  a = 2.0E+00 * volume / ( double ) ( n + 2 );
  b = ( double ) ( 4 - n ) * volume / 2.0
    / ( double ) ( ( n + 2 ) * ( n + 2 ) );
  c = volume / ( double ) ( ( n + 2 ) * ( n + 2 ) );

  r = std::sqrt ( ( double ) ( n + 2 ) / 2.0 );
  s = std::sqrt ( ( double ) ( n + 2 ) / 4.0 );

  webbur::r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = a;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - r;
    w[k] = b;
    k = k + 1;
    x[i+k*n] = + r;
    w[k] = b;
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - s;
      x[j+k*n] = - s;
      w[k] = c;
      k = k + 1;
      x[i+k*n] = - s;
      x[j+k*n] = + s;
      w[k] = c;
      k = k + 1;
      x[i+k*n] = + s;
      x[j+k*n] = - s;
      w[k] = c;
      k = k + 1;
      x[i+k*n] = + s;
      x[j+k*n] = + s;
      w[k] = c;
    }
  }
  return;
}
//****************************************************************************80

int en_her_05_2_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_05_2_SIZE sizes the Stroud rule 5.2 for region EN_HER.
//
//  Discussion:
//
//    The rule has order O = 2 * N^2 + 1.
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_HER_01_1_SIZE, the order.
//
{
  int o;

  o = 2 * n * n + 1;

  return o;
}
//****************************************************************************80

double en_her_monomial_integral ( int n, int alpha[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_HER_MONOMIAL_INTEGRAL evaluates monomial integrals in EN_HER.
//
//  Discussion:
//
//    ALPHA is the set of polynomial exponents.
//
//    EN_HER is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
//
//    The integral to be evaluated is
//
//      value = integral ( EN ) x(1)^alpha(1) * x(2)^alpha(2) * ... 
//        * x(n)^alpha(n) * w(x) dx
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
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int ALPHA[N], the polynomial exponents.
//    0 <= ALPHA[*].
//
//    Output, double EN_HER_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double arg;
  int i;
  double value;

  for ( i = 0; i < n; i++ )
  {
    if ( alpha[i] < 0 )
    {
      std::cerr << "\n";
      std::cerr << "EN_HER_MONOMIAL_INTEGRAL - Fatal error//\n";
      std::cerr << "  ALPHA[" << i << "] < 0.\n";
      std::exit ( 1 );
    }
  }

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    if ( ( alpha[i] % 2 == 1 ) )
    {
      value = 0.0;
      break;
    }
    else
    {
      arg = ( ( double ) ( alpha[i] + 1 ) ) / 2.0;
      value = value * webbur::r8_gamma ( arg );
    }
  }

  return value;
}
//****************************************************************************80

double ep1_glg_monomial_integral ( int expon, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EP1_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EP1.
//
//  Discussion:
//
//    EP1_GLG is the interval [0,+oo) with generalized Laguerre weight function:
//
//      w(alpha;x) = x^alpha exp ( - x )
//
//    value = integral ( 0 <= x < +oo ) x^expon x^alpha exp ( - x ) dx
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
//  Parameters:
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, double EP1_GLG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double arg;
  double exact;

  arg = alpha + ( double ) ( expon + 1 );

  exact = webbur::r8_gamma ( arg );

  return exact;
}
//****************************************************************************80

double ep1_lag_monomial_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    EP1_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EP1.
//
//  Discussion:
//
//    EP1 is the interval [0,+oo) with exponential or Laguerre weight function:
//
//      w(x) = exp ( - x )
//
//    value = integral ( 0 <= x < oo ) x^expon exp ( - x ) dx
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
//  Parameters:
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double EP1_LAG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double value;

  value = webbur::r8_factorial ( expon );

  return value;
}
//****************************************************************************80

void epn_glg_01_1 ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_01_1 implements a precision 1 rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int i;
  int k;
  double value2;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "EPN_GLG_01_1 - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  expon = 0;
  volume_1d = webbur::ep1_glg_monomial_integral ( expon, alpha );
  volume = std::pow ( volume_1d, n );

  expon = 1;
  value2 = webbur::ep1_glg_monomial_integral ( expon, alpha );

  webbur::r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = value2 / volume_1d;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int epn_glg_01_1_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_01_1_SIZE sizes a precision 1 rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, int EPN_GLG_01_1_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "EPN_GLG_01_1_SIZE - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void epn_glg_02_xiu ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_02_XIU implements the Xiu precision 2 rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "EPN_GLG_02_XIU - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = webbur::r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = - 1.0;
  delta0 = alpha + 1.0;
  c1 = - alpha - 1.0;

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( std::sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = webbur::ep1_glg_monomial_integral ( expon, alpha );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int epn_glg_02_xiu_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_02_XIU_SIZE sizes the Xiu precision 2 rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
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
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, int EPN_GLG_02_XIU_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "EPN_GLG_02_XIUI_SIZE - Fatal error!\n";
    std::cerr << "  ALPHA <= -1.0\n";
    std::exit ( 1 );
  }

  o = n + 1;

  return o;
}
//****************************************************************************80

double epn_glg_monomial_integral ( int n, int expon[], double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EPN.
//
//  Discussion:
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//    value = integral ( EPN ) 
//      product ( 1 <= i <= n ) x(I)^expon(i) x(i)^alpha exp ( - x(i) ) dx(i)
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int EXPON[N], the exponents.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, double EPN_GLG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = webbur::ep1_glg_monomial_integral ( expon[i], alpha );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void epn_lag_01_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_01_1 implements a precision 1 rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int i;
  int k;
  double value2;
  double volume;
  double volume_1d;

  expon = 0;
  volume_1d = webbur::ep1_lag_monomial_integral ( expon );
  volume = std::pow ( volume_1d, n );

  expon = 1;
  value2 = webbur::ep1_lag_monomial_integral ( expon );

  webbur::r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = value2 / volume_1d;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int epn_lag_01_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_01_1_SIZE sizes a precision 1 rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EPN_LOG_01_1_SIZE, the order.
//
{
  int o;

  o = 1;

  return o;
}
//****************************************************************************80

void epn_lag_02_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_02_XIU implements the Xiu precision 2 rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = webbur::r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = - 1.0;
  delta0 = 1.0;
  c1 = - 1.0;

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( std::sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = webbur::ep1_lag_monomial_integral ( expon );
  volume = std::pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int epn_lag_02_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_02_XIU_SIZE sizes the Xiu precision 2 rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
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
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EPN_LAG_02_XIU_SIZE, the order.
//
{
  int o;

  o = n + 1;

  return o;
}
//****************************************************************************80

double epn_lag_monomial_integral ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EPN.
//
//  Discussion:
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//    value = integral ( EPN ) 
//      product ( 1 <= i <= n ) x(I)^expon(i) exp ( -x(i) ) dx(i)
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
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int EXPON(N), the exponents.
//
//    Output, double EPN_LAG_MONOMIAL_VALUE, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = webbur::ep1_lag_monomial_integral ( expon[i] );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void gw_02_xiu ( int n, int o, double gamma0, double delta0, double c1, 
  double volume_1d, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    GW_02_XIU implements the Golub-Welsch version of the Xiu rule.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    It is assumed that the integral is over an N-dimensional region,
//    and has the form
//
//      Integral f(x) w(x) dx
//
//    where w(x) is separable into identical and independent components:
//
//      w(x) = v(x1) * v(x2) * ... * v(xn)
//
//    Associated with the weight function v(x), we assume there is a
//    family of orthogonal polynomials satisfying a three-term recurrence
//    of the form:
//
//      x P(n,x) = An * P(n+1,x) + Bn * P(n,x) + Cn * P(n-1,x)
//
//    with P(0,x) = 1, and P(-1,x) = 0.
//
//    This routine can construct the desired quadrature rule by knowing
//    the values of C1, used in the definition of P2, the values
//    GAMMA0 = 1/A0 and DELTA0 = - B0/A0, for which it is the case that
//    P(1,X) = GAMMA0 * X + DELTA0, and the value of VOLUME_1D, that is,
//    the 1D integral of v(x) over the region.
//
//    Note the values for the following standard polynomial families:
//
//    Chebyshev Type 1
//      V(X) =      1 / sqrt ( 1 - X^2 )
//      Interval =  [-1,+1]
//      GAMMA0 =    1.0
//      DELTA0 =    0.0
//      C1 =        1/2
//      VOLUME_1D = PI
//
//    Chebyshev Type 2
//      V(X) =      sqrt ( 1 - X^2 )
//      Interval =  [-1,+1]
//      GAMMA0 =    2.0
//      DELTA0 =    0.0
//      C1 =        1/2
//      VOLUME_1D = PI / 2
//
//    Gegenbauer
//      V(X) =      ( 1 - X^2 )^A
//      Interval =  [-1,+1]
//      GAMMA0 =    2 * A + 1
//      DELTA0 =    0.0
//      C1 =        ( 2 * A + 1 ) / ( 2 A + 3 )
//      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
//
//    Gegenbauer* (Removes singularity at ALPHA = -0.5):
//      V(X) =      ( 1 - X^2 )^A
//      Interval =  [-1,+1]
//      GAMMA0 =    1
//      DELTA0 =    0.0
//      C1 =        1 / ( 2 A + 3 )
//      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
//
//    Generalized Hermite
//      V(X) = |x|^A exp ( - x^2 )
//      Interval = (-oo,+oo)
//      GAMMA0 =    2
//      DELTA0 =    0
//      C1 =        2+2A
//      VOLUME_1D = Gamma((A+1)/2)
//
//    Generalized Laguerre
//      V(X) =       x^A exp ( - x )
//      Interval =  [0,+oo)
//      GAMMA0 =    -1.0
//      DELTA0 =     A+1.0
//      C1 =        -A-1.0
//      VOLUME_1D =  Gamma(A+1)
//
//    Hermite (physicist)
//      V(X) =      exp ( - x^2 )
//      Interval =  (-oo,+oo)
//      GAMMA0 =    2.0
//      DELTA0 =    0.0
//      C1 =        1.0
//      VOLUME_1D = sqrt ( PI )
//
//    Hermite (probabilist)
//      V(X) =      exp ( - x^2 / 2 )
//      Interval =  (-oo,+oo)
//      GAMMA0 =    1.0
//      DELTA0 =    0.0
//      C1 =        1.0
//      VOLUME_1D = sqrt ( 2 PI )
//
//    Jacobi
//      V(X) =      (1-x)^A (1+x)^B
//      Interval =  [-1,+1]
//      GAMMA0 =    (A+B+2)/2  
//      DELTA0 =    (A-B)/2
//      C1 =        2(A+1)(B+1)/(A+B+3)/(A+B+2)
//      VOLUME_1D = 2^(A+B+1) * Gamma(A+1) * Gamma(B+1) / ( A+B+1) / Gamma(A+B+1)
//
//    Laguerre
//      V(X) =       exp ( - x )
//      Interval =  [0,+oo)
//      GAMMA0 =    -1.0
//      DELTA0 =     1.0
//      C1 =        -1.0
//      VOLUME_1D =  1.0
//
//    Legendre
//      V(X) =      1.0
//      Interval =  [-1,+1]
//      GAMMA0 =    1.0
//      DELTA0 =    0.0
//      C1 =        1/3
//      VOLUME_1D = 2.0
//                                  
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Input, double GAMMA0, the ratio 1 / A0.
//
//    Input, double DELTA0, the ratio B0 / A0.
//
//    Input, double C1, the coefficient of P(0,X) in the definition of P(2,X).
//
//    Input, double VOLUME_1D, the 1D integral of V(X).
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );
      x[i+j*n] = std::sqrt ( 2.0 ) * std::cos ( arg );
      i = i + 1;
      x[i+j*n] = std::sqrt ( 2.0 ) * std::sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = webbur::r8_mop ( j );
      i = i + 1;
    }
  }
//
//  Adjust for the GW rule.
//
  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( std::sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }
//
//  The weights are equal.
//
  for ( j = 0; j < o; j++ )
  {
    w[j] = std::pow ( volume_1d, n ) / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int gw_02_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    GW_02_XIU_SIZE sizes the Golub Welsch version of the Xiu rule.
//
//  Discussion:
//
//    The rule has order O = N + 1;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int GW_02_XIU_SIZE, the order.
//
{
  int o;

  o = n + 1;

  return o;
}
//****************************************************************************80

double *monomial_value ( int dim_num, int point_num, double x[], int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points at which the
//    monomial is to be evaluated.
//
//    Input, double X[DIM_NUM*POINT_NUM], the point coordinates.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Output, double MONOMIAL_VALUE[POINT_NUM], the value of the monomial.
//
{
  int dim;
  int point;
  double *value;

  value = new double[point_num];

  for ( point = 0; point < point_num; point++ )
  {
    value[point] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0 != expon[dim] )
    {
      for ( point = 0; point < point_num; point++ )
      {
        value[point] = value[point] * std::pow ( x[dim+point*dim_num], expon[dim] );
      }
    }
  }

  return value;
}
}
