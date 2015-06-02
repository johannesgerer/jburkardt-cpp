# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "quadmom.hpp"
# include "toms655.hpp"

int main ( );
void quadmom_prb01 ( );
void quadmom_prb02 ( );
void quadmom_prb03 ( );
void quadmom_prb04 ( );
void quadmom_prb05 ( );
void quadmom_prb06 ( );
void quadmom_prb07 ( );
void quadmom_prb08 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QUADMOM_PRB.
//
//  Discussion:
//
//    QUADMOM_PRB tests the QUADMOM library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "QUADMOM_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the QUADMOM library.\n";

  quadmom_prb01 ( );
  quadmom_prb02 ( );
  quadmom_prb03 ( );
  quadmom_prb04 ( );
  quadmom_prb05 ( );
  quadmom_prb06 ( );
  quadmom_prb07 ( );
  quadmom_prb08 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QUADMOM_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

 return 0;
}
//****************************************************************************80

void quadmom_prb01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUADMOM_PRB01 tests the QUADMOM procedure for the Legendre weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, John Welsch,
//    Calculation of Gaussian Quadrature Rules,
//    Mathematics of Computation,
//    Volume 23, Number 106, April 1969, pages 221-230.
//
{
  double a;
  double alpha2;
  double b;
  double beta2;
  int kind;
  int lo;
  int m;
  double *moment;
  int n;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "QUADMOM_PRB01:\n";
  cout << "  Compute the points and weights of a quadrature rule\n";
  cout << "  for the Legendre weight, rho(x)=1, over [-1,+1],\n";
  cout << "  using Golub and Welsch's moment method.\n";
  cout << "  Compare with a standard calculation.\n";
//
//  N is the order of the rule we want to compute.
//
  n = 5;
//
//  Compute M = 2*N+1 moments for the Legendre weight on [-1,+1].
//
  m = 2 * n + 1;
  a = -1.0;
  b = 1.0;

  moment = moments_legendre ( m, a, b );
//
//  Compute the points and weights by the method of moments.
//
  x1 = new double[n];
  w1 = new double[n];

  moment_method ( n, moment, x1, w1 );
//
//  Compute the points and weights the usual way.
//
  kind = 1;
  alpha2 = 0.0;
  beta2 = 0.0;
  a = -1.0;
  b = +1.0;
  lo = 0;
  x2 = new double[n];
  w2 = new double[n];

  cgqf ( n, kind, alpha2, beta2, a, b, lo, x2, w2 );
//
//  Compare the results.
//
  r8vec2_print ( n, x1, x2, 
    "  Points from GW moment and orthogonal polynomial methods:" );

  r8vec2_print ( n, w1, w2, 
    "  Weights from GW moment and orthogonal polynomial methods:" );
//
//  Free memory.
//
  delete [] moment;
  delete [] w1;
  delete [] w2;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void quadmom_prb02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUADMOM_PRB02 tests the QUADMOM procedure for the standard Gaussian weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, John Welsch,
//    Calculation of Gaussian Quadrature Rules,
//    Mathematics of Computation,
//    Volume 23, Number 106, April 1969, pages 221-230.
//
{
  double a;
  double alpha2;
  double b;
  double beta2;
  int i;
  int kind;
  int lo;
  int m;
  double *moment;
  int n;
  const double pi = 3.141592653589793;
  double sigma;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "QUADMOM_PRB02:\n";
  cout << "  Compute the points and weights of a quadrature rule for\n";
  cout << "  the standard Gaussian weight, rho(x)=exp(-x^2/2)/sqrt(2pi),\n";
  cout << "  over (-oo,+oo), using Golub and Welsch's moment method.\n";
  cout << "  Compare with a standard calculation.\n";
//
//  N is the order of the quadrature rule.
//
  n = 5;
//
//  Compute the M = 2 * N + 1 moments for the standard Gaussian weight on (-oo,+oo).
//
  m = 2 * n + 1;
  moment = moments_normal_01 ( m );
//
//  Compute the points and weights by the method of moments.
//
  x1 = new double[n];
  w1 = new double[n];

  moment_method ( n, moment, x1, w1 );
//
//  Compute the points and weights the usual way.
//
  kind = 6;
  alpha2 = 0.0;
  beta2 = 0.0;
  a = 0.0;
  b = +0.5;
  lo = 0;
  x2 = new double[n];
  w2 = new double[n];

  cgqf ( n, kind, alpha2, beta2, a, b, lo, x2, w2 );
//
//  The CGQF weights need to be normalized by sigma * sqrt ( 2 * pi )
//  because they don't divide the Gaussian PDF by that factor.
//
  sigma = 1.0;
  for ( i = 0; i < n; i++ )
  {
    w2[i] = w2[i] / sigma / sqrt ( 2.0 * pi );
  }
//
//  Compare the results.
//
  r8vec2_print ( n, x1, x2, 
    "  Points from GW moment and orthogonal polynomial methods:" );

  r8vec2_print ( n, w1, w2, 
    "  Weights from GW moment and orthogonal polynomial methods:" );
//
//  Free memory.
//
  delete [] moment;
  delete [] w1;
  delete [] w2;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void quadmom_prb03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUADMOM_PRB03 tests the QUADMOM procedure for the general Gaussian weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, John Welsch,
//    Calculation of Gaussian Quadrature Rules,
//    Mathematics of Computation,
//    Volume 23, Number 106, April 1969, pages 221-230.
//
{
  double a;
  double alpha2;
  double b;
  double beta2;
  int i;
  int kind;
  int lo;
  int m;
  double *moment;
  double mu;
  int n;
  const double pi = 3.141592653589793;
  double sigma;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "QUADMOM_PRB03:\n";
  cout << "  Compute the points and weights of a quadrature rule for\n";
  cout << "  a general Gaussian weight,\n";
  cout << "  rho(mu,s;x)=exp(-((x-mu)/sigma)^2/2)/sigma^2/sqrt(2pi),\n";
  cout << "  over (-oo,+oo), using Golub and Welsch''s moment method.\n";
  cout << "  Compare with a standard calculation.\n";
//
//  N is the order of the quadrature rule.
//
  n = 5;
//
//  Compute the M = 2 * N + 1 moments for a general Gaussian weight on (-oo,+oo).
//
  m = 2 * n + 1;
  mu = 1.0;
  sigma = 2.0;

  cout << "\n";
  cout << "  MU = " << mu << "\n";
  cout << "  SIGMA = " << sigma << "\n";

  moment = moments_normal ( m, mu, sigma );
//
//  Compute the points and weights by the method of moments.
//
  x1 = new double[n];
  w1 = new double[n];

  moment_method ( n, moment, x1, w1 );
//
//  Compute the points and weights the usual way.
//
  kind = 6;
  alpha2 = 0.0;
  beta2 = 0.0;
  a = 1.0;
  b = 0.5 / sigma / sigma;
  lo = 0;

  x2 = new double[n];
  w2 = new double[n];

  cgqf ( n, kind, alpha2, beta2, a, b, lo, x2, w2 );
//
//  The CGQF weights need to be normalized by sigma * sqrt ( 2 * pi )
//  because they don't divide the Gaussian PDF by that factor.
//
  for ( i = 0; i < n; i++ )
  {
    w2[i] = w2[i] / sigma / sqrt ( 2.0 * pi );
  }
//
//  Compare the results.
//
  r8vec2_print ( n, x1, x2, 
    "  Points from GW moment and orthogonal polynomial methods:" );

  r8vec2_print ( n, w1, w2, 
    "  Weights from GW moment and orthogonal polynomial methods:" );
//
//  Free memory.
//
  delete [] moment;
  delete [] w1;
  delete [] w2;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void quadmom_prb04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUADMOM_PRB04 tests the QUADMOM procedure for the Laguerre weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, John Welsch,
//    Calculation of Gaussian Quadrature Rules,
//    Mathematics of Computation,
//    Volume 23, Number 106, April 1969, pages 221-230.
//
{
  double a;
  double alpha2;
  double b;
  double beta2;
  int kind;
  int lo;
  int m;
  double *moment;
  int n;
  const double pi = 3.141592653589793;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "QUADMOM_PRB04:\n";
  cout << "  Compute the points and weights of a quadrature rule for\n";
  cout << "  the Laguerre weight, rho(x)=exp(-x),\n";
  cout << "  over [0,+oo), using Golub and Welsch's moment method.\n";
  cout << "  Compare with a standard calculation.\n";
//
//  N is the order of the quadrature rule.
//
  n = 5;
//
//  Compute the M = 2 * N + 1 moments for the Laguerre weight on [0,+oo).
//
  m = 2 * n + 1;
  moment = moments_laguerre ( m );
//
//  Compute the points and weights by the method of moments.
//
  x1 = new double[n];
  w1 = new double[n];

  moment_method ( n, moment, x1, w1 );
//
//  Compute the points and weights the usual way.
//
  kind = 5;
  alpha2 = 0.0;
  beta2 = 0.0;
  a = 0.0;
  b = +1.0;
  lo = 0;

  x2 = new double[n];
  w2 = new double[n];

  cgqf ( n, kind, alpha2, beta2, a, b, lo, x2, w2 );
//
//  Compare the results.
//
  r8vec2_print ( n, x1, x2, 
    "  Points from GW moment and orthogonal polynomial methods:" );

  r8vec2_print ( n, w1, w2, 
    "  Weights from GW moment and orthogonal polynomial methods:" );
//
//  Free memory.
//
  delete [] moment;
  delete [] w1;
  delete [] w2;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void quadmom_prb05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUADMOM_PRB05 tests the QUADMOM procedure for the truncated normal weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, John Welsch,
//    Calculation of Gaussian Quadrature Rules,
//    Mathematics of Computation,
//    Volume 23, Number 106, April 1969, pages 221-230.
//
{
  double a;
  double b;
  int m;
  double *moment;
  double mu;
  int n;
  int order;
  double sigma;
  double *w1;
  double *x1;

  cout << "\n";
  cout << "QUADMOM_PRB05:\n";
  cout << "  Compute the points and weights of a quadrature rule for\n";
  cout << "  a truncated normal weight,\n";
  cout << "  rho(mu,s;x)=exp(-((x-mu)/sigma)^2/2)/sigma^2/sqrt(2pi),\n";
  cout << "  over [a,b], using Golub and Welsch's moment method.\n";
//
//  N is the order of the quadrature rule.
//
  n = 5;
//
//  Compute the M = 2 * N + 1 moments.
//
  m = 2 * n + 1;
  mu = 100.0;
  sigma = 25.0;
  a = 50.0;
  b = 150.0;

  cout << "\n";
  cout << "  MU = " << mu << "\n";
  cout << "  SIGMA = " << sigma << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";

  moment = moments_truncated_normal_ab ( m, mu, sigma, a, b );
//
//  Compute the points and weights by the method of moments.
//
  x1 = new double[n];
  w1 = new double[n];

  moment_method ( n, moment, x1, w1 );
//
//  Print the results.
//
  r8vec_print ( n, x1, "  Points from GW moment method:" );

  r8vec_print ( n, w1, "  Weights from GW moment method:" );
//
//  Free memory.
//
  delete [] moment;
  delete [] w1;
  delete [] x1;

  return;
}
//****************************************************************************80

void quadmom_prb06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUADMOM_PRB06 tests QUADMOM for the lower truncated normal weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, John Welsch,
//    Calculation of Gaussian Quadrature Rules,
//    Mathematics of Computation,
//    Volume 23, Number 106, April 1969, pages 221-230.
//
{
  double a;
  int m;
  double *moment;
  double mu;
  int n;
  int order;
  double sigma;
  double *w1;
  double *x1;

  cout << "\n";
  cout << "QUADMOM_PRB06:\n";
  cout << "  Compute the points and weights of a quadrature rule for\n";
  cout << "  a truncated normal weight,\n";
  cout << "  rho(mu,s;x)=exp(-((x-mu)/sigma)^2/2)/sigma^2/sqrt(2pi),\n";
  cout << "  over [a,+oo), using Golub and Welsch's moment method.\n";
//
//  N is the order of the quadrature rule.
//
  n = 9;
//
//  Compute the M = 2 * N + 1 moments.
//
  m = 2 * n + 1;
  mu = 2.0;
  sigma = 0.5;
  a = 0.0;

  cout << "\n";
  cout << "  MU = " << mu << "\n";
  cout << "  SIGMA = " << sigma << "\n";
  cout << "  A = " << a << "\n";

  moment = moments_truncated_normal_a ( m, mu, sigma, a );
//
//  Compute the points and weights by the method of moments.
//
  x1 = new double[n];
  w1 = new double[n];

  moment_method ( n, moment, x1, w1 );
//
//  Print the results.
//
  r8vec_print ( n, x1, "  Points from GW moment method:" );

  r8vec_print ( n, w1, "  Weights from GW moment method:" );
//
//  Free memory.
//
  delete [] moment;
  delete [] w1;
  delete [] x1;

  return;
}
//****************************************************************************80

void quadmom_prb07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUADMOM_PRB07 tests QUADMOM for the upper truncated normal weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, John Welsch,
//    Calculation of Gaussian Quadrature Rules,
//    Mathematics of Computation,
//    Volume 23, Number 106, April 1969, pages 221-230.
//
{
  double b;
  int m;
  double *moment;
  double mu;
  int n;
  int order;
  double sigma;
  double *w1;
  double *x1;

  cout << "\n";
  cout << "QUADMOM_PRB07:\n";
  cout << "  Compute the points and weights of a quadrature rule for\n";
  cout << "  a truncated normal weight,\n";
  cout << "  rho(mu,s;x)=exp(-((x-mu)/sigma)^2/2)/sigma^2/sqrt(2pi),\n";
  cout << "  over (-oo,b], using Golub and Welsch's moment method.\n";
//
//  N is the order of the quadrature rule.
//
  n = 9;
//
//  Compute the M = 2 * N + 1 moments.
//
  m = 2 * n + 1;
  mu = 2.0;
  sigma = 0.5;
  b = 3.0;

  cout << "\n";
  cout << "  MU = " << mu << "\n";
  cout << "  SIGMA = " << sigma << "\n";
  cout << "  B = " << b << "\n";

  moment = moments_truncated_normal_b ( m, mu, sigma, b );
//
//  Compute the points and weights by the method of moments.
//
  x1 = new double[n];
  w1 = new double[n];

  moment_method ( n, moment, x1, w1 );
//
//  Print the results.
//
  r8vec_print ( n, x1, "  Points from GW moment method:" );

  r8vec_print ( n, w1, "  Weights from GW moment method:" );
//
//  Free memory.
//
  delete [] moment;
  delete [] w1;
  delete [] x1;

  return;
}
//****************************************************************************80

void quadmom_prb08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUADMOM_PRB08 integrates sin(x) against a lower truncated normal weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gene Golub, John Welsch,
//    Calculation of Gaussian Quadrature Rules,
//    Mathematics of Computation,
//    Volume 23, Number 106, April 1969, pages 221-230.
//
{
  double a;
  int i;
  int m;
  double *moment;
  double mu;
  int n;
  int order;
  double q;
  double sigma;
  double *w1;
  double *x1;

  cout << "\n";
  cout << "QUADMOM_PRB08:\n";
  cout << "  Integrate sin(x) against a lower truncated normal weight.\n";

  mu = 0.0;
  sigma = 1.0;
  a = -3.0;

  cout << "\n";
  cout << "  MU = " << mu << "\n";
  cout << "  SIGMA = " << sigma << "\n";
  cout << "  A = " << a << "\n";
  cout << "\n";
  cout << "   N   Estimate\n";
  cout << "\n";
//
//  N is the order of the quadrature rule.
//
  for ( n = 1; n <= 9; n++ )
  {
//
//  Compute the M = 2 * N + 1 moments.
//
    m = 2 * n + 1;

    moment = moments_truncated_normal_a ( m, mu, sigma, a );
//
//  Compute the points and weights by the method of moments.
//
    x1 = new double[n];
    w1 = new double[n];

    moment_method ( n, moment, x1, w1 );

    q = 0.0;
    for ( i = 0; i < n; i++ )
    {
      q = q + w1[i] * sin ( x1[i] );
    }
    cout << "  " << setw(2) << n
         << "  " << setw(14) << q << "\n";
//
//  Free memory.
//
    delete [] moment;
    delete [] w1;
    delete [] x1;
  }

  return;
}
