# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "black_scholes.hpp"

int main ( );
void asset_path_test ( );
void binomial_test ( );
void bsf_test ( );
void forward_test ( );
void mc_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLACK_SCHOLES_PRB.
//
//  Discussion:
//
//    BLACK_SCHOLES_PRB tests the BLACK_SCHOLES library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BLACK_SCHOLES_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BLACK_SCHOLES library.\n";

  asset_path_test ( );
  binomial_test ( );
  bsf_test ( );
  forward_test ( );
  mc_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BLACK_SCHOLES_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void asset_path_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ASSET_PATH_TEST tests ASSET_PATH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double mu;
  int n = 100;
  string output_filename;
  double *s;
  double s0;
  int seed;
  double sigma;
  double t1;

  cout << "\n";
  cout << "ASSET_PATH_TEST:\n";
  cout << "  Demonstrate the simulated of an asset price path.\n";

  s0 = 2.0;
  mu = 0.1;
  sigma = 0.3;
  t1 = 1.0;
  seed = 123456789;

  cout << "\n";
  cout << "  The asset price at time 0      S0    = " << s0 << "\n";
  cout << "  The asset expected growth rate MU    = " << mu << "\n";
  cout << "  The asset volatility           SIGMA = " << sigma << "\n";
  cout << "  The expiry date                T1    = " << t1 << "\n";
  cout << "  The number of time steps       N     = " << n << "\n";
  cout << "  The random number seed was     SEED  = " << seed << "\n";

  s = asset_path ( s0, mu, sigma, t1, n, &seed );

  r8vec_print_part ( n + 1, s, 10, "  Partial results:" );

  output_filename = "asset_path.txt";
  r8vec_write ( output_filename, n + 1, s );

  cout << "\n";
  cout << "  Full results written to \"" << output_filename << "\".\n";

  delete [] s;

  return;
}
//****************************************************************************80

void binomial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_TEST tests BINOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  double e;
  int m;
  double r;
  double s0;
  double sigma;
  double t1;

  cout << "\n";
  cout << "BINOMIAL_TEST:\n";
  cout << "  A demonstration of the binomial method\n";
  cout << "  for option valuation.\n";

  s0 = 2.0;
  e = 1.0;
  r = 0.05;
  sigma = 0.25;
  t1 = 3.0;
  m = 256;

  cout << "\n";
  cout << "  The asset price at time 0 S0    = " << s0 << "\n";
  cout << "  The exercise price        E     = " << e << "\n";
  cout << "  The interest rate         R     = " << r << "\n";
  cout << "  The asset volatility      SIGMA = " << sigma << "\n";
  cout << "  The expiry date           T1    = " << t1 << "\n";
  cout << "  The number of intervals   M     = " << m << "\n";

  c = binomial ( s0, e, r, sigma, t1, m );

  cout << "\n";
  cout << "  The option value is " << c << "\n";

  return;
}
//****************************************************************************80

void bsf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BSF_TEST tests BSF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  double e;
  double r;
  double s0;
  double sigma;
  double t0;
  double t1;

  cout << "\n";
  cout << "BSF_TEST:\n";
  cout << "  A demonstration of the Black-Scholes formula\n";
  cout << "  for option valuation.\n";

  s0 = 2.0;
  t0 = 0.0;
  e = 1.0;
  r = 0.05;
  sigma = 0.25;
  t1 = 3.0;

  cout << "\n";
  cout << "  The asset price at time T0 S0    = " << s0 << "\n";
  cout << "  The time                   T0    = " << t0 << "\n";
  cout << "  The exercise price         E     = " << e << "\n";
  cout << "  The interest rate          R     = " << r << "\n";
  cout << "  The asset volatility       SIGMA = " << sigma << "\n";
  cout << "  The expiry date            T1    = " << t1 << "\n";

  c = bsf ( s0, t0, e, r, sigma, t1 );

  cout << "\n";
  cout << "  The option value C = " << c << "\n";

  return;
}
//****************************************************************************80

void forward_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FORWARD_TEST tests FORWARD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double e;
  int i;
  int nt;
  int nx;
  double r;
  double s;
  double sigma;
  double smax;
  double smin;
  double t1;
  double *u;

  cout << "\n";
  cout << "FORWARD_TEST:\n";
  cout << "  A demonstration of the forward difference method\n";
  cout << "  for option valuation.\n";

  e = 4.0;
  r = 0.03;
  sigma = 0.50;
  t1 = 1.0;
  nx = 11;
  nt = 29;
  smax = 10.0;

  cout << "\n";
  cout << "  The exercise price        E =     " << e << "\n";
  cout << "  The interest rate         R =     " << r << "\n";
  cout << "  The asset volatility      SIGMA = " << sigma << "\n";
  cout << "  The expiry date           T1 =    " << t1 << "\n";
  cout << "  The number of space steps NX =    " << nx << "\n";
  cout << "  The number of time steps  NT =    " << nt << "\n";
  cout << "  The value of              SMAX =  " << smax << "\n";

  u = forward ( e, r, sigma, t1, nx, nt, smax );

  cout << "\n";
  cout << "         Initial          Option\n";
  cout << "           Value           Value\n"; 
  cout << "\n";

  smin = 0.0;
  for ( i = 0; i < nx - 1; i++ )
  {
    s = ( ( nx - i - 2 ) * smin +  ( i + 1 ) * smax ) / ( double ) ( nx - 1 );
    cout << "  " << setw(14) << s
         << "  " << setw(14) << u[i+nt*(nx-1)] << "\n";
  }

  delete [] u;

  return;
}
//****************************************************************************80

void mc_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MC_TEST tests MC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *conf;
  double e;
  int m;
  double r;
  double s0;
  int seed;
  double sigma;
  double t1;

  cout << "\n";
  cout << "MC_TEST:\n";
  cout << "  A demonstration of the Monte Carlo method\n";
  cout << "  for option valuation.\n";

  s0 = 2.0;
  e = 1.0;
  r = 0.05;
  sigma = 0.25;
  t1 = 3.0;
  m = 1000000;
  seed = 123456789;

  cout << "\n";
  cout << "  The asset price at time 0, S0    = " << s0 << "\n";
  cout << "  The exercise price         E     = " << e << "\n";
  cout << "  The interest rate          R     = " << r << "\n";
  cout << "  The asset volatility       SIGMA = " << sigma << "\n";
  cout << "  The expiry date            T1    = " << t1 << "\n";
  cout << "  The number of simulations  M     = " << m << "\n";

  conf = mc ( s0, e, r, sigma, t1, m, &seed );

  cout << "\n";
  cout << "  The confidence interval is [" << conf[0] << ", " << conf[1] << "].\n";

  return;
}
