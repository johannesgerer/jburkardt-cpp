# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "ou.hpp"

int main ( );
void ou_euler_test ( );
void ou_euler_maruyama_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for OU_PRB.
//
//  Discussion:
//
//    OU_PRB tests the OU library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "OU_PRB:\n";
  cout << "  C++ version.\n";
  cout << "  Test the OU library.\n";

  ou_euler_test ( );
  ou_euler_maruyama_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "OU_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void ou_euler_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    OU_EULER_TEST tests OU_EULER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double mu;
  int n;
  int seed;
  double sigma;
  double theta;
  double tmax;
  double x0;

  cout << "\n";
  cout << "OU_EULER_TEST:\n";
  cout << "  Estimate a solution to the Ornstein-Uhlenbeck equation\n";
  cout << "  using the Euler method for stochastic differential equations.\n";
  cout << "\n";

  theta = 2.0;
  cout << "  Using decay rate THETA = " << theta << "\n";
  mu = 1.0;
  cout << "  Using mean MU = " << mu << "\n";
  sigma = 0.15;
  cout << "  Using variance SIGMA = " << sigma << "\n";
  x0 = 2.0;
  cout << "  Using initial value X0 = " << x0 << "\n";
  tmax = 3.0;
  cout << "  Using final time TMAX = " << tmax << "\n";
  n = 10000;
  cout << "  Using number of timesteps N = " << n << "\n";
  seed = 123456789;
  cout << "  Using value of random SEED = " << seed << "\n";

  ou_euler ( theta, mu, sigma, x0, tmax, n, seed );

  return;
}
//****************************************************************************80

void ou_euler_maruyama_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    OU_EULER_MARUYAMA_TEST tests OU_EULER_MARUYAMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double mu;
  int n;
  int r;
  int seed;
  double sigma;
  double theta;
  double tmax;
  double x0;

  cout << "\n";
  cout << "OU_EULER_MARUYAMA_TEST:\n";
  cout << "  Estimate a solution to the Ornstein-Uhlenbeck equation\n";
  cout << "  using the Euler-Maruyama method for stochastic \n";
  cout << "  differential equations.\n";
  cout << "\n";

  theta = 2.0;
  cout << "  Using decay rate THETA = " << theta << "\n";
  mu = 1.0;
  cout << "  Using mean MU = " << mu << "\n";
  sigma = 0.15;
  cout << "  Using variance SIGMA = " << sigma << "\n";
  x0 = 2.0;
  cout << "  Using initial value X0 = " << x0 << "\n";
  tmax = 3.0;
  cout << "  Using final time TMAX = " << tmax << "\n";
  n = 10000;
  cout << "  Using number of large timesteps N = " << n << "\n";
  r = 16;
  cout << "  Using number small time steps per one large time step R = " << r << "\n";
  seed = 123456789;
  cout << "  Using value of random SEED = " << seed << "\n";

  ou_euler_maruyama ( theta, mu, sigma, x0, tmax, n, r, seed );

  return;
}

