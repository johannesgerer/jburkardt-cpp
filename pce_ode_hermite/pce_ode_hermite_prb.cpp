# include <cstdlib>
# include <iostream>
# include <cmath>
# include <iomanip>
# include <ctime>

using namespace std;

# include "pce_ode_hermite.hpp"

int main ( );
void pce_ode_hermite_test01 ( );
void pce_ode_hermite_test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PCE_ODE_HERMITE_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PCE_ODE_HERMITE_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test PCE_ODE_HERMITE.\n";

  pce_ode_hermite_test01 ( );
  pce_ode_hermite_test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PCE_ODE_HERMITE_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void pce_ode_hermite_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    PCE_ODE_HERMITE_TEST01 runs a test problem with PCE_ODE_HERMITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double alpha_mu;
  double alpha_sigma;
  int i;
  int np = 4;
  int nt = 200;
  double *t;
  double tf;
  double ti;
  double *u;
  double *uex;
  double ui;

  t = new double[nt+1];
  u = new double[(nt+1)*(np+1)];
  uex = new double[nt+1];

  cout << "\n";
  cout << "PCE_ODE_HERMITE_TEST01:\n";
  cout << "  Call PCE_ODE_HERMITE to compute a polynomial chaos expansion\n";
  cout << "  for the ODE:\n";
  cout << "\n";
  cout << "    u' = - alpha * u,\n";
  cout << "    u(0) = 1.\n";

  ti = 0.0;
  tf = 2.0;
  ui = 1.0;
  alpha_mu = 0.0;
  alpha_sigma = 1.0;

  cout << "\n";
  cout << "  Initial time         TI = " << ti << "\n";
  cout << "  Final time           TF = " << tf << "\n";
  cout << "  Number of time steps NT = " << nt << "\n";
  cout << "  Initial condition    UI = " << ui << "\n";
  cout << "  Expansion degree     NP = " << np << "\n";
  cout << "  E(ALPHA)       ALPHA_MU = " << alpha_mu << "\n";
  cout << "  STD(ALPHA)  ALPHA_SIGMA = " << alpha_sigma << "\n";

  pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u );
//
//  Evaluate the exact expected value function.
//
  for ( i = 0; i <= nt; i++ )
  {
    uex[i] = ui * exp ( t[i] * t[i] / 2.0 );
  }
//
//  Compare the first computed component against the exact expected value.
//
  cout << "\n";
  cout << " i  T(i)  E(U(T(i)))    U(T(i),0)\n";
  cout << "\n";
  for ( i = 0; i <= nt; i = i + 10 )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(6) << t[i]
         << "  " << setw(14) << uex[i]
         << "  " << setw(14) << u[i+0*(nt+1)]
         << "  " << setw(14) << r8_abs ( uex[i] - u[i+0*(nt+1)] ) << "\n";
  }

  delete [] t;
  delete [] u;
  delete [] uex;

  return;
}
//****************************************************************************80

void pce_ode_hermite_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    PCE_ODE_HERMITE_TEST02 looks at convergence behavior for a fixed time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double alpha_mu;
  double alpha_sigma;
  double ep[6];
  int i;
  int np;
  int nt = 2000;
  double *t;
  double tf;
  double ti;
  double *u;
  double uexf;
  double ui;

  t = new double[nt+1];

  cout << "\n";
  cout << "PCE_ODE_HERMITE_TEST02:\n";
  cout << "  Examine convergence behavior as the PCE degree increases:\n";
  cout << "\n";
  cout << "    u' = - alpha * u,\n";
  cout << "    u(0) = 1.\n";

  ti = 0.0;
  tf = 2.0;
  ui = 1.0;
  alpha_mu = 0.0;
  alpha_sigma = 1.0;

  cout << "\n";
  cout << "  Initial time         TI = " << ti << "\n";
  cout << "  Final time           TF = " << tf << "\n";
  cout << "  Number of time steps NT = " << nt << "\n";
  cout << "  Initial condition    UI = " << ui << "\n";
  cout << "  E(ALPHA)       ALPHA_MU = " << alpha_mu << "\n";
  cout << "  STD(ALPHA)  ALPHA_SIGMA = " << alpha_sigma << "\n";

  uexf = ui * exp ( tf * tf / 2.0 );

  for ( np = 0; np <= 5; np++ )
  {
    u = new double[(nt+1)*(np+1)];

    pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u );

    ep[np] = r8_abs ( uexf - u[nt+0*(nt+1)] );

    delete [] u;
  }
//
//  Plot error in expected value as a function of the PCE degree.
//
  cout << "\n";
  cout << "    NP     Error(NP)     Log(Error(NP))\n";
  cout << "\n";
  for ( np = 0; np <= 5; np++ )
  {
    cout << "  " << setw(4) << np
         << "  " << setw(14) << ep[np]
         << "  " << setw(14) << log ( ep[np] ) << "\n";
  }

  delete [] t;

  return;
}
