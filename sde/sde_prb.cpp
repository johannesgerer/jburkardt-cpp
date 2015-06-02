# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sde.hpp"

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

//****************************************************************************

int main ( )

//****************************************************************************
//
//  Purpose:
//
//    MAIN is the main program for SDE_PRB.
//
//  Discussion:
//
//    SDE_PRB tests the SDE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SDE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SDE library.\n";

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
//
//  Terminate.
//
  cout << "\n";
  cout << "SDE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************

void test01 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST01 tests BPATH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n = 500;
  int seed;
  double *w;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  BPATH generates a sample Brownian motion path.\n";

  seed = 123456789;

  w = bpath ( seed, n );

  bpath_gnuplot ( n, w );

  delete [] w;

  return;
}
//****************************************************************************

void test02 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST02 tests BPATH_AVERAGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  int m = 1000;
  int n = 500;
  int seed;
  double *u;
  double *umean;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  BPATH_AVERAGE generates many Brownian paths\n";
  cout << "  and averages them.\n";

  seed = 123456789;
  u = new double[m*(n+1)];
  umean = new double[n+1];

  bpath_average ( seed, m, n, u, umean, error );

  bpath_average_gnuplot ( m, n, u, umean );

  delete [] u;
  delete [] umean;

  return;
}
//****************************************************************************

void test03 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST03 tests CHAIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  int n = 200;
  int seed;
  double *vem;
  double *xem;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  CHAIN solves a stochastic differential equation for\n";
  cout << "  a function of a stochastic variable X.\n";
  cout << "  We can solve for X(t), and then evaluate V(X(t)).\n";
  cout << "  Or, we can apply the stochastic chain rule to derive an\n";
  cout << "  an SDE for V, and solve that.\n";

  seed = 123456789;
  xem = new double[n+1];
  vem = new double[n+1];
  chain ( seed, n, xem, vem, diff );

  cout << "\n";
  cout << "  Maximum | Sqrt(X) - V | = " << diff << "\n";

  chain_gnuplot ( n, xem, vem );

  delete [] vem;
  delete [] xem;

  return;
}
//****************************************************************************

void test04 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST04 tests EM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  int n = 256;
  int seed;
  double *t;
  double *t2;
  double *xem;
  double *xtrue;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  EM solves a stochastic differential equation\n";
  cout << "  using the Euler-Maruyama method.\n";

  seed = 123456789;

  t = new double[n+1];
  t2 = new double[1+n/4];
  xem = new double[1+n/4];
  xtrue = new double[n+1];

  em ( seed, n, t, xtrue, t2, xem, diff );

  cout << "\n";
  cout << "  | Exact X(T) - EM X(T) | = " << diff << "\n";

  em_gnuplot ( n, t, xtrue, t2, xem );

  delete [] t;
  delete [] t2;
  delete [] xem;
  delete [] xtrue;

  return;
}
//****************************************************************************

void test05 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST05 tests EMSTRONG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *dtvals;
  int m = 100;
  int n = 512;
  int p_max = 6;
  int seed;
  double *xerr;

  dtvals = new double[p_max];
  xerr = new double[p_max];

  cout << "\n";
  cout << "TEST05:\n";
  cout << "  EMSTRONG investigates the strong convergence\n";
  cout << "  of the Euler-Maruyama method.\n";

  seed = 123456789;

  emstrong ( seed, m, n, p_max, dtvals, xerr );

  emstrong_gnuplot ( p_max, dtvals, xerr );

  delete [] dtvals;
  delete [] xerr;

  return;
}
//****************************************************************************

void test06 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST06 tests EMWEAK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *dtvals;
  int m = 50000;
  int method;
  int p_max = 5;
  int seed;
  double *xerr;

  dtvals = new double[p_max];
  xerr = new double[p_max];

  cout << "\n";
  cout << "TEST06:\n";
  cout << "  EMWEAK investigates the weak convergence\n";
  cout << "  of the Euler-Maruyama method.\n";

  seed = 123456789;
  method = 0;

  emweak ( seed, method, m, p_max, dtvals, xerr );

  emweak_gnuplot ( p_max, dtvals, xerr, method );

  seed = 123456789;
  method = 1;

  emweak ( seed, method, m, p_max, dtvals, xerr );

  emweak_gnuplot ( p_max, dtvals, xerr, method );

  delete [] dtvals;
  delete [] xerr;

  return;
}
//****************************************************************************

void test07 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST07 tests MILSTRONG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *dtvals;
  int p_max = 4;
  int seed;
  double *xerr;

  dtvals = new double[p_max];
  xerr = new double[p_max];

  cout << "\n";
  cout << "TEST07:\n";
  cout << "  MILSTRONG investigates the strong convergence\n";
  cout << "  of the Milstein method.\n";

  seed = 123456789;

  milstrong ( seed, p_max, dtvals, xerr );

  milstrong_gnuplot ( p_max, dtvals, xerr );

  delete [] dtvals;
  delete [] xerr;

  return;
}
//****************************************************************************

void test08 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST08 tests STAB_ASYMPTOTIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n = 1000;
  int p_max = 3;
  int seed;

  cout << "\n";
  cout << "TEST08:\n";
  cout << "  STAB_ASYMPTOTIC investigates the asymptotic\n";
  cout << "  stability of the Euler-Maruyama method.\n";
  cout << "\n";
  cout << "  For technical reasons, the plotting is done\n";
  cout << "  in the same routine as the computations.\n";

  seed = 123456789;

  stab_asymptotic ( seed, n, p_max );

  return;
}
//****************************************************************************

void test09 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST09 tests STAB_MEANSQUARE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  int seed;

  cout << "\n";
  cout << "TEST09:\n";
  cout << "  STAB_MEANSQUARE investigates the mean square\n";
  cout << "  stability of the Euler-Maruyama method.\n";
  cout << "\n";
  cout << "  For technical reasons, the plotting is done\n";
  cout << "  in the same routine as the computations.\n";

  seed = 123456789;

  stab_meansquare ( seed );

  return;
}
//****************************************************************************

void test10 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST10 tests STOCHASTIC_INTEGRAL_ITO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  int i;
  int n;
  int seed;

  cout << "\n";
  cout << "TEST10:\n";
  cout << "  Estimate the Ito integral of W(t) dW over [0,1].\n";
  cout << "\n";
  cout << "                                                 Abs          Rel\n";
  cout << "         N        Exact        Estimate          Error        Error\n";
  cout << "\n";

  n = 100;
  seed = 123456789;

  for ( i = 1; i <= 7; i++ )
  {
    stochastic_integral_ito ( n, seed, estimate, exact, error );

    cout << "  " << setw(8) << n
         << "  " << setw(16) << exact
         << "  " << setw(16) << estimate
         << "  " << setw(16) << error
         << "  " << setw(10) << error / exact << "\n";

    n = n * 4;
  }
  return;
}
//****************************************************************************

void test11 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST11 tests STOCHASTIC_INTEGRAL_STRAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  int i;
  int n;
  int seed;

  cout << "\n";
  cout << "TEST11:\n";
  cout << "  Estimate the Stratonovich integral of W(t) dW over [0,1].\n";
  cout << "\n";
  cout << "                                                 Abs          Rel\n";
  cout << "         N        Exact        Estimate          Error        Error\n";
  cout << "\n";

  n = 100;
  seed = 123456789;

  for ( i = 1; i <= 7; i++ )
  {
    stochastic_integral_strat ( n, seed, estimate, exact, error );

    cout << "  " << setw(8) << n
         << "  " << setw(16) << exact
         << "  " << setw(16) << estimate
         << "  " << setw(16) << error
         << "  " << setw(10) << error / exact << "\n";

    n = n * 4;
  }
  return;
}
