# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "stochastic_rk.hpp"

int main ( );
void test01 ( );
double fi ( double x );
double gi ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for STOCHASTIC_RK_PRB.
//
//  Discussion:
//
//    STOCHASTIC_RK_PRB tests the STOCHASTIC_RK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "STOCHASTIC_RK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the STOCHASTIC_RK library.\n";
 
  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "STOCHASTIC_RK_PRB\n";
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
//    TEST01 tests RK1_TI_STEP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double h;
  int i;
  int n;
  double q;
  int seed;
  double t;
  double t0 = 0.0;
  double tn = 1.0;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  RK1_TI_STEP uses a first order RK method\n";
  cout << "  for a problem whose right hand side does not\n";
  cout << "  depend explicitly on time.\n";

  n = 10;
  x = new double[n+1];
  h = ( tn - t0 ) / ( double ) ( n );
  q = 1.0;
  seed = 123456789;

  i = 0;
  t = t0;
  x[i] = 0.0;

  cout << "\n";
  cout << "         I           T             X\n";
  cout << "\n";
  cout << "  " << setw(8) << i
       << "  " << setw(14) << t
       << "  " << setw(14) << x[i] << "\n";

  for ( i = 1; i <= n; i++ )
  {
    t = ( ( double ) ( n - i ) * t0   
        + ( double ) (     i ) * tn )
        / ( double ) ( n     );

    x[i] = rk1_ti_step ( x[i-1], t, h, q, fi, gi, &seed );

    cout << "  " << setw(8) << i
         << "  " << setw(14) << t
         << "  " << setw(14) << x[i] << "\n";
  }

  delete [] x;

  return;
}
//****************************************************************************80

double fi ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FI is a time invariant deterministic right hand side.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double FI, the value.
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double gi ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    GI is a time invariant stochastic right hand side.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double GI, the value.
//
{
  double value;

  value = 1.0;

  return value;
}
