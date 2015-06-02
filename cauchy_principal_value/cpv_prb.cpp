# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "cpv.hpp"

int main ( );
void cpv_test01 ( );
double f01 ( double t );
void cpv_test02 ( );
double f02 ( double t );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    CPV_PRB tests the CPV library.
//
//  Location:
//
//    http://people.sc.fsu.edu/~jburkardt/cpp_src/cauchy_principal_value/cpv_prb.cpp
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CPV_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CPV library.\n";

  cpv_test01 ( );
  cpv_test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CPV_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void cpv_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    CPV_TEST01 seeks the CPV of Integral ( -1 <= t <= 1 ) exp(t) / t dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double exact;
  int n;
  double value;

  cout << "\n";
  cout << "CPV_TEST01:\n";
  cout << "  CPV of Integral ( -1 <= t <= 1 ) exp(t) / t dt\n";

  cout << "\n";
  cout << "   N           Estimate             Error\n";
  cout << "\n";

  exact = 2.11450175075;
  a = -1.0;
  b = +1.0;
  for ( n = 2; n <= 8; n = n + 2 )
  {
    value = cpv ( f01, a, b, n );
    cout << "  " << setw(2) << n
         << "  " << setw(24) << value
         << "  " << setw(14) <<  fabs ( value - exact ) << "\n";
  }

  return;
}
//****************************************************************************80

double f01 ( double t )

//****************************************************************************80
//
//  Purpose:
//
//    F01 evaluates the integrand of Integral ( -1 <= t <= 1 ) exp(t) / t dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real T, the argument.
//
//    Output, real VALUE, the value of the integrand.
//
{
  double value;

  value = exp ( t );

  return value;
}
//****************************************************************************80

void cpv_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    CPV_TEST02 is another test.
//
// Discussion:
//
//    We seek
//      CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1-t)^3 dt )
//    which we must rewrite as
//      CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1+t+t^2) 1/(1-t) dt )
//    so that our "integrand" is 1/(1+t+t^2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double delta;
  double exact;
  int k;
  int n;
  double r1;
  double r2;
  double r3;
  double r4;
  double value;

  cout << "\n";
  cout << "CPV_TEST02:\n";
  cout << "  Compute CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1-t)^3 dt )\n";
  cout << "  Try this for delta = 1, 1/2, 1/4.\n";
  cout << "\n";
  cout << "   N          Estimate                  Exact                  Error\n";
  delta = 1.0;
  for ( k = 1; k <= 3; k++ )
  {
    cout << "\n";
    r1 = pow (   delta + 1.5, 2 ) + 0.75;
    r2 = pow ( - delta + 1.5, 2 ) + 0.75;
    r3 = atan ( sqrt ( 0.75 ) / (   delta + 1.5 ) );
    r4 = atan ( sqrt ( 0.75 ) / ( - delta + 1.5 ) );
    exact = - log ( r1 / r2 ) / 6.0 + ( r3 - r4 ) / sqrt ( 3.0 );
    for ( n = 2; n <= 8; n = n + 2 )
    {
      a = 1.0 - delta;
      b = 1.0 + delta;
      value = cpv ( f02, a, b, n );
      cout << "  " << setw(2) << n
           << "  " << setw(24) << value
           << "  " << setw(24) << exact
           << "  " << setw(14) <<  fabs ( value - exact ) << "\n";
    }
    delta = delta / 2.0;
  }

  return;
}
//****************************************************************************80

double f02 ( double t )

//****************************************************************************80
//
//  Purpose:
//
//    F02: integrand of Integral ( 1-delta <= t <= 1+delta ) 1/(1-t^3) dt
//
//  Discussion:
//
//    1/(1-t^3) = 1/(1+t+t^2) * 1/(1-t)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the evaluation point.
//
//    Output, double F02, the value of the integrand at T.
//
{
  double value;

  value = 1.0 / ( 1.0 + t + t * t );

  return value;
}
