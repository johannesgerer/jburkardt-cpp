# include <cstdlib>

using namespace std;

# include "rk4.hpp"

//****************************************************************************80

double rk4 ( double t0, double u0, double dt, double f ( double t, double u ) )

//****************************************************************************80
//
//  Purpose:
// 
//    RK4 takes one Runge-Kutta step.
//
//  Discussion:
//
//    It is assumed that an initial value problem, of the form
//
//      du/dt = f ( t, u )
//      u(t0) = u0
//
//    is being solved.
//
//    If the user can supply current values of t, u, a stepsize dt, and a
//    function to evaluate the derivative, this function can compute the
//    fourth-order Runge Kutta estimate to the solution at time t+dt.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T0, the current time.
//
//    Input, double U0, the solution estimate at the current time.
//
//    Input, double DT, the time step.
//
//    Input, double F ( double T, double U ), a function which evaluates
//    the derivative, or right hand side of the problem.
//
//    Output, double RK4, the fourth-order Runge-Kutta solution estimate
//    at time T0+DT.
//
{
  double f1;
  double f2;
  double f3;
  double f4;
  double u1;
//
//  Get four sample values of the derivative.
//
  f1 = f ( t0,          u0 );
  f2 = f ( t0 + dt / 2, u0 + dt * f1 / 2 );
  f3 = f ( t0 + dt / 2, u0 + dt * f2 / 2 );
  f4 = f ( t0 + dt,     u0 + dt * f3 );
//
//  Combine them to estimate the solution U1 at time T1 = T0 + DT.
//
  u1 = u0 + dt * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0;

  return u1;
}
