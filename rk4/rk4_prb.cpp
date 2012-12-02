# include <cstdlib>
# include <iostream>
# include <cmath>

using namespace std;

# include "rk4.hpp"

double f3 ( double t, double u );

//****************************************************************************80

int main ( ) 

//****************************************************************************80
//
//  Purpose:
// 
//    RK4_PRB demonstrates the use of the RK4 one-step ODE solver.
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
//    Local, double DT, the time step.
//
//    Local, double T0, the time at which the solution is known.
//
//    Local, double TMAX, the maximum time at which a solution is desired.
//
//    Local, double U0, the estimated solution at time T0.
//
{
  double dt = 0.1;
  double pi = 3.14159265;
  double t0 = 0.0;
  double t1;
  double tmax = 12.0 * pi;
  double u0 = 0.5;
  double u1;

  while ( true )
  {
//
//  Print (T0,U0).
//
    cout << "  " << t0 << "  " << u0 << "\n";
//
//  Stop if we've exceeded TMAX.
//
    if ( tmax <= t0 )
    {
      break;
    }
//
//  Otherwise, advance to time T1, and have RK4 estimate 
//  the solution U1 there.
//
    t1 = t0 + dt;
    u1 = rk4 ( t0, u0, dt, f3 );
//
//  Shift the data to prepare for another step.
//
    t0 = t1;
    u0 = u1;
  }
  return 0;
}
//****************************************************************************80

double f3 ( double t, double u )

//****************************************************************************80
//
//  Purpose:
// 
//    F3 evaluates the right hand side of a particular ODE.
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
//    Input, double T, the current time.
//
//    Input, double U, the current solution value.
//
//    Output, double F3, the value of the derivative, dU/dT.
//
{
  double dudt;
  
  dudt = u * cos ( t );
  
  return dudt;
}  
