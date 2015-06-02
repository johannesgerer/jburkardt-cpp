# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "rk4.hpp"

int main ( );
void test01 ( );
double test01_f ( double t, double u );
void test02 ( );
double *test02_f ( double t, int n, double u[] );

//****************************************************************************80

int main ( ) 

//****************************************************************************80
//
//  Purpose:
// 
//    MAIN is the main program for RK4_PRB.
//
//  Discussion:
//
//    RK4_PRB tests the RK4 library.
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
{
  timestamp ( );
  printf ( "\n" );
  printf ( "RK4_PRB\n" );
  printf ( "  C++ version\n" );
  printf ( "  Test the RK4 library.\n" );

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "RK4_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( ) 

//****************************************************************************80
//
//  Purpose:
// 
//    TEST01 demonstrates the use of RK4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2013
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

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  RK4 takes one Runge Kutta step for a scalar ODE.\n" );

  printf ( "\n" );
  printf ( "          T          U[T]\n" );
  printf ( "\n" );

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
    u1 = rk4 ( t0, u0, dt, test01_f );
//
//  Shift the data to prepare for another step.
//
    t0 = t1;
    u0 = u1;
  }
  return;
}
//****************************************************************************80

double test01_f ( double t, double u )

//****************************************************************************80
//
//  Purpose:
// 
//    TEST01_F evaluates the right hand side of a particular ODE.
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
//    Output, double TEST01_F, the value of the derivative, dU/dT.
//
{
  double dudt;
  
  dudt = u * cos ( t );
  
  return dudt;
}  
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests the RK4 routine for a vector ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double dt = 0.2;
  int i;
  int n = 2;
  double t0;
  double t1;
  double tmax = 12.0 * 3.141592653589793;
  double *u0;
  double *u1;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  RK4VEC takes a Runge Kutta step for a vector ODE.\n";

  cout << "\n";
  cout << "       T       U[0]       U[1]\n";
  cout << "\n";
  t0 = 0.0;

  u0 = new double[n];
  u0[0] = 0.0;
  u0[1] = 1.0;

  for ( ; ; )
  {
//
//  Print (T0,U0).
//
    cout << "  " << setw(14) << t0
         << "  " << setw(14) << u0[0]
         << "  " << setw(14) << u0[1] << "\n";
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
    u1 = rk4vec ( t0, n, u0, dt, test02_f );
//
//  Shift the data to prepare for another step.
//
    t0 = t1;
    for ( i = 0; i < n; i++ )
    {
      u0[i] = u1[i];
    }
    delete [] u1;
  }
  return;
}
//****************************************************************************80

double *test02_f ( double t, int n, double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02_F evaluates the right hand side of a vector ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the current time.
//
//    Input, int N, the dimension of the system.
//
//    Input, double U[N], the current solution value.
//
//    Output, double TEST02_F[N], the value of the derivative, dU/dT.
//
{
  double *uprime;

  uprime = new double[n];

  uprime[0] =   u[1];
  uprime[1] = - u[0];
 
  return uprime;
}
