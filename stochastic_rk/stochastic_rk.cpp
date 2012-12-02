# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "stochastic_rk.hpp"

//****************************************************************************80

double rk1_ti_step ( double x, double t, double h, double q, 
  double fi ( double x ), double gi ( double x ), int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RK1_TI_STEP takes one step of a stochastic Runge Kutta scheme.
//
//  Discussion:
//
//    The Runge-Kutta scheme is first-order, and suitable for time-invariant
//    systems in which F and G do not depend explicitly on time.
//
//    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
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
//  Reference:
//
//    Jeremy Kasdin,
//    Runge-Kutta algorithm for the numerical integration of
//    stochastic differential equations,
//    Journal of Guidance, Control, and Dynamics,
//    Volume 18, Number 1, January-February 1995, pages 114-120.
//
//    Jeremy Kasdin,
//    Discrete Simulation of Colored Noise and Stochastic Processes
//    and 1/f^a Power Law Noise Generation,
//    Proceedings of the IEEE,
//    Volume 83, Number 5, 1995, pages 802-827.
//
//  Parameters:
//
//    Input, double X, the value at the current time.
//
//    Input, double T, the current time.
//
//    Input, double H, the time step.
//
//    Input, double Q, the spectral density of the input white noise.
//
//    Input, double FI ( double X ), the name of the deterministic
//    right hand side function.
//
//    Input, double GI ( double X ), the name of the stochastic
//    right hand side function.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double RK1_TI_STEP, the value at time T+H.
//
{
  double a21;
  double k1;
  double q1;
  double t1;
  double w1;
  double x1;
  double xstar;

  a21 = 1.0;

  q1 = 1.0;

  t1 = t;
  x1 = x;
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1;

  xstar = x1 + a21 * k1;

  return xstar;
}
//****************************************************************************80

double rk2_ti_step ( double x, double t, double h, double q, 
  double fi ( double x ), double gi ( double x ), int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RK2_TI_STEP takes one step of a stochastic Runge Kutta scheme.
//
//  Discussion:
//
//    The Runge-Kutta scheme is second-order, and suitable for time-invariant
//    systems.
//
//    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
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
//  Reference:
//
//    Jeremy Kasdin,
//    Runge-Kutta algorithm for the numerical integration of
//    stochastic differential equations,
//    Journal of Guidance, Control, and Dynamics,
//    Volume 18, Number 1, January-February 1995, pages 114-120.
//
//    Jeremy Kasdin,
//    Discrete Simulation of Colored Noise and Stochastic Processes
//    and 1/f^a Power Law Noise Generation,
//    Proceedings of the IEEE,
//    Volume 83, Number 5, 1995, pages 802-827.
//
//  Parameters:
//
//    Input, double X, the value at the current time.
//
//    Input, double T, the current time.
//
//    Input, double H, the time step.
//
//    Input, double Q, the spectral density of the input white noise.
//
//    Input, double FI ( double X ), the name of the deterministic
//    right hand side function.
//
//    Input, double GI ( double X ), the name of the stochastic
//    right hand side function.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double RK2_TI_STEP, the value at time T+H.
//
{
  double a21;
  double a31;
  double a32;
  double k1;
  double k2;
  double q1;
  double q2;
  double t1;
  double t2;
  double w1;
  double w2;
  double x1;
  double x2;
  double xstar;

  a21 =   1.0;
  a31 =   0.5;
  a32 =   0.5;

  q1 = 2.0;
  q2 = 2.0;

  t1 = t;
  x1 = x;
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1;

  t2 = t1 + a21 * h;
  x2 = x1 + a21 * k1;
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2;

  xstar = x1 + a31 * k1 + a32 * k2;

  return xstar;
}
//****************************************************************************80

double rk3_ti_step ( double x, double t, double h, double q, 
  double fi ( double x ), double gi ( double x ), int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RK3_TI_STEP takes one step of a stochastic Runge Kutta scheme.
//
//  Discussion:
//
//    The Runge-Kutta scheme is third-order, and suitable for time-invariant
//    systems in which F and G do not depend explicitly on time.
//
//    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
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
//  Reference:
//
//    Jeremy Kasdin,
//    Runge-Kutta algorithm for the numerical integration of
//    stochastic differential equations,
//    Journal of Guidance, Control, and Dynamics,
//    Volume 18, Number 1, January-February 1995, pages 114-120.
//
//    Jeremy Kasdin,
//    Discrete Simulation of Colored Noise and Stochastic Processes
//    and 1/f^a Power Law Noise Generation,
//    Proceedings of the IEEE,
//    Volume 83, Number 5, 1995, pages 802-827.
//
//  Parameters:
//
//    Input, double X, the value at the current time.
//
//    Input, double T, the current time.
//
//    Input, double H, the time step.
//
//    Input, double Q, the spectral density of the input white noise.
//
//    Input, double FI ( double X ), the name of the deterministic
//    right hand side function.
//
//    Input, double GI ( double X ), the name of the stochastic
//    right hand side function.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double RK3_TI_STEP, the value at time T+H.
//
{
  double a21;
  double a31;
  double a32;
  double a41;
  double a42;
  double a43;
  double k1;
  double k2;
  double k3;
  double q1;
  double q2;
  double q3;
  double t1;
  double t2;
  double t3;
  double w1;
  double w2;
  double w3;
  double x1;
  double x2;
  double x3;
  double xstar;

  a21 =   1.52880952525675;
  a31 =   0.0;
  a32 =   0.51578733443615;
  a41 =   0.53289582961739;
  a42 =   0.25574324768195;
  a43 =   0.21136092270067;

  q1 = 1.87653936176981;
  q2 = 3.91017166264989;
  q3 = 4.73124353935667;

  t1 = t;
  x1 = x;
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1;

  t2 = t1 + a21 * h;
  x2 = x1 + a21 * k1;
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2;

  t3 = t1 + a31 * h  + a32 * h;
  x3 = x1 + a31 * k1 + a32 * k2;
  w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h );
  k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3;

  xstar = x1 + a41 * k1 + a42 * k2 + a43 * k3;

  return xstar;
}
//****************************************************************************80

double rk4_ti_step ( double x, double t, double h, double q, 
  double fi ( double x ), double gi ( double x ), int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RK4_TI_STEP takes one step of a stochastic Runge Kutta scheme.
//
//  Discussion:
//
//    The Runge-Kutta scheme is fourth-order, and suitable for time-invariant
//    systems in which F and G do not depend explicitly on time.
//
//    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
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
//  Reference:
//
//    Jeremy Kasdin,
//    Runge-Kutta algorithm for the numerical integration of
//    stochastic differential equations,
//    Journal of Guidance, Control, and Dynamics,
//    Volume 18, Number 1, January-February 1995, pages 114-120.
//
//    Jeremy Kasdin,
//    Discrete Simulation of Colored Noise and Stochastic Processes
//    and 1/f^a Power Law Noise Generation,
//    Proceedings of the IEEE,
//    Volume 83, Number 5, 1995, pages 802-827.
//
//  Parameters:
//
//    Input, double X, the value at the current time.
//
//    Input, double T, the current time.
//
//    Input, double H, the time step.
//
//    Input, double Q, the spectral density of the input white noise.
//
//    Input, double FI ( double X ), the name of the deterministic
//    right hand side function.
//
//    Input, double GI ( double X ), the name of the stochastic
//    right hand side function.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double RK4_TI_STEP, the value at time T+H.
//
{
  double a21;
  double a31;
  double a32;
  double a41;
  double a42;
  double a43;
  double a51;
  double a52;
  double a53;
  double a54;
  double k1;
  double k2;
  double k3;
  double k4;
  double q1;
  double q2;
  double q3;
  double q4;
  double t1;
  double t2;
  double t3;
  double t4;
  double w1;
  double w2;
  double w3;
  double w4;
  double x1;
  double x2;
  double x3;
  double x4;
  double xstar;

  a21 =   2.71644396264860;
  a31 = - 6.95653259006152;
  a32 =   0.78313689457981;
  a41 =   0.0;
  a42 =   0.48257353309214;
  a43 =   0.26171080165848;
  a51 =   0.47012396888046;
  a52 =   0.36597075368373;
  a53 =   0.08906615686702;
  a54 =   0.07483912056879;

  q1 =   2.12709852335625;
  q2 =   2.73245878238737;
  q3 =  11.22760917474960;
  q4 =  13.36199560336697;

  t1 = t;
  x1 = x;
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1;

  t2 = t1 + a21 * h;
  x2 = x1 + a21 * k1;
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2;

  t3 = t1 + a31 * h  + a32 * h;
  x3 = x1 + a31 * k1 + a32 * k2;
  w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h );
  k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3;

  t4 = t1 + a41 * h  + a42 * h + a43 * h;
  x4 = x1 + a41 * k1 + a42 * k2;
  w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h );
  k4 = h * fi ( x4 ) + h * gi ( x4 ) * w4;

  xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;

  return xstar;
}
//****************************************************************************80

double rk1_tv_step ( double x, double t, double h, double q, 
  double fv ( double t, double x ), double gv ( double t, double x ), 
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RK1_TV_STEP takes one step of a stochastic Runge Kutta scheme.
//
//  Discussion:
//
//    The Runge-Kutta scheme is first-order, and suitable for time-varying
//    systems.
//
//    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
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
//  Reference:
//
//    Jeremy Kasdin,
//    Runge-Kutta algorithm for the numerical integration of
//    stochastic differential equations,
//    Journal of Guidance, Control, and Dynamics,
//    Volume 18, Number 1, January-February 1995, pages 114-120.
//
//    Jeremy Kasdin,
//    Discrete Simulation of Colored Noise and Stochastic Processes
//    and 1/f^a Power Law Noise Generation,
//    Proceedings of the IEEE,
//    Volume 83, Number 5, 1995, pages 802-827.
//
//  Parameters:
//
//    Input, double X, the value at the current time.
//
//    Input, double T, the current time.
//
//    Input, double H, the time step.
//
//    Input, double Q, the spectral density of the input white noise.
//
//    Input, double FV ( double T, double X ), the name of the deterministic
//    right hand side function.
//
//    Input, double GV ( double T, double X ), the name of the stochastic
//    right hand side function.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double RK1_TV_STEP the value at time T+H.
//
{
  double a21;
  double k1;
  double q1;
  double t1;
  double w1;
  double x1;
  double xstar;

  a21 = 1.0;

  q1 = 1.0;

  t1 = t;
  x1 = x;
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1;

  xstar = x1 + a21 * k1;

  return xstar;
}
//****************************************************************************80

double rk2_tv_step ( double x, double t, double h, double q, 
  double fv ( double t, double x ), double gv ( double t, double x ), 
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RK2_TV_STEP takes one step of a stochastic Runge Kutta scheme.
//
//  Discussion:
//
//    The Runge-Kutta scheme is second-order, and suitable for time-varying
//    systems.
//
//    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
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
//  Reference:
//
//    Jeremy Kasdin,
//    Runge-Kutta algorithm for the numerical integration of
//    stochastic differential equations,
//    Journal of Guidance, Control, and Dynamics,
//    Volume 18, Number 1, January-February 1995, pages 114-120.
//
//    Jeremy Kasdin,
//    Discrete Simulation of Colored Noise and Stochastic Processes
//    and 1/f^a Power Law Noise Generation,
//    Proceedings of the IEEE,
//    Volume 83, Number 5, 1995, pages 802-827.
//
//  Parameters:
//
//    Input, double X, the value at the current time.
//
//    Input, double T, the current time.
//
//    Input, double H, the time step.
//
//    Input, double Q, the spectral density of the input white noise.
//
//    Input, double FV ( double T, double X ), the name of the deterministic
//    right hand side function.
//
//    Input, double GV ( double T, double X ), the name of the stochastic
//    right hand side function.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double RK2_TV_STEP, the value at time T+H.
//
{
  double a21;
  double a31;
  double a32;
  double k1;
  double k2;
  double q1;
  double q2;
  double t1;
  double t2;
  double w1;
  double w2;
  double x1;
  double x2;
  double xstar;

  a21 =   1.0;
  a31 =   0.5;
  a32 =   0.5;

  q1 = 2.0;
  q2 = 2.0;

  t1 = t;
  x1 = x;
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1;

  t2 = t1 + a21 * h;
  x2 = x1 + a21 * k1;
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
  k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2;

  xstar = x1 + a31 * k1 + a32 * k2;

  return xstar;
}
//****************************************************************************80

double rk4_tv_step ( double x, double t, double h, double q, 
  double fv ( double t, double x ), double gv ( double t, double x ), 
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RK4_TV_STEP takes one step of a stochastic Runge Kutta scheme.
//
//  Discussion:
//
//    The Runge-Kutta scheme is fourth-order, and suitable for time-varying
//    systems.
//
//    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
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
//  Reference:
//
//    Jeremy Kasdin,
//    Runge-Kutta algorithm for the numerical integration of
//    stochastic differential equations,
//    Journal of Guidance, Control, and Dynamics,
//    Volume 18, Number 1, January-February 1995, pages 114-120.
//
//    Jeremy Kasdin,
//    Discrete Simulation of Colored Noise and Stochastic Processes
//    and 1/f^a Power Law Noise Generation,
//    Proceedings of the IEEE,
//    Volume 83, Number 5, 1995, pages 802-827.
//
//  Parameters:
//
//    Input, double X, the value at the current time.
//
//    Input, double T, the current time.
//
//    Input, double H, the time step.
//
//    Input, double Q, the spectral density of the input white noise.
//
//    Input, double FV ( double T, double X ), the name of the deterministic
//    right hand side function.
//
//    Input, double GV ( double T, double X ), the name of the stochastic
//    right hand side function.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double RK4_TV_STEP, the value at time T+H.
//
{
  double a21;
  double a31;
  double a32;
  double a41;
  double a42;
  double a43;
  double a51;
  double a52;
  double a53;
  double a54;
  double k1;
  double k2;
  double k3;
  double k4;
  double q1;
  double q2;
  double q3;
  double q4;
  double t1;
  double t2;
  double t3;
  double t4;
  double w1;
  double w2;
  double w3;
  double w4;
  double x1;
  double x2;
  double x3;
  double x4;
  double xstar;

  a21 =   0.66667754298442;
  a31 =   0.63493935027993;
  a32 =   0.00342761715422;
  a41 = - 2.32428921184321;
  a42 =   2.69723745129487;
  a43 =   0.29093673271592;
  a51 =   0.25001351164789;
  a52 =   0.67428574806272;
  a53 = - 0.00831795169360;
  a54 =   0.08401868181222;

  q1 = 3.99956364361748;
  q2 = 1.64524970733585;
  q3 = 1.59330355118722;
  q4 = 0.26330006501868;

  t1 = t;
  x1 = x;
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1;

  t2 = t1 + a21 * h;
  x2 = x1 + a21 * k1;
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
  k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2;

  t3 = t1 + a31 * h  + a32 * h;
  x3 = x1 + a31 * k1 + a32 * k2;
  w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h );
  k3 = h * fv ( t3, x3 ) + h * gv ( t3, x3 ) * w3;

  t4 = t1 + a41 * h  + a42 * h  + a43 * h;
  x4 = x1 + a41 * k1 + a42 * k2 + a43 * k3;
  w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h );
  k4 = h * fv ( t4, x4 ) + h * gv ( t4, x4 ) * w4;

  xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;

  return xstar;
}
//****************************************************************************80

double r8_normal_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01 returns a unit pseudonormal R8.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has 
//    mean 0 and standard deviation 1.
//
//    Because this routine uses the Box Muller method, it requires pairs
//    of uniform random values to generate a pair of normal random values.
//    This means that on every other call, the code can use the second
//    value that it calculated.
//
//    However, if the user has changed the SEED value between calls,
//    the routine automatically resets itself and discards the saved data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_NORMAL_01, a normally distributed random value.
//
{
# define R8_PI 3.141592653589793

  double r1;
  double r2;
  static int seed1 = 0;
  static int seed2 = 0;
  static int seed3 = 0;
  static int used = 0;
  double v1;
  static double v2 = 0.0;
//
//  If USED is odd, but the input SEED does not match
//  the output SEED on the previous call, then the user has changed
//  the seed.  Wipe out internal memory.
//
  if ( ( used % 2 ) == 1 )
  {
    if ( *seed != seed2 )
    {
      used = 0;
      seed1 = 0;
      seed2 = 0;
      seed3 = 0;
      v2 = 0.0;
    }
  }
//
//  If USED is even, generate two uniforms, create two normals,
//  return the first normal and its corresponding seed.
//
  if ( ( used % 2 ) == 0 )
  {
    seed1 = *seed;

    r1 = r8_uniform_01 ( seed );

    if ( r1 == 0.0 )
    {
      cerr << "\n";
      cerr << "R8_NORMAL_01 - Fatal error!\n";
      cerr << "  R8_UNIFORM_01 returned a value of 0.\n";
      exit ( 1 );
    }

    seed2 = *seed;
    r2 = r8_uniform_01 ( seed );
    seed3 = *seed;
    *seed = seed2;

    v1 = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * R8_PI * r2 );
    v2 = sqrt ( - 2.0 * log ( r1 ) ) * sin ( 2.0 * R8_PI * r2 );
  }
//
//  If USED is odd (and the input SEED matched the output value from
//  the previous call), return the second normal and its corresponding seed.
//
  else
  {
    v1 = v2;
    *seed = seed3;
  }

  used = used + 1;

  return v1;
# undef R8_PI
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

