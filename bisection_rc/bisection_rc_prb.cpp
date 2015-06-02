# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "bisection_rc.hpp"

int main ( );
void test01 ( );
double f01 ( double x );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BISECTION_RC_PRB.
//
//  Discussion:
//
//    BISECTION_RC_TEST tests the BISECTION_RC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BISECTION_RC_PRB:\n";
  cout << "  C++ version.\n";
  cout << "  Test the BISECTION_RC library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BISECTION_RC_PRB:\n";
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
//    TEST01 tests BISECTION_RC, evaluating the function in a separate routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double dx;
  double dx_tol;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Demonstrate BISECTION_RC on a simple example.\n";
  cout << "  The function is evaluated in a separate routine.\n";

  fx_tol = 1.0E-08;
  dx_tol = 1.0E-06;
  it = 0;
  it_max = 30;

  a = 0.0;
  b = 1.0;
  fx = 0.0;
  job = 0;

  cout << "\n";
  cout << "     I      X               FX              DX\n";
  cout << "\n";

  for ( ; ; )
  {
    x = bisection_rc ( a, b, fx, job );

    if ( job < 0 )
    {
      cout << "\n";
      cout << "  Error return.\n";
      break;
    }

    it = it + 1;

    fx = f01 ( x );

    if ( it <= 2 )
    {
      dx = fabs ( b - a );
    }
    else
    {
      dx = 0.5 * fabs ( b - a );
    }

    cout << "  " << setw(4) << it
         << "  " << setw(14) << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << dx << "\n";

    if ( fabs ( fx ) <= fx_tol )
    {
      cout << "\n";
      cout << "  Function is small.\n";
      break;
    }

    if ( dx <= dx_tol )
    {
      cout << "\n";
      cout << "  Interval is tiny.\n";
      break;
    }

    if ( it_max <= it )
    {
      cout << "\n";
      cout << "  Reached iteration limit.\n";
      break;
    }

  }

  cout << "\n";
  cout << "  A = " << setw(14) << a << " F(A) = " << f01 ( a ) << "\n";
  cout << "  X = " << setw(14) << x << " F(X) = " << f01 ( x ) << "\n";
  cout << "  B = " << setw(14) << b << " F(B) = " << f01 ( b ) << "\n";

  return;
}
//****************************************************************************80

double f01 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F01 evaluates the function f(x) = cos ( x ) - x which is zero around 0.74
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double F01, the function value.
//
{
  double value;

  value = cos ( x ) - x;

  return value;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BISECTION_RC, evaluating the function within the routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double dx;
  double dx_tol;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Demonstrate BISECTION_RC on a simple example.\n";
  cout << "  The function is evaluated within this routine.\n";

  fx_tol = 1.0E-09;
  dx_tol = 1.0E-09;
  it = 0;
  it_max = 30;

  a = 0.0;
  b = 1.0;
  fx = 0.0;
  job = 0;

  cout << "\n";
  cout << "     I      X               FX              DX\n";
  cout << "\n";

  for ( ; ; )
  {
    x = bisection_rc ( a, b, fx, job );

    if ( job < 0 )
    {
      cout << "\n";
      cout << "  Error return.\n";
      break;
    }

    it = it + 1;

    fx = cos ( 100.0 * x ) - 4.0 * erf ( 30.0 * x - 10.0 );

    if ( it <= 2 )
    {
      dx = fabs ( b - a );
    }
    else
    {
      dx = 0.5 * fabs ( b - a );
    }

    cout << "  " << setw(4) << it
         << "  " << setw(14) << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << dx << "\n";

    if ( fabs ( fx ) <= fx_tol )
    {
      cout << "\n";
      cout << "  Function is small.\n";
      break;
    }

    if ( dx <= dx_tol )
    {
      cout << "\n";
      cout << "  Interval is tiny.\n";
      break;
    }

    if ( it_max <= it )
    {
      cout << "\n";
      cout << "  Reached iteration limit.\n";
      break;
    }

  }

  cout << "\n";
  fx = cos ( 100.0 * a ) - 4.0 * erf ( 30.0 * a - 10.0 );
  cout << "  A = " << setw(14) << a << ", F(A) = " << setw(14) << fx << "\n";
  fx = cos ( 100.0 * x ) - 4.0 * erf ( 30.0 * x - 10.0 );
  cout << "  X = " << setw(14) << x << ", F(X) = " << setw(14) << fx << "\n";
  fx = cos ( 100.0 * b ) - 4.0 * erf ( 30.0 * b - 10.0 );
  cout << "  B = " << setw(14) << b << ", F(B) = " << setw(14) << fx << "\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests BISECTION_RC, to invert the cardiod CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double alpha = 0.0;
  double b;
  double beta = 0.25;
  double cdf;
  double dx;
  double dx_tol;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double r8_pi = 3.141592653589793;
  double x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Demonstrate BISECTION_RC on a probability example.\n";
  cout << "\n";
  cout << "  The cardioid probability density function has a\n";
  cout << "  cumulative density function of the form:\n";
  cout << "    CDF(X) = ( pi + x - alpha + 2 beta * sin ( x - alpha ) ) / ( 2 * pi )\n";
  cout << "  where alpha and beta are parameters, and x is a value\n";
  cout << "  in the range -pi <= x <= +pi.\n";
  cout << "\n";
  cout << "  CDF(X) is the probability that a random sample will have\n";
  cout << "  a value less than or equal to X.\n";
  cout << "\n";
  cout << "  As X moves from -pi to +pi,\n";
  cout << "  the CDF rises from 0 (no probability)\n";
  cout << "  to 1 (certain probability).\n";
  cout << "\n";
  cout << "  Assuming that:\n";
  cout << "  * ALPHA = " << alpha << "\n";
  cout << "  * BETA =  " << beta << "\n";
  cout << "  determine the value X where the Cardioid CDF is exactly 0.75.\n";

  fx_tol = 1.0E-05;
  dx_tol = 1.0E-08;
  it = 0;
  it_max = 30;

  job = 0;
  a = - r8_pi;
  b = + r8_pi;

  fx = 0.0;

  cout << "\n";
  cout << "     I      X               FX              DX\n";
  cout << "\n";

  for ( ; ; )
  {
    x = bisection_rc ( a, b, fx, job );

    if ( job < 0 )
    {
      cout << "\n";
      cout << "  Error return.\n";
      break;
    }

    it = it + 1;

    cdf = ( r8_pi + x - alpha + 2.0 * beta * sin ( x - alpha ) ) / ( 2.0 * r8_pi );
    fx = cdf - 0.75;

    if ( it <= 2 )
    {
      dx = fabs ( b - a );
    }
    else
    {
      dx = 0.5 * fabs ( b - a );
    }

    cout << "  " << setw(4) << it
         << "  " << setw(14) << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << dx << "\n";

    if ( fabs ( fx ) <= fx_tol )
    {
      cout << "\n";
      cout << "  Function is small.\n";
      break;
    }

    if ( dx <= dx_tol )
    {
      cout << "\n";
      cout << "  Interval is tiny.\n";
      break;
    }

    if ( it_max <= it )
    {
      cout << "\n";
      cout << "  Reached iteration limit.\n";
      break;
    }

  }

  cout << "\n";
  cdf = ( r8_pi + a - alpha + 2.0 * beta * sin ( a - alpha ) ) / ( 2.0 * r8_pi );
  fx = cdf - 0.75;
  cout << "  A = " << setw(14) << a << ", F(A) = " << setw(14) << fx << "\n";
  cdf = ( r8_pi + x - alpha + 2.0 * beta * sin ( x - alpha ) ) / ( 2.0 * r8_pi );
  fx = cdf - 0.75;
  cout << "  X = " << setw(14) << x << ", F(X) = " << setw(14) << fx << "\n";
  cdf = ( r8_pi + b - alpha + 2.0 * beta * sin ( b - alpha ) ) / ( 2.0 * r8_pi );
  fx = cdf - 0.75;
  cout << "  B = " << setw(14) << b << ", F(B) = " << setw(14) << fx << "\n";

  cout << "\n";
  cout << "  Look at the actual cardioid CDF value now:\n";
  cout << "\n";
  cdf = ( r8_pi + x - alpha + 2.0 * beta * sin ( x - alpha ) ) / ( 2.0 * r8_pi );
  cout << "  Cardioid(" << x << ") = " << cdf << "\n";

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests BISECTION_RC for the pipe freezing problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Cleve Moler,
//    Numerical Computing with MATLAB,
//    SIAM, 2004,
//    ISBN13: 978-0-898716-60-3,
//    LC: QA297.M625,
//    ebook: http://www.mathworks.com/moler/chapters.html
//
{
  double a;
  double alpha;
  double b;
  double dx;
  double dx_tol;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double t;
  double tc;
  double ti;
  double x;

  cout << "\n";
  cout << "BISECTION_RC_TEST04\n";
  cout << "  The freezing pipe problem.\n";
  cout << "\n";
  cout << "  At the beginning of a cold spell, the soil is at a uniform\n";
  cout << "  temperature of Ti.  The cold spell applies a uniform air\n";
  cout << "  temperature of Tc, which begins to cool the soil.\n";
  cout << "  As a function of depth x and time t, the soil temperature\n";
  cout << "  will now cool down as:\n";
  cout << "    ( T(x,t) - Tc ) / ( Ti - Tc ) = erf ( 0.5 * x / sqrt ( alpha * t ) ).\n";
  cout << "  where:\n";
  cout << "    Ti =  20 degrees centigrade,\n";
  cout << "    Tc = -15 degrees centigrade,\n";
  cout << "    alpha = 0.000000138 meter^2 / second, thermal conductivity;\n";
  cout << "    and erf() is the error function.\n";
  cout << "  Water freezes at 0 degrees centigrade.\n";
  cout << "\n";
  cout << "  What depth x in meters must a water pipe be buried so that it will\n";
  cout << "  not freeze even if this cold snap lasts for 60 days?\n";
//
//  Problem parameters.
//
  ti = 20.0;
  tc = -15.0;
  t = 60.0 * 24.0 * 60.0 * 60.0;
  alpha = 0.000000138;
//
//  Iteration parameters.
//
  fx_tol = 1.0E-09;
  dx_tol = 1.0E-09;
  it = 0;
  it_max = 30;
  job = 0;
  fx = 0.0;
//
//  Initial guess for interval.
//
  a = 0.0;
  b = 1000.0;

  cout << "\n";
  cout << "     I      X               FX              DX\n";
  cout << "\n";

  for ( ; ; )
  {
    x = bisection_rc ( a, b, fx, job );

    if ( job < 0 )
    {
      cout << "\n";
      cout << "  Error return.\n";
      break;
    }

    it = it + 1;

    fx = tc + ( ti - tc ) * erf ( 0.5 * x / sqrt ( alpha * t ) );

    if ( it <= 2 )
    {
      dx = fabs ( b - a );
    }
    else
    {
      dx = 0.5 * fabs ( b - a );
    }

    cout << "  " << setw(4) << it
         << "  " << setw(14) << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << dx << "\n";

    if ( fabs ( fx ) <= fx_tol )
    {
      cout << "\n";
      cout << "  Function is small.\n";
      break;
    }

    if ( dx <= dx_tol )
    {
      cout << "\n";
      cout << "  Interval is tiny.\n";
      break;
    }

    if ( it_max <= it )
    {
      cout << "\n";
      cout << "  Reached iteration limit.\n";
      break;
    }

  }

  cout << "\n";
  fx = tc + ( ti - tc ) * erf ( 0.5 * a / sqrt ( alpha * t ) );
  cout << "  A = " << a << ", F(A) = " << fx << "\n";
  fx = tc + ( ti - tc ) * erf ( 0.5 * x / sqrt ( alpha * t ) );
  cout << "  X = " << x << ", F(X) = " << fx << "\n";
  fx = tc + ( ti - tc ) * erf ( 0.5 * b / sqrt ( alpha * t ) );
  cout << "  B = " << b << ", F(B) = " << fx << "\n";

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests BISECTION_RC for Kepler's problem.
//
//  Discussion:
//
//    Kepler's equation has the form:
//
//      X = M + E * sin ( X )
//
//    X represents the eccentric anomaly of a planet, the angle between the
//    perihelion (the point on the orbit nearest to the sun) through the sun 
//    to the center of the ellipse, and the line from the center of the ellipse
//    to the planet.
//
//    There are two parameters, E and M:
//
//    * E is the eccentricity of the orbit, which should be between 0 and 1.0;
//
//    * M is the angle from the perihelion made by a fictitious planet traveling
//      on a circular orbit centered at the sun, and traveling at a constant
//      angular velocity equal to the average angular velocity of the true
//      planet.  M is usually between 0 and 180 degrees, but can have any value.
//
//    For convenience, X and M are measured in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Cleve Moler,
//    Numerical Computing with MATLAB,
//    SIAM, 2004,
//    ISBN13: 978-0-898716-60-3,
//    LC: QA297.M625,
//    ebook: http://www.mathworks.com/moler/chapters.html
//
{
  double ad;
  double ar;
  double bd;
  double br;
  double dx;
  double dx_tol;
  double e;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double md;
  double mr;
  const double r8_pi = 3.141592653589793;
  double xd;
  double xr;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  The Kepler equation.\n";
  cout << "\n";
  cout << "  Kepler's equation has the form\n";
  cout << "\n";
  cout << "    X = M + E * sin ( X )\n";
  cout << "\n";
  cout << "  X represents the eccentric anomaly of a planet, the angle between the\n";
  cout << "  perihelion (the point on the orbit nearest to the sun) through the sun\n";
  cout << "  to the center of the ellipse, and the line from the center of the ellipse\n";
  cout << "  to the planet.\n";
  cout << "\n";
  cout << "  There are two parameters, E and M:\n";
  cout << "\n";
  cout << "  * E is the eccentricity of the orbit, which should be between 0 and 1.0;\n";
  cout << "\n";
  cout << "  * M is the angle from the perihelion made by a fictitious planet traveling\n";
  cout << "    on a circular orbit centered at the sun, and traveling at a constant\n";
  cout << "    angular velocity equal to the average angular velocity of the true\n";
  cout << "    planet.  M is usually between 0 and 180 degrees, but can have any value.\n";
  cout << "\n";
  cout << "  For convenience, X and M are measured in degrees.\n";
//
//  Problem parameters.
//
  md = 24.851090;
  mr = md * r8_pi / 180.0;
  e = 0.1;

  cout << "\n";
  cout << "  Given eccentricity E = " << e << "\n";
  cout << "  Given angle M = " << md << " (degrees)\n";
  cout << "                = " << mr << " (radians)\n";
  cout << "\n";
  cout << "  Given E and M, find corresponding X.\n";
//
//  Iteration parameters.
//
  fx_tol = 1.0E-09;
  dx_tol = 1.0E-09;
  it = 0;;
  it_max = 30;
  job = 0;
  fx = 0.0;
//
//  Initial guess for interval.
//
  ad = 0.0;
  bd = 180.0;

  ar = ad * r8_pi / 180.0;
  br = bd * r8_pi / 180.0;

  cout << "\n";
  cout << "     I      X               FX              DX\n";
  cout << "\n";

  for ( ; ; )
  {
    xr = bisection_rc ( ar, br, fx, job );

    if ( job < 0 )
    {
      cout << "\n";
      cout << "  Error return.\n";
      break;
    }

    it = it + 1;

    fx = xr - mr - e * sin ( xr );

    if ( it <= 2 )
    {
      dx = fabs ( br - ar );
    }
    else
    {
      dx = 0.5 * fabs ( br - ar );
    }

    cout << "  " << setw(4) << it
         << "  " << setw(14) << xr
         << "  " << setw(14) << fx
         << "  " << setw(14) << dx << "\n";
 
    if ( fabs ( fx ) <= fx_tol )
    {
      cout << "\n";
      cout << "  Function is small.\n";
      break;
    }

    if ( dx <= dx_tol )
    {
      cout << "\n";
      cout << "  Interval is tiny.\n";
      break;
    }

    if ( it_max <= it )
    {
      cout << "\n";
      cout << "  Reached iteration limit.\n";
      break;
    }

  }

  cout << "\n";
  cout << "  In Radians:\n";
  cout << "\n";
  fx = ar - mr - e * sin ( ar );
  cout << "  A = " << ar << ", F(A) = " << fx << "\n";
  fx = xr - mr - e * sin ( xr );
  cout << "  X = " << xr << ", F(X) = " << fx << "\n";
  fx = br - mr - e * sin ( br );
  cout << "  B = " << br << ", F(B) = " << fx << "\n";

  ad = ar * 180.0 / r8_pi;
  xd = xr * 180.0 / r8_pi;
  bd = br * 180.0 / r8_pi;

  cout << "\n";
  cout << "  In Degrees:\n";
  cout << "\n";
  fx = ( ad - md ) * r8_pi / 180.0 - e * sin ( ad * r8_pi / 180.0 );
  cout << "  A = " << ad << ", F(A) = " << fx << "\n";
  fx = ( xd - md ) * r8_pi / 180.0 - e * sin ( xd * r8_pi / 180.0 );
  cout << "  X = " << xd << ", F(X) = " << fx << "\n";
  fx = ( bd - md ) * r8_pi / 180.0 - e * sin ( bd * r8_pi / 180.0);
  cout << "  B = " << bd << ", F(B) = " << fx << "\n";

  return;
}

