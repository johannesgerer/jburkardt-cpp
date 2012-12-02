# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "test_int_2d.hpp"

//****************************************************************************80

void legendre_dr_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_DR_COMPUTE: Gauss-Legendre quadrature by Davis-Rabinowitz method.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2007
//
//  Author:
//
//    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be greater than 0.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pi = 3.141592653589793;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "LEGENDRE_DR_COMPUTE - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    exit ( 1 );
  }

  e1 = ( double ) ( n * ( n + 1 ) );

  m = ( n + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * n + 2 );

    x0 = cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( n ) ) 
      / ( double ) ( 8 * n * n ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= n; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( n ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
//
//  Initial approximation H:
//
    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );
//
//  Refine H using one step of Newton's method:
//
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    x[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    w[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx );
  }

  if ( ( n % 2 ) == 1 )
  {
    x[0] = 0.0;
  }
//
//  Shift the data up.
//
  nmove = ( n + 1 ) / 2;
  ncopy = n - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = n + 1 - i;
    x[iback-1] = x[iback-ncopy-1];
    w[iback-1] = w[iback-ncopy-1];
  }
//
//  Reflect values for the negative abscissas.
//
  for ( i = 1; i <= n - nmove; i++ )
  {
    x[i-1] = - x[n-i];
    w[i-1] = w[n-i];
  }

  return;
}
//****************************************************************************80

double p00_exact ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_EXACT returns the exact integral for any problem.
//
//  Discussion:
//
//    This routine provides a "generic" interface to the exact integral
//    routines for the various problems, and allows a problem to be called
//    by index rather than by name.
//
//    In some cases, the "exact" value of the integral is in fact
//    merely a respectable approximation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Output, double P00_EXACT, the exact value of the integral.
//
{
  double exact;

  if ( problem == 1 )
  {
    exact = p01_exact ( );
  }
  else if ( problem == 2 )
  {
    exact = p02_exact ( );
  }
  else if ( problem == 3 )
  {
    exact = p03_exact ( );
  }
  else if ( problem == 4 )
  {
    exact = p04_exact ( );
  }
  else if ( problem == 5 )
  {
    exact = p05_exact ( );
  }
  else if ( problem == 6 )
  {
    exact = p06_exact ( );
  }
  else if ( problem == 7 )
  {
    exact = p07_exact ( );
  }
  else if ( problem == 8 )
  {
    exact = p08_exact ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_EXACT - Fatal error!\n";
    cerr << "  Illegal problem index = " << problem << "\n";
    exit ( 1 );
  }
  return exact;
}
//****************************************************************************80

void p00_fun ( int problem, int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_FUN evaluates the integrand for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  if ( problem == 1 )
  {
    p01_fun ( n, x, fx );
  }
  else if ( problem == 2 )
  {
    p02_fun ( n, x, fx );
  }
  else if ( problem == 3 )
  {
    p03_fun ( n, x, fx );
  }
  else if ( problem == 4 )
  {
    p04_fun ( n, x, fx );
  }
  else if ( problem == 5 )
  {
    p05_fun ( n, x, fx );
  }
  else if ( problem == 6 )
  {
    p06_fun ( n, x, fx );
  }
  else if ( problem == 7 )
  {
    p07_fun ( n, x, fx );
  }
  else if ( problem == 8 )
  {
    p08_fun ( n, x, fx );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_FUN - Fatal error!\n";
    cerr << "  Illegal problem index = " << problem << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void p00_lim ( int problem, double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_LIM returns the integration limits for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  if ( problem == 1 )
  {
    p01_lim ( a, b );
  }
  else if ( problem == 2 )
  {
    p02_lim ( a, b );
  }
  else if ( problem == 3 )
  {
    p03_lim ( a, b );
  }
  else if ( problem == 4 )
  {
    p04_lim ( a, b );
  }
  else if ( problem == 5 )
  {
    p05_lim ( a, b );
  }
  else if ( problem == 6 )
  {
    p06_lim ( a, b );
  }
  else if ( problem == 7 )
  {
    p07_lim ( a, b );
  }
  else if ( problem == 8 )
  {
    p08_lim ( a, b );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_LIM - Fatal error!\n";
    cerr << "  Illegal problem index = " << problem << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int p00_problem_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROBLEM_NUM returns the number of test integration problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P00_PROBLEM_NUM, the number of problems.
//
{
  int problem_num;

  problem_num = 8;

  return problem_num;
}
//****************************************************************************80

double p01_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_EXACT returns the exact integral for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P01_EXACT, the value of the integral.
//
{
  double exact;
  double pi = 3.141592653589793;

  exact = pi * pi / 6.0;

  return exact;
}
//****************************************************************************80

void p01_fun ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_FUN evaluates the integrand for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993,
//    ISBN: 047193898X,
//    LC: QA299.3E93.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the value of the integrand at X.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    fx[j] = 1.0 / ( 1.0 - x[0+j*2] * x[1+j*2] );
  }
  return;
}
//****************************************************************************80

void p01_lim ( double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_LIM returns the integration limits for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  a[0] = 0.0;
  a[1] = 0.0;

  b[0] = 1.0;
  b[1] = 1.0;

  return;
}
//****************************************************************************80

double p02_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_EXACT returns the exact integral for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P02_EXACT, the value of the integral.
//
{
  double exact;
  double pi = 3.141592653589793;

  exact = 2.0 * pi * log ( 2.0 );

  return exact;
}
//****************************************************************************80

void p02_fun ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_FUN evaluates the integrand for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993,
//    ISBN: 047193898X,
//    LC: QA299.3E93.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the value of the integrand at X.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    fx[j] = 1.0 / sqrt ( 1.0 - pow ( x[0+j*2] * x[1+j*2], 2 ) );
  }
  return;
}
//****************************************************************************80

void p02_lim ( double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_LIM returns the integration limits for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  a[0] = -1.0;
  a[1] = -1.0;

  b[0] = 1.0;
  b[1] = 1.0;

  return;
}
//****************************************************************************80

double p03_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_EXACT returns the exact integral for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P03_EXACT, the value of the integral.
//
{
  double exact;

  exact = ( 16.0 / 3.0 ) * ( 2.0 - sqrt ( 2.0 ) );

  return exact;
}
//****************************************************************************80

void p03_fun ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_FUN evaluates the integrand for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993,
//    ISBN: 047193898X,
//    LC: QA299.3E93.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the value of the integrand at X.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    fx[j] = 1.0 / sqrt ( 2.0 - x[0+j*2] - x[1+j*2] );
  }
  return;
}
//****************************************************************************80

void p03_lim ( double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_LIM returns the integration limits for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  a[0] = -1.0;
  a[1] = -1.0;

  b[0] = 1.0;
  b[1] = 1.0;

  return;
}
//****************************************************************************80

double p04_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_EXACT returns the exact integral for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P04_EXACT, the value of the integral.
//
{
  double exact;

  exact = ( sqrt ( 32.0 ) / 3.0 ) * ( sqrt ( 27.0 ) - sqrt ( 8.0 ) - 1.0 );

  return exact;
}
//****************************************************************************80

void p04_fun ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_FUN evaluates the integrand for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993,
//    ISBN: 047193898X,
//    LC: QA299.3E93.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the value of the integrand at X.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    fx[j] = 1.0 / sqrt ( 3.0 - x[0+j*2] - 2.0 * x[1+j*2] );
  }
  return;
}
//****************************************************************************80

void p04_lim ( double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_LIM returns the integration limits for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  a[0] = -1.0;
  a[1] = -1.0;

  b[0] = 1.0;
  b[1] = 1.0;

  return;
}
//****************************************************************************80

double p05_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_EXACT returns the exact integral for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P05_EXACT, the value of the integral.
//
{
  double exact;

  exact = 4.0 / 9.0;

  return exact;
}
//****************************************************************************80

void p05_fun ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_FUN evaluates the integrand for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993,
//    ISBN: 047193898X,
//    LC: QA299.3E93.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the value of the integrand at X.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    fx[j] = sqrt ( x[0+j*2] * x[1+j*2] );
  }
  return;
}
//****************************************************************************80

void p05_lim ( double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_LIM returns the integration limits for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  a[0] = 0.0;
  a[1] = 0.0;

  b[0] = 1.0;
  b[1] = 1.0;

  return;
}
//****************************************************************************80

double p06_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_EXACT returns the exact integral for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P06_EXACT, the value of the integral.
//
{
  double exact;
  double pi = 3.141592653589793;

  exact = ( 5.0 / 3.0 ) + ( pi / 16.0 );

  return exact;
}
//****************************************************************************80

void p06_fun ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_FUN evaluates the integrand for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993,
//    ISBN: 047193898X,
//    LC: QA299.3E93.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the value of the integrand at X.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    fx[j] = r8_abs ( pow ( x[0+j*2], 2 ) + pow ( x[1+j*2], 2 ) - 0.25 );
  }
  return;
}
//****************************************************************************80

void p06_lim ( double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_LIM returns the integration limits for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  a[0] = -1.0;
  a[1] = -1.0;

  b[0] = 1.0;
  b[1] = 1.0;

  return;
}
//****************************************************************************80

double p07_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_EXACT returns the exact integral for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P07_EXACT, the value of the integral.
//
{
  double exact;

  exact = 8.0 / 15.0;

  return exact;
}
//****************************************************************************80

void p07_fun ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_FUN evaluates the integrand for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993,
//    ISBN: 047193898X,
//    LC: QA299.3E93.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the value of the integrand at X.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    fx[j] = sqrt ( r8_abs ( x[0+j*2] - x[1+j*2] ) );
  }
  return;
}
//****************************************************************************80

void p07_lim ( double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_LIM returns the integration limits for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  a[0] = 0.0;
  a[1] = 0.0;

  b[0] = 1.0;
  b[1] = 1.0;

  return;
}
//****************************************************************************80

double p08_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_EXACT returns the exact integral for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P08_EXACT, the value of the integral.
//
{
  double exact;
  static double pi = 3.141592653589793;

  exact = 0.25 * pi * pow ( r8_erf ( 1.0 ) + r8_erf ( 4.0 ), 2 );

  return exact;
}
//****************************************************************************80

void p08_fun ( int n, double x[], double fx[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_FUN evaluates the integrand for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[2*N], the evaluation points.
//
//    Output, double FX[N], the value of the integrand at X.
//
{
  double arg;
  int j;

  for ( j = 0; j < n; j++ )
  {
    arg = pow ( x[0+j*2] - 4.0, 2 ) + pow ( x[1+j*2] - 1.0, 2 );
    fx[j] = exp ( - arg );
  }
  return;
}
//****************************************************************************80

void p08_lim ( double a[2], double b[2] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_LIM returns the integration limits for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double A[2], B[2], the lower and upper limits of integration.
//
{
  a[0] = 0.0;
  a[1] = 0.0;

  b[0] = 5.0;
  b[1] = 5.0;

  return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

 double r8_csevl ( double x, double a[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CSEVL evaluates a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double CS[N], the Chebyshev coefficients.
//
//    Input, int N, the number of Chebyshev coefficients.
//
//    Output, double R8_CSEVL, the Chebyshev series evaluated at X.
//
{
  double b0;
  double b1;
  double b2;
  int i;
  double twox;
  double value;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  Number of terms <= 0.\n";
    exit ( 1 );
  }

  if ( 1000 < n )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  Number of terms greater than 1000.\n";
    exit ( 1 );
 }

  if ( x < -1.1 || 1.1 < x )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  X outside (-1,+1).\n";
    exit ( 1 );
  }

  twox = 2.0 * x;
  b1 = 0.0;
  b0 = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[i];
  }

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double r8_erf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERF evaluates the error function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ERF, the error function of X.
//
{
  static double erfcs[21] = {
    -0.49046121234691808039984544033376E-01,
    -0.14226120510371364237824741899631,
    +0.10035582187599795575754676712933E-01,
    -0.57687646997674847650827025509167E-03,
    +0.27419931252196061034422160791471E-04,
    -0.11043175507344507604135381295905E-05,
    +0.38488755420345036949961311498174E-07,
    -0.11808582533875466969631751801581E-08,
    +0.32334215826050909646402930953354E-10,
    -0.79910159470045487581607374708595E-12,
    +0.17990725113961455611967245486634E-13,
    -0.37186354878186926382316828209493E-15,
    +0.71035990037142529711689908394666E-17,
    -0.12612455119155225832495424853333E-18,
    +0.20916406941769294369170500266666E-20,
    -0.32539731029314072982364160000000E-22,
    +0.47668672097976748332373333333333E-24,
    -0.65980120782851343155199999999999E-26,
    +0.86550114699637626197333333333333E-28,
    -0.10788925177498064213333333333333E-29,
    +0.12811883993017002666666666666666E-31 };
  static int nterf = 0;
  static double sqeps = 0.0;
  static double sqrtpi = 1.77245385090551602729816748334115;
  double value;
  static double xbig = 0.0;
  double y;

  if ( nterf == 0 )
  {
    nterf = r8_inits ( erfcs, 21, 0.1 * r8_mach ( 3 ) );
    xbig = sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) );
    sqeps = sqrt ( 2.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= sqeps )
  {
    value = 2.0 * x / sqrtpi;
  }
  else if ( y <= 1.0 )
  {
    value = x * ( 1.0 + r8_csevl ( 2.0 * x * x - 1.0, erfcs, nterf ) );
  }
  else if ( y <= xbig )
  {
    value = 1.0 - r8_erfc ( y );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  else
  {
    value = 1.0;
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

double r8_erfc ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERFC evaluates the co-error function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ERFC, the co-error function of X.
//
{
  static double erc2cs[49] = {
    -0.6960134660230950112739150826197E-01,
    -0.4110133936262089348982212084666E-01,
    +0.3914495866689626881561143705244E-02,
    -0.4906395650548979161280935450774E-03,
    +0.7157479001377036380760894141825E-04,
    -0.1153071634131232833808232847912E-04,
    +0.1994670590201997635052314867709E-05,
    -0.3642666471599222873936118430711E-06,
    +0.6944372610005012589931277214633E-07,
    -0.1371220902104366019534605141210E-07,
    +0.2788389661007137131963860348087E-08,
    -0.5814164724331161551864791050316E-09,
    +0.1238920491752753181180168817950E-09,
    -0.2690639145306743432390424937889E-10,
    +0.5942614350847910982444709683840E-11,
    -0.1332386735758119579287754420570E-11,
    +0.3028046806177132017173697243304E-12,
    -0.6966648814941032588795867588954E-13,
    +0.1620854541053922969812893227628E-13,
    -0.3809934465250491999876913057729E-14,
    +0.9040487815978831149368971012975E-15,
    -0.2164006195089607347809812047003E-15,
    +0.5222102233995854984607980244172E-16,
    -0.1269729602364555336372415527780E-16,
    +0.3109145504276197583836227412951E-17,
    -0.7663762920320385524009566714811E-18,
    +0.1900819251362745202536929733290E-18,
    -0.4742207279069039545225655999965E-19,
    +0.1189649200076528382880683078451E-19,
    -0.3000035590325780256845271313066E-20,
    +0.7602993453043246173019385277098E-21,
    -0.1935909447606872881569811049130E-21,
    +0.4951399124773337881000042386773E-22,
    -0.1271807481336371879608621989888E-22,
    +0.3280049600469513043315841652053E-23,
    -0.8492320176822896568924792422399E-24,
    +0.2206917892807560223519879987199E-24,
    -0.5755617245696528498312819507199E-25,
    +0.1506191533639234250354144051199E-25,
    -0.3954502959018796953104285695999E-26,
    +0.1041529704151500979984645051733E-26,
    -0.2751487795278765079450178901333E-27,
    +0.7290058205497557408997703680000E-28,
    -0.1936939645915947804077501098666E-28,
    +0.5160357112051487298370054826666E-29,
    -0.1378419322193094099389644800000E-29,
    +0.3691326793107069042251093333333E-30,
    -0.9909389590624365420653226666666E-31,
    +0.2666491705195388413323946666666E-31 };
  static double erfccs[59] = {
    +0.715179310202924774503697709496E-01,
    -0.265324343376067157558893386681E-01,
    +0.171115397792085588332699194606E-02,
    -0.163751663458517884163746404749E-03,
    +0.198712935005520364995974806758E-04,
    -0.284371241276655508750175183152E-05,
    +0.460616130896313036969379968464E-06,
    -0.822775302587920842057766536366E-07,
    +0.159214187277090112989358340826E-07,
    -0.329507136225284321486631665072E-08,
    +0.722343976040055546581261153890E-09,
    -0.166485581339872959344695966886E-09,
    +0.401039258823766482077671768814E-10,
    -0.100481621442573113272170176283E-10,
    +0.260827591330033380859341009439E-11,
    -0.699111056040402486557697812476E-12,
    +0.192949233326170708624205749803E-12,
    -0.547013118875433106490125085271E-13,
    +0.158966330976269744839084032762E-13,
    -0.472689398019755483920369584290E-14,
    +0.143587337678498478672873997840E-14,
    -0.444951056181735839417250062829E-15,
    +0.140481088476823343737305537466E-15,
    -0.451381838776421089625963281623E-16,
    +0.147452154104513307787018713262E-16,
    -0.489262140694577615436841552532E-17,
    +0.164761214141064673895301522827E-17,
    -0.562681717632940809299928521323E-18,
    +0.194744338223207851429197867821E-18,
    -0.682630564294842072956664144723E-19,
    +0.242198888729864924018301125438E-19,
    -0.869341413350307042563800861857E-20,
    +0.315518034622808557122363401262E-20,
    -0.115737232404960874261239486742E-20,
    +0.428894716160565394623737097442E-21,
    -0.160503074205761685005737770964E-21,
    +0.606329875745380264495069923027E-22,
    -0.231140425169795849098840801367E-22,
    +0.888877854066188552554702955697E-23,
    -0.344726057665137652230718495566E-23,
    +0.134786546020696506827582774181E-23,
    -0.531179407112502173645873201807E-24,
    +0.210934105861978316828954734537E-24,
    -0.843836558792378911598133256738E-25,
    +0.339998252494520890627359576337E-25,
    -0.137945238807324209002238377110E-25,
    +0.563449031183325261513392634811E-26,
    -0.231649043447706544823427752700E-26,
    +0.958446284460181015263158381226E-27,
    -0.399072288033010972624224850193E-27,
    +0.167212922594447736017228709669E-27,
    -0.704599152276601385638803782587E-28,
    +0.297976840286420635412357989444E-28,
    -0.126252246646061929722422632994E-28,
    +0.539543870454248793985299653154E-29,
    -0.238099288253145918675346190062E-29,
    +0.109905283010276157359726683750E-29,
    -0.486771374164496572732518677435E-30,
    +0.152587726411035756763200828211E-30 };
  static double erfcs[21] = {
    -0.49046121234691808039984544033376E-01,
    -0.14226120510371364237824741899631,
    +0.10035582187599795575754676712933E-01,
    -0.57687646997674847650827025509167E-03,
    +0.27419931252196061034422160791471E-04,
    -0.11043175507344507604135381295905E-05,
    +0.38488755420345036949961311498174E-07,
    -0.11808582533875466969631751801581E-08,
    +0.32334215826050909646402930953354E-10,
    -0.79910159470045487581607374708595E-12,
    +0.17990725113961455611967245486634E-13,
    -0.37186354878186926382316828209493E-15,
    +0.71035990037142529711689908394666E-17,
    -0.12612455119155225832495424853333E-18,
    +0.20916406941769294369170500266666E-20,
    -0.32539731029314072982364160000000E-22,
    +0.47668672097976748332373333333333E-24,
    -0.65980120782851343155199999999999E-26,
    +0.86550114699637626197333333333333E-28,
    -0.10788925177498064213333333333333E-29,
    +0.12811883993017002666666666666666E-31 };
  double eta;
  static int nterc2 = 0;
  static int nterf = 0;
  static int nterfc = 0;
  static double sqeps = 0.0;
  static double sqrtpi = 1.77245385090551602729816748334115;
  double value;
  static double xmax = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nterf == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nterf = r8_inits ( erfcs, 21, eta );
    nterfc = r8_inits ( erfccs, 59, eta );
    nterc2 = r8_inits ( erc2cs, 49, eta );

    xsml = - sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) );
    xmax = sqrt (- log ( sqrtpi * r8_mach ( 1 ) ) );
    xmax = xmax - 0.5 * log ( xmax ) / xmax - 0.01;
    sqeps = sqrt ( 2.0 * r8_mach ( 3 ) );
  }

  if ( x <= xsml )
  {
    value = 2.0;
    return value;
  }

  if ( xmax < x )
  {
    cerr << "\n";
    cerr << "R8_ERFC - Warning!\n";
    cerr << "  X so big that ERFC underflows.\n";
    value = 0.0;
    return value;
  }

  y = r8_abs ( x );

  if ( y < sqeps )
  {
    value = 1.0 - 2.0 * x / sqrtpi;
    return value;
  }
  else if ( y <= 1.0 )
  {
    value = 1.0 - x * ( 1.0 
      + r8_csevl ( 2.0 * x * x - 1.0, erfcs, nterf ) );
    return value;
  }

  y = y * y;

  if ( y <= 4.0 )
  {
    value = exp ( - y ) / r8_abs ( x ) * ( 0.5 
      + r8_csevl ( ( 8.0 / y - 5.0 ) / 3.0, erc2cs, nterc2 ) );
  }
  else 
  {
    value = exp ( - y ) / r8_abs ( x ) * ( 0.5 
      + r8_csevl ( 8.0 / y - 1.0, erfccs, nterfc ) );
  }

  if ( x < 0.0 )
  {
    value = 2.0 - value;
  }

  return value;
}
//****************************************************************************80

int r8_inits ( double dos[], int nos, double eta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_INITS initializes a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double DOS[NOS], the Chebyshev coefficients.
//
//    Input, int NOS, the number of coefficients.
//
//    Input, double ETA, the desired accuracy.
//
//    Output, int R8_INITS, the number of terms of the series needed
//    to ensure the requested accuracy.
//
{
  double err;
  int i;
  int value;

  if ( nos < 1 )
  {
    cerr << "\n";
    cerr << "R8_INITS - Fatal error!\n";
    cerr << "  Number of coefficients < 1.\n";
    exit ( 1 );
  }

  err = 0.0;

  for ( i = nos - 1; 0 <= i; i-- )
  {
    err = err + r8_abs ( dos[i] );
    if ( eta < err )
    {
      value = i + 1;
      return value;
    }
  }

  value = i;
  cerr << "\n";
  cerr << "R8_INITS - Warning!\n";
  cerr << "  ETA may be too small.\n";

  return value;
}
//****************************************************************************80

double r8_mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MACH returns double precision real machine constants.
//
//  Discussion:
//
//    Assuming that the internal representation of a double precision real
//    number is in base B, with T the number of base-B digits in the mantissa,
//    and EMIN the smallest possible exponent and EMAX the largest possible 
//    exponent, then
//
//      R8_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
//      R8_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
//      R8_MACH(3) = B^(-T), the smallest relative spacing.
//      R8_MACH(4) = B^(1-T), the largest relative spacing.
//      R8_MACH(5) = log10(B).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, Norman Schryer,
//    Algorithm 528:
//    Framework for a Portable Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, Number 2, June 1978, page 176-188.
//
//  Parameters:
//
//    Input, int I, chooses the parameter to be returned.
//    1 <= I <= 5.
//
//    Output, double R8_MACH, the value of the chosen parameter.
//
{
  double value;

  if ( i == 1 )
  {
    value = 4.450147717014403E-308;
  }
  else if ( i == 2 )
  {
    value = 8.988465674311579E+307;
  }
  else if ( i == 3 )
  {
    value = 1.110223024625157E-016;
  }
  else if ( i == 4 )
  {
    value = 2.220446049250313E-016;
  }
  else if ( i == 5 )
  {
    value = 0.301029995663981E+000;
  }
  else if ( 5 < i )
  {
    cerr << "\n";
    cerr << "R8_MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 5.\n";
    cerr << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

double *r8mat_uniform_01 ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2005
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
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }

  return value;
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
//    24 September 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
