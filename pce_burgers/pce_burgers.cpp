# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

int main ( );
double he_double_product_integral ( int i, int j );
double he_triple_product_integral ( int i, int j, int k );
double r8_factorial ( int n );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PCE_BURGERS.
//
//  Discussion:
//
//    The time-dependent viscous Burgers equation to be solved is:
//
//      du/dt = - d ( u*(1/2-u)) /dx + nu d2u/dx2
//
//    with boundary conditions
//
//      u(-3.0) = 0.0, u(+3.0) = 1.0.
//
//    The viscosity nu is assumed to be an uncertain quantity with
//    normal distribution of known mean and variance.
//
//    A polynomial chaos expansion is to be used, with Hermite polynomial
//    basis functions h(i,x), 0 <= i <= n.
//
//    Because the first two Hermite polynomials are simply 1 and x, 
//    we have that 
//
//      nu = nu_mean * h(0,x) + nu_variance * h(1,x).
//
//    We replace the time derivative by an explicit Euler approximation,
//    so that the equation now describes the value of U(x,t+dt) in terms
//    of known data at time t.
//
//    Now assume that the solution U(x,t) can be approximated
//    by the truncated expansion:
//
//      U(x,t) = sum ( 0 <= i <= n ) c(i,t) * h(i,x)
//
//    In the equation, we replace U by its expansion, and then multiply
//    successively by each of the basis functions h(*,x) to get a set of
//    n+1 equations that can be used to determine the values of c(i,t+dt).
//
//    This process is repeated until the desired final time is reached.
//
//    At any time, the coefficients c(0,t) contain information definining
//    the expected value of u(x,t) at that time, while the higher order coefficients
//    can be used to deterimine higher moments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2012
//
//  Author:
//
//    The original FORTRAN90 version of this program was written by Gianluca Iaccarino.
//    This C++ version is by John Burkardt.
//
//  Local parameters:
//
//    Local, double DT, the timestep.
//
//    Local, double DX, the spacing between grid points.
//
//    Local, int N, the number of intervals in the spatial domain.
//
//    Local, double NUMEAN, the mean of viscosity.
//
//    Local, double NUVARIANCE, the variance of viscosity.
//
//    Local, int P, the order of the PC expansion.
//
//    Local, double T, the current time.
//
//    Local, double TF, the final integration time.
//
//    Local, double U1[(N+1)*(P+1)], the PCE representation at the current time.
//
//    Local, double U2[(N+1)*(P+1)], the PCE representation for the next time.
//
//    Local, double X[N+1], the grid points.
//
{
  double conv;
  double dp;
  double dt;
  double dx;
  int i;
  int it;
  int ix;
  int j;
  int k;
  int n;
  int nt;
  double numean;
  double nuvariance;
  ofstream output;
  string output_filename;
  int p;
  double t1;
  double t2;
  double term1;
  double term2;
  double tf;
  double ti;
  double tp;
  double *u1;
  double *u2;
  double *umean;
  double *uvariance;
  double visc[2];
  double *x;

  p = 5 ;
  n = 32;
  nt = 2000;
  ti = 0.0;
  tf = 2.0;
  dt = ( tf - ti ) / ( double ) ( nt );
  numean = 0.25;
  nuvariance = 0.08;

  timestamp ( );
  cout << "\n";
  cout << "PCE_BURGERS:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Polynomial Chaos Solution\n";
  cout << "  1D Burgers equation\n";
  cout << "  Original version by Gianluca Iaccarino\n";
  cout << "\n";
  cout << "  PCE order       = " << p << "\n";
  cout << "  Number of cells = " << n << "\n";
  cout << "  Time step       = " << dt << "\n";
  cout << "  Initial time    = " << ti << "\n";
  cout << "  Final time      = " << tf << "\n";
  cout << "  Viscosity Mean  = " << numean << "\n";
  cout << "  Viscosity Var   = " << nuvariance << "\n";
  cout << "\n";

  u1 = new double[(n+1)*(p+1)];
  u2 = new double[(n+1)*(p+1)];
  x = new double[n+1];
//
//  Define some numerical parameters
//
  dx = 6.0 / ( double ) ( n );
  conv = dt / ( 2.0 * dx );
//
//  The expansion for viscosity stops at the linear term.
//
  visc[0] = numean * dt / ( dx * dx );
  visc[1] = nuvariance * dt / ( dx * dx );
//
//  Define a uniform grid
//
  for ( i = 0; i <= n; i++ )
  {
    x[i] = ( ( double ) ( n - i ) * ( -3.0 )   
           + ( double ) (     i ) * ( +3.0 ) ) 
           / ( double ) ( n     );
  }
//
//  Set the initial conditions
//
  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= n; i++ )
    {
      u1[i+j*(n+1)] = 0.0;
    }
  }

  for ( i = 0; i <= n; i++ )
  {
    u1[i+0*(n+1)] = 0.5 + x[i] / 6.0;
  }
//
//  Write the current solution.
//
  output_filename = "burgers.history.txt";

  output.open ( output_filename.c_str ( ) );

  output << "----------\n";
  output << "T = " << t1 << "\n";
  for ( i = 0; i <= n; i++ )
  {
    output << "  " << setw(14) << x[i];
    for ( j = 0; j <= p; j++ )
    {
      output << "  " << setw(14) << u1[i+j*(n+1)];
    }
    output << "\n";
  }
//
//  Time integration
//
  t1 = ti;

  for ( it = 1; it <= nt; it++ )
  {
    t2 = ( ( double ) ( nt - it ) * ti   
         + ( double ) (      it ) * tf ) 
         / ( double ) ( nt      );
//
//  Boundary conditions.
//
    for ( j = 0; j <= p; j++ )
    {
      u2[0+j*(n+1)] = 0.0;
    }
    u2[n+0*(n+1)] = 1.0;
    for ( j = 1; j <= p; j++ )
    {
      u2[n+j*(n+1)] = 0.0;
    }

    for ( k = 0; k <= p; k++ )
    {
      dp = he_double_product_integral ( k, k );

      for ( ix = 1; ix < n; ix++ )
      {
//
//  Viscous term.
//
        term1 = visc[0] * ( u1[ix+1+k*(n+1)] - 2.0 * u1[ix+k*(n+1)] + u1[ix-1+k*(n+1)] );
        i = 1;
        for ( j = 0; j <= p; j++ )
        {
          tp = he_triple_product_integral ( i, j, k );
          term1 = term1 + visc[i] 
            * ( u1[ix+1+j*(n+1)] - 2.0 * u1[ix+j*(n+1)] + u1[ix-1+j*(n+1)] ) * tp / dp;
        }
//
//  Convective term.
//
        term2 = - conv * 0.5 * ( u1[ix+1+k*(n+1)] - u1[ix-1+k*(n+1)] );
        for ( j = 0; j <= p; j++ )
        {
          for ( i = 0; i <= p; i++ )
          {
            tp = he_triple_product_integral ( i, j, k );
            term2 = term2 + ( conv 
              * u1[ix+i*(n+1)] * ( u1[ix+1+j*(n+1)] - u1[ix-1+j*(n+1)] ) * tp ) / dp;
          }
        }
        u2[ix+k*(n+1)] = u1[ix+k*(n+1)] + term1 + term2;
      }
    }

    t1 = t2;
    for ( j = 0; j <= p; j++ )
    {
      for ( i = 0; i <= n; i++ )
      {
        u1[i+j*(n+1)] = u2[i+j*(n+1)];
      }
    }
//
//  Print solution every 100 time steps.
//
    if ( ( it % 100 ) == 0 )
    {
      output << "----------\n";
      output << "T = " << t1 << "\n";
      for ( i = 0; i <= n; i++ )
      {
        output << "  " << setw(14) << x[i];
        for ( j = 0; j <= p; j++ )
        {
          output << "  " << setw(14) << u1[i+j*(n+1)];
        }
        output << "\n";
      }
    }
  }
  output.close ( );
  cout << "  Time history in \"" << output_filename << "\".\n";
//
//  Compute the mean and variance.
//
  umean = new double[n+1];
  uvariance = new double[n+1];

  for ( i = 0; i <= n; i++ )
  {
    umean[i] = u1[i+0*(n+1)];
  }

  for ( i = 0; i <= n; i++ )
  {
    uvariance[i] = 0.0;
    for ( j = 1; j <= p; j++ )
    {
      dp = he_double_product_integral ( j, j );
      uvariance[i] = uvariance[i] + pow ( u1[i+j*(n+1)], 2 ) * dp;
    }
  }
//
//  Write data about the solution at the final time.
//
  output_filename = "burgers.moments.txt";
  output.open ( output_filename.c_str ( ) );
  output << "X E[U] Var[U]\n";
  for ( i = 0; i <= n; i++ )
  {
    output << "  " << setw(18) << x[i]
           << "  " << setw(18) << umean[i]
           << "  " << setw(18) << uvariance[i] << "\n";
  }
  output.close ( );
  cout << "  Moments written to \"" << output_filename << "\".\n";

  output_filename = "burgers.modes.txt";
  output.open ( output_filename.c_str ( ) );
  output << "X U_0 ... U_P\n";
  for ( i = 0; i <= n; i++ )
  {
    output << "  " << setw(20) << x[i];
    for ( j = 0; j <= p; j++ )
    {
      output << "  " << setw(20) << u1[i+j*(n+1)];
    }
    output << "\n";
  }
  output.close ( );
  cout << "  Final modes written to \"" << output_filename << "\".\n";
//
//  Free memory.
//
  delete [] u1;
  delete [] u2;
  delete [] umean;
  delete [] uvariance;
  delete [] x;
//
//  Terminate.
//
  cout << "\n";
  cout << "PCE_BURGERS:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double he_double_product_integral ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    HE_DOUBLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*e^(-x^2/2).
//
//  Discussion:
//
//    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
//    Princeton, 2010,
//    ISBN13: 978-0-691-14212-8,
//    LC: QA274.23.X58.
//
//  Parameters:
//
//    Input, int I, J, the polynomial indices.
//
//    Output, double HE_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  double value;

  if ( i == j )
  {
    value = r8_factorial ( i );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double he_triple_product_integral ( int i, int j, int k )

//****************************************************************************80
//
//  Purpose:
//
//    HE_TRIPLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
//
//  Discussion:
//
//    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
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
//  Reference:
//
//    Dongbin Xiu,
//    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
//    Princeton, 2010,
//    ISBN13: 978-0-691-14212-8,
//    LC: QA274.23.X58.
//
//  Parameters:
//
//    Input, int I, J, K, the polynomial indices.
//
//    Output, double HE_TRIPLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  int s;
  double value;

  s = ( i + j + k ) / 2;

  if ( s < i || s < j || s < k )
  {
    value = 0.0;
  }
  else if ( ( ( i + j + k ) % 2 ) != 0 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_factorial ( i ) / r8_factorial ( s - i ) 
          * r8_factorial ( j ) / r8_factorial ( s - j ) 
          * r8_factorial ( k ) / r8_factorial ( s - k );
  }

  return value;
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
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
