# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

# include "black_scholes.hpp"

//****************************************************************************80

double *asset_path ( double s0, double mu, double sigma, double t1, int n, 
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    ASSET_PATH simulates the behavior of an asset price over time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    Original MATLAB version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    Black-Scholes for Scientific Computing Students,
//    Computing in Science and Engineering,
//    November/December 2004, Volume 6, Number 6, pages 72-79.
//
//  Parameters:
//
//    Input, double S0, the asset price at time 0.
//
//    Input, double MU, the expected growth rate.
//
//    Input, double SIGMA, the volatility of the asset.
//
//    Input, double T1, the expiry date.
//
//    Input, integer N, the number of steps to take between 0 and T1.
//
//    Input/output, int SEED, a seed for the random number generator.
//
//    Output, double ASSET_PATH[N+1], the option values from time 0 to T1 
//    in equal steps.
//
{
  double dt;
  int i;
  double p;
  double *r;
  double *s;

  dt = t1 / ( double ) ( n );

  r = r8vec_normal_01_new ( n, seed );

  s = new double[n+1];

  s[0] = s0;
  p = s0;
  for ( i = 1; i <= n; i++ )
  {
    p = p * exp ( ( mu - sigma * sigma ) * dt + sigma * sqrt ( dt ) * r[i-1] );
    s[i] = p;
  }

  delete [] r;

  return s;
}
//****************************************************************************80

double binomial ( double s0, double e, double r, double sigma, double t1, 
  int m )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL uses the binomial method for a European call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    Original MATLAB version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    Black-Scholes for Scientific Computing Students,
//    Computing in Science and Engineering,
//    November/December 2004, Volume 6, Number 6, pages 72-79.
//
//  Parameters:
//
//    Input, double S0, the asset price at time 0.
//
//    Input, double E, the exercise price.
//
//    Input, double R, the interest rate.
//
//    Input, double SIGMA, the volatility of the asset.
//
//    Input, double T1, the expiry date.
//
//    Input, int M, the number of steps to take 
//    between 0 and T1.
//
//    Output, double BIONOMIAL, the option value at time 0.
//
{
  double a;
  double b;
  double c;
  double d;
  double dt;
  int i;
  int n;
  double p;
  double u;
  double *w;
//
//  Time stepsize.
//
  dt = t1 / ( double ) ( m );

  a = 0.5 * ( exp ( - r * dt ) + exp ( ( r + sigma * sigma ) * dt ) );

  d = a - sqrt ( a * a - 1.0 );
  u = a + sqrt ( a * a - 1.0 );

  p = ( exp ( r * dt ) - d ) / ( u - d );

  w = new double[m+1];

  for ( i = 0; i <= m; i++ )
  {
    w[i] = r8_max ( s0 * pow ( d, m - i ) * pow ( u, i ) - e, 0.0 );
  }
//
//  Trace backwards to get the option value at time 0.
//
  for ( n = m - 1; 0 <= n; n-- )
  {
    for ( i = 0; i <= n; i++ )
    {
      w[i] = ( 1.0 - p ) * w[i] + p * w[i+1];
    }
  }

  for ( i = 0; i < m + 1; i++ )
  {
    w[i] = exp ( - r * t1 ) * w[i];
  }
  c = w[0];

  delete [] w;

  return c;
}
//****************************************************************************80

double bsf ( double s0, double t0, double e, double r, double sigma, double t1 )

//****************************************************************************80
//
//  Purpose:
//
//    BSF evaluates the Black-Scholes formula for a European call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    Original MATLAB version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    Black-Scholes for Scientific Computing Students,
//    Computing in Science and Engineering,
//    November/December 2004, Volume 6, Number 6, pages 72-79.
//
//  Parameters:
//
//    Input, double S0, the asset price at time T0.
//
//    Input, double T0, the time at which the asset price is known.
//
//    Input, double E, the exercise price.
//
//    Input, double R, the interest rate.
//
//    Input, double SIGMA, the volatility of the asset.
//
//    Input, double T1, the expiry date.
//
//    Output, double BSF, the value of the call option.
//
{
  double c;
  double d1;
  double d2;
  double n1;
  double n2;
  double tau;

  tau = t1 - t0;

  if ( 0.0 < tau )
  {
    d1 = ( log ( s0 / e ) + ( r + 0.5 * sigma * sigma ) * tau ) 
      / ( sigma * sqrt ( tau ) );

    d2 = d1 - sigma * sqrt ( tau );

    n1 = 0.5 * ( 1.0 + erf ( d1 / sqrt ( 2.0 ) ) );
    n2 = 0.5 * ( 1.0 + erf ( d2 / sqrt ( 2.0 ) ) );

    c = s0 * n1 - e * exp ( - r * tau ) * n2;
  }
  else
  {
    c = r8_max ( s0 - e, 0.0 );
  }
  return c;
}
//****************************************************************************80

double *forward ( double e, double r, double sigma, double t1, int nx, 
  int nt, double smax )

//****************************************************************************80
//
//  Purpose:
//
//    FORWARD uses the forward difference method to value a European call option.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    Original MATLAB version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    Black-Scholes for Scientific Computing Students,
//    Computing in Science and Engineering,
//    November/December 2004, Volume 6, Number 6, pages 72-79.
//
//  Parameters:
//
//    Input, double E, the exercise price.
//
//    Input, double R, the interest rate.
//
//    Input, double SIGMA, the volatility of the asset.
//
//    Input, double T1, the expiry date.
//
//    Input, int NX, the number of "space" steps used to 
//    divide the interval [0,L].
//
//    Input, int NT, the number of time steps.
//
//    Input, double SMAX, the maximum value of S to consider.
//
//    Output, double U[(NX-1)*(NT+1)], the value of the European 
//    call option.
//
{
  double *a;
  double *b;
  double *c;
  double dt;
  double dx;
  int i;
  int j;
  double p;
  double t;
  double *u;
  double u0;

  dt = t1 / ( double ) ( nt );
  dx = smax / ( double ) ( nx );

  a = new double[nx-1];
  b = new double[nx-1];
  c = new double[nx-1];

  for ( i = 0; i < nx - 1; i++ )
  {
    b[i] = 1.0 - r * dt - dt * pow ( sigma * ( i + 1 ), 2 );
  }

  for ( i = 0; i < nx - 2; i++ )
  {
    c[i] = 0.5 * dt * pow ( sigma * ( i + 1 ), 2 ) + 0.5 * dt * r * ( i + 1 );
  }

  for ( i = 1; i < nx - 1; i++ )
  {
    a[i] = 0.5 * dt * pow ( sigma * ( i + 1 ), 2 ) - 0.5 * dt * r * ( i + 1 );
  }

  u = new double[(nx-1)*(nt+1)];

  u0 = 0.0;
  for ( i = 0; i < nx - 1; i++ )
  {
    u0 = u0 + dx;
    u[i+0*(nx-1)] = r8_max ( u0 - e, 0.0 );
  }
  
  for ( j = 0; j < nt; j++ )
  {
    t = ( double ) ( j ) * t1 / ( double ) ( nt );

    p = 0.5 * dt * ( nx - 1 ) * ( sigma * sigma * ( nx - 1 ) + r ) 
      * ( smax - e * exp ( - r * t ) );

    for ( i = 0; i < nx - 1; i++ )
    {
      u[i+(j+1)*(nx-1)] = b[i] * u[i+j*(nx-1)];
    }
    for ( i = 0; i < nx - 2; i++ )
    {
      u[i+(j+1)*(nx-1)] = u[i+(j+1)*(nx-1)] + c[i] * u[i+1+j*(nx-1)];
    }
    for ( i = 1; i < nx - 1; i++ )
    {
      u[i+(j+1)*(nx-1)] = u[i+(j+1)*(nx-1)] + a[i] * u[i-1+j*(nx-1)];
    }
    u[nx-2+(j+1)*(nx-1)] = u[nx-2+(j+1)*(nx-1)] + p;
  }

  delete [] a;
  delete [] b;
  delete [] c;

  return u;
}
//****************************************************************************80

double *mc ( double s0, double e, double r, double sigma, double t1, int m, 
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    MC uses Monte Carlo valuation on a European call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    Original MATLAB version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    Black-Scholes for Scientific Computing Students,
//    Computing in Science and Engineering,
//    November/December 2004, Volume 6, Number 6, pages 72-79.
//
//  Parameters:
//
//    Input, double S0, the asset price at time 0.
//
//    Input, double E, the exercise price.
//
//    Input, double R, the interest rate.
//
//    Input, double SIGMA, the volatility of the asset.
//
//    Input, double T1, the expiry date.
//
//    Input, int M, the number of simulations.
//
//    Input/output, int SEED, a seed for the random
//    number generator.
//
//    Output, double MC[2], the estimated range of the valuation.
//
{
  double *conf;
  int i;
  double pmean;
  double *pvals;
  double std;
  double *svals;
  double *u;
  double width;

  u = r8vec_normal_01_new ( m, seed );

  svals = new double[m];

  for ( i = 0; i < m; i++ )
  {
    svals[i] = s0 * exp ( ( r - 0.5 * sigma * sigma ) * t1 
      + sigma * sqrt ( t1 ) * u[i] );
  }

  pvals = new double[m];

  for ( i = 0; i < m; i++ )
  {
    pvals[i] = exp ( - r * t1 ) * r8_max ( svals[i] - e, 0.0 );
  }

  pmean = 0.0;
  for ( i = 0; i < m; i++ )
  {
    pmean = pmean + pvals[i];
  }
  pmean = pmean / ( double ) ( m );

  std = 0.0;
  for ( i = 0; i < m; i++ )
  {
    std = std + pow ( pvals[i] - pmean, 2 );
  }
  std = sqrt ( std / ( double ) ( m - 1 ) );

  width = 1.96 * std / sqrt ( ( double ) ( m ) );

  conf = new double[2];

  conf[0] = pmean - width;
  conf[1] = pmean + width;

  delete [] pvals;
  delete [] svals;
  delete [] u;

  return conf;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double *r8vec_normal_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, double R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, double Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
  int i;
  int m;
  static int made = 0;
  double pi = 3.141592653589793;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;
//
//  I'd like to allow the user to reset the internal data.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }

  x = new double[n];
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
//
//  Maybe we don't need any more values.
//
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );
    y =         sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * pi * r[1] );

    saved = 1;

    made = made + 2;

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
    y           = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
}
//****************************************************************************80

void r8vec_print_part ( int n, double a[], int max_print, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_PART prints "part" of an R8VEC.
//
//  Discussion:
//
//    The user specifies MAX_PRINT, the maximum number of lines to print.
//
//    If N, the size of the vector, is no more than MAX_PRINT, then
//    the entire vector is printed, one entry per line.
//
//    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
//    followed by a line of periods suggesting an omission,
//    and the last entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, int MAX_PRINT, the maximum number of lines
//    to print.
//
//    Input, string TITLE, a title.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i] << "\n";
    }
    cout << "  ........  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i] << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i] << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i]
         << "  " << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
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
//    19 August 2004
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
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void r8vec_write ( string output_filename, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_WRITE writes an R8VEC file.
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
//    10 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the data.
//
{
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8VEC_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    output << "  " << setw(24) << setprecision(16) << x[j] << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
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

