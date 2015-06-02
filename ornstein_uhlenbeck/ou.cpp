# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "ou.hpp"

//****************************************************************************80

void ou_euler ( double theta, double mu, double sigma, double x0, double tmax, 
  int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    OU_EULER applies the Euler method to the Ornstein-Uhlenbeck SDE.
//
//  Discussion:
//
//    The stochastic differential equation (SDE) is:
//
//      dx(t) = theta * ( mu - x(t) ) dt + sigma dW,   
//      x(0) = x0.
//
//    The discretized Brownian path uses a constant stepsize.
//
//    For an SDE of the form:
//
//      dx = f(x(t)) dt + g(x(t)) dW(t),
//
//    the Euler method has the form:
//
//      x(j) = x(j-1) + f(x(j-1)) * dt + g(x(j-1)) * dW(j-1)
//
//    Note that if SIGMA is zero, the problem becomes deterministic,
//    with solution:
//
//      x(t) = mu + ( x0 - mu ) * exp ( - theta * t )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546
//
//  Parameters:
//
//    Input, double THETA, MU, SIGMA, the value of problem parameters.
//
//    Input, double X0, the initial condition.  When studying many
//    realizations of this problem, it is usual for X0 to be chosen
//    from a normal distribution.
//
//    Input, double TMAX, the final time.
//
//    Input, int N, the number of time steps.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  string command_filename = "ou_euler_commands.txt";
  ofstream command_unit;
  string data_filename = "ou_euler_data.txt";
  ofstream data_unit;
  double dt;
  double *dw;
  int j;
  double *t;
  double *x;

  cout << "\n";
  cout << "OU_EULER:\n";
  cout << "  C++ version\n";
  cout << "  Use an Euler method to approximate the solution of\n";
  cout << "  the Ornstein-Uhlenbeck stochastic differential equation:\n";
  cout << "\n";
  cout << "    d x(t) = theta * ( mu - x(t) ) dt + sigma dW\n";
  cout << "\n";
  cout << "  with initial condition x(0) = x0.\n";
//
//  Set the discrete time stepsize.
//
  dt = tmax / ( double ) ( n );
//
//  Compute the Brownian increments.
//
  dw = r8vec_normal_01_new ( n, seed );
  for ( j = 0; j < n; j++ )
  {
    dw[j] = dw[j] * sqrt ( dt );
  }
//
//  Carry out the Euler approximate integration process.
//
  t = r8vec_linspace_new ( n + 1, 0.0, tmax );

  x = new double[n+1];
  x[0] = x0;
  for ( j = 1; j <= n; j++ )
  {
    x[j] = x[j-1] + dt * theta * ( mu - x[j-1] ) + sigma * dw[j-1];
  }
//
//  Create the plot data file.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j <= n; j++ )
  {
    data_unit << "  " << setw(14) << t[j]
              << "  " << setw(14) << x[j] << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";
//
//  Create the plot command file.
//
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'ou_euler.png'\n";
  command_unit << "set xlabel '<--- T --->'\n";
  command_unit << "set ylabel '<--- X(T) --->'\n";
  command_unit << "set title 'Euler Solution of Ornstein-Uhlenbeck SDE'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 lw 3 linecolor rgb 'blue'\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file '" << command_filename << "'\n";
//
//  Free memory.
//
  delete [] dw;
  delete [] t;
  delete [] x;

  return;
}
//****************************************************************************80

void ou_euler_maruyama ( double theta, double mu, double sigma, double x0, 
  double tmax, int n, int r, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    OU_EULER_MARUYAMA applies Euler-Maruyama to the Ornstein-Uhlenbeck SDE.
//
//  Discussion:
//
//    The stochastic differential equation (SDE) is:
//
//      dx = theta * ( mu - x(t) ) dt + sigma dW,   
//      x(0) = x0,
//
//    The discretized Brownian path uses a constant stepsize.
//
//    A "large" time step DT_LARGE is used for the smooth portion of the
//    equation, while a smaller time step DT_SMALL is used for the
//    discretized Brownian path.  We take R small steps to equal one 
//    large step, so that:
//
//      dt_large = r * dt_small = tmax / n
//
//    For an SDE of the form:
//
//      dx = f(x(t)) dt + g(x(t)) dW(t)
//
//    the Euler-Maruyama method has the form:
//
//      x(j) = x(j-1) + f(X(j-1)) * dt_large + g(X(j-1)) * dW(j-1)
//
//    where dW(j-1) is approximated by the sum of R normal random values
//    multiplied by the square root of DT_SMALL.
//
//    Note that if SIGMA is zero, the problem becomes deterministic,
//    with solution
//
//      x(t) = mu + ( x0 - mu ) * exp ( - theta * t )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546
//
//  Parameters:
//
//    Input, double THETA, MU, SIGMA, the value of problem parameters.
//
//    Input, double X0, the initial condition.  When studying many
//    realizations of this problem, it is usual for X0 to be chosen
//    from a normal distribution.
//
//    Input, double TMAX, the final time.
//
//    Input, int N, the number of large scale time steps.
//
//    Input, int R, the number of small scale time steps per single
//    large scale time step.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  string command_filename = "ou_euler_maruyama_commands.txt";
  ofstream command_unit;
  string data_filename = "ou_euler_maruyama_data.txt";
  ofstream data_unit;
  double dt_large;
  double dt_small;
  double *dw;
  double dw_sum;
  int i;
  int j;
  double *t;
  double *x;

  cout << "\n";
  cout << "OU_EULER_MARUYAMA:\n";
  cout << "  C++ version\n";
  cout << "  Use an Euler-Maruyama method to approximate the solution of\n";
  cout << "  the Ornstein-Uhlenbeck stochastic differential equation:\n";
  cout << "\n";
  cout << "    d x(t) = theta * ( mu - x(t) ) dt + sigma dW\n";
  cout << "\n";
  cout << "  with initial condition x(0) = x0.\n";
//
//  Set time steps.
//
  dt_large = tmax / ( double ) ( n );
  dt_small = tmax / ( double ) ( n ) / ( double ) ( r );
//
//  Carry out the Euler-Maruyama approximate integration process.
//
  t = r8vec_linspace_new ( n + 1, 0.0, tmax );

  x = new double[n+1];
  x[0] = x0;
  for ( j = 1; j <= n; j++ )
  {
    dw = r8vec_normal_01_new ( r, seed );
    dw_sum = 0.0;
    for ( i = 0; i < r; i++ )
    {
      dw_sum = dw_sum + dw[i];
    }
    dw_sum = dw_sum * sqrt ( dt_small );
    x[j] = x[j-1] + dt_large * theta * ( mu - x[j-1] ) + sigma * dw_sum;
    delete [] dw;
  }
//
//  Plot the approximate solution.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j <= n; j++ )
  {
    data_unit << "  " << setw(14) << t[j]
              << "  " << setw(14) << x[j] << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";

  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'ou_euler_maruyama.png'\n";
  command_unit << "set xlabel '<--- T --->'\n";
  command_unit << "set ylabel '<--- X(T) --->'\n";
  command_unit << "set title 'Euler-Maruyama Solution of Ornstein-Uhlenbeck SDE'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 lw 3 linecolor rgb 'blue'\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file '" << command_filename << "'\n";
//
//  Free memory.
//
  delete [] t;
  delete [] x;

  return;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

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
//    09 April 2012
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
//
//    In other words, the interval is divided into N-1 even subintervals,
//    and the endpoints of intervals are used as the points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

double *r8vec_normal_01_new ( int n, int &seed )

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
//    Input/output, int &SEED, a seed for the random number generator.
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

double *r8vec_uniform_01_new ( int n, int &seed )

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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

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
