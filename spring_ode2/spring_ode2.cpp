# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

int main ( );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPRING_ODE2.
//
//  Discussion:
//
//    This is a revision of the SPRING_ODE code.
//
//    In this revision of the program, we want to use vectors (C arrays) to 
//    store the data, and we want to write the data out to a file in a form 
//    that Gnuplot (or other plotting programs) can use.
//
//    Hooke's law for a spring observes that the restoring force is
//    proportional to the displacement: F = - k x
//
//    Newton's law relates the force to acceleration: F = m a
//
//    Putting these together, we have
//
//      m * d^2 x/dt^2 = - k * x
//
//    We can add a damping force with coefficient c:
//
//      m * d^2 x/dt^2 = - k * x - c * dx/dt
//
//    If we write this as a pair of first order equations for (x,v), we have
//
//          dx/dt = v
//      m * dv/dt = - k * x - c * v
//
//    and now we can approximate these values for small time steps.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2012
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
  float c;
  float dt;
  int i;
  float k;
  float m;
  int n = 101;
  float t[101];
  float t_final;
  float v[101];
  float x[101];

  timestamp ( );
  cerr << "\n";
  cerr << "SPRING_ODE2\n";
  cerr << "  C++ version\n";
  cerr << "  Approximate the solution of a spring equation.\n";
  cerr << "  Write data to a file for use by gnuplot.\n";
  cerr << "\n";
//
//  Data
//
  m = 1.0;
  k = 1.0;
  c = 0.3;
  t_final = 20.0;
  dt = t_final / ( float ) ( n - 1 );
//
//  Store the initial conditions in entry 0.
//
  t[0] = 0.0;
  x[0] = 1.0;
  v[0] = 0.0;
//
//  Compute the approximate solution at equally spaced times in entries 1 through N-1.
//
  for ( i = 1; i < n; i++ )
  {
    t[i] = ( float ) ( i ) * t_final / ( float ) ( n - 1 );
    x[i] = x[i-1] + dt * v[i-1];
    v[i] = v[i-1] + ( dt / m ) * ( - k * x[i-1] - c * v[i-1] );
  }
//
//  Write the data to a file for plotting, possibly by Gnuplot.
//  Gnuplot expects T, X and V to be columns of output.
//
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << t[i]
         << "  " << x[i]
         << "  " << v[i] << "\n";
  }
//
//  Terminate.
//
  cerr << "\n";
  cerr << "SPRING_ODE2:\n";
  cerr << "  Normal end of execution.\n";
  cerr << "\n";
  timestamp ( );

  return 0;
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

  std::cerr << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
