# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

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
//    09 October 2013
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
  double c;
  string command_filename = "spring_ode2_commands.txt";
  ofstream command_unit;
  string data_filename = "spring_ode2_data.txt";
  ofstream data_unit;
  double dt;
  int i;
  double k;
  double m;
  int n = 101;
  double t[101];
  double t_final;
  double v[101];
  double x[101];

  timestamp ( );
  cout << "\n";
  cout << "SPRING_ODE2\n";
  cout << "  C++ version\n";
  cout << "  Approximate the solution of a spring equation.\n";
  cout << "  Write data to a file for use by gnuplot.\n";
//
//  Data
//
  m = 1.0;
  k = 1.0;
  c = 0.3;
  t_final = 20.0;
  dt = t_final / ( double ) ( n - 1 );
//
//  Store the initial conditions in entry 0.
//
  t[0] = 0.0;
  x[0] = 1.0;
  v[0] = 0.0;
//
//  Compute the approximate solution at equally spaced times 
//  in entries 1 through N-1.
//
  for ( i = 1; i < n; i++ )
  {
    t[i] = ( double ) ( i ) * t_final / ( double ) ( n - 1 );
    x[i] = x[i-1] + dt * v[i-1];
    v[i] = v[i-1] + ( dt / m ) * ( - k * x[i-1] - c * v[i-1] );
  }
//
//  Create the plot data file.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << setw(14) << t[i]
              << "  " << setw(14) << x[i]
              << "  " << setw(14) << v[i] << "\n";
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
  command_unit << "set output 'xv_time.png'\n";
  command_unit << "set xlabel '<--- T --->'\n";
  command_unit << "set ylabel '<--- X(T), V(T) --->'\n";
  command_unit << "set title 'Position and Velocity versus Time'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 lw 3 linecolor rgb 'blue',"
               << " '' using 1:3 lw 3 linecolor rgb 'red'\n";
  command_unit << "set output 'xv_phase.png'\n";
  command_unit << "set xlabel '<--- X(T) --->'\n";
  command_unit << "set ylabel '<--- V(T) --->'\n";
  command_unit << "set title 'Position versus Velocity'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename 
               << "' using 2:3 lw 3 linecolor rgb 'green'\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file '" << command_filename << "'\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "SPRING_ODE2:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
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

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
