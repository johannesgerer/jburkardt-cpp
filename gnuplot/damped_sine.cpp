# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <fstream>

using namespace std;

int main ( );
double *correlation_damped_sine ( int n, double rho[], double rho0 );
void correlation_plot ( int n, double rho[], double c[], string header, 
  string title );
double *r8vec_linspace_new ( int n, double a, double b );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    DAMPED_SINE evaluates and plots the damped sine correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int n = 101;
  double *rho;
  double rho0;

  cout << "\n";
  cout << "DAMPED_SINE\n";
  cout << "  C++ version\n";
  cout << "  Demonstrating how a correlation function can be\n";
  cout << "  evaluated and plotted using GNUPLOT.\n";

  rho0 = 1.0;
  rho = r8vec_linspace_new ( n, -12.0, 12.0 );
  c = correlation_damped_sine ( n, rho, rho0 );
  correlation_plot ( n, rho, c, "damped_sine", "Damped sine correlation" );

  delete [] rho;
  delete [] c;
//
//  Terminate.
//
  cout << "\n";
  cout << "DAMPED_SINE\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

double *correlation_damped_sine ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_DAMPED_SINE evaluates the damped sine correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double CORRELATION_DAMPED_SINE[N], the correlations.
//
{
  double *c;
  int i;
  double rhohat;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( rho[i] == 0.0 )
    {
      c[i] = 1.0;
    }
    else
    {
      rhohat = fabs ( rho[i] ) / rho0;
      c[i] = sin ( rhohat ) / rhohat;
    }
  }

  return c;
}
//****************************************************************************80

void correlation_plot ( int n, double rho[], double c[], string header, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_PLOT makes a plot of a correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double C[N], the correlations.
//
//    Input, string HEADER, an identifier for the files.
//
//    Input, string TITLE, a title for the plot.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  double rho0;

  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << rho[i]
              << "  " << c[i] << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";

  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output \"" << header << ".png\"\n";
  command_unit << "set xlabel 'Distance Rho'\n";
  command_unit << "set ylabel 'Correlation C(Rho)'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot \"" << data_filename << "\" using 1:2 lw 3 linecolor rgb \"blue\"\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file \"" << command_filename << "\"\n";

  return;
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
