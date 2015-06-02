# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

int main ( );
double *r8vec_linspace_new ( int n, double a, double b );
void timestamp ( );
double *trisolve ( int n, double a[], double b[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_ADVECTION_DIFFUSION_STEADY solves steady advection diffusion equation.
//
//  Discussion:
//
//    The steady advection diffusion equation has the form:
//
//      v ux - k * uxx = 0
//
//    where V (the advection velocity) and K (the diffusivity) are positive 
//    constants, posed in the region
//
//      a = 0 < x < 1 = b
//
//    with boundary conditions
//
//      u(0) = 0, u(1) = 1.
//
//    The discrete solution is unreliable when dx > 2 * k / v / ( b - a ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *a3;
  double b;
  string command_filename = "fd1d_advection_diffusion_steady_commands.txt";
  ofstream command_unit;
  string data_filename = "fd1d_advection_diffusion_steady_data.txt";
  ofstream data_unit;
  double dx;
  double *f;
  int i;
  int j;
  double k;
  int nx;
  double r;
  double *u;
  double v;
  double *w;
  double *x;

  timestamp ( );
  cout << "\n";
  cout << "FD1D_ADVECTION_DIFFUSION_STEADY:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Solve the 1D steady advection diffusion equation:,\n";
  cout << "    v du/dx - k d2u/dx2 = 0\n";
  cout << "  with constant, positive velocity V and diffusivity K\n";
  cout << "  over the interval:\n";
  cout << "    0.0 <= x <= 1.0\n";
  cout << "  with boundary conditions:\n";
  cout << "    u(0) = 0, u(1) = 1.\n";
  cout << "\n";
  cout << "  Use finite differences\n";
  cout << "   d u/dx  = (u(t,x+dx)-u(t,x-dx))/2/dx\n";
  cout << "   d2u/dx2 = (u(x+dx)-2u(x)+u(x-dx))/dx^2\n";
//
//  Physical constants.
//
  v = 1.0;
  k = 0.05;
  cout << "\n";
  cout << "  Diffusivity K = " << k << "\n";
  cout << "  Velocity V    = " << v << "\n";
//
//  Spatial discretization.
//
  nx = 101;
  a = 0.0;
  b = 1.0;
  dx = ( b - a ) / ( double ) ( nx - 1 );
  x = r8vec_linspace_new ( nx, a, b );

  cout << "  Number of nodes NX = " << nx << "\n";
  cout << "  DX = " << dx << "\n";
  cout << "  Maximum safe DX is " << 2.0 * k / v / ( b - a ) << "\n";
//
//  Set up the tridiagonal linear system corresponding to the boundary 
//  conditions and advection-diffusion equation.
//
  a3 = new double[nx*3];
  f = new double[nx];

  a3[0+1*nx] = 1.0;
  f[0] = 0.0;

  for ( i = 1; i < nx - 1; i++ )
  {
    a3[i+0*nx] = - v / dx / 2.0 -           k / dx / dx;
    a3[i+1*nx] =                    + 2.0 * k / dx / dx;
    a3[i+2*nx] = + v / dx / 2.0 -           k / dx / dx;
    f[i] = 0.0;
  }

  a3[nx-1+1*nx] = 1.0;
  f[nx-1] = 1.0;

  u = trisolve ( nx, a3, f );
//
//  The exact solution to the differential equation is known.
//
  r = v * ( b - a ) / k;

  w = new double[nx];

  for ( i = 0; i < nx; i++ )
  {
    w[i] = ( 1.0 - exp ( r * x[i] ) ) / ( 1.0 - exp ( r ) );
  }
//
//  Write data file.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j < nx; j++ )
  {
    data_unit << x[j] << "  " 
              << u[j] << "  " 
              << w[j] << "\n";
  }
  data_unit.close ( );

  cout << "\n";
  cout << "  Gnuplot data written to file '" << data_filename << "'.\n";
//
//  Write command file.
//
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "set term png\n";
  command_unit << "set output 'fd1d_advection_diffusion_steady.png'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "unset key\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---U(X)--->'\n";
  command_unit << "set title 'Exact: green line, Approx: red dots'\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 with points pt 7 ps 2,\\\n";
  command_unit << "'' using 1:3 with lines lw 3\n";
  command_unit << "quit\n";

  command_unit.close ( );

  cout << "  Gnuplot commands written to '" << command_filename << "'\n";
//
//  Free memory.
//
  delete [] a3;
  delete [] f;
  delete [] u;
  delete [] w;
  delete [] x;
//
//  Terminate.
//
  cout << "\n";
  cout << "FD1D_ADVECTION_DIFFUSION_STEADY\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
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
//****************************************************************************80

double *trisolve ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRISOLVE factors and solves a tridiagonal system.
//
//  Discussion:
//
//    The three nonzero diagonals of the N by N matrix are stored as 3
//    columns of an N by 3 matrix.
//
//  Example:
//
//    Here is how a tridiagonal matrix of order 5 would be stored:
//
//       *  A11 A12
//      A21 A22 A23
//      A32 A33 A34
//      A43 A44 A45
//      A54 A55  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, double A[N*3].
//    On input, the tridiagonal matrix.
//    On output, the data in these vectors has been overwritten
//    by factorization information.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double TRISOLVE[N], the solution of the linear system.
//
{
  int i;
  double *x;
  double xmult;
//
//  The diagonal entries can't be zero.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[i+1*n] == 0.0 )
    {
      cerr << "\n";
      cerr << "TRISOLVE - Fatal error!\n";
      cerr << "  A(" << i << ",2) = 0.\n";
      exit ( 1 );
    }
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( i = 1; i < n; i++ )
  {
    xmult = a[i+0*n] / a[i-1+1*n];
    a[i+1*n] = a[i+1*n] - xmult * a[i-1+2*n];
    x[i]   = x[i]   - xmult * x[i-1];
  }

  x[n-1] = x[n-1] / a[n-1+1*n];
  for ( i = n - 2; 0 <= i; i-- )
  {
    x[i] = ( x[i] - a[i+2*n] * x[i+1] ) / a[i+1*n];
  }

  return x;
}
