# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

using namespace std;

# include "fd1d_heat_explicit.hpp"

//****************************************************************************80

double *fd1d_heat_explicit ( int x_num, double x[], double t, double dt, 
  double cfl, double *rhs ( int x_num, double x[], double t ), 
  void bc ( int x_num, double x[], double t, double h[] ), double h[] )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_HEAT_EXPLICIT: Finite difference solution of 1D heat equation.
//
//  Discussion:
//
//    This program takes one time step to solve the 1D heat equation 
//    with an explicit method.
//
//    This program solves
//
//      dUdT - k * d2UdX2 = F(X,T)
//
//    over the interval [A,B] with boundary conditions
//
//      U(A,T) = UA(T),
//      U(B,T) = UB(T),
//
//    over the time interval [T0,T1] with initial conditions
//
//      U(X,T0) = U0(X)
//
//    The code uses the finite difference method to approximate the
//    second derivative in space, and an explicit forward Euler approximation
//    to the first derivative in time.
//
//    The finite difference form can be written as
//
//      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
//      ------------------  = F(X,T) + k *  ------------------------------------
//               dt                                   dx * dx
//
//    or, assuming we have solved for all values of U at time T, we have
//
//      U(X,T+dt) = U(X,T) 
//        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
//
//    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
//
//      cfl = k * dt / dx / dx
//
//    In order for accurate results to be computed by this explicit method,
//    the CFL coefficient must be less than 0.5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of points to use in the 
//    spatial dimension.
//
//    Input, double X(X_NUM), the coordinates of the nodes.
//
//    Input, double T, the current time.
//
//    Input, double DT, the size of the time step.
//
//    Input, double CFL, the Courant-Friedrichs-Loewy coefficient,
//    computed by FD1D_HEAT_EXPLICIT_CFL.
//
//    Input, double H[X_NUM], the solution at the current time.
//
//    Input, double *RHS ( int x_num, double x[], double t ), the function 
//    which evaluates the right hand side.
//
//    Input, void BC ( int x_num, double x[], double t, double h[] ), 
//    the function which evaluates the boundary conditions.
//
//    Output, double FD1D_HEAT_EXPLICIT[X_NUM)], the solution at time T+DT.
//
{
  double *f;
  double *h_new;
  int j;

  f = rhs ( x_num, x, t );

  h_new = new double[x_num];

  h_new[0] = 0.0;

  for ( j = 1; j < x_num - 1; j++ )
  {
    h_new[j] = h[j] + dt * f[j] 
      + cfl * (         h[j-1]
                - 2.0 * h[j]
                +       h[j+1] );
  }
  h_new[x_num-1] = 0.0;

  bc ( x_num, x, t + dt, h_new );

  delete [] f;

  return h_new;
}
//****************************************************************************80

double fd1d_heat_explicit_cfl ( double k, int t_num, double t_min, double t_max, 
  int x_num, double x_min, double x_max )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_HEAT_EXPLICIT_CFL: compute the Courant-Friedrichs-Loewy coefficient.
//
//  Discussion:
//
//    The equation to be solved has the form:
//
//      dUdT - k * d2UdX2 = F(X,T)
//
//    over the interval [X_MIN,X_MAX] with boundary conditions
//
//      U(X_MIN,T) = U_X_MIN(T),
//      U(X_MIN,T) = U_X_MAX(T),
//
//    over the time interval [T_MIN,T_MAX] with initial conditions
//
//      U(X,T_MIN) = U_T_MIN(X)
//
//    The code uses the finite difference method to approximate the
//    second derivative in space, and an explicit forward Euler approximation
//    to the first derivative in time.
//
//    The finite difference form can be written as
//
//      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
//      ------------------  = F(X,T) + k *  ------------------------------------
//               dt                                   dx * dx
//
//    or, assuming we have solved for all values of U at time T, we have
//
//      U(X,T+dt) = U(X,T) 
//        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
//
//    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
//
//      cfl = k * dt / dx / dx
//
//    In order for accurate results to be computed by this explicit method,
//    the CFL coefficient must be less than 0.5!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Lindfield, John Penny,
//    Numerical Methods Using MATLAB,
//    Second Edition,
//    Prentice Hall, 1999,
//    ISBN: 0-13-012641-1,
//    LC: QA297.P45.
//
//  Parameters:
//
//    Input, double K, the heat conductivity coefficient.
//
//    Input, int T_NUM, the number of time values, including 
//    the initial value.
//
//    Input, double T_MIN, T_MAX, the minimum and maximum times.
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X_MIN, X_MAX, the minimum and maximum spatial 
//    coordinates.
//
//    Output, double FD1D_HEAT_EXPLICIT_CFL, the Courant-Friedrichs-Loewy coefficient.
//
{
  double cfl;
  double dx;
  double dt;

  dx = ( x_max - x_min ) / ( double ) ( x_num - 1 );
  dt = ( t_max - t_min ) / ( double ) ( t_num - 1 );
//
//  Check the CFL condition, print out its value, and quit if it is too large.
//
  cfl = k * dt / dx / dx;

  cout << "\n";
  cout << "  CFL stability criterion value = " << cfl << "\n";

  if ( 0.5 <= cfl )
  {
    cerr << "\n";
    cerr << "FD1D_HEAT_EXPLICIT_CFL - Fatal error!\n";
    cerr << "  CFL condition failed.\n";
    cerr << "  0.5 <= K * dT / dX / dX = CFL.\n";
    exit ( 1 );
  }

  return cfl;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
