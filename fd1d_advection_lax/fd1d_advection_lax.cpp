# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
double *initial_condition ( int nx, double x[] );
double *r8vec_linspace_new ( int n, double a, double b );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_ADVECTION_LAX solves the advection equation using the Lax method.
//
//  Discussion:
//
//    The Lax method is stable for the advection problem, if the time step
//    satisifies the Courant-Friedrichs-Levy (CFL) condition:
//
//      dt <= dx / c
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  string command_filename = "advection_commands.txt";
  ofstream command_unit;
  string data_filename = "advection_data.txt";
  ofstream data_unit;
  double dt;
  double dx;
  int i;
  int j;
  int jm1;
  int jp1;
  int nx;
  int nt;
  int nt_step;
  int plotstep;
  double t;
  double *u;
  double *unew;
  double *x;

  timestamp ( );
  cout << "\n";
  cout << "FD1D_ADVECTION_LAX:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Solve the constant-velocity advection equation in 1D,\n";
  cout << "    du/dt = - c du/dx\n";
  cout << "  over the interval:\n";
  cout << "    0.0 <= x <= 1.0\n";
  cout << "  with periodic boundary conditions, and\n";
  cout << "  with a given initial condition\n";
  cout << "    u(0,x) = (10x-4)^2 (6-10x)^2 for 0.4 <= x <= 0.6\n";
  cout << "           = 0 elsewhere.\n";
  cout << "\n";
  cout << "  We modify the FTCS method using the Lax method:\n";
  cout << "    du/dt = (u(t+dt,x)-0.5*u(t,x-dx)-0.5*u(t,x+dx))/dt\n";
  cout << "    du/dx = (u(t,x+dx)-u(t,x-dx))/2/dx\n";

  nx = 101;
  dx = 1.0 / ( double ) ( nx - 1 );
  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( nx, a, b );
  nt = 1000;
  dt = 1.0 / ( double ) ( nt );
  c = 1.0;

  u = initial_condition ( nx, x );
//
//  Open data file, and write solutions as they are computed.
//
  data_unit.open ( data_filename.c_str ( ) );

  t = 0.0;
  data_unit << "  " << x[0]
            << "  " << t
            << "  " << u[0] << "\n";
  for ( j = 0; j < nx; j++ )
  {
    data_unit << "  " << x[j]
              << "  " << t
              << "  " << u[j] << "\n";
  }
  data_unit << "\n";

  nt_step = 100;

  cout << "\n";
  cout << "  Number of nodes NX = " << nx << "\n";
  cout << "  Number of time steps NT = " << nt << "\n";
  cout << "  Constant velocity C = " << c << "\n";
  cout << "  CFL condition: dt (" << dt << ") <= dx / c (" << dx / c << ")\n";

  unew = new double[nx];

  for ( i = 0; i < nt; i++ )
  {
    for ( j = 0; j < nx; j++ )
    {
      jm1 = i4_wrap ( j - 1, 0, nx - 1 );
      jp1 = i4_wrap ( j + 1, 0, nx - 1 );
      unew[j] = 0.5 * u[jp1] + 0.5 * u[jm1] 
        - c * dt / dx / 2.0 * ( u[jp1] - u[jm1] );
    }
    for ( j = 0; j < nx; j++ )
    {
      u[j] = unew[j];
    }
    if ( i == nt_step - 1 )
    {
      t = ( double ) ( i ) * dt;
      for ( j = 0; j < nx; j++ )
      {
        data_unit << "  " << x[j]
                  << "  " << t
                  << "  " << u[j] << "\n";
      }
      data_unit << "\n";
      nt_step = nt_step + 100;
    }
  }
//
//  Close the data file once the computation is done.
//
  data_unit.close ( );

  cout << "\n";
  cout << "  Plot data written to the file \"" << data_filename << "\"\n";
//
//  Write gnuplot command file.
//
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "set term png\n";
  command_unit << "set output 'advection_lax.png'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "unset key\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Time--->'\n";
  command_unit << "splot '" << data_filename << "' using 1:2:3 with lines\n";
  command_unit << "quit\n";

  command_unit.close ( );

  cout << "  Gnuplot command data written to the file \"" << command_filename << "\"\n";
//
//  Free memory.
//
  delete [] u;
  delete [] unew;
  delete [] x;
//
//  Terminate.
//
  cout << "\n";
  cout << "FD1D_ADVECTION_LAX\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 December 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  if ( ilo <= ihi )
  {
    jlo = ilo;
    jhi = ihi;
  }
  else
  {
    jlo = ihi;
    jhi = ilo;
  }

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

double *initial_condition ( int nx, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    INITIAL_CONDITION sets the initial condition.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 December 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, the number of nodes.
//
//    Input, double X[NX], the coordinates of the nodes.
//
//    Output, double INITIAL_CONDITION[NX], the value of the initial condition.
//
{
  int i;
  double *u;

  u = new double[nx];

  for ( i = 0; i < nx; i++ )
  {
    if  ( 0.4 <= x[i] && x[i] <= 0.6 )
    {
      u[i] = pow ( 10.0 * x[i] - 4.0, 2 )
           * pow ( 6.0 - 10.0 * x[i], 2 );
    }
    else
    {
      u[i] = 0.0;
    }
  }
  return u;
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
