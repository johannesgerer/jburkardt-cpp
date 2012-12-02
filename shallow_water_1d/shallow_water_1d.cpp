# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void boundary_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] );
void initial_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double *r8vec_linspace_new ( int n, double a_first, double a_last );
void r8vec_write ( string output_filename, int n, double x[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SHALLOW_WATER_1D.
//
//  Discussion:
//
//    SHALLOW_WATER_1D approximates the 1D shallow water equations.
//
//    This code can be considered a 1D version of Cleve Moler's shallow
//    water equation solver.
//
//    The version of the shallow water equations being solved here is in
//    conservative form, and omits the Coriolis force.  The state variables
//    are H (the height) and UH (the mass velocity).
//
//    The equations have the form
//
//      dH/dt + d UH/dx = 0
//
//      d UH/dt + d ( U^2 H + 1/2 g H^2 )/dx = 0
//
//    Here U is the ordinary velocity, U = UH/H, and g is the gravitational
//    acceleration.
//
//    The initial conditions are used to specify ( H, UH ) at an equally
//    spaced set of points, and then the Lax-Wendroff method is used to advance
//    the solution through a number of equally spaced points in time, with 
//    boundary conditions supplying the first and last spatial values.
//
//
//    Some input values will result in an unstable calculation that
//    quickly blows up.  This is related to the Courant-Friedrichs-Lewy
//    condition, which requires that DT be small enough, relative to DX and
//    the velocity, that information cannot cross an entire cell.
//
//    A "reasonable" set of input quantities is
//
//      shallow_water_1d 41 100 1.0 0.2 9.8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Cleve Moler,
//    "The Shallow Water Equations",
//    Experiments with MATLAB.
//
//  Parameters:
//
//    Input, integer NX, the number of spatial nodes.
//
//    Input, integer NT, the number of times steps.
//
//    Input, real X_LENGTH, the length of the region.
//
//    Input, real T_LENGTH, the time extent.
//
//    Input, real G, the gravity constant.  G = 9.8 meters per second**2.
//
//    Output, real H_ARRAY[NX*(NT+1)], the height for all space and time points.
//
//    Output, real UH_ARRAY[NX*(NT+1], the mass velocity for all space and time points.
//
//    Output, real X[NX], the X coordinates.
//
//    Output, real T[NT+1], the T coordinates.
//
{
  double dx;
  double dt;
  string filename_h;
  string filename_t;
  string filename_uh;
  string filename_x;
  double g;
  double *h;
  double *h_array;
  double *hm;
  int i;
  int it;
  int nt;
  int nx;
  double *t;
  double t_length;
  double *uh;
  double *uh_array;
  double *uhm;
  double *x;
  double x_length;

  timestamp ( );
  cout << "\n";
  cout << "SHALLOW_WATER_1D\n";
  cout << "  C++ version\n";
  cout << "\n";
//
//  Get the quadrature file root name:
//
  if ( argc <= 1 )
  {
    nx = 41;
  }
  else
  {
    nx = atoi ( argv[1] );
  }
  cout << "  NX = " << nx << "\n";

  if ( argc <= 2 )
  {
    nt = 100;
  }
  else
  {
    nt = atoi ( argv[2] );
  }
  cout << "  NT = " << nt << "\n";

  if ( argc <= 3 )
  {
    x_length = 1.0;
  }
  else
  {
    x_length = atof ( argv[3] );
  }
  cout << "  X_LENGTH = " << x_length << "\n";

  if ( argc <= 4 )
  {
    t_length = 0.2;
  }
  else
  {
    t_length = atof ( argv[4] );
  }
  cout << "  T_LENGTH = " << t_length << "\n";

  if ( argc <= 5 )
  {
    g = 9.8;
  }
  else
  {
    g = atof ( argv[5] );
  }
  cout << "  G = " << g << "\n";
//
//  Allocate space.
//
  h = new double[nx];
  h_array = new double[nx*(nt+1)];
  hm = new double[nx-1];
  t = new double[nt+1];
  uh = new double[nx];
  uh_array = new double[nx*(nt+1)];
  uhm = new double[nx-1];
  x = new double[nx];
//
//  Define the locations of the nodes and time steps and the spacing.
//
  x = r8vec_linspace_new ( nx, 0.0, x_length );
  t = r8vec_linspace_new ( nt + 1, 0.0, t_length );

  dx = x_length / ( double ) ( nx - 1 );
  dt = t_length / ( double ) ( nt );
//
//  Apply the initial conditions.
//
  initial_conditions ( nx, nt, x, t[0], h, uh );
//
//  Apply the boundary conditions.
//
  boundary_conditions ( nx, nt, x, t[0], h, uh );
//
//  Store the first time step into H_ARRAY and UH_ARRAY.
//
  for ( i = 0; i < nx; i++ )
  {
    h_array[i+0*nx] = h[i];
    uh_array[i+0*nx] = uh[i];
  }
//
//  Take NT more time steps.
//
  for ( it = 1; it <= nt; it++ )
  {
//
//  Take a half time step, estimating H and UH at the NX-1 spatial midpoints.
//
    for ( i = 0; i < nx - 1; i++ )
    {
      hm[i] = ( h[i] + h[i+1] ) / 2.0 
        - ( dt / 2.0 ) * ( uh[i+1] - uh[i] ) / dx;
    }
    for ( i = 0; i < nx - 1; i++ )
    {
      uhm[i] = ( uh[i] + uh[i+1] ) / 2.0 
        - ( dt / 2.0 ) * ( 
          uh[i+1] * uh[i+1] / h[i+1] + 0.5 * g * h[i+1] * h[i+1]
        - uh[i] * uh[i]  / h[i] - 0.5 * g * h[i] * h[i] ) / dx;
    }
//
//  Take a full time step, evaluating the derivative at the half time step,
//  to estimate the solution at the NX-2 nodes.
//
    for ( i = 1; i < nx - 1; i++ )
    {
      h[i] = h[i] 
        - dt * ( uhm[i] - uhm[i-1] ) / dx;
    }
    for ( i = 1; i < nx; i++ )
    {
      uh[i] = uh[i] 
        - dt * ( 
          uhm[i] * uhm[i]  / hm[i] + 0.5 * g * hm[i] * hm[i]
        - uhm[i-1] * uhm[i-1]  / hm[i-1] - 0.5 * g * hm[i-1] * hm[i-1] ) / dx;
    }
//
//  Update the boundary conditions.
//
    boundary_conditions ( nx, nt, x, t[it], h, uh );
//
//  Copy data into the big arrays.
//
    for ( i = 0; i < nx; i++ )
    {
      h_array[i+it*nx] = h[i];
      uh_array[i+it*nx] = uh[i];
    } 
  }
//
//  Write data to files.
//
  filename_x = "sw1d_x.txt";
  filename_t = "sw1d_t.txt";
  filename_h = "sw1d_h.txt";
  filename_uh = "sw1d_uh.txt";

  r8vec_write ( filename_x, nx, x );
  r8vec_write ( filename_t, nt + 1, t );
  r8mat_write ( filename_h, nx, nt + 1, h_array );
  r8mat_write ( filename_uh, nx, nt + 1, uh_array );

  cout << "\n";
  cout << "  X  values saved in file \"" << filename_x << "\".\n";
  cout << "  T  values saved in file \"" << filename_t << "\".\n";
  cout << "  H  values saved in file \"" << filename_h << "\".\n";
  cout << "  UH values saved in file \"" << filename_uh << "\".\n";
//
//  Free memory.
//
  delete [] h;
  delete [] h_array;
  delete [] hm;
  delete [] t;
  delete [] uh;
  delete [] uh_array;
  delete [] uhm;
  delete [] x;
//
//  Terminate.
//
  cout << "\n";
  cout << "SHALLOW_WATER_1D:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void boundary_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] )

//****************************************************************************80
//
//  Purpose:
//
//    INITIAL_CONDITIONS sets the initial conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, the number of spatial nodes.
//
//    Input, int NT, the number of times steps.
//
//    Input, double X[NX], the coordinates of the nodes.
//
//    Input, double T, the current time.
//
//    Input/output, double H[NX], the height, with H(1) and H(NX) 
//    adjusted for boundary conditions.
//
//    Input/output, double UH[NX], the mass velocity, with UH(1) 
//    and UH(NX) adjusted for boundary conditions.
//
{
  int bc;

  bc = 1;
//
//  Periodic boundary conditions on H and UH.
//
  if ( bc == 1 )
  {
    h[0]     = h[nx-2];
    h[nx-1]  = h[1];
    uh[0]    = uh[nx-2];
    uh[nx-1] = uh[1];
  }
//
//  Free boundary conditions on H and UH.
//
  else if ( bc == 2 )
  {
    h[0]     = h[1];
    h[nx-1]  = h[nx-2];
    uh[0]    = uh[1];
    uh[nx-1] = uh[nx-2];
  }
//
//  Reflective boundary conditions on UH, free boundary conditions on H.
//
  else if ( bc == 3 )
  {
    h[0]     =   h[1];
    h[nx-1]  =   h[nx-2];
    uh[0]    = - uh[1];
    uh[nx-1] = - uh[nx-2];
  }
  return;
}
//****************************************************************************80

void initial_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] )

//****************************************************************************80
//
//  Purpose:
//
//    INITIAL_CONDITIONS sets the initial conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, the number of spatial nodes.
//
//    Input, int NT, the number of times steps.
//
//    Input, double X[NX], the coordinates of the nodes.
//
//    Input, double T, the current time.
//
//    Output, double H[NX], the initial height for all space.
//
//    Output, double UH[NX], the initial mass velocity for all space.
//
{
  int i;
  double pi = 3.141592653589793;

  for ( i = 0; i < nx; i++ )
  {
    h[i] = 2.0 + sin ( 2.0 * pi * x[i] );
  }
  for ( i = 0; i < nx; i++ )
  {
    uh[i] = 0.0;
  }
  return;
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
