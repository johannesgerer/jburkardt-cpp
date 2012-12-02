# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

using namespace std;

# include "fd1d_heat_steady.hpp"

//****************************************************************************80

double *fd1d_heat_steady ( int n, double a, double b, double ua, double ub,
  double k ( double x ), double f ( double x ), double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_HEAT_STEADY solves the steady 1D heat equation.
//
//  Discussion:
//
//    This program seeks a solution of the steady heat equation:
//
//      - d/dx ( K(X) dUdx ) = F(X)
//
//    over the interval [A,B] with boundary conditions
//
//      U(A) = UA,
//      U(B) = UB.
//
//    The code uses the finite difference method to approximate the
//    second derivative in space.  This results in a sparse linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of grid points.
//
//    Input, double A, B, the interval endpoints.
//
//    Input, double UA, UB, the values prescribed for U at the endpoints.
//
//    Input, double K ( double X ), evaluates the thermal conductance at the N
//    points X.  Set K(X) = 1 if you don't care about this coefficient.
//
//    Input, double F ( double X ), evaluates the heat source term at the N
//    points X.  Set F(X) = 0 if you don't want any heat sources.
//
//    Input, double X[N], the grid points.
//
//    Output, double FD1D_HEAT_STEADY[N], the approximation to the solution
//    at the grid points.
//
{
  double dx;
  int i;
  double *rhs;
  double *tri;
  double *u;
  double xm;
  double xp;

  cout << "\n";
  cout << "FD1D_HEAT_STEADY\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Finite difference solution of\n";
  cout << "  the steady 1D heat equation\n";
  cout << "\n";
  cout << "    - d/dx ( k(x) dUdx ) = F(x)\n";
  cout << "\n";
  cout << "  for space interval A <= X <= B with boundary conditions\n";
  cout << "\n";
  cout << "    U(A) = UA\n";
  cout << "    U(B) = UB\n";
  cout << "\n";
  cout << "  A second order difference approximation is used.\n";
//
//  Set the spacing.
//
  dx = ( b - a ) / ( double ) ( n - 1 );
//
//  Set up the tridiagonal matrix.
//
  tri = new double[3*n];
  rhs = new double[n];

  tri[0+0*3] = 0.0;
  tri[1+0*3] = 1.0;
  tri[2+0*3] = 0.0;
  rhs[0] = ua;

  for ( i = 1; i < n - 1; i++ )
  {
    xm = ( x[i-1] + x[i] ) / 2.0;
    xp = ( x[i] + x[i+1] ) / 2.0;

    tri[0+i*3] = - k ( xm )              / dx / dx;
    tri[1+i*3] = ( k ( xm ) + k ( xp ) ) / dx / dx;
    tri[2+i*3] =            - k ( xp )   / dx / dx;

    rhs[i] = f ( x[i] );
  }

  tri[0+(n-1)*3] = 0.0;
  tri[1+(n-1)*3] = 1.0;
  tri[2+(n-1)*3] = 0.0;
  rhs[n-1] = ub;
//
//  Solve the linear system.
//
  u = r83np_fs ( n, tri, rhs );

  delete [] tri;
  delete [] rhs;

  return u;
}
//****************************************************************************80

double *r83np_fs ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83NP_FS factors and solves an R83NP system.
//
//  Discussion:
//
//    The R83NP storage format is used for a tridiagonal matrix.
//    The subdiagonal   is in entries (0,1:N-1),
//    the diagonal      is in entries (1,0:N-1),
//    the superdiagonal is in entries (2,0:N-2).
//
//    This algorithm requires that each diagonal entry be nonzero.
//    It does not use pivoting, and so can fail on systems that
//    are actually nonsingular.
//
//    The "R83NP" format used for this routine is different from the R83 format.
//    Here, we insist that the nonzero entries
//    for a given row now appear in the corresponding column of the
//    packed array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A21 A32 A43 A54
//      A11 A22 A33 A44 A55
//      A12 A23 A34 A45  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, double A[3*N].
//    On input, the nonzero diagonals of the linear system.
//    On output, the data in these vectors has been overwritten
//    by factorization information.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R83NP_FS[N], the solution of the linear system.
//
{
  int i;
  double *x;
//
//  Check.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      cerr << "\n";
      cerr << "R83NP_FS - Fatal error!\n";
      cerr << "  A[1+" << i << "*3] = 0.\n";
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
    a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3] / a[1+(i-1)*3];
    x[i]     = x[i]     - x[i-1]       * a[0+i*3] / a[1+(i-1)*3];
  }

  x[n-1] = x[n-1] / a[1+(n-1)*3];
  for ( i = n-2; 0 <= i; i-- )
  {
    x[i] = ( x[i] - a[2+i*3] * x[i+1] ) / a[1+i*3];
  }

  return x;
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

double *r8vec_even ( int n, double alo, double ahi )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN returns an R8VEC of values evenly spaced between ALO and AHI.
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
//    18 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, double ALO, AHI, the low and high values.
//
//    Output, double R8VEC_EVEN[N], N evenly spaced values.
//    Normally, A[0] = ALO and A[N-1] = AHI.
//    However, if N = 1, then A[0] = 0.5*(ALO+AHI).
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - i - 1 ) * alo
             + ( double ) (     i     ) * ahi )
             / ( double ) ( n     - 1 );
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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
