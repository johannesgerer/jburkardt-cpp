# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "fd1d_bvp.hpp"

//****************************************************************************80

double *fd1d_bvp  ( int n, double a ( double x ), double aprime ( double x ),
  double c ( double x ), double f ( double x ), double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_BVP solves a two point boundary value problem.
//
//  Discussion:
//
//    The program uses the finite difference method to solve a BVP
//    (boundary value problem) in one dimension.
//
//    The problem is defined on the region X[0] <= x <= X[N-1].
//
//    The following differential equation is imposed in the region:
//
//      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
//
//    where a(x), c(x), and f(x) are given functions.  We write out
//    the equation in full as
//
//      - a(x) * u''(x) - a'(x) * u'(x) + c(x) * u(x) = f(x)
//
//    At the boundaries, the following conditions are applied:
//
//      u(X[0]) = 0.0
//      u(X[N-1]) = 0.0
//
//    We replace the function U(X) by a vector of N values U associated
//    with the nodes.
//
//    The first and last values of U are determined by the boundary conditions.
//
//    At each interior node I, we write an equation to help us determine
//    U(I).  We do this by approximating the derivatives of U(X) by
//    finite differences.  Let us write XL, XM, and XR for X(I-1), X(I) and X(I+1).
//    Similarly we have UL, UM, and UR.  Other quantities to be evaluated at
//    X(I) = XM will also be labeled with an M:
//
//      - AM * ( UL - 2 UM + UR ) / DX^2 - A'M * ( UL - UR ) / ( 2 * DX ) = FM
//
//    These N-2 linear equations for the unknown coefficients complete the
//    linear system and allow us to compute the finite difference approximation
//    to the solution of the BVP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double A ( double x ), evaluates a(x);
//
//    Input, double APRIME ( double x ), evaluates a'(x);
//
//    Input, double C ( double x ), evaluates c(x);
//
//    Input, double F ( double x ), evaluates f(x);
//
//    Input, double X[N], the mesh points, which may be nonuniformly spaced.
//
//    Output, double FD1D_BVP[N], the value of the finite difference
//    approximation to the solution.
//
{
  double am;
  double apm;
  double cm;
  double fm;
  int i;
  double *rhs;
  double *tri;
  double *u;
  double xm;
//
//  Equation 1 is the left boundary condition, U(X[0]) = 0.0;
//
  tri = new double[3*n];
  rhs = new double[n];

  tri[0+0*3] = 0.0;
  tri[1+0*3] = 1.0;
  tri[2+0*3] = 0.0;
  rhs[0] = 0.0;
//
//  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1),
//  and so on.
//
  for ( i = 1; i < n - 1; i++ )
  {
    xm  = x[i];
    am  = a ( xm );
    apm = aprime ( xm );
    cm  = c ( xm );
    fm  = f ( xm );

    tri[0+i*3] = - 2.0 * am / ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] )
      + apm / ( x[i+1] - x[i-1] );

    tri[1+i*3] = + 2.0 * am / ( x[i] - x[i-1] ) / ( x[i+1] - x[i] )
      + cm;

    tri[2+i*3] = - 2.0 * am / ( x[i+1] - x[i] ) / ( x[i+1] - x[i-1] )
      - apm / ( x[i+1] - x[i-1] );

    rhs[i]   = fm;
  }
//
//  Equation N is the right boundary condition, U(X[N-1]) = 0.0;
//
  tri[0+(n-1)*3] = 0.0;
  tri[1+(n-1)*3] = 1.0;
  tri[2+(n-1)*3] = 0.0;
  rhs[n-1] = 0.0;
//
//  Solve the linear system.
//
  u = r83np_fs ( n, tri, rhs );

  delete [] rhs;
  delete [] tri;

  return u;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = - x;
  }
  return value;
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
//    R8VEC_EVEN returns N real values, evenly spaced between ALO and AHI.
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
//    17 February 2004
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
//    Normally, A(1) = ALO and A(N) = AHI.
//    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
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
    for ( i = 1; i <= n; i++ )
    {
      a[i-1] = ( ( double ) ( n - i     ) * alo
               + ( double ) (     i - 1 ) * ahi )
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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
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
