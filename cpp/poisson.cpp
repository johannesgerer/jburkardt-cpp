# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>

using namespace std;

int main ( int argc, char *argv[] );
void init_serial ( int nx, int ny, double **f, double **a );
void r8mat_delete ( double **a, int m, int n );
double **r8mat_new ( int m, int n );
double r8_max ( double x, double y );
double **sweep_serial ( int nx, int ny, double dx, double dy, double **f, 
  double **u );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POISSON_SERIAL.
//
//  Discussion:
//
//    POISSON is a program for solving the Poisson problem.
//
//    The Poisson equation
//
//      DEL^2 U(X,Y) = F(X,Y)
//
//    is solved on the unit square [0,1] x [0,1] using a grid of [0,NX+1] by
//    [0,NY+1] evenly spaced points.  
//
//    The boundary conditions and F are set so that the exact solution is
//
//      U(X,Y) = sin ( X * Y )
//
//    The Jacobi iteration is repeatedly applied until convergence is detected.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2007
//
//  Author:
//
//    John Burkardt
//
{
  int nx = 9;
  int ny = 9;

  double **a;
  double anorm;
  double **b;
  double bnorm;
  bool converged;
  double diff;
  double dx;
  double dy;
  double error;
  double **f;
  int i;
  int it;
  int it_max = 100;
  int j;
  double tolerance = 0.001;
  double u_exact;
  double x;
  double y;
//
//  Print a message.
//
  cout << "\n";
  cout << "POISSON:\n";
  cout << "  C++ version\n";
  cout << "  A program for solving the Poisson equation.\n";
  cout << "\n";
  cout << "  The number of interior X grid lines is " << nx << "\n";
  cout << "  The number of interior Y grid lines is " << ny << "\n";
//
//  Initialize the data.
//
  dx = 1.0 / ( double ) ( nx + 1 );
  dy = 1.0 / ( double ) ( ny + 1 );

  a = r8mat_new ( nx + 2, ny + 2 );
  f = r8mat_new ( nx + 2, ny + 2 );

  init_serial ( nx, ny, f, a );

  cout << "\n";
  cout << "  The X grid spacing is " << dx << "\n";
  cout << "  The Y grid spacing is " << dy << "\n";
//
//  Do the iteration.
//
  converged = false;

  cout << "\n";
  cout << "Step    ||U||         ||Unew||     ||Unew-U||     ||Unew-Exact||\n";
  cout << "\n";

  for ( it = 1; it <= it_max; it++ )
  {
//
//  Perform one Jacobi sweep, computing B from A.
//
    b = sweep_serial ( nx, ny, dx, dy, f, a );
//
//  Perform a second Jacobi sweep, computing A from B.
//
    r8mat_delete ( a, nx + 2, ny + 2 );

    a = sweep_serial ( nx, ny, dx, dy, f, b );
//
//  Check for convergence.
//
    anorm = 0.0;
    bnorm = 0.0;
    diff = 0.0;
    for ( j = 1; j <= ny; j++ )
    {
      for ( i = 1; i <= nx; i++ )
      {
        anorm = r8_max ( anorm, fabs ( a[i][j] ) );
        bnorm = r8_max ( bnorm, fabs ( b[i][j] ) );
        diff = r8_max ( diff, fabs ( a[i][j] - b[i][j] ) );
      }
    }
    r8mat_delete ( b, nx + 2, ny + 2 );

    error = 0.0;
    for ( j = 1; j <= ny; j++ )
    {
      y = ( double ) ( j ) / ( double ) ( ny + 1 );
      for ( i = 1; i <= nx; i++ )
      {
        x = ( double ) ( i ) / ( double ) ( nx + 1 );
        u_exact = sin ( x * y );
        error = r8_max ( error, fabs ( u_exact - a[i][j] ) );
      }
    }

    cout << "  " << setw(3)  << it
         << "  " << setw(12) << bnorm
         << "  " << setw(12) << anorm
         << "  " << setw(12) << diff
         << "  " << setw(12) << error << "\n";

    if ( diff <= tolerance )
    {
      converged = true;
      break;
    }

  }

  if ( converged )
  {
    cout << "\n";
    cout << "POISSON:\n";
    cout << "  The iteration has converged\n";
  }
  else
  {
    cout << "\n";
    cout << "POISSON:\n";
    cout << "  The iteration has NOT converged\n";
  }

  r8mat_delete ( a, nx + 2, ny + 2 );
  r8mat_delete ( f, nx + 2, ny + 2 );
//
//  Terminate.
//
  cout << "\n";
  cout << "POISSON:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void init_serial ( int nx, int ny, double **f, double **a )

//****************************************************************************80
//
//  Purpose:
//
//    INIT_SERIAL initializes the arrays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the X and Y grid dimensions.
//
//    Output, double F[0:NX+1][0:NY+1], the initialized right hand side data.
//
//    Output, double A[0:NX+1][0:NY+1], the initialized solution estimate.
//
{
  double fnorm;
  int i;
  int j;
  double u;
  double unorm;
  double x;
  double y;
//
//  Set the initial guesses to zero.
//
  for ( i = 0; i <= nx+1; i++ )
  {
    for ( j = 0; j <= ny+1; j++ )
    {
      a[i][j] = 0.0;
    }
  }
//
//  Just for debugging, compute the max norm of the exact solution 
//  on the interior nodes.
//
  unorm = 0.0;
  for ( j = 1; j <= ny; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny + 1 );
    for ( i = 1; i <= nx; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( nx + 1 );
      u = sin ( x * y );
      unorm = r8_max ( unorm, fabs ( u ) );
    }
  }
  cout << "\n";
  cout << "INIT_SERIAL:\n";
  cout << "  Max norm of exact solution U at interior nodes = " 
       << unorm << "\n";;
//
//  The "boundary" entries of F will store the boundary values of the solution.
//
  fnorm = 0.0;

  x = 0.0;
  for ( j = 0; j <= ny+1; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny + 1 );
    f[0][j] = sin ( x * y );
    a[0][j] = f[0][j];
    fnorm = r8_max ( fnorm, f[0][j] );
  }

  x = 1.0;
  for ( j = 0; j <= ny+1; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny + 1 );
    f[nx+1][j] = sin ( x * y );
    a[nx+1][j] = f[nx+1][j];
    fnorm = r8_max ( fnorm, f[nx+1][j] );
  }

  y = 0.0;
  for ( i = 0; i <= nx + 1; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( nx + 1 );
    f[i][0] = sin ( x * y );
    a[i][0] = f[i][0];
    fnorm = r8_max ( fnorm, f[i][0] );
  }

  y = 1.0;
  for ( i = 0; i <= nx + 1; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( nx + 1 );
    f[i][ny+1] = sin ( x * y );
    a[i][ny+1] = f[i][ny+1];
    fnorm = r8_max ( fnorm, f[i][ny+1] );
  }

  cout << "\n";
  cout << "INIT_SERIAL:\n";
  cout << "  Max norm of boundary values = " << fnorm << "\n";
//
//  The "interior" entries of F will store the right hand sides 
//  of the Poisson equation.
//
  fnorm = 0.0;
  for ( j = 1; j <= ny; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny + 1 );
    for ( i = 1; i <= nx; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( nx + 1 );
      f[i][j] = - ( x * x + y * y ) * sin ( x * y );
      fnorm = r8_max ( fnorm, fabs ( f[i][j] ) );
    }
  }

  cout << "\n";
  cout << "INIT_SERIAL:\n";
  cout << "  Max norm of right hand side F at interior nodes = "
       << fnorm << "\n";

  return;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}
//****************************************************************************80

void r8mat_delete ( double **a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DELETE frees memory associated with an R8MAT.
//
//  Discussion:
//
//    This function releases the memory associated with an array that was 
//    created by a command like
//      double **a;
//      a = r8mat_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double **A, the pointer to the array.
//
//    Input, int M, N, the number of rows and columns in the array.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}
//****************************************************************************80

double **r8mat_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NEW allocates a new R8MAT.
//
//  Discussion:
//
//    A declaration of the form
//      double **a;
//    is necesary.  Then an assignment of the form:
//      a = r8mat_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17.0;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, double R8MAT_NEW[M][N], a new matrix.
//
{
  double **a;
  int i;

  a = new double *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8MAT_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new double[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R8MAT_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}
//****************************************************************************80

double **sweep_serial ( int nx, int ny, double dx, double dy, double **f, 
  double **u )

//****************************************************************************80
//
//  Purpose:
//
//   SWEEP_SERIAL carries out one step of the Jacobi iteration.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the X and Y grid dimensions.
//
//    Input, double DX, DY, the spacing between grid points.
//
//    Input, double F[0:NX+1][0:NY+1], the right hand side data.
//
//    Input, double U[0:NX+1][0:NY+1], the previous solution estimate.
//
//    Output, double SWEEP_SERIAL[0:NX+1][0:NY+1], the updated solution
//    estimate.
//
{
  int i;
  int j;
  double **unew;

  unew = r8mat_new ( nx + 2, ny + 2 );

  for ( j = 0; j <= ny + 1; j++ )
  {
    for ( i = 0; i <= nx + 1; i++ )
    {
      if ( i == 0 || j == 0 || i == nx + 1 || j == ny + 1 )
      {
        unew[i][j] = u[i][j];
      }
      else
      { 
        unew[i][j] = ( u[i-1][j] + u[i][j+1] + u[i][j-1] + u[i+1][j] ) / 4.0 
          - f[i][j] * dx * dy;
      }
    }
  }
  return unew;
}
