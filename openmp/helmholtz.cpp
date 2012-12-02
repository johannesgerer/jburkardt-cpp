# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );
void driver ( int m, int n, int it_max, double alpha, double omega, double tol );
void error_check ( int m, int n, double alpha, double u[], double f[] );
void jacobi ( int m, int n, double alpha, double omega, double u[], double f[], 
  double tol, int it_max );
double *rhs_set ( int m, int n, double alpha );
double u_exact ( double x, double y );
double uxx_exact ( double x, double y );
double uyy_exact ( double x, double y );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HELMHOLTZ.
//
//  Discussion:
//
//    HELMHOLTZ solves a discretized Helmholtz equation.
//
//    The two dimensional region given is:
//
//      -1 <= X <= +1
//      -1 <= Y <= +1
//
//    The region is discretized by a set of M by N nodes:
//
//      P(I,J) = ( X(I), Y(J) )
//
//    where, for 0 <= I <= M-1, 0 <= J <= N - 1, (C/C++ convention)
//
//      X(I) = ( 2 * I - M + 1 ) / ( M - 1 )
//      Y(J) = ( 2 * J - N + 1 ) / ( N - 1 )
//
//    The Helmholtz equation for the scalar function U(X,Y) is
//
//      - Uxx(X,Y) -Uyy(X,Y) + ALPHA * U(X,Y) = F(X,Y)
//
//    where ALPHA is a positive constant.  We suppose that Dirichlet
//    boundary conditions are specified, that is, that the value of
//    U(X,Y) is given for all points along the boundary.
//
//    We suppose that the right hand side function F(X,Y) is specified in 
//    such a way that the exact solution is
//
//      U(X,Y) = ( 1 - X^2 ) * ( 1 - Y^2 )
//
//    Using standard finite difference techniques, the second derivatives
//    of U can be approximated by linear combinations of the values
//    of U at neighboring points.  Using this fact, the discretized
//    differential equation becomes a set of linear equations of the form:
//
//      A * U = F
//
//    These linear equations are then solved using a form of the Jacobi 
//    iterative method with a relaxation factor.
//
//    Directives are used in this code to achieve parallelism.
//    All do loops are parallized with default 'static' scheduling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2009.
//
//  Author:
//
//    Original FORTRAN77 version by Joseph Robicheaux, Sanjiv Shah.
//    C++ version by John Burkardt.
//
{
  double alpha = 0.25;
  int it_max = 100;
  int m = 500;
  int n = 500;
  double omega = 1.1;
  double tol = 1.0E-08;
  double wtime;

  cout << "\n";
  cout << "HELMHOLTZ\n";
  cout << "  C++/OpenMP version\n";
  cout << "\n";
  cout << "  A program which solves the 2D Helmholtz equation.\n";

  cout << "\n";
  cout << "  This program is being run in parallel.\n";

  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  cout << "\n";
  cout << "  The region is [-1,1] x [-1,1].\n";
  cout << "  The number of nodes in the X direction is M = " << m << "\n";
  cout << "  The number of nodes in the Y direction is N = " << n << "\n";
  cout << "  Number of variables in linear system M * N  = " << m * n << "\n";
  cout << "  The scalar coefficient in the Helmholtz equation is ALPHA = " 
       << alpha << "\n";
  cout << "  The relaxation value is OMEGA = " << omega << "\n";
  cout << "  The error tolerance is TOL = " << tol << "\n";
  cout << "  The maximum number of Jacobi iterations is IT_MAX = " 
       << it_max << "\n";
//
//  Call the driver routine.
//
  wtime = omp_get_wtime ( );

  driver ( m, n, it_max, alpha, omega, tol );

  wtime = omp_get_wtime ( ) - wtime;

  cout << "\n";
  cout << "  Wall clock time in seconds = " << wtime << "\n";
//
//  Termiante.
//
  cout << "\n";
  cout << "HELMHOLTZ\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void driver ( int m, int n, int it_max, double alpha, double omega, double tol )

//****************************************************************************80
//
//  Purpose:
//
//    DRIVER allocates arrays and solves the problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2007
//
//  Author:
//
//    Original FORTRAN77 version by Joseph Robicheaux, Sanjiv Shah.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int M, N, the number of grid points in the 
//    X and Y directions.
//
//    Input, int IT_MAX, the maximum number of Jacobi 
//    iterations allowed.
//
//    Input, double ALPHA, the scalar coefficient in the
//    Helmholtz equation.
//
//    Input, double OMEGA, the relaxation parameter, which
//    should be strictly between 0 and 2.  For a pure Jacobi method,
//    use OMEGA = 1.
//
//    Input, double TOL, an error tolerance for the linear
//    equation solver.
//
{
  double *f;
  int i;
  int j;
  double *u;
//
//  Initialize the data.
//
  f = rhs_set ( m, n, alpha );

  u = new double[m*n];

# pragma omp parallel \
  shared ( m, n, u ) \
  private ( i, j )

# pragma omp for

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      u[i+j*m] = 0.0;
    }
  }
//
//  Solve the Helmholtz equation.
//
  jacobi ( m, n, alpha, omega, u, f, tol, it_max );
//
//  Determine the error.
//
  error_check ( m, n, alpha, u, f );

  delete [] f;
  delete [] u;

  return;
}
//****************************************************************************80

void error_check ( int m, int n, double alpha, double u[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ERROR_CHECK determines the error in the numerical solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2007
//
//  Author:
//
//    Original FORTRAN77 version by Joseph Robicheaux, Sanjiv Shah.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int M, N, the number of grid points in the 
//    X and Y directions.
//
//    Input, double ALPHA, the scalar coefficient in the
//    Helmholtz equation.  ALPHA should be positive.
//
//    Input, double U[M*N], the solution of the Helmholtz equation 
//    at the grid points.
//
//    Input, double F[M*N], values of the right hand side function 
//    for the Helmholtz equation at the grid points.
//
{
  double error_norm;
  int i;
  int j;
  double u_norm;
  double u_true;
  double u_true_norm;
  double x;
  double y;

  u_norm = 0.0;

# pragma omp parallel \
  shared ( m, n, u ) \
  private ( i, j )

# pragma omp for reduction ( + : u_norm )

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      u_norm = u_norm + u[i+j*m] * u[i+j*m];
    }
  }

  u_norm = sqrt ( u_norm );

  u_true_norm = 0.0;
  error_norm = 0.0;

# pragma omp parallel \
  shared ( m, n, u ) \
  private ( i, j, u_true, x, y )

# pragma omp for reduction ( + : error_norm, u_true_norm )

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x = ( double ) ( 2 * i - m + 1 ) / ( double ) ( m - 1 );
      y = ( double ) ( 2 * j - n + 1 ) / ( double ) ( n - 1 );
      u_true = u_exact ( x, y );
      error_norm = error_norm + ( u[i+j*m] - u_true ) * ( u[i+j*m] - u_true );
      u_true_norm = u_true_norm + u_true * u_true;
    }
  }

  error_norm = sqrt ( error_norm );
  u_true_norm = sqrt ( u_true_norm );

  cout << "\n";
  cout << "  Computed U l2 norm :       " << u_norm << "\n";
  cout << "  Computed U_EXACT l2 norm : " << u_true_norm << "\n";
  cout << "  Error l2 norm:             " << error_norm << "\n";

  return;
}
//****************************************************************************80

void jacobi ( int m, int n, double alpha, double omega, double u[], double f[], 
  double tol, int it_max )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI applies the Jacobi iterative method to solve the linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2007
//
//  Author:
//
//    Original FORTRAN77 version by Joseph Robicheaux, Sanjiv Shah.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int M, N, the number of grid points in the 
//    X and Y directions.
//
//    Input, double ALPHA, the scalar coefficient in the
//    Helmholtz equation.  ALPHA should be positive.
//
//    Input, double OMEGA, the relaxation parameter, which
//    should be strictly between 0 and 2.  For a pure Jacobi method,
//    use OMEGA = 1.
//
//    Input/output, double U(M,N), the solution of the Helmholtz
//    equation at the grid points.
//
//    Input, double F(M,N), values of the right hand side function 
//    for the Helmholtz equation at the grid points.
//
//    Input, double TOL, an error tolerance for the linear
//    equation solver.
//
//    Input, int IT_MAX, the maximum number of Jacobi 
//    iterations allowed.
//
{
  double ax;
  double ay;
  double b;
  double dx;
  double dy;
  double error;
  double error_norm;
  int i;
  int it;
  int j;
  double *u_old;
//
//  Initialize the coefficients.
//
  dx = 2.0 / ( double ) ( m - 1 );
  dy = 2.0 / ( double ) ( n - 1 );

  ax = - 1.0 / dx / dx;
  ay = - 1.0 / dy / dy;
  b  = + 2.0 / dx / dx + 2.0 / dy / dy + alpha;

  u_old = new double[m*n];

  for ( it = 1; it <= it_max; it++ )
  {
    error_norm = 0.0;
//
//  Copy new solution into old.
//
# pragma omp parallel \
  shared ( m, n, u, u_old ) \
  private ( i, j )

# pragma omp for

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        u_old[i+m*j] = u[i+m*j];
      }
    }
//
//  Compute stencil, residual, and update.
//
# pragma omp parallel \
  shared ( ax, ay, b, f, m, n, omega, u, u_old ) \
  private ( error, i, j )

# pragma omp for reduction ( + : error_norm )

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
//
//  Evaluate the residual.
//
        if ( i == 0 || i == m - 1 || j == 0 || j == n - 1 )
        {
          error = u_old[i+j*m] - f[i+j*m];
        }
        else
        {
          error = ( ax * ( u_old[i-1+j*m] + u_old[i+1+j*m] ) 
            + ay * ( u_old[i+(j-1)*m] + u_old[i+(j+1)*m] ) 
            + b * u_old[i+j*m] - f[i+j*m] ) / b;
        }
//
//  Update the solution.
//
        u[i+j*m] = u_old[i+j*m] - omega * error;
//
//  Accumulate the residual error.
//
        error_norm = error_norm + error * error;
      }
    }
//
//  Error check.
//
    error_norm = sqrt ( error_norm ) / ( double ) ( m * n );

    cout << "  " << setw(4) << it
         << "  Residual RMS " << error_norm << "\n";

    if ( error_norm <= tol )
    {
      break;
    }

  }

  cout << "\n";
  cout << "  Total number of iterations " << it << "\n";

  delete [] u_old;

  return;
}
//****************************************************************************80

double *rhs_set ( int m, int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    RHS_SET sets the right hand side F(X,Y).
//
//  Discussion:
//
//    The routine assumes that the exact solution and its second
//    derivatives are given by the routine EXACT.
//
//    The appropriate Dirichlet boundary conditions are determined
//    by getting the value of U returned by EXACT.
//
//    The appropriate right hand side function is determined by
//    having EXACT return the values of U, UXX and UYY, and setting
//
//      F = -UXX - UYY + ALPHA * U
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2007
//
//  Author:
//
//    Original FORTRAN77 version by Joseph Robicheaux, Sanjiv Shah.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int M, N, the number of grid points in the 
//    X and Y directions.
//
//    Input, double ALPHA, the scalar coefficient in the
//    Helmholtz equation.  ALPHA should be positive.
//
//    Output, double RHS[M*N], values of the right hand side function 
//    for the Helmholtz equation at the grid points.
//
{
  double *f;
  double f_norm;
  int i;
  int j;
  double x;
  double y;

  f = new double[m*n];

# pragma omp parallel \
  shared ( alpha, f, m, n ) \
  private ( i, j, x, y )
  {
# pragma for

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        f[i+j*m] = 0.0;
      }
    }
//
//  Set the boundary conditions.
//
# pragma omp for

    for ( i = 0; i < m; i++ )
    {
      j = 0;
      y = ( double ) ( 2 * j - n + 1 ) / ( double ) ( n - 1 );
      x = ( double ) ( 2 * i - m + 1 ) / ( double ) ( m - 1 );
      f[i+j*m] = u_exact ( x, y );
    }

# pragma omp for

    for ( i = 0; i < m; i++ )
    {
      j = n - 1;
      y = ( double ) ( 2 * j - n + 1 ) / ( double ) ( n - 1 );
      x = ( double ) ( 2 * i - m + 1 ) / ( double ) ( m - 1 );
      f[i+j*m] = u_exact ( x, y );
    }

# pragma omp for

    for ( j = 0; j < n; j++ )
    {
      i = 0;
      x = ( double ) ( 2 * i - m + 1 ) / ( double ) ( m - 1 );
      y = ( double ) ( 2 * j - n + 1 ) / ( double ) ( n - 1 );
      f[i+j*m] = u_exact ( x, y );
    }

# pragma omp for

    for ( j = 0; j < n; j++ )
    {
      i = m - 1;
      x = ( double ) ( 2 * i - m + 1 ) / ( double ) ( m - 1 );
      y = ( double ) ( 2 * j - n + 1 ) / ( double ) ( n - 1 );
      f[i+j*m] = u_exact ( x, y );
    }
//
//  Set the right hand side F.
//

# pragma omp for

    for ( j = 1; j < n - 1; j++ )
    {
      for ( i = 1; i < m - 1; i++ )
      {
        x = ( double ) ( 2 * i - m + 1 ) / ( double ) ( m - 1 );
        y = ( double ) ( 2 * j - n + 1 ) / ( double ) ( n - 1 );
        f[i+j*m] = - uxx_exact ( x, y ) - uyy_exact ( x, y ) + alpha * u_exact ( x, y );
      }
    }
  }

  f_norm = 0.0;

# pragma omp parallel \
  shared ( f, m, n ) \
  private ( i, j )

# pragma omp for reduction ( + : f_norm )

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      f_norm = f_norm + f[i+j*m] * f[i+j*m];
    }
  }
  f_norm = sqrt ( f_norm );

  cout << "\n";
  cout << "  Right hand side l2 norm = " << f_norm << "\n";

  return f;
}
//****************************************************************************80

double u_exact ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    U_EXACT returns the exact value of U(X,Y).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the point at which the values are needed.
//
//    Output, double U_EXACT, the value of the exact solution.
//
{
  double value;

  value = ( 1.0 - x * x ) * ( 1.0 - y * y );

  return value;
}
//****************************************************************************80

double uxx_exact ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    UXX_EXACT returns the exact second X derivative of the solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the point at which the values are needed.
//
//    Output, double UXX_EXACT, the exact second X derivative.
//
{
  double value;

  value = -2.0 * ( 1.0 + y ) * ( 1.0 - y );

  return value;
}
//****************************************************************************80

double uyy_exact ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    UYY_EXACT returns the exact second Y derivative of the solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the point at which the values are needed.
//
//    Output, double UYY_EXACT, the exact second Y derivative.
//
{
  double value;

  value = -2.0 * ( 1.0 + x ) * ( 1.0 - x );

  return value;
}
