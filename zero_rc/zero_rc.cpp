# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "zero_rc.hpp"

//****************************************************************************80

double root_rc ( double x, double fx, double &ferr, double &xerr, double q[9] )

//****************************************************************************80
//
//  Purpose:
//
//    ROOT_RC solves a single nonlinear equation using reverse communication.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2013
//
//  Author:
//
//    Original FORTRAN77 version by Gaston Gonnet.
//    C++ version by John Burkardt.
//
//  Reference:
// 
//    Gaston Gonnet,
//    On the Structure of Zero Finders,
//    BIT Numerical Mathematics,
//    Volume 17, Number 2, June 1977, pages 170-183.
//
//  Parameters:
//
//    Input, double X, an estimate for the root.  On the first
//    call, this must be a value chosen by the user.  Thereafter, it may
//    be a value chosen by the user, or the value of ROOT returned on the
//    previous call to the function.
//
//    Input, double FX, the value of the function at X.
//
//    Output, double &FERR, the smallest value of F encountered.
//
//    Output, double &XERR, the width of the change-in-sign interval,
//    if one was encountered.
//
//    Input/output, double Q[9], storage needed by the function.
//    Before the first call, the user must set Q(1) to 0.
//
//    Output, double ROOT_RC, an improved estimate for the root.
//
{
  double d;
  double decr;
  int i;
  double p;
  double r;
  double u;
  double v;
  double w;
  double xnew;
  double z;
//
//  If we found an exact zero, there is nothing more to do.
//
  if ( fx == 0.0 )
  {
    ferr = 0.0;
    xerr = 0.0;
    xnew = x;
    return xnew;
  }

  ferr = r8_abs ( fx );
//
//  If this is the first time, initialize, estimate the first root, and exit.
//
  if ( q[0] == 0.0 )
  {
    q[0] = fx;
    q[1] = x;
    for ( i = 2; i < 9; i++ )
    {
      q[i] = 0.0;
    }
    xnew = x + fx;
    xerr = r8_huge ( );
    return xnew;
  }
//
//  This is not the first call.
//
  q[8] = q[8] + 1.0;
//
//  Check for too many iterations.
//
  if ( 80.0 < q[8] )
  {
    cout << "\n";
    cout << "ROOT_RC - Fatal error!\n";
    cout << "  Number of iterations = " << ( int ) q[8] << "\n";
    exit ( 1 );
  }
//
//  Check for a repeated X value.
//
  if ( ( 2.0 <= q[8] && x == q[3] ) || x == q[1] )
  {
    cout << "\n";
    cout << "ROOT_RC - Fatal error!\n";
    cout << "  Value of X has been input before.\n";
    exit ( 1 );
  }
//
//  Push X -> A -> B -> C
//
  for ( i = 5; 2 <= i; i-- )
  {
    q[i] = q[i-2];
  }
  q[0] = fx;
  q[1] = x;
//
//  If we have a change-in-sign interval, store the opposite value.
//
  if ( r8_sign ( q[0] ) != r8_sign ( q[2] ) )
  {
    q[6] = q[2];
    q[7] = q[3];
  }
//
//  Calculate XERR.
//
  if ( q[6] != 0.0 )
  {
    xerr = r8_abs ( q[7] - q[1] );
  }
  else
  {
    xerr = r8_huge ( );
  }
//
//  If more than 30 iterations, and we have change-in-sign interval, bisect.
//
  if ( 30.0 < q[8] && q[6] != 0.0 )
  {
    xnew = q[1] + ( q[7] - q[1] ) / 2.0;
    return xnew;
  }

  v = ( q[2] - q[0] ) / ( q[3] - q[1] );
//
//  If 3 or more points, try Muller.
//
  if ( q[4] != 0.0 )
  {
    u = ( q[4] - q[2] ) / ( q[5] - q[3] );
    w = q[3] - q[1];
    z = ( q[5] - q[1] ) / w;
    r = ( z + 1.0 ) * v - u;

    if ( r != 0.0 )
    {
      p = 2.0 * z * q[0] / r;
      d = 2.0 * p / ( w * r ) * ( v - u );
      if ( -1.0 <= d )
      {
        xnew = q[1] - p / ( 1.0 + sqrt ( 1.0 + d ) );
        if ( q[6] == 0.0 || 
             ( q[1] < xnew && xnew < q[7] ) || 
             ( q[7] < xnew && xnew < q[1] ) )
        {
          return xnew;
        }
      }
    }
  }
//
//  Try the secant step.
//
  if ( q[0] != q[2] || q[6] == 0.0 )
  {
    if ( q[0] == q[2] )
    {
      cout << "\n";
      cout << "ROOT_RC - Fatal error!\n";
      cout << "  Cannot apply any method.\n";
      exit ( 1 );
    }
    decr = q[0] / v;
    if ( r8_abs ( decr ) * 4.6E+18 < r8_abs ( q[1] ) )
    {
      decr = 1.74E-18 * r8_abs ( q[1] ) * r8_sign ( decr );
    }
    xnew = q[1] - decr;
    if ( q[6] == 0.0 || 
        ( q[1] < xnew && xnew < q[7] ) || 
        ( q[7] < xnew && xnew < q[1] ) )
    {
      return xnew;
    }
  }
//
//  Apply bisection.
//
  xnew = q[1] + ( q[7] - q[1] ) / 2.0;

  return xnew;
}
//****************************************************************************80

void roots_rc ( int n, double x[], double fx[], double &ferr, double xnew[], 
  double q[] )

//****************************************************************************80
//
//  Purpose:
//
//    ROOTS_RC solves a system of nonlinear equations using reverse communication.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2013
//
//  Author:
//
//    Original FORTRAN77 version by Gaston Gonnet.
//    C++ version by John Burkardt.
//
//  Reference:
//     
//    Gaston Gonnet,
//    On the Structure of Zero Finders,
//    BIT Numerical Mathematics,
//    Volume 17, Number 2, June 1977, pages 170-183.
//
//  Parameters:
//
//    Input, int N, the number of equations.
//
//    Input, double X[N].  Before the first call, the user should
//    set X to an initial guess or estimate for the root.  Thereafter, the input
//    value of X should be the output value of XNEW from the previous call.
//
//    Input, double FX[N], the value of the function at XNEW.
//
//    Output, double &FERR, the function error, that is, the sum of
//    the absolute values of the most recently computed function vector.
//
//    Output, double XNEW[N], a new point at which a function 
//    value is requested.
//
//    Workspace, double Q[(2*N+2)*(N+2)].  Before the first call 
//    for a given problem, the user must set Q(2*N+1,1) to 0.0.
//
{
  double damp;
  int i;
  int j;
  int jsma;
  int jsus;
  int lda;
  double sump;
  double t;

  lda = 2 * n + 2;

  ferr = 0.0;
  for ( i = 0; i < n; i++ )
  {
    ferr = ferr + r8_abs ( fx[i] );
  }
//
//  Initialization if Q(2*N+1,1) = 0.0.
//
  if ( q[2*n+1+0*lda] == 0.0 )
  {
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= n + 1; j++ )
      {
        q[i-1+(j-1)*lda] = 0.0;
        q[i+(j-1)*lda] = 0.0;
      }
      q[i-1+(i-1)*lda] = 100.0;
      q[i+n-1+(i-1)*lda] = 1.0;
    }

    for ( j = 1; j <= n; j++ )
    {
      q[2*n+(j-1)*lda] = r8_huge ( );
    }

    for ( j = 1; j <= n; j++ )
    {
      q[2*n+1+(j-1)*lda] = ( double ) ( n );
    }

    for ( i = 1; i <= n; i++ )
    {
      q[i+n-1+n*lda] = x[i-1];
    }

    for ( i = 1; i <= n; i++ )
    {
      q[i-1+n*lda] = fx[i-1];
    }

    q[2*n+n*lda] = ferr;
    q[2*n+1+n*lda] = 0.0;
    damp = 0.99;
  }
  else
  {
    jsus = 1;
    for ( i = 2; i <= n + 1; i++ )
    {
      if ( ( double ) ( 2 * n ) <= q[2*n+1+(i-1)*lda] )
      {
        q[2*n+(i-1)*lda] = r8_huge ( );
      }
      if ( q[2*n+1+(jsus-1)*lda] < ( n + 3 ) / 2 )
      {
        jsus = i;
      }
      if ( ( n + 3 ) / 2 <= q[2*n+1+(i-1)*lda] && 
        q[2*n+(jsus-1)*lda] < q[2*n+(i-1)*lda] ) 
      {
        jsus = i;
      }
    }

    for ( i = 1; i <= n; i++ )
    {
      q[i+n-1+(jsus-1)*lda] = x[i-1];
      q[i-1+(jsus-1)*lda] = fx[i-1];
    }

    q[2*n+(jsus-1)*lda] = ferr;
    q[2*n+1+(jsus-1)*lda] = 0;
    jsma = 1;
    damp = 0.0;

    for ( j = 1; j <= n + 1; j++ )
    {
      if ( r8_huge ( ) / 10.0 < q[2*n+(j-1)*lda] )
      {
        damp = 0.99;
      }
      if ( q[2*n+(j-1)*lda] < q[2*n+(jsma-1)*lda] )
      {
        jsma = j;
      }
    }

    if ( jsma != n + 1 )
    {
      for ( i = 1; i <= 2 * n + 2; i++ )
      {
        t = q[i-1+(jsma-1)*lda];
        q[i-1+(jsma-1)*lda] = q[i-1+n*lda];
        q[i-1+n*lda] = t;
      }
    }

  }

  for ( i = 1; i <= n; i++ )
  {
    q[i-1+(n+1)*lda] = q[i-1+n*lda];
  }
//
//  Call the linear equation solver, which should not destroy the matrix 
//  in Q(1:N,1:N), and should overwrite the solution into Q(1:N,N+2).
//
  r8mat_fs ( lda, n, q, q+(n+1)*lda );

  sump = 0.0;
  for ( i = 1; i <= n; i++ )
  {
    sump = sump + q[i-1+(n+1)*lda];
  }

  if ( r8_abs ( 1.0 - sump ) <= 1.0E-10 )
  {
    cerr << "\n";
    cerr << "ROOTS_RC - Fatal error!\n";
    cerr << "  SUMP almost exactly 1.\n";
    cerr << "  SUMP = " << sump << "\n";
    exit ( 1 );
  }

  for ( i = 1; i <= n; i++ )
  {
    xnew[i-1] = q[i+n-1+n*lda];
    for ( j = 1; j <= n; j++ )
    {
      xnew[i-1] = xnew[i-1] - q[i+n-1+(j-1)*lda] * q[j-1+(n+1)*lda];
    }
//
//  If system not complete, damp the solution.
//
    xnew[i-1] = xnew[i-1] / ( 1.0 - sump ) * ( 1.0 - damp ) + q[i+n-1+n*lda] * damp;
  }

  for ( j = 1; j <= n + 1; j++ )
  {
    q[2*n+1+(j-1)*lda] = q[2*n+1+(j-1)*lda] + 1.0;
  }

  return;
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
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

void r8mat_fs ( int lda, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_FS factors and solves a system with one right hand side.
//
//  Discussion:
//
//    This routine differs from R8MAT_FSS in two ways:
//    * only one right hand side is allowed;
//    * the input matrix A is not modified.
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LDA, the leading dimension of the matrix.
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N*N], the coefficient matrix of the linear system.
//    The matrix is stored in an LDAxN array.
//
//    Input/output, double X[N], on input, the right hand side of the
//    linear system.  On output, the solution of the linear system.
//
{
  double *a2;
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  a2 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*lda];
    }
  }

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = r8_abs ( a2[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a2[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a2[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cout << "\n";
      cout << "R8MAT_FS - Fatal error!\n";
      cout << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                  = a2[jcol-1+(j-1)*n];
        a2[jcol-1+(j-1)*n] = a2[ipiv-1+(j-1)*n];
        a2[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }
//
//  Scale the pivot row.
//
    t = a2[jcol-1+(jcol-1)*n];
    a2[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a2[jcol-1+(j-1)*n] = a2[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a2[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a2[i-1+(jcol-1)*n];
        a2[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a2[i-1+(j-1)*n] = a2[i-1+(j-1)*n] + t * a2[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a2[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  delete [] a2;

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
