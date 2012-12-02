# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "multigrid_poisson_1d.hpp"

//****************************************************************************80

void ctof ( int nc, double uc[], int nf, double uf[] )

//****************************************************************************80
//
//  Purpose:
//
//    CTOF transfers data from a coarse to a finer grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, int NC, the number of coarse nodes.
//
//    Input, double UC[NC], the coarse correction data.
//
//    Input, int NF, the number of fine nodes.
//
//    Input/output, double UF[NF], on input, the fine grid data.
//    On output, the data has been updated with prolonged coarse 
//    correction data.
//
{
  int ic;
  int iff;

  for ( ic = 0; ic < nc; ic++ )
  {
    iff = 2 * ic;
    uf[iff] = uf[iff] + uc[ic];
  }

  for ( ic = 0; ic < nc - 1; ic++ )
  {
    iff = 2 * ic + 1;
    uf[iff] = uf[iff] + 0.5 * ( uc[ic] + uc[ic+1] );
  }

  return;
}
//****************************************************************************80

void ftoc ( int nf, double uf[], double rf[], int nc, double uc[], 
  double rc[] )

//****************************************************************************80
//
//  Purpose:
//
//    FTOC transfers data from a fine grid to a coarser grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, int NF, the number of fine nodes.
//
//    Input, double UF[NF], the fine data.
//
//    Input, double RF[NF], the right hand side for the fine grid.
//
//    Input, int NC, the number of coarse nodes.
//
//    Output, double UC[NC], the coarse grid data, set to zero.
//
//    Output, double RC[NC], the right hand side for the coarse grid.
//
{
  int ic;
  int iff;

  for ( ic = 0; ic < nc; ic++ )
  {
    uc[ic] = 0.0;
  }

  rc[0] = 0.0;
  for ( ic = 1; ic < nc - 1; ic++ )
  {
    iff = 2 * ic;
    rc[ic] = 4.0 * ( rf[iff] + uf[iff-1] - 2.0 * uf[iff] + uf[iff+1] );
  }
  rc[nc-1] = 0.0;

  return;
}
//****************************************************************************80

void gauss_seidel ( int n, double r[], double u[], double &dif_l1 )

//****************************************************************************80
//
//  Purpose:
//
//    GAUSS_SEIDEL carries out one step of a Gauss-Seidel iteration.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, int N, the number of unknowns.
//
//    Input, double R[N], the right hand side.
//
//    Input/output, double U[N], the estimated solution.
//
//    Output, double &DIF_L1, the L1 norm of the difference between the
//    input and output solution estimates.
//
{
  int i;
  double u_old;

  dif_l1 = 0.0;

  for ( i = 1; i < n - 1; i++ )
  {
    u_old = u[i];
    u[i] = 0.5 * ( u[i-1] + u[i+1] + r[i] );
    dif_l1 = dif_l1 + r8_abs ( u[i] - u_old );
  }

  return;
}
//****************************************************************************80

int i4_log_2 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
//
//  Example:
//
//        I  I4_LOG_10
//    -----  --------
//        0    0
//        1    0
//        2    1
//        3    1
//        4    2
//        5    2
//        7    2
//        8    3
//        9    3
//     1000    9
//     1024   10
//
//  Discussion:
//
//    I4_LOG_2 ( I ) + 1 is the number of binary digits in I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number whose logarithm base 2 is desired.
//
//    Output, int I4_LOG_2, the integer part of the logarithm base 2 of
//    the absolute value of X.
//
{
  int i_abs;
  int two_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    two_pow = 2;

    i_abs = abs ( i );

    while ( two_pow <= i_abs )
    {
      value = value + 1;
      two_pow = two_pow * 2;
    }
  }

  return value;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

void monogrid_poisson_1d ( int n, double force ( double x ), 
  double exact ( double x ), int &it_num, double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOGRID_POISSON_1D solves a 1D PDE, using the Gauss-Seidel method.
//
//  Discussion:
//
//    This routine solves a 1D boundary value problem of the form
//
//      - U''(X) = F(X) for A < X < B,
//
//    with boundary conditions U(A) = UA, U(B) = UB.
//
//    The Gauss-Seidel method is used. 
//
//    This routine is provided primarily for comparison with the
//    multigrid solver.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, int N, the number of intervals.
//
//    Input, double FORCE ( double x ), the name of the function 
//    which evaluates the right hand side.
//
//    Input, double EXACT ( double x ), the name of the function 
//    which evaluates the exact solution.
//
//    Output, int &IT_NUM, the number of iterations.
//
//    Output, double U[N+1], the computed solution.
//
{
  double d1;
  double h;
  int i;
  double *r;
  double tol;
  double x;
//
//  Initialization.
//
  tol = 0.0001;

  r = new double[n+1];

  r[0] = 0.0;
  h = 1.0 / ( double ) ( n );
  for ( i = 1; i < n; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( n );
    r[i] = h * h * force ( x );
  }
  r[n] = 0.0;

  for ( i = 0; i <= n; i++ )
  {
    u[i] = 0.0;
  }
  it_num = 0;
//
//  Gauss-Seidel iteration.
//
  for ( ; ; )
  {
    it_num = it_num + 1;

    gauss_seidel ( n + 1, r, u, d1 );

    if ( d1 <= tol )
    {
      break;
    }
  }

  delete [] r;

  return;
}
//****************************************************************************80

void multigrid_poisson_1d ( int n, double force ( double x ), 
  double exact ( double x ), int &it_num, double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_POISSON_1D solves a 1D PDE using the multigrid method.
//
//  Discussion:
//
//    This routine solves a 1D boundary value problem of the form
//
//      - U''(X) = F(X) for A < X < B,
//
//    with boundary conditions U(A) = UA, U(B) = UB.
//
//    The multigrid method is used. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    Original FORTRAN77 version by William Hager.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, int N, the number of intervals.
//    N must be a power of 2.
//
//    Input, double FORCE ( double x ), the name of the function 
//    which evaluates the right hand side.
//
//    Input, double EXACT ( double x ), the name of the function 
//    which evaluates the exact solution.
//
//    Output, int &IT_NUM, the number of iterations.
//
//    Output, double U[N+1], the computed solution.
//
{
  double d0;
  double d1;
  double h;
  int i;
  int it;
  int j;
  int k;
  int l;
  int ll;
  int m;
  int nl;
  double *r;
  double s;
  double tol;
  double utol;
  double *uu;
  double x;
//
//  Determine if we have enough storage.
//
  k = i4_log_2 ( n );

  if ( n != i4_power ( 2, k ) )
  {
    cout << "\n";
    cout << "MULTIGRID_POISSON_1D - Fatal error!\n";
    cout << "  N is not a power of 2.\n";
    exit ( 1 );
  }

  nl = n + n + k - 2;
//
//  Initialization.
//
  it = 4;
  it_num = 0;
  tol = 0.0001;
  utol = 0.7;
  m = n;
//
//  Set the right hand side.
//
  r = new double[nl];
  r[0] = 0.0;
  h = 1.0 / ( double ) ( n );
  for ( i = 1; i < n; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( n );
    r[i] = h * h * force ( x );
  }
  r[n] = 0.0;

  uu = new double[nl];

  for ( i = 0; i < nl; i++ )
  {
    uu[i] = 0.0;
  }
//
//  L points to first entry of solution
//  LL points to penultimate entry.
//
  l = 0;
  ll = n - 1;
//
//  Gauss-Seidel iteration
//
  d1 = 0.0;
  j = 0;

  for ( ; ; )
  {
    d0 = d1;
    j = j + 1;
    gauss_seidel ( n + 1, r + l, uu + l, d1 );
    it_num = it_num + 1;
//
//  Do at least 4 iterations at each level.
//
    if ( j < it )
    {
      continue;
    }
//
//  Enough iterations, satisfactory decrease, on finest grid, exit.
//
    else if ( d1 < tol && n == m )
    {
      break;
    }
//
//  Enough iterations, satisfactory convergence, go finer.
//
    else if ( d1 < tol )
    {
      ctof ( n + 1, uu + l, n + n + 1, uu + (l-1-n-n) );

      n = n + n;
      ll = l - 2;
      l = l - 1 - n;
      j = 0;
    }
//
//  Enough iterations, slow convergence, 2 < N, go coarser.
//
    else if ( utol * d0 <= d1 && 2 < n )
    {
      ftoc ( n + 1, uu + l, r + l, (n/2)+1, uu+(l+n+1), r+(l+n+1) );

      n = n / 2;
      l = ll + 2;
      ll = ll + n + 1;
      j = 0;
    }
  }

  for ( i = 0; i < n + 1; i++ )
  {
    u[i] = uu[i];
  }
  delete [] r;
  delete [] uu;

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
//    18 August 2004
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
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
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
