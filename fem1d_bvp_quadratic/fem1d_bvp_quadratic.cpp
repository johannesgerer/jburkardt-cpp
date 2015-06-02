# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "fem1d_bvp_quadratic.hpp"

//****************************************************************************80

double *fem1d_bvp_quadratic ( int n, double a ( double x ), double c ( double x ), 
  double f ( double x ), double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_QUADRATIC solves a two point boundary value problem.
//
//  Discussion:
//
//    The program uses the finite element method, with piecewise quadratic basis
//    functions to solve a boundary value problem in one dimension.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double A ( double X ), evaluates a(x);
//
//    Input, double C ( double X ), evaluates c(x);
//
//    Input, double F ( double X ), evaluates f(x);
//
//    Input, double X[N], the mesh points.
//
//    Output, double FEM1D_BVP_QUADRATIC[N], the finite element coefficients, 
//    which are also the value of the computed solution at the mesh points.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956,
     0.000000000000000000000000000000,
     0.774596669241483377035853079956 };
  double al;
  double am;
  double ar;
  double *amat;
  double axq;
  double *b;
  double bm;
  double cxq;
  int e;
  int e_num;
  double fxq;
  int i;
  int ierror;
  int j;
  int l;
  int m;
  int q;
  int quad_num = QUAD_NUM;
  int r;
  double *u;
  double weight[QUAD_NUM] = { 
    0.555555555555555555555555555556,
    0.888888888888888888888888888889,
    0.555555555555555555555555555556 };
  double wq;
  double vl;
  double vlp;
  double vm;
  double vmp;
  double vr;
  double vrp;
  double xl;
  double xm;
  double xq;
  double xr;
//
//  Zero out the matrix and right hand side.
//
  amat = r8mat_zero_new ( n, n );
  b = r8vec_zero_new ( n );
//
//  Integrate over element E.
//
  e_num = ( n - 1 ) / 2;

  for ( e = 0; e < e_num; e++ )
  {
//
//  Element E uses nodes
//    L = 2 * E
//    M = 2 * E + 1
//    R = 2 * E + 2
//
    l = 2 * e;
    m = 2 * e + 1;
    r = 2 * e + 2;

    xl = x[l];
    xm = x[m];
    xr = x[r];

    for ( q = 0; q < quad_num; q++ )
    {

      xq = ( ( 1.0 - abscissa[q] ) * xl 
           + ( 1.0 + abscissa[q] ) * xr ) 
           /   2.0;

      wq = weight[q] * ( xr - xl ) / 2.0;

      axq = a ( xq );
      cxq = c ( xq );
      fxq = f ( xq );

      vl = ( ( xq - xm ) / ( xl - xm ) ) 
         * ( ( xq - xr ) / ( xl - xr ) );

      vm = ( ( xq - xl ) / ( xm - xl ) ) 
         * ( ( xq - xr ) / ( xm - xr ) );

      vr = ( ( xq - xl ) / ( xr - xl ) ) 
         * ( ( xq - xm ) / ( xr - xm ) );

      vlp = (         1.0 / ( xl - xm ) ) 
          * ( ( xq - xr ) / ( xl - xr ) ) 
          + ( ( xq - xm ) / ( xl - xm ) ) 
          * (         1.0 / ( xl - xr ) );

      vmp = (         1.0 / ( xm - xl ) ) 
          * ( ( xq - xr ) / ( xm - xr ) ) 
          + ( ( xq - xl ) / ( xm - xl ) ) 
          * (         1.0 / ( xm - xr ) );

      vrp = (         1.0 / ( xr - xl ) ) 
          * ( ( xq - xm ) / ( xr - xm ) ) 
          + ( ( xq - xl ) / ( xr - xl ) ) 
          * (         1.0 / ( xr - xm ) );

      amat[l+l*n] = amat[l+l*n] + wq * ( vlp * axq * vlp + vl * cxq * vl );
      amat[l+m*n] = amat[l+m*n] + wq * ( vlp * axq * vmp + vl * cxq * vm );
      amat[l+r*n] = amat[l+r*n] + wq * ( vlp * axq * vrp + vl * cxq * vr );
      b[l]   = b[l]   + wq * ( vl * fxq );

      amat[m+l*n] = amat[m+l*n] + wq * ( vmp * axq * vlp + vm * cxq * vl );
      amat[m+m*n] = amat[m+m*n] + wq * ( vmp * axq * vmp + vm * cxq * vm );
      amat[m+r*n] = amat[m+r*n] + wq * ( vmp * axq * vrp + vm * cxq * vr );
      b[m] =   b[m]   + wq * ( vm * fxq );

      amat[r+l*n] = amat[r+l*n] + wq * ( vrp * axq * vlp + vr * cxq * vl );
      amat[r+m*n] = amat[r+m*n] + wq * ( vrp * axq * vmp + vr * cxq * vm );
      amat[r+r*n] = amat[r+r*n] + wq * ( vrp * axq * vrp + vr * cxq * vr );
      b[r] =   b[r]   + wq * ( vr * fxq );
    }
  }
//
//  Equation 0 is the left boundary condition, U(0.0) = 0.0;
//
  i = 0;
  for ( j = 0; j < n; j++ )
  {
    amat[i+j*n] = 0.0;
  }
  amat[i+i*n] = 1.0;
  b[i] = 0.0;
//
//  Equation N-1 is the right boundary condition, U(1.0) = 0.0;
//
  i = n - 1;
  for ( j = 0; j < n; j++ )
  {
    amat[i+j*n] = 0.0;
  }
  amat[i+i*n] = 1.0;
  b[i] = 0.0;
//
//  Solve the linear system.
//
  u = r8mat_solve2 ( n, amat, b, &ierror );

  delete [] amat;
  delete [] b;

  return u;
# undef QUAD_NUM
}
//****************************************************************************80

double h1s_error_quadratic ( int n, double x[], double u[], 
  double exact_ux ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    H1S_ERROR_QUADRATIC estimates the seminorm error of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise quadratic elements used for the basis.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact derivative is known.
//
//    This function estimates the seminorm of the error:
//
//      SEMINORM = Integral ( A <= X <= B ) ( dU(X)/dx - EXACT_UX(X) )^2 dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT_UX ( X ), returns the value of the exact
//    derivative at the point X.
//
//    Output, double H1S_ERROR_QUADRATIC, the estimated seminorm of 
//    the error.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956,
     0.000000000000000000000000000000,
     0.774596669241483377035853079956 };
  int e;
  int e_num;
  double exq;
  double h1s;
  int i;
  int l;
  int m;
  int q;
  int quad_num = QUAD_NUM;
  int r;
  double uxq;
  double vlp;
  double vmp;
  double vrp;
  double weight[QUAD_NUM] = { 
    0.555555555555555555555555555556,
    0.888888888888888888888888888889,
    0.555555555555555555555555555556 };
  double wq;
  double xl;
  double xm;
  double xq;
  double xr;

  h1s = 0.0;
//
//  Integrate over element E.
//
  e_num = ( n - 1 ) / 2;

  for ( e = 0; e < e_num; e++ )
  {
    l = 2 * e;
    m = 2 * e + 1;
    r = 2 * e + 2;

    xl = x[l];
    xm = x[m];
    xr = x[r];

    for ( q = 0; q < quad_num; q++ )
    {

      xq = ( ( 1.0 - abscissa[q] ) * xl 
           + ( 1.0 + abscissa[q] ) * xr ) 
           /   2.0;

      wq = weight[q] * ( xr - xl ) / 2.0;

      vlp = (         1.0 / ( xl - xm ) ) 
          * ( ( xq - xr ) / ( xl - xr ) ) 
          + ( ( xq - xm ) / ( xl - xm ) ) 
          * (         1.0 / ( xl - xr ) );

      vmp = (         1.0 / ( xm - xl ) ) 
          * ( ( xq - xr ) / ( xm - xr ) ) 
          + ( ( xq - xl ) / ( xm - xl ) ) 
          * (         1.0 / ( xm - xr ) );

      vrp = (         1.0 / ( xr - xl ) ) 
          * ( ( xq - xm ) / ( xr - xm ) ) 
          + ( ( xq - xl ) / ( xr - xl ) ) 
          * (         1.0 / ( xr - xm ) );

      uxq = u[l] * vlp + u[m] * vmp + u[r] * vrp;

      exq = exact_ux ( xq );
 
      h1s = h1s + wq * pow ( uxq - exq, 2 );
    }
  }
  h1s = sqrt ( h1s );

  return h1s;
# undef QUAD_NUM
}
//****************************************************************************80

int *i4vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}
//****************************************************************************80

double l1_error ( int n, double x[], double u[], double exact ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    L1_ERROR estimates the L1 error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates the little l1 norm of the error:
//      L1_NORM = sum ( 1 <= I <= N ) abs ( U(i) - EXACT(X(i)) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT ( X ), returns the value of the exact
//    solution at the point X.
//
//    Output, double L1_ERROR, the estimated L2 norm of the error.
//
{
  int i;
  double e1;

  e1 = 0.0;

  for ( i = 0; i < n; i++ )
  {
    e1 = e1 + fabs ( u[i] - exact ( x[i] ) );
  }

  e1 = e1 / ( double ) ( n );

  return e1;
}
//****************************************************************************80

double l2_error_quadratic ( int n, double x[], double u[], 
  double exact ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    L2_ERROR_QUADRATIC estimates the L2 error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise quadratic elements used for the basis.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates the L2 norm of the error:
//
//      L2_NORM = Integral ( A <= X <= B ) ( U(X) - EXACT(X) )^2 dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT ( X ), returns the value of the exact
//    solution at the point X.
//
//    Output, double L2_ERROR_QUADRATIC, the estimated L2 norm of the error.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956,
     0.000000000000000000000000000000,
     0.774596669241483377035853079956 };
  int e;
  int e_num;
  double e2;
  double eq;
  int i;
  int l;
  int m;
  int q;
  int quad_num = QUAD_NUM;
  int r;
  double uq;
  double vl;
  double vm;
  double vr;
  double weight[QUAD_NUM] = { 
    0.555555555555555555555555555556,
    0.888888888888888888888888888889,
    0.555555555555555555555555555556 };
  double wq;
  double xl;
  double xm;
  double xq;
  double xr;

  e2 = 0.0;
//
//  Integrate over element E.
//
  e_num = ( n - 1 ) / 2;

  for ( e = 0; e < e_num; e++ )
  {
    l = 2 * e;
    m = 2 * e + 1;
    r = 2 * e + 2;

    xl = x[l];
    xm = x[m];
    xr = x[r];

    for ( q = 0; q < quad_num; q++ )
    {

      xq = ( ( 1.0 - abscissa[q] ) * xl 
           + ( 1.0 + abscissa[q] ) * xr ) 
           /   2.0;

      wq = weight[q] * ( xr - xl ) / 2.0;

      vl = ( ( xq - xm ) / ( xl - xm ) ) 
         * ( ( xq - xr ) / ( xl - xr ) );

      vm = ( ( xq - xl ) / ( xm - xl ) ) 
         * ( ( xq - xr ) / ( xm - xr ) );

      vr = ( ( xq - xl ) / ( xr - xl ) ) 
         * ( ( xq - xm ) / ( xr - xm ) );

      uq = u[l] * vl + u[m] * vm + u[r] * vr;
      eq = exact ( xq );

      e2 = e2 + wq * pow ( uq - eq, 2 );
    }
  }
  e2 = sqrt ( e2 );

  return e2;
# undef QUAD_NUM
}
//****************************************************************************80

double *r8mat_solve2 ( int n, double a[], double b[], int *ierror )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE2 computes the solution of an N by N linear system.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
//
//    The linear system may be represented as
//
//      A*X = B
//
//    If the linear system is singular, but consistent, then the routine will
//    still produce a solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of equations.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix to be inverted.
//    On output, A has been overwritten.
//
//    Input/output, double B[N].
//    On input, B is the right hand side of the system.
//    On output, B has been overwritten.
//
//    Output, int *IERROR.
//    0, no error detected.
//    1, consistent singularity.
//    2, inconsistent singularity.
//
//    Output, double R8MAT_SOLVE2[N], the solution of the linear system.
//
{
  double amax;
  int i;
  int imax;
  int j;
  int k;
  int *piv;
  double *x;

  *ierror = 0;

  piv = i4vec_zero_new ( n );
  x = r8vec_zero_new ( n );
//
//  Process the matrix.
//
  for ( k = 1; k <= n; k++ )
  {
//
//  In column K:
//    Seek the row IMAX with the properties that:
//      IMAX has not already been used as a pivot;
//      A(IMAX,K) is larger in magnitude than any other candidate.
//
    amax = 0.0;
    imax = 0;
    for ( i = 1; i <= n; i++ )
    {
      if ( piv[i-1] == 0 )
      {
        if ( amax < fabs ( a[i-1+(k-1)*n] ) )
        {
          imax = i;
          amax = fabs ( a[i-1+(k-1)*n] );
        }
      }
    }
//
//  If you found a pivot row IMAX, then,
//    eliminate the K-th entry in all rows that have not been used for pivoting.
//
    if ( imax != 0 )
    {
      piv[imax-1] = k;
      for ( j = k+1; j <= n; j++ )
      {
        a[imax-1+(j-1)*n] = a[imax-1+(j-1)*n] / a[imax-1+(k-1)*n];
      }
      b[imax-1] = b[imax-1] / a[imax-1+(k-1)*n];
      a[imax-1+(k-1)*n] = 1.0;

      for ( i = 1; i <= n; i++ )
      {
        if ( piv[i-1] == 0 )
        {
          for ( j = k+1; j <= n; j++ )
          {
            a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - a[i-1+(k-1)*n] * a[imax-1+(j-1)*n];
          }
          b[i-1] = b[i-1] - a[i-1+(k-1)*n] * b[imax-1];
          a[i-1+(k-1)*n] = 0.0;
        }
      }
    }
  }
//
//  Now, every row with nonzero IPIV begins with a 1, and
//  all other rows are all zero.  Begin solution.
//
  for ( j = n; 1 <= j; j-- )
  {
    imax = 0;
    for ( k = 1; k <= n; k++ )
    {
      if ( piv[k-1] == j )
      {
        imax = k;
      }
    }

    if ( imax == 0 )
    {
      x[j-1] = 0.0;

      if ( b[j-1] == 0.0 )
      {
        *ierror = 1;
        cout << "\n";
        cout << "R8MAT_SOLVE2 - Warning:\n";
        cout << "  Consistent singularity, equation = " << j << "\n";
      }
      else
      {
        *ierror = 2;
        cout << "\n";
        cout << "R8MAT_SOLVE2 - Warning:\n";
        cout << "  Inconsistent singularity, equation = " << j << "\n";
      }
    }
    else
    {
      x[j-1] = b[imax-1];

      for ( i = 1; i <= n; i++ )
      {
        if ( i != imax )
        {
          b[i-1] = b[i-1] - a[i-1+(j-1)*n] * x[j-1];
        }
      }
    }
  }

  delete [] piv;

  return x;
}
//****************************************************************************80

double *r8mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ZERO_NEW returns a new zeroed R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, double R8MAT_ZERO[M*N], the new zeroed matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

double *r8vec_even_new ( int n, double alo, double ahi )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN_NEW returns an R8VEC of values evenly spaced between ALO and AHI.
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
//    Output, double R8VEC_EVEN_NEW[N], N evenly spaced values.
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

double *r8vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO_NEW creates and zeroes an R8VEC.
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
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
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
