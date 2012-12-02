# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "fem1d_bvp_linear.hpp"

//****************************************************************************80

double compute_l2_error ( int n, double x[], double u[], 
  double exact ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    COMPUTE_L2_ERROR estimates the L2 error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise linear elements used for the basis.
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
//    17 February 2012
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
//    Output, double COMPUTE_L2_ERROR, the estimated L2 norm of the error.
//
{
  double abscissa[2];
  double eq;
  int i;
  double l2_error;
  int q;
  int quad_num;
  double ul;
  double ur;
  double uq;
  double weight[2];
  double wq;
  double xl;
  double xq;
  double xr;

  l2_error = 0.0;
//
//  Quadrature definitions.
//
  quad_num = 2;
  abscissa[0] = -0.577350269189625764509148780502;
  abscissa[1] = +0.577350269189625764509148780502;
  weight[0] = 1.0;
  weight[1] = 1.0;
//
//  Integrate over each interval.
//
  for ( i = 0; i < n - 1; i++ )
  {
    xl = x[i];
    xr = x[i+1];
    ul = u[i];
    ur = u[i+1];

    for ( q = 0; q < quad_num; q++ )
    {
      xq = ( ( 1.0 - abscissa[q] ) * xl   
           + ( 1.0 + abscissa[q] ) * xr ) 
           /   2.0;

      wq = weight[q] * ( xr - xl ) / 2.0;
//
//  Use the fact that U is a linear combination of piecewise linears.
//
      uq = ( ( xr - xq      ) * ul 
           + (      xq - xl ) * ur ) 
           / ( xr      - xl );

      eq = exact ( xq );

      l2_error = l2_error + wq * pow ( uq - eq, 2 );
    }
  }
  l2_error = sqrt ( l2_error );

  return l2_error;
}
//****************************************************************************80

double compute_seminorm_error ( int n, double x[], double u[], 
  double exact_ux ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    COMPUTE_SEMINORM_ERROR estimates the seminorm error of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise linear elements used for the basis.
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
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X(N), the mesh points.
//
//    Input, double U(N), the finite element coefficients.
//
//    Input, function EQ = EXACT_UX ( X ), returns the value of the exact
//    derivative at the point X.
//
//    Output, double COMPUTE_SEMINORM_ERROR, the estimated seminorm of 
//    the error.
//
{
  double abscissa[2];
  double exq;
  int i;
  int q;
  int quad_num;
  double seminorm_error;
  double ul;
  double ur;
  double uxq;
  double weight[2];
  double wq;
  double xl;
  double xq;
  double xr;

  seminorm_error = 0.0;
//
//  Quadrature definitions.
//
  quad_num = 2;
  abscissa[0] = -0.577350269189625764509148780502;
  abscissa[1] = +0.577350269189625764509148780502;
  weight[0] = 1.0;
  weight[1] = 1.0;
//
//  Integrate over each interval.
//
  for ( i = 0; i < n - 1; i++ )
  {
    xl = x[i];
    xr = x[i+1];
    ul = u[i];
    ur = u[i+1];

    for ( q = 0; q < quad_num; q++ )
    {
      xq = ( ( 1.0 - abscissa[q] ) * xl   
           + ( 1.0 + abscissa[q] ) * xr ) 
           /   2.0;

      wq = weight[q] * ( xr - xl ) / 2.0;
//
//  The piecewise linear derivative is a constant in the interval.
//
      uxq = ( ur - ul ) / ( xr - xl );

      exq = exact_ux ( xq );
 
      seminorm_error = seminorm_error + wq * pow ( uxq - exq, 2);
    }
  }
  seminorm_error = sqrt ( seminorm_error );

  return seminorm_error;
}
//****************************************************************************80

double *fem1d_bvp_linear ( int n, double a ( double x ), double c ( double x ), 
  double f ( double x ), double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR solves a two point boundary value problem.
//
//  Discussion:
//
//    The program uses the finite element method, with piecewise linear basis
//    functions to solve a boundary value problem in one dimension.
//
//    The problem is defined on the region 0 <= x <= 1.
//
//    The following differential equation is imposed between 0 and 1:
//
//      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
//
//    where a(x), c(x), and f(x) are given functions.
//
//    At the boundaries, the following conditions are applied:
//
//      u(0.0) = 0.0
//      u(1.0) = 0.0
//
//    A set of N equally spaced nodes is defined on this
//    interval, with 0 = X(1) < X(2) < ... < X(N) = 1.0.
//
//    At each node I, we associate a piecewise linear basis function V(I,X),
//    which is 0 at all nodes except node I.  This implies that V(I,X) is
//    everywhere 0 except that
//
//    for X(I-1) <= X <= X(I):
//
//      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) ) 
//
//    for X(I) <= X <= X(I+1):
//
//      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )
//
//    We now assume that the solution U(X) can be written as a linear
//    sum of these basis functions:
//
//      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)
//
//    where U(X) on the left is the function of X, but on the right,
//    is meant to indicate the coefficients of the basis functions.
//
//    To determine the coefficient U(J), we multiply the original
//    differential equation by the basis function V(J,X), and use
//    integration by parts, to arrive at the I-th finite element equation:
//
//        Integral A(X) * U'(X) * V'(I,X) + C(X) * U(X) * V(I,X) dx 
//      = Integral F(X) * V(I,X) dx
//
//    We note that the functions U(X) and U'(X) can be replaced by
//    the finite element form involving the linear sum of basis functions,
//    but we also note that the resulting integrand will only be nonzero
//    for terms where J = I - 1, I, or I + 1.
//
//    By writing this equation for basis functions I = 2 through N - 1,
//    and using the boundary conditions, we have N linear equations
//    for the N unknown coefficients U(1) through U(N), which can
//    be easily solved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
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
//    Output, double FEM1D_BVP_LINEAR[N], the finite element coefficients, 
//    which are also the value of the computed solution at the mesh points.
//
{
# define QUAD_NUM 2

  double abscissa[QUAD_NUM] = {
    -0.577350269189625764509148780502,
    +0.577350269189625764509148780502 };
  double al;
  double am;
  double ar;
  double *amat;
  double axq;
  double *b;
  double bm;
  double cxq;
  double fxq;
  double h;
  int i;
  int ierror;
  int q;
  int quad_num = QUAD_NUM;
  double *u;
  double weight[QUAD_NUM] = { 1.0, 1.0 };
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
//  Equation 1 is the left boundary condition, U(0.0) = 0.0;
//
  amat[0+0*n] = 1.0;
  b[0] = 0.0;
//
//  Equation I involves the basis function at node I.
//  This basis function is nonzero from X(I-1) to X(I+1).
//  Equation I looks like this:
//
//    Integral A(X) U'(X) V'(I,X) 
//           + C(X) * U(X) V(I,X) dx 
//  = Integral F(X) V(I,X) dx
//
//  Then, we realize that U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X), 
//  (U(X) means the function; U(J) is the coefficient of V(J,X) ).
//
//  The only V functions that are nonzero when V(I,X) is nonzero are
//  V(I-1,X) and V(I+1,X). 
//
//  Let's use the shorthand 
//
//    VL(X) = V(I-1,X)
//    VM(X) = V(I,X)
//    VR(X) = V(I+1,X)
//
//  So our equation becomes
//
//    Integral A(X) [ VL'(X) U(I-1) + VM'(X) U(I) + VR'(X) U(I+1) ] * VM'(X)
//           + C(X) [ VL(X)  U(I-1) + VM(X)  U(I) + VR(X)  U(I+1) ] * VM(X) dx
//  = Integral F(X) VM(X) dx.
//
//  
//
//  This is actually a set of N-2 linear equations for the N coefficients U.
//
//  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1), 
//  and so on.
//
  for ( i = 1; i < n - 1; i++ )
  {
//
//  Get the left, right and middle coordinates.
//
    xl = x[i-1];
    xm = x[i];
    xr = x[i+1];
//
//  Make temporary variables for A(I,I-1), A(I,I), A(I,I+1) and B(I).
//
    al = 0.0;
    am = 0.0;
    ar = 0.0;
    bm = 0.0;
//
//  We approximate the integrals by using a weighted sum of
//  the integrand values at quadrature points.
//
    for ( q = 0; q < quad_num; q++ )
    {
//
//  Integrate over the LEFT interval, between XL and XM, where:
//
//  VL(X) = ( XM - X       ) / ( XM - XL )
//  VM(X) = (      X  - XL ) / ( XM - XL )
//  VR(X) = 0
//
//  VL'(X) =             - 1 / ( XM - XL )
//  VM'(X) =             + 1 / ( XM - XL ) 
//  VR'(X) = 0
//
      xq = ( ( 1.0 - abscissa[q] ) * xl 
           + ( 1.0 + abscissa[q] ) * xm ) 
           /   2.0;

      wq = weight[q] * ( xm - xl ) / 2.0;

      vl =  ( xm - xq ) / ( xm - xl );
      vlp =      - 1.0  / ( xm - xl );

      vm =  ( xq - xl ) / ( xm - xl );
      vmp =      + 1.0  / ( xm - xl );

      vr =  0.0;
      vrp = 0.0;

      axq = a ( xq );
      cxq = c ( xq );
      fxq = f ( xq );

      al = al + wq * ( axq * vlp * vmp + cxq * vl * vm );
      am = am + wq * ( axq * vmp * vmp + cxq * vm * vm );
      ar = ar + wq * ( axq * vrp * vmp + cxq * vr * vm );
      bm = bm + wq * ( fxq * vm );
//
//  Integrate over the RIGHT interval, between XM and XR, where:
//
//  VL(X) = 0
//  VM(X) = ( XR - X       ) / ( XR - XM )
//  VR(X) = (      X  - XM ) / ( XR - XM )
//
//  VL'(X) = 0
//  VM'(X) =             - 1 / ( XR - XM )
//  VR'(X) =             + 1 / ( XR - XM ) 
//
      xq = ( ( 1.0 - abscissa[q] ) * xm 
           + ( 1.0 + abscissa[q] ) * xr ) 
           /   2.0;

      wq = weight[q] * ( xr - xm ) / 2.0;

      vl = 0.0;
      vlp = 0.0;

      vm = ( xr - xq ) / ( xr - xm );
      vmp =     - 1.0  / ( xr - xm );

      vr = ( xq - xm ) / ( xr - xm );
      vrp =      1.0   / ( xr - xm );

      axq = a ( xq );
      cxq = c ( xq );
      fxq = f ( xq );

      al = al + wq * ( axq * vlp * vmp + cxq * vl * vm );
      am = am + wq * ( axq * vmp * vmp + cxq * vm * vm );
      ar = ar + wq * ( axq * vrp * vmp + cxq * vr * vm );
      bm = bm + wq * ( fxq * vm );
    }
    amat[i+(i-1)*n] = al;
    amat[i+ i   *n] = am;
    amat[i+(i+1)*n] = ar;

    b[i] = bm;
  }
//
//  Equation N is the right boundary condition, U(1.0) = 0.0;
//
  amat[n-1+(n-1)*n] = 1.0;
  b[n-1] = 0.0;
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
        if ( amax < r8_abs ( a[i-1+(k-1)*n] ) )
        {
          imax = i;
          amax = r8_abs ( a[i-1+(k-1)*n] );
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
